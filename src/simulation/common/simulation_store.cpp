#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include <jsoncons/json.hpp>

#include <h5.hpp>
#include <md.hpp>

#include "h5_traits.hpp"
#include "json_traits.hpp"
#include "particle_data.hpp"
#include "simulation_store.hpp"



namespace
{
    inline md::scalar
    quantize(md::scalar value, int bits)
    {
        // Similar to HDF5's builtin scaleoffset filter, but with binary scaling
        // rather than decimal. This fills the lower mantissa bits with zeros,
        // making the value more compressible.
        int exp;
        auto const mant = std::frexp(value, &exp);
        auto const scaled_mant = std::nearbyint(std::ldexp(mant, bits));
        return std::ldexp(scaled_mant, exp - bits);
    }
}


JSONCONS_ALL_MEMBER_TRAITS(
    simulation_context,

    time,
    wall_semiaxes,
    core_scale,
    bond_scale,
    mean_energy,
    wall_energy
)


simulation_store::simulation_store(std::string const& filename)
: _store(filename, "r+")
{
}


void
simulation_store::set_stage(std::string const& name)
{
    _stage = name;
}


simulation_config
simulation_store::load_config()
{
    std::string config_text;
    _store.dataset<h5::str>("/metadata/config").read(config_text);
    return parse_simulation_config(config_text);
}


metaphase_design
simulation_store::load_metaphase_design()
{
    return {
        .seed              = load_seed(":metaphase"),
        .chains            = load_chains(":metaphase"),
        .kinetochore_beads = load_kinetochore_beads(),
    };
}


interphase_design
simulation_store::load_interphase_design()
{
    return {
        .seed            = load_seed(":interphase"),
        .particles       = load_interphase_particles(),
        .chains          = load_chains(),
        .nucleolar_bonds = load_nucleolar_bonds(),
    };
}


void
simulation_store::save_context(md::step step, simulation_context const& context)
{
    std::string context_json;
    jsoncons::encode_json(context, context_json);
    _store.dataset<h5::str>(make_stage_path(step, "context")).write(context_json);
}


simulation_context
simulation_store::load_context(md::step step)
{
    std::string context_json;
    _store.dataset<h5::str>(make_stage_path(step, "context")).read(context_json);
    return jsoncons::decode_json<simulation_context>(context_json);
}


void
simulation_store::append_frame(md::step step)
{
    std::vector<std::string> frame_index;

    auto dataset = _store.dataset<h5::str, 1>(make_stage_path(".steps"));
    if (dataset) {
        dataset.read_fit(frame_index);
    }
    frame_index.push_back(std::to_string(step));

    dataset.write(frame_index);
}


void
simulation_store::save_positions(
    md::step step, md::array_view<md::point const> positions
)
{
    // Quantize coordinate values for better compression. Five significant
    // digits ought to be sufficient for genome-wide simulation, so use
    // 16 bits.
    constexpr int fraction_bits = 16;
    constexpr int compression_level = 6;

    std::vector<std::array<md::scalar, 3>> positions_array;
    for (auto const& position : positions) {
        positions_array.push_back({
            quantize(position.x, fraction_bits),
            quantize(position.y, fraction_bits),
            quantize(position.z, fraction_bits),
        });
    }

    _store.dataset<h5::f32, 2>(make_stage_path(step, "positions")).write(
        positions_array,
        {.compression = compression_level, .scaleoffset = {}}
    );
}


std::vector<md::point>
simulation_store::load_positions(md::step step)
{
    std::vector<std::array<md::scalar, 3>> coords_array;
    _store.dataset<h5::f32, 2>(make_stage_path(step, "positions")).read_fit(coords_array);

    std::vector<md::point> positions;
    positions.reserve(coords_array.size());
    for (auto const& coords : coords_array) {
        positions.push_back({coords[0], coords[1], coords[2]});
    }

    return positions;
}


void
simulation_store::save_contacts(
    md::step step, std::vector<std::array<std::uint32_t, 3>> const& contacts
)
{
    if (contacts.empty()) {
        return;
    }

    _store.dataset<h5::i32, 2>(make_stage_path(step, "contacts")).write(
        contacts, {.compression = 4, .scaleoffset = 0}
    );
}


std::uint64_t
simulation_store::load_seed(std::string suffix)
{
    std::uint64_t value;
    _store.dataset<h5::u64>("/metadata/seed" + suffix).read(value);
    return value;
}


std::vector<chain_range>
simulation_store::load_chains(std::string suffix)
{
    std::vector<std::string> chain_names;
    std::vector<std::array<md::index, 2>> chain_ranges;
    _store.dataset<h5::str, 1>("/metadata/chain_names").read_fit(chain_names);
    _store.dataset<h5::u32, 2>("/metadata/chain_ranges" + suffix).read_fit(chain_ranges);

    std::vector<chain_range> chains;
    for (std::size_t i = 0; i < chain_ranges.size(); i++) {
        chains.push_back({
            .name  = chain_names[i],
            .start = chain_ranges[i][0],
            .end   = chain_ranges[i][1],
        });
    }

    return chains;
}


std::vector<md::index>
simulation_store::load_kinetochore_beads()
{
    std::vector<md::index> indices;
    _store.dataset<h5::u32, 1>("/metadata/kinetochore_beads:metaphase").read_fit(indices);
    return indices;
}


std::vector<nucleolar_bond>
simulation_store::load_nucleolar_bonds()
{
    std::vector<std::array<md::index, 2>> index_pairs;
    _store.dataset<h5::u32, 2>("/metadata/nucleolar_bonds").read_fit(index_pairs);

    std::vector<nucleolar_bond> bonds;
    for (auto const& pair : index_pairs) {
        bonds.push_back({
            .nor_index = pair[0],
            .nuc_index = pair[1],
        });
    }

    return bonds;
}


std::vector<particle_data>
simulation_store::load_interphase_particles()
{
    std::vector<std::array<md::scalar, 2>> ab_factors;
    _store.dataset<h5::f32, 2>("/metadata/ab_factors").read_fit(ab_factors);

    std::vector<particle_data> particles;
    for (auto const& ab : ab_factors) {
        particles.push_back({
            .a_factor = ab[0],
            .b_factor = ab[1],
        });
    }

    return particles;
}


std::string
simulation_store::make_stage_path(std::string sub_path) const
{
    return "/stages/" + _stage + "/" + sub_path;
}


std::string
simulation_store::make_stage_path(md::step step, std::string sub_path) const
{
    return make_stage_path(std::to_string(step) + "/" + sub_path);
}
