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
    interphase_context,

    time,
    wall_semiaxes,
    core_scale,
    bond_scale,
    mean_energy,
    wall_energy
)


JSONCONS_ALL_MEMBER_TRAITS(
    prometaphase_context::microtubule_record,

    length,
    oscillation_phase
)


JSONCONS_ALL_MEMBER_TRAITS(
    prometaphase_context,

    time,
    microtubules
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


anatelophase_design
simulation_store::load_anatelophase_design()
{
    // Anaphase and telophase share the same design.
    std::string const stage = "anaphase";
    return {
        .seed   = load_seed(stage),
        .chains = load_chains(stage),
    };
}


interphase_design
simulation_store::load_interphase_design()
{
    std::string const stage = "interphase";
    return {
        .seed            = load_seed(stage),
        .particles       = load_interphase_particles(stage),
        .chains          = load_chains(stage),
        .nucleolar_bonds = load_nucleolar_bonds(stage),
    };
}


prometaphase_design
simulation_store::load_prometaphase_design()
{
    std::string const stage = "prometaphase";
    return {
        .seed              = load_seed(stage),
        .chains            = load_chains(stage),
        .sister_chromatids = load_sister_chromatids(stage),
        .pole_positions    = load_pole_positions(stage),
    };
}


void
simulation_store::save_interphase_context(md::step step, interphase_context const& context)
{
    std::string context_json;
    jsoncons::encode_json(context, context_json);
    _store.dataset<h5::str>(locate_data(step, "context")).write(context_json);
}


interphase_context
simulation_store::load_interphase_context(md::step step)
{
    std::string context_json;
    _store.dataset<h5::str>(locate_data(step, "context")).read(context_json);
    return jsoncons::decode_json<interphase_context>(context_json);
}


void
simulation_store::save_prometaphase_context(md::step step, prometaphase_context const& context)
{
    std::string context_json;
    jsoncons::encode_json(context, context_json);
    _store.dataset<h5::str>(locate_data(step, "context")).write(context_json);
}


prometaphase_context
simulation_store::load_prometaphase_context(md::step step)
{
    std::string context_json;
    _store.dataset<h5::str>(locate_data(step, "context")).read(context_json);
    return jsoncons::decode_json<prometaphase_context>(context_json);
}


bool
simulation_store::check_positions(md::step step)
{
    return bool(_store.dataset<h5::f32, 2>(locate_data(step, "positions")));
}


void
simulation_store::clear_frames()
{
    if (auto dataset = _store.dataset<h5::str, 1>(locate_data(".steps"))) {
        std::vector<std::string> empty;
        dataset.write(empty);
    }
}


void
simulation_store::append_frame(md::step step)
{
    std::vector<std::string> frame_index;

    auto dataset = _store.dataset<h5::str, 1>(locate_data(".steps"));
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
    // log2(10^5) ~ 16 bits.
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

    _store.dataset<h5::f32, 2>(locate_data(step, "positions")).write(
        positions_array,
        {.compression = compression_level, .scaleoffset = {}}
    );
}


std::vector<md::step>
simulation_store::load_steps()
{
    std::vector<md::step> steps;

    // FIXME: Why strings?
    if (auto dataset = _store.dataset<h5::str, 1>(locate_data(".steps"))) {
        std::vector<std::string> step_keys;
        dataset.read_fit(step_keys);
        for (auto const& step_key : step_keys) {
            steps.push_back(std::stol(step_key));
        }
    }

    return steps;
}


std::vector<md::point>
simulation_store::load_positions(md::step step)
{
    std::vector<std::array<md::scalar, 3>> coords_array;
    _store.dataset<h5::f32, 2>(locate_data(step, "positions")).read_fit(coords_array);

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

    // Integral values compresses pretty good with the adaptive (zero)
    // scaleoffset filter.
    _store.dataset<h5::i32, 2>(locate_data(step, "contacts")).write(
        contacts, {.compression = 4, .scaleoffset = 0}
    );
}


std::uint64_t
simulation_store::load_seed(std::string const& stage)
{
    std::uint64_t value;
    _store.dataset<h5::u64>("/stages/" + stage + "/metadata/seed").read(value);
    return value;
}


std::vector<chain_range>
simulation_store::load_chains(std::string const& stage)
{
    std::vector<std::string> chain_names;
    std::vector<std::array<md::index, 2>> chain_ranges;
    _store.dataset<h5::str, 1>(locate_metadata_in(stage, "chain_names")).read_fit(chain_names);
    _store.dataset<h5::u32, 2>(locate_metadata_in(stage, "chain_ranges")).read_fit(chain_ranges);

    std::vector<chain_range> chains;
    for (std::size_t i = 0; i < chain_ranges.size(); i++) {
        chains.push_back({
            .name  = chain_names[i],
            .start = chain_ranges[i][0],
            .end   = chain_ranges[i][1],
        });
    }

    auto kinetochore_dataset = _store.dataset<h5::u32, 1>(locate_metadata_in(stage, "kinetochore_beads"));
    if (kinetochore_dataset) {
        std::vector<md::index> kinetochore_beads;
        kinetochore_dataset.read_fit(kinetochore_beads);

        if (kinetochore_beads.size() != chains.size()) {
            throw simulation_store_error("chains and kinetochore_beads datasets mismatch");
        }

        for (md::index i = 0; i < chains.size(); i++) {
            chains[i].kinetochore = kinetochore_beads[i];
        }
    }

    return chains;
}


std::vector<nucleolar_bond>
simulation_store::load_nucleolar_bonds(std::string const& stage)
{
    std::vector<std::array<md::index, 2>> index_pairs;
    _store.dataset<h5::u32, 2>(locate_metadata_in(stage, "nucleolar_bonds")).read_fit(index_pairs);

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
simulation_store::load_interphase_particles(std::string const& stage)
{
    std::vector<std::array<md::scalar, 2>> ab_factors;
    _store.dataset<h5::f32, 2>(locate_metadata_in(stage, "ab_factors")).read_fit(ab_factors);

    std::vector<particle_data> particles;
    for (auto const& ab : ab_factors) {
        particles.push_back({
            .a_factor = ab[0],
            .b_factor = ab[1],
        });
    }

    return particles;
}


std::vector<std::array<std::size_t, 2>>
simulation_store::load_sister_chromatids(std::string const& stage)
{
    std::vector<std::array<std::size_t, 2>> sister_chromatids;
    _store.dataset<h5::u32, 2>(locate_metadata_in(stage, "sister_chromatids")).read_fit(sister_chromatids);
    return sister_chromatids;
}


std::array<md::point, 2>
simulation_store::load_pole_positions(std::string const& stage)
{
    std::vector<std::array<md::scalar, 3>> pole_positions;
    _store.dataset<h5::f32, 2>(locate_metadata_in(stage, "pole_positions")).read_fit(pole_positions);
    if (pole_positions.size() != 2) {
        throw simulation_store_error("unexpected pole_positions shape");
    }
    auto const coords_to_point = [](std::array<md::scalar, 3> const& coords) -> md::point {
        return { coords[0], coords[1], coords[2] };
    };
    return {
        coords_to_point(pole_positions[0]),
        coords_to_point(pole_positions[1]),
    };
}


std::string
simulation_store::locate_metadata_in(std::string const& stage, std::string const& key) const
{
    return "/stages/" + stage + "/metadata/" + key;
}


std::string
simulation_store::locate_data(std::string const& key) const
{
    return "/stages/" + _stage + "/" + key;
}


std::string
simulation_store::locate_data(md::step step, std::string const& key) const
{
    return locate_data(std::to_string(step) + "/" + key);
}
