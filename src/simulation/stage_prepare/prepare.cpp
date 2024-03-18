#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include <h5.hpp>

#include "../common/h5_traits.hpp"
#include "prepare.hpp"


// Numerical code used in the "/metadata/particle_types" dataset.
enum class particle_type
{
    unknown    = 0,
    a          = 1,
    b          = 2,
    u          = 3,
    centromere = 4,
    active_nor = 5,
    silent_nor = 6,
    nucleolus  = 7,
};


// Numerical code used in the "/metadata/particle_types:metaphase" dataset.
enum class metaphase_particle_type
{
    unknown    = 0,
    arm        = 1,
    centromere = 2,
};


// Structure holding the parameters of a simulated particle.
struct particle_data
{
    double        a_factor = 0;
    double        b_factor = 0;
    particle_type type     = particle_type::unknown;
};


// Structure holding the region in a one-dimensional array of particles that
// constitutes a single chromosome in a simulation.
struct chain_assignment
{
    std::string name;
    std::size_t start = 0;
    std::size_t end   = 0;
};


// Trajectory file preparation pipeline.
class preparation_pipeline
{
public:
    preparation_pipeline(
        h5::file& store,
        simulation_config const& config,
        chain_definitions const& chains,
        std::uint32_t master_seed
    );
    void run();

private:
    void define_chains();
    void define_nucleolar_particles();
    void define_metaphase_chains();

    void write_inputs();
    void write_interphase_particles();
    void write_metaphase_particles();
    void write_seeds();
    void write_stages();

private:
    h5::file&                               _store;
    simulation_config const                 _config;
    chain_definitions const                 _chains;
    std::uint32_t const                     _master_seed;
    std::vector<particle_data>              _particles;
    std::vector<chain_assignment>           _chain_assignments;
    std::vector<std::size_t>                _nor_indices;
    std::vector<std::array<std::size_t, 2>> _nucleolar_bonds;
    std::vector<metaphase_particle_type>    _metaphase_particle_types;
    std::vector<chain_assignment>           _metaphase_chain_assignments;
    std::vector<std::size_t>                _metaphase_kinetochore_indices;
};


preparation_pipeline::preparation_pipeline(
    h5::file& store,
    simulation_config const& config,
    chain_definitions const& chains,
    std::uint32_t master_seed
)
: _store(store)
, _config(config)
, _chains(chains)
, _master_seed(master_seed)
{
}


void
preparation_pipeline::run()
{
    define_chains();
    define_nucleolar_particles();
    define_metaphase_chains();

    write_inputs();
    write_interphase_particles();
    write_metaphase_particles();
    write_seeds();
    write_stages();
}


namespace
{
    // Checks if comma-delimited string `tags` contains a field with given
    // string `tag`.
    bool
    check_tag(std::string const& tags, std::string const& tag)
    {
        constexpr char delimiter = ',';

        for (auto it = tags.begin(); it != tags.end(); ) {
            auto const break_point = std::find(it, tags.end(), delimiter);
            if (std::equal(it, break_point, tag.begin(), tag.end())) {
                return true;
            }
            it = break_point;
            if (it != tags.end()) {
                ++it;
            }
        }
        return false;
    }
}


// Defines chromosomal chains in interphase simulation.
void
preparation_pipeline::define_chains()
{
    // Input:  (none)
    // Output: _particles, _nor_indices, _chain_assignments

    std::vector<std::pair<std::string, particle_type>> const tag_type_map = {
        {"anor", particle_type::active_nor},
        {"bnor", particle_type::silent_nor},
        {"cen", particle_type::centromere},
        {"A", particle_type::a},
        {"B", particle_type::b},
        {"u", particle_type::u},
    };

    for (auto const& chain : _chains.chains) {
        auto const chain_start = _particles.size();
        auto const chain_end = chain_start + chain.beads.size();

        for (auto const& bead : chain.beads) {
            auto const bead_index = _particles.size();
            auto bead_type = particle_type::unknown;

            for (auto const& [tag, type] : tag_type_map) {
                if (check_tag(bead.tags, tag)) {
                    bead_type = type;
                    break;
                }
            }

            if (bead_type == particle_type::active_nor) {
                _nor_indices.push_back(bead_index);
            }

            _particles.push_back({
                .a_factor = bead.a_factor,
                .b_factor = bead.b_factor,
                .type     = bead_type,
            });
        }

        _chain_assignments.push_back({
            .name  = chain.name,
            .start = chain_start,
            .end   = chain_end,
        });
    }
}


// Defines nucleolar particles in interphase simulation.
void
preparation_pipeline::define_nucleolar_particles()
{
    // Input:  _particles, _nor_indices
    // Output: _particles, _nucleolar_bonds

    for (auto const& nor_index : _nor_indices) {
        for (std::size_t i = 0; i < _config.interphase.nucleolus_bead_count; i++) {
            auto const nucleolus_index = _particles.size();
            _particles.push_back({
                .a_factor = _config.interphase.nucleolus_ab_factor.a,
                .b_factor = _config.interphase.nucleolus_ab_factor.b,
                .type     = particle_type::nucleolus,
            });
            _nucleolar_bonds.push_back({ nor_index, nucleolus_index });
        }
    }
}


namespace
{
    // Determines the half-open window [start, end) of given size around
    // `center`, clipping the lower and upper ends to given bounds.
    std::array<std::size_t, 2>
    define_clipped_window(
        std::size_t center,
        std::size_t window,
        std::size_t lower_bound,
        std::size_t upper_bound
    )
    {
        auto const lower_window = window / 2;
        auto const upper_window = window - lower_window;

        auto const lower =
            (lower_bound + lower_window > center)
                ? lower_bound
                : center - lower_window;

        auto const upper =
            (center + upper_window > upper_bound)
                ? upper_bound
                : center + upper_window;

        return { lower, upper };
    }
}


// Defines coarse-grained chromosomal chains in metaphase simulation.
void
preparation_pipeline::define_metaphase_chains()
{
    // Input:  _particles, _chain_assignments
    // Output: _metaphase_particle_types, _metaphase_chain_assignments, _metaphase_kinetochore_indices

    auto const coarse_graining = _config.anatelophase.coarse_graining;
    auto const dragged_beads = _config.anatelophase.dragged_beads;

    // Identify the centromeric region [start, end) of each chain.
    std::vector<std::array<std::size_t, 2>> centromere_ranges;

    for (auto const& assign : _chain_assignments) {
        std::size_t start = assign.start;
        std::size_t end = assign.end;
        bool seen = false;

        for (std::size_t i = assign.start; i < assign.end; i++) {
            if (_particles[i].type == particle_type::centromere) {
                if (!seen) {
                    start = i;
                    seen = true;
                }
                end = i + 1;
            }
        }

        if (!seen) {
            // No annotation. Proceed anyway, marking the whole chromosome
            // as being centromeric.
            std::clog << "No centromere found on " << assign.name << "\n";
        }

        centromere_ranges.push_back({ start, end });
    }

    // Define metaphase chains by coarse-graining the interphase ones.

    for (std::size_t chain_index = 0; chain_index < _chain_assignments.size(); chain_index++) {
        auto const& assign = _chain_assignments[chain_index];
        auto const length = assign.end - assign.start;

        auto const coarse_length = length / coarse_graining;
        auto const coarse_chain_start = _metaphase_particle_types.size();
        auto const coarse_chain_end = coarse_chain_start + coarse_length;

        // Map centromere to `dragged_beads` coarse-grained beads around the
        // midpoint of the centromeric region of the original chain.
        auto const [cen_start, cen_end] = centromere_ranges[chain_index];
        auto const centromere_midpoint = (cen_start + cen_end) / 2;

        auto const kinetochore_center = (centromere_midpoint - assign.start) / coarse_graining;
        auto const [kinetochore_start, kinetochore_end] = define_clipped_window(
            kinetochore_center, dragged_beads, 0, coarse_length
        );

        for (std::size_t bin = 0; bin < coarse_length; bin++) {
            auto const bead_index = _metaphase_particle_types.size();

            auto type = metaphase_particle_type::arm;
            if (bin >= kinetochore_start && bin < kinetochore_end) {
                type = metaphase_particle_type::centromere;
                _metaphase_kinetochore_indices.push_back(bead_index);
            }

            _metaphase_particle_types.push_back(type);
        }

        _metaphase_chain_assignments.push_back({
            .name  = assign.name,
            .start = coarse_chain_start,
            .end   = coarse_chain_end,
        });
    }
}


// Writes input parameters under "/metadata" hierarchy.
void
preparation_pipeline::write_inputs()
{
    // "/metadata/config" is the serialization of the actual parameters used in
    // the simulation.
    _store.dataset<h5::u32>("/metadata/master_seed").write(_master_seed);
    _store.dataset<h5::str>("/metadata/config").write(format_simulation_config(_config));
    _store.dataset<h5::str>("/metadata/config_source").write(_config.source);
    _store.dataset<h5::str>("/metadata/chains_source").write(_chains.source);
}


// Writes paramerters for the interphase particles.
void
preparation_pipeline::write_interphase_particles()
{
    // Store metadata datasets that define interphase particles and topology:
    // - particle_types  (*) enum
    // - ab_factors  (*) float
    // - chain_names  (*) str
    // - chain_ranges  (*, 2) int
    // - nucleolar_bonds  (*, 2) int

    h5::enums<int> const particle_types_enum = {
        {"unknown", int(particle_type::unknown)},
        {"a", int(particle_type::a)},
        {"b", int(particle_type::b)},
        {"u", int(particle_type::u)},
        {"centromere", int(particle_type::centromere)},
        {"active_nor", int(particle_type::active_nor)},
        {"silent_nor", int(particle_type::silent_nor)},
        {"nucleolus", int(particle_type::nucleolus)},
    };

    std::vector<int> particle_types;
    std::vector<std::array<double, 2>> ab_factors;
    for (auto const& particle : _particles) {
        particle_types.push_back(int(particle.type));
        ab_factors.push_back({ particle.a_factor, particle.b_factor });
    }

    std::vector<std::string> chain_names;
    std::vector<std::array<std::size_t, 2>> chain_ranges;
    for (auto const& assign : _chain_assignments) {
        chain_names.push_back(assign.name);
        chain_ranges.push_back({ assign.start, assign.end });
    }

    _store.dataset<h5::i32, 1>("/metadata/particle_types", particle_types_enum).write(particle_types);
    _store.dataset<h5::f32, 2>("/metadata/ab_factors").write(ab_factors);
    _store.dataset<h5::str, 1>("/metadata/chain_names").write(chain_names);
    _store.dataset<h5::i32, 2>("/metadata/chain_ranges").write(chain_ranges);
    _store.dataset<h5::i32, 2>("/metadata/nucleolar_bonds").write(_nucleolar_bonds);
}


// Writes parameters for the metaphase particles.
void
preparation_pipeline::write_metaphase_particles()
{
    // Store metadata datasets that define interphase particles and topology:
    // - particle_types:metaphase  (*) enum
    // - chain_ranges:metaphase  (*, 2) int
    // - kinetochore_beads:metaphase  (*) int

    h5::enums<int> const particle_types_enum = {
        {"unknown", int(metaphase_particle_type::unknown)},
        {"arm", int(metaphase_particle_type::arm)},
        {"centromere", int(metaphase_particle_type::centromere)},
    };

    std::vector<int> particle_types;
    for (auto const& type : _metaphase_particle_types) {
        particle_types.push_back(int(type));
    }

    std::vector<std::array<std::size_t, 2>> chain_ranges;
    for (auto const& assign : _metaphase_chain_assignments) {
        chain_ranges.push_back({ assign.start, assign.end });
    }

    _store.dataset<h5::i32, 1>("/metadata/particle_types:metaphase", particle_types_enum).write(particle_types);
    _store.dataset<h5::i32, 2>("/metadata/chain_ranges:metaphase").write(chain_ranges);
    _store.dataset<h5::i32, 1>("/metadata/kinetochore_beads:metaphase").write(_metaphase_kinetochore_indices);
}


// Derives and stores the random seed used in each simulation stage.
void
preparation_pipeline::write_seeds()
{
    std::seed_seq seed_seq = { _master_seed };
    std::array<std::uint32_t, 2> seed_values;
    seed_seq.generate(seed_values.begin(), seed_values.end());

    _store.dataset<h5::u32>("/metadata/seed:metaphase").write(seed_values[0]);
    _store.dataset<h5::u32>("/metadata/seed:interphase").write(seed_values[1]);
}


namespace
{
    void
    link_dataset(h5::file& store, std::string const& dest, std::string const& name)
    {
        // Allow intermediate groups to be automatically created.
        h5::unique_hid<H5Pclose> link_props = H5Pcreate(H5P_LINK_CREATE);
        if (link_props < 0) {
            throw std::runtime_error("failed to create link props");
        }
        if (H5Pset_create_intermediate_group(link_props, 1) < 0) {
            throw std::runtime_error("failed to configure link props");
        }

        auto const ret = H5Lcreate_soft(
            dest.c_str(), store.handle(), name.c_str(), link_props, H5P_DEFAULT
        );
        if (ret < 0) {
            throw std::runtime_error("failed to create link " + name + " -> " + dest);
        }
    }
}


// Initializes simulation stages.
void
preparation_pipeline::write_stages()
{
    // Set up metadata datasets in each simulation stage (/stages/*/metadata).
    // These are used in analysis scripts.

    auto do_link = [this](
        std::string const& stage,
        std::string const& key,
        std::string const& suffix = ""
    ) {
        link_dataset(
            _store,
            "/metadata/" + key + suffix,
            "/stages/" + stage + "/metadata/" + key
        );
    };

    std::initializer_list<std::string> const metaphase_stages = {"anaphase", "telophase"};
    for (auto const& stage : metaphase_stages) {
        do_link(stage, "particle_types", ":metaphase");
        do_link(stage, "chain_names");
        do_link(stage, "chain_ranges", ":metaphase");
    }
    do_link("anaphase", "kinetochore_beads", ":metaphase");

    std::initializer_list<std::string> const standard_stages = {"relaxation", "interphase"};
    for (auto const& stage : standard_stages) {
        do_link(stage, "particle_types");
        do_link(stage, "chain_names");
        do_link(stage, "chain_ranges");
        do_link(stage, "ab_factor");
        do_link(stage, "nucleolar_bonds");
    }
}


// ---------------------------------------------------------------------------

void
prepare_simulation_store(
    h5::file& store,
    simulation_config const& config,
    chain_definitions const& chains,
    std::uint32_t master_seed
)
{
    preparation_pipeline(store, config, chains, master_seed).run();
}
