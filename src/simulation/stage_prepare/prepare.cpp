#include <algorithm>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

#include <h5.hpp>

#include "../common/h5_traits.hpp"
#include "h5_misc.hpp"
#include "prepare.hpp"


// Numerical code used in the "/metadata/particle_types" dataset.
enum class interphase_particle_type
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


// Numerical code used in the "/metadata/particle_types:mitotic_phase" dataset.
enum class mitotic_particle_type
{
    unknown     = 0,
    arm         = 1,
    kinetochore = 2,
};


// Structure holding the parameters of a simulated particle.
struct interphase_particle_data
{
    double                   a_factor = 0;
    double                   b_factor = 0;
    interphase_particle_type type     = interphase_particle_type::unknown;
};


// Structure holding the region in a one-dimensional array of particles that
// constitutes a single chromosome in a simulation.
struct chain_assignment
{
    std::string                name;
    std::size_t                start = 0;
    std::size_t                end   = 0;
    std::optional<std::size_t> kinetochore;
};


// Trajectory file preparation pipeline.
class init_pipeline
{
public:
    init_pipeline(
        h5::file& store,
        simulation_config const& config,
        chain_definitions const& chains,
        std::uint32_t master_seed
    );
    void run();

private:
    void define_chains();
    void define_nucleolar_particles();
    void define_anatelophase_chains();
    void define_prometaphase_chains();

    void write_inputs();
    void write_interphase_particles();
    void write_anatelophase_particles();
    void write_prometaphase_particles();
    void write_seeds();

private:
    struct interphase_preparation
    {
        std::vector<interphase_particle_data>   particles;
        std::vector<chain_assignment>           chains;
        std::vector<std::size_t>                nor_indices;
        std::vector<std::array<std::size_t, 2>> nucleolar_bonds;
    };

    struct anatelophase_preparation
    {
        std::vector<mitotic_particle_type> particles;
        std::vector<chain_assignment>      chains;
    };

    struct prometaphase_preparation
    {
        std::vector<mitotic_particle_type>      particles;
        std::vector<chain_assignment>           chains;
        std::vector<std::array<std::size_t, 2>> sister_chromatids;
        std::array<md::point, 2>                pole_positions;
    };

    h5::file&                _store;
    simulation_config const  _config;
    chain_definitions const  _chains;
    std::uint32_t const      _master_seed;
    interphase_preparation   _interphase;
    anatelophase_preparation _anatelophase;
    prometaphase_preparation _prometaphase;
};


init_pipeline::init_pipeline(
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
init_pipeline::run()
{
    define_chains();
    define_nucleolar_particles();
    define_anatelophase_chains();
    define_prometaphase_chains();

    write_inputs();
    write_interphase_particles();
    write_anatelophase_particles();
    write_prometaphase_particles();
    write_seeds();
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
init_pipeline::define_chains()
{
    // Input:  (none)
    // Output: _interphase.{particles, nor_indices, chains}

    std::vector<std::pair<std::string, interphase_particle_type>> const tag_type_map = {
        {"anor", interphase_particle_type::active_nor},
        {"bnor", interphase_particle_type::silent_nor},
        {"cen", interphase_particle_type::centromere},
        {"A", interphase_particle_type::a},
        {"B", interphase_particle_type::b},
        {"u", interphase_particle_type::u},
    };

    for (auto const& chain : _chains.chains) {
        auto const chain_start = _interphase.particles.size();
        auto const chain_end = chain_start + chain.beads.size();

        for (auto const& bead : chain.beads) {
            auto const bead_index = _interphase.particles.size();
            auto bead_type = interphase_particle_type::unknown;

            for (auto const& [tag, type] : tag_type_map) {
                if (check_tag(bead.tags, tag)) {
                    bead_type = type;
                    break;
                }
            }

            if (bead_type == interphase_particle_type::active_nor) {
                _interphase.nor_indices.push_back(bead_index);
            }

            _interphase.particles.push_back({
                .a_factor = bead.a_factor,
                .b_factor = bead.b_factor,
                .type     = bead_type,
            });
        }

        _interphase.chains.push_back({
            .name        = chain.name,
            .start       = chain_start,
            .end         = chain_end,
            .kinetochore = {},
        });
    }
}


// Defines nucleolar particles in interphase simulation.
void
init_pipeline::define_nucleolar_particles()
{
    // Input:  _interphase.{particles, nor_indices}
    // Output: _interphase.{particles, nucleolar_bonds}

    for (auto const& nor_index : _interphase.nor_indices) {
        for (std::size_t i = 0; i < _config.interphase.nucleolus_bead_count; i++) {
            auto const nucleolus_index = _interphase.particles.size();
            _interphase.particles.push_back({
                .a_factor = _config.interphase.nucleolus_ab_factor.a,
                .b_factor = _config.interphase.nucleolus_ab_factor.b,
                .type     = interphase_particle_type::nucleolus,
            });
            _interphase.nucleolar_bonds.push_back({ nor_index, nucleolus_index });
        }
    }
}


// Defines coarse-grained chromosomal chains in anatelophase simulation.
void
init_pipeline::define_anatelophase_chains()
{
    // Input:  _interphase.{particles, chains}
    // Output: _anatelophase.{particle_types, chains}

    auto const coarse_graining = _config.mitotic_phase.coarse_graining;

    // Identify the centromeric region [start, end) of each chain.
    std::vector<std::array<std::size_t, 2>> centromere_ranges;

    for (auto const& assign : _interphase.chains) {
        std::size_t start = assign.start;
        std::size_t end = assign.end;
        bool seen = false;

        for (std::size_t i = assign.start; i < assign.end; i++) {
            if (_interphase.particles[i].type == interphase_particle_type::centromere) {
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

    // Define anatelophase chains by coarse-graining the interphase ones.

    for (std::size_t chain_index = 0; chain_index < _interphase.chains.size(); chain_index++) {
        auto const& assign = _interphase.chains[chain_index];
        auto const length = assign.end - assign.start;

        auto const coarse_length = length / coarse_graining;
        auto const coarse_chain_start = _anatelophase.particles.size();
        auto const coarse_chain_end = coarse_chain_start + coarse_length;

        // Model kinetochore-attached chromosomal region as a single bead at
        // the midpoint of the centromeric region of the original chain.
        auto const [cen_start, cen_end] = centromere_ranges[chain_index];
        auto const centromere_midpoint = (cen_start + cen_end) / 2;
        auto const kinetochore_offset = (centromere_midpoint - assign.start) / coarse_graining;

        std::optional<std::size_t> kinetochore_index;

        for (std::size_t bin = 0; bin < coarse_length; bin++) {
            auto const bead_index = _anatelophase.particles.size(); // next index

            auto type = mitotic_particle_type::arm;
            if (bin == kinetochore_offset) {
                type = mitotic_particle_type::kinetochore;
                kinetochore_index = bead_index;
            }

            _anatelophase.particles.push_back(type);
        }

        _anatelophase.chains.push_back({
            .name        = assign.name,
            .start       = coarse_chain_start,
            .end         = coarse_chain_end,
            .kinetochore = kinetochore_index,
        });
    }
}


// Defines coarse-grained chromosomal chains in anatelophase simulation.
void
init_pipeline::define_prometaphase_chains()
{
    // Input:  _anatelophase.{particles, chains}
    // Output: _prometaphase.{particles, chains, sister_chromatids, pole_positions}

    std::size_t const anatelo_chain_count = _anatelophase.chains.size();
    for (std::size_t i = 0; i < anatelo_chain_count; i++) {
        auto const target_chromatid = 2 * i;
        auto const sister_chromatid = target_chromatid + 1;
        _prometaphase.sister_chromatids.push_back({ target_chromatid, sister_chromatid });
    }

    std::string const sister_suffix = "-copy";
    for (auto const& anatelo_assign : _anatelophase.chains) {
        auto const chain_length = anatelo_assign.end - anatelo_assign.start;
        auto const kinetochore_offset = *anatelo_assign.kinetochore - anatelo_assign.start;

        auto const target_start = anatelo_assign.start * 2;
        auto const target_end = target_start + chain_length;
        auto const target_kinetochore = target_start + kinetochore_offset;

        auto const sister_start = target_end;
        auto const sister_end = sister_start + chain_length;
        auto const sister_kinetochore = sister_start + kinetochore_offset;

        _prometaphase.chains.push_back({
            .name        = anatelo_assign.name,
            .start       = target_start,
            .end         = target_end,
            .kinetochore = target_kinetochore,
        });

        _prometaphase.chains.push_back({
            .name        = anatelo_assign.name + sister_suffix,
            .start       = sister_start,
            .end         = sister_end,
            .kinetochore = sister_kinetochore,
        });

        using std::ptrdiff_t;
        auto const anatelo_beg = _anatelophase.particles.begin() + ptrdiff_t(anatelo_assign.start);
        auto const anatelo_end = _anatelophase.particles.begin() + ptrdiff_t(anatelo_assign.end);
        std::copy(anatelo_beg, anatelo_end, std::back_inserter(_prometaphase.particles));
        std::copy(anatelo_beg, anatelo_end, std::back_inserter(_prometaphase.particles));
    }

    md::point const origin = {};
    md::vector const pole_axis = _config.mitotic_phase.spindle_axis;
    md::point const target_pole = origin - pole_axis;
    md::point const sister_pole = origin + pole_axis;
    _prometaphase.pole_positions = { target_pole, sister_pole };
}


// Writes input parameters under "/metadata" hierarchy.
void
init_pipeline::write_inputs()
{
    // Note: "/metadata/config" is the serialization of the actual parameters
    // used in the simulation.
    _store.dataset<h5::u32>("/metadata/master_seed").write(_master_seed);
    _store.dataset<h5::str>("/metadata/config").write(format_simulation_config(_config));
    _store.dataset<h5::str>("/metadata/config_source").write(_config.source);
    _store.dataset<h5::str>("/metadata/chains_source").write(_chains.source);
}


// Writes paramerters for the interphase particles.
void
init_pipeline::write_interphase_particles()
{
    // Store metadata datasets that define interphase particles and topology:
    // - particle_types  (*) enum
    // - ab_factors  (*) float
    // - chain_names  (*) str
    // - chain_ranges  (*, 2) int
    // - nucleolar_bonds  (*, 2) int
    // under the metadata groups for the interphase and relaxation stages.

    h5::enums<int> const particle_types_enum = {
        {"unknown", int(interphase_particle_type::unknown)},
        {"a", int(interphase_particle_type::a)},
        {"b", int(interphase_particle_type::b)},
        {"u", int(interphase_particle_type::u)},
        {"centromere", int(interphase_particle_type::centromere)},
        {"active_nor", int(interphase_particle_type::active_nor)},
        {"silent_nor", int(interphase_particle_type::silent_nor)},
        {"nucleolus", int(interphase_particle_type::nucleolus)},
    };

    std::vector<int> particle_types;
    std::vector<std::array<double, 2>> ab_factors;
    for (auto const& particle : _interphase.particles) {
        particle_types.push_back(int(particle.type));
        ab_factors.push_back({ particle.a_factor, particle.b_factor });
    }

    std::vector<std::string> chain_names;
    std::vector<std::array<std::size_t, 2>> chain_ranges;
    for (auto const& assign : _interphase.chains) {
        chain_names.push_back(assign.name);
        chain_ranges.push_back({ assign.start, assign.end });
    }

    auto const make_interphase_path = [](std::string const& key) {
        return "/stages/interphase/metadata/" + key;
    };
    auto const make_relaxation_path = [](std::string const& key) {
        return "/stages/relaxation/metadata/" + key;
    };
    _store.dataset<h5::i32, 1>(make_interphase_path("particle_types"), particle_types_enum).write(particle_types);
    _store.dataset<h5::f32, 2>(make_interphase_path("ab_factors")).write(ab_factors);
    _store.dataset<h5::str, 1>(make_interphase_path("chain_names")).write(chain_names);
    _store.dataset<h5::i32, 2>(make_interphase_path("chain_ranges")).write(chain_ranges);
    _store.dataset<h5::i32, 2>(make_interphase_path("nucleolar_bonds")).write(_interphase.nucleolar_bonds);

    // The relaxation and interphase stage share some metadata.
    auto const share_metadata = [&, this](std::string const& key) {
        h5_link_path(_store, make_interphase_path(key), make_relaxation_path(key));
    };
    share_metadata("particle_types");
    share_metadata("ab_factors");
    share_metadata("chain_names");
    share_metadata("chain_ranges");
    share_metadata("nucleolar_bonds");
}


// Writes parameters for the anatelophase particles.
void
init_pipeline::write_anatelophase_particles()
{
    // Store anatelophase datasets that define anatelophase particles and topology:
    // - particle_types (*) enum
    // - chain_names (*) str
    // - chain_ranges (*, 2) int
    // - kinetochore_beads (*) int
    // under the metadata groups for the anaphase and telophase stages.

    h5::enums<int> const particle_types_enum = {
        {"unknown", int(mitotic_particle_type::unknown)},
        {"arm", int(mitotic_particle_type::arm)},
        {"kinetochore", int(mitotic_particle_type::kinetochore)},
    };

    std::vector<int> particle_types;
    for (auto const& type : _anatelophase.particles) {
        particle_types.push_back(int(type));
    }

    std::vector<std::string> chain_names;
    std::vector<std::array<std::size_t, 2>> chain_ranges;
    std::vector<std::size_t> kinetochore_beads;
    for (auto const& assign : _anatelophase.chains) {
        chain_names.push_back(assign.name);
        chain_ranges.push_back({ assign.start, assign.end });
        kinetochore_beads.push_back(assign.kinetochore.value_or(-1));
    }

    auto const make_anaphase_path = [](std::string const& key) {
        return "/stages/anaphase/metadata/" + key;
    };
    auto const make_telophase_path = [](std::string const& key) {
        return "/stages/telophase/metadata/" + key;
    };
    _store.dataset<h5::i32, 1>(make_anaphase_path("particle_types"), particle_types_enum).write(particle_types);
    _store.dataset<h5::str, 1>(make_anaphase_path("chain_names")).write(chain_names);
    _store.dataset<h5::i32, 2>(make_anaphase_path("chain_ranges")).write(chain_ranges);
    _store.dataset<h5::i32, 1>(make_anaphase_path("kinetochore_beads")).write(kinetochore_beads);

    // The anaphase and telophase stage share some metadata.
    auto const share_metadata = [&, this](std::string const& key) {
        h5_link_path(_store, make_anaphase_path(key), make_telophase_path(key));
    };
    share_metadata("particle_types");
    share_metadata("chain_names");
    share_metadata("chain_ranges");
}

// Writes parameters for the prometaphase particles.
void
init_pipeline::write_prometaphase_particles()
{
    // Store prometaphase datasets that define prometaphase particles and topology:
    // - particle_types (*) enum
    // - chain_names (*) str
    // - chain_ranges (*, 2) int
    // - kinetochore_beads (*) int
    // - sister_chromatids (*) int
    // - pole_positions (2, 3) float
    // under the metadata group for the prometaphase stage.

    h5::enums<int> const particle_types_enum = {
        {"unknown", int(mitotic_particle_type::unknown)},
        {"arm", int(mitotic_particle_type::arm)},
        {"kinetochore", int(mitotic_particle_type::kinetochore)},
    };

    std::vector<int> particle_types;
    for (auto const& type : _prometaphase.particles) {
        particle_types.push_back(int(type));
    }

    std::vector<std::string> chain_names;
    std::vector<std::array<std::size_t, 2>> chain_ranges;
    std::vector<std::size_t> kinetochore_beads;
    for (auto const& assign : _prometaphase.chains) {
        chain_names.push_back(assign.name);
        chain_ranges.push_back({ assign.start, assign.end });
        kinetochore_beads.push_back(assign.kinetochore.value_or(-1));
    }

    std::vector<std::array<md::scalar, 3>> pole_positions;
    for (auto const& position : _prometaphase.pole_positions) {
        pole_positions.push_back({ position.x, position.y, position.z });
    }

    auto const make_path = [](std::string const& key) {
        return "/stages/prometaphase/metadata/" + key;
    };
    _store.dataset<h5::i32, 1>(make_path("particle_types"), particle_types_enum).write(particle_types);
    _store.dataset<h5::str, 1>(make_path("chain_names")).write(chain_names);
    _store.dataset<h5::i32, 2>(make_path("chain_ranges")).write(chain_ranges);
    _store.dataset<h5::i32, 1>(make_path("kinetochore_beads")).write(kinetochore_beads);
    _store.dataset<h5::i32, 2>(make_path("sister_chromatids")).write(_prometaphase.sister_chromatids);
    _store.dataset<h5::f32, 2>(make_path("pole_positions")).write(pole_positions);
}


// Derives and stores the random seed used in each simulation stage.
void
init_pipeline::write_seeds()
{
    std::seed_seq seed_seq = { _master_seed };
    std::array<std::uint32_t, 3> seed_values;
    seed_seq.generate(seed_values.begin(), seed_values.end());

    auto const make_seed_path = [](std::string const& stage) {
        return "/stages/" + stage + "/metadata/seed";
    };
    _store.dataset<h5::u32>(make_seed_path("anaphase")).write(seed_values[0]);
    _store.dataset<h5::u32>(make_seed_path("interphase")).write(seed_values[1]);
    _store.dataset<h5::u32>(make_seed_path("prometaphase")).write(seed_values[2]);
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
    init_pipeline(store, config, chains, master_seed).run();
}
