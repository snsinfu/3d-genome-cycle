#include "simulation_driver.hpp"

#include <algorithm>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <vector>

#include <md.hpp>

#include "../common/forcefield/kinetochore_fiber_forcefield.hpp"


namespace
{
    template<typename T>
    std::shared_ptr<T>
    copy_shared(T const& obj)
    {
        return std::make_shared<T>(obj);
    }


    template<typename RNG>
    md::vector
    normal_vector(RNG& random)
    {
        std::normal_distribution<md::scalar> normal;
        return {normal(random), normal(random), normal(random)};
    }
}


simulation_driver::simulation_driver(simulation_store& store)
: _store(store)
, _config(store.load_config().mitotic_phase)
, _design(store.load_anatelophase_design())
, _random(_design.seed)
{
    setup();
}


void
simulation_driver::setup()
{
    setup_particles();
    setup_forcefield();
}


void
simulation_driver::setup_particles()
{
    for (auto const& chain : _design.chains) {
        for (md::index i = chain.start; i < chain.end; i++) {
            _system.add_particle({
                .mobility = _config.core_mobility,
            });
        }
    }
}


void
simulation_driver::setup_forcefield()
{
    setup_repulsive_forcefield();
    setup_connectivity_forcefield();
    setup_dragging_forcefield();
    setup_packing_forcefield();
}


void
simulation_driver::setup_repulsive_forcefield()
{
    // General repulsion for avoiding chain crossings.

    _system.add_forcefield(
        md::make_neighbor_pairwise_forcefield(
            md::softcore_potential<2, 3> {
                .energy   = _config.core_repulsion,
                .diameter = _config.core_diameter,
            }
        )
        .set_neighbor_distance(_config.core_diameter)
    );
}


void
simulation_driver::setup_connectivity_forcefield()
{
    // Spring bonds and bending cost.

    auto make_bonding_forcefield = [this](md::scalar bond_spring) {
        auto bonds = md::make_bonded_pairwise_forcefield(
            md::semispring_potential {
                .spring_constant      = bond_spring,
                .equilibrium_distance = _config.bond_length,
            }
        );
        for (auto const& chain : _design.chains) {
            bonds.add_bonded_range(chain.start, chain.end);
        }
        return copy_shared(bonds);
    };

    _anaphase_bonding_forcefield = make_bonding_forcefield(
        _config.bond_spring
    );
    _telophase_bonding_forcefield = make_bonding_forcefield(
        _config.bond_spring * _config.telophase_bond_spring_multiplier
    );

    auto make_bending_forcefield = [this](md::scalar bending_energy) {
        auto bends = md::make_bonded_triplewise_forcefield(
            md::cosine_bending_potential {
                .bending_energy = bending_energy,
            }
        );
        for (auto const& chain : _design.chains) {
            if (_config.penalize_centromere_bending) {
                bends.add_bonded_range(chain.start, chain.end);
            } else {
                bends.add_bonded_range(chain.start, chain.kinetochore);
                bends.add_bonded_range(chain.kinetochore + 1, chain.end);
            }
        }
        return copy_shared(bends);
    };

    _anaphase_bending_forcefield = make_bending_forcefield(
        _config.bending_energy
    );
    _telophase_bending_forcefield = make_bending_forcefield(
        _config.bending_energy * _config.telophase_bending_energy_multiplier
    );

    // Forcefield will switch when transitioning from anaphase to telpohase.
    _system.add_forcefield(_anaphase_bonding_forcefield);
    _system.add_forcefield(_anaphase_bending_forcefield);
}


void
simulation_driver::setup_dragging_forcefield()
{
    kinetochore_fiber_forcefield kfiber_forcefield;

    md::point const origin = {0, 0, 0};
    md::point const pole_position = origin + _config.anaphase_spindle_shift;
    kfiber_forcefield.set_pole_position(pole_position);

    auto const compute_mobility = [](md::scalar bead_mobility, chain_range const& chain) {
        return bead_mobility / md::scalar(chain.end - chain.start);
    };

    for (auto const& chain : _design.chains) {
        kfiber_forcefield.add_kinetochore({
            .particle_index    = chain.kinetochore,
            .mobility          = compute_mobility(_config.core_mobility, chain),
            .decay_rate        = _config.kfiber_decay_rate_anaphase,
            .stationary_length = _config.kfiber_length_anaphase,
        });
    }

    _dragging_forcefield = copy_shared(kfiber_forcefield);
}


void
simulation_driver::setup_packing_forcefield()
{
    // Weak harmonic well potential prevents open diffusion.

    _packing_forcefield = copy_shared(
        md::make_point_source_forcefield(
            md::semispring_potential {
                .spring_constant      = _config.telophase_packing_spring,
                .equilibrium_distance = _config.telophase_packing_radius,
            }
        )
        .set_point_source({0, 0, 0})
    );
}


void
simulation_driver::run()
{
    _store.clear_frames();
    _store.set_stage("anaphase");

    run_initialization();
    run_dragging_stage();
    run_packing_stage();

    std::clog << "Finished.\n";
}


void
simulation_driver::run_initialization()
{
    auto positions = _system.view_positions();

    // Initial structure may be given.
    if (_store.check_positions(0)) {
        auto const init_positions = _store.load_positions(0);
        if (init_positions.size() != positions.size()) {
            throw std::runtime_error("initial structure size mismatch");
        }
        std::copy(init_positions.begin(), init_positions.end(), positions.begin());
        return;
    }

    // Initialize chains as randomly-directed rods.
    md::point const origin = {};
    md::point const start_center = origin - _config.spindle_axis;

    for (auto const& chain : _design.chains) {
        auto const centroid =
            start_center + _config.anaphase_start_stddev * normal_vector(_random);
        auto const step =
            _config.bond_length * md::normalize(normal_vector(_random));

        auto pos = centroid - step * md::scalar(chain.end - chain.start) / 2;
        for (md::index i = chain.start; i < chain.end; i++) {
            positions[i] = pos;
            pos += step;
        }
    }
}


void
simulation_driver::run_dragging_stage()
{
    // Enable dragging force.
    _system.add_forcefield(_dragging_forcefield);

    // Save snapshots under /stages/telophase hierarchy.
    _store.set_stage("anaphase");

    auto const callback = [this](md::step step) {
        if (step % _config.sampling_interval == 0) {
            _store.save_positions(step, _system.view_positions());
            _store.append_frame(step);
        }

        if (step % _config.logging_interval == 0) {
            print_progress("anaphase", step);
        }
    };

    callback(0);

    md::simulate_brownian_dynamics(_system, {
        .temperature = _config.temperature,
        .timestep    = _config.timestep,
        .steps       = _config.anaphase_steps,
        .seed        = _random(),
        .callback    = callback,
    });
}


void
simulation_driver::run_packing_stage()
{
    // Microtubules detach from kinetochores and nuclear membrane forms. Turn
    // the dragging force off and enable packing force.
    _system.remove_forcefield(_dragging_forcefield);
    _system.add_forcefield(_packing_forcefield);

    // The rigidity of chain would change because of dissociation of condensin
    // from chromosomes during telophase.
    _system.remove_forcefield(_anaphase_bonding_forcefield);
    _system.remove_forcefield(_anaphase_bending_forcefield);
    _system.add_forcefield(_telophase_bonding_forcefield);
    _system.add_forcefield(_telophase_bending_forcefield);

    // Save snapshots under /stages/telophase hierarchy.
    _store.set_stage("telophase");

    auto const callback = [this](md::step step) {
        if (step % _config.sampling_interval == 0) {
            _store.save_positions(step, _system.view_positions());
            _store.append_frame(step);
        }

        if (step % _config.logging_interval == 0) {
            print_progress("telophase", step);
        }
    };

    callback(0);

    md::simulate_brownian_dynamics(_system, {
        .temperature = _config.temperature,
        .timestep    = _config.timestep,
        .steps       = _config.telophase_steps,
        .seed        = _random(),
        .callback    = callback,
    });
}


void
simulation_driver::print_progress(std::string const& stage, md::step step)
{
    auto const wallclock_time = std::time(nullptr);

    std::clog
        << "[" + stage + "] "
        << std::put_time(std::localtime(&wallclock_time), "%F %T")
        << '\t'
        << step
        << '\t'
        << "E: "
        << _system.compute_energy() / md::scalar(_system.particle_count())
        << '\n';
}
