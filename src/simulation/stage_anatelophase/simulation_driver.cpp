#include <ctime>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <vector>

#include <md.hpp>

#include "simulation_driver.hpp"


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
, _config(store.load_config().anatelophase)
, _design(store.load_metaphase_design())
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

    auto bonds = _system.add_forcefield(
        md::make_bonded_pairwise_forcefield(
            md::semispring_potential {
                .spring_constant      = _config.bond_spring,
                .equilibrium_distance = _config.bond_length,
            }
        )
    );

    auto bends = _system.add_forcefield(
        md::make_bonded_triplewise_forcefield(
            md::cosine_bending_potential {
                .bending_energy = _config.bending_energy,
            }
        )
    );

    for (auto const& chain : _design.chains) {
        bonds->add_bonded_range(chain.start, chain.end);
        bends->add_bonded_range(chain.start, chain.end);
    }
}


void
simulation_driver::setup_dragging_forcefield()
{
    // Spindle core attracts centromeres.

    _dragging_forcefield = copy_shared(
        md::make_point_source_forcefield(
            md::harmonic_potential {
                .spring_constant = _config.dragging_spring,
            }
        )
        .set_point_source({0, 0, 0})
        .set_point_source_targets(_design.kinetochore_beads)
    );
}


void
simulation_driver::setup_packing_forcefield()
{
    // Weak harmonic well potential prevents open diffusion.

    _packing_forcefield = copy_shared(
        md::make_point_source_forcefield(
            md::semispring_potential {
                .spring_constant      = _config.packing_spring,
                .equilibrium_distance = _config.packing_radius,
            }
        )
        .set_point_source({0, 0, 0})
    );
}


void
simulation_driver::run()
{
    run_initialization();
    run_dragging_stage();
    run_packing_stage();
    run_refinement();

    std::clog << "Finished.\n";
}


void
simulation_driver::run_initialization()
{
    // Initialize chains as randomly-directed rods.
    auto positions = _system.view_positions();

    for (auto const& chain : _design.chains) {
        auto const centroid =
            _config.start_center
            + _config.start_stddev * normal_vector(_random);
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
        .steps       = _config.dragging_steps,
        .seed        = _random(),
        .callback    = callback,
    });
}


void
simulation_driver::run_packing_stage()
{
    // Enable packing force.
    _system.add_forcefield(_packing_forcefield);

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
        .steps       = _config.packing_steps,
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
