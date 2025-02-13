#include "simulation_driver.hpp"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>

#include "../common/forcefield/kinetochore_fiber_forcefield.hpp"
#include "../common/potentials/force_flux_potential.hpp"


simulation_driver::simulation_driver(simulation_store& store)
: _store(store)
, _config(store.load_config().mitotic_phase)
, _design(store.load_prometaphase_design())
{
    _random.seed(_design.seed);

    setup();
}


void
simulation_driver::setup()
{
    _store.set_stage("prometaphase");

    setup_particles();
    setup_repulsive_forcefield();
    setup_connectivity_forcefield();
    setup_sister_forcefield();
    setup_kinetochore_forcefield();
    setup_polar_ejection_forcefield();
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
simulation_driver::setup_repulsive_forcefield()
{
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

        if (_config.penalize_centromere_bending) {
            bends->add_bonded_range(chain.start, chain.end);
        } else {
            // Exclude kinetochore (centromere) bead.
            bends->add_bonded_range(chain.start, chain.kinetochore);
            bends->add_bonded_range(chain.kinetochore + 1, chain.end);
        }
    }
}


void
simulation_driver::setup_sister_forcefield()
{
    auto bonds = _system.add_forcefield(
        md::make_bonded_pairwise_forcefield(
            md::semispring_potential {
                .spring_constant      = _config.bond_spring,
                .equilibrium_distance = _config.sister_separation,
            }
        )
    );

    for (auto const& [target_index, sister_index] : _design.sister_chromatids) {
        bonds->add_bonded_pair(
            _design.chains[target_index].kinetochore,
            _design.chains[sister_index].kinetochore
        );
    }
}


void
simulation_driver::setup_kinetochore_forcefield()
{
    md::point const target_pole = _design.pole_positions[0];
    md::point const sister_pole = _design.pole_positions[1];

    auto& target_forcefield = *_system.add_forcefield(
        kinetochore_fiber_forcefield()
        .set_pole_position(target_pole)
    );

    auto& sister_forcefield = *_system.add_forcefield(
        kinetochore_fiber_forcefield()
        .set_pole_position(sister_pole)
    );

    auto const compute_mobility = [](md::scalar bead_mobility, chain_range const& chain) {
        return bead_mobility / md::scalar(chain.end - chain.start);
    };

    for (auto const [target_index, sister_index] : _design.sister_chromatids) {
        chain_range const& target_chain = _design.chains[target_index];
        chain_range const& sister_chain = _design.chains[sister_index];

        target_forcefield.add_kinetochore({
            .particle_index    = target_chain.kinetochore,
            .mobility          = compute_mobility(_config.core_mobility, target_chain),
            .decay_rate        = _config.kfiber_decay_rate_prometaphase,
            .stationary_length = _config.kfiber_length_prometaphase,
        });

        sister_forcefield.add_kinetochore({
            .particle_index    = sister_chain.kinetochore,
            .mobility          = compute_mobility(_config.core_mobility, sister_chain),
            .decay_rate        = _config.kfiber_decay_rate_prometaphase,
            .stationary_length = _config.kfiber_length_prometaphase,
        });
    }
}


void
simulation_driver::setup_polar_ejection_forcefield()
{
    force_flux_potential const potential = {
        .constant_force    = _config.polar_ejection_force,
        .reactive_distance = std::sqrt(_config.polar_ejection_cross_section),
    };

    md::point const target_pole = _design.pole_positions[0];
    md::point const sister_pole = _design.pole_positions[1];

    _system.add_forcefield(
        md::make_point_source_forcefield(potential)
        .set_point_source(target_pole)
    );

    _system.add_forcefield(
        md::make_point_source_forcefield(potential)
        .set_point_source(sister_pole)
    );
}


void
simulation_driver::run()
{
    _store.set_stage("prometaphase");
    _store.clear_frames();

    run_initialization();
    run_sampling();
}


void
simulation_driver::run_initialization()
{
    if (!_store.check_positions(0)) {
        throw std::runtime_error("no initial structure is given");
    }

    auto const init_positions = _store.load_positions(0);
    auto positions = _system.view_positions();
    if (init_positions.size() != positions.size()) {
        throw std::runtime_error("initial structure size mismatch");
    }

    std::copy(init_positions.begin(), init_positions.end(), positions.begin());
}


void
simulation_driver::run_sampling()
{
    auto const callback = [&](md::step step) {
        update_context(step);

        if (step % _config.logging_interval == 0) {
            print_progress(step);
        }

        if (step % _config.sampling_interval == 0) {
            save_sample(step);
        }
    };

    callback(0);

    md::simulate_brownian_dynamics(_system, {
        .temperature = _config.temperature,
        .timestep    = _config.timestep,
        .steps       = _config.prometaphase_steps,
        .seed        = _random(),
        .callback    = callback,
    });
}


void
simulation_driver::update_context(md::step step)
{
    _context.time = double(step) * _config.timestep;
}


void
simulation_driver::print_progress(md::step step)
{
    auto const wallclock_time = std::time(nullptr);
    auto const energy = _system.compute_potential_energy();

    std::clog
        << "[prometaphase] "
        << std::put_time(std::localtime(&wallclock_time), "%F %T")
        << '\t'
        << step
        << '\t'
        << "E: " << energy
        << '\n';
}


void
simulation_driver::save_sample(md::step step)
{
    _store.save_positions(step, _system.view_positions());
    _store.append_frame(step);
}
