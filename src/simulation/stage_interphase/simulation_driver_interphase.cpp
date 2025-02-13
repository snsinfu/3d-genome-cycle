#include <cmath>

#include <md.hpp>

#include "simulation_driver.hpp"


void
simulation_driver::run_simulation()
{
    _store.set_stage("interphase");
    _store.clear_frames();

    // Note: No 

    auto callback = [this](md::step step) {
        _context.time = md::scalar(step) * _config.timestep;

        // Calculate energy only when it is required.
        auto const with_logging = step % _config.logging_interval == 0;
        auto const with_sampling = step % _config.sampling_interval == 0;
        auto const sample_frame = step / _config.sampling_interval;

        if (with_logging || with_sampling) {
            _context.mean_energy =
                _system.compute_energy() / md::scalar(_system.particle_count());
        }

        if (with_logging) {
            print_progress("interphase", step);
        }

        if (with_sampling) {
            _store.save_positions(step, _system.view_positions());
            _store.save_interphase_context(step, _context);
        }

        if (step % _config.contactmap_update_interval == 0) {
            _contact_map.update(_system.view_positions());
        }

        if (with_sampling && sample_frame % _config.contactmap_output_window == 0) {
            _store.save_contacts(step, _contact_map.accumulate());
            _contact_map.clear();
        }

        if (with_sampling) {
            _store.append_frame(step);
        }

        update_core_scale();
        update_wall_semiaxes();
    };

    callback(0);

    md::simulate_brownian_dynamics(_system, {
        .temperature = _config.temperature,
        .timestep    = _config.timestep,
        .steps       = _config.steps,
        .seed        = _random(),
        .callback    = callback,
    });
}


void
simulation_driver::update_core_scale()
{
    auto const t_bead = _context.time / _config.core_scale_tau;
    auto const t_bond = _context.time / _config.bond_scale_tau;
    _context.core_scale = 1 - (1 - _config.core_scale_init) * std::exp(-t_bead);
    _context.bond_scale = 1 - (1 - _config.bond_scale_init) * std::exp(-t_bond);

    _contact_map.set_contact_distance(_config.contactmap_distance * _context.core_scale);
}


void
simulation_driver::update_wall_semiaxes()
{
    auto& semiaxes = _context.wall_semiaxes;

    md::vector net_force;
    net_force += _compute_packing_reaction();
    net_force -= _config.wall_semiaxes_spring.hadamard(semiaxes);

    // Simulate ad-hoc overdamped motion of the wall.
    semiaxes += _config.timestep * _config.wall_mobility * net_force;
}
