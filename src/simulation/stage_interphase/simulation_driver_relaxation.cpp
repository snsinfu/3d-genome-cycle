#include <algorithm>

#include <md.hpp>

#include "simulation_driver.hpp"


void
simulation_driver::run_relaxation()
{
    _store.set_stage("relaxation");

    // Load initial structure.
    {
        auto const init_positions = _store.load_positions(0);
        auto positions = _system.view_positions();
        std::copy(
            init_positions.begin(),
            init_positions.end(),
            positions.begin()
        );
    }

    auto callback = [this](md::step step) {
        // Calculating energy is expensive. So update stats only when needed.
        auto const with_logging = step % _config.relaxation_logging_interval == 0;
        auto const with_sampling = step % _config.relaxation_sampling_interval == 0;

        if (with_logging || with_sampling) {
            _context.mean_energy =
                _system.compute_energy() / md::scalar(_system.particle_count());
        }

        if (with_logging) {
            print_progress("relaxation", step);
        }

        if (with_sampling) {
            _store.save_positions(step, _system.view_positions());
            _store.save_context(step, _context);
            _store.append_frame(step);
        }
    };

    callback(0);

    md::simulate_brownian_dynamics(_system, {
        .temperature = _config.temperature,
        .timestep    = _config.timestep,
        .spacestep   = _config.relaxation_spacestep,
        .steps       = _config.relaxation_steps,
        .seed        = _random(),
        .callback    = callback,
    });
}
