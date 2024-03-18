#include <algorithm>
#include <cassert>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <random>

#include <md.hpp>

#include "../common/simulation_store.hpp"

#include "simulation_driver.hpp"


simulation_driver::simulation_driver(simulation_store& store)
: _store(store)
, _config(store.load_config().interphase)
, _design(store.load_interphase_design())
, _random(_design.seed)
{
    setup();
}


void
simulation_driver::setup()
{
    setup_particles();
    setup_forcefield();
    setup_context();
}


void
simulation_driver::setup_context()
{
    _context = {
        .time          = 0,
        .wall_semiaxes = _config.wall_semiaxes_init,
        .core_scale    = _config.core_scale_init,
        .bond_scale    = _config.bond_scale_init,
    };
}


void
simulation_driver::run()
{
    run_relaxation();
    run_simulation();
}


void
simulation_driver::print_progress(std::string const& phase, md::step step)
{
    auto const wallclock_time = std::time(nullptr);
    auto const effective_radius = std::cbrt(
        _context.wall_semiaxes.x *
        _context.wall_semiaxes.y *
        _context.wall_semiaxes.z
    );

    std::clog
        << "[" + phase + "] "
        << std::put_time(std::localtime(&wallclock_time), "%F %T")
        << '\t'
        << step
        << '\t'
        << "t: "
        << _context.time
        << '\t'
        << "R: "
        << effective_radius
        << '\t'
        << "E: "
        << _context.mean_energy
        << '\n';
}
