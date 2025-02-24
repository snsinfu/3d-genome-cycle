#pragma once

#include <functional>
#include <random>

#include <md.hpp>

#include "../common/simulation_config.hpp"
#include "../common/simulation_context.hpp"
#include "../common/simulation_store.hpp"

#include "contact_map.hpp"


class simulation_driver
{
public:
    explicit simulation_driver(simulation_store& store);
    void     run();

private:
    void setup();
    void setup_particles();
    void setup_forcefield();
    void setup_repulsive_forcefield();
    void setup_connectivity_forcefield();
    void setup_loop_forcefield();
    void setup_nucleolus_forcefield();
    void setup_membrane_forcefield();
    void setup_context();

    void run_relaxation();
    void run_simulation();

    void print_progress(std::string const& phase, md::step step);

    void update_core_scale();
    void update_wall_semiaxes();

private:
    simulation_store&  _store;
    interphase_config  _config;
    interphase_design  _design;
    interphase_context _context;
    contact_map        _contact_map;
    md::system         _system;
    std::mt19937_64    _random;

    std::function<md::vector()> _compute_packing_reaction;
};
