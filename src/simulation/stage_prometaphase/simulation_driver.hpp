#pragma once

#include <memory>
#include <random>
#include <string>

#include <md.hpp>

#include "../common/simulation_config.hpp"
#include "../common/simulation_store.hpp"


class simulation_driver
{
public:
    explicit simulation_driver(simulation_store& store);
    void run();

private:
    void setup();
    void setup_particles();
    void setup_repulsive_forcefield();
    void setup_connectivity_forcefield();
    void setup_sister_forcefield();
    void setup_kinetochore_forcefield();
    void setup_polar_ejection_forcefield();

    void run_initialization();
    void run_sampling();

    void update_context(md::step step);
    void print_progress(md::step step);
    void save_sample(md::step step);

    struct context_type
    {
        double time                   = 0;
        double kinetochore_attachment = 0;
    };

private:
    simulation_store&    _store;
    mitotic_phase_config _config;
    prometaphase_design  _design;
    md::system           _system;
    std::mt19937_64      _random;
    context_type         _context;
};
