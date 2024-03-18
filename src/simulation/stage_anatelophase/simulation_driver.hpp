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
    void setup_forcefield();
    void setup_repulsive_forcefield();
    void setup_connectivity_forcefield();
    void setup_dragging_forcefield();
    void setup_packing_forcefield();

    void run_initialization();
    void run_dragging_stage();
    void run_packing_stage();
    void run_refinement();

    void print_progress(std::string const& stage, md::step step);

private:
    simulation_store&               _store;
    anatelophase_config             _config;
    metaphase_design                _design;
    md::system                      _system;
    std::mt19937_64                 _random;
    std::shared_ptr<md::forcefield> _dragging_forcefield;
    std::shared_ptr<md::forcefield> _packing_forcefield;
};
