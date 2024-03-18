#pragma once

#include <string>

#include <md.hpp>


struct ab_factor_config
{
    md::scalar a = 0;
    md::scalar b = 0;
};


struct anatelophase_config
{
    // Extra coarse-graining in the anatelophase simulation
    md::index  coarse_graining = 100;

    // Overdamped Langevin dynamics
    md::scalar temperature       = 1;
    md::scalar timestep          = 1e-4;
    md::step   dragging_steps    = 200000;
    md::step   packing_steps     = 50000;
    md::step   sampling_interval = 1000;
    md::step   logging_interval  = 10000;

    // Initialization
    md::point  start_center = {0, 0, 5};
    md::scalar start_stddev = 1;

    // Repulsive polymer
    md::scalar core_diameter  = 0.3;
    md::scalar core_repulsion = 2;
    md::scalar bond_length    = 0.3;
    md::scalar bond_spring    = 1000;
    md::scalar bending_energy = 0;
    md::scalar core_mobility  = 1;

    // Dragging
    md::index  dragged_beads   = 3;
    md::scalar dragging_spring = 1.5;
    md::scalar packing_radius  = 1.2;
    md::scalar packing_spring  = 10;
};


struct interphase_config
{
    // Overdamped Langevin dynamics
    md::scalar temperature                  = 1;
    md::scalar timestep                     = 1e-5;
    md::step   steps                        = 700000;
    md::step   sampling_interval            = 1000;
    md::step   logging_interval             = 1000;
    md::scalar relaxation_spacestep         = 0.001;
    md::step   relaxation_steps             = 10000;
    md::step   relaxation_sampling_interval = 1000;
    md::step   relaxation_logging_interval  = 100;

    // Contact map
    md::scalar contactmap_distance        = 0.24;
    md::step   contactmap_update_interval = 20;
    md::step   contactmap_output_window   = 10;

    // Repulsive copolymer
    md::scalar a_core_diameter        = 0.30;
    md::scalar b_core_diameter        = 0.24;
    md::scalar a_core_repulsion       = 2.5;
    md::scalar b_core_repulsion       = 2.5;
    md::scalar a_core_bond_spring     = 100;
    md::scalar b_core_bond_spring     = 50;
    md::scalar a_core_bond_length     = 0;
    md::scalar b_core_bond_length     = 0;
    md::scalar a_core_2nd_bond_spring = 0;
    md::scalar b_core_2nd_bond_spring = 0;
    md::scalar a_core_mobility        = 1;
    md::scalar b_core_mobility        = 1;

    // Scheduled expansion
    md::scalar core_scale_init = 0.5;
    md::scalar core_scale_tau  = 0.5;
    md::scalar bond_scale_init = 0.5;
    md::scalar bond_scale_tau  = 0.5;

    // Nucleolar particles
    md::index        nucleolus_bead_count     = 2;
    ab_factor_config nucleolus_ab_factor      = {0, 10};
    md::scalar       nucleolus_bond_spring    = 10;
    md::scalar       nucleolus_bond_length    = 0;
    md::scalar       nucleolus_droplet_energy = 0.3;
    md::scalar       nucleolus_droplet_decay  = 0.2;
    md::scalar       nucleolus_droplet_cutoff = 0.4;
    md::scalar       nucleolus_mobility       = 1;

    // Ellipsoidal, moving wall
    md::vector       wall_semiaxes_init   = {2, 2, 2};
    md::vector       wall_semiaxes_spring = {3e4, 3e4, 3e4};
    md::scalar       wall_packing_spring  = 1000;
    ab_factor_config wall_ab_factor       = {0, 10};
    md::scalar       wall_mobility        = 2e-4;
};


struct simulation_config
{
    anatelophase_config anatelophase;
    interphase_config   interphase;
    std::string         source;
};


/** Parses JSON representation of `simulation_config` structure. */
simulation_config parse_simulation_config(std::string const& text);

/** Formats `simulation_config` structure as a JSON string. */
std::string format_simulation_config(simulation_config const& config);
