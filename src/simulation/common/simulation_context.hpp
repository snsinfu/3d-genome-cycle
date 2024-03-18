#pragma once

#include <md.hpp>


// Aggregate of simulation-global variables and statistics.
struct simulation_context
{
    // Global variables
    md::scalar time          = 0;
    md::vector wall_semiaxes = {0, 0, 0};
    md::scalar core_scale    = 1;
    md::scalar bond_scale    = 1;

    // Statistics
    md::scalar mean_energy = 0;
    md::scalar wall_energy = 0;
};
