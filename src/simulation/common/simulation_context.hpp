#pragma once

#include <md.hpp>


// Aggregate of global state and statistics in an interphase simulation.
struct interphase_context
{
    // State
    md::scalar time          = 0;
    md::vector wall_semiaxes = {0, 0, 0};
    md::scalar core_scale    = 1;
    md::scalar bond_scale    = 1;

    // Statistics
    md::scalar mean_energy = 0;
    md::scalar wall_energy = 0;
};


// Aggregate of global state and statistics in a prometaphase simulation.
struct prometaphase_context
{
    struct microtubule_record
    {
        md::scalar length            = 0;
        md::scalar oscillation_phase = 0;
    };

    // State
    md::scalar                      time = 0;
    std::vector<microtubule_record> microtubules;
};
