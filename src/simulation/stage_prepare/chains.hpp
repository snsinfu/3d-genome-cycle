#pragma once

#include <cstdint>
#include <string>
#include <vector>


struct bead_definition
{
    std::uint32_t bin_start = 0;
    std::uint32_t bin_end   = 0;
    double        a_factor  = 0;
    double        b_factor  = 0;
    std::string   tags;
};


struct chain_definition
{
    std::string                  name;
    std::vector<bead_definition> beads;
};


struct chain_definitions
{
    std::vector<chain_definition> chains;
    std::string                   source;
};


chain_definitions load_chains(std::string const& filename);
