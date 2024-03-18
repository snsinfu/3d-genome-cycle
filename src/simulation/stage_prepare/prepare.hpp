#pragma once

#include <cstdint>

#include <h5.hpp>

#include "../common/simulation_config.hpp"
#include "chains.hpp"


void prepare_simulation_store(
    h5::file& store,
    simulation_config const& config,
    chain_definitions const& chains,
    std::uint32_t master_seed
);
