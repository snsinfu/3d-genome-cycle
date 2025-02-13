#pragma once

#include "../common/simulation_store.hpp"


void transition_interphase(simulation_store& store);
void transition_prometaphase(simulation_store& store);
void transition_cycle(simulation_store& prev, simulation_store& next);
