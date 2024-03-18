#include <md.hpp>

#include "../common/particle_data.hpp"

#include "simulation_driver.hpp"


void
simulation_driver::setup_particles()
{
    // Set particle data.
    _system.add_attribute(particle_data_attribute);

    for (md::index i = 0; i < _design.particles.size(); i++) {
        auto part = _system.add_particle();
        part.view(particle_data_attribute) = _design.particles[i];
    }

    // Mobility varies between chromatin and nucleolar particles.
    auto mobilities = _system.view_mobilities();

    for (auto const& chain : _design.chains) {
        // Chromatin
        for (md::index i = chain.start; i < chain.end; i++) {
            auto const& data = _design.particles[i];
            mobilities[i] = data.a_factor >= data.b_factor
                ? _config.a_core_mobility
                : _config.b_core_mobility;
        }
    }

    for (auto const& [nor, nuc] : _design.nucleolar_bonds) {
        mobilities[nuc] = _config.nucleolus_mobility;
    }
}
