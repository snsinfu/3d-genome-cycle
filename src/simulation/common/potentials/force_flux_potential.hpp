#pragma once

#include <cmath>

#include <md.hpp>


// Potential energy of the form
//
//    u(r) = b f arctan(b / r)
//
// This reduces to the inverse-square law, or Coulomb, potential when r >> b.
// The parameters describe
//
//    b:  Square root of the product of the reaction constant and the reaction
//        cross-section for the flux and a particle;
//    f:  Force excerted by a unit amount of flux.
//
struct force_flux_potential
{
    md::scalar constant_force    = 0;
    md::scalar reactive_distance = 1;

    inline md::scalar evaluate_energy(md::vector const& r) const
    {
        md::scalar const r1 = r.norm();
        return constant_force * reactive_distance * std::atan2(reactive_distance, r1);
    }

    inline md::vector evaluate_force(md::vector const& r) const
    {
        md::scalar const r2 = r.squared_norm();
        md::scalar const r1 = std::sqrt(r2);
        md::scalar const r3 = r1 * r2;
        md::scalar const b2 = reactive_distance * reactive_distance;
        return constant_force * b2 / (b2 * r1 + r3) * r;
    }
};
