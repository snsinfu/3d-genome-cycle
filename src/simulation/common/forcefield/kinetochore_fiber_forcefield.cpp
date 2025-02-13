#include "kinetochore_fiber_forcefield.hpp"


using this_type = kinetochore_fiber_forcefield;
using kinetochore_spec = kinetochore_fiber_forcefield::kinetochore_spec;


this_type&
kinetochore_fiber_forcefield::set_pole_position(md::point const& pos)
{
    _pole_position = pos;
    return *this;
}


void
kinetochore_fiber_forcefield::add_kinetochore(kinetochore_spec const& spec)
{
    _kinetochores.push_back(spec);
}


namespace
{
    // The time evolution of microtubule length l follows the equation
    //
    //   dl/dt = a - kl
    //
    // where a and k are the polymerization and depolymerization rates,
    // respectively. When a particle is attached to the plus end of the
    // microtubule, its one-dimensional position relative to the minus
    // end x is constrained at x = l so that
    //
    //   dx/dt = -k (x - a/k) ,
    //
    // or, the particle effectively obeys the overdamped dynamics under a
    // spring potential. Hence, given mobility μ, the effective potential
    // for the particle is
    //
    //   u(x) = K/2 (x - b) ,
    //   K = k / μ ,
    //   b = a / k .
    //
    // b is the stationary length of the microtubule.
    //
    static inline md::spring_potential
    make_potential(kinetochore_spec const& spec)
    {
        return md::spring_potential {
            .spring_constant      = spec.decay_rate / spec.mobility,
            .equilibrium_distance = spec.stationary_length,
        };
    }
}


md::scalar
kinetochore_fiber_forcefield::compute_energy(md::system const& system)
{
    md::array_view<md::point const> const positions = system.view_positions();
    md::scalar energy = 0;

    for (kinetochore_spec const& spec : _kinetochores) {
        md::vector const r = positions[spec.particle_index] - _pole_position;
        energy += make_potential(spec).evaluate_energy(r);
    }

    return energy;
}


void
kinetochore_fiber_forcefield::compute_force(
    md::system const& system,
    md::array_view<md::vector> forces
)
{
    md::array_view<md::point const> const positions = system.view_positions();

    for (kinetochore_spec const& spec : _kinetochores) {
        md::vector const r = positions[spec.particle_index] - _pole_position;
        forces[spec.particle_index] += make_potential(spec).evaluate_force(r);
    }
}
