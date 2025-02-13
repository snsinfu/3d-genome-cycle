#pragma once

#include <vector>

#include <md.hpp>


class kinetochore_fiber_forcefield : public md::forcefield
{
public:
    using this_type = kinetochore_fiber_forcefield;

    struct kinetochore_spec
    {
        md::index  particle_index    = 0;
        md::scalar mobility          = 1;
        md::scalar decay_rate        = 0;
        md::scalar stationary_length = 0;
    };

    this_type& set_pole_position(md::point const& pos);
    void       add_kinetochore(kinetochore_spec const& spec);

    md::scalar compute_energy(md::system const& system) override;
    void       compute_force(md::system const& system, md::array_view<md::vector> forces) override;

private:
    md::point                     _pole_position;
    std::vector<kinetochore_spec> _kinetochores;
};
