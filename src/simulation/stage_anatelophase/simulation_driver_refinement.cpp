#include <iostream>
#include <vector>

#include <md.hpp>
#include <spline.hpp>

#include "simulation_driver.hpp"


namespace
{
    // Computes the number of particles in an interphase simulation.
    std::size_t
    compute_particle_count(interphase_design const& design)
    {
        std::size_t sum = 0;
        for (auto const& chain : design.chains) {
            sum += chain.end - chain.start;
        }
        sum += design.nucleolar_bonds.size();
        return sum;
    }


    // Fits a cubic spline curve to `chain` and resample points along the
    // curve to fill `new_chain`.
    void
    resample_chain(
        md::array_view<md::point const> chain,
        md::array_view<md::point> new_chain
    )
    {
        constexpr auto spline_bc = cubic_spline::not_a_knot;

        std::vector<double> ts, xs, ys, zs;
        for (auto const& point : chain) {
            auto const t = (0.5 + double(ts.size())) / double(chain.size());
            ts.push_back(t);
            xs.push_back(point.x);
            ys.push_back(point.y);
            zs.push_back(point.z);
        }

        cubic_spline x_spline(ts, xs, spline_bc);
        cubic_spline y_spline(ts, ys, spline_bc);
        cubic_spline z_spline(ts, zs, spline_bc);

        for (std::size_t i = 0; i < new_chain.size(); i++) {
            auto const t = (0.5 + double(i)) / double(new_chain.size());
            new_chain[i] = { x_spline(t), y_spline(t), z_spline(t) };
        }
    }
}


void
simulation_driver::run_refinement()
{
    std::clog << "Refining the structure\n";

    auto const new_design = _store.load_interphase_design();
    auto const new_particle_count = compute_particle_count(new_design);
    std::vector<md::point> new_positions(new_particle_count);

    auto const positions_view = _system.view_positions();
    auto const new_positions_view = md::array_view<md::point>(new_positions);

    for (std::size_t i = 0; i < _design.chains.size(); i++) {
        resample_chain(
            positions_view.subview(
                _design.chains[i].start,
                _design.chains[i].end - _design.chains[i].start
            ),
            new_positions_view.subview(
                new_design.chains[i].start,
                new_design.chains[i].end - new_design.chains[i].start
            )
        );
    }

    for (auto const& [nor, nuc] : new_design.nucleolar_bonds) {
        new_positions[nuc] = new_positions[nor];
    }

    _store.set_stage("relaxation");
    _store.save_positions(0, new_positions);
}
