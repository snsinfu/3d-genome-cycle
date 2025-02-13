#include <iostream>
#include <vector>

#include <md.hpp>
#include <spline.hpp>

#include "misc.hpp"
#include "transition.hpp"


namespace
{
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

        for (md::index i = 0; i < new_chain.size(); i++) {
            auto const t = (0.5 + double(i)) / double(new_chain.size());
            new_chain[i] = { x_spline(t), y_spline(t), z_spline(t) };
        }
    }
}


// Refines coarse telophase structure into interphase one.
void
transition_interphase(simulation_store& store)
{
    std::clog << "Refining structure... ";

    std::string const prev_stage = "telophase";
    std::string const next_stage = "relaxation";

    auto const telophase_design = store.load_anatelophase_design();
    auto const interphase_design = store.load_interphase_design();

    store.set_stage(prev_stage);
    std::vector<md::point> const telophase_positions =
        store.load_positions(store.load_steps().back());
    std::vector<md::point> interphase_positions(interphase_design.particles.size());

    for (md::index chain_index = 0; chain_index < telophase_design.chains.size(); chain_index++) {
        resample_chain(
            view_slice(
                telophase_positions,
                telophase_design.chains[chain_index].start,
                telophase_design.chains[chain_index].end
            ),
            view_slice(
                interphase_positions,
                interphase_design.chains[chain_index].start,
                interphase_design.chains[chain_index].end
            )
        );
    }

    for (auto const& [nor, nuc] : interphase_design.nucleolar_bonds) {
        interphase_positions[nuc] = interphase_positions[nor];
    }

    store.set_stage(next_stage);
    store.save_positions(0, interphase_positions);

    std::clog << "OK\n";
}
