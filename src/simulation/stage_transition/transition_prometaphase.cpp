#include <algorithm>
#include <iostream>
#include <vector>

#include <md.hpp>

#include "misc.hpp"
#include "transition.hpp"


namespace
{
    // Returns the number of particles in a simulation.
    md::index
    count_particles(prometaphase_design const& design)
    {
        md::index beads_count = 0;
        for (auto const& chain : design.chains) {
            beads_count += chain.end - chain.start;
        }
        return beads_count;
    }


    // Computes the centroid of a group of points.
    md::point
    compute_centroid(md::array_view<md::point const> const& points)
    {
        md::point const origin = {};
        md::vector coords = {};
        md::scalar count = 0;

        for (md::point const& point : points) {
            coords += point - origin;
            count++;
        }
        coords /= count;

        return origin + coords;
    }
}


void
transition_prometaphase(simulation_store& store)
{
    // Coarse-grain interphase structure and duplicate sister chromatids to
    // generate the initial structure for a prometaphase simulation.
    std::clog << "Coarse-graining structure... ";

    std::string const prev_stage = "interphase";
    std::string const next_stage = "prometaphase";

    auto const config = store.load_config();
    auto const interphase_design = store.load_interphase_design();
    auto const prometaphase_design = store.load_prometaphase_design();

    store.set_stage(prev_stage);
    std::vector<md::point> const interphase_positions =
        store.load_positions(store.load_steps().back());
    std::vector<md::point> prometaphase_positions(count_particles(prometaphase_design));

    // Spindle poles and sister chromatids are positioned as follows:
    //
    //            spindle_axis vector
    //            ------->
    //   o====[s]:[t]====o
    //
    //   o spindle poles
    //   [s] sister chromatid
    //   [t] target chromatid
    //   ==== microtubules
    //
    // An interphase (G1 phase) chromosome is coarse-grained into a chromatid,
    // and its displaced replica becomes the corresponding sister one. This
    // process models the S, G2, and prophase.
    md::vector const sister_displacement =
        -config.mitotic_phase.sister_separation * config.mitotic_phase.spindle_axis.normalize();
    std::size_t const coarse_graining = config.mitotic_phase.coarse_graining;

    auto const chrom_count = interphase_design.chains.size();
    for (std::size_t chrom_index = 0; chrom_index < chrom_count; chrom_index++) {
        auto const source_chain = interphase_design.chains[chrom_index];
        auto const [target_chain_index, sister_chain_index] = prometaphase_design.sister_chromatids[chrom_index];
        auto const target_chain = prometaphase_design.chains[target_chain_index];
        auto const sister_chain = prometaphase_design.chains[sister_chain_index];
        auto const coarse_length = target_chain.end - target_chain.start;
        auto const source_length = source_chain.end - source_chain.start;

        for (std::size_t offset = 0; offset < coarse_length; offset++) {
            auto const source_start = source_chain.start + coarse_graining * offset;
            auto const source_end = std::min(source_start + coarse_graining, source_start + source_length);
            auto const target_index = target_chain.start + offset;
            auto const sister_index = sister_chain.start + offset;
            auto const centroid = compute_centroid(view_slice(interphase_positions, source_start, source_end));
            prometaphase_positions[target_index] = centroid;
            prometaphase_positions[sister_index] = centroid + sister_displacement;
        }
    }

    store.set_stage(next_stage);
    store.save_positions(0, prometaphase_positions);

    std::clog << "OK\n";
}
