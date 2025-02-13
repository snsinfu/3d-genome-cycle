#include <iostream>
#include <string>

#include <md.hpp>

#include "misc.hpp"
#include "transition.hpp"


namespace
{
    // Returns the number of particles in a simulation.
    md::index
    count_particles(anatelophase_design const& design)
    {
        md::index beads_count = 0;
        for (auto const& chain : design.chains) {
            beads_count += chain.end - chain.start;
        }
        return beads_count;
    }
}


void
transition_cycle(simulation_store& prev, simulation_store& next)
{
    std::clog << "Copying into a daughter cell... ";

    std::string const prev_stage = "prometaphase";
    std::string const next_stage = "anaphase";

    auto const metaphase_design = prev.load_prometaphase_design();
    auto const anaphase_design = next.load_anatelophase_design();
    auto const config = next.load_config();

    prev.set_stage(prev_stage);
    std::vector<md::point> const metaphase_positions =
        prev.load_positions(prev.load_steps().back());
    std::vector<md::point> anaphase_positions(count_particles(anaphase_design));

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
    // We take the target chromatid. Its associated spindle pole becomes the
    // new origin, so let us displace the target chromatid by the spindle_axis
    // vector.
    md::vector const displacement = -config.mitotic_phase.spindle_axis;

    for (md::index chrom_index = 0; chrom_index < anaphase_design.chains.size(); chrom_index++) {

        auto const [target_index, sister_index] = metaphase_design.sister_chromatids[chrom_index];
        auto const metaphase_chain = metaphase_design.chains[target_index];
        auto const anaphase_chain = anaphase_design.chains[chrom_index];
        auto const chain_length = metaphase_chain.end - metaphase_chain.start;

        for (md::index offset = 0; offset < chain_length; offset++) {
            anaphase_positions[anaphase_chain.start + offset] =
                metaphase_positions[metaphase_chain.start + offset]
                + displacement;
        }
    }

    next.set_stage(next_stage);
    next.save_positions(0, anaphase_positions);

    std::clog << "OK\n";
}
