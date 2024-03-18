#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include <h5.hpp>
#include <md.hpp>

#include "particle_data.hpp"
#include "simulation_config.hpp"
#include "simulation_context.hpp"


struct chain_range
{
    std::string name;
    md::index   start = 0;
    md::index   end   = 0;
};


struct nucleolar_bond
{
    md::index nor_index = 0;
    md::index nuc_index = 0;
};


struct index_range
{
    md::index begin = 0;
    md::index end   = 0;
};


struct metaphase_design
{
    std::uint64_t            seed = 0;
    std::vector<chain_range> chains;
    std::vector<md::index>   kinetochore_beads;
};


struct interphase_design
{
    std::uint64_t               seed = 0;
    std::vector<particle_data>  particles;
    std::vector<chain_range>    chains;
    std::vector<nucleolar_bond> nucleolar_bonds;
};


class simulation_store
{
public:
    // Constructor takes the filename of the HDF5 file to operate on and opens
    // it in read-write mode.
    explicit simulation_store(std::string const& filename);

    // Metadata
    simulation_config  load_config();
    metaphase_design   load_metaphase_design();
    interphase_design  load_interphase_design();

    // Set the HDF5 hierarchy to save snapshots to.
    void set_stage(std::string const& name);

    // Snapshot
    void append_frame(md::step step);
    void save_positions(md::step step, md::array_view<md::point const> positions);
    void save_context(md::step step, simulation_context const& context);
    void save_contacts(md::step step, std::vector<std::array<std::uint32_t, 3>> const& contacts);

    std::vector<md::point> load_positions(md::step step);
    simulation_context     load_context(md::step step);

private:
    std::uint64_t               load_seed(std::string suffix = "");
    std::vector<chain_range>    load_chains(std::string suffix = "");
    std::vector<md::index>      load_kinetochore_beads();
    std::vector<nucleolar_bond> load_nucleolar_bonds();
    std::vector<particle_data>  load_interphase_particles();

    std::string make_stage_path(std::string sub_path) const;
    std::string make_stage_path(md::step step, std::string sub_path) const;

private:
    h5::file    _store;
    std::string _stage = "unknown";
};
