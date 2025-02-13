#pragma once

#include <array>
#include <cstdint>
#include <stdexcept>
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
    md::index   start       = 0;
    md::index   end         = 0;
    md::index   kinetochore = 0;    // TODO
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


struct anatelophase_design
{
    std::uint64_t            seed = 0;
    std::vector<chain_range> chains;
};


struct interphase_design
{
    std::uint64_t               seed = 0;
    std::vector<particle_data>  particles;
    std::vector<chain_range>    chains;
    std::vector<nucleolar_bond> nucleolar_bonds;
};


struct prometaphase_design
{
    std::uint64_t                         seed = 0;
    std::vector<chain_range>              chains;
    std::vector<std::array<md::index, 2>> sister_chromatids;
    std::array<md::point, 2>              pole_positions;
};


class simulation_store
{
public:
    // Constructor takes the filename of the HDF5 file to operate on and opens
    // it in read-write mode.
    explicit simulation_store(std::string const& filename);

    // Metadata
    simulation_config   load_config();
    anatelophase_design load_anatelophase_design();
    interphase_design   load_interphase_design();
    prometaphase_design load_prometaphase_design();

    // Set the HDF5 hierarchy to save snapshots to.
    void set_stage(std::string const& name);

    // Snapshot
    void clear_frames();
    void append_frame(md::step step);
    void save_positions(md::step step, md::array_view<md::point const> positions);
    void save_interphase_context(md::step step, interphase_context const& context);
    void save_prometaphase_context(md::step step, prometaphase_context const& context);
    void save_contacts(md::step step, std::vector<std::array<std::uint32_t, 3>> const& contacts);

    std::vector<md::step>  load_steps();
    std::vector<md::point> load_positions(md::step step);
    interphase_context     load_interphase_context(md::step step);
    prometaphase_context   load_prometaphase_context(md::step step);
    bool                   check_positions(md::step step);

private:
    std::uint64_t               load_seed(std::string const& stage);
    std::vector<chain_range>    load_chains(std::string const& stage);
    std::vector<nucleolar_bond> load_nucleolar_bonds(std::string const& stage);
    std::vector<particle_data>  load_interphase_particles(std::string const& stage);

    std::vector<std::array<std::size_t, 2>> load_sister_chromatids(std::string const& stage);
    std::array<md::point, 2>                load_pole_positions(std::string const& stage);

    std::string locate_metadata_in(std::string const& stage, std::string const& key) const;
    std::string locate_data(std::string const& key) const;
    std::string locate_data(md::step step, std::string const& key) const;

private:
    h5::file    _store;
    std::string _stage = "unknown";
};


class simulation_store_error : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;
};
