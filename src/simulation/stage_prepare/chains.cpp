#include <sstream>
#include <stdexcept>
#include <utility>

#include <tsv.hpp>

#include "chains.hpp"
#include "io.hpp"


using namespace std::string_literals;


chain_definitions
load_chains(std::string const& filename)
{
    struct record_type
    {
        std::string   chain;
        std::uint32_t start;
        std::uint32_t end;
        double        a;
        double        b;
        std::string   tags;
    };

    chain_definitions defs;
    std::vector<record_type> records;

    try {
        defs.source = load_text(filename);
        records = tsv::load<record_type>(std::istringstream(defs.source));
    } catch (std::exception const& ex) {
        throw std::runtime_error("cannot load chain definitions: "s + ex.what());
    }

    // Group contiguous records that have the same chain name.
    chain_definition current_chain;

    auto append_chain = [&] {
        if (current_chain.beads.size() > 0) {
            defs.chains.push_back(std::move(current_chain));
        }
    };

    for (auto const& record : records) {
        if (record.chain != current_chain.name) {
            append_chain();
            current_chain = {.name = record.chain, .beads = {}};
        }

        current_chain.beads.push_back({
            .bin_start = record.start,
            .bin_end   = record.end,
            .a_factor  = record.a,
            .b_factor  = record.b,
            .tags      = record.tags,
        });
    }
    append_chain();

    return defs;
}
