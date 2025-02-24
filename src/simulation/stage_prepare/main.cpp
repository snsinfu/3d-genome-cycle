#include <cstddef>
#include <cstdint>
#include <exception>
#include <iostream>
#include <map>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include <getopt.hpp>
#include <h5.hpp>
#include <tsv.hpp>

#include "../common/simulation_config.hpp"
#include "chains.hpp"
#include "io.hpp"
#include "prepare.hpp"


struct program_options
{
    std::string                  trajectory_filename;
    std::string                  config_filename;
    std::string                  chains_filename;
    std::optional<std::uint32_t> seed;
    bool                         help = false;
};


class usage_error : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;
};


static program_options    parse_options(int argc, char** argv);
static void               show_short_usage(std::ostream& out);
static void               show_usage(std::ostream& out);
static simulation_config  load_config(std::string const& filename);
static std::uint32_t      make_master_seed(program_options const& options);


int
main(int argc, char** argv)
{
    try {
        auto const options = parse_options(argc, argv);

        if (options.help) {
            show_usage(std::cout);
            return 0;
        }

        h5::file store(options.trajectory_filename, "w");
        auto const config = load_config(options.config_filename);
        auto const chains = load_chains(options.chains_filename);
        auto const master_seed = make_master_seed(options);
        prepare_simulation_store(store, config, chains, master_seed);

        return 0;
    } catch (usage_error const& err) {
        show_short_usage(std::cerr);
        return 1;
    } catch (std::exception const& err) {
        std::cerr << "error: " << err.what() << '\n';
        return 1;
    }
}


void
show_short_usage(std::ostream& out)
{
    out << R"(usage:
  prepare init [-s seed] -o <trajectory.h5> <config.json> <chains.tsv>
  prepare interphase <trajectory.h5>
  prepare prometaphase <trajectory.h5>
)";
}


void
show_usage(std::ostream& out)
{
    show_short_usage(out);
    out << R"(
Create and prepare a trajectory file.

init options:
  -s seed           Specify random seed (default: random)
  -o trajectory.h5  Trajectory file to create (required)
  -h                Show usage message and exit

init positional arguments:
  config.json       Simulation parameters
  chains.tsv        TSV file specifying chromosome chains to be simulated

interphase/prometaphase positional arguments:
  trajectory.h5     Existing trajectory file to use
)";
}


program_options
parse_options(int argc, char** argv)
{
    program_options options;
    cxx::getopt getopt;

    for (int ch; (ch = getopt(argc, argv, "s:o:h")) != -1; ) {
        switch (ch) {
        case 'h':
            // Do not parse the remaining command-line if help is requested.
            options.help = true;
            return options;

        case 's':
            options.seed = std::uint32_t(std::stol(getopt.optarg));
            break;

        case 'o':
            options.trajectory_filename = getopt.optarg;
            break;

        default:
            throw usage_error("bad option");
        }
    }

    argc -= getopt.optind;
    argv += getopt.optind;

    if (argc != 2) {
        throw usage_error("config and chain definition files must be specified");
    }

    options.config_filename = argv[0];
    options.chains_filename = argv[1];

    return options;
}


simulation_config
load_config(std::string const& filename)
{
    return parse_simulation_config(load_text(filename));
}


std::uint32_t
make_master_seed(program_options const& options)
{
    if (options.seed) {
        return *options.seed;
    }

    std::random_device random_source;
    return random_source();
}
