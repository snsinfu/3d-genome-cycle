#include <exception>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

#include <getopt.hpp>

#include "../common/simulation_store.hpp"
#include "transition.hpp"


enum class program_mode
{
    help,
    interphase,
    prometaphase,
    cycle,
};


struct program_options
{
    program_mode mode = program_mode::help;
    std::string  target_filename;
    std::string  source_filename;   // used in 'cycle' mode
};


class usage_error : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;
};


static void            show_short_usage(std::ostream& out);
static void            show_usage(std::ostream& out);
static program_options parse_options(int argc, char** argv);


int
main(int argc, char** argv)
{
    try {
        auto const options = parse_options(argc, argv);

        switch (options.mode) {
        case program_mode::help:
            show_usage(std::cout);
            break;

        case program_mode::interphase:
            {
                simulation_store store(options.target_filename);
                transition_interphase(store);
            }
            break;

        case program_mode::prometaphase:
            {
                simulation_store store(options.target_filename);
                transition_prometaphase(store);
            }
            break;

        case program_mode::cycle:
            {
                simulation_store prev(options.source_filename);
                simulation_store next(options.target_filename);
                transition_cycle(prev, next);
            }
            break;
        }

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
  transition interphase <trajectory.h5>
  transition prometaphase <trajectory.h5>
  transition cycle <prev.h5> <next.h5>
)";
}


void
show_usage(std::ostream& out)
{
    show_short_usage(out);
    out << R"(
Convert and transfer simulation state to the next stage.

interphase:
  Refine the final 'telophase' structure in the given trajectory
  file into the initial interphase 'relaxation' structure.

prometaphase:
  Coarse-grain the final 'interphase' structure in the given
  trajectory file into the initial 'prometaphase' structure with
  generated sister chromatids.

cycle:
  Copy main set of chromosomes in the final 'prometaphase' step
  into the initial 'anaphase' structure in another trajectory.
)";
}


program_options
parse_options(int argc, char** argv)
{
    program_options options;

    // Parse mode argument
    if (argc <= 1) {
        throw usage_error("mode must be specified");
    }
    std::string const mode_name = argv[1];

    std::map<std::string, program_mode> const mode_map = {
        {"help", program_mode::help},
        {"-h", program_mode::help},
        {"interphase", program_mode::interphase},
        {"prometaphase", program_mode::prometaphase},
        {"cycle", program_mode::cycle},
    };

    if (auto const mode_it = mode_map.find(mode_name); mode_it != mode_map.end()) {
        options.mode = mode_it->second;
    } else {
        throw usage_error("unrecognized mode");
    }

    // Parse mode-specific options (now mode is argv[0], the program name)
    argc--;
    argv++;

    switch (options.mode) {
    case program_mode::help:
        return options;

    case program_mode::interphase:
    case program_mode::prometaphase:
        if (argc != 2) {
            throw usage_error("single trajectory file must be specified");
        }
        options.target_filename = argv[1];
        break;

    case program_mode::cycle:
        if (argc != 3) {
            throw usage_error("prev and next trajectory files must be specified");
        }
        options.source_filename = argv[1];
        options.target_filename = argv[2];
        break;
    }

    return options;
}
