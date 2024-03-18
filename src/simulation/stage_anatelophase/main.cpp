#include <iostream>

#include "../common/simulation_store.hpp"
#include "simulation_driver.hpp"


int
main(int argc, char** argv)
{
    if (argc != 2) {
        std::cerr << "usage: anatelophase <trajectory.h5>\n";
        return 1;
    }

    simulation_store store(argv[1]);
    simulation_driver driver(store);
    driver.run();

    return 0;
}
