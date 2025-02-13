#pragma once

#include <string>

#include <h5.hpp>


// Creates a soft link at name pointing to target.
inline void
h5_link_path(h5::file& store, std::string const& target, std::string const& name)
{
    // Allow intermediate groups to be automatically created.
    h5::unique_hid<H5Pclose> link_props = H5Pcreate(H5P_LINK_CREATE);
    if (link_props < 0) {
        throw h5::exception("failed to create link props");
    }
    if (H5Pset_create_intermediate_group(link_props, 1) < 0) {
        throw h5::exception("failed to configure link props");
    }

    auto const ret = H5Lcreate_soft(
        target.c_str(), store.handle(), name.c_str(), link_props, H5P_DEFAULT
    );
    if (ret < 0) {
        throw h5::exception("failed to create link " + name + " -> " + target);
    }
}
