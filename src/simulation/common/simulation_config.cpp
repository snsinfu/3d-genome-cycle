#include <cstddef>
#include <array>
#include <string>
#include <vector>

#include <jsoncons/json.hpp>
#include <md.hpp>

#include "json_traits.hpp"
#include "simulation_config.hpp"


template<class Json>
struct jsoncons::json_type_traits<Json, ab_factor_config>
{
    using value_type = ab_factor_config;
    using scalar_type = md::scalar;
    using allocator_type = typename Json::allocator_type;

    static bool is(const Json& j) noexcept
    {
        return j.is_array() && j.size() == 2;
    }

    static value_type as(const Json& j)
    {
        auto const coords = j.template as<std::array<scalar_type, 2>>();
        return {coords[0], coords[1]};
    }

    static Json to_json(value_type const& value, allocator_type alloc = {})
    {
        Json j{jsoncons::json_array_arg_t{}, jsoncons::semantic_tag::none, alloc};
        j.push_back(value.a);
        j.push_back(value.b);
        return j;
    }
};



JSONCONS_N_MEMBER_TRAITS(
    mitotic_phase_config,

    // Required fields
    0,

    // Optional fields
    temperature,
    timestep,
    anaphase_steps,
    telophase_steps,
    prometaphase_steps,
    sampling_interval,
    logging_interval,

    anaphase_start_stddev,

    coarse_graining,
    core_diameter,
    core_repulsion,
    bond_length,
    bond_spring,
    bending_energy,
    penalize_centromere_bending,
    core_mobility,

    sister_separation,
    sister_spring,

    spindle_axis,
    kfiber_decay_rate_prometaphase,
    kfiber_decay_rate_anaphase,
    kfiber_length_prometaphase,
    kfiber_length_anaphase,
    polar_ejection_force,
    polar_ejection_cross_section,

    anaphase_spindle_shift,
    telophase_packing_radius,
    telophase_packing_spring,
    telophase_bond_spring_multiplier,
    telophase_bending_energy_multiplier
)



JSONCONS_N_MEMBER_TRAITS(
    interphase_config,

    // Required fields
    0,

    // Optional fields
    temperature,
    timestep,
    steps,
    sampling_interval,
    logging_interval,
    relaxation_spacestep,
    relaxation_steps,
    relaxation_sampling_interval,
    relaxation_logging_interval,

    contactmap_distance,
    contactmap_update_interval,
    contactmap_output_window,

    a_core_diameter,
    b_core_diameter,
    a_core_repulsion,
    b_core_repulsion,
    a_core_bond_spring,
    b_core_bond_spring,
    a_core_bond_length,
    b_core_bond_length,
    a_core_mobility,
    b_core_mobility,

    core_scale_init,
    core_scale_tau,
    bond_scale_init,
    bond_scale_tau,

    nucleolus_bead_count,
    nucleolus_ab_factor,
    nucleolus_bond_spring,
    nucleolus_bond_length,
    nucleolus_droplet_energy,
    nucleolus_droplet_decay,
    nucleolus_droplet_cutoff,
    nucleolus_mobility,

    wall_semiaxes_init,
    wall_semiaxes_spring,
    wall_packing_spring,
    wall_ab_factor,
    wall_mobility
)


JSONCONS_ALL_MEMBER_TRAITS(
    simulation_config,

    // Required fields
    mitotic_phase,
    interphase
)


simulation_config parse_simulation_config(std::string const& text)
{
    auto config = jsoncons::decode_json<simulation_config>(text);
    config.source = text;
    return config;
}


std::string format_simulation_config(simulation_config const& config)
{
    std::string text;
    jsoncons::encode_json(config, text);
    return text;
}
