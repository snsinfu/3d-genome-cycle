#pragma once

#include <cstddef>
#include <array>
#include <vector>

#include <jsoncons/json.hpp>

#include <md.hpp>


template<typename Json, typename T>
struct coords_json_type_traits
{
    using value_type = T;
    using scalar_type = md::scalar;
    using allocator_type = typename Json::allocator_type;

    static bool is(const Json& j) noexcept
    {
        return j.is_array() && j.size() == 3;
    }

    static value_type as(const Json& j)
    {
        auto const coords = j.template as<std::array<scalar_type, 3>>();
        return {coords[0], coords[1], coords[2]};
    }

    static Json to_json(value_type const& value, allocator_type alloc = {})
    {
        Json j(jsoncons::json_array_arg_t(), jsoncons::semantic_tag::none, alloc);
        j.push_back(value.x);
        j.push_back(value.y);
        j.push_back(value.z);
        return j;
    }
};


template<class Json>
struct jsoncons::json_type_traits<Json, md::point> : coords_json_type_traits<Json, md::point>
{
};


template<class Json>
struct jsoncons::json_type_traits<Json, md::vector> : coords_json_type_traits<Json, md::vector>
{
};
