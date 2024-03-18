#pragma once

#include <cstddef>
#include <array>
#include <stdexcept>
#include <vector>

#include <h5.hpp>
#include <md.hpp>


template<typename T, std::size_t N>
struct h5::buffer_traits<std::vector<std::array<T, N>>>
{
    using buffer_type = std::vector<std::array<T, N>>;
    using value_type = T;
    static constexpr std::size_t rank = 2;
    static constexpr std::size_t dimension = N;

    static void reshape(buffer_type& buffer, h5::shape<rank> const& shape)
    {
        if (shape.dims[1] != dimension) {
            throw std::invalid_argument("invalid dataset size");
        }
        buffer.resize(shape.dims[0]);
    }

    static h5::shape<rank> shape(buffer_type const& buffer)
    {
        return {buffer.size(), dimension};
    }

    static value_type* data(buffer_type& buffer)
    {
        return buffer.data()->data();
    }

    static value_type const* data(buffer_type const& buffer)
    {
        return buffer.data()->data();
    }
};
