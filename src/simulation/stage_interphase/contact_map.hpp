#pragma once

#include <array>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include <md.hpp>


struct contact_pair
{
    std::uint32_t i = 0;
    std::uint32_t j = 0;
};


inline bool
operator==(contact_pair const& a, contact_pair const& b) noexcept
{
    return a.i == b.i && a.j == b.j;
}


inline bool
operator!=(contact_pair const& a, contact_pair const& b) noexcept
{
    return !(a == b);
}


template<>
struct std::hash<contact_pair>
{
    inline std::size_t operator()(contact_pair const& pair) const noexcept
    {
        std::size_t value;
        value = std::size_t(pair.i) << 32;
        value |= std::size_t(pair.j);
        return value;
    }
};


// Accumulates time-integrated contact map of moving points.
class contact_map
{
public:
    // Sets contact distance used in following updates.
    void set_contact_distance(md::scalar dist);

    // Returns the current contact distance.
    md::scalar contact_distance() const;

    // Clears contact map in-place.
    void clear();

    // Computes contact map of given points and adds to the ensemble.
    void update(md::array_view<const md::point> points);

    // Returns (i,j,v)-style contact map where i and j are indices and v is the
    // number of contacts.
    std::vector<std::array<std::uint32_t, 3>> accumulate();

private:
    md::scalar                                      _contact_distance = 0;
    std::unordered_map<contact_pair, std::uint32_t> _contacts;
};
