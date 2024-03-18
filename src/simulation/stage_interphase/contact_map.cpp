#include <algorithm>
#include <array>
#include <cstdint>
#include <utility>
#include <vector>

#include <md.hpp>

#include "contact_map.hpp"


void
contact_map::set_contact_distance(md::scalar dist)
{
    _contact_distance = dist;
}


md::scalar
contact_map::contact_distance() const
{
    return _contact_distance;
}


void
contact_map::clear()
{
    _contacts.clear();
}


void
contact_map::update(md::array_view<const md::point> points)
{
    struct contact_output_iterator
    {
        std::unordered_map<contact_pair, std::uint32_t>& contacts;

        contact_output_iterator operator++(int)
        {
            return *this;
        }

        contact_output_iterator& operator*()
        {
            return *this;
        }

        void operator=(std::pair<md::index, md::index> const& pair)
        {
            auto const cpair = contact_pair {
                .i = std::uint32_t(pair.first),
                .j = std::uint32_t(pair.second),
            };
            contacts[cpair] += 1;
        }
    };

    md::neighbor_searcher<md::open_box> searcher({}, _contact_distance);
    searcher.set_points(points);
    searcher.search(contact_output_iterator(_contacts));
}


std::vector<std::array<std::uint32_t, 3>>
contact_map::accumulate()
{
    std::vector<std::array<std::uint32_t, 3>> contacts;
    contacts.reserve(_contacts.size());
    for (auto const& [pair, count] : _contacts) {
        contacts.push_back({ pair.i, pair.j, count });
    }

    // The contact map compresses much better when sorted.
    auto encode_ij = [](std::uint32_t i, std::uint32_t j) {
        return (std::uint64_t(i) << 32) | j;
    };
    auto compare_ij = [=](auto const& a, auto const& b) {
        return encode_ij(a[0], a[1]) < encode_ij(b[0], b[1]);
    };
    std::sort(contacts.begin(), contacts.end(), compare_ij);

    return contacts;
}
