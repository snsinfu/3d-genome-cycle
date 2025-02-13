#pragma once


// Creates an array_view into a subspan [start,end) of a vector.
template<typename T>
md::array_view<T>
view_slice(std::vector<T>& vec, std::size_t start, std::size_t end)
{
    return md::array_view<T>(vec).subview(start, end - start);
}


template<typename T>
md::array_view<T const>
view_slice(std::vector<T> const& vec, std::size_t start, std::size_t end)
{
    return md::array_view<T const>(vec).subview(start, end - start);
}
