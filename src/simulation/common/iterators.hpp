#pragma once


template<typename Iterator>
class iterator_range
{
public:
    using iterator = Iterator;

    inline iterator_range(iterator const& begin, iterator const& end)
    : _begin(begin), _end(end)
    {
    }

    inline iterator const& begin() const noexcept { return _begin; }
    inline iterator const&   end() const noexcept { return _end; }

private:
    iterator _begin;
    iterator _end;
};
