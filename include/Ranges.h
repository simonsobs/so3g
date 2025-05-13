#pragma once

#include <stdint.h>
#include <vector>
#include <utility>

#include <nanobind/nanobind.h>

using namespace std;

namespace nb = nanobind;


// Template class for working with Ranges

template <typename T>
class Ranges {
public:
    T count;
    T reference;
    vector<pair<T,T>> segments;

    // Construction
    Ranges() : count{0}, reference(0) {}
    Ranges(T count) : count{count}, reference(0) {}
    Ranges(T count, T reference) : count{count}, reference(reference) {}

    // Basic ops
    Ranges<T>& merge(const Ranges<T> &src);
    Ranges<T>& intersect(const Ranges<T> &src);
    Ranges<T>& add_interval(const T start, const T end);
    Ranges<T>& append_interval_no_check(const T start, const T end);
    Ranges<T>& buffer(const T buff);
    Ranges<T>& close_gaps(const T gap);
    Ranges<T> buffered(const T buff) const;
    Ranges<T> complement() const;
    Ranges<T> zeros_like() const;
    Ranges<T> ones_like() const;

    void cleanup();

    void safe_set_count(T count_);

    // Operators.
    Ranges<T> operator~() const;
    Ranges<T> operator-() const;
    void operator+=(const Ranges<T> &src);
    void operator*=(const Ranges<T> &src);
    Ranges<T> operator+(const Ranges<T> &src) const;
    Ranges<T> operator*(const Ranges<T> &src) const;

    string Description() const;
};


typedef Ranges<int32_t> RangesInt32;


void register_ranges(nb::module_ & m);
