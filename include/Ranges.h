#pragma once

#include <stdint.h>

#include <nanobind/nanobind.h>

#include "numpy_assist.h"

using namespace std;

namespace nb = nanobind;


// Template class for working with intervals -- pairs of objects of
// the same (well-ordered) type, with operations defined that support
// intersection, union, and subtraction.

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

    static Ranges<T> * from_array(const nb::object &src, const nb::object &count);

    // Basic ops
    Ranges<T>& merge(const Ranges<T> &src);
    Ranges<T>& intersect(const Ranges<T> &src);
    Ranges<T>& add_interval(const T start, const T end);
    Ranges<T>& _add_interval_numpysafe(const nb::object start,
                                       const nb::object end);
    Ranges<T>& append_interval_no_check(const T start, const T end);
    Ranges<T>& buffer(const T buff);
    Ranges<T>& close_gaps(const T gap);
    Ranges<T> buffered(const T buff) const;
    Ranges<T> complement() const;
    Ranges<T> zeros_like() const;
    Ranges<T> ones_like() const;

    void cleanup();

    nb::object ranges() const;

    Ranges<T> getitem(nb::object indices);
    nb::object shape();
    void safe_set_count(T count_);

    // Operators.
    Ranges<T> operator~() const;
    Ranges<T> operator-() const;
    void operator+=(const Ranges<T> &src);
    void operator*=(const Ranges<T> &src);
    Ranges<T> operator+(const Ranges<T> &src) const;
    Ranges<T> operator*(const Ranges<T> &src) const;

    // Special conversions.
    static nb::object from_bitmask(const nb::object &src, int n_bits);
    static nb::object bitmask(const nb::list &ivlist, int n_bits);
    static nb::object from_mask(const nb::object &src);
    nb::object mask();

    string Description() const;
};


// Support for working with RangesMatrix, which is basically just a list of Ranges
template <typename T>
vector<Ranges<T>> extract_ranges(const nb::object & ival_list) {
    vector<Ranges<T>> v(nb::len(ival_list));
    for (int i=0; i<nb::len(ival_list); i++)
        v[i] = nb::cast<Ranges<T>>(ival_list[i]);
    return v;
}

typedef Ranges<int32_t> RangesInt32;


void register_ranges(nb::module_ & m);
