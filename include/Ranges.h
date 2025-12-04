#pragma once

#include <stdint.h>

#include <pybind11/pybind11.h>

#include "numpy_assist.h"

using namespace std;

namespace py = pybind11;


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

    static Ranges<T> * from_array(const py::object &src, const T count);

    // Basic ops
    Ranges<T>& merge(const Ranges<T> &src);
    Ranges<T>& intersect(const Ranges<T> &src);
    Ranges<T>& add_interval(const T start, const T end);
    // Ranges<T>& _add_interval_numpysafe(const py::object start,
    //                                    const py::object end);
    Ranges<T>& append_interval_no_check(const T start, const T end);
    Ranges<T>& buffer(const T buff);
    Ranges<T>& close_gaps(const T gap);
    Ranges<T> buffered(const T buff) const;
    Ranges<T> complement() const;
    Ranges<T> zeros_like() const;
    Ranges<T> ones_like() const;

    void cleanup();

    py::object ranges() const;

    Ranges<T> getitem(py::object indices);
    py::object shape();
    void safe_set_count(T count_);

    // Operators.
    Ranges<T> operator~() const;
    Ranges<T> operator-() const;
    void operator+=(const Ranges<T> &src);
    void operator*=(const Ranges<T> &src);
    Ranges<T> operator+(const Ranges<T> &src) const;
    Ranges<T> operator*(const Ranges<T> &src) const;

    // Special conversions.
    static py::object from_bitmask(const py::object &src, int n_bits);
    static py::object bitmask(const py::list &ivlist, int n_bits);
    static py::object from_mask(const py::object &src);
    py::object mask();

    string Description() const;
};


// Support for working with RangesMatrix, which is basically just a list of Ranges
template <typename T>
vector<Ranges<T>> extract_ranges(const py::object & ival_obj) {
    py::list ival_list = py::cast<py::list>(ival_obj);
    vector<Ranges<T>> v(py::len(ival_list));
    for (int i=0; i<py::len(ival_list); i++)
        v[i] = py::cast<Ranges<T>>(ival_list[i]);
    return v;
}

typedef Ranges<int32_t> RangesInt32;


void register_ranges(py::module_ & m);
