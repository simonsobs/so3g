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
class Intervals {
public:
    pair<T,T> domain;
    vector<pair<T,T>> segments;

    // Construction
    Intervals();
    Intervals(pair<T,T> domain) : domain{domain} {}
    Intervals(T start, T end) : Intervals(make_pair(start,end)) {}
    Intervals(Intervals const & other);

    static Intervals<T> * from_array(const py::object &src);

    // Basic ops
    Intervals<T>& merge(const Intervals<T> &src);
    Intervals<T>& intersect(const Intervals<T> &src);
    Intervals<T>& add_interval(const T start, const T end);
    Intervals<T>& append_interval_no_check(const T start, const T end);
    Intervals<T> complement() const;

    void set_domain(T start, T end);
    pair<T,T> get_domain();

    void cleanup();

    py::object array() const;

    Intervals<T> getitem(py::object indices);

    // Operators.
    Intervals<T> operator~() const;
    Intervals<T> operator-() const;
    void operator+=(const Intervals<T> &src);
    void operator-=(const Intervals<T> &src);
    Intervals<T> operator+(const Intervals<T> &src) const;
    Intervals<T> operator-(const Intervals<T> &src) const;
    Intervals<T> operator*(const Intervals<T> &src) const;

    // Special conversions.
    static py::object from_mask(const py::object &src, int n_bits);
    static py::object mask(const py::list &ivlist, int n_bits);

    string Description() const;
};


typedef Intervals<double> IntervalsDouble;
typedef Intervals<int64_t> IntervalsInt;
typedef Intervals<int32_t> IntervalsInt32;


void register_intervals(py::module_ & m);
