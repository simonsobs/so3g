#pragma once

#include <stdint.h>
#include <vector>
#include <utility>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

using namespace std;

namespace nb = nanobind;

// G3Time is just an int64_t
typedef int64_t G3Time;


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

    //static Intervals<T> from_array(const bp::numpy::ndarray &src);
    static Intervals<T> from_array(const bp::object &src);

    // Basic ops
    Intervals<T>& merge(const Intervals<T> &src);
    Intervals<T>& intersect(const Intervals<T> &src);
    Intervals<T>& add_interval(const T start, const T end);
    Intervals<T>& append_interval_no_check(const T start, const T end);
    Intervals<T> complement() const;

    void set_domain(T start, T end);
    pair<T,T> get_domain();

    void cleanup();

    nb::ndarray<T> array() const;

    Intervals<T> getitem(bp::object indices);

    // Operators.
    Intervals<T> operator~() const;
    Intervals<T> operator-() const;
    void operator+=(const Intervals<T> &src);
    void operator-=(const Intervals<T> &src);
    Intervals<T> operator+(const Intervals<T> &src) const;
    Intervals<T> operator-(const Intervals<T> &src) const;
    Intervals<T> operator*(const Intervals<T> &src) const;

    // Special conversions.
    static bp::object from_mask(const bp::object &src, int n_bits);
    static bp::object mask(const bp::list &ivlist, int n_bits);

    string Description() const;
};


typedef Intervals<double> IntervalsDouble;
typedef Intervals<int64_t> IntervalsInt;
typedef Intervals<int32_t> IntervalsInt32;
typedef Intervals<G3Time> IntervalsTime;


void register_intervals(nb::module_ & m);

