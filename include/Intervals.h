#pragma once

#include <stdint.h>
#include <vector>
#include <utility>

#include <nanobind/nanobind.h>

using namespace std;

namespace nb = nanobind;


// Template class for working with intervals -- pairs of objects of
// the same (well-ordered) type, with operations defined that support
// intersection, union, and subtraction.

template <typename T>
class Intervals {
public:
    typedef T value_type;
    pair<T,T> domain;
    vector<pair<T,T>> segments;

    // Construction
    Intervals();
    Intervals(pair<T,T> domain) : domain{domain} {}
    Intervals(T start, T end) : Intervals(make_pair(start,end)) {}

    // Basic ops
    Intervals<T>& merge(const Intervals<T> &src);
    Intervals<T>& intersect(const Intervals<T> &src);
    Intervals<T>& add_interval(const T start, const T end);
    Intervals<T>& append_interval_no_check(const T start, const T end);
    Intervals<T> complement() const;

    void set_domain(T start, T end);
    pair<T,T> get_domain();

    void cleanup();

    // Operators.
    Intervals<T> operator~() const;
    Intervals<T> operator-() const;
    void operator+=(const Intervals<T> &src);
    void operator-=(const Intervals<T> &src);
    Intervals<T> operator+(const Intervals<T> &src) const;
    Intervals<T> operator-(const Intervals<T> &src) const;
    Intervals<T> operator*(const Intervals<T> &src) const;

    string Description() const;
};


typedef Intervals<double> IntervalsDouble;
typedef Intervals<int64_t> IntervalsInt;
typedef Intervals<int32_t> IntervalsInt32;


void register_intervals(nb::module_ & m);

