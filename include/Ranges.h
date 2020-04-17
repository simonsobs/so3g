#pragma once

#include <G3Frame.h>
#include <G3Map.h>
#include <G3TimeStamp.h>

#include <stdint.h>

using namespace std;
namespace bp = boost::python;

// Template class for working with intervals -- pairs of objects of
// the same (well-ordered) type, with operations defined that support
// intersection, union, and subtraction.

template <typename T>
class Ranges : public G3FrameObject {
public:
    T count;
    T reference;
    vector<pair<T,T>> segments;

    // Construction
    Ranges() : count{0}, reference(0) {}
    Ranges(T count) : count{count}, reference(0) {}
    Ranges(T count, T reference) : count{count}, reference(reference) {}

    static Ranges<T> from_array(const bp::object &src, int count);

    // Basic ops
    Ranges<T>& merge(const Ranges<T> &src);
    Ranges<T>& intersect(const Ranges<T> &src);
    Ranges<T>& add_interval(const T start, const T end);
    Ranges<T>& append_interval_no_check(const T start, const T end);
    Ranges<T>& buffer(const T buff);
    Ranges<T> buffered(const T buff);
    Ranges<T> complement() const;

    void cleanup();

    bp::object ranges() const;

    Ranges<T> getitem(bp::object indices);
    bp::object shape();
    void safe_set_count(T count_);

    // Operators.
    Ranges<T> operator~() const;
    Ranges<T> operator-() const;
    void operator+=(const Ranges<T> &src);
    void operator*=(const Ranges<T> &src);
    Ranges<T> operator+(const Ranges<T> &src) const;
    Ranges<T> operator*(const Ranges<T> &src) const;

    // Special conversions.
    static bp::object from_bitmask(const bp::object &src, int n_bits);
    static bp::object bitmask(const bp::list &ivlist, int n_bits);
    static bp::object from_mask(const bp::object &src);
    bp::object mask();

    // Required for G3FrameObjects.
    string Description() const;
    template <class A> void serialize(A &ar, unsigned v);
};


typedef Ranges<int32_t> RangesInt32;

G3_SERIALIZABLE(RangesInt32, 0);

G3MAP_OF(std::string, RangesInt32, MapRangesInt32);

