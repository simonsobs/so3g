#pragma once

#include <stdint.h>
#include <limits>
#include <stdexcept>

#include "numpy_assist.h"

using namespace std;
namespace bp = boost::python;

// Template class for working with intervals -- pairs of objects of
// the same (well-ordered) type, with operations defined that support
// intersection, union, and subtraction.

template <typename T>
class Ranges {
public:
    T count;
    T reference;
    vector<pair<T,T>> segments;

    inline void validate_count(int64_t count_) const {
        if (count_ < 0 || count_ > static_cast<int64_t>(std::numeric_limits<T>::max())) {
            throw std::overflow_error("Count exceeds the maximum expressible value for type T");
        }
    }

    // Construction
    Ranges() : count{0}, reference(0) {}
    Ranges(int64_t count_) : reference(0) {
        validate_count(count_);
        count = count_;
    }
    
    Ranges(int64_t count_, T reference_) : reference(reference_) {
        validate_count(count_);
        count = count_;
    }

    static Ranges<T> from_array(const bp::object &src, const bp::object &count);

    // Basic ops
    Ranges<T>& merge(const Ranges<T> &src);
    Ranges<T>& intersect(const Ranges<T> &src);
    Ranges<T>& add_interval(const T start, const T end);
    Ranges<T>& _add_interval_numpysafe(const bp::object start,
                                       const bp::object end);
    Ranges<T>& append_interval_no_check(const T start, const T end);
    Ranges<T>& buffer(const T buff);
    Ranges<T>& close_gaps(const T gap);
    Ranges<T> buffered(const T buff) const;
    Ranges<T> complement() const;
    Ranges<T> zeros_like() const;
    Ranges<T> ones_like() const;
    
    void cleanup();

    bp::object ranges() const;

    Ranges<T> getitem(bp::object indices);
    bp::object shape();
    void safe_set_count(int64_t count_) {
        validate_count(count_);
        count = count_;
        cleanup();
    }

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

    string Description() const;
};


// Support for working with RangesMatrix, which is basically just a list of Ranges
template <typename T>
vector<Ranges<T>> extract_ranges(const bp::object & ival_list) {
    const int N = bp::len(ival_list);
    vector<Ranges<T>> v(N);
    for (int i=0; i<N ; i++)
        v[i] = bp::extract<Ranges<T>>(ival_list[i])();
    return v;
}

typedef Ranges<int32_t> RangesInt32;
