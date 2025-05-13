
#include <iostream>
#include <sstream>
#include <limits>
#include <type_traits>
#include <cstring>
#include <tuple>

// Needed to access self in bindings
#include <nanobind/operators.h>

// Needed to work with numpy arrays in bindings
#include <nanobind/ndarray.h>

// Needed to access std::tuple return from slice.compute()
#include <nanobind/stl/tuple.h>

#include <bindings.h>
#include <exceptions.h>
#include <Intervals.h>


//
// Default constructors, explicitly defined for each type, to set a
// sensible (perhaps) default domain.
//

template <>
Intervals<double>::Intervals() {
    domain = make_pair(-std::numeric_limits<double>::infinity(),
                       +std::numeric_limits<double>::infinity());
}

template <>
Intervals<int64_t>::Intervals() {
    domain = make_pair(INT64_MIN, INT64_MAX);
}

template <>
Intervals<int32_t>::Intervals() {
    domain = make_pair(INT32_MIN, INT32_MAX);
}

//
// Some support templates for Description() -- these are broadly
// applicable so consider having them live more publicly.
//

// _ival_type_name() returns the standard G3 data type suffix.

template <typename T>
const char *_ival_type_name();
template <>
const char *_ival_type_name<int32_t>() { return "Int32"; }
template <>
const char *_ival_type_name<int64_t>() { return "Int"; }
template <>
const char *_ival_type_name<double> () { return "Double"; }

// _ival_cute_lim() allows standard limits (such as INT32_MAX) to be printed as such.

template <typename T>
std::string _ival_cute_lim(T val) {
    std::ostringstream s;
    s << val;
    return s.str();
}
template <>
std::string _ival_cute_lim(int32_t val) {
    std::ostringstream s;
    if (val == INT32_MIN)
        s << "INT32_MIN";
    else if (val == INT32_MAX)
        s << "INT32_MAX";
    else
        s << val;
    return s.str();
}
template <>
std::string _ival_cute_lim(int64_t val) {
    std::ostringstream s;
    if (val == INT64_MIN)
        s << "INT64_MIN";
    else if (val == INT64_MAX)
        s << "INT64_MAX";
    else
        s << val;
    return s.str();
}


template <typename T>
std::string Intervals<T>::Description() const
{
	std::ostringstream s;
	s << "Intervals" << _ival_type_name<T>() << "("
          << "domain=(" << _ival_cute_lim<T>(domain.first)
          << "," << _ival_cute_lim<T>(domain.second) << "), "
          << "ivals=" << segments.size() << ")";
	return s.str();
}

template <typename T>
pair<T,T> merge_pair(pair<T,T> a, pair<T,T> b) {
    pair<T,T> out = a;
    out.first = min(out.first, b.first);
    out.second = max(out.second, b.second);
    return out;
}

template <typename T>
void Intervals<T>::cleanup()
{
    auto p = segments.begin();
    while (p != segments.end()) {
        // Truncate intervals that cross the domain boundary.
        if (p->first < domain.first)
            p->first = domain.first;
        if (p->second > domain.second)
            p->second = domain.second;
        // Delete empty or negative intervals.
        if (p->first >= p->second) {
            segments.erase(p);
            continue;
        }
        // Check for overlap with next interval, if any.
        auto q = p+1;
        if (q == segments.end())
            break;
        if (p->second >= q->first) {
            p->second = max(q->second, p->second);
            segments.erase(q);
        } else
            ++p;
    }
}

template <typename T>
Intervals<T>& Intervals<T>::add_interval(const T start, const T end)
{
    // We can optimize this later.  For now just do something that is
    // obviously correct.
    auto p = lower_bound(segments.begin(), segments.end(), make_pair(start, end));
    segments.insert(p, make_pair(start, end));
    cleanup();

    return *this;
}

template <typename T>
Intervals<T>& Intervals<T>::append_interval_no_check(const T start, const T end)
{
    segments.push_back(make_pair(start, end));

    return *this;
}

template <typename T>
Intervals<T>& Intervals<T>::merge(const Intervals<T> &src)
{
    domain.first = max(domain.first, src.domain.first);
    domain.second = min(domain.second, src.domain.second);
    domain.second = max(domain.first, domain.second); // enforce first <= second.

    auto p0 = this->segments.begin();
    auto p1 = src.segments.begin();
    while (p0 != this->segments.end() && p1 != src.segments.end()) {
        if (p1->second < p0->first) {
            // p1 is entirely before p0.  Insert it and advance both iters.
            p0 = segments.insert(p0, *(p1++)) + 1;
        } else if (p0->second < p1->first) {
            // p0 is entirely before p1.
            p0++;
        } else {
            // The two intervals overlap, so merge them into p0.
            p0->first = min(p0->first, p1->first);
            p0->second = max(p0->second, p1->second);
            p1++;
        }
    }
    // Any trailing intervals in p1?
    while (p1 != src.segments.end())
        this->segments.push_back(*(p1++));

    cleanup();
    return *this;
}

//
// Machinery for converting between Interval vector<pair>
// representation and buffers (such as numpy arrays).
//

template <typename T>
static inline
pair<T,T> interval_pair(char *p1, char *p2) {
    return make_pair(*reinterpret_cast<T*>(p1),
                     *reinterpret_cast<T*>(p2));
}

//
// Implementation of the algebra
//

template <typename T>
Intervals<T>& Intervals<T>::intersect(const Intervals<T> &src)
{
    auto output = -(*this);
    output -= src;
    *this = output.complement();
    return *this;
}

template <typename T>
void Intervals<T>::set_domain(T start, T end)
{
    domain.first = start;
    domain.second = max(start, end);
    cleanup();
}

template <typename T>
pair<T,T> Intervals<T>::get_domain()
{
    return make_pair(domain.first, domain.second);
}

template <typename T>
Intervals<T> Intervals<T>::complement() const
{
    Intervals<T> output(domain.first, domain.second);

    T next_start = domain.first;
    for (auto p: segments) {
        output.segments.push_back(make_pair(next_start, p.first));
        next_start = p.second;
    }
    output.segments.push_back(make_pair(next_start, domain.second));
    output.cleanup();
    return output;
}

//
// Operators
//

template <typename T>
Intervals<T> Intervals<T>::operator~() const
{
    return complement();
}

template <typename T>
Intervals<T> Intervals<T>::operator-() const
{
    return complement();
}

template <typename T>
void Intervals<T>::operator+=(const Intervals<T> &src)
{
    merge(src);
}

template <typename T>
void Intervals<T>::operator-=(const Intervals<T> &src)
{
    merge(-src);
}

template <typename T>
Intervals<T> Intervals<T>::operator+(const Intervals<T> &src) const
{
    Intervals<T> output = *this;
    output += src;
    return output;
}

template <typename T>
Intervals<T> Intervals<T>::operator-(const Intervals<T> &src) const
{
    Intervals<T> output = *this;
    output -= src;
    return output;
}

template <typename T>
Intervals<T> Intervals<T>::operator*(const Intervals<T> &src) const
{
    auto output = *this;
    output.intersect(src);
    return output;
}


// Conversions between Intervals and mask arrays.  These are templated on both the
// intervals type and the mask array type.

template <
    typename C,
    typename M,
    typename std::enable_if_t<std::is_integral<C>::value, bool>* = nullptr
>
nb::list intervals_from_mask(
    nb::ndarray_view<M, 1, 'C'> & view, int n_bits
) {
    // Fast-access view of the array
    auto p = view.data();
    if (n_bits < 0) {
        // Use all bits
        n_bits = 8 * sizeof(M);
    }
    C count = view.shape(0);

    vector<Intervals<C>> output;
    vector<C> start;
    for (int bit = 0; bit < n_bits; ++bit) {
        output.push_back(Intervals<C>(0, count));
        start.push_back(-1);
    }

    auto last = p[0];
    last = 0;
    for (C i = 0; i < count; i++) {
        auto d = p[i] ^ last;
        for (int bit = 0; bit < n_bits; ++bit) {
            if (d & (1 << bit)) {
                if (start[bit] >= 0) {
                    output[bit].segments.push_back(
                        interval_pair<C>((char*)&start[bit], (char*)&i)
                    );
                    start[bit] = -1;
                } else {
                    start[bit] = i;
                }
            }
        }
        last = p[i];
    }
    for (int bit = 0; bit < n_bits; ++bit) {
        if (start[bit] >= 0) {
            output[bit].segments.push_back(
                interval_pair<C>((char*)&start[bit], (char*)&count)
            );
        }
    }

    // Once added to the list, we can't modify further.
    nb::list bits;
    for (auto i: output) {
        bits.append(i);
    }
    return bits;
}

template <
    typename C,
    typename M,
    typename std::enable_if_t<! std::is_integral<C>::value, bool>* = nullptr
>
nb::list intervals_from_mask(
    nb::ndarray_view<M, 1, 'C'> & view, int n_bits
) {
    throw general_agreement_exception("from_mask() is only valid for integral types");
    return nb::list();
}

template <
    typename C,
    typename M,
    typename std::enable_if_t<std::is_integral<C>::value, bool>* = nullptr
>
nb::object intervals_to_mask(nb::list ivlist, int n_bits) {
    if (n_bits < 0) {
        n_bits = ivlist.size();
    }
    // Ensure that the mask array type has sufficient number of bits.
    if (n_bits > nb::dtype<M>().bits) {
        std::ostringstream err;
        err << "Mask array type has " << nb::dtype<M>().bits <<
            " bits, which cannot hold the requested " << n_bits << " bits.";
        throw general_agreement_exception(err.str().c_str());
    }

    vector<Intervals<C>> ivals;
    vector<int> indexes;
    pair<C, C> domain;

    for (size_t i = 0; i < ivlist.size(); i++) {
        // Convert list element to a C++ object.  Note that this does make a copy
        // when pushed onto the std::vector.
        ivals.push_back(nb::cast<Intervals<C>>(ivlist[i]));
        indexes.push_back(0);
        if (i == 0) {
            domain = ivals[i].domain;
        } else if (domain != ivals[i].domain) {
            throw agreement_exception(
                "ivlist[0]", "all other ivlist[i]", "domain"
            );
        }
    }

    size_t n = domain.second - domain.first;

    // Create the returned array
    size_t shape[1] = {n};
    auto array = create_ndarray<M, 1>(shape);
    auto raw = array.data();
    auto n_byte = array.itemsize();
    uint8_t * ptr = reinterpret_cast<uint8_t*>(raw);

    // Fill the array
    for (size_t bit = 0; bit < ivals.size(); ++bit) {
        for (auto & p: ivals[bit].segments) {
            for (
                int i = p.first - domain.first;
                i < p.second - domain.first;
                i++
            ) {
                ptr[i*n_byte + bit/8] |= (1<<(bit%8));
            }
        }
    }
    return array.cast();
}

template <
    typename C,
    typename M,
    typename std::enable_if_t<! std::is_integral<C>::value, bool>* = nullptr
>
nb::object intervals_to_mask(nb::list ivlist, int n_bits) {
    throw general_agreement_exception("mask() is only valid for integral types");
    return nb::object();
}


// Helper function to register an Intervals class for a concrete type.

template <typename C>
void intervals_bindings(nb::module_ & m, char const * name) {

    nb::class_<Intervals<C>>(m, name)
        .def(nb::init<C, C>())
        .def("__str__", &Intervals<C>::Description)
        .def("add_interval", &Intervals<C>::add_interval, nb::rv_policy::none)
        .def(
            "append_interval_no_check",
            &Intervals<C>::append_interval_no_check,
            nb::rv_policy::none
        )
        .def("merge", &Intervals<C>::merge, nb::rv_policy::none)
        .def("intersect", &Intervals<C>::intersect, nb::rv_policy::none)
        .def_prop_rw("domain",
            [](Intervals<C> & slf) {
                return nb::make_tuple(slf.get_domain());
            },
            [](Intervals<C> & slf, nb::object value) {
                if (nb::isinstance<nb::list>(value)) {
                    auto v = nb::cast<nb::list>(value);
                    if (v.size() != 2) {
                        throw shape_exception("domain", "!= 2");
                    }
                    slf.set_domain(nb::cast<C>(v[0]), nb::cast<C>(v[1]));
                } else if (nb::isinstance<nb::tuple>(value)) {
                    auto v = nb::cast<nb::tuple>(value);
                    if (v.size() != 2) {
                        throw shape_exception("domain", "!= 2");
                    }
                    slf.set_domain(nb::cast<C>(v[0]), nb::cast<C>(v[1]));
                } else {
                    throw general_agreement_exception(
                        "Only list or tuple values can be used to set domain"
                    );
                }
            }
        )
        .def("complement", &Intervals<C>::complement, nb::rv_policy::take_ownership)
        .def(
            "copy",
            [](Intervals<C> & slf) {return Intervals<C>(slf);},
            nb::rv_policy::take_ownership
        )
        .def_static("from_array", [](nb::ndarray<C, nb::ndim<2>, nb::c_contig> input) {
            Intervals<C> output;
            // Fast-access view of the array
            auto v = input.view();
            size_t n_seg = v.shape(0);
            if (v.shape(1) != 2) {
                throw shape_exception("input array shape[1]", "!= 2");
            }
            auto data = v.data();
            for (size_t i = 0; i < n_seg; ++i) {
                output.add_interval(data[2*i], data[2*i + 1]);
            }
            return output;
        }, nb::rv_policy::take_ownership)
        .def("array", [](Intervals<C> & slf) {
            auto n_seg = slf.segments.size();
            size_t shape[2] = {n_seg, 2};
            auto array = create_ndarray<C, 2>(shape);
            auto raw = array.data();

            // Fill the array
            size_t seg = 0;
            for (auto const & p : slf.segments) {
                raw[2 * seg] = p.first;
                raw[2 * seg + 1] = p.second;
                seg++;
            }
            return array;
        })
        .def("__getitem__", [](Intervals<C> & slf, nb::slice slc) {
            if (! std::is_integral<C>::value) {
                throw ValueError_exception("Intervals __getitem__: type not integral");
            }
            auto count = slf.domain.second - slf.domain.first;
            auto slc_par = slc.compute(count);
            auto start = slc_par.template get<0>();
            auto stop = slc_par.template get<1>();
            auto step = slc_par.template get<2>();
            auto slen = slc_par.template get<3>();

            // Now interpret start, stop in domain reference.
            start += slf.domain.first;
            stop += slf.domain.first;

            if (start > slf.domain.second) {
                return Intervals<C>(slf.domain.second, slf.domain.second);
            }
            if (stop < slf.domain.first) {
                return Intervals<C>(slf.domain.first, slf.domain.first);
            }

            if (start < slf.domain.first) {
                start = slf.domain.first;
            }
            if (stop > slf.domain.second) {
                stop = slf.domain.second;
            }

            auto output = slf; // copy
            output.set_domain(start, stop);
            return output;
        }, nb::rv_policy::take_ownership)
        .def_static("from_mask", [](
            nb::ndarray<nb::ndim<1>, nb::c_contig> input, int n_bits
        ) {
            // Do runtime creation of a view with a specific type, and then
            // dispatch that to the proper function.  See:
            // https://nanobind.readthedocs.io/en/latest/ndarray.html#specializing-views-at-runtime

            if (input.dtype() == nb::dtype<uint8_t>()) {
                auto view = input.view<uint8_t>();
                return intervals_from_mask<C, uint8_t>(view, n_bits);
            } else if (input.dtype() == nb::dtype<uint16_t>()) {
                auto view = input.view<uint16_t>();
                return intervals_from_mask<C, uint16_t>(view, n_bits);
            } else if (input.dtype() == nb::dtype<uint32_t>()) {
                auto view = input.view<uint32_t>();
                return intervals_from_mask<C, uint32_t>(view, n_bits);
            } else if (input.dtype() == nb::dtype<uint64_t>()) {
                auto view = input.view<uint64_t>();
                return intervals_from_mask<C, uint64_t>(view, n_bits);
            } else {
                std::ostringstream err;
                err << "Mask array type is not one of the supported types"
                    << " (uint8, uint16, uint32, uint64)";
                throw general_agreement_exception(err.str().c_str());
            }
        }, nb::rv_policy::take_ownership)
        .def_static("mask", [](nb::list ivlist, int n_bits){
            if (n_bits < 0) {
                // Use the length of the interval list
                n_bits = ivlist.size();
            }
            if (n_bits <= 8) {
                return intervals_to_mask<C, uint8_t>(ivlist, n_bits);
            } else if (n_bits <= 16) {
                return intervals_to_mask<C, uint16_t>(ivlist, n_bits);
            } else if (n_bits <= 32) {
                return intervals_to_mask<C, uint32_t>(ivlist, n_bits);
            } else if (n_bits <= 64) {
                return intervals_to_mask<C, uint64_t>(ivlist, n_bits);
            } else {
                throw general_agreement_exception(
                    "Number of requested bits must be <= 64"
                );
            }
        }, nb::rv_policy::take_ownership)
        .def(-nb::self)
        .def(~nb::self)
        .def(nb::self += nb::self)
        .def(nb::self -= nb::self)
        .def(nb::self + nb::self)
        .def(nb::self - nb::self)
        .def(nb::self * nb::self);

    return;
}


void register_intervals(nb::module_ & m) {
    // Concrete intervals types
    intervals_bindings<double>(m, "IntervalsDouble");
    intervals_bindings<int64_t>(m, "IntervalsInt");
    intervals_bindings<int32_t>(m, "IntervalsInt32");
    intervals_bindings<int64_t>(m, "IntervalsTime");
    return;
}



