#define NO_IMPORT_ARRAY

#include <iostream>
#include <limits>
#include <type_traits>

#include <nanobind/nanobind.h>
#include <nanobind/operators.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/string.h>
#include <nanobind/ndarray.h>

#include "so3g_numpy.h"

#include "Intervals.h"
#include "exceptions.h"
#include <exception>

namespace nb = nanobind;


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

template <typename T>
Intervals<T>::Intervals(Intervals<T> const & other) {
    std::cerr << "Copy constructor input = " << other.Description() << std::endl;
    domain = other.domain;
    segments = other.segments;
    std::cerr << "Copy constructor output = " << Description() << std::endl;
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
    ostringstream dbg;
    dbg << "DBG add_interval (" << start << "," << end << ")";
    std::cerr << dbg.str() << std::endl;
    // We can optimize this later.  For now just do something that is
    // obviously correct.
    auto p = lower_bound(segments.begin(), segments.end(), make_pair(start, end));
    segments.insert(p, make_pair(start, end));
    cleanup();

    dbg.str("");
    dbg << "DBG segments size now " << segments.size();
    std::cerr << dbg.str() << std::endl;

    return *this;
}

template <typename T>
Intervals<T>& Intervals<T>::append_interval_no_check(const T start, const T end)
{
    ostringstream dbg;
    dbg << "DBG add_interval_no_check (" << start << "," << end << ")";
    std::cerr << dbg.str() << std::endl;
    segments.push_back(make_pair(start, end));

    dbg.str("");
    dbg << "DBG segments size now " << segments.size();
    std::cerr << dbg.str() << std::endl;

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


template <typename T>
static inline
int interval_extract(const std::pair<T,T> *src, char *dest) {
    auto Tdest = reinterpret_cast<T*>(dest);
    *(Tdest) = src->first;
    *(Tdest+1) = src->second;
    return 2 * sizeof(*Tdest);
}

template <typename T>
static inline int get_dtype() {
    return NPY_NOTYPE;
}

template <>
inline int get_dtype<std::int64_t>() {
    return NPY_INT64;
}

template <>
inline int get_dtype<std::int32_t>() {
    return NPY_INT32;
}

template <>
inline int get_dtype<double>() {
    return NPY_FLOAT64;
}

template <typename T>
static int format_to_dtype(const BufferWrapper<T> &view)
{
    if (strcmp(view->format, "b") == 0 ||
        strcmp(view->format, "h") == 0 ||
        strcmp(view->format, "i") == 0 ||
        strcmp(view->format, "l") == 0 ||
        strcmp(view->format, "q") == 0) {
        switch(view->itemsize) {
        case 1:
            return NPY_INT8;
        case 2:
            return NPY_INT16;
        case 4:
            return NPY_INT32;
        case 8:
            return NPY_INT64;
        }
    } else if (strcmp(view->format, "c") == 0 ||
               strcmp(view->format, "B") == 0 ||
               strcmp(view->format, "H") == 0 ||
               strcmp(view->format, "I") == 0 ||
               strcmp(view->format, "L") == 0 ||
               strcmp(view->format, "Q") == 0) {
        switch(view->itemsize) {
        case 1:
            return NPY_UINT8;
        case 2:
            return NPY_UINT16;
        case 4:
            return NPY_UINT32;
        case 8:
            return NPY_UINT64;
        }
    } else if (strcmp(view->format, "f") == 0 ||
               strcmp(view->format, "d") == 0) {
        switch(view->itemsize) {
        case 4:
            return NPY_FLOAT32;
        case 8:
            return NPY_FLOAT64;
        }
    }

    return NPY_NOTYPE;
}

template <typename T>
Intervals<T> * Intervals<T>::from_array(const nb::object & src)
{
    Intervals<T> * output = new Intervals<T>();

    BufferWrapper<T> buf("src", src, false, vector<int>{-1, 2});
    char *d = (char*)buf->buf;
    size_t n_seg = buf->shape[0];

    std::cerr << "Intervals from_array shape = " << buf->shape[0] << "," << buf->shape[1] << std::endl;
    std::cerr << "Intervals from_array strides = " << buf->strides[0] << "," << buf->strides[1] << std::endl;

    for (size_t i = 0; i < n_seg; ++i) {
        std::pair<T,T> pr = interval_pair<T>(d, d+(buf->strides[1]));
        std::cerr << "Intervals from_array seg " << i << " pair = " << pr.first << "," << pr.second << std::endl;
        output->segments.push_back(pr);
        //output.segments.push_back(interval_pair<T>(d, d+buf->strides[1]));
        d += buf->strides[0];
        std::cerr << "Intervals from_array segments size now " << output->segments.size() << std::endl;
    }

    return output;
}

template <typename T>
nb::object Intervals<T>::array() const
{
    npy_intp dims[2];
    std::cerr << "Intervals segments.size() = " << segments.size() << std::endl;
    dims[0] = (npy_intp)segments.size();
    dims[1] = 2;
    int dtype = get_dtype<T>();
    if (dtype == NPY_NOTYPE)
        throw general_agreement_exception("array() not implemented for this domain dtype.");
    PyObject *v = PyArray_SimpleNew(2, dims, dtype);
    if (v == NULL) {
        ostringstream dstr;
        dstr << "Failed to allocate Intervals numpy array of size (";
        dstr << dims[0] << ", " << dims[1] << ")";
        throw RuntimeError_exception(dstr.str().c_str());
    }
    char *ptr = reinterpret_cast<char*>((PyArray_DATA((PyArrayObject*)v)));
    for (auto p = segments.begin(); p != segments.end(); ++p) {
        ptr += interval_extract((&*p), ptr);
    }
    return nb::steal<nb::object>(v);
}


/**
 * Bit-mask conversion - convert between list<IntervalsInt> and
 * ndarray bit-masks.
 *
 * Uses SFINAE of templates in order to have separate implementations
 * for integer and non-integer types.  (Templating the class member
 * function directly is a lot more difficult.)  The .from_mask(...)
 * call is implemented in the inline function from_mask_.  The
 * .mask(...) call is implemented in the inline function mask_.
 */

// from_mask_() definitions, for .from_mask().
//
// intType is the type of the Interval, which should be a simple
// signed integer type (e.g. int64_t).  numpyType is a simple unsigned
// type carried in the ndarray (e.g. uint8_t).  The n_bits argument is
// used to specify how many of the LSBs should be processed; if
// negative this will default to match the numpyType (e.g. 8 for
// uint8_t).

template <typename intType, typename numpyType,
          typename std::enable_if<!std::is_integral<intType>::value,
                                  int>::type* = nullptr>
static inline nb::object from_mask_(void *buf, intType count, int n_bits)
{
    throw dtype_exception("target", "Interval<> over integral type.");
    return nb::object();
}

template <typename intType, typename numpyType,
          typename std::enable_if<std::is_integral<intType>::value,
                                  int>::type* = nullptr>
static inline nb::object from_mask_(void *buf, intType count, int n_bits)
{
    if (n_bits < 0)
        n_bits = 8*sizeof(numpyType);

    auto p = (numpyType*)buf;

    vector<Intervals<intType>> output;;
    vector<intType> start;
    for (int bit=0; bit<n_bits; ++bit) {
        output.push_back(Intervals<intType>(0, count));
        start.push_back(-1);
    }

    numpyType last = 0;
    for (intType i=0; i<count; i++) {
        numpyType d = p[i] ^ last;
        for (int bit=0; bit<n_bits; ++bit) {
            if (d & (1 << bit)) {
                if (start[bit] >= 0) {
                    output[bit].segments.push_back(
                        interval_pair<intType>((char*)&start[bit], (char*)&i));
                    start[bit] = -1;
                } else {
                    start[bit] = i;
                }
            }
        }
        last = p[i];
    }
    for (int bit=0; bit<n_bits; ++bit) {
        if (start[bit] >= 0)
            output[bit].segments.push_back(
                interval_pair<intType>((char*)&start[bit], (char*)&count));
    }

    // Once added to the list, we can't modify further.
    nb::list bits;
    for (auto i: output)
        bits.append(i);
    return bits;
}

template <typename T>
nb::object Intervals<T>::from_mask(const nb::object &src, int n_bits)
{
    BufferWrapper<T> buf("src", src, false);

    if (buf->ndim != 1)
        throw shape_exception("src", "must be 1-d");

    int p_count = buf->shape[0];
    void *p = buf->buf;

    int dtype = format_to_dtype(buf);
    switch(dtype) {
    case NPY_UINT8:
    case NPY_INT8:
        return from_mask_<T,uint8_t>(p, p_count, n_bits);
    case NPY_UINT16:
    case NPY_INT16:
        return from_mask_<T,uint16_t>(p, p_count, n_bits);
    case NPY_UINT32:
    case NPY_INT32:
        return from_mask_<T,uint32_t>(p, p_count, n_bits);
    case NPY_UINT64:
    case NPY_INT64:
        return from_mask_<T,uint64_t>(p, p_count, n_bits);
    }

    throw dtype_exception("src", "integer type");
    return nb::object();
}


// mask_() definitions, for .mask().
//
// intType is the type of the Interval, which should be a simple
// signed integer type (e.g. int64_t).  The numpy array data type is
// determined based on the number of bits requested, either explicitly
// through the n_bits argument (which must be large enough to handle
// the list) or implicitly through the length of the list of
// Interval<T> objects.

template <typename intType,typename std::enable_if<!std::is_integral<intType>::value,
                                                   int>::type* = nullptr>
static inline nb::object mask_(const nb::list &ivlist, int n_bits)
{
    intType x;
    throw dtype_exception("ivlist", "Interval<> over integral type.");
    return nb::object();
}

template <typename intType, typename std::enable_if<std::is_integral<intType>::value,
                                                    int>::type* = nullptr>
static inline nb::object mask_(const nb::list &ivlist, int n_bits)
{
    vector<Intervals<intType>> ivals;
    vector<int> indexes;

    pair<intType,intType> domain;

    for (long i=0; i<nb::len(ivlist); i++) {
        std::cerr << "Intervals mask processing list item " << i << std::endl;
        indexes.push_back(0);
        ivals.push_back(nb::cast<Intervals<intType>>(ivlist[i]));
        if (i==0) {
            domain = ivals[i].domain;
        } else if (domain != ivals[i].domain) {
            throw agreement_exception("ivlist[0]", "all other ivlist[i]", "domain");
        }
    }
    std::cerr << "Intervals mask done processing list" << std::endl;

    // Determine the output mask size based on n_bits... which may be unspecified.
    int npy_type = NPY_UINT8;
    if (n_bits < 0)
        n_bits = ivals.size();
    else if (n_bits < ivals.size())
        throw general_agreement_exception("Input list has more items than the "
                                          "output mask size (n_bits).");
    std::cerr << "Intervals mask using n_bits = " << n_bits << std::endl;

    if (n_bits <= 8)
        npy_type = NPY_UINT8;
    else if (n_bits <= 16)
        npy_type = NPY_UINT16;
    else if (n_bits <= 32)
        npy_type = NPY_UINT32;
    else if (n_bits <= 64)
        npy_type = NPY_UINT64;
    else {
	std::ostringstream err;
        err << "No integer type is available to host the " << n_bits
            << " requested to encode this mask.";
        throw general_agreement_exception(err.str());
    }
    std::cerr << "Intervals mask use npy_type = " << npy_type << std::endl;

    int n = domain.second - domain.first;
    npy_intp dims[1];
    dims[0] = n;

    std::cerr << "Intervals mask allocate 1D array of len " << n << std::endl;
    PyObject *v = PyArray_SimpleNew(1, dims, npy_type);
    if (v == NULL) {
        ostringstream dstr;
        dstr << "Failed to allocate Intervals mask array of size (";
        dstr << dims[0] << ",)";
        throw RuntimeError_exception(dstr.str().c_str());
    }

    // Assumes little-endian.
    int n_byte = PyArray_ITEMSIZE((PyArrayObject*)v);
    std::cerr << "Intervals mask n_byte = " << n_byte << std::endl;

    uint8_t *ptr = reinterpret_cast<uint8_t*>((PyArray_DATA((PyArrayObject*)v)));
    memset(ptr, 0, n*n_byte);
    for (long bit=0; bit<ivals.size(); ++bit) {
        for (auto p: ivals[bit].segments) {
            for (int i=p.first - domain.first; i<p.second - domain.first; i++) {
                std::cerr << "Intervals mask bit " << bit << ", i " << i << ": ptr[" << i*n_byte + bit/8 << "] |= " << (1<<(bit%8)) << std::endl;
                ptr[i*n_byte + bit/8] |= (1<<(bit%8));
            }
        }
    }

    std::cerr << "Intervals mask return" << std::endl;
    return nb::steal<nb::object>(v);
}

template <typename T>
nb::object Intervals<T>::mask(const nb::list &ivlist, int n_bits)
{
    return mask_<T>(ivlist, n_bits);
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


// Slicing!
template <typename T,
          typename std::enable_if<!std::is_integral<T>::value,
                                  int>::type* = nullptr>
static inline Intervals<T> _getitem_(Intervals<T> &src, nb::object indices)
{
    throw dtype_exception("target", "Interval<> over integral type.");
    return Intervals<T>();
}

template <typename objType, typename T>
static inline T extract_or_default(objType src, T default_)
{
    T result;
    if (nb::try_cast<T>(src, result)) {
        // Successful cast
        return result;
    } else {
        return default_;
    }
}

template <typename T,
          typename std::enable_if<std::is_integral<T>::value,
                                  int>::type* = nullptr>
static inline Intervals<T> _getitem_(Intervals<T> &src, nb::object indices)
{
    if (nb::isinstance<nb::slice>(indices)) {
        nb::slice sl = nb::cast<nb::slice>(indices);

        T count = src.domain.second - src.domain.first;
        auto slc_par = sl.compute(count);

        T start = slc_par.template get<0>();
        T stop = slc_par.template get<1>();
        T step = slc_par.template get<2>();

        assert(step == 1);
        if (start < 0)
            start = count + start;
        if (stop < 0)
            stop = count + stop;
        if (stop < start)
            stop = start;

        // Now interpret start, stop in domain reference.
        start += src.domain.first;
        stop += src.domain.first;

        if (start > src.domain.second)
            return Intervals<T>(src.domain.second, src.domain.second);
        if (stop < src.domain.first)
            return Intervals<T>(src.domain.first, src.domain.first);

        if (start < src.domain.first)
            start = src.domain.first;
        if (stop > src.domain.second)
            stop = src.domain.second;

        auto output = src; // copy
        output.set_domain(start, stop);
        return output;
    }
    return Intervals<T>();
}

template <typename T>
Intervals<T> Intervals<T>::getitem(nb::object indices)
{
    return _getitem_(*this, indices);
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


// Helper function to register an Intervals class for a concrete type.

template <typename C>
void intervals_bindings(nb::module_ & m, char const * name) {

    nb::class_<Intervals<C>>(m, name)
        .def(nb::init<>())
        .def(nb::init<C, C>(),
            R"(
            A finite series of non-overlapping semi-open intervals
            )"
        )
        .def("__str__", &Intervals<C>::Description)
        .def("add_interval", &Intervals<C>::add_interval, nb::rv_policy::none,
            nb::arg("start"),
            nb::arg("end"),
            R"(
            Merge an interval into the set.
            )"
        )
        .def("append_interval_no_check", &Intervals<C>::append_interval_no_check,
            nb::rv_policy::none,
            nb::arg("start"),
            nb::arg("end"),
            R"(
            Append an interval to the set without checking for overlap or sequence.
            )"
        )
        .def("merge", &Intervals<C>::merge, nb::rv_policy::none,
            R"(
            Merge an Intervals into the set.
            )"
        )
        .def("intersect", &Intervals<C>::intersect, nb::rv_policy::none,
            nb::arg("source"),
            R"(
            Intersect another Intervals object with this one.
            )"
        )
        .def_prop_rw("domain",
            [](Intervals<C> & slf) {
                auto dom = slf.get_domain();
                return nb::make_tuple(dom.first, dom.second);
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
            },
            R"(
            Interval set domain (settable, with consequences).
            )"
        )
        .def("complement", &Intervals<C>::complement, nb::rv_policy::take_ownership,
            R"(
            Return the complement (over domain).
            )"
        )
        .def(
            "copy",
            [](Intervals<C> & slf) {
                auto obj = Intervals<C>(slf);
                std::cerr << "DBG bindings: " << obj.Description() << std::endl;
                return obj;
            }, nb::rv_policy::move,
            R"(
            Get a new object with a copy of the data.
            )"
        )
        .def_static("from_array", &Intervals<C>::from_array,
            nb::rv_policy::take_ownership,
            nb::arg("input_array"),
            R"(
            Return an Intervals object based on an (n,2) ndarray.
            )"
        )
        .def("array", &Intervals<C>::array, nb::rv_policy::take_ownership,
            R"(
            Return the intervals as a 2-d numpy array.
            )"
        )
        .def("__getitem__", &Intervals<C>::getitem, nb::rv_policy::take_ownership)
        .def_static("from_mask", &Intervals<C>::from_mask,
            nb::rv_policy::take_ownership,
            nb::arg("input_array"),
            nb::arg("n_bits"),
            R"(
            Return a list Intervals.

            The Intervals are extracted from the first n_bits of the input_array
            (a 1-D array of integral type).
            )"
        )
        .def_static("mask", &Intervals<C>::mask, nb::rv_policy::take_ownership,
            nb::arg("intervals_list"),
            nb::arg("n_bits"),
            R"(
            Return an ndarray bitmask from a list of Intervals.

            The dtype will be the smallest available to hold n_bits.
            )"
        )
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
    return;
}
