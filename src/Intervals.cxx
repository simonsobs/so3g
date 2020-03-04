#define NO_IMPORT_ARRAY

#include <pybindings.h>

#include <iostream>
#include <limits>
#include <type_traits>

#include <boost/python.hpp>
#include <cereal/types/utility.hpp>
#include <container_pybindings.h>

#include "so3g_numpy.h"

#include "Intervals.h"
#include "exceptions.h"
#include <exception>

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

// The G3Time internal encoding is an int64 with the number of 100 MHz
// ticks since the unix epoch.  Make our default domain span from a
// while ago to a while from now.

#define G3TIME_LO                   0LL  // Jan 1 1970
#define G3TIME_HI  725811840000000000LL  // Jan 1 2200

template <>
Intervals<G3Time>::Intervals() {
    domain = make_pair(G3Time(G3TIME_LO), G3Time(G3TIME_HI));
}


template <typename T>
std::string Intervals<T>::Description() const
{
	std::ostringstream s;
	s << "Intervals over domain [" << domain.first << "," << domain.second << ")";
	return s.str();
}

template <typename T>
template <class A> void Intervals<T>::serialize(A &ar, unsigned v)
{
	using namespace cereal;
        // v is the version code!

	ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));
	ar & make_nvp("domain", domain);
	ar & make_nvp("segments", segments);
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

template <>
inline
pair<G3Time,G3Time> interval_pair(char *p1, char *p2) {
    return make_pair(G3Time(*reinterpret_cast<G3TimeStamp*>(p1)),
                     G3Time(*reinterpret_cast<G3TimeStamp*>(p2)));
}


template <typename T>
static inline
int interval_extract(const std::pair<T,T> *src, char *dest) {
    auto Tdest = reinterpret_cast<T*>(dest);
    *(Tdest) = src->first;
    *(Tdest+1) = src->second;
    return 2 * sizeof(*Tdest);
}

template <>
inline
int interval_extract(const std::pair<G3Time,G3Time> *src, char *dest) {
    auto Tdest = reinterpret_cast<G3TimeStamp*>(dest);
    *(Tdest) = src->first.time;
    *(Tdest+1) = src->second.time;
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

template <>
inline int get_dtype<G3Time>() {
    return NPY_INT64;
}

static int format_to_dtype(const Py_buffer &view)
{
    if (strcmp(view.format, "b") == 0 ||
        strcmp(view.format, "h") == 0 ||
        strcmp(view.format, "i") == 0 ||
        strcmp(view.format, "l") == 0 ||
        strcmp(view.format, "q") == 0) {
        switch(view.itemsize) {
        case 1:
            return NPY_INT8;
        case 2:
            return NPY_INT16;
        case 4:
            return NPY_INT32;
        case 8:
            return NPY_INT64;
        }
    } else if (strcmp(view.format, "c") == 0 ||
               strcmp(view.format, "B") == 0 ||
               strcmp(view.format, "H") == 0 ||
               strcmp(view.format, "I") == 0 ||
               strcmp(view.format, "L") == 0 ||
               strcmp(view.format, "Q") == 0) {
        switch(view.itemsize) {
        case 1:
            return NPY_UINT8;
        case 2:
            return NPY_UINT16;
        case 4:
            return NPY_UINT32;
        case 8:
            return NPY_UINT64;
        }
    } else if (strcmp(view.format, "f") == 0 ||
               strcmp(view.format, "d") == 0) {
        switch(view.itemsize) {
        case 4:
            return NPY_FLOAT32;
        case 8:
            return NPY_FLOAT64;
        }
    } 

    return NPY_NOTYPE;
}


template <typename T>
Intervals<T> Intervals<T>::from_array(const bp::object &src)
{
    Intervals<T> output;

    // Get a view...
    BufferWrapper buf;
    if (PyObject_GetBuffer(src.ptr(), &buf.view,
                           PyBUF_FORMAT | PyBUF_ANY_CONTIGUOUS) == -1) {
        PyErr_Clear();
        throw buffer_exception("src");
    } 

    if (buf.view.ndim != 2 || buf.view.shape[1] != 2)
        throw shape_exception("src", "must have shape (n,2)");

    int dtype = format_to_dtype(buf.view);
    if (dtype != get_dtype<T>())
        throw dtype_exception("src", "matching Interval class");

    char *d = (char*)buf.view.buf;
    int n_seg = buf.view.shape[0];
    for (int i=0; i<n_seg; ++i) {
        output.segments.push_back(interval_pair<T>(d, d+buf.view.strides[1]));
        d += buf.view.strides[0];
    }
    
    return output;
}

template <typename T>
bp::object Intervals<T>::array() const
{
    npy_intp dims[2] = {0, 2};
    dims[0] = (npy_intp)segments.size();
    int dtype = get_dtype<T>();
    if (dtype == NPY_NOTYPE)
        throw general_agreement_exception("array() not implemented for this domain dtype.");

    PyObject *v = PyArray_SimpleNew(2, dims, dtype);
    char *ptr = reinterpret_cast<char*>((PyArray_DATA((PyArrayObject*)v)));
    for (auto p = segments.begin(); p != segments.end(); ++p) {
        ptr += interval_extract((&*p), ptr);
    }
    return bp::object(bp::handle<>(v));
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
static inline bp::object from_mask_(void *buf, intType count, int n_bits)
{
    throw dtype_exception("target", "Interval<> over integral type.");
    return bp::object();
}

template <typename intType, typename numpyType,
          typename std::enable_if<std::is_integral<intType>::value,
                                  int>::type* = nullptr>
static inline bp::object from_mask_(void *buf, intType count, int n_bits)
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
    bp::list bits;
    for (auto i: output)
        bits.append(i);
    return bits;
}

template <typename T>
bp::object Intervals<T>::from_mask(const bp::object &src, int n_bits)
{
    BufferWrapper buf;
    if (PyObject_GetBuffer(src.ptr(), &buf.view,
                           PyBUF_FORMAT | PyBUF_ANY_CONTIGUOUS) == -1) {
        PyErr_Clear();
        throw buffer_exception("src");
    }

    if (buf.view.ndim != 1)
        throw shape_exception("src", "must be 1-d");

    int p_count = buf.view.shape[0];
    void *p = buf.view.buf;

    int dtype = format_to_dtype(buf.view);
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
    return bp::object();
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
static inline bp::object mask_(const bp::list &ivlist, int n_bits)
{
    intType x;
    throw dtype_exception("ivlist", "Interval<> over integral type.");
    return bp::object();
}

template <typename intType, typename std::enable_if<std::is_integral<intType>::value,
                                                    int>::type* = nullptr>
static inline bp::object mask_(const bp::list &ivlist, int n_bits)
{
    vector<Intervals<intType>> ivals;
    vector<int> indexes;

    pair<intType,intType> domain;

    for (long i=0; i<bp::len(ivlist); i++) {
        indexes.push_back(0);
        ivals.push_back(bp::extract<Intervals<intType>>(ivlist[i]));
        if (i==0) {
            domain = ivals[i].domain;
        } else if (domain != ivals[i].domain) {
            throw agreement_exception("ivlist[0]", "all other ivlist[i]", "domain");
        }
    }

    // Determine the output mask size based on n_bits... which may be unspecified.
    int npy_type = NPY_UINT8;
    if (n_bits < 0)
        n_bits = ivals.size();
    else if (n_bits < ivals.size())
        throw general_agreement_exception("Input list has more items than the "
                                          "output mask size (n_bits).");

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

    int n = domain.second - domain.first;
    npy_intp dims[1] = {n};
    PyObject *v = PyArray_SimpleNew(1, dims, npy_type);

    // Assumes little-endian.
    int n_byte = PyArray_ITEMSIZE((PyArrayObject*)v);
    uint8_t *ptr = reinterpret_cast<uint8_t*>((PyArray_DATA((PyArrayObject*)v)));
    memset(ptr, 0, n*n_byte);
    for (long bit=0; bit<ivals.size(); ++bit) {
        for (auto p: ivals[bit].segments) {
            for (int i=p.first - domain.first; i<p.second - domain.first; i++)
                ptr[i*n_byte + bit/8] |= (1<<(bit%8));
        }
    }

    return bp::object(bp::handle<>(v));
}

template <typename T>
bp::object Intervals<T>::mask(const bp::list &ivlist, int n_bits)
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
static inline Intervals<T> _getitem_(Intervals<T> &src, bp::object indices)
{
    throw dtype_exception("target", "Interval<> over integral type.");
    return Intervals<T>();
}

template <typename objType, typename T>
static inline T extract_or_default(objType src, T default_)
{
    bp::extract<T> ex(src);
    if (ex.check())
        return ex();
    return default_;
}

template <typename T,
          typename std::enable_if<std::is_integral<T>::value,
                                  int>::type* = nullptr>
static inline Intervals<T> _getitem_(Intervals<T> &src, bp::object indices)
{
    bp::extract<bp::slice> ex(indices);
    if (ex.check()) {
        T count = src.domain.second - src.domain.first;

        auto sl = ex();
        T start = extract_or_default(sl.start(), 0);
        T stop = extract_or_default(sl.stop(), count);
        T step = extract_or_default(sl.step(), 1);

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
Intervals<T> Intervals<T>::getitem(bp::object indices)
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

//
// boost-python registration.
//

using namespace boost::python;

#define EXPORT_INTERVALS(DOMAIN_TYPE, CLASSNAME) \
    EXPORT_FRAMEOBJECT(CLASSNAME, init<>(), \
   "A finite series of non-overlapping semi-open intervals on a domain " \
   "of type: " #DOMAIN_TYPE ".") \
    .def(init<const DOMAIN_TYPE&, const DOMAIN_TYPE&>("Initialize with domain.")) \
    .def("add_interval", &CLASSNAME::add_interval, \
         return_internal_reference<>(), \
         "Merge an interval into the set.") \
    .def("append_interval_no_check", &CLASSNAME::append_interval_no_check, \
         return_internal_reference<>(), \
         "Append an interval to the set without checking for overlap or sequence.") \
    .def("merge", &CLASSNAME::merge, \
         return_internal_reference<>(), \
         "Merge an Intervals into the set.") \
    .def("intersect", &CLASSNAME::intersect, \
         return_internal_reference<>(), \
         "Intersect another " #DOMAIN_TYPE "with this one.") \
    .add_property(                                                 \
        "domain",                                                  \
        +[](const CLASSNAME& A) {                                  \
             return make_tuple( A.domain.first, A.domain.second ); \
         },                                                        \
        +[](CLASSNAME& A, object _domain) {                        \
             A.set_domain(extract<DOMAIN_TYPE>(_domain[0]),        \
                          extract<DOMAIN_TYPE>(_domain[1]));       \
         },                                                        \
        "Interval set domain (settable, with consequences).")      \
    .def("complement", &CLASSNAME::complement, \
         "Return the complement (over domain).") \
    .def("array", &CLASSNAME::array, \
         "Return the intervals as a 2-d numpy array.") \
    .def("from_array", &CLASSNAME::from_array,              \
         "Return a " #DOMAIN_TYPE " based on an (n,2) ndarray.") \
    .staticmethod("from_array")                                  \
    .def("from_mask", &CLASSNAME::from_mask,                     \
         "Return a list of " #CLASSNAME " extracted from an ndarray encoding a bitmask.") \
    .staticmethod("from_mask")                                   \
    .def("mask", &CLASSNAME::mask,                                      \
         "Return an ndarray bitmask from a list of" #CLASSNAME ".") \
    .staticmethod("mask")                                               \
    .def("copy", \
         +[](CLASSNAME& A) { \
              return CLASSNAME(A); \
          }, \
         "Get a new object with a copy of the data.") \
    .def("__getitem__", &CLASSNAME::getitem) \
    .def(-self) \
    .def(~self) \
    .def(self += self) \
    .def(self -= self) \
    .def(self + self) \
    .def(self - self) \
    .def(self * self); \
    register_g3map<Map ## CLASSNAME>("Map" #CLASSNAME, "Mapping from "  \
        "strings to Intervals over " #DOMAIN_TYPE ".")


G3_SERIALIZABLE_CODE(IntervalsDouble);
G3_SERIALIZABLE_CODE(IntervalsInt);
G3_SERIALIZABLE_CODE(IntervalsInt32);
G3_SERIALIZABLE_CODE(IntervalsTime);

G3_SERIALIZABLE_CODE(MapIntervalsDouble);
G3_SERIALIZABLE_CODE(MapIntervalsInt);
G3_SERIALIZABLE_CODE(MapIntervalsInt32);
G3_SERIALIZABLE_CODE(MapIntervalsTime);

PYBINDINGS("so3g")
{
    EXPORT_INTERVALS(double,  IntervalsDouble);
    EXPORT_INTERVALS(int64_t, IntervalsInt);
    EXPORT_INTERVALS(int32_t, IntervalsInt32);
    EXPORT_INTERVALS(G3Time,  IntervalsTime);
}
