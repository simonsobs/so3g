#define NO_IMPORT_ARRAY

#include <pybindings.h>

#include <iostream>
#include <limits>

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
inline
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
inline
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
inline int get_dtype() {
    return NPY_NOTYPE;
}

template <>
inline int get_dtype<std::int64_t>() {
    return NPY_INT64;
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
    PyObject *v = PyArray_SimpleNew(2, dims, dtype);
    char *ptr = reinterpret_cast<char*>((PyArray_DATA((PyArrayObject*)v)));
    for (auto p = segments.begin(); p != segments.end(); ++p) {
        ptr += interval_extract((&*p), ptr);
    }
    return bp::object(bp::handle<>(v));
}


//
// Bit-mask conversion - create list<IntervalsInt> from ndarray bit-masks.
//

template <typename T>
bp::object Intervals<T>::from_mask(const bp::object &src, int n_bits)
{
     throw "Not implemented for non-index Intervals types.";
     bp::list bits;
     Intervals<T> output;
     bits.append(output);
     return bits;
}

template <>
bp::object Intervals<int64_t>::from_mask(const bp::object &src, int n_bits)
{
    // Get a view...
    BufferWrapper buf;
    if (PyObject_GetBuffer(src.ptr(), &buf.view,
                           PyBUF_FORMAT | PyBUF_ANY_CONTIGUOUS) == -1) {
        PyErr_Clear();
        throw buffer_exception("src");
    }

    if (buf.view.ndim != 1)
        throw shape_exception("src", "must be 1-d");

    int dtype = format_to_dtype(buf.view);
    if (dtype != NPY_UINT8)
        throw dtype_exception("src", "uint8");

    int64_t n = buf.view.shape[0];

    if (n_bits < 0)
        n_bits = 8*sizeof(uint8_t);

    vector<Intervals<int64_t>> output;;
    vector<int64_t> start;
    for (int i=0; i<n_bits; i++) {
        output.push_back(Intervals<int64_t>(0, n));
        start.push_back(-1);
    }

    uint8_t *p = (uint8_t*)buf.view.buf;

    uint8_t last = 0;
    for (int64_t i=0; i<n; i++) {
        uint8_t d = p[i] ^ last;
        for (int bit=0; bit<n_bits; ++bit) {
            if (d & (1 << bit)) {
                if (start[bit] >= 0) {
                    output[bit].segments.push_back(
                        interval_pair<int64_t>((char*)&start[bit], (char*)&i));
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
                interval_pair<int64_t>((char*)&start[bit], (char*)&n));
    }

    // Once added to the list, we can't modify further.
    bp::list bits;
    for (auto i: output)
        bits.append(i);
    return bits;
}


// Convert a list of Intervals<T> to a ndarray bitmask vector.

template <typename T>
bp::object Intervals<T>::mask(const bp::list &ivlist, int n_bits)
{
    npy_intp dims[1] = {0};
    int dtype = get_dtype<T>();
    PyObject *v = PyArray_SimpleNew(1, dims, NPY_UINT8);
    return bp::object(bp::handle<>(v));
}

template <>
bp::object Intervals<int64_t>::mask(const bp::list &ivlist, int n_bits)
{
    assert(n_bits == 8);
    vector<Intervals<int64_t>> ivals;
    vector<int> indexes;

    pair<int64_t,int64_t> domain;

    for (long i=0; i<bp::len(ivlist); i++) {
        indexes.push_back(0);
        ivals.push_back(bp::extract<Intervals<int64_t>>(ivlist[i]));
        if (i==0) {
            domain = ivals[i].domain;
        } else if (domain != ivals[i].domain) {
            throw agreement_exception("ivlist[0]", "all other ivlist[i]", "domain");
        }
    }

    int n = domain.second - domain.first;
    npy_intp dims[1] = {n};
    PyObject *v = PyArray_SimpleNew(1, dims, NPY_UINT8);

    uint8_t *ptr = reinterpret_cast<uint8_t*>((PyArray_DATA((PyArrayObject*)v)));
    memset(ptr, 0, n);
    for (long bit=0; bit<ivals.size(); ++bit) {
        for (auto p: ivals[bit].segments) {
            for (int i=p.first - domain.first; i<p.second - domain.first; i++)
                ptr[i] |= (1<<bit);
        }
    }

    return bp::object(bp::handle<>(v));
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
    .def("mask", &CLASSNAME::mask,                     \
         "Return an ndarray bitmask from a list of " #DOMAIN_TYPE ".") \
    .staticmethod("mask")                                        \
    .def("from_mask", &CLASSNAME::from_mask,                     \
         "Return a " #DOMAIN_TYPE " based on an ndarray of uint8.") \
    .staticmethod("from_mask")                                   \
    .def("copy", \
         +[](CLASSNAME& A) { \
              return CLASSNAME(A); \
          }, \
         "Get a new object with a copy of the data.") \
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
G3_SERIALIZABLE_CODE(IntervalsTime);

G3_SERIALIZABLE_CODE(MapIntervalsDouble);
G3_SERIALIZABLE_CODE(MapIntervalsInt);
G3_SERIALIZABLE_CODE(MapIntervalsTime);

PYBINDINGS("so3g")
{
    EXPORT_INTERVALS(double,  IntervalsDouble);
    EXPORT_INTERVALS(int64_t, IntervalsInt);
    EXPORT_INTERVALS(G3Time,  IntervalsTime);
}
