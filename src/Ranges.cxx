#define NO_IMPORT_ARRAY

#include <pybindings.h>

#include <iostream>
#include <limits>
#include <type_traits>

#include <boost/python.hpp>

#include "so3g_numpy.h"

#include "Ranges.h"
#include "exceptions.h"
#include <exception>


template <typename T>
std::string Ranges<T>::Description() const
{
	std::ostringstream s;
	s << "Ranges(n=" << count
          << ":rngs=" << segments.size() << ")";
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
void Ranges<T>::cleanup()
{
    auto p = segments.begin();
    while (p != segments.end()) {
        // Truncate intervals that cross the domain boundary.
        if (p->first < 0)
            p->first = 0;
        if (p->second > count)
            p->second = count;
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
Ranges<T>& Ranges<T>::_add_interval_numpysafe(
    const bp::object start_obj, const bp::object end_obj)
{
    int start = numpysafe_extract_int(start_obj, "start");
    int end = numpysafe_extract_int(end_obj, "end");
    return add_interval(start, end);
}

template <typename T>
Ranges<T>& Ranges<T>::add_interval(const T start, const T end)
{
    // We can optimize this later.  For now just do something that is
    // obviously correct.
    auto p = lower_bound(segments.begin(), segments.end(), make_pair(start, end));
    segments.insert(p, make_pair(start, end));
    cleanup();

    return *this;
}

template <typename T>
Ranges<T>& Ranges<T>::append_interval_no_check(const T start, const T end)
{
    segments.push_back(make_pair(start, end));

    return *this;
}

template <typename T>
Ranges<T>& Ranges<T>::merge(const Ranges<T> &src)
{
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

//Buffer Range in place
template <typename T>
Ranges<T>& Ranges<T>::buffer(const T buff)
{
    auto p0 = this->segments.begin();
    while (p0 != this->segments.end()) {
        p0->first -= buff;
        p0->second += buff;
        p0++;
    }
    cleanup();
    return *this;
}

//Return newly buffered range
template <typename T>
Ranges<T> Ranges<T>::buffered(const T buff) const
{
    Ranges<T> output(count, reference);

    for (auto p: segments) {
        output.segments.push_back(make_pair(p.first - buff, p.second + buff));
    }
    output.cleanup();
    return output;
}

//Close gaps between Ranges if they are leq gap
template <typename T>
Ranges<T>& Ranges<T>::close_gaps(const T gap)
{
    auto p = segments.begin();
    while (p != segments.end()) {
        // Check for distance from the front
        if (p->first <= gap){
            p->first = 0;
        }
        if (p->second >= count-gap) {
            p->second = count;
        }
        // Check for distances from the next interval.
        auto q = p+1;
        if (q == segments.end())
            break;
        // if distance is leq gap, close interval
        if (q->first - p->second <= gap) {
            p->second = q->second;
            segments.erase(q);
        } else{
            p++;
        }
    }
    
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
inline int get_dtype<int64_t>() {
    return NPY_INT64;
}

template <>
inline int get_dtype<int32_t>() {
    return NPY_INT32;
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
    } 

    return NPY_NOTYPE;
}


template <typename T>
Ranges<T> Ranges<T>::from_array(const bp::object &src, const bp::object &count)
{
    Ranges<T> output;
    BufferWrapper<T> buf("src", src, false, vector<int>{-1, 2});

    char *d = (char*)buf->buf;
    int n_seg = buf->shape[0];
    for (int i=0; i<n_seg; ++i) {
        output.segments.push_back(interval_pair<T>(d, d+buf->strides[1]));
        d += buf->strides[0];
    }
    output.count = numpysafe_extract_int(count, "count");

    output.cleanup();
    return output;
}

template <typename T>
bp::object Ranges<T>::ranges() const
{
    npy_intp dims[2] = {0, 2};
    dims[0] = (npy_intp)segments.size();
    int dtype = get_dtype<T>();
    if (dtype == NPY_NOTYPE)
        throw general_agreement_exception("ranges() not implemented for this domain dtype.");

    PyObject *v = PyArray_SimpleNew(2, dims, dtype);
    char *ptr = reinterpret_cast<char*>((PyArray_DATA((PyArrayObject*)v)));
    for (auto p = segments.begin(); p != segments.end(); ++p) {
        ptr += interval_extract((&*p), ptr);
    }
    return bp::object(bp::handle<>(v));
}


// 
// Bit-mask conversion - convert between list<RangesInt> and
// ndarray bit-masks.
// 
// intType is the type of the Interval, which should be a simple
// signed integer type (e.g. int64_t).  numpyType is a simple unsigned
// type carried in the ndarray (e.g. uint8_t).  The n_bits argument is
// used to specify how many of the LSBs should be processed; if
// negative this will default to match the numpyType (e.g. 8 for
// uint8_t).

template <typename intType, typename numpyType>
static inline bp::object from_bitmask_(void *buf, intType count, int n_bits)
{
    bool return_singleton = (n_bits == 0);

    if (n_bits < 0)
        n_bits = 8*sizeof(numpyType);
    else if (return_singleton)
        n_bits = 1;

    auto p = (numpyType*)buf;

    vector<Ranges<intType>> output;
    vector<intType> start;
    for (int bit=0; bit<n_bits; ++bit) {
        output.push_back(Ranges<intType>(count));
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

    if (return_singleton)
        return bp::object(output[0]);

    // Once added to the list, we can't modify further.
    bp::list bits;
    for (auto i: output)
        bits.append(i);
    return bits;
}

template <typename T>
bp::object Ranges<T>::from_bitmask(const bp::object &src, int n_bits)
{
    // Buffer protocol doesn't work directly on bool arrays, so if
    // what we've been passed is definitely a bool array, get a view
    // of it as a uint8 array and work with that.  (Wrap it with
    // bp::object so references are counted properly.)
    bp::object target(src);
    PyObject* obj_ptr = target.ptr();
    if (PyArray_Check(obj_ptr))
        if (PyArray_ISBOOL((PyArrayObject *)obj_ptr)) {
            obj_ptr = PyArray_Cast((PyArrayObject *)obj_ptr, NPY_UINT8);
            target = bp::object(bp::handle<>(obj_ptr));
        }

    BufferWrapper<T> buf("src", target, false);

    if (buf->ndim != 1)
        throw shape_exception("src", "must be 1-d");

    int p_count = buf->shape[0];
    void *p = buf->buf;

    int dtype = format_to_dtype(buf);
    switch(dtype) {
    case NPY_BOOL:
    case NPY_UINT8:
    case NPY_INT8:
        return from_bitmask_<T,uint8_t>(p, p_count, n_bits);
    case NPY_UINT16:
    case NPY_INT16:
        return from_bitmask_<T,uint16_t>(p, p_count, n_bits);
    case NPY_UINT32:
    case NPY_INT32:
        return from_bitmask_<T,uint32_t>(p, p_count, n_bits);
    case NPY_UINT64:
    case NPY_INT64:
        return from_bitmask_<T,uint64_t>(p, p_count, n_bits);
    }

    throw dtype_exception("src", "integer type");
    return bp::object();
}

template <typename T>
bp::object Ranges<T>::from_mask(const bp::object &src)
{
    return from_bitmask(src, 0);
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
static inline bp::object mask_(vector<Ranges<intType>> ivals, int n_bits)
{
    intType x;
    throw dtype_exception("ivlist", "Interval<> over integral type.");
    return bp::object();
}

template <typename intType, typename std::enable_if<std::is_integral<intType>::value,
                                                    int>::type* = nullptr>
static inline bp::object mask_(vector<Ranges<intType>> ivals, int n_bits)
{
    vector<int> indexes;
    int count = 0;
    
    for (long i=0; i<ivals.size(); i++) {
        indexes.push_back(0);
        if (i==0) {
            count = ivals[i].count;
        } else if (count != ivals[i].count) {
            throw agreement_exception("ranges[0]", "all other ranges[i]", "count");
        }
    }

    // Determine the output mask size based on n_bits... which may be unspecified.
    int npy_type = NPY_NOTYPE;

    if (n_bits < 0)
        n_bits = ivals.size();
    else if (n_bits == 0 && ivals.size() == 1)
        npy_type = NPY_BOOL;  // This will return a boolean array.
    else if (n_bits < ivals.size())
        throw general_agreement_exception("Input list has more items than the "
                                          "output mask size (n_bits).");

    if (npy_type == NPY_BOOL) { }
    else if (n_bits <= 8)
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

    npy_intp dims[1] = {count};
    PyObject *v = PyArray_SimpleNew(1, dims, npy_type);

    // Assumes little-endian.
    int n_byte = PyArray_ITEMSIZE((PyArrayObject*)v);
    uint8_t *ptr = reinterpret_cast<uint8_t*>((PyArray_DATA((PyArrayObject*)v)));
    memset(ptr, 0, count*n_byte);
    for (long bit=0; bit<ivals.size(); ++bit) {
        for (auto p: ivals[bit].segments) {
            for (int i=p.first; i<p.second; i++)
                ptr[i*n_byte + bit/8] |= (1<<(bit%8));
        }
    }

    return bp::object(bp::handle<>(v));
}

template <typename intType>
static inline bp::object mask_(const bp::list &ivlist, int n_bits)
{
    vector<Ranges<intType>> ivals;
    for (long i=0; i<bp::len(ivlist); i++)
        ivals.push_back(bp::extract<Ranges<intType>>(ivlist[i]));
    return mask_(ivals, n_bits);
}


template <typename T>
bp::object Ranges<T>::bitmask(const bp::list &ivlist, int n_bits)
{
    return mask_<T>(ivlist, n_bits);
}

template <typename T>
bp::object Ranges<T>::mask()
{
    vector<Ranges<T>> temp;
    temp.push_back(*this);
    return mask_<T>(temp, 0);
}



//
// Implementation of the algebra
//
 
template <typename T>
Ranges<T>& Ranges<T>::intersect(const Ranges<T> &src)
{
    *this = (this->complement() + src.complement()).complement();
    return *this;
}
    
template <typename T>
Ranges<T> Ranges<T>::complement() const
{
    Ranges<T> output(count, reference);

    T next_start = 0;
    for (auto p: segments) {
        output.segments.push_back(make_pair(next_start, p.first));
        next_start = p.second;
    }
    output.segments.push_back(make_pair(next_start, count));
    output.cleanup();
    return output;
}

// Make empty range to match this range
template <typename T>
Ranges<T> Ranges<T>::zeros_like() const
{
    Ranges<T> output(count, reference);
    return output;
    
}

//make "full" range to match this range
template <typename T>
Ranges<T> Ranges<T>::ones_like() const
{
    Ranges<T> output(count, reference);
    output.add_interval(0, count);
    return output;
    
}


// Slicing!
template <typename T,
          typename std::enable_if<!std::is_integral<T>::value,
                                  int>::type* = nullptr>
static inline Ranges<T> _getitem_(Ranges<T> &src, bp::object indices)
{
    throw dtype_exception("target", "Interval<> over integral type.");
    return Ranges<T>();
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
static inline Ranges<T> _getitem_(Ranges<T> &src, bp::object indices)
{
    bp::object target = indices;

    // Allow user to pass in a length-one tuple of slices.
    bp::extract<bp::tuple> tu_ex(target);
    if (tu_ex.check()) {
        if (bp::len(tu_ex()) == 0)
            target = bp::slice();
        else {
            assert(bp::len(tu_ex()) == 1);
            target = tu_ex()[0];
        }
    }

    bp::extract<bp::slice> ex(target);
    if (ex.check()) {
        T n = src.count;

        auto sl = ex();
        T start = extract_or_default(sl.start(), 0);
        T stop = extract_or_default(sl.stop(), n);
        T step = extract_or_default(sl.step(), 1);

        assert(step == 1);
        if (start < 0)
            start = n + start;
        if (stop < 0)
            stop = n + stop;
        if (stop < start)
            stop = start;
        if (start > n)
            start = n;
        if (stop > n)
            stop = n;

        auto output = Ranges<T>(stop - start, src.reference - start);
        for (auto p: src.segments)
            if (p.second > start && p.first < stop)
                output.segments.push_back(make_pair(p.first - start, p.second - start));
        output.cleanup();

        return output;
    }
    return Ranges<T>();
}

template <typename T>
Ranges<T> Ranges<T>::getitem(bp::object indices)
{
    return _getitem_(*this, indices);
}

template <typename T>
bp::object Ranges<T>::shape()
{
    vector<T> temp = {count};
    return bp::tuple(temp);
}

template <typename T>
void Ranges<T>::safe_set_count(T count_)
{
    count = count_;
    cleanup();
}


//
// Operators
//

template <typename T>
Ranges<T> Ranges<T>::operator~() const
{
    return complement();
}

template <typename T>
void Ranges<T>::operator+=(const Ranges<T> &src)
{
    merge(src);
}

template <typename T>
void Ranges<T>::operator*=(const Ranges<T> &src)
{
    intersect(src);
}

template <typename T>
Ranges<T> Ranges<T>::operator+(const Ranges<T> &src) const
{
    Ranges<T> output = *this;
    output += src;
    return output;
}

template <typename T>
Ranges<T> Ranges<T>::operator*(const Ranges<T> &src) const
{
    auto output = *this;
    output *= src;
    return output;
}


//
// boost-python registration.
//

using namespace boost::python;


#define EXPORT_RANGES(DOMAIN_TYPE, CLASSNAME)                           \
    bp::class_<CLASSNAME>(#CLASSNAME,                                   \
    "A finite series of non-overlapping semi-open intervals on a domain\n" \
    "of type: " #DOMAIN_TYPE ".\n\n"                                    \
    "To create an empty object, instantiate with just a sample count:\n" \
    "``" #CLASSNAME "(count)``.\n"                                      \
    "\n"                                                                \
    "Alternately, consider convenience methods such as ``from_mask``,\n" \
    "``from_array``, and ``from_bitmask``; see below.\n"                \
    "\n"                                                                \
    "In addition to the methods explained below, note the that following\n" \
    "operators have been defined and perform as follows (where ``r1`` and\n" \
    "``r2`` are objects of this class:\n"                               \
    "\n"                                                                \
    "- ``~r1`` is equivalent to ``r1.complement()``\n"                  \
    "- ``r1 *= r2`` is equivalent to ``r1.intersect(r2)``\n"            \
    "- ``r1 += r2`` is equivalent to ``r1.merge(r2)``\n"                \
    "- ``r1 * r2`` and ``r1 + r2`` behave as you might expect, returning a\n" \
    "  new object and leaving ``r1`` and ``r2`` unmodified.\n"          \
    "\n"                                                                \
    "The object also supports slicing.  For example, if ``r1`` has\n"   \
    "count = 100 then r1[10:-5] returns a new object (not a view)\n"    \
    "that has count = 85.  A data member ``reference`` keeps track\n"   \
    "of the history of shifts in the first sample; for example if\n"    \
    "r1.reference = 0 then r1[10:-5].reference will be -10.  This\n"    \
    "variable can be interpreted as giving the logical index, in\n"     \
    "the new index system, of where index=0 of the original object\n"   \
    "would be found.  This is useful for bookkeeping in some cases.\n") \
    .def(init<const DOMAIN_TYPE&>("Initialize with count."))            \
    .def(init<const DOMAIN_TYPE&, const DOMAIN_TYPE&>("Initialize with count and reference.")) \
    .def("__str__", &CLASSNAME::Description) \
    .add_property("count", &CLASSNAME::count, &CLASSNAME::safe_set_count) \
    .add_property("reference", &CLASSNAME::reference)                   \
    .def("add_interval", &CLASSNAME::_add_interval_numpysafe,           \
         return_internal_reference<>(),                                 \
         args("self", "start", "end"),                                  \
         "Merge an interval into the set.")                             \
    .def("append_interval_no_check", &CLASSNAME::append_interval_no_check, \
         return_internal_reference<>(),                                 \
         args("self", "start", "end"),                                  \
         "Append an interval to the set without checking for overlap or sequence.") \
    .def("merge", &CLASSNAME::merge,                                    \
         return_internal_reference<>(),                                 \
         args("self", "src"),                                           \
         "Merge ranges from another " #CLASSNAME " into this one.")     \
    .def("buffer", &CLASSNAME::buffer,                                  \
        return_internal_reference<>(),                                  \
        args("self", "buff"),                                           \
        "Buffer each interval by an amount specified by buff")          \
    .def("buffered", &CLASSNAME::buffered,                              \
        args("self", "buff"),                                           \
        "Return an interval buffered by buff")                          \
    .def("close_gaps", &CLASSNAME::close_gaps,                          \
        return_internal_reference<>(),                                  \
        args("self", "gap"),                                            \
        "Remove gaps between ranges less than gap")                     \
    .def("intersect", &CLASSNAME::intersect,                            \
         return_internal_reference<>(),                                 \
         args("self", "src"),                                           \
         "Intersect another " #CLASSNAME " with this one.")             \
    .def("complement", &CLASSNAME::complement,                          \
         "Return the complement (over domain).")                        \
    .def("zeros_like", &CLASSNAME::zeros_like,                          \
         "Return range of same length but no intervals")                \
    .def("ones_like", &CLASSNAME::ones_like,                            \
         "Return range of same length and interval spanning count")     \
    .def("ranges", &CLASSNAME::ranges,                                  \
         "Return the intervals as a 2-d numpy array of ranges.")        \
    .def("from_array", &CLASSNAME::from_array,                          \
         args("data", "count"),                                         \
         "The input data must be an (n,2) shape ndarray of int32. "     \
         "The integer count sets the domain of the object.")            \
    .staticmethod("from_array")                                         \
    .def("from_bitmask", &CLASSNAME::from_mask,                         \
         args("bitmask_array"),                                         \
         "Return a list of " #CLASSNAME " extracted from an ndarray encoding a bitmask.") \
    .staticmethod("from_bitmask")                                       \
    .def("from_mask", &CLASSNAME::from_mask,                            \
         args("bool_array"),                                            \
         "Return a list of " #CLASSNAME " extracted from an ndarray of bool.") \
    .staticmethod("from_mask")                                          \
    .def("bitmask", &CLASSNAME::bitmask,                                \
         args("ranges_list", "n_bits"),                                 \
         "Return an ndarray bitmask from a list of" #CLASSNAME ".\n"    \
         "n_bits determines the output integer type.  Bits are assigned from \n" \
         "LSB onwards; use None in the list to skip a bit.")            \
    .staticmethod("bitmask")                                            \
    .def("mask", &CLASSNAME::mask,                                      \
         "Return a boolean mask from this Ranges object.")              \
    .def("copy",                                                        \
         +[](CLASSNAME& A) {                                            \
             return CLASSNAME(A);                                       \
         },                                                             \
         "Get a new object with a copy of the data.")                   \
    .def("__getitem__", &CLASSNAME::getitem)                            \
    .add_property("shape", &CLASSNAME::shape)                           \
    .def(~self)                                                         \
    .def(self += self)                                                  \
    .def(self *= self)                                                  \
    .def(self + self)                                                   \
    .def(self * self);


PYBINDINGS("so3g")
{
    docstring_options local_docstring_options(true, true, false);
    EXPORT_RANGES(int32_t, RangesInt32);
}
