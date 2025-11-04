#define NO_IMPORT_ARRAY

#include <iostream>
#include <limits>
#include <type_traits>

#include <nanobind/nanobind.h>
#include <nanobind/operators.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/string.h>

#include "so3g_numpy.h"

#include "Ranges.h"
#include "exceptions.h"
#include <exception>

namespace nb = nanobind;


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
    const nb::object start_obj, const nb::object end_obj)
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
Ranges<T> Ranges<T>::from_array(const nb::object &src, const nb::object &count)
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
nb::object Ranges<T>::ranges() const
{
    npy_intp dims[2] = {0, 2};
    dims[0] = (npy_intp)segments.size();
    int dtype = get_dtype<T>();
    if (dtype == NPY_NOTYPE)
        throw general_agreement_exception("ranges() not implemented for this domain dtype.");

    PyObject *v = PyArray_SimpleNew(2, dims, dtype);
    if (v == NULL) {
        ostringstream dstr;
        dstr << "Failed to allocate Ranges numpy array of size (";
        dstr << dims[0] << ", " << dims[1] << ")";
        throw RuntimeError_exception(dstr.str().c_str());
    }
    char *ptr = reinterpret_cast<char*>((PyArray_DATA((PyArrayObject*)v)));
    for (auto p = segments.begin(); p != segments.end(); ++p) {
        ptr += interval_extract((&*p), ptr);
    }
    return nb::steal<nb::object>(v);
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
static inline nb::object from_bitmask_(void *buf, intType count, int n_bits)
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
        return nb::cast(output[0]);

    // Once added to the list, we can't modify further.
    nb::list bits;
    for (auto i: output)
        bits.append(i);
    return bits;
}

template <typename T>
nb::object Ranges<T>::from_bitmask(const nb::object &src, int n_bits)
{
    // Buffer protocol doesn't work directly on bool arrays, so if
    // what we've been passed is definitely a bool array, get a view
    // of it as a uint8 array and work with that.  (Wrap it with
    // nb::object so references are counted properly.)
    nb::object target(src);
    PyObject* obj_ptr = target.ptr();
    if (PyArray_Check(obj_ptr))
        if (PyArray_ISBOOL((PyArrayObject *)obj_ptr)) {
            obj_ptr = PyArray_Cast((PyArrayObject *)obj_ptr, NPY_UINT8);
            target = nb::steal<nb::object>(obj_ptr);
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
    return nb::object();
}

template <typename T>
nb::object Ranges<T>::from_mask(const nb::object &src)
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
static inline nb::object mask_(vector<Ranges<intType>> ivals, int n_bits)
{
    intType x;
    throw dtype_exception("ivlist", "Interval<> over integral type.");
    return nb::object();
}

template <typename intType, typename std::enable_if<std::is_integral<intType>::value,
                                                    int>::type* = nullptr>
static inline nb::object mask_(vector<Ranges<intType>> ivals, int n_bits)
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
    if (v == NULL) {
        ostringstream dstr;
        dstr << "Failed to allocate Ranges mask array of size (";
        dstr << dims[0] << ",)";
        throw RuntimeError_exception(dstr.str().c_str());
    }

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

    return nb::steal<nb::object>(v);
}

template <typename intType>
static inline nb::object mask_(const nb::list &ivlist, int n_bits)
{
    vector<Ranges<intType>> ivals;
    for (long i=0; i<nb::len(ivlist); i++)
        ivals.push_back(nb::cast<Ranges<intType>>(ivlist[i]));
    return mask_(ivals, n_bits);
}


template <typename T>
nb::object Ranges<T>::bitmask(const nb::list &ivlist, int n_bits)
{
    return mask_<T>(ivlist, n_bits);
}

template <typename T>
nb::object Ranges<T>::mask()
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
static inline Ranges<T> _getitem_(Ranges<T> &src, nb::object indices)
{
    throw dtype_exception("target", "Interval<> over integral type.");
    return Ranges<T>();
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
static inline Ranges<T> _getitem_(Ranges<T> &src, nb::object indices)
{
    nb::object target = indices;

    // Allow user to pass in a length-one tuple of slices.
    nb::slice slc;
    if (nb::isinstance<nb::tuple>(indices)) {
        nb::tuple temp = nb::cast<nb::tuple>(indices);
        if (nb::len(temp) != 1) {
            throw general_agreement_exception("Ranges __getitem__ got tuple of len >1");
        }
        slc = nb::cast<nb::slice>(temp[0]);
    } else if (nb::isinstance<nb::slice>(indices)) {
        slc = nb::cast<nb::slice>(indices);
    } else {
        return Ranges<T>();
    }

    T n = src.count;
    auto slc_par = slc.compute(n);
    T start = slc_par.template get<0>();
    T stop = slc_par.template get<1>();
    T step = slc_par.template get<2>();

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

template <typename T>
Ranges<T> Ranges<T>::getitem(nb::object indices)
{
    return _getitem_(*this, indices);
}

template <typename T>
nb::object Ranges<T>::shape()
{
    return nb::make_tuple(count);
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


// Helper function to register an Ranges class for a concrete type.

template <typename C>
void ranges_bindings(nb::module_ & m, char const * name) {

    nb::class_<Ranges<C>>(m, name,
            R"(
            A finite series of non-overlapping semi-open intervals on a domain.

            To create an empty object, instantiate with just a sample count.
            Alternately, consider convenience methods such as ``from_mask``,
            ``from_array``, and ``from_bitmask``; see below.

            In addition to the methods explained below, note the that following
            operators have been defined and perform as follows (where ``r1`` and
            ``r2`` are objects of this class:

            - ``~r1`` is equivalent to ``r1.complement()``
            - ``r1 *= r2`` is equivalent to ``r1.intersect(r2)``
            - ``r1 += r2`` is equivalent to ``r1.merge(r2)``
            - ``r1 * r2`` and ``r1 + r2`` behave as you might expect, returning a
                new object and leaving ``r1`` and ``r2`` unmodified.

            The object also supports slicing.  For example, if ``r1`` has
            count = 100 then r1[10:-5] returns a new object (not a view)
            that has count = 85.  A data member ``reference`` keeps track
            of the history of shifts in the first sample; for example if
            r1.reference = 0 then r1[10:-5].reference will be -10.  This
            variable can be interpreted as giving the logical index, in
            the new index system, of where index=0 of the original object
            would be found.  This is useful for bookkeeping in some cases.
            )"
        )
        .def(nb::init<>())
        .def(nb::init<C>(),
            R"(
            Initialize with count.
            )"
        )
        .def(nb::init<C, C>(),
            R"(
            Initialize with count and reference.
            )"
        )
        .def("__str__", &Ranges<C>::Description)
        .def_prop_rw("count",
            [](Ranges<C> & slf){
                return slf.count;
            },
            &Ranges<C>::safe_set_count
        )
        .def_prop_ro("reference",
            [](Ranges<C> & slf){
                return slf.reference;
            }
        )
        .def_prop_ro("shape", &Ranges<C>::shape)
        .def("add_interval", &Ranges<C>::add_interval, nb::rv_policy::none,
            nb::arg("start"),
            nb::arg("end"),
            R"(
            Merge an interval into the set.
            )"
        )
        .def("append_interval_no_check", &Ranges<C>::append_interval_no_check,
            nb::rv_policy::none,
            nb::arg("start"),
            nb::arg("end"),
            R"(
            Append an interval to the set without checking for overlap or sequence.
            )"
        )
        .def("merge", &Ranges<C>::merge, nb::rv_policy::none,
            R"(
            Merge ranges from another object into this one.
            )"
        )
        .def("intersect", &Ranges<C>::intersect, nb::rv_policy::none,
            nb::arg("source"),
            R"(
            Intersect another Ranges object with this one.
            )"
        )
        .def("complement", &Ranges<C>::complement, nb::rv_policy::take_ownership,
            R"(
            Return the complement (over domain).
            )"
        )
        .def(
            "copy",
            [](Ranges<C> & slf) {return Ranges<C>(slf);},
            nb::rv_policy::take_ownership,
            R"(
            Get a new object with a copy of the data.
            )"
        )
        .def("buffer", &Ranges<C>::buffer, nb::rv_policy::none,
            nb::arg("buff"),
            R"(
            Buffer each interval by an amount specified by buff
            )"
        )
        .def("buffered", &Ranges<C>::buffered, nb::rv_policy::take_ownership,
            nb::arg("buff"),
            R"(
            Return an interval buffered by buff
            )"
        )
        .def("close_gaps", &Ranges<C>::close_gaps, nb::rv_policy::none,
            nb::arg("gap"),
            R"(
            Remove gaps between ranges less than gap
            )"
        )
        .def("zeros_like", &Ranges<C>::zeros_like, nb::rv_policy::take_ownership,
            R"(
            Return range of same length but no intervals
            )"
        )
        .def("ones_like", &Ranges<C>::ones_like, nb::rv_policy::take_ownership,
            R"(
            Return range of same length and interval spanning count
            )"
        )
        .def("ranges", &Ranges<C>::ranges, nb::rv_policy::take_ownership,
            R"(
            Return the intervals as a 2-d numpy array of ranges.
            )"
        )
        .def_static("from_array", &Ranges<C>::from_array,
            nb::rv_policy::take_ownership,
            nb::arg("src"),
            nb::arg("count"),
            R"(
            The input data must be an (n,2) shape ndarray of int32.
            The integer count sets the domain of the object.
            )"
        )
        .def("__getitem__", &Ranges<C>::getitem, nb::rv_policy::take_ownership)
        .def_static("from_bitmask", &Ranges<C>::from_mask,
            nb::rv_policy::take_ownership,
            nb::arg("bitmask_array"),
            R"(
            Return a list of Ranges extracted from an ndarray encoding a bitmask.
            )"
        )
        .def_static("from_mask", &Ranges<C>::from_mask,
            nb::rv_policy::take_ownership,
            nb::arg("bool_array"),
            R"(
            Return a list of Ranges extracted from an ndarray of bool.
            )"
        )
        .def_static("bitmask", &Ranges<C>::bitmask, nb::rv_policy::take_ownership,
            nb::arg("ranges_list"),
            nb::arg("n_bits"),
            R"(
            Return an ndarray bitmask from a list of Ranges.  n_bits determines
            the output integer type.  Bits are assigned from LSB onwards; use None
            in the list to skip a bit.
            )"
        )
        .def("mask", &Ranges<C>::mask, nb::rv_policy::take_ownership,
            R"(
            Return a boolean mask from this Ranges object.
            )"
        )
        .def(~nb::self)
        .def(nb::self += nb::self)
        .def(nb::self *= nb::self)
        .def(nb::self + nb::self)
        .def(nb::self * nb::self);

    return;
}


void register_ranges(nb::module_ & m) {
    // Concrete Ranges types
    ranges_bindings<int32_t>(m, "RangesInt32");
    return;
}
