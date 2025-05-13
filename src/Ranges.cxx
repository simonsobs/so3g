
#include <iostream>
#include <limits>
#include <type_traits>

// Needed to access self in bindings
#include <nanobind/operators.h>

// Needed to work with numpy arrays in bindings
#include <nanobind/ndarray.h>

// Needed to access std::tuple return from slice.compute()
#include <nanobind/stl/tuple.h>

#include <bindings.h>
#include <exceptions.h>
#include <Ranges.h>


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


// Conversions between Ranges and bit mask arrays.  These are templated on both the
// Ranges type and the mask array type.

template <
    typename C,
    typename M,
    typename std::enable_if_t<std::is_integral<C>::value, bool>* = nullptr
>
nb::object ranges_from_mask(
    nb::ndarray_view<M, 1, 'C'> & view, int n_bits
) {
    bool return_singleton = (n_bits == 0);

    auto p = view.data();
    if (n_bits < 0) {
        // Use all bits
        n_bits = 8 * sizeof(M);
    } else if (n_bits == 0) {
        // Single mask
        n_bits = 1;
    }
    C count = view.shape(0);

    vector<Ranges<C>> output;
    vector<C> start;
    for (int bit = 0; bit < n_bits; ++bit) {
        output.push_back(Ranges<C>(count));
        start.push_back(-1);
    }

    M last = 0;
    for (C i = 0; i < count; i++) {
        M d = p[i] ^ last;
        for (int bit = 0; bit < n_bits; ++bit) {
            if (d & (1 << bit)) {
                if (start[bit] >= 0) {
                    output[bit].segments.push_back(make_pair(start[bit], i));
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
            output[bit].segments.push_back(make_pair(start[bit], count));
        }
    }

    if (return_singleton) {
        return nb::cast(output[0]);
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
nb::object ranges_from_mask(
    nb::ndarray_view<M, 1, 'C'> & view, int n_bits
) {
    throw general_agreement_exception("from_mask() is only valid for integral types");
    return nb::object();
}


template <
    typename C,
    typename M,
    typename std::enable_if_t<std::is_integral<C>::value, bool>* = nullptr
>
nb::object ranges_to_mask(nb::list rlist, int n_bits) {
    if (n_bits < 0) {
        n_bits = rlist.size();
    } else if (n_bits == 0) {
        n_bits = 1;
    }
    // Ensure that the mask array type has sufficient number of bits.
    if (n_bits > nb::dtype<M>().bits) {
        std::ostringstream err;
        err << "Mask array type has " << nb::dtype<M>().bits <<
            " bits, which cannot hold the requested " << n_bits << " bits.";
        throw general_agreement_exception(err.str().c_str());
    }

    vector<Ranges<C>> rvals;
    vector<int> indexes;
    C count;

    for (size_t i = 0; i < rlist.size(); i++) {
        // Convert list element to a C++ object.  Note that this does make a copy
        // when pushed onto the std::vector.
        rvals.push_back(nb::cast<Ranges<C>>(rlist[i]));
        indexes.push_back(0);
        if (i == 0) {
            count = rvals[i].count;
        } else if (count != rvals[i].count) {
            throw agreement_exception("ranges[0]", "all other ranges[i]", "count");
        }
    }

    // Create the returned array
    size_t shape[1] = {(size_t)count};
    auto array = create_ndarray<M, 1>(shape);
    auto raw = array.data();
    auto n_byte = array.itemsize();
    uint8_t * ptr = reinterpret_cast<uint8_t*>(raw);

    // Fill the array
    for (size_t bit = 0; bit < rvals.size(); ++bit) {
        for (auto & p : rvals[bit].segments) {
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
nb::object ranges_to_mask(nb::list rlist, int n_bits) {
    throw general_agreement_exception("mask() is only valid for integral types");
    return nb::object();
}


// Helper function to register a Ranges class for a concrete type.

template <typename C>
void ranges_bindings(nb::module_ & m, char const * name) {

    nb::class_<Ranges<C>>(m, name)
        .def(nb::init<>())
        .def(nb::init<C>())
        .def(nb::init<C, C>())
        .def("__str__", &Ranges<C>::Description)
        .def("add_interval", [](Intervals<C> & slf, C start, C end){
            slf.add_interval(start, end);
        }, nb::rv_policy::none)
        .def("add_interval", [](Intervals<C> & slf, nb::object start, nb::object end){
            slf.add_interval(nb::cast<C>(start), nb::cast<C>(end));
        }, nb::rv_policy::none)
        .def(
            "append_interval_no_check",
            &Ranges<C>::append_interval_no_check,
            nb::rv_policy::none
        )
        .def("merge", &Ranges<C>::merge, nb::rv_policy::none)
        .def("intersect", &Ranges<C>::intersect, nb::rv_policy::none)
        .def("buffer", &Ranges<C>::buffer, nb::rv_policy::none)
        .def("close_gaps", &Ranges<C>::close_gaps, nb::rv_policy::none)
        .def("buffered", &Ranges<C>::buffered, nb::rv_policy::take_ownership)
        .def("complement", &Ranges<C>::complement, nb::rv_policy::take_ownership)
        .def("zeros_like", &Ranges<C>::zeros_like, nb::rv_policy::take_ownership)
        .def("ones_like", &Ranges<C>::ones_like, nb::rv_policy::take_ownership)
        .def("shape", [](Ranges<C> & slf) {
            return nb::make_tuple({slf.count});
        }, nb::rv_policy::take_ownership)
        .def_static("from_array", [](
            nb::ndarray<C, nb::ndim<2>, nb::c_contig> input, nb::object count
        ) {
            Ranges<C> output;
            // Fast-access view of the array
            auto v = input.view();
            size_t n_seg = v.shape(0);
            if (v.shape(1) != 2) {
                throw shape_exception("input array shape[1]", "!= 2");
            }
            auto data = v.data();
            for (size_t i = 0; i < n_seg; ++i) {
                output.append_interval_nocheck(data[2*i], data[2*i + 1]);
            }
            output.count = nb::cast<C>(count);
            output.cleanup();
            return output;
        }, nb::rv_policy::take_ownership)
        .def("ranges", [](Ranges<C> & slf) {
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
        }, nb::rv_policy::take_ownership)
        .def("__getitem__", [](Ranges<C> & slf, nb::object indices) {
            if (! std::is_integral<C>::value) {
                throw ValueError_exception("Ranges __getitem__: type not integral");
            }
            // From original code- Accept either a slice or a single-element tuple
            // with a slice.  Not sure why...
            nb::slice slc;
            if (nb::isinstance<nb::tuple>(indices)) {
                slc = nb::cast<nb::slice>(nb::cast<nb::tuple>(indices)[0]);
            } else {
                slc = nb::cast<nb::slice>(indices);
            }

            C n = slf.count;
            auto slc_par = slc.compute(n);
            auto start = slc_par.template get<0>();
            auto stop = slc_par.template get<1>();
            auto step = slc_par.template get<2>();
            auto slen = slc_par.template get<3>();

            if (step != 1) {
                throw general_agreement_exception("Ranges __getitem__: step size != 1");
            }
            if (start < 0) {
                start = n + start;
            }
            if (stop < 0) {
                stop = n + stop;
            }
            if (stop < start) {
                stop = start;
            }
            if (start > n) {
                start = n;
            }
            if (stop > n) {
                stop = n;
            }

            auto output = Ranges<T>(stop - start, slf.reference - start);
            for (auto const & p : slf.segments) {
                if (p.second > start && p.first < stop) {
                    output.segments.push_back(
                        make_pair(p.first - start, p.second - start)
                    );
                }
            }
            output.cleanup();
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
                return ranges_from_mask<C, uint8_t>(view, n_bits);
            } else if (input.dtype() == nb::dtype<uint16_t>()) {
                auto view = input.view<uint16_t>();
                return ranges_from_mask<C, uint16_t>(view, n_bits);
            } else if (input.dtype() == nb::dtype<uint32_t>()) {
                auto view = input.view<uint32_t>();
                return ranges_from_mask<C, uint32_t>(view, n_bits);
            } else if (input.dtype() == nb::dtype<uint64_t>()) {
                auto view = input.view<uint64_t>();
                return ranges_from_mask<C, uint64_t>(view, n_bits);
            } else {
                std::ostringstream err;
                err << "Mask array type is not one of the supported types"
                    << " (uint8, uint16, uint32, uint64)";
                throw general_agreement_exception(err.str().c_str());
            }
        }, nb::rv_policy::take_ownership)
        .def_static("mask", [](nb::list rlist, int n_bits){
            if (n_bits < 0) {
                // Use the length of the interval list
                n_bits = rlist.size();
            }
            if (n_bits <= 8) {
                return ranges_to_mask<C, uint8_t>(rlist, n_bits);
            } else if (n_bits <= 16) {
                return ranges_to_mask<C, uint16_t>(rlist, n_bits);
            } else if (n_bits <= 32) {
                return ranges_to_mask<C, uint32_t>(rlist, n_bits);
            } else if (n_bits <= 64) {
                return ranges_to_mask<C, uint64_t>(rlist, n_bits);
            } else {
                throw general_agreement_exception(
                    "Number of requested bits must be <= 64"
                );
            }
        }, nb::rv_policy::take_ownership)
        .def(-nb::self)
        .def(~nb::self)
        .def(nb::self += nb::self)
        .def(nb::self *= nb::self)
        .def(nb::self + nb::self)
        .def(nb::self * nb::self);




// Support for working with RangesMatrix, which is basically just a list of Ranges
template <typename T>
vector<Ranges<T>> extract_ranges(const bp::object & ival_list) {
    vector<Ranges<T>> v(bp::len(ival_list));
    for (int i=0; i<bp::len(ival_list); i++)
        v[i] = bp::extract<Ranges<T>>(ival_list[i])();
    return v;
}

    return;
}



void register_ranges(nb::module_ & m) {


    return;
}


// //
// // boost-python registration.
// //

// using namespace boost::python;


// #define EXPORT_RANGES(DOMAIN_TYPE, CLASSNAME)                           \
//     bp::class_<CLASSNAME>(#CLASSNAME,                                   \
//     "A finite series of non-overlapping semi-open intervals on a domain\n" \
//     "of type: " #DOMAIN_TYPE ".\n\n"                                    \
//     "To create an empty object, instantiate with just a sample count:\n" \
//     "``" #CLASSNAME "(count)``.\n"                                      \
//     "\n"                                                                \
//     "Alternately, consider convenience methods such as ``from_mask``,\n" \
//     "``from_array``, and ``from_bitmask``; see below.\n"                \
//     "\n"                                                                \
//     "In addition to the methods explained below, note the that following\n" \
//     "operators have been defined and perform as follows (where ``r1`` and\n" \
//     "``r2`` are objects of this class:\n"                               \
//     "\n"                                                                \
//     "- ``~r1`` is equivalent to ``r1.complement()``\n"                  \
//     "- ``r1 *= r2`` is equivalent to ``r1.intersect(r2)``\n"            \
//     "- ``r1 += r2`` is equivalent to ``r1.merge(r2)``\n"                \
//     "- ``r1 * r2`` and ``r1 + r2`` behave as you might expect, returning a\n" \
//     "  new object and leaving ``r1`` and ``r2`` unmodified.\n"          \
//     "\n"                                                                \
//     "The object also supports slicing.  For example, if ``r1`` has\n"   \
//     "count = 100 then r1[10:-5] returns a new object (not a view)\n"    \
//     "that has count = 85.  A data member ``reference`` keeps track\n"   \
//     "of the history of shifts in the first sample; for example if\n"    \
//     "r1.reference = 0 then r1[10:-5].reference will be -10.  This\n"    \
//     "variable can be interpreted as giving the logical index, in\n"     \
//     "the new index system, of where index=0 of the original object\n"   \
//     "would be found.  This is useful for bookkeeping in some cases.\n") \
//     .def(init<const DOMAIN_TYPE&>("Initialize with count."))            \
//     .def(init<const DOMAIN_TYPE&, const DOMAIN_TYPE&>("Initialize with count and reference.")) \
//     .def("__str__", &CLASSNAME::Description) \
//     .add_property("count", &CLASSNAME::count, &CLASSNAME::safe_set_count) \
//     .add_property("reference", &CLASSNAME::reference)                   \
//     .def("add_interval", &CLASSNAME::_add_interval_numpysafe,           \
//          return_internal_reference<>(),                                 \
//          args("self", "start", "end"),                                  \
//          "Merge an interval into the set.")                             \
//     .def("append_interval_no_check", &CLASSNAME::append_interval_no_check, \
//          return_internal_reference<>(),                                 \
//          args("self", "start", "end"),                                  \
//          "Append an interval to the set without checking for overlap or sequence.") \
//     .def("merge", &CLASSNAME::merge,                                    \
//          return_internal_reference<>(),                                 \
//          args("self", "src"),                                           \
//          "Merge ranges from another " #CLASSNAME " into this one.")     \
//     .def("buffer", &CLASSNAME::buffer,                                  \
//         return_internal_reference<>(),                                  \
//         args("self", "buff"),                                           \
//         "Buffer each interval by an amount specified by buff")          \
//     .def("buffered", &CLASSNAME::buffered,                              \
//         args("self", "buff"),                                           \
//         "Return an interval buffered by buff")                          \
//     .def("close_gaps", &CLASSNAME::close_gaps,                          \
//         return_internal_reference<>(),                                  \
//         args("self", "gap"),                                            \
//         "Remove gaps between ranges less than gap")                     \
//     .def("intersect", &CLASSNAME::intersect,                            \
//          return_internal_reference<>(),                                 \
//          args("self", "src"),                                           \
//          "Intersect another " #CLASSNAME " with this one.")             \
//     .def("complement", &CLASSNAME::complement,                          \
//          "Return the complement (over domain).")                        \
//     .def("zeros_like", &CLASSNAME::zeros_like,                          \
//          "Return range of same length but no intervals")                \
//     .def("ones_like", &CLASSNAME::ones_like,                            \
//          "Return range of same length and interval spanning count")     \
//     .def("ranges", &CLASSNAME::ranges,                                  \
//          "Return the intervals as a 2-d numpy array of ranges.")        \
//     .def("from_array", &CLASSNAME::from_array,                          \
//          args("data", "count"),                                         \
//          "The input data must be an (n,2) shape ndarray of int32. "     \
//          "The integer count sets the domain of the object.")            \
//     .staticmethod("from_array")                                         \
//     .def("from_bitmask", &CLASSNAME::from_mask,                         \
//          args("bitmask_array"),                                         \
//          "Return a list of " #CLASSNAME " extracted from an ndarray encoding a bitmask.") \
//     .staticmethod("from_bitmask")                                       \
//     .def("from_mask", &CLASSNAME::from_mask,                            \
//          args("bool_array"),                                            \
//          "Return a list of " #CLASSNAME " extracted from an ndarray of bool.") \
//     .staticmethod("from_mask")                                          \
//     .def("bitmask", &CLASSNAME::bitmask,                                \
//          args("ranges_list", "n_bits"),                                 \
//          "Return an ndarray bitmask from a list of" #CLASSNAME ".\n"    \
//          "n_bits determines the output integer type.  Bits are assigned from \n" \
//          "LSB onwards; use None in the list to skip a bit.")            \
//     .staticmethod("bitmask")                                            \
//     .def("mask", &CLASSNAME::mask,                                      \
//          "Return a boolean mask from this Ranges object.")              \
//     .def("copy",                                                        \
//          +[](CLASSNAME& A) {                                            \
//              return CLASSNAME(A);                                       \
//          },                                                             \
//          "Get a new object with a copy of the data.")                   \
//     .def("__getitem__", &CLASSNAME::getitem)                            \
//     .add_property("shape", &CLASSNAME::shape)                           \
//     .def(~self)                                                         \
//     .def(self += self)                                                  \
//     .def(self *= self)                                                  \
//     .def(self + self)                                                   \
//     .def(self * self);


// PYBINDINGS("so3g")
// {
//     docstring_options local_docstring_options(true, true, false);
//     EXPORT_RANGES(int32_t, RangesInt32);
// }
