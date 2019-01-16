#include <pybindings.h>

#include <iostream>
#include <Intervals.h>

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

#include <container_pybindings.h>
// Need cereal utility.hpp to encode pair<int,int>.
#include <cereal/types/utility.hpp>

//
// Default constructors, explicitly defined for each type, to set a
// sensible (perhaps) default domain.
//

// The G3Time internal encoding is an int64 with the number of 100 MHz
// ticks since the unix epoch.  Therefore...
#define THE_FUTURE_INT (4000000000L * 100000000L) // Sometime in the year 2096, in G3Time.

template <>
Intervals<double>::Intervals() {
    domain = make_pair(0., (double)THE_FUTURE_INT);
}

template <>
Intervals<int64_t>::Intervals() {
    domain = make_pair(0, THE_FUTURE_INT);
}

template <>
Intervals<G3Time>::Intervals() {
    domain = make_pair(G3Time(0), G3Time(THE_FUTURE_INT));
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
void Intervals<T>::cleanup()
{
    auto p = segments.begin();
    while (p != segments.end()) {
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
    // Expand the global domain, perhaps.
    domain.first  = min(domain.first,  start);
    domain.second = max(domain.second, end  );

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
    auto p0 = this->segments.begin();
    auto p1 = src.segments.begin();
    while (p0 != this->segments.end() && p1 != src.segments.end()) {
        // Which interval starts first?
        if (p0->first <= p1->first) {
            // p0 starts first.  Do they overlap?
            if (p1->first <= p0->second) {
                // Yes.  Increment p0 until we find the end of overlap.
                auto p0seg = p0+1;
                while (p0seg != segments.end() && p0seg->first <= p1->second)
                    ++p0seg;
                // Erase all but one of the overlapping intervals;
                // modify the remaining one.
                segments.erase(p0+1, p0seg);
                p0->second = max(p0->second, p1->second);
                ++p1;
            } else {
                // No. Move on, looking for the insertion point.
                ++p0;
            }
        } else {
            // p1 starts first.  Do they overlap?
            if (p0->first <= p1->second) {
                // Yes.
                p0->first = p1->first;
                p0->second = max(p0->second, p1->second);
            } else {
                // No.
                this->segments.insert(p0++, *(p1++));
            }
        }
    }
    while (p1 != src.segments.end())
        this->segments.push_back(*(p1++));

    return *this;
}

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
Intervals<T> Intervals<T>::from_array(const bp::numpy::ndarray &src)
{
    Intervals<T> output;
    assert(output.get_dtype() == src.get_dtype());

    char* d = src.get_data();
    auto strides = src.get_strides();
    auto shape = src.get_shape();
    assert(src.get_nd() == 2);
    for (int i=0; i < shape[0]; ++i) {
        output.segments.push_back(interval_pair<T>(d, d+strides[1]));
        d += strides[0];
    }
    return output;
}

template <typename T>
Intervals<T>& Intervals<T>::intersect(const Intervals<T> &src)
{
    auto output = -(*this);
    output -= src;
    *this = output.complement();
    return *this;
}
    
template <typename T>
void Intervals<T>::trim_to(T start, T end)
{
    // Once again, just be correct here... optimize later.
    auto p0 = segments.begin();
    while (p0 != segments.end()) {
        if (p0->second < start)
            p0->first = max(p0->first, start);
        p0++;
    }

    auto p1 = p0;
    while (p1 != segments.end()) {
        if (p1->second > end) {
            if (p1->first > end)
                break;
            (p1++)->second = end;
            break;
        }
    }

    segments.erase(p1, segments.end());
    segments.erase(segments.begin(), p0);
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

// Expose

template <typename T>
bp::numpy::dtype Intervals<T>::get_dtype() const
{
    return bp::numpy::dtype::get_builtin<T>();
}

template <>
bp::numpy::dtype Intervals<G3Time>::get_dtype() const
{
    return bp::numpy::dtype::get_builtin<G3TimeStamp>();
}


template <typename T>
bp::numpy::ndarray Intervals<T>::array() const
{
    bp::tuple shape = bp::make_tuple(segments.size(), 2);
    bp::numpy::dtype dt = get_dtype();
    bp::numpy::ndarray v = bp::numpy::empty(shape,dt);
    T* d = reinterpret_cast<T*>(v.get_data());
    for (auto p = segments.begin(); p != segments.end(); ++p) {
        *(d++) = p->first;
        *(d++) = p->second;
    }
    return v;
}

template <>
bp::numpy::ndarray Intervals<G3Time>::array() const
{
    bp::tuple shape = bp::make_tuple(segments.size(), 2);
    bp::numpy::dtype dt = get_dtype();
    bp::numpy::ndarray v = bp::numpy::empty(shape,dt);
    G3TimeStamp* d = reinterpret_cast<G3TimeStamp*>(v.get_data());
    for (auto p = segments.begin(); p != segments.end(); ++p) {
        *(d++) = p->first.time;
        *(d++) = p->second.time;
    }
    return v;
}


//
// Operators
//

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
    output.intersection(src);
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
    .def("add_interval", &CLASSNAME::add_interval, \
         return_internal_reference<>(), \
         "Merge an interval into the set.") \
    .def("merge", &CLASSNAME::merge, \
         return_internal_reference<>(), \
         "Merge an Intervals into the set.") \
    .def("intersect", &CLASSNAME::intersect, \
         return_internal_reference<>(), \
         "Intersect another " #DOMAIN_TYPE "with this one.") \
    .def("trim_to", &CLASSNAME::trim_to, \
         "Trim content to the specified range.") \
    .def("complement", &CLASSNAME::complement, \
         "Return the complement (over domain).") \
    .def("array", &CLASSNAME::array, \
         "Return the intervals as a 2-d numpy array.") \
    .def("from_array", &CLASSNAME::from_array, \
         "Return a " #DOMAIN_TYPE " based on an (n,2) ndarray.") \
    .staticmethod("from_array") \
    .def(-self) \
    .def(self += self) \
    .def(self -= self) \
    .def(self + self) \
    .def(self - self); \
    register_g3map<Map ## CLASSNAME>("Map" #CLASSNAME, "Mapping from " \
        "strings to Intervals over " #DOMAIN_TYPE ".")


G3_SERIALIZABLE_CODE(IntervalsFloat);
G3_SERIALIZABLE_CODE(IntervalsInt);
G3_SERIALIZABLE_CODE(IntervalsTime);

G3_SERIALIZABLE_CODE(MapIntervalsFloat);
G3_SERIALIZABLE_CODE(MapIntervalsInt);
G3_SERIALIZABLE_CODE(MapIntervalsTime);

PYBINDINGS("so3g")
{
    EXPORT_INTERVALS(double,  IntervalsFloat);
    EXPORT_INTERVALS(int64_t, IntervalsInt);
    EXPORT_INTERVALS(G3Time,  IntervalsTime);
}
