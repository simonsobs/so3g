#pragma once

#include <G3Frame.h>
#include <G3Map.h>
#include <G3TimeStamp.h>
#include <boost/python/numpy.hpp>
#include <numpy/arrayobject.h>

using namespace std;
namespace bp = boost::python;

class G3Ndarray : public G3FrameObject {
public:
    G3Ndarray();
    G3Ndarray(const G3Ndarray &);
    G3Ndarray(const bp::numpy::ndarray &);
    ~G3Ndarray();
    template <class A> void save(A &ar, unsigned v) const;
    template <class A> void load(A &ar, unsigned v);
    std::string Description() const;
    bp::numpy::ndarray to_array() const;
    PyArrayObject *data;
};

namespace cereal {
    template <class A> struct specialize<A, G3Ndarray, cereal::specialization::member_load_save> {};
}

G3_SERIALIZABLE(G3Ndarray, 0);
G3_POINTERS(G3Ndarray);
