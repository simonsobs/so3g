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
       G3Ndarray(G3Ndarray &);
       G3Ndarray(bp::numpy::ndarray &);
       ~G3Ndarray();
       template <class A> void save(A &ar, unsigned v) const;
       template <class A> void load(A &ar, unsigned v);
       std::string Description() const;
       PyArrayObject * data;
};

G3_SERIALIZABLE(G3Ndarray, 0);
