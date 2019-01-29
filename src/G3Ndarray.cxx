#include <pybindings.h>

#include <iostream>
#include <G3Ndarray.h>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <numpy/arrayobject.h>
#include <container_pybindings.h>
// Need cereal utility.hpp to encode pair<int,int>.
#include <cereal/types/utility.hpp>

G3Ndarray::G3Ndarray():data(NULL) {}
G3Ndarray::G3Ndarray(G3Ndarray & src) {
	data = (PyArrayObject*) PyArray_FromAny((PyObject*)src.data, NULL, 0, 0, NPY_ARRAY_ENSUREARRAY, NULL);
}
G3Ndarray::G3Ndarray(bp::numpy::ndarray & barray) {
	data = (PyArrayObject*) PyArray_FromAny((PyObject*)barray.ptr(), NULL, 0, 0, NPY_ARRAY_ENSUREARRAY, NULL);
}
G3Ndarray::~G3Ndarray() { Py_XDECREF((PyObject*)data); }

std::string G3Ndarray::Description() const
{
	return "G3Ndarray(dummy)";
}

// Serialization: save dtype type-int, ndim, shape, strides and then data
// loading: load type, ndim, shape, strides. Update our properties. Then load data
// But how does the update part look?

template <class A> void G3Ndarray::save(A &ar, unsigned v) const // v is the version code
{
	using namespace cereal;
	// Save the basic frame object
	ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));
	// Then handle the actual ndarray.
	ar & make_nvp("ndim",  PyArray_NDIM(data));
	ar & make_nvp("dtype", PyArray_DESCR(data)->type);
	ar & make_nvp("shape", binary_data(PyArray_DIMS(data), PyArray_NDIM(data)));
	long int size = 1;
	for(int i = 0; i < PyArray_NDIM(data); i++) size *= PyArray_DIM(data, i);
	// Copy over data into contiguous structure. It seems like the easiest way
	// to do this is to make a whole new numpy array
	PyArrayObject * contig = (PyArrayObject*) PyArray_NewCopy(data, NPY_CORDER);
	ar & make_nvp("data", binary_data((char*)PyArray_DATA(contig), size*PyArray_DESCR(data)->elsize));
	Py_DECREF((PyObject*)contig);
}

template <class A> void G3Ndarray::load(A &ar, unsigned v) {
	using namespace cereal;
	// Load the basic frame object
	ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));
	// Load the information necessary to reconstruct numpy array
	npy_intp ndim, typenum;
	ar & make_nvp("ndim",  ndim);
	ar & make_nvp("dtype", typenum);
	vector<npy_intp> shape(ndim), stride(ndim, 1);
	ar & make_nvp("shape", cereal::binary_data(&shape[0], ndim));
	long int size = 1;
	for(int i = 0; i < ndim; i++) size *= shape[i];
	// Make a new PyArrayObject with these properties
	Py_XDECREF(data);
	data = (PyArrayObject*) PyArray_SimpleNew(ndim, &shape[0], typenum);
	ar & make_nvp("data", binary_data((char*)PyArray_DATA(data), size*PyArray_DESCR(data)->elsize));
}

G3_SPLIT_SERIALIZABLE_CODE(G3Ndarray);

using namespace boost::python;

//PYBINDINGS("so3g")
//{
//	EXPORT_FRAMEOBJECT("
//
//    EXPORT_INTERVALS(double,  IntervalsFloat);
//    EXPORT_INTERVALS(int64_t, IntervalsInt);
//    EXPORT_INTERVALS(G3Time,  IntervalsTime);
//}

//#endif
//
//
//template <class A> void G3Ndarray::save(A &ar, unsigned v) const
//{
//}
//
//template <class A> void G3Ndarray::load(A &ar, unsigned v)
//{
//}
//
////template void G3Ndarray::save(cereal::PortableBinaryOutputArchive &, unsigned) const;
////template void G3Ndarray::load(cereal::PortableBinaryInputArchive  &, unsigned);
//
//
PYBINDINGS("so3g")
{
	EXPORT_FRAMEOBJECT(G3Ndarray, init<>(), "G3Ndarray cow");
}
