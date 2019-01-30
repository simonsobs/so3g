#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL Py_Array_API_SO3G

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
G3Ndarray::G3Ndarray(const G3Ndarray & src) {
    data = (PyArrayObject*) PyArray_FromAny((PyObject*)src.data, NULL, 0, 0, NPY_ARRAY_ENSUREARRAY, NULL);
}

G3Ndarray::G3Ndarray(bp::numpy::ndarray & barray)
{
    // Store a reference to the PyArrayObject.
    PyObject *p = barray.ptr();
    Py_INCREF(p);
    data = reinterpret_cast<PyArrayObject*>(p);
}

G3Ndarray::~G3Ndarray() {
    Py_XDECREF((PyObject*)data);
}

std::string G3Ndarray::Description() const
{
    if (data == NULL)
	return "G3Ndarray(dummy)";
    else {
	std::ostringstream s;
        double *x = reinterpret_cast<double*>(PyArray_DATA(data));
	s << "G3Ndarray(" << PyArray_NDIM(data) << ":" << x[0] <<  ")";
	return s.str();
    }
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
        npy_intp ndim = PyArray_NDIM(data);
        npy_intp dtype = PyArray_TYPE(data);
	ar & make_nvp("ndim",  ndim);
	ar & make_nvp("dtype", dtype);
        auto dims = PyArray_DIMS(data);
	ar & make_nvp("shape", binary_data(dims, PyArray_NDIM(data)*sizeof(*dims)));
	long int size = 1;
	for(int i = 0; i < PyArray_NDIM(data); i++) size *= PyArray_DIM(data, i);
	// Copy over data into contiguous structure. It seems like the easiest way
	// to do this is to make a whole new numpy array
	PyArrayObject *contig = (PyArrayObject*) PyArray_NewCopy(data, NPY_CORDER);
        ar & make_nvp("data", binary_data((char*)PyArray_DATA(contig),
                                          size*PyArray_DESCR(contig)->elsize));
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
	npy_intp shape[32];
	ar & make_nvp("shape", cereal::binary_data(&shape[0], ndim*sizeof(shape[0])));
	long int size = 1;
	for(int i = 0; i < ndim; i++) size *= shape[i];
        cout << ndim << ":" << shape[0] << " of " << typenum << " i.e. " << size << endl;
	// Make a new PyArrayObject with these properties
	Py_XDECREF(data);
	data = (PyArrayObject*) PyArray_SimpleNew(ndim, shape, typenum);
	ar & make_nvp("data", binary_data((char*)PyArray_DATA(data), size*PyArray_DESCR(data)->elsize));
}

G3_SPLIT_SERIALIZABLE_CODE(G3Ndarray);

using namespace boost::python;

PYBINDINGS("so3g")
{
    EXPORT_FRAMEOBJECT(G3Ndarray, init<>(), "G3Ndarray cow")
        .def(bp::init<numpy::ndarray&>("Capture ndarray as a G3Ndarray."))
        ;
}
