#define NO_IMPORT_ARRAY

#include <pybindings.h>
#include <iostream>
#include <G3Ndarray.h>
#include <boost/python.hpp>
#include <container_pybindings.h>
#include <cereal/types/utility.hpp>

G3Ndarray::G3Ndarray():data(NULL) {}
G3Ndarray::G3Ndarray(const G3Ndarray & src) {
    data = (PyArrayObject*) PyArray_FromAny((PyObject*)src.data, NULL, 0, 0, NPY_ARRAY_ENSUREARRAY, NULL);
}

// Construct G3Ndarray based on an array-like object.
G3Ndarray::G3Ndarray(const bp::object &bobject)
{
    PyObject *ob = PyArray_FromAny(bobject.ptr(), NULL, 0, 0, 0, NULL);
    if (ob == NULL)
        throw exception();

    data = reinterpret_cast<PyArrayObject*>(ob);
}

G3Ndarray::~G3Ndarray() {
    Py_XDECREF((PyObject*)data);
}

std::string G3Ndarray::Description() const
{
    if (data == NULL)
        return "G3Ndarray()";
    else {
        std::ostringstream s;
        s << "G3Ndarray(shape=(";
        for (auto i=0; i<PyArray_NDIM(data); i++) {
            if (i!=0)
                s << ",";
            s << PyArray_DIMS(data)[i];
        }
        s << "))";
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
    // SKN: I think this returns the type char (for example 'f'), which does not uniquely identify the
    // array type. There's something else called type_int (for example 7) which does uniquely identify
    // it. I think this is what simple_new wants.
    npy_intp type_num = PyArray_TYPE(data);
    ar & make_nvp("ndim", ndim);
    ar & make_nvp("type_num", type_num);
    auto dims = PyArray_DIMS(data);
    ar & make_nvp("shape", binary_data(dims, PyArray_NDIM(data)*sizeof(*dims)));
    long int size = 1;
    for(int i = 0; i < PyArray_NDIM(data); i++) size *= PyArray_DIM(data, i);
    // Copy over data into contiguous structure. It seems like the easiest way
    // to do this is to make a whole new numpy array
    PyArrayObject *contig = (PyArrayObject*) PyArray_NewCopy(data, NPY_CORDER);
    ar & make_nvp("data", binary_data((char*)PyArray_DATA(contig), size*PyArray_ITEMSIZE(contig)));
    Py_DECREF((PyObject*)contig);
}

template <class A> void G3Ndarray::load(A &ar, unsigned v) {
    using namespace cereal;
    // Load the basic frame object
    ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));
    // Load the information necessary to reconstruct numpy array
    npy_intp ndim, type_num;
    ar & make_nvp("ndim", ndim);
    ar & make_nvp("type_num", type_num);
    vector<npy_intp> shape(ndim);
    ar & make_nvp("shape", cereal::binary_data(&shape[0], ndim*sizeof(npy_intp)));
    long int size = 1;
    for(int i = 0; i < ndim; i++) size *= shape[i];
    // Make a new PyArrayObject with these properties
    Py_XDECREF(data);
    data = (PyArrayObject*) PyArray_SimpleNew(ndim, &shape[0], type_num);
    ar & make_nvp("data", binary_data((char*)PyArray_DATA(data), size*PyArray_ITEMSIZE(data)));
}

bp::object G3Ndarray::to_array() const {
    return bp::object(bp::handle<>(bp::borrowed(reinterpret_cast<PyObject*>(data))));
}

G3_SPLIT_SERIALIZABLE_CODE(G3Ndarray);

PYBINDINGS("so3g")
{
    using namespace boost::python;

    EXPORT_FRAMEOBJECT(G3Ndarray, init<>(), "G3Ndarray default constructor")
    .def(init<const bp::object&>("Construct G3Ndarray from numpy array"))
    .def("to_array", &G3Ndarray::to_array, "Get the wrapped numpy array")
    ;
}
