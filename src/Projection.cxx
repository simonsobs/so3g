#define NO_IMPORT_ARRAY

// debug
#include <iostream>
using namespace std;

#include <pybindings.h>

#include <assert.h>
#include <math.h>

#include <container_pybindings.h>

#include "so3g_numpy.h"
//#include "Python.h"

#include <Projection.h>
#include "exceptions.h"

GnomonicGridder::GnomonicGridder()
{
    for (int i=0; i<2; ++i) {
        naxis[i] = 256;
        crpix[i] = 0;
        crval[i] = 0; 
        cdelt[i] = .0001;
    }
}

bp::object GnomonicGridder::zeros(
    bp::object shape)
{
    int size = 1;
    int dimi = 0;
    npy_intp dims[32];
    if (shape == bp::object()) {
    } else {
        // add shape...
    }
    dims[dimi++] = naxis[1];
    dims[dimi++] = naxis[0];
    size *= naxis[0] * naxis[1];
    int dtype = NPY_FLOAT64;
    PyObject *v = PyArray_SimpleNew(dimi, dims, dtype);
    memset(PyArray_DATA((PyArrayObject*)v), 0, size * PyArray_ITEMSIZE((PyArrayObject*)v));
    return bp::object(bp::handle<>(v));
}

bp::object GnomonicGridder::to_map(
    bp::object map, bp::object qpoint, bp::object signal, bp::object weights)
{
    if (map.ptr() == Py_None) {
        // New zero'd map.
        map = zeros(bp::object());
    } 

    if (qpoint.ptr() == Py_None)
        return map;

    // qpoint is just X,Y believe it or not.
    BufferWrapper qbuf;
    if (PyObject_GetBuffer(qpoint.ptr(), &qbuf.view,
                           PyBUF_FORMAT | PyBUF_ANY_CONTIGUOUS) == -1) {
        PyErr_Clear();
        throw buffer_exception("qpoint");
    } 

    // Buffer access to the map.
    BufferWrapper mapbuf;
    if (PyObject_GetBuffer(map.ptr(), &mapbuf.view,
                           PyBUF_FORMAT | PyBUF_ANY_CONTIGUOUS) == -1) {
        PyErr_Clear();
        throw buffer_exception("map");
    } 

    const double *xy = (double*)qbuf.view.buf;
    const double *x = xy;
    const double *y = xy + qbuf.view.strides[0] / sizeof(*y);
    
    // Update the map...
    double *mapd = (double*)mapbuf.view.buf;

    for (int i=0; i<qbuf.view.shape[1]; i++) {
        double ix = (x[i] - crval[1]) / cdelt[1] + crpix[1] + 0.5;
        if (ix < 0 || ix >= naxis[0])
            continue;
        double iy = (y[i] - crval[0]) / cdelt[0] + crpix[0] + 0.5;
        if (iy < 0 || iy >= naxis[1])
            continue;
        cout << ix << " " << iy << endl;
        mapd[(int)ix +(int)iy*naxis[0]] += 1.;
    }
    
    return map;
}


PYBINDINGS("so3g")
{
    bp::class_<GnomonicGridder>("GnomonicGridder")
        .def("zeros", &GnomonicGridder::zeros)
        .def("to_map", &GnomonicGridder::to_map);
}
