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

Spin0Weightor::~Spin0Weightor() {
    if (_handle != nullptr)
        Py_DecRef(_handle);
}

bool Spin0Weightor::TestInputs(bp::object map, bp::object weight)
{
    // Weights are always 1.  Map must be simple.
    BufferWrapper mapbuf;
    if (PyObject_GetBuffer(map.ptr(), &mapbuf.view,
                           PyBUF_FORMAT | PyBUF_ANY_CONTIGUOUS) == -1) {
        PyErr_Clear();
        throw buffer_exception("map");
    } 
    if (mapbuf.view.ndim != 3)
        throw shape_exception("map", "must have shape (n_y,n_x,n_map)");
    
    if (mapbuf.view.shape[2] != 1)
        throw shape_exception("map", "must have shape (n_y,n_x,1)");

    // Insist that user passed in None for the weights.
    if (weight.ptr() != Py_None) {
        throw shape_exception("weight", "must be None");
    }
    // The fully abstracted weights are shape (n_det n_time, n_map).
    // But in this specialization, we know n_map = 1, and the all
    // other elements should answer with weight 1.
    
    return true;
}

inline
bool Spin0Weightor::GetWeights(
    BufferWrapper &inline_weightbuf, double x, double y, double phi)
{
    if (_handle == nullptr) {
        // This is sketch-town.
        npy_intp dims[3] = {1,1,1};
        int dtype = NPY_FLOAT64;
        PyObject *v = PyArray_SimpleNew(3, dims, dtype);
        cout << "buffer: " << PyArray_DATA((PyArrayObject*)v) << endl;
        double *d = (double*)PyArray_DATA((PyArrayObject*)v);
        *d = 1.;
        PyObject_GetBuffer(v, &inline_weightbuf.view, PyBUF_FORMAT | PyBUF_ANY_CONTIGUOUS);
        _handle = v;
        inline_weightbuf.view.strides[0] = 0;
        inline_weightbuf.view.strides[1] = 0;
        inline_weightbuf.view.strides[2] = 0;
    }
    return true;
}

template<typename W>
GnomonicGridder<W>::GnomonicGridder()
{
    for (int i=0; i<2; ++i) {
        naxis[i] = 256;
        crpix[i] = 0;
        crval[i] = 0; 
        cdelt[i] = .0001;
    }
}

template<typename W>
bp::object GnomonicGridder<W>::zeros(
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

/** to_map(map, qpoint, signal, weights)
 *
 *  Full dimensionalities are:
 *     map:      (ny, nx, n_map)
 *     qpoint:   (n_det, n_t, n_coord)
 *     signal:   (n_sig, n_det, n_t)
 *     weight:   (n_sig, n_det, n_map)
 *
 *  Template over classes:
 *
 *  - Quaternion coords.
 *  - Quaternion boresight + offsets
 */

template<typename W>
bp::object GnomonicGridder<W>::to_map(
    bp::object map, bp::object qpoint, bp::object signal, bp::object weight)
{
    BufferWrapper mapbuf;
    if (PyObject_GetBuffer(map.ptr(), &mapbuf.view,
                           PyBUF_FORMAT | PyBUF_ANY_CONTIGUOUS) == -1) {
        PyErr_Clear();
        throw buffer_exception("map");
    } 

    if (mapbuf.view.ndim != 3)
        throw shape_exception("map", "must have shape (n_y,n_x,n_map)");

    BufferWrapper qpointbuf;
    if (PyObject_GetBuffer(qpoint.ptr(), &qpointbuf.view,
                           PyBUF_FORMAT | PyBUF_ANY_CONTIGUOUS) == -1) {
        PyErr_Clear();
        throw buffer_exception("qpoint");
    } 

    if (qpointbuf.view.ndim != 3)
        throw shape_exception("qpoint", "must have shape (n_det,n_t,n_coord)");

    BufferWrapper signalbuf;
    if (PyObject_GetBuffer(signal.ptr(), &signalbuf.view,
                           PyBUF_FORMAT | PyBUF_ANY_CONTIGUOUS) == -1) {
        PyErr_Clear();
        throw buffer_exception("signal");
    } 

    if (signalbuf.view.ndim != 3)
        throw shape_exception("signalbuf", "must have shape (n_sig,n_det,n_t)");

    BufferWrapper weightbuf;
    if (PyObject_GetBuffer(weight.ptr(), &weightbuf.view,
                           PyBUF_FORMAT | PyBUF_ANY_CONTIGUOUS) == -1) {
        PyErr_Clear();
        throw buffer_exception("weight");
    } 

    if (weightbuf.view.ndim != 3)
        throw shape_exception("weightbuf", "must have shape (n_sig,n_det,n_map)");

    int ny = mapbuf.view.shape[0];
    int nx = mapbuf.view.shape[1];
    int nmap = mapbuf.view.shape[2];
    int ndet = qpointbuf.view.shape[0];
    int nt = qpointbuf.view.shape[1];
    int ncoord = qpointbuf.view.shape[2];
    int nsig = signalbuf.view.shape[0];
    
    // Check that everything agrees...
    // ...
    
    // Update the map...
    double *mapd = (double*)mapbuf.view.buf;
    for (int idet=0; idet<ndet; ++idet) {
        const void *xy = ((char*)qpointbuf.view.buf + qpointbuf.view.strides[0] * idet);
        for (int it=0; it<nt; ++it) {
            double *x = (double*)((char*)xy + qpointbuf.view.strides[1] * it);
            double *y = (double*)((char*)xy + qpointbuf.view.strides[1] * it + qpointbuf.view.strides[2]);

            // Compute map pixel.
            double ix = (*x - crval[1]) / cdelt[1] + crpix[1] + 0.5;
            if (ix < 0 || ix >= naxis[1])
                continue;
            double iy = (*y - crval[0]) / cdelt[0] + crpix[0] + 0.5;
            if (iy < 0 || iy >= naxis[0])
                continue;
            
            int pix = int(ix);
            int piy = int(iy);
            for (int imap=0; imap<nmap; ++imap) {
                for (int isig=0; isig<nsig; ++isig) {
                    double sig = *(double*)((char*)signalbuf.view.buf +
                                            signalbuf.view.strides[0]*isig +
                                            signalbuf.view.strides[1]*idet +
                                            signalbuf.view.strides[2]*it);
                    double wt = *(double*)((char*)weightbuf.view.buf +
                                           weightbuf.view.strides[0]*isig +
                                           weightbuf.view.strides[1]*idet +
                                           weightbuf.view.strides[2]*imap);
                    //cout << pix << " " << piy << " " << sig << " " << wt << endl;
                    *(double*)((char*)mapbuf.view.buf +
                               mapbuf.view.strides[0]*piy +
                               mapbuf.view.strides[1]*pix +
                               mapbuf.view.strides[2]*imap) = sig*wt;
                }
            }
        }
    }

    return map;
}


/** to_map2(map, qpoint, qofs, signal, weights)
 *
 *  Full dimensionalities are:
 *     map:      (ny, nx, n_map)
 *     qbore:    (n_t, n_coord)
 *     qofs:     (n_det, n_coord)
 *     signal:   (n_sig, n_det, n_t)
 *     weight:   (n_sig, n_det, n_map)
 *
 *  Template over classes:
 *
 *  - Quaternion coords.
 *  - Quaternion boresight + offsets
 */

template<typename W>
bp::object GnomonicGridder<W>::to_map2(
    bp::object map, bp::object qbore, bp::object qofs, bp::object signal, bp::object weight)
{
    BufferWrapper mapbuf;
    if (PyObject_GetBuffer(map.ptr(), &mapbuf.view,
                           PyBUF_FORMAT | PyBUF_ANY_CONTIGUOUS) == -1) {
        PyErr_Clear();
        throw buffer_exception("map");
    } 

    if (mapbuf.view.ndim != 3)
        throw shape_exception("map", "must have shape (n_y,n_x,n_map)");

    BufferWrapper qborebuf;
    if (PyObject_GetBuffer(qbore.ptr(), &qborebuf.view,
                           PyBUF_FORMAT | PyBUF_ANY_CONTIGUOUS) == -1) {
        PyErr_Clear();
        throw buffer_exception("qbore");
    } 

    if (qborebuf.view.ndim != 2)
        throw shape_exception("qbore", "must have shape (n_t,n_coord)");

    BufferWrapper qofsbuf;
    if (PyObject_GetBuffer(qofs.ptr(), &qofsbuf.view,
                           PyBUF_FORMAT | PyBUF_ANY_CONTIGUOUS) == -1) {
        PyErr_Clear();
        throw buffer_exception("qofs");
    } 

    if (qofsbuf.view.ndim != 2)
        throw shape_exception("qofs", "must have shape (n_det,n_coord)");

    BufferWrapper signalbuf;
    if (PyObject_GetBuffer(signal.ptr(), &signalbuf.view,
                           PyBUF_FORMAT | PyBUF_ANY_CONTIGUOUS) == -1) {
        PyErr_Clear();
        throw buffer_exception("map");
    } 

    if (signalbuf.view.ndim != 3)
        throw shape_exception("signalbuf", "must have shape (n_sig,n_det,n_t)");

    // BufferWrapper weightbuf;
    // if (PyObject_GetBuffer(weight.ptr(), &weightbuf.view,
    //                        PyBUF_FORMAT | PyBUF_ANY_CONTIGUOUS) == -1) {
    //     PyErr_Clear();
    //     throw buffer_exception("map");
    // } 

    // if (weightbuf.view.ndim != 3)
    //     throw shape_exception("weightbuf", "must have shape (n_sig,n_det,n_map)");

    int ny = mapbuf.view.shape[0];
    int nx = mapbuf.view.shape[1];
    int nmap = mapbuf.view.shape[2];
    int nt = qborebuf.view.shape[0];
    int ncoord = qborebuf.view.shape[1];
    int ndet = qofsbuf.view.shape[0];
    int nsig = signalbuf.view.shape[0];
    
    // Check that everything agrees...
    // ...

    //Instantiate the weight-computing class.
    auto weightor = W();
    
    //Initialize it / check inputs.
    weightor.TestInputs(map, weight);

    BufferWrapper inline_weightbuf;
    
    // Update the map...
    double *mapd = (double*)mapbuf.view.buf;
    for (int idet=0; idet<ndet; ++idet) {
        const double dx = *(double*)((char*)qofsbuf.view.buf +
                                     qofsbuf.view.strides[0] * idet);
        const double dy = *(double*)((char*)qofsbuf.view.buf +
                                     qofsbuf.view.strides[0] * idet + 
                                     qofsbuf.view.strides[1]);
        for (int it=0; it<nt; ++it) {
            double *x = (double*)((char*)qborebuf.view.buf +
                                  qborebuf.view.strides[0] * it);
            double *y = (double*)((char*)qborebuf.view.buf +
                                  qborebuf.view.strides[0] * it +
                                  qborebuf.view.strides[1]);

            // Compute map pixel.
            double ix = (*x + dx - crval[1]) / cdelt[1] + crpix[1] + 0.5;
            if (ix < 0 || ix >= naxis[1])
                continue;
            double iy = (*y + dy - crval[0]) / cdelt[0] + crpix[0] + 0.5;
            if (iy < 0 || iy >= naxis[0])
                continue;
            
            int pix = int(ix);
            int piy = int(iy);

            // Now the Weightor class will give us the right weights to use.
            weightor.GetWeights(inline_weightbuf, ix, iy, 0.);
            
            for (int imap=0; imap<nmap; ++imap) {
                for (int isig=0; isig<nsig; ++isig) {
                    double sig = *(double*)((char*)signalbuf.view.buf +
                                            signalbuf.view.strides[0]*isig +
                                            signalbuf.view.strides[1]*idet +
                                            signalbuf.view.strides[2]*it);
                    // hot-wire simple case.
                    //double wt = *(double*)((char*)inline_weightbuf.view.buf);
                    // if (it == 0) {
                    //     cout << "bufferX: " << inline_weightbuf.view.buf << endl;
                    //     cout << wt << endl;
                    // }
                    double wt = *(double*)((char*)inline_weightbuf.view.buf +
                                           inline_weightbuf.view.strides[0]*isig +
                                           inline_weightbuf.view.strides[1]*idet +
                                           inline_weightbuf.view.strides[2]*imap);
                    //cout << pix << " " << piy << " " << sig << " " << wt << endl;
                    *(double*)((char*)mapbuf.view.buf +
                               mapbuf.view.strides[0]*piy +
                               mapbuf.view.strides[1]*pix +
                               mapbuf.view.strides[2]*imap) = sig*wt;
                }
            }
        }
    }

    return map;
}

typedef GnomonicGridder<Weightor> GnomonicGridderX;
typedef GnomonicGridder<Spin0Weightor> GnomonicGridder0;

PYBINDINGS("so3g")
{
    bp::class_<GnomonicGridderX>("GnomonicGridderBase");
    bp::class_<GnomonicGridder0>("GnomonicGridder")
        .def("zeros", &GnomonicGridder0::zeros)
        .def("to_map", &GnomonicGridder0::to_map)
        .def("to_map2", &GnomonicGridder0::to_map2);
}
