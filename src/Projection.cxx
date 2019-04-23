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

inline bool isNone(bp::object &pyo)
{
    return (pyo.ptr() == Py_None);
}

void Pointer::Init(BufferWrapper &qborebuf, BufferWrapper &qofsbuf)
{
    _qborebuf = &qborebuf;
    _qofsbuf = &qofsbuf;
}

inline
void Pointer::InitPerDet(int idet)
{
    _dx = *(double*)((char*)_qofsbuf->view.buf +
                                 _qofsbuf->view.strides[0] * idet);
    _dy = *(double*)((char*)_qofsbuf->view.buf +
                                 _qofsbuf->view.strides[0] * idet + 
                                 _qofsbuf->view.strides[1]);
}

inline
void Pointer::GetCoords(int idet, int it, double *coords)
{
    coords[0] = _dx + *(double*)((char*)_qborebuf->view.buf +
                                 _qborebuf->view.strides[0] * it);
    coords[1] = _dy + *(double*)((char*)_qborebuf->view.buf +
                                 _qborebuf->view.strides[0] * it +
                                 _qborebuf->view.strides[1]);
}


Pixelizor::Pixelizor()
{
    for (int i=0; i<2; ++i) {
        naxis[i] = 256;
        crpix[i] = 0;
        crval[i] = 0; 
        cdelt[i] = .0001;
    }
}

void Pixelizor::Init(BufferWrapper &mapbuf)
{
    _mapbuf = &mapbuf;
}


bp::object Pixelizor::zeros(bp::object shape)
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


inline
int Pixelizor::GetPixel(int i_det, int i_t, const double *coords)
{
    double ix = (coords[0] - crval[1]) / cdelt[1] + crpix[1] + 0.5;
    if (ix < 0 || ix >= naxis[1])
        return -1;

    double iy = (coords[1] - crval[0]) / cdelt[0] + crpix[0] + 0.5;
    if (iy < 0 || iy >= naxis[0])
        return -1;
            
    int pixel_offset = _mapbuf->view.strides[0]*int(iy) +
        _mapbuf->view.strides[1]*int(ix);

    return pixel_offset;
}


/** Accumulator - transfer signal from map domain to time domain.
 *
 */

AccumulatorSpin0::~AccumulatorSpin0() {
    if (_handle != nullptr)
        Py_DecRef(_handle);
}

bool AccumulatorSpin0::TestInputs(bp::object map, bp::object signal,
                                  bp::object weight)
{
    // Check that map and signal have 1d output.
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
    if (!isNone(weight)) {
        throw shape_exception("weight", "must be None");
    }
    // The fully abstracted weights are shape (n_det n_time, n_map).
    // But in this specialization, we know n_map = 1, and the all
    // other elements should answer with weight 1.
    return true;
}

void AccumulatorSpin0::Init(BufferWrapper &inline_weightbuf)
{
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

inline
void AccumulatorSpin0::Forward(const BufferWrapper &inline_weightbuf,
                               const BufferWrapper &signalbuf,
                               BufferWrapper &mapbuf,
                               const int idet,
                               const int it,
                               const double* coords,
                               const int pixel_offset)
{
    // double wt = *(double*)((char*)inline_weightbuf.view.buf +
    //                        inline_weightbuf.view.strides[1]*idet);
    double sig = *(double*)((char*)signalbuf.view.buf +
                            signalbuf.view.strides[1]*idet +
                            signalbuf.view.strides[2]*it);
    *(double*)((char*)mapbuf.view.buf + pixel_offset) += sig;
}

AccumulatorSpin2::~AccumulatorSpin2() {
    if (_handle != nullptr)
        Py_DecRef(_handle);
}

bool AccumulatorSpin2::TestInputs(bp::object map, bp::object signal,
                                  bp::object weight)
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
    
    if (mapbuf.view.shape[2] != 3)
        throw shape_exception("map", "must have shape (n_y,n_x,3)");

    // Insist that user passed in None for the weights.
    if (weight.ptr() != Py_None) {
        throw shape_exception("weight", "must be None");
    }

    // Weights will be computed from input coordinates, esp 2phi.
    
    return true;
}

void AccumulatorSpin2::Init(BufferWrapper &inline_weightbuf)
{
    if (_handle == nullptr) {
        // This is sketch-town.
        npy_intp dims[3] = {1,1,3};
        int dtype = NPY_FLOAT64;
        PyObject *v = PyArray_SimpleNew(3, dims, dtype);
        cout << "buffer: " << PyArray_DATA((PyArrayObject*)v) << endl;
        double *d = (double*)PyArray_DATA((PyArrayObject*)v);
        d[0] = 1.;
        PyObject_GetBuffer(v, &inline_weightbuf.view, PyBUF_FORMAT | PyBUF_ANY_CONTIGUOUS);
        _handle = v;
        inline_weightbuf.view.strides[0] = 0;
        inline_weightbuf.view.strides[1] = 0;
        inline_weightbuf.view.strides[2] = sizeof(*d);
    }
}

inline
void AccumulatorSpin2::Forward(const BufferWrapper &inline_weightbuf,
                               const BufferWrapper &signalbuf,
                               BufferWrapper &mapbuf,
                               const int idet,
                               const int it,
                               const double* coords,
                               const int pixel_offset)
{
    double sig = *(double*)((char*)signalbuf.view.buf +
                            signalbuf.view.strides[1]*idet +
                            signalbuf.view.strides[2]*it);
    for (int imap=0; imap<3; ++imap) {
        const double wt = *(double*)((char*)inline_weightbuf.view.buf +
                               inline_weightbuf.view.strides[1]*idet +
                               inline_weightbuf.view.strides[2]*imap);
        *(double*)((char*)mapbuf.view.buf +
            mapbuf.view.strides[2]*imap +
            pixel_offset) += sig * wt;
    }
}


template<typename P, typename Z, typename A>
bp::object ProjectionEngine<P,Z,A>::zeros(
    bp::object shape)
{
    Pixelizor p;
    return p.zeros(shape);
}


/** to_map(map, qpoint, qofs, signal, weights)
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

template<typename P, typename Z, typename A>
bp::object ProjectionEngine<P,Z,A>::to_map(
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
    auto pointer = P();
    auto pixelizor = Z();
    auto accumulator = A();

    //Initialize it / check inputs.
    //pointer.TestInputs(...);
    //pixelizor.TestInputs(...);
    accumulator.TestInputs(map, signal, weight);
    
    BufferWrapper inline_weightbuf;
    pointer.Init(qborebuf, qofsbuf);
    pixelizor.Init(mapbuf);
    accumulator.Init(inline_weightbuf);
    
    // Update the map...
    double *mapd = (double*)mapbuf.view.buf;
    double coords[4];

    for (int idet=0; idet<ndet; ++idet) {
        pointer.InitPerDet(idet);
        for (int it=0; it<nt; ++it) {
            pointer.GetCoords(idet, it, (double*)coords);
            int pixel_offset = pixelizor.GetPixel(idet, it, (double*)coords);
            if (pixel_offset < 0)
                continue;

            // Now the Weightor class will give us the right weights to use.
            accumulator.Forward(inline_weightbuf, signalbuf, mapbuf,
                                idet, it, coords, pixel_offset);
        }
    }

    return map;
}

typedef ProjectionEngine<Pointer,Pixelizor,AccumulatorSpin0> ProjectionEngine0;
typedef ProjectionEngine<Pointer,Pixelizor,AccumulatorSpin2> ProjectionEngine2;

PYBINDINGS("so3g")
{
    bp::class_<ProjectionEngine0>("ProjectionEngine")
        .def("zeros", &ProjectionEngine0::zeros)
        .def("to_map", &ProjectionEngine0::to_map);
    bp::class_<ProjectionEngine2>("ProjectionEngine2")
        .def("zeros", &ProjectionEngine2::zeros)
        .def("to_map", &ProjectionEngine2::to_map);
}
