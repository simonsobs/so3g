#define NO_IMPORT_ARRAY

// debug
#include <iostream>
using namespace std;

#include <pybindings.h>

#include <assert.h>
#include <math.h>

#include <container_pybindings.h>

#include "so3g_numpy.h"
#include <Projection.h>
#include "exceptions.h"

inline bool isNone(bp::object &pyo)
{
    return (pyo.ptr() == Py_None);
}

bool PointerP_Simple_Flat::TestInputs(
    bp::object &map, bp::object &pbore, bp::object &pdet,
    bp::object &signal, bp::object &weight)
{
    // Boresight and Detector must present and inter-compatible.

    if (PyObject_GetBuffer(pbore.ptr(), &_pborebuf.view,
                           PyBUF_RECORDS) == -1) {
        PyErr_Clear();
        throw buffer_exception("pbore");
    }
    if (PyObject_GetBuffer(pdet.ptr(), &_pdetbuf.view,
                           PyBUF_RECORDS) == -1) {
        PyErr_Clear();
        throw buffer_exception("pdet");
    }

    if (_pborebuf.view.ndim != 2)
        throw shape_exception("pbore", "must have shape (n_t,n_coord)");
    if (_pdetbuf.view.ndim != 2)
        throw shape_exception("pdet", "must have shape (n_det,n_coord)");

    if (_pborebuf.view.shape[1] != 4)
        throw shape_exception("pbore", "must have shape (n_t,4)");
    if (_pdetbuf.view.shape[1] != 4)
        throw shape_exception("pdet", "must have shape (n_det,4)");

    n_time = _pborebuf.view.shape[0];
    n_det = _pdetbuf.view.shape[0];

    return true;
}

inline
void PointerP_Simple_Flat::InitPerDet(int i_det)
{
    const char *det = (char*)_pdetbuf.view.buf
        + _pdetbuf.view.strides[0] * i_det;
    for (int ic = 0; ic < 4; ++ic)
        _coords[ic] = *(double*)(det + _pdetbuf.view.strides[1] * ic);
}

inline
void PointerP_Simple_Flat::GetCoords(int i_det, int i_time, double *coords)
{
    for (int ic=0; ic<4; ic++)
        coords[ic] = *(double*)((char*)_pborebuf.view.buf +
                                _pborebuf.view.strides[0] * i_time +
                                _pborebuf.view.strides[1] * ic);
    coords[0] += _coords[0];
    coords[1] += _coords[1];
    const double coords_2_ = coords[2];
    coords[2] = coords[2] * _coords[2] - coords[3] * _coords[3];
    coords[3] = coords[3] * _coords[2] + coords_2_ * _coords[3];
}

Pixelizor2_Flat::Pixelizor2_Flat(
    int ny, int nx,
    double dy, double dx,
    double y0, double x0,
    double iy0, double ix0)
{
    naxis[0] = ny;
    naxis[1] = nx;
    cdelt[0] = dy;
    cdelt[1] = dx;
    crval[0] = y0;
    crval[1] = x0;
    crpix[0] = iy0;
    crpix[1] = ix0;

    // These will be set in context. 
    strides[0] = 0;
    strides[1] = 0;
}

bool Pixelizor2_Flat::TestInputs(bp::object &map, bp::object &pbore, bp::object &pdet,
                             bp::object &signal, bp::object &weight)
{
    // Flat pixelizor.  Requires:
    // 1. Map has the right shape -- but map can be none.
    // 2. Coords have the right format (unchecked currently).

    if (!isNone(map)) {
        BufferWrapper mapbuf;
        if (PyObject_GetBuffer(map.ptr(), &mapbuf.view,
                               PyBUF_RECORDS) == -1) {
            PyErr_Clear();
            throw buffer_exception("map");
        }
        int ndim = mapbuf.view.ndim;
        if (mapbuf.view.ndim < 2)
            throw shape_exception("map", "must have shape (...,n_y,n_x)");
        if (mapbuf.view.shape[ndim-2] != naxis[0])
            throw shape_exception("map", "dimension -2 must match naxis[0]");
        if (mapbuf.view.shape[ndim-1] != naxis[1])
            throw shape_exception("map", "dimension -1 must match naxis[1]");

        // Note these are byte offsets, not index.
        strides[0] = mapbuf.view.strides[ndim-2];
        strides[1] = mapbuf.view.strides[ndim-1];
    } else {
        // Set it up to return naive C-ordered pixel indices.
        strides[0] = naxis[1];
        strides[1] = 1;
    }

    // [Check/declare coords format.]
    return true;
}

bp::object Pixelizor2_Flat::zeros(int count)
{
    int size = 1;
    int dimi = 0;
    npy_intp dims[32];

    if (count >= 0) {
        dims[dimi++] = count;
        size *= count;
    }

    dims[dimi++] = naxis[0];
    dims[dimi++] = naxis[1];
    size *= naxis[0] * naxis[1];

    int dtype = NPY_FLOAT64;
    PyObject *v = PyArray_ZEROS(dimi, dims, dtype, 0);
    return bp::object(bp::handle<>(v));
}


inline
int Pixelizor2_Flat::GetPixel(int i_det, int i_time, const double *coords)
{
    double ix = (coords[0] - crval[1]) / cdelt[1] + crpix[1] + 0.5;
    if (ix < 0 || ix >= naxis[1])
        return -1;

    double iy = (coords[1] - crval[0]) / cdelt[0] + crpix[0] + 0.5;
    if (iy < 0 || iy >= naxis[0])
        return -1;
            
    int pixel_offset = strides[0]*int(iy) + strides[1]*int(ix);

    return pixel_offset;
}


/** Accumulator - transfer signal from map domain to time domain.
 *
 */

template <typename SpinClass>
bool Accumulator<SpinClass>::TestInputs(
    bp::object &map, bp::object &pbore, bp::object &pdet,
    bp::object &signal, bp::object &weight)
{
    const int N = SpinClass::comp_count;
    if (need_map) {
        if (PyObject_GetBuffer(map.ptr(), &_mapbuf.view,
                               PyBUF_RECORDS) == -1) {
            PyErr_Clear();
            throw buffer_exception("map");
        }
        if (_mapbuf.view.ndim < 2)
            throw shape_exception("map", "must have shape (n_map,n_axis0,...)");

        if (_mapbuf.view.shape[0] != N)
            throw shape_exception("map", "must have shape (n_comp,n_axis0,...)");
    } else if (need_weight_map) {
        if (PyObject_GetBuffer(map.ptr(), &_mapbuf.view,
                               PyBUF_RECORDS) == -1) {
            PyErr_Clear();
            throw buffer_exception("map");
        }
        if (_mapbuf.view.ndim < 3)
            throw shape_exception("map", "must have shape (n_map,n_axis0,...)");
        if (_mapbuf.view.shape[0] != N ||
            _mapbuf.view.shape[1] != N)
            throw shape_exception("map", "must have shape (n_comp,n_comp,n_axis0,...)");
    }

    if (need_signal) {
        if (PyObject_GetBuffer(signal.ptr(), &_signalbuf.view,
                               PyBUF_RECORDS) == -1) {
            PyErr_Clear();
            throw buffer_exception("signal");
        }
        if (_signalbuf.view.ndim != 3)
            throw shape_exception("signal", "must have shape (n_sig,n_det,n_t)");
    }

    // Insist that user passed in None for the weights.
    if (weight.ptr() != Py_None) {
        throw shape_exception("weight", "must be None");
    }

    return true;
}

template <>
inline
void Accumulator<SpinT>::PixelWeight(
    const double* coords, double *pwt)
{
    pwt[0] = 1;
}

template <>
inline
void Accumulator<SpinQU>::PixelWeight(
    const double* coords, double *pwt)
{
    const double c = coords[2];
    const double s = coords[3];
    pwt[0] = c*c - s*s;
    pwt[1] = 2*c*s;
}

template<>
inline
void Accumulator<SpinTQU>::PixelWeight(
    const double* coords, double *pwt)
{
    const double c = coords[2];
    const double s = coords[3];
    pwt[0] = 1.;
    pwt[1] = c*c - s*s;
    pwt[2] = 2*c*s;
}

template <typename SpinClass>
inline
void Accumulator<SpinClass>::Forward(
    const int i_det, const int i_time,
    const int pixel_offset, const double* coords, const double* weights)
{
    if (pixel_offset < 0) return;
    const int N = SpinClass::comp_count;
    double sig = *(double*)((char*)_signalbuf.view.buf +
                            _signalbuf.view.strides[1]*i_det +
                            _signalbuf.view.strides[2]*i_time);
    double wt[N];
    PixelWeight(coords, wt);
    for (int imap=0; imap<N; ++imap) {
        *(double*)((char*)_mapbuf.view.buf +
                   _mapbuf.view.strides[0]*imap +
                   pixel_offset) += sig * wt[imap];
    }
}

template <typename SpinClass>
inline
void Accumulator<SpinClass>::ForwardWeight(
    const int i_det, const int i_time,
    const int pixel_offset, const double* coords, const double* weights)
{
    if (pixel_offset < 0) return;
    const int N = SpinClass::comp_count;
    double wt[N];
    PixelWeight(coords, wt);
    for (int imap=0; imap<N; ++imap) {
        for (int jmap=imap; jmap<N; ++jmap) {
            *(double*)((char*)_mapbuf.view.buf +
                       _mapbuf.view.strides[0]*imap +
                       _mapbuf.view.strides[1]*jmap +
                       pixel_offset) += wt[imap] * wt[jmap];
        }
    }
}

template <typename SpinClass>
inline
void Accumulator<SpinClass>::Reverse(
    const int i_det, const int i_time,
    const int pixel_offset, const double* coords, const double* weights)
{
    if (pixel_offset < 0) return;
    const int N = SpinClass::comp_count;
    double wt[N];
    PixelWeight(coords, wt);
    double _sig = 0.;
    for (int imap=0; imap<N; ++imap) {
        _sig += *(double*)((char*)_mapbuf.view.buf +
                           _mapbuf.view.strides[0]*imap +
                           pixel_offset) * wt[imap];
    }
    double *sig = (double*)((char*)_signalbuf.view.buf +
                            _signalbuf.view.strides[1]*i_det +
                            _signalbuf.view.strides[2]*i_time);
    *sig += _sig;
}


/** to_map(map, qpoint, pofs, signal, weights)
 *
 *  Each argument is an ndarray.  In the general case the dimensionalities are:
 *
 *     map:      (n_map, ny, nx, ...)
 *     pbore:    (n_t, n_coord)
 *     pofs:     (n_det, n_coord)
 *     signal:   (n_det, n_t)
 *     weight:   (n_sig, n_det, n_map)
 *
 *  Notes:
 *
 *  - The map dimensions (ny, nx, ...) are meant to be general enough
 *    to also capture 1-dimensional pixelization systems (such as
 *    healpix), or higher dimensional grids, e.g., (nz, ny, nx).
 *
 *  - Normally we will have n_coord=4.  In the special case of a Flat
 *    (2-d cartesian) projection space with no detector orientation
 *    information, n_coord=2.
 *
 */

template<typename P, typename Z, typename A>
ProjectionEngine<P,Z,A>::ProjectionEngine(Z pixelizor)
{
    _pixelizor = pixelizor;
}

template<typename P, typename Z, typename A>
bp::object ProjectionEngine<P,Z,A>::to_map(
    bp::object map, bp::object pbore, bp::object pofs, bp::object signal, bp::object weight)
{
    auto pointer = P();
    auto accumulator = A(true, true, false);

    //Initialize it / check inputs.
    pointer.TestInputs(map, pbore, pofs, signal, weight);
    int n_det = pointer.DetCount();
    int n_time = pointer.TimeCount();

    //Do we need a map?  Now is the time.
    if (isNone(map)) {
        int n_comp = accumulator.ComponentCount();
        map = _pixelizor.zeros(n_comp);
    }

    _pixelizor.TestInputs(map, pbore, pofs, signal, weight);
    accumulator.TestInputs(map, pbore, pofs, signal, weight);

    for (int i_det = 0; i_det < n_det; ++i_det) {
        pointer.InitPerDet(i_det);
        for (int i_time = 0; i_time < n_time; ++i_time) {
            double coords[4];
            double weights[4];
            int pixel_offset;
            pointer.GetCoords(i_det, i_time, (double*)coords);
            pixel_offset = _pixelizor.GetPixel(i_det, i_time, (double*)coords);
            accumulator.Forward(i_det, i_time, pixel_offset, coords, weights);
        }
    }

    return map;
}

template<typename P, typename Z, typename A>
bp::object ProjectionEngine<P,Z,A>::to_weight_map(
    bp::object map, bp::object pbore, bp::object pofs, bp::object signal, bp::object weight)
{
    auto pointer = P();
    auto accumulator = A(false, false, true);

    //Initialize it / check inputs.
    pointer.TestInputs(map, pbore, pofs, signal, weight);
    int n_det = pointer.DetCount();
    int n_time = pointer.TimeCount();

    //Do we need a map?  Now is the time.
    if (isNone(map)) {
        int n_comp = accumulator.ComponentCount();
        map = _pixelizor.zeros(n_comp * n_comp);
        auto v0 = (PyArrayObject*)map.ptr();
        npy_intp dims[32] = {n_comp, n_comp};
        int dimi = 2;
        for (int d=1; d<PyArray_NDIM(v0); d++)
            dims[dimi++] = PyArray_DIM(v0, d);
        PyArray_Dims padims = {dims, dimi};
        PyObject *v1 = PyArray_Newshape(v0, &padims,  NPY_ANYORDER);
        map = bp::object(bp::handle<>(v1));
    }

    _pixelizor.TestInputs(map, pbore, pofs, signal, weight);
    accumulator.TestInputs(map, pbore, pofs, signal, weight);

    for (int i_det = 0; i_det < n_det; ++i_det) {
        pointer.InitPerDet(i_det);
        for (int i_time = 0; i_time < n_time; ++i_time) {
            double coords[4];
            double weights[4];
            int pixel_offset;
            pointer.GetCoords(i_det, i_time, (double*)coords);
            pixel_offset = _pixelizor.GetPixel(i_det, i_time, (double*)coords);
            accumulator.ForwardWeight(i_det, i_time, pixel_offset, coords, weights);
        }
    }

    return map;
}

template<typename P, typename Z, typename A>
bp::object ProjectionEngine<P,Z,A>::from_map(
    bp::object map, bp::object pbore, bp::object pofs, bp::object signal, bp::object weight)
{
    auto pointer = P();
    auto accumulator = A(true, true, false);

    // Initialize pointer and _pixelizor.
    pointer.TestInputs(map, pbore, pofs, signal, weight);

    int n_det = pointer.DetCount();
    int n_time = pointer.TimeCount();
    
    // Do we have a signal array?  Now is the time.
    if (isNone(signal)) {
        npy_intp dims[3] = {1, n_det, n_time};
        PyObject *v = PyArray_ZEROS(3, dims, NPY_FLOAT64, 0);
        signal = bp::object(bp::handle<>(v));
    }

    // Initialize accumulator.
    _pixelizor.TestInputs(map, pbore, pofs, signal, weight);
    accumulator.TestInputs(map, pbore, pofs, signal, weight);
    
    for (int i_det = 0; i_det < n_det; ++i_det) {
        pointer.InitPerDet(i_det);
        for (int i_time = 0; i_time < n_time; ++i_time) {
            double coords[4];
            double weights[4];
            int pixel_offset;
            pointer.GetCoords(i_det, i_time, (double*)coords);
            pixel_offset = _pixelizor.GetPixel(i_det, i_time, (double*)coords);
            accumulator.Reverse(i_det, i_time, pixel_offset, coords, weights);
        }
    }

    return signal;
}

template<typename P, typename Z, typename A>
bp::object ProjectionEngine<P,Z,A>::coords(
    bp::object pbore, bp::object pofs, bp::object coord)
{

    BufferWrapper coordbuf;
    if (PyObject_GetBuffer(coord.ptr(), &coordbuf.view,
                           PyBUF_RECORDS) == -1) {
        PyErr_Clear();
        throw buffer_exception("coord");
    }

    if (coordbuf.view.ndim != 3)
        throw shape_exception("coord", "must have shape (n_det,n_t,n_coord)");

    int n_det = coordbuf.view.shape[0];
    int n_time = coordbuf.view.shape[1];

    auto pointer = P();
    auto _none = bp::object();
    pointer.TestInputs(_none, pbore, pofs, _none, _none);

    auto coords_out = (char*)coordbuf.view.buf;

    for (int i_det = 0; i_det < n_det; ++i_det) {
        pointer.InitPerDet(i_det);
        for (int i_time = 0; i_time < n_time; ++i_time) {
            double coords[4];
            pointer.GetCoords(i_det, i_time, (double*)coords);
            for (int ic=0; ic<4; ic++) {
                *(double*)(coords_out
                  + coordbuf.view.strides[0] * i_det
                  + coordbuf.view.strides[1] * i_time
                  + coordbuf.view.strides[2] * ic) = coords[ic];
            }
        }
    }

    return coord;
}

template<typename P, typename Z, typename A>
bp::object ProjectionEngine<P,Z,A>::pixels(
    bp::object pbore, bp::object pofs, bp::object pixel)
{
    BufferWrapper pixelbuf;
    if (PyObject_GetBuffer(pixel.ptr(), &pixelbuf.view,
                           PyBUF_RECORDS) == -1) {
        PyErr_Clear();
        throw buffer_exception("pixel");
    }

    if (pixelbuf.view.ndim != 2)
        throw shape_exception("pixel", "must have shape (n_det,n_t)");

    int n_det = pixelbuf.view.shape[0];
    int n_time = pixelbuf.view.shape[1];

    auto pointer = P();
    auto _none = bp::object();

    pointer.TestInputs(_none, pbore, pofs, _none, _none);
    _pixelizor.TestInputs(_none, _none, _none, _none, _none);

    for (int i_det = 0; i_det < n_det; ++i_det) {
        pointer.InitPerDet(i_det);
        for (int i_time = 0; i_time < n_time; ++i_time) {
            double coords[4];
            pointer.GetCoords(i_det, i_time, (double*)coords);
            int pixel_offset = _pixelizor.GetPixel(i_det, i_time, (double*)coords);

            *(int*)((char*)pixelbuf.view.buf
                    + pixelbuf.view.strides[0] * i_det
                    + pixelbuf.view.strides[1] * i_time) = pixel_offset;
        }
    }

    return pixel;
}

typedef ProjectionEngine<PointerP_Simple_Flat,Pixelizor2_Flat,Accumulator<SpinT>>
  ProjectionEngine0;
typedef ProjectionEngine<PointerP_Simple_Flat,Pixelizor2_Flat,Accumulator<SpinQU>>
  ProjectionEngine1;
typedef ProjectionEngine<PointerP_Simple_Flat,Pixelizor2_Flat,Accumulator<SpinTQU>>
  ProjectionEngine2;

#define EXPORT_ENGINE(CLASSNAME)                                        \
    bp::class_<CLASSNAME>(#CLASSNAME, bp::init<Pixelizor2_Flat>())      \
    .def("to_map", &CLASSNAME::to_map)                                  \
    .def("to_weight_map", &CLASSNAME::to_weight_map)                    \
    .def("from_map", &CLASSNAME::from_map)                              \
    .def("coords", &CLASSNAME::coords)                                  \
    .def("pixels", &CLASSNAME::pixels);

PYBINDINGS("so3g")
{
    EXPORT_ENGINE(ProjectionEngine0);
    EXPORT_ENGINE(ProjectionEngine1);
    EXPORT_ENGINE(ProjectionEngine2);
    bp::class_<Pixelizor2_Flat>("Pixelizor2_Flat", bp::init<int,int,double,double,
                          double,double,double,double>())
        .def("zeros", &Pixelizor2_Flat::zeros);
}