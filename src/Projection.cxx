#define NO_IMPORT_ARRAY

// debug
#include <iostream>
using namespace std;

#include <pybindings.h>

#include <assert.h>
#include <math.h>

#include <omp.h>

#include <container_pybindings.h>

#include "so3g_numpy.h"
#include <Projection.h>
#include <Intervals.h>
#include "exceptions.h"

#include <boost/math/quaternion.hpp>

typedef boost::math::quaternion<double> quatd;

inline bool isNone(bp::object &pyo)
{
    return (pyo.ptr() == Py_None);
}

template <typename CoordSys>
bool Pointer<CoordSys>::TestInputs(
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

template <typename CoordSys>
inline
void Pointer<CoordSys>::InitPerDet(int i_det, double *dofs)
{
    const char *det = (char*)_pdetbuf.view.buf
        + _pdetbuf.view.strides[0] * i_det;
    for (int ic = 0; ic < 4; ++ic)
        dofs[ic] = *(double*)(det + _pdetbuf.view.strides[1] * ic);
}

/* ProjQuat: Not a projection -- returns the quaternion rotation
 * components.
 */

template <>
inline
void Pointer<ProjQuat>::GetCoords(int i_det, int i_time,
                                  const double *dofs, double *coords)
{
    double _qbore[4];
    for (int ic=0; ic<4; ic++)
        _qbore[ic] = *(double*)((char*)_pborebuf.view.buf +
                            _pborebuf.view.strides[0] * i_time +
                            _pborebuf.view.strides[1] * ic);

    // What could possibly go wrong.
    const quatd *qbore = reinterpret_cast<const quatd*>(_qbore);
    const quatd *qofs = reinterpret_cast<const quatd*>(dofs);
    quatd *qdet = reinterpret_cast<quatd*>(coords);
    *qdet =(*qbore) * (*qofs);
}

/* ProjFlat: Not a spherical projection -- assumes flat space (as in
 * FITS X,Y type coordinates). */

template <>
inline
void Pointer<ProjFlat>::GetCoords(int i_det, int i_time,
                                  const double *dofs, double *coords)
{
    for (int ic=0; ic<4; ic++)
        coords[ic] = *(double*)((char*)_pborebuf.view.buf +
                                _pborebuf.view.strides[0] * i_time +
                                _pborebuf.view.strides[1] * ic);
    coords[0] += dofs[0];
    coords[1] += dofs[1];
    const double coords_2_ = coords[2];
    coords[2] = coords[2] * dofs[2] - coords[3] * dofs[3];
    coords[3] = coords[3] * dofs[2] + coords_2_ * dofs[3];
}

/* ProjARC: the zenithal equidistant projection.
 *
 * The first two coordinates are R(lat)*sin(lon) and -R(lat)*cos(lon)
 * [see FITS-II].  Then cos and sin of parallactic angle.  For ARC,
 * R(lat) = 90 - lat = theta.
 */

template <>
inline
void Pointer<ProjARC>::GetCoords(int i_det, int i_time,
                                 const double *dofs, double *coords)
{
    double _qbore[4];
    for (int ic=0; ic<4; ic++)
        _qbore[ic] = *(double*)((char*)_pborebuf.view.buf +
                            _pborebuf.view.strides[0] * i_time +
                            _pborebuf.view.strides[1] * ic);

    // What could possibly go wrong.
    const quatd *qbore = reinterpret_cast<const quatd*>(_qbore);
    const quatd *qofs = reinterpret_cast<const quatd*>(dofs);
    quatd qdet = (*qbore) * (*qofs);

    const double a = qdet.R_component_1();
    const double b = qdet.R_component_2();
    const double c = qdet.R_component_3();
    const double d = qdet.R_component_4();

    const double cos_theta2_sq = a*a + d*d;

    const double sc = c*a + d*b;
    const double ss = a*b - c*d;
    const double half_sin_theta = sqrt(sc*sc + ss*ss);
    double R_factor;
    if (half_sin_theta < 1e-8)
        R_factor = 2 + 1.33333333333*half_sin_theta*half_sin_theta;
    else
        R_factor = asin(half_sin_theta*2) / half_sin_theta;

    coords[0] = ss * R_factor;
    coords[1] = sc * R_factor;
    coords[2] = (a*a - d*d) / cos_theta2_sq;
    coords[3] = (2*a*d) / cos_theta2_sq;
}

/* ProjTAN: the tangent plane (gnomonic) projection.  It is zenithal.
 *
 * The first two coordinates are R(lat)*sin(lon) and -R(lat)*cos(lon)
 * [see FITS-II].  Then cos and sin of parallactic angle.  For TAN,
 * R(lat) = tan(lat) = cot(theta).
 */

template <>
inline
void Pointer<ProjTAN>::GetCoords(int i_det, int i_time,
                                 const double *dofs, double *coords)
{
    double _qbore[4];
    for (int ic=0; ic<4; ic++)
        _qbore[ic] = *(double*)((char*)_pborebuf.view.buf +
                            _pborebuf.view.strides[0] * i_time +
                            _pborebuf.view.strides[1] * ic);

    // What could possibly go wrong.
    const quatd *qbore = reinterpret_cast<const quatd*>(_qbore);
    const quatd *qofs = reinterpret_cast<const quatd*>(dofs);
    quatd qdet = (*qbore) * (*qofs);

    const double a = qdet.R_component_1();
    const double b = qdet.R_component_2();
    const double c = qdet.R_component_3();
    const double d = qdet.R_component_4();

    const double cos_theta2_sq = a*a + d*d;
    const double cos_theta = 2*cos_theta2_sq - 1;

    coords[0] = (a*b - c*d) * 2 / cos_theta;
    coords[1] = (a*c + b*d) * 2 / cos_theta;
    coords[2] = (a*a - d*d) / cos_theta2_sq;
    coords[3] = (2*a*d) / cos_theta2_sq;
}

/* ProjZEA: the zenithal equal area projection.
 *
 * The first two coordinates are R(lat)*sin(lon) and -R(lat)*cos(lon)
 * [see FITS-II].  Then cos and sin of parallactic angle.  For ZEA,
 * R(lat) = sqrt(2(1-cos(lat))) = sqrt(2(1-sin(theta))).
 */

template <>
inline
void Pointer<ProjZEA>::GetCoords(int i_det, int i_time,
                                 const double *dofs, double *coords)
{
    double _qbore[4];
    for (int ic=0; ic<4; ic++)
        _qbore[ic] = *(double*)((char*)_pborebuf.view.buf +
                            _pborebuf.view.strides[0] * i_time +
                            _pborebuf.view.strides[1] * ic);

    // What could possibly go wrong.
    const quatd *qbore = reinterpret_cast<const quatd*>(_qbore);
    const quatd *qofs = reinterpret_cast<const quatd*>(dofs);
    quatd qdet = (*qbore) * (*qofs);

    const double a = qdet.R_component_1();
    const double b = qdet.R_component_2();
    const double c = qdet.R_component_3();
    const double d = qdet.R_component_4();

    const double cos_theta2_sq = a*a + d*d;
    const double cos_theta = 2*cos_theta2_sq - 1;
    const double cos_theta2 = sqrt(cos_theta2_sq);

    coords[0] = (a*b - c*d) * 2 / cos_theta2;
    coords[1] = (a*c + b*d) * 2 / cos_theta2;
    coords[2] = (a*a - d*d) / cos_theta2_sq;
    coords[3] = (2*a*d) / cos_theta2_sq;
}

/* ProjCEA: Cylindrical projection.
 *
 * First two coordinates are lon (in radians) and sin(lat).  Then cos
 * and sin of parallactic angle.
 */

template <>
inline
void Pointer<ProjCEA>::GetCoords(int i_det, int i_time,
                                 const double *dofs, double *coords)
{
    double _qbore[4];
    for (int ic=0; ic<4; ic++)
        _qbore[ic] = *(double*)((char*)_pborebuf.view.buf +
                            _pborebuf.view.strides[0] * i_time +
                            _pborebuf.view.strides[1] * ic);

    // What could possibly go wrong.
    const quatd *qbore = reinterpret_cast<const quatd*>(_qbore);
    const quatd *qofs = reinterpret_cast<const quatd*>(dofs);
    quatd qdet = (*qbore) * (*qofs);

    const double a = qdet.R_component_1();
    const double b = qdet.R_component_2();
    const double c = qdet.R_component_3();
    const double d = qdet.R_component_4();

    const double cos_theta = a*a - b*b - c*c + d*d;
    const double half_sin_theta = 0.5 * sqrt(1 - cos_theta*cos_theta);

    coords[0] = atan2(c*d - a*b, c*a + d*b);
    coords[1] = cos_theta; // Yes, cos(theta) = sin(lat).
    coords[2] = (a*c - b*d) / half_sin_theta;
    coords[3] = (c*d + a*b) / half_sin_theta;
}

/* ProjCAR: Cylindrical projection.
 *
 * First two coordinates are lon and lat (in radians).  Then cos and
 * sin of parallactic angle.
 */

template <>
inline
void Pointer<ProjCAR>::GetCoords(int i_det, int i_time,
                                      const double *dofs, double *coords)
{
    double _qbore[4];
    for (int ic=0; ic<4; ic++)
        _qbore[ic] = *(double*)((char*)_pborebuf.view.buf +
                            _pborebuf.view.strides[0] * i_time +
                            _pborebuf.view.strides[1] * ic);

    // What could possibly go wrong.
    const quatd *qbore = reinterpret_cast<const quatd*>(_qbore);
    const quatd *qofs = reinterpret_cast<const quatd*>(dofs);
    quatd qdet = (*qbore) * (*qofs);

    const double a = qdet.R_component_1();
    const double b = qdet.R_component_2();
    const double c = qdet.R_component_3();
    const double d = qdet.R_component_4();

    const double cos_theta = a*a - b*b - c*c + d*d;
    const double half_sin_theta = 0.5 * sqrt(1 - cos_theta*cos_theta);

    coords[0] = atan2(c*d - a*b, c*a + d*b);
    coords[1] = asin(cos_theta);   // Yes, cos(theta) = sin(lat).
    coords[2] = (a*c - b*d) / half_sin_theta;
    coords[3] = (c*d + a*b) / half_sin_theta;
}

Pixelizor2_Flat::Pixelizor2_Flat(
    int ny, int nx,
    double dy, double dx,
    double iy0, double ix0)
{
    naxis[0] = ny;
    naxis[1] = nx;
    cdelt[0] = dy;
    cdelt[1] = dx;
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
    double ix = coords[0] / cdelt[1] + crpix[1] + 0.5;
    if (ix < 0 || ix >= naxis[1])
        return -1;

    double iy = coords[1] / cdelt[0] + crpix[0] + 0.5;
    if (iy < 0 || iy >= naxis[0])
        return -1;

    int pixel_offset = strides[0]*int(iy) + strides[1]*int(ix);

    return pixel_offset;
}

std::pair<int,int> Pixelizor2_Flat::IndexRange()
{
    return make_pair(0,naxis[0]*naxis[1]);
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
        _signalspace = new SignalSpace<FSIGNAL>(
            signal, "signal", FSIGNAL_NPY_TYPE, n_det, n_time);
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
    const double* coords, FSIGNAL *pwt)
{
    pwt[0] = 1;
}

template <>
inline
void Accumulator<SpinQU>::PixelWeight(
    const double* coords, FSIGNAL *pwt)
{
    const double c = coords[2];
    const double s = coords[3];
    pwt[0] = c*c - s*s;
    pwt[1] = 2*c*s;
}

template<>
inline
void Accumulator<SpinTQU>::PixelWeight(
    const double* coords, FSIGNAL *pwt)
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
    const int pixel_offset, const double* coords, const FSIGNAL* weights)
{
    if (pixel_offset < 0) return;
    const FSIGNAL sig = *(_signalspace->data_ptr[i_det] +
                          _signalspace->steps[0]*i_time);
    const int N = SpinClass::comp_count;
    FSIGNAL wt[N];
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
    const int pixel_offset, const double* coords, const FSIGNAL* weights)
{
    if (pixel_offset < 0) return;
    const int N = SpinClass::comp_count;
    FSIGNAL wt[N];
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
    const int pixel_offset, const double* coords, const FSIGNAL* weights)
{
    if (pixel_offset < 0) return;
    const int N = SpinClass::comp_count;
    FSIGNAL wt[N];
    PixelWeight(coords, wt);
    double _sig = 0.;
    for (int imap=0; imap<N; ++imap) {
        _sig += *(double*)((char*)_mapbuf.view.buf +
                           _mapbuf.view.strides[0]*imap +
                           pixel_offset) * wt[imap];
    }
    FSIGNAL *sig = (_signalspace->data_ptr[i_det] +
                    _signalspace->steps[0]*i_time);
    *sig += _sig;
}


template <typename DTYPE>
bool SignalSpace<DTYPE>::_Validate(bp::object input, std::string var_name,
                                   int dtype, std::vector<int> dims)
{
    // The first axis is special; identify it with detector count.
    int n_det = dims[0];

    // We want a list of arrays here.
    bp::list sig_list;
    auto list_extractor = bp::extract<bp::list>(input);
    if (isNone(input)) {
        npy_intp _dims[dims.size()];
        for (int d=0; d<dims.size(); ++d)
            _dims[d] = dims[d];
        for (int i=0; i<n_det; ++i) {
            PyObject *v = PyArray_ZEROS(dims.size()-1, _dims+1, dtype, 0);
            sig_list.append(bp::object(bp::handle<>(v)));
        }
    } else if (list_extractor.check()) {
        sig_list = list_extractor();
    } else {
        // Probably an array... listify it.
        for (int i=0; i<bp::len(input); ++i)
            sig_list.append(input[i]);
    }
    ret_val = sig_list;

    if (bp::len(sig_list) != n_det)
        throw shape_exception(var_name, "must contain (n_det) vectors"); 

    // Now extract each list member into our vector of BufferWrappers..
    data_ptr = (DTYPE**)calloc(n_det, sizeof(*data_ptr));

    bw.reserve(n_det);
    for (int i=0; i<n_det; i++) {
        bw.push_back(BufferWrapper());
        BufferWrapper &_bw = bw[i];
        bp::object item = bp::extract<bp::object>(sig_list[i])();
        if (PyObject_GetBuffer(item.ptr(),
                               &_bw.view, PyBUF_RECORDS) == -1) {
            PyErr_Clear();
            throw buffer_exception(var_name);
        }
        if (_bw.view.ndim != dims.size() - 1)
            throw shape_exception(var_name, "must have the right number of dimsons");
        if (strcmp(bw[i].view.format, bw[0].view.format) != 0)
            throw dtype_exception(var_name, "[all elements must have same type]");
        for (int d=0; d<dims.size() - 1; ++d) {
            if (bw[i].view.shape[d] != dims[d+1])
                throw shape_exception(var_name, "must have right shape in all dimensions");
            if (bw[i].view.strides[d] != bw[0].view.strides[d])
                throw shape_exception(var_name, "[all elements must have same stride]");
        }
        data_ptr[i] = (DTYPE*)bw[i].view.buf;
    }
    // Check the dtype
    // FIXME; this does not check the dtype, it only checks the sizes.
    if (bw[0].view.itemsize != sizeof(DTYPE))
        throw dtype_exception(var_name, "[itemsize does not match expectation]");

    // Check the stride and store the step in units of the itemsize.
    //steps.empty();
    for (int d=0; d<dims.size()-1; d++) {
        if (bw[0].view.strides[d] % bw[0].view.itemsize != 0)
            throw shape_exception(var_name, "stride is non-integral; realign.");
        //steps.push_back(bw[0].view.strides[d] / bw[0].view.itemsize);
        steps[d] = bw[0].view.strides[d] / bw[0].view.itemsize;
    }
    return true;
}

template <typename DTYPE>
SignalSpace<DTYPE>::SignalSpace(
    bp::object input, std::string var_name, int dtype, int n_det, int n_time)
{
    vector<int> dims = {n_det, n_time};
    _Validate(input, var_name, dtype, dims);
}

template <typename DTYPE>
SignalSpace<DTYPE>::SignalSpace(
    bp::object input, std::string var_name, int dtype, int n_det, int n_time,
    int n_thirdaxis)
{
    vector<int> dims = {n_det, n_time, n_thirdaxis};
    _Validate(input, var_name, dtype, dims);
}


/** to_map(map, qpoint, pofs, signal, weights)
 *
 *  Each argument is an ndarray.  In the general case the dimensionalities are:
 *
 *     map:      (n_map, ny, nx, ...)
 *     pbore:    (n_t, n_coord)
 *     pofs:     (n_det, n_coord)
 *     signal:   (n_det, n_t)
 *     weight:   (n_det, n_map)
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
    //Initialize it / check inputs.
    auto pointer = P();
    pointer.TestInputs(map, pbore, pofs, signal, weight);
    int n_det = pointer.DetCount();
    int n_time = pointer.TimeCount();

    auto accumulator = A(true, true, false, n_det, n_time);

    //Do we need a map?  Now is the time.
    if (isNone(map)) {
        int n_comp = accumulator.ComponentCount();
        map = _pixelizor.zeros(n_comp);
    }

    _pixelizor.TestInputs(map, pbore, pofs, signal, weight);
    accumulator.TestInputs(map, pbore, pofs, signal, weight);

    for (int i_det = 0; i_det < n_det; ++i_det) {
        double dofs[4];
        pointer.InitPerDet(i_det, dofs);
        for (int i_time = 0; i_time < n_time; ++i_time) {
            double coords[4];
            FSIGNAL weights[4];
            int pixel_offset;
            pointer.GetCoords(i_det, i_time, (double*)dofs, (double*)coords);
            pixel_offset = _pixelizor.GetPixel(i_det, i_time, (double*)coords);
            accumulator.Forward(i_det, i_time, pixel_offset, coords, weights);
        }
    }

    return map;
}

template<typename P, typename Z, typename A>
bp::object ProjectionEngine<P,Z,A>::to_map_omp(
    bp::object map, bp::object pbore, bp::object pofs, bp::object signal, bp::object weight,
    bp::object thread_intervals)
{
    auto _none = bp::object();

    //Initialize it / check inputs.
    auto pointer = P();
    pointer.TestInputs(map, pbore, pofs, signal, weight);
    int n_det = pointer.DetCount();
    int n_time = pointer.TimeCount();

    auto accumulator = A(true, true, false, n_det, n_time);

    //Do we need a map?  Now is the time.
    if (isNone(map)) {
        int n_comp = accumulator.ComponentCount();
        map = _pixelizor.zeros(n_comp);
    }

    _pixelizor.TestInputs(map, pbore, pofs, signal, weight);
    accumulator.TestInputs(map, pbore, pofs, signal, weight);

    // Indexed by i_thread, i_det.
    vector<vector<IntervalsInt32>> ivals;

    auto ival_list_list = bp::extract<bp::list>(thread_intervals)();
    for (int i=0; i<bp::len(ival_list_list); i++) {
        auto ival_list =  bp::extract<bp::list>(ival_list_list[i])();
        vector<IntervalsInt32> v(bp::len(ival_list));
        for (int j=0; j<bp::len(ival_list); j++)
            v[j] = bp::extract<IntervalsInt32>(ival_list[j])();
        ivals.push_back(v);
    }

#pragma omp parallel
    {
        // The principle here is that all threads loop over all
        // detectors, but the sample ranges encoded in ivals are
        // disjoint.
        const int i_thread = omp_get_thread_num();
        for (int i_det = 0; i_det < n_det; ++i_det) {
            double dofs[4];
            pointer.InitPerDet(i_det, dofs);
            for (auto rng: ivals[i_thread][i_det].segments) {
                for (int i_time = rng.first; i_time < rng.second; ++i_time) {
                    double coords[4];
                    FSIGNAL weights[4];
                    int pixel_offset;
                    pointer.GetCoords(i_det, i_time, (double*)dofs, (double*)coords);
                    pixel_offset = _pixelizor.GetPixel(i_det, i_time, (double*)coords);
                    accumulator.Forward(i_det, i_time, pixel_offset, coords, weights);
                }
            }
        }
    }
    return map;
}

template<typename P, typename Z, typename A>
bp::object ProjectionEngine<P,Z,A>::to_weight_map(
    bp::object map, bp::object pbore, bp::object pofs, bp::object signal, bp::object weight)
{
    //Initialize it / check inputs.
    auto pointer = P();
    pointer.TestInputs(map, pbore, pofs, signal, weight);
    int n_det = pointer.DetCount();
    int n_time = pointer.TimeCount();

    auto accumulator = A(false, false, true, n_det, n_time);

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
        // pointer.InitPerDet(i_det);
        double dofs[4];
        pointer.InitPerDet(i_det, dofs);
        for (int i_time = 0; i_time < n_time; ++i_time) {
            double coords[4];
            FSIGNAL weights[4];
            int pixel_offset;
            pointer.GetCoords(i_det, i_time, (double*)dofs, (double*)coords);
            pixel_offset = _pixelizor.GetPixel(i_det, i_time, (double*)coords);
            accumulator.ForwardWeight(i_det, i_time, pixel_offset, coords, weights);
        }
    }

    return map;
}

template<typename P, typename Z, typename A>
bp::object ProjectionEngine<P,Z,A>::to_weight_map_omp(
    bp::object map, bp::object pbore, bp::object pofs, bp::object signal, bp::object weight,
    bp::object thread_intervals)
{
    auto _none = bp::object();

    //Initialize it / check inputs.
    auto pointer = P();
    pointer.TestInputs(map, pbore, pofs, signal, weight);
    int n_det = pointer.DetCount();
    int n_time = pointer.TimeCount();

    auto accumulator = A(false, false, true, n_det, n_time);

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

    // Indexed by i_thread, i_det.
    vector<vector<IntervalsInt32>> ivals;

    auto ival_list_list = bp::extract<bp::list>(thread_intervals)();
    for (int i=0; i<bp::len(ival_list_list); i++) {
        auto ival_list =  bp::extract<bp::list>(ival_list_list[i])();
        vector<IntervalsInt32> v(bp::len(ival_list));
        for (int j=0; j<bp::len(ival_list); j++)
            v[j] = bp::extract<IntervalsInt32>(ival_list[j])();
        ivals.push_back(v);
    }

#pragma omp parallel
    {
        // The principle here is that all threads loop over all
        // detectors, but the sample ranges encoded in ivals are
        // disjoint.
        const int i_thread = omp_get_thread_num();
        for (int i_det = 0; i_det < n_det; ++i_det) {
            double dofs[4];
            pointer.InitPerDet(i_det, dofs);
            for (auto rng: ivals[i_thread][i_det].segments) {
                for (int i_time = rng.first; i_time < rng.second; ++i_time) {
                    double coords[4];
                    FSIGNAL weights[4];
                    int pixel_offset;
                    pointer.GetCoords(i_det, i_time, (double*)dofs, (double*)coords);
                    pixel_offset = _pixelizor.GetPixel(i_det, i_time, (double*)coords);
                    accumulator.ForwardWeight(i_det, i_time, pixel_offset, coords, weights);
                }
            }
        }
    }

    return map;
}

template<typename P, typename Z, typename A>
bp::object ProjectionEngine<P,Z,A>::from_map(
    bp::object map, bp::object pbore, bp::object pofs, bp::object signal, bp::object weight)
{
    // Initialize pointer and _pixelizor.
    auto pointer = P();
    pointer.TestInputs(map, pbore, pofs, signal, weight);
    int n_det = pointer.DetCount();
    int n_time = pointer.TimeCount();

    // Initialize accumulator -- create signal if it DNE.
    auto accumulator = A(true, true, false, n_det, n_time);
    accumulator.TestInputs(map, pbore, pofs, signal, weight);

    _pixelizor.TestInputs(map, pbore, pofs, signal, weight);

#pragma omp parallel for
    for (int i_det = 0; i_det < n_det; ++i_det) {
        double dofs[4];
        pointer.InitPerDet(i_det, dofs);
        for (int i_time = 0; i_time < n_time; ++i_time) {
            double coords[4];
            FSIGNAL weights[4];
            int pixel_offset;
            pointer.GetCoords(i_det, i_time, (double*)dofs, (double*)coords);
            pixel_offset = _pixelizor.GetPixel(i_det, i_time, (double*)coords);
            accumulator.Reverse(i_det, i_time, pixel_offset, coords, weights);
        }
    }

    return accumulator._signalspace->ret_val;
}

template<typename P, typename Z, typename A>
bp::object ProjectionEngine<P,Z,A>::coords(
    bp::object pbore, bp::object pofs, bp::object coord)
{
    auto _none = bp::object();
    auto pointer = P();
    pointer.TestInputs(_none, pbore, pofs, _none, _none);

    int n_det = pointer.DetCount();
    int n_time = pointer.TimeCount();

    const int n_coord = 4;
    auto coord_buf_man = SignalSpace<double>(
        coord, "coord", NPY_FLOAT64, n_det, n_time, n_coord);

#pragma omp parallel for
    for (int i_det = 0; i_det < n_det; ++i_det) {
        double dofs[4];
        pointer.InitPerDet(i_det, dofs);

        double* const coords_det = coord_buf_man.data_ptr[i_det];
        const int step0 = coord_buf_man.steps[0];
        const int step1 = coord_buf_man.steps[1];

        for (int i_time = 0; i_time < n_time; ++i_time) {
            double coords[4];
            pointer.GetCoords(i_det, i_time, (double*)dofs, (double*)coords);
            for (int ic=0; ic<4; ic++)
                *(coords_det + step0 * i_time + step1 * ic) = coords[ic];
        }
    }

    return coord_buf_man.ret_val;
}

template<typename P, typename Z, typename A>
bp::object ProjectionEngine<P,Z,A>::pixels(
    bp::object pbore, bp::object pofs, bp::object pixel)
{
    auto _none = bp::object();

    auto pointer = P();
    pointer.TestInputs(_none, pbore, pofs, _none, _none);
    _pixelizor.TestInputs(_none, _none, _none, _none, _none);

    int n_det = pointer.DetCount();
    int n_time = pointer.TimeCount();

    auto pixel_buf_man = SignalSpace<int32_t>(
        pixel, "pixel", NPY_INT32, n_det, n_time);

#pragma omp parallel for
    for (int i_det = 0; i_det < n_det; ++i_det) {
        double dofs[4];
        pointer.InitPerDet(i_det, dofs);
        int* const pix_buf = pixel_buf_man.data_ptr[i_det];
        const int step = pixel_buf_man.steps[0];
        for (int i_time = 0; i_time < n_time; ++i_time) {
            double coords[4];
            pointer.GetCoords(i_det, i_time, (double*)dofs, (double*)coords);
            int pixel_offset = _pixelizor.GetPixel(i_det, i_time, (double*)coords);
            pix_buf[i_time * step] = pixel_offset;
            // pix_buf[i_time * pixel_buf_man.steps[0]] = pixel_offset;
        }
    }

    return pixel_buf_man.ret_val;
}


template<typename P, typename Z, typename A>
bp::object ProjectionEngine<P,Z,A>::pixel_ranges(
    bp::object pbore, bp::object pofs)
{
    auto pointer = P();
    auto _none = bp::object();

    pointer.TestInputs(_none, pbore, pofs, _none, _none);
    _pixelizor.TestInputs(_none, _none, _none, _none, _none);

    int n_det = pointer.DetCount();
    int n_time = pointer.TimeCount();

    auto pix_range = _pixelizor.IndexRange();

    vector<vector<IntervalsInt32>> ranges;

#pragma omp parallel
    {
        int n_domain = omp_get_num_threads();
#pragma omp single
        {
            for (int i=0; i<n_domain; ++i) {
                vector<IntervalsInt32> v(n_det);
                ranges.push_back(v);
            }
        }

        int pix_lo = pix_range.first;
        int pix_step = (pix_range.second - pix_range.first +
                        n_domain - 1) / n_domain;

#pragma omp for
        for (int i_det = 0; i_det < n_det; ++i_det) {
            double dofs[4];
            pointer.InitPerDet(i_det, dofs);
            int last_slice = -1;
            int slice_start = 0;
            for (int i_time = 0; i_time < n_time; ++i_time) {
                double coords[4];
                pointer.GetCoords(i_det, i_time, (double*)dofs, (double*)coords);
                int pixel_offset = _pixelizor.GetPixel(i_det, i_time, (double*)coords);
                int this_slice = -1;
                if (pixel_offset >= 0)
                    this_slice = (pixel_offset - pix_lo) / pix_step;
                if (this_slice != last_slice) {
                    if (last_slice >= 0)
                        ranges[last_slice][i_det].append_interval_no_check(
                            slice_start, i_time);
                    slice_start = i_time;
                    last_slice = this_slice;
                }
            }
            if (last_slice >= 0)
                ranges[last_slice][i_det].append_interval_no_check(
                    slice_start, n_time);
        }
    }

    // Convert super vector to a list and return
    auto ivals_out = bp::list();
    for (int j=0; j<ranges.size(); j++) {
        auto ivals = bp::list();
        for (int i_det=0; i_det<n_det; i_det++) {
            auto iv = ranges[j][i_det];
            ivals.append(bp::object(iv));
        }
        ivals_out.append(bp::extract<bp::object>(ivals)());
    }
    return bp::extract<bp::object>(ivals_out);
}

//Flat.
typedef ProjectionEngine<Pointer<ProjFlat>,Pixelizor2_Flat,Accumulator<SpinT>>
  ProjEng_Flat_T;
typedef ProjectionEngine<Pointer<ProjFlat>,Pixelizor2_Flat,Accumulator<SpinQU>>
  ProjEng_Flat_QU;
typedef ProjectionEngine<Pointer<ProjFlat>,Pixelizor2_Flat,Accumulator<SpinTQU>>
  ProjEng_Flat_TQU;
typedef ProjectionEngine<Pointer<ProjCEA>,Pixelizor2_Flat,Accumulator<SpinT>>

//Cylindrical.
  ProjEng_CEA_T;
typedef ProjectionEngine<Pointer<ProjCEA>,Pixelizor2_Flat,Accumulator<SpinQU>>
  ProjEng_CEA_QU;
typedef ProjectionEngine<Pointer<ProjCEA>,Pixelizor2_Flat,Accumulator<SpinTQU>>
  ProjEng_CEA_TQU;
typedef ProjectionEngine<Pointer<ProjCAR>,Pixelizor2_Flat,Accumulator<SpinT>>
  ProjEng_CAR_T;
typedef ProjectionEngine<Pointer<ProjCAR>,Pixelizor2_Flat,Accumulator<SpinQU>>
  ProjEng_CAR_QU;
typedef ProjectionEngine<Pointer<ProjCAR>,Pixelizor2_Flat,Accumulator<SpinTQU>>
  ProjEng_CAR_TQU;

//Zenithal.
typedef ProjectionEngine<Pointer<ProjARC>,Pixelizor2_Flat,Accumulator<SpinT>>
  ProjEng_ARC_T;
typedef ProjectionEngine<Pointer<ProjARC>,Pixelizor2_Flat,Accumulator<SpinQU>>
  ProjEng_ARC_QU;
typedef ProjectionEngine<Pointer<ProjARC>,Pixelizor2_Flat,Accumulator<SpinTQU>>
  ProjEng_ARC_TQU;
typedef ProjectionEngine<Pointer<ProjTAN>,Pixelizor2_Flat,Accumulator<SpinT>>
  ProjEng_TAN_T;
typedef ProjectionEngine<Pointer<ProjTAN>,Pixelizor2_Flat,Accumulator<SpinQU>>
  ProjEng_TAN_QU;
typedef ProjectionEngine<Pointer<ProjTAN>,Pixelizor2_Flat,Accumulator<SpinTQU>>
  ProjEng_TAN_TQU;
typedef ProjectionEngine<Pointer<ProjZEA>,Pixelizor2_Flat,Accumulator<SpinT>>
  ProjEng_ZEA_T;
typedef ProjectionEngine<Pointer<ProjZEA>,Pixelizor2_Flat,Accumulator<SpinQU>>
  ProjEng_ZEA_QU;
typedef ProjectionEngine<Pointer<ProjZEA>,Pixelizor2_Flat,Accumulator<SpinTQU>>
  ProjEng_ZEA_TQU;

#define EXPORT_ENGINE(CLASSNAME)                                        \
    bp::class_<CLASSNAME>(#CLASSNAME, bp::init<Pixelizor2_Flat>())      \
    .def("to_map", &CLASSNAME::to_map)                                  \
    .def("to_map_omp", &CLASSNAME::to_map_omp)                          \
    .def("to_weight_map", &CLASSNAME::to_weight_map)                    \
    .def("to_weight_map_omp", &CLASSNAME::to_weight_map_omp)            \
    .def("from_map", &CLASSNAME::from_map)                              \
    .def("coords", &CLASSNAME::coords)                                  \
    .def("pixels", &CLASSNAME::pixels)                                  \
    .def("pixel_ranges", &CLASSNAME::pixel_ranges);

PYBINDINGS("so3g")
{
    EXPORT_ENGINE(ProjEng_Flat_T);
    EXPORT_ENGINE(ProjEng_Flat_QU);
    EXPORT_ENGINE(ProjEng_Flat_TQU);
    EXPORT_ENGINE(ProjEng_CAR_T);
    EXPORT_ENGINE(ProjEng_CAR_QU);
    EXPORT_ENGINE(ProjEng_CAR_TQU);
    EXPORT_ENGINE(ProjEng_CEA_T);
    EXPORT_ENGINE(ProjEng_CEA_QU);
    EXPORT_ENGINE(ProjEng_CEA_TQU);
    EXPORT_ENGINE(ProjEng_ARC_T);
    EXPORT_ENGINE(ProjEng_ARC_QU);
    EXPORT_ENGINE(ProjEng_ARC_TQU);
    EXPORT_ENGINE(ProjEng_TAN_T);
    EXPORT_ENGINE(ProjEng_TAN_QU);
    EXPORT_ENGINE(ProjEng_TAN_TQU);
    EXPORT_ENGINE(ProjEng_ZEA_T);
    EXPORT_ENGINE(ProjEng_ZEA_QU);
    EXPORT_ENGINE(ProjEng_ZEA_TQU);
    bp::class_<Pixelizor2_Flat>("Pixelizor2_Flat", bp::init<int,int,double,double,
                          double,double>())
        .def("zeros", &Pixelizor2_Flat::zeros);
}
