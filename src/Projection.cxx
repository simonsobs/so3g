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
#include <Ranges.h>
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
    bp::object &signal, bp::object &det_weights)
{
    // Boresight and Detector must present and inter-compatible.
    _pborebuf = BufferWrapper<double>("boresight", pbore, false,
                                      vector<int>{-1, 4});
    _pdetbuf = BufferWrapper<double>("detectors", pdet, false,
                                     vector<int>{-1, 4});
    n_time = _pborebuf->shape[0];
    n_det = _pdetbuf->shape[0];
    return true;
}

template <typename CoordSys>
inline
void Pointer<CoordSys>::InitPerDet(int i_det, double *dofs)
{
    const char *det = (char*)_pdetbuf->buf
        + _pdetbuf->strides[0] * i_det;
    for (int ic = 0; ic < 4; ++ic)
        dofs[ic] = *(double*)(det + _pdetbuf->strides[1] * ic);
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
        _qbore[ic] = *(double*)((char*)_pborebuf->buf +
                            _pborebuf->strides[0] * i_time +
                            _pborebuf->strides[1] * ic);

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
        coords[ic] = *(double*)((char*)_pborebuf->buf +
                                _pborebuf->strides[0] * i_time +
                                _pborebuf->strides[1] * ic);
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
        _qbore[ic] = *(double*)((char*)_pborebuf->buf +
                            _pborebuf->strides[0] * i_time +
                            _pborebuf->strides[1] * ic);

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
        _qbore[ic] = *(double*)((char*)_pborebuf->buf +
                            _pborebuf->strides[0] * i_time +
                            _pborebuf->strides[1] * ic);

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
        _qbore[ic] = *(double*)((char*)_pborebuf->buf +
                            _pborebuf->strides[0] * i_time +
                            _pborebuf->strides[1] * ic);

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
        _qbore[ic] = *(double*)((char*)_pborebuf->buf +
                            _pborebuf->strides[0] * i_time +
                            _pborebuf->strides[1] * ic);

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
        _qbore[ic] = *(double*)((char*)_pborebuf->buf +
                            _pborebuf->strides[0] * i_time +
                            _pborebuf->strides[1] * ic);

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
                             bp::object &signal, bp::object &det_weights)
{
    // Flat pixelizor.  Requires:
    // 1. Map has the right shape -- but map can be none.
    // 2. Coords have the right format (unchecked currently).

    // If map is provided, it can have any number of leading
    // dimensions but the last two dimensions must match naxis.
    BufferWrapper<double> mapbuf("map", map, true,
                                 vector<int>{-2,naxis[0],naxis[1]});
    if (mapbuf->buf != NULL) {
        // Note these are byte offsets, not index.
        int ndim = mapbuf->ndim;
        strides[0] = mapbuf->strides[ndim-2];
        strides[1] = mapbuf->strides[ndim-1];
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

class Tiled {
public:
    vector<BufferWrapper<double>> mapbufs;
    vector<int> shape;
    vector<int> tile_shape;
    
    bool TestInputs(bp::object &map, bool need_map, bool need_weight_map, int comp_count) {
        std::cout << "Tiled::TestInputs\n";
        vector<int> map_shape_req;
        if (need_map) {
            // The map is mandatory, and the leading axis must match the
            // component count.  It can have 1+ other dimensions.
            map_shape_req = {comp_count,-1,-3};
        } else if (need_weight_map) {
            // The map is mandatory, and the two leading axes must match
            // the component count.  It can have 1+ other dimensions.
            map_shape_req = {comp_count,comp_count,-1,-3};
        }
        
        if (map_shape_req.size())
            mapbufs.push_back(
                BufferWrapper<double>("map", map, false, map_shape_req));

        return true;
    }        
    double *pix(int imap, int pixel_offset) {
        return (double*)((char*)mapbufs[0]->buf +
                         mapbufs[0]->strides[0]*imap +
                         pixel_offset);
    }
    double *wpix(int imap, int jmap, int pixel_offset) {
        return (double*)((char*)mapbufs[0]->buf +
                         mapbufs[0]->strides[0]*imap +
                         mapbufs[0]->strides[1]*jmap +
                         pixel_offset);
    }
};

class NonTiled {
public:
    BufferWrapper<double> mapbuf;
    bool TestInputs(bp::object &map, bool need_map, bool need_weight_map, int comp_count) {
        std::cout << "NonTiled::TestInputs\n";
        if (need_map) {
            // The map is mandatory, and the leading axis must match the
            // component count.  It can have 1+ other dimensions.
            mapbuf = BufferWrapper<double>("map", map, false,
                                           vector<int>{comp_count,-1,-3});
        } else if (need_weight_map) {
            // The map is mandatory, and the two leading axes must match
            // the component count.  It can have 1+ other dimensions.
            mapbuf = BufferWrapper<double>("map", map, false,
                                           vector<int>{comp_count,comp_count,-1,-3});
        }
        return true;
    }        
    double *pix(int imap, int pixel_offset) {
        return (double*)((char*)mapbuf->buf +
                         mapbuf->strides[0]*imap +
                         pixel_offset);
    }
    double *wpix(int imap, int jmap, int pixel_offset) {
        return (double*)((char*)mapbuf->buf +
                         mapbuf->strides[0]*imap +
                         mapbuf->strides[1]*jmap +
                         pixel_offset);
    }
};

template <typename SpinClass, typename TilingSystem>
bool Accumulator<SpinClass,TilingSystem>::TestInputs(
    bp::object &map, bp::object &pbore, bp::object &pdet,
    bp::object &signal, bp::object &det_weights)
{
    tiling.TestInputs(map, need_map, need_weight_map, SpinClass::comp_count);

    if (need_signal) {
        _signalspace = new SignalSpace<FSIGNAL>(
            signal, "signal", FSIGNAL_NPY_TYPE, n_det, n_time);
    }

     _det_weights = BufferWrapper<FSIGNAL>(
         "det_weights", det_weights, true, vector<int>{n_det});

    return true;
}

// template <typename TilingSystem>
// inline
// void Accumulator<SpinT,TilingSystem>::SpinProjFactors(
//     const double* coords, FSIGNAL *projfacs)
// {
//     projfacs[0] = 1;
// }

// template <typename TilingSystem>
// inline
// void Accumulator<SpinQU,TilingSystem>::SpinProjFactors(
//     const double* coords, FSIGNAL *projfacs)
// {
//     const double c = coords[2];
//     const double s = coords[3];
//     projfacs[0] = c*c - s*s;
//     projfacs[1] = 2*c*s;
// }

// template <typename TilingSystem>
// inline
// void Accumulator<SpinTQU,TilingSystem>::SpinProjFactors(
//     const double* coords, FSIGNAL *projfacs)
// {
//     const double c = coords[2];
//     const double s = coords[3];
//     projfacs[0] = 1.;
//     projfacs[1] = c*c - s*s;
//     projfacs[2] = 2*c*s;
// }

template <typename SpinClass,typename TilingSystem>
inline
void Accumulator<SpinClass,TilingSystem>::Forward(
    const int i_det, const int i_time,
    const int pixel_offset, const double* coords, const FSIGNAL* weights)
{
    if (pixel_offset < 0) return;
    const FSIGNAL sig = *(_signalspace->data_ptr[i_det] +
                          _signalspace->steps[0]*i_time);

    FSIGNAL det_wt = 1.;
    if (_det_weights->obj != NULL)
        det_wt = *(FSIGNAL*)((char*)_det_weights->buf +
                             _det_weights->strides[0]*i_det);

    const int N = SpinClass::comp_count;
    FSIGNAL pf[N];
    SpinClass::ProjFactors(coords, pf);
    for (int imap=0; imap<N; ++imap) {
        // *(double*)((char*)tiling.mapbuf->buf +
        //            tiling.mapbuf->strides[0]*imap +
        //            pixel_offset) += sig * pf[imap] * det_wt;
        *tiling.pix(imap, pixel_offset) += sig * pf[imap] * det_wt;
    }
}

template <typename SpinClass,typename TilingSystem>
inline
void Accumulator<SpinClass,TilingSystem>::ForwardWeight(
    const int i_det, const int i_time,
    const int pixel_offset, const double* coords, const FSIGNAL* weights)
{
    if (pixel_offset < 0) return;

    FSIGNAL det_wt = 1.;
    if (_det_weights->obj != NULL)
        det_wt = *(FSIGNAL*)((char*)_det_weights->buf + _det_weights->strides[0]*i_det);

    const int N = SpinClass::comp_count;
    FSIGNAL pf[N];
    SpinClass::ProjFactors(coords, pf);
    for (int imap=0; imap<N; ++imap) {
        for (int jmap=imap; jmap<N; ++jmap) {
            // *(double*)((char*)tiling.mapbuf->buf +
            //            tiling.mapbuf->strides[0]*imap +
            //            tiling.mapbuf->strides[1]*jmap +
            //            pixel_offset) += pf[imap] * pf[jmap] * det_wt;
            *tiling.wpix(imap, jmap, pixel_offset) += pf[imap] * pf[jmap] * det_wt;
        }
    }
}

template <typename SpinClass,typename TilingSystem>
inline
void Accumulator<SpinClass,TilingSystem>::Reverse(
    const int i_det, const int i_time,
    const int pixel_offset, const double* coords, const FSIGNAL* weights)
{
    if (pixel_offset < 0) return;
    const int N = SpinClass::comp_count;
    FSIGNAL pf[N];
    SpinClass::ProjFactors(coords, pf);
    double _sig = 0.;
    for (int imap=0; imap<N; ++imap) {
        // _sig += *(double*)((char*)tiling.mapbuf->buf +
        //                    tiling.mapbuf->strides[0]*imap +
        //                    pixel_offset) * pf[imap];
        _sig += *tiling.pix(imap, pixel_offset) * pf[imap];
    }
    FSIGNAL *sig = (_signalspace->data_ptr[i_det] +
                    _signalspace->steps[0]*i_time);
    *sig += _sig;
}


template <typename DTYPE>
bool SignalSpace<DTYPE>::_Validate(bp::object input, std::string var_name,
                                   int dtype)
{
    // We want a list of arrays here.
    bp::list sig_list;
    auto list_extractor = bp::extract<bp::list>(input);
    if (isNone(input)) {
        npy_intp _dims[dims.size()];
        for (int d=0; d<dims.size(); ++d) {
            if (dims[d] < 0)
                throw shape_exception(var_name, "Cannot create space with wildcard dimensons.");
            _dims[d] = dims[d];
        }
        for (int i=0; i<dims[0]; ++i) {
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

    if (dims[0] == -1) {
        dims[0] = bp::len(sig_list);
        if (dims[0] == 0)
            throw shape_exception(var_name, "has not been tested on shape 0 objects");
    } else if (bp::len(sig_list) != dims[0])
        throw shape_exception(var_name, "must contain (n_det) vectors"); 

    const int n_det = dims[0];

    // Now extract each list member into our vector of BufferWrappers..
    data_ptr = (DTYPE**)calloc(n_det, sizeof(*data_ptr));

    bw.reserve(n_det);

    // Copy dims[1:] into sub_dims; potentially update sub_dims during
    // extraction, then copy those back into dims.
    vector<int> sub_dims(dims.begin()+1, dims.end());

    for (int i=0; i<n_det; i++) {
        bp::object item = bp::extract<bp::object>(sig_list[i])();
        bw.push_back(BufferWrapper<DTYPE>(var_name, item, false, sub_dims));
        if (i == 0) {
            sub_dims.clear();
            for (int d=0; d<bw[0]->ndim; d++)
                sub_dims.push_back(bw[0]->shape[d]);
        } else {
            for (int d=0; d<sub_dims.size(); ++d) {
                if (bw[i]->strides[d] != bw[0]->strides[d])
                    throw shape_exception(var_name, "[all elements must have same stride]");
            }
        }
        data_ptr[i] = (DTYPE*)bw[i]->buf;
    }

    // Store the step in units of the itemsize; update dims from sub_dims.
    for (int d=0; d<dims.size()-1; d++) {
        dims[d+1] = sub_dims[d];
        if (bw[0]->strides[d] % bw[0]->itemsize != 0)
            throw shape_exception(var_name, "stride is non-integral; realign.");
        steps[d] = bw[0]->strides[d] / bw[0]->itemsize;
    }
    return true;
}

template <typename DTYPE>
SignalSpace<DTYPE>::SignalSpace(
    bp::object input, std::string var_name, int dtype, int n_det, int n_time)
{
    dims = {n_det, n_time};
    _Validate(input, var_name, dtype);
}

template <typename DTYPE>
SignalSpace<DTYPE>::SignalSpace(
    bp::object input, std::string var_name, int dtype, int n_det, int n_time,
    int n_thirdaxis)
{
    dims = {n_det, n_time, n_thirdaxis};
    _Validate(input, var_name, dtype);
}


/** to_map(map, qpoint, pofs, signal, weights)
 *
 *  Each argument is an ndarray.  In the general case the dimensionalities are:
 *
 *     map:      (n_map, ny, nx, ...)
 *     pbore:    (n_t, n_coord)
 *     pofs:     (n_det, n_coord)
 *     signal:   (n_det, n_t)
 *     det_weights: (n_det)
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
    bp::object map, bp::object pbore, bp::object pofs, bp::object signal, bp::object det_weights)
{
    //Initialize it / check inputs.
    auto pointer = P();
    pointer.TestInputs(map, pbore, pofs, signal, det_weights);
    int n_det = pointer.DetCount();
    int n_time = pointer.TimeCount();

    auto accumulator = A(true, true, false, n_det, n_time);

    //Do we need a map?  Now is the time.
    if (isNone(map)) {
        int n_comp = accumulator.ComponentCount();
        map = _pixelizor.zeros(n_comp);
    }

    _pixelizor.TestInputs(map, pbore, pofs, signal, det_weights);
    accumulator.TestInputs(map, pbore, pofs, signal, det_weights);

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
    bp::object map, bp::object pbore, bp::object pofs, bp::object signal, bp::object det_weights,
    bp::object thread_intervals)
{
    auto _none = bp::object();

    //Initialize it / check inputs.
    auto pointer = P();
    pointer.TestInputs(map, pbore, pofs, signal, det_weights);
    int n_det = pointer.DetCount();
    int n_time = pointer.TimeCount();

    auto accumulator = A(true, true, false, n_det, n_time);

    //Do we need a map?  Now is the time.
    if (isNone(map)) {
        int n_comp = accumulator.ComponentCount();
        map = _pixelizor.zeros(n_comp);
    }

    _pixelizor.TestInputs(map, pbore, pofs, signal, det_weights);
    accumulator.TestInputs(map, pbore, pofs, signal, det_weights);

    // Indexed by i_thread, i_det.
    vector<vector<RangesInt32>> ivals;

    // Descend two levels.. don't assume it's a list, just that it has
    // len and [].
    for (int i=0; i<bp::len(thread_intervals); i++) {
        bp::object ival_list = thread_intervals[i];
        vector<RangesInt32> v(bp::len(ival_list));
        for (int j=0; j<bp::len(ival_list); j++)
            v[j] = bp::extract<RangesInt32>(ival_list[j])();
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
            for (auto const &rng: ivals[i_thread][i_det].segments) {
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
    bp::object map, bp::object pbore, bp::object pofs, bp::object signal, bp::object det_weights)
{
    //Initialize it / check inputs.
    auto pointer = P();
    pointer.TestInputs(map, pbore, pofs, signal, det_weights);
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

    _pixelizor.TestInputs(map, pbore, pofs, signal, det_weights);
    accumulator.TestInputs(map, pbore, pofs, signal, det_weights);

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
    bp::object map, bp::object pbore, bp::object pofs, bp::object signal, bp::object det_weights,
    bp::object thread_intervals)
{
    auto _none = bp::object();

    //Initialize it / check inputs.
    auto pointer = P();
    pointer.TestInputs(map, pbore, pofs, signal, det_weights);
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

    _pixelizor.TestInputs(map, pbore, pofs, signal, det_weights);
    accumulator.TestInputs(map, pbore, pofs, signal, det_weights);

    // Indexed by i_thread, i_det.
    vector<vector<RangesInt32>> ivals;

    // Descend two levels.. don't assume it's a list, just that it has
    // len and [].
    for (int i=0; i<bp::len(thread_intervals); i++) {
        bp::object ival_list = thread_intervals[i];
        vector<RangesInt32> v(bp::len(ival_list));
        for (int j=0; j<bp::len(ival_list); j++)
            v[j] = bp::extract<RangesInt32>(ival_list[j])();
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
            for (auto const &rng: ivals[i_thread][i_det].segments) {
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
    bp::object map, bp::object pbore, bp::object pofs, bp::object signal, bp::object det_weights)
{
    // Initialize pointer and _pixelizor.
    auto pointer = P();
    pointer.TestInputs(map, pbore, pofs, signal, det_weights);
    int n_det = pointer.DetCount();
    int n_time = pointer.TimeCount();

    // Initialize accumulator -- create signal if it DNE.
    auto accumulator = A(true, true, false, n_det, n_time);
    accumulator.TestInputs(map, pbore, pofs, signal, det_weights);

    _pixelizor.TestInputs(map, pbore, pofs, signal, det_weights);

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
bp::object ProjectionEngine<P,Z,A>::pointing_matrix(
    bp::object pbore, bp::object pofs, bp::object pixel, bp::object proj)
{
    auto _none = bp::object();

    auto pointer = P();
    pointer.TestInputs(_none, pbore, pofs, _none, _none);
    _pixelizor.TestInputs(_none, _none, _none, _none, _none);

    int n_det = pointer.DetCount();
    int n_time = pointer.TimeCount();

    auto pixel_buf_man = SignalSpace<int32_t>(
        pixel, "pixel", NPY_INT32, n_det, n_time);

    auto accumulator = A(false, false, false, n_det, n_time);
    accumulator.TestInputs(_none, pbore, pofs, _none, _none);

    const int n_spin = accumulator.ComponentCount();
    auto proj_buf_man = SignalSpace<FSIGNAL>(
        proj, "proj", FSIGNAL_NPY_TYPE, n_det, n_time, n_spin);

#pragma omp parallel for
    for (int i_det = 0; i_det < n_det; ++i_det) {
        double dofs[4];
        pointer.InitPerDet(i_det, dofs);
        int* const pix_buf = pixel_buf_man.data_ptr[i_det];
        FSIGNAL* const proj_buf = proj_buf_man.data_ptr[i_det];
        const int step = pixel_buf_man.steps[0];
        for (int i_time = 0; i_time < n_time; ++i_time) {
            double coords[4];
            FSIGNAL pf[n_spin];
            pointer.GetCoords(i_det, i_time, (double*)dofs, (double*)coords);
            int pixel_offset = _pixelizor.GetPixel(i_det, i_time, (double*)coords);
            accumulator.SpinProjFactors(coords, pf);

            pix_buf[i_time * step] = pixel_offset;
            for (int i_spin = 0; i_spin < n_spin; i_spin++)
                proj_buf[i_time * proj_buf_man.steps[0] +
                         i_spin * proj_buf_man.steps[1]] = pf[i_spin];
        }
    }

    return bp::make_tuple(pixel_buf_man.ret_val,
                          proj_buf_man.ret_val);
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

    vector<vector<RangesInt32>> ranges;

#pragma omp parallel
    {
        int n_domain = omp_get_num_threads();
#pragma omp single
        {
            for (int i=0; i<n_domain; ++i) {
                vector<RangesInt32> v(n_det);
                for (auto &_v: v)
                    _v.count = n_time;
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
typedef ProjectionEngine<Pointer<ProjFlat>,Pixelizor2_Flat,Accumulator<SpinT,NonTiled>>
  ProjEng_Flat_T;
typedef ProjectionEngine<Pointer<ProjFlat>,Pixelizor2_Flat,Accumulator<SpinQU,NonTiled>>
  ProjEng_Flat_QU;
typedef ProjectionEngine<Pointer<ProjFlat>,Pixelizor2_Flat,Accumulator<SpinTQU,NonTiled>>
  ProjEng_Flat_TQU;

typedef ProjectionEngine<Pointer<ProjFlat>,Pixelizor2_Flat,Accumulator<SpinT,Tiled>>
  ProjEng_Flat_Tiled_T;
typedef ProjectionEngine<Pointer<ProjFlat>,Pixelizor2_Flat,Accumulator<SpinQU,Tiled>>
  ProjEng_Flat_Tiled_QU;
typedef ProjectionEngine<Pointer<ProjFlat>,Pixelizor2_Flat,Accumulator<SpinTQU,Tiled>>
  ProjEng_Flat_Tiled_TQU;

/*
//Cylindrical.
typedef ProjectionEngine<Pointer<ProjCEA>,Pixelizor2_Flat,Accumulator<SpinT,NonTiled>>
  ProjEng_CEA_T;
typedef ProjectionEngine<Pointer<ProjCEA>,Pixelizor2_Flat,Accumulator<SpinQU,NonTiled>>
  ProjEng_CEA_QU;
typedef ProjectionEngine<Pointer<ProjCEA>,Pixelizor2_Flat,Accumulator<SpinTQU,NonTiled>>
  ProjEng_CEA_TQU;
typedef ProjectionEngine<Pointer<ProjCAR>,Pixelizor2_Flat,Accumulator<SpinT,NonTiled>>
  ProjEng_CAR_T;
typedef ProjectionEngine<Pointer<ProjCAR>,Pixelizor2_Flat,Accumulator<SpinQU,NonTiled>>
  ProjEng_CAR_QU;
typedef ProjectionEngine<Pointer<ProjCAR>,Pixelizor2_Flat,Accumulator<SpinTQU,NonTiled>>
  ProjEng_CAR_TQU;

typedef ProjectionEngine<Pointer<ProjCAR>,Pixelizor2_Flat,Accumulator<SpinT,Tiled>>
  ProjEng_CAR_Tiled_T;
typedef ProjectionEngine<Pointer<ProjCAR>,Pixelizor2_Flat,Accumulator<SpinQU,Tiled>>
  ProjEng_CAR_Tiled_QU;
typedef ProjectionEngine<Pointer<ProjCAR>,Pixelizor2_Flat,Accumulator<SpinTQU,Tiled>>
  ProjEng_CAR_Tiled_TQU;

//Zenithal.
typedef ProjectionEngine<Pointer<ProjARC>,Pixelizor2_Flat,Accumulator<SpinT,NonTiled>>
  ProjEng_ARC_T;
typedef ProjectionEngine<Pointer<ProjARC>,Pixelizor2_Flat,Accumulator<SpinQU,NonTiled>>
  ProjEng_ARC_QU;
typedef ProjectionEngine<Pointer<ProjARC>,Pixelizor2_Flat,Accumulator<SpinTQU,NonTiled>>
  ProjEng_ARC_TQU;
typedef ProjectionEngine<Pointer<ProjTAN>,Pixelizor2_Flat,Accumulator<SpinT,NonTiled>>
  ProjEng_TAN_T;
typedef ProjectionEngine<Pointer<ProjTAN>,Pixelizor2_Flat,Accumulator<SpinQU,NonTiled>>
  ProjEng_TAN_QU;
typedef ProjectionEngine<Pointer<ProjTAN>,Pixelizor2_Flat,Accumulator<SpinTQU,NonTiled>>
  ProjEng_TAN_TQU;
typedef ProjectionEngine<Pointer<ProjZEA>,Pixelizor2_Flat,Accumulator<SpinT,NonTiled>>
  ProjEng_ZEA_T;
typedef ProjectionEngine<Pointer<ProjZEA>,Pixelizor2_Flat,Accumulator<SpinQU,NonTiled>>
  ProjEng_ZEA_QU;
typedef ProjectionEngine<Pointer<ProjZEA>,Pixelizor2_Flat,Accumulator<SpinTQU,NonTiled>>
  ProjEng_ZEA_TQU;
*/

/*
 * We also have a generic ProjEng_Precomp, which is basically what you
 * can use if you have precomputed pixel index and spin projection
 * factors for every sample.  This is agnostic of coordinate system,
 * and so not crazily templated.
 */

bp::object ProjEng_Precomp::to_map(
    bp::object map, bp::object pixel_index, bp::object spin_proj,
    bp::object signal, bp::object det_weights)
{
    /* You won't get far without pixel_index, so use that to nail down
     * the n_time and n_det. */
    auto pixel_buf_man = SignalSpace<int32_t>(
        pixel_index, "pixel_index", NPY_INT32, -1, -1);

    int n_det = pixel_buf_man.dims[0];
    int n_time = pixel_buf_man.dims[1];

    auto signal_man = SignalSpace<FSIGNAL>(
        signal, "signal", FSIGNAL_NPY_TYPE, n_det, n_time);

    BufferWrapper<double> map_buf("map", map, false,
                                  vector<int>{-1,-3});
    int n_spin = map_buf->shape[0];

    auto spin_proj_man = SignalSpace<FSIGNAL>(
        spin_proj, "spin_proj", FSIGNAL_NPY_TYPE, n_det, n_time, n_spin);

    BufferWrapper<FSIGNAL> _det_weights("det_weights", det_weights, true,
                                        vector<int>{n_det});

    for (int i_det = 0; i_det < n_det; ++i_det) {
        FSIGNAL det_wt = 1.;
        if (_det_weights->obj != NULL)
            det_wt = *(FSIGNAL*)((char*)_det_weights->buf + _det_weights->strides[0]*i_det);

        for (int i_time = 0; i_time < n_time; ++i_time) {
            const int pixel_index = *(pixel_buf_man.data_ptr[i_det] +
                                      pixel_buf_man.steps[0]*i_time);
            if (pixel_index < 0)
                continue;
            const FSIGNAL sig = *(signal_man.data_ptr[i_det] +
                                  signal_man.steps[0]*i_time);
            for (int i_spin = 0; i_spin < n_spin; ++i_spin) {
                FSIGNAL sp = *(spin_proj_man.data_ptr[i_det] +
                               spin_proj_man.steps[0] * i_time +
                               spin_proj_man.steps[1] * i_spin);
                ((double*)map_buf->buf)[pixel_index +
                                             map_buf->strides[0] * i_spin / sizeof(double)]
                    += det_wt * sig * sp;
            }
        }
    }

    return map;
}


bp::object ProjEng_Precomp::to_weight_map(
    bp::object map, bp::object pixel_index, bp::object spin_proj,
    bp::object signal, bp::object det_weights)
{
    /* You won't get far without pixel_index, so use that to nail down
     * the n_time and n_det. */
    auto pixel_buf_man = SignalSpace<int32_t>(
        pixel_index, "pixel_index", NPY_INT32, -1, -1);

    int n_det = pixel_buf_man.dims[0];
    int n_time = pixel_buf_man.dims[1];

    auto signal_man = SignalSpace<FSIGNAL>(
        signal, "signal", FSIGNAL_NPY_TYPE, n_det, n_time);

    BufferWrapper<double> map_buf("map", map, false,
                                  vector<int>{-1,-3});
    int n_spin = map_buf->shape[0];

    auto spin_proj_man = SignalSpace<FSIGNAL>(
        spin_proj, "spin_proj", FSIGNAL_NPY_TYPE, n_det, n_time, n_spin);

    BufferWrapper<FSIGNAL> _det_weights("det_weights", det_weights, true,
                                        vector<int>{n_det});

    for (int i_det = 0; i_det < n_det; ++i_det) {
        FSIGNAL det_wt = 1.;
        if (_det_weights->obj != NULL)
            det_wt = *(FSIGNAL*)((char*)_det_weights->buf + _det_weights->strides[0]*i_det);

        for (int i_time = 0; i_time < n_time; ++i_time) {
            const int pixel_index = *(pixel_buf_man.data_ptr[i_det] +
                                      pixel_buf_man.steps[0]*i_time);
            if (pixel_index < 0)
                continue;
            const FSIGNAL sig = *(signal_man.data_ptr[i_det] +
                                  signal_man.steps[0]*i_time);
            FSIGNAL *sp = (spin_proj_man.data_ptr[i_det] +
                           spin_proj_man.steps[0] * i_time);
            for (int i_spin = 0; i_spin < n_spin; ++i_spin) {
                for (int j_spin = i_spin; j_spin < n_spin; ++j_spin) {
                    ((double*)map_buf->buf)[pixel_index +
                                                (map_buf->strides[0] * i_spin +
                                                 map_buf->strides[1] * j_spin) / sizeof(double)]
                        += det_wt * sp[spin_proj_man.steps[1] * i_spin] *
                        sp[spin_proj_man.steps[1] * j_spin];
                }
            }
        }
    }

    return map;
}


bp::object ProjEng_Precomp::from_map(
    bp::object map, bp::object pixel_index, bp::object spin_proj,
    bp::object signal, bp::object det_weights)
{
    /* You won't get far without pixel_index, so use that to nail down
     * the n_time and n_det. */
    auto pixel_buf_man = SignalSpace<int32_t>(
        pixel_index, "pixel_index", NPY_INT32, -1, -1);

    int n_det = pixel_buf_man.dims[0];
    int n_time = pixel_buf_man.dims[1];

    auto signal_man = SignalSpace<FSIGNAL>(
        signal, "signal", FSIGNAL_NPY_TYPE, n_det, n_time);

    BufferWrapper<double> map_buf("map", map, false,
                                  vector<int>{-1,-3});
    int n_spin = map_buf->shape[0];

    auto spin_proj_man = SignalSpace<FSIGNAL>(
        spin_proj, "spin_proj", FSIGNAL_NPY_TYPE, n_det, n_time, n_spin);

    BufferWrapper<FSIGNAL> _det_weights("det_weights", det_weights, true,
                                        vector<int>{n_det});

#pragma omp parallel for
    for (int i_det = 0; i_det < n_det; ++i_det) {
        FSIGNAL det_wt = 1.;
        if (_det_weights->obj != NULL)
            det_wt = *(FSIGNAL*)((char*)_det_weights->buf + _det_weights->strides[0]*i_det);

        for (int i_time = 0; i_time < n_time; ++i_time) {
            const int pixel_index = *(pixel_buf_man.data_ptr[i_det] +
                                      pixel_buf_man.steps[0]*i_time);
            if (pixel_index < 0)
                continue;
            FSIGNAL sig = 0;
            for (int i_spin = 0; i_spin < n_spin; ++i_spin) {
                FSIGNAL sp = *(spin_proj_man.data_ptr[i_det] +
                               spin_proj_man.steps[0] * i_time +
                               spin_proj_man.steps[1] * i_spin);
                sig += sp * ((double*)map_buf->buf)[
                    pixel_index + map_buf->strides[0] * i_spin / sizeof(double)];
            }
            *(signal_man.data_ptr[i_det] + signal_man.steps[0]*i_time) += sig;
        }
    }

    return signal_man.ret_val;
}



#define EXPORT_ENGINE(CLASSNAME)                                        \
    bp::class_<CLASSNAME>(#CLASSNAME, bp::init<Pixelizor2_Flat>())      \
    .def("to_map", &CLASSNAME::to_map)                                  \
    .def("to_map_omp", &CLASSNAME::to_map_omp)                          \
    .def("to_weight_map", &CLASSNAME::to_weight_map)                    \
    .def("to_weight_map_omp", &CLASSNAME::to_weight_map_omp)            \
    .def("from_map", &CLASSNAME::from_map)                              \
    .def("coords", &CLASSNAME::coords)                                  \
    .def("pixels", &CLASSNAME::pixels)                                  \
    .def("pointing_matrix", &CLASSNAME::pointing_matrix)                \
    .def("pixel_ranges", &CLASSNAME::pixel_ranges);

#define EXPORT_PRECOMP(CLASSNAME)                                       \
    bp::class_<CLASSNAME>(#CLASSNAME)                                   \
    .def("to_map", &CLASSNAME::to_map)                                  \
    .def("to_weight_map", &CLASSNAME::to_weight_map)                    \
    .def("from_map", &CLASSNAME::from_map);


PYBINDINGS("so3g")
{
    EXPORT_ENGINE(ProjEng_Flat_T);
    EXPORT_ENGINE(ProjEng_Flat_QU);
    EXPORT_ENGINE(ProjEng_Flat_TQU);
    EXPORT_ENGINE(ProjEng_Flat_Tiled_T);
    EXPORT_ENGINE(ProjEng_Flat_Tiled_QU);
    EXPORT_ENGINE(ProjEng_Flat_Tiled_TQU);
//    EXPORT_ENGINE(ProjEng_CAR_T);
//    EXPORT_ENGINE(ProjEng_CAR_QU);
//    EXPORT_ENGINE(ProjEng_CAR_TQU);
//    EXPORT_ENGINE(ProjEng_CEA_T);
//    EXPORT_ENGINE(ProjEng_CEA_QU);
//    EXPORT_ENGINE(ProjEng_CEA_TQU);
//    EXPORT_ENGINE(ProjEng_ARC_T);
//    EXPORT_ENGINE(ProjEng_ARC_QU);
//    EXPORT_ENGINE(ProjEng_ARC_TQU);
//    EXPORT_ENGINE(ProjEng_TAN_T);
//    EXPORT_ENGINE(ProjEng_TAN_QU);
//    EXPORT_ENGINE(ProjEng_TAN_TQU);
//    EXPORT_ENGINE(ProjEng_ZEA_T);
//    EXPORT_ENGINE(ProjEng_ZEA_QU);
//    EXPORT_ENGINE(ProjEng_ZEA_TQU);
//
    EXPORT_PRECOMP(ProjEng_Precomp);

    bp::class_<Pixelizor2_Flat>("Pixelizor2_Flat", bp::init<int,int,double,double,
                          double,double>())
        .def("zeros", &Pixelizor2_Flat::zeros);
}
