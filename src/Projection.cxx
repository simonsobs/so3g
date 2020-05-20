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

inline bool isNone(const bp::object &pyo)
{
    return (pyo.ptr() == Py_None);
}



// Proto-templates.
template <int N>
class SpinClass {
public:
    static const int comp_count = N;
};



// template<CoordSys> options
class ProjQuat;
class ProjFlat;
class ProjCEA;
class ProjCAR;
class ProjARC;
class ProjTAN;
class ProjZEA;

// template<PixelSys> options
// class Pixelizor2_Flat;
// class Pixelizor2_Flat_Tiled;

// template<TilingSys> options
class Tiled;
class NonTiled;

// template<SpinSys> options
class SpinT : public SpinClass<1> {};
class SpinQU : public SpinClass<2> {};
class SpinTQU : public SpinClass<3> {};


template <typename CoordSys>
class Pointer {
public:
    bool TestInputs(bp::object &map, bp::object &pbore, bp::object &pdet,
                    bp::object &signal, bp::object &det_weights);
    void InitPerDet(int i_det, double *dofs);
    int DetCount() { return n_det; }
    int TimeCount() { return n_time; }
    void GetCoords(int i_det, int i_time, const double *dofs, double *coords);
private:
    BufferWrapper<double> _pborebuf;
    BufferWrapper<double> _pdetbuf;
    int n_det;
    int n_time;
};


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


//! Pixelizors
//
template <typename TilingSys>
class Pixelizor2_Flat;

template <>
class Pixelizor2_Flat<NonTiled> {
public:
    static const int index_count = 2;
    Pixelizor2_Flat(int ny, int nx,
                    double dy=1., double dx=1.,
                    double iy0=0., double ix0=0.) {
        std::cout << ny << " "
                  << nx << "\n";

        naxis[0] = ny;
        naxis[1] = nx;
        cdelt[0] = dy;
        cdelt[1] = dx;
        crpix[0] = iy0;
        crpix[1] = ix0;
    };
    Pixelizor2_Flat() : naxis{1,1} {};
    Pixelizor2_Flat(bp::object args) {
        bp::tuple args_tuple = bp::extract<bp::tuple>(args);
        naxis[0] = bp::extract<int>(args_tuple[0])();
        naxis[1] = bp::extract<int>(args_tuple[1])();
        cdelt[0] = bp::extract<double>(args_tuple[2])();
        cdelt[1] = bp::extract<double>(args_tuple[3])();
        crpix[0] = bp::extract<double>(args_tuple[4])();
        crpix[1] = bp::extract<double>(args_tuple[5])();
    }
    ~Pixelizor2_Flat() {};
    // bool TestInputs(bp::object &map, bp::object &pbore, bp::object &pdet,
    //                 bp::object &signal, bp::object &det_weights) {
    //     // Flat pixelizor.  Requires:
    //     // 1. Map has the right shape -- but map can be none.
    //     // 2. Coords have the right format (unchecked currently).

    //     // If map is provided, it can have any number of leading
    //     // dimensions but the last two dimensions must match naxis.
    //     BufferWrapper<double> mapbuf("map", map, true,
    //                                  vector<int>{-2,naxis[0],naxis[1]});
    //     if (mapbuf->buf != NULL) {
    //         // Note these are byte offsets, not index.
    //         // int ndim = mapbuf->ndim;
    //         // strides[0] = mapbuf->strides[ndim-2];
    //         // strides[1] = mapbuf->strides[ndim-1];
    //     } else {
    //         // Set it up to return naive C-ordered pixel indices.
    //         // strides[0] = naxis[1];
    //         // strides[1] = 1;
    //     }

    //     // [Check/declare coords format.]
    //     return true;
    // }

    bp::object zeros(vector<int> shape)
        {
            shape.push_back(naxis[0]);
            shape.push_back(naxis[1]);
            int ndim = 0;
            npy_intp dims[32];
            for (auto d: shape)
                dims[ndim++] = d;
            // while (ndim <
            // for (int i; dimi < shape.size(); dimi++)
            //     dims[dimi] = shape[dimi];
            //     size *= shape[dimi];
            // }

            // dims[dimi++] = naxis[0];
            // dims[dimi++] = naxis[1];

            int dtype = NPY_FLOAT64;
            PyObject *v = PyArray_ZEROS(ndim, dims, dtype, 0);
            return bp::object(bp::handle<>(v));
        }

    inline
    void GetPixel(int i_det, int i_time, const double *coords, int *pixel_index)
        {
            double ix = coords[0] / cdelt[1] + crpix[1] + 0.5;
            if (ix < 0 || ix >= naxis[1]) {
                *pixel_index = -1;
                return;
            }

            double iy = coords[1] / cdelt[0] + crpix[0] + 0.5;
            if (iy < 0 || iy >= naxis[0]) {
                *pixel_index = -1;
                return;
            }

            pixel_index[0] = int(iy);
            pixel_index[1] = int(ix);
            // *pixel_index =  strides[0]*int(iy) + strides[1]*int(ix);
        }

    int crpix[2];
    double cdelt[2];
    int naxis[2];

    // From Tiled/NonTiled
    BufferWrapper<double> mapbuf;
    bool TestInputs(bp::object &map, bool need_map, bool need_weight_map, int comp_count) {
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
    double *pix(int imap, const int pixel_index[]) {
        return (double*)((char*)mapbuf->buf +
                         mapbuf->strides[0]*imap +
                         mapbuf->strides[1]*pixel_index[0] +
                         mapbuf->strides[2]*pixel_index[1]);
    }
    double *wpix(int imap, int jmap, const int pixel_index[]) {
        return (double*)((char*)mapbuf->buf +
                         mapbuf->strides[0]*imap +
                         mapbuf->strides[1]*jmap +
                         mapbuf->strides[2]*pixel_index[0] +
                         mapbuf->strides[3]*pixel_index[1]);
    }
    int stripe(const int pixel_index[], int thread_count) {
        if (pixel_index[0] < 0)
            return -1;
        return pixel_index[0] * thread_count / (mapbuf->shape[0] / mapbuf->itemsize);
    }
};

template <>
class Pixelizor2_Flat<Tiled> {
public:
    static const int index_count = 3;
    Pixelizor2_Flat() {};
    ~Pixelizor2_Flat() {};
    Pixelizor2_Flat(int ny, int nx, double dy, double dx,
                    double iy0, double ix0, int tiley, int tilex) {
        parent_pix = Pixelizor2_Flat<NonTiled>(ny, nx, dy, dx, iy0, ix0);
        tile_shape[0] = tiley;
        tile_shape[1] = tilex;
    };
    Pixelizor2_Flat(Pixelizor2_Flat<NonTiled> _parent_pix, int tiley, int tilex) {
        parent_pix = _parent_pix;
        tile_shape[0] = tiley;
        tile_shape[1] = tilex;
    };
    Pixelizor2_Flat(bp::object args) {
        bp::tuple args_tuple = bp::extract<bp::tuple>(args);
        Pixelizor2_Flat(bp::extract<int>(args_tuple[0]),
                        bp::extract<int>(args_tuple[1]),
                        bp::extract<double>(args_tuple[2]),
                        bp::extract<double>(args_tuple[3]),
                        bp::extract<double>(args_tuple[4]),
                        bp::extract<double>(args_tuple[5]),
                        bp::extract<int>(args_tuple[6]),
                        bp::extract<int>(args_tuple[7]));
    }

    bp::object zeros(int count) {
        // If a mapping routine tries to generate an empty map, on the
        // fly, don't let it.
        throw shape_exception("Flat_Tiled", "Tiled output must be pre-initialized.");
    }

    bp::object zeros_tiled(int count, bp::object active_tiles) {
        int base_size = 1;
        int base_dimi = 0;
        npy_intp dims[32];
        int dtype = NPY_FLOAT64;

        if (count >= 0) {
            dims[base_dimi++] = count;
            base_size *= count;
        }

        int n_ty = (parent_pix.naxis[0] + tile_shape[0] - 1) / tile_shape[0];
        int n_tx = (parent_pix.naxis[1] + tile_shape[1] - 1) / tile_shape[1];
        vector<bool> populate(n_ty * n_tx, false);
        if (isNone(active_tiles)) {
            populate = vector<bool>(populate.size(), true);
        } else {
            for (int i=0; i<bp::len(active_tiles); i++) {
                // We're using C API instead of boost because this:
                //   bp::extract<int>(active_tiles[i])
                // does not work on an ndarray, nor even on a list of
                // elements extracted from an array, unless one carefully
                // casts them to int first.
                int idx = PyLong_AsLong(bp::object(active_tiles[i]).ptr());
                if (idx >= 0 && idx < n_tx*n_ty)
                    populate[idx] = true;
            }
        }

        bp::list maps_out;
        for (int i_ty = 0; i_ty < n_ty; i_ty++) {
            dims[base_dimi] = min(tile_shape[0], parent_pix.naxis[0] - i_ty * tile_shape[0]);
            for (int i_tx = 0; i_tx < n_tx; i_tx++) {
                dims[base_dimi+1] = min(tile_shape[1], parent_pix.naxis[1] - i_tx * tile_shape[1]);
                if (populate[i_ty * n_tx + i_tx]) {
                    PyObject *v = PyArray_ZEROS(base_dimi+2, dims, dtype, 0);
                    maps_out.append(bp::handle<>(v));
                } else
                    maps_out.append(bp::object());
            }
        }
        return maps_out;
    }

    inline
    void GetPixel(int i_det, int i_time, const double *coords, int *pixel_index) {
        double ix = coords[0] / parent_pix.cdelt[1] + parent_pix.crpix[1] + 0.5;
        if (ix < 0 || ix >= parent_pix.naxis[1]) {
            pixel_index[0] = -1;
            return;
        }

        double iy = coords[1] / parent_pix.cdelt[0] + parent_pix.crpix[0] + 0.5;
        if (iy < 0 || iy >= parent_pix.naxis[0]) {
            *pixel_index = -1;
            return;
        }

        int sub_y = int(iy) / tile_shape[0];
        int sub_x = int(ix) / tile_shape[1];
        pixel_index[0] = sub_y * ((parent_pix.naxis[1] + tile_shape[1] - 1) / tile_shape[1]) + sub_x;
        //pixel_index[1] = strides[0]*int(iy) + strides[1]*int(ix);
        pixel_index[1] = int(iy) - sub_y * tile_shape[0];
        pixel_index[2] = int(ix) - sub_x * tile_shape[1];
    }

    Pixelizor2_Flat<NonTiled> parent_pix;

    // From Tiled/NonTiled
    vector<BufferWrapper<double>> mapbufs;
    vector<int> shape;
    vector<int> tile_shape;

    bool TestInputs(bp::object &map, bool need_map, bool need_weight_map, int comp_count) {
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
        if (map_shape_req.size() == 0)
            return true;

        // The user must pass in the following, as a list or tuple:
        //     [(shape, tile_shape), [map0, ..., mapN-1]]
        //
        // The shape and tile_shape are both two-tuples of ints.  The
        // list of maps must be complete (one entry per tile) but can
        // have null entries.

        bp::object tiling_info = map[0];
        bp::object map_list = map[1];
        shape.clear();
        shape.push_back(bp::extract<int>(tiling_info[0][0]));
        shape.push_back(bp::extract<int>(tiling_info[0][1]));
        tile_shape.clear();
        tile_shape.push_back(bp::extract<int>(tiling_info[1][0]));
        tile_shape.push_back(bp::extract<int>(tiling_info[1][1]));
        
        // int ntile0 = (shape[0] + tile_shape[0] - 1) / tile_shape[0];
        // int ntile1 = (shape[1] + tile_shape[1] - 1) / tile_shape[1];
        // for (int i0; i0<ntile0; i++) { }
        mapbufs.clear();
        for (int i_tile = 0; i_tile < bp::len(map_list); i_tile++) {
            if (isNone(map_list[i_tile])) {
                mapbufs.push_back(BufferWrapper<double>());
            } else {
                mapbufs.push_back(
                    BufferWrapper<double>("map", map_list[i_tile], false, map_shape_req));
            }
        }

        return true;
    }
    double *pix(int imap, const int pixel_index[]) {
        BufferWrapper<double> mapbuf = mapbufs[pixel_index[0]];
        return (double*)((char*)mapbuf->buf +
                         mapbuf->strides[0]*imap +
                         mapbuf->strides[1]*pixel_index[1] +
                         mapbuf->strides[2]*pixel_index[2]);
    }
    double *wpix(int imap, int jmap, const int pixel_index[]) {
        // Expensive copy !
        BufferWrapper<double> mapbuf = mapbufs[pixel_index[0]];
        return (double*)((char*)mapbuf->buf +
                         mapbuf->strides[0]*imap +
                         mapbuf->strides[1]*jmap +
                         mapbuf->strides[2]*pixel_index[1] +
                         mapbuf->strides[3]*pixel_index[2]);
    }
    int stripe(const int pixel_index[], int thread_count) {
        return -1;
    }
};

/*
// Tiled Pixelizor

Pixelizor2_Flat_Tiled::Pixelizor2_Flat_Tiled(
    int ny, int nx,
    double dy, double dx,
    double iy0, double ix0,
    int tiley, int tilex)
{
    parent_pix = Pixelizor2_Flat(ny, nx, dy, dx, iy0, ix0);
    tile_shape[0] = tiley;
    tile_shape[1] = tilex;
}

Pixelizor2_Flat_Tiled::Pixelizor2_Flat_Tiled(
    Pixelizor2_Flat _parent_pix,
    int tiley, int tilex)
{
    parent_pix = _parent_pix;
    tile_shape[0] = tiley;
    tile_shape[1] = tilex;
}

bool Pixelizor2_Flat_Tiled::TestInputs(bp::object &map, bp::object &pbore, bp::object &pdet,
                             bp::object &signal, bp::object &det_weights)
{
    // Flat pixelizor.  Requires:
    // - If map is not None, it must be a list of the right length
    //   with the right shapes.

    if (isNone(map))
        return true;

    // [Check/declare coords format.]
    return true;
}

bp::object Pixelizor2_Flat_Tiled::zeros(int count)
{
    // If a mapping routine tries to generate an empty map, on the
    // fly, don't let it.
    throw shape_exception("Flat_Tiled", "Tiled output must be pre-initialized.");
}

bp::object Pixelizor2_Flat_Tiled::zeros_tiled(int count, bp::object active_tiles)
{
    int base_size = 1;
    int base_dimi = 0;
    npy_intp dims[32];
    int dtype = NPY_FLOAT64;
    
    if (count >= 0) {
        dims[base_dimi++] = count;
        base_size *= count;
    }

    int n_ty = (parent_pix.naxis[0] + tile_shape[0] - 1) / tile_shape[0];
    int n_tx = (parent_pix.naxis[1] + tile_shape[1] - 1) / tile_shape[1];
    vector<bool> populate(n_ty * n_tx, false);
    if (isNone(active_tiles)) {
        populate = vector<bool>(populate.size(), true);
    } else {
        for (int i=0; i<bp::len(active_tiles); i++) {
            // We're using C API instead of boost because this:
            //    bp::extract<int>(active_tiles[i])
            // does not work on an ndarray, nor even on a list of
            // elements extracted from an array, unless one carefully
            // casts them to int first.
            int idx = PyLong_AsLong(bp::object(active_tiles[i]).ptr());
            if (idx >= 0 && idx < n_tx*n_ty)
                populate[idx] = true;
        }
    }
        
    bp::list maps_out;
    for (int i_ty = 0; i_ty < n_ty; i_ty++) {
        dims[base_dimi] = min(tile_shape[0], parent_pix.naxis[0] - i_ty * tile_shape[0]);
        for (int i_tx = 0; i_tx < n_tx; i_tx++) {
            dims[base_dimi+1] = min(tile_shape[1], parent_pix.naxis[1] - i_tx * tile_shape[1]);
            if (populate[i_ty * n_tx + i_tx]) {
                PyObject *v = PyArray_ZEROS(base_dimi+2, dims, dtype, 0);
                maps_out.append(bp::handle<>(v));
            } else
                maps_out.append(bp::object());
        }
    }
    return maps_out;
}


inline
void Pixelizor2_Flat_Tiled::GetPixel(int i_det, int i_time, const double *coords, int *pixel_index)
{
    double ix = coords[0] / parent_pix.cdelt[1] + parent_pix.crpix[1] + 0.5;
    if (ix < 0 || ix >= parent_pix.naxis[1]) {
        pixel_index[0] = -1;
        return;
    }

    double iy = coords[1] / parent_pix.cdelt[0] + parent_pix.crpix[0] + 0.5;
    if (iy < 0 || iy >= parent_pix.naxis[0]) {
        *pixel_index = -1;
        return;
    }

    int sub_y = int(iy) / tile_shape[0];
    int sub_x = int(ix) / tile_shape[1];
    pixel_index[0] = sub_y * ((parent_pix.naxis[1] + tile_shape[1] - 1) / tile_shape[1]) + sub_x;
    //pixel_index[1] = strides[0]*int(iy) + strides[1]*int(ix);
    pixel_index[1] = int(iy) - sub_y * tile_shape[0];
    pixel_index[2] = int(ix) - sub_x * tile_shape[1];
}
*/

/** Accumulator - transfer signal from map domain to time domain.
 *
 */

class Tiled {
public:
    vector<BufferWrapper<double>> mapbufs;
    vector<int> shape;
    vector<int> tile_shape;
    
    bool TestInputs(bp::object &map, bool need_map, bool need_weight_map, int comp_count) {
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
        if (map_shape_req.size() == 0)
            return true;

        // The user must pass in the following, as a list or tuple:
        //     [(shape, tile_shape), [map0, ..., mapN-1]]
        //
        // The shape and tile_shape are both two-tuples of ints.  The
        // list of maps must be complete (one entry per tile) but can
        // have null entries.

        bp::object tiling_info = map[0];
        bp::object map_list = map[1];
        shape.clear();
        shape.push_back(bp::extract<int>(tiling_info[0][0]));
        shape.push_back(bp::extract<int>(tiling_info[0][1]));
        tile_shape.clear();
        tile_shape.push_back(bp::extract<int>(tiling_info[1][0]));
        tile_shape.push_back(bp::extract<int>(tiling_info[1][1]));
        
        // int ntile0 = (shape[0] + tile_shape[0] - 1) / tile_shape[0];
        // int ntile1 = (shape[1] + tile_shape[1] - 1) / tile_shape[1];
        // for (int i0; i0<ntile0; i++) { }
        mapbufs.clear();
        for (int i_tile = 0; i_tile < bp::len(map_list); i_tile++) {
            if (isNone(map_list[i_tile])) {
                mapbufs.push_back(BufferWrapper<double>());
            } else {
                mapbufs.push_back(
                    BufferWrapper<double>("map", map_list[i_tile], false, map_shape_req));
            }
        }

        return true;
    }        
    double *pix(int imap, const int pixel_index[]) {
        BufferWrapper<double> mapbuf = mapbufs[pixel_index[0]];
        return (double*)((char*)mapbuf->buf +
                         mapbuf->strides[0]*imap +
                         mapbuf->strides[1]*pixel_index[1] +
                         mapbuf->strides[2]*pixel_index[2]);
    }
    double *wpix(int imap, int jmap, const int pixel_index[]) {
        // Expensive copy !
        BufferWrapper<double> mapbuf = mapbufs[pixel_index[0]];
        return (double*)((char*)mapbuf->buf +
                         mapbuf->strides[0]*imap +
                         mapbuf->strides[1]*jmap +
                         mapbuf->strides[2]*pixel_index[1] +
                         mapbuf->strides[3]*pixel_index[2]);
    }
    int stripe(const int pixel_index[], int thread_count) {
        if (pixel_index[0] < 0)
            return -1;
        return pixel_index[0] % thread_count;
    }
};


template <typename SpinSys>
inline
void spin_proj_factors(const double* coords, FSIGNAL *projfacs);


template <>
inline void spin_proj_factors<SpinT>(const double* coords, FSIGNAL *projfacs)
{
     projfacs[0] = 1.;
}

template <>
inline void spin_proj_factors<SpinQU>(const double* coords, FSIGNAL *projfacs)
{
    const double c = coords[2];
    const double s = coords[3];
    projfacs[0] = c*c - s*s;
    projfacs[1] = 2*c*s;
}

template <>
inline void spin_proj_factors<SpinTQU>(const double* coords, FSIGNAL *projfacs)
{
     const double c = coords[2];
     const double s = coords[3];
     projfacs[0] = 1.;
     projfacs[1] = c*c - s*s;
     projfacs[2] = 2*c*s;
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

template<typename C, typename P, typename T, typename S>
ProjectionEngine<C,P,T,S>::ProjectionEngine(bp::object pix_args)
{
    _pixelizor = P(pix_args);
}


/*
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
                    int pixel_offset[Z::index_count];
                    pointer.GetCoords(i_det, i_time, (double*)dofs, (double*)coords);
                    _pixelizor.GetPixel(i_det, i_time, (double*)coords, pixel_offset);
                    accumulator.Forward(i_det, i_time, pixel_offset, coords, weights);
                }
            }
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
        map = _pixelizor.zeros(std::vector<int>{n_comp,n_comp});
        // auto v0 = (PyArrayObject*)map.ptr();
        // npy_intp dims[32] = {n_comp, n_comp};
        // int dimi = 2;
        // for (int d=1; d<PyArray_NDIM(v0); d++)
        //     dims[dimi++] = PyArray_DIM(v0, d);
        // PyArray_Dims padims = {dims, dimi};
        // PyObject *v1 = PyArray_Newshape(v0, &padims,  NPY_ANYORDER);
        // map = bp::object(bp::handle<>(v1));
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
                    int pixel_offset[Z::index_count];
                    pointer.GetCoords(i_det, i_time, (double*)dofs, (double*)coords);
                    _pixelizor.GetPixel(i_det, i_time, (double*)coords, pixel_offset);
                    accumulator.ForwardWeight(i_det, i_time, pixel_offset, coords, weights);
                }
            }
        }
    }

    return map;
}
*/

template<typename C, typename P, typename T, typename S>
int ProjectionEngine<C,P,T,S>::index_count() const {
    return P::index_count;
}

template<typename C, typename P, typename T, typename S>
int ProjectionEngine<C,P,T,S>::comp_count() const {
    return S::comp_count;
}

template<typename C, typename P, typename T, typename S>
bp::object ProjectionEngine<C,P,T,S>::coords(
    bp::object pbore, bp::object pofs, bp::object coord)
{
    auto _none = bp::object();
    auto pointer = Pointer<C>();
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

template<typename C, typename P, typename T, typename S>
bp::object ProjectionEngine<C,P,T,S>::pixels(
    bp::object pbore, bp::object pofs, bp::object pixel)
{
    auto _none = bp::object();

    auto pointer = Pointer<C>();
    pointer.TestInputs(_none, pbore, pofs, _none, _none);
    //_pixelizor.TestInputs(_none, _none, _none, _none, _none);

    int n_det = pointer.DetCount();
    int n_time = pointer.TimeCount();

    auto pixel_buf_man = SignalSpace<int32_t>(
        pixel, "pixel", NPY_INT32, n_det, n_time, P::index_count);

#pragma omp parallel for
    for (int i_det = 0; i_det < n_det; ++i_det) {
        double dofs[4];
        pointer.InitPerDet(i_det, dofs);
        int* const pix_buf = pixel_buf_man.data_ptr[i_det];
        int pixel_offset[P::index_count];
        for (int i_time = 0; i_time < n_time; ++i_time) {
            double coords[4];
            pointer.GetCoords(i_det, i_time, (double*)dofs, (double*)coords);
            _pixelizor.GetPixel(i_det, i_time, (double*)coords, pixel_offset);
            // Fix me -- multi-dim!
            for (int i_dim = 0; i_dim < P::index_count; i_dim++)
                pix_buf[i_time * pixel_buf_man.steps[0] +
                        i_dim * pixel_buf_man.steps[1]] = pixel_offset[i_dim];
            // pix_buf[i_time * pixel_buf_man.steps[0]] = pixel_offset;
        }
    }

    return pixel_buf_man.ret_val;
}

template<typename C, typename P, typename T, typename S>
bp::object ProjectionEngine<C,P,T,S>::pointing_matrix(
    bp::object pbore, bp::object pofs, bp::object pixel, bp::object proj)
{
    auto _none = bp::object();

    auto pointer = Pointer<C>();
    pointer.TestInputs(_none, pbore, pofs, _none, _none);
    //_pixelizor.TestInputs(_none, _none, _none, _none, _none);

    int n_det = pointer.DetCount();
    int n_time = pointer.TimeCount();

    auto pixel_buf_man = SignalSpace<int32_t>(
        pixel, "pixel", NPY_INT32, n_det, n_time, P::index_count);

    auto proj_buf_man = SignalSpace<FSIGNAL>(
        proj, "proj", FSIGNAL_NPY_TYPE, n_det, n_time, S::comp_count);

#pragma omp parallel for
    for (int i_det = 0; i_det < n_det; ++i_det) {
        double dofs[4];
        pointer.InitPerDet(i_det, dofs);
        int* const pix_buf = pixel_buf_man.data_ptr[i_det];
        FSIGNAL* const proj_buf = proj_buf_man.data_ptr[i_det];
        const int step = pixel_buf_man.steps[0];
        int pixel_offset[P::index_count];
        for (int i_time = 0; i_time < n_time; ++i_time) {
            double coords[4];
            FSIGNAL pf[S::comp_count];
            pointer.GetCoords(i_det, i_time, (double*)dofs, (double*)coords);
            _pixelizor.GetPixel(i_det, i_time, (double*)coords, pixel_offset);
            spin_proj_factors<S>(coords, pf);

            for (int i_dim = 0; i_dim < P::index_count; i_dim++)
                pix_buf[i_time * pixel_buf_man.steps[0] +
                        i_dim * pixel_buf_man.steps[1]] = pixel_offset[i_dim];
            for (int i_spin = 0; i_spin < S::comp_count; i_spin++)
                proj_buf[i_time * proj_buf_man.steps[0] +
                         i_spin * proj_buf_man.steps[1]] = pf[i_spin];
        }
    }

    return bp::make_tuple(pixel_buf_man.ret_val,
                          proj_buf_man.ret_val);
}

/*
template<typename P, typename Z, typename A>
bp::object ProjectionEngine<P,Z,A>::pixel_ranges(
    bp::object pbore, bp::object pofs)
{
    auto pointer = P();
    auto accumulator = A();
    auto _none = bp::object();

    pointer.TestInputs(_none, pbore, pofs, _none, _none);
    _pixelizor.TestInputs(_none, _none, _none, _none, _none);
    //accumulator.TestInputs(_none, pbore, pofs, _none, _none);

    int n_det = pointer.DetCount();
    int n_time = pointer.TimeCount();

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

#pragma omp for
        for (int i_det = 0; i_det < n_det; ++i_det) {
            double dofs[4];
            pointer.InitPerDet(i_det, dofs);
            int last_slice = -1;
            int slice_start = 0;
            int pixel_offset[Z::index_count];
            for (int i_time = 0; i_time < n_time; ++i_time) {
                double coords[4];
                pointer.GetCoords(i_det, i_time, (double*)dofs, (double*)coords);
                _pixelizor.GetPixel(i_det, i_time, (double*)coords, pixel_offset);
                int this_slice = accumulator.Stripe(pixel_offset, n_domain);
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
*/

template<typename C, typename P, typename T, typename S>
bp::object ProjectionEngine<C,P,T,S>::zeros(bp::object shape)
{
    vector<int> dims;
    bp::extract<int> int_ex(shape);
    if (int_ex.check()) {
        dims.push_back(int_ex());
        return _pixelizor.zeros(dims);
    }

    bp::extract<bp::tuple> tuple_ex(shape);
    if (tuple_ex.check()) {
        auto tuple = tuple_ex();
        for (int i=0; i<bp::len(tuple); i++)
            dims.push_back(bp::extract<int>(tuple[i])());
        return _pixelizor.zeros(dims);
    }
    
    return bp::object();  //None on fall-through
}

template<typename C, typename P, typename T, typename S>
bp::object ProjectionEngine<C,P,T,S>::from_map(
    bp::object map, bp::object pbore, bp::object pofs, bp::object signal, bp::object det_weights)
{
    // Initialize pointer and _pixelizor.
    auto pointer = Pointer<C>();
    pointer.TestInputs(map, pbore, pofs, signal, det_weights);
    int n_det = pointer.DetCount();
    int n_time = pointer.TimeCount();

    // Confirm that map has the right meta-shape.
    _pixelizor.TestInputs(map, true, false, S::comp_count);

    // Get pointers to the signal and (optional) per-det weights.
    auto _signalspace = new SignalSpace<FSIGNAL>(
            signal, "signal", FSIGNAL_NPY_TYPE, n_det, n_time);
    auto _det_weights = BufferWrapper<FSIGNAL>(
         "det_weights", det_weights, true, vector<int>{n_det});

#pragma omp parallel for
    for (int i_det = 0; i_det < n_det; ++i_det) {
        double dofs[4];
        pointer.InitPerDet(i_det, dofs);
        int pixel_offset[P::index_count];
        FSIGNAL pf[S::comp_count]; 
        for (int i_time = 0; i_time < n_time; ++i_time) {
            double coords[4];
            pointer.GetCoords(i_det, i_time, (double*)dofs, (double*)coords);
            _pixelizor.GetPixel(i_det, i_time, (double*)coords, pixel_offset);
            if (pixel_offset[0] < 0) continue;
            spin_proj_factors<S>(coords, pf);
            FSIGNAL *sig = (_signalspace->data_ptr[i_det] +
                            _signalspace->steps[0]*i_time);
            for (int imap=0; imap<S::comp_count; ++imap)
                *sig += *_pixelizor.pix(imap, pixel_offset) * pf[imap];
        }
    }

    return _signalspace->ret_val;
}


template<typename C, typename P, typename T, typename S>
bp::object ProjectionEngine<C,P,T,S>::to_map(
    bp::object map, bp::object pbore, bp::object pofs, bp::object signal, bp::object det_weights)
{
    //Initialize it / check inputs.
    auto pointer = Pointer<C>();
    pointer.TestInputs(map, pbore, pofs, signal, det_weights);
    int n_det = pointer.DetCount();
    int n_time = pointer.TimeCount();

    //Do we need a map?  Now is the time.
    if (isNone(map))
        map = _pixelizor.zeros(vector<int>{S::comp_count});

    // Confirm that map has the right meta-shape.
    _pixelizor.TestInputs(map, true, false, S::comp_count);

    // Get pointers to the signal and (optional) per-det weights.
    auto _signalspace = new SignalSpace<FSIGNAL>(
            signal, "signal", FSIGNAL_NPY_TYPE, n_det, n_time);
    auto _det_weights = BufferWrapper<FSIGNAL>(
         "det_weights", det_weights, true, vector<int>{n_det});

    for (int i_det = 0; i_det < n_det; ++i_det) {
        double dofs[4];
        pointer.InitPerDet(i_det, dofs);
        for (int i_time = 0; i_time < n_time; ++i_time) {
            double coords[4];
            int pixel_offset[P::index_count];
            FSIGNAL pf[S::comp_count];
            FSIGNAL det_wt = 1.;

            pointer.GetCoords(i_det, i_time, (double*)dofs, (double*)coords);
            _pixelizor.GetPixel(i_det, i_time, (double*)coords, pixel_offset);

            if (pixel_offset[0] < 0) continue;
            const FSIGNAL sig = *(_signalspace->data_ptr[i_det] +
                                  _signalspace->steps[0]*i_time);
            if (_det_weights->obj != NULL)
                det_wt = *(FSIGNAL*)((char*)_det_weights->buf +
                                     _det_weights->strides[0]*i_det);
            spin_proj_factors<S>(coords, pf);
            for (int imap=0; imap<S::comp_count; ++imap)
                *_pixelizor.pix(imap, pixel_offset) += sig * pf[imap] * det_wt;
        }
    }

    return map;
}

template<typename C, typename P, typename T, typename S>
bp::object ProjectionEngine<C,P,T,S>::to_weight_map(
    bp::object map, bp::object pbore, bp::object pofs, bp::object signal, bp::object det_weights)
{
    //Initialize it / check inputs.
    auto pointer = Pointer<C>();
    pointer.TestInputs(map, pbore, pofs, signal, det_weights);
    int n_det = pointer.DetCount();
    int n_time = pointer.TimeCount();

    //Do we need a map?  Now is the time.
    if (isNone(map))
        map = _pixelizor.zeros(vector<int>{S::comp_count,S::comp_count});

    // Confirm that map has the right meta-shape.
    _pixelizor.TestInputs(map, false, true, S::comp_count);

    // Get pointer to (optional) per-det weights.
    auto _det_weights = BufferWrapper<FSIGNAL>(
         "det_weights", det_weights, true, vector<int>{n_det});

    for (int i_det = 0; i_det < n_det; ++i_det) {
        double dofs[4];
        pointer.InitPerDet(i_det, dofs);
        for (int i_time = 0; i_time < n_time; ++i_time) {
            double coords[4];
            int pixel_offset[P::index_count];
            FSIGNAL pf[S::comp_count];
            FSIGNAL det_wt = 1.;

            pointer.GetCoords(i_det, i_time, (double*)dofs, (double*)coords);
            _pixelizor.GetPixel(i_det, i_time, (double*)coords, pixel_offset);

            if (pixel_offset[0] < 0) continue;
            if (_det_weights->obj != NULL)
                det_wt = *(FSIGNAL*)((char*)_det_weights->buf + _det_weights->strides[0]*i_det);

            spin_proj_factors<S>(coords, pf);
            for (int imap=0; imap<S::comp_count; ++imap)
                for (int jmap=imap; jmap<S::comp_count; ++jmap)
                    *_pixelizor.wpix(imap, jmap, pixel_offset) += pf[imap] * pf[jmap] * det_wt;
        }
    }

    return map;
}


////Flat.
//typedef ProjectionEngine<Pointer<ProjFlat>,Pixelizor2_Flat,Accumulator<SpinT,NonTiled>>
//  ProjEng_Flat_T;
//typedef ProjectionEngine<Pointer<ProjFlat>,Pixelizor2_Flat,Accumulator<SpinQU,NonTiled>>
//  ProjEng_Flat_QU;
//typedef ProjectionEngine<Pointer<ProjFlat>,Pixelizor2_Flat,Accumulator<SpinTQU,NonTiled>>
//  ProjEng_Flat_TQU;
//
//typedef ProjectionEngine<Pointer<ProjFlat>,Pixelizor2_Flat_Tiled,Accumulator<SpinT,Tiled>>
//  ProjEng_Flat_Tiled_T;
//typedef ProjectionEngine<Pointer<ProjFlat>,Pixelizor2_Flat_Tiled,Accumulator<SpinQU,Tiled>>
//  ProjEng_Flat_Tiled_QU;
//typedef ProjectionEngine<Pointer<ProjFlat>,Pixelizor2_Flat_Tiled,Accumulator<SpinTQU,Tiled>>
//  ProjEng_Flat_Tiled_TQU;
//
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
/*
template<typename Z, typename T>
class ProjEng_Precomp {
public:
    ProjEng_Precomp() {};
    bp::object to_map(bp::object map, bp::object pixel_index, bp::object spin_proj,
                      bp::object signal, bp::object weights);
    bp::object to_weight_map(bp::object map, bp::object pixel_index, bp::object spin_proj,
                             bp::object signal, bp::object weights);
    bp::object from_map(bp::object map, bp::object pixel_index, bp::object spin_proj,
                        bp::object signal, bp::object weights);
};

template <typename Z, typename T>
bp::object ProjEng_Precomp<Z,T>::to_map(
    bp::object map, bp::object pixel_index, bp::object spin_proj,
    bp::object signal, bp::object det_weights)
{
    // You won't get far without pixel_index, so use that to nail down
    // the n_time and n_det.
    auto pixel_buf_man = SignalSpace<int32_t>(
        pixel_index, "pixel_index", NPY_INT32, -1, -1, Z::index_count);
    int n_det = pixel_buf_man.dims[0];
    int n_time = pixel_buf_man.dims[1];

    // Similarly, spin_proj tells you the number of components.
    auto spin_proj_man = SignalSpace<FSIGNAL>(
        spin_proj, "spin_proj", FSIGNAL_NPY_TYPE, n_det, n_time, -1);
    int n_spin = spin_proj_man.dims[2];

    auto tiling = T();
    tiling.TestInputs(map, true, false, n_spin);

    // Check the signal.
    auto signal_man = SignalSpace<FSIGNAL>(
        signal, "signal", FSIGNAL_NPY_TYPE, n_det, n_time);

    BufferWrapper<FSIGNAL> _det_weights("det_weights", det_weights, true,
                                        vector<int>{n_det});

    // Below we assume the pixel sub-indices are close-packed.
    if (pixel_buf_man.steps[1] != 1)
        throw shape_exception("pixel_index",
                              "Fast dimension of pixel indices must be close-packed.");

    for (int i_det = 0; i_det < n_det; ++i_det) {
        FSIGNAL det_wt = 1.;
        if (_det_weights->obj != NULL)
            det_wt = *(FSIGNAL*)((char*)_det_weights->buf + _det_weights->strides[0]*i_det);

        for (int i_time = 0; i_time < n_time; ++i_time) {
            const int *pixel_ofs = pixel_buf_man.data_ptr[i_det] +
                pixel_buf_man.steps[0]*i_time;
            if (pixel_ofs[0] < 0)
                continue;
            const FSIGNAL sig = *(signal_man.data_ptr[i_det] +
                                  signal_man.steps[0]*i_time);
            for (int i_spin = 0; i_spin < n_spin; ++i_spin) {
                const FSIGNAL sp = *(spin_proj_man.data_ptr[i_det] +
                                     spin_proj_man.steps[0] * i_time +
                                     spin_proj_man.steps[1] * i_spin);
                *tiling.pix(i_spin, pixel_ofs) += det_wt * sig * sp;
            }
        }
    }

    return map;
}

template <typename Z, typename T>
bp::object ProjEng_Precomp<Z,T>::to_weight_map(
    bp::object map, bp::object pixel_index, bp::object spin_proj,
    bp::object signal, bp::object det_weights)
{
    // You won't get far without pixel_index, so use that to nail down
    // the n_time and n_det.
    auto pixel_buf_man = SignalSpace<int32_t>(
        pixel_index, "pixel_index", NPY_INT32, -1, -1, Z::index_count);
    int n_det = pixel_buf_man.dims[0];
    int n_time = pixel_buf_man.dims[1];

    // Similarly, spin_proj tells you the number of components.
    auto spin_proj_man = SignalSpace<FSIGNAL>(
        spin_proj, "spin_proj", FSIGNAL_NPY_TYPE, n_det, n_time, -1);
    int n_spin = spin_proj_man.dims[2];

    auto tiling = T();
    tiling.TestInputs(map, false, true, n_spin);

    BufferWrapper<FSIGNAL> _det_weights("det_weights", det_weights, true,
                                        vector<int>{n_det});

    // Below we assume the pixel sub-indices are close-packed.
    if (pixel_buf_man.steps[1] != 1)
        throw shape_exception("pixel_index",
                              "Fast dimension of pixel indices must be close-packed.");

    for (int i_det = 0; i_det < n_det; ++i_det) {
        FSIGNAL det_wt = 1.;
        if (_det_weights->obj != NULL)
            det_wt = *(FSIGNAL*)((char*)_det_weights->buf + _det_weights->strides[0]*i_det);

        for (int i_time = 0; i_time < n_time; ++i_time) {
            const int *pixel_ofs = pixel_buf_man.data_ptr[i_det] +
                pixel_buf_man.steps[0]*i_time;
            if (pixel_ofs[0] < 0)
                continue;
            const FSIGNAL *sp = (spin_proj_man.data_ptr[i_det] +
                                 spin_proj_man.steps[0] * i_time);
            for (int i_spin = 0; i_spin < n_spin; ++i_spin) {
                for (int j_spin = i_spin; j_spin < n_spin; ++j_spin) {
                    *tiling.wpix(i_spin, j_spin, pixel_ofs) += det_wt *
                        sp[spin_proj_man.steps[1] * i_spin] *
                        sp[spin_proj_man.steps[1] * j_spin];
                }
            }
        }
    }

    return map;
}

template <typename Z, typename T>
bp::object ProjEng_Precomp<Z,T>::from_map(
    bp::object map, bp::object pixel_index, bp::object spin_proj,
    bp::object signal, bp::object det_weights)
{
    // You won't get far without pixel_index, so use that to nail down
    // the n_time and n_det.
    auto pixel_buf_man = SignalSpace<int32_t>(
        pixel_index, "pixel_index", NPY_INT32, -1, -1, Z::index_count);
    int n_det = pixel_buf_man.dims[0];
    int n_time = pixel_buf_man.dims[1];

    // Similarly, spin_proj tells you the number of components.
    auto spin_proj_man = SignalSpace<FSIGNAL>(
        spin_proj, "spin_proj", FSIGNAL_NPY_TYPE, n_det, n_time, -1);
    int n_spin = spin_proj_man.dims[2];

    auto tiling = T();
    tiling.TestInputs(map, true, false, n_spin);

    // Check the signal.
    auto signal_man = SignalSpace<FSIGNAL>(
        signal, "signal", FSIGNAL_NPY_TYPE, n_det, n_time);

    BufferWrapper<FSIGNAL> _det_weights("det_weights", det_weights, true,
                                        vector<int>{n_det});

#pragma omp parallel for
    for (int i_det = 0; i_det < n_det; ++i_det) {
        FSIGNAL det_wt = 1.;
        if (_det_weights->obj != NULL)
            det_wt = *(FSIGNAL*)((char*)_det_weights->buf + _det_weights->strides[0]*i_det);

        for (int i_time = 0; i_time < n_time; ++i_time) {
            const int *pixel_ofs = pixel_buf_man.data_ptr[i_det] +
                pixel_buf_man.steps[0]*i_time;
            if (pixel_ofs[0] < 0)
                continue;
            FSIGNAL sig = 0;
            for (int i_spin = 0; i_spin < n_spin; ++i_spin) {
                FSIGNAL sp = *(spin_proj_man.data_ptr[i_det] +
                               spin_proj_man.steps[0] * i_time +
                               spin_proj_man.steps[1] * i_spin);
                sig += sp * *tiling.pix(i_spin, pixel_ofs);
            }
            *(signal_man.data_ptr[i_det] + signal_man.steps[0]*i_time) += sig;
        }
    }

    return signal_man.ret_val;
}
*/


#define TYPEDEF_BATCH(PIX)                                              \
    typedef ProjectionEngine<Proj ## PIX, Pixelizor2_Flat<NonTiled>, NonTiled, SpinTQU> ProjEng_##PIX##_TQU_NonTiled;

TYPEDEF_BATCH(Flat)
// TYPEDEF_BATCH(CAR)
// TYPEDEF_BATCH(TAN)

#define EXPORT_ENGINE(CLASSNAME, INIT_CLASS)                            \
    bp::class_<CLASSNAME>(#CLASSNAME, bp::init<bp::object>())           \
    .add_property("index_count", &CLASSNAME::index_count)               \
    .add_property("comp_count", &CLASSNAME::comp_count)                 \
    .def("coords", &CLASSNAME::coords)                                  \
    .def("pixels", &CLASSNAME::pixels)                                  \
    .def("pointing_matrix", &CLASSNAME::pointing_matrix)                \
    .def("zeros", &CLASSNAME::zeros)                                    \
    .def("from_map", &CLASSNAME::from_map)                              \
    .def("to_map", &CLASSNAME::to_map)                                  \
    .def("to_weight_map", &CLASSNAME::to_weight_map)                    \
    ;

/*
    .def("to_map_omp", &CLASSNAME::to_map_omp)                          \
    .def("to_weight_map_omp", &CLASSNAME::to_weight_map_omp)            \
    .def("pixel_ranges", &CLASSNAME::pixel_ranges)                      \
*/

// typedef ProjEng_Precomp<Pixelizor2_Flat,NonTiled> ProjEng_Precomp_;
// typedef ProjEng_Precomp<Pixelizor2_Flat_Tiled,Tiled> ProjEng_Precomp_Tiled;

#define EXPORT_PRECOMP(CLASSNAME)                                       \
    bp::class_<CLASSNAME>(#CLASSNAME)                                   \
    .def("to_map", &CLASSNAME::to_map)                                  \
    .def("to_weight_map", &CLASSNAME::to_weight_map)                    \
    .def("from_map", &CLASSNAME::from_map);


#define EXPORT_BATCH(PIX)                                     \
    EXPORT_ENGINE(PIX ## _TQU_NonTiled,   Flat);

/*
    EXPORT_ENGINE(PIX ## _QU,  Flat);                         \
    EXPORT_ENGINE(PIX ## _TQU, Flat);                         \
    EXPORT_ENGINE(PIX ## _Tiled_T,   Flat_Tiled);             \
    EXPORT_ENGINE(PIX ## _Tiled_QU,  Flat_Tiled);             \
    EXPORT_ENGINE(PIX ## _Tiled_TQU, Flat_Tiled);
*/

// There are probably better ways to do this..
template<typename T>
inline
int _index_count(const T &) { return T::index_count; }

PYBINDINGS("so3g")
{
    EXPORT_BATCH(ProjEng_Flat);
    // EXPORT_BATCH(ProjEng_CAR);
    // EXPORT_BATCH(ProjEng_TAN);
    // EXPORT_ENGINE(ProjEng_Flat_T, Flat);
    // EXPORT_ENGINE(ProjEng_Flat_QU, Flat);
    // EXPORT_ENGINE(ProjEng_Flat_TQU, Flat);
    // EXPORT_ENGINE(ProjEng_Flat_Tiled_T,   Flat_Tiled);
    // EXPORT_ENGINE(ProjEng_Flat_Tiled_QU,  Flat_Tiled);
    // EXPORT_ENGINE(ProjEng_Flat_Tiled_TQU, Flat_Tiled);
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
    // EXPORT_PRECOMP(ProjEng_Precomp_);
    // EXPORT_PRECOMP(ProjEng_Precomp_Tiled);
}
