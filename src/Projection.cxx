#define NO_IMPORT_ARRAY

#include <iostream>
#include <assert.h>
#include <cmath>
#include <cstdio>
#include <typeinfo>

#ifdef _OPENMP
# include <omp.h>
#endif // ifdef _OPENMP

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "quaternion.h"
#include "so3g_numpy.h"
#include "Projection.h"
#include "Ranges.h"
#include "exceptions.h"
#include "so_linterp.h"

using namespace std;

namespace py = pybind11;

#include "healpix_bare.c"


// TRIG_TABLE_SIZE
//
// Set this macro to enable trig interpolation tables.  Studies show
// that a table sizeof 16k gives errors below .1 arcsecond.  To use
// math library asin and atan2, leave the macro undefined.
//
#define TRIG_TABLE_SIZE (1 << 14)

#ifdef TRIG_TABLE_SIZE
static asinTable asin_lookup(TRIG_TABLE_SIZE);
static atan2Table atan2_lookup(TRIG_TABLE_SIZE);
#  define ASIN asin_lookup.get
#  define ATAN2 atan2_lookup.get
#else
#  define ASIN asin
#  define ATAN2 atan2
#endif


typedef quaternion<double> quatd;

inline bool isNone(const py::object &pyo)
{
    return pyo.is_none();
}


// Convert an n-dimensional array into a list of array slices.
py::list list_of_arrays(py::object input) {
    if (input.is_none()) {
        return py::list();
    }
    if (py::isinstance<py::list>(input)) {
        // Already a list, pass-through
        return input;
    }

    // The returned list
    py::list output;

    auto in_arr = py::cast<py::array>(input);
    int ndim = in_arr.ndim();
    if (ndim == 1) {
        output.append(in_arr);
    } else {
        for (int i = 0; i < in_arr.shape(0); ++i) {
            auto slc = py::make_tuple(py::int_(i), py::ellipsis());
            py::array in_arr_row = py::cast<py::array>(in_arr[slc]).squeeze();
            output.append(in_arr_row);
        }
    }
    return output;
}


// ProjEng template system
//
// ProjEng classes will be templated like:
//
//     template<typename C, typename P, typename S>
//
// Parameters C, P, S are shorthands for:
//
// - C: CoordSys, determines how the boresight position and detector
//   offsets are interpreted and combined.  The most obvious cases
//   (e.g. ProjTAN, ProjCAR), labeled with FITS projection codes,
//   imply that boresight and detector offsets should be treated as
//   rotation quaternions, multiplied together, and interpreted as
//   celestial positions within a particular projection.  Other cases
//   include ProjFlat (where boresight and detector offsets live in a
//   flat 2-d space) and ProjQuat (where requests for coordinates will
//   simply return the composed quaternion, without celestial
//   interpretation).  See Pointer<C> implementation for more details.
//
// - P: PixelSys, determining how computed coordinates (in a
//   projection plane, say) should be turned into pixel indices.  This
//   is the bit that actually needs to know the footprint of a
//   specific map -- the equivalent of CRVAL, CRPIX and CDELT.
//
// - S: SpinSys, determining the number of spin-components and the
//   equations for their weights.  This provides the abstraction for
//   T-only, TQU, and QU-only mapping, by specifying how to compute
//   "pointing matrix" projection factors for each active map
//   component.  The input here is the polarized detector's position
//   angle, gamma, relative to the reference axis of the active
//   projection.  This system can support mapping of other spin fields
//   quite easily; the most obvious candidate being spin-1, to carry
//   information such as in-plane pointing displacements.

// Bilinear mapmaking (and interpolated mapmaking in general)
//
// Originally only nearest-neighbor mapmaking was supported. Here,
// a sample takes the value of the nearest pixel, regardless of where
// inside this pixel it hits. This means each sample only needs to
// care about a single pixel, which is computed by Pixelizor's GetPixel.
// This simple relationship makes it elatively easy to avoid
// clobbering when using threads in the tod2map operation, as one could
// assign responsibility for different areas of pixels to different threads,
// and translate this directly into sample ranges for those threads.
//
// In interpolated mapmaking like bilinear mapmaking, the value of a sample
// is given by some linear combination of several of the nearest pixels.
// For bilinear mapmaking, these are the four pixels above/below and left/right
// of the sample, which are each weighted based on how close each pixel is to
// the sample. The pixels and corresponding weights are calculated using
// Pixelizor's GetPixels function, which in turn uses the Pixelizor's Interpol
// template argument to determine which type of interpolation should be used
// (including the "no interpolation" nearest neighbor case).
//
// However, since each sample now hits multiple pixels, it is more difficult
// build a set of sample ranges that touch disjoint sets of pixels and so can
// be iterated over in parallel without locking in the tod2map operation.
// There will always be some samples that hit the edge of each thread's pixel
// groups, resulting in it affecting multiple threads' pixels. We handle this
// by generalizing the sample ranges from ranges = vector<vector<RangesInt32>>, which
// is [nthread][ndet][nrange] in shape, to ranges = vector<vector<vector<RangesInt32>>>
// which is [nbunch][nthread][ndet][nrange], and where we iterate over the
// outermost "bunch" level in series, and then in parallel over the next level.
// This general framework should work for all parallelization schemes, but what
// we've implemented so far is a simple case where all pixels that unambiguously
// fall inside a single threads' pixel domain are done as before, and then thew few samples
// that fall in multiple domains are done by a single thread in the end. This
// means that ranges[0] has shape [nthread][ndet][nrange] while ranges[1] has
// shape [1][ndet][nrange]. More complicated schemes could be used, but this
// seems to work well enough.

// Arbitrary detector response support
// -----------------------------------
// It's useful to support detectors with polarization response < 1. It's also
// sometimes useful with detectors with 0 total intensity response. This is
// implemented as follows.
//
// At the top level, the functions where the responsivity is relevant, like
// from_map, to_map etc. take a py::object reponse argument, representing
// the [ndet,2] detector responsivities.
//
// This is then converted into a BufferWrapper<FSIGNAL>, before each detector's
// responsivity is extracted using get_response(BufferWrapper<FSIGNAL> & buf, int i_det),
// which returns a Response { FSIGNAL T, P; } struct.
//
// Finally, this Response struct is passed to spin_proj_factors, e.g.
//
//  void spin_proj_factors<SpinTQU>(const double* coords, const Response & response, FSIGNAL *projfacs)
//  {
//       const double c = coords[2];
//       const double s = coords[3];
//       projfacs[0] = response.T;
//       projfacs[1] = response.P*(c*c - s*s);
//       projfacs[2] = response.P*(2*c*s);
//  }

// State the template<CoordSys> options
class ProjQuat;
class ProjFlat;
class ProjCEA;
class ProjCAR;
class ProjARC;
class ProjSIN;
class ProjTAN;
class ProjZEA;

// State the template<TilingSys> options
class Tiled;
class NonTiled;

// State the pointing matrix interpolation options
struct NearestNeighbor { static const int interp_count = 1; };
struct Bilinear        { static const int interp_count = 4; };

// Helper template for SpinSys...
template <int N>
class SpinClass {
public:
    static const int comp_count = N;
};

// State the template<SpinSys> options
class SpinT : public SpinClass<1> {};
class SpinQU : public SpinClass<2> {};
class SpinTQU : public SpinClass<3> {};


template <typename CoordSys>
class Pointer {
public:
    bool TestInputs(py::object &map, py::object &pbore, py::object &pdet,
                    py::object &signal, py::object &det_weights);
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
    py::object &map, py::object &pbore, py::object &pdet,
    py::object &signal, py::object &det_weights)
{
    // Boresight and Detector must present and inter-compatible.
    _pborebuf = BufferWrapper<double>("boresight", pbore, false,
                                      vector<int>{-1, 4});
    _pdetbuf = BufferWrapper<double>("detectors", pdet, false,
                                     vector<int>{-1, 4});
    n_time = _pborebuf->shape[0];
    n_det = _pdetbuf->shape[0];

    // Check the detectors for NaNs.
    for (int i=0; i<n_det; i++) {
        for (int ic = 0; ic < 4; ++ic) {
            const char *x = (char*)_pdetbuf->buf
                + _pdetbuf->strides[0] * i
                + _pdetbuf->strides[1] * ic;
            if (std::isnan(*(double*)x)) {
                std::ostringstream err;
                err << "Pointing offset error: nan found at index " << i << ".";
                throw ValueError_exception(err.str());
            }
        }
    }
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

// ProjQuat: Not a projection -- returns the quaternion rotation
// components.

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

// ProjFlat: Not a spherical projection -- assumes flat space (as in
// FITS X,Y type coordinates).

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

// ProjARC: the zenithal equidistant projection.
//
// The first two coordinates are R(lat)*sin(lon) and -R(lat)*cos(lon)
// [see FITS-II].  Then cos and sin of parallactic angle.  For ARC,
// R(lat) = 90 - lat = theta.

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
    const double sin_theta2_sq = b*b + c*c;

    const double sc = c*a + d*b;
    const double ss = a*b - c*d;
    const double half_sin_theta = sqrt(sc*sc + ss*ss);

    double R_factor;
    if (half_sin_theta < 1e-8)
        R_factor = 2 + 1.33333333333*half_sin_theta*half_sin_theta;
    else
        R_factor = ASIN(half_sin_theta*2) / half_sin_theta;

    // Allow for theta in [pi/2, pi).
    if (cos_theta2_sq < sin_theta2_sq)
        R_factor = M_PI / half_sin_theta - R_factor;

    coords[0] = ss * R_factor;
    coords[1] = sc * R_factor;
    coords[2] = (a*a - d*d) / cos_theta2_sq;
    coords[3] = (2*a*d) / cos_theta2_sq;
}

// ProjSIN: the orthographic projection
//
// This is the projection plane most naturally associated with our
// "xieta" system -- pixels are on lines of constant (xi, eta).
//
// The first two coordinates are R(lat)*sin(lon) and -R(lat)*cos(lon)
// [see FITS-II].  Then cos and sin of parallactic angle.  For SIN,
// R(lat) = sin(90 - lat) = sin(theta).

template <>
inline
void Pointer<ProjSIN>::GetCoords(int i_det, int i_time,
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

    coords[0] = 2 * (a*b - c*d);
    coords[1] = 2 * (c*a + d*b);
    coords[2] = (a*a - d*d) / cos_theta2_sq;
    coords[3] = (2*a*d) / cos_theta2_sq;
}

// ProjTAN: the tangent plane (gnomonic) projection.  It is zenithal.
//
// The first two coordinates are R(lat)*sin(lon) and -R(lat)*cos(lon)
// [see FITS-II].  Then cos and sin of parallactic angle.  For TAN,
// R(lat) = tan(lat) = cot(theta).

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

// ProjZEA: the zenithal equal area projection.
//
// The first two coordinates are R(lat)*sin(lon) and -R(lat)*cos(lon)
// [see FITS-II].  Then cos and sin of parallactic angle.  For ZEA,
// R(lat) = sqrt(2(1-cos(lat))) = sqrt(2(1-sin(theta))).

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

// ProjCEA: Cylindrical projection.
//
// First two coordinates are lon (in radians) and sin(lat).  Then cos
// and sin of parallactic angle.

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

// ProjCAR: Cylindrical projection.
//
// First two coordinates are lon and lat (in radians).  Then cos and
// sin of parallactic angle.

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

    coords[0] = ATAN2(c*d - a*b, c*a + d*b);
    coords[1] = ASIN(cos_theta);   // Yes, cos(theta) = sin(lat).
    coords[2] = (a*c - b*d) / half_sin_theta;
    coords[3] = (c*d + a*b) / half_sin_theta;
}


//! Pixelizors
//

const int NSIDE_MAX=8192; // int32 should work up to nside8192
template <typename TilingSys>
class Pixelizor_Healpix;

// Pixelizor_Healpix<NonTiled>
//
// Make [npix] full-sky maps in HEALPIX NEST or RING ordering
// int32 for pixel indexes will work up to NSIDE=8192. Assuming int is int32

template<>
class Pixelizor_Healpix<NonTiled> {
public:
    static const int index_count = 1;
    static const int interp_count = 1;
    Pixelizor_Healpix(int nside, bool nest=true){
        this->nside = nside;
        this->nest = nest;
        npix = nside2npix(nside);
        check_nside(nside);
    };
    Pixelizor_Healpix() {};
    Pixelizor_Healpix(py::object args) {
        py::tuple args_tuple = py::cast<py::tuple>(args);
        // args[0]: int nside
        nside = py::cast<int>(args_tuple[0]);
        // args[1]: bool isnest
        nest = py::cast<bool>(args_tuple[1]);
        npix = nside2npix(nside);
        check_nside(nside);
    }
    ~Pixelizor_Healpix() {};

    py::object zeros(vector<int> shape) {
        std::cerr << "PixHP non-tile input shape has len " << shape.size() << std::endl;
        shape.push_back(npix);
        int ndim = 0;
        npy_intp dims[32];
        for (auto d: shape) {
            dims[ndim++] = d;
        }
        std::cerr << "PixHP non-tile ndim = " << ndim << " (";
        for (int i = 0; i < ndim; ++i) {
            std::cerr << dims[i] << ", ";
        }
        std::cerr << ")" << std::endl;

        int dtype = NPY_FLOAT64;
        PyObject *v = PyArray_ZEROS(ndim, dims, dtype, 0);
        return py::reinterpret_steal<py::object>(v);
    }

    inline
    void GetPixel(int i_det, int i_time, const double *coords, int *pixel_index) {
        double phi = coords[0]; // lon, -2pi->2pi
        double theta = M_PI/2 - coords[1]; // colat, 0->pi
        t_ang ang = (t_ang) {theta, phi};
        int ix = (nest) ? ang2nest(nside, ang) : ang2ring(nside, ang);
        pixel_index[0] = ix;
    }
    inline
    int GetPixels(int i_det, int i_time, const double *coords, int pixinds[interp_count][index_count], FSIGNAL pixweights[interp_count]){
        double phi = coords[0]; // lon, -2pi->2pi
        double theta = M_PI/2 - coords[1]; // colat, 0->pi
        t_ang ang = (t_ang) {theta, phi};
        int ix = (nest) ? ang2nest(nside, ang) : ang2ring(nside, ang);
        pixinds[0][0] = ix;
        pixweights[0] = 1;
        return 1;
    }

    bool TestInputs(py::object &map, bool need_map, bool need_weight_map, int comp_count) {
        if (need_map) {
            // The map is mandatory, and the leading axis must match the
            // component count.  Then the one healpix pixel axis.
            mapbuf = BufferWrapper<double>("map", map, false,
                                           vector<int>{comp_count,npix});
        } else if (need_weight_map) {
            // The map is mandatory, and the two leading axes must match
            // the component count.  Then the one healpix pixel axis.
            mapbuf = BufferWrapper<double>("map", map, false,
                                           vector<int>{comp_count,comp_count,npix});
        }
        return true;
    }
    double *pix(int imap, const int pixel_index[]) {
        return (double*)((char*)mapbuf->buf +
                         mapbuf->strides[0]*imap +
                         mapbuf->strides[1]*pixel_index[0]);
    }
    double *wpix(int imap, int jmap, const int pixel_index[]) {
        return (double*)((char*)mapbuf->buf +
                         mapbuf->strides[0]*imap +
                         mapbuf->strides[1]*jmap +
                         mapbuf->strides[2]*pixel_index[0]);
    }
    int stripe(const int pixel_index[], int thread_count) {
        return pixel_index[0] * thread_count / npix;
    }
    int tile_count() {
        return -1;
    }
    void check_nside(int nside){
        bool isok = (nside > 0) & (nside < NSIDE_MAX) & (nside & (nside-1)) == 0; // check nside is a power of 2
        if (! isok){
          std::ostringstream err;
          err << "Invalid nside " << nside;
          throw ValueError_exception(err.str());
        }
    }

    int nside;
    int npix;
    BufferWrapper<double> mapbuf;
    bool nest;
};

// Pixelizor_Healpix<Tiled>
//
// Make a len(nTile) list of [npix/nTile] partial-sky maps in HEALPIX NEST ordering
// inactive tiles will contain None
// int32 for pixel indexes will work up to NSIDE=8192. Assuming int is int32

template<>
class Pixelizor_Healpix<Tiled> {
public:
    static const int index_count = 2;
    static const int interp_count = 1;
    Pixelizor_Healpix(int nside){
        this->nside = nside;
        this->nest = true; // we only tile maps in nest
        populate = vector<bool>(1, true);
        check_nside(nside);
    };
    Pixelizor_Healpix() {};
    Pixelizor_Healpix(py::object args) {
        // args[0]: int nside
        py::tuple args_tuple = py::cast<py::tuple>(args);
        nside = py::cast<int>(args_tuple[0]);
        // args[1]: bool nestin: MUST BE true ; check that
        bool nestin = py::cast<bool>(args_tuple[1]); // "nest" argument, not used for tiled maps
        if (! nestin){
            std::ostringstream err;
            err << "RING not supported for tiled maps";
            throw ValueError_exception(err.str());
        }
        // args[2] int nside_tile; nside defining the tiling
        int nside_tile = py::cast<int>(args_tuple[2]);
        ntiles = nside2npix(nside_tile);
        if (py::len(args) >= 4) {
            // args[3] list(int) of indexes for active tiles
            py::list active_tiles = list_of_arrays(args_tuple[3]);
            if (! isNone(active_tiles)) {
                populate = vector<bool>(ntiles, false);
                for (int i=0; i < py::len(active_tiles); i++) {
                    int idx = py::cast<int>(active_tiles[i]);
                    if (idx >= 0 && idx < ntiles) {
                        populate[idx] = true;
                    }
                }
            }
        }
        int npix = nside2npix(nside);
        npix_per_tile = npix / ntiles;
        check_nside(nside); // check validity
        check_nside(nside_tile);
        if (nside_tile > nside) {
            std::ostringstream err;
            err << "Invalid nside_tile " << nside_tile << " > nside " << nside;
            throw ValueError_exception(err.str());
        }
    }
  ~Pixelizor_Healpix() {};

    py::object zeros(vector<int> shape) {
        std::cerr << "PixHP tiled input shape has len " << shape.size() << std::endl;
        int dtype = NPY_FLOAT64;
        int ndim = 0;
        npy_intp dims[32];
        for (auto d: shape) {
            dims[ndim++] = d;
        }
        ndim += 1;

        if (populate.size() == 0)
            throw RuntimeError_exception("Cannot create blank tiled map unless "
                                  "user has specified what tiles to populate.");

        py::list maps_out;
        auto pop_iter = populate.begin();
        for (int itile = 0; itile < ntiles; itile++){
            std::cerr << "PixHP zeros, tile " << itile;
            bool pop_this = (pop_iter != populate.end()) && *(pop_iter++);
            if (pop_this) {
                dims[ndim-1] = npix_per_tile;
                std::cerr << "PixHP zeros, pop_this = true, append ndim = " << ndim << " (";
                for (int i = 0; i < ndim; ++i) {
                    std::cerr << dims[i] << ",";
                }
                std::cerr << ")" << std::endl;
                PyObject *v = PyArray_ZEROS(ndim, dims, dtype, 0);
                maps_out.append(py::reinterpret_steal<py::object>(v));
            } else {
                std::cerr << "PixHP zeros, pop_this = false, append None" << std::endl;
                maps_out.append(py::none());
            }
        }

        std::cerr << "PixHP tiled ndim = " << ndim << " (";
        for (int i = 0; i < ndim; ++i) {
            std::cerr << dims[i] << ", ";
        }
        std::cerr << ")" << std::endl;
        std::cerr << "PixHP final zeros has " << py::len(maps_out) << " tiles" << std:: endl;
        return maps_out;
    }

    inline
    void GetPixel(int i_det, int i_time, const double *coords, int *pixel_index) {
        double phi = coords[0]; // lon, -2pi->2pi
        double theta = M_PI/2 - coords[1]; // colat, 0->pi
        t_ang ang = (t_ang) {theta, phi};
        int ix = ang2nest(nside, ang);
        //int ix = (nest) ? ang2nest(nside, ang) : ang2ring(nside, ang);
        pixel_index[0] = ix / npix_per_tile;
        pixel_index[1] = ix % npix_per_tile;
    }
    inline
    int GetPixels(int i_det, int i_time, const double *coords, int pixinds[interp_count][index_count], FSIGNAL pixweights[interp_count]){
        double phi = coords[0]; // lon, -2pi->2pi
        double theta = M_PI/2 - coords[1]; // colat, 0->pi
        t_ang ang = (t_ang) {theta, phi};
        int ix = ang2nest(nside, ang);
        //int ix = (nest) ? ang2nest(nside, ang) : ang2ring(nside, ang);
        pixinds[0][0] = ix / npix_per_tile;
        pixinds[0][1] = ix % npix_per_tile;
        pixweights[0] = 1;
        return 1;
    }

    bool TestInputs(py::object & map_obj, bool need_map, bool need_weight_map, int comp_count) {
        // Get list of input maps
        py::list map = list_of_arrays(map_obj);

        vector<int> map_shape_req;
        if (need_map) {
            // The map is mandatory, and the leading axis must match the
            // component count, second axis the npix_per_tile.
            map_shape_req = {comp_count,npix_per_tile};
        } else if (need_weight_map) {
            // The map is mandatory, and the two leading axes must match
            // the component count.
            map_shape_req = {comp_count,comp_count,npix_per_tile};
        }
        if (map_shape_req.size() == 0)
            return true;

        mapbufs.clear();
        // Number of tiles
        int n_tile = py::len(map);

        std::cerr << "healpix tiled TestInputs got map with n_tile = " << n_tile << std::endl;
        for (int i_tile = 0; i_tile < n_tile; i_tile++) {
            if (! isNone(map[i_tile])) {
                py::array mp = py::cast<py::array>(map[i_tile]);
                std::cerr << "healpix tiled TestInputs tile " << i_tile << " map shape: ";
                for (int im = 0; im < mp.ndim(); ++im) {
                    std::cerr << mp.shape(im) << ",";
                }
                std::cerr << std::endl;
            }
        }

        for (int i_tile = 0; i_tile < n_tile; i_tile++) {
            if (isNone(map[i_tile])) {
                if (populate[i_tile])
                    throw tiling_exception(i_tile, "Projector expects tile but it is missing.");
                mapbufs.push_back(BufferWrapper<double>());
            } else {
                mapbufs.push_back(
                    BufferWrapper<double>("map", map[i_tile], false, map_shape_req)
                );
            }
        }

        return true;
    }

    double *pix(int imap, const int pixel_index[]) {
        const BufferWrapper<double> &mapbuf = mapbufs[pixel_index[0]];
        if (mapbuf->buf == nullptr)
            throw tiling_exception(pixel_index[0],
                                   "Attempted pointing operation on non-instantiated tile.");
        return (double*)((char*)mapbuf->buf +
                         mapbuf->strides[0]*imap +
                         mapbuf->strides[1]*pixel_index[1]);
    }
    double *wpix(int imap, int jmap, const int pixel_index[]) {
        // Expensive shared_ptr copy?
        const BufferWrapper<double> &mapbuf = mapbufs[pixel_index[0]];
        if (mapbuf->buf == nullptr)
            throw tiling_exception(pixel_index[0],
                                   "Attempted pointing operation on non-instantiated tile.");
        return (double*)((char*)mapbuf->buf +
                         mapbuf->strides[0]*imap +
                         mapbuf->strides[1]*jmap +
                         mapbuf->strides[2]*pixel_index[1]);
    }

    int stripe(const int pixel_index[], int thread_count) {
        int ipix = pixel_index[0] * npix_per_tile + pixel_index[1];
        int npix = npix_per_tile * ntiles;
        return ipix * thread_count / npix;
    }
    int tile_count() {
        return ntiles;
    }

    void check_nside(int nside){
        bool isok = (nside > 0) & (nside < NSIDE_MAX) & (nside & (nside-1)) == 0; // check power of 2
        if (! isok){
            std::ostringstream err;
            err << "Invalid nside " << nside;
            throw ValueError_exception(err.str());
        }
    }

    int nside;
    int ntiles;
    int npix_per_tile;
    vector<bool> populate;
    vector<BufferWrapper<double>> mapbufs;
    bool nest;
};

template <typename TilingSys, typename Interpol = NearestNeighbor>
class Pixelizor2_Flat;

template <typename Interpol>
class Pixelizor2_Flat<NonTiled, Interpol> {
public:
    static const int index_count = 2;
    static const int interp_count = Interpol::interp_count;
    Pixelizor2_Flat(int ny, int nx,
                    double dy=1., double dx=1.,
                    double iy0=0., double ix0=0.) {
        naxis[0] = ny;
        naxis[1] = nx;
        cdelt[0] = dy;
        cdelt[1] = dx;
        crpix[0] = iy0;
        crpix[1] = ix0;
        std::cerr << "pixelizor2 non-tile ctor naxis = " << naxis[0] << "," << naxis[1] << std::endl;
    };
    Pixelizor2_Flat() : naxis{1,1} {};
    Pixelizor2_Flat(py::object args) {
        py::tuple args_tuple = py::cast<py::tuple>(args);
        naxis[0] = py::cast<int>(args_tuple[0]);
        naxis[1] = py::cast<int>(args_tuple[1]);
        cdelt[0] = py::cast<double>(args_tuple[2]);
        cdelt[1] = py::cast<double>(args_tuple[3]);
        crpix[0] = py::cast<double>(args_tuple[4]);
        crpix[1] = py::cast<double>(args_tuple[5]);
        std::cerr << "pixelizor2 non-tile pyctor naxis = " << naxis[0] << "," << naxis[1] << std::endl;
    }
    ~Pixelizor2_Flat() {};

    py::object zeros(vector<int> shape) {
        shape.push_back(naxis[0]);
        shape.push_back(naxis[1]);
        int ndim = 0;
        npy_intp dims[32];
        for (auto d: shape)
            dims[ndim++] = d;

        int dtype = NPY_FLOAT64;
        PyObject *v = PyArray_ZEROS(ndim, dims, dtype, 0);
        return py::reinterpret_steal<py::object>(v);
    }

    inline
    void GetPixel(int i_det, int i_time, const double *coords, int *pixel_index) {
        //double ix = coords[0] / cdelt[1] + crpix[1] + 0.5;
        double ix = coords[0] / cdelt[1] + crpix[1] - 1 + 0.5;
        if (ix < 0 || ix >= naxis[1]) {
            *pixel_index = -1;
            return;
        }

        //double iy = coords[1] / cdelt[0] + crpix[0] + 0.5;
        double iy = coords[1] / cdelt[0] + crpix[0] - 1 + 0.5;
        if (iy < 0 || iy >= naxis[0]) {
            *pixel_index = -1;
            return;
        }

        pixel_index[0] = int(iy);
        pixel_index[1] = int(ix);
    }
    inline
    int GetPixels(int i_det, int i_time, const double *coords, int pixinds[interp_count][index_count], FSIGNAL pixweights[interp_count]);

    bool TestInputs(py::object &map, bool need_map, bool need_weight_map, int comp_count) {
        py::array mp = py::cast<py::array>(map);
        std::cerr << "pixelizer2 non-tiled TestInputs map shape: ";
        for (int im = 0; im < mp.ndim(); ++im) {
            std::cerr << mp.shape(im) << ",";
        }
        std::cerr << std::endl;

        if (need_map) {
            // The map is mandatory, and the leading axis must match the
            // component count.  And then 2 celestial axes.
            mapbuf = BufferWrapper<double>("map", map, false,
                                           vector<int>{comp_count,-1,-1});
        } else if (need_weight_map) {
            // The map is mandatory, and the two leading axes must match
            // the component count.  And then 2 celestial axes.
            mapbuf = BufferWrapper<double>("map", map, false,
                                           vector<int>{comp_count,comp_count,-1,-1});
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
        return pixel_index[1] * thread_count / naxis[1];
    }
    int tile_count() {
        return -1;
    }

    double crpix[2];
    double cdelt[2];
    int naxis[2];
    BufferWrapper<double> mapbuf;
};

// Take care of the parts that depend on the interpolation order.
// First NearestNeighbor interpolation

template<>
inline int Pixelizor2_Flat<NonTiled, NearestNeighbor>::GetPixels(int i_det, int i_time, const double *coords, int pixinds[interp_count][index_count], FSIGNAL pixweights[interp_count]) {
    // For nearest neighbor this basically reduces to GetPixel
    double x = coords[0] / cdelt[1] + crpix[1] - 1 + 0.5;
    double y = coords[1] / cdelt[0] + crpix[0] - 1 + 0.5;
    if (x < 0 || x >= naxis[1]) return 0;
    if (y < 0 || y >= naxis[0]) return 0;
    // SKN: On my laptop int(y)-int(y<0) takes 0.63 ns vs. 0.84 ns for int(floor(y)).
    // Both are very fast, so maybe it's better to go for the latter, since it's clearer.
    // For comparison the plain cast was 0.45 ns.
    pixinds[0][0] = int(y)-int(y<0);
    pixinds[0][1] = int(x)-int(x<0);
    pixweights[0] = 1;
    return 1;
}

// Then Bilinear interpolation

template<>
inline int Pixelizor2_Flat<NonTiled, Bilinear>::GetPixels(int i_det, int i_time, const double *coords, int pixinds[interp_count][index_count], FSIGNAL pixweights[interp_count]) {
    // For bilinear mapmaking we need to visit the four bounding pixels
    double x  = coords[0] / cdelt[1] + crpix[1] - 1 + 0.5;
    double y  = coords[1] / cdelt[0] + crpix[0] - 1 + 0.5;
    int    x1 = int(x)-int(x<0);
    int    y1 = int(y)-int(y<0);
    double wx[2] = {x-x1, 1-(x-x1)};
    double wy[2] = {y-y1, 1-(y-y1)};
    // Loop through the our cases
    int iout = 0;
    for(int iy = y1; iy < y1+2; iy++) {
        if(iy < 0 || iy >= naxis[0]) continue;
        for(int ix = x1; ix < x1+2; ix++) {
            if(ix < 0 || ix >= naxis[1]) continue;
            pixinds[iout][0] = iy;
            pixinds[iout][1] = ix;
            pixweights[iout] = wy[iy-y1]*wx[ix-x1];
            iout++;
        }
    }
    return iout;
}

template <typename Interpol>
class Pixelizor2_Flat<Tiled, Interpol> {
public:
    static const int index_count = 3;
    static const int interp_count = Interpol::interp_count;
    Pixelizor2_Flat() {};
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
    Pixelizor2_Flat(py::object args) {
        parent_pix = Pixelizor2_Flat<NonTiled>(args); // first 6...
        py::tuple args_tuple = py::cast<py::tuple>(args);

        tile_shape[0] = py::cast<int>(args_tuple[6]);
        tile_shape[1] = py::cast<int>(args_tuple[7]);
        // Check for tile list as arg[8].
        if (py::len(args) >= 9) {
            int n_ty = (parent_pix.naxis[0] + tile_shape[0] - 1) / tile_shape[0];
            int n_tx = (parent_pix.naxis[1] + tile_shape[1] - 1) / tile_shape[1];
            populate = vector<bool>(n_ty * n_tx, false);
            py::list active_tiles = list_of_arrays(args_tuple[8]);
            if (! isNone(active_tiles)) {
                for (int i=0; i < py::len(active_tiles); i++) {
                    int idx = py::cast<int>(active_tiles[i]);
                    if (idx >= 0 && idx < n_tx*n_ty) {
                        populate[idx] = true;
                    }
                }
            }
        }
    }
    ~Pixelizor2_Flat() {};

    py::object zeros(vector<int> shape) {
        int dtype = NPY_FLOAT64;
        int ndim = 0;
        npy_intp dims[32];
        for (auto d: shape)
            dims[ndim++] = d;
        ndim += 2;

        int n_ty = (parent_pix.naxis[0] + tile_shape[0] - 1) / tile_shape[0];
        int n_tx = (parent_pix.naxis[1] + tile_shape[1] - 1) / tile_shape[1];

        if (populate.size() == 0)
            throw shape_exception("zeros", "Cannot create blank tiled map unless "
                                  "user has specified what tiles to populate.");

        py::list maps_out;
        auto pop_iter = populate.begin();
        for (int i_ty = 0; i_ty < n_ty; i_ty++) {
            for (int i_tx = 0; i_tx < n_tx; i_tx++) {
                bool pop_this = (pop_iter != populate.end()) && *(pop_iter++);
                if (pop_this) {
                    dims[ndim-2] = min(tile_shape[0], parent_pix.naxis[0] - i_ty * tile_shape[0]);
                    dims[ndim-1] = min(tile_shape[1], parent_pix.naxis[1] - i_tx * tile_shape[1]);
                    PyObject *v = PyArray_ZEROS(ndim, dims, dtype, 0);
                    maps_out.append(py::reinterpret_steal<py::object>(v));
                } else
                    maps_out.append(py::none());
            }
        }
        return maps_out;
    }

    inline
    void GetPixel(int i_det, int i_time, const double *coords, int *pixel_index) {
        //double ix = coords[0] / parent_pix.cdelt[1] + parent_pix.crpix[1] + 0.5;
        double ix = coords[0] / parent_pix.cdelt[1] + parent_pix.crpix[1] - 1 + 0.5;
        if (ix < 0 || ix >= parent_pix.naxis[1]) {
            pixel_index[0] = -1;
            return;
        }

        //double iy = coords[1] / parent_pix.cdelt[0] + parent_pix.crpix[0] + 0.5;
        double iy = coords[1] / parent_pix.cdelt[0] + parent_pix.crpix[0] - 1 + 0.5;
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
    inline
    int GetPixels(int i_det, int i_time, const double *coords, int pixinds[interp_count][index_count], FSIGNAL pixweights[interp_count]);

    bool TestInputs(py::object & map_obj, bool need_map, bool need_weight_map, int comp_count) {
        // List of input maps
        py::list map = list_of_arrays(map_obj);

        vector<int> map_shape_req;
        if (need_map) {
            // The map is mandatory, and the leading axis must match the
            // component count.  And then 2 celestial axes.
            map_shape_req = {comp_count,-1,-1};
        } else if (need_weight_map) {
            // The map is mandatory, and the two leading axes must match
            // the component count.  And then 2 celestial axes.
            map_shape_req = {comp_count,comp_count,-1,-1};
        }
        if (map_shape_req.size() == 0)
            return true;

        mapbufs.clear();
        // Number of tiles
        int n_tile = py::len(map);
        std::cerr << "pixelizer2 tiled TestInputs got map with n_tile = " << n_tile << std::endl;
        for (int i_tile = 0; i_tile < n_tile; i_tile++) {
            if (! isNone(map[i_tile])) {
                py::array mp = py::cast<py::array>(map[i_tile]);
                std::cerr << "pixelizer2 tiled TestInputs tile " << i_tile << " map shape: ";
                for (int im = 0; im < mp.ndim(); ++im) {
                    std::cerr << mp.shape(im) << ",";
                }
                std::cerr << std::endl;
            }
        }

        std::cerr << "pixelizer2 tiled TestInputs: map_shape_req = ";
        for (int im = 0; im < map_shape_req.size(); ++im) {
            std::cerr << map_shape_req[im] << ",";
        }
        std::cerr << std::endl;

        for (int i_tile = 0; i_tile < n_tile; i_tile++) {
            if (isNone(map[i_tile])) {
                if (populate[i_tile])
                    throw tiling_exception(i_tile, "Projector expects tile but it is missing.");
                mapbufs.push_back(BufferWrapper<double>());
            } else {
                mapbufs.push_back(
                    BufferWrapper<double>("map", map[i_tile], false, map_shape_req)
                );
            }
        }

        return true;
    }
    double *pix(int imap, const int pixel_index[]) {
        const BufferWrapper<double> &mapbuf = mapbufs[pixel_index[0]];
        if (mapbuf->buf == nullptr)
            throw tiling_exception(pixel_index[0],
                                   "Attempted pointing operation on non-instantiated tile.");
        return (double*)((char*)mapbuf->buf +
                         mapbuf->strides[0]*imap +
                         mapbuf->strides[1]*pixel_index[1] +
                         mapbuf->strides[2]*pixel_index[2]);
    }
    double *wpix(int imap, int jmap, const int pixel_index[]) {
        // Expensive shared_ptr copy?
        const BufferWrapper<double> &mapbuf = mapbufs[pixel_index[0]];
        if (mapbuf->buf == nullptr)
            throw tiling_exception(pixel_index[0],
                                   "Attempted pointing operation on non-instantiated tile.");
        return (double*)((char*)mapbuf->buf +
                         mapbuf->strides[0]*imap +
                         mapbuf->strides[1]*jmap +
                         mapbuf->strides[2]*pixel_index[1] +
                         mapbuf->strides[3]*pixel_index[2]);
    }
    int stripe(const int pixel_index[], int thread_count) {
        return pixel_index[0] % thread_count;
    }
    int tile_count() {
        int n_ty = (parent_pix.naxis[0] + tile_shape[0] - 1) / tile_shape[0];
        int n_tx = (parent_pix.naxis[1] + tile_shape[1] - 1) / tile_shape[1];
        return n_tx * n_ty;
    }

    Pixelizor2_Flat<NonTiled> parent_pix;
    vector<bool> populate;
    int tile_shape[2];

    vector<BufferWrapper<double>> mapbufs;
};

// Take care of the parts that depend on the interpolation order.
// First NearestNeighbor interpolation
template<>
inline int Pixelizor2_Flat<Tiled, NearestNeighbor>::GetPixels(int i_det, int i_time, const double *coords, int pixinds[interp_count][index_count], FSIGNAL pixweights[interp_count]) {
    int ix = int(coords[0] / parent_pix.cdelt[1] + parent_pix.crpix[1] - 1 + 0.5);
    if (ix < 0 || ix >= parent_pix.naxis[1]) return 0;
    int iy = int(coords[1] / parent_pix.cdelt[0] + parent_pix.crpix[0] - 1 + 0.5);
    if (iy < 0 || iy >= parent_pix.naxis[0]) return 0;
    // Get the tile and tile offsets
    int sub_y = iy / tile_shape[0];
    int sub_x = ix / tile_shape[1];
    pixinds[0][0] = sub_y * ((parent_pix.naxis[1] + tile_shape[1] - 1) / tile_shape[1]) + sub_x;
    pixinds[0][1] = iy - sub_y * tile_shape[0];
    pixinds[0][2] = ix - sub_x * tile_shape[1];
    pixweights[0] = 1;
    return 1;
}

// Then Bilinear interpolation
template<>
inline int Pixelizor2_Flat<Tiled, Bilinear>::GetPixels(int i_det, int i_time, const double *coords, int pixinds[interp_count][index_count], FSIGNAL pixweights[interp_count]) {
    // For bilinear mapmaking we need to visit the four bounding pixels
    double x  = coords[0] / parent_pix.cdelt[1] + parent_pix.crpix[1] - 1 + 0.5;
    double y  = coords[1] / parent_pix.cdelt[0] + parent_pix.crpix[0] - 1 + 0.5;
    int    x1 = int(x);
    int    y1 = int(y);
    double wx[2] = {x-x1, 1-(x-x1)};
    double wy[2] = {y-y1, 1-(y-y1)};
    // Loop through the our cases
    int iout = 0;
    for(int iy = y1; iy < y1+2; iy++) {
        if(iy < 0 || iy >= parent_pix.naxis[0]) continue;
        int sub_y = iy / tile_shape[0];
        for(int ix = x1; ix < x1+2; ix++) {
            if(ix < 0 || ix >= parent_pix.naxis[1]) continue;
            int sub_x = ix / tile_shape[1];
            pixinds[iout][0] = sub_y * ((parent_pix.naxis[1] + tile_shape[1] - 1) / tile_shape[1]) + sub_x;
            pixinds[iout][1] = iy - sub_y * tile_shape[0];
            pixinds[iout][2] = ix - sub_x * tile_shape[1];
            pixweights[iout] = wy[iy-y1]*wx[ix-x1];
            iout++;
        }
    }
    return iout;
}

// This struct is used for representing the total intensity and polarization
// response of a single detector. Originally this was just a length-2 array,
// but those can't easily be returned from functions. If I could ensure
// the buffer had stride-1 in the last dimension, then one could simply pass
// in the address to the first element for a detector, BufferWrapper doesn't make
// that straightforward, so I think this approach is best.
struct Response { FSIGNAL T, P; };
Response get_response(BufferWrapper<FSIGNAL> & buf, int i_det) {
    Response r = {*buf.ptr_2d(i_det,0), *buf.ptr_2d(i_det,1)};
    return r;
}

template <typename SpinSys>
static inline
void spin_proj_factors(const double* coords, const Response & response, FSIGNAL *projfacs);

template <>
inline void spin_proj_factors<SpinT>(const double* coords, const Response & response, FSIGNAL *projfacs)
{
    projfacs[0] = response.T;
}

template <>
inline
void spin_proj_factors<SpinQU>(const double* coords, const Response & response, FSIGNAL *projfacs)
{
    const double c = coords[2];
    const double s = coords[3];
    projfacs[0] = response.P*(c*c - s*s);
    projfacs[1] = response.P*(2*c*s);
}

template <>
inline
void spin_proj_factors<SpinTQU>(const double* coords, const Response & response, FSIGNAL *projfacs)
{
     const double c = coords[2];
     const double s = coords[3];
     projfacs[0] = response.T;
     projfacs[1] = response.P*(c*c - s*s);
     projfacs[2] = response.P*(2*c*s);
}


template <typename DTYPE>
constexpr bool SignalSpace<DTYPE>::_Validate(py::object input, std::string var_name,
                                   int dtype)
{
    // We want a list of arrays here.
    py::list sig_list;
    if (isNone(input)) {
        std::cerr << "SignalSpace.Validate(None," << var_name << "," << dtype << ")" << std::endl;
        npy_intp _dims[dims.size()];
        for (int d=0; d<dims.size(); ++d) {
            if (dims[d] < 0)
                throw shape_exception(var_name, "Cannot create space with wildcard dimensons.");
            _dims[d] = dims[d];
        }
        std::ostringstream dbg;
        for (int i = 1; i < dims.size(); ++i) {
            dbg << _dims[i] << ",";
        }
        for (int i=0; i<dims[0]; ++i) {
            PyObject *v = PyArray_ZEROS(dims.size()-1, _dims+1, dtype, 0);
            std::cerr << "SignalSpace.Validate append row " << i << " zeros(" << dims.size()-1 << ",(" << dbg.str() << "))" << std::endl;
            sig_list.append(py::reinterpret_steal<py::object>(v));
        }
    } else {
        sig_list = list_of_arrays(input);
    }
    ret_val = sig_list;

    if (dims[0] == -1) {
        dims[0] = py::len(sig_list);
        if (dims[0] == 0) {
            throw shape_exception(var_name, "has not been tested on shape 0 objects");
        }
    } else if (py::len(sig_list) != dims[0]) {
        throw shape_exception(var_name, "must contain (n_det) vectors");
    }

    const int n_det = dims[0];

    // Now extract each list member into our vector of BufferWrappers..
    // data_ptr is freed in the destructor.
    data_ptr = (DTYPE**)calloc(n_det, sizeof(*data_ptr));

    std::cerr << "SignalSpace.Validate: reserve bw size " << n_det << std::endl;
    bw.reserve(n_det);

    // Copy dims[1:] into sub_dims; potentially update sub_dims during
    // extraction, then copy those back into dims.
    vector<int> sub_dims(dims.begin()+1, dims.end());

    for (int i=0; i<n_det; i++) {
        py::object item = py::cast<py::object>(sig_list[i]);
        std::cerr << "SignalSpace.Validate: det " << i << " push back " << var_name << ", " << item << ", false, ";
        for (int d=0; d<sub_dims.size(); d++) {
            std::cerr << sub_dims[d] << ":";
        }
        std::cerr << std::endl;
        bw.push_back(BufferWrapper<DTYPE>(var_name, item, false, sub_dims));
        std::cerr << "SignalSpace.Validate: det " << i << " buffer wrapper now has " << bw.size() << " elements" << std::endl;
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
        if (bw[0]->strides[d] % bw[0]->itemsize != 0) {
            throw shape_exception(var_name, "stride is non-integral; realign.");
        }
        steps[d] = bw[0]->strides[d] / bw[0]->itemsize;
    }
    return true;
}

template <typename DTYPE>
SignalSpace<DTYPE>::SignalSpace(
    py::object input, std::string var_name, int dtype, int n_det, int n_time)
{
    dims = {n_det, n_time};
    _Validate(input, var_name, dtype);
}

template <typename DTYPE>
SignalSpace<DTYPE>::SignalSpace(
    py::object input, std::string var_name, int dtype, int n_det, int n_time,
    int n_thirdaxis)
{
    dims = {n_det, n_time, n_thirdaxis};
    _Validate(input, var_name, dtype);
}

template<typename C, typename P, typename S>
ProjectionEngine<C,P,S>::ProjectionEngine(py::object pix_args)
{
    _pixelizor = P(pix_args);
}

template<typename C, typename P, typename S>
int ProjectionEngine<C,P,S>::index_count() const {
    return P::index_count;
}

template<typename C, typename P, typename S>
int ProjectionEngine<C,P,S>::comp_count() const {
    return S::comp_count;
}

template<typename C, typename P, typename S>
py::object ProjectionEngine<C,P,S>::coords(
    py::object pbore, py::object pofs, py::object coord)
{
    auto _none = py::object();
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

template<typename C, typename P, typename S>
py::object ProjectionEngine<C,P,S>::pixels(
    py::object pbore, py::object pofs, py::object pixel)
{
    auto _none = py::object();

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
        int pixel_offset[P::index_count] = {-1};
        for (int i_time = 0; i_time < n_time; ++i_time) {
            double coords[4];
            pointer.GetCoords(i_det, i_time, (double*)dofs, (double*)coords);
            _pixelizor.GetPixel(i_det, i_time, (double*)coords, pixel_offset);
            for (int i_dim = 0; i_dim < P::index_count; i_dim++)
                pix_buf[i_time * pixel_buf_man.steps[0] +
                        i_dim * pixel_buf_man.steps[1]] = pixel_offset[i_dim];
            // pix_buf[i_time * pixel_buf_man.steps[0]] = pixel_offset;
        }
    }

    return pixel_buf_man.ret_val;
}

// NB: This function assumes nearest neigbor. It doesn't look like
// it can be generalized with its current interface, since it returns
// an [ndet,ntime,{y,x,...}] array which can't handle multiple pixels per
// sample
template<typename C, typename P, typename S>
py::object ProjectionEngine<C,P,S>::pointing_matrix(
    py::object pbore, py::object pofs, py::object response, py::object pixel, py::object proj)
{
    auto _none = py::object();

    auto pointer = Pointer<C>();
    pointer.TestInputs(_none, pbore, pofs, _none, _none);
    //_pixelizor.TestInputs(_none, _none, _none, _none, _none);

    int n_det = pointer.DetCount();
    int n_time = pointer.TimeCount();

    auto pixel_buf_man = SignalSpace<int32_t>(
        pixel, "pixel", NPY_INT32, n_det, n_time, P::index_count);
    auto proj_buf_man = SignalSpace<FSIGNAL>(
        proj, "proj", FSIGNAL_NPY_TYPE, n_det, n_time, S::comp_count);
    auto _response = BufferWrapper<FSIGNAL>(
         "response", response, false, vector<int>{n_det,2});

#pragma omp parallel for
    for (int i_det = 0; i_det < n_det; ++i_det) {
        double dofs[4];
        pointer.InitPerDet(i_det, dofs);
        int* const pix_buf = pixel_buf_man.data_ptr[i_det];
        FSIGNAL* const proj_buf = proj_buf_man.data_ptr[i_det];
        const int step = pixel_buf_man.steps[0];
        int pixel_offset[P::index_count] = {-1};
        Response resp = get_response(_response, i_det);
        for (int i_time = 0; i_time < n_time; ++i_time) {
            double coords[4];
            FSIGNAL pf[S::comp_count];
            pointer.GetCoords(i_det, i_time, (double*)dofs, (double*)coords);
            _pixelizor.GetPixel(i_det, i_time, (double*)coords, pixel_offset);
            spin_proj_factors<S>(coords, resp, pf);

            for (int i_dim = 0; i_dim < P::index_count; i_dim++)
                pix_buf[i_time * pixel_buf_man.steps[0] +
                        i_dim * pixel_buf_man.steps[1]] = pixel_offset[i_dim];
            for (int i_spin = 0; i_spin < S::comp_count; i_spin++)
                proj_buf[i_time * proj_buf_man.steps[0] +
                         i_spin * proj_buf_man.steps[1]] = pf[i_spin];
        }
    }

    return py::make_tuple(pixel_buf_man.ret_val,
                          proj_buf_man.ret_val);
}

template<typename C, typename P, typename S>
py::object ProjectionEngine<C,P,S>::pixel_ranges(
    py::object pbore, py::object pofs, py::object map, int n_domain)
{
    auto _none = py::object();

    auto pointer = Pointer<C>();
    pointer.TestInputs(map, pbore, pofs, _none, _none);
    int n_det = pointer.DetCount();
    int n_time = pointer.TimeCount();

    bool from_map = !isNone(map);
    if (from_map)
        _pixelizor.TestInputs(map, true, false, S::comp_count);

    // Default value of the number of domains to set up is one per available
    // thread
    if (n_domain <= 0) {
        n_domain = 1;
        #ifdef _OPENMP
        n_domain = omp_get_max_threads();
        #endif
    }

    // ranges is [nbunch,ndomain,ndet][nrange]
    // we loop in series over bunches, then in parallel over
    // the domains in each bunch. This could be set up multiple
    // ways. Here, we use 2 bunches. The first has nthread domains,
    // and takes care of the samples that unambiguously hit a domain.
    // The second bunch contains just one pseudo-domain for all the samples
    // that straddle a domain edge. These will be done in series.
    vector<vector<vector<RangesInt32>>> ranges(2);
    auto & par_ranges = ranges[0];
    auto & ser_ranges = ranges[1];
    auto empty_range = vector<RangesInt32>(n_det, RangesInt32(n_time));
    for (int i=0; i<n_domain; i++)
        par_ranges.push_back(empty_range);
    ser_ranges.push_back(empty_range);

    // Find the domain for each det-sample and build ranges of consecutive
    // samples with the same domain
    #pragma omp parallel for
    for (int i_det = 0; i_det < n_det; ++i_det) {
        double dofs[4];
        pointer.InitPerDet(i_det, dofs);
        int last_domain = -1;
        int domain_start = 0;
        int pixinds[P::interp_count][P::index_count] = {-1};
        FSIGNAL pixweights[P::interp_count];
        for (int i_time = 0; i_time < n_time; ++i_time) {
            double coords[4];
            pointer.GetCoords(i_det, i_time, (double*)dofs, (double*)coords);
            // Get all the pixels this sample hits
            int n_point = _pixelizor.GetPixels(i_det, i_time, coords, pixinds, pixweights);
            // Assign a domain to this sample. If all pixels in this sample hit the same domain,
            // then this will be assigned to samp_domain. If the domain is ambiguous, then it will
            // be set to n_domain.
            int samp_domain = -1;
            for(int i_point = 0; i_point < n_point; ++i_point) {
                int point_domain = from_map ? _pixelizor.pix(0, pixinds[i_point])[0] : _pixelizor.stripe(pixinds[i_point], n_domain);
                if(i_point > 0 && point_domain != samp_domain) {
                    // This is a mixed-domain sample. Put it in the serial catagory
                    samp_domain = n_domain;
                    break;
                }
                samp_domain = point_domain;
            }
            // Did the domain status change?
            if (samp_domain != last_domain) {
                // Avoid triggering on the first sample, before we have established any domain
                if (last_domain >= 0) {
                    // Assign to either a parallel domain or the serial misc category
                    auto & target_ranges = last_domain < n_domain ? par_ranges[last_domain] : ser_ranges[0];
                    target_ranges[i_det].append_interval_no_check(domain_start, i_time);
                }
                domain_start = i_time;
                last_domain  = samp_domain;
            }
        }
        // Close the last interval
        if (last_domain >= 0) {
            auto & target_ranges = last_domain < n_domain ? par_ranges[last_domain] : ser_ranges[0];
            target_ranges[i_det].append_interval_no_check(domain_start, n_time);
        }
    }

    // Convert super vector to a list and return
    auto ivals = py::list();
    for (int i=0; i<ranges.size(); i++) {
        auto bunches = py::list();
        for (int j=0; j<ranges[i].size(); j++) {
            auto domains = py::list();
            for (int i_det=0; i_det<n_det; i_det++) {
                auto iv = py::cast(ranges[i][j][i_det]);
                domains.append(py::object(iv));
            }
            bunches.append(py::cast<py::object>(domains));
        }
        ivals.append(py::cast<py::object>(bunches));
    }
    return py::cast<py::object>(ivals);
}

template<typename C, typename P, typename S>
vector<int> ProjectionEngine<C,P,S>::tile_hits(
    py::object pbore, py::object pofs)
{
    auto _none = py::object();

    auto pointer = Pointer<C>();
    pointer.TestInputs(_none, pbore, pofs, _none, _none);
    int n_det = pointer.DetCount();
    int n_time = pointer.TimeCount();

    int n_tile = _pixelizor.tile_count();
    if (n_tile < 0)
        throw RuntimeError_exception("No tiles in this pixelization.");

    vector<int> hits(n_tile);
    vector<vector<int>> temp;

#pragma omp parallel
    {
        int n_domain = 1;
        #ifdef _OPENMP
        n_domain = omp_get_num_threads();
        #endif
#pragma omp single
        {
            for (int i=0; i<n_domain; i++)
                temp.push_back(vector<int>(n_tile));
        }

#pragma omp for
        for (int i_det = 0; i_det < n_det; ++i_det) {
            double dofs[4];
            pointer.InitPerDet(i_det, dofs);
            int pixinds[P::interp_count][P::index_count] = {-1};
            FSIGNAL pixweights[P::interp_count];
            int thread = 0;
            #ifdef _OPENMP
            thread = omp_get_thread_num();
            #endif
            for (int i_time = 0; i_time < n_time; ++i_time) {
                double coords[4];
                pointer.GetCoords(i_det, i_time, (double*)dofs, (double*)coords);
                int n_point = _pixelizor.GetPixels(i_det, i_time, coords, pixinds, pixweights);
                for(int i_point = 0; i_point < n_point; ++i_point)
                    temp[thread][pixinds[i_point][0]]++;
            }
        }
#pragma omp single
        {
            for (int i=0; i<n_domain; i++) {
                for (int j=0; j<n_tile; j++)
                    hits[j] += temp[i][j];
            }
        }
    }

    return hits;
}

//tile_ranges: create RangesInt32 information for n threads, such that
//each thread is active only on certain tiles.  tile_map should be a
//list of lists of tiles.
template<typename C, typename P, typename S>
py::object ProjectionEngine<C,P,S>::tile_ranges(
    py::object pbore, py::object pofs, py::list tile_lists)
{
    auto _none = py::object();

    auto pointer = Pointer<C>();
    pointer.TestInputs(_none, pbore, pofs, _none, _none);
    int n_det = pointer.DetCount();
    int n_time = pointer.TimeCount();

    int n_tile = _pixelizor.tile_count();
    if (n_tile < 0)
        throw RuntimeError_exception("No tiles in this pixelization.");
    int n_domain = py::len(tile_lists);

    // Make a vector that maps tile into thread.
    vector<int> thread_idx(n_tile, -1);
    for (int i=0; i<py::len(tile_lists); i++) {
        py::list tile_list = tile_lists[i];
        for (int j=0; j<py::len(tile_list); j++) {
            int tile_idx = PyLong_AsLong(py::object(tile_list[j]).ptr());
            thread_idx[tile_idx] = i;
       }
    }

    // ranges is [nbunch,ndomain,ndet][nrange]
    // we loop in series over bunches, then in parallel over
    // the domains in each bunch. This could be set up multiple
    // ways. Here, we use 2 bunches. The first has nthread domains,
    // and takes care of the samples that unambiguously hit a domain.
    // The second bunch contains just one pseudo-domain for all the samples
    // that straddle a domain edge. These will be done in series.
    vector<vector<vector<RangesInt32>>> ranges(2);
    auto & par_ranges = ranges[0];
    auto & ser_ranges = ranges[1];
    auto empty_range = vector<RangesInt32>(n_det, RangesInt32(n_time));
    for (int i=0; i<n_domain; i++)
        par_ranges.push_back(empty_range);
    ser_ranges.push_back(empty_range);

    #pragma omp parallel for
    for (int i_det = 0; i_det < n_det; ++i_det) {
        double dofs[4];
        pointer.InitPerDet(i_det, dofs);
        int last_domain = -1;
        int domain_start = 0;
        int pixinds[P::interp_count][P::index_count] = {-1};
        FSIGNAL pixweights[P::interp_count];
        for (int i_time = 0; i_time < n_time; ++i_time) {
            double coords[4];
            pointer.GetCoords(i_det, i_time, (double*)dofs, (double*)coords);
            // set samp_domain from thread_idx if all points for this sample agree, otherwise
            // set it to n_domain
            int samp_domain = -1;
            int n_point = _pixelizor.GetPixels(i_det, i_time, coords, pixinds, pixweights);
            for(int i_point = 0; i_point < n_point; ++i_point) {
                int point_domain = thread_idx[pixinds[i_point][0]];
                if(i_point > 0 && point_domain != samp_domain) {
                    // This is a mixed-domain sample. Put it in the serial catagory
                    samp_domain = n_domain;
                    break;
                }
                samp_domain = point_domain;
            }
            // Did the domain status change?
            if (samp_domain != last_domain) {
                // Avoid triggering on the first sample, before we have established any domain
                if (last_domain >= 0) {
                    // Assign to either a parallel domain or the serial misc category
                    auto & target_ranges = last_domain < n_domain ? par_ranges[last_domain] : ser_ranges[0];
                    target_ranges[i_det].append_interval_no_check(domain_start, i_time);
                }
                domain_start = i_time;
                last_domain  = samp_domain;
            }
        }
        // Close the last interval
        if (last_domain >= 0) {
            auto & target_ranges = last_domain < n_domain ? par_ranges[last_domain] : ser_ranges[0];
            target_ranges[i_det].append_interval_no_check(domain_start, n_time);
        }
    }

    // Convert super vector to a list and return
    auto ivals = py::list();
    for (int i=0; i<ranges.size(); i++) {
        auto bunches = py::list();
        for (int j=0; j<ranges[i].size(); j++) {
            auto domains = py::list();
            for (int i_det=0; i_det<n_det; i_det++) {
                auto iv = py::cast(ranges[i][j][i_det]);
                domains.append(py::object(iv));
            }
            bunches.append(py::cast<py::object>(domains));
        }
        ivals.append(py::cast<py::object>(bunches));
    }
    return py::cast<py::object>(ivals);
}

template<typename C, typename P, typename S>
py::object ProjectionEngine<C,P,S>::zeros(py::object shape)
{
    vector<int> dims;
    int dim;
    if (py::isinstance<int>(shape)) {
        // scalar
        dims.push_back(py::cast<int>(shape));
        return _pixelizor.zeros(dims);
    } else if (py::isinstance<py::tuple>(shape)) {
        // tuple
        py::tuple shp = py::cast<py::tuple>(shape);
        for (int i = 0; i < py::len(shp); ++i) {
            dims.push_back(py::cast<int>(shp[i]));
        }
        return _pixelizor.zeros(dims);
    } else {
        // invalid
        return py::object();
    }
}

template<typename C, typename P, typename S>
py::object ProjectionEngine<C,P,S>::from_map(
    py::object map, py::object pbore, py::object pofs, py::object response, py::object signal)
{
    auto _none = py::object();

    // Initialize pointer and _pixelizor.
    auto pointer = Pointer<C>();
    pointer.TestInputs(map, pbore, pofs, signal, _none);
    int n_det = pointer.DetCount();
    int n_time = pointer.TimeCount();

    // Confirm that map has the right meta-shape.
    _pixelizor.TestInputs(map, true, false, S::comp_count);

    // Get pointers to the signal and (optional) per-det weights.
    auto _signalspace = SignalSpace<FSIGNAL>(
            signal, "signal", FSIGNAL_NPY_TYPE, n_det, n_time);
    auto _response = BufferWrapper<FSIGNAL>(
         "response", response, false, vector<int>{n_det,2});

    #pragma omp parallel for
    for (int i_det = 0; i_det < n_det; ++i_det) {
        double dofs[4];
        pointer.InitPerDet(i_det, dofs);
        int pixinds[P::interp_count][P::index_count] = {-1};
        FSIGNAL pixweights[P::interp_count];
        FSIGNAL pf[S::comp_count];
        Response resp = get_response(_response, i_det);
        for (int i_time = 0; i_time < n_time; ++i_time) {
            double coords[4];
            pointer.GetCoords(i_det, i_time, (double*)dofs, (double*)coords);
            spin_proj_factors<S>(coords, resp, pf);
            FSIGNAL *sig = (_signalspace.data_ptr[i_det] + _signalspace.steps[0]*i_time);
            int n_point = _pixelizor.GetPixels(i_det, i_time, coords, pixinds, pixweights);
            for(int i_point = 0; i_point < n_point; ++i_point)
                for (int imap=0; imap<S::comp_count; ++imap)
                    *sig += *_pixelizor.pix(imap, pixinds[i_point]) * pf[imap] * pixweights[i_point];
        }
    }

    return _signalspace.ret_val;
}

template<typename C, typename P, typename S>
static
void to_map_single_thread(Pointer<C> &pointer,
                          BufferWrapper<FSIGNAL> & _response,
                          P &_pixelizor,
                          const vector<RangesInt32> & ivals,
                          BufferWrapper<FSIGNAL> &_det_weights,
                          SignalSpace<FSIGNAL> *_signalspace)
{
    int n_det = pointer.DetCount();
    for (int i_det = 0; i_det < n_det; ++i_det) {
        FSIGNAL det_wt = 1.;
        if (_det_weights->obj != NULL){
            det_wt = *(FSIGNAL*)((char*)_det_weights->buf +
                                 _det_weights->strides[0]*i_det);
            if (det_wt == 0.) continue;
        }
        double dofs[4];
        double coords[4];
        FSIGNAL pf[S::comp_count];
        pointer.InitPerDet(i_det, dofs);
        Response resp = get_response(_response, i_det);
        // Pointing matrix interpolation stuff
        int pixinds[P::interp_count][P::index_count] = {-1};
        FSIGNAL pixweights[P::interp_count] = {0};
        for (auto const &rng: ivals[i_det].segments) {
            for (int i_time = rng.first; i_time < rng.second; ++i_time) {
                pointer.GetCoords(i_det, i_time, (double*)dofs, (double*)coords);
                spin_proj_factors<S>(coords, resp, pf);
                const FSIGNAL sig = _signalspace->data_ptr[i_det][_signalspace->steps[0]*i_time];
                // In interpolated mapmaking like bilinear mampamking, each sample hits multipe
                // pixels, each with its own weight.
                int n_point = _pixelizor.GetPixels(i_det, i_time, coords, pixinds, pixweights);
                for(int i_point = 0; i_point < n_point; ++i_point)
                    for (int i_map=0; i_map<S::comp_count; ++i_map)
                        *_pixelizor.pix(i_map, pixinds[i_point]) += sig * pf[i_map] * pixweights[i_point] * det_wt;
            }
        }
    }
}

template<typename C, typename P, typename S>
static
void to_weight_map_single_thread(Pointer<C> &pointer,
                                 BufferWrapper<FSIGNAL> & _response,
                                 P &_pixelizor,
                                 vector<RangesInt32> ivals,
                                 BufferWrapper<FSIGNAL> &_det_weights)
{
    int n_det = pointer.DetCount();
    for (int i_det = 0; i_det < n_det; ++i_det) {
        FSIGNAL det_wt = 1.;
        if (_det_weights->obj != NULL){
            det_wt = *(FSIGNAL*)((char*)_det_weights->buf +
                                 _det_weights->strides[0]*i_det);
            if (det_wt == 0.) continue;
        }
        double dofs[4];
        double coords[4];
        FSIGNAL pf[S::comp_count];
        pointer.InitPerDet(i_det, dofs);
        Response resp = get_response(_response, i_det);
        int pixinds[P::interp_count][P::index_count] = {-1};
        FSIGNAL pixweights[P::interp_count] = {0};
        for (auto const &rng: ivals[i_det].segments) {
            for (int i_time = rng.first; i_time < rng.second; ++i_time) {
                pointer.GetCoords(i_det, i_time, (double*)dofs, (double*)coords);
                spin_proj_factors<S>(coords, resp, pf);

                int n_point = _pixelizor.GetPixels(i_det, i_time, coords, pixinds, pixweights);
                // Is this enough? Do we need a double loop over i_point?
                // The weight map is supposed to be the diagonal of P'WP
                // P_pi W_ii P_ip = P_pi**2 * W_ii. So we just need to square the pixweight.
                // The i_point,j_point stuff would be off-diagonal elements which are outside the
                // scope of this function
                for(int i_point = 0; i_point < n_point; ++i_point)
                    for (int imap=0; imap<S::comp_count; ++imap)
                        for (int jmap=imap; jmap<S::comp_count; ++jmap)
                            *_pixelizor.wpix(imap, jmap, pixinds[i_point]) += pf[imap] * pf[jmap] * pixweights[i_point] * pixweights[i_point] * det_wt;
            }
        }
    }
}

static
vector<vector<vector<RangesInt32>>> derive_ranges(
    py::object thread_intervals, int n_det, int n_time,
    std::string arg_name)
{
    // The first index of the returned object should correspond to
    // (OMP) execution thread; the second index is over detectors.
    // The input object can be any of the following:
    //
    // - None: Result will be full coverage for all detectors in a
    //   single thread.
    // - List of Ranges: This will be promoted to a single thread.
    // - List containing N lists of Ranges: N threads, with each list
    //   used by corresponding thread.
    // - List containing M lists of N[M] lists of Ranges: M bunches
    //   which will be looped over in serial, each containing N[M]
    //   threads-lists that will be done in parallel.
    vector<vector<vector<RangesInt32>>> ivals;

    if (isNone(thread_intervals)) {
        // It's None. Generate a single bunch with a single-thread covering all samples
        auto r = RangesInt32(n_time).add_interval(0, n_time);
        vector<vector<RangesInt32>> v(1, vector<RangesInt32>(n_det, r));
        ivals.push_back(v);
    } else if (py::isinstance<RangesInt32>(thread_intervals[0])) {
        // It's a RangesMatrix (ndet,nranges). Promote to single thread, single bunch
        ivals.push_back(vector<vector<RangesInt32>>(1, extract_ranges<int32_t>(thread_intervals)));
    } else if (py::isinstance<RangesInt32>(thread_intervals[0][0])) {
        // It's a per-thread RangesMatrix (nthread,ndet,nranges). Promote to single bunch
        vector<vector<RangesInt32>> bunch;
        py::list leading = py::cast<py::list>(thread_intervals);
        for (int i=0; i<py::len(leading); i++)
            bunch.push_back(extract_ranges<int32_t>(leading[i]));
        ivals.push_back(bunch);
    } else if (py::isinstance<RangesInt32>(thread_intervals[0][0][0])) {
        // It's a full multi-bunch (nbunch,nthread,ndet,nranges) thing.
        py::list leading_i = py::cast<py::list>(thread_intervals);
        for (int i=0; i<py::len(leading_i); i++) {
            vector<vector<RangesInt32>> bunch;
            py::list leading_j = py::cast<py::list>(leading_i[i]);
            for (int j=0; j<py::len(leading_j); j++)
                bunch.push_back(extract_ranges<int32_t>(leading_j[j]));
            ivals.push_back(bunch);
        }
    } else {
        // This should not happen
        assert(false);
    }
    // Check that these all have the right shape. Maybe consider a standard
    // for loop instead of foreach to give more useful error messages using
    // the index
    for(auto const & bunch: ivals)
    for(auto const & ival: bunch) {
        if (ival.size() != n_det) {
            std::ostringstream err;
            err << "Expected RangesMatrix with n_det=" << n_det
                << " but encountered n_det=" << ival.size();
            throw shape_exception(arg_name, err.str());
        }
        for (auto const &ri: ival) {
            if (ri.count != n_time) {
                std::ostringstream err;
                err << "Expected RangesMatrix with n_time=" << n_time
                    << " but encountered n_time=" << ri.count;
                throw shape_exception(arg_name, err.str());
            }
        }
    }
    return ivals;
}

template<typename C, typename P, typename S>
py::object ProjectionEngine<C,P,S>::to_map(
    py::object map, py::object pbore, py::object pofs, py::object response,
    py::object signal, py::object det_weights, py::object thread_intervals)
{
    //Initialize it / check inputs.
    auto pointer = Pointer<C>();
    std::cerr << "::to_map calling pointer.TestInputs" << std::endl;
    pointer.TestInputs(map, pbore, pofs, signal, det_weights);
    std::cerr << "::to_map calling pointer.DetCount" << std::endl;
    int n_det = pointer.DetCount();
    std::cerr << "::to_map calling pointer.TimeCount" << std::endl;
    int n_time = pointer.TimeCount();

    //Do we need a map?  Now is the time.
    if (isNone(map)) {
        std::cerr << "::to_map calling pixelizer.zeros to init map" << std::endl;
        std::cerr << "::to_map zeros " << S::comp_count << " components" << std::endl;
        map = _pixelizor.zeros(vector<int>{S::comp_count});
    }

    // Confirm that map has the right meta-shape.
    std::cerr << "::to_map calling pixelizer.TestInputs" << std::endl;
    _pixelizor.TestInputs(map, true, false, S::comp_count);

    // Get pointers to the signal and (optional) per-det weights.
    std::cerr << "::to_map create SignalSpace" << std::endl;
    auto _signalspace = SignalSpace<FSIGNAL>(
            signal, "signal", FSIGNAL_NPY_TYPE, n_det, n_time);
    std::cerr << "::to_map create det_weights" << std::endl;
    auto _det_weights = BufferWrapper<FSIGNAL>(
         "det_weights", det_weights, true, vector<int>{n_det});
    std::cerr << "::to_map create response" << std::endl;
    auto _response = BufferWrapper<FSIGNAL>(
         "response", response, false, vector<int>{n_det,2});

    // For multi-threading, the principle here is that we loop serially
    // over bunches, and then inside each block all threads loop over
    // all detectors in parallel, but the sample ranges are pixel-disjoint.
    std::cerr << "::to_map call derive_ranges" << std::endl;
    auto bunches = derive_ranges(thread_intervals, n_det, n_time,
                               "thread_intervals");

    // First loop over serial bunches
    for(int i_bunch = 0; i_bunch < bunches.size(); i_bunch++) {
        const auto & ivals = bunches[i_bunch];
        // Then loop over parallel bunches. This works even if omp is not enabled,
        // or if ivals.size() == 1
        #pragma omp parallel for
        for (int i_thread = 0; i_thread < ivals.size(); i_thread++) {
            to_map_single_thread<C,P,S>(pointer, _response, _pixelizor, ivals[i_thread], _det_weights, &_signalspace);
        }
    }
    std::cerr << "::to_map finished loop over bunches" << std::endl;
    return map;
}

template<typename C, typename P, typename S>
py::object ProjectionEngine<C,P,S>::to_weight_map(
    py::object map, py::object pbore, py::object pofs, py::object response,
    py::object det_weights, py::object thread_intervals)
{
    auto _none = py::object();

    //Initialize it / check inputs.
    auto pointer = Pointer<C>();
    pointer.TestInputs(map, pbore, pofs, _none, det_weights);
    int n_det = pointer.DetCount();
    int n_time = pointer.TimeCount();

    //Do we need a map?  Now is the time.
    if (isNone(map)) {
        map = _pixelizor.zeros(vector<int>{S::comp_count,S::comp_count});
    }

    // Confirm that map has the right meta-shape.
    _pixelizor.TestInputs(map, false, true, S::comp_count);

    // Get pointer to (optional) per-det weights.
    auto _det_weights = BufferWrapper<FSIGNAL>(
         "det_weights", det_weights, true, vector<int>{n_det});
    auto _response = BufferWrapper<FSIGNAL>(
         "response", response, false, vector<int>{n_det,2});

    // For multi-threading, the principle here is that we loop serially
    // over bunches, and then inside each block all threads loop over
    // all detectors in parallel, but the sample ranges are pixel-disjoint.
    auto bunches = derive_ranges(thread_intervals, n_det, n_time,
                               "thread_intervals");

    // First loop over serial bunches
    for(int i_bunch = 0; i_bunch < bunches.size(); i_bunch++) {
        const auto & ivals = bunches[i_bunch];
        // Then loop over parallel bunches. This works even if omp is not enabled,
        // or if ivals.size() == 1
        #pragma omp parallel for
        for (int i_thread = 0; i_thread < ivals.size(); i_thread++) {
            to_weight_map_single_thread<C,P,S>(pointer, _response, _pixelizor, ivals[i_thread], _det_weights);
        }
    }
    return map;
}

// ProjEng_Precomp
//
// The ProjEng_Precomp may be used for very fast projection
// operations, provided you have precomputed pixel index and spin
// projection factors for every sample.  This is agnostic of
// coordinate system, and so not crazily templated. It is
// hard-coded to nearest neighbor interpolation though
// (though one could emulate interpolation by calling it repeatedly
// with different pixel indices and weights).
//

template<typename TilingSys>
class ProjEng_Precomp {
public:
    ProjEng_Precomp() {};
    py::object from_map(py::object map, py::object pixel_index, py::object spin_proj,
                        py::object signal);
    py::object to_map(py::object map, py::object pixel_index, py::object spin_proj,
                      py::object signal, py::object weights, py::object thread_intervals);
    py::object to_weight_map(py::object map, py::object pixel_index, py::object spin_proj,
                             py::object weights, py::object thread_intervals);
};


template<typename TilingSys>
py::object ProjEng_Precomp<TilingSys>::from_map(
    py::object map, py::object pixel_index, py::object spin_proj,
    py::object signal)
{
    // You won't get far without pixel_index, so use that to nail down
    // the n_time and n_det.
    auto pixel_buf_man = SignalSpace<int32_t>(
        pixel_index, "pixel_index", NPY_INT32, -1, -1, -1);
    int n_det = pixel_buf_man.dims[0];
    int n_time = pixel_buf_man.dims[1];
    int n_pix = pixel_buf_man.dims[2];

    // Similarly, spin_proj tells you the number of components.
    auto spin_proj_man = SignalSpace<FSIGNAL>(
        spin_proj, "spin_proj", FSIGNAL_NPY_TYPE, n_det, n_time, -1);
    int n_spin = spin_proj_man.dims[2];

    auto tiling = Pixelizor2_Flat<TilingSys>();
    tiling.TestInputs(map, true, false, n_spin);

    // Check the signal.
    auto signal_man = SignalSpace<FSIGNAL>(
        signal, "signal", FSIGNAL_NPY_TYPE, n_det, n_time);

    // Below we assume the pixel sub-indices are close-packed.
    if (pixel_buf_man.steps[1] != 1)
        throw shape_exception("pixel_index",
                              "Fast dimension of pixel indices must be close-packed.");

#pragma omp parallel for
    for (int i_det = 0; i_det < n_det; ++i_det) {
        for (int i_time = 0; i_time < n_time; ++i_time) {
            const int *pixel_ofs = pixel_buf_man.data_ptr[i_det] +
                pixel_buf_man.steps[0]*i_time;
            if (pixel_ofs[0] < 0)
                continue;
            FSIGNAL sig = 0;
            for (int i_spin = 0; i_spin < n_spin; ++i_spin) {
                FSIGNAL pf = *(spin_proj_man.data_ptr[i_det] +
                               spin_proj_man.steps[0] * i_time +
                               spin_proj_man.steps[1] * i_spin);
                sig += pf * *tiling.pix(i_spin, pixel_ofs);
            }
            *(signal_man.data_ptr[i_det] + signal_man.steps[0]*i_time) += sig;
        }
    }

    return signal_man.ret_val;
}

template<typename TilingSys>
static
void precomp_to_map_single_thread(Pixelizor2_Flat<TilingSys> &tiling,
                                  const SignalSpace<int32_t> &pixel_buf_man,
                                  const SignalSpace<FSIGNAL> &spin_proj_man,
                                  vector<RangesInt32> ivals,
                                  BufferWrapper<FSIGNAL> &_det_weights,
                                  SignalSpace<FSIGNAL> *_signalspace)
{
    int n_det = pixel_buf_man.dims[0];
    int n_spin = spin_proj_man.dims[2];
    assert(spin_proj_man.steps[1] == 1); // assumed below

    for (int i_det = 0; i_det < n_det; ++i_det) {
        FSIGNAL det_wt = 1.;
        if (_det_weights->obj != NULL){
            det_wt = *(FSIGNAL*)((char*)_det_weights->buf + _det_weights->strides[0]*i_det);
            if (det_wt == 0.) continue;
        }
        for (auto const &rng: ivals[i_det].segments) {
            for (int i_time = rng.first; i_time < rng.second; ++i_time) {
                const int *pixel_offset = pixel_buf_man.data_ptr[i_det] +
                    pixel_buf_man.steps[0]*i_time;
                if (pixel_offset[0] < 0)
                    continue;
                const FSIGNAL sig = *(_signalspace->data_ptr[i_det] +
                                      _signalspace->steps[0]*i_time);
                const FSIGNAL *pf = (spin_proj_man.data_ptr[i_det] +
                                     spin_proj_man.steps[0] * i_time);
                for (int i_spin = 0; i_spin < n_spin; ++i_spin) {
                    *tiling.pix(i_spin, pixel_offset) += det_wt * sig * pf[i_spin];
                }
            }
        }
    }
}

template<typename TilingSys>
static
void precomp_to_weight_map_single_thread(Pixelizor2_Flat<TilingSys> &tiling,
                                         const SignalSpace<int32_t> &pixel_buf_man,
                                         const SignalSpace<FSIGNAL> &spin_proj_man,
                                         vector<RangesInt32> ivals,
                                         BufferWrapper<FSIGNAL> &_det_weights)
{
    int n_det = pixel_buf_man.dims[0];
    int n_spin = spin_proj_man.dims[2];
    assert(spin_proj_man.steps[1] == 1); // assumed below
    for (int i_det = 0; i_det < n_det; ++i_det) {
        FSIGNAL det_wt = 1.;
        if (_det_weights->obj != NULL){
            det_wt = *(FSIGNAL*)((char*)_det_weights->buf + _det_weights->strides[0]*i_det);
            if (det_wt == 0.) continue;
        }
        for (auto const &rng: ivals[i_det].segments) {
            for (int i_time = rng.first; i_time < rng.second; ++i_time) {
                const int *pixel_offset = pixel_buf_man.data_ptr[i_det] +
                    pixel_buf_man.steps[0]*i_time;
                if (pixel_offset[0] < 0)
                    continue;
                const FSIGNAL *pf = (spin_proj_man.data_ptr[i_det] +
                                     spin_proj_man.steps[0] * i_time);
                for (int imap=0; imap < n_spin; ++imap)
                    for (int jmap=imap; jmap < n_spin; ++jmap)
                        *tiling.wpix(imap, jmap, pixel_offset) += pf[imap] * pf[jmap] * det_wt;
            }
        }
    }
}

template<typename TilingSys>
py::object ProjEng_Precomp<TilingSys>::to_map(
    py::object map, py::object pixel_index, py::object spin_proj,
    py::object signal, py::object det_weights, py::object thread_intervals)
{
    // You won't get far without pixel_index, so use that to nail down
    // the n_time and n_det.
    auto pixel_buf_man = SignalSpace<int32_t>(
        pixel_index, "pixel_index", NPY_INT32, -1, -1, -1);
    int n_det = pixel_buf_man.dims[0];
    int n_time = pixel_buf_man.dims[1];
    int n_pix = pixel_buf_man.dims[2];
    std::cerr << "BIND to_map: n_det = " << n_det << ", n_time = " << n_time << ", n_pix = " << n_pix << std::endl;

    // Similarly, spin_proj tells you the number of components.
    auto spin_proj_man = SignalSpace<FSIGNAL>(
        spin_proj, "spin_proj", FSIGNAL_NPY_TYPE, n_det, n_time, -1);
    int n_spin = spin_proj_man.dims[2];
    std::cerr << "BIND to_map: n_spin / n_components = " << n_spin << std::endl;

    // Unlike the on-the-fly class, we aren't able to make the map
    // here, because there's no initialized pixelizor.
    auto pixelizor = Pixelizor2_Flat<TilingSys>();
    pixelizor.TestInputs(map, true, false, n_spin);

    // Check the signal.
    auto _signalspace = SignalSpace<FSIGNAL>(
        signal, "signal", FSIGNAL_NPY_TYPE, n_det, n_time);

    BufferWrapper<FSIGNAL> _det_weights("det_weights", det_weights, true,
                                        vector<int>{n_det});

    // Below we assume the pixel sub-indices are close-packed.
    if (pixel_buf_man.steps[1] != 1)
        throw shape_exception("pixel_index",
                              "Fast dimension of pixel indices must be close-packed.");

    auto bunches = derive_ranges(thread_intervals, n_det, n_time, "thread_intervals");

    // First loop over serial bunches
    for(int i_bunch = 0; i_bunch < bunches.size(); i_bunch++) {
        const auto & ivals = bunches[i_bunch];
        // Then loop over parallel bunches. This works even if omp is not enabled,
        // or if ivals.size() == 1
        #pragma omp parallel for
        for (int i_thread = 0; i_thread < ivals.size(); i_thread++) {
            precomp_to_map_single_thread<TilingSys>(
                pixelizor, pixel_buf_man, spin_proj_man,
                ivals[i_thread], _det_weights, &_signalspace);
        }
    }
    return map;
}

template<typename TilingSys>
py::object ProjEng_Precomp<TilingSys>::to_weight_map(
    py::object map, py::object pixel_index, py::object spin_proj,
    py::object det_weights, py::object thread_intervals)
{
    // You won't get far without pixel_index, so use that to nail down
    // the n_time and n_det.
    auto pixel_buf_man = SignalSpace<int32_t>(
        pixel_index, "pixel_index", NPY_INT32, -1, -1, -1);
    int n_det = pixel_buf_man.dims[0];
    int n_time = pixel_buf_man.dims[1];
    int n_pix = pixel_buf_man.dims[2];

    // Similarly, spin_proj tells you the number of components.
    auto spin_proj_man = SignalSpace<FSIGNAL>(
        spin_proj, "spin_proj", FSIGNAL_NPY_TYPE, n_det, n_time, -1);
    int n_spin = spin_proj_man.dims[2];

    // Unlike the on-the-fly class, we aren't able to make the map
    // here, because there's no initialized pixelizor.
    auto pixelizor = Pixelizor2_Flat<TilingSys>();
    pixelizor.TestInputs(map, false, true, n_spin);

    BufferWrapper<FSIGNAL> _det_weights("det_weights", det_weights, true,
                                        vector<int>{n_det});

    // Below we assume the pixel sub-indices are close-packed.
    if (pixel_buf_man.steps[1] != 1)
        throw shape_exception("pixel_index",
                              "Fast dimension of pixel indices must be close-packed.");

    auto bunches = derive_ranges(thread_intervals, n_det, n_time, "thread_intervals");

    // First loop over serial bunches
    for(int i_bunch = 0; i_bunch < bunches.size(); i_bunch++) {
        const auto & ivals = bunches[i_bunch];
        // Then loop over parallel bunches. This works even if omp is not enabled,
        // or if ivals.size() == 1
        #pragma omp parallel for
        for (int i_thread = 0; i_thread < ivals.size(); i_thread++) {
            precomp_to_weight_map_single_thread<TilingSys>(
                pixelizor, pixel_buf_man, spin_proj_man,
                ivals[i_thread], _det_weights);
        }
    }
    return map;
}

typedef ProjEng_Precomp<NonTiled> ProjEng_Precomp_NonTiled;
typedef ProjEng_Precomp<Tiled> ProjEng_Precomp_Tiled;

#define PROJENG(PIX, SPIN, TILING) ProjEng_##PIX##_##SPIN##_##TILING
#define PROJENG_INTERP(PIX, SPIN, TILING, INTERP) ProjEng_##PIX##_##SPIN##_##TILING##_##INTERP

#define TYPEDEF_TILING(PIX, SPIN, TILING) \
    typedef ProjectionEngine<Proj ## PIX,Pixelizor2_Flat<TILING>,Spin##SPIN> \
        PROJENG(PIX, SPIN, TILING); \
    typedef ProjectionEngine<Proj ## PIX,Pixelizor2_Flat<TILING,Bilinear>,Spin##SPIN> \
        PROJENG_INTERP(PIX, SPIN, TILING, Bilinear);

typedef ProjectionEngine<ProjCAR, Pixelizor_Healpix<NonTiled>, SpinT> ProjEng_HP_T_NonTiled;
typedef ProjectionEngine<ProjCAR, Pixelizor_Healpix<NonTiled>, SpinQU> ProjEng_HP_QU_NonTiled;
typedef ProjectionEngine<ProjCAR, Pixelizor_Healpix<NonTiled>, SpinTQU> ProjEng_HP_TQU_NonTiled;
typedef ProjectionEngine<ProjCAR, Pixelizor_Healpix<Tiled>, SpinT> ProjEng_HP_T_Tiled;
typedef ProjectionEngine<ProjCAR, Pixelizor_Healpix<Tiled>, SpinQU> ProjEng_HP_QU_Tiled;
typedef ProjectionEngine<ProjCAR, Pixelizor_Healpix<Tiled>, SpinTQU> ProjEng_HP_TQU_Tiled;

#define TYPEDEF_SPIN(PIX, SPIN)                 \
    TYPEDEF_TILING(PIX, SPIN, Tiled)            \
    TYPEDEF_TILING(PIX, SPIN, NonTiled)

#define TYPEDEF_PIX(PIX)                         \
    TYPEDEF_SPIN(PIX, T)                         \
    TYPEDEF_SPIN(PIX, QU)                        \
    TYPEDEF_SPIN(PIX, TQU)

TYPEDEF_PIX(Flat)
TYPEDEF_PIX(Quat)
TYPEDEF_PIX(CAR)
TYPEDEF_PIX(CEA)
TYPEDEF_PIX(ARC)
TYPEDEF_PIX(SIN)
TYPEDEF_PIX(TAN)
TYPEDEF_PIX(ZEA)

#define STRINGIFY(X) #X

#define EXPORT_ENGINE(CLASSNAME)                                        \
    py::class_<CLASSNAME>(m, STRINGIFY(CLASSNAME))                      \
    .def(py::init<py::object>())                                        \
    .def_property_readonly("index_count", &CLASSNAME::index_count)               \
    .def_property_readonly("comp_count", &CLASSNAME::comp_count)                 \
    .def("coords", &CLASSNAME::coords)                                  \
    .def("pixels", &CLASSNAME::pixels)                                  \
    .def("tile_hits", &CLASSNAME::tile_hits)                            \
    .def("tile_ranges", &CLASSNAME::tile_ranges)                        \
    .def("pointing_matrix", &CLASSNAME::pointing_matrix)                \
    .def("pixel_ranges", &CLASSNAME::pixel_ranges)                      \
    .def("zeros", &CLASSNAME::zeros)                                    \
    .def("from_map", &CLASSNAME::from_map)                              \
    .def("to_map", &CLASSNAME::to_map)                                  \
    .def("to_weight_map", &CLASSNAME::to_weight_map)                    \
    ;


#define EXPORT_TILING(PIX, SPIN, TILING)        \
    EXPORT_ENGINE(PROJENG(PIX, SPIN, TILING)) \
    EXPORT_ENGINE(PROJENG_INTERP(PIX, SPIN, TILING, Bilinear))

#define EXPORT_SPIN(PIX, SPIN)                                          \
    EXPORT_TILING(PIX, SPIN, Tiled)                                     \
    EXPORT_TILING(PIX, SPIN, NonTiled)

#define EXPORT_PIX(PIX) \
    EXPORT_SPIN(PIX, T)                                                 \
    EXPORT_SPIN(PIX, QU)                                                \
    EXPORT_SPIN(PIX, TQU)

#define EXPORT_PRECOMP(CLASSNAME)                                       \
    py::class_<CLASSNAME>(m, #CLASSNAME)                                   \
    .def("from_map", &CLASSNAME::from_map)                              \
    .def("to_map", &CLASSNAME::to_map)                                  \
    .def("to_weight_map", &CLASSNAME::to_weight_map)                    \
;


// There are probably better ways to do this..
template<typename T>
inline
int _index_count(const T &) { return T::index_count; }


void register_projection(py::module_ & m) {
    EXPORT_PIX(Flat);
    EXPORT_PIX(Quat);
    EXPORT_PIX(CAR);
    EXPORT_PIX(CEA);
    EXPORT_PIX(ARC);
    EXPORT_PIX(SIN);
    EXPORT_PIX(TAN);
    EXPORT_PIX(ZEA);

    EXPORT_PRECOMP(ProjEng_Precomp_NonTiled);
    EXPORT_PRECOMP(ProjEng_Precomp_Tiled);

    EXPORT_ENGINE(ProjEng_HP_T_NonTiled);
    EXPORT_ENGINE(ProjEng_HP_QU_NonTiled);
    EXPORT_ENGINE(ProjEng_HP_TQU_NonTiled);
    EXPORT_ENGINE(ProjEng_HP_T_Tiled);
    EXPORT_ENGINE(ProjEng_HP_QU_Tiled);
    EXPORT_ENGINE(ProjEng_HP_TQU_Tiled);

    return;
}
