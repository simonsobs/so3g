#include <boost/python.hpp>
#include "exceptions.h"
#include "numpy_assist.h"

namespace bp = boost::python;

/* For detector timestreams, a.k.a. "signal", float32 is sufficient in
 * most cases.  We don't want to template this, but let's leave our
 * options open by using this typedef. */
typedef float FSIGNAL;
#define FSIGNAL_NPY_TYPE NPY_FLOAT32
#define FSIGNAL_BUFFER_FORMAT "f"


class Pixelizor2_Flat {
public:
    static const int index_count = 2;
    Pixelizor2_Flat() {};
    Pixelizor2_Flat(int ny, int nx,
                    double dy=1., double dx=1.,
                    double iy0=0., double ix0=0.);
    ~Pixelizor2_Flat() {};
    bool TestInputs(bp::object &map, bp::object &pbore, bp::object &pdet,
                    bp::object &signal, bp::object &det_weights);
    bp::object zeros(int count);
    void GetPixel(int i_det, int i_time, const double *coords, int *pixel_index);
    int crpix[2];
    double cdelt[2];
    int naxis[2];
};

class Pixelizor2_Flat_Tiled {
public:
    static const int index_count = 3;
    Pixelizor2_Flat_Tiled() {};
    Pixelizor2_Flat_Tiled(Pixelizor2_Flat _parent_pix,
                          int tiley, int tilex);
    Pixelizor2_Flat_Tiled(int ny, int nx,
                          double dy=1., double dx=1.,
                          double iy0=0., double ix0=0.,
                          int tiley=0, int tilex=0);
    ~Pixelizor2_Flat_Tiled() {};
    bool TestInputs(bp::object &map, bp::object &pbore, bp::object &pdet,
                    bp::object &signal, bp::object &det_weights);
    bp::object zeros(int count);
    bp::object zeros_tiled(int count, bp::object active_tiles);
    void GetPixel(int i_det, int i_time, const double *coords, int *pixel_index);
private:
    Pixelizor2_Flat parent_pix;
    int tile_shape[2];
};

template <typename DTYPE>
class SignalSpace {
public:
    SignalSpace(bp::object input, std::string var_name,
                int dtype, int n_det, int n_time);
    SignalSpace(bp::object input, std::string var_name,
                int dtype, int n_det, int n_time, int n_thirdaxis);
    ~SignalSpace() { if (data_ptr) free(data_ptr); };

    // Use arrays rather than vectors for anything that will be in
    // inner loop (even if it's inlined, the compiler might worry).
    DTYPE **data_ptr = NULL;
    int steps[64];
    vector<int> dims;

    vector<BufferWrapper<DTYPE>> bw;
    bp::object ret_val;

private:
    bool _Validate(bp::object input, std::string var_name,
                   int dtype);
};


// template<typename C, typename P, typename T, typename S>
template<typename CoordSys, typename PixelSys, typename TilingSys, typename SpinSys>
class ProjectionEngine {
public:
    ProjectionEngine(PixelSys pixelizor);
    bp::object coords(bp::object pbore, bp::object pofs,
                      bp::object coord);
    bp::object pixels(bp::object pbore, bp::object pofs, bp::object pixel);
    bp::object pointing_matrix(bp::object pbore, bp::object pofs,
                               bp::object pixel, bp::object proj);
    bp::object from_map(bp::object map, bp::object pbore, bp::object pofs,
                        bp::object signal, bp::object weights);
    // bp::object to_map(bp::object map, bp::object pbore, bp::object pofs,
    //                   bp::object signal, bp::object weights);
    // bp::object to_map_omp(bp::object map, bp::object pbore, bp::object pofs,
    //                       bp::object signal, bp::object weights,
    //                       bp::object thread_intervals);
    // bp::object to_weight_map(bp::object map, bp::object pbore, bp::object pofs,
    //                          bp::object signal, bp::object weights);
    // bp::object to_weight_map_omp(bp::object map, bp::object pbore, bp::object pofs,
    //                              bp::object signal, bp::object weights,
    //                              bp::object thread_intervals);
    // bp::object pixel_ranges(bp::object pbore, bp::object pofs);
private:
    PixelSys _pixelizor;
};
