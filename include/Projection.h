#pragma once

#include <vector>

#include <nanobind/nanobind.h>

#include "exceptions.h"
#include "numpy_assist.h"

using namespace std;

namespace nb = nanobind;


// For detector timestreams, a.k.a. "signal", float32 is sufficient in
// most cases.  We don't want to template this, but let's leave our
// options open by using this typedef.

typedef float FSIGNAL;
#define FSIGNAL_NPY_TYPE NPY_FLOAT32
#define FSIGNAL_BUFFER_FORMAT "f"


template <typename DTYPE>
class SignalSpace {
public:
    SignalSpace(nb::object input, std::string var_name,
                int dtype, int n_det, int n_time);
    SignalSpace(nb::object input, std::string var_name,
                int dtype, int n_det, int n_time, int n_thirdaxis);
    ~SignalSpace() { if (data_ptr) free(data_ptr); };

    // Use arrays rather than vectors for anything that will be in
    // inner loop (even if it's inlined, the compiler might worry).
    DTYPE **data_ptr = NULL;
    int steps[64];
    vector<int> dims;

    vector<BufferWrapper<DTYPE>> bw;
    nb::object ret_val;

private:
    bool _Validate(nb::object input, std::string var_name,
                   int dtype);
};


// template<typename C, typename P, typename S>
template<typename CoordSys, typename PixelSys, typename SpinSys>
class ProjectionEngine {
public:
    //ProjectionEngine(PixelSys pixelizor);
    ProjectionEngine(nb::object pix_args);
    nb::object coords(nb::object pbore, nb::object pofs,
                      nb::object coord);
    nb::object pixels(nb::object pbore, nb::object pofs, nb::object pixel);
    vector<int> tile_hits(nb::object pbore, nb::object pofs);
    nb::object tile_ranges(nb::object pbore, nb::object pofs, nb::list tile_lists);
    nb::object pointing_matrix(nb::object pbore, nb::object pofs, nb::object response,
                               nb::object pixel, nb::object proj);
    nb::object zeros(nb::object shape);
    nb::object pixel_ranges(nb::object pbore, nb::object pofs, nb::object map, int n_domain=-1);
    nb::object from_map(nb::object map, nb::object pbore, nb::object pofs,
                        nb::object response, nb::object signal);
    nb::object to_map(nb::object map, nb::object pbore, nb::object pofs, nb::object response,
                      nb::object signal, nb::object det_weights,
                      nb::object thread_intervals);
    nb::object to_weight_map(nb::object map, nb::object pbore, nb::object pofs,
                  nb::object response, nb::object det_weights, nb::object thread_intervals);

    int comp_count() const;
    int index_count() const;

private:
    PixelSys _pixelizor;
};


void register_projection(nb::module_ & m);
