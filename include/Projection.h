#pragma once

#include <vector>

#include <pybind11/pybind11.h>

#include "exceptions.h"
#include "numpy_assist.h"

using namespace std;

namespace py = pybind11;


// For detector timestreams, a.k.a. "signal", float32 is sufficient in
// most cases.  We don't want to template this, but let's leave our
// options open by using this typedef.

typedef float FSIGNAL;
#define FSIGNAL_NPY_TYPE NPY_FLOAT32
#define FSIGNAL_BUFFER_FORMAT "f"


template <typename DTYPE>
class SignalSpace {
public:
    SignalSpace(py::object input, std::string var_name,
                int dtype, int n_det, int n_time);
    SignalSpace(py::object input, std::string var_name,
                int dtype, int n_det, int n_time, int n_thirdaxis);
    ~SignalSpace() { if (data_ptr) free(data_ptr); };

    // Use arrays rather than vectors for anything that will be in
    // inner loop (even if it's inlined, the compiler might worry).
    DTYPE **data_ptr = NULL;
    int steps[64];
    vector<int> dims;

    vector<BufferWrapper<DTYPE>> bw;
    py::object ret_val;

private:
    constexpr bool _Validate(py::object input, std::string var_name,
                   int dtype);
};


// template<typename C, typename P, typename S>
template<typename CoordSys, typename PixelSys, typename SpinSys>
class ProjectionEngine {
public:
    //ProjectionEngine(PixelSys pixelizor);
    ProjectionEngine(py::object pix_args);
    py::object coords(py::object pbore, py::object pofs,
                      py::object coord);
    py::object pixels(py::object pbore, py::object pofs, py::object pixel);
    vector<int> tile_hits(py::object pbore, py::object pofs);
    py::object tile_ranges(py::object pbore, py::object pofs, py::list tile_lists);
    py::object pointing_matrix(py::object pbore, py::object pofs, py::object response,
                               py::object pixel, py::object proj);
    py::object zeros(py::object shape);
    py::object pixel_ranges(py::object pbore, py::object pofs, py::object map, int n_domain=-1);
    py::object from_map(py::object map, py::object pbore, py::object pofs,
                        py::object response, py::object signal);
    py::object to_map(py::object map, py::object pbore, py::object pofs, py::object response,
                      py::object signal, py::object det_weights,
                      py::object thread_intervals);
    py::object to_weight_map(py::object map, py::object pbore, py::object pofs,
                  py::object response, py::object det_weights, py::object thread_intervals);

    int comp_count() const;
    int index_count() const;

private:
    PixelSys _pixelizor;
};


void register_projection(py::module_ & m);
