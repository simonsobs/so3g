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


// template<typename C, typename P, typename S>
template<typename CoordSys, typename PixelSys, typename SpinSys>
class ProjectionEngine {
public:
    //ProjectionEngine(PixelSys pixelizor);
    ProjectionEngine(bp::object pix_args);
    bp::object coords(bp::object pbore, bp::object pofs,
                      bp::object coord);
    bp::object pixels(bp::object pbore, bp::object pofs, bp::object pixel);
    bp::object pointing_matrix(bp::object pbore, bp::object pofs,
                               bp::object pixel, bp::object proj);
    bp::object zeros(bp::object shape);
    bp::object from_map(bp::object map, bp::object pbore, bp::object pofs,
                        bp::object signal, bp::object weights);
    bp::object to_map(bp::object map, bp::object pbore, bp::object pofs,
                      bp::object signal, bp::object weights);
    bp::object to_weight_map(bp::object map, bp::object pbore, bp::object pofs,
                             bp::object signal, bp::object weights);
    // bp::object to_map_omp(bp::object map, bp::object pbore, bp::object pofs,
    //                       bp::object signal, bp::object weights,
    //                       bp::object thread_intervals);
    // bp::object to_weight_map_omp(bp::object map, bp::object pbore, bp::object pofs,
    //                              bp::object signal, bp::object weights,
    //                              bp::object thread_intervals);
    // bp::object pixel_ranges(bp::object pbore, bp::object pofs);

    int comp_count() const;
    int index_count() const;

private:
    PixelSys _pixelizor;
};
