#ifdef _OPENMP
# include <omp.h>
#endif // ifdef _OPENMP

#include <string>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

// See this header file for discussion of numpy compilation issues.
#include "so3g_numpy.h"

// Note _version.h is supposed to be auto-generated during build.  If
// that breaks at some point, you can replace it with a single
// definition:
//   #define SO3G_VERSION_STRING "unknown"
#include "_version.h"

// Include headers with registration functions
#include "exceptions.h"
#include "hkagg.h"
#include "so_linterp.h"
#include "Butterworth.h"
#include "Intervals.h"
#include "Ranges.h"
#include "Projection.h"

// Declaration here, since there is no header file for array_ops.
void register_array_ops(nb::module_ & m);


namespace nb = nanobind;


const std::string version()
{
    return SO3G_VERSION_STRING;
}

nb::object useful_info() {
    int omp_num_threads = 1;
#pragma omp parallel
    {
        #ifdef _OPENMP
        if (omp_get_thread_num() == 0)
            omp_num_threads = omp_get_num_threads();
        #endif
    }
    nb::dict output;
    output["omp_num_threads"] = omp_num_threads;
    output["version"] = version();
    return output;
}


static void* _so3g_import_array() {
    import_array();
    return NULL;
}


NB_MODULE(libso3g, m) {
    _so3g_import_array();

    m.def("version", &version);
    m.def("useful_info", &useful_info);

    register_exceptions(m);
    register_hkagg(m);
    register_butterworth(m);
    register_so_linterp(m);
    register_intervals(m);
    register_ranges(m);
    register_array_ops(m);
    register_projection(m);
}
