
#include <string>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#ifdef _OPENMP
# include <omp.h>
#endif // ifdef _OPENMP

#include <hkagg.h>
#include <exceptions.h>
#include <Butterworth.h>
#include <so_linterp.h>
//#include <Intervals.h>

// Note _version.h is auto-generated during build.  You can
// also override it here if needed:
//   #define SO3G_VERSION_STRING "unknown"
#include "_version.h"

namespace nb = nanobind;


const std::string version() {
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


NB_MODULE(libso3g, m) {
    m.def("version", &version);
    m.def("useful_info", &useful_info);

    register_hkagg(m);
    register_exeptions(m);
    register_butterworth(m);
    register_so_linterp(m);
    //register_intervals(m);
}
