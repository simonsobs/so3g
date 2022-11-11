#include <boost/python.hpp>
#ifdef _OPENMP
# include <omp.h>
#endif // ifdef _OPENMP

// See this header file for discussion of numpy compilation issues.
#include "so3g_numpy.h"

#include <container_pybindings.h>

#include <string>

// Note _version.h is supposed to be auto-generated during build.  If
// that breaks at some point, you can replace it with a single
// definition:
//   #define SO3G_VERSION_STRING "unknown"
#include "_version.h"

namespace bp = boost::python;

const std::string version()
{
    return SO3G_VERSION_STRING;
}

bp::object useful_info() {
    int omp_num_threads = 1;
#pragma omp parallel
    {
        #ifdef _OPENMP
        if (omp_get_thread_num() == 0)
            omp_num_threads = omp_get_num_threads();
        #endif
    }
    bp::dict output;
    output["omp_num_threads"] = omp_num_threads;
    output["version"] = version();
    return output;
}




PYBINDINGS("so3g") {
    bp::def("version", version);
    bp::def("useful_info", useful_info);
}

static void* _so3g_import_array() {
    import_array();
    return NULL;
}

BOOST_PYTHON_MODULE(so3g) {
    _so3g_import_array();
    G3ModuleRegistrator::CallRegistrarsFor("so3g");
}
