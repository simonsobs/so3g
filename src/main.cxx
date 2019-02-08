#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

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
namespace np = boost::python::numpy;

const std::string version()
{
    return SO3G_VERSION_STRING;
}

PYBINDINGS("so3g") {
    bp::def("version", version);
}

static void* _so3g_import_array() {
    import_array();
    return NULL;
}

BOOST_PYTHON_MODULE(so3g) {
    _so3g_import_array();
    np::initialize();
    G3ModuleRegistrator::CallRegistrarsFor("so3g");
}
