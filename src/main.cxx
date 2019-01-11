#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

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

BOOST_PYTHON_MODULE(so3g) {
    np::initialize();
    G3ModuleRegistrator::CallRegistrarsFor("so3g");
}
