#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

// Some macro defines are necessary when working with numpy C api
// across multiple translation units (i.e. source files):
//
// - https://sourceforge.net/p/numpy/mailman/message/5700519/
// - https://stackoverflow.com/questions/38003707/trouble-with-numpy-c-api
//
// As a result, independently of the boost numpy stuff, we need to
// tell numpy to drop in a single copy of the initialized API function
// pointers.  All files using the API should define the
// PY_ARRAY_UNIQUE_SYMBOL variable to the same value.  Then we call
// import_array() in one of those source files (this one), and define
// NO_IMPORT_ARRAY in the others (both of these before including
// arrayobject.h).

#define PY_ARRAY_UNIQUE_SYMBOL Py_Array_API_SO3G
#include <numpy/arrayobject.h>

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
