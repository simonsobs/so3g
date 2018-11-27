#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

#include <container_pybindings.h>

namespace np = boost::python::numpy;

BOOST_PYTHON_MODULE(so3g) {
    np::initialize();
    G3ModuleRegistrator::CallRegistrarsFor("so3g");
}
