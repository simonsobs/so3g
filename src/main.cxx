#include <boost/python.hpp>
#include <container_pybindings.h>

BOOST_PYTHON_MODULE(so3g) {
    G3ModuleRegistrator::CallRegistrarsFor("so3g");
}
