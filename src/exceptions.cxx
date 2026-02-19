
#include "exceptions.h"

namespace py = pybind11;


void register_exceptions(py::module_ & m) {
    py::exception<RuntimeError_exception>(m, "SO3G_RuntimeError", PyExc_RuntimeError);
    py::exception<TypeError_exception>(m, "SO3G_TypeError", PyExc_TypeError);
    py::exception<ValueError_exception>(m, "SO3G_ValueError", PyExc_ValueError);
}
