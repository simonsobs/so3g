
#include "exceptions.h"

namespace nb = nanobind;


void register_exceptions(nb::module_ & m) {
    nb::exception<RuntimeError_exception>(m, "SO3G_RuntimeError", PyExc_RuntimeError);
    nb::exception<TypeError_exception>(m, "SO3G_TypeError", PyExc_TypeError);
    nb::exception<ValueError_exception>(m, "SO3G_ValueError", PyExc_ValueError);
}
