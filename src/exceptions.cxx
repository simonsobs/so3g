
#include "exceptions.h"

namespace py = pybind11;


void register_exceptions(py::module_ & m) {
    py::exception<value_exception>(m, "SO3G_value_exception", PyExc_ValueError);
    py::exception<buffer_exception>(m, "SO3G_buffer_exception", PyExc_RuntimeError);
    py::exception<shape_exception>(m, "SO3G_shape_exception", PyExc_RuntimeError);
    py::exception<dtype_exception>(m, "SO3G_dtype_exception", PyExc_ValueError);
    py::exception<agreement_exception>(
        m, "SO3G_agreement_exception", PyExc_RuntimeError
    );
    py::exception<tiling_exception>(m, "SO3G_tiling_exception", PyExc_RuntimeError);
    py::exception<general_agreement_exception>(
        m, "SO3G_general_agreement_exception", PyExc_ValueError
    );
    py::exception<alloc_exception>(m, "SO3G_alloc_exception", PyExc_RuntimeError);
}
