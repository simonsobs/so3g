#include <boost/python/exception_translator.hpp>

#include <container_pybindings.h>
#include "exceptions.h"

// Here we define the map from C++ exceptions defined in exceptions.h
// to Python exceptions.  Currently we use built-in python exceptions,
// with informative error messages.

static void translate_RuntimeError(so3g_exception const& e)
{
    PyErr_SetString(PyExc_RuntimeError, e.msg_for_python().c_str());
}

static void translate_TypeError(so3g_exception const& e)
{
    PyErr_SetString(PyExc_TypeError, e.msg_for_python().c_str());
}

static void translate_ValueError(so3g_exception const& e)
{
    PyErr_SetString(PyExc_ValueError, e.msg_for_python().c_str());
}

static void translate_KeyError(so3g_exception const& e)
{
    PyErr_SetString(PyExc_KeyError, e.msg_for_python().c_str());
}


namespace bp = boost::python;
PYBINDINGS("so3g")
{
    bp::register_exception_translator<buffer_exception>    (&translate_TypeError);
    bp::register_exception_translator<dtype_exception>     (&translate_ValueError);
    bp::register_exception_translator<shape_exception>     (&translate_RuntimeError);
    bp::register_exception_translator<agreement_exception> (&translate_RuntimeError);
    bp::register_exception_translator<general_agreement_exception> (&translate_ValueError);
    bp::register_exception_translator<key_crawling_exception>(&translate_KeyError);
    bp::register_exception_translator<value_decoding_exception>(&translate_ValueError);
}
