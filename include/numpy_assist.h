#pragma once

#include <exception>
#include <memory>
#include <vector>
#include <iostream>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "exceptions.h"

namespace py = pybind11;


// check_buffer_type<T>(const Py_buffer &view)
//
// This function checks whether the view is consistent with the
// requested type, by checking view.format and view.itemsize.  Returns
// true if so, and false otherwise.

template <typename T>
static bool _check_buffer_helper(const Py_buffer &view, std::string opts) {
    bool code_ok = false;
    if (view.format == NULL || *view.format == 0)
        return false;
    for (auto c: opts)
        code_ok = code_ok || (view.format[0] == c);
    return code_ok && view.itemsize == sizeof(T);
}

template <typename T,
          typename std::enable_if<std::is_integral<T>::value,
                                  int>::type* = nullptr>
static bool check_buffer_type(const Py_buffer &view) {
    return _check_buffer_helper<T>(view, "bhilq");
}

template <typename T,
          typename std::enable_if<!std::is_integral<T>::value,
                                  int>::type* = nullptr>
static bool check_buffer_type(const Py_buffer &view) {
    return false;
}

template <>
bool check_buffer_type<float>(const Py_buffer &view) {
    return _check_buffer_helper<float>(view, "f");
}

template <>
bool check_buffer_type<double>(const Py_buffer &view) {
    return _check_buffer_helper<double>(view, "d");
}


// type_name<T>()
//
// This function returns a string description of type T, for reporting
// in exceptions.  The strings here should probably look like
// numpy.dtype names; i.e. "int32" not "int" or "int32_t".

template <typename T>
static std::string type_name() {
    return "unknown";
}

template <>
std::string type_name<int32_t>() {
    return "int32";
}

template <>
std::string type_name<int64_t>() {
    return "int64";
}

template <>
std::string type_name<float>() {
    return "float32";
}

template <>
std::string type_name<double>() {
    return "float64";
}


static std::string shape_string(std::vector<int> shape)
{
    std::ostringstream s;
    s << "(";
    for (int i=0; i<shape.size(); i++) {
        if (i > 0)
            s << ", ";
        if (shape[i] >= 0)
            s << shape[i];
        else if (shape[i] == -1)
            s << "*";
        else if (shape[i] == -2)
            s << "...";
        else
            s << "!error";
    }
    s << ")";
    return s.str();
}


// class BufferWrapper
//
// A wrapper for Py_buffer pointer, templated for a particular data
// type.  The underlying Py_buffer is privately held, and
// reference-counted with a shared_ptr so it is freed when all copies
// of the parent of gone out of scope.  Various constructors help to
// eliminate boilerplate when wrapping objects that must be of a
// certain shape / type.

template <typename T>
class BufferWrapper {
public:
    // Through the -> operator you can access the fields of the
    // Py_buffer.
    Py_buffer *operator->() const {
        return view.get();
    }

    BufferWrapper() {
        auto p = (Py_buffer*)calloc(1, sizeof(Py_buffer));
        view = std::shared_ptr<Py_buffer>(p, PyBuffer_Release);
    }

    // Constructor with no shape or type checking.
    BufferWrapper(std::string name, const py::object &src, bool optional)
        : BufferWrapper() {
        if (optional && (src.is_none()))
            return;
        //std::cerr << "BufferWrapper ctor " << name << " src = " << src.ptr() << ")" << std::endl;
        if (PyObject_GetBuffer(src.ptr(), view.get(),
                               PyBUF_RECORDS) == -1) {
            PyErr_Clear();
            throw buffer_exception(name);
        }
        // std::cerr << "BufferWrapper ctor " << name << " view now = " << view.get() << std::endl;
        //int ndim = view.get()->ndim;
        // std::cerr << "BufferWrapper ctor " << name << " view ndim = " << ndim << std::endl;
        // std::cerr << "BufferWrapper ctor " << name << " (";
        // for (int i = 0; i < ndim; ++i) {
        //     std::cerr << view.get()->shape[i] << ",";
        // }
        //std::cerr << ")" << std::endl;
    }

    // Constructor with shape and type checking.
    BufferWrapper(std::string name, const py::object &src, bool optional,
                  std::vector<int> shape)
        : BufferWrapper(name, src, optional) {

        std::cerr << "BufferWrapper start 4-arg ctor" << std::endl;

        // "optional" items will cause the parent constructor to
        // succeed, but will leave buffer pointer unset.
        if (view->buf == NULL) {
            std::cerr << "BufferWrapper: view->buf == NULL, return" << std::endl;
            return;
        }

        if (!check_buffer_type<T>(*(view.get()))) {
            std::cerr << "BufferWrapper: view invalid type = " << type_name<T>() << std::endl;
            throw dtype_exception(name, type_name<T>());
        }

        std::vector<int> vshape;
        std::cerr << "BufferWrapper: vshape = ";
        for (int i=0; i<view->ndim; i++) {
            std::cerr << view->shape[i] << ",";
            vshape.push_back(view->shape[i]);
        }
        std::cerr << std::endl;

        // Note special value -1 is as in numpy -- matches a single
        // axis.  Special value -2 is treated as an ellipsis -- can be
        // specified at most once, and matches 0 or more axes.
        int i=0, j=0;
        int ellipsis_count = 0;
        while (i < shape.size()) {
            if (shape[i] == -2) {
                if (ellipsis_count++) {
                    std::ostringstream s;
                    s << "Invalid shape specifier " << shape_string(shape) << " (multiple elipses).";
                    throw shape_exception(name, s.str());
                }
                // Ignore 0 or more leading entries.
                j = vshape.size() - (shape.size() - i) + 1;
            } else if (j >= vshape.size()) {
                break;
            } else if (shape[i] == -1 || shape[i] == vshape[j]) {
                // Match.
                j++;
            } else
                break;
            i++;
        }
        if (i != shape.size() || j != vshape.size()) {
            std::ostringstream s;
            s << "Expected " << shape_string(shape) << " but got " <<
                shape_string(vshape) << ".";
            std::cerr << s.str() << std::endl;
            throw shape_exception(name, s.str());
        }
    }

    // These accessors are for convenience but you should probably
    // roll your own for inner loops.
    inline
    T* ptr_1d(int i0) {
        return (T*)((char*)view->buf + view->strides[0] * i0);
    }

    inline
    T* ptr_2d(int i0, int i1) {
        return (T*)((char*)view->buf + view->strides[0] * i0 + view->strides[1] * i1);
    }

    inline
    bool test() {
        return view->obj != nullptr;
    }

private:
    std::shared_ptr<Py_buffer> view;
};


// Convert an n-dimensional array into a list of array slices.
py::list list_of_arrays(py::object input);
