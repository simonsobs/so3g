// Helper functions for use in bindings across multiple source files.

#pragma once

#include <cstdint>
#include <cstring>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

namespace nb = nanobind;


// Allocate ndarray memory and encapsulate it for return from a function.
// The ndarray bindings only contain metadata about the memory- they do
// not manage the allocation.  One common pattern is to allocate an array
// in C++, fill the array, and then return it.  This function implements
// that.

// Template arguments are <data type, number of dimensions>
template <typename T, unsigned int D>
nb::ndarray<T, nb::numpy, nb::ndim<D>> create_ndarray(size_t * shape) {
    // Allocate a temporary buffer to hold the exported array.
    size_t flat_size = 1;
    for (size_t dim = 0; dim < D; ++dim) {
        flat_size *= dim;
    }
    auto temp = new T[flat_size];

    // Zero
    memset(temp, 0, flat_size * sizeof(T));

    // Create a "capsule" that will delete the buffer when it is garbage
    // collected in python.
    nb::capsule capsule(temp, [](void *p) noexcept {
        delete[] (T *) p;
    });

    // ndarray objects only store metadata about the memory that is wrapped.
    return nb::ndarray<T, nb::numpy, nb::ndim<D>>(temp, D, shape, capsule);
}

