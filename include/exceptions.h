#pragma once

#include <stdexcept>
#include <sstream>

#include <pybind11/pybind11.h>

namespace py = pybind11;


// Current C++ guidance seems to recommend a "wide" versus "deep"
// exception hierarchy.  It is also important to not store a string
// as a member of a custom exception class, since the associated
// copy constructor must be non-throwing (and copying the member
// string cannot guarantee that).  Instead, we use std::runtime_error
// with has a built-in string storage that meets that requirement.
// We derive all custom exceptions from runtime_error.

class value_exception : public std::runtime_error
{
public:
    value_exception(std::string text) : std::runtime_error(text) {}
};

class buffer_exception : public std::runtime_error
{
public:
    buffer_exception(std::string var_name) : std::runtime_error(
        std::string("Argument '") +
        var_name +
        std::string("' does not expose buffer protocol, ") +
        std::string("is not contiguous, or does not export a format.")
    ) {}
};

class shape_exception : public std::runtime_error
{
public:
    shape_exception(
        std::string var_name, std::string detail
    ) : std::runtime_error(
        std::string("Buffer '") +
        var_name +
        std::string("' has incompatible shape: ") +
        detail +
        std::string(".")
    ) {}
};

class dtype_exception : public std::runtime_error
{
public:
    dtype_exception(
        std::string var_name, std::string type_str
    ) : std::runtime_error(
        std::string("Expected buffer '") +
        var_name +
        std::string("' to contain items of type ") +
        type_str +
        std::string(".")
    ) {}
};

class agreement_exception : public std::runtime_error
{
public:
    agreement_exception(
        std::string var1, std::string var2, std::string prop
    ) : std::runtime_error(
        std::string("Expected buffers '") +
        var1 +
        std::string("' and '") +
        var2 +
        std::string("' to have the same ") +
        prop +
        std::string(".")
    ) {}
};

class tiling_exception : public std::runtime_error
{
public:
    tiling_exception(
        int tile_idx, std::string msg
    ) : std::runtime_error(
        std::string("Tiling problem (index ") +
        std::to_string(tile_idx) +
        std::string("): ") +
        msg
    ) {}
};

class alloc_exception : public std::runtime_error
{
public:
    alloc_exception(
        std::string msg
    ) : std::runtime_error(
        std::string("Failed allocation: ") +
        msg
    ) {}
};

class general_agreement_exception : public std::runtime_error
{
public:
    general_agreement_exception(std::string text) : std::runtime_error(text) {}
};


void register_exceptions(py::module_ & m);
