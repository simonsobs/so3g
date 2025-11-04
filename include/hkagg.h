#pragma once

#include <stdint.h>

#include <pybind11/pybind11.h>
#include <pybind11/native_enum.h>

namespace py = pybind11;

using namespace std;


enum HKFrameType {
     session = 0,
     status = 1,
     data = 2,
};


void register_hkagg(py::module_ & m);
