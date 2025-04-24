#pragma once

#include <stdint.h>

#include <nanobind/nanobind.h>

using namespace std;

namespace nb = nanobind;


enum HKFrameType {
     session = 0,
     status = 1,
     data = 2,
};

void register_hkagg(nb::module_ & m);
