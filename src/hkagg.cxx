
#include <iostream>

#include <nanobind/nanobind.h>

#include <hkagg.h>


namespace nb = nanobind;

NB_MODULE(libso3g, m) {
    nb::enum_<HKFrameType>()
    .value("session", HKFrameType::session)
    .value("status", HKFrameType::status)
    .value("data", HKFrameType::data)
    .export_values();
}

