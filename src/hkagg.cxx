
#include <iostream>

#include <hkagg.h>

namespace nb = nanobind;


void register_hkagg(nb::module_ & m) {
    nb::enum_<HKFrameType>(m, "HKFrameType")
    .value("session", HKFrameType::session)
    .value("status", HKFrameType::status)
    .value("data", HKFrameType::data)
    .export_values();
}

