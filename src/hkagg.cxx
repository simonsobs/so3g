
#include <iostream>

#include <hkagg.h>
#include "exceptions.h"

namespace py = pybind11;

int hk_frame_type_int(HKFrameType typ) {
    // Manually convert enum to int for python.
    if (typ == HKFrameType::session) {
        return 0;
    } else if (typ == HKFrameType::status) {
        return 1;
    } else if (typ == HKFrameType::data) {
        return 2;
    } else {
        throw RuntimeError_exception("Invalid HKFrameType enum value");
    }
    return -1;
}

void register_hkagg(py::module_ & m) {
    py::enum_<HKFrameType>(m, "HKFrameType", "Identifier for generic HK streams.")
    .value("session", HKFrameType::session)
    .value("status", HKFrameType::status)
    .value("data", HKFrameType::data);

    m.def("hk_frame_type_int", &hk_frame_type_int);
}
