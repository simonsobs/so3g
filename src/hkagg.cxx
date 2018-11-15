#include <pybindings.h>

#include <iostream>
#include <boost/python.hpp>

#include <container_pybindings.h>
#include <test.h>

std::string HKInfo::Description() const
{
	std::ostringstream s;
	s << "hk_source '" << hk_source << "' @session_id:" << session_id;
	return s.str();
}

std::string HKInfo::Summary() const
{
    return Description();
}

template <class A> void HKInfo::serialize(A &ar, unsigned v)
{
	using namespace cereal;
        // v is the version code!

	ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));
	ar & make_nvp("hk_source", hk_source);
	ar & make_nvp("session_id", session_id);
}

namespace bp = boost::python;

PYBINDINGS("so3g")
{
    EXPORT_FRAMEOBJECT(HKInfo, init<>(),
    "Housekeeping Info Frame.  Carries session description information.  "
    "Subsequent Housekeeping Data Frames carry a matching session_id.")
    .def_readwrite("session_id", &HKInfo::session_id,
    "Identifier associated with aggregation software session.")
    .def_readwrite("hk_source", &HKInfo::hk_source,
    "Source name, or something.")
    ;
}
