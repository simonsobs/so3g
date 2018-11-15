#include <pybindings.h>

#include <iostream>
#include <boost/python.hpp>

#include <container_pybindings.h>
#include <hkagg.h>


/* IrregBlockDouble */

std::string IrregBlockDouble::Description() const
{
	std::ostringstream s;
	s << "Double data (" << data.size() << " vectors) with timestamp.";
	return s.str();
}

std::string IrregBlockDouble::Summary() const
{
    return Description();
}

template <class A> void IrregBlockDouble::serialize(A &ar, unsigned v)
{
	using namespace cereal;
        // v is the version code!

	ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));
	ar & make_nvp("prefix", prefix);
	ar & make_nvp("t", t);
	ar & make_nvp("data", data);
}


G3_SERIALIZABLE_CODE(IrregBlockDouble);


namespace bp = boost::python;

PYBINDINGS("so3g")
{
    EXPORT_FRAMEOBJECT(IrregBlockDouble, init<>(),
    "Data block for irregularly sampled data.")
    .def_readwrite("prefix", &IrregBlockDouble::prefix,
    "Prefix for field names.")
    .def_readwrite("data", &IrregBlockDouble::data,
    "Map to HK data vectors.")
    .def_readwrite("t", &IrregBlockDouble::t,
    "Timestamp vector.")
    ;

    bp::enum_<HKFrameType>("HKFrameType",
                           "Identifier for generic HK streams.")
        .value("session", HKFrameType::session)
        .value("status",  HKFrameType::status)
        .value("data",    HKFrameType::data)
        ;

}
