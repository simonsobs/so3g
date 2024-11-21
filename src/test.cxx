#include <pybindings.h>

#include <iostream>
#include <test.h>

#include <boost/python.hpp>

#include <container_pybindings.h>

void TestClass::runme() {
    cout << "I guess it works." << endl;
}


std::string TestFrame::Description() const
{
	std::ostringstream s;
	s << "data_source '" << data_source << "' @session_id:" << session_id;
	return s.str();
}

std::string TestFrame::Summary() const
{
    return Description();
}

template <class A> void TestFrame::serialize(A &ar, unsigned v)
{
	using namespace cereal;
        // v is the version code!

	ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));
	ar & make_nvp("data_source", data_source);
	ar & make_nvp("session_id", session_id);
}

void greet() {
    cout << "test 1 complete." << endl;
}

namespace bp = boost::python;

PYBINDINGS("_libso3g")
{
    bp::def("greet", greet);
    bp::class_<TestClass>("TestClass")
        .def("runme", &TestClass::runme);

    EXPORT_FRAMEOBJECT(TestFrame, init<>(),
    "TestFrame for demonstration.")
    .def_readwrite("session_id", &TestFrame::session_id,
    "Example integer.")
    .def_readwrite("data_source", &TestFrame::data_source,
    "Example string.")
    ;

    //	register_g3map<DfMuxWiringMap>("DfMuxWiringMap", "Mapping from "
    //	   "logical detector ID string (same as used in timestreams) to wiring "
    //	   "information (the board, module, and channel to which a given "
    //	   "detector is connected)");
}
