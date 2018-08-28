#include <pybindings.h>

#include <iostream>
#include <test.h>

#include <boost/python.hpp>

#include <container_pybindings.h>

void TestClass::runme() {
    cout << "I guess it works." << endl;
}



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

//G3_SERIALIZABLE_CODE(HKInfo);

void greet() {
    cout << "test 1 complete." << endl;
}

using namespace boost::python;

//BOOST_PYTHON_MODULE(libso3g)

PYBINDINGS("so3g")
{
    class_<TestClass>("TestClass")
        .def("runme", &TestClass::runme);

    EXPORT_FRAMEOBJECT(HKInfo, init<>(),
    "Bolometer wiring information. Module and channel IDs are stored "
    "zero-indexed, but be aware that they often printed one-indexed "
    "for compatibility with pydfmux.")
    .def_readwrite("session_id", &HKInfo::session_id,
    "IP Address of the board, encoded as an int using struct")
    .def_readwrite("hk_source", &HKInfo::hk_source,
    "Serial number of the readout board to which this channel is "
    "attached.")
    ;

    //	register_g3map<DfMuxWiringMap>("DfMuxWiringMap", "Mapping from "
    //	   "logical detector ID string (same as used in timestreams) to wiring "
    //	   "information (the board, module, and channel to which a given "
    //	   "detector is connected)");
}

BOOST_PYTHON_MODULE(so3g) {
    def("greet", greet);
    G3ModuleRegistrator::CallRegistrarsFor("so3g");
}
