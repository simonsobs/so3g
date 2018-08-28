#include <pybindings.h>

#include <iostream>
#include <SampleFlags.h>

#include <boost/python.hpp>

#include <container_pybindings.h>

std::string SampleFlagSet::Description() const
{
	std::ostringstream s;
	s << "SampleFlagSet " << name << " for " << flags.size() << " flags.";
	return s.str();
}

std::string SampleFlagSet::Summary() const
{
    return Description();
}

template <class A> void SampleFlagSet::serialize(A &ar, unsigned v)
{
	using namespace cereal;
        // v is the version code!

	ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));
	ar & make_nvp("name", name);
}

using namespace boost::python;

PYBINDINGS("so3g")
{
   EXPORT_FRAMEOBJECT(SampleFlagSet, init<>(),
   "Per-sample flag set, stored efficiently.")
    // register_g3map<SampleFlagSet>("SampleFlagSet", "Mapping from "
    //                                "strings to arrays of integers.")
    .def_readwrite("name", &SampleFlagSet::name,
                   "IP Address of the board, encoded as an int using struct")
    .def_readwrite("flags", &SampleFlagSet::flags,
                   "IP Address of the board, encoded as an int using struct")
    ;
}

// BOOST_PYTHON_MODULE(so3g) {
//     //    def("greet", greet);
//     G3ModuleRegistrator::CallRegistrarsFor("so3g");
// }
