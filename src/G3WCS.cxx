#include <pybindings.h>
#include <G3WCS.h>
#include <boost/python.hpp>
#include <container_pybindings.h>

G3WCS::G3WCS() {}
G3WCS::G3WCS(const G3WCS & other):header(other.header) {}
G3WCS::G3WCS(const string & hdr):header(hdr) {}

string G3WCS::Description() const
{
    return "G3WCS("+header+")";
}

template <class A> void G3WCS::serialize(A &ar, unsigned v)
{
    using namespace cereal;
    ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));
    ar & make_nvp("header", header);
}

using namespace boost::python;

G3_SERIALIZABLE_CODE(G3WCS);

PYBINDINGS("so3g")
{
    EXPORT_FRAMEOBJECT(G3WCS, init<>(), "G3WCS default constructor")
    .def(init<const string &>("Construct G3Ndarray from string"))
    .def_readwrite("header", &G3WCS::header)
    ;
}
