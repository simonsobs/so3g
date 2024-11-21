#define NO_IMPORT_ARRAY

#include <pybindings.h>
#include <G3Ndmap.h>
#include <boost/python.hpp>
#include <container_pybindings.h>

G3Ndmap::G3Ndmap() {}
G3Ndmap::G3Ndmap(const G3Ndmap & other):data(other.data),wcs(other.wcs) {}
G3Ndmap::G3Ndmap(const G3Ndarray & data_, const G3WCS & wcs_):data(data_), wcs(wcs_) {}
G3Ndmap::G3Ndmap(const bp::object & barray, const string & header):data(barray), wcs(header) {}

string G3Ndmap::Description() const
{
    return "G3Ndmap("+data.Description()+", " + wcs.Description() + ")";
}

template <class A> void G3Ndmap::serialize(A &ar, unsigned v)
{
    using namespace cereal;
    ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));
    ar & make_nvp("data", data);
    ar & make_nvp("wcs", wcs);
}

using namespace boost::python;

G3_SERIALIZABLE_CODE(G3Ndmap);

PYBINDINGS("so3g")
{
    EXPORT_FRAMEOBJECT(G3Ndmap, init<>(), "G3Ndmap default constructor")
    .def(init<const G3Ndarray &, const G3WCS &>("Construct G3Nmap from a G3Ndarray and G3WCS"))
    .def(init<const bp::object &, const string &>("Construct G3Nmap from a numpy array and a header string"))
    .def_readwrite("data", &G3Ndmap::data)
    .def_readwrite("wcs", &G3Ndmap::wcs)
    ;
}
