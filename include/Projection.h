#include <vector>
#include <cstdint>

#include <boost/python.hpp>

namespace bp = boost::python;

class GnomonicGridder {
public:
    int crpix[2];
    double crval[2];
    double cdelt[2];
    int naxis[2];

    GnomonicGridder();
    bp::object zeros(bp::object shape);
    bp::object to_map(bp::object map, bp::object qpoint,
                      bp::object signal, bp::object weights);
};
