#include <vector>
#include <cstdint>

#include <boost/python.hpp>

namespace bp = boost::python;

class BufferWrapper;

/* To-do -- can we make this a properly abstract base class?  Initial
 * efforts on this lead to run-time symbol errors, because it somehow
 * is still trying to find Weightor even if only Spin0Weightor is in
 * use, through a python export.  */

class Weightor {
public:
    Weightor() {};
    ~Weightor() {};
    bool TestInputs(bp::object map, bp::object weight) {return false;};
    bool GetWeights(BufferWrapper &inline_weightbuf, double x, double y, double phi)
        { return false; };
};

class Spin0Weightor : public Weightor {
public:
    Spin0Weightor() {_handle = nullptr; };
    ~Spin0Weightor();
    bool TestInputs(bp::object map, bp::object weight);
    bool GetWeights(BufferWrapper &inline_weightbuf, double x, double y, double phi);
private:
    PyObject *_handle;
};

template<typename W>
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
    bp::object to_map2(bp::object map, bp::object qbore, bp::object qofs,
                      bp::object signal, bp::object weights);
};

