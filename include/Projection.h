#include <boost/python.hpp>

namespace bp = boost::python;

class BufferWrapper;

/** ProjectionOptimizer is the base class for our optimization model.
 *
 *  Though what we get from this is quite unclear, since interface
 *  definitions don't seem to work well with boost::python in the mix.
 */

class ProjectionOptimizer {
public:
    ProjectionOptimizer() {};
    ~ProjectionOptimizer() {};
    bool TestInputs(bp::object &map, bp::object &pbore, bp::object &pdet,
                    bp::object &signal, bp::object &weight);
    inline
    void InitPerDet(int idet) {};
    inline
    void InitPerSample() {};
};


    
class Pointer : public ProjectionOptimizer {
public:
    //Pointer();
    //~Pointer() {};
    void Init(BufferWrapper &qborebuf, BufferWrapper &qofsbuf);
    void InitPerDet(int idet);
    void GetCoords(int i_det, int i_t, double *coords);
private:
    double _dx;
    double _dy;
    BufferWrapper *_qborebuf;
    BufferWrapper *_qofsbuf;
};


    
class Pixelizor : public ProjectionOptimizer {
public:
    Pixelizor();
    ~Pixelizor() {};
    //bool TestInputs(bp::object map, bp::object weight) {return false;};
    bp::object zeros(bp::object shape);
    void Init(BufferWrapper &mapbuf);
    int GetPixel(int i_det, int i_t, const double *coords);
private:
    int crpix[2];
    double crval[2];
    double cdelt[2];
    int naxis[2];
    BufferWrapper *_mapbuf;
};


/* Accumulator classes. */

class Accumulator : public ProjectionOptimizer {
public:
    Accumulator() {};
    ~Accumulator() {};
//    bool TestInputs(bp::object map, bp::object signal) {return false;};
    void Init(BufferWrapper &inline_weightbuf) {};
    void Forward(const BufferWrapper &inline_weightbuf,
                 const BufferWrapper &signalbuf,
                 BufferWrapper &mapbuf,
                 const int idet,
                 const int it,
                 const double* coords,
                 const int pixel_index) {};
};

class AccumulatorSpin0 : public Accumulator {
public:
    AccumulatorSpin0() {_handle = nullptr; };
    ~AccumulatorSpin0();
    bool TestInputs(bp::object map, bp::object signal, bp::object weight);
    void Init(BufferWrapper &inline_weightbuf);
    void Forward(const BufferWrapper &inline_weightbuf,
                 const BufferWrapper &signalbuf,
                 BufferWrapper &mapbuf,
                 const int idet,
                 const int it,
                 const double* coords,
                 const int pixel_index);
private:
    PyObject *_handle;
};


class AccumulatorSpin2 : public Accumulator {
public:
    AccumulatorSpin2() {_handle = nullptr; };
    ~AccumulatorSpin2();
    bool TestInputs(bp::object map, bp::object signal, bp::object weight);
    void Init(BufferWrapper &inline_weightbuf);
    void Forward(const BufferWrapper &inline_weightbuf,
                 const BufferWrapper &signalbuf,
                 BufferWrapper &mapbuf,
                 const int idet,
                 const int it,
                 const double* coords,
                 const int pixel_index);
private:
    PyObject *_handle;
};


template<typename P, typename Z, typename A>
class ProjectionEngine {
public:
    bp::object zeros(bp::object shape);
    bp::object to_map(bp::object map, bp::object qbore, bp::object qofs,
                      bp::object signal, bp::object weights);
};

