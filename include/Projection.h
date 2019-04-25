#include <boost/python.hpp>
#include "exceptions.h"

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
    virtual bool TestInputs(bp::object &map, bp::object &pbore, bp::object &pdet,
                            bp::object &signal, bp::object &weight) = 0;
    inline
    void InitPerDet(int i_det) {};
    inline
    void InitPerSample() {};
};
    
class Pointer : public ProjectionOptimizer {
public:
    bool TestInputs(bp::object &map, bp::object &pbore, bp::object &pdet,
                    bp::object &signal, bp::object &weight) { return false; };
    void InitPerDet(int i_det);
    void GetCoords(int i_det, int i_time, double *coords);
};

class PointerP_Simple_Flat : public Pointer {
public:
    bool TestInputs(bp::object &map, bp::object &pbore, bp::object &pdet,
                    bp::object &signal, bp::object &weight);
    void InitPerDet(int i_det);
    int DetCount() { return n_det; }
    int TimeCount() { return n_time; }
    void GetCoords(int i_det, int i_timeime, double *coords);
private:
    double _coords[4];
    BufferWrapper _pborebuf;
    BufferWrapper _pdetbuf;
    int n_det;
    int n_time;
};
    
class Pixelizor2_Flat : public ProjectionOptimizer {
public:
    Pixelizor2_Flat() {};
    Pixelizor2_Flat(int nx, int ny,
              double dx=1., double dy=1.,
              double x0=0., double y0=0.,
              double ix0=0., double iy0=0.);
    ~Pixelizor2_Flat() {};
    bool TestInputs(bp::object &map, bp::object &pbore, bp::object &pdet,
                    bp::object &signal, bp::object &weight);
    bp::object zeros(int count);
    int GetPixel(int i_det, int i_time, const double *coords);
private:
    int crpix[2];
    double crval[2];
    double cdelt[2];
    int naxis[2];
    int strides[2];
};


/* Accumulator classes. */

class Accumulator : public ProjectionOptimizer {
public:
    Accumulator() {};
    Accumulator(bool _need_map, bool _need_signal, bool _need_weight_map):
        need_map(_need_map), need_signal(_need_signal),
        need_weight_map(_need_weight_map) {};
    ~Accumulator() {};
    bool TestInputs(bp::object &map, bp::object &pbore, bp::object &pdet,
                    bp::object &signal, bp::object &weight) { return false; };
    int ComponentCount() {return 0;}
    void Forward(const int i_det,
                 const int i_time,
                 const int pixel_index,
                 const double* coords,
                 const double* weights) {};
    void Reverse(const int i_det,
                 const int i_time,
                 const int pixel_index,
                 const double* coords,
                 const double* weights) {};
protected:
    bool need_map = true;
    bool need_signal = true;
    bool need_weight_map = false;
};

class AccumulatorT_Flat : public Accumulator {
public:
    AccumulatorT_Flat(bool _need_map, bool _need_signal, bool _need_weight_map):
        Accumulator(_need_map, _need_signal, _need_weight_map) {};
    bool TestInputs(bp::object &map, bp::object &pbore, bp::object &pdet,
                    bp::object &signal, bp::object &weight);
    int ComponentCount() {return 1;}
    void Forward(const int i_det,
                 const int i_time,
                 const int pixel_index,
                 const double* coords,
                 const double* weights);
    void ForwardWeight(const int i_det,
                       const int i_time,
                       const int pixel_index,
                       const double* coords,
                       const double* weights);
    void Reverse(const int i_det,
                 const int i_time,
                 const int pixel_index,
                 const double* coords,
                 const double* weights);
private:
    BufferWrapper _mapbuf;
    BufferWrapper _signalbuf;
};

class AccumulatorTQU_Flat : public Accumulator {
public:
    AccumulatorTQU_Flat(bool _need_map, bool _need_signal, bool _need_weight_map):
        Accumulator(_need_map, _need_signal, _need_weight_map) {};
    bool TestInputs(bp::object &map, bp::object &pbore, bp::object &pdet,
                    bp::object &signal, bp::object &weight);
    int ComponentCount() {return 3;}
    void Forward(const int i_det,
                 const int i_time,
                 const int pixel_index,
                 const double* coords,
                 const double* weights);
    void ForwardWeight(const int i_det,
                       const int i_time,
                       const int pixel_index,
                       const double* coords,
                       const double* weights);
    void Reverse(const int i_det,
                 const int i_time,
                 const int pixel_index,
                 const double* coords,
                 const double* weights);
private:
    BufferWrapper _mapbuf;
    BufferWrapper _signalbuf;
};


template<typename P, typename Z, typename A>
class ProjectionEngine {
public:
    ProjectionEngine(Z pixelizor);
    bp::object to_map(bp::object map, bp::object pbore, bp::object pofs,
                      bp::object signal, bp::object weights);
    bp::object to_weight_map(bp::object map, bp::object pbore, bp::object pofs,
                             bp::object signal, bp::object weights);
    bp::object from_map(bp::object map, bp::object pbore, bp::object pofs,
                        bp::object signal, bp::object weights);
    bp::object coords(bp::object pbore, bp::object pofs,
                      bp::object coord);
    bp::object pixels(bp::object pbore, bp::object pofs, bp::object pixel);
private:
    Z _pixelizor;
};
