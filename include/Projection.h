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
    
class CoordQuatCyl;
class CoordQuatZen;
class CoordFlat;

template <typename CoordSys>
class Pointer : public ProjectionOptimizer {
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
              double ix0=0., double iy0=0.);
    ~Pixelizor2_Flat() {};
    bool TestInputs(bp::object &map, bp::object &pbore, bp::object &pdet,
                    bp::object &signal, bp::object &weight);
    bp::object zeros(int count);
    int GetPixel(int i_det, int i_time, const double *coords);
    std::pair<int,int> IndexRange();
private:
    int crpix[2];
    double cdelt[2];
    int naxis[2];
    int strides[2];
};


/** Accumulator class templates.
 *
 * The SpinClass is used primarily to activate appropriate PixelWeight
 * specializations depending on the map components.  For example,
 * SpinT implies weight = 1, while SpinTQU implies weight triplet (1,
 * cos(2 phi), sin(2 phi)).  SpinClass also provides the number of
 * basic map components.
 */

template <int N>
class SpinClass {
public:
    static const int comp_count = N;
};

class SpinT : public SpinClass<1> {};
class SpinQU : public SpinClass<2> {};
class SpinTQU : public SpinClass<3> {};

template <typename SpinClass>
class Accumulator : public ProjectionOptimizer {
public:
    Accumulator<SpinClass>() {};
    Accumulator<SpinClass>(bool _need_map, bool _need_signal, bool _need_weight_map):
        need_map(_need_map), need_signal(_need_signal),
        need_weight_map(_need_weight_map) {};
    ~Accumulator<SpinClass>() {};
    inline int ComponentCount() {return SpinClass::comp_count;}
    void PixelWeight(const double *coords, double *wt);
    bool TestInputs(bp::object &map, bp::object &pbore, bp::object &pdet,
                    bp::object &signal, bp::object &weight);
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
protected:
    bool need_map = true;
    bool need_signal = true;
    bool need_weight_map = false;
    BufferWrapper _mapbuf;
    BufferWrapper _signalbuf;
};

template<typename P, typename Z, typename A>
class ProjectionEngine {
public:
    ProjectionEngine(Z pixelizor);
    bp::object to_map(bp::object map, bp::object pbore, bp::object pofs,
                      bp::object signal, bp::object weights);
    bp::object to_map_omp(bp::object map, bp::object pbore, bp::object pofs,
                          bp::object signal, bp::object weights,
                          bp::object thread_intervals);
    bp::object to_weight_map(bp::object map, bp::object pbore, bp::object pofs,
                             bp::object signal, bp::object weights);
    bp::object to_weight_map_omp(bp::object map, bp::object pbore, bp::object pofs,
                                 bp::object signal, bp::object weights,
                                 bp::object thread_intervals);
    bp::object from_map(bp::object map, bp::object pbore, bp::object pofs,
                        bp::object signal, bp::object weights);
    bp::object coords(bp::object pbore, bp::object pofs,
                      bp::object coord);
    bp::object pixels(bp::object pbore, bp::object pofs, bp::object pixel);
    bp::object pixel_ranges(bp::object pbore, bp::object pofs);
private:
    Z _pixelizor;
};
