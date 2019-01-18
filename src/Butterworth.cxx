#include <assert.h>
#include <math.h>

#include <pybindings.h>
#include <container_pybindings.h>

#include "Butterworth.h"

BFilterBank::BFilterBank(const BFilterBank& a) {
    // Copy the parameters but reset the accumulators... that's probably evil.
    for (auto p: a.par)
        par.push_back(p);
    if (a.w.size() > 0)
        init(a.w[0].size());
}


BFilterBank& BFilterBank::add(BFilterParams bp) {
    par.push_back(bp);
    return *this;
}

BFilterBank& BFilterBank::init(int n_chan) {
    assert(n_chan == 1); // Implementation assumes shape = (n_samp).
    // Populate accumulators based on filter bank.
    w.clear();
    for (auto bp: par) {
        auto v = vector<array<int64_t,2>>(n_chan);
        w.push_back(v);
    }
    return *this;
}

void BFilterBank::apply_numpy(np::ndarray& input, np::ndarray& output)
{
    // Record the type and shape, then check that input and output
    // match and are acceptable.
    np::dtype dtype = input.get_dtype();
    int n_dims = input.get_nd();
    int n_samp = input.shape(0);

    assert(output.get_flags() & np::ndarray::C_CONTIGUOUS);
    assert(input.get_flags() & np::ndarray::C_CONTIGUOUS);
    
    assert(output.get_nd() == n_dims);
    assert(output.shape(0) == n_samp);
    assert(output.get_dtype() == dtype);
    
    if (dtype == np::dtype::get_builtin<int32_t>()) {
        int32_t *in = reinterpret_cast<int32_t*>(input.get_data());
        int32_t *out = reinterpret_cast<int32_t*>(output.get_data());
        apply(in, out, n_samp);
    } else if (dtype == np::dtype::get_builtin<float>()) {
        float *in = reinterpret_cast<float*>(input.get_data());
        float *out = reinterpret_cast<float*>(output.get_data());
        apply_to_float(in, out, 1., n_samp);
    } else {
        assert(0); // dtype unknown.
    }
}


void BFilterBank::apply(int32_t* input, int32_t* output, int n_samp)
{
    int n_filt = w.size();
    assert(n_filt == par.size());
    int n_chan = w[0].size();

#pragma omp parallel for shared(input, output)
    for (int ic=0; ic<n_chan; ic++) {
        for (int i=0; i<n_samp; i++) {
            int32_t x = input[i];
            for (int ib=0; ib<n_filt; ib++) {
                const BFilterParams &_p = par[ib];
                auto &_w = w[ib][ic];
                int64_t c = (-_w[1] * _p.b0 + _w[0] * _p.b1) >> _p.b_bits;
                int64_t W = (x << _p.p_bits) - c;
                _w[0] = _w[1];
                _w[1] = W;
                x = (_w[0] + (_w[1]<<1) + W) >> _p.shift;
            }
            output[i] = x;
        }
    }
}

void BFilterBank::apply_to_float(float *input, float *output, float unit, int n_samp)
{
    int n_filt = w.size();
    assert(n_filt == par.size());
    int n_chan = w[0].size();

#pragma omp parallel for shared(input, output)
    for (int ic=0; ic<n_chan; ic++) {
        float *src = input  + ic*n_samp;
        float *dst = output + ic*n_samp;
        for (int i=0; i<n_samp; i++) {
            int32_t x = roundf(*(src++) / unit);
            for (int ib=0; ib<n_filt; ib++) {
                const BFilterParams &_p = par[ib];
                auto &_w = w[ib][ic];
                int64_t c = (-_w[1] * _p.b0 + _w[0] * _p.b1) >> _p.b_bits;
                int64_t W = (x << _p.p_bits) - c;
                _w[0] = _w[1];
                _w[1] = W;
                x = (_w[0] + (_w[1]<<1) + W) >> _p.shift;
            }
            *(dst++) = x * unit;
        }
    }
}


PYBINDINGS("so3g")
{
    bp::class_<BFilterParams>("BFilterParams",
                              bp::init<int32_t, int32_t, int, int, int>() );

    bp::class_<BFilterBank>("BFilterBank")
        .def("add", &BFilterBank::add,
             bp::return_internal_reference<>() )
        .def("init", &BFilterBank::init,
             bp::return_internal_reference<>() )
        .def("apply", &BFilterBank::apply_numpy);
}

