#include <assert.h>
#include <math.h>

#include <exceptions.h>
#include <Butterworth.h>

// Needed to work with numpy arrays in bindings
#include <nanobind/ndarray.h>

using namespace std;

namespace nb = nanobind;


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


void register_butterworth(nb::module_ & m) {
    nb::class_<BFilterParams>(m, "BFilterParams")
    .def(nb::init<int32_t, int32_t, int, int, int>())
    .def_rw("b0", &BFilterParams::b0)
    .def_rw("b1", &BFilterParams::b1)
    .def_rw("b_bits", &BFilterParams::b_bits)
    .def_rw("p_bits", &BFilterParams::p_bits)
    .def_rw("shift", &BFilterParams::shift);

    nb::class_<BFilterBank>(m, "BFilterBank")
    .def(nb::init<>())
    .def("add", &BFilterBank::add, nb::rv_policy::none)
    .def("init", &BFilterBank::init, nb::rv_policy::none)
    .def("apply", [](
        BFilterBank & slf,
        nb::ndarray<int32_t, nb::ndim<1>, nb::c_contig> input, 
        nb::ndarray<int32_t, nb::ndim<1>, nb::c_contig> output
    ) {
        auto v_in = input.view();
        auto v_out = output.view();
        auto n_samp = input.shape(0);
        if (v_out.shape(0) != n_samp) {
            throw agreement_exception("input", "output", "shape");
        }
        slf.apply(&v_in(0), &v_out(0), (int)n_samp);
    })
    .def("apply", [](
        BFilterBank & slf,
        nb::ndarray<float, nb::ndim<1>, nb::c_contig> input, 
        nb::ndarray<float, nb::ndim<1>, nb::c_contig> output
    ) {
        auto v_in = input.view();
        auto v_out = output.view();
        auto n_samp = input.shape(0);
        if (v_out.shape(0) != n_samp) {
            throw agreement_exception("input", "output", "shape");
        }
        slf.apply_to_float(&v_in(0), &v_out(0), 1.0, (int)n_samp);
    });

    return;
}

