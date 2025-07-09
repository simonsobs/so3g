#define NO_IMPORT_ARRAY

#include <complex>
#include <stdint.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
extern "C" {
    #include <cblas.h>
    // Additional prototypes for Fortran LAPACK routines.
    // dposv: solve Ax = b for A positive definite.
    void dposv_(const char* uplo, int* n, int* nrhs, double* a, int* lda,
                double* b, int* ldb, int* info );
}

#include <boost/python.hpp>
#ifdef _OPENMP
# include <omp.h>
#endif // ifdef _OPENMP

#include <gsl/gsl_spline.h>
#include <gsl/gsl_statistics.h>

#include <pybindings.h>
#include "so3g_numpy.h"
#include "numpy_assist.h"
#include "Ranges.h"

// TODO: Generalize to double precision too.
// This implements Jon's noise model for ACT. It takes in
// * ft[ndet,nfreq]        the fourier transform of the time-ordered data
// * bins[nbin,{from,to}]  the start and end of each bin
// * iD[nbin,ndet]         the inverse uncorrelated variance for each detector per bin
// * iV[nbin,ndet,nvec]    matrix representing the scaled eivenvectors per bin
void nmat_detvecs_apply(const bp::object & ft, const bp::object & bins, const bp::object & iD, const bp::object & iV, float s, float norm) {
    // Should pass in this too
    BufferWrapper<float>               ft_buf  ("ft",   ft,   false, std::vector<int>{-1,-1});
    BufferWrapper<int32_t>             bins_buf("bins", bins, false, std::vector<int>{-1, 2});
    int ndet = ft_buf->shape[0], nmode = ft_buf->shape[1], nbin = bins_buf->shape[0];
    BufferWrapper<float>               iD_buf  ("iD",   iD,   false, std::vector<int>{nbin,ndet});
    BufferWrapper<float>               iV_buf  ("iV",   iV,   false, std::vector<int>{nbin,ndet,-1});
    int nvec = iV_buf->shape[2];
    if (ft_buf->strides[1] != ft_buf->itemsize || ft_buf->strides[0] != ft_buf->itemsize*nmode)
        throw buffer_exception("ft must be C-contiguous along last axis");
    if (bins_buf->strides[1] != bins_buf->itemsize || bins_buf->strides[0] != bins_buf->itemsize*2)
        throw buffer_exception("bins must be C-contiguous along last axis");
    if (iD_buf->strides[1] != iD_buf->itemsize || iD_buf->strides[0] != iD_buf->itemsize*ndet)
        throw buffer_exception("iD must be C-contiguous along last axis");
    if (iV_buf->strides[2] != iV_buf->itemsize || iV_buf->strides[1] != iV_buf->itemsize*nvec || iV_buf->strides[0] != iV_buf->itemsize*nvec*ndet)
        throw buffer_exception("iV must be C-contiguous along last axis");
    // Internally we work with a real view of ft, with twice as many elements to compensate
    //int nmode = 2*nfreq;
    float   * ft_   = (float*)   ft_buf->buf;
    int32_t * bins_ = (int32_t*) bins_buf->buf;
    float   * iD_   = (float*)   iD_buf->buf;
    float   * iV_   = (float*)   iV_buf->buf;

    // Ok, actually do the work
    for(int bi = 0; bi < nbin; bi++) {
        int b1 = min(2*bins_[2*bi+0],nmode-1);
        int b2 = min(2*bins_[2*bi+1],nmode);
        int nm = b2-b1;
        float * biD = iD_ + bi*ndet;
        float * biV = iV_ + bi*ndet*nvec;

        // what I want to do
        // ft    = ftod[:,b[0]:b[1]]
        // iD    = self.iD[bi]/norm
        // iV    = self.iV[bi]/norm**0.5
        // ft[:] = iD[:,None]*ft + self.s*iV.dot(iV.T.dot(ft))
        // So first do iV.T [nvec,ndet] dot ft [ndet,nm] -> Q [nvec,nm]
        float * Q = new float[nvec*nm];
        cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, nvec, nm, ndet, 1.0f, biV, nvec, ft_+b1, nmode, 0.0f, Q, nm);
        // Handle the uncorrelated part
        //#pragma omp parallel for
        for(int di = 0; di < ndet; di++)
            for(int i = b1; i < b2; i++)
                ft_[di*nmode+i] *= biD[di]/norm;
        // Do ft += s*iV[ndet,nvec] dot Q [nvec,nm]
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, ndet, nm, nvec, s/norm, biV, nvec, Q, nm, 1.0f, ft_+b1, nmode);
        delete [] Q;
    }
}

// Support of maximum-liklihood sample cut handling. This got a bit long, so it should
// probably be moved into its own file.

// Forward declarations of helper functions
int get_dtype(const bp::object &);
int pcut_full_measure_helper(const vector<RangesInt32> &);
template <typename T> void pcut_full_tod2vals_helper(const vector<RangesInt32> &, T *, int, T *);
template <typename T> void pcut_full_vals2tod_helper(const vector<RangesInt32> &, T *, int, T *);
template <typename T> void pcut_full_translate_helper(const vector<RangesInt32> &, const vector<RangesInt32> &, int, T *, int, T*);
int pcut_poly_measure_helper(const vector<RangesInt32> &, int, int nmax);
template <typename T> void pcut_poly_tod2vals_helper(const vector<RangesInt32> &, int, int, T *, int, T *);
template <typename T> void pcut_poly_vals2tod_helper(const vector<RangesInt32> &, int, int, T *, int, T *);
template <typename T> void pcut_clear_helper(const vector<RangesInt32> &, T *, int);
template <typename T> void pcut_poly_translate_helper(const vector<RangesInt32> &, const vector<RangesInt32> &, int, int, int, T *, int, T *);

// The main cuts processing function. In Maximum-likelihood map-making cuts are handled
// as special degrees of freedom, and are part of the pointing matrix. This function
// provides the operations corresponding to the cuts part of the pointing matrix P and
// P'. In particular, operation == "insert" corresponds to the model-to-data matrix P,
// while operation == "extract" corresponds to the data-to-model matrix P'. Additionally,
// we provide operation == "measure", which doesn't touch any of the arrays, and just returns
// how long the vals argument should be.
//
// Note that "extract" does not in general estimate the coefficients of the model. It just does
// P'd, while the coefficients would (in the absence of any sky signal etc) be (P'P)'P'd.
//
// Arguments:
//  range_matrix: the .ranges member of a python RangesMatrix
//  operation: "measure", "insert" or "extract"
//  model: The type of model used
//   * "full": One degree of freedom per sample
//   * "poly": Model is a legendre polynomial in each range.
//        Order of poly determined by params["resolution"] (samples per order) and params["nmax"] (max order)
//  tod:  numpy array with shape [ndet,nsamp] of dtype float32 or float64
//  vals: numpy array with shape [:] of same dtype as tod. Holds the model degrees of freedom.
//
// TODO: To be able to process cuts in parallel, we need a lookup table for where in vals each
// cut range starts. This will be fast enough to build on the fly. Would pass this as an extra
// argument to the helper functions.

int process_cuts(const bp::object & range_matrix, const std::string & operation, const std::string & model, const bp::dict & params, const bp::object & tod, const bp::object & vals) {
    auto ranges = extract_ranges<int32_t>(range_matrix);
    // Decoding these up here lets us avoid some duplication later
    int resolution, nmax;
    if     (model == "full") {}
    else if(model == "poly") {
        resolution = bp::extract<int>(params.get("resolution"));
        nmax       = bp::extract<int>(params.get("nmax"));
    } else throw ValueError_exception("process_cuts model can only be 'full' or 'poly'");

    if(operation == "measure") {
        if     (model == "full") return pcut_full_measure_helper(ranges);
        else if(model == "poly") return pcut_poly_measure_helper(ranges, resolution, nmax);
    } else {
        int dtype = get_dtype(tod);
        if(dtype == NPY_FLOAT) {
            BufferWrapper<float> tod_buf  ("tod",  tod,  false, std::vector<int>{-1,-1});
            BufferWrapper<float> vals_buf ("vals", vals, false, std::vector<int>{-1});
            int ndet = tod_buf->shape[0], nsamp = tod_buf->shape[1];
            if(operation == "insert") {
                if     (model == "full")
                    pcut_full_vals2tod_helper(ranges, (float*)tod_buf->buf, nsamp, (float*) vals_buf->buf);
                else if(model == "poly")
                    pcut_poly_vals2tod_helper(ranges, resolution, nmax, (float*)tod_buf->buf, nsamp, (float*) vals_buf->buf);
            } else if(operation == "extract") {
                if     (model == "full")
                    pcut_full_tod2vals_helper(ranges, (float*)tod_buf->buf, nsamp, (float*) vals_buf->buf);
                else if(model == "poly")
                    pcut_poly_tod2vals_helper(ranges, resolution, nmax, (float*)tod_buf->buf, nsamp, (float*) vals_buf->buf);
            } else if(operation == "clear") {
                pcut_clear_helper(ranges, (float*)tod_buf->buf, nsamp);
            } else throw ValueError_exception("process_cuts operation can only be 'measure', 'insert' or 'extract'");
        } else if(dtype == NPY_DOUBLE) {
            BufferWrapper<double> tod_buf  ("tod",  tod,  false, std::vector<int>{-1,-1});
            BufferWrapper<double> vals_buf ("vals", vals, false, std::vector<int>{-1});
            int ndet = tod_buf->shape[0], nsamp = tod_buf->shape[1];
            if(operation == "insert") {
                if     (model == "full")
                    pcut_full_vals2tod_helper(ranges, (double*)tod_buf->buf, nsamp, (double*) vals_buf->buf);
                else if(model == "poly")
                    pcut_poly_vals2tod_helper(ranges, resolution, nmax, (double*)tod_buf->buf, nsamp, (double*) vals_buf->buf);
            } else if(operation == "extract") {
                if     (model == "full")
                    pcut_full_tod2vals_helper(ranges, (double*)tod_buf->buf, nsamp, (double*) vals_buf->buf);
                else if(model == "poly")
                    pcut_poly_tod2vals_helper(ranges, resolution, nmax, (double*)tod_buf->buf, nsamp, (double*) vals_buf->buf);
            } else if(operation == "clear") {
                pcut_clear_helper(ranges, (float*)tod_buf->buf, nsamp);
            } else throw ValueError_exception("process_cuts operation can only be 'measure', 'insert' or 'extract'");
        } else throw TypeError_exception("process_cuts only supports float32 and float64");
    }
    return 0;
}

void translate_cuts(const bp::object & irange_matrix, const bp::object & orange_matrix, const std::string & model, const bp::dict & params, const bp::object & ivals, bp::object & ovals) {
    // Decoding these up here lets us avoid some duplication later
    int resolution, nmax;
    if       (model == "full") {
        // nothing to do - res and nmax not used here
    } else if(model == "poly") {
        resolution = bp::extract<int>(params.get("resolution"));
        nmax       = bp::extract<int>(params.get("nmax"));
    } else {
        throw ValueError_exception("process_cuts model can only be 'full' or 'poly'");
    }
    auto iranges = extract_ranges<int32_t>(irange_matrix);
    auto oranges = extract_ranges<int32_t>(orange_matrix);
    int  insamp  = iranges[0].count;
    int  onsamp  = oranges[0].count;
    int  dtype   = get_dtype(ivals);
    if(dtype == NPY_FLOAT) {
        BufferWrapper<float> ivals_buf ("ivals", ivals, false, std::vector<int>{-1});
        BufferWrapper<float> ovals_buf ("ovals", ovals, false, std::vector<int>{-1});
        if(model == "full")
            pcut_full_translate_helper(iranges, oranges, insamp, (float*)ivals_buf->buf, onsamp, (float*)ovals_buf->buf);
        else if(model == "poly")
            pcut_poly_translate_helper(iranges, oranges, resolution, nmax, insamp, (float*)ivals_buf->buf, onsamp, (float*)ovals_buf->buf);
    } else if(dtype == NPY_DOUBLE) {
        BufferWrapper<double> ivals_buf ("ivals", ivals, false, std::vector<int>{-1});
        BufferWrapper<double> ovals_buf ("ovals", ovals, false, std::vector<int>{-1});
        if(model == "full")
            pcut_full_translate_helper(iranges, oranges, insamp, (double*)ivals_buf->buf, onsamp, (double*)ovals_buf->buf);
        else if(model == "poly")
            pcut_poly_translate_helper(iranges, oranges, resolution, nmax, insamp, (double*)ivals_buf->buf, onsamp, (double*)ovals_buf->buf);
    }
}

// Helpers for the cuts

int get_dtype(const bp::object & arr) {
    PyObject *ob = PyArray_FromAny(arr.ptr(), NULL, 0, 0, 0, NULL);
    if (ob == NULL) throw exception();
    PyArrayObject * a = reinterpret_cast<PyArrayObject*>(ob);
    int res = PyArray_TYPE(a);
    Py_DECREF(ob);
    return res;
}

// This is all from Jon's work for ACT.
int get_npoly(int gapsize, int resolution, int nmax) {
    if(nmax < 1) nmax = 1;
    if     (gapsize <=  1) return min<int>(1, nmax);
    else if(gapsize <=  3) return min<int>(2, nmax);
    else if(gapsize <=  6) return min<int>(3, nmax);
    else if(gapsize <= 20) return min<int>(4, nmax);
    else return min<int>(5 + gapsize/resolution, nmax);
}

// Full cut treatment, with one degree of freedom for each sample in the cut
int pcut_full_measure_helper(const vector<RangesInt32> & rangemat) {
    int n = 0;
    for(int di = 0; di < rangemat.size(); di++)
        for (auto const &r: rangemat[di].segments)
            n += r.second-r.first;
    return n;
}
template <typename T>
void pcut_full_tod2vals_helper(const vector<RangesInt32> & rangemat, T * tod, int nsamp, T * vals) {
    int i = 0;
    for(int di = 0; di < rangemat.size(); di++)
        for (auto const &r: rangemat[di].segments)
            for(int j = r.first; j < r.second; j++, i++)
                vals[i] = tod[di*nsamp+j];
}
template <typename T>
void pcut_full_vals2tod_helper(const vector<RangesInt32> & rangemat, T * tod, int nsamp, T * vals) {
    int i = 0;
    for(int di = 0; di < rangemat.size(); di++)
        for (auto const &r: rangemat[di].segments)
            for(int j = r.first; j < r.second; j++, i++)
                tod[di*nsamp+j] = vals[i];
}

template <typename T>
void pcut_full_translate_helper(const vector<RangesInt32> & iranges, const vector<RangesInt32> & oranges, int insamp, T * ivals, int onsamp, T * ovals) {
    // Translate cut degrees of freedom from one sample rate to another.
    // Simple interpolation for upscaling, averaging when downscaling
    int ioff = 0, ooff = 0;
    for(int di = 0; di < iranges.size(); di++) {
        for(int ri = 0; ri < iranges[di].segments.size(); ri++) {
            auto const & irange = iranges[di].segments[ri];
            auto const & orange = oranges[di].segments[ri];
            if(onsamp >= insamp) {
                for(int64_t osamp = orange.first, oind = ooff; osamp < orange.second; osamp++, oind++) {
                    int64_t isamp = osamp * insamp / onsamp;
                    int64_t iind  = ioff + isamp - irange.first;
                    // Nearest neighbor should be good enough
                    ovals[oind] = ivals[iind];
                }
            } else {
                vector<int> ohits(orange.second-orange.first);
                // Accumulate
                for(int64_t isamp = irange.first, iind = ioff; isamp < irange.second; isamp++, iind++) {
                    int64_t osamp = isamp * onsamp / insamp;
                    int64_t oind  = ooff + osamp - orange.first;
                    ovals[oind] += ivals[iind];
                    ohits[oind-ooff]++;
                }
                // average
                for(int osamp = orange.first, oind = ooff; osamp < orange.second; osamp++, oind++)
                    ovals[oind] /= ohits[oind-ooff];
            }
            ioff += irange.second - irange.first;
            ooff += orange.second - orange.first;
        }
    }
}

// Polynomial cut treatment
int pcut_poly_measure_helper(const vector<RangesInt32> & rangemat, int resolution, int nmax) {
    int n = 0;
    for(int di = 0; di < rangemat.size(); di++)
        for (auto const &r: rangemat[di].segments)
            n += get_npoly(r.second-r.first, resolution, nmax);
    return n;
}
template <typename T>
void pcut_poly_tod2vals_helper(const vector<RangesInt32> & rangemat, int resolution, int nmax, T * tod, int nsamp, T * vals) {
    int i = 0;
    for(int di = 0; di < rangemat.size(); di++) {
        for (auto const &r: rangemat[di].segments) {
            int np = get_npoly(r.second-r.first, resolution, nmax);
            if(np <= 1) {
                for(int s = r.first; s < r.second; s++)
                    vals[i] += tod[di*nsamp+s];
                i++;
            } else {
                for(int p = 0; p < np; p++) vals[i+p] = 0;
                for(int s = r.first; s < r.second; s++) {
                    T x = -1 + 2*(s-r.first)/T(r.second-r.first-1);
                    T t = tod[di*nsamp+s];
                    vals[i] += t;
                    if(np > 1) vals[i+1] += t*x;
                    if(np > 2) {
                        T Pa = x, Pb = 1, Pc = 0;
                        for(int p = 2; p < np; p++) {
                            Pc = Pb; Pb = Pa; Pa = ((2*p-1)*x*Pb-(p-1)*Pc)/p;
                            vals[i+p] += t*Pa;
                        }
                    }
                    i += np;
                }
            }
        }
    }
}
template <typename T>
void pcut_poly_vals2tod_helper(const vector<RangesInt32> & rangemat, int resolution, int nmax, T * tod, int nsamp, T * vals) {
    int i = 0;
    for(int di = 0; di < rangemat.size(); di++) {
        for (auto const &r: rangemat[di].segments) {
            int np = get_npoly(r.second-r.first, resolution, nmax);
            if(np <= 1) {
                for(int s = r.first; s < r.second; s++)
                    tod[di*nsamp+s] = vals[i];
                i++;
            } else {
                for(int s = r.first; s < r.second; s++) {
                    T x = -1 + 2*(s-r.first)/T(r.second-r.first-1);
                    T t = vals[i];
                    if(np > 1) t += x*vals[i+1];
                    if(np > 2) {
                        T Pa = x, Pb = 1, Pc = 0;
                        for(int p = 2; p < np; p++) {
                            Pc = Pb; Pb = Pa; Pa = ((2*p-1)*x*Pb-(p-1)*Pc)/p;
                            t += Pa*vals[i+p];
                        }
                    }
                    tod[di*nsamp+s] = t;
                    i += np;
                }
            }
        }
    }
}
template <typename T>
void pcut_clear_helper(const vector<RangesInt32> & rangemat, T * tod, int nsamp) {
    #pragma omp parallel for
    for(int di = 0; di < rangemat.size(); di++)
        for (auto const &r: rangemat[di].segments)
            for(int s = r.first; s < r.second; s++)
                tod[di*nsamp+s] = 0;
}

template <typename T>
void pcut_poly_translate_helper(const vector<RangesInt32> & iranges, const vector<RangesInt32> & oranges, int resolution, int nmax, int insamp, T * ivals, int onsamp, T * ovals) {
    // Translate cut degrees of freedom from one sample rate to another.
    // Zero-pad poly coeffs if upsampling, truncate if downsampling
    int ioff = 0, ooff = 0;
    for(int di = 0; di < iranges.size(); di++) {
        for(int ri = 0; ri < iranges[di].segments.size(); ri++) {
            auto const & irange = iranges[di].segments[ri];
            auto const & orange = oranges[di].segments[ri];
            int inp = get_npoly(irange.second-irange.first, resolution, nmax);
            int onp = get_npoly(orange.second-orange.first, resolution, nmax);
            int i;
            for(i = 0; i < min(inp,onp); i++)
                ovals[ooff+i] = ivals[ioff+i];
            for(; i < onp; i++)
                ovals[ooff+i] = 0;
            ioff += inp;
            ooff += onp;
        }
    }
}

// get_gap_fill_poly_single
//
// Single detector processor for get_gap_fill_poly.  Fit polynomials
// to each stretch of data in "gaps".  Requires square matrix a and
// vector b to be preallocated for ncoeff.
//
// The arguments (inplace, extract) behave like this:
//
//   (true,  valid):   The data are replaced with the model, and the
//                     replaced samples are copied into *extract.
//   (false, valid):   The data are left unchanged and the model
//                     values are placed in *extract.
//   (true,  nullptr): The data are replaced with the model, and
//                     the original data samples are discarded.
//   (false, nullptr): The data is left unchanged and the model
//                     is discarded (probably not what you want...).
//
// When extract is non-null, it needs to be the right size... this
// should be checked by get_gap_fill_poly.

template <typename T>
void get_gap_fill_poly_single(const RangesInt32 &gaps, T *data,
                              double *a, double *b,
                              int buffer, int ncoeff,
                              bool inplace, T *extract)
{
    // Generate Ranges corresponding to samples near the edges of
    // intervals we want to fill.
    RangesInt32 rsegs = gaps.buffered((int32_t)buffer);
    rsegs.intersect(gaps.complement());

    // We are guaranteed that there is one or two rseg intervals
    // between each of our gaps, and zero or one rseg intervals
    // before the first gap and after the second gap.

    // Index of the rseg we're working with.
    int model_i = 0;

    for (auto const &gap: gaps.segments) {
        // std::cout << "GAP\n";
        int contrib_samps = 0;
        int contrib_segs = 0;
        memset(a, 0, ncoeff*ncoeff*sizeof(*a));
        memset(b, 0, ncoeff*sizeof(*b));
        double x0 = gap.first;

        for (; model_i < rsegs.segments.size(); model_i++) {
            // Advance until just before this gap.
            auto const &ival = rsegs.segments[model_i];
            if (ival.second + 1 < gap.first)
                continue;
            // Include this interval in the fit.
            for (int i=ival.first; i<ival.second; i++) {
                double xx = 1.;
                double yxx = data[i];
                double x = i - x0;
                for (int j=0; j<ncoeff; j++) {
                    b[j] += yxx;
                    yxx *= x;
                    a[j] += xx;
                    xx *= x;
                }
                // Now down...
                for (int k=2; k<ncoeff+1; k++) {
                    a[k*ncoeff - 1] += xx;
                    xx *= x;
                }
            }
            contrib_samps += ival.second - ival.first;
            contrib_segs += 1;
            // If this was the right side interval, bail
            // now (so we can possibly re-use this
            // interval for next gap).
            if (ival.first > gap.first)
                break;
        }

        // Restrict order based on number of contributing samples.
        int n_keep = std::min(contrib_samps / 10 + 1, ncoeff);
        if (contrib_samps > 0) {
            // Fill in a's interior.
            for (int r=1; r<ncoeff; r++)
                for (int c=0; c<ncoeff-1; c++)
                    a[r*ncoeff + c] = a[r*ncoeff + c - ncoeff + 1];

            // Re-organize a if the order has changed...
            if (n_keep < ncoeff) {
                for (int r=1; r<n_keep; r++)
                    for (int c=0; c<n_keep; c++)
                        a[r*n_keep+c] = a[r*ncoeff+c];
            }

            // Solve the system...
            int one = 1;
            int err = 0;
            dposv_("Upper", &n_keep, &one, a, &n_keep, b, &n_keep, &err);
        }
        T *write_to = nullptr;
        T *save_data = nullptr;
        if (inplace) {
            // Copy original data to extract, write model to data.
            save_data = extract;
            write_to = data + gap.first;
        } else {
            // Write results to extract, do not touch data.
            write_to = extract;
        }

        if (save_data != nullptr) {
            for (int i=gap.first; i<gap.second; i++, save_data++)
                *save_data = data[i];
        }
        if (write_to != nullptr) {
            for (int i=gap.first; i<gap.second; i++, write_to++) {
                double xx = 1.;
                *write_to = 0.;
                for (int j=0; j<n_keep; j++) {
                    *write_to += xx * b[j];
                    xx *= (i - gap.first);
                }
            }
        }
        if (extract != nullptr)
            extract += (gap.second - gap.first);
    }
}

template <typename T>
void get_gap_fill_poly(const bp::object ranges,
                       const bp::object tod,
                       int buffer,
                       int order,
                       bool inplace,
                       const bp::object ex)
{
    // As a test, copy data from rangemat into segment.
    auto rangemat = extract_ranges<int32_t>(ranges);
    int ndet = rangemat.size();

    BufferWrapper<T> tod_buf  ("tod",  tod,  false, std::vector<int>{ndet,-1});
    int nsamp = tod_buf->shape[1];

    int ncoeff = order + 1; // Let us not speak of order again.
    double *a = (double*)malloc(ncoeff*(ncoeff+1)*sizeof(*a));
    double *b = a + ncoeff*ncoeff;

    T *ex_data = nullptr;
    std::vector<int> ex_offsets;

    if (ex.ptr() != Py_None) {
        // Compute offsets of each detector into ex.
        int n = 0;
        for (auto const &r: rangemat) {
            ex_offsets.push_back(n);
            for (auto const &r: r.segments)
                n += (r.second - r.first);
        }
        BufferWrapper<T> ex_buf("ex", ex, false, std::vector<int>{n});
        ex_data = (T*)ex_buf->buf;
    }

    for (int di=0; di < rangemat.size(); di++) {
        T* data = (T*)((char*)tod_buf->buf + di*tod_buf->strides[0]);
        T* _ex = ex_data;
        if (_ex != nullptr)
            _ex += ex_offsets[di];
        get_gap_fill_poly_single(rangemat[di], data, a, b, buffer, ncoeff,
                                 inplace, _ex);
    }

    free(a);
}


void test_buffer_wrapper(const bp::object array,
                         const bp::object dims)
{
    std::vector<int> _dims(bp::len(dims));
    for (int i=0; i<bp::len(dims); i++)
        _dims[i] = bp::extract<double>(dims[i]);
    BufferWrapper<double> array_buf  ("array",  array,  false, _dims);

}


template <typename T>
void _moment(T* data, T* output, int moment, bool central, int start, int stop)
{
    int bsize = stop - start;
    // Could replace the loops with boost accumulators?
    T center = 0.0;
    if(central || moment == 1) {
        for(int si = start; si < stop; si++) {
            center = center + data[si];
        }
        center = center / bsize;
    }
    T val = 0;
    if(moment == 1) {
        val = center;
    }
    else {
        for(int si = start; si < stop; si++) {
            val = val + pow(data[si] - center, moment);
        }
        val = val / bsize;
    }
    for(int si = start; si < stop; si++) {
        output[si] = val;
    }

}

template <typename T>
void _block_moment(T* tod_data, T* output, int bsize, int moment, bool central, int ndet, int nsamp, int shift)
{
    int nblock = ((nsamp - shift)/bsize) + 1;
    #pragma omp parallel for
    for(int di = 0; di < ndet; di++)
    {
        int ioff = di*nsamp;
        // do the the pre-shift portion
        if(shift > 0){
            _moment(tod_data, output, moment, central, ioff, ioff+shift);
        }

        for(int bi = 0; bi < nblock; bi++) {
            int start =  (bi * bsize) + shift;
            int stop = min(start + bsize, nsamp);
            _moment(tod_data, output, moment, central, ioff+start, ioff+stop);
        }
    }
}

template <typename T>
void block_moment(const bp::object & tod, const bp::object & out, int bsize, int moment, bool central, int shift)
{
    BufferWrapper<T> tod_buf  ("tod",  tod,  false, std::vector<int>{-1, -1});
    int ndet = tod_buf->shape[0];
    int nsamp = tod_buf->shape[1];
    T* tod_data = (T*)tod_buf->buf;
    if (tod_buf->strides[1] != tod_buf->itemsize || tod_buf->strides[0] != tod_buf->itemsize*nsamp)
        throw buffer_exception("tod must be C-contiguous along last axis");
    BufferWrapper<T> out_buf  ("out",  out,  false, std::vector<int>{ndet, nsamp});
    if (out_buf->strides[1] != out_buf->itemsize || out_buf->strides[0] != out_buf->itemsize*nsamp)
        throw buffer_exception("out must be C-contiguous along last axis");
    T* output = (T*)out_buf->buf;
    _block_moment(tod_data, output, bsize, moment, central, ndet, nsamp, shift);
}

template <typename T>
void _minmax(T* data, T* output, int mode, int start, int stop)
{
    T val;
    if(mode == 0){ // get the min
        val = *(std::min_element(data+start, data+stop));
    }
    else if(mode == 1){ // get the max
        val = *(std::max_element(data+start, data+stop));
    }
    else{ // get the peak to peak
        auto min = std::min_element(data+start, data+stop);
        auto max = std::max_element(data+start, data+stop);
        val = *max - *min;
    }
    for(int si = start; si < stop; si++) {
        output[si] = val;
    }

}

template <typename T>
void _block_minmax(T* tod_data, T* output, int bsize, int mode, int ndet, int nsamp, int shift)
{
    int nblock = ((nsamp - shift)/bsize) + 1;
    #pragma omp parallel for
    for(int di = 0; di < ndet; di++)
    {
        int ioff = di*nsamp;
        // do the the pre-shift portion
        if(shift > 0){
            _minmax(tod_data, output, mode, ioff, ioff+shift);
        }

        for(int bi = 0; bi < nblock; bi++) {
            int start =  (bi * bsize) + shift;
            int stop = min(start + bsize, nsamp);
            _minmax(tod_data, output, mode, ioff+start, ioff+stop);
        }
    }
}

template <typename T>
void block_minmax(const bp::object & tod, const bp::object & out, int bsize, int mode, int shift)
{
    BufferWrapper<T> tod_buf  ("tod",  tod,  false, std::vector<int>{-1, -1});
    int ndet = tod_buf->shape[0];
    int nsamp = tod_buf->shape[1];
    T* tod_data = (T*)tod_buf->buf;
    if (tod_buf->strides[1] != tod_buf->itemsize || tod_buf->strides[0] != tod_buf->itemsize*nsamp)
        throw buffer_exception("tod must be C-contiguous along last axis");
    BufferWrapper<T> out_buf  ("out",  out,  false, std::vector<int>{ndet, nsamp});
    if (out_buf->strides[1] != out_buf->itemsize || out_buf->strides[0] != out_buf->itemsize*nsamp)
        throw buffer_exception("out must be C-contiguous along last axis");
    T* output = (T*)out_buf->buf;
    _block_minmax(tod_data, output, bsize, mode, ndet, nsamp, shift);
}


void _clean_flag(int* flag_data, int width, int ndet, int nsamp)
{
    #pragma omp parallel for
    for(int di = 0; di < ndet; di++) {
        int ioff = di*nsamp;
        int* det_flag = flag_data + ioff;
        int count = 0;
        for(int si = 0; si < nsamp; si++) {
            // If we run into a 0
            if(det_flag[si]==0) {
                // If this block was too small
                if(count<width) {
                    for(int i = si - count; i < si; i++){
                        det_flag[i] = 0;
                    }
                }
                // Reset count
                count = 0;
            }
            else {
                count = count + 1;
            }
        }
    }
}

void clean_flag(const bp::object & flag, int width)
{
    BufferWrapper<int> flag_buf  ("flag", flag, false, std::vector<int>{-1, -1});
    int ndet = flag_buf->shape[0];
    int nsamp = flag_buf->shape[1];
    if (flag_buf->strides[1] != flag_buf->itemsize || flag_buf->strides[0] != flag_buf->itemsize*nsamp)
        throw buffer_exception("flag must be C-contiguous along last axis");
    int* flag_data = (int*)flag_buf->buf;
    _clean_flag(flag_data, width, ndet, nsamp);

}

template <typename T>
void _jumps_thresh_on_mfilt(T* mfilt, int* flag, T* size, int bsize, int shift, T prefac, bool samp_uncert, bool check_flag, int ndet, int nsamp)
{
    // The sensivity of the filter: s = (n_left * n_right)/(n_left + n_right)
    // For a jump of height j, the filter responce will be s*j
    // We use this to kill features smaller than we expect from the min jump size
    // Because s = 0 at the window edges we skip those indices
    # pragma omp parallel for
    for(int di = 0; di < ndet; di++){
        int ioff = di*nsamp;
        for(int si = 0; si < nsamp; si++){
            if(si < shift){
                flag[ioff+si] = 0;
                continue;
            }
            int _si = si - shift;
            T n_left = (_si % bsize) + 1;
            T n_right = (min(nsamp - shift, bsize*(1 + (_si/bsize))) - _si) - 1;
            T thresh = prefac * size[di] * (n_left * n_right) / (n_left + n_right);
            if((n_left+.5<(float)bsize / 4.)||(n_right+.5<(float) bsize / 4.)){
                flag[ioff+si] = 0;
                continue;
            }
            if(samp_uncert){
                thresh = thresh * (1 - ((abs(n_left - n_right) - .5)/(2*n_left*n_right)));
            }
            int flagged = (abs(mfilt[ioff+si]) > thresh);
            if(check_flag){
                flagged = flagged && (flag[ioff+si] == 1);
            }
            flag[ioff+si] = flagged;
        }
    }
}

template <typename T>
void _jumps_matched_filter(T* tod_data, T* output, int bsize, int shift, int ndet, int nsamp)
{
    // Get the matched filter, this is basically convolving with a step
    _block_moment(tod_data, output, bsize, 1, 0, ndet, nsamp, shift);
    #pragma omp parallel for
    for(int di = 0; di < ndet; di++) {
        int ioff = di*nsamp;
        T val = 0;
        for(int si = 0; si < nsamp; si++) {
            int i = ioff + si;
            val = val + tod_data[i] - output[i];
            output[i] = val;
        }
    }
}

template <typename T>
void matched_jumps(const bp::object & tod, const bp::object & out, const bp::object & min_size, int bsize)
{
    BufferWrapper<T> tod_buf  ("tod",  tod,  false, std::vector<int>{-1, -1});
    int ndet = tod_buf->shape[0];
    int nsamp = tod_buf->shape[1];
    if (tod_buf->strides[1] != tod_buf->itemsize || tod_buf->strides[0] != tod_buf->itemsize*nsamp)
        throw buffer_exception("tod must be C-contiguous along last axis");
    T* tod_data = (T*)tod_buf->buf;
    BufferWrapper<int> out_buf  ("out",  out,  false, std::vector<int>{ndet, nsamp});
    if (out_buf->strides[1] != out_buf->itemsize || out_buf->strides[0] != out_buf->itemsize*nsamp)
        throw buffer_exception("out must be C-contiguous along last axis");
    int* output = (int*)out_buf->buf;
    BufferWrapper<T> size_buf ("min_size",  min_size,  false, std::vector<int>{ndet});
    if (size_buf->strides[0] != size_buf->itemsize)
        throw buffer_exception("min_size must be C-contiguous along last axis");
    T* size = (T*)size_buf->buf;
    T* buffer = new T[ndet * nsamp];
    
    int half_win = bsize / 2;
    int quarter_win = bsize / 4;

    // Get the first matched filter, this is basically convolving with a step
    _jumps_matched_filter(tod_data, buffer, bsize, 0, ndet, nsamp);
    // For this first round of cleaning we use min_size/2
    // Note that after this filtering we are left with at least win_size/4 width
    _jumps_thresh_on_mfilt(buffer, output, size, bsize, 0, (T).5, false, false, ndet, nsamp);
    // Clean spurs 
    _clean_flag(output, quarter_win, ndet, nsamp);
    // Recall that we set _min_size to be half the actual peak min above
    // We allow for .5 samples worth of uncertainty here
    _jumps_thresh_on_mfilt(buffer, output, size, bsize, 0, (T)1., true, true, ndet, nsamp);
    
    // Now do the shifted filter
    _jumps_matched_filter(tod_data, buffer, bsize, half_win, ndet, nsamp);
    int* shift_flag = new int[ndet * nsamp];
    _jumps_thresh_on_mfilt(buffer, shift_flag, size, bsize, half_win, (T).5, false, false, ndet, nsamp);
    _clean_flag(shift_flag, quarter_win, ndet, nsamp);
    _jumps_thresh_on_mfilt(buffer, shift_flag, size, bsize, half_win, (T)1., true, true, ndet, nsamp);
    delete buffer;

    // Now we combine
    #pragma omp parallel for
    for(int di = 0; di < ndet; di++) {
        int ioff = di*nsamp;
        for(int si = 0; si < nsamp; si++) {
            int i = ioff + si;
            output[i] = output[i] || shift_flag[i]; 
        }
    }
    delete shift_flag;
}

template <typename T>
void find_quantized_jumps(const bp::object & tod, const bp::object & out, const bp::object & atol, int win_size, T scale)
{
    BufferWrapper<T> tod_buf  ("tod",  tod,  false, std::vector<int>{-1, -1});
    int ndet = tod_buf->shape[0];
    int nsamp = tod_buf->shape[1];
    if (tod_buf->strides[1] != tod_buf->itemsize || tod_buf->strides[0] != tod_buf->itemsize*nsamp)
        throw buffer_exception("tod must be C-contiguous along last axis");
    T* tod_data = (T*)tod_buf->buf;
    BufferWrapper<T> out_buf  ("out",  out,  false, std::vector<int>{ndet, nsamp});
    if (out_buf->strides[1] != out_buf->itemsize || out_buf->strides[0] != out_buf->itemsize*nsamp)
        throw buffer_exception("out must be C-contiguous along last axis");
    T* output = (T*)out_buf->buf;
    BufferWrapper<T> tol_buf  ("atol",  atol,  false, std::vector<int>{ndet});
    if (tol_buf->strides[0] != tol_buf->itemsize)
        throw buffer_exception("atol must be C-contiguous along last axis");
    T* tol = (T*)tol_buf->buf;

    #pragma omp parallel for
    for(int di = 0; di < ndet; di++) {
        int ioff = di*nsamp;
        T* det_data = tod_data + ioff;
        T* det_out = output + ioff;
        for(int si = 0; si < nsamp; si++) {
            T val;
            int idx = 0;
            if(si > win_size) {
                idx = si - win_size;
            }
            T ratio = (det_data[si] - det_data[idx])/scale;
            if(abs(ratio) < .5) {
                det_out[si] = 0;
                continue;
            }
            T rounded = (T)round(ratio);
            if(abs(ratio - rounded) <= tol[di]) {
                det_out[si] = rounded * scale;
            }
            else {
                det_out[si] = 0;
            }
        }
    }
}

template <typename T>
void subtract_jump_heights(const bp::object & tod, const bp::object & out, const bp::object & heights, const bp::object & jumps) {
    BufferWrapper<T> tod_buf  ("tod",  tod,  false, std::vector<int>{-1, -1});
    int ndet = tod_buf->shape[0];
    int nsamp = tod_buf->shape[1];
    if (tod_buf->strides[1] != tod_buf->itemsize || tod_buf->strides[0] != tod_buf->itemsize*nsamp)
        throw buffer_exception("tod must be C-contiguous along last axis");
    T* tod_data = (T*)tod_buf->buf;
    BufferWrapper<T> out_buf  ("out",  out,  false, std::vector<int>{ndet, nsamp});
    if (out_buf->strides[1] != out_buf->itemsize || out_buf->strides[0] != out_buf->itemsize*nsamp)
        throw buffer_exception("out must be C-contiguous along last axis");
    T* output = (T*)out_buf->buf;
    BufferWrapper<T> h_buf  ("heights",  heights,  false, std::vector<int>{ndet, nsamp});
    if (h_buf->strides[1] != h_buf->itemsize || h_buf->strides[0] != h_buf->itemsize*nsamp)
        throw buffer_exception("heights must be C-contiguous along last axis");
    T* jump_heights = (T*)h_buf->buf;
    auto ranges = extract_ranges<int32_t>(jumps);

    #pragma omp parallel for
    for(int di = 0; di < ranges.size(); di++) {
        int start = 0;
        int stop = 0;
        T min_h;
        T max_h;
        T height;
        T to_sub = 0;
        for (auto const &r: ranges[di].segments) {
            start = di*nsamp + r.first;
            for(int j = stop; j < start && to_sub != 0; j++) {
                output[j] = tod_data[j] - to_sub;
            }
            stop = di*nsamp + r.second;
            min_h = *(std::min_element(jump_heights+start, jump_heights+stop));
            max_h = *(std::max_element(jump_heights+start, jump_heights+stop));
            // Decide whether this is a negative or positive jump.
            height = (abs(min_h) > abs(max_h)) ? min_h : max_h;
            to_sub = to_sub + height;
            for(int j = start; j < stop && to_sub != 0; j++) {
                output[j] = tod_data[j] - to_sub;
            }
        }
        for(int j = stop; j < di*nsamp + nsamp && to_sub != 0; j++) {
            output[j] = tod_data[j] - to_sub;
        }
    }
}

template <typename T>
using _interp_func_pointer = void (*)(const double* x, const double* y,
                                      const double* x_interp, T* y_interp,
                                      const int n_x, const int n_x_interp,
                                      gsl_spline* spline, gsl_interp_accel* acc);

template <typename T>
void _linear_interp(const double* x, const double* y, const double* x_interp,
                    T* y_interp, const int n_x, const int n_x_interp,
                    gsl_spline* spline, gsl_interp_accel* acc)
{
    // Re-initialize for each row
    gsl_spline_init(spline, x, y, n_x);

    double x_step_left = x[1] - x[0];
    double x_step_right = x[n_x - 1] - x[n_x - 2];

    double x_min = x[0];
    double x_max = x[n_x - 1];

    double slope_left = (y[1] - y[0]) / x_step_left;
    double slope_right = (y[n_x - 1] -  y[n_x - 2]) / x_step_right;

    for (int si = 0; si < n_x_interp; ++si) {
        // Points below minimum value
        if (x_interp[si] < x_min) {
            y_interp[si] = y[0] + slope_left * (x_interp[si] - x_min);
        }
        // Points above maximum value
        else if (x_interp[si] >= x_max) {
            y_interp[si] = y[n_x - 1] + slope_right * (x_interp[si] - x_max);
        } 
        else {
            y_interp[si] = gsl_spline_eval(spline, x_interp[si], acc);
        }
    }
}

template <typename T>
void _interp1d(const bp::object & x, const bp::object & y, const bp::object & x_interp,
               bp::object & y_interp, const gsl_interp_type* interp_type,
               _interp_func_pointer<T> interp_func)
{
    BufferWrapper<T> y_buf  ("y",  y,  false, std::vector<int>{-1, -1});
    if (y_buf->strides[1] != y_buf->itemsize)
        throw ValueError_exception("Argument 'y' must be contiguous in last axis.");
    T* y_data = (T*)y_buf->buf;
    const int n_rows = y_buf->shape[0];
    const int n_x = y_buf->shape[1];

    BufferWrapper<T> x_buf  ("x",  x,  false, std::vector<int>{n_x});
    if (x_buf->strides[0] != x_buf->itemsize)
        throw ValueError_exception("Argument 'x' must be a C-contiguous 1d array");
    T* x_data = (T*)x_buf->buf;

    BufferWrapper<T> y_interp_buf  ("y_interp",  y_interp,  false, std::vector<int>{n_rows, -1});
    if (y_interp_buf->strides[1] != y_interp_buf->itemsize)
        throw ValueError_exception("Argument 'y_interp' must be contiguous in last axis.");
    T* y_interp_data = (T*)y_interp_buf->buf;
    const int n_x_interp = y_interp_buf->shape[1];

    BufferWrapper<T> x_interp_buf  ("x_interp",  x_interp,  false, std::vector<int>{n_x_interp});
    if (x_interp_buf->strides[0] != x_interp_buf->itemsize)
        throw ValueError_exception("Argument 'x_interp' must be a C-contiguous 1d array");
    T* x_interp_data = (T*)x_interp_buf->buf;

    if constexpr (std::is_same<T, double>::value) {
        // Strides for non-contiguous rows
        int y_data_stride = y_buf->strides[0] / sizeof(double);
        int y_interp_data_stride = y_interp_buf->strides[0] / sizeof(double);

        #pragma omp parallel
        {
            // Create one accel and spline per thread
            gsl_interp_accel* acc = gsl_interp_accel_alloc();
            gsl_spline* spline = gsl_spline_alloc(interp_type, n_x);

            #pragma omp for
            for (int row = 0; row < n_rows; ++row) {

                int y_row_start = row * y_data_stride;
                int y_row_end = y_row_start + n_x;
                int y_interp_row_start = row * y_interp_data_stride;

                T* y_row = y_data + y_row_start;
                T* y_interp_row = y_interp_data + y_interp_row_start;

                interp_func(x_data, y_row, x_interp_data, y_interp_row,
                            n_x, n_x_interp, spline, acc);
            }

            // Free gsl objects
            gsl_spline_free(spline);
            gsl_interp_accel_free(acc);
        }
    }
    else if constexpr (std::is_same<T, float>::value) {
        // Strides for non-contiguous rows
        int y_data_stride = y_buf->strides[0] / sizeof(float);
        int y_interp_data_stride = y_interp_buf->strides[0] / sizeof(float);

        // Transform x and x_interp to double arrays for gsl
        double x_dbl[n_x], x_interp_dbl[n_x_interp];

        std::transform(x_data, x_data + n_x, x_dbl, 
                       [](float value) { return static_cast<double>(value); });

        std::transform(x_interp_data, x_interp_data + n_x_interp, x_interp_dbl,
                       [](float value) { return static_cast<double>(value); });

        #pragma omp parallel
        {
            // Create one accel and spline per thread
            gsl_interp_accel* acc = gsl_interp_accel_alloc();
            gsl_spline* spline = gsl_spline_alloc(interp_type, n_x);

            #pragma omp for
            for (int row = 0; row < n_rows; ++row) {

                int y_row_start = row * y_data_stride;
                int y_row_end = y_row_start + n_x;
                int y_interp_row_start = row * y_interp_data_stride;

                // Transform y row to double array for gsl
                double y_dbl[n_x];

                std::transform(y_data + y_row_start, y_data + y_row_end, y_dbl,
                           [](float value) { return static_cast<double>(value); });

                T* y_interp_row = y_interp_data + y_interp_row_start;

                // Don't copy y_interp to doubles as it is cast during assignment
                interp_func(x_dbl, y_dbl, x_interp_dbl, y_interp_row,
                            n_x, n_x_interp, spline, acc);
            }

            // Free gsl objects
            gsl_spline_free(spline);
            gsl_interp_accel_free(acc);
        }
    }
}

void interp1d_linear(const bp::object & x, const bp::object & y,
                     const bp::object & x_interp, bp::object & y_interp)
{
    // Get data type
    int dtype = get_dtype(y);

    if (dtype == NPY_FLOAT) {
        // GSL interpolation type
        const gsl_interp_type* interp_type = gsl_interp_linear;
        // Pointer to interpolation function
        _interp_func_pointer<float> interp_func = &_linear_interp<float>;

        _interp1d<float>(x, y, x_interp, y_interp, interp_type, interp_func);
    }
    else if (dtype == NPY_DOUBLE) {
        // GSL interpolation type
        const gsl_interp_type* interp_type = gsl_interp_linear;
        // Pointer to interpolation function
        _interp_func_pointer<double> interp_func = &_linear_interp<double>;

        _interp1d<double>(x, y, x_interp, y_interp, interp_type, interp_func);
    }
    else {
        throw TypeError_exception("Only float32 or float64 arrays are supported.");
    }
}

template <typename T>
T _calculate_median(const T* data, const int n)
{
    // Copy to prevent overwriting input with gsl median
    // Explicitly cast to double here due to gsl
    std::vector<double> data_copy(n);
    std::transform(data, data + n, data_copy.begin(), [](double val) {
        return static_cast<double>(val);
    });

    // GSL is much faster than a naive std::sort implementation
    return gsl_stats_median(data_copy.data(), 1, n);
}

template <typename T>
void _detrend(T* data, const int ndets, const int nsamps, const int row_stride,
              const std::string & method, const int linear_ncount,
              const int nthreads)
{
    if (method == "mean") {
        #pragma omp parallel for num_threads(nthreads)
        for (int i = 0; i < ndets; ++i) {
            int ioff = i * row_stride;

            T* data_row = data + ioff;

            // This is significantly faster than gsl_stats_mean
            T det_mean = 0.;
            for (int si = 0; si < nsamps; ++si) {
                det_mean += data_row[si];
            }

            det_mean /= nsamps;

            for (int si = 0; si < nsamps; ++si) {
                data_row[si] -= det_mean;
            }
        }
    }
    else if (method == "median") {
        #pragma omp parallel for num_threads(nthreads)
        for (int i = 0; i < ndets; ++i) {
            int ioff = i * row_stride;

            T* data_row = data + ioff;

            T det_median = _calculate_median(data_row, nsamps);

            for (int si = 0; si < nsamps; ++si) {
                data_row[si] -= det_median;
            }
        }
    }
    else if (method == "linear") {
        // Default ncount
        int ncount = linear_ncount;
        if (ncount == -1) {
            ncount = nsamps / 2;
        }

        T x[nsamps];
        T step = 1.0 / (nsamps - 1);

        // Equivalent to np.linspace(0.,1.,nsamp)
        for (int si = 0; si < nsamps; ++si) {
            x[si] = si * step;
        }

        ncount = std::max(1, std::min(ncount, nsamps / 2));

        int last_offset = nsamps - ncount;

        #pragma omp parallel for num_threads(nthreads)
        for (int i = 0; i < ndets; ++i) {
            int ioff = i * row_stride;

            T* data_row = data + ioff;

            // Mean of first and last ncount samples
            T det_mean_first = 0.;
            T det_mean_last = 0.;

            for (int si = 0; si < ncount; ++si) {
                det_mean_last += data_row[si + last_offset];
                det_mean_first += data_row[si];
            }

            T slope = (det_mean_last - det_mean_first) / ncount;

            T det_mean = 0.;
            for (int si = 0; si < nsamps; ++si) {
                data_row[si] -= slope * x[si];
                det_mean += data_row[si];
            }

            det_mean /= nsamps;
            for (int si = 0; si < nsamps; ++si) {
                data_row[si] -= det_mean;
            }
        }
    }
    else {
        throw ValueError_exception("Unupported detrend method. Supported methods "
                                   "are 'mean', 'median', and 'linear'");
    }
}

template <typename T>
void _detrend_buffer(bp::object & tod, const std::string & method,
                     const int linear_ncount)
{
    BufferWrapper<T> tod_buf  ("tod",  tod,  false, std::vector<int>{-1, -1});
    if (tod_buf->strides[1] != tod_buf->itemsize)
        throw ValueError_exception("Argument 'tod' must be contiguous in last axis.");
    T* tod_data = (T*)tod_buf->buf;
    const int ndets = tod_buf->shape[0];
    const int nsamps = tod_buf->shape[1];

    int row_stride = tod_buf->strides[0] / sizeof(T);

    // _detrend may be called from within a parallel loop internally, so manage
    // parallelization explicitly
    int nthreads = 1;
     #pragma omp parallel
    {
        #ifdef _OPENMP
        if (omp_get_thread_num() == 0)
            nthreads = omp_get_num_threads();
        #endif
    }

    // We want _detrend to accept C++ types so it can be used internally
    // for Welch psd calculations, hence the hierarchical function calls
    _detrend<T>(tod_data, ndets, nsamps, row_stride, method, linear_ncount, nthreads);
}

void detrend(bp::object & tod, const std::string & method, const int linear_ncount)
{
    // Get data type
    int dtype = get_dtype(tod);

    if (dtype == NPY_FLOAT) {
        _detrend_buffer<float>(tod, method, linear_ncount);
    }
    else if (dtype == NPY_DOUBLE) {
        _detrend_buffer<double>(tod, method, linear_ncount);
    }
    else {
        throw TypeError_exception("Only float32 or float64 arrays are supported.");
    }
}

template <typename T>
int _find_bin_index(const T* bin_edges, T value, int nbins) {
    int left = 0;
    int right = nbins;

    // Mimic np.clip and assign out-of-bounds points to
    // first and last bin
    if (value < bin_edges[left]) {
        return 0;
    }
    else if (value >= bin_edges[right]) {
        return right - 1;
    }

    while (left < right) {
        int mid = left + (right - left) / 2;

        if (value >= bin_edges[mid] && value < bin_edges[mid + 1]) {
            return mid;
        }
        else if (value >= bin_edges[mid + 1]) {
            left = mid + 1;
        }
        else {
            right = mid;
        }
    }

    return -1;
}

template <typename T>
void _bin_signal(const bp::object & bin_by, const bp::object & signal,
                 const bp::object & weight, bp::object & binned_sig,
                 bp::object & binned_sig_sigma, bp::object & bin_counts,
                 const bp::object & bin_edges, const T lower, const T upper,
                 const int* flags_data=nullptr,
                 const std::vector<int> & flags_shape={0,0},
                 const int flags_stride=0)
{
    // signal
    BufferWrapper<T> signal_buf  ("signal",  signal,  false, std::vector<int>{-1, -1});
    if (signal_buf->strides[1] != signal_buf->itemsize)
        throw ValueError_exception("Argument 'signal' must be contiguous in last axis.");
    const int ndets = signal_buf->shape[0];
    const int nsamps = signal_buf->shape[1];
    T* signal_data = (T*)signal_buf->buf;

    // Check flags shape
    bool is_flags_2d = false;
    if (flags_data) {
        if (flags_shape.size() == 2) {
            is_flags_2d = true;
            if (flags_shape[0] != ndets || flags_shape[1] != nsamps) {
                throw ValueError_exception("2D 'flags' array has incorrect shape");
            }
        }
        else if (flags_shape[0] != nsamps) {
            throw ValueError_exception("1D 'flags' array has incorrect shape");
        }
    }

    // bin_by
    BufferWrapper<T> bin_by_buf  ("bin_by",  bin_by,  false, std::vector<int>{nsamps});
    if (bin_by_buf->strides[0] != bin_by_buf->itemsize)
        throw ValueError_exception("Argument 'bin_by' must be a C-contiguous 1d array");
    T* bin_by_data = (T*)bin_by_buf->buf;

    // weight
    BufferWrapper<T> weight_buf_temp  ("weight",  weight, true);
    if (weight_buf_temp->ndim == 1 && weight_buf_temp->strides[0] != weight_buf_temp->itemsize)
        throw ValueError_exception("Argument 'weight' must be a C-contiguous 1d array");
    else if (weight_buf_temp->ndim == 2 && weight_buf_temp->strides[1] != weight_buf_temp->itemsize)
        throw ValueError_exception("Argument 'weight' must be contiguous in last axis.");

    // Check weight dimensions
    bool is_weight_2d = false;
    std::vector<int> weight_dims;
    if (weight_buf_temp->ndim == 2) {
        weight_dims.push_back(ndets);
        is_weight_2d = true;
    }
    if (weight_buf_temp->ndim >= 1) {
        weight_dims.push_back(nsamps);
    }

    BufferWrapper<T> weight_buf  ("weight",  weight, false, weight_dims);
    T* weight_data = (T*)weight_buf->buf;

    // bin_edges
    BufferWrapper<T> bin_edges_buf  ("bin_edges",  bin_edges,  false, std::vector<int>{-1});
    if (bin_edges_buf->strides[0] != bin_edges_buf->itemsize)
        throw ValueError_exception("Argument 'bin_edges' must be a C-contiguous 1d array");
    const int nbins = bin_edges_buf->shape[0] - 1;
    T* bin_edges_data = (T*)bin_edges_buf->buf;

    // binned_sig
    BufferWrapper<T> binned_sig_buf  ("binned_sig",  binned_sig,  false, std::vector<int>{ndets, nbins});
    if (binned_sig_buf->strides[1] != binned_sig_buf->itemsize)
        throw ValueError_exception("Argument 'binned_sig' must be contiguous in last axis.");
    T* binned_sig_data = (T*)binned_sig_buf->buf;

    // binned_sig_sigma
    BufferWrapper<T> binned_sig_sigma_buf  ("binned_sig_sigma",  binned_sig_sigma,  false, std::vector<int>{ndets, nbins});
    if (binned_sig_sigma_buf->strides[1] != binned_sig_sigma_buf->itemsize)
        throw ValueError_exception("Argument 'binned_sig_sigma' must be contiguous in last axis.");
    T* binned_sig_sigma_data = (T*)binned_sig_sigma_buf->buf;

    // bin_counts
    BufferWrapper<int> bin_counts_buf  ("bin_counts",  bin_counts,  false, std::vector<int>{ndets, nbins});
    if (bin_counts_buf->strides[1] != bin_counts_buf->itemsize)
        throw ValueError_exception("Argument 'bin_counts' must be contiguous in last axis.");
    int* bin_counts_data = (int*)bin_counts_buf->buf;

    // Strides
    int signal_stride = signal_buf->strides[0] / sizeof(T);
    int weight_stride = 0;
    if (is_weight_2d) {
        weight_stride = weight_buf->strides[0] / sizeof(T);
    }
    int binned_sig_stride = binned_sig_buf->strides[0] / sizeof(T);
    int binned_sig_sigma_stride = binned_sig_sigma_buf->strides[0] / sizeof(T);
    int bin_counts_stride = bin_counts_buf->strides[0] / sizeof(int);

    // Map from data column to bin index
    int* bin_indices = (int*) malloc(nsamps * sizeof(int));
    for (int i = 0; i < nsamps; ++i) {
        bin_indices[i] = _find_bin_index(bin_edges_data, bin_by_data[i], nbins);
    }

    if (!flags_data) {
        for (int i = 0; i < nbins; ++i) {
            bin_counts_data[i] = 0;
        }
        for (int i = 0; i < nsamps; ++i) {
            if (bin_by[i] < lower || bin_by[i] > upper)
                continue;
            int bin = bin_indices[i];
            bin_counts_data[bin] += weight_data[i];
        }

        // Set all other detectors to first det bins if no flags
        #pragma omp parallel for
        for (int i = 1; i < ndets; ++i) {
            int binned_ioff = i * bin_counts_stride;
            int* bin_counts_row = bin_counts_data + binned_ioff;
            for (int j = 0; j < nbins; ++j) {
                bin_counts_row[j] = bin_counts_data[j];
            }
        }
    }

    T* binned_sig_sq_mean = (T*) malloc(nbins * ndets * sizeof(T));

    #pragma omp parallel for
    for (int i = 0; i < ndets; ++i) {
        int ioff = i * signal_stride;
        int binned_ioff = i * bin_counts_stride;
        int weight_ioff = i * weight_stride;

        T* signal_row = signal_data + ioff;
        T* binned_sig_row = binned_sig_data + (i * binned_sig_stride);
        T* binned_sig_sq_mean_row = binned_sig_sq_mean + (i * nbins);
        T* binned_sig_sigma_row = binned_sig_sigma_data + (i * binned_sig_sigma_stride);
        T* weight_row = weight_data + weight_stride;
        int* bin_counts_row = bin_counts_data + binned_ioff;

        // Zero out binned data
        for (int j = 0; j < nbins; ++j) {
            binned_sig_row[j] = 0;
            binned_sig_sq_mean_row[j] = 0;

            if (flags_data) {
                bin_counts_row[j] = 0;
            }
        }

        // Populate binned data
        for (int j = 0; j < nsamps; ++j) {
            bool samp_flagged = false;
            if (flags_data) {
                int flags_ioff = i * flags_stride;
                samp_flagged = flags_data[flags_ioff + j];
            }
            if (!samp_flagged) {
                int bin = bin_indices[j];
                if (flags_data) {
                    bin_counts_row[bin] += weight_row[j];
                }

                if (bin_counts_row[bin] > 0) {
                    T ws = weight_row[j] * signal_row[j];
                    binned_sig_row[bin] += ws;
                    binned_sig_sq_mean_row[bin] += ws * ws;
                }
            }
        }
        // Normalize
        for (int j = 0; j < nbins; ++j) {
            if (bin_counts_row[j] > 0) {
                binned_sig_row[j] /= bin_counts_row[j];
                binned_sig_sq_mean_row[j] /= bin_counts_row[j];
                binned_sig_sigma_row[j] =
                    std::sqrt(std::abs(binned_sig_sq_mean_row[j] -
                                       binned_sig_row[j] * binned_sig_row[j])) /
                                       std::sqrt(bin_counts_row[j]);
            }
            else {
                binned_sig_row[j] = std::numeric_limits<T>::quiet_NaN();
                binned_sig_sigma_row[j] = std::numeric_limits<T>::quiet_NaN();
            }
        }
    }

    free(bin_indices);
    free(binned_sig_sq_mean);
}

void bin_signal(const bp::object & bin_by, const bp::object & signal,
                const bp::object & weight, bp::object & binned_sig,
                bp::object & binned_sig_sigma, bp::object & bin_counts,
                const bp::object & bin_edges, const double lower,
                const double upper)
{
    // Get data type
    int dtype = get_dtype(signal);

    if (dtype == NPY_FLOAT) {
        _bin_signal<float>(bin_by, signal, weight, binned_sig, binned_sig_sigma,
                           bin_counts, bin_edges, (float)lower, (float)upper);
    }
    else if (dtype == NPY_DOUBLE) {
        _bin_signal<double>(bin_by, signal, weight, binned_sig, binned_sig_sigma,
                            bin_counts, bin_edges, (double)lower, (double)upper);
    }
    else {
        throw TypeError_exception("Only float32 or float64 arrays are supported.");
    }
}

void bin_flagged_signal(const bp::object & bin_by, const bp::object & signal,
                        const bp::object & weight, bp::object & binned_sig,
                        bp::object & binned_sig_sigma, bp::object & bin_counts,
                        const bp::object & bin_edges, const double lower,
                        const double upper, const bp::object & flags)
{
    // Get data type
    int dtype = get_dtype(signal);

    // flags
    BufferWrapper<int> flags_buf  ("flags",  flags,  false);
    if (flags_buf->ndim == 1 && flags_buf->strides[0] != flags_buf->itemsize)
        throw ValueError_exception("Argument 'flags' must be a C-contiguous 1d array");
    else if (flags_buf->ndim == 2 && flags_buf->strides[1] != flags_buf->itemsize)
        throw ValueError_exception("Argument 'flags' must be contiguous in last axis.");

    std::vector<int> flags_shape;
    flags_shape.push_back(flags_buf->shape[0]);

    int flags_stride = 0;
    if (flags_buf->ndim == 2) {
        flags_shape.push_back(flags_buf->shape[1]);
        flags_stride = flags_buf->strides[0] / sizeof(int);
    }

    int* flags_data = (int*)flags_buf->buf;

    if (dtype == NPY_FLOAT) {
        _bin_signal<float>(bin_by, signal, weight, binned_sig, binned_sig_sigma,
                           bin_counts, bin_edges, (float)lower, (float)upper,
                           flags_data, flags_shape, flags_stride);
    }
    else if (dtype == NPY_DOUBLE) {
        _bin_signal<double>(bin_by, signal, weight, binned_sig, binned_sig_sigma,
                            bin_counts, bin_edges, (double)lower, (double)upper,
                            flags_data, flags_shape, flags_stride);
    }
    else {
        throw TypeError_exception("Only float32 or float64 arrays are supported.");
    }
}


PYBINDINGS("so3g")
{
    bp::def("nmat_detvecs_apply", nmat_detvecs_apply);
    bp::def("process_cuts",  process_cuts);
    bp::def("translate_cuts", translate_cuts);
    bp::def("get_gap_fill_poly",  get_gap_fill_poly<float>,
            "get_gap_fill_poly(ranges, signal, buffer, order, extract)\n"
            "\n"
            "Do polynomial gap-filling on a float32 array.\n"
            "\n"
            "Args:\n"
            "  ranges: RangesMatrix with shape (ndet, nsamp)\n"
            "  signal: data array (float32) with shape (ndet, nsamp)\n"
            "  buffer: integer stating max number of samples to use on each end\n"
            "  order: order of polynomial to use (1 means linear)\n"
            "  inplace: whether to overwrite data array with the model\n"
            "  extract: array to write the original data samples (inplace)\n"
            "    or the model (!inplace) into.\n");
    bp::def("get_gap_fill_poly64",  get_gap_fill_poly<double>,
            "get_gap_fill_poly64(ranges, signal, buffer, order, extract)\n"
            "\n"
            "Do polynomial gap-filling for float64 data.\n"
            "\n"
            "See details in docstring for get_gap_fill_poly.\n");
    bp::def("test_buffer_wrapper", test_buffer_wrapper,
            "Pass array and list of dims to match against its shape.");
    bp::def("block_moment", block_moment<float>,
            "block_moment(tod, out, bsize, moment, central, shift)\n"
            "\n"
            "Compute the nth moment in blocks on a float32 array.\n"
            "\n"
            "Args:\n"
            "  tod: data array (float32) with shape (ndet, nsamp)\n"
            "  out: output array (float32) with shape (ndet, nsamp)\n"
            "       can be the same as tod\n"
            "  bsize: number of samples in each block\n"
            "  moment: moment to compute, should be >= 1\n"
            "  central: whether to compute the central moment in each block\n"
            "  shift: sample to start block at, used in each row\n");
    bp::def("block_moment64", block_moment<double>,
            "block_moment64(tod, out, bsize, moment, central, shift)\n"
            "\n"
            "Compute the nth moment in blocks on a float32 array.\n"
            "\n"
            "See details in docstring for block_moment.\n");
    bp::def("block_minmax", block_minmax<float>,
            "block_minmax(tod, out, bsize, mode, shift)\n"
            "\n"
            "Compute the minimum, maximum, or peak to peak in blocks on a float32 array.\n"
            "\n"
            "Args:\n"
            "  tod: data array (float32) with shape (ndet, nsamp)\n"
            "  out: output array (float32) with shape (ndet, nsamp)\n"
            "       can be the same as tod\n"
            "  bsize: number of samples in each block\n"
            "  mode: if 0 compute the block minimum, if 1 the maximum, anything else will compute the peak to peak\n"
            "  shift: sample to start block at, used in each row\n");
    bp::def("block_minmax64", block_minmax<double>,
            "block_minmax64(tod, out, bsize, mode, shift)\n"
            "\n"
            "Compute the minimum, maximum, or peak to peak in blocks on a float64 array.\n"
            "\n"
            "See details in docstring for block_minmax.\n");
    bp::def("matched_jumps", matched_jumps<float>,
            "matched_jumps(tod, out, min_size, bsize)\n"
            "\n"
            "Flag jumps with the matched filter for a unit jump in a float32 array.\n"
            "\n"
            "Args:\n"
            "  tod: data array (float32) with shape (ndet, nsamp)\n"
            "  out: output array (int32) with shape (ndet, nsamp)\n"
            "  min_size: minimum jump size for each det, shape (ndet,)\n"
            "  bsize: number of samples in each block\n");
    bp::def("matched_jumps64", matched_jumps<double>,
            "matched_jumps64(tod, out, min_size, bsize)\n"
            "\n"
            "Flag jumps with the matched filter for a unit jump in a float64 array.\n"
            "\n"
            "See details in docstring for matched_jumps.\n");
    bp::def("find_quantized_jumps", find_quantized_jumps<float>,
            "find_quantized_jumps(tod, out, atol, win_size, scale)"
            "\n"
            "Search for jumps that are a multiple of a known value in a float32 array.\n"
            "Output will be 0 where jumps are not found and the assumed jump height where jumps are found.\n"
            "\n"
            "Args:\n"
            "  tod: data array (float32) with shape (ndet, nsamp)\n"
            "  out: output array (float32) with shape (ndet, nsamp)\n"
            "  atol: how close to the multiple of scale a value needs to be to be a jump in the same units as the signal\n"
            "        should be an array (float32) with shape (ndet,)\n"
            "  win_size: size of window to use as buffer when differencing\n"
            "  scale: the scale of jumps to look for\n");
    bp::def("find_quantized_jumps64", find_quantized_jumps<double>,
            "find_quantized_jumps64(tod, out, atol, win_size, scale)\n"
            "\n"
            "Search for jumps that are a multiple of a known value in a float64 array.\n"
            "Output will be 0 where jumps are not found and the assumed jump height where jumps are found.\n"
            "\n"
            "See details in docstring for find_quantized_jumps.\n");
    bp::def("subtract_jump_heights", subtract_jump_heights<float>,
            "subtract_jump_heights(tod, out, heights, jumps)"
            "\n"
            "For each detector, compute the cumulative effect of the jumps identified by the array 'heights' and the RangesMatrix 'jumps'."
            "For each range in 'jumps', the values from 'heights' are checked and the size of the jump is either the largest positive"
            "or the largest negative number (whichever has the largest absolute value)."
            "The 'output' value is the difference of 'tod' and the accumulated jump vector."
            "\n"
            "Args:\n"
            "  tod: data array (float32) with shape (ndet, nsamp)\n"
            "  out: output array (float32) with shape (ndet, nsamp)\n"
            "       can be the same as tod\n"
            "  heights: the height of the jump at each samples\n"
            "           should be an array (float32) with shape (ndet, nsamp)\n"
            "  jumps: RangesMatrix with the jump locations and shape (ndet, nsamp).\n");
    bp::def("subtract_jump_heights64", subtract_jump_heights<double>,
            "subtract_jump_heights64(tod, out, heights, jumps)"
            "\n"
            "Subtract fit jump heights from known jump locatations in a float64 array."
            "If multiple samples in a jump have different heights, the largest height is used.\n"
            "\n"
            "See details in docstring for subtract_jump_heights.\n");
    bp::def("clean_flag", clean_flag,
            "clean_flag(flag, width)"
            "\n"
            "Clean a flag inplace by unflagging regions without enough contiguous flagged values.\n"
            "\n"
            "Args:\n"
            "  flag: flag array (int) with shape (ndet, nsamp)\n"
            "  width: the minimum number of contiguous flagged samples\n");
    bp::def("interp1d_linear", interp1d_linear,
            "interp1d_linear(x, y, x_interp, y_interp)"
            "\n"
            "Perform linear interpolation over rows of float32 or float64 array with GSL.\n"
            "This function uses OMP to parallelize over the dets (rows) axis.\n"
            "Vector x must be strictly increasing. Values for x_interp beyond the "
            "domain of x will be computed based on extrapolation."
            "\n"
            "Args:\n"
            "  x: independent variable (float32/float64) of data with shape (nsamp,)\n"
            "  y: data array (float32/float64) with shape (ndet, nsamp)\n"
            "  x_interp: independent variable (float32/float64) for interpolated data "
            "            with shape (nsamp_interp,)\n"
            "  y_interp: interpolated data array (float32/float64) output buffer to be modified "
            "            with shape (ndet, nsamp_interp)\n");
    bp::def("detrend", detrend,
            "detrend(tod, method, ncount)"
            "\n"
            "Detrend each row of an array (float32/float64). This function uses OMP to parallelize "
            "over the dets (rows) axis."
            "\n"
            "Args:\n"
            "  tod: input array (float32/float64) buffer with shape (ndet, nsamp) that is to be detrended. "
            "       The data is modified in place.\n"
            "  method: how to detrend data.  Options are 'mean', 'median', and 'linear'. Linear calculates "
            "          and subtracts the slope from either end of each row as determined from 'linear_ncount'.\n"
            "  linear_ncount: Number (int) of samples to use on each end, when measuring mean level for 'linear'"
            "                 detrend. Must be a positive integer or -1.  If -1, nsamps / 2 will be used. Values "
            "                 larger than 1 suppress the influence of white noise.\n");
    bp::def("bin_signal", bin_signal,
            "bin_signal(bin_by, signal, weight, binned_sig, binned_sig_sigma, bin_counts, bin_edges, lower, upper)"
            "\n"
            "Bin time-ordered data by ``bin_by`` and return the binned signal and its standard deviation.\n"
            "This function uses OMP to parallelize over the dets (rows) axis. Supports unequal bin widths.\n"
            "Args:\n"
            "  bin_by: the array (float32/float64) by which signal is binned with shape (nsamp)\n"
            "  signal: the signal array (float32/float64) to be binned with shape (ndet,nsamp)\n"
            "  weight: array (float32/float64) of weights for the signal values.  May have shapes\n"
            "          of (nsamps) or (ndets, nsamps)\n"
            "  binned_sig: binned signal array (float32/float64) with shape (ndet,nbin).\n"
            "              Modified in place.\n"
            "  binned_sig_sigma: estimated sigma of binned signal (float32/float64) with shape (ndet,nbin).\n"
            "                    Modified in place."
            "  bin_counts: counts of binned samples (int32) with shape (ndet,nbin).  Modified in place.\n"
            "  bin_edges: array (float32/float64) of bin edges with length=nbins+1.  Must be monotonically increasing but\n"
            "             but may have different widths.\n"
            "  lower: lower bin range (float64)\n"
            "  upper: upper bin range (float64)\n");
    bp::def("bin_flagged_signal", bin_flagged_signal,
            "bin_signal(bin_by, signal, weight, binned_sig, binned_sig_sigma, bin_counts, bin_edges, lower, upper, flags)\n"
            "\n"
            "Bin time-ordered data by ``bin_by`` and return the binned signal and its standard deviation.\n"
            "This function uses OMP to parallelize over the dets (rows) axis. Supports unequal bin widths.\n"
            "Args:\n"
            "  bin_by: the array (float32/float64) by which signal is binned with shape (nsamp)\n"
            "  signal: the signal array (float32/float64) to be binned with shape (ndet,nsamp)\n"
            "  weight: array (float32/float64) of weights for the signal values.  May have shapes\n"
            "          of (nsamps) or (ndets, nsamps)\n"
            "  binned_sig: binned signal array (float32/float64) with shape (ndet,nbin).\n"
            "              Modified in place.\n"
            "  binned_sig_sigma: estimated sigma of binned signal (float32/float64) with shape (ndet,nbin).\n"
            "                    Modified in place.\n"
            "  bin_counts: counts of binned samples (int32) with shape (ndet,nbin).  Modified in place.\n"
            "  bin_edges: array (float32/float64) of bin edges with length=nbins+1.  Must be monotonically increasing but\n"
            "             but may have different widths.\n"
            "  lower: lower bin range. Data points falling outside this range will be ignored. (float64)\n"
            "  upper: upper bin range. Data points falling outside this range will be ignored.  (float64)\n"
            "  flags: array (int32) indicating whether to exclude flagged samples when binning the signal.\n"
            "         Can be of shape (nsamp) or (ndet,nsamp).\n");
}
