#define NO_IMPORT_ARRAY

#include <complex>
#include <stdint.h>
#include <stdio.h>
#include <vector>
#include <string>
extern "C" {
    #include <cblas.h>
    // Additional prototypes for Fortran LAPACK routines.
    // dposv: solve Ax = b for A positive definite.
    void dposv_(const char* uplo, int* n, int* nrhs, double* a, int* lda,
                double* b, int* ldb, int* info );
}

#include <boost/python.hpp>

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
    Py_DECREF(a);
    return PyArray_TYPE(a);
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
}
