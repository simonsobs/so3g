#define NO_IMPORT_ARRAY

#include <pybindings.h>
#include <complex>
#include <stdint.h>
#include <stdio.h>
#include "array_ops.h"
#include "numpy_assist.h"
#include "Ranges.h"
#include <vector>
extern "C" {
  #include <cblas.h>
}
#include <string>
#include <cblas.h>
#include "so3g_numpy.h"

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
template <typename T> int pcut_full_tod2vals_helper(const vector<RangesInt32> &, T *, int, int, T *);
template <typename T> int pcut_full_vals2tod_helper(const vector<RangesInt32> &, T *, int, int, T *);
int pcut_poly_measure_helper(const vector<RangesInt32> &, int, int nmax);
template <typename T> void pcut_poly_tod2vals_helper(const vector<RangesInt32> &, int, int, T *, int, int, T *);
template <typename T> void pcut_poly_vals2tod_helper(const vector<RangesInt32> &, int, int, T *, int, int, T *);
template <typename T> int pcut_clear_helper(const vector<RangesInt32> &, T *, int, int);

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
	} else throw general_exception("process_cuts model can only be 'full' or 'poly'");

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
					pcut_full_vals2tod_helper(ranges, (float*)tod_buf->buf, ndet, nsamp, (float*) vals_buf->buf);
				else if(model == "poly")
					pcut_poly_vals2tod_helper(ranges, resolution, nmax, (float*)tod_buf->buf, ndet, nsamp, (float*) vals_buf->buf);
			} else if(operation == "extract") {
				if     (model == "full")
					pcut_full_tod2vals_helper(ranges, (float*)tod_buf->buf, ndet, nsamp, (float*) vals_buf->buf);
				else if(model == "poly")
					pcut_poly_tod2vals_helper(ranges, resolution, nmax, (float*)tod_buf->buf, ndet, nsamp, (float*) vals_buf->buf);
			} else if(operation == "clear") {
				pcut_clear_helper(ranges, (float*)tod_buf->buf, ndet, nsamp);
			} else throw general_exception("process_cuts operation can only be 'measure', 'insert' or 'extract'");
		} else if(dtype == NPY_DOUBLE) {
			BufferWrapper<double> tod_buf  ("tod",  tod,  false, std::vector<int>{-1,-1});
			BufferWrapper<double> vals_buf ("vals", vals, false, std::vector<int>{-1});
			int ndet = tod_buf->shape[0], nsamp = tod_buf->shape[1];
			if(operation == "insert") {
				if     (model == "full")
					pcut_full_vals2tod_helper(ranges, (double*)tod_buf->buf, ndet, nsamp, (double*) vals_buf->buf);
				else if(model == "poly")
					pcut_poly_vals2tod_helper(ranges, resolution, nmax, (double*)tod_buf->buf, ndet, nsamp, (double*) vals_buf->buf);
			} else if(operation == "extract") {
				if     (model == "full")
					pcut_full_tod2vals_helper(ranges, (double*)tod_buf->buf, ndet, nsamp, (double*) vals_buf->buf);
				else if(model == "poly")
					pcut_poly_tod2vals_helper(ranges, resolution, nmax, (double*)tod_buf->buf, ndet, nsamp, (double*) vals_buf->buf);
			} else if(operation == "clear") {
				pcut_clear_helper(ranges, (float*)tod_buf->buf, ndet, nsamp);
			} else throw general_exception("process_cuts operation can only be 'measure', 'insert' or 'extract'");
		} else throw general_exception("process_cuts only supports float32 and float64");
	}
	return 0;
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
int pcut_full_tod2vals_helper(const vector<RangesInt32> & rangemat, T * tod, int ndet, int nsamp, T * vals) {
	int i = 0;
	for(int di = 0; di < rangemat.size(); di++)
		for (auto const &r: rangemat[di].segments)
			for(int j = r.first; j < r.second; j++, i++)
				vals[i] = tod[di*nsamp+j];
}
template <typename T>
int pcut_full_vals2tod_helper(const vector<RangesInt32> & rangemat, T * tod, int ndet, int nsamp, T * vals) {
	int i = 0;
	for(int di = 0; di < rangemat.size(); di++)
		for (auto const &r: rangemat[di].segments)
			for(int j = r.first; j < r.second; j++, i++)
				tod[di*nsamp+j] = vals[i];
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
void pcut_poly_tod2vals_helper(const vector<RangesInt32> & rangemat, int resolution, int nmax, T * tod, int ndet, int nsamp, T * vals) {
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
void pcut_poly_vals2tod_helper(const vector<RangesInt32> & rangemat, int resolution, int nmax, T * tod, int ndet, int nsamp, T * vals) {
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
int pcut_clear_helper(const vector<RangesInt32> & rangemat, T * tod, int ndet, int nsamp) {
	#pragma omp parallel for
	for(int di = 0; di < rangemat.size(); di++)
		for (auto const &r: rangemat[di].segments)
			for(int s = r.first; s < r.second; s++)
				tod[di*nsamp+s] = 0;
}


PYBINDINGS("so3g")
{
	bp::def("nmat_detvecs_apply", nmat_detvecs_apply);
	bp::def("process_cuts",  process_cuts);
}
