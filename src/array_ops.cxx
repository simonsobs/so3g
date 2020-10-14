#define NO_IMPORT_ARRAY

#include <pybindings.h>
#include <complex>
#include <stdint.h>
#include <stdio.h>
#include "array_ops.h"
#include "numpy_assist.h"
#include <vector>
#include <cblas.h>

#if 0
enum {CblasRowMajor=101, CblasColMajor=102};
enum {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113};
void cblas_sgemm(int layout, int TransA, int TransB, const int M, const int N, const int K, const float alpha, const float * A, const int lda, const float * B, const int ldb, const float beta, float * C, const int ldc );
#endif

template <typename T>
const T & min(const T & a, const T & b) { return a < b ? a : b; }

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

PYBINDINGS("so3g")
{
	bp::def("nmat_detvecs_apply", nmat_detvecs_apply);
}
