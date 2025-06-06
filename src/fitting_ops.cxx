#define NO_IMPORT_ARRAY
#define GLOG_USE_GLOG_EXPORT

#include <stdint.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <cmath>

#include <boost/python.hpp>
#include <omp.h>

#include <gsl/gsl_statistics.h>
#include <glog/logging.h>
#include <ceres/ceres.h>

#include <pybindings.h>
#include "so3g_numpy.h"
#include "numpy_assist.h"
#include "Ranges.h"
#include "array_ops.h"
#include "fitting_ops.h"


template <typename CostFunction>
void _least_squares(const double* x, const double* y, double* p, const int n,
                    const int nthreads=1, const double* lower_bounds=nullptr,
                    const double* upper_bounds=nullptr, double* c=nullptr)
{
    // Set up ceres problem
    ceres::Problem problem = CostFunction::create(n, x, y, p);

    ceres::Solver::Options options;
    options.linear_solver_type = ceres::LinearSolverType::DENSE_QR;
    options.logging_type = ceres::LoggingType::SILENT;
    options.minimizer_progress_to_stdout = false;
    options.num_threads = nthreads;

    // Set bounds
    if (lower_bounds) {
        for (int i = 0; i < CostFunction::model::nparams; ++i) {
            if (!std::isnan(lower_bounds[i])) {
                problem.SetParameterLowerBound(p, i, lower_bounds[i]);
            }
        }
    }
    if (upper_bounds) {
        for (int i = 0; i < CostFunction::model::nparams; ++i) {
            if (!std::isnan(upper_bounds[i])) {
                problem.SetParameterUpperBound(p, i, upper_bounds[i]);
            }
        }
    }

    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    if (summary.IsSolutionUsable()) {
        if (c) {
            ceres::Covariance::Options covariance_options;
            covariance_options.sparse_linear_algebra_library_type =
                ceres::SparseLinearAlgebraLibraryType::EIGEN_SPARSE;
            covariance_options.algorithm_type =
                ceres::CovarianceAlgorithmType::DENSE_SVD;
            covariance_options.null_space_rank = -1;

            ceres::Covariance covariance(covariance_options);
            std::vector<std::pair<const double*, const double*>> covariance_blocks;
            covariance_blocks.emplace_back(p, p);
            if (covariance.Compute(covariance_blocks, &problem)) {
                covariance.GetCovarianceBlock(p, p, c);
            }
            else {
                for (int i = 0; i < CostFunction::model::nparams; ++i) {
                    p[i] = std::numeric_limits<double>::quiet_NaN();
                    c[i] = std::numeric_limits<double>::quiet_NaN();
                }
            }
        }
    }
    else {
        for (int i = 0; i < CostFunction::model::nparams; ++i) {
            p[i] = std::numeric_limits<double>::quiet_NaN();
            if (c) {
                c[i] = std::numeric_limits<double>::quiet_NaN();
            }
        }
    }
}

template <typename T>
bool _invert_matrix(const T* matrix, T* inverse, const int n) {
    // Calculate inverse with Gaussian elimination (faster than GSL
    // for small numbers of parameters)
    std::unique_ptr<T[]> aug_matrix = std::make_unique<T[]>(n * 2 * n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            aug_matrix[i * 2 * n + j] = matrix[i * n + j];
        }
        for (int j = 0; j < n; ++j) {
            aug_matrix[i * 2 * n + (n + j)] = (i == j) ? 1.0 : 0.0;
        }
    }

    for (int i = 0; i < n; ++i) {
        T pivot = aug_matrix[i * 2 * n + i];
        if (pivot == 0.0) {
            return false;
        }

        for (int j = 0; j < 2 * n; ++j) {
            aug_matrix[i * 2 * n + j] /= pivot;
        }

        for (int k = 0; k < n; ++k) {
            if (k == i) continue;
            T factor = aug_matrix[k * 2 * n + i];
            for (int j = 0; j < 2 * n; ++j) {
                aug_matrix[k * 2 * n + j] -= factor * aug_matrix[i * 2 * n + j];
            }
        }
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            inverse[i * n + j] = aug_matrix[i * 2 * n + (n + j)];
        }
    }

    return true;
}

// Compute the Hessian matrix using finite differences
template <typename Likelihood>
void _compute_hessian(ceres::GradientProblem& problem, double* p,
                      double* hessian, const double epsilon=1e-5) {

    double gradient_plus[Likelihood::model::nparams];
    double gradient_minus[Likelihood::model::nparams];
    double perturbed_parameters[Likelihood::model::nparams];
    double cost;

    // Loop over each parameter
    for (int i = 0; i < Likelihood::model::nparams; ++i) {
        // Perturb in the positive direction
        std::copy(p, p + Likelihood::model::nparams, perturbed_parameters);
        perturbed_parameters[i] += epsilon;
        problem.Evaluate(perturbed_parameters, &cost, gradient_plus);

        // Perturb in the negative direction
        std::copy(p, p + Likelihood::model::nparams, perturbed_parameters);
        perturbed_parameters[i] -= epsilon;
        problem.Evaluate(perturbed_parameters, &cost, gradient_minus);

        // Second derivative approximation
        for (int j = 0; j < Likelihood::model::nparams; ++j) {
            hessian[i * Likelihood::model::nparams + j] =
                (gradient_plus[j] - gradient_minus[j]) / (2 * epsilon);
        }
    }
}

// General minimization function
template <typename Likelihood>
void _minimize(const double* x, const double* y, double* p, const int n,
               const double tol, const int niter, double* c=nullptr,
               const double epsilon=1e-5)
{
    // Create problem
    ceres::GradientProblem problem(Likelihood::create(n, x, y));

    // Options for solving
    ceres::GradientProblemSolver::Options options;
    options.logging_type = ceres::LoggingType::SILENT;
    options.minimizer_progress_to_stdout = false;
    options.line_search_direction_type = ceres::BFGS;
    options.max_num_iterations = niter;
    options.function_tolerance = tol;

    ceres::GradientProblemSolver::Summary summary;

    // Run the solver
    ceres::Solve(options, problem, p, &summary);

    // Calculate uncertainties.  Ceres does not include a way to calculate
    // uncertainties for gradient problems but this seems fast enough.
    if (summary.IsSolutionUsable()) {
        if (c) {
            double hessian[Likelihood::model::nparams * Likelihood::model::nparams];
            double inv_hessian[Likelihood::model::nparams * Likelihood::model::nparams];

            _compute_hessian<Likelihood>(problem, p, hessian, epsilon);
            bool not_singular = _invert_matrix(hessian, inv_hessian,
                                               Likelihood::model::nparams);

            if (not_singular) {
                for (int i = 0; i < Likelihood::model::nparams; ++i) {
                    c[i] = std::sqrt(inv_hessian[i * Likelihood::model::nparams + i]);
                }
            }
            else {
                for (int i = 0; i < Likelihood::model::nparams; ++i) {
                    c[i] = std::numeric_limits<double>::quiet_NaN();
                }
            }
        }
    }
    else {
        for (int i = 0; i < Likelihood::model::nparams; ++i) {
            p[i] = std::numeric_limits<double>::quiet_NaN();
            if (c) {
                c[i] = std::numeric_limits<double>::quiet_NaN();
            }
        }
    }
}

// Get indices in an array corresponding to input values.  Assumes both
// are monotonically increasing.
auto _get_array_indices(const double* x, const std::vector<double>& vals,
                        const int nsamps)
{
    std::vector<int> indices(vals.size(), -1);

    int j = 0;
    for (int i = 0; i < vals.size(); ++i) {
        while (j < nsamps && x[j] <= vals[i]) {
            ++j;
        }
        indices[i] = j - 1;
    }

    return indices;
}

// Get indices corresponding to lower freq and white noise
// limits.
auto _get_frequency_limits(const double* f, const double lowf,
                           const double fwhite_lower,
                           const double fwhite_upper,
                           const int nsamps)
{
    if (fwhite_lower < lowf) {
        throw ValueError_exception("fwhite lower < lower freq.");
    }

    if (fwhite_lower >= fwhite_upper) {
        throw ValueError_exception("fwhite lower >= fwhite upper.");
    }

    std::vector<int> f_indx = _get_array_indices(f, {lowf, fwhite_lower,
                                                     fwhite_upper}, nsamps);
    int fwhite_size = f_indx[2] - f_indx[1] + 1;

    return std::make_tuple(f_indx[0], std::vector<int>{f_indx[1], f_indx[2]},
                           fwhite_size);
}

template <typename CostFunc, typename Likelihood, typename T>
void _fit_noise(const double* f, const double* log_f, const double* pxx,
                T* p,  T* c, const int ndets, const int nsamps,
                const int lowf_i, const std::vector<int> fwhite_i,
                const int fwhite_size, const double tol, const int niter,
                const double epsilon)
{
    double log_pxx[nsamps];
    for (int i = 0; i < nsamps; ++i) {
        log_pxx[i] = std::log10(pxx[i]);
    }

    // Estimate of white noise
    double wnest = _calculate_median(pxx + fwhite_i[0], fwhite_size);

    // Fit 1/f to line in logspace
    double pfit[2] = {1., -1.};

    // Some weak bounds (w > 0, alpha < 0)
    double lower_bounds[2] = {0., std::numeric_limits<double>::quiet_NaN()};
    double upper_bounds[2] = {std::numeric_limits<double>::quiet_NaN(), 0.};

    _least_squares<CostFunc>(log_f, log_pxx, pfit, lowf_i + 1, 1,
                             lower_bounds, upper_bounds);

    // Don't run minimization if least squares fit failed
    if (std::isnan(pfit[0]) || std::isnan(pfit[1])) {
        // Implict cast from double to T
        for (int i = 0; i < Likelihood::model::nparams; ++i) {
            p[i] = std::numeric_limits<double>::quiet_NaN();
            c[i] = std::numeric_limits<double>::quiet_NaN();
        }
    }
    else {
        // Find guess for fknee
        int fi = 0;
        double min_f = std::numeric_limits<double>::max();
        for (int i = 0; i < nsamps; ++i) {
            double model_val = CostFunc::model::eval(log_f[i], pfit);
            double diff = std::abs(std::pow(10., model_val) - wnest);
            if (diff <= min_f) {
                min_f = diff;
                fi = i;
            }
        }

        // Initial parameters
        double p0[Likelihood::model::nparams] = {f[fi], wnest, -pfit[1]};
        double c0[Likelihood::model::nparams];

        // Minimize
        _minimize<Likelihood>(f, pxx, p0, nsamps, tol, niter, c0, epsilon);

        // Implict cast from double to T
        for (int i = 0; i < Likelihood::model::nparams; ++i) {
            p[i] = p0[i];
            c[i] = c0[i];
        }
    }
}

template <typename T>
void _fit_noise_buffer(const bp::object & f, const bp::object & pxx,
                       bp::object & p, bp::object & c, const double lowf,
                       const double fwhite_lower, const double fwhite_upper,
                       const double tol, const int niter, const double epsilon)
{
    using Likelihood = NegLogLikelihood<NoiseModel>;
    using CostFunc = CostFunction<PolynomialModel<1>>;

    // Disable google logging to prevent failed fit warnings for each detector
    google::InitGoogleLogging("logger");
    FLAGS_minloglevel = google::FATAL + 1;
    FLAGS_logtostderr = false;

    BufferWrapper<T> pxx_buf  ("pxx",  pxx,  false, std::vector<int>{-1, -1});
    if (pxx_buf->strides[1] != pxx_buf->itemsize)
        throw ValueError_exception("Argument 'pxx' must be contiguous in last axis.");
    T* pxx_data = (T*)pxx_buf->buf;
    const int ndets = pxx_buf->shape[0];
    const int nsamps = pxx_buf->shape[1];
    const int row_stride = pxx_buf->strides[0] / sizeof(T);

    BufferWrapper<T> f_buf  ("f",  f,  false, std::vector<int>{nsamps});
    if (f_buf->strides[0] != f_buf->itemsize)
        throw ValueError_exception("Argument 'f' must be a C-contiguous 1d array.");
    T* f_data = (T*)f_buf->buf;

    BufferWrapper<T> p_buf  ("p",  p,  false, std::vector<int>{ndets, Likelihood::model::nparams});
    if (p_buf->strides[1] != p_buf->itemsize)
        throw ValueError_exception("Argument 'p' must be contiguous in last axis.");
    T* p_data = (T*)p_buf->buf;
    const int p_stride = p_buf->strides[0] / sizeof(T);

    BufferWrapper<T> c_buf  ("c",  c,  false, std::vector<int>{ndets, Likelihood::model::nparams});
    if (c_buf->strides[1] != c_buf->itemsize)
        throw ValueError_exception("Argument 'c' must be contiguous in last axis.");
    T* c_data = (T*)c_buf->buf;
    const int c_stride = c_buf->strides[0] / sizeof(T);

    if constexpr (std::is_same<T, float>::value) {
        // Copy f to double
        double f_double[nsamps];

        std::transform(f_data, f_data + nsamps, f_double,
                       [](float value) { return static_cast<double>(value); });

        // Get frequency bounds
        auto [lowf_i, fwhite_i, fwhite_size] =
            _get_frequency_limits(f_double, lowf, fwhite_lower, fwhite_upper, nsamps);

        // Fit in logspace
        double log_f[nsamps];
        for (int i = 0; i < nsamps; ++i) {
            log_f[i] = std::log10(f_double[i]);
        }

        #pragma omp parallel for
        for (int i = 0; i < ndets; ++i) {
            int ioff = i * row_stride;
            int poff = i * p_stride;
            int coff = i * c_stride;

            // Copy pxx row to double
            double pxx_det[nsamps];

            std::transform(pxx_data + ioff, pxx_data + ioff + nsamps, pxx_det,
                       [](float value) { return static_cast<double>(value); });

            // Cast implicitly on assignment
            T* p_det = p_data + poff;
            T* c_det = c_data + coff;

            _fit_noise<CostFunc, Likelihood>(f_double, log_f, pxx_det, p_det,
                                             c_det, ndets, nsamps, lowf_i, fwhite_i,
                                             fwhite_size, tol, niter, epsilon);
        }
    }
    else if constexpr (std::is_same<T, double>::value) {
        // Get frequency bounds
        auto [lowf_i, fwhite_i, fwhite_size] =
            _get_frequency_limits(f_data, lowf, fwhite_lower, fwhite_upper, nsamps);

        // Fit in logspace
        double log_f[nsamps];
        for (int i = 0; i < nsamps; ++i) {
            log_f[i] = std::log10(f_data[i]);
        }

        #pragma omp parallel for
        for (int i = 0; i < ndets; ++i) {
            int ioff = i * row_stride;
            int poff = i * p_stride;
            int coff = i * c_stride;

            double* pxx_det = pxx_data + ioff;
            T* p_det = p_data + poff;
            T* c_det = c_data + coff;

            _fit_noise<CostFunc, Likelihood>(f_data, log_f, pxx_det, p_det, c_det,
                                             ndets, nsamps, lowf_i, fwhite_i, fwhite_size,
                                             tol, niter, epsilon);
        }
    }
    google::ShutdownGoogleLogging();
}

void fit_noise(const bp::object & f, const bp::object & pxx, bp::object & p, bp::object & c,
               const double lowf, const double fwhite_lower, const double fwhite_upper,
               const double tol, const int niter, const double epsilon)
{
    // Get data type
    int dtype = get_dtype(pxx);

    if (dtype == NPY_FLOAT) {
        _fit_noise_buffer<float>(f, pxx, p, c, lowf, fwhite_lower, fwhite_upper, tol, niter, epsilon);
    }
    else if (dtype == NPY_DOUBLE) {
        _fit_noise_buffer<double>(f, pxx, p, c, lowf, fwhite_lower, fwhite_upper, tol, niter, epsilon);
    }
    else {
        throw TypeError_exception("Only float32 or float64 arrays are supported.");
    }
}


PYBINDINGS("so3g")
{
     bp::def("fit_noise", fit_noise,
             "fit_noise(f, pxx, p, c, lowf, fwhite_lower, fwhite_upper, tol, niter, epsilon)"
             "\n"
             "Fits noise model with white and 1/f components to the PSD of signal. Uses a MLE\n"
             "method that minimizes a log-likelihood. OMP is used to parallelize across dets (rows)."
             "\n"
             "Args:\n"
             "  f: frequency array (float32/64) with dimensions (nsamps,).\n"
             "     Should be positive definite and strictly increasing.\n"
             "  pxx: PSD array (float32/64) with dimensions (ndets, nsamps).\n"
             "  p: output parameter array (float32/64) with dimensions (ndets, nparams).\n"
             "     This is modified in place and input values are ignored.\n"
             "  c: output uncertaintiy array (float32/64) with dimensions (ndets, nparams).\n"
             "     This is modified in place and input values are ignored.\n"
             "  lowf: Frequency below which the 1/f noise index and fknee are estimated for initial\n"
             "        guess passed to least_squares fit (float64).\n"
             "  fwhite_lower: Lower frequency used to estimate white noise for initial guess passed to\n"
             "                least_squares fit (float64).  Should be < fwhite_upper and >= lowf.\n"
             "  fwhite_upper: Upper frequency used to estimate white noise for initial guess passed to\n"
             "                least_squares fit (float64).  Should be > fwhite_lower and lowf.\n"
             "  tol: absolute tolerance for minimization (float64).\n"
             "  niter: total number of iterations to run minimization for (int).\n"
             "  epsilon: Value to perturb gradients by when calculating uncertainties with the inverse\n"
             "           Hessian matrix (float64). Affects minimization only.\n");
}