#pragma once

#include <ceres/ceres.h>

template <int Degree>
struct PolynomialModel
{
    // Ceres requires number of params
    // to be known at compile time
    static constexpr int nparams = Degree + 1;

    template <typename T>
    static T eval(T x, const T* params)
    {
        const T p0 = params[0];
        T result = p0;
        for (int i = 1; i < nparams; ++i) {
            const T p = params[i];
            result += p * ceres::pow(x, T(i));
        }

        return result;
    }
    // Not needed for least squares as ceres 
    // supports boundaries
    template <typename T>
    static bool check_bounds(const T* params)
    {
        return true;
    }
};

struct NoiseModel
{
    // Ceres requires number of params
    // to be known at compile time
    static constexpr int nparams = 3;

    template <typename T>
    static T eval(T f, const T* params)
    {
        const T fknee = params[0];
        const T w = params[1];
        const T alpha = params[2];

        return w * (1.0 + ceres::pow(fknee / f, alpha));
    }

    // Slightly hacky way of bounds checking but is
    // suggested by Ceres to ensure it never goes
    // out of bounds
    template <typename T>
    static bool check_bounds(const T* params)
    {
        const T w = params[1];
        if (w <= 0.0) {
            return false;
        }
        return true;
    }
};

// Model independent cost function for least-squares fitting
template <typename Model>
struct CostFunction
{
    using model = Model;

    CostFunction(int n, const double* x_data, const double* y_data)
        : n_pts(n), x(x_data), y(y_data) {}

    template <typename T>
    bool operator()(const T* const params, T* residual) const {
        for (int i = 0; i < n_pts; ++i) {
            T model = Model::eval(T(x[i]), params);
            residual[i] = T(y[i]) - model;
        }
        return true;
    }

    static ceres::Problem create(const int n, const double* xx,
                                 const double* yy, double* p)
    {
        ceres::Problem problem;

        problem.AddResidualBlock(
            new ceres::AutoDiffCostFunction<CostFunction<Model>, 
            ceres::DYNAMIC, Model::nparams>(
                new CostFunction<Model>(n, xx, yy), n), nullptr, p);

        return problem;
    }

private:
    const int n_pts;
    const double* x;
    const double* y;
};

// Model independent Negative Log Likelihood for generalized 
// unconstrained minimization
template <typename Model>
struct NegLogLikelihood
{
    using model = Model;

    NegLogLikelihood(int n, const double* x_data, const double* y_data)
        : n_pts(n), x(x_data), y(y_data) {}

    template <typename T>
    bool operator()(const T* const params, T* cost) const
    {
        // Check bounds (saves a lot of time)
        if (!model::check_bounds(params)) {
            return false;
        }

        cost[0] = T(0.);
        for (int i = 0; i < n_pts; ++i) {
            T model = Model::eval(T(x[i]), params);
            cost[0] += ceres::log(model) + T(y[i]) / model;
        }

        return true;
    }

    static ceres::FirstOrderFunction* create(int n, const double* xx,
                                             const double* yy)
    {
        // Ceres takes ownership of pointers so no cleanup is required
        return new ceres::AutoDiffFirstOrderFunction<NegLogLikelihood<Model>,
            Model::nparams>(new NegLogLikelihood<Model>(n, xx, yy));
    }

private:
    const int n_pts;
    const double* x;
    const double* y;
};
