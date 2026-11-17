#pragma once
//
// -| --------- The base of the fused MetaOrdDTA functions on Stan math only ------------------------------------------
//
// The fused functions of the MetaOrdDTA models (the likelihood, the induced Dirichlet of the cutpoints, the Se / Sp
// priors at chosen cutoffs, the Dirichlet cutpoint medians and the partially centered cutpoints, each in one
// autodiff node with its hand-written gradient) use a few elementwise helpers on buffers of doubles: k_exp(p, n),
// k_log(p, n), ... in place, k_log_into(x, n, out), k_lgamma_and_digamma(x, n, lgamma, digamma), one scratch buffer
// per thread, and some constants. This file gives each of them with the standard library and Stan math only: no
// SIMD instructions, no -m flags, any CPU (x86-64 without AVX, ARM, ...). The fused functions themselves are the
// scalar code of the fused versions, unchanged apart from these helpers and their elementary calls (log, exp,
// Phi, ...), which are std:: or stan::math:: calls on doubles.
//
// Everything is in the namespace fused_stan_math; the Stan-facing functions keep the names the Stan files declare
// (stanc --allow-undefined) and are put in the global namespace with a using-declaration at the end of each file.
//
// ------------------------------------------------------------------------------------------------------------------

#include <stan/math/rev.hpp>
#include <stan/math/prim.hpp>

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <ostream>
#include <vector>

#ifndef ALWAYS_INLINE
#define ALWAYS_INLINE inline
#endif
#ifndef VECTORCALL
#define VECTORCALL
#endif
#ifndef ALIGN64
#define ALIGN64 alignas(64)
#endif
#ifndef ALIGN32
#define ALIGN32 alignas(32)
#endif

//// one value at a time (the fused functions loop over their values in blocks of BMVP_LANES)
#ifndef BMVP_LANES
#define BMVP_LANES 1
#endif

namespace fused_stan_math {

using stan::math::var;
using stan::math::vari;
using VecVar = Eigen::Matrix<var, -1, 1>;

// ---- constants
static constexpr double INV_SQRT_2PI = 0.398942280401432677939946059934;
static constexpr double SQRT_2PI     = 2.506628274631000502415765284811;
static constexpr double LN2          = 0.693147180559945309417232121458;

// ---- the normal tail of the probit link: below -2.5 the log of the derivative of Phi^-1 is taken from the
//      continued fraction of the Mills ratio (40 terms), which keeps its precision far in the tail
namespace reference {

static constexpr double NORMAL_TAIL_THRESH = 2.5;
static constexpr int    MILLS_CF_DEPTH     = 40;

}  // namespace reference

// ---- one value: the buffer helpers of the fused functions read and write doubles
using simd_t = double;
inline simd_t simd_loadu(const double* p) { return *p; }
inline void   simd_storeu(double* p, simd_t v) { *p = v; }
inline simd_t simd_load(const double* p) { return *p; }
inline void   simd_store(double* p, simd_t v) { *p = v; }

template <typename K>
inline void apply_kernel(double* p, const int n, K k) {

    for (int i = 0; i < n; ++i) p[i] = k(p[i]);

}

// ---- the sum of n doubles in index order
inline double sum_in_a_fixed_order(const double* p, const Eigen::Index n) {

    double sum = 0.0;
    for (Eigen::Index i = 0; i < n; ++i) sum += p[i];
    return sum;

}

inline void copy_n(const double* x, int n, double* out) { for (int i = 0; i < n; ++i) out[i] = x[i]; }




// ---- the elementwise functions on a buffer of n doubles, in place, with Stan math on doubles
inline void k_Phi(double* p, int n)             { for (int i = 0; i < n; ++i) p[i] = stan::math::Phi(p[i]); }
inline void k_log_Phi(double* p, int n)         { for (int i = 0; i < n; ++i) p[i] = stan::math::std_normal_lcdf(p[i]); }
inline void k_inv_Phi(double* p, int n)         { for (int i = 0; i < n; ++i) p[i] = stan::math::inv_Phi(p[i]); }
inline void k_exp(double* p, int n)             { for (int i = 0; i < n; ++i) p[i] = std::exp(p[i]); }
inline void k_log(double* p, int n)             { for (int i = 0; i < n; ++i) p[i] = std::log(p[i]); }
inline void k_log1m(double* p, int n)           { for (int i = 0; i < n; ++i) p[i] = stan::math::log1m(p[i]); }
inline void k_log1p(double* p, int n)           { for (int i = 0; i < n; ++i) p[i] = std::log1p(p[i]); }
inline void k_log1p_exp(double* p, int n)       { for (int i = 0; i < n; ++i) p[i] = stan::math::log1p_exp(p[i]); }
inline void k_log1m_exp(double* p, int n)       { for (int i = 0; i < n; ++i) p[i] = stan::math::log1m_exp(p[i]); }
inline void k_inv_logit(double* p, int n)       { for (int i = 0; i < n; ++i) p[i] = stan::math::inv_logit(p[i]); }
inline void k_log_inv_logit(double* p, int n)   { for (int i = 0; i < n; ++i) p[i] = stan::math::log_inv_logit(p[i]); }
inline void k_log1m_inv_logit(double* p, int n) { for (int i = 0; i < n; ++i) p[i] = stan::math::log1m_inv_logit(p[i]); }
inline void k_logit(double* p, int n)           { for (int i = 0; i < n; ++i) p[i] = stan::math::logit(p[i]); }
inline void k_expm1(double* p, int n)           { for (int i = 0; i < n; ++i) p[i] = std::expm1(p[i]); }
inline void k_tanh(double* p, int n)            { for (int i = 0; i < n; ++i) p[i] = std::tanh(p[i]); }
inline void k_square(double* p, int n)          { for (int i = 0; i < n; ++i) p[i] = p[i] * p[i]; }
inline void k_sqrt(double* p, int n)            { for (int i = 0; i < n; ++i) p[i] = std::sqrt(p[i]); }




// ---- lgamma and digamma of n values (Stan math's lgamma and digamma on doubles); with_digamma = false: lgamma only
template <bool with_digamma>
inline void k_lgamma_and_digamma(const double* x,
                                 const int n,
                                 double* lgamma_of_x,
                                 double* digamma_of_x) {

    for (int i = 0; i < n; ++i) {
        const double x_i = x[i];
        lgamma_of_x[i] = stan::math::lgamma(x_i);
        if constexpr (with_digamma) digamma_of_x[i] = stan::math::digamma(x_i);
    }

}

inline void k_lgamma(double* p,
                     const int n) {

    k_lgamma_and_digamma<false>(p, n, p, nullptr);

}

// ---- log(x) and log1m(y) = log1p(-y) of n values into out
inline void k_log_into(const double* x,
                       const int n,
                       double* out) {

    for (int i = 0; i < n; ++i) out[i] = std::log(x[i]);

}

inline void k_log1m_into(const double* y,
                         const int n,
                         double* out) {

    for (int i = 0; i < n; ++i) out[i] = std::log1p(-y[i]);

}

inline void k_log_and_log1m_into(const double* y,
                                 const int n,
                                 double* log_y,
                                 double* log1m_y) {

    for (int i = 0; i < n; ++i) {
        log_y[i] = std::log(y[i]);
        log1m_y[i] = std::log1p(-y[i]);
    }

}

// ---- the logs of x_1 .. x_n (log_x) and of 1 - y_1 .. 1 - y_m (log1m_y); either n or m may be 0
inline void k_log_and_log1m_short(const double* x,
                                  const int n,
                                  const double* y,
                                  const int m,
                                  double* log_x,
                                  double* log1m_y) {

    for (int i = 0; i < n; ++i) log_x[i] = std::log(x[i]);
    for (int j = 0; j < m; ++j) log1m_y[j] = std::log1p(-y[j]);

}




// ---- one scratch buffer of doubles per thread for the values of a call (no allocation once it is large enough);
//      a function takes it once and calls no other function that takes it
inline double* nico_scratch_of_doubles(const std::size_t n) {

    static thread_local std::vector<double> scratch;
    if (scratch.size() < n) scratch.resize(((n + 7) / 8) * 8);
    return scratch.data();

}

}  // namespace fused_stan_math
