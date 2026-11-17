#pragma once
// ---- Stan math only: a copy of NicoSIMD/stan_external_binomial_lpmf.hpp (md5 94d2d9abd86550a7fc52e95bf5753895) with its scalar code and
//      its hand-written gradients, the SIMD branches left out, and log, exp, Phi, ... of doubles from
//      std:: and stan::math:: (the buffer helpers k_exp, k_log, ... are in
//      stan_external_functions_on_Stan_math_only.hpp). Where the notes below speak of SIMD or NicoSIMD
//      kernels, the same values come from Stan math here: no -m flags, any CPU.
//
// -| --------- NicoSIMD binomial log probability, a Stan external function -------------------------------------
//
//     fused_binomial_lpmf(n | N, theta) = sum_i [ n_i log(theta_i) + (N_i - n_i) log1m(theta_i)
//                                               + log choose(N_i, n_i) ],
//
// the log probability of Stan's binomial_lpmf(n | N, theta), with n and N arrays of int (data) and theta a vector
// of probabilities. A general function: any binomial counts with their own probabilities (e.g. the counts above
// each threshold of ordinal and continuation-ratio models, dose-response data, aggregated Bernoulli data).
// ONE autodiff node for the whole vector, with the BayesMVP log kernel and the double-precision log1m over all
// its elements (Stan's version: one node too, but std::log, log1m and three lgamma() per element, at every call),
// and:
//   - log choose(N_i, n_i) = log N_i! - log n_i! - log (N_i - n_i)!, from a table of log k! for k < 4096 (made
//     once, with std::lgamma(k + 1); std::lgamma above it), left out with propto, as in Stan;
//   - n_i log(theta_i) is 0 when n_i = 0, and (N_i - n_i) log1m(theta_i) is 0 when n_i = N_i, as in Stan (so
//     theta_i = 0 or 1 is allowed where its counts allow it);
//   - d / d theta_i = n_i / theta_i - (N_i - n_i) / (1 - theta_i), with the same two cases.
// The arguments are checked with Stan's own check functions (and so the same exceptions) only when the pass over
// the values finds a bad one (n_i outside 0, ..., N_i, or theta_i outside [0, 1] or NaN), so that valid
// arguments cost no extra pass.
//
// Declared in the Stan functions block (no body; stanc --allow-undefined) as
//
//     real fused_binomial_lpmf(array[] int n, array[] int N, vector theta);
//
// which also gives  fused_binomial_lupmf(n | N, theta)  and  n ~ fused_binomial(N, theta).
//
// It uses the namespace fused_stan_math, BMVP_LANES, the kernel k_log() and the double-precision log1p of
// stan_external_functions.hpp, which it includes (#pragma once), so that a user header may include this file
// alone. Its helpers are in their own namespace, fused_stan_math::binomial.
//
// -----------------------------------------------------------------------------------------------------------

#include "stan_external_functions_on_Stan_math_only.hpp"

#include <stan/math/rev.hpp>
#include <stan/math/prim.hpp>

#include <Eigen/Dense>
#include <cmath>
#include <ostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace fused_stan_math {
namespace binomial {

static constexpr const char* function_name = "fused_binomial_lpmf";
static constexpr int size_of_the_table_of_log_factorials = 4096;




// ---- log k! for k = 0, ..., 4095 (a function-local static: made once, by the first call, thread-safe)
inline const double* table_of_log_factorials() {

    static const std::vector<double> table = [] {
        std::vector<double> log_factorials(size_of_the_table_of_log_factorials);
        for (int k = 0; k < size_of_the_table_of_log_factorials; ++k) {
            log_factorials[k] = std::lgamma(k + 1.0);
        }
        return log_factorials;
    }();
    return table.data();

}

inline double log_factorial(const int k, const double* table) {

    return (k < size_of_the_table_of_log_factorials) ? table[k] : std::lgamma(k + 1.0);

}




// ---- log(1 - p) in place, in double precision
inline void log1m_in_double_precision(double* p,
                                      const int n) {

    for (int i = 0; i < n; ++i) {
        p[i] = std::log1p(-p[i]);
    }

}




// ---- the sum of the log probabilities, from the values of theta (M elements). log_theta and log1m_theta:
//      buffers of M doubles; with compute_partials, d / d theta_i is written to partials[i]. Returns false,
//      and nothing else, when the pass over the values finds a bad argument (the caller then runs Stan's
//      checks).
template <bool include_log_choose, bool compute_partials>
inline bool sum_of_log_probabilities(const int* n, const int* N, const double* theta, const int M,
                                     double* log_theta, double* log1m_theta, double* partials, double& sum) {

    int n_bad = 0;
    for (int i = 0; i < M; ++i) {
        const double theta_i = theta[i];
        n_bad += (n[i] < 0) | (n[i] > N[i]) | !(theta_i >= 0.0 && theta_i <= 1.0);
        log_theta[i] = theta_i;
        log1m_theta[i] = theta_i;
    }
    if (n_bad > 0) return false;
    k_log(log_theta, M);
    log1m_in_double_precision(log1m_theta, M);
    double total = 0.0;
    for (int i = 0; i < M; ++i) {
        const double successes = n[i];
        const double failures = N[i] - n[i];
        total += ((successes > 0.0) ? successes * log_theta[i] : 0.0)
                 + ((failures > 0.0) ? failures * log1m_theta[i] : 0.0);
        if constexpr (compute_partials) {
            partials[i] = ((successes > 0.0) ? successes / theta[i] : 0.0)
                          - ((failures > 0.0) ? failures / (1.0 - theta[i]) : 0.0);
        }
    }
    if constexpr (include_log_choose) {
        const double* table = table_of_log_factorials();
        for (int i = 0; i < M; ++i) {
            total += log_factorial(N[i], table) - log_factorial(n[i], table) - log_factorial(N[i] - n[i], table);
        }
    }
    sum = total;
    return true;

}




// ---- Stan's checks of the arguments (they throw; reached only when the pass over the values found a bad one)
template <typename T_n, typename T_N, typename T_theta_values>
inline void run_the_checks_of_Stan(const T_n& n, const T_N& N, const T_theta_values& theta_values) {

    stan::math::check_nonnegative(function_name, "Population size parameter", N);
    stan::math::check_bounded(function_name, "Successes variable", n, 0, N);
    stan::math::check_bounded(function_name, "Probability parameter", theta_values, 0.0, 1.0);
    throw std::domain_error(std::string(function_name) + ": an argument is out of its range");

}

}  // namespace binomial
}  // namespace fused_stan_math




// ---- the Stan-facing function (global scope; stanc calls it as  fused_binomial_lpmf<propto__>(n, N, theta,
//      pstream__), with n and N arrays of int and theta a vector or vector expression, of double or of var)
namespace fused_stan_math {

template <bool propto, typename T_n, typename T_N, typename T_theta,
          stan::require_all_std_vector_vt<std::is_integral, T_n, T_N>* = nullptr,
          stan::require_eigen_col_vector_t<T_theta>* = nullptr>
inline stan::return_type_t<T_theta> fused_binomial_lpmf(const T_n& n, const T_N& N, const T_theta& theta,
                                                       std::ostream* /*pstream__*/) {

    namespace binomial = fused_stan_math::binomial;
    using stan::math::var;
    stan::math::check_consistent_sizes(binomial::function_name, "Successes variable", n,
                                       "Population size parameter", N, "Probability parameter", theta);
    const int M = static_cast<int>(n.size());
    if (M == 0) return 0.0;
    const auto& theta_ref = stan::math::to_ref(theta);
    const Eigen::Matrix<double, -1, 1> theta_values = stan::math::value_of(theta_ref);
    if constexpr (stan::is_constant_all<T_theta>::value) {
        // (nothing depends on a parameter: with propto, Stan checks the arguments, then returns 0)
        std::vector<double> log_theta(M);
        std::vector<double> log1m_theta(M);
        double sum = 0.0;
        if (!binomial::sum_of_log_probabilities<true, false>(n.data(), N.data(), theta_values.data(), M,
                                                             log_theta.data(), log1m_theta.data(), nullptr,
                                                             sum)) {

            binomial::run_the_checks_of_Stan(n, N, theta_values);

        }
        return propto ? 0.0 : sum;
    } else {
        auto& arena = stan::math::ChainableStack::instance_->memalloc_;
        double* log_theta = arena.alloc_array<double>(M);
        double* log1m_theta = arena.alloc_array<double>(M);
        double* partials = arena.alloc_array<double>(M);
        double sum = 0.0;
        if (!binomial::sum_of_log_probabilities<!propto, true>(n.data(), N.data(), theta_values.data(), M,
                                                               log_theta, log1m_theta, partials, sum)) {

            binomial::run_the_checks_of_Stan(n, N, theta_values);

        }
        stan::arena_t<Eigen::Matrix<var, Eigen::Dynamic, 1>> theta_vars = theta_ref;
        return stan::math::make_callback_var(sum, [theta_vars, partials, M](auto& vi) mutable {
            const double adjoint = vi.adj();
            for (int i = 0; i < M; ++i) theta_vars.coeffRef(i).adj() += adjoint * partials[i];
        });
    }

}

}  // namespace fused_stan_math

using fused_stan_math::fused_binomial_lpmf;

























// ---- Stratified bivariate binomial cumulative-link likelihood (fuses linear predictor + link + binomial in ONE node)
namespace fused_stan_math {
namespace strat_bivariate_binomial {

static constexpr const char* function_name = "fused_strat_bivariate_binomial_lpmf";

template <bool include_log_choose>
inline bool forward_pass(
        const int* n,
        const int* N,
        const double* latent,
        const int link,
        const double bound,
        const int K,
        double* p,
        double* log_p,
        double* log1m_p,
        double* p_passes,
        double& sum) {

    int n_bad = 0;
    for (int t = 0; t < K; ++t) {
        n_bad += (n[t] < 0) | (n[t] > N[t]);
    }
    if (n_bad > 0) return false;
    copy_n(latent, K, p);
    if (link == 1) {
        k_Phi(p, K);
    } else {
        k_inv_logit(p, K);
    }
    const double upper = 1.0 - bound;
    for (int t = 0; t < K; ++t) {
        const double ratio = p[t];
        const double after_fmin = (ratio < upper) ? ratio : upper;
        const double after_fmax = (bound > after_fmin) ? bound : after_fmin;
        p[t] = after_fmax;
        p_passes[t] = ((ratio < upper) && !(bound > ratio)) ? 1.0 : 0.0;
        log_p[t] = after_fmax;
        log1m_p[t] = after_fmax;
    }
    k_log(log_p, K);
    binomial::log1m_in_double_precision(log1m_p, K);
    double total = 0.0;
    for (int t = 0; t < K; ++t) {
        const double successes = n[t];
        const double failures = N[t] - n[t];
        total += ((successes > 0.0) ? successes * log_p[t] : 0.0)
                 + ((failures > 0.0) ? failures * log1m_p[t] : 0.0);
    }
    if constexpr (include_log_choose) {
        const double* table = binomial::table_of_log_factorials();
        for (int t = 0; t < K; ++t) {
            total += binomial::log_factorial(N[t], table)
                     - binomial::log_factorial(n[t], table)
                     - binomial::log_factorial(N[t] - n[t], table);
        }
    }
    sum = total;
    return true;
}

}  // namespace strat_bivariate_binomial

template <bool propto,
          typename T_n, typename T_N,
          typename T_k, typename T_p, typename T_c,
          typename T_mu, typename T_sigma, typename T_L, typename T_z,
          stan::require_all_std_vector_vt<std::is_integral, T_n, T_N, T_k, T_p, T_c>* = nullptr,
          stan::require_eigen_col_vector_t<T_z>* = nullptr>
inline stan::return_type_t<T_mu, T_sigma, T_L, T_z>
fused_strat_bivariate_binomial_lpmf(
        const T_n& n,
        const T_N& N,
        const T_k& term_k,
        const T_p& term_pair,
        const T_c& term_c,
        const std::vector<Eigen::Matrix<T_mu, -1, 1>>& mu,
        const std::vector<Eigen::Matrix<T_sigma, -1, 1>>& sigma,
        const std::vector<Eigen::Matrix<T_L, -1, -1>>& L_Omega,
        const T_z& z,
        const int link,
        const double bound,
        std::ostream* /*pstream__*/) {

    namespace sb = fused_stan_math::strat_bivariate_binomial;
    using stan::math::var;

    const int K = static_cast<int>(n.size());
    if (static_cast<int>(N.size()) != K ||
        static_cast<int>(term_k.size()) != K ||
        static_cast<int>(term_pair.size()) != K ||
        static_cast<int>(term_c.size()) != K) {
        throw std::invalid_argument(std::string(sb::function_name) + ": inconsistent term array sizes");
    }
    const int n_thr = static_cast<int>(mu.size());
    if (static_cast<int>(sigma.size()) != n_thr || static_cast<int>(L_Omega.size()) != n_thr) {
        throw std::invalid_argument(std::string(sb::function_name) + ": mu, sigma, L_Omega must have length n_thr");
    }
    for (int k = 0; k < n_thr; ++k) {
        if (mu[k].size() != 2 || sigma[k].size() != 2 || L_Omega[k].rows() != 2 || L_Omega[k].cols() != 2) {
            throw std::invalid_argument(std::string(sb::function_name) + ": invalid per-threshold dimensions");
        }
    }
    const auto& z_ref = stan::math::to_ref(z);
    if ((z_ref.size() % 2) != 0) {
        throw std::invalid_argument(std::string(sb::function_name) + ": z must have even length 2 * n_pairs");
    }
    const int n_pairs = static_cast<int>(z_ref.size()) / 2;
    if (K == 0) return 0.0;

    constexpr bool is_mu_var    = !stan::is_constant_all<T_mu>::value;
    constexpr bool is_sigma_var = !stan::is_constant_all<T_sigma>::value;
    constexpr bool is_L_var     = !stan::is_constant_all<T_L>::value;
    constexpr bool is_z_var     = !stan::is_constant_all<T_z>::value;
    constexpr bool is_any_var   = is_mu_var || is_sigma_var || is_L_var || is_z_var;

    const int n_z = 2 * n_pairs;
    if constexpr (!is_any_var) {
        std::vector<double> thr_buf(6 * static_cast<std::size_t>(n_thr));
        std::vector<double> z_val(n_z);
        std::vector<double> d_buf(5 * static_cast<std::size_t>(K));
        double* mu0  = thr_buf.data();
        double* mu1  = mu0 + n_thr;
        double* sig0 = mu1 + n_thr;
        double* sig1 = sig0 + n_thr;
        double* l21  = sig1 + n_thr;
        double* l22  = l21 + n_thr;
        for (int k = 0; k < n_thr; ++k) {
            mu0[k]  = stan::math::value_of(mu[k].coeff(0));
            mu1[k]  = stan::math::value_of(mu[k].coeff(1));
            sig0[k] = stan::math::value_of(sigma[k].coeff(0));
            sig1[k] = stan::math::value_of(sigma[k].coeff(1));
            l21[k]  = stan::math::value_of(L_Omega[k].coeff(1, 0));
            l22[k]  = stan::math::value_of(L_Omega[k].coeff(1, 1));
        }
        for (int j = 0; j < n_z; ++j) z_val[j] = stan::math::value_of(z_ref.coeff(j));
        double* latent   = d_buf.data();
        double* p_buf    = latent + K;
        double* log_p    = p_buf + K;
        double* log1m_p  = log_p + K;
        double* p_passes = log1m_p + K;
        for (int t = 0; t < K; ++t) {
            const int k  = term_k[t] - 1;
            const int pr = term_pair[t] - 1;
            const int c  = term_c[t];
            if (k < 0 || k >= n_thr || pr < 0 || pr >= n_pairs || (c != 1 && c != 2)) {
                throw std::domain_error(std::string(sb::function_name) + ": invalid term index");
            }
            const int zb = 2 * pr;
            latent[t] = (c == 1)
                ? (mu0[k] + sig0[k] * z_val[zb])
                : (mu1[k] + sig1[k] * (l21[k] * z_val[zb] + l22[k] * z_val[zb + 1]));
        }
        double sum = 0.0;
        if (!sb::forward_pass<true>(n.data(), N.data(), latent, link, bound, K,
                                    p_buf, log_p, log1m_p, p_passes, sum)) {
            binomial::run_the_checks_of_Stan(n, N, Eigen::Matrix<double, -1, 1>::Constant(K, 0.5));
        }
        return propto ? 0.0 : sum;
    } else {
        auto& arena = stan::math::ChainableStack::instance_->memalloc_;
        int* i_buf      = arena.alloc_array<int>(5 * K);
        double* d_buf   = arena.alloc_array<double>(5 * K + 6 * n_thr + n_z);
        int* n_kept     = i_buf;
        int* N_kept     = n_kept + K;
        int* k_kept     = N_kept + K;
        int* zbase_kept = k_kept + K;
        int* c_kept     = zbase_kept + K;

        double* latent   = d_buf;
        double* p_buf    = latent + K;
        double* log_p    = p_buf + K;
        double* log1m_p  = log_p + K;
        double* p_passes = log1m_p + K;
        double* mu0      = p_passes + K;
        double* mu1      = mu0 + n_thr;
        double* sig0     = mu1 + n_thr;
        double* sig1     = sig0 + n_thr;
        double* l21      = sig1 + n_thr;
        double* l22      = l21 + n_thr;
        double* z_val    = l22 + n_thr;

        for (int k = 0; k < n_thr; ++k) {
            mu0[k]  = stan::math::value_of(mu[k].coeff(0));
            mu1[k]  = stan::math::value_of(mu[k].coeff(1));
            sig0[k] = stan::math::value_of(sigma[k].coeff(0));
            sig1[k] = stan::math::value_of(sigma[k].coeff(1));
            l21[k]  = stan::math::value_of(L_Omega[k].coeff(1, 0));
            l22[k]  = stan::math::value_of(L_Omega[k].coeff(1, 1));
        }
        for (int j = 0; j < n_z; ++j) z_val[j] = stan::math::value_of(z_ref.coeff(j));

        for (int t = 0; t < K; ++t) {
            const int k  = term_k[t] - 1;
            const int pr = term_pair[t] - 1;
            const int c  = term_c[t];
            if (k < 0 || k >= n_thr || pr < 0 || pr >= n_pairs || (c != 1 && c != 2)) {
                throw std::domain_error(std::string(sb::function_name) + ": invalid term index");
            }
            const int zb = 2 * pr;
            n_kept[t]     = n[t];
            N_kept[t]     = N[t];
            k_kept[t]     = k;
            zbase_kept[t] = zb;
            c_kept[t]     = c;
            latent[t] = (c == 1)
                ? (mu0[k] + sig0[k] * z_val[zb])
                : (mu1[k] + sig1[k] * (l21[k] * z_val[zb] + l22[k] * z_val[zb + 1]));
        }

        double sum = 0.0;
        if (!sb::forward_pass<!propto>(n_kept, N_kept, latent, link, bound, K,
                                       p_buf, log_p, log1m_p, p_passes, sum)) {
            binomial::run_the_checks_of_Stan(n, N, Eigen::Matrix<double, -1, 1>::Constant(K, 0.5));
        }

        stan::arena_t<Eigen::Matrix<var, -1, 1>> mu0_vars(is_mu_var ? n_thr : 0);
        stan::arena_t<Eigen::Matrix<var, -1, 1>> mu1_vars(is_mu_var ? n_thr : 0);
        stan::arena_t<Eigen::Matrix<var, -1, 1>> sig0_vars(is_sigma_var ? n_thr : 0);
        stan::arena_t<Eigen::Matrix<var, -1, 1>> sig1_vars(is_sigma_var ? n_thr : 0);
        stan::arena_t<Eigen::Matrix<var, -1, 1>> l21_vars(is_L_var ? n_thr : 0);
        stan::arena_t<Eigen::Matrix<var, -1, 1>> l22_vars(is_L_var ? n_thr : 0);
        stan::arena_t<Eigen::Matrix<var, -1, 1>> z_vars;
        if constexpr (is_mu_var) {
            for (int k = 0; k < n_thr; ++k) {
                mu0_vars.coeffRef(k) = mu[k].coeff(0);
                mu1_vars.coeffRef(k) = mu[k].coeff(1);
            }
        }
        if constexpr (is_sigma_var) {
            for (int k = 0; k < n_thr; ++k) {
                sig0_vars.coeffRef(k) = sigma[k].coeff(0);
                sig1_vars.coeffRef(k) = sigma[k].coeff(1);
            }
        }
        if constexpr (is_L_var) {
            for (int k = 0; k < n_thr; ++k) {
                l21_vars.coeffRef(k) = L_Omega[k].coeff(1, 0);
                l22_vars.coeffRef(k) = L_Omega[k].coeff(1, 1);
            }
        }
        if constexpr (is_z_var) {
            z_vars = z_ref;
        }

        return stan::math::make_callback_var(
            sum,
            [mu0_vars, mu1_vars, sig0_vars, sig1_vars, l21_vars, l22_vars, z_vars,
             n_kept, N_kept, k_kept, zbase_kept, c_kept,
             latent, p_buf, p_passes, log_p, log1m_p,
             sig0, sig1, l21, l22, z_val, link, K, n_thr](auto& vi) mutable {
                const double adjoint = vi.adj();
                double* d_S = log_p;
                double* d_latent = log1m_p;
                for (int t = 0; t < K; ++t) {
                    const double successes = n_kept[t];
                    const double failures = N_kept[t] - n_kept[t];
                    d_S[t] = adjoint * p_passes[t]
                             * (((successes > 0.0) ? successes / p_buf[t] : 0.0)
                                - ((failures > 0.0) ? failures / (1.0 - p_buf[t]) : 0.0));
                }
                if (link == 1) {
                    for (int t = 0; t < K; ++t) d_latent[t] = -0.5 * latent[t] * latent[t];
                    k_exp(d_latent, K);
                    for (int t = 0; t < K; ++t) d_latent[t] *= INV_SQRT_2PI * d_S[t];
                } else {
                    for (int t = 0; t < K; ++t) d_latent[t] = d_S[t] * p_buf[t] * (1.0 - p_buf[t]);
                }

                static constexpr int STACK_MAX_THR = 256;
                double stack_acc[6 * STACK_MAX_THR];
                std::vector<double> heap_acc;
                double* acc_base = (n_thr <= STACK_MAX_THR) ? stack_acc : (heap_acc.resize(6 * n_thr), heap_acc.data());
                std::fill_n(acc_base, 6 * n_thr, 0.0);
                double* adj_mu0  = acc_base;
                double* adj_mu1  = adj_mu0 + n_thr;
                double* adj_sig0 = adj_mu1 + n_thr;
                double* adj_sig1 = adj_sig0 + n_thr;
                double* adj_l21  = adj_sig1 + n_thr;
                double* adj_l22  = adj_l21 + n_thr;

                for (int t = 0; t < K; ++t) {
                    const double dl = d_latent[t];
                    const int k  = k_kept[t];
                    const int zb = zbase_kept[t];
                    const double z0 = z_val[zb];
                    if (c_kept[t] == 1) {
                        if constexpr (is_mu_var)    adj_mu0[k]  += dl;
                        if constexpr (is_sigma_var) adj_sig0[k] += dl * z0;
                        if constexpr (is_z_var)     z_vars.coeffRef(zb).adj() += dl * sig0[k];
                    } else {
                        const double z1 = z_val[zb + 1];
                        const double z_lin = l21[k] * z0 + l22[k] * z1;
                        if constexpr (is_mu_var)    adj_mu1[k]  += dl;
                        if constexpr (is_sigma_var) adj_sig1[k] += dl * z_lin;
                        const double dz = dl * sig1[k];
                        if constexpr (is_L_var) {
                            adj_l21[k] += dz * z0;
                            adj_l22[k] += dz * z1;
                        }
                        if constexpr (is_z_var) {
                            z_vars.coeffRef(zb).adj()     += dz * l21[k];
                            z_vars.coeffRef(zb + 1).adj() += dz * l22[k];
                        }
                    }
                }

                for (int k = 0; k < n_thr; ++k) {
                    if constexpr (is_mu_var) {
                        mu0_vars.coeffRef(k).adj() += adj_mu0[k];
                        mu1_vars.coeffRef(k).adj() += adj_mu1[k];
                    }
                    if constexpr (is_sigma_var) {
                        sig0_vars.coeffRef(k).adj() += adj_sig0[k];
                        sig1_vars.coeffRef(k).adj() += adj_sig1[k];
                    }
                    if constexpr (is_L_var) {
                        l21_vars.coeffRef(k).adj() += adj_l21[k];
                        l22_vars.coeffRef(k).adj() += adj_l22[k];
                    }
                }
            });
    }
}

}  // namespace fused_stan_math

using fused_stan_math::fused_strat_bivariate_binomial_lpmf;
