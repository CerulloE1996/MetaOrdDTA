#pragma once
// ---- Stan math only: a copy of NicoSIMD/stan_external_binomial_counts_above_successive_cutpoints_of_a_cumulative_link.hpp (md5 81cc2a4a5180a7838270f674fa3a9cc7) with its scalar code and
//      its hand-written gradients, the SIMD branches left out, and log, exp, Phi, ... of doubles from
//      std:: and stan::math:: (the buffer helpers k_exp, k_log, ... are in
//      stan_external_functions_on_Stan_math_only.hpp). Where the notes below speak of SIMD or NicoSIMD
//      kernels, the same values come from Stan math here: no -m flags, any CPU.
//
// -| --------- NicoSIMD binomial counts above successive cutpoints of a cumulative link, a Stan external fn ----
//
//     nico_binomial_counts_above_cuts_cumul_link_lpmf(
//           n | n_previous, latent, is_first, link, bound)
//       = sum_k binomial_lpmf(n_k | n_previous_k, p_k),
//
//     S_k  = F(latent_k)       (F = Phi when link = 1, inv_logit otherwise: the probability above cutpoint k),
//     S0_k = 1 when is_first_k = 1 (the first reported cutpoint of its group), S_(k-1) otherwise,
//     p_k  = fmax(bound, fmin(1 - bound, S_k / S0_k))      (Stan's fmax and fmin, NaN included),
//
// with n and n_previous arrays of int (data: the count above cutpoint k, and the count above the group's previous
// reported cutpoint or the group's total), latent a vector (the location of the group minus cutpoint k) and
// is_first a vector of 0 and 1 (data). The probabilities above successive cutpoints of one group and its counts
// above them: the conditional binomial (continuation) form of the multinomial over its categories. A general
// function: any cumulative-link (probit or logit) ordinal model of counts reported at some or all of the
// cutpoints of each group (ordinal meta-analysis of test accuracy, interval-censored grouped data, ordinal
// regression with aggregated counts per category). The same values and derivatives as the Stan statements
//
//     vector[K] S  = (link == 1) ? Phi(latent) : inv_logit(latent);
//     vector[K] S0 = append_row(1.0, S[1:(K - 1)]) .* (1 - is_first) + is_first;
//     target += binomial_lpmf(n | n_previous, fmax(bound, fmin(1 - bound, S ./ S0)));
//
// with every constant (log choose(n_previous_k, n_k)) unless propto, as binomial_lpmf, in ONE autodiff node in
// place of about eight: the forward pass and the reverse pass over all K elements, with BayesMVP's SIMD kernels
// for Phi, inv_logit, exp and log, the double-precision log1m (the same kernels as nico_Phi, nico_inv_logit and
// nico_binomial_lpmf, so the same values), and the element-wise arithmetic in branch-free loops over the elements
// (each element reads only its own values and its predecessor's S, so there is no recurrence and the compiler
// vectorises them).
// The counts are checked as in nico_binomial_lpmf: Stan's checks run only when the pass over the values finds a
// bad one.
//
// Declared in the Stan functions block (no body; stanc --allow-undefined) as
//
//     real nico_binomial_counts_above_cuts_cumul_link_lpmf(
//           array[] int n, array[] int n_previous, vector latent, data vector is_first, data int link,
//           data real bound);
//
// It uses the namespace fused_stan_math, the kernels of stan_external_functions.hpp and the table of log k! of
// stan_external_binomial_lpmf.hpp, which it includes (#pragma once). Its helpers are in their own namespace,
// fused_stan_math::binomial_counts_above_successive_cutpoints.
//
// -----------------------------------------------------------------------------------------------------------

#include "stan_external_functions_on_Stan_math_only.hpp"
#include "stan_external_binomial_lpmf.hpp"

#include <stan/math/rev.hpp>
#include <stan/math/prim.hpp>

#include <Eigen/Dense>
#include <cmath>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace fused_stan_math {
namespace binomial_counts_above_successive_cutpoints {

static constexpr const char* function_name =
    "nico_binomial_counts_above_cuts_cumul_link_lpmf";




// ---- the forward pass over K elements, into buffers of K doubles each:
//        S: F(latent); p: the bounded conditional probabilities; log_p, log1m_p: their logs;
//        p_passes: 1 where p_k is S_k / S0_k (the bounds of fmax and fmin not taken, not NaN), 0 otherwise.
//      Returns false, and nothing else, when a count is out of range (the caller then runs Stan's checks).
template <bool include_log_choose>
inline bool forward_pass(const int* n, const int* n_previous, const double* latent, const double* is_first,
                         const int link, const double bound, const int K, double* S, double* p,
                         double* log_p, double* log1m_p, double* p_passes, double& sum) {

    int n_bad = 0;
    for (int k = 0; k < K; ++k) {
        n_bad += (n[k] < 0) | (n[k] > n_previous[k]);
    }
    if (n_bad > 0) return false;
    copy_n(latent, K, S);
    if (link == 1) {
        k_Phi(S, K);
    } else {
        k_inv_logit(S, K);
    }
    const double upper = 1.0 - bound;
    for (int k = 0; k < K; ++k) {
        //// (S0: 1 at the first reported cutpoint of its group, the previous S otherwise, as
        ////  append_row(1.0, S[1:(K - 1)]) .* (1 - is_first) + is_first)
        const double S_previous = (k == 0) ? 1.0 : S[k - 1];
        const double S0 = S_previous * (1.0 - is_first[k]) + is_first[k];
        const double ratio = S[k] / S0;
        //// (Stan's fmin(1 - bound, ratio): ratio when ratio < 1 - bound, else the bound, NaN giving the bound;
        ////  then fmax(bound, .): the bound when bound > value, else the value)
        const double after_fmin = (ratio < upper) ? ratio : upper;
        const double after_fmax = (bound > after_fmin) ? bound : after_fmin;
        p[k] = after_fmax;
        p_passes[k] = ((ratio < upper) && !(bound > ratio)) ? 1.0 : 0.0;
        log_p[k] = after_fmax;
        log1m_p[k] = after_fmax;
    }
    k_log(log_p, K);
    binomial::log1m_in_double_precision(log1m_p, K);
    double total = 0.0;
    for (int k = 0; k < K; ++k) {
        const double successes = n[k];
        const double failures = n_previous[k] - n[k];
        total += ((successes > 0.0) ? successes * log_p[k] : 0.0)
                 + ((failures > 0.0) ? failures * log1m_p[k] : 0.0);
    }
    if constexpr (include_log_choose) {
        const double* table = binomial::table_of_log_factorials();
        for (int k = 0; k < K; ++k) {
            total += binomial::log_factorial(n_previous[k], table) - binomial::log_factorial(n[k], table)
                     - binomial::log_factorial(n_previous[k] - n[k], table);
        }
    }
    sum = total;
    return true;

}




// ---- the reverse pass: d / d latent_k, from the forward pass's buffers (adjoint: the node's adjoint). The
//      chain of the Stan statements above: d lp / d p_k = n_k / p_k - (n_previous_k - n_k) / (1 - p_k) where p_k
//      passes the bounds (0 otherwise); p_k = S_k / S0_k gives S_k the share 1 / S0_k and, where S0_k is S_(k-1),
//      S_(k-1) the share - S_k / S0_k^2 (as Stan's division: - adjoint * (S_k / S0_k) / S0_k); then
//      d S_k / d latent_k = phi(latent_k) (probit) or S_k (1 - S_k) (logit).
inline void reverse_pass(const int* n, const int* n_previous, const double* latent, const double* is_first,
                         const int link, const int K, const double* S, const double* p, const double* p_passes,
                         const double adjoint, double* d_S, double* d_latent) {

    for (int k = 0; k < K; ++k) d_S[k] = 0.0;
    for (int k = 0; k < K; ++k) {
        const double successes = n[k];
        const double failures = n_previous[k] - n[k];
        const double d_p = adjoint * p_passes[k]
                           * (((successes > 0.0) ? successes / p[k] : 0.0)
                              - ((failures > 0.0) ? failures / (1.0 - p[k]) : 0.0));
        const double S_previous = (k == 0) ? 1.0 : S[k - 1];
        const double S0 = S_previous * (1.0 - is_first[k]) + is_first[k];
        d_S[k] += d_p / S0;
        if (k > 0) d_S[k - 1] -= d_p * (S[k] / S0) / S0 * (1.0 - is_first[k]);
    }
    if (link == 1) {
        //// (phi(latent) = exp(-latent^2 / 2) / sqrt(2 pi), with the exp kernel, as nico_Phi's derivative)
        for (int k = 0; k < K; ++k) d_latent[k] = -0.5 * latent[k] * latent[k];
        k_exp(d_latent, K);
        for (int k = 0; k < K; ++k) d_latent[k] *= INV_SQRT_2PI * d_S[k];
    } else {
        for (int k = 0; k < K; ++k) d_latent[k] = d_S[k] * S[k] * (1.0 - S[k]);
    }

}




// ---- Stan's checks of the counts (they throw; reached only when the pass over the values found a bad one)
template <typename T_n, typename T_n_previous>
inline void run_the_checks_of_Stan(const T_n& n, const T_n_previous& n_previous) {

    stan::math::check_nonnegative(function_name, "Population size parameter", n_previous);
    stan::math::check_bounded(function_name, "Successes variable", n, 0, n_previous);
    throw std::domain_error(std::string(function_name) + ": an argument is out of its range");

}

}  // namespace binomial_counts_above_successive_cutpoints
}  // namespace fused_stan_math




// ---- the Stan-facing function (global scope; stanc calls it as
//      nico_binomial_counts_above_cuts_cumul_link_lpmf<propto__>(n, n_previous,
//      latent, is_first, link, bound, pstream__))
namespace fused_stan_math {

template <bool propto, typename T_n, typename T_n_previous, typename T_latent, typename T_is_first,
          stan::require_all_std_vector_vt<std::is_integral, T_n, T_n_previous>* = nullptr,
          stan::require_eigen_col_vector_t<T_latent>* = nullptr,
          stan::require_eigen_col_vector_vt<std::is_arithmetic, T_is_first>* = nullptr>
inline stan::return_type_t<T_latent> nico_binomial_counts_above_cuts_cumul_link_lpmf(
        const T_n& n, const T_n_previous& n_previous, const T_latent& latent, const T_is_first& is_first,
        const int link, const double bound, std::ostream* /*pstream__*/) {

    namespace counts = fused_stan_math::binomial_counts_above_successive_cutpoints;
    using stan::math::var;
    stan::math::check_consistent_sizes(counts::function_name, "Successes variable", n,
                                       "Population size parameter", n_previous, "Latent values", latent,
                                       "Indicators of the first cutpoints", is_first);
    const int K = static_cast<int>(n.size());
    if (K == 0) return 0.0;
    const auto& latent_ref = stan::math::to_ref(latent);
    const Eigen::Matrix<double, -1, 1> latent_values = stan::math::value_of(latent_ref);
    const Eigen::Matrix<double, -1, 1> is_first_values = stan::math::value_of(is_first);
    if constexpr (stan::is_constant_all<T_latent>::value) {
        std::vector<double> buffers(5 * static_cast<std::size_t>(K));
        double sum = 0.0;
        if (!counts::forward_pass<true>(n.data(), n_previous.data(), latent_values.data(), is_first_values.data(),
                                        link, bound, K, buffers.data(), buffers.data() + K,
                                        buffers.data() + 2 * K, buffers.data() + 3 * K, buffers.data() + 4 * K,
                                        sum)) {

            counts::run_the_checks_of_Stan(n, n_previous);

        }
        return propto ? 0.0 : sum;
    } else {
        auto& arena = stan::math::ChainableStack::instance_->memalloc_;
        //// (on the arena, kept for the reverse pass: the data, latent's values, S, p and where p passes)
        int* n_kept = arena.alloc_array<int>(K);
        int* n_previous_kept = arena.alloc_array<int>(K);
        double* latent_kept = arena.alloc_array<double>(K);
        double* is_first_kept = arena.alloc_array<double>(K);
        double* S = arena.alloc_array<double>(K);
        double* p = arena.alloc_array<double>(K);
        double* p_passes = arena.alloc_array<double>(K);
        double* log_p = arena.alloc_array<double>(K);
        double* log1m_p = arena.alloc_array<double>(K);
        for (int k = 0; k < K; ++k) {
            n_kept[k] = n[k];
            n_previous_kept[k] = n_previous[k];
            latent_kept[k] = latent_values(k);
            is_first_kept[k] = is_first_values(k);
        }
        double sum = 0.0;
        if (!counts::forward_pass<!propto>(n_kept, n_previous_kept, latent_kept, is_first_kept, link, bound, K,
                                           S, p, log_p, log1m_p, p_passes, sum)) {

            counts::run_the_checks_of_Stan(n, n_previous);

        }
        stan::arena_t<Eigen::Matrix<var, Eigen::Dynamic, 1>> latent_vars = latent_ref;
        //// (log_p and log1m_p are not needed after the forward pass: their arena space holds d_S and d_latent)
        return stan::math::make_callback_var(sum, [latent_vars, n_kept, n_previous_kept, latent_kept,
                                                   is_first_kept, link, K, S, p, p_passes, log_p,
                                                   log1m_p](auto& vi) mutable {
            double* d_S = log_p;
            double* d_latent = log1m_p;
            counts::reverse_pass(n_kept, n_previous_kept, latent_kept, is_first_kept, link, K, S, p, p_passes,
                                 vi.adj(), d_S, d_latent);
            for (int k = 0; k < K; ++k) latent_vars.coeffRef(k).adj() += d_latent[k];
        });
    }

}

}  // namespace fused_stan_math

using fused_stan_math::nico_binomial_counts_above_cuts_cumul_link_lpmf;

// ---- Indexed overloads: fuse gathering of (x_2, n, Xbeta - C, ...) directly without Stan tape slicing/n_terms ----
namespace fused_stan_math {

// 1) Two-class random-cutpoint ordinal bivariate MA: latent[c, s, cut_i] = Xbeta(s, c) - C[c](position_of_each_C[s, cut_i] - 1)
template <bool propto, typename T_loc, typename T_C,
          stan::require_eigen_t<T_loc>* = nullptr>
inline stan::return_type_t<T_loc, T_C>
nico_binomial_counts_above_cuts_cumul_link_lpmf(
        const std::vector<std::vector<std::vector<int>>>& x_2,
        const std::vector<std::vector<std::vector<int>>>& n_prev,
        const T_loc& Xbeta,
        const std::vector<Eigen::Matrix<T_C, -1, 1>>& C,
        const std::vector<std::vector<int>>& position_of_each_C,
        const std::vector<int>& n_obs_cutpoints,
        const std::vector<int>& K_fold_CV_indicator,
        const int link,
        const double bound,
        std::ostream* /*pstream__*/) {

    namespace counts = fused_stan_math::binomial_counts_above_successive_cutpoints;
    using stan::math::var;
    const int n_studies = static_cast<int>(n_obs_cutpoints.size());
    if (static_cast<int>(C.size()) != 2) {
        throw std::domain_error(std::string(counts::function_name) + ": C must have 2 classes");
    }

    int K = 0;
    for (int s = 0; s < n_studies; ++s) {
        if (K_fold_CV_indicator[s] == 1) K += 2 * n_obs_cutpoints[s];
    }
    if (K == 0) return 0.0;

    constexpr bool is_loc_var = !stan::is_constant_all<T_loc>::value;
    constexpr bool is_C_var   = !stan::is_constant_all<T_C>::value;
    constexpr bool is_any_var = is_loc_var || is_C_var;

    const auto& Xbeta_ref = stan::math::to_ref(Xbeta);

    if constexpr (!is_any_var) {
        std::vector<int> n_kept(K), n_prev_kept(K);
        std::vector<double> latent_val(K), is_first_val(K);
        std::vector<double> buffers(5 * static_cast<std::size_t>(K));
        int t = 0;
        for (int c = 0; c < 2; ++c) {
            for (int s = 0; s < n_studies; ++s) {
                if (K_fold_CV_indicator[s] != 1) continue;
                const double loc_sc = stan::math::value_of(Xbeta_ref.coeff(s, c));
                for (int cut_i = 0; cut_i < n_obs_cutpoints[s]; ++cut_i) {
                    const int p = position_of_each_C[s][cut_i] - 1;
                    if (p < 0 || p >= C[c].size()) {
                        throw std::domain_error(std::string(counts::function_name) + ": invalid position_of_each_C");
                    }
                    n_kept[t]       = x_2[c][s][cut_i];
                    n_prev_kept[t]  = n_prev[c][s][cut_i];
                    is_first_val[t] = (cut_i == 0) ? 1.0 : 0.0;
                    latent_val[t]   = loc_sc - stan::math::value_of(C[c].coeff(p));
                    ++t;
                }
            }
        }
        double sum = 0.0;
        if (!counts::forward_pass<true>(n_kept.data(), n_prev_kept.data(), latent_val.data(), is_first_val.data(),
                                        link, bound, K, buffers.data(), buffers.data() + K,
                                        buffers.data() + 2 * K, buffers.data() + 3 * K, buffers.data() + 4 * K,
                                        sum)) {
            counts::run_the_checks_of_Stan(n_kept, n_prev_kept);
        }
        return propto ? 0.0 : sum;
    } else {
        auto& arena = stan::math::ChainableStack::instance_->memalloc_;
        int* i_buf           = arena.alloc_array<int>(5 * K);
        double* d_buf        = arena.alloc_array<double>(7 * K);
        int* n_kept          = i_buf;
        int* n_prev_kept     = i_buf + K;
        int* loc_s           = i_buf + 2 * K;
        int* loc_c           = i_buf + 3 * K;
        int* C_pos           = i_buf + 4 * K;
        double* latent_kept  = d_buf;
        double* is_first     = d_buf + K;
        double* S            = d_buf + 2 * K;
        double* p_buf        = d_buf + 3 * K;
        double* p_passes     = d_buf + 4 * K;
        double* log_p        = d_buf + 5 * K;
        double* log1m_p      = d_buf + 6 * K;

        int t = 0;
        for (int c = 0; c < 2; ++c) {
            for (int s = 0; s < n_studies; ++s) {
                if (K_fold_CV_indicator[s] != 1) continue;
                const double loc_sc = stan::math::value_of(Xbeta_ref.coeff(s, c));
                for (int cut_i = 0; cut_i < n_obs_cutpoints[s]; ++cut_i) {
                    const int p = position_of_each_C[s][cut_i] - 1;
                    if (p < 0 || p >= C[c].size()) {
                        throw std::domain_error(std::string(counts::function_name) + ": invalid position_of_each_C");
                    }
                    n_kept[t]      = x_2[c][s][cut_i];
                    n_prev_kept[t] = n_prev[c][s][cut_i];
                    loc_s[t]       = s;
                    loc_c[t]       = c;
                    C_pos[t]       = p;
                    is_first[t]    = (cut_i == 0) ? 1.0 : 0.0;
                    latent_kept[t] = loc_sc - stan::math::value_of(C[c].coeff(p));
                    ++t;
                }
            }
        }
        double sum = 0.0;
        if (!counts::forward_pass<!propto>(n_kept, n_prev_kept, latent_kept, is_first, link, bound, K,
                                           S, p_buf, log_p, log1m_p, p_passes, sum)) {
            counts::run_the_checks_of_Stan(
                std::vector<int>(n_kept, n_kept + K),
                std::vector<int>(n_prev_kept, n_prev_kept + K));
        }

        stan::arena_t<Eigen::Matrix<var, -1, -1>> Xbeta_vars;
        stan::arena_t<Eigen::Matrix<var, -1, 1>> C_vars_0;
        stan::arena_t<Eigen::Matrix<var, -1, 1>> C_vars_1;
        if constexpr (is_loc_var) Xbeta_vars = Xbeta_ref;
        if constexpr (is_C_var) {
            C_vars_0 = C[0];
            C_vars_1 = C[1];
        }
        return stan::math::make_callback_var(
            sum,
            [Xbeta_vars, C_vars_0, C_vars_1, n_kept, n_prev_kept, loc_s, loc_c, C_pos,
             latent_kept, is_first, link, K, S, p_buf, p_passes, log_p, log1m_p](auto& vi) mutable {
                double* d_S      = log_p;
                double* d_latent = log1m_p;
                counts::reverse_pass(n_kept, n_prev_kept, latent_kept, is_first, link, K,
                                     S, p_buf, p_passes, vi.adj(), d_S, d_latent);
                for (int t = 0; t < K; ++t) {
                    const double dl = d_latent[t];
                    const int c = loc_c[t];
                    if constexpr (is_loc_var) {
                        Xbeta_vars.coeffRef(loc_s[t], c).adj() += dl;
                    }
                    if constexpr (is_C_var) {
                        if (c == 0) C_vars_0.coeffRef(C_pos[t]).adj() -= dl;
                        else        C_vars_1.coeffRef(C_pos[t]).adj() -= dl;
                    }
                }
            });
    }
}

// 2) Two-class fixed-cutpoint ordinal bivariate MA: uses 3D cutpoint_index[c][s][cut_i + 1]
template <bool propto, typename T_loc, typename T_C,
          stan::require_eigen_t<T_loc>* = nullptr>
inline stan::return_type_t<T_loc, T_C>
nico_binomial_counts_above_cuts_cumul_link_lpmf(
        const std::vector<std::vector<std::vector<int>>>& x_2,
        const std::vector<std::vector<std::vector<int>>>& n_prev,
        const T_loc& Xbeta,
        const std::vector<Eigen::Matrix<T_C, -1, 1>>& C,
        const std::vector<std::vector<std::vector<int>>>& cutpoint_index,
        const std::vector<int>& n_obs_cutpoints,
        const std::vector<int>& K_fold_CV_indicator,
        const int link,
        const double bound,
        std::ostream* /*pstream__*/) {

    namespace counts = fused_stan_math::binomial_counts_above_successive_cutpoints;
    using stan::math::var;
    const int n_studies = static_cast<int>(n_obs_cutpoints.size());
    if (static_cast<int>(C.size()) != 2) {
        throw std::domain_error(std::string(counts::function_name) + ": C must have 2 classes");
    }

    int K = 0;
    for (int s = 0; s < n_studies; ++s) {
        if (K_fold_CV_indicator[s] == 1) K += 2 * n_obs_cutpoints[s];
    }
    if (K == 0) return 0.0;

    constexpr bool is_loc_var = !stan::is_constant_all<T_loc>::value;
    constexpr bool is_C_var   = !stan::is_constant_all<T_C>::value;
    constexpr bool is_any_var = is_loc_var || is_C_var;

    const auto& Xbeta_ref = stan::math::to_ref(Xbeta);

    if constexpr (!is_any_var) {
        std::vector<int> n_kept(K), n_prev_kept(K);
        std::vector<double> latent_val(K), is_first_val(K);
        std::vector<double> buffers(5 * static_cast<std::size_t>(K));
        int t = 0;
        for (int c = 0; c < 2; ++c) {
            for (int s = 0; s < n_studies; ++s) {
                if (K_fold_CV_indicator[s] != 1) continue;
                const double loc_sc = stan::math::value_of(Xbeta_ref.coeff(s, c));
                for (int cut_i = 0; cut_i < n_obs_cutpoints[s]; ++cut_i) {
                    const int k = cutpoint_index[c][s][cut_i + 1] - 1;
                    if (k < 0 || k >= C[c].size()) {
                        throw std::domain_error(std::string(counts::function_name) + ": invalid cutpoint_index");
                    }
                    n_kept[t]       = x_2[c][s][cut_i];
                    n_prev_kept[t]  = n_prev[c][s][cut_i];
                    is_first_val[t] = (cut_i == 0) ? 1.0 : 0.0;
                    latent_val[t]   = loc_sc - stan::math::value_of(C[c].coeff(k));
                    ++t;
                }
            }
        }
        double sum = 0.0;
        if (!counts::forward_pass<true>(n_kept.data(), n_prev_kept.data(), latent_val.data(), is_first_val.data(),
                                        link, bound, K, buffers.data(), buffers.data() + K,
                                        buffers.data() + 2 * K, buffers.data() + 3 * K, buffers.data() + 4 * K,
                                        sum)) {
            counts::run_the_checks_of_Stan(n_kept, n_prev_kept);
        }
        return propto ? 0.0 : sum;
    } else {
        auto& arena = stan::math::ChainableStack::instance_->memalloc_;
        int* i_buf           = arena.alloc_array<int>(5 * K);
        double* d_buf        = arena.alloc_array<double>(7 * K);
        int* n_kept          = i_buf;
        int* n_prev_kept     = i_buf + K;
        int* loc_s           = i_buf + 2 * K;
        int* loc_c           = i_buf + 3 * K;
        int* C_pos           = i_buf + 4 * K;
        double* latent_kept  = d_buf;
        double* is_first     = d_buf + K;
        double* S            = d_buf + 2 * K;
        double* p_buf        = d_buf + 3 * K;
        double* p_passes     = d_buf + 4 * K;
        double* log_p        = d_buf + 5 * K;
        double* log1m_p      = d_buf + 6 * K;

        int t = 0;
        for (int c = 0; c < 2; ++c) {
            for (int s = 0; s < n_studies; ++s) {
                if (K_fold_CV_indicator[s] != 1) continue;
                const double loc_sc = stan::math::value_of(Xbeta_ref.coeff(s, c));
                for (int cut_i = 0; cut_i < n_obs_cutpoints[s]; ++cut_i) {
                    const int k = cutpoint_index[c][s][cut_i + 1] - 1;
                    if (k < 0 || k >= C[c].size()) {
                        throw std::domain_error(std::string(counts::function_name) + ": invalid cutpoint_index");
                    }
                    n_kept[t]      = x_2[c][s][cut_i];
                    n_prev_kept[t] = n_prev[c][s][cut_i];
                    loc_s[t]       = s;
                    loc_c[t]       = c;
                    C_pos[t]       = k;
                    is_first[t]    = (cut_i == 0) ? 1.0 : 0.0;
                    latent_kept[t] = loc_sc - stan::math::value_of(C[c].coeff(k));
                    ++t;
                }
            }
        }
        double sum = 0.0;
        if (!counts::forward_pass<!propto>(n_kept, n_prev_kept, latent_kept, is_first, link, bound, K,
                                           S, p_buf, log_p, log1m_p, p_passes, sum)) {
            counts::run_the_checks_of_Stan(
                std::vector<int>(n_kept, n_kept + K),
                std::vector<int>(n_prev_kept, n_prev_kept + K));
        }

        stan::arena_t<Eigen::Matrix<var, -1, -1>> Xbeta_vars;
        stan::arena_t<Eigen::Matrix<var, -1, 1>> C_vars_0;
        stan::arena_t<Eigen::Matrix<var, -1, 1>> C_vars_1;
        if constexpr (is_loc_var) Xbeta_vars = Xbeta_ref;
        if constexpr (is_C_var) {
            C_vars_0 = C[0];
            C_vars_1 = C[1];
        }
        return stan::math::make_callback_var(
            sum,
            [Xbeta_vars, C_vars_0, C_vars_1, n_kept, n_prev_kept, loc_s, loc_c, C_pos,
             latent_kept, is_first, link, K, S, p_buf, p_passes, log_p, log1m_p](auto& vi) mutable {
                double* d_S      = log_p;
                double* d_latent = log1m_p;
                counts::reverse_pass(n_kept, n_prev_kept, latent_kept, is_first, link, K,
                                     S, p_buf, p_passes, vi.adj(), d_S, d_latent);
                for (int t = 0; t < K; ++t) {
                    const double dl = d_latent[t];
                    const int c = loc_c[t];
                    if constexpr (is_loc_var) {
                        Xbeta_vars.coeffRef(loc_s[t], c).adj() += dl;
                    }
                    if constexpr (is_C_var) {
                        if (c == 0) C_vars_0.coeffRef(C_pos[t]).adj() -= dl;
                        else        C_vars_1.coeffRef(C_pos[t]).adj() -= dl;
                    }
                }
            });
    }
}

// 3) Ordinal HSROC (shared C vector + Xlocations + Xscale): latent[c, s, cut_i] = -(C[p] + Xlocations[s, c]) * Xscale[s, c]
template <bool propto, typename T_loc, typename T_scale, typename T_C,
          stan::require_eigen_t<T_loc>* = nullptr,
          stan::require_eigen_t<T_scale>* = nullptr,
          stan::require_eigen_col_vector_t<T_C>* = nullptr>
inline stan::return_type_t<T_loc, T_scale, T_C>
nico_binomial_counts_above_cuts_cumul_link_lpmf(
        const std::vector<std::vector<std::vector<int>>>& x_2,
        const std::vector<std::vector<std::vector<int>>>& n_prev,
        const T_loc& Xlocations,
        const T_scale& Xscale,
        const T_C& C,
        const std::vector<std::vector<int>>& position_of_each_C,
        const std::vector<int>& n_obs_cutpoints,
        const std::vector<int>& K_fold_CV_indicator,
        const int link,
        const double bound,
        std::ostream* /*pstream__*/) {

    namespace counts = fused_stan_math::binomial_counts_above_successive_cutpoints;
    using stan::math::var;
    const int n_studies = static_cast<int>(n_obs_cutpoints.size());

    int K = 0;
    for (int s = 0; s < n_studies; ++s) {
        if (K_fold_CV_indicator[s] == 1) K += 2 * n_obs_cutpoints[s];
    }
    if (K == 0) return 0.0;

    constexpr bool is_loc_var   = !stan::is_constant_all<T_loc>::value;
    constexpr bool is_scale_var = !stan::is_constant_all<T_scale>::value;
    constexpr bool is_C_var     = !stan::is_constant_all<T_C>::value;
    constexpr bool is_any_var   = is_loc_var || is_scale_var || is_C_var;

    const auto& Xloc_ref   = stan::math::to_ref(Xlocations);
    const auto& Xscale_ref = stan::math::to_ref(Xscale);
    const auto& C_ref      = stan::math::to_ref(C);

    if constexpr (!is_any_var) {
        std::vector<int> n_kept(K), n_prev_kept(K);
        std::vector<double> latent_val(K), is_first_val(K);
        std::vector<double> buffers(5 * static_cast<std::size_t>(K));
        int t = 0;
        for (int c = 0; c < 2; ++c) {
            for (int s = 0; s < n_studies; ++s) {
                if (K_fold_CV_indicator[s] != 1) continue;
                const double loc_sc = stan::math::value_of(Xloc_ref.coeff(s, c));
                const double scl_sc = stan::math::value_of(Xscale_ref.coeff(s, c));
                for (int cut_i = 0; cut_i < n_obs_cutpoints[s]; ++cut_i) {
                    const int p = position_of_each_C[s][cut_i] - 1;
                    if (p < 0 || p >= C_ref.size()) {
                        throw std::domain_error(std::string(counts::function_name) + ": invalid position_of_each_C");
                    }
                    n_kept[t]       = x_2[c][s][cut_i];
                    n_prev_kept[t]  = n_prev[c][s][cut_i];
                    is_first_val[t] = (cut_i == 0) ? 1.0 : 0.0;
                    latent_val[t]   = -(stan::math::value_of(C_ref.coeff(p)) + loc_sc) * scl_sc;
                    ++t;
                }
            }
        }
        double sum = 0.0;
        if (!counts::forward_pass<true>(n_kept.data(), n_prev_kept.data(), latent_val.data(), is_first_val.data(),
                                        link, bound, K, buffers.data(), buffers.data() + K,
                                        buffers.data() + 2 * K, buffers.data() + 3 * K, buffers.data() + 4 * K,
                                        sum)) {
            counts::run_the_checks_of_Stan(n_kept, n_prev_kept);
        }
        return propto ? 0.0 : sum;
    } else {
        auto& arena = stan::math::ChainableStack::instance_->memalloc_;
        int* n_kept          = arena.alloc_array<int>(K);
        int* n_prev_kept     = arena.alloc_array<int>(K);
        int* loc_s           = arena.alloc_array<int>(K);
        int* loc_c           = arena.alloc_array<int>(K);
        int* C_pos           = arena.alloc_array<int>(K);
        double* latent_kept  = arena.alloc_array<double>(K);
        double* is_first     = arena.alloc_array<double>(K);
        double* sum_c_loc    = arena.alloc_array<double>(K);
        double* scale_kept   = arena.alloc_array<double>(K);
        double* S            = arena.alloc_array<double>(K);
        double* p_buf        = arena.alloc_array<double>(K);
        double* p_passes     = arena.alloc_array<double>(K);
        double* log_p        = arena.alloc_array<double>(K);
        double* log1m_p      = arena.alloc_array<double>(K);

        int t = 0;
        for (int c = 0; c < 2; ++c) {
            for (int s = 0; s < n_studies; ++s) {
                if (K_fold_CV_indicator[s] != 1) continue;
                const double loc_sc = stan::math::value_of(Xloc_ref.coeff(s, c));
                const double scl_sc = stan::math::value_of(Xscale_ref.coeff(s, c));
                for (int cut_i = 0; cut_i < n_obs_cutpoints[s]; ++cut_i) {
                    const int p = position_of_each_C[s][cut_i] - 1;
                    if (p < 0 || p >= C_ref.size()) {
                        throw std::domain_error(std::string(counts::function_name) + ": invalid position_of_each_C");
                    }
                    const double Cp_loc = stan::math::value_of(C_ref.coeff(p)) + loc_sc;
                    n_kept[t]      = x_2[c][s][cut_i];
                    n_prev_kept[t] = n_prev[c][s][cut_i];
                    loc_s[t]       = s;
                    loc_c[t]       = c;
                    C_pos[t]       = p;
                    sum_c_loc[t]   = Cp_loc;
                    scale_kept[t]  = scl_sc;
                    is_first[t]    = (cut_i == 0) ? 1.0 : 0.0;
                    latent_kept[t] = -Cp_loc * scl_sc;
                    ++t;
                }
            }
        }
        double sum = 0.0;
        if (!counts::forward_pass<!propto>(n_kept, n_prev_kept, latent_kept, is_first, link, bound, K,
                                           S, p_buf, log_p, log1m_p, p_passes, sum)) {
            counts::run_the_checks_of_Stan(
                std::vector<int>(n_kept, n_kept + K),
                std::vector<int>(n_prev_kept, n_prev_kept + K));
        }

        stan::arena_t<Eigen::Matrix<var, -1, -1>> Xloc_vars;
        stan::arena_t<Eigen::Matrix<var, -1, -1>> Xscale_vars;
        stan::arena_t<Eigen::Matrix<var, -1, 1>>  C_vars;
        if constexpr (is_loc_var)   Xloc_vars   = Xloc_ref;
        if constexpr (is_scale_var) Xscale_vars = Xscale_ref;
        if constexpr (is_C_var)     C_vars      = C_ref;
        return stan::math::make_callback_var(
            sum,
            [Xloc_vars, Xscale_vars, C_vars, n_kept, n_prev_kept, loc_s, loc_c, C_pos,
             sum_c_loc, scale_kept, latent_kept, is_first, link, K, S, p_buf, p_passes, log_p, log1m_p](auto& vi) mutable {
                double* d_S      = log_p;
                double* d_latent = log1m_p;
                counts::reverse_pass(n_kept, n_prev_kept, latent_kept, is_first, link, K,
                                     S, p_buf, p_passes, vi.adj(), d_S, d_latent);
                for (int t = 0; t < K; ++t) {
                    const double dl = d_latent[t];
                    const double d_sum = -dl * scale_kept[t];
                    if constexpr (is_loc_var) {
                        Xloc_vars.coeffRef(loc_s[t], loc_c[t]).adj() += d_sum;
                    }
                    if constexpr (is_C_var) {
                        C_vars.coeffRef(C_pos[t]).adj() += d_sum;
                    }
                    if constexpr (is_scale_var) {
                        Xscale_vars.coeffRef(loc_s[t], loc_c[t]).adj() += -dl * sum_c_loc[t];
                    }
                }
            });
    }
}

}  // namespace fused_stan_math

using fused_stan_math::nico_binomial_counts_above_cuts_cumul_link_lpmf;

























namespace fused_stan_math {

// 4) Two-class Box-Cox / fixed cutpoints + study-class scales (Jones):
//    latent[c, s, cut_i] = (Xbeta(s, c) - C[c](cutpoint_index[c][s][cut_i + 1] - 1)) / Xscale(s, c)
template <bool propto, typename T_loc, typename T_scale, typename T_C,
          stan::require_eigen_t<T_loc>* = nullptr,
          stan::require_eigen_t<T_scale>* = nullptr>
inline stan::return_type_t<T_loc, T_scale, T_C>
nico_binomial_counts_above_cuts_cumul_link_lpmf(
        const std::vector<std::vector<std::vector<int>>>& x_2,
        const std::vector<std::vector<std::vector<int>>>& n_prev,
        const T_loc& Xbeta,
        const T_scale& Xscale,
        const std::vector<Eigen::Matrix<T_C, -1, 1>>& C,
        const std::vector<std::vector<std::vector<int>>>& cutpoint_index,
        const std::vector<int>& n_obs_cutpoints,
        const std::vector<int>& K_fold_CV_indicator,
        const int link,
        const double bound,
        std::ostream* /*pstream__*/) {

    namespace counts = fused_stan_math::binomial_counts_above_successive_cutpoints;
    using stan::math::var;
    const int n_studies = static_cast<int>(n_obs_cutpoints.size());
    if (static_cast<int>(C.size()) != 2) {
        throw std::domain_error(std::string(counts::function_name) + ": C must have 2 classes");
    }

    int K = 0;
    for (int s = 0; s < n_studies; ++s) {
        if (K_fold_CV_indicator[s] == 1) K += 2 * n_obs_cutpoints[s];
    }
    if (K == 0) return 0.0;

    constexpr bool is_loc_var   = !stan::is_constant_all<T_loc>::value;
    constexpr bool is_scale_var = !stan::is_constant_all<T_scale>::value;
    constexpr bool is_C_var     = !stan::is_constant_all<T_C>::value;
    constexpr bool is_any_var   = is_loc_var || is_scale_var || is_C_var;

    const auto& Xloc_ref   = stan::math::to_ref(Xbeta);
    const auto& Xscale_ref = stan::math::to_ref(Xscale);
    const int n_cs = 2 * n_studies;

    if constexpr (!is_any_var) {
        std::vector<int> n_kept(K), n_prev_kept(K);
        std::vector<double> latent_val(K), is_first_val(K);
        std::vector<double> buffers(5 * static_cast<std::size_t>(K));
        int t = 0;
        for (int c = 0; c < 2; ++c) {
            const Eigen::Index C_sz = C[c].size();
            for (int s = 0; s < n_studies; ++s) {
                if (K_fold_CV_indicator[s] != 1) continue;
                const int m_s = n_obs_cutpoints[s];
                const double loc_sc = stan::math::value_of(Xloc_ref.coeff(s, c));
                const double inv_s  = 1.0 / stan::math::value_of(Xscale_ref.coeff(s, c));
                for (int cut_i = 0; cut_i < m_s; ++cut_i) {
                    const int k = cutpoint_index[c][s][cut_i + 1] - 1;
                    if (k < 0 || k >= C_sz) {
                        throw std::domain_error(std::string(counts::function_name) + ": invalid cutpoint_index");
                    }
                    n_kept[t]       = x_2[c][s][cut_i];
                    n_prev_kept[t]  = n_prev[c][s][cut_i];
                    is_first_val[t] = (cut_i == 0) ? 1.0 : 0.0;
                    latent_val[t]   = (loc_sc - stan::math::value_of(C[c].coeff(k))) * inv_s;
                    ++t;
                }
            }
        }
        double sum = 0.0;
        if (!counts::forward_pass<true>(n_kept.data(), n_prev_kept.data(), latent_val.data(), is_first_val.data(),
                                        link, bound, K, buffers.data(), buffers.data() + K,
                                        buffers.data() + 2 * K, buffers.data() + 3 * K, buffers.data() + 4 * K,
                                        sum)) {
            counts::run_the_checks_of_Stan(n_kept, n_prev_kept);
        }
        return propto ? 0.0 : sum;
    } else {
        auto& arena = stan::math::ChainableStack::instance_->memalloc_;
        int* i_buf           = arena.alloc_array<int>(3 * K + n_cs);
        double* d_buf        = arena.alloc_array<double>(7 * K + n_cs);
        int* n_kept          = i_buf;
        int* n_prev_kept     = i_buf + K;
        int* C_pos           = i_buf + 2 * K;
        int* cs_len          = i_buf + 3 * K;
        double* latent_kept  = d_buf;
        double* is_first     = d_buf + K;
        double* S            = d_buf + 2 * K;
        double* p_buf        = d_buf + 3 * K;
        double* p_passes     = d_buf + 4 * K;
        double* log_p        = d_buf + 5 * K;
        double* log1m_p      = d_buf + 6 * K;
        double* inv_scl      = d_buf + 7 * K;

        int t = 0;
        for (int c = 0; c < 2; ++c) {
            const Eigen::Index C_sz = C[c].size();
            for (int s = 0; s < n_studies; ++s) {
                const int cs = c * n_studies + s;
                if (K_fold_CV_indicator[s] != 1) {
                    cs_len[cs] = 0;
                    inv_scl[cs] = 0.0;
                    continue;
                }
                const int m_s = n_obs_cutpoints[s];
                cs_len[cs] = m_s;
                const double loc_sc = stan::math::value_of(Xloc_ref.coeff(s, c));
                const double inv_s  = 1.0 / stan::math::value_of(Xscale_ref.coeff(s, c));
                inv_scl[cs] = inv_s;
                for (int cut_i = 0; cut_i < m_s; ++cut_i) {
                    const int k = cutpoint_index[c][s][cut_i + 1] - 1;
                    if (k < 0 || k >= C_sz) {
                        throw std::domain_error(std::string(counts::function_name) + ": invalid cutpoint_index");
                    }
                    n_kept[t]      = x_2[c][s][cut_i];
                    n_prev_kept[t] = n_prev[c][s][cut_i];
                    C_pos[t]       = k;
                    is_first[t]    = (cut_i == 0) ? 1.0 : 0.0;
                    latent_kept[t] = (loc_sc - stan::math::value_of(C[c].coeff(k))) * inv_s;
                    ++t;
                }
            }
        }
        double sum = 0.0;
        if (!counts::forward_pass<!propto>(n_kept, n_prev_kept, latent_kept, is_first, link, bound, K,
                                           S, p_buf, log_p, log1m_p, p_passes, sum)) {
            counts::run_the_checks_of_Stan(
                std::vector<int>(n_kept, n_kept + K),
                std::vector<int>(n_prev_kept, n_prev_kept + K));
        }

        stan::arena_t<Eigen::Matrix<var, -1, -1>> Xloc_vars;
        stan::arena_t<Eigen::Matrix<var, -1, -1>> Xscale_vars;
        stan::arena_t<Eigen::Matrix<var, -1, 1>>  C_vars_0;
        stan::arena_t<Eigen::Matrix<var, -1, 1>>  C_vars_1;
        if constexpr (is_loc_var)   Xloc_vars   = Xloc_ref;
        if constexpr (is_scale_var) Xscale_vars = Xscale_ref;
        if constexpr (is_C_var) {
            C_vars_0 = C[0];
            C_vars_1 = C[1];
        }
        return stan::math::make_callback_var(
            sum,
            [Xloc_vars, Xscale_vars, C_vars_0, C_vars_1, n_studies, n_kept, n_prev_kept, C_pos, cs_len,
             inv_scl, latent_kept, is_first, link, K, S, p_buf, p_passes, log_p, log1m_p](auto& vi) mutable {
                double* d_S      = log_p;
                double* d_latent = log1m_p;
                counts::reverse_pass(n_kept, n_prev_kept, latent_kept, is_first, link, K,
                                     S, p_buf, p_passes, vi.adj(), d_S, d_latent);
                int n_thr_0 = 0;
                int n_thr_1 = 0;
                double stack_dC[256];
                std::vector<double> heap_dC;
                double* d_C0 = nullptr;
                double* d_C1 = nullptr;
                if constexpr (is_C_var) {
                    n_thr_0 = static_cast<int>(C_vars_0.size());
                    n_thr_1 = static_cast<int>(C_vars_1.size());
                    const int tot_C = n_thr_0 + n_thr_1;
                    if (tot_C <= 256) {
                        d_C0 = stack_dC;
                        std::fill_n(d_C0, tot_C, 0.0);
                    } else {
                        heap_dC.assign(tot_C, 0.0);
                        d_C0 = heap_dC.data();
                    }
                    d_C1 = d_C0 + n_thr_0;
                }
                int pos = 0;
                for (int c = 0; c < 2; ++c) {
                    double* d_Cc = (c == 0) ? d_C0 : d_C1;
                    for (int s = 0; s < n_studies; ++s) {
                        const int cs = c * n_studies + s;
                        const int m_s = cs_len[cs];
                        if (m_s == 0) continue;
                        const double inv_s = inv_scl[cs];
                        double acc_loc = 0.0;
                        double acc_scl = 0.0;
                        for (int cut_i = 0; cut_i < m_s; ++cut_i, ++pos) {
                            const double dl = d_latent[pos];
                            acc_loc += dl;
                            if constexpr (is_scale_var) {
                                acc_scl += dl * latent_kept[pos];
                            }
                            if constexpr (is_C_var) {
                                d_Cc[C_pos[pos]] -= dl * inv_s;
                            }
                        }
                        if constexpr (is_loc_var) {
                            Xloc_vars.coeffRef(s, c).adj() += acc_loc * inv_s;
                        }
                        if constexpr (is_scale_var) {
                            Xscale_vars.coeffRef(s, c).adj() -= acc_scl * inv_s;
                        }
                    }
                }
                if constexpr (is_C_var) {
                    for (int k = 0; k < n_thr_0; ++k) C_vars_0.coeffRef(k).adj() += d_C0[k];
                    for (int k = 0; k < n_thr_1; ++k) C_vars_1.coeffRef(k).adj() += d_C1[k];
                }
            });
    }
}

}  // namespace fused_stan_math

using fused_stan_math::nico_binomial_counts_above_cuts_cumul_link_lpmf;
