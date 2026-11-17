#pragma once
// ---- Stan math only: a copy of NicoSIMD/stan_external_induced_dirichlet_prior_of_ordered_cutpoints.hpp (md5 7036303ebeb7f0a5b49c51855904839c) with its scalar code and
//      its hand-written gradients, the SIMD branches left out, and log, exp, Phi, ... of doubles from
//      std:: and stan::math:: (the buffer helpers k_exp, k_log, ... are in
//      stan_external_functions_on_Stan_math_only.hpp). Where the notes below speak of SIMD or NicoSIMD
//      kernels, the same values come from Stan math here: no -m flags, any CPU.
//
// -| --------- NicoSIMD induced-Dirichlet prior of ordered cutpoints, a Stan external function ----------------
//
//     fused_induced_dirichlet_prior_of_ordered_cutpoints_lpdf(C | alpha, link)
//       = J(C) + dirichlet_lpdf(p | alpha),
//     fused_induced_dirichlet_prior_of_ordered_cutpoints_lpdf(C | alpha, link, location)
//       = the same at the latent location: sum over the rows r of C of J(C_r - location) +
//         dirichlet_lpdf(p_r | alpha), with p_r the category probabilities of C_r - location,
//
// with, for each set of K ordered cutpoints c (a vector C, or each row of a matrix C: one set per group, e.g.
// the random cutpoints of each study),
//
//     F_k = Phi(c_k) when link = 1, inv_logit(c_k) otherwise          (k = 1, ..., K),
//     p   = (F_1, F_2 - F_1, ..., F_K - F_(K-1), 1 - F_K)             (the K + 1 category probabilities),
//     J(c) = sum_k log f(c_k): std_normal_lpdf(c) when link = 1 (its -0.5 log(2 pi) per value left out with
//            propto), sum(log_inv_logit(c) + log1m_inv_logit(c)) otherwise (always kept),
//
// and alpha a vector of K + 1 Dirichlet weights, data (a prior) or parameters (a between-group model of the
// cutpoints, e.g. alpha = 0.01 + phi kappa). The induced Dirichlet prior (or between-group model) of the
// cutpoints of an ordinal model (Betancourt's): a Dirichlet on the category probabilities at a latent location
// (0 in the first form), with the log Jacobian of the map from the cutpoints to the cumulative probabilities,
// J. A general function: any cumulative-link (probit or logit) ordinal model with an induced Dirichlet on its
// cutpoints, fixed or random over groups (ordinal regression, ordinal meta-analysis and network meta-analysis of
// test accuracy, ordinal latent class models). The same values and derivatives as the Stan statements, for each
// set c of cutpoints,
//
//     vector[K] F = (link == 1) ? Phi(c - location) : inv_logit(c - location);
//     vector[K + 1] p = cumul_probs_to_ord_probs(F);
//     target += (link == 1) ? std_normal_lupdf(c - location)
//                           : sum(log_inv_logit(c - location) + log1m_inv_logit(c - location));
//     target += dirichlet_lupdf(p | alpha);
//
// (with propto; every constant, lgamma(sum(alpha)) - sum(lgamma(alpha)) of a data alpha and -0.5 log(2 pi) per
// value, without it, as Stan's _lpdf; the lgamma terms of a parameter alpha are always kept, as Stan's), in ONE
// autodiff node for all the sets in place of about K + 6 per set: BayesMVP's SIMD kernels for Phi, inv_logit,
// exp, log, log_inv_logit and log1m_inv_logit (the same kernels as fused_Phi, fused_inv_logit,
// fused_log_inv_logit and fused_log1m_inv_logit), and the element-wise arithmetic in branch-free loops (each
// category reads only its own cumulative probability and its predecessor's, so there is no recurrence).
// d / d c_k = (d lp / d F_k) f(c_k) + d J / d c_k, with d lp / d F_k = (alpha_k - 1) / p_k - (alpha_(k+1) - 1)
// / p_(k+1), f = phi or F (1 - F), and d J / d c_k = - c_k (probit) or 1 - 2 F_k (logit) (c at the location);
// d / d alpha_j = sum over the sets of log p_j, plus (number of sets) (digamma(sum(alpha)) - digamma(alpha_j))
// (the K + 1 digamma with Stan's scalar function, once per call).
//
// Declared in the Stan functions block (no body; stanc --allow-undefined) as any of
//
//     real fused_induced_dirichlet_prior_of_ordered_cutpoints_lpdf(vector C, data vector alpha, data int link);
//     real fused_induced_dirichlet_prior_of_ordered_cutpoints_lpdf(vector C, vector alpha, data int link,
//                                                                 data real location);
//     real fused_induced_dirichlet_prior_of_ordered_cutpoints_lpdf(row_vector C, vector alpha, data int link,
//                                                                 data real location);
//     real fused_induced_dirichlet_prior_of_ordered_cutpoints_lpdf(matrix C, vector alpha, data int link,
//                                                                 data real location);
//
// It uses the namespace fused_stan_math and the kernels of stan_external_functions.hpp, which it includes
// (#pragma once). Its helpers are in their own namespace, fused_stan_math::induced_dirichlet_prior.
//
// -----------------------------------------------------------------------------------------------------------

#include "stan_external_functions_on_Stan_math_only.hpp"

#include <stan/math/rev.hpp>
#include <stan/math/prim.hpp>

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <ostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace fused_stan_math {
namespace induced_dirichlet_prior {

static constexpr const char* function_name = "fused_induced_dirichlet_prior_of_ordered_cutpoints_lpdf";
static constexpr double half_log_of_two_pi = 0.91893853320467274178;




// ---- the forward pass of one set of K cutpoints x (the cutpoints minus the location), into buffers: lF (K: log F(x),
//      log Phi(x) or log_inv_logit(x)), lS (K + 1, the first K used: log(1 - F(x)), from the same kernel at -x),
//      log_p (2 K + 1: the log of the K + 1 category probabilities, then K values of work space), J_terms (K: log
//      f(x_k) without the constant); its log probability, without the lgamma terms of alpha (added once by the
//      caller), into sum. Every category probability is computed on the LOG scale, with BayesMVP's SIMD kernels
//      (fused_log_Phi with its log-scale tails, log_inv_logit, log1m_exp), no clamps:
//          log p_1     = log F(x_1),      log p_(K+1) = log(1 - F(x_K)),
//          log p_k     = log F(x_k) + log1m_exp(log F(x_(k-1)) - log F(x_k))          (x_(k-1) < 0),
//                      = log S(x_(k-1)) + log1m_exp(log S(x_k) - log S(x_(k-1)))      (x_(k-1) >= 0, S = 1 - F),
//      so that a category probability far in a tail, or between two close cutpoints, is never rounded to 0 (as
//      Phi(x_k) - Phi(x_(k-1)) in double precision was, which threw and left the chain without a gradient); two
//      equal cutpoints give log p = -inf, i.e. a log density of -inf (a rejected proposal), never an error. A term
//      with alpha_j = 1 is 0 and is left out (0 x log p). Returns false, and nothing else, only when the cutpoints
//      are not in increasing order or a value is not a number (the caller then runs Stan's checks).
template <bool include_the_constant_of_the_normal>
inline bool forward_pass(const double* x, const double* alpha, const int link, const int K, double* lF,
                         double* lS, double* log_p, double* J_terms, double& sum) {

    copy_n(x, K, lF);
    for (int k = 0; k < K; ++k) lS[k] = -x[k];
    if (link == 1) {
        k_log_Phi(lF, K);
        k_log_Phi(lS, K);
    } else {
        k_log_inv_logit(lF, K);
        k_log_inv_logit(lS, K);
    }
    log_p[0] = lF[0];
    log_p[K] = lS[K - 1];
    double* work = log_p + (K + 1);
    for (int k = 1; k < K; ++k) {
        work[k - 1] = (x[k - 1] >= 0.0) ? (lS[k] - lS[k - 1]) : (lF[k - 1] - lF[k]);
        if (!(work[k - 1] <= 0.0)) return false;
    }
    if (K > 1) k_log1m_exp(work, K - 1);
    for (int k = 1; k < K; ++k) log_p[k] = ((x[k - 1] >= 0.0) ? lS[k - 1] : lF[k]) + work[k - 1];
    double total = 0.0;
    for (int j = 0; j <= K; ++j) {
        if (std::isnan(log_p[j])) return false;
        if (alpha[j] != 1.0) total += (alpha[j] - 1.0) * log_p[j];
    }
    //// (J: the log density of F at each cutpoint: -x^2 / 2 (probit), log F + log(1 - F) (logit))
    if (link == 1) {
        for (int k = 0; k < K; ++k) J_terms[k] = -0.5 * x[k] * x[k];
        if constexpr (include_the_constant_of_the_normal) total -= K * half_log_of_two_pi;
    } else {
        for (int k = 0; k < K; ++k) J_terms[k] = lF[k] + lS[k];
    }
    for (int k = 0; k < K; ++k) total += J_terms[k];
    sum = total;
    return true;

}




// ---- the reverse pass of one set: d / d x_k (adjoint: the node's adjoint), from the forward pass's buffers, with
//      f(x_k) / p_j = exp(log f(x_k) - log p_j) (the SIMD exp kernel), so that a tiny p_j never divides 0 by 0:
//          d / d x_k = (alpha_k - 1) f(x_k) / p_k - (alpha_(k+1) - 1) f(x_k) / p_(k+1) + d / d x_k log f(x_k);
//      ratio_to_own_category and ratio_to_next_category: K values of work space each, made once by the caller.
inline void reverse_pass(const double* x, const double* alpha, const int link, const int K, const double* lF,
                         const double* lS, const double* log_p, const double adjoint, double* d_x,
                         double* ratio_to_own_category, double* ratio_to_next_category) {

    for (int k = 0; k < K; ++k) {
        const double log_f = (link == 1) ? (-0.5 * x[k] * x[k] - half_log_of_two_pi) : (lF[k] + lS[k]);
        ratio_to_own_category[k] = log_f - log_p[k];
        ratio_to_next_category[k] = log_f - log_p[k + 1];
    }
    k_exp(ratio_to_own_category, K);
    k_exp(ratio_to_next_category, K);
    for (int k = 0; k < K; ++k) {
        double d_F = 0.0;
        if (alpha[k] != 1.0) d_F += (alpha[k] - 1.0) * ratio_to_own_category[k];
        if (alpha[k + 1] != 1.0) d_F -= (alpha[k + 1] - 1.0) * ratio_to_next_category[k];
        const double d_J = (link == 1) ? -x[k] : (std::exp(lS[k]) - std::exp(lF[k]));
        d_x[k] = adjoint * (d_F + d_J);
    }

}




// ---- Stan's checks (they throw; reached only when a category probability was not in (0, 1])
template <typename T_p, typename T_alpha>
inline void run_the_checks_of_Stan(const T_p& p_values, const T_alpha& alpha) {

    stan::math::check_positive(function_name, "prior size parameter", alpha);
    stan::math::check_simplex(function_name, "probabilities", p_values);
    throw std::domain_error(std::string(function_name) + ": a category probability is not positive");

}

}  // namespace induced_dirichlet_prior
}  // namespace fused_stan_math




// ---- the Stan-facing function, at a latent location (global scope; stanc calls it as
//      fused_induced_dirichlet_prior_of_ordered_cutpoints_lpdf<propto__>(C, alpha, link, location, pstream__)):
//      C a vector or row_vector (one set of cutpoints) or a matrix (one set per row); alpha data or parameters
namespace fused_stan_math {

template <bool propto, typename T_C, typename T_alpha, stan::require_eigen_t<T_C>* = nullptr,
          stan::require_eigen_col_vector_t<T_alpha>* = nullptr>
inline stan::return_type_t<T_C, T_alpha> fused_induced_dirichlet_prior_of_ordered_cutpoints_lpdf(
        const T_C& C, const T_alpha& alpha, const int link, const double location, std::ostream* /*pstream__*/) {

    namespace prior = fused_stan_math::induced_dirichlet_prior;
    using stan::math::var;
    constexpr bool C_is_one_set = stan::is_eigen_vector<T_C>::value;
    constexpr bool C_is_constant = stan::is_constant_all<T_C>::value;
    constexpr bool alpha_is_constant = stan::is_constant_all<T_alpha>::value;
    //// (the lgamma terms of alpha: left out with propto only when alpha is data, as Stan's dirichlet_lupdf)
    constexpr bool include_the_lgamma_terms = !propto || !alpha_is_constant;
    constexpr bool include_the_constant_of_the_normal = !propto;
    const int R = C_is_one_set ? 1 : static_cast<int>(C.rows());
    const int K = C_is_one_set ? static_cast<int>(C.size()) : static_cast<int>(C.cols());
    stan::math::check_size_match(prior::function_name, "Number of prior weights", alpha.size(),
                                 "Number of cutpoints plus 1", K + 1);
    if ((K == 0) || (R == 0)) return 0.0;
    if constexpr (propto && C_is_constant && alpha_is_constant) return 0.0;
    const auto& C_ref = stan::math::to_ref(C);
    const auto& alpha_ref = stan::math::to_ref(alpha);
    //// (the buffers: x = C - location (row by row), alpha, F (log F), p (log(1 - F)) and log p of every set, kept for the reverse
    ////  pass on the arena when there is one (J_terms for the forward pass only, its space reused for d_x), then
    ////  the reverse pass's two K-value work spaces, used by each set in turn; in a std::vector when nothing is
    ////  differentiated, so that no evaluation of the log density alone leaves memory on the arena)
    const std::size_t RK = static_cast<std::size_t>(R) * K;
    const std::size_t n_buffer_values = 3 * RK + (K + 1) + static_cast<std::size_t>(R) * (K + 1) +
                                        static_cast<std::size_t>(R) * (2 * K + 1) +
                                        2 * static_cast<std::size_t>(K);
    std::vector<double> buffer_without_autodiff;
    double* buffers = nullptr;
    if constexpr (C_is_constant && alpha_is_constant) {
        buffer_without_autodiff.resize(n_buffer_values);
        buffers = buffer_without_autodiff.data();
    } else {
        buffers = stan::math::ChainableStack::instance_->memalloc_.alloc_array<double>(n_buffer_values);
    }
    double* x = buffers;
    double* alpha_kept = x + RK;
    double* F = alpha_kept + (K + 1);
    double* p = F + RK;
    double* log_p = p + static_cast<std::size_t>(R) * (K + 1);
    double* J_terms = log_p + static_cast<std::size_t>(R) * (2 * K + 1);
    double* ratios = J_terms + RK;
    for (int r = 0; r < R; ++r) {
        for (int k = 0; k < K; ++k) {
            const double value_of_C = C_is_one_set ? stan::math::value_of(C_ref.coeff(k))
                                                   : stan::math::value_of(C_ref.coeff(r, k));
            x[static_cast<std::size_t>(r) * K + k] = value_of_C - location;
        }
    }
    for (int j = 0; j <= K; ++j) alpha_kept[j] = stan::math::value_of(alpha_ref.coeff(j));
    double sum = 0.0;
    for (int r = 0; r < R; ++r) {
        double sum_of_set = 0.0;
        double* p_r = p + static_cast<std::size_t>(r) * (K + 1);
        const std::size_t start_of_set = static_cast<std::size_t>(r) * K;
        if (!prior::forward_pass<include_the_constant_of_the_normal>(
                    x + start_of_set, alpha_kept, link, K, F + start_of_set, p_r,
                    log_p + static_cast<std::size_t>(r) * (2 * K + 1), J_terms + start_of_set, sum_of_set)) {

            Eigen::Matrix<double, -1, 1> p_of_the_set(K + 1);
            const double* log_p_r = log_p + static_cast<std::size_t>(r) * (2 * K + 1);
            for (int j = 0; j <= K; ++j) p_of_the_set(j) = std::exp(log_p_r[j]);
            prior::run_the_checks_of_Stan(p_of_the_set, Eigen::Map<const Eigen::Matrix<double, -1, 1>>(alpha_kept,
                                                                                                       K + 1));

        }
        sum += sum_of_set;
    }
    if constexpr (include_the_lgamma_terms) {
        double sum_of_alpha = 0.0;
        double sum_of_lgamma_of_alpha = 0.0;
        for (int j = 0; j <= K; ++j) {
            sum_of_alpha += alpha_kept[j];
            sum_of_lgamma_of_alpha += std::lgamma(alpha_kept[j]);
        }
        sum += R * (std::lgamma(sum_of_alpha) - sum_of_lgamma_of_alpha);
    }
    if constexpr (C_is_constant && alpha_is_constant) {
        return sum;
    } else {
        using C_vars_type = std::conditional_t<C_is_constant,
                                               Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>,
                                               stan::plain_type_t<stan::math::promote_scalar_t<var, T_C>>>;
        using alpha_vars_type = Eigen::Matrix<var, Eigen::Dynamic, 1>;
        stan::arena_t<C_vars_type> C_vars;
        stan::arena_t<alpha_vars_type> alpha_vars;
        if constexpr (!C_is_constant) C_vars = C_ref;
        if constexpr (!alpha_is_constant) alpha_vars = alpha_ref;
        return stan::math::make_callback_var(sum, [C_vars, alpha_vars, x, alpha_kept, link, R, K, F, p, log_p,
                                                   J_terms, ratios](auto& vi) mutable {
            if constexpr (!C_is_constant) {
                for (int r = 0; r < R; ++r) {
                    double* d_x = J_terms + static_cast<std::size_t>(r) * K;
                    prior::reverse_pass(x + static_cast<std::size_t>(r) * K, alpha_kept, link, K,
                                        F + static_cast<std::size_t>(r) * K,
                                        p + static_cast<std::size_t>(r) * (K + 1),
                                        log_p + static_cast<std::size_t>(r) * (2 * K + 1), vi.adj(), d_x,
                                        ratios, ratios + K);
                    for (int k = 0; k < K; ++k) {
                        if constexpr (C_is_one_set) {
                            C_vars.coeffRef(k).adj() += d_x[k];
                        } else {
                            C_vars.coeffRef(r, k).adj() += d_x[k];
                        }
                    }
                }
            }
            if constexpr (!alpha_is_constant) {
                double sum_of_alpha = 0.0;
                for (int j = 0; j <= K; ++j) sum_of_alpha += alpha_kept[j];
                const double digamma_of_the_sum = stan::math::digamma(sum_of_alpha);
                for (int j = 0; j <= K; ++j) {
                    double sum_of_log_p = 0.0;
                    for (int r = 0; r < R; ++r) {
                        sum_of_log_p += log_p[static_cast<std::size_t>(r) * (2 * K + 1) + j];
                    }
                    const double d_alpha = sum_of_log_p +
                                           R * (digamma_of_the_sum - stan::math::digamma(alpha_kept[j]));
                    alpha_vars.coeffRef(j).adj() += vi.adj() * d_alpha;
                }
            }
        });
    }

}

}  // namespace fused_stan_math

using fused_stan_math::fused_induced_dirichlet_prior_of_ordered_cutpoints_lpdf;




// ---- the Stan-facing function at the latent location 0 (stanc calls it as
//      fused_induced_dirichlet_prior_of_ordered_cutpoints_lpdf<propto__>(C, alpha, link, pstream__))
namespace fused_stan_math {

template <bool propto, typename T_C, typename T_alpha, stan::require_eigen_t<T_C>* = nullptr,
          stan::require_eigen_col_vector_t<T_alpha>* = nullptr>
inline stan::return_type_t<T_C, T_alpha> fused_induced_dirichlet_prior_of_ordered_cutpoints_lpdf(
        const T_C& C, const T_alpha& alpha, const int link, std::ostream* pstream__) {

    return fused_induced_dirichlet_prior_of_ordered_cutpoints_lpdf<propto>(C, alpha, link, 0.0, pstream__);

}

// ---- Ragged multi-group (varying cutpoints per group) induced-Dirichlet in ONE autodiff node ----
template <bool propto, typename T_C, typename T_alpha, typename T_loc,
          stan::require_eigen_col_vector_t<T_C>* = nullptr,
          stan::require_eigen_col_vector_t<T_alpha>* = nullptr>
inline stan::return_type_t<T_C, T_alpha>
fused_induced_dirichlet_prior_of_ordered_cutpoints_ragged_lpdf(
        const T_C& C,
        const T_alpha& alpha,
        const std::vector<int>& n_cutpoints_per_group,
        const std::vector<int>& cutpoint_k_of_each_C,
        const int link,
        const T_loc& location_arg,
        std::ostream* /*pstream__*/) {

    const double location = stan::math::value_of(location_arg);

    namespace prior = fused_stan_math::induced_dirichlet_prior;
    using stan::math::var;
    static constexpr const char* fn_name = "fused_induced_dirichlet_prior_of_ordered_cutpoints_ragged_lpdf";
    constexpr bool C_is_constant     = stan::is_constant_all<T_C>::value;
    constexpr bool alpha_is_constant = stan::is_constant_all<T_alpha>::value;
    constexpr bool include_the_lgamma_terms = !propto || !alpha_is_constant;
    constexpr bool include_the_constant_of_the_normal = !propto;

    const int n_cat = static_cast<int>(alpha.size());
    const int n_C   = static_cast<int>(C.size());
    const int R     = static_cast<int>(n_cutpoints_per_group.size());
    if (n_cat < 2) {
        throw std::domain_error(std::string(fn_name) + ": alpha must have length >= 2");
    }
    stan::math::check_size_match(fn_name, "cutpoint_k_of_each_C", cutpoint_k_of_each_C.size(), "C", n_C);

    int sum_K = 0;
    int R_nonempty = 0;
    for (int r = 0; r < R; ++r) {
        const int Kr = n_cutpoints_per_group[r];
        if (Kr < 0 || Kr > n_cat - 1) {
            throw std::domain_error(std::string(fn_name) + ": invalid n_cutpoints_per_group entry");
        }
        if (Kr > 0) {
            ++R_nonempty;
            for (int m = 0; m < Kr; ++m) {
                const int k = cutpoint_k_of_each_C[sum_K + m];
                if (k < 1 || k > n_cat - 1 || (m > 0 && k <= cutpoint_k_of_each_C[sum_K + m - 1])) {
                    throw std::domain_error(std::string(fn_name) + ": cutpoint_k_of_each_C must rise in 1..n_thr");
                }
            }
            sum_K += Kr;
        }
    }
    if (sum_K != n_C) {
        throw std::domain_error(std::string(fn_name) + ": sum(n_cutpoints_per_group) != size(C)");
    }
    if (n_C == 0 || R_nonempty == 0) return 0.0;
    if constexpr (propto && C_is_constant && alpha_is_constant) return 0.0;

    const auto& C_ref     = stan::math::to_ref(C);
    const auto& alpha_ref = stan::math::to_ref(alpha);

    const int n_total_segs = n_C + R_nonempty;
    const int n_interior   = n_C - R_nonempty;

    // Double buffers: x (n_C), lF (n_C), lS (n_C), log_p (n_total_segs), seg_alpha (n_total_segs),
    //                 u_seg_alpha (n_total_segs), u_lgamma (n_total_segs), u_digamma (n_total_segs),
    //                 alpha_kept (n_cat), cum_alpha (n_cat+1), ratio_own (n_C), ratio_next (n_C), work (max(1, n_interior))
    const std::size_t n_dbl_bufs = static_cast<std::size_t>(5) * n_C + static_cast<std::size_t>(5) * n_total_segs +
                                   static_cast<std::size_t>(2 * n_cat + 1) + static_cast<std::size_t>(std::max(1, n_interior));
    const std::size_t n_int_bufs = static_cast<std::size_t>(3) * n_total_segs + static_cast<std::size_t>(2) * n_C;
    std::vector<double> dbl_prim;
    std::vector<int> int_prim;
    double* d_buf = nullptr;
    int* i_buf = nullptr;
    if constexpr (C_is_constant && alpha_is_constant) {
        dbl_prim.resize(n_dbl_bufs);
        int_prim.resize(n_int_bufs);
        d_buf = dbl_prim.data();
        i_buf = int_prim.data();
    } else {
        auto& arena = stan::math::ChainableStack::instance_->memalloc_;
        d_buf = arena.alloc_array<double>(n_dbl_bufs);
        i_buf = arena.alloc_array<int>(n_int_bufs);
    }

    double* x           = d_buf;
    double* lF          = x + n_C;
    double* lS          = lF + n_C;
    double* log_p       = lS + n_C;
    double* seg_alpha   = log_p + n_total_segs;
    double* u_seg_alpha = seg_alpha + n_total_segs;
    double* u_lgamma    = u_seg_alpha + n_total_segs;
    double* u_digamma   = u_lgamma + n_total_segs;
    double* alpha_kept  = u_digamma + n_total_segs;
    double* cum_alpha   = alpha_kept + n_cat;
    double* ratio_own   = cum_alpha + (n_cat + 1);
    double* ratio_next  = ratio_own + n_C;
    double* work        = ratio_next + n_C;

    int* seg_lo        = i_buf;
    int* seg_hi        = seg_lo + n_total_segs;
    int* seg_uid       = seg_hi + n_total_segs;
    int* seg_own_of_C  = seg_uid + n_total_segs;
    int* seg_next_of_C = seg_own_of_C + n_C;

    cum_alpha[0] = 0.0;
    for (int j = 0; j < n_cat; ++j) {
        const double aj = stan::math::value_of(alpha_ref.coeff(j));
        alpha_kept[j] = aj;
        if (!(aj > 0.0) || !std::isfinite(aj)) {
            stan::math::check_positive(fn_name, "prior size parameter",
                                       Eigen::Map<const Eigen::Matrix<double, -1, 1>>(alpha_kept, n_cat));
        }
        cum_alpha[j + 1] = cum_alpha[j] + aj;
    }
    const double alpha_0 = cum_alpha[n_cat];

    for (int p = 0; p < n_C; ++p) {
        const double xp = stan::math::value_of(C_ref.coeff(p)) - location;
        x[p]  = xp;
        lF[p] = xp;
        lS[p] = -xp;
    }
    if (link == 1) {
        k_log_Phi(lF, n_C);
        k_log_Phi(lS, n_C);
    } else {
        k_log_inv_logit(lF, n_C);
        k_log_inv_logit(lS, n_C);
    }

    int w_idx = 0;
    int off = 0;
    for (int r = 0; r < R; ++r) {
        const int Kr = n_cutpoints_per_group[r];
        if (Kr == 0) continue;
        for (int m = 1; m < Kr; ++m) {
            const int p = off + m;
            const double d = (x[p - 1] >= 0.0) ? (lS[p] - lS[p - 1]) : (lF[p - 1] - lF[p]);
            if (!(d <= 0.0)) {
                throw std::domain_error(std::string(fn_name) + ": cutpoints are not strictly increasing");
            }
            work[w_idx++] = d;
        }
        off += Kr;
    }
    if (n_interior > 0) {
        k_log1m_exp(work, n_interior);
    }

    // Stack lookup table for deduplicating unique (lo, hi) segments (handles n_cat <= 128 with zero heap allocs)
    static constexpr int STACK_MAX_NCAT = 128;
    int16_t seg_map_stack[(STACK_MAX_NCAT + 1) * (STACK_MAX_NCAT + 1)];
    const bool use_seg_dedup = (n_cat <= STACK_MAX_NCAT);
    const int stride = n_cat + 1;
    if (use_seg_dedup) {
        std::fill_n(seg_map_stack, stride * stride, static_cast<int16_t>(-1));
    }

    double sum = 0.0;
    int n_u = 0;
    w_idx = 0;
    off = 0;
    int s_idx = 0;
    for (int r = 0; r < R; ++r) {
        const int Kr = n_cutpoints_per_group[r];
        if (Kr == 0) continue;
        for (int s = 0; s <= Kr; ++s) {
            const int lo = (s == 0) ? 0 : cutpoint_k_of_each_C[off + s - 1];
            const int hi = (s == Kr) ? n_cat : cutpoint_k_of_each_C[off + s];
            double lp;
            if (s == 0) {
                lp = lF[off];
            } else if (s == Kr) {
                lp = lS[off + Kr - 1];
            } else {
                const int p = off + s;
                lp = ((x[p - 1] >= 0.0) ? lS[p - 1] : lF[p]) + work[w_idx++];
            }
            if (std::isnan(lp)) {
                throw std::domain_error(std::string(fn_name) + ": category log probability is NaN");
            }
            if (s < Kr) {
                seg_own_of_C[off + s]  = s_idx;
                seg_next_of_C[off + s] = s_idx + 1;
            }
            const double sa = cum_alpha[hi] - cum_alpha[lo];
            seg_lo[s_idx]    = lo;
            seg_hi[s_idx]    = hi;
            seg_alpha[s_idx] = sa;
            log_p[s_idx]     = lp;
            if (sa != 1.0) sum += (sa - 1.0) * lp;

            int uid;
            if (use_seg_dedup) {
                int16_t& cell = seg_map_stack[lo * stride + hi];
                if (cell < 0) {
                    uid = n_u++;
                    cell = static_cast<int16_t>(uid);
                    u_seg_alpha[uid] = sa;
                } else {
                    uid = cell;
                }
            } else {
                uid = n_u++;
                u_seg_alpha[uid] = sa;
            }
            seg_uid[s_idx] = uid;
            ++s_idx;
        }
        off += Kr;
    }

    if constexpr (include_the_lgamma_terms || !alpha_is_constant) {
        u_seg_alpha[n_u] = alpha_0;
        k_lgamma_and_digamma<!alpha_is_constant>(u_seg_alpha, n_u + 1, u_lgamma, u_digamma);
        if constexpr (include_the_lgamma_terms) {
            double sum_lgamma_seg = 0.0;
            for (int s_i = 0; s_i < n_total_segs; ++s_i) {
                sum_lgamma_seg += u_lgamma[seg_uid[s_i]];
            }
            sum += static_cast<double>(R_nonempty) * u_lgamma[n_u] - sum_lgamma_seg;
        }
    }
    if (link == 1) {
        double sum_sq = 0.0;
        for (int p = 0; p < n_C; ++p) sum_sq += x[p] * x[p];
        sum += -0.5 * sum_sq;
        if constexpr (include_the_constant_of_the_normal) {
            sum -= static_cast<double>(n_C) * prior::half_log_of_two_pi;
        }
    } else {
        for (int p = 0; p < n_C; ++p) sum += lF[p] + lS[p];
    }

    if constexpr (C_is_constant && alpha_is_constant) {
        return sum;
    } else {
        stan::arena_t<Eigen::Matrix<var, -1, 1>> C_vars;
        stan::arena_t<Eigen::Matrix<var, -1, 1>> alpha_vars;
        if constexpr (!C_is_constant) C_vars = C_ref;
        if constexpr (!alpha_is_constant) alpha_vars = alpha_ref;
        return stan::math::make_callback_var(
            sum,
            [C_vars, alpha_vars, x, lF, lS, log_p, seg_alpha, u_digamma,
             ratio_own, ratio_next, seg_lo, seg_hi, seg_uid, seg_own_of_C, seg_next_of_C,
             link, n_C, n_cat, n_total_segs, n_u](auto& vi) mutable {
                const double adj = vi.adj();
                if constexpr (!C_is_constant) {
                    for (int p = 0; p < n_C; ++p) {
                        const double log_f = (link == 1)
                            ? (-0.5 * x[p] * x[p] - prior::half_log_of_two_pi)
                            : (lF[p] + lS[p]);
                        ratio_own[p]  = log_f - log_p[seg_own_of_C[p]];
                        ratio_next[p] = log_f - log_p[seg_next_of_C[p]];
                    }
                    k_exp(ratio_own, n_C);
                    k_exp(ratio_next, n_C);
                    if (link != 1) {
                        k_exp(lF, n_C);
                        k_exp(lS, n_C);
                    }
                    for (int p = 0; p < n_C; ++p) {
                        const double a_own  = seg_alpha[seg_own_of_C[p]];
                        const double a_next = seg_alpha[seg_next_of_C[p]];
                        double d_F = 0.0;
                        if (a_own  != 1.0) d_F += (a_own  - 1.0) * ratio_own[p];
                        if (a_next != 1.0) d_F -= (a_next - 1.0) * ratio_next[p];
                        const double d_J = (link == 1) ? (-x[p]) : (lS[p] - lF[p]);
                        C_vars.coeffRef(p).adj() += adj * (d_F + d_J);
                    }
                }
                if constexpr (!alpha_is_constant) {
                    const double dg_sum = u_digamma[n_u];
                    double diff_stack[STACK_MAX_NCAT + 2];
                    std::vector<double> diff_heap;
                    double* diff = (n_cat <= STACK_MAX_NCAT) ? diff_stack : (diff_heap.resize(n_cat + 1), diff_heap.data());
                    std::fill_n(diff, n_cat + 1, 0.0);
                    for (int s_i = 0; s_i < n_total_segs; ++s_i) {
                        const double w = log_p[s_i] + dg_sum - u_digamma[seg_uid[s_i]];
                        diff[seg_lo[s_i]] += w;
                        diff[seg_hi[s_i]] -= w;
                    }
                    double running = 0.0;
                    for (int c = 0; c < n_cat; ++c) {
                        running += diff[c];
                        alpha_vars.coeffRef(c).adj() += adj * running;
                    }
                }
            });
    }
}

}  // namespace fused_stan_math

using fused_stan_math::fused_induced_dirichlet_prior_of_ordered_cutpoints_lpdf;
using fused_stan_math::fused_induced_dirichlet_prior_of_ordered_cutpoints_ragged_lpdf;

























// ---- Fused 2-class ragged Induced-Dirichlet + Binomial cumulative-link likelihood (precomputed segments, single C copy) ----
namespace fused_stan_math {

template <bool propto,
          typename T_loc,
          typename T_C,
          typename T_alpha,
          typename T_blp,
          stan::require_eigen_t<T_loc>* = nullptr,
          stan::require_eigen_col_vector_t<T_blp>* = nullptr>
inline stan::return_type_t<T_loc, T_C, T_alpha>
fused_ord_biv_randthr_ind_dir_and_binomial_lpmf(
        const std::vector<std::vector<std::vector<int>>>& x_2,
        const std::vector<std::vector<std::vector<int>>>& n_prev,
        const T_loc& Xbeta,
        const std::vector<Eigen::Matrix<T_C, -1, 1>>& C,
        const std::vector<Eigen::Matrix<T_alpha, -1, 1>>& alpha,
        const std::vector<std::vector<int>>& position_of_each_C,
        const std::vector<int>& n_obs_cutpoints,
        const std::vector<int>& n_reported_cutpoints,
        const std::vector<int>& id_u_seg_lo,
        const std::vector<int>& id_u_seg_hi,
        const std::vector<int>& id_seg_uid,
        const std::vector<int>& id_seg_own_of_C,
        const std::vector<int>& id_seg_next_of_C,
        const std::vector<int>& K_fold_CV_indicator,
        const int link,
        const T_blp& baseline_location_prior,
        const double bound,
        std::ostream* /*pstream__*/) {
    namespace prior  = fused_stan_math::induced_dirichlet_prior;
    namespace counts = fused_stan_math::binomial_counts_above_successive_cutpoints;
    using stan::math::var;
    static constexpr const char* fn_name = "fused_ord_biv_randthr_ind_dir_and_binomial_lpmf";

    if (static_cast<int>(C.size()) != 2 || static_cast<int>(alpha.size()) != 2 ||
        baseline_location_prior.size() != 2) {
        throw std::domain_error(std::string(fn_name) + ": C, alpha, baseline_location_prior must have length 2");
    }
    const int n_studies = static_cast<int>(n_reported_cutpoints.size());
    const int n_cat = static_cast<int>(alpha[0].size());
    const int n_C = static_cast<int>(C[0].size());
    if (static_cast<int>(alpha[1].size()) != n_cat || static_cast<int>(C[1].size()) != n_C) {
        throw std::domain_error(std::string(fn_name) + ": C[0]/C[1] or alpha[0]/alpha[1] size mismatch");
    }

    int R_nonempty = 0;
    for (int r = 0; r < n_studies; ++r) {
        if (n_reported_cutpoints[r] > 0) ++R_nonempty;
    }
    const int n_total_segs = static_cast<int>(id_seg_uid.size());
    const int n_u_segs = static_cast<int>(id_u_seg_lo.size());
    const int n_interior = n_C - R_nonempty;

    int K_bin = 0;
    for (int s = 0; s < static_cast<int>(n_obs_cutpoints.size()); ++s) {
        if (K_fold_CV_indicator[s] == 1) K_bin += 2 * n_obs_cutpoints[s];
    }

    constexpr bool C_is_const     = stan::is_constant_all<T_C>::value;
    constexpr bool alpha_is_const = stan::is_constant_all<T_alpha>::value;
    constexpr bool loc_is_const   = stan::is_constant_all<T_loc>::value;
    constexpr bool is_any_var     = !C_is_const || !alpha_is_const || !loc_is_const;
    constexpr bool include_lgamma = !propto || !alpha_is_const;
    constexpr bool include_norm_const = !propto;

    if constexpr (propto && !is_any_var) return 0.0;

    const auto& Xbeta_ref  = stan::math::to_ref(Xbeta);
    const auto& C0_ref     = stan::math::to_ref(C[0]);
    const auto& C1_ref     = stan::math::to_ref(C[1]);
    const auto& alpha0_ref = stan::math::to_ref(alpha[0]);
    const auto& alpha1_ref = stan::math::to_ref(alpha[1]);
    const auto& blp_ref    = stan::math::to_ref(baseline_location_prior);
    const double blp[2] = { stan::math::value_of(blp_ref.coeff(0)),
                            stan::math::value_of(blp_ref.coeff(1)) };

    const std::size_t exact_id_dbl =
        static_cast<std::size_t>(14) * n_C
        + static_cast<std::size_t>(std::max(1, 2 * n_interior))
        + static_cast<std::size_t>(2) * n_total_segs
        + static_cast<std::size_t>(6) * (n_u_segs + 1)
        + static_cast<std::size_t>(2) * (n_cat + 1);
    const std::size_t bin_dbl_count = static_cast<std::size_t>(7) * K_bin;
    const std::size_t tot_dbl = exact_id_dbl + bin_dbl_count;

    const std::size_t id_int_count =
        static_cast<std::size_t>(2) * n_u_segs + n_total_segs + 2 * n_C;
    const std::size_t bin_int_count = static_cast<std::size_t>(5) * K_bin;
    const std::size_t tot_int = id_int_count + bin_int_count;

    std::vector<double> d_prim;
    std::vector<int> i_prim;
    double* d_buf = nullptr;
    int* i_buf = nullptr;
    if constexpr (!is_any_var) {
        d_prim.resize(std::max<std::size_t>(1, tot_dbl));
        i_prim.resize(std::max<std::size_t>(1, tot_int));
        d_buf = d_prim.data();
        i_buf = i_prim.data();
    } else {
        auto& arena = stan::math::ChainableStack::instance_->memalloc_;
        d_buf = arena.alloc_array<double>(std::max<std::size_t>(1, tot_dbl));
        i_buf = arena.alloc_array<int>(std::max<std::size_t>(1, tot_int));
    }

    double* C_val       = d_buf;
    double* x           = C_val + 2 * n_C;
    double* lF          = x + 2 * n_C;
    double* lS          = lF + 2 * n_C;
    double* ratio_own   = lS + 2 * n_C;
    double* ratio_next  = ratio_own + 2 * n_C;
    double* dC_buf      = ratio_next + 2 * n_C;
    double* log_p       = dC_buf + 2 * n_C;
    double* u_seg_alpha = log_p + 2 * n_total_segs;
    double* u_lgamma    = u_seg_alpha + 2 * (n_u_segs + 1);
    double* u_digamma   = u_lgamma + 2 * (n_u_segs + 1);
    double* cum_alpha   = u_digamma + 2 * (n_u_segs + 1);
    double* work        = cum_alpha + 2 * (n_cat + 1);
    double* bin_d_base  = work + std::max(1, 2 * n_interior);

    int* seg_lo       = i_buf;                                 // n_u_segs
    int* seg_hi       = seg_lo + n_u_segs;                     // n_u_segs
    int* seg_uid      = seg_hi + n_u_segs;                     // n_total_segs
    int* seg_own      = seg_uid + n_total_segs;                // n_C
    int* seg_next     = seg_own + n_C;                         // n_C
    int* bin_i_base   = seg_next + n_C;                        // 5 * K_bin

    for (int u = 0; u < n_u_segs; ++u) {
        seg_lo[u] = id_u_seg_lo[u];
        seg_hi[u] = id_u_seg_hi[u];
    }
    for (int s_i = 0; s_i < n_total_segs; ++s_i) seg_uid[s_i] = id_seg_uid[s_i];
    for (int p = 0; p < n_C; ++p) {
        seg_own[p] = id_seg_own_of_C[p];
        seg_next[p] = id_seg_next_of_C[p];
    }

    double sum_total = 0.0;

    // 1) Read alpha[0..1], C[0..1], compute C - blp and batched lF/lS
    if (n_C > 0 && R_nonempty > 0) {
        for (int c = 0; c < 2; ++c) {
            const auto& a_ref = (c == 0) ? alpha0_ref : alpha1_ref;
            double* ca = cum_alpha + c * (n_cat + 1);
            ca[0] = 0.0;
            for (int j = 0; j < n_cat; ++j) {
                const double aj = stan::math::value_of(a_ref.coeff(j));
                if (!(aj > 0.0) || !std::isfinite(aj)) {
                    stan::math::check_positive(fn_name, "prior size parameter", aj);
                }
                ca[j + 1] = ca[j] + aj;
            }
            const double* ca_c = ca;
            double* usa = u_seg_alpha + c * (n_u_segs + 1);
            for (int u = 0; u < n_u_segs; ++u) {
                usa[u] = ca_c[seg_hi[u]] - ca_c[seg_lo[u]];
            }
            usa[n_u_segs] = ca_c[n_cat];

            const auto& c_ref = (c == 0) ? C0_ref : C1_ref;
            double* Cv = C_val + c * n_C;
            double* xc = x + c * n_C;
            double* lFc = lF + c * n_C;
            double* lSc = lS + c * n_C;
            const double loc_c = blp[c];
            for (int p = 0; p < n_C; ++p) {
                const double cp = stan::math::value_of(c_ref.coeff(p));
                const double xp = cp - loc_c;
                Cv[p] = cp;
                xc[p] = xp;
                lFc[p] = xp;
                lSc[p] = -xp;
            }
        }
        // 2 batched SIMD passes over ALL 4 * n_C values (lF[0..1], lS[0..1])
        if (link == 1) {
            k_log_Phi(lF, 4 * n_C); // lF and lS are contiguous (x->C_val precedes, lF is 2*n_C followed immediately by lS 2*n_C)
        } else {
            k_log_inv_logit(lF, 4 * n_C);
        }

        int w_idx = 0;
        for (int c = 0; c < 2; ++c) {
            const double* xc = x + c * n_C;
            const double* lFc = lF + c * n_C;
            const double* lSc = lS + c * n_C;
            int off = 0;
            for (int r = 0; r < n_studies; ++r) {
                const int Kr = n_reported_cutpoints[r];
                if (Kr == 0) continue;
                for (int m = 1; m < Kr; ++m) {
                    const int p = off + m;
                    const double d = (xc[p - 1] >= 0.0) ? (lSc[p] - lSc[p - 1]) : (lFc[p - 1] - lFc[p]);
                    if (!(d <= 0.0)) {
                        throw std::domain_error(std::string(fn_name) + ": cutpoints are not strictly increasing");
                    }
                    work[w_idx++] = d;
                }
                off += Kr;
            }
        }
        if (n_interior > 0) {
            k_log1m_exp(work, 2 * n_interior);
        }

        if constexpr (include_lgamma || !alpha_is_const) {
            k_lgamma_and_digamma<!alpha_is_const>(
                u_seg_alpha, 2 * (n_u_segs + 1), u_lgamma, u_digamma);
        }

        w_idx = 0;
        for (int c = 0; c < 2; ++c) {
            const double* xc = x + c * n_C;
            const double* lFc = lF + c * n_C;
            const double* lSc = lS + c * n_C;
            const double* usa = u_seg_alpha + c * (n_u_segs + 1);
            const double* ulg = u_lgamma + c * (n_u_segs + 1);
            double* lp_c = log_p + c * n_total_segs;

            int off = 0;
            int s_idx = 0;
            for (int r = 0; r < n_studies; ++r) {
                const int Kr = n_reported_cutpoints[r];
                if (Kr == 0) continue;
                for (int s = 0; s <= Kr; ++s) {
                    double lp;
                    if (s == 0) lp = lFc[off];
                    else if (s == Kr) lp = lSc[off + Kr - 1];
                    else {
                        const int p = off + s;
                        lp = ((xc[p - 1] >= 0.0) ? lSc[p - 1] : lFc[p]) + work[w_idx++];
                    }
                    if (std::isnan(lp)) {
                        throw std::domain_error(std::string(fn_name) + ": category log probability is NaN");
                    }
                    lp_c[s_idx] = lp;
                    const double sa = usa[seg_uid[s_idx]];
                    if (sa != 1.0) sum_total += (sa - 1.0) * lp;
                    ++s_idx;
                }
                off += Kr;
            }
            if constexpr (include_lgamma) {
                double sum_lg_seg = 0.0;
                for (int s_i = 0; s_i < n_total_segs; ++s_i) {
                    sum_lg_seg += ulg[seg_uid[s_i]];
                }
                sum_total += static_cast<double>(R_nonempty) * ulg[n_u_segs] - sum_lg_seg;
            }
        }
        if (link == 1) {
            double sum_sq = 0.0;
            for (int p = 0; p < 2 * n_C; ++p) sum_sq += x[p] * x[p];
            sum_total += -0.5 * sum_sq;
            if constexpr (include_norm_const) {
                sum_total -= static_cast<double>(2 * n_C) * prior::half_log_of_two_pi;
            }
        } else {
            for (int p = 0; p < 2 * n_C; ++p) sum_total += lF[p] + lS[p];
        }
    }

    // 2) Binomial counts forward pass reusing already-read C_val
    int* n_kept      = bin_i_base;
    int* n_prev_kept = bin_i_base + K_bin;
    int* loc_s       = bin_i_base + 2 * K_bin;
    int* loc_c       = bin_i_base + 3 * K_bin;
    int* C_pos       = bin_i_base + 4 * K_bin;
    double* latent_kept = bin_d_base;
    double* is_first    = bin_d_base + K_bin;
    double* S           = bin_d_base + 2 * K_bin;
    double* p_buf       = bin_d_base + 3 * K_bin;
    double* p_passes    = bin_d_base + 4 * K_bin;
    double* b_log_p     = bin_d_base + 5 * K_bin;
    double* b_log1m_p   = bin_d_base + 6 * K_bin;

    if (K_bin > 0) {
        int t = 0;
        const int n_st_obs = static_cast<int>(n_obs_cutpoints.size());
        for (int c = 0; c < 2; ++c) {
            const double* Cv = C_val + c * n_C;
            for (int s = 0; s < n_st_obs; ++s) {
                if (K_fold_CV_indicator[s] != 1) continue;
                const double loc_sc = stan::math::value_of(Xbeta_ref.coeff(s, c));
                for (int cut_i = 0; cut_i < n_obs_cutpoints[s]; ++cut_i) {
                    const int p = position_of_each_C[s][cut_i] - 1;
                    if (p < 0 || p >= n_C) {
                        throw std::domain_error(std::string(fn_name) + ": invalid position_of_each_C");
                    }
                    n_kept[t]      = x_2[c][s][cut_i];
                    n_prev_kept[t] = n_prev[c][s][cut_i];
                    loc_s[t]       = s;
                    loc_c[t]       = c;
                    C_pos[t]       = p;
                    is_first[t]    = (cut_i == 0) ? 1.0 : 0.0;
                    latent_kept[t] = loc_sc - Cv[p];
                    ++t;
                }
            }
        }
        double sum_bin = 0.0;
        if (!counts::forward_pass<!propto>(n_kept, n_prev_kept, latent_kept, is_first, link, bound, K_bin,
                                           S, p_buf, b_log_p, b_log1m_p, p_passes, sum_bin)) {
            counts::run_the_checks_of_Stan(
                std::vector<int>(n_kept, n_kept + K_bin),
                std::vector<int>(n_prev_kept, n_prev_kept + K_bin));
        }
        sum_total += sum_bin;
    }

    if constexpr (!is_any_var) {
        return sum_total;
    } else {
        stan::arena_t<Eigen::Matrix<var, -1, -1>> Xbeta_vars;
        stan::arena_t<Eigen::Matrix<var, -1, 1>> C_vars_0, C_vars_1;
        stan::arena_t<Eigen::Matrix<var, -1, 1>> alpha_vars_0, alpha_vars_1;
        if constexpr (!loc_is_const) Xbeta_vars = Xbeta_ref;
        if constexpr (!C_is_const) {
            C_vars_0 = C0_ref;
            C_vars_1 = C1_ref;
        }
        if constexpr (!alpha_is_const) {
            alpha_vars_0 = alpha0_ref;
            alpha_vars_1 = alpha1_ref;
        }
        return stan::math::make_callback_var(
            sum_total,
            [Xbeta_vars, C_vars_0, C_vars_1, alpha_vars_0, alpha_vars_1,
             x, lF, lS, ratio_own, ratio_next, dC_buf, log_p, u_seg_alpha, u_digamma,
             seg_lo, seg_hi, seg_uid, seg_own, seg_next,
             n_kept, n_prev_kept, loc_s, loc_c, C_pos, latent_kept, is_first,
             S, p_buf, p_passes, b_log_p, b_log1m_p,
             link, n_C, n_cat, n_total_segs, n_u_segs, R_nonempty, K_bin](auto& vi) mutable {
                const double adj = vi.adj();
                if constexpr (!C_is_const) {
                    std::fill_n(dC_buf, 2 * n_C, 0.0);
                }
                if (K_bin > 0) {
                    double* d_S = b_log_p;
                    double* d_latent = b_log1m_p;
                    counts::reverse_pass(n_kept, n_prev_kept, latent_kept, is_first, link, K_bin,
                                         S, p_buf, p_passes, adj, d_S, d_latent);
                    for (int t = 0; t < K_bin; ++t) {
                        const double dl = d_latent[t];
                        const int c = loc_c[t];
                        if constexpr (!loc_is_const) {
                            Xbeta_vars.coeffRef(loc_s[t], c).adj() += dl;
                        }
                        if constexpr (!C_is_const) {
                            dC_buf[c * n_C + C_pos[t]] -= dl;
                        }
                    }
                }
                if (n_C > 0 && R_nonempty > 0) {
                    if constexpr (!C_is_const) {
                        for (int c = 0; c < 2; ++c) {
                            const double* xc = x + c * n_C;
                            const double* lFc = lF + c * n_C;
                            const double* lSc = lS + c * n_C;
                            const double* lp_c = log_p + c * n_total_segs;
                            double* ro = ratio_own + c * n_C;
                            double* rn = ratio_next + c * n_C;
                            for (int p = 0; p < n_C; ++p) {
                                const double log_f = (link == 1)
                                    ? (-0.5 * xc[p] * xc[p] - prior::half_log_of_two_pi)
                                    : (lFc[p] + lSc[p]);
                                ro[p] = log_f - lp_c[seg_own[p]];
                                rn[p] = log_f - lp_c[seg_next[p]];
                            }
                        }
                        // Batched SIMD k_exp over BOTH classes at once (4*n_C in 2 contiguous arrays)
                        k_exp(ratio_own, 4 * n_C); // ratio_own (2*n_C) followed immediately by ratio_next (2*n_C)!
                        if (link != 1) {
                            k_exp(lF, 4 * n_C);
                        }
                        for (int c = 0; c < 2; ++c) {
                            const double* lFc = lF + c * n_C;
                            const double* lSc = lS + c * n_C;
                            const double* xc = x + c * n_C;
                            const double* ro = ratio_own + c * n_C;
                            const double* rn = ratio_next + c * n_C;
                            const double* usa = u_seg_alpha + c * (n_u_segs + 1);
                            double* dCc = dC_buf + c * n_C;
                            for (int p = 0; p < n_C; ++p) {
                                const double a_own  = usa[seg_uid[seg_own[p]]];
                                const double a_next = usa[seg_uid[seg_next[p]]];
                                double d_F = 0.0;
                                if (a_own  != 1.0) d_F += (a_own  - 1.0) * ro[p];
                                if (a_next != 1.0) d_F -= (a_next - 1.0) * rn[p];
                                const double d_J = (link == 1) ? (-xc[p]) : (lSc[p] - lFc[p]);
                                dCc[p] += adj * (d_F + d_J);
                            }
                        }
                        for (int p = 0; p < n_C; ++p) C_vars_0.coeffRef(p).adj() += dC_buf[p];
                        for (int p = 0; p < n_C; ++p) C_vars_1.coeffRef(p).adj() += dC_buf[n_C + p];
                    }
                    if constexpr (!alpha_is_const) {
                        double diff_stack[256];
                        std::vector<double> diff_heap;
                        double* diff = (n_cat + 1 <= 256) ? diff_stack : (diff_heap.resize(n_cat + 1), diff_heap.data());
                        for (int c = 0; c < 2; ++c) {
                            const double* lp_c = log_p + c * n_total_segs;
                            const double* udg = u_digamma + c * (n_u_segs + 1);
                            const double dg_sum = udg[n_u_segs];
                            std::fill_n(diff, n_cat + 1, 0.0);
                            for (int s_i = 0; s_i < n_total_segs; ++s_i) {
                                const int u = seg_uid[s_i];
                                const double w = lp_c[s_i] + dg_sum - udg[u];
                                diff[seg_lo[u]] += w;
                                diff[seg_hi[u]] -= w;
                            }
                            auto& av = (c == 0) ? alpha_vars_0 : alpha_vars_1;
                            double running = 0.0;
                            for (int k = 0; k < n_cat; ++k) {
                                running += diff[k];
                                av.coeffRef(k).adj() += adj * running;
                            }
                        }
                    }
                } else if constexpr (!C_is_const) {
                    for (int p = 0; p < n_C; ++p) C_vars_0.coeffRef(p).adj() += dC_buf[p];
                    for (int p = 0; p < n_C; ++p) C_vars_1.coeffRef(p).adj() += dC_buf[n_C + p];
                }
            });
    }
}

}  // namespace fused_stan_math

using fused_stan_math::fused_ord_biv_randthr_ind_dir_and_binomial_lpmf;


// ---- Fused 1-class (shared cutpoints C, alpha) ragged Induced-Dirichlet + HSROC Binomial cumulative-link likelihood ----
namespace fused_stan_math {

template <bool propto,
          typename T_loc,
          typename T_scale,
          typename T_C,
          typename T_alpha,
          stan::require_eigen_t<T_loc>* = nullptr,
          stan::require_eigen_t<T_scale>* = nullptr>
inline stan::return_type_t<T_loc, T_scale, T_C, T_alpha>
fused_ord_hsroc_randthr_ind_dir_and_binomial_lpmf(
        const std::vector<std::vector<std::vector<int>>>& x_2,
        const std::vector<std::vector<std::vector<int>>>& n_prev,
        const T_loc& Xlocations,
        const T_scale& Xscale,
        const Eigen::Matrix<T_C, -1, 1>& C,
        const Eigen::Matrix<T_alpha, -1, 1>& alpha,
        const std::vector<std::vector<int>>& position_of_each_C,
        const std::vector<int>& n_obs_cutpoints,
        const std::vector<int>& n_reported_cutpoints,
        const std::vector<int>& id_u_seg_lo,
        const std::vector<int>& id_u_seg_hi,
        const std::vector<int>& id_seg_uid,
        const std::vector<int>& id_seg_own_of_C,
        const std::vector<int>& id_seg_next_of_C,
        const std::vector<int>& K_fold_CV_indicator,
        const int link,
        const double baseline_location_prior,
        const double bound,
        std::ostream* /*pstream__*/) {
    namespace prior  = fused_stan_math::induced_dirichlet_prior;
    namespace counts = fused_stan_math::binomial_counts_above_successive_cutpoints;
    using stan::math::var;
    static constexpr const char* fn_name = "fused_ord_hsroc_randthr_ind_dir_and_binomial_lpmf";

    const int n_studies = static_cast<int>(n_reported_cutpoints.size());
    const int n_cat = static_cast<int>(alpha.size());
    const int n_C = static_cast<int>(C.size());

    int R_nonempty = 0;
    for (int r = 0; r < n_studies; ++r) {
        if (n_reported_cutpoints[r] > 0) ++R_nonempty;
    }
    const int n_total_segs = static_cast<int>(id_seg_uid.size());
    const int n_u_segs = static_cast<int>(id_u_seg_lo.size());
    const int n_interior = n_C - R_nonempty;

    int K_bin = 0;
    for (int s = 0; s < static_cast<int>(n_obs_cutpoints.size()); ++s) {
        if (K_fold_CV_indicator[s] == 1) K_bin += 2 * n_obs_cutpoints[s];
    }

    constexpr bool C_is_const     = stan::is_constant_all<T_C>::value;
    constexpr bool alpha_is_const = stan::is_constant_all<T_alpha>::value;
    constexpr bool loc_is_const   = stan::is_constant_all<T_loc>::value;
    constexpr bool scale_is_const = stan::is_constant_all<T_scale>::value;
    constexpr bool is_any_var     = !C_is_const || !alpha_is_const || !loc_is_const || !scale_is_const;
    constexpr bool include_lgamma = !propto || !alpha_is_const;
    constexpr bool include_norm_const = !propto;

    if constexpr (propto && !is_any_var) return 0.0;

    const auto& Xloc_ref   = stan::math::to_ref(Xlocations);
    const auto& Xscale_ref = stan::math::to_ref(Xscale);
    const auto& C_ref      = stan::math::to_ref(C);
    const auto& alpha_ref  = stan::math::to_ref(alpha);

    const std::size_t id_dbl_count =
        static_cast<std::size_t>(7) * n_C
        + static_cast<std::size_t>(std::max(1, n_interior))
        + static_cast<std::size_t>(n_total_segs)
        + static_cast<std::size_t>(3) * (n_u_segs + 1)
        + static_cast<std::size_t>(n_cat + 1);
    const std::size_t bin_dbl_count = static_cast<std::size_t>(9) * K_bin;
    const std::size_t tot_dbl = id_dbl_count + bin_dbl_count;

    const std::size_t id_int_count =
        static_cast<std::size_t>(2) * n_u_segs + n_total_segs + 2 * n_C;
    const std::size_t bin_int_count = static_cast<std::size_t>(5) * K_bin;
    const std::size_t tot_int = id_int_count + bin_int_count;

    std::vector<double> d_prim;
    std::vector<int> i_prim;
    double* d_buf = nullptr;
    int* i_buf = nullptr;
    if constexpr (!is_any_var) {
        d_prim.resize(std::max<std::size_t>(1, tot_dbl));
        i_prim.resize(std::max<std::size_t>(1, tot_int));
        d_buf = d_prim.data();
        i_buf = i_prim.data();
    } else {
        auto& arena = stan::math::ChainableStack::instance_->memalloc_;
        d_buf = arena.alloc_array<double>(std::max<std::size_t>(1, tot_dbl));
        i_buf = arena.alloc_array<int>(std::max<std::size_t>(1, tot_int));
    }

    double* C_val       = d_buf;
    double* x           = C_val + n_C;
    double* lF          = x + n_C;
    double* lS          = lF + n_C;
    double* ratio_own   = lS + n_C;
    double* ratio_next  = ratio_own + n_C;
    double* dC_buf      = ratio_next + n_C;
    double* log_p       = dC_buf + n_C;
    double* u_seg_alpha = log_p + n_total_segs;
    double* u_lgamma    = u_seg_alpha + (n_u_segs + 1);
    double* u_digamma   = u_lgamma + (n_u_segs + 1);
    double* cum_alpha   = u_digamma + (n_u_segs + 1);
    double* work        = cum_alpha + (n_cat + 1);
    double* bin_d_base  = work + std::max(1, n_interior);

    int* seg_lo       = i_buf;
    int* seg_hi       = seg_lo + n_u_segs;
    int* seg_uid      = seg_hi + n_u_segs;
    int* seg_own      = seg_uid + n_total_segs;
    int* seg_next     = seg_own + n_C;
    int* bin_i_base   = seg_next + n_C;

    for (int u = 0; u < n_u_segs; ++u) {
        seg_lo[u] = id_u_seg_lo[u];
        seg_hi[u] = id_u_seg_hi[u];
    }
    for (int s_i = 0; s_i < n_total_segs; ++s_i) seg_uid[s_i] = id_seg_uid[s_i];
    for (int p = 0; p < n_C; ++p) {
        seg_own[p] = id_seg_own_of_C[p];
        seg_next[p] = id_seg_next_of_C[p];
    }

    double sum_total = 0.0;

    if (n_C > 0 && R_nonempty > 0) {
        cum_alpha[0] = 0.0;
        for (int j = 0; j < n_cat; ++j) {
            const double aj = stan::math::value_of(alpha_ref.coeff(j));
            if (!(aj > 0.0) || !std::isfinite(aj)) {
                stan::math::check_positive(fn_name, "prior size parameter", aj);
            }
            cum_alpha[j + 1] = cum_alpha[j] + aj;
        }
        for (int u = 0; u < n_u_segs; ++u) {
            u_seg_alpha[u] = cum_alpha[seg_hi[u]] - cum_alpha[seg_lo[u]];
        }
        u_seg_alpha[n_u_segs] = cum_alpha[n_cat];

        for (int p = 0; p < n_C; ++p) {
            const double cp = stan::math::value_of(C_ref.coeff(p));
            const double xp = cp - baseline_location_prior;
            C_val[p] = cp;
            x[p]     = xp;
            lF[p]    = xp;
            lS[p]    = -xp;
        }
        if (link == 1) {
            k_log_Phi(lF, 2 * n_C);
        } else {
            k_log_inv_logit(lF, 2 * n_C);
        }

        int w_idx = 0;
        int off = 0;
        for (int r = 0; r < n_studies; ++r) {
            const int Kr = n_reported_cutpoints[r];
            if (Kr == 0) continue;
            for (int m = 1; m < Kr; ++m) {
                const int p = off + m;
                const double d = (x[p - 1] >= 0.0) ? (lS[p] - lS[p - 1]) : (lF[p - 1] - lF[p]);
                if (!(d <= 0.0)) {
                    throw std::domain_error(std::string(fn_name) + ": cutpoints are not strictly increasing");
                }
                work[w_idx++] = d;
            }
            off += Kr;
        }
        if (n_interior > 0) {
            k_log1m_exp(work, n_interior);
        }

        if constexpr (include_lgamma || !alpha_is_const) {
            k_lgamma_and_digamma<!alpha_is_const>(
                u_seg_alpha, n_u_segs + 1, u_lgamma, u_digamma);
        }

        w_idx = 0;
        off = 0;
        int s_idx = 0;
        for (int r = 0; r < n_studies; ++r) {
            const int Kr = n_reported_cutpoints[r];
            if (Kr == 0) continue;
            for (int s = 0; s <= Kr; ++s) {
                double lp;
                if (s == 0) lp = lF[off];
                else if (s == Kr) lp = lS[off + Kr - 1];
                else {
                    const int p = off + s;
                    lp = ((x[p - 1] >= 0.0) ? lS[p - 1] : lF[p]) + work[w_idx++];
                }
                if (std::isnan(lp)) {
                    throw std::domain_error(std::string(fn_name) + ": category log probability is NaN");
                }
                log_p[s_idx] = lp;
                const double sa = u_seg_alpha[seg_uid[s_idx]];
                if (sa != 1.0) sum_total += (sa - 1.0) * lp;
                ++s_idx;
            }
            off += Kr;
        }
        if constexpr (include_lgamma) {
            double sum_lg_seg = 0.0;
            for (int s_i = 0; s_i < n_total_segs; ++s_i) {
                sum_lg_seg += u_lgamma[seg_uid[s_i]];
            }
            sum_total += static_cast<double>(R_nonempty) * u_lgamma[n_u_segs] - sum_lg_seg;
        }
        if (link == 1) {
            double sum_sq = 0.0;
            for (int p = 0; p < n_C; ++p) sum_sq += x[p] * x[p];
            sum_total += -0.5 * sum_sq;
            if constexpr (include_norm_const) {
                sum_total -= static_cast<double>(n_C) * prior::half_log_of_two_pi;
            }
        } else {
            for (int p = 0; p < n_C; ++p) sum_total += lF[p] + lS[p];
        }
    }

    int* n_kept       = bin_i_base;
    int* n_prev_kept  = bin_i_base + K_bin;
    int* loc_s        = bin_i_base + 2 * K_bin;
    int* loc_c        = bin_i_base + 3 * K_bin;
    int* C_pos        = bin_i_base + 4 * K_bin;
    double* latent_kept = bin_d_base;
    double* is_first    = bin_d_base + K_bin;
    double* S           = bin_d_base + 2 * K_bin;
    double* p_buf       = bin_d_base + 3 * K_bin;
    double* p_passes    = bin_d_base + 4 * K_bin;
    double* b_log_p     = bin_d_base + 5 * K_bin;
    double* b_log1m_p   = bin_d_base + 6 * K_bin;
    double* scl_kept    = bin_d_base + 7 * K_bin;
    double* neg_CpL     = bin_d_base + 8 * K_bin;

    if (K_bin > 0) {
        int t = 0;
        const int n_st_obs = static_cast<int>(n_obs_cutpoints.size());
        for (int c = 0; c < 2; ++c) {
            for (int s = 0; s < n_st_obs; ++s) {
                if (K_fold_CV_indicator[s] != 1) continue;
                const double loc_sc = stan::math::value_of(Xloc_ref.coeff(s, c));
                const double scl_sc = stan::math::value_of(Xscale_ref.coeff(s, c));
                for (int cut_i = 0; cut_i < n_obs_cutpoints[s]; ++cut_i) {
                    const int p = position_of_each_C[s][cut_i] - 1;
                    if (p < 0 || p >= n_C) {
                        throw std::domain_error(std::string(fn_name) + ": invalid position_of_each_C");
                    }
                    const double cp = (n_C > 0 && R_nonempty > 0)
                        ? C_val[p] : stan::math::value_of(C_ref.coeff(p));
                    n_kept[t]      = x_2[c][s][cut_i];
                    n_prev_kept[t] = n_prev[c][s][cut_i];
                    loc_s[t]       = s;
                    loc_c[t]       = c;
                    C_pos[t]       = p;
                    is_first[t]    = (cut_i == 0) ? 1.0 : 0.0;
                    scl_kept[t]    = scl_sc;
                    neg_CpL[t]     = -(cp + loc_sc);
                    latent_kept[t] = neg_CpL[t] * scl_sc;
                    ++t;
                }
            }
        }
        double sum_bin = 0.0;
        if (!counts::forward_pass<!propto>(n_kept, n_prev_kept, latent_kept, is_first, link, bound, K_bin,
                                           S, p_buf, b_log_p, b_log1m_p, p_passes, sum_bin)) {
            counts::run_the_checks_of_Stan(
                std::vector<int>(n_kept, n_kept + K_bin),
                std::vector<int>(n_prev_kept, n_prev_kept + K_bin));
        }
        sum_total += sum_bin;
    }

    if constexpr (!is_any_var) {
        return sum_total;
    } else {
        stan::arena_t<Eigen::Matrix<var, -1, -1>> Xloc_vars, Xscale_vars;
        stan::arena_t<Eigen::Matrix<var, -1, 1>> C_vars, alpha_vars;
        if constexpr (!loc_is_const)   Xloc_vars   = Xloc_ref;
        if constexpr (!scale_is_const) Xscale_vars = Xscale_ref;
        if constexpr (!C_is_const)     C_vars      = C_ref;
        if constexpr (!alpha_is_const) alpha_vars  = alpha_ref;
        return stan::math::make_callback_var(
            sum_total,
            [Xloc_vars, Xscale_vars, C_vars, alpha_vars,
             x, lF, lS, ratio_own, ratio_next, dC_buf, log_p, u_seg_alpha, u_digamma,
             seg_lo, seg_hi, seg_uid, seg_own, seg_next,
             n_kept, n_prev_kept, loc_s, loc_c, C_pos, latent_kept, is_first,
             S, p_buf, p_passes, b_log_p, b_log1m_p, scl_kept, neg_CpL,
             link, n_C, n_cat, n_total_segs, n_u_segs, R_nonempty, K_bin](auto& vi) mutable {
                const double adj = vi.adj();
                if constexpr (!C_is_const) {
                    std::fill_n(dC_buf, n_C, 0.0);
                }
                if (K_bin > 0) {
                    double* d_S = b_log_p;
                    double* d_latent = b_log1m_p;
                    counts::reverse_pass(n_kept, n_prev_kept, latent_kept, is_first, link, K_bin,
                                         S, p_buf, p_passes, adj, d_S, d_latent);
                    for (int t = 0; t < K_bin; ++t) {
                        const double dl = d_latent[t];
                        const int s = loc_s[t];
                        const int c = loc_c[t];
                        const double dl_neg_scl = -dl * scl_kept[t];
                        if constexpr (!loc_is_const)   Xloc_vars.coeffRef(s, c).adj()   += dl_neg_scl;
                        if constexpr (!C_is_const)     dC_buf[C_pos[t]]                 += dl_neg_scl;
                        if constexpr (!scale_is_const) Xscale_vars.coeffRef(s, c).adj() += dl * neg_CpL[t];
                    }
                }
                if (n_C > 0 && R_nonempty > 0) {
                    if constexpr (!C_is_const) {
                        for (int p = 0; p < n_C; ++p) {
                            const double log_f = (link == 1)
                                ? (-0.5 * x[p] * x[p] - prior::half_log_of_two_pi)
                                : (lF[p] + lS[p]);
                            ratio_own[p]  = log_f - log_p[seg_own[p]];
                            ratio_next[p] = log_f - log_p[seg_next[p]];
                        }
                        k_exp(ratio_own, 2 * n_C);
                        if (link != 1) {
                            k_exp(lF, 2 * n_C);
                        }
                        for (int p = 0; p < n_C; ++p) {
                            const double a_own  = u_seg_alpha[seg_uid[seg_own[p]]];
                            const double a_next = u_seg_alpha[seg_uid[seg_next[p]]];
                            double d_F = 0.0;
                            if (a_own  != 1.0) d_F += (a_own  - 1.0) * ratio_own[p];
                            if (a_next != 1.0) d_F -= (a_next - 1.0) * ratio_next[p];
                            const double d_J = (link == 1) ? (-x[p]) : (lS[p] - lF[p]);
                            dC_buf[p] += adj * (d_F + d_J);
                        }
                        for (int p = 0; p < n_C; ++p) C_vars.coeffRef(p).adj() += dC_buf[p];
                    }
                    if constexpr (!alpha_is_const) {
                        double diff_stack[256];
                        std::vector<double> diff_heap;
                        double* diff = (n_cat + 1 <= 256) ? diff_stack : (diff_heap.resize(n_cat + 1), diff_heap.data());
                        const double dg_sum = u_digamma[n_u_segs];
                        std::fill_n(diff, n_cat + 1, 0.0);
                        for (int s_i = 0; s_i < n_total_segs; ++s_i) {
                            const int u = seg_uid[s_i];
                            const double w = log_p[s_i] + dg_sum - u_digamma[u];
                            diff[seg_lo[u]] += w;
                            diff[seg_hi[u]] -= w;
                        }
                        double running = 0.0;
                        for (int k = 0; k < n_cat; ++k) {
                            running += diff[k];
                            alpha_vars.coeffRef(k).adj() += adj * running;
                        }
                    }
                } else if constexpr (!C_is_const) {
                    for (int p = 0; p < n_C; ++p) C_vars.coeffRef(p).adj() += dC_buf[p];
                }
            });
    }
}

}  // namespace fused_stan_math

using fused_stan_math::fused_ord_hsroc_randthr_ind_dir_and_binomial_lpmf;


// ---- Precomputed-index 2-class ragged Induced-Dirichlet + SIMD-reverse Binomial cumulative-link likelihood ----
namespace fused_stan_math {
namespace obiv_rc_detail {

inline void counts_reverse_pass_simd(const int* n,
                                     const int* n_previous,
                                     const double* latent,
                                     const double* is_first,
                                     const int link,
                                     const int K,
                                     const double* S,
                                     const double* p,
                                     const double* p_passes,
                                     const double adjoint,
                                     double* d_S,
                                     double* d_latent,
                                     double* t2_tmp) {
    //// (the log-scale reverse pass of the binomial counts above successive cutpoints, the same for every
    ////  caller: finite where Phi underflows; t2_tmp is no longer needed)
    (void) t2_tmp;
    fused_stan_math::binomial_counts_above_successive_cutpoints::reverse_pass(
        n, n_previous, latent, is_first, link, K, S, p, p_passes, adjoint, d_S, d_latent);
    //// (the linear-scale pass it replaces, kept for reference:)
    // int k = 0;
    // for (; k < K; ++k) {
    //     const double succ = static_cast<double>(n[k]);
    //     const double fail = static_cast<double>(n_previous[k] - n[k]);
    //     const double d_p = adjoint * p_passes[k]
    //         * (((succ > 0.0) ? (succ / p[k]) : 0.0) - ((fail > 0.0) ? (fail / (1.0 - p[k])) : 0.0));
    //     const double S_prev = (k == 0) ? 1.0 : S[k - 1];
    //     const double S0 = S_prev * (1.0 - is_first[k]) + is_first[k];
    //     d_S[k] = d_p / S0;
    //     t2_tmp[k] = (k > 0) ? (d_p * (S[k] / S0) / S0 * (1.0 - is_first[k])) : 0.0;
    // }
    // for (int i = 0; i < K - 1; ++i) {
    //     d_S[i] -= t2_tmp[i + 1];
    // }
    // if (link == 1) {
    //     for (int i = 0; i < K; ++i) d_latent[i] = -0.5 * latent[i] * latent[i];
    //     fused_stan_math::k_exp(d_latent, K);
    //     for (int i = 0; i < K; ++i) d_latent[i] *= INV_SQRT_2PI * d_S[i];
    // } else {
    //     for (int i = 0; i < K; ++i) d_latent[i] = d_S[i] * S[i] * (1.0 - S[i]);
    // }
}

}  // namespace obiv_rc_detail

template <bool propto,
          typename T_loc,
          typename T_C,
          typename T_alpha,
          typename T_blp,
          stan::require_eigen_t<T_loc>* = nullptr,
          stan::require_eigen_col_vector_t<T_blp>* = nullptr>
inline stan::return_type_t<T_loc, T_C, T_alpha>
fused_ord_biv_randthr_ind_dir_and_binomial_precomputed_lpmf(
        const std::vector<int>& bin_n_in,
        const std::vector<int>& bin_n_prev_in,
        const std::vector<int>& bin_sc_in,
        const std::vector<int>& bin_Cpos_in,
        const Eigen::VectorXd& bin_is_first_in,
        const int K_bin,
        const int bin_is_seq,
        const T_loc& Xbeta,
        const std::vector<Eigen::Matrix<T_C, -1, 1>>& C,
        const std::vector<Eigen::Matrix<T_alpha, -1, 1>>& alpha,
        const std::vector<int>& id_u_seg_lo,
        const std::vector<int>& id_u_seg_hi,
        const std::vector<int>& id_u_seg_cnt,
        const std::vector<int>& id_seg_uid,
        const std::vector<int>& id_seg_own_of_C,
        const std::vector<int>& id_seg_next_of_C,
        const std::vector<int>& id_int_p,
        const std::vector<int>& id_seg_p,
        const std::vector<int>& id_seg_w,
        const int n_interior,
        const int R_nonempty,
        const int link,
        const T_blp& baseline_location_prior,
        const double bound,
        std::ostream* /*pstream__*/) {
    namespace prior = fused_stan_math::induced_dirichlet_prior;
    namespace counts = fused_stan_math::binomial_counts_above_successive_cutpoints;
    using stan::math::var;
    static constexpr const char* fn_name = "fused_ord_biv_randthr_ind_dir_and_binomial_precomputed_lpmf";

    const int n_studies = static_cast<int>(Xbeta.rows());
    const int n_cat = static_cast<int>(alpha[0].size());
    const int n_C = static_cast<int>(C[0].size());
    const int n_total_segs = n_C + R_nonempty;
    const int n_u_segs = static_cast<int>(id_u_seg_lo.size());

    constexpr bool C_is_const = stan::is_constant_all<T_C>::value;
    constexpr bool alpha_is_const = stan::is_constant_all<T_alpha>::value;
    constexpr bool loc_is_const = stan::is_constant_all<T_loc>::value;
    constexpr bool is_any_var = !C_is_const || !alpha_is_const || !loc_is_const;
    constexpr bool include_lgamma = !propto || !alpha_is_const;
    constexpr bool include_norm_const = !propto;

    if constexpr (propto && !is_any_var) return 0.0;

    const auto& Xbeta_ref = stan::math::to_ref(Xbeta);
    const auto& C0_ref = stan::math::to_ref(C[0]);
    const auto& C1_ref = stan::math::to_ref(C[1]);
    const auto& alpha0_ref = stan::math::to_ref(alpha[0]);
    const auto& alpha1_ref = stan::math::to_ref(alpha[1]);
    const auto& blp_ref = stan::math::to_ref(baseline_location_prior);
    const double blp[2] = {stan::math::value_of(blp_ref.coeff(0)), stan::math::value_of(blp_ref.coeff(1))};

    const std::size_t exact_id_dbl =
        static_cast<std::size_t>(14) * n_C
        + static_cast<std::size_t>(std::max(1, 2 * n_interior))
        + static_cast<std::size_t>(2) * n_total_segs
        + static_cast<std::size_t>(7) * (n_u_segs + 1)
        + static_cast<std::size_t>(2) * (n_cat + 1)
        + static_cast<std::size_t>(2) * std::max(1, n_studies);
    const std::size_t bin_dbl_count = static_cast<std::size_t>(8) * K_bin;
    const std::size_t tot_dbl = exact_id_dbl + bin_dbl_count;

    const std::size_t id_int_count =
        static_cast<std::size_t>(3) * n_u_segs + n_total_segs + 2 * n_C;
    const std::size_t bin_int_count = static_cast<std::size_t>(4) * K_bin;
    const std::size_t tot_int = id_int_count + bin_int_count;

    std::vector<double> d_prim;
    std::vector<int> i_prim;
    double* d_buf = nullptr;
    int* i_buf = nullptr;
    if constexpr (!is_any_var) {
        d_prim.resize(std::max<std::size_t>(1, tot_dbl));
        i_prim.resize(std::max<std::size_t>(1, tot_int));
        d_buf = d_prim.data();
        i_buf = i_prim.data();
    } else {
        auto& arena = stan::math::ChainableStack::instance_->memalloc_;
        d_buf = arena.alloc_array<double>(std::max<std::size_t>(1, tot_dbl));
        i_buf = arena.alloc_array<int>(std::max<std::size_t>(1, tot_int));
    }

    double* C_val = d_buf;
    double* x = C_val + 2 * n_C;
    double* lF = x + 2 * n_C;
    double* lS = lF + 2 * n_C;
    double* ratio_own = lS + 2 * n_C;
    double* ratio_next = ratio_own + 2 * n_C;
    double* dC_buf = ratio_next + 2 * n_C;
    double* log_p = dC_buf + 2 * n_C;
    double* u_seg_alpha = log_p + 2 * n_total_segs;
    double* u_lgamma = u_seg_alpha + 2 * (n_u_segs + 1);
    double* u_digamma = u_lgamma + 2 * (n_u_segs + 1);
    double* u_lp_sum = u_digamma + 2 * (n_u_segs + 1);
    double* cum_alpha = u_lp_sum + (n_u_segs + 1);
    double* d_Xb_tmp = cum_alpha + 2 * (n_cat + 1);
    double* work = d_Xb_tmp + 2 * std::max(1, n_studies);
    double* bin_d_base = work + std::max(1, 2 * n_interior);

    int* seg_lo = i_buf;
    int* seg_hi = seg_lo + n_u_segs;
    int* u_seg_cnt = seg_hi + n_u_segs;
    int* seg_uid = u_seg_cnt + n_u_segs;
    int* seg_own = seg_uid + n_total_segs;
    int* seg_next = seg_own + n_C;
    int* bin_i_base = seg_next + n_C;

    if (n_u_segs > 0) {
        std::memcpy(seg_lo, id_u_seg_lo.data(), static_cast<std::size_t>(n_u_segs) * sizeof(int));
        std::memcpy(seg_hi, id_u_seg_hi.data(), static_cast<std::size_t>(n_u_segs) * sizeof(int));
        std::memcpy(u_seg_cnt, id_u_seg_cnt.data(), static_cast<std::size_t>(n_u_segs) * sizeof(int));
    }
    if (n_total_segs > 0) {
        std::memcpy(seg_uid, id_seg_uid.data(), static_cast<std::size_t>(n_total_segs) * sizeof(int));
    }
    if (n_C > 0) {
        std::memcpy(seg_own, id_seg_own_of_C.data(), static_cast<std::size_t>(n_C) * sizeof(int));
        std::memcpy(seg_next, id_seg_next_of_C.data(), static_cast<std::size_t>(n_C) * sizeof(int));
    }

    double sum_total = 0.0;

    if (n_C > 0 && R_nonempty > 0) {
        for (int c = 0; c < 2; ++c) {
            const auto& a_ref = (c == 0) ? alpha0_ref : alpha1_ref;
            double* ca = cum_alpha + c * (n_cat + 1);
            ca[0] = 0.0;
            for (int j = 0; j < n_cat; ++j) {
                const double aj = stan::math::value_of(a_ref.coeff(j));
                if (!(aj > 0.0) || !std::isfinite(aj)) {
                    stan::math::check_positive(fn_name, "prior size parameter", aj);
                }
                ca[j + 1] = ca[j] + aj;
            }
            double* usa = u_seg_alpha + c * (n_u_segs + 1);
            for (int u = 0; u < n_u_segs; ++u) usa[u] = ca[seg_hi[u]] - ca[seg_lo[u]];
            usa[n_u_segs] = ca[n_cat];

            const auto& c_ref = (c == 0) ? C0_ref : C1_ref;
            double* Cv = C_val + c * n_C;
            double* xc = x + c * n_C;
            double* lFc = lF + c * n_C;
            double* lSc = lS + c * n_C;
            const double loc_c = blp[c];
            for (int p = 0; p < n_C; ++p) {
                const double cp = stan::math::value_of(c_ref.coeff(p));
                const double xp = cp - loc_c;
                Cv[p] = cp; xc[p] = xp; lFc[p] = xp; lSc[p] = -xp;
            }
        }
        if (link == 1) k_log_Phi(lF, 4 * n_C);
        else           k_log_inv_logit(lF, 4 * n_C);

        if (n_interior > 0) {
            const int* int_p = id_int_p.data();
            bool bad_ord = false;
            for (int c = 0; c < 2; ++c) {
                const double* xc = x + c * n_C;
                const double* lFc = lF + c * n_C;
                const double* lSc = lS + c * n_C;
                double* wk = work + c * n_interior;
                for (int wi = 0; wi < n_interior; ++wi) {
                    const int p = int_p[wi];
                    const double d = (xc[p - 1] >= 0.0) ? (lSc[p] - lSc[p - 1]) : (lFc[p - 1] - lFc[p]);
                    if (!(d <= 0.0)) bad_ord = true;
                    wk[wi] = d;
                }
            }
            if (bad_ord) throw std::domain_error(std::string(fn_name) + ": cutpoints are not strictly increasing");
            k_log1m_exp(work, 2 * n_interior);
        }
        if constexpr (include_lgamma || !alpha_is_const) {
            k_lgamma_and_digamma<!alpha_is_const>(u_seg_alpha, 2 * (n_u_segs + 1), u_lgamma, u_digamma);
        }

        const int* sp_arr = id_seg_p.data();
        const int* sw_arr = id_seg_w.data();
        bool bad_lp = false;
        for (int c = 0; c < 2; ++c) {
            const double* xc = x + c * n_C;
            const double* lFc = lF + c * n_C;
            const double* lSc = lS + c * n_C;
            const double* wk = work + c * n_interior;
            const double* usa = u_seg_alpha + c * (n_u_segs + 1);
            const double* ulg = u_lgamma + c * (n_u_segs + 1);
            double* lp_c = log_p + c * n_total_segs;
            for (int si = 0; si < n_total_segs; ++si) {
                const int p = sp_arr[si];
                const int wi = sw_arr[si];
                double lp;
                if (wi == -1)      lp = lFc[p];
                else if (wi == -2) lp = lSc[p];
                else               lp = ((xc[p - 1] >= 0.0) ? lSc[p - 1] : lFc[p]) + wk[wi];
                if (std::isnan(lp)) bad_lp = true;
                lp_c[si] = lp;
                const double sa = usa[seg_uid[si]];
                if (sa != 1.0) sum_total += (sa - 1.0) * lp;
            }
            if constexpr (include_lgamma) {
                double sum_lg_seg = 0.0;
                for (int u = 0; u < n_u_segs; ++u) sum_lg_seg += static_cast<double>(u_seg_cnt[u]) * ulg[u];
                sum_total += static_cast<double>(R_nonempty) * ulg[n_u_segs] - sum_lg_seg;
            }
        }
        if (bad_lp) throw std::domain_error(std::string(fn_name) + ": category log probability is NaN");

        if (link == 1) {
            double sum_sq = 0.0;
            for (int p = 0; p < 2 * n_C; ++p) sum_sq += x[p] * x[p];
            sum_total += -0.5 * sum_sq;
            if constexpr (include_norm_const) sum_total -= static_cast<double>(2 * n_C) * prior::half_log_of_two_pi;
        } else {
            for (int p = 0; p < 2 * n_C; ++p) sum_total += lF[p] + lS[p];
        }
    }

    int* n_kept = bin_i_base;
    int* n_prev_kept = bin_i_base + K_bin;
    int* bin_sc = bin_i_base + 2 * K_bin;
    int* C_pos = bin_i_base + 3 * K_bin;
    double* latent_kept = bin_d_base;
    double* is_first = bin_d_base + K_bin;
    double* S = bin_d_base + 2 * K_bin;
    double* p_buf = bin_d_base + 3 * K_bin;
    double* p_passes = bin_d_base + 4 * K_bin;
    double* b_log_p = bin_d_base + 5 * K_bin;
    double* b_log1m_p = bin_d_base + 6 * K_bin;
    double* t2_tmp = bin_d_base + 7 * K_bin;

    if (K_bin > 0) {
        std::memcpy(n_kept, bin_n_in.data(), static_cast<std::size_t>(K_bin) * sizeof(int));
        std::memcpy(n_prev_kept, bin_n_prev_in.data(), static_cast<std::size_t>(K_bin) * sizeof(int));
        std::memcpy(bin_sc, bin_sc_in.data(), static_cast<std::size_t>(K_bin) * sizeof(int));
        std::memcpy(C_pos, bin_Cpos_in.data(), static_cast<std::size_t>(K_bin) * sizeof(int));
        std::memcpy(is_first, bin_is_first_in.data(), static_cast<std::size_t>(K_bin) * sizeof(double));

        double xb_stack[512];
        std::vector<double> xb_heap;
        double* xb_val = (2 * n_studies <= 512) ? xb_stack : (xb_heap.resize(2 * n_studies), xb_heap.data());
        for (int c = 0; c < 2; ++c) {
            for (int s = 0; s < n_studies; ++s) xb_val[c * n_studies + s] = stan::math::value_of(Xbeta_ref.coeff(s, c));
        }
        if (bin_is_seq == 1) {
            for (int t = 0; t < K_bin; ++t) latent_kept[t] = xb_val[bin_sc[t]] - C_val[t];
        } else {
            for (int t = 0; t < K_bin; ++t) latent_kept[t] = xb_val[bin_sc[t]] - C_val[C_pos[t]];
        }
        double sum_bin = 0.0;
        if (!counts::forward_pass<!propto>(n_kept, n_prev_kept, latent_kept, is_first, link, bound, K_bin,
                                           S, p_buf, b_log_p, b_log1m_p, p_passes, sum_bin)) {
            counts::run_the_checks_of_Stan(
                std::vector<int>(n_kept, n_kept + K_bin),
                std::vector<int>(n_prev_kept, n_prev_kept + K_bin));
        }
        sum_total += sum_bin;
    }

    if constexpr (!is_any_var) {
        return sum_total;
    } else {
        stan::arena_t<Eigen::Matrix<var, -1, -1>> Xbeta_vars;
        stan::arena_t<Eigen::Matrix<var, -1, 1>> C_vars_0, C_vars_1;
        stan::arena_t<Eigen::Matrix<var, -1, 1>> alpha_vars_0, alpha_vars_1;
        if constexpr (!loc_is_const) Xbeta_vars = Xbeta_ref;
        if constexpr (!C_is_const) { C_vars_0 = C0_ref; C_vars_1 = C1_ref; }
        if constexpr (!alpha_is_const) { alpha_vars_0 = alpha0_ref; alpha_vars_1 = alpha1_ref; }
        return stan::math::make_callback_var(
            sum_total,
            [Xbeta_vars, C_vars_0, C_vars_1, alpha_vars_0, alpha_vars_1,
             x, lF, lS, ratio_own, ratio_next, dC_buf, log_p, u_seg_alpha, u_digamma, u_lp_sum, d_Xb_tmp,
             seg_lo, seg_hi, u_seg_cnt, seg_uid, seg_own, seg_next,
             n_kept, n_prev_kept, bin_sc, C_pos, latent_kept, is_first,
             S, p_buf, p_passes, b_log_p, b_log1m_p, t2_tmp,
             link, n_studies, n_C, n_cat, n_total_segs, n_u_segs, R_nonempty, K_bin, bin_is_seq](auto& vi) mutable {
                const double adj = vi.adj();
                if constexpr (!loc_is_const) std::fill_n(d_Xb_tmp, 2 * n_studies, 0.0);
                if (K_bin > 0) {
                    double* d_S = b_log_p;
                    double* d_latent = b_log1m_p;
                    obiv_rc_detail::counts_reverse_pass_simd(
                        n_kept, n_prev_kept, latent_kept, is_first, link, K_bin,
                        S, p_buf, p_passes, adj, d_S, d_latent, t2_tmp);
                    if constexpr (!loc_is_const) {
                        for (int t = 0; t < K_bin; ++t) d_Xb_tmp[bin_sc[t]] += d_latent[t];
                        for (int c = 0; c < 2; ++c) {
                            const int off = c * n_studies;
                            for (int s = 0; s < n_studies; ++s) {
                                Xbeta_vars.coeffRef(s, c).adj() += d_Xb_tmp[off + s];
                            }
                        }
                    }
                    if constexpr (!C_is_const) {
                        if (bin_is_seq == 1) {
                            for (int t = 0; t < 2 * n_C; ++t) dC_buf[t] = -d_latent[t];
                        } else {
                            std::fill_n(dC_buf, 2 * n_C, 0.0);
                            for (int t = 0; t < K_bin; ++t) dC_buf[C_pos[t]] -= d_latent[t];
                        }
                    }
                } else if constexpr (!C_is_const) {
                    std::fill_n(dC_buf, 2 * n_C, 0.0);
                }

                if (n_C > 0 && R_nonempty > 0) {
                    if constexpr (!C_is_const) {
                        for (int c = 0; c < 2; ++c) {
                            const double* xc = x + c * n_C;
                            const double* lFc = lF + c * n_C;
                            const double* lSc = lS + c * n_C;
                            const double* lp_c = log_p + c * n_total_segs;
                            double* ro = ratio_own + c * n_C;
                            double* rn = ratio_next + c * n_C;
                            for (int p = 0; p < n_C; ++p) {
                                const double log_f = (link == 1)
                                    ? (-0.5 * xc[p] * xc[p] - prior::half_log_of_two_pi)
                                    : (lFc[p] + lSc[p]);
                                ro[p] = log_f - lp_c[seg_own[p]];
                                rn[p] = log_f - lp_c[seg_next[p]];
                            }
                        }
                        k_exp(ratio_own, 4 * n_C);
                        if (link != 1) k_exp(lF, 4 * n_C);
                        for (int c = 0; c < 2; ++c) {
                            const double* lFc = lF + c * n_C;
                            const double* lSc = lS + c * n_C;
                            const double* xc = x + c * n_C;
                            const double* ro = ratio_own + c * n_C;
                            const double* rn = ratio_next + c * n_C;
                            const double* usa = u_seg_alpha + c * (n_u_segs + 1);
                            double* dCc = dC_buf + c * n_C;
                            for (int p = 0; p < n_C; ++p) {
                                const double a_own = usa[seg_uid[seg_own[p]]];
                                const double a_next = usa[seg_uid[seg_next[p]]];
                                double d_F = 0.0;
                                if (a_own != 1.0) d_F += (a_own - 1.0) * ro[p];
                                if (a_next != 1.0) d_F -= (a_next - 1.0) * rn[p];
                                const double d_J = (link == 1) ? (-xc[p]) : (lSc[p] - lFc[p]);
                                dCc[p] += adj * (d_F + d_J);
                            }
                        }
                        for (int p = 0; p < n_C; ++p) C_vars_0.coeffRef(p).adj() += dC_buf[p];
                        for (int p = 0; p < n_C; ++p) C_vars_1.coeffRef(p).adj() += dC_buf[n_C + p];
                    }
                    if constexpr (!alpha_is_const) {
                        double diff_stack[256];
                        std::vector<double> diff_heap;
                        double* diff = (n_cat + 1 <= 256) ? diff_stack : (diff_heap.resize(n_cat + 1), diff_heap.data());
                        for (int c = 0; c < 2; ++c) {
                            const double* lp_c = log_p + c * n_total_segs;
                            const double* udg = u_digamma + c * (n_u_segs + 1);
                            const double dg_sum = udg[n_u_segs];
                            std::fill_n(u_lp_sum, n_u_segs, 0.0);
                            for (int s_i = 0; s_i < n_total_segs; ++s_i) u_lp_sum[seg_uid[s_i]] += lp_c[s_i];
                            std::fill_n(diff, n_cat + 1, 0.0);
                            for (int u = 0; u < n_u_segs; ++u) {
                                const double w = u_lp_sum[u] + static_cast<double>(u_seg_cnt[u]) * (dg_sum - udg[u]);
                                diff[seg_lo[u]] += w;
                                diff[seg_hi[u]] -= w;
                            }
                            auto& av = (c == 0) ? alpha_vars_0 : alpha_vars_1;
                            double running = 0.0;
                            for (int k = 0; k < n_cat; ++k) {
                                running += diff[k];
                                av.coeffRef(k).adj() += adj * running;
                            }
                        }
                    }
                } else if constexpr (!C_is_const) {
                    for (int p = 0; p < n_C; ++p) C_vars_0.coeffRef(p).adj() += dC_buf[p];
                    for (int p = 0; p < n_C; ++p) C_vars_1.coeffRef(p).adj() += dC_buf[n_C + p];
                }
            });
    }
}

}  // namespace fused_stan_math

using fused_stan_math::fused_ord_biv_randthr_ind_dir_and_binomial_precomputed_lpmf;
