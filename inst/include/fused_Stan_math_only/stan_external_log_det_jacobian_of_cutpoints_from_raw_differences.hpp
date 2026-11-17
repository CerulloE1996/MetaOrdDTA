#pragma once
// ---- Stan math only: a copy of NicoSIMD/stan_external_log_det_jacobian_of_cutpoints_from_raw_differences.hpp (md5 d19a410e199ab26d7e974e171e426ffd) with its scalar code and
//      its hand-written gradients, the SIMD branches left out, and log, exp, Phi, ... of doubles from
//      std:: and stan::math:: (the buffer helpers k_exp, k_log, ... are in
//      stan_external_functions_on_Stan_math_only.hpp). Where the notes below speak of SIMD or NicoSIMD
//      kernels, the same values come from Stan math here: no -m flags, any CPU.
//
// -| --------- NicoSIMD log determinant of the Jacobian of ordered cutpoints built from raw differences ---------
//
//     nico_log_det_jac_of_cuts_from_raw_diffs(C_raw, softplus)
//
// For a set of K ordered cutpoints c built from K unconstrained values r (a vector or row_vector C_raw, or each
// row of a matrix C_raw: one set per group, e.g. the random cutpoints of each study) as
//
//     c_1 = r_1,   c_k = c_(k-1) + exp(r_k)         (softplus = 0),
//     c_1 = r_1,   c_k = c_(k-1) + log1p_exp(r_k)   (softplus = 1),          k = 2, ..., K,
//
// the log determinant of the Jacobian of the map from r to c (lower triangular, so the sum of the logs of its
// diagonal: d c_k / d r_k), summed over every set:
//
//     softplus = 0:   sum over the sets of sum_(k = 2..K) r_k,
//     softplus = 1:   sum over the sets of sum_(k = 2..K) log_inv_logit(r_k).
//
// A general function: any ordinal model whose ordered cutpoints are the first one plus the exp (or softplus) of
// unconstrained raw differences (ordinal regression, ordinal meta-analysis and network meta-analysis of test
// accuracy, ordinal latent class models), with fixed cutpoints (one set) or random cutpoints (one set per group).
// The same value and derivatives as MetaOrdDTA's raw_C_to_C_log_det_J_lp (Stan_fns_Jacobian.stan, with the
// Jacobians of construct_C(): softplus -> log_inv_logit(r), exp -> r),
//
//     if (softplus == 1)  log_det_J = sum(log_inv_logit(C_raw[2:K]));     (C_raw[, 2:K] for a matrix)
//     else                log_det_J = sum(C_raw[2:K]);
//
// in ONE autodiff node for all the sets, in place of one node per element plus a sum per set: the sum with
// Eigen's SIMD reduction (AVX2 / AVX-512 packets), and in the softplus case BayesMVP's SIMD kernels for
// log_inv_logit (the value, the same kernel as nico_log_inv_logit) and inv_logit (the derivative,
// d / d r log_inv_logit(r) = inv_logit(-r), the same kernel as nico_inv_logit), computed in the forward pass and
// kept on the arena, so that the reverse pass is the adjoint update alone. d / d r_1 = 0 (the first cutpoint is
// r_1 itself); d / d r_k = 1 (softplus = 0) or inv_logit(-r_k) (softplus = 1), k = 2, ..., K.
//
// Declared in the Stan functions block (no body; stanc --allow-undefined) as any of
//
//     real nico_log_det_jac_of_cuts_from_raw_diffs(vector C_raw, data int softplus);
//     real nico_log_det_jac_of_cuts_from_raw_diffs(row_vector C_raw, data int softplus);
//     real nico_log_det_jac_of_cuts_from_raw_diffs(matrix C_raw, data int softplus);
//
// and used as   target += nico_log_det_jac_of_cuts_from_raw_diffs(C_raw, softplus);
// The same function is also given under the name of MetaOrdDTA's Stan function with nico_ in front,
//
//     real nico_raw_C_to_C_log_det_J_lp(vector C_raw, data int softplus);   (also row_vector and matrix)
//
// (an _lp function, so stanc passes it lp__ and lp_accum__, which it does not use).
// It uses the namespace fused_stan_math and the kernels of stan_external_functions.hpp, which it includes
// (#pragma once).
//
// -----------------------------------------------------------------------------------------------------------

#include "stan_external_functions_on_Stan_math_only.hpp"

#include <stan/math/rev.hpp>
#include <stan/math/prim.hpp>

#include <Eigen/Dense>
#include <cstddef>
#include <ostream>
#include <type_traits>
#include <vector>

// ---- the Stan-facing function (global scope; stanc calls it as
//      nico_log_det_jac_of_cuts_from_raw_diffs(C_raw, softplus, pstream__)): C_raw a vector or
//      row_vector (one set of cutpoints) or a matrix (one set per row)
namespace fused_stan_math {

template <typename T_C_raw, stan::require_eigen_t<T_C_raw>* = nullptr>
inline stan::return_type_t<T_C_raw> nico_log_det_jac_of_cuts_from_raw_diffs(
        const T_C_raw& C_raw, const int softplus, std::ostream* /*pstream__*/) {

    using stan::math::var;
    constexpr bool C_raw_is_one_set  = stan::is_eigen_vector<T_C_raw>::value;
    constexpr bool C_raw_is_constant = stan::is_constant_all<T_C_raw>::value;
    const int R = C_raw_is_one_set ? 1 : static_cast<int>(C_raw.rows());
    const int K = C_raw_is_one_set ? static_cast<int>(C_raw.size()) : static_cast<int>(C_raw.cols());
    if ((K < 2) || (R == 0)) return 0.0;
    const auto& C_raw_ref = stan::math::to_ref(C_raw);
    //// (the K - 1 raw differences of every set, set after set: r_2, ..., r_K of set 1, then of set 2, ...)
    const int n_differences_per_set = K - 1;
    const std::size_t n_differences = static_cast<std::size_t>(R) * n_differences_per_set;
    //// (the buffers: the raw differences, then their derivatives in the softplus case; on the arena when the
    ////  derivatives are kept for the reverse pass, in a std::vector otherwise, so that no evaluation of the log
    ////  density alone, and no exp case, leaves memory on the arena)
    const bool keep_the_derivatives = !C_raw_is_constant && (softplus == 1);
    std::vector<double> buffer_without_autodiff;
    double* r = nullptr;
    if (keep_the_derivatives) {
        r = stan::math::ChainableStack::instance_->memalloc_.alloc_array<double>(2 * n_differences);
    } else {
        buffer_without_autodiff.resize(n_differences);
        r = buffer_without_autodiff.data();
    }
    for (int s = 0; s < R; ++s) {
        double* r_s = r + static_cast<std::size_t>(s) * n_differences_per_set;
        for (int k = 1; k < K; ++k) {
            r_s[k - 1] = C_raw_is_one_set ? stan::math::value_of(C_raw_ref.coeff(k))
                                          : stan::math::value_of(C_raw_ref.coeff(s, k));
        }
    }
    double log_det_J = 0.0;
    double* d_r = keep_the_derivatives ? r + n_differences : nullptr;
    if (softplus == 1) {
        //// (log_inv_logit(r) and, for the reverse pass, inv_logit(-r), with BayesMVP's SIMD kernels)
        std::vector<double> log_inv_logit_of_r(r, r + n_differences);
        fused_stan_math::k_log_inv_logit(log_inv_logit_of_r.data(), static_cast<int>(n_differences));
        //// (added in index order, so that the sum does not depend on where the buffer sits in memory)
        for (std::size_t i = 0; i < n_differences; ++i) log_det_J += log_inv_logit_of_r[i];
        if (keep_the_derivatives) {
            for (std::size_t i = 0; i < n_differences; ++i) d_r[i] = -r[i];
            fused_stan_math::k_inv_logit(d_r, static_cast<int>(n_differences));
        }
    } else {
        //// (added in index order, so that the sum does not depend on where the buffer sits in memory)
        for (std::size_t i = 0; i < n_differences; ++i) log_det_J += r[i];
    }
    if constexpr (C_raw_is_constant) {
        return log_det_J;
    } else {
        using C_raw_vars_type = stan::plain_type_t<stan::math::promote_scalar_t<var, T_C_raw>>;
        stan::arena_t<C_raw_vars_type> C_raw_vars = C_raw_ref;
        const bool softplus_case = (softplus == 1);
        return stan::math::make_callback_var(log_det_J, [C_raw_vars, d_r, softplus_case, R, K,
                                                         n_differences_per_set](auto& vi) mutable {
            const double adjoint = vi.adj();
            for (int s = 0; s < R; ++s) {
                const double* d_r_s = softplus_case ? d_r + static_cast<std::size_t>(s) * n_differences_per_set
                                                    : nullptr;
                for (int k = 1; k < K; ++k) {
                    const double d = softplus_case ? adjoint * d_r_s[k - 1] : adjoint;
                    if constexpr (C_raw_is_one_set) {
                        C_raw_vars.coeffRef(k).adj() += d;
                    } else {
                        C_raw_vars.coeffRef(s, k).adj() += d;
                    }
                }
            }
        });
    }

}

}  // namespace fused_stan_math

using fused_stan_math::nico_log_det_jac_of_cuts_from_raw_diffs;
//
// ---- the same function under the name the MetaOrdDTA NMA Stan files declare:
//          real nico_log_det_jacobian_of_cutpoints_from_raw_differences(vector C_raw, data int softplus);
//          real nico_log_det_jacobian_of_cutpoints_from_raw_differences(matrix C_raw, data int softplus);
namespace fused_stan_math {

template <typename T_C_raw, stan::require_eigen_t<T_C_raw>* = nullptr>
inline stan::return_type_t<T_C_raw> nico_log_det_jacobian_of_cutpoints_from_raw_differences(
        const T_C_raw& C_raw, const int softplus, std::ostream* pstream__) {

    return nico_log_det_jac_of_cuts_from_raw_diffs(C_raw, softplus, pstream__);

}

}  // namespace fused_stan_math

using fused_stan_math::nico_log_det_jacobian_of_cutpoints_from_raw_differences;
//
// ---- the same function under the name of MetaOrdDTA's Stan function with nico_ in front (an _lp function: stanc
//      calls it as nico_raw_C_to_C_log_det_J_lp<propto__>(C_raw, softplus, lp__, lp_accum__, pstream__))
namespace fused_stan_math {

template <bool propto__, typename T_C_raw, typename T_lp, typename T_lp_accum,
          stan::require_eigen_t<T_C_raw>* = nullptr>
inline stan::return_type_t<T_C_raw> nico_raw_C_to_C_log_det_J_lp(const T_C_raw& C_raw, const int softplus,
                                                                 T_lp& /*lp__*/, T_lp_accum& /*lp_accum__*/,
                                                                 std::ostream* pstream__) {

    return nico_log_det_jac_of_cuts_from_raw_diffs(C_raw, softplus, pstream__);

}

}  // namespace fused_stan_math

using fused_stan_math::nico_raw_C_to_C_log_det_J_lp;
























// ===========================================================================================================
// ---- Multi-width NicoSIMD::k_* kernels for ordered fixed-cutpoint bivariate model log-posterior and cutpoints ----
//      (supports AVX-512, AVX2, ARM64 NEON via SIMDe, and Scalar)
// ===========================================================================================================

#include "stan_external_binomial_counts_above_successive_cutpoints_of_a_cumulative_link.hpp"
#include "stan_external_induced_dirichlet_prior_of_ordered_cutpoints.hpp"
#include "stan_external_induced_dirichlet_with_Se_Sp_prior.hpp"
#include "stan_external_lkj_corr_cholesky_lpdf.hpp"

namespace fused_stan_math {
namespace ord_biv_fixedthr_detail {

static constexpr double HALF_LOG_TWO_PI = 0.91893853320467274178;

template <bool include_const>
inline double normal_term_val_and_adj(const double y, const double mu, const double sd, double& d_y) {
    if (sd <= 0.0 || !std::isfinite(sd)) {
        stan::math::check_positive_finite("nico_ord_biv_fixedthr_model_lp", "Scale parameter", sd);
    }
    const double inv_sd = 1.0 / sd;
    const double z = (y - mu) * inv_sd;
    d_y = -z * inv_sd;
    double lp = -0.5 * z * z;
    if constexpr (include_const) lp -= std::log(sd) + HALF_LOG_TWO_PI;
    return lp;
}

template <bool include_const>
inline double lkj_2x2_val_and_adj(const double L00, const double L01, const double L11,
                                  const double eta, double& d_L11) {
    stan::math::check_positive("nico_ord_biv_fixedthr_model_lp", "Shape parameter", eta);
    if (!(L00 > 0.0) || !(L11 > 0.0) || (L01 != 0.0)) {
        Eigen::Matrix<double, 2, 2> bad_L;
        bad_L << L00, L01, 0.0, L11;
        stan::math::check_cholesky_factor("nico_ord_biv_fixedthr_model_lp", "Random variable", bad_L);
    }
    double lp = 0.0;
    if constexpr (include_const) lp += nico_lkj_constant_from_data(eta, 2U, nullptr);
    d_L11 = 0.0;
    if (eta != 1.0) {
        const double coef = 2.0 * eta - 2.0;
        lp += coef * std::log(L11);
        d_L11 = coef / L11;
    }
    return lp;
}

template <typename T_mean, typename T_sd>
inline void check_and_add_normal_block(const char* name, const double* vals, int n,
                                       const T_mean& mean_mat, const T_sd& sd_mat, int row,
                                       bool eval_consts, double& total, double* d_out) {
    for (int j = 0; j < n; ++j) {
        const double mu = stan::math::value_of(mean_mat.coeff(row, j));
        const double sd = stan::math::value_of(sd_mat.coeff(row, j));
        double d_j = 0.0;
        total += eval_consts ? normal_term_val_and_adj<true>(vals[j], mu, sd, d_j)
                             : normal_term_val_and_adj<false>(vals[j], mu, sd, d_j);
        if (d_out) d_out[j] += d_j;
        (void)name;
    }
}

struct FCClassBuffers {
    int K = 0;
    int m = 0;
    double* x_shift = nullptr;      // K
    double* alpha = nullptr;        // K + 1
    double* F = nullptr;            // K
    double* p = nullptr;            // K + 1
    double* log_p = nullptr;        // 2 * K + 1
    double* work_id = nullptr;      // K (forward J_terms / reverse d_x)
    int* chosen_idx = nullptr;      // m (0-based)
    double* d_x_chosen = nullptr;   // m
    bool add_se_sp_to_C = false;
    bool add_se_sp_to_beta_mu = false;
    double sum_d_x_chosen = 0.0;
};

template <bool include_consts>
inline void eval_class_cutpoint_prior(
        const int c_idx,
        const double* C_val,
        const double baseline_loc_prior,
        const Eigen::VectorXd& alpha_vec,
        const int link,
        const int use_Se_Sp_prior,
        const int Se_Sp_prior_type,
        const int m_c,
        const std::vector<int>& chosen_c,
        const Eigen::VectorXd& alpha_seg_c,
        const Eigen::VectorXd& dir_A_c,
        const Eigen::MatrixXd& beta_td_c,
        const Eigen::MatrixXd& link_ms_td_c,
        const double* beta_mu_row,
        const double* baseline_case,
        int n_cov,
        bool is_var_C,
        bool is_var_beta_mu,
        double* scratch_buf,
        FCClassBuffers& cb,
        double& total) {
    namespace prior = fused_stan_math::induced_dirichlet_prior;
    namespace Se_Sp = fused_stan_math::Se_Sp_prior_at_chosen_cutoffs;
    const int K = static_cast<int>(alpha_vec.size()) - 1;
    cb.K = K;
    cb.m = m_c;
    cb.add_se_sp_to_C = false;
    cb.add_se_sp_to_beta_mu = false;
    cb.sum_d_x_chosen = 0.0;

    for (int k = 0; k < K; ++k) {
        cb.x_shift[k] = C_val[k] - baseline_loc_prior;
    }
    for (int j = 0; j <= K; ++j) cb.alpha[j] = alpha_vec.coeff(j);

    double id_sum = 0.0;
    if (!prior::forward_pass<include_consts>(cb.x_shift, cb.alpha, link, K,
                                             cb.F, cb.p, cb.log_p, cb.work_id, id_sum)) {
        Eigen::Matrix<double, -1, 1> p_set(K + 1);
        for (int j = 0; j <= K; ++j) p_set(j) = std::exp(cb.log_p[j]);
        prior::run_the_checks_of_Stan(p_set, alpha_vec);
    }
    total += id_sum;

    if constexpr (include_consts) {
        double sum_a = 0.0, sum_lg = 0.0;
        for (int j = 0; j <= K; ++j) {
            sum_a += cb.alpha[j];
            sum_lg += std::lgamma(cb.alpha[j]);
        }
        total += std::lgamma(sum_a) - sum_lg;
    }

    if (use_Se_Sp_prior >= 1 && m_c > 0) {
        for (int j = 0; j < m_c; ++j) {
            const int p0 = chosen_c[j] - 1;
            if (p0 < 0 || p0 >= K) {
                throw std::domain_error("nico_ord_biv_fixedthr_model_lp: chosen_cutoffs out of bounds");
            }
            cb.chosen_idx[j] = p0;
            cb.d_x_chosen[j] = 0.0;
        }
        Se_Sp::buffers_of_the_Se_Sp_prior sp_buf;
        double* x_ch = scratch_buf;
        sp_buf.place(x_ch + m_c, m_c);

        if (use_Se_Sp_prior == 2) {
            for (int j = 0; j < m_c; ++j) x_ch[j] = cb.x_shift[cb.chosen_idx[j]];
            Se_Sp::values_of_the_hyperparameters(
                "nico_ord_biv_fixedthr_model_lp",
                alpha_seg_c.head(m_c + 1), Se_Sp_prior_type,
                dir_A_c.head(m_c + 1), beta_td_c.topRows(m_c), link_ms_td_c.topRows(m_c),
                sp_buf);
            if (!Se_Sp::values_of_one_set(x_ch, link, sp_buf)) {
                throw std::domain_error("nico_ord_biv_fixedthr_model_lp: cutpoints at chosen cutoffs not increasing");
            }
            total += Se_Sp::dirichlet_of_the_segments(sp_buf.segments, sp_buf.log_q, sp_buf.log_f, m_c,
                                                      include_consts, -1.0, cb.d_x_chosen, nullptr);
            total += Se_Sp::Se_Sp_prior(x_ch, link, Se_Sp_prior_type, sp_buf,
                                        include_consts, false, cb.d_x_chosen);
            cb.add_se_sp_to_C = is_var_C;
        } else {
            // use_Se_Sp_prior == 1: log_f + Se_Sp_prior at x_chosen = C[chosen] - Xbeta_baseline
            double xbeta_base = 0.0;
            for (int i = 0; i < n_cov; ++i) xbeta_base += baseline_case[i] * beta_mu_row[i];
            for (int j = 0; j < m_c; ++j) x_ch[j] = C_val[cb.chosen_idx[j]] - xbeta_base;
            const int type_of_H = (Se_Sp_prior_type == 1 || Se_Sp_prior_type == 2) ? Se_Sp_prior_type : 4;
            Eigen::VectorXd ones_unused = Eigen::VectorXd::Ones(m_c + 1);
            Se_Sp::values_of_the_hyperparameters(
                "nico_ord_biv_fixedthr_model_lp",
                ones_unused, type_of_H,
                dir_A_c.head(m_c + 1), beta_td_c.topRows(m_c), link_ms_td_c.topRows(m_c),
                sp_buf);
            if (type_of_H == 1) {
                if (!Se_Sp::values_of_one_set(x_ch, link, sp_buf)) {
                    throw std::domain_error("nico_ord_biv_fixedthr_model_lp: x_chosen is not increasing");
                }
            } else {
                Se_Sp::values_of_one_set_in_any_order(x_ch, link, sp_buf);
            }
            for (int j = 0; j < m_c; ++j) {
                total += sp_buf.log_f[j];
                cb.d_x_chosen[j] += (link == 1) ? -x_ch[j] : (std::exp(sp_buf.lS[j]) - std::exp(sp_buf.lF[j]));
            }
            total += Se_Sp::Se_Sp_prior(x_ch, link, type_of_H, sp_buf,
                                        include_consts, false, cb.d_x_chosen);
            cb.add_se_sp_to_C = is_var_C;
            cb.add_se_sp_to_beta_mu = is_var_beta_mu;
            if (is_var_beta_mu) {
                double s = 0.0;
                for (int j = 0; j < m_c; ++j) s += cb.d_x_chosen[j];
                cb.sum_d_x_chosen = s;
            }
        }
    }
    (void)c_idx;
}

}  // namespace ord_biv_fixedthr_detail

// ---- 1) Fast ordered cutpoint construction C[k] from raw differences C_raw with one reverse_pass_callback:
//         C[0] = C_raw[0]; C[k] = C[k - 1] + (softplus == 1 ? log1p_exp(C_raw[k]) : exp(C_raw[k]))
template <typename T_C_raw, stan::require_eigen_col_vector_t<T_C_raw>* = nullptr>
inline Eigen::Matrix<stan::return_type_t<T_C_raw>, -1, 1>
nico_construct_C(const T_C_raw& C_raw, const int softplus, std::ostream* /*pstream__*/) {
    using stan::math::var;
    const int K = static_cast<int>(C_raw.size());
    const auto& C_raw_ref = stan::math::to_ref(C_raw);
    if constexpr (stan::is_constant_all<T_C_raw>::value) {
        Eigen::Matrix<double, -1, 1> C_out(K);
        if (K == 0) return C_out;
        C_out.coeffRef(0) = C_raw_ref.coeff(0);
        if (K > 1) {
            std::vector<double> inc(K - 1);
            for (int j = 0; j < K - 1; ++j) inc[j] = C_raw_ref.coeff(j + 1);
            if (softplus == 1) k_log1p_exp(inc.data(), K - 1);
            else               k_exp(inc.data(), K - 1);
            for (int k = 1; k < K; ++k) C_out.coeffRef(k) = C_out.coeffRef(k - 1) + inc[k - 1];
        }
        return C_out;
    } else {
        Eigen::Matrix<var, -1, 1> C_out(K);
        if (K == 0) return C_out;
        if (K == 1) {
            C_out.coeffRef(0) = C_raw_ref.coeff(0);
            return C_out;
        }
        auto& arena = stan::math::ChainableStack::instance_->memalloc_;
        double* deriv = arena.alloc_array<double>(K - 1);
        std::vector<double> inc(K - 1);
        for (int j = 0; j < K - 1; ++j) {
            const double rk = stan::math::value_of(C_raw_ref.coeff(j + 1));
            inc[j] = rk;
            if (softplus == 1) deriv[j] = rk;
        }
        if (softplus == 1) {
            k_log1p_exp(inc.data(), K - 1);
            k_inv_logit(deriv, K - 1);
        } else {
            k_exp(inc.data(), K - 1);
            for (int j = 0; j < K - 1; ++j) deriv[j] = inc[j];
        }
        double c_val = stan::math::value_of(C_raw_ref.coeff(0));
        C_out.coeffRef(0) = var(new stan::math::vari(c_val, false));
        for (int k = 1; k < K; ++k) {
            c_val += inc[k - 1];
            C_out.coeffRef(k) = var(new stan::math::vari(c_val, false));
        }
        stan::arena_t<Eigen::Matrix<var, -1, 1>> C_raw_vars = C_raw_ref;
        stan::arena_t<Eigen::Matrix<var, -1, 1>> C_out_vars = C_out;
        stan::math::reverse_pass_callback([C_raw_vars, C_out_vars, deriv, K]() mutable {
            double suffix_adj = 0.0;
            for (int k = K - 1; k >= 1; --k) {
                suffix_adj += C_out_vars.coeff(k).adj();
                C_raw_vars.coeffRef(k).adj() += suffix_adj * deriv[k - 1];
            }
            suffix_adj += C_out_vars.coeff(0).adj();
            C_raw_vars.coeffRef(0).adj() += suffix_adj;
        });
        return C_out;
    }
}

// ---- 2) Single fused O_biv_FC model_lp kernel:
//         Fuses beta_mu normal + beta_SD normal + 2x2 LKJ-Cholesky + beta_z std_normal +
//         raw_C_to_C log-det-J (both classes) + Induced-Dirichlet (+ Se/Sp prior for use_Se_Sp_prior 1/2) +
//         flattened Binomial cumulative-link likelihood into ONE autodiff callback node.
template <bool propto__,
          typename T_bmu, typename T_bsd, typename T_L, typename T_bz,
          typename T_Craw1, typename T_Craw2, typename T_C1, typename T_C2, typename T_Xbeta,
          typename T_lp, typename T_lp_accum>
inline stan::return_type_t<T_bmu, T_bsd, T_L, T_bz, T_Craw1, T_Craw2, T_C1, T_C2, T_Xbeta>
nico_ord_biv_fixedthr_model_lp(
        const T_bmu& beta_mu,
        const T_bsd& beta_SD,
        const T_L& beta_L_Omega,
        const T_bz& beta_z,
        const T_Craw1& C_raw_1,
        const T_Craw2& C_raw_2,
        const T_C1& C_1,
        const T_C2& C_2,
        const T_Xbeta& Xbeta,
        const int n_covariates_nd,
        const int n_covariates_d,
        const Eigen::MatrixXd& prior_beta_mu_mean,
        const Eigen::MatrixXd& prior_beta_mu_SD,
        const Eigen::VectorXd& prior_beta_SD_mean,
        const Eigen::VectorXd& prior_beta_SD_SD,
        const double prior_beta_corr_LKJ,
        const int softplus,
        const int use_probit_link,
        const std::vector<Eigen::VectorXd>& prior_dirichlet_alpha,
        const Eigen::VectorXd& baseline_location_prior,
        const Eigen::VectorXd& baseline_case_nd,
        const Eigen::VectorXd& baseline_case_d,
        const int use_Se_Sp_prior,
        const int Se_Sp_prior_type,
        const std::vector<int>& n_chosen_cutoffs,
        const std::vector<std::vector<int>>& chosen_cutoffs,
        const std::vector<Eigen::VectorXd>& alpha_of_the_segments,
        const std::vector<Eigen::VectorXd>& prior_Se_Sp_dirichlet_A,
        const std::vector<Eigen::MatrixXd>& prior_beta_td,
        const std::vector<Eigen::MatrixXd>& prior_link_mean_SD_td,
        const std::vector<int>& term_n,
        const std::vector<int>& term_n_prev,
        const std::vector<int>& term_s,
        const std::vector<int>& term_c,
        const std::vector<int>& term_k,
        const Eigen::VectorXd& term_is_first,
        const double bound,
        T_lp& /*lp__*/,
        T_lp_accum& /*lp_accum__*/,
        std::ostream* /*pstream__*/) {

    namespace det = fused_stan_math::ord_biv_fixedthr_detail;
    namespace counts = fused_stan_math::binomial_counts_above_successive_cutpoints;
    namespace prior = fused_stan_math::induced_dirichlet_prior;
    namespace Se_Sp = fused_stan_math::Se_Sp_prior_at_chosen_cutoffs;
    using stan::math::var;

    constexpr bool bmu_var   = !stan::is_constant_all<T_bmu>::value;
    constexpr bool bsd_var   = !stan::is_constant_all<T_bsd>::value;
    constexpr bool L_var     = !stan::is_constant_all<T_L>::value;
    constexpr bool bz_var    = !stan::is_constant_all<T_bz>::value;
    constexpr bool craw1_var = !stan::is_constant_all<T_Craw1>::value;
    constexpr bool craw2_var = !stan::is_constant_all<T_Craw2>::value;
    constexpr bool C1_var    = !stan::is_constant_all<T_C1>::value;
    constexpr bool C2_var    = !stan::is_constant_all<T_C2>::value;
    constexpr bool xb_var    = !stan::is_constant_all<T_Xbeta>::value;
    constexpr bool any_var   = bmu_var || bsd_var || L_var || bz_var || craw1_var || craw2_var || C1_var || C2_var || xb_var;

    const auto& bmu_ref   = stan::math::to_ref(beta_mu);
    const auto& bsd_ref   = stan::math::to_ref(beta_SD);
    const auto& L_ref     = stan::math::to_ref(beta_L_Omega);
    const auto& bz_ref    = stan::math::to_ref(beta_z);
    const auto& craw1_ref = stan::math::to_ref(C_raw_1);
    const auto& craw2_ref = stan::math::to_ref(C_raw_2);
    const auto& C1_ref    = stan::math::to_ref(C_1);
    const auto& C2_ref    = stan::math::to_ref(C_2);
    const auto& xb_ref    = stan::math::to_ref(Xbeta);

    const int n_thr     = static_cast<int>(C1_ref.size());
    const int n_studies = static_cast<int>(bz_ref.rows());
    const int L_terms   = static_cast<int>(term_n.size());
    const int n_cov[2]  = {n_covariates_nd, n_covariates_d};
    const int max_cov   = std::max(n_covariates_nd, n_covariates_d);
    const int n_diffs   = (n_thr > 1) ? (2 * (n_thr - 1)) : 0;
    const int m0        = n_chosen_cutoffs[0];
    const int m1        = n_chosen_cutoffs[1];
    const int m_max     = std::max(m0, m1);

    const std::size_t class_buf_sz = static_cast<std::size_t>(4 * n_thr + 2 + (2 * n_thr + 1) + 2 * m_max);
    const std::size_t sp_scratch_sz = static_cast<std::size_t>(m_max + Se_Sp::buffers_of_the_Se_Sp_prior::n_values(m_max));
    const std::size_t n_doubles =
        static_cast<std::size_t>(n_cov[0] + n_cov[1])  // d_bmu
        + 2U                                           // d_bsd
        + static_cast<std::size_t>(2 * n_studies)      // bz_val
        + static_cast<std::size_t>(n_diffs)            // d_craw
        + 2U * static_cast<std::size_t>(n_thr)         // C_val
        + 2U * class_buf_sz                            // Induced-Dirichlet + Se/Sp per class
        + static_cast<std::size_t>(2 * n_thr)          // id reverse ratios
        + sp_scratch_sz                                // scratch for Se/Sp eval
        + static_cast<std::size_t>(7 * L_terms)        // binomial forward/reverse
        + static_cast<std::size_t>(max_cov + n_diffs); // temp bmu row vals & diff vals

    std::vector<double> host_d;
    int* chosen_mem = nullptr;
    std::vector<int> host_i;
    double* buf = nullptr;
    if constexpr (!any_var) {
        host_d.resize(n_doubles);
        host_i.resize(static_cast<std::size_t>(m0 + m1));
        buf = host_d.data();
        chosen_mem = host_i.data();
    } else {
        auto& arena = stan::math::ChainableStack::instance_->memalloc_;
        buf = arena.alloc_array<double>(n_doubles);
        chosen_mem = arena.alloc_array<int>(m0 + m1);
    }

    double* p = buf;
    double* d_bmu0 = p; p += n_cov[0];
    double* d_bmu1 = p; p += n_cov[1];
    double* d_bsd  = p; p += 2;
    double* bz_val = p; p += 2 * n_studies;
    double* d_craw = p; p += n_diffs;
    double* C1_val = p; p += n_thr;
    double* C2_val = p; p += n_thr;

    det::FCClassBuffers cb[2];
    for (int c = 0; c < 2; ++c) {
        cb[c].x_shift    = p; p += n_thr;
        cb[c].alpha      = p; p += n_thr + 1;
        cb[c].F          = p; p += n_thr;
        cb[c].p          = p; p += n_thr + 1;
        cb[c].log_p      = p; p += 2 * n_thr + 1;
        cb[c].work_id    = p; p += n_thr;
        cb[c].chosen_idx = (c == 0) ? chosen_mem : (chosen_mem + m0);
        cb[c].d_x_chosen = p; p += m_max;
    }
    double* id_ratios   = p; p += 2 * n_thr;
    double* sp_scratch  = p; p += sp_scratch_sz;
    double* latent      = p; p += L_terms;
    double* binom_S     = p; p += L_terms;
    double* binom_p     = p; p += L_terms;
    double* binom_ppass = p; p += L_terms;
    double* binom_lp    = p; p += L_terms;
    double* binom_l1mp  = p; p += L_terms;
    double* bmu_row_tmp = p; p += max_cov;
    double* diff_tmp    = p;

    for (int j = 0; j < n_cov[0]; ++j) d_bmu0[j] = 0.0;
    for (int j = 0; j < n_cov[1]; ++j) d_bmu1[j] = 0.0;

    double total = 0.0;

    // 1. beta_mu priors
    for (int c = 0; c < 2; ++c) {
        double* d_bmu_c = (c == 0) ? d_bmu0 : d_bmu1;
        for (int j = 0; j < n_cov[c]; ++j) {
            const double val = stan::math::value_of(bmu_ref.coeff(c, j));
            const double mu  = prior_beta_mu_mean.coeff(c, j);
            const double sd  = prior_beta_mu_SD.coeff(c, j);
            double d_j = 0.0;
            total += (!propto__) ? det::normal_term_val_and_adj<true>(val, mu, sd, d_j)
                                 : det::normal_term_val_and_adj<false>(val, mu, sd, d_j);
            d_bmu_c[j] = d_j;
        }
    }

    // 2. beta_SD priors
    for (int c = 0; c < 2; ++c) {
        const double val = stan::math::value_of(bsd_ref.coeff(c));
        const double mu  = prior_beta_SD_mean.coeff(c);
        const double sd  = prior_beta_SD_SD.coeff(c);
        double d_j = 0.0;
        total += (!propto__) ? det::normal_term_val_and_adj<true>(val, mu, sd, d_j)
                             : det::normal_term_val_and_adj<false>(val, mu, sd, d_j);
        d_bsd[c] = d_j;
    }

    // 3. 2x2 LKJ-corr-Cholesky prior
    double d_L11 = 0.0;
    const double L00 = stan::math::value_of(L_ref.coeff(0, 0));
    const double L01 = stan::math::value_of(L_ref.coeff(0, 1));
    const double L11 = stan::math::value_of(L_ref.coeff(1, 1));
    total += (!propto__) ? det::lkj_2x2_val_and_adj<true>(L00, L01, L11, prior_beta_corr_LKJ, d_L11)
                         : det::lkj_2x2_val_and_adj<false>(L00, L01, L11, prior_beta_corr_LKJ, d_L11);

    // 4. beta_z standard Normal prior
    {
        double sum_sq = 0.0;
        for (int c = 0; c < 2; ++c) {
            for (int s = 0; s < n_studies; ++s) {
                const double zval = stan::math::value_of(bz_ref.coeff(s, c));
                bz_val[c * n_studies + s] = zval;
                sum_sq += zval * zval;
            }
        }
        total += -0.5 * sum_sq;
        if constexpr (!propto__) total -= static_cast<double>(2 * n_studies) * det::HALF_LOG_TWO_PI;
    }

    // 5. Cutpoint raw-differences log-det Jacobian (both classes in one SIMD batch)
    if (n_diffs > 0) {
        const int half = n_thr - 1;
        for (int k = 1; k < n_thr; ++k) {
            diff_tmp[k - 1]        = stan::math::value_of(craw1_ref.coeff(k));
            diff_tmp[half + k - 1] = stan::math::value_of(craw2_ref.coeff(k));
        }
        if (softplus == 1) {
            for (int i = 0; i < n_diffs; ++i) d_craw[i] = -diff_tmp[i];
            k_log_inv_logit(diff_tmp, n_diffs);
            for (int i = 0; i < n_diffs; ++i) total += diff_tmp[i];
            if constexpr (craw1_var || craw2_var) k_inv_logit(d_craw, n_diffs);
        } else {
            for (int i = 0; i < n_diffs; ++i) {
                total += diff_tmp[i];
                d_craw[i] = 1.0;
            }
        }
    }

    // 6. Induced-Dirichlet + optional Se/Sp priors
    for (int k = 0; k < n_thr; ++k) {
        C1_val[k] = stan::math::value_of(C1_ref.coeff(k));
        C2_val[k] = stan::math::value_of(C2_ref.coeff(k));
    }
    for (int c = 0; c < 2; ++c) {
        const double* C_val_c = (c == 0) ? C1_val : C2_val;
        const bool is_C_v = (c == 0) ? C1_var : C2_var;
        for (int j = 0; j < n_cov[c]; ++j) bmu_row_tmp[j] = stan::math::value_of(bmu_ref.coeff(c, j));
        const Eigen::VectorXd& base_vec = (c == 0) ? baseline_case_nd : baseline_case_d;
        det::eval_class_cutpoint_prior<!propto__>(
            c, C_val_c, baseline_location_prior.coeff(c), prior_dirichlet_alpha[c],
            use_probit_link, use_Se_Sp_prior, Se_Sp_prior_type,
            n_chosen_cutoffs[c], chosen_cutoffs[c], alpha_of_the_segments[c],
            prior_Se_Sp_dirichlet_A[c], prior_beta_td[c], prior_link_mean_SD_td[c],
            bmu_row_tmp, base_vec.data(), n_cov[c], is_C_v, bmu_var,
            sp_scratch, cb[c], total);
    }

    // 7. Flattened Binomial cumulative-link likelihood across both classes
    if (L_terms > 0) {
        Eigen::MatrixXd xb_val(n_studies, 2);
        for (int c = 0; c < 2; ++c) {
            for (int s = 0; s < n_studies; ++s) {
                xb_val.coeffRef(s, c) = stan::math::value_of(xb_ref.coeff(s, c));
            }
        }
        for (int t = 0; t < L_terms; ++t) {
            const int s0 = term_s[t] - 1;
            const int c0 = term_c[t] - 1;
            const int k0 = term_k[t] - 1;
            const double* C_val_c = (c0 == 0) ? C1_val : C2_val;
            latent[t] = xb_val.coeff(s0, c0) - C_val_c[k0];
        }
        double binom_sum = 0.0;
        if (!counts::forward_pass<!propto__>(
                term_n.data(), term_n_prev.data(), latent, term_is_first.data(),
                use_probit_link, bound, L_terms, binom_S, binom_p, binom_lp, binom_l1mp, binom_ppass, binom_sum)) {
            counts::run_the_checks_of_Stan(term_n, term_n_prev);
        }
        total += binom_sum;
    }

    if constexpr (!any_var) {
        return total;
    } else {
        stan::arena_t<Eigen::Matrix<var, -1, -1>> bmu_vars;
        stan::arena_t<Eigen::Matrix<var, 1, -1>>  bsd_vars;
        stan::arena_t<Eigen::Matrix<var, -1, -1>> L_vars;
        stan::arena_t<Eigen::Matrix<var, -1, -1>> bz_vars;
        stan::arena_t<Eigen::Matrix<var, -1, 1>>  craw1_vars;
        stan::arena_t<Eigen::Matrix<var, -1, 1>>  craw2_vars;
        stan::arena_t<Eigen::Matrix<var, -1, 1>>  C1_vars;
        stan::arena_t<Eigen::Matrix<var, -1, 1>>  C2_vars;
        stan::arena_t<Eigen::Matrix<var, -1, -1>> xb_vars;

        if constexpr (bmu_var)   bmu_vars   = bmu_ref;
        if constexpr (bsd_var)   bsd_vars   = bsd_ref;
        if constexpr (L_var)     L_vars     = L_ref;
        if constexpr (bz_var)    bz_vars    = bz_ref;
        if constexpr (craw1_var) craw1_vars = craw1_ref;
        if constexpr (craw2_var) craw2_vars = craw2_ref;
        if constexpr (C1_var)    C1_vars    = C1_ref;
        if constexpr (C2_var)    C2_vars    = C2_ref;
        if constexpr (xb_var)    xb_vars    = xb_ref;

        int* term_s_mem = nullptr;
        int* term_c_mem = nullptr;
        int* term_k_mem = nullptr;
        int* term_n_mem = nullptr;
        int* term_np_mem = nullptr;
        double* term_if_mem = nullptr;
        if (L_terms > 0) {
            auto& arena = stan::math::ChainableStack::instance_->memalloc_;
            term_s_mem  = arena.alloc_array<int>(L_terms);
            term_c_mem  = arena.alloc_array<int>(L_terms);
            term_k_mem  = arena.alloc_array<int>(L_terms);
            term_n_mem  = arena.alloc_array<int>(L_terms);
            term_np_mem = arena.alloc_array<int>(L_terms);
            term_if_mem = arena.alloc_array<double>(L_terms);
            for (int t = 0; t < L_terms; ++t) {
                term_s_mem[t]  = term_s[t] - 1;
                term_c_mem[t]  = term_c[t] - 1;
                term_k_mem[t]  = term_k[t] - 1;
                term_n_mem[t]  = term_n[t];
                term_np_mem[t] = term_n_prev[t];
                term_if_mem[t] = term_is_first.coeff(t);
            }
        }
        double* base_nd_mem = nullptr;
        double* base_d_mem  = nullptr;
        {
            auto& arena = stan::math::ChainableStack::instance_->memalloc_;
            base_nd_mem = arena.alloc_array<double>(n_covariates_nd);
            base_d_mem  = arena.alloc_array<double>(n_covariates_d);
            for (int i = 0; i < n_covariates_nd; ++i) base_nd_mem[i] = baseline_case_nd.coeff(i);
            for (int i = 0; i < n_covariates_d; ++i)  base_d_mem[i]  = baseline_case_d.coeff(i);
        }

        return stan::math::make_callback_var(
            total,
            [bmu_vars, bsd_vars, L_vars, bz_vars, craw1_vars, craw2_vars, C1_vars, C2_vars, xb_vars,
             n_cov_0 = n_cov[0], n_cov_1 = n_cov[1], d_bmu0, d_bmu1, d_bsd, d_L11, bz_val,
             n_studies, n_thr, n_diffs, d_craw, cb0 = cb[0], cb1 = cb[1], id_ratios,
             use_probit_link, base_nd_mem, base_d_mem,
             L_terms, term_n_mem, term_np_mem, term_s_mem, term_c_mem, term_k_mem, term_if_mem,
             latent, binom_S, binom_p, binom_ppass, binom_lp, binom_l1mp](auto& vi) mutable {
                const double adj = vi.adj();

                if constexpr (bmu_var) {
                    for (int j = 0; j < n_cov_0; ++j) bmu_vars.coeffRef(0, j).adj() += adj * d_bmu0[j];
                    for (int j = 0; j < n_cov_1; ++j) bmu_vars.coeffRef(1, j).adj() += adj * d_bmu1[j];
                    if (cb0.add_se_sp_to_beta_mu) {
                        const double sg = adj * cb0.sum_d_x_chosen;
                        for (int j = 0; j < n_cov_0; ++j) bmu_vars.coeffRef(0, j).adj() -= sg * base_nd_mem[j];
                    }
                    if (cb1.add_se_sp_to_beta_mu) {
                        const double sg = adj * cb1.sum_d_x_chosen;
                        for (int j = 0; j < n_cov_1; ++j) bmu_vars.coeffRef(1, j).adj() -= sg * base_d_mem[j];
                    }
                }
                if constexpr (bsd_var) {
                    bsd_vars.coeffRef(0).adj() += adj * d_bsd[0];
                    bsd_vars.coeffRef(1).adj() += adj * d_bsd[1];
                }
                if constexpr (L_var) {
                    if (d_L11 != 0.0) L_vars.coeffRef(1, 1).adj() += adj * d_L11;
                }
                if constexpr (bz_var) {
                    for (int c = 0; c < 2; ++c) {
                        const double* bzc = bz_val + c * n_studies;
                        for (int s = 0; s < n_studies; ++s) {
                            bz_vars.coeffRef(s, c).adj() -= adj * bzc[s];
                        }
                    }
                }
                if (n_diffs > 0) {
                    const int half = n_thr - 1;
                    if constexpr (craw1_var) {
                        for (int k = 1; k < n_thr; ++k) craw1_vars.coeffRef(k).adj() += adj * d_craw[k - 1];
                    }
                    if constexpr (craw2_var) {
                        for (int k = 1; k < n_thr; ++k) craw2_vars.coeffRef(k).adj() += adj * d_craw[half + k - 1];
                    }
                }
                if constexpr (C1_var) {
                    double* d_x = cb0.work_id;
                    prior::reverse_pass(cb0.x_shift, cb0.alpha, use_probit_link, n_thr,
                                        cb0.F, cb0.p, cb0.log_p, adj, d_x, id_ratios, id_ratios + n_thr);
                    if (cb0.add_se_sp_to_C) {
                        for (int j = 0; j < cb0.m; ++j) d_x[cb0.chosen_idx[j]] += adj * cb0.d_x_chosen[j];
                    }
                    for (int k = 0; k < n_thr; ++k) C1_vars.coeffRef(k).adj() += d_x[k];
                }
                if constexpr (C2_var) {
                    double* d_x = cb1.work_id;
                    prior::reverse_pass(cb1.x_shift, cb1.alpha, use_probit_link, n_thr,
                                        cb1.F, cb1.p, cb1.log_p, adj, d_x, id_ratios, id_ratios + n_thr);
                    if (cb1.add_se_sp_to_C) {
                        for (int j = 0; j < cb1.m; ++j) d_x[cb1.chosen_idx[j]] += adj * cb1.d_x_chosen[j];
                    }
                    for (int k = 0; k < n_thr; ++k) C2_vars.coeffRef(k).adj() += d_x[k];
                }
                if (L_terms > 0) {
                    double* d_S = binom_lp;
                    double* d_latent = binom_l1mp;
                    counts::reverse_pass(term_n_mem, term_np_mem, latent, term_if_mem,
                                         use_probit_link, L_terms, binom_S, binom_p, binom_ppass,
                                         adj, d_S, d_latent);
                    for (int t = 0; t < L_terms; ++t) {
                        const double dl = d_latent[t];
                        const int c0 = term_c_mem[t];
                        if constexpr (xb_var) xb_vars.coeffRef(term_s_mem[t], c0).adj() += dl;
                        if (c0 == 0) {
                            if constexpr (C1_var) C1_vars.coeffRef(term_k_mem[t]).adj() -= dl;
                        } else {
                            if constexpr (C2_var) C2_vars.coeffRef(term_k_mem[t]).adj() -= dl;
                        }
                    }
                }
            });
    }
}

}  // namespace fused_stan_math

using fused_stan_math::nico_construct_C;
using fused_stan_math::nico_ord_biv_fixedthr_model_lp;
