#pragma once
// ---- Stan math only: a copy of NicoSIMD/stan_external_set0_marginal_at_chosen_cutoffs.hpp (md5 ab51b130ff6c95011fa76ee307777023) with its scalar code and
//      its hand-written gradients, the SIMD branches left out, and log, exp, Phi, ... of doubles from
//      std:: and stan::math:: (the buffer helpers k_exp, k_log, ... are in
//      stan_external_functions_on_Stan_math_only.hpp). Where the notes below speak of SIMD or NicoSIMD
//      kernels, the same values come from Stan math here: no -m flags, any CPU.
// ---- NicoSIMD exact Set-0 marginal log-density at chosen cutoffs via 24-node Laplace-centered Gauss-Hermite ----
//
// Evaluates the exact normalized marginal log-density of x_S = C_S - Xbeta_baseline under Set 0:
//
//   log p_Set0(x_S) = log integral_{-inf}^{inf} N(b | 0, sqrt(v0))
//                       * Dirichlet(Delta F(x_S + b) | alpha_of_the_segments)
//                       * prod_{j=1}^m f(x_j + b) db
//
// where F = Phi (use_probit_link == 1) or inv_logit (use_probit_link == 0), f = F', and v0 > 0 is the
// prior variance of Xbeta_baseline under Set 0. Because g(b) is strictly concave in b, 8 Newton steps
// locate the exact mode b* and curvature scale sigma* = (-g''(b*))^{-1/2} to machine precision, and a
// 24-node probabilist Gauss-Hermite quadrature centered at (b*, sigma*) matches adaptive integration to ~1e-15
// across the entire parameter space (including extreme tails and K = 63 thresholds).
// Uses NicoSIMD's k_log_Phi / k_log_inv_logit / k_log1m_exp / k_exp kernels (AVX-512 on 8-lane hosts,
// AVX2 on 4-lane hosts such as ALICE with zero AVX-512 instructions) and one reverse_pass_callback node.

#include "stan_external_functions_on_Stan_math_only.hpp"

#include <stan/math/rev.hpp>
#include <stan/math/prim.hpp>
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <ostream>
#include <stdexcept>
#include <vector>

namespace fused_stan_math::set0_marginal_detail {

static constexpr double HALF_LOG_TWO_PI = 0.918938533204672741780329736406;

static constexpr double GH24_T[24] = {
    -8.50780351919527433e+00,
    -7.43789066602167104e+00,
    -6.54167500509862876e+00,
    -5.73274717525120181e+00,
    -4.97804137463911367e+00,
    -4.26038360501990354e+00,
    -3.56930676407355918e+00,
    -2.89772864322331181e+00,
    -2.24046785169175244e+00,
    -1.59348042981641846e+00,
    -9.53421922932108146e-01,
    -3.17370096629455034e-01,
     3.17370096629453258e-01,
     9.53421922932111698e-01,
     1.59348042981641846e+00,
     2.24046785169175156e+00,
     2.89772864322331269e+00,
     3.56930676407356007e+00,
     4.26038360501990709e+00,
     4.97804137463911989e+00,
     5.73274717525120003e+00,
     6.54167500509863320e+00,
     7.43789066602166216e+00,
     8.50780351919525657e+00
};

static constexpr double GH24_LW[24] = {
    -3.69042806607304854e+01,
    -2.86212344883108756e+01,
    -2.24843031451895179e+01,
    -1.76020197436324572e+01,
    -1.36185797828409356e+01,
    -1.03468101270751411e+01,
    -7.67407824125704963e+00,
    -5.52732585071174043e+00,
    -3.85723446710271611e+00,
    -2.63012623460007200e+00,
    -1.82350086207048112e+00,
    -1.42349743041213506e+00,
    -1.42349743041215793e+00,
    -1.82350086207048467e+00,
    -2.63012623460007866e+00,
    -3.85723446710268059e+00,
    -5.52732585071170757e+00,
    -7.67407824125703453e+00,
    -1.03468101270751269e+01,
    -1.36185797828409267e+01,
    -1.76020197436324182e+01,
    -2.24843031451895001e+01,
    -2.86212344883107441e+01,
    -3.69042806607303291e+01
};

inline void eval_log_integrand(const double b,
                               const double* x,
                               const double* w_seg,
                               const int m,
                               const double inv_v0,
                               const double lg_const,
                               const int link,
                               double* z_buf,
                               double* lF_buf,
                               double* lS_buf,
                               double* lf_buf,
                               double* log_q_buf,
                               double* work_buf,
                               double* r_own_buf,
                               double* r_nxt_buf,
                               double* dg_dz_out,
                               double& val_out,
                               double& dg_db_out,
                               double& d2g_db2_out) {
    for (int j = 0; j < m; ++j) {
        const double zj = x[j] + b;
        z_buf[j]  = zj;
        lF_buf[j] = zj;
        lS_buf[j] = -zj;
    }
    if (link == 1) {
        fused_stan_math::k_log_Phi(lF_buf, m);
        fused_stan_math::k_log_Phi(lS_buf, m);
        for (int j = 0; j < m; ++j) {
            lf_buf[j] = -0.5 * z_buf[j] * z_buf[j] - HALF_LOG_TWO_PI;
        }
    } else {
        fused_stan_math::k_log_inv_logit(lF_buf, m);
        fused_stan_math::k_log_inv_logit(lS_buf, m);
        for (int j = 0; j < m; ++j) {
            lf_buf[j] = lF_buf[j] + lS_buf[j];
        }
    }

    log_q_buf[0] = lF_buf[0];
    log_q_buf[m] = lS_buf[m - 1];
    if (m > 1) {
        for (int j = 1; j < m; ++j) {
            work_buf[j - 1] = (z_buf[j - 1] >= 0.0) ? (lS_buf[j] - lS_buf[j - 1]) : (lF_buf[j - 1] - lF_buf[j]);
        }
        fused_stan_math::k_log1m_exp(work_buf, m - 1);
        for (int j = 1; j < m; ++j) {
            log_q_buf[j] = ((z_buf[j - 1] >= 0.0) ? lS_buf[j - 1] : lF_buf[j]) + work_buf[j - 1];
        }
    }

    double val = -0.5 * b * b * inv_v0 + lg_const;
    for (int i = 0; i <= m; ++i) {
        if (w_seg[i] != 1.0) val += (w_seg[i] - 1.0) * log_q_buf[i];
    }
    for (int j = 0; j < m; ++j) {
        val += lf_buf[j];
        r_own_buf[j] = lf_buf[j] - log_q_buf[j];
        r_nxt_buf[j] = lf_buf[j] - log_q_buf[j + 1];
    }
    fused_stan_math::k_exp(r_own_buf, m);
    fused_stan_math::k_exp(r_nxt_buf, m);

    double sum_dg_dz = 0.0;
    for (int j = 0; j < m; ++j) {
        const double d_f = (link == 1) ? (-z_buf[j]) : (std::exp(lS_buf[j]) - std::exp(lF_buf[j]));
        const double g_j = (w_seg[j] - 1.0) * r_own_buf[j] - (w_seg[j + 1] - 1.0) * r_nxt_buf[j] + d_f;
        dg_dz_out[j] = g_j;
        sum_dg_dz += g_j;
    }
    dg_db_out = -b * inv_v0 + sum_dg_dz;

    double sum_d2 = -inv_v0;
    for (int j = 0; j < m; ++j) {
        if (link == 1) {
            sum_d2 -= 1.0;
        } else {
            sum_d2 -= 2.0 * std::exp(lf_buf[j]);
        }
    }
    auto dlogf_dz = [&](const int idx) -> double {
        return (link == 1) ? (-z_buf[idx]) : (std::exp(lS_buf[idx]) - std::exp(lF_buf[idx]));
    };
    if (w_seg[0] != 1.0) {
        const double u0 = r_own_buf[0];
        sum_d2 += (w_seg[0] - 1.0) * (dlogf_dz(0) * u0 - u0 * u0);
    }
    if (w_seg[m] != 1.0) {
        const double um = r_nxt_buf[m - 1];
        sum_d2 += (w_seg[m] - 1.0) * (-dlogf_dz(m - 1) * um - um * um);
    }
    for (int i = 1; i < m; ++i) {
        if (w_seg[i] != 1.0) {
            const double r_hi = r_own_buf[i];
            const double r_lo = r_nxt_buf[i - 1];
            const double ui   = r_hi - r_lo;
            sum_d2 += (w_seg[i] - 1.0) * (dlogf_dz(i) * r_hi - dlogf_dz(i - 1) * r_lo - ui * ui);
        }
    }
    val_out = val;
    d2g_db2_out = sum_d2;
}

}  // namespace fused_stan_math::set0_marginal_detail

template <typename T_x, typename T_seg, typename T_v0>
inline stan::return_type_t<T_x>
fused_log_Set0_marginal_at_chosen_cutoffs(const T_x& x_chosen,
                                         const T_seg& alpha_of_the_segments,
                                         const T_v0& v0_arg,
                                         const int use_probit_link,
                                         std::ostream* pstream__) {
    (void) pstream__;
    namespace detail = fused_stan_math::set0_marginal_detail;
    using T = stan::return_type_t<T_x>;
    const double v0 = stan::math::value_of(v0_arg);
    const int m = static_cast<int>(x_chosen.size());
    if (m == 0) return 0.0;
    if (alpha_of_the_segments.size() != m + 1) {
        throw std::invalid_argument("fused_log_Set0_marginal_at_chosen_cutoffs: alpha_of_the_segments size != m + 1");
    }
    if (!(v0 > 0.0)) {
        throw std::domain_error("fused_log_Set0_marginal_at_chosen_cutoffs: v0 must be positive");
    }

    const auto& x_ref   = stan::math::to_ref(x_chosen);
    const auto& seg_ref = stan::math::to_ref(alpha_of_the_segments);

    static constexpr int STACK_MAX_M = 16;
    // Newton buffers: x_val(m), w_seg(m+1), z/lF/lS/lf/work/r_own/r_nxt/dg_dz (8m+4 + m)
    // Batched GH24 buffers (overlapping reuse after Newton finishes!): 24 * (9*m + 1) doubles
    const int newton_doubles = 10 * m + 5;
    const int gh_doubles     = 24 * (9 * m + 1);
    const int total_scratch  = std::max(newton_doubles, gh_doubles) + 2 * m + 1;
    double stack_buf[std::max(10 * STACK_MAX_M + 5, 24 * (9 * STACK_MAX_M + 1)) + 2 * STACK_MAX_M + 1];
    std::vector<double> heap_buf;
    double* scratch = (m <= STACK_MAX_M) ? stack_buf : (heap_buf.resize(total_scratch), heap_buf.data());

    double* x_val    = scratch;
    double* w_seg    = x_val + m;
    double* work_all = w_seg + (m + 1);

    for (int j = 0; j < m; ++j) {
        x_val[j] = stan::math::value_of(x_ref.coeff(j));
        if (!std::isfinite(x_val[j])) {
            throw std::domain_error("fused_log_Set0_marginal_at_chosen_cutoffs: x_chosen is not finite");
        }
        if (j > 0 && !(x_val[j] > x_val[j - 1])) {
            throw std::domain_error("fused_log_Set0_marginal_at_chosen_cutoffs: x_chosen is not strictly increasing");
        }
    }

    double sum_seg = 0.0;
    double sum_lgamma_seg = 0.0;
    for (int i = 0; i <= m; ++i) {
        w_seg[i] = stan::math::value_of(seg_ref.coeff(i));
        sum_seg += w_seg[i];
        sum_lgamma_seg += std::lgamma(w_seg[i]);
    }
    const double lg_const = std::lgamma(sum_seg) - sum_lgamma_seg;
    const double inv_v0   = 1.0 / v0;

    // Initial bracketed guess from Beta/Dirichlet mean probabilities
    double b_star = 0.0;
    double cum_w  = 0.0;
    for (int j = 0; j < m; ++j) {
        cum_w += w_seg[j];
        const double p_mean = cum_w / sum_seg;
        const double q_mean = (use_probit_link == 1)
                              ? stan::math::inv_Phi(p_mean)
                              : std::log(p_mean / (1.0 - p_mean));
        b_star += (q_mean - x_val[j]);
    }
    b_star /= static_cast<double>(m);

    double* z_buf     = work_all;
    double* lF_buf    = z_buf + m;
    double* lS_buf    = lF_buf + m;
    double* lf_buf    = lS_buf + m;
    double* log_q_buf = lf_buf + m;
    double* work_buf  = log_q_buf + (m + 1);
    double* r_own_buf = work_buf + m;
    double* r_nxt_buf = r_own_buf + m;
    double* dg_dz_tmp = r_nxt_buf + m;

    double val = 0.0, dg_db = 0.0, d2g_db2 = -inv_v0;
    for (int it = 0; it < 10; ++it) {
        detail::eval_log_integrand(b_star, x_val, w_seg, m, inv_v0, lg_const, use_probit_link,
                                   z_buf, lF_buf, lS_buf, lf_buf, log_q_buf, work_buf, r_own_buf, r_nxt_buf,
                                   dg_dz_tmp, val, dg_db, d2g_db2);
        double step = -dg_db / d2g_db2;
        if (step >  2.0) step =  2.0;
        if (step < -2.0) step = -2.0;
        b_star += step;
        if (std::abs(step) <= 1e-14 * (1.0 + std::abs(b_star))) break;
    }
    detail::eval_log_integrand(b_star, x_val, w_seg, m, inv_v0, lg_const, use_probit_link,
                               z_buf, lF_buf, lS_buf, lf_buf, log_q_buf, work_buf, r_own_buf, r_nxt_buf,
                               dg_dz_tmp, val, dg_db, d2g_db2);
    const double sig_star = 1.0 / std::sqrt(-d2g_db2);

    // Batched SIMD evaluation of ALL 24 Gauss-Hermite nodes at once (24*m elements per kernel call)
    const int N24 = 24 * m;
    double* z_24      = work_all;
    double* lF_24     = z_24 + N24;
    double* lS_24     = lF_24 + N24;
    double* lf_24     = lS_24 + N24;
    double* log_q_24  = lf_24 + N24;
    double* work_24   = log_q_24 + 24 * (m + 1);
    double* r_own_24  = work_24 + N24;
    double* r_nxt_24  = r_own_24 + N24;
    double* grad_nodes = r_nxt_24 + N24;

    double b_24[24];
    for (int i = 0; i < 24; ++i) {
        const double b_i = b_star + sig_star * detail::GH24_T[i];
        b_24[i] = b_i;
        for (int j = 0; j < m; ++j) {
            const double zij = x_val[j] + b_i;
            const int idx = i * m + j;
            z_24[idx]  = zij;
            lF_24[idx] = zij;
            lS_24[idx] = -zij;
        }
    }
    if (use_probit_link == 1) {
        fused_stan_math::k_log_Phi(lF_24, N24);
        fused_stan_math::k_log_Phi(lS_24, N24);
        for (int idx = 0; idx < N24; ++idx) {
            lf_24[idx] = -0.5 * z_24[idx] * z_24[idx] - detail::HALF_LOG_TWO_PI;
        }
    } else {
        fused_stan_math::k_log_inv_logit(lF_24, N24);
        fused_stan_math::k_log_inv_logit(lS_24, N24);
        for (int idx = 0; idx < N24; ++idx) {
            lf_24[idx] = lF_24[idx] + lS_24[idx];
        }
    }
    if (m > 1) {
        int wi = 0;
        for (int i = 0; i < 24; ++i) {
            const int base = i * m;
            for (int j = 1; j < m; ++j) {
                work_24[wi++] = (z_24[base + j - 1] >= 0.0)
                    ? (lS_24[base + j] - lS_24[base + j - 1])
                    : (lF_24[base + j - 1] - lF_24[base + j]);
            }
        }
        fused_stan_math::k_log1m_exp(work_24, 24 * (m - 1));
        wi = 0;
        for (int i = 0; i < 24; ++i) {
            const int base = i * m;
            const int qbase = i * (m + 1);
            log_q_24[qbase]     = lF_24[base];
            log_q_24[qbase + m] = lS_24[base + m - 1];
            for (int j = 1; j < m; ++j) {
                log_q_24[qbase + j] = ((z_24[base + j - 1] >= 0.0) ? lS_24[base + j - 1] : lF_24[base + j])
                                      + work_24[wi++];
            }
        }
    } else {
        for (int i = 0; i < 24; ++i) {
            log_q_24[2 * i]     = lF_24[i];
            log_q_24[2 * i + 1] = lS_24[i];
        }
    }
    for (int i = 0; i < 24; ++i) {
        const int base  = i * m;
        const int qbase = i * (m + 1);
        for (int j = 0; j < m; ++j) {
            r_own_24[base + j] = lf_24[base + j] - log_q_24[qbase + j];
            r_nxt_24[base + j] = lf_24[base + j] - log_q_24[qbase + j + 1];
        }
    }
    fused_stan_math::k_exp(r_own_24, N24);
    fused_stan_math::k_exp(r_nxt_24, N24);
    if constexpr (stan::is_var<T>::value) {
        if (use_probit_link != 1) {
            fused_stan_math::k_exp(lF_24, N24);
            fused_stan_math::k_exp(lS_24, N24);
        }
    }

    double lw_terms[24];
    double max_lw = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < 24; ++i) {
        const double t_i = detail::GH24_T[i];
        const double b_i = b_24[i];
        const int base  = i * m;
        const int qbase = i * (m + 1);
        double val_i = -0.5 * b_i * b_i * inv_v0 + lg_const;
        for (int s_i = 0; s_i <= m; ++s_i) {
            if (w_seg[s_i] != 1.0) val_i += (w_seg[s_i] - 1.0) * log_q_24[qbase + s_i];
        }
        for (int j = 0; j < m; ++j) {
            val_i += lf_24[base + j];
            if constexpr (stan::is_var<T>::value) {
                const double d_f = (use_probit_link == 1) ? (-z_24[base + j]) : (lS_24[base + j] - lF_24[base + j]);
                grad_nodes[base + j] = (w_seg[j] - 1.0) * r_own_24[base + j]
                                       - (w_seg[j + 1] - 1.0) * r_nxt_24[base + j] + d_f;
            }
        }
        const double lw_i = detail::GH24_LW[i] + val_i + 0.5 * t_i * t_i;
        lw_terms[i] = lw_i;
        if (lw_i > max_lw) max_lw = lw_i;
    }

    double sum_exp = 0.0;
    double w_norm[24];
    for (int i = 0; i < 24; ++i) {
        w_norm[i] = lw_terms[i] - max_lw;
    }
    fused_stan_math::k_exp(w_norm, 24);
    for (int i = 0; i < 24; ++i) {
        sum_exp += w_norm[i];
    }
    const double inv_sum = 1.0 / sum_exp;
    const double log_marginal = -0.5 * std::log(v0) + std::log(sig_star) + max_lw + std::log(sum_exp);

    if constexpr (!stan::is_var<T>::value) {
        return log_marginal;
    } else {
        stan::arena_t<Eigen::Matrix<stan::math::var, -1, 1>> arena_x = x_ref;
        double* arena_grad = stan::math::ChainableStack::instance_->memalloc_.alloc_array<double>(m);
        for (int j = 0; j < m; ++j) {
            double g_j = 0.0;
            for (int i = 0; i < 24; ++i) {
                g_j += (w_norm[i] * inv_sum) * grad_nodes[i * m + j];
            }
            arena_grad[j] = g_j;
        }
        stan::math::var res_var(new stan::math::vari(log_marginal, false));
        stan::math::reverse_pass_callback([arena_x, arena_grad, res_var, m]() mutable {
            const double adj = res_var.adj();
            for (int j = 0; j < m; ++j) {
                arena_x.coeffRef(j).adj() += adj * arena_grad[j];
            }
        });
        return res_var;
    }
}


// ---- Batched multi-width Set-0 Gaussian convolution marginal (shared by O_biv_FC & O_biv_RC) ----
namespace fused_stan_math::set0_marginal_detail {

// Scratch size required by eval_set0_marginal_fast_raw:
//   [z_24 (N24), lF_24 (N24), lS_24 (N24), lf_24 (N24), log_q_24 (24*(m+1)),
//    work_24 (N24), r_own_24 (N24), r_nxt_24 (N24), grad_nodes (N24)] = 24*(9*m+1) + 2*m
inline constexpr int raw_scratch_size(const int m) {
    return std::max(10 * m + 5, 24 * (9 * m + 1)) + 2 * m + 4;
}

// Fast batchedSet0 marginal evaluator directly on double arrays (zero Stan vars or allocations)
// Notice [lF_buf (m), lS_buf (m)] contiguous in Newton, and [lF_24 (N24), lS_24 (N24)] and
// [r_own_24 (N24), r_nxt_24 (N24)] contiguous in GH24 -> single SIMD call each!
template <bool need_grad>
inline double eval_set0_marginal_fast_raw(const double* x_val,
                                          const double* w_seg,
                                          const int m,
                                          const double v0,
                                          const int use_probit_link,
                                          const bool all_unit,
                                          double* grad_x_out,
                                          double* scratch) {
    if (m == 0) return 0.0;
    for (int j = 0; j < m; ++j) {
        if (!std::isfinite(x_val[j])) {
            stan::math::check_finite("fused_log_Set0_marginal_at_chosen_cutoffs", "x_chosen", x_val[j]);
        }
        if (j > 0 && !(x_val[j] > x_val[j - 1])) {
            throw std::domain_error("fused_log_Set0_marginal_at_chosen_cutoffs: x_chosen is not strictly increasing");
        }
    }
    if (use_probit_link == 1 && all_unit) {
        const double inv_denom = 1.0 / (1.0 + static_cast<double>(m) * v0);
        const double c = v0 * inv_denom;
        double sx = 0.0, sx2 = 0.0;
        for (int j = 0; j < m; ++j) {
            sx += x_val[j];
            sx2 += x_val[j] * x_val[j];
        }
        if constexpr (need_grad) {
            const double shift = c * sx;
            for (int j = 0; j < m; ++j) grad_x_out[j] = -(x_val[j] - shift);
        }
        return std::lgamma(static_cast<double>(m) + 1.0)
             - static_cast<double>(m) * HALF_LOG_TWO_PI
             + 0.5 * std::log(inv_denom)
             - 0.5 * (sx2 - c * sx * sx);
    }

    double sum_seg = 0.0;
    double sum_lgamma_seg = 0.0;
    for (int i = 0; i <= m; ++i) {
        sum_seg += w_seg[i];
        sum_lgamma_seg += std::lgamma(w_seg[i]);
    }
    const double lg_const = std::lgamma(sum_seg) - sum_lgamma_seg;
    const double inv_v0 = 1.0 / v0;

    double b_star = 0.0;
    double cum_w = 0.0;
    for (int j = 0; j < m; ++j) {
        cum_w += w_seg[j];
        const double p_mean = cum_w / sum_seg;
        const double q_mean = (use_probit_link == 1)
            ? stan::math::inv_Phi(p_mean)
            : std::log(p_mean / (1.0 - p_mean));
        b_star += (q_mean - x_val[j]);
    }
    b_star /= static_cast<double>(m);

    // Contiguous Newton layout: [z(m), lFS(2m), lf(m), log_q(m+1), work(m), r_both(2m), dg_dz(m)]
    double* z_buf = scratch;
    double* lF_buf = z_buf + m;
    double* lS_buf = lF_buf + m;
    double* lf_buf = lS_buf + m;
    double* log_q_buf = lf_buf + m;
    double* work_buf = log_q_buf + (m + 1);
    double* r_own_buf = work_buf + m;
    double* r_nxt_buf = r_own_buf + m;

    auto eval_newton = [&](const double b, double& dg_db_out, double& d2g_db2_out) {
        for (int j = 0; j < m; ++j) {
            const double zj = x_val[j] + b;
            z_buf[j] = zj;
            lF_buf[j] = zj;
            lS_buf[j] = -zj;
        }
        if (use_probit_link == 1) {
            fused_stan_math::k_log_Phi(lF_buf, 2 * m);  // batched [lF; lS]
            for (int j = 0; j < m; ++j) lf_buf[j] = -0.5 * z_buf[j] * z_buf[j] - HALF_LOG_TWO_PI;
        } else {
            fused_stan_math::k_log_inv_logit(lF_buf, 2 * m);
            for (int j = 0; j < m; ++j) lf_buf[j] = lF_buf[j] + lS_buf[j];
        }
        log_q_buf[0] = lF_buf[0];
        log_q_buf[m] = lS_buf[m - 1];
        if (m > 1) {
            for (int j = 1; j < m; ++j) {
                work_buf[j - 1] = (z_buf[j - 1] >= 0.0) ? (lS_buf[j] - lS_buf[j - 1]) : (lF_buf[j - 1] - lF_buf[j]);
            }
            fused_stan_math::k_log1m_exp(work_buf, m - 1);
            for (int j = 1; j < m; ++j) {
                log_q_buf[j] = ((z_buf[j - 1] >= 0.0) ? lS_buf[j - 1] : lF_buf[j]) + work_buf[j - 1];
            }
        }
        for (int j = 0; j < m; ++j) {
            r_own_buf[j] = lf_buf[j] - log_q_buf[j];
            r_nxt_buf[j] = lf_buf[j] - log_q_buf[j + 1];
        }
        fused_stan_math::k_exp(r_own_buf, 2 * m);  // batched [r_own; r_nxt]
        double sum_dg_dz = 0.0;
        double sum_d2 = -inv_v0;
        if (use_probit_link == 1) {
            for (int j = 0; j < m; ++j) {
                sum_dg_dz += (w_seg[j] - 1.0) * r_own_buf[j] - (w_seg[j + 1] - 1.0) * r_nxt_buf[j] - z_buf[j];
            }
            sum_d2 -= static_cast<double>(m);
            if (w_seg[0] != 1.0) {
                const double u0 = r_own_buf[0];
                sum_d2 += (w_seg[0] - 1.0) * (-z_buf[0] * u0 - u0 * u0);
            }
            if (w_seg[m] != 1.0) {
                const double um = r_nxt_buf[m - 1];
                sum_d2 += (w_seg[m] - 1.0) * (z_buf[m - 1] * um - um * um);
            }
            for (int i = 1; i < m; ++i) {
                if (w_seg[i] != 1.0) {
                    const double r_hi = r_own_buf[i];
                    const double r_lo = r_nxt_buf[i - 1];
                    const double ui = r_hi - r_lo;
                    sum_d2 += (w_seg[i] - 1.0) * (-z_buf[i] * r_hi + z_buf[i - 1] * r_lo - ui * ui);
                }
            }
        } else {
            for (int j = 0; j < m; ++j) {
                const double df = std::exp(lS_buf[j]) - std::exp(lF_buf[j]);
                sum_dg_dz += (w_seg[j] - 1.0) * r_own_buf[j] - (w_seg[j + 1] - 1.0) * r_nxt_buf[j] + df;
                sum_d2 -= 2.0 * std::exp(lf_buf[j]);
            }
            auto dlogf_l = [&](int idx) { return std::exp(lS_buf[idx]) - std::exp(lF_buf[idx]); };
            if (w_seg[0] != 1.0) {
                const double u0 = r_own_buf[0];
                sum_d2 += (w_seg[0] - 1.0) * (dlogf_l(0) * u0 - u0 * u0);
            }
            if (w_seg[m] != 1.0) {
                const double um = r_nxt_buf[m - 1];
                sum_d2 += (w_seg[m] - 1.0) * (-dlogf_l(m - 1) * um - um * um);
            }
            for (int i = 1; i < m; ++i) {
                if (w_seg[i] != 1.0) {
                    const double r_hi = r_own_buf[i];
                    const double r_lo = r_nxt_buf[i - 1];
                    const double ui = r_hi - r_lo;
                    sum_d2 += (w_seg[i] - 1.0) * (dlogf_l(i) * r_hi - dlogf_l(i - 1) * r_lo - ui * ui);
                }
            }
        }
        dg_db_out = -b * inv_v0 + sum_dg_dz;
        d2g_db2_out = sum_d2;
    };

    double dg_db = 0.0, d2g_db2 = -inv_v0;
    for (int it = 0; it < 10; ++it) {
        eval_newton(b_star, dg_db, d2g_db2);
        double step = -dg_db / d2g_db2;
        if (step > 2.0) step = 2.0;
        if (step < -2.0) step = -2.0;
        b_star += step;
        if (std::abs(step) <= 1e-13 * (1.0 + std::abs(b_star))) break;
    }
    eval_newton(b_star, dg_db, d2g_db2);
    const double sig_star = 1.0 / std::sqrt(-d2g_db2);

    // Batched 24-pt Gauss-Hermite with contiguous [lF_24; lS_24] and [r_own_24; r_nxt_24]
    const int N24 = 24 * m;
    double* z_24 = scratch;
    double* lF_24 = z_24 + N24;
    double* lS_24 = lF_24 + N24;
    double* lf_24 = lS_24 + N24;
    double* log_q_24 = lf_24 + N24;
    double* work_24 = log_q_24 + 24 * (m + 1);
    double* r_own_24 = work_24 + N24;
    double* r_nxt_24 = r_own_24 + N24;
    double* grad_nodes = r_nxt_24 + N24;

    double b_24[24];
    for (int i = 0; i < 24; ++i) {
        const double bi = b_star + sig_star * GH24_T[i];
        b_24[i] = bi;
        const int base = i * m;
        for (int j = 0; j < m; ++j) {
            const double zij = x_val[j] + bi;
            z_24[base + j] = zij;
            lF_24[base + j] = zij;
            lS_24[base + j] = -zij;
        }
    }
    if (use_probit_link == 1) {
        fused_stan_math::k_log_Phi(lF_24, 2 * N24);
        for (int idx = 0; idx < N24; ++idx) {
            lf_24[idx] = -0.5 * z_24[idx] * z_24[idx] - HALF_LOG_TWO_PI;
        }
    } else {
        fused_stan_math::k_log_inv_logit(lF_24, 2 * N24);
        for (int idx = 0; idx < N24; ++idx) {
            lf_24[idx] = lF_24[idx] + lS_24[idx];
        }
    }
    if (m > 1) {
        int wi = 0;
        for (int i = 0; i < 24; ++i) {
            const int base = i * m;
            for (int j = 1; j < m; ++j) {
                work_24[wi++] = (z_24[base + j - 1] >= 0.0)
                    ? (lS_24[base + j] - lS_24[base + j - 1])
                    : (lF_24[base + j - 1] - lF_24[base + j]);
            }
        }
        fused_stan_math::k_log1m_exp(work_24, 24 * (m - 1));
        wi = 0;
        for (int i = 0; i < 24; ++i) {
            const int base = i * m;
            const int qbase = i * (m + 1);
            log_q_24[qbase] = lF_24[base];
            log_q_24[qbase + m] = lS_24[base + m - 1];
            for (int j = 1; j < m; ++j) {
                log_q_24[qbase + j] = ((z_24[base + j - 1] >= 0.0) ? lS_24[base + j - 1] : lF_24[base + j])
                                    + work_24[wi++];
            }
        }
    } else {
        for (int i = 0; i < 24; ++i) {
            log_q_24[2 * i] = lF_24[i];
            log_q_24[2 * i + 1] = lS_24[i];
        }
    }
    for (int i = 0; i < 24; ++i) {
        const int base = i * m;
        const int qbase = i * (m + 1);
        for (int j = 0; j < m; ++j) {
            r_own_24[base + j] = lf_24[base + j] - log_q_24[qbase + j];
            r_nxt_24[base + j] = lf_24[base + j] - log_q_24[qbase + j + 1];
        }
    }
    fused_stan_math::k_exp(r_own_24, 2 * N24);
    if constexpr (need_grad) {
        if (use_probit_link != 1) {
            fused_stan_math::k_exp(lF_24, 2 * N24);
        }
    }

    double lw_terms[24];
    double max_lw = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < 24; ++i) {
        const double t_i = GH24_T[i];
        const double b_i = b_24[i];
        const int base = i * m;
        const int qbase = i * (m + 1);
        double val_i = -0.5 * b_i * b_i * inv_v0 + lg_const;
        for (int s_i = 0; s_i <= m; ++s_i) {
            if (w_seg[s_i] != 1.0) val_i += (w_seg[s_i] - 1.0) * log_q_24[qbase + s_i];
        }
        for (int j = 0; j < m; ++j) {
            val_i += lf_24[base + j];
            if constexpr (need_grad) {
                const double d_f = (use_probit_link == 1) ? (-z_24[base + j]) : (lS_24[base + j] - lF_24[base + j]);
                grad_nodes[base + j] = (w_seg[j] - 1.0) * r_own_24[base + j]
                                     - (w_seg[j + 1] - 1.0) * r_nxt_24[base + j] + d_f;
            }
        }
        const double lw_i = GH24_LW[i] + val_i + 0.5 * t_i * t_i;
        lw_terms[i] = lw_i;
        if (lw_i > max_lw) max_lw = lw_i;
    }

    double w_norm[24];
    for (int i = 0; i < 24; ++i) w_norm[i] = lw_terms[i] - max_lw;
    fused_stan_math::k_exp(w_norm, 24);
    double sum_exp = 0.0;
    for (int i = 0; i < 24; ++i) sum_exp += w_norm[i];
    const double inv_sum = 1.0 / sum_exp;
    if constexpr (need_grad) {
        for (int j = 0; j < m; ++j) {
            double g_j = 0.0;
            for (int i = 0; i < 24; ++i) g_j += (w_norm[i] * inv_sum) * grad_nodes[i * m + j];
            grad_x_out[j] = g_j;
        }
    }
    return -0.5 * std::log(v0) + std::log(sig_star) + max_lw + std::log(sum_exp);
}

}  // namespace fused_stan_math::set0_marginal_detail

namespace fused_stan_math {

template <typename T_x, typename T_w, typename T_v0,
          stan::require_eigen_col_vector_t<T_x>* = nullptr,
          stan::require_eigen_col_vector_t<T_w>* = nullptr>
inline stan::return_type_t<T_x, T_w, T_v0>
fused_log_Set0_marginal_at_chosen_cutoffs_fast(const T_x& x_chosen,
                                              const T_w& alpha_of_the_segments,
                                              const T_v0& v0_arg,
                                              const int use_probit_link,
                                              std::ostream* pstream__) {
    (void) pstream__;
    namespace detail = fused_stan_math::set0_marginal_detail;
    using T = stan::return_type_t<T_x, T_w, T_v0>;
    const int m = static_cast<int>(x_chosen.size());
    if (m == 0) return 0.0;
    if (static_cast<int>(alpha_of_the_segments.size()) != m + 1) {
        throw std::invalid_argument(
            "fused_log_Set0_marginal_at_chosen_cutoffs_fast: alpha_of_the_segments must have length m + 1");
    }
    if constexpr (!stan::is_constant_all<T_w, T_v0>::value) {
        return fused_log_Set0_marginal_at_chosen_cutoffs(
            x_chosen, alpha_of_the_segments, v0_arg, use_probit_link, pstream__);
    } else {
        const double v0 = stan::math::value_of(v0_arg);
        if (!(v0 > 0.0) || !std::isfinite(v0)) {
            throw std::domain_error("fused_log_Set0_marginal_at_chosen_cutoffs_fast: v0 must be positive");
        }
        const auto& x_ref = stan::math::to_ref(x_chosen);
        const auto& w_ref = stan::math::to_ref(alpha_of_the_segments);
        static constexpr int STACK_MAX_M = 16;
        double stack_buf[detail::raw_scratch_size(STACK_MAX_M) + 2 * STACK_MAX_M + 2];
        std::vector<double> heap_buf;
        const std::size_t need = static_cast<std::size_t>(detail::raw_scratch_size(m) + 2 * m + 2);
        double* buf = (m <= STACK_MAX_M) ? stack_buf : (heap_buf.resize(need), heap_buf.data());
        double* x_val = buf;
        double* w_seg = x_val + m;
        double* scratch = w_seg + (m + 1);
        for (int j = 0; j < m; ++j) x_val[j] = stan::math::value_of(x_ref.coeff(j));
        bool all_unit = true;
        for (int i = 0; i <= m; ++i) {
            w_seg[i] = stan::math::value_of(w_ref.coeff(i));
            if (w_seg[i] != 1.0) all_unit = false;
        }
        if constexpr (stan::is_constant_all<T_x>::value) {
            return detail::eval_set0_marginal_fast_raw<false>(
                x_val, w_seg, m, v0, use_probit_link, all_unit, nullptr, scratch);
        } else {
            double* arena_grad = stan::math::ChainableStack::instance_->memalloc_.alloc_array<double>(m);
            const double val = detail::eval_set0_marginal_fast_raw<true>(
                x_val, w_seg, m, v0, use_probit_link, all_unit, arena_grad, scratch);
            stan::arena_t<Eigen::Matrix<stan::math::var, -1, 1>> arena_x = x_ref;
            stan::math::var res_var(new stan::math::vari(val, false));
            stan::math::reverse_pass_callback([arena_x, arena_grad, res_var, m]() mutable {
                const double adj = res_var.adj();
                for (int j = 0; j < m; ++j) arena_x.coeffRef(j).adj() += adj * arena_grad[j];
            });
            return T(res_var);
        }
    }
}

}  // namespace fused_stan_math

using fused_stan_math::fused_log_Set0_marginal_at_chosen_cutoffs_fast;
