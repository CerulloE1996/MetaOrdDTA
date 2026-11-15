#pragma once
// ---- NicoSIMD fused Sean Pinkney partially-centered study cutpoints and total log-Jacobian ----
#if __has_include("stan_external_functions.hpp")
#include "stan_external_functions.hpp"
#include "stan_external_dirichlet_cutpoint_median.hpp"
#elif __has_include("NicoSIMD/stan_external_functions.hpp")
#include "NicoSIMD/stan_external_functions.hpp"
#include "NicoSIMD/stan_external_dirichlet_cutpoint_median.hpp"
#endif

#include <stan/math/rev.hpp>
#include <stan/math/prim.hpp>
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

template <typename T_mu,
          typename T_tau,
          typename T_z,
          typename T_w,
          typename T_s,
          typename T_l,
          typename T_wd>
inline Eigen::Matrix<stan::return_type_t<T_mu, T_tau, T_z>, -1, 1>
nico_pinkney_C_and_log_J_from_pairs(const T_mu& logit_centre_of_pair,
                                    const T_tau& spread_free_of_pair,
                                    const T_z& C_raw_c,
                                    const T_w& Pinkney_w_c,
                                    const T_s& logit_start,
                                    const T_l& C_raw_lower_c,
                                    const T_wd& C_raw_width,
                                    const double sum_log_width,
                                    const std::vector<int>& pair_of_each_C,
                                    const std::vector<int>& is_increment_of_each_C,
                                    const int softplus,
                                    std::ostream* pstream__) {
    (void) pstream__;
    using T = stan::return_type_t<T_mu, T_tau, T_z>;
    using T_mu_scalar  = stan::value_type_t<T_mu>;
    using T_tau_scalar = stan::value_type_t<T_tau>;
    using T_z_scalar   = stan::value_type_t<T_z>;

    const Eigen::Index n_C = C_raw_c.size();
    const Eigen::Index n_pairs = logit_centre_of_pair.size();
    if (n_C == 0) {
        Eigen::Matrix<T, -1, 1> empty_res(1);
        empty_res[0] = 0.0;
        return empty_res;
    }

    const auto& mu_ref  = stan::math::to_ref(logit_centre_of_pair);
    const auto& tau_ref = stan::math::to_ref(spread_free_of_pair);
    const auto& z_ref   = stan::math::to_ref(C_raw_c);
    const auto& w_ref   = stan::math::to_ref(Pinkney_w_c);
    const auto& s_ref   = stan::math::to_ref(logit_start);
    const auto& l_ref   = stan::math::to_ref(C_raw_lower_c);
    const auto& wd_ref  = stan::math::to_ref(C_raw_width);

    Eigen::VectorXd mu_val(n_pairs);
    Eigen::VectorXd tau_val(n_pairs);
    Eigen::VectorXd log_tau_val(n_pairs);
    for (Eigen::Index u = 0; u < n_pairs; ++u) {
        mu_val[u]  = stan::math::value_of(mu_ref.coeff(u));
        tau_val[u] = stan::math::value_of(tau_ref.coeff(u));
    }
    log_tau_val = tau_val;
    nico::math::k_log(log_tau_val.data(), static_cast<int>(n_pairs));

    if constexpr (!stan::is_var<T>::value) {
        Eigen::VectorXd theta(n_C);
        Eigen::VectorXd abs_theta(n_C);
        Eigen::VectorXd D_vec(n_C);
        double sum_log_tau_over_D = 0.0;
        for (Eigen::Index p = 0; p < n_C; ++p) {
            const int u = pair_of_each_C[p] - 1;
            const double w_p = stan::math::value_of(w_ref.coeff(p));
            const double tau_u = tau_val[u];
            const double D_p = tau_u - w_p * (tau_u - 1.0);
            const double eta_sh = (1.0 - w_p) * stan::math::value_of(s_ref.coeff(p))
                                  + stan::math::value_of(z_ref.coeff(p));
            const double th = (w_p * mu_val[u] + tau_u * eta_sh) / D_p;
            theta[p] = th;
            abs_theta[p] = std::abs(th);
            D_vec[p] = D_p;
            sum_log_tau_over_D += log_tau_val[u];
        }
        Eigen::VectorXd log_D = D_vec;
        nico::math::k_log(log_D.data(), static_cast<int>(n_C));
        for (Eigen::Index p = 0; p < n_C; ++p) sum_log_tau_over_D -= log_D[p];

        Eigen::VectorXd q = theta;
        nico::math::k_inv_logit(q.data(), static_cast<int>(n_C));
        Eigen::VectorXd log_inv_abs = abs_theta;
        nico::math::k_log_inv_logit(log_inv_abs.data(), static_cast<int>(n_C));

        Eigen::VectorXd r_vec(n_C);
        double log_J = sum_log_width + sum_log_tau_over_D;
        for (Eigen::Index p = 0; p < n_C; ++p) {
            r_vec[p] = stan::math::value_of(l_ref.coeff(p))
                       + stan::math::value_of(wd_ref.coeff(p)) * q[p];
            log_J += 2.0 * log_inv_abs[p] - abs_theta[p];
        }

        Eigen::VectorXd step_val = r_vec;
        if (softplus == 1) {
            nico::math::k_log1p_exp(step_val.data(), static_cast<int>(n_C));
            Eigen::VectorXd lil_r = r_vec;
            nico::math::k_log_inv_logit(lil_r.data(), static_cast<int>(n_C));
            for (Eigen::Index p = 0; p < n_C; ++p) {
                if (is_increment_of_each_C[p] == 1) log_J += lil_r[p];
            }
        } else {
            nico::math::k_exp(step_val.data(), static_cast<int>(n_C));
            for (Eigen::Index p = 0; p < n_C; ++p) {
                if (is_increment_of_each_C[p] == 1) log_J += r_vec[p];
            }
        }

        Eigen::Matrix<double, -1, 1> result(n_C + 1);
        for (Eigen::Index p = 0; p < n_C; ++p) {
            result[p] = (is_increment_of_each_C[p] == 1) ? (result[p - 1] + step_val[p]) : r_vec[p];
        }
        result[n_C] = log_J;
        return result;
    } else {
        stan::arena_t<Eigen::Matrix<T_mu_scalar,  -1, 1>> arena_mu  = mu_ref;
        stan::arena_t<Eigen::Matrix<T_tau_scalar, -1, 1>> arena_tau = tau_ref;
        stan::arena_t<Eigen::Matrix<T_z_scalar,   -1, 1>> arena_z   = z_ref;
        stan::arena_t<Eigen::Matrix<stan::math::var, -1, 1>> arena_result(n_C + 1);
        stan::arena_t<Eigen::Matrix<double, 8, -1>> arena_buf(8, n_C);
        int* arena_pair = stan::math::ChainableStack::instance_->memalloc_.alloc_array<int>(n_C);
        int* arena_inc  = stan::math::ChainableStack::instance_->memalloc_.alloc_array<int>(n_C);

        Eigen::VectorXd theta(n_C);
        Eigen::VectorXd abs_theta(n_C);
        Eigen::VectorXd D_vec(n_C);
        double sum_log_tau_over_D = 0.0;

        for (Eigen::Index p = 0; p < n_C; ++p) {
            const int u = pair_of_each_C[p] - 1;
            arena_pair[p] = u;
            arena_inc[p]  = is_increment_of_each_C[p];
            const double w_p = stan::math::value_of(w_ref.coeff(p));
            const double tau_u = tau_val[u];
            const double D_p = tau_u - w_p * (tau_u - 1.0);
            const double inv_D_p = 1.0 / D_p;
            const double eta_sh = (1.0 - w_p) * stan::math::value_of(s_ref.coeff(p))
                                  + stan::math::value_of(arena_z.coeff(p));
            const double th = (w_p * mu_val[u] + tau_u * eta_sh) * inv_D_p;
            theta[p] = th;
            abs_theta[p] = std::abs(th);
            D_vec[p] = D_p;
            sum_log_tau_over_D += log_tau_val[u];
            arena_buf.coeffRef(0, p) = inv_D_p;
            arena_buf.coeffRef(1, p) = w_p;
            arena_buf.coeffRef(2, p) = tau_u;
            arena_buf.coeffRef(3, p) = eta_sh - (1.0 - w_p) * th;
            arena_buf.coeffRef(7, p) = 1.0 / tau_u;
        }

        Eigen::VectorXd log_D = D_vec;
        nico::math::k_log(log_D.data(), static_cast<int>(n_C));
        for (Eigen::Index p = 0; p < n_C; ++p) sum_log_tau_over_D -= log_D[p];

        Eigen::VectorXd q = theta;
        nico::math::k_inv_logit(q.data(), static_cast<int>(n_C));
        Eigen::VectorXd one_m_q = -theta;
        nico::math::k_inv_logit(one_m_q.data(), static_cast<int>(n_C));
        Eigen::VectorXd log_inv_abs = abs_theta;
        nico::math::k_log_inv_logit(log_inv_abs.data(), static_cast<int>(n_C));

        Eigen::VectorXd r_vec(n_C);
        double log_J = sum_log_width + sum_log_tau_over_D;
        for (Eigen::Index p = 0; p < n_C; ++p) {
            const double wd_p = stan::math::value_of(wd_ref.coeff(p));
            const double q_p = q[p];
            const double omq_p = one_m_q[p];
            r_vec[p] = stan::math::value_of(l_ref.coeff(p)) + wd_p * q_p;
            arena_buf.coeffRef(4, p) = q_p - omq_p;
            arena_buf.coeffRef(5, p) = wd_p * q_p * omq_p;
            log_J += 2.0 * log_inv_abs[p] - abs_theta[p];
        }

        Eigen::VectorXd step_val = r_vec;
        if (softplus == 1) {
            nico::math::k_log1p_exp(step_val.data(), static_cast<int>(n_C));
            Eigen::VectorXd s_r = r_vec;
            nico::math::k_inv_logit(s_r.data(), static_cast<int>(n_C));
            Eigen::VectorXd lil_r = r_vec;
            nico::math::k_log_inv_logit(lil_r.data(), static_cast<int>(n_C));
            for (Eigen::Index p = 0; p < n_C; ++p) {
                arena_buf.coeffRef(6, p) = s_r[p];
                if (arena_inc[p] == 1) log_J += lil_r[p];
            }
        } else {
            nico::math::k_exp(step_val.data(), static_cast<int>(n_C));
            for (Eigen::Index p = 0; p < n_C; ++p) {
                arena_buf.coeffRef(6, p) = step_val[p];
                if (arena_inc[p] == 1) log_J += r_vec[p];
            }
        }

        double prev_C = 0.0;
        for (Eigen::Index p = 0; p < n_C; ++p) {
            const double C_p = (arena_inc[p] == 1) ? (prev_C + step_val[p]) : r_vec[p];
            prev_C = C_p;
            arena_result.coeffRef(p) = stan::math::var(new stan::math::vari(C_p, false));
        }
        arena_result.coeffRef(n_C) = stan::math::var(new stan::math::vari(log_J, false));

        const bool softplus_flag = (softplus == 1);
        stan::math::reverse_pass_callback(
            [arena_mu, arena_tau, arena_z, arena_result, arena_buf, arena_pair, arena_inc, n_C, softplus_flag]() mutable {
                const double adj_J = arena_result.coeff(n_C).adj();
                double running_adj_C = 0.0;
                for (Eigen::Index p = n_C - 1; p >= 0; --p) {
                    const double adj_C_here = arena_result.coeff(p).adj() + running_adj_C;
                    const bool is_inc = (arena_inc[p] == 1);
                    running_adj_C = is_inc ? adj_C_here : 0.0;
                    double adj_r;
                    if (!is_inc) {
                        adj_r = adj_C_here;
                    } else if (!softplus_flag) {
                        adj_r = adj_C_here * arena_buf.coeff(6, p) + adj_J;
                    } else {
                        const double s_p = arena_buf.coeff(6, p);
                        adj_r = adj_C_here * s_p + adj_J * (1.0 - s_p);
                    }
                    const double adj_theta = adj_r * arena_buf.coeff(5, p) - adj_J * arena_buf.coeff(4, p);
                    const double inv_D_p = arena_buf.coeff(0, p);
                    const double w_p     = arena_buf.coeff(1, p);
                    const double tau_u   = arena_buf.coeff(2, p);
                    const int u          = arena_pair[p];
                    if constexpr (stan::is_var<T_z_scalar>::value) {
                        arena_z.coeffRef(p).adj() += adj_theta * (tau_u * inv_D_p);
                    }
                    if constexpr (stan::is_var<T_mu_scalar>::value) {
                        arena_mu.coeffRef(u).adj() += adj_theta * (w_p * inv_D_p);
                    }
                    if constexpr (stan::is_var<T_tau_scalar>::value) {
                        arena_tau.coeffRef(u).adj() += adj_theta * (arena_buf.coeff(3, p) * inv_D_p)
                                                       + adj_J * (w_p * arena_buf.coeff(7, p) * inv_D_p);
                    }
                }
            });
        return Eigen::Matrix<T, -1, 1>(arena_result);
    }
}

// ---- Fully fused Pinkney partially-centered study cutpoints and total log-Jacobian directly from alpha ----
template <typename T_alpha,
          typename T_z,
          typename T_loc,
          typename T_w,
          typename T_s,
          typename T_l,
          typename T_u>
inline Eigen::Matrix<stan::return_type_t<T_alpha, T_z>, -1, 1>
nico_pinkney_C_and_log_J_from_alpha(const T_alpha& alpha,
                                    const T_z& C_raw_c,
                                    const T_loc& baseline_location_prior_arg,
                                    const T_w& Pinkney_w_c,
                                    const T_s& logit_start,
                                    const T_l& C_raw_lower_c,
                                    const T_u& C_raw_upper_c,
                                    const std::vector<int>& k_of_each_C,
                                    const std::vector<int>& k_previous_of_each_C,
                                    const int use_probit_link,
                                    const int softplus,
                                    std::ostream* pstream__) {
    (void) pstream__;
    const double baseline_location_prior = stan::math::value_of(baseline_location_prior_arg);
    using T = stan::return_type_t<T_alpha, T_z>;
    using T_alpha_scalar = stan::value_type_t<T_alpha>;
    using T_z_scalar     = stan::value_type_t<T_z>;
    static constexpr const char* function_name = "nico_pinkney_C_and_log_J_from_alpha";

    const Eigen::Index n_cat = alpha.size();
    const Eigen::Index K     = n_cat - 1;
    const Eigen::Index n_C   = C_raw_c.size();
    if (K < 1) {
        throw std::domain_error(std::string(function_name) + ": alpha must have length >= 2");
    }
    stan::math::check_size_match(function_name, "Pinkney_w_c", Pinkney_w_c.size(), "C_raw_c", n_C);
    stan::math::check_size_match(function_name, "logit_start", logit_start.size(), "C_raw_c", n_C);
    stan::math::check_size_match(function_name, "C_raw_lower_c", C_raw_lower_c.size(), "C_raw_c", n_C);
    stan::math::check_size_match(function_name, "C_raw_upper_c", C_raw_upper_c.size(), "C_raw_c", n_C);
    stan::math::check_size_match(function_name, "k_of_each_C", k_of_each_C.size(), "C_raw_c", n_C);
    stan::math::check_size_match(function_name, "k_previous_of_each_C", k_previous_of_each_C.size(), "C_raw_c", n_C);

    if (n_C == 0) {
        Eigen::Matrix<T, -1, 1> empty_res(1);
        empty_res[0] = 0.0;
        return empty_res;
    }

    const auto& alpha_ref = stan::math::to_ref(alpha);
    const auto& z_ref     = stan::math::to_ref(C_raw_c);
    const auto& w_ref     = stan::math::to_ref(Pinkney_w_c);
    const auto& s_ref     = stan::math::to_ref(logit_start);
    const auto& l_ref     = stan::math::to_ref(C_raw_lower_c);
    const auto& u_ref     = stan::math::to_ref(C_raw_upper_c);

    constexpr bool is_var_call = stan::is_var<T>::value;
    constexpr bool is_alpha_var = stan::is_var<T_alpha_scalar>::value;

    // 1) Deduplicate unique (k, kp, lb, wd) cutpoint pairs in O(n_C) with zero heap allocs
    static constexpr int STACK_MAX_K = 128;
    int16_t pair_map_stack[STACK_MAX_K * (STACK_MAX_K + 1)];
    const bool use_pair_dedup = (K <= STACK_MAX_K);
    const int p_stride = static_cast<int>(K) + 1;
    if (use_pair_dedup) {
        std::fill_n(pair_map_stack, static_cast<int>(K) * p_stride, static_cast<int16_t>(-1));
    }

    const std::size_t n_cut_bufs = 6 * static_cast<std::size_t>(K) + static_cast<std::size_t>(n_cat);
    const std::size_t n_dbl_tot  = n_cut_bufs + 21 * static_cast<std::size_t>(n_C);
    const std::size_t n_int_tot  = 4 * static_cast<std::size_t>(n_C);
    std::vector<double> dbl_prim;
    std::vector<int> int_prim;
    std::vector<unsigned char> cf_prim;
    double* buf_base = nullptr;
    int* ibuf_base   = nullptr;
    unsigned char* u_clamp = nullptr;
    if constexpr (!is_var_call) {
        dbl_prim.resize(n_dbl_tot);
        int_prim.resize(n_int_tot);
        cf_prim.resize(n_C);
        buf_base  = dbl_prim.data();
        ibuf_base = int_prim.data();
        u_clamp   = cf_prim.data();
    } else {
        auto& arena = stan::math::ChainableStack::instance_->memalloc_;
        buf_base  = arena.alloc_array<double>(n_dbl_tot);
        ibuf_base = arena.alloc_array<int>(n_int_tot);
        u_clamp   = arena.alloc_array<unsigned char>(n_C);
    }

    double* P_prob     = buf_base;
    double* U_prob     = P_prob + K;
    double* centre_cut = U_prob + K;
    double* density    = centre_cut + K;
    double* var_fac    = density + K;
    double* tmp_cut    = var_fac + K;
    double* alpha_val  = tmp_cut + K;

    double* p_base     = alpha_val + n_cat;
    double* u_step     = p_base;
    double* u_vs       = u_step + n_C;
    double* u_em1      = u_vs + n_C;
    double* u_es       = u_em1 + n_C;
    double* u_place    = u_es + n_C;
    double* u_mu       = u_place + n_C;
    double* u_tau      = u_mu + n_C;
    double* u_log_tau  = u_tau + n_C;
    double* u_asig     = u_log_tau + n_C;
    double* u_bsig     = u_asig + n_C;
    double* u_spread   = u_bsig + n_C;
    double* u_lb       = u_spread + n_C;
    double* u_wd       = u_lb + n_C;

    double* c_base     = u_wd + n_C;
    double* theta_tmp  = c_base;
    double* abs_th_tmp = theta_tmp + n_C;
    double* log_D_tmp  = abs_th_tmp + n_C;
    double* r_tmp      = log_D_tmp + n_C;
    double* step_v_tmp = r_tmp + n_C;
    double* q_buf      = step_v_tmp + n_C;
    double* omq_buf    = q_buf + n_C;
    double* sr_buf     = omq_buf + n_C;

    int* pair_of_C = ibuf_base;
    int* is_inc_C  = pair_of_C + n_C;
    int* u_k       = is_inc_C + n_C;
    int* u_kp      = u_k + n_C;

    for (Eigen::Index j = 0; j < n_cat; ++j) {
        alpha_val[j] = stan::math::value_of(alpha_ref.coeff(j));
    }
    double running_S = 0.0;
    for (Eigen::Index k = 0; k < K; ++k) {
        running_S += alpha_val[k];
        P_prob[k] = running_S;
    }
    const double alpha_0   = running_S + alpha_val[K];
    const double inv_a0    = 1.0 / alpha_0;
    const double inv_a0_p1 = 1.0 / (alpha_0 + 1.0);
    double running_U = alpha_val[K];
    for (Eigen::Index k = K - 1; k >= 0; --k) {
        U_prob[k] = running_U * inv_a0;
        P_prob[k] *= inv_a0;
        if (k > 0) running_U += alpha_val[k];
    }

    if (use_probit_link == 1) {
        nico::math::copy_n(P_prob, static_cast<int>(K), centre_cut);
        nico::math::k_inv_Phi(centre_cut, static_cast<int>(K));
        nico::math::copy_n(U_prob, static_cast<int>(K), tmp_cut);
        nico::math::k_inv_Phi(tmp_cut, static_cast<int>(K));
        for (Eigen::Index k = 0; k < K; ++k) {
            if (P_prob[k] > 0.5) centre_cut[k] = -tmp_cut[k];
            density[k] = -0.5 * centre_cut[k] * centre_cut[k];
        }
        nico::math::k_exp(density, static_cast<int>(K));
        for (Eigen::Index k = 0; k < K; ++k) {
            density[k] *= nico::math::INV_SQRT_2PI;
        }
    } else {
        nico::math::copy_n(P_prob, static_cast<int>(K), centre_cut);
        nico::math::k_log(centre_cut, static_cast<int>(K));
        nico::math::copy_n(U_prob, static_cast<int>(K), tmp_cut);
        nico::math::k_log(tmp_cut, static_cast<int>(K));
        for (Eigen::Index k = 0; k < K; ++k) {
            centre_cut[k] -= tmp_cut[k];
            density[k] = P_prob[k] * U_prob[k];
        }
    }
    for (Eigen::Index k = 0; k < K; ++k) {
        var_fac[k] = (P_prob[k] * U_prob[k]) / (density[k] * density[k]);
    }

    int n_pairs = 0;
    double sum_log_width = 0.0;
    for (Eigen::Index p = 0; p < n_C; ++p) {
        const int k  = k_of_each_C[p] - 1;
        const int kp = k_previous_of_each_C[p] - 1;
        if (k < 0 || k >= K || (kp >= 0 && (kp >= k || p == 0)) || kp < -1) {
            throw std::domain_error(std::string(function_name) + ": invalid k_of_each_C or k_previous_of_each_C");
        }
        is_inc_C[p] = (kp >= 0) ? 1 : 0;
        const double lb_p = stan::math::value_of(l_ref.coeff(p));
        const double wd_p = stan::math::value_of(u_ref.coeff(p)) - lb_p;
        sum_log_width += std::log(wd_p);

        int u = -1;
        if (use_pair_dedup) {
            const int16_t cached = pair_map_stack[k * p_stride + (kp + 1)];
            if (cached >= 0 && u_lb[cached] == lb_p && u_wd[cached] == wd_p) {
                u = cached;
            }
        }
        if (u < 0) {
            u = n_pairs++;
            if (use_pair_dedup && pair_map_stack[k * p_stride + (kp + 1)] < 0) {
                pair_map_stack[k * p_stride + (kp + 1)] = static_cast<int16_t>(u);
            }
            u_k[u]  = k;
            u_kp[u] = kp;
            u_lb[u] = lb_p;
            u_wd[u] = wd_p;

            double centre_u;
            double spread_u;
            unsigned char flags = 0;
            if (kp < 0) {
                centre_u = centre_cut[k] + baseline_location_prior;
                spread_u = std::sqrt(var_fac[k] * inv_a0_p1);
                u_step[u] = 0.0;
                u_vs[u]   = 0.0;
                u_em1[u]  = 0.0;
                u_es[u]   = 0.0;
            } else {
                const double step   = centre_cut[k] - centre_cut[kp];
                const double inv_dp = 1.0 / (density[k] * density[kp]);
                double raw_vs = (var_fac[k] + var_fac[kp] - 2.0 * P_prob[kp] * U_prob[k] * inv_dp) * inv_a0_p1;
                if (raw_vs < 1e-14) {
                    raw_vs = 1e-14;
                    flags |= 1u;
                }
                const double sqrt_vs = std::sqrt(raw_vs);
                u_step[u] = step;
                u_vs[u]   = raw_vs;
                if (softplus == 1) {
                    const double es = std::exp(step);
                    double em1 = es - 1.0;
                    if (em1 < 1e-12) {
                        em1 = 1e-12;
                        flags |= 2u;
                    }
                    u_es[u]  = es;
                    u_em1[u] = em1;
                    centre_u = std::log(em1);
                    spread_u = (sqrt_vs * es) / em1;
                } else {
                    u_em1[u] = 0.0;
                    u_es[u]  = 0.0;
                    centre_u = std::log(step);
                    spread_u = sqrt_vs / step;
                }
            }
            u_spread[u] = spread_u;
            const double place_raw = (centre_u - lb_p) / wd_p;
            const double a_raw = 500.0 * (place_raw - 0.001);
            const double b_raw = 500.0 * (place_raw - 0.999);
            if (b_raw > 40.0) {
                flags |= 4u;
            } else if (a_raw < -40.0) {
                flags |= 8u;
            }
            u_asig[u]  = (flags & 4u) ? 40.0 : ((flags & 8u) ? -40.0 : a_raw);
            u_bsig[u]  = (flags & 4u) ? 40.0 : ((flags & 8u) ? -40.0 : b_raw);
            u_clamp[u] = flags;
        }
        pair_of_C[p] = u;
    }

    // Batch SIMD ONLY over unique pairs (reuse r_tmp / step_v_tmp as scratch -> zero heap allocs)
    double* l1pe_a = r_tmp;
    double* l1pe_b = step_v_tmp;
    nico::math::copy_n(u_asig, n_pairs, l1pe_a);
    nico::math::copy_n(u_bsig, n_pairs, l1pe_b);
    nico::math::k_log1p_exp(l1pe_a, n_pairs);
    nico::math::k_log1p_exp(l1pe_b, n_pairs);
    nico::math::k_inv_logit(u_asig, n_pairs);
    nico::math::k_inv_logit(u_bsig, n_pairs);

    for (int u = 0; u < n_pairs; ++u) {
        const unsigned char flags = u_clamp[u];
        double place_u;
        if (flags & 4u) {
            place_u = 0.999;
            u_asig[u] = 1.0;
            u_bsig[u] = 1.0;
        } else if (flags & 8u) {
            place_u = 0.001;
            u_asig[u] = 0.0;
            u_bsig[u] = 0.0;
        } else {
            place_u = 0.001 + (l1pe_a[u] - l1pe_b[u]) / 500.0;
        }
        u_place[u] = place_u;
        u_mu[u]    = place_u;
        const double A_u   = u_wd[u] * place_u * (1.0 - place_u);
        const double sp_u  = u_spread[u];
        const double tau_u = sp_u / std::sqrt(A_u * A_u + sp_u * sp_u);
        u_tau[u]     = tau_u;
        u_log_tau[u] = tau_u;
    }
    nico::math::k_logit(u_mu, n_pairs);
    nico::math::k_log(u_log_tau, n_pairs);

    // Pinkney partial-centering per element
    double sum_log_tau_over_D = 0.0;
    for (Eigen::Index p = 0; p < n_C; ++p) {
        const int u         = pair_of_C[p];
        const double w_p    = stan::math::value_of(w_ref.coeff(p));
        const double tau_u  = u_tau[u];
        const double D_p    = tau_u - w_p * (tau_u - 1.0);
        const double eta_sh = (1.0 - w_p) * stan::math::value_of(s_ref.coeff(p))
                              + stan::math::value_of(z_ref.coeff(p));
        const double th     = (w_p * u_mu[u] + tau_u * eta_sh) / D_p;
        theta_tmp[p]  = th;
        abs_th_tmp[p] = std::abs(th);
        log_D_tmp[p]  = D_p;
        sum_log_tau_over_D += u_log_tau[u];
    }
    nico::math::k_log(log_D_tmp, static_cast<int>(n_C));
    for (Eigen::Index p = 0; p < n_C; ++p) sum_log_tau_over_D -= log_D_tmp[p];

    for (Eigen::Index p = 0; p < n_C; ++p) {
        q_buf[p]   = theta_tmp[p];
        omq_buf[p] = -theta_tmp[p];
    }
    nico::math::k_inv_logit(q_buf, static_cast<int>(n_C));
    nico::math::k_inv_logit(omq_buf, static_cast<int>(n_C));
    nico::math::k_log_inv_logit(abs_th_tmp, static_cast<int>(n_C));

    double log_J = sum_log_width + sum_log_tau_over_D;
    for (Eigen::Index p = 0; p < n_C; ++p) {
        const int u = pair_of_C[p];
        r_tmp[p] = u_lb[u] + u_wd[u] * q_buf[p];
        log_J += 2.0 * abs_th_tmp[p] - std::abs(theta_tmp[p]);
        step_v_tmp[p] = r_tmp[p];
    }

    if (softplus == 1) {
        nico::math::k_log1p_exp(step_v_tmp, static_cast<int>(n_C));
        nico::math::copy_n(r_tmp, static_cast<int>(n_C), sr_buf);
        nico::math::k_inv_logit(sr_buf, static_cast<int>(n_C));
        nico::math::copy_n(r_tmp, static_cast<int>(n_C), abs_th_tmp);
        nico::math::k_log_inv_logit(abs_th_tmp, static_cast<int>(n_C));
        for (Eigen::Index p = 0; p < n_C; ++p) {
            if (is_inc_C[p] == 1) log_J += abs_th_tmp[p];
        }
    } else {
        nico::math::k_exp(step_v_tmp, static_cast<int>(n_C));
        for (Eigen::Index p = 0; p < n_C; ++p) {
            if (is_inc_C[p] == 1) log_J += r_tmp[p];
        }
    }

    if constexpr (!is_var_call) {
        Eigen::Matrix<double, -1, 1> result(n_C + 1);
        for (Eigen::Index p = 0; p < n_C; ++p) {
            result[p] = (is_inc_C[p] == 1) ? (result[p - 1] + step_v_tmp[p]) : r_tmp[p];
        }
        result[n_C] = log_J;
        return result;
    } else {
        stan::arena_t<Eigen::Matrix<T_alpha_scalar, -1, 1>> arena_alpha = alpha_ref;
        stan::arena_t<Eigen::Matrix<T_z_scalar, -1, 1>>     arena_z     = z_ref;
        stan::arena_t<Eigen::Matrix<stan::math::var, -1, 1>> arena_result(n_C + 1);
        stan::arena_t<Eigen::Matrix<double, 7, -1>> arena_buf(7, n_C);

        double prev_C = 0.0;
        for (Eigen::Index p = 0; p < n_C; ++p) {
            const int u = pair_of_C[p];
            const double C_p = (is_inc_C[p] == 1) ? (prev_C + step_v_tmp[p]) : r_tmp[p];
            prev_C = C_p;
            arena_result.coeffRef(p) = stan::math::var(new stan::math::vari(C_p, false));

            const double w_p   = stan::math::value_of(w_ref.coeff(p));
            const double tau_u = u_tau[u];
            const double D_p   = tau_u - w_p * (tau_u - 1.0);
            const double inv_D = 1.0 / D_p;
            const double eta_sh = (1.0 - w_p) * stan::math::value_of(s_ref.coeff(p))
                                  + stan::math::value_of(arena_z.coeff(p));
            arena_buf.coeffRef(0, p) = inv_D;
            arena_buf.coeffRef(1, p) = w_p;
            arena_buf.coeffRef(2, p) = tau_u;
            arena_buf.coeffRef(3, p) = eta_sh - (1.0 - w_p) * theta_tmp[p];
            arena_buf.coeffRef(4, p) = q_buf[p] - omq_buf[p];
            arena_buf.coeffRef(5, p) = u_wd[u] * q_buf[p] * omq_buf[p];
            arena_buf.coeffRef(6, p) = (softplus == 1) ? sr_buf[p] : step_v_tmp[p];
        }
        arena_result.coeffRef(n_C) = stan::math::var(new stan::math::vari(log_J, false));

        const bool softplus_flag = (softplus == 1);
        const bool probit_flag   = (use_probit_link == 1);
        stan::math::reverse_pass_callback(
            [arena_alpha, arena_z, arena_result, arena_buf,
             P_prob, U_prob, centre_cut, density, var_fac,
             u_step, u_vs, u_em1, u_es, u_place, u_tau,
             u_asig, u_bsig, u_spread, u_wd,
             pair_of_C, is_inc_C, u_k, u_kp, u_clamp,
             K, n_cat, n_C, n_pairs, inv_a0, inv_a0_p1, softplus_flag, probit_flag]() mutable {
                const double adj_J = arena_result.coeff(n_C).adj();
                double running_adj_C = 0.0;

                static constexpr int STACK_MAX_PAIRS = 4096;
                double stack_pair_adj[2 * STACK_MAX_PAIRS];
                std::vector<double> heap_pair_adj;
                double* adj_mu_u  = nullptr;
                double* adj_tau_u = nullptr;
                if constexpr (is_alpha_var) {
                    if (n_pairs <= STACK_MAX_PAIRS) {
                        adj_mu_u = stack_pair_adj;
                    } else {
                        heap_pair_adj.resize(2 * n_pairs);
                        adj_mu_u = heap_pair_adj.data();
                    }
                    adj_tau_u = adj_mu_u + n_pairs;
                    std::fill_n(adj_mu_u, 2 * n_pairs, 0.0);
                }

                for (Eigen::Index p = n_C - 1; p >= 0; --p) {
                    const double adj_C_here = arena_result.coeff(p).adj() + running_adj_C;
                    const bool is_inc = (is_inc_C[p] == 1);
                    running_adj_C = is_inc ? adj_C_here : 0.0;

                    double adj_r;
                    if (!is_inc) {
                        adj_r = adj_C_here;
                    } else if (!softplus_flag) {
                        adj_r = adj_C_here * arena_buf.coeff(6, p) + adj_J;
                    } else {
                        const double s_p = arena_buf.coeff(6, p);
                        adj_r = adj_C_here * s_p + adj_J * (1.0 - s_p);
                    }
                    const double adj_theta = adj_r * arena_buf.coeff(5, p) - adj_J * arena_buf.coeff(4, p);
                    const double inv_D_p   = arena_buf.coeff(0, p);
                    const double w_p       = arena_buf.coeff(1, p);
                    const double tau_u     = arena_buf.coeff(2, p);

                    if constexpr (stan::is_var<T_z_scalar>::value) {
                        arena_z.coeffRef(p).adj() += adj_theta * (tau_u * inv_D_p);
                    }
                    if constexpr (is_alpha_var) {
                        const int u = pair_of_C[p];
                        adj_mu_u[u]  += adj_theta * (w_p * inv_D_p);
                        adj_tau_u[u] += adj_theta * (arena_buf.coeff(3, p) * inv_D_p)
                                        + adj_J * (w_p * inv_D_p / tau_u);
                    }
                }

                if constexpr (is_alpha_var) {
                    double stack_cut_adj[5 * STACK_MAX_K];
                    std::vector<double> heap_cut_adj;
                    double* adj_centre = (K <= STACK_MAX_K) ? stack_cut_adj : (heap_cut_adj.resize(5 * K), heap_cut_adj.data());
                    std::fill_n(adj_centre, 5 * K, 0.0);
                    double* adj_vf   = adj_centre + K;
                    double* adj_P    = adj_vf + K;
                    double* adj_U    = adj_P + K;
                    double* adj_dens = adj_U + K;
                    double adj_alpha0 = 0.0;

                    for (int u = 0; u < n_pairs; ++u) {
                        const double amu  = adj_mu_u[u];
                        const double atau = adj_tau_u[u];
                        if (amu == 0.0 && atau == 0.0) continue;

                        const double wd_u    = u_wd[u];
                        const double place_u = u_place[u];
                        const double sp_u    = u_spread[u];
                        const double A_u     = wd_u * place_u * (1.0 - place_u);
                        const double denom2  = A_u * A_u + sp_u * sp_u;
                        const double denom   = std::sqrt(denom2);
                        const double inv_d3  = 1.0 / (denom2 * denom);

                        const double adj_place = amu / (place_u * (1.0 - place_u))
                            + atau * (-sp_u * A_u * inv_d3) * wd_u * (1.0 - 2.0 * place_u);
                        const double adj_spread = atau * (A_u * A_u * inv_d3);

                        const unsigned char flags = u_clamp[u];
                        const double d_clamp = ((flags & 12u) != 0) ? 0.0 : (u_asig[u] - u_bsig[u]);
                        const double adj_centre_u = adj_place * d_clamp / wd_u;

                        const int k  = u_k[u];
                        const int kp = u_kp[u];
                        if (kp < 0) {
                            adj_centre[k] += adj_centre_u;
                            adj_vf[k]     += adj_spread * 0.5 * inv_a0_p1 / sp_u;
                            adj_alpha0    -= adj_spread * 0.5 * sp_u * inv_a0_p1;
                        } else {
                            const double step = u_step[u];
                            const double vs   = u_vs[u];
                            const bool clamped_vs  = (flags & 1u) != 0;
                            const bool clamped_em1 = (flags & 2u) != 0;
                            double adj_step = 0.0;
                            double adj_vs   = 0.0;
                            if (!softplus_flag) {
                                adj_step = (adj_centre_u - adj_spread * sp_u) / step;
                                if (!clamped_vs) {
                                    adj_vs = adj_spread * 0.5 / (std::sqrt(vs) * step);
                                }
                            } else {
                                const double es  = u_es[u];
                                const double em1 = u_em1[u];
                                if (!clamped_em1) {
                                    adj_step = adj_centre_u * (es / em1) - adj_spread * (sp_u / em1);
                                } else {
                                    adj_step = adj_spread * sp_u;
                                }
                                if (!clamped_vs) {
                                    adj_vs = adj_spread * 0.5 * sp_u / vs;
                                }
                            }
                            adj_centre[k]  += adj_step;
                            adj_centre[kp] -= adj_step;
                            if (adj_vs != 0.0) {
                                const double s_vs = adj_vs * inv_a0_p1;
                                const double inv_dp = 1.0 / (density[k] * density[kp]);
                                const double cross  = s_vs * 2.0 * P_prob[kp] * U_prob[k] * inv_dp;
                                adj_alpha0   -= adj_vs * vs * inv_a0_p1;
                                adj_vf[k]    += s_vs;
                                adj_vf[kp]   += s_vs;
                                adj_P[kp]    -= s_vs * 2.0 * U_prob[k] * inv_dp;
                                adj_U[k]     -= s_vs * 2.0 * P_prob[kp] * inv_dp;
                                adj_dens[k]  += cross / density[k];
                                adj_dens[kp] += cross / density[kp];
                            }
                        }
                    }

                    double dot_PU = 0.0;
                    for (Eigen::Index k = 0; k < K; ++k) {
                        const double dk = density[k];
                        const double pk = P_prob[k];
                        const double uk = U_prob[k];
                        const double av = adj_vf[k];
                        const double inv_dk2 = 1.0 / (dk * dk);
                        adj_P[k]    += av * uk * inv_dk2;
                        adj_U[k]    += av * pk * inv_dk2;
                        adj_dens[k] -= av * 2.0 * var_fac[k] / dk;

                        if (probit_flag) {
                            const double ac = adj_centre[k] - adj_dens[k] * dk * centre_cut[k];
                            if (pk <= 0.5) {
                                adj_P[k] += ac / dk;
                            } else {
                                adj_U[k] -= ac / dk;
                            }
                        } else {
                            adj_P[k] += adj_dens[k] * uk + adj_centre[k] / pk;
                            adj_U[k] += adj_dens[k] * pk - adj_centre[k] / uk;
                        }
                        dot_PU += adj_P[k] * pk + adj_U[k] * uk;
                    }

                    double stack_suf[STACK_MAX_K + 2];
                    std::vector<double> heap_suf;
                    double* suf_P = (n_cat <= STACK_MAX_K + 1) ? stack_suf : (heap_suf.resize(n_cat), heap_suf.data());
                    suf_P[K] = 0.0;
                    for (Eigen::Index j = K - 1; j >= 0; --j) {
                        suf_P[j] = suf_P[j + 1] + adj_P[j];
                    }
                    double pre_U = 0.0;
                    for (Eigen::Index j = 0; j < n_cat; ++j) {
                        if (j > 0) pre_U += adj_U[j - 1];
                        const double da = (suf_P[j] + pre_U - dot_PU) * inv_a0 + adj_alpha0;
                        arena_alpha.coeffRef(j).adj() += da;
                    }
                }
            });
        return Eigen::Matrix<T, -1, 1>(arena_result);
    }
}


// ---- Precomputed-pair Pinkney C & total log-J directly from alpha (zero hash dedup, zero std::log(wd_p) during HMC) ----
namespace nico::math {

template <typename T_alpha,
          typename T_z,
          typename T_loc_prior,
          typename T_w,
          typename T_s,
          typename T_l,
          typename T_wd,
          typename T_ulb,
          typename T_uwd>
inline Eigen::Matrix<stan::return_type_t<T_alpha, T_z, T_loc_prior>, -1, 1>
nico_pinkney_C_and_log_J_from_alpha_precomputed(
        const T_alpha& alpha,
        const T_z& C_raw_c,
        const T_loc_prior& baseline_location_prior_arg,
        const T_w& Pinkney_w_c,
        const T_s& logit_start,
        const T_l& C_raw_lower_c,
        const T_wd& C_raw_width_c,
        const double sum_log_width_c,
        const std::vector<int>& pair_of_each_C,
        const std::vector<int>& is_increment_of_each_C,
        const std::vector<int>& u_k_in,
        const std::vector<int>& u_kp_in,
        const T_ulb& u_lb_in,
        const T_uwd& u_wd_in,
        const int use_probit_link,
        const int softplus,
        std::ostream* pstream__) {
    static constexpr const char* function_name = "nico_pinkney_C_and_log_J_from_alpha_precomputed";
    using T = stan::return_type_t<T_alpha, T_z, T_loc_prior>;
    using T_alpha_scalar = stan::value_type_t<T_alpha>;
    using T_z_scalar = stan::value_type_t<T_z>;
    (void) pstream__;

    const double baseline_location_prior = stan::math::value_of(baseline_location_prior_arg);
    const Eigen::Index n_cat = alpha.size();
    const Eigen::Index K = n_cat - 1;
    const Eigen::Index n_C = C_raw_c.size();
    const int n_pairs = static_cast<int>(u_k_in.size());
    if (K < 1) throw std::domain_error(std::string(function_name) + ": alpha must have length >= 2");
    if (n_C == 0) {
        Eigen::Matrix<T, -1, 1> empty_res(1);
        empty_res[0] = 0.0;
        return empty_res;
    }

    const auto& alpha_ref = stan::math::to_ref(alpha);
    const auto& z_ref     = stan::math::to_ref(C_raw_c);
    const auto& w_ref     = stan::math::to_ref(Pinkney_w_c);
    const auto& s_ref     = stan::math::to_ref(logit_start);
    const auto& ulb_ref   = stan::math::to_ref(u_lb_in);
    const auto& uwd_ref   = stan::math::to_ref(u_wd_in);
    (void) C_raw_lower_c;
    (void) C_raw_width_c;

    constexpr bool is_var_call = stan::is_var<T>::value;
    constexpr bool is_alpha_var = stan::is_var<T_alpha_scalar>::value;

    const std::size_t n_cut_bufs = 6 * static_cast<std::size_t>(K) + static_cast<std::size_t>(n_cat);
    const std::size_t n_dbl_tot  = n_cut_bufs
                                 + 13 * static_cast<std::size_t>(n_pairs)
                                 + 8 * static_cast<std::size_t>(n_C);
    const std::size_t n_int_tot  = 2 * static_cast<std::size_t>(n_C)
                                 + 2 * static_cast<std::size_t>(n_pairs);

    std::vector<double> dbl_prim;
    std::vector<int> int_prim;
    std::vector<unsigned char> cf_prim;
    double* buf_base = nullptr;
    int* ibuf_base   = nullptr;
    unsigned char* u_clamp = nullptr;
    if constexpr (!is_var_call) {
        dbl_prim.resize(n_dbl_tot);
        int_prim.resize(n_int_tot);
        cf_prim.resize(n_pairs);
        buf_base  = dbl_prim.data();
        ibuf_base = int_prim.data();
        u_clamp   = cf_prim.data();
    } else {
        auto& arena = stan::math::ChainableStack::instance_->memalloc_;
        buf_base  = arena.alloc_array<double>(n_dbl_tot);
        ibuf_base = arena.alloc_array<int>(n_int_tot);
        u_clamp   = arena.alloc_array<unsigned char>(n_pairs);
    }

    double* P_prob     = buf_base;
    double* U_prob     = P_prob + K;
    double* centre_cut = U_prob + K;
    double* density    = centre_cut + K;
    double* var_fac    = density + K;
    double* tmp_cut    = var_fac + K;
    double* alpha_val  = tmp_cut + K;

    double* u_step     = alpha_val + n_cat;
    double* u_vs       = u_step + n_pairs;
    double* u_em1      = u_vs + n_pairs;
    double* u_es       = u_em1 + n_pairs;
    double* u_place    = u_es + n_pairs;
    double* u_mu       = u_place + n_pairs;
    double* u_tau      = u_mu + n_pairs;
    double* u_log_tau  = u_tau + n_pairs;
    double* u_asig     = u_log_tau + n_pairs;
    double* u_bsig     = u_asig + n_pairs;
    double* u_spread   = u_bsig + n_pairs;
    double* u_lb       = u_spread + n_pairs;
    double* u_wd       = u_lb + n_pairs;

    double* theta_tmp  = u_wd + n_pairs;
    double* abs_th_tmp = theta_tmp + n_C;
    double* log_D_tmp  = abs_th_tmp + n_C;
    double* r_tmp      = log_D_tmp + n_C;
    double* step_v_tmp = r_tmp + n_C;
    double* q_buf      = step_v_tmp + n_C;
    double* omq_buf    = q_buf + n_C;
    double* sr_buf     = omq_buf + n_C;

    int* pair_of_C = ibuf_base;
    int* is_inc_C  = pair_of_C + n_C;
    int* u_k       = is_inc_C + n_C;
    int* u_kp      = u_k + n_pairs;

    for (Eigen::Index j = 0; j < n_cat; ++j) alpha_val[j] = stan::math::value_of(alpha_ref.coeff(j));
    double running_S = 0.0;
    for (Eigen::Index k = 0; k < K; ++k) {
        running_S += alpha_val[k];
        P_prob[k] = running_S;
    }
    const double alpha_0   = running_S + alpha_val[K];
    const double inv_a0    = 1.0 / alpha_0;
    const double inv_a0_p1 = 1.0 / (alpha_0 + 1.0);
    double running_U = alpha_val[K];
    for (Eigen::Index k = K - 1; k >= 0; --k) {
        U_prob[k] = running_U * inv_a0;
        P_prob[k] *= inv_a0;
        if (k > 0) running_U += alpha_val[k];
    }

    if (use_probit_link == 1) {
        nico::math::copy_n(P_prob, static_cast<int>(K), centre_cut);
        nico::math::k_inv_Phi(centre_cut, static_cast<int>(K));
        nico::math::copy_n(U_prob, static_cast<int>(K), tmp_cut);
        nico::math::k_inv_Phi(tmp_cut, static_cast<int>(K));
        for (Eigen::Index k = 0; k < K; ++k) {
            if (P_prob[k] > 0.5) centre_cut[k] = -tmp_cut[k];
            density[k] = -0.5 * centre_cut[k] * centre_cut[k];
        }
        nico::math::k_exp(density, static_cast<int>(K));
        for (Eigen::Index k = 0; k < K; ++k) density[k] *= nico::math::INV_SQRT_2PI;
    } else {
        nico::math::copy_n(P_prob, static_cast<int>(K), centre_cut);
        nico::math::k_log(centre_cut, static_cast<int>(K));
        nico::math::copy_n(U_prob, static_cast<int>(K), tmp_cut);
        nico::math::k_log(tmp_cut, static_cast<int>(K));
        for (Eigen::Index k = 0; k < K; ++k) {
            centre_cut[k] -= tmp_cut[k];
            density[k] = P_prob[k] * U_prob[k];
        }
    }
    for (Eigen::Index k = 0; k < K; ++k) {
        var_fac[k] = (P_prob[k] * U_prob[k]) / (density[k] * density[k]);
    }

    // Unique pairs only (n_pairs ~ O(150) vs n_C ~ O(1485))
    for (int u = 0; u < n_pairs; ++u) {
        const int k  = u_k_in[u] - 1;
        const int kp = u_kp_in[u] - 1;
        u_k[u]  = k;
        u_kp[u] = kp;
        const double lb_p = stan::math::value_of(ulb_ref.coeff(u));
        const double wd_p = stan::math::value_of(uwd_ref.coeff(u));
        u_lb[u] = lb_p;
        u_wd[u] = wd_p;

        double centre_u;
        double spread_u;
        unsigned char flags = 0;
        if (kp < 0) {
            centre_u = centre_cut[k] + baseline_location_prior;
            spread_u = std::sqrt(var_fac[k] * inv_a0_p1);
            u_step[u] = 0.0;
            u_vs[u]   = 0.0;
            u_em1[u]  = 0.0;
            u_es[u]   = 0.0;
        } else {
            const double step   = centre_cut[k] - centre_cut[kp];
            const double inv_dp = 1.0 / (density[k] * density[kp]);
            double raw_vs = (var_fac[k] + var_fac[kp] - 2.0 * P_prob[kp] * U_prob[k] * inv_dp) * inv_a0_p1;
            if (raw_vs < 1e-14) {
                raw_vs = 1e-14;
                flags |= 1u;
            }
            const double sqrt_vs = std::sqrt(raw_vs);
            u_step[u] = step;
            u_vs[u]   = raw_vs;
            if (softplus == 1) {
                const double es = std::exp(step);
                double em1 = es - 1.0;
                if (em1 < 1e-12) {
                    em1 = 1e-12;
                    flags |= 2u;
                }
                u_es[u]  = es;
                u_em1[u] = em1;
                centre_u = std::log(em1);
                spread_u = (sqrt_vs * es) / em1;
            } else {
                u_em1[u] = 0.0;
                u_es[u]  = 0.0;
                centre_u = std::log(step);
                spread_u = sqrt_vs / step;
            }
        }
        u_spread[u] = spread_u;
        const double place_raw = (centre_u - lb_p) / wd_p;
        const double a_raw = 500.0 * (place_raw - 0.001);
        const double b_raw = 500.0 * (place_raw - 0.999);
        if (b_raw > 40.0) flags |= 4u;
        else if (a_raw < -40.0) flags |= 8u;
        u_asig[u]  = (flags & 4u) ? 40.0 : ((flags & 8u) ? -40.0 : a_raw);
        u_bsig[u]  = (flags & 4u) ? 40.0 : ((flags & 8u) ? -40.0 : b_raw);
        u_clamp[u] = flags;
    }

    double* l1pe_a = r_tmp;
    double* l1pe_b = step_v_tmp;
    nico::math::copy_n(u_asig, n_pairs, l1pe_a);
    nico::math::copy_n(u_bsig, n_pairs, l1pe_b);
    nico::math::k_log1p_exp(l1pe_a, n_pairs);
    nico::math::k_log1p_exp(l1pe_b, n_pairs);
    nico::math::k_inv_logit(u_asig, n_pairs);
    nico::math::k_inv_logit(u_bsig, n_pairs);

    for (int u = 0; u < n_pairs; ++u) {
        const unsigned char flags = u_clamp[u];
        double place_u;
        if (flags & 4u) {
            place_u = 0.999;
            u_asig[u] = 1.0;
            u_bsig[u] = 1.0;
        } else if (flags & 8u) {
            place_u = 0.001;
            u_asig[u] = 0.0;
            u_bsig[u] = 0.0;
        } else {
            place_u = 0.001 + (l1pe_a[u] - l1pe_b[u]) / 500.0;
        }
        u_place[u] = place_u;
        u_mu[u]    = place_u;
        const double A_u   = u_wd[u] * place_u * (1.0 - place_u);
        const double sp_u  = u_spread[u];
        const double tau_u = sp_u / std::sqrt(A_u * A_u + sp_u * sp_u);
        u_tau[u]     = tau_u;
        u_log_tau[u] = tau_u;
    }
    nico::math::k_logit(u_mu, n_pairs);
    nico::math::k_log(u_log_tau, n_pairs);

    double sum_log_tau_over_D = 0.0;
    for (Eigen::Index p = 0; p < n_C; ++p) {
        const int u = pair_of_each_C[p] - 1;
        pair_of_C[p] = u;
        is_inc_C[p] = is_increment_of_each_C[p];
        const double w_p    = stan::math::value_of(w_ref.coeff(p));
        const double tau_u  = u_tau[u];
        const double D_p    = tau_u - w_p * (tau_u - 1.0);
        const double eta_sh = (1.0 - w_p) * stan::math::value_of(s_ref.coeff(p))
                              + stan::math::value_of(z_ref.coeff(p));
        const double th     = (w_p * u_mu[u] + tau_u * eta_sh) / D_p;
        theta_tmp[p]  = th;
        abs_th_tmp[p] = std::abs(th);
        log_D_tmp[p]  = D_p;
        sum_log_tau_over_D += u_log_tau[u];
    }
    nico::math::k_log(log_D_tmp, static_cast<int>(n_C));
    for (Eigen::Index p = 0; p < n_C; ++p) sum_log_tau_over_D -= log_D_tmp[p];

    for (Eigen::Index p = 0; p < n_C; ++p) {
        q_buf[p]   = theta_tmp[p];
        omq_buf[p] = -theta_tmp[p];
    }
    nico::math::k_inv_logit(q_buf, static_cast<int>(n_C));
    nico::math::k_inv_logit(omq_buf, static_cast<int>(n_C));
    nico::math::k_log_inv_logit(abs_th_tmp, static_cast<int>(n_C));

    double log_J = sum_log_width_c + sum_log_tau_over_D;
    for (Eigen::Index p = 0; p < n_C; ++p) {
        const int u = pair_of_C[p];
        r_tmp[p] = u_lb[u] + u_wd[u] * q_buf[p];
        log_J += 2.0 * abs_th_tmp[p] - std::abs(theta_tmp[p]);
        step_v_tmp[p] = r_tmp[p];
    }

    if (softplus == 1) {
        nico::math::k_log1p_exp(step_v_tmp, static_cast<int>(n_C));
        nico::math::copy_n(r_tmp, static_cast<int>(n_C), sr_buf);
        nico::math::k_inv_logit(sr_buf, static_cast<int>(n_C));
        nico::math::copy_n(r_tmp, static_cast<int>(n_C), abs_th_tmp);
        nico::math::k_log_inv_logit(abs_th_tmp, static_cast<int>(n_C));
        for (Eigen::Index p = 0; p < n_C; ++p) {
            if (is_inc_C[p] == 1) log_J += abs_th_tmp[p];
        }
    } else {
        nico::math::k_exp(step_v_tmp, static_cast<int>(n_C));
        for (Eigen::Index p = 0; p < n_C; ++p) {
            if (is_inc_C[p] == 1) log_J += r_tmp[p];
        }
    }

    if constexpr (!is_var_call) {
        Eigen::Matrix<double, -1, 1> result(n_C + 1);
        for (Eigen::Index p = 0; p < n_C; ++p) {
            result[p] = (is_inc_C[p] == 1) ? (result[p - 1] + step_v_tmp[p]) : r_tmp[p];
        }
        result[n_C] = log_J;
        return result;
    } else {
        stan::arena_t<Eigen::Matrix<T_alpha_scalar, -1, 1>> arena_alpha = alpha_ref;
        stan::arena_t<Eigen::Matrix<T_z_scalar, -1, 1>>     arena_z     = z_ref;
        stan::arena_t<Eigen::Matrix<stan::math::var, -1, 1>> arena_result(n_C + 1);
        stan::arena_t<Eigen::Matrix<double, 7, -1>> arena_buf(7, n_C);

        double prev_C = 0.0;
        for (Eigen::Index p = 0; p < n_C; ++p) {
            const int u = pair_of_C[p];
            const double C_p = (is_inc_C[p] == 1) ? (prev_C + step_v_tmp[p]) : r_tmp[p];
            prev_C = C_p;
            arena_result.coeffRef(p) = stan::math::var(new stan::math::vari(C_p, false));

            const double w_p   = stan::math::value_of(w_ref.coeff(p));
            const double tau_u = u_tau[u];
            const double D_p   = tau_u - w_p * (tau_u - 1.0);
            const double inv_D = 1.0 / D_p;
            const double eta_sh = (1.0 - w_p) * stan::math::value_of(s_ref.coeff(p))
                                  + stan::math::value_of(arena_z.coeff(p));
            arena_buf.coeffRef(0, p) = inv_D;
            arena_buf.coeffRef(1, p) = w_p;
            arena_buf.coeffRef(2, p) = tau_u;
            arena_buf.coeffRef(3, p) = eta_sh - (1.0 - w_p) * theta_tmp[p];
            arena_buf.coeffRef(4, p) = q_buf[p] - omq_buf[p];
            arena_buf.coeffRef(5, p) = u_wd[u] * q_buf[p] * omq_buf[p];
            arena_buf.coeffRef(6, p) = (softplus == 1) ? sr_buf[p] : step_v_tmp[p];
        }
        arena_result.coeffRef(n_C) = stan::math::var(new stan::math::vari(log_J, false));

        const bool softplus_flag = (softplus == 1);
        const bool probit_flag   = (use_probit_link == 1);
        stan::math::reverse_pass_callback(
            [arena_alpha, arena_z, arena_result, arena_buf,
             P_prob, U_prob, centre_cut, density, var_fac,
             u_step, u_vs, u_em1, u_es, u_place, u_tau,
             u_asig, u_bsig, u_spread, u_wd,
             pair_of_C, is_inc_C, u_k, u_kp, u_clamp,
             K, n_cat, n_C, n_pairs, inv_a0, inv_a0_p1, softplus_flag, probit_flag]() mutable {
                const double adj_J = arena_result.coeff(n_C).adj();
                double running_adj_C = 0.0;

                static constexpr int STACK_MAX_PAIRS = 4096;
                static constexpr int STACK_MAX_K = 128;
                double stack_pair_adj[2 * STACK_MAX_PAIRS];
                std::vector<double> heap_pair_adj;
                double* adj_mu_u  = nullptr;
                double* adj_tau_u = nullptr;
                if constexpr (is_alpha_var) {
                    if (n_pairs <= STACK_MAX_PAIRS) {
                        adj_mu_u = stack_pair_adj;
                    } else {
                        heap_pair_adj.resize(2 * n_pairs);
                        adj_mu_u = heap_pair_adj.data();
                    }
                    adj_tau_u = adj_mu_u + n_pairs;
                    std::fill_n(adj_mu_u, 2 * n_pairs, 0.0);
                }

                for (Eigen::Index p = n_C - 1; p >= 0; --p) {
                    const double adj_C_here = arena_result.coeff(p).adj() + running_adj_C;
                    const bool is_inc = (is_inc_C[p] == 1);
                    running_adj_C = is_inc ? adj_C_here : 0.0;

                    double adj_r;
                    if (!is_inc) {
                        adj_r = adj_C_here;
                    } else if (!softplus_flag) {
                        adj_r = adj_C_here * arena_buf.coeff(6, p) + adj_J;
                    } else {
                        const double s_p = arena_buf.coeff(6, p);
                        adj_r = adj_C_here * s_p + adj_J * (1.0 - s_p);
                    }
                    const double adj_theta = adj_r * arena_buf.coeff(5, p) - adj_J * arena_buf.coeff(4, p);
                    const double inv_D_p   = arena_buf.coeff(0, p);
                    const double w_p       = arena_buf.coeff(1, p);
                    const double tau_u     = arena_buf.coeff(2, p);

                    if constexpr (stan::is_var<T_z_scalar>::value) {
                        arena_z.coeffRef(p).adj() += adj_theta * (tau_u * inv_D_p);
                    }
                    if constexpr (is_alpha_var) {
                        const int u = pair_of_C[p];
                        adj_mu_u[u]  += adj_theta * (w_p * inv_D_p);
                        adj_tau_u[u] += adj_theta * (arena_buf.coeff(3, p) * inv_D_p)
                                        + adj_J * (w_p * inv_D_p / tau_u);
                    }
                }

                if constexpr (is_alpha_var) {
                    double stack_cut_adj[5 * STACK_MAX_K];
                    std::vector<double> heap_cut_adj;
                    double* adj_centre = (K <= STACK_MAX_K) ? stack_cut_adj : (heap_cut_adj.resize(5 * K), heap_cut_adj.data());
                    std::fill_n(adj_centre, 5 * K, 0.0);
                    double* adj_vf   = adj_centre + K;
                    double* adj_P    = adj_vf + K;
                    double* adj_U    = adj_P + K;
                    double* adj_dens = adj_U + K;
                    double adj_alpha0 = 0.0;

                    for (int u = 0; u < n_pairs; ++u) {
                        const double amu  = adj_mu_u[u];
                        const double atau = adj_tau_u[u];
                        if (amu == 0.0 && atau == 0.0) continue;

                        const double wd_u    = u_wd[u];
                        const double place_u = u_place[u];
                        const double sp_u    = u_spread[u];
                        const double A_u     = wd_u * place_u * (1.0 - place_u);
                        const double denom2  = A_u * A_u + sp_u * sp_u;
                        const double denom   = std::sqrt(denom2);
                        const double inv_d3  = 1.0 / (denom2 * denom);

                        const double adj_place = amu / (place_u * (1.0 - place_u))
                            + atau * (-sp_u * A_u * inv_d3) * wd_u * (1.0 - 2.0 * place_u);
                        const double adj_spread = atau * (A_u * A_u * inv_d3);

                        const unsigned char flags = u_clamp[u];
                        const double d_clamp = ((flags & 12u) != 0) ? 0.0 : (u_asig[u] - u_bsig[u]);
                        const double adj_centre_u = adj_place * d_clamp / wd_u;

                        const int k  = u_k[u];
                        const int kp = u_kp[u];
                        if (kp < 0) {
                            adj_centre[k] += adj_centre_u;
                            adj_vf[k]     += adj_spread * 0.5 * inv_a0_p1 / sp_u;
                            adj_alpha0    -= adj_spread * 0.5 * sp_u * inv_a0_p1;
                        } else {
                            const double step = u_step[u];
                            const double vs   = u_vs[u];
                            const bool clamped_vs  = (flags & 1u) != 0;
                            const bool clamped_em1 = (flags & 2u) != 0;
                            double adj_step = 0.0;
                            double adj_vs   = 0.0;
                            if (!softplus_flag) {
                                adj_step = (adj_centre_u - adj_spread * sp_u) / step;
                                if (!clamped_vs) adj_vs = adj_spread * 0.5 / (std::sqrt(vs) * step);
                            } else {
                                const double es  = u_es[u];
                                const double em1 = u_em1[u];
                                if (!clamped_em1) {
                                    adj_step = adj_centre_u * (es / em1) - adj_spread * (sp_u / em1);
                                } else {
                                    adj_step = adj_spread * sp_u;
                                }
                                if (!clamped_vs) adj_vs = adj_spread * 0.5 * sp_u / vs;
                            }
                            adj_centre[k]  += adj_step;
                            adj_centre[kp] -= adj_step;
                            if (adj_vs != 0.0) {
                                const double s_vs = adj_vs * inv_a0_p1;
                                const double inv_dp = 1.0 / (density[k] * density[kp]);
                                const double cross  = s_vs * 2.0 * P_prob[kp] * U_prob[k] * inv_dp;
                                adj_alpha0   -= adj_vs * vs * inv_a0_p1;
                                adj_vf[k]    += s_vs;
                                adj_vf[kp]   += s_vs;
                                adj_P[kp]    -= s_vs * 2.0 * U_prob[k] * inv_dp;
                                adj_U[k]     -= s_vs * 2.0 * P_prob[kp] * inv_dp;
                                adj_dens[k]  += cross / density[k];
                                adj_dens[kp] += cross / density[kp];
                            }
                        }
                    }

                    double dot_PU = 0.0;
                    for (Eigen::Index k = 0; k < K; ++k) {
                        const double dk = density[k];
                        const double pk = P_prob[k];
                        const double uk = U_prob[k];
                        const double av = adj_vf[k];
                        const double inv_dk2 = 1.0 / (dk * dk);
                        adj_P[k]    += av * uk * inv_dk2;
                        adj_U[k]    += av * pk * inv_dk2;
                        adj_dens[k] -= av * 2.0 * var_fac[k] / dk;

                        if (probit_flag) {
                            const double ac = adj_centre[k] - adj_dens[k] * dk * centre_cut[k];
                            if (pk <= 0.5) adj_P[k] += ac / dk;
                            else           adj_U[k] -= ac / dk;
                        } else {
                            adj_P[k] += adj_dens[k] * uk + adj_centre[k] / pk;
                            adj_U[k] += adj_dens[k] * pk - adj_centre[k] / uk;
                        }
                        dot_PU += adj_P[k] * pk + adj_U[k] * uk;
                    }

                    double stack_suf[STACK_MAX_K + 2];
                    std::vector<double> heap_suf;
                    double* suf_P = (n_cat <= STACK_MAX_K + 1) ? stack_suf : (heap_suf.resize(n_cat), heap_suf.data());
                    suf_P[K] = 0.0;
                    for (Eigen::Index j = K - 1; j >= 0; --j) {
                        suf_P[j] = suf_P[j + 1] + adj_P[j];
                    }
                    double pre_U = 0.0;
                    for (Eigen::Index j = 0; j < n_cat; ++j) {
                        if (j > 0) pre_U += adj_U[j - 1];
                        const double da = (suf_P[j] + pre_U - dot_PU) * inv_a0 + adj_alpha0;
                        arena_alpha.coeffRef(j).adj() += da;
                    }
                }
            });
        return Eigen::Matrix<T, -1, 1>(arena_result);
    }
}

}  // namespace nico::math

using nico::math::nico_pinkney_C_and_log_J_from_alpha_precomputed;


// ---- Dual-class vectorized precomputed-pair Pinkney C & log-J + Fused O_biv_RC priors node ----
#if __has_include("stan_external_dirichlet_cutpoint_median.hpp")
#include "stan_external_dirichlet_cutpoint_median.hpp"
#include "stan_external_induced_dirichlet_with_Se_Sp_prior.hpp"
#include "stan_external_lkj_corr_cholesky_lpdf.hpp"
#include "stan_external_set0_marginal_at_chosen_cutoffs.hpp"
#elif __has_include("NicoSIMD/stan_external_dirichlet_cutpoint_median.hpp")
#include "NicoSIMD/stan_external_dirichlet_cutpoint_median.hpp"
#include "NicoSIMD/stan_external_induced_dirichlet_with_Se_Sp_prior.hpp"
#include "NicoSIMD/stan_external_lkj_corr_cholesky_lpdf.hpp"
#include "NicoSIMD/stan_external_set0_marginal_at_chosen_cutoffs.hpp"
#endif

namespace nico::math {
namespace pinkney_fast_detail {

template <int n_cls,
          typename T_alpha,
          typename T_z,
          typename T_w,
          typename T_s,
          typename T_ulb,
          typename T_uwd>
inline Eigen::Matrix<stan::return_type_t<T_alpha, T_z>, -1, 1>
pinkney_from_alpha_precomputed_impl(
        const T_alpha* const* alpha_ptrs,
        const T_z* const* z_ptrs,
        const double* blp_cls,
        const T_w* const* w_ptrs,
        const T_s* const* s_ptrs,
        const double* sum_log_wd_cls,
        const std::vector<int>& pair_of_each_C,
        const std::vector<int>& is_increment_of_each_C,
        const std::vector<int>& u_k_in,
        const std::vector<int>& u_kp_in,
        const T_ulb& u_lb_in,
        const T_uwd& u_wd_in,
        const int use_probit_link,
        const int softplus) {
    static constexpr const char* fn_name = "nico_pinkney_C_and_log_J_from_alpha_precomputed";
    using T = stan::return_type_t<T_alpha, T_z>;
    using T_alpha_scalar = stan::value_type_t<T_alpha>;
    using T_z_scalar = stan::value_type_t<T_z>;

    const Eigen::Index n_cat = alpha_ptrs[0]->size();
    const Eigen::Index K = n_cat - 1;
    const Eigen::Index n_C_one = z_ptrs[0]->size();
    const Eigen::Index n_C_tot = static_cast<Eigen::Index>(n_cls) * n_C_one;
    const int n_pairs_one = static_cast<int>(u_k_in.size());
    const int n_pairs_tot = n_cls * n_pairs_one;
    const int K_tot = n_cls * static_cast<int>(K);

    if (K < 1) throw std::domain_error(std::string(fn_name) + ": alpha must have length >= 2");
    if (n_C_one == 0) {
        Eigen::Matrix<T, -1, 1> empty_res(n_cls);
        for (int c = 0; c < n_cls; ++c) empty_res[c] = 0.0;
        return empty_res;
    }

    const auto& ulb_ref = stan::math::to_ref(u_lb_in);
    const auto& uwd_ref = stan::math::to_ref(u_wd_in);

    constexpr bool is_var_call = stan::is_var<T>::value;
    constexpr bool is_alpha_var = stan::is_var<T_alpha_scalar>::value;

    const std::size_t rev_scratch_sz = (is_var_call && is_alpha_var)
        ? static_cast<std::size_t>(2 * n_pairs_one + 5 * static_cast<int>(K) + static_cast<int>(n_cat))
        : 0U;
    const std::size_t n_cut_bufs = 9 * static_cast<std::size_t>(K_tot) + static_cast<std::size_t>(n_cls * n_cat);
    const std::size_t n_pair_bufs = 17 * static_cast<std::size_t>(n_pairs_tot) + 3 * static_cast<std::size_t>(n_pairs_one);
    const std::size_t n_C_bufs = 12 * static_cast<std::size_t>(n_C_tot);
    const std::size_t n_dbl_tot = n_cut_bufs + n_pair_bufs + n_C_bufs + rev_scratch_sz;
    const std::size_t n_int_tot = 2 * static_cast<std::size_t>(n_C_one) + 2 * static_cast<std::size_t>(n_pairs_one);

    std::vector<double> dbl_prim;
    std::vector<int> int_prim;
    std::vector<unsigned char> cf_prim;
    double* buf_base = nullptr;
    int* ibuf_base = nullptr;
    unsigned char* u_clamp = nullptr;
    if constexpr (!is_var_call) {
        dbl_prim.resize(n_dbl_tot);
        int_prim.resize(n_int_tot);
        cf_prim.resize(n_pairs_tot);
        buf_base = dbl_prim.data();
        ibuf_base = int_prim.data();
        u_clamp = cf_prim.data();
    } else {
        auto& arena = stan::math::ChainableStack::instance_->memalloc_;
        buf_base = arena.alloc_array<double>(n_dbl_tot);
        ibuf_base = arena.alloc_array<int>(n_int_tot);
        u_clamp = arena.alloc_array<unsigned char>(n_pairs_tot);
    }

    // Cutpoint buffers (contiguous across classes for SIMD / ivdep batching)
    double* P_prob = buf_base;                              // K_tot
    double* U_prob = P_prob + K_tot;                        // K_tot
    double* centre_cut = U_prob + K_tot;                    // K_tot
    double* density = centre_cut + K_tot;                   // K_tot
    double* inv_dd = density + K_tot;                       // K_tot
    double* inv_P = inv_dd + K_tot;                         // K_tot
    double* inv_U = inv_P + K_tot;                          // K_tot
    double* var_fac = inv_U + K_tot;                        // K_tot
    double* tmp_cut = var_fac + K_tot;                      // K_tot
    double* alpha_val = tmp_cut + K_tot;                    // n_cls * n_cat

    // Pair buffers
    double* u_step = alpha_val + n_cls * n_cat;
    double* u_vs = u_step + n_pairs_tot;
    double* u_em1 = u_vs + n_pairs_tot;
    double* u_es = u_em1 + n_pairs_tot;
    double* u_place = u_es + n_pairs_tot;
    double* u_mu = u_place + n_pairs_tot;
    double* u_tau = u_mu + n_pairs_tot;
    double* u_log_tau = u_tau + n_pairs_tot;
    double* u_asig = u_log_tau + n_pairs_tot;
    double* u_bsig = u_asig + n_pairs_tot;
    double* u_spread = u_bsig + n_pairs_tot;
    double* u_inv_tau = u_spread + n_pairs_tot;
    double* u_inv_pom = u_inv_tau + n_pairs_tot;
    double* u_inv_d3 = u_inv_pom + n_pairs_tot;
    double* u_inv_sp = u_inv_d3 + n_pairs_tot;
    double* u_inv_step = u_inv_sp + n_pairs_tot;
    double* u_lb_one = u_inv_step + n_pairs_tot;
    double* u_wd_one = u_lb_one + n_pairs_one;
    double* u_inv_wd_one = u_wd_one + n_pairs_one;

    // C buffers (12 * n_C_tot)
    double* theta_tmp = u_inv_wd_one + n_pairs_one;
    double* abs_th_tmp = theta_tmp + n_C_tot;
    double* log_D_tmp = abs_th_tmp + n_C_tot;
    double* inv_D_tmp = log_D_tmp + n_C_tot;
    double* r_tmp = inv_D_tmp + n_C_tot;
    double* step_v_tmp = r_tmp + n_C_tot;
    double* q_buf = step_v_tmp + n_C_tot;
    double* dth_dz = q_buf + n_C_tot;
    double* dth_dmu = dth_dz + n_C_tot;
    double* dth_dtau_a = dth_dmu + n_C_tot;
    double* dj_dtau_a = dth_dtau_a + n_C_tot;
    double* uwd_qomq_a = dj_dtau_a + n_C_tot;
    double* rev_work = uwd_qomq_a + n_C_tot;
    double* l1pe_pair_a = dth_dz;

    int* pair_of_C = ibuf_base;
    int* is_inc_C = pair_of_C + n_C_one;
    int* u_k = is_inc_C + n_C_one;
    int* u_kp = u_k + n_pairs_one;

    // 1) Read alpha and compute P_prob, U_prob for all classes
    double inv_a0[2], inv_a0_p1[2];
    for (int c = 0; c < n_cls; ++c) {
        const auto& aref = *alpha_ptrs[c];
        double* av = alpha_val + c * n_cat;
        double* Pc = P_prob + c * K;
        double* Uc = U_prob + c * K;
        double runS = 0.0;
        for (Eigen::Index j = 0; j < n_cat; ++j) av[j] = stan::math::value_of(aref.coeff(j));
        for (Eigen::Index k = 0; k < K; ++k) {
            runS += av[k];
            Pc[k] = runS;
        }
        const double a0 = runS + av[K];
        inv_a0[c] = 1.0 / a0;
        inv_a0_p1[c] = 1.0 / (a0 + 1.0);
        double runU = av[K];
        for (Eigen::Index k = K - 1; k >= 0; --k) {
            Uc[k] = runU * inv_a0[c];
            Pc[k] *= inv_a0[c];
            if (k > 0) runU += av[k];
        }
    }

    if (use_probit_link == 1) {
        nico::math::copy_n(P_prob, K_tot, centre_cut);
        nico::math::copy_n(U_prob, K_tot, tmp_cut);
        nico::math::k_inv_Phi(centre_cut, K_tot);
        nico::math::k_inv_Phi(tmp_cut, K_tot);
        for (int k = 0; k < K_tot; ++k) {
            if (P_prob[k] > 0.5) centre_cut[k] = -tmp_cut[k];
            density[k] = -0.5 * centre_cut[k] * centre_cut[k];
        }
        nico::math::k_exp(density, K_tot);
        for (int k = 0; k < K_tot; ++k) density[k] *= nico::math::INV_SQRT_2PI;
    } else {
        nico::math::copy_n(P_prob, K_tot, centre_cut);
        nico::math::copy_n(U_prob, K_tot, tmp_cut);
        nico::math::k_log(centre_cut, K_tot);
        nico::math::k_log(tmp_cut, K_tot);
        for (int k = 0; k < K_tot; ++k) {
            centre_cut[k] -= tmp_cut[k];
            density[k] = P_prob[k] * U_prob[k];
        }
        if constexpr (is_var_call && is_alpha_var) {
            for (int k = 0; k < K_tot; ++k) {
                inv_P[k] = 1.0 / P_prob[k];
                inv_U[k] = 1.0 / U_prob[k];
            }
        }
    }
    for (int k = 0; k < K_tot; ++k) {
        const double id = 1.0 / density[k];
        inv_dd[k] = id;
        var_fac[k] = (P_prob[k] * U_prob[k]) * (id * id);
    }

    // 2) Vectorized unique-pair computations across all classes (zero raw intrinsics)
    for (int u = 0; u < n_pairs_one; ++u) {
        u_k[u] = u_k_in[u] - 1;
        u_kp[u] = u_kp_in[u] - 1;
        const double lbp = stan::math::value_of(ulb_ref.coeff(u));
        const double wdp = stan::math::value_of(uwd_ref.coeff(u));
        u_lb_one[u] = lbp;
        u_wd_one[u] = wdp;
        u_inv_wd_one[u] = 1.0 / wdp;
    }

    double* sqrt_work = u_asig;
    double* log_work = u_place;
    for (int c = 0; c < n_cls; ++c) {
        const double* cc = centre_cut + c * K;
        const double* pp = P_prob + c * K;
        const double* up = U_prob + c * K;
        const double* vf = var_fac + c * K;
        const double* idd = inv_dd + c * K;
        const double ia0p1 = inv_a0_p1[c];
        for (int u = 0; u < n_pairs_one; ++u) {
            const int gu = c * n_pairs_one + u;
            const int k = u_k[u];
            const int kp = u_kp[u];
            if (kp < 0) {
                u_step[gu] = 0.0;
                const double vs0 = vf[k] * ia0p1;
                u_vs[gu] = vs0;
                sqrt_work[gu] = vs0;
                u_clamp[gu] = 0;
                log_work[gu] = 1.0;
            } else {
                const double step = cc[k] - cc[kp];
                const double inv_dp = idd[k] * idd[kp];
                double raw_vs = (vf[k] + vf[kp] - 2.0 * pp[kp] * up[k] * inv_dp) * ia0p1;
                unsigned char fl = 0;
                if (raw_vs < 1e-14) { raw_vs = 1e-14; fl |= 1u; }
                u_step[gu] = step;
                u_vs[gu] = raw_vs;
                sqrt_work[gu] = raw_vs;
                u_clamp[gu] = fl;
                log_work[gu] = step;
            }
        }
    }
    nico::math::k_sqrt(sqrt_work, n_pairs_tot);
    if (softplus == 1) {
        nico::math::copy_n(u_step, n_pairs_tot, u_em1);
        nico::math::k_expm1(u_em1, n_pairs_tot);
        for (int c = 0; c < n_cls; ++c) {
            for (int u = 0; u < n_pairs_one; ++u) {
                const int gu = c * n_pairs_one + u;
                if (u_kp[u] >= 0) {
                    const double es = u_em1[gu] + 1.0;
                    u_es[gu] = es;
                    if (u_em1[gu] < 1e-12) { u_em1[gu] = 1e-12; u_clamp[gu] |= 2u; }
                    log_work[gu] = u_em1[gu];
                } else {
                    u_es[gu] = 0.0;
                    u_em1[gu] = 1.0;
                }
            }
        }
    }
    nico::math::k_log(log_work, n_pairs_tot);
    for (int c = 0; c < n_cls; ++c) {
        const double* cc = centre_cut + c * K;
        const double blp = blp_cls[c];
        for (int u = 0; u < n_pairs_one; ++u) {
            const int gu = c * n_pairs_one + u;
            const int k = u_k[u];
            const int kp = u_kp[u];
            const double sq_vs = sqrt_work[gu];
            const double sp_denom = (kp < 0) ? 1.0 : ((softplus == 1) ? u_em1[gu] : u_step[gu]);
            const double inv_sd = 1.0 / sp_denom;
            u_inv_step[gu] = inv_sd;
            const double centre_u = (kp < 0) ? (cc[k] + blp) : log_work[gu];
            const double spread_u = (kp < 0) ? sq_vs : ((softplus == 1) ? (sq_vs * u_es[gu] * inv_sd) : (sq_vs * inv_sd));
            u_spread[gu] = spread_u;
            u_inv_sp[gu] = 1.0 / spread_u;
            const double place_raw = (centre_u - u_lb_one[u]) * u_inv_wd_one[u];
            const double a_raw = 500.0 * (place_raw - 0.001);
            const double b_raw = 500.0 * (place_raw - 0.999);
            unsigned char fl = u_clamp[gu];
            if (b_raw > 40.0) fl |= 4u;
            else if (a_raw < -40.0) fl |= 8u;
            u_asig[gu] = (fl & 4u) ? 40.0 : ((fl & 8u) ? -40.0 : a_raw);
            u_bsig[gu] = (fl & 4u) ? 40.0 : ((fl & 8u) ? -40.0 : b_raw);
            u_clamp[gu] = fl;
        }
    }

    nico::math::copy_n(u_asig, 2 * n_pairs_tot, l1pe_pair_a);
    nico::math::k_log1p_exp(l1pe_pair_a, 2 * n_pairs_tot);
    nico::math::k_inv_logit(u_asig, 2 * n_pairs_tot);
    const double* l1pe_pair_b = l1pe_pair_a + n_pairs_tot;

    for (int c = 0; c < n_cls; ++c) {
        for (int u = 0; u < n_pairs_one; ++u) {
            const int gu = c * n_pairs_one + u;
            const unsigned char fl = u_clamp[gu];
            double place_u;
            if (fl & 4u) {
                place_u = 0.999; u_asig[gu] = 1.0; u_bsig[gu] = 1.0;
            } else if (fl & 8u) {
                place_u = 0.001; u_asig[gu] = 0.0; u_bsig[gu] = 0.0;
            } else {
                place_u = 0.001 + (l1pe_pair_a[gu] - l1pe_pair_b[gu]) * 0.002;
            }
            u_place[gu] = place_u;
            u_mu[gu] = place_u;
            const double A_u = u_wd_one[u] * place_u * (1.0 - place_u);
            const double sp_u = u_spread[gu];
            u_inv_d3[gu] = A_u * A_u + sp_u * sp_u;
            u_inv_pom[gu] = 1.0 / (place_u * (1.0 - place_u));
        }
    }
    nico::math::k_logit(u_mu, n_pairs_tot);
    nico::math::k_sqrt(u_inv_d3, n_pairs_tot);
    for (int gu = 0; gu < n_pairs_tot; ++gu) {
        const double sp_u = u_spread[gu];
        const double denom = u_inv_d3[gu];
        const double inv_den = 1.0 / denom;
        const double tau_u = sp_u * inv_den;
        u_tau[gu] = tau_u;
        u_inv_tau[gu] = denom * u_inv_sp[gu];
        u_log_tau[gu] = tau_u;
        u_inv_d3[gu] = inv_den / (denom * denom);
    }
    nico::math::k_log(u_log_tau, n_pairs_tot);

    // 3) Evaluate all C cutpoints across classes in vectorized batches
    double sum_log_tau_over_D[2] = {0.0, 0.0};
    for (Eigen::Index p = 0; p < n_C_one; ++p) {
        pair_of_C[p] = pair_of_each_C[p] - 1;
        is_inc_C[p] = is_increment_of_each_C[p];
    }
    for (int c = 0; c < n_cls; ++c) {
        const auto& w_ref = *w_ptrs[c];
        const auto& s_ref = *s_ptrs[c];
        const auto& z_ref = *z_ptrs[c];
        const int u_off = c * n_pairs_one;
        const Eigen::Index p_off = static_cast<Eigen::Index>(c) * n_C_one;
        double slt = 0.0;
        for (Eigen::Index p = 0; p < n_C_one; ++p) {
            const Eigen::Index gp = p_off + p;
            const int gu = u_off + pair_of_C[p];
            const double w_p = stan::math::value_of(w_ref.coeff(p));
            const double tau_u = u_tau[gu];
            const double D_p = tau_u - w_p * (tau_u - 1.0);
            const double eta_sh = (1.0 - w_p) * stan::math::value_of(s_ref.coeff(p))
                                + stan::math::value_of(z_ref.coeff(p));
            log_D_tmp[gp] = D_p;
            r_tmp[gp] = w_p * u_mu[gu] + tau_u * eta_sh;
            theta_tmp[gp] = w_p;
            step_v_tmp[gp] = eta_sh;
            slt += u_log_tau[gu];
        }
        sum_log_tau_over_D[c] = slt;
    }
    for (Eigen::Index gp = 0; gp < n_C_tot; ++gp) {
        inv_D_tmp[gp] = 1.0 / log_D_tmp[gp];
    }
    for (int c = 0; c < n_cls; ++c) {
        const int u_off = c * n_pairs_one;
        const Eigen::Index p_off = static_cast<Eigen::Index>(c) * n_C_one;
        for (Eigen::Index p = 0; p < n_C_one; ++p) {
            const Eigen::Index gp = p_off + p;
            const int gu = u_off + pair_of_C[p];
            const double inv_D = inv_D_tmp[gp];
            const double th = r_tmp[gp] * inv_D;
            const double w_p = theta_tmp[gp];
            const double eta_sh = step_v_tmp[gp];
            const double tau_u = u_tau[gu];
            theta_tmp[gp] = th;
            abs_th_tmp[gp] = std::abs(th);
            q_buf[gp] = th;
            if constexpr (is_var_call) {
                dth_dz[gp] = tau_u * inv_D;
                dth_dmu[gp] = w_p * inv_D;
                dth_dtau_a[gp] = (eta_sh - (1.0 - w_p) * th) * inv_D;
                dj_dtau_a[gp] = w_p * inv_D * u_inv_tau[gu];
            }
        }
    }

    nico::math::k_log(log_D_tmp, static_cast<int>(n_C_tot));
    nico::math::k_inv_logit(q_buf, static_cast<int>(n_C_tot));
    nico::math::k_log_inv_logit(abs_th_tmp, static_cast<int>(n_C_tot));

    double log_J_cls[2] = {0.0, 0.0};
    for (int c = 0; c < n_cls; ++c) {
        const Eigen::Index p_off = static_cast<Eigen::Index>(c) * n_C_one;
        double lj = sum_log_wd_cls[c] + sum_log_tau_over_D[c];
        for (Eigen::Index p = 0; p < n_C_one; ++p) {
            const Eigen::Index gp = p_off + p;
            const int u = pair_of_C[p];
            const double q = q_buf[gp];
            const double rv = u_lb_one[u] + u_wd_one[u] * q;
            r_tmp[gp] = rv;
            step_v_tmp[gp] = rv;
            lj += 2.0 * abs_th_tmp[gp] - std::abs(theta_tmp[gp]) - log_D_tmp[gp];
            if constexpr (is_var_call) {
                const double omq = 1.0 - q;
                uwd_qomq_a[gp] = u_wd_one[u] * q * omq;
                theta_tmp[gp] = q - omq;
            }
        }
        log_J_cls[c] = lj;
    }

    if (softplus == 1) {
        nico::math::k_log1p_exp(step_v_tmp, static_cast<int>(n_C_tot));
        nico::math::copy_n(r_tmp, static_cast<int>(n_C_tot), q_buf);
        nico::math::k_inv_logit(q_buf, static_cast<int>(n_C_tot));
    } else {
        nico::math::k_exp(step_v_tmp, static_cast<int>(n_C_tot));
        if constexpr (is_var_call) nico::math::copy_n(step_v_tmp, static_cast<int>(n_C_tot), q_buf);
    }

    if constexpr (!is_var_call) {
        Eigen::Matrix<double, -1, 1> result(n_C_tot + n_cls);
        for (int c = 0; c < n_cls; ++c) {
            const Eigen::Index p_off = static_cast<Eigen::Index>(c) * n_C_one;
            double lj_extra = 0.0;
            for (Eigen::Index p = 0; p < n_C_one; ++p) {
                const Eigen::Index gp = p_off + p;
                const int inc = is_inc_C[p];
                result[gp] = (inc == 1) ? (result[gp - 1] + step_v_tmp[gp]) : r_tmp[gp];
                if (inc == 1) lj_extra += (softplus == 1) ? (r_tmp[gp] - step_v_tmp[gp]) : r_tmp[gp];
            }
            result[n_C_tot + c] = log_J_cls[c] + lj_extra;
        }
        return result;
    } else {
        stan::arena_t<Eigen::Matrix<T_alpha_scalar, -1, 1>> arena_a0, arena_a1;
        stan::arena_t<Eigen::Matrix<T_z_scalar, -1, 1>> arena_z0, arena_z1;
        if constexpr (is_alpha_var) {
            arena_a0 = stan::math::to_ref(*alpha_ptrs[0]);
            if constexpr (n_cls == 2) arena_a1 = stan::math::to_ref(*alpha_ptrs[1]);
        }
        if constexpr (stan::is_var<T_z_scalar>::value) {
            arena_z0 = stan::math::to_ref(*z_ptrs[0]);
            if constexpr (n_cls == 2) arena_z1 = stan::math::to_ref(*z_ptrs[1]);
        }

        stan::arena_t<Eigen::Matrix<stan::math::var, -1, 1>> arena_result(n_C_tot + n_cls);
        for (int c = 0; c < n_cls; ++c) {
            const Eigen::Index p_off = static_cast<Eigen::Index>(c) * n_C_one;
            double prev_C = 0.0;
            double lj_extra = 0.0;
            for (Eigen::Index p = 0; p < n_C_one; ++p) {
                const Eigen::Index gp = p_off + p;
                const int inc = is_inc_C[p];
                const double C_p = (inc == 1) ? (prev_C + step_v_tmp[gp]) : r_tmp[gp];
                prev_C = C_p;
                arena_result.coeffRef(gp) = stan::math::var(new stan::math::vari(C_p, false));
                if (inc == 1) lj_extra += (softplus == 1) ? (r_tmp[gp] - step_v_tmp[gp]) : r_tmp[gp];
            }
            arena_result.coeffRef(n_C_tot + c) = stan::math::var(new stan::math::vari(log_J_cls[c] + lj_extra, false));
        }

        const bool softplus_flag = (softplus == 1);
        const bool probit_flag = (use_probit_link == 1);
        const double inv_a0_0 = inv_a0[0], inv_a0_p1_0 = inv_a0_p1[0];
        const double inv_a0_1 = (n_cls == 2) ? inv_a0[1] : 0.0;
        const double inv_a0_p1_1 = (n_cls == 2) ? inv_a0_p1[1] : 0.0;

        stan::math::reverse_pass_callback(
            [arena_a0, arena_a1, arena_z0, arena_z1, arena_result,
             P_prob, U_prob, centre_cut, density, inv_dd, inv_P, inv_U, var_fac,
             u_vs, u_es, u_place, u_asig, u_bsig, u_spread, u_inv_pom, u_inv_d3,
             u_inv_sp, u_inv_step, u_wd_one, u_inv_wd_one,
             dth_dz, dth_dmu, dth_dtau_a, dj_dtau_a, uwd_qomq_a, theta_tmp, q_buf,
             pair_of_C, is_inc_C, u_k, u_kp, u_clamp, rev_work,
             K, n_cat, n_C_one, n_C_tot, n_pairs_one,
             inv_a0_0, inv_a0_p1_0, inv_a0_1, inv_a0_p1_1,
             softplus_flag, probit_flag]() mutable {
                for (int c = 0; c < n_cls; ++c) {
                    const double adj_J = arena_result.coeff(n_C_tot + c).adj();
                    const Eigen::Index p_off = static_cast<Eigen::Index>(c) * n_C_one;
                    const int u_off = c * n_pairs_one;

                    double* adj_mu_u = rev_work;
                    double* adj_tau_u = adj_mu_u + n_pairs_one;
                    if constexpr (is_alpha_var) {
                        std::fill_n(adj_mu_u, 2 * n_pairs_one, 0.0);
                    }

                    auto& az = (c == 0) ? arena_z0 : arena_z1;
                    double running_adj_C = 0.0;
                    for (Eigen::Index p = n_C_one - 1; p >= 0; --p) {
                        const Eigen::Index gp = p_off + p;
                        const double adj_C_here = arena_result.coeff(gp).adj() + running_adj_C;
                        const bool is_inc = (is_inc_C[p] == 1);
                        running_adj_C = is_inc ? adj_C_here : 0.0;

                        double adj_r;
                        if (!is_inc) {
                            adj_r = adj_C_here;
                        } else if (!softplus_flag) {
                            adj_r = adj_C_here * q_buf[gp] + adj_J;
                        } else {
                            const double s_p = q_buf[gp];
                            adj_r = adj_C_here * s_p + adj_J * (1.0 - s_p);
                        }
                        const double adj_theta = adj_r * uwd_qomq_a[gp] - adj_J * theta_tmp[gp];
                        if constexpr (stan::is_var<T_z_scalar>::value) {
                            az.coeffRef(p).adj() += adj_theta * dth_dz[gp];
                        }
                        if constexpr (is_alpha_var) {
                            const int u = pair_of_C[p];
                            adj_mu_u[u] += adj_theta * dth_dmu[gp];
                            adj_tau_u[u] += adj_theta * dth_dtau_a[gp] + adj_J * dj_dtau_a[gp];
                        }
                    }

                    if constexpr (is_alpha_var) {
                        const double inv_a0 = (c == 0) ? inv_a0_0 : inv_a0_1;
                        const double inv_a0_p1 = (c == 0) ? inv_a0_p1_0 : inv_a0_p1_1;
                        const double* Pc = P_prob + c * K;
                        const double* Uc = U_prob + c * K;
                        const double* cc = centre_cut + c * K;
                        const double* dd = density + c * K;
                        const double* idd = inv_dd + c * K;
                        const double* vf = var_fac + c * K;

                        double* adj_centre = adj_tau_u + n_pairs_one;
                        std::fill_n(adj_centre, 5 * K, 0.0);
                        double* adj_vf = adj_centre + K;
                        double* adj_P = adj_vf + K;
                        double* adj_U = adj_P + K;
                        double* adj_dens = adj_U + K;
                        double adj_alpha0 = 0.0;

                        for (int u = 0; u < n_pairs_one; ++u) {
                            const double amu = adj_mu_u[u];
                            const double atau = adj_tau_u[u];
                            if (amu == 0.0 && atau == 0.0) continue;
                            const int gu = u_off + u;
                            const double wd_u = u_wd_one[u];
                            const double inv_wd_u = u_inv_wd_one[u];
                            const double place_u = u_place[gu];
                            const double sp_u = u_spread[gu];
                            const double inv_sp_u = u_inv_sp[gu];
                            const double inv_step_u = u_inv_step[gu];
                            const double A_u = wd_u * place_u * (1.0 - place_u);
                            const double inv_d3 = u_inv_d3[gu];

                            const double adj_place = amu * u_inv_pom[gu]
                                + atau * (-sp_u * A_u * inv_d3) * wd_u * (1.0 - 2.0 * place_u);
                            const double adj_spread = atau * (A_u * A_u * inv_d3);

                            const unsigned char fl = u_clamp[gu];
                            const double d_clamp = ((fl & 12u) != 0) ? 0.0 : (u_asig[gu] - u_bsig[gu]);
                            const double adj_centre_u = adj_place * d_clamp * inv_wd_u;

                            const int k = u_k[u];
                            const int kp = u_kp[u];
                            if (kp < 0) {
                                adj_centre[k] += adj_centre_u;
                                adj_vf[k] += adj_spread * 0.5 * inv_a0_p1 * inv_sp_u;
                                adj_alpha0 -= adj_spread * 0.5 * sp_u * inv_a0_p1;
                            } else {
                                const double vs = u_vs[gu];
                                const bool clamped_vs = (fl & 1u) != 0;
                                const bool clamped_em1 = (fl & 2u) != 0;
                                double adj_step = 0.0;
                                double adj_vs = 0.0;
                                if (!softplus_flag) {
                                    adj_step = (adj_centre_u - adj_spread * sp_u) * inv_step_u;
                                    if (!clamped_vs) adj_vs = adj_spread * 0.5 * inv_sp_u * inv_step_u * inv_step_u;
                                } else {
                                    const double es = u_es[gu];
                                    if (!clamped_em1) adj_step = (adj_centre_u * es - adj_spread * sp_u) * inv_step_u;
                                    else              adj_step = adj_spread * sp_u;
                                    if (!clamped_vs) adj_vs = adj_spread * 0.5 * sp_u / vs;
                                }
                                adj_centre[k] += adj_step;
                                adj_centre[kp] -= adj_step;
                                if (adj_vs != 0.0) {
                                    const double s_vs = adj_vs * inv_a0_p1;
                                    const double inv_dp = idd[k] * idd[kp];
                                    const double cross = s_vs * 2.0 * Pc[kp] * Uc[k] * inv_dp;
                                    adj_alpha0 -= adj_vs * vs * inv_a0_p1;
                                    adj_vf[k] += s_vs;
                                    adj_vf[kp] += s_vs;
                                    adj_P[kp] -= s_vs * 2.0 * Uc[k] * inv_dp;
                                    adj_U[k] -= s_vs * 2.0 * Pc[kp] * inv_dp;
                                    adj_dens[k] += cross * idd[k];
                                    adj_dens[kp] += cross * idd[kp];
                                }
                            }
                        }

                        double dot_PU = 0.0;
                        const double* iPc = probit_flag ? nullptr : (inv_P + c * K);
                        const double* iUc = probit_flag ? nullptr : (inv_U + c * K);
                        for (Eigen::Index k = 0; k < K; ++k) {
                            const double dk = dd[k];
                            const double idk = idd[k];
                            const double pk = Pc[k];
                            const double uk = Uc[k];
                            const double av = adj_vf[k];
                            const double inv_dk2 = idk * idk;
                            adj_P[k] += av * uk * inv_dk2;
                            adj_U[k] += av * pk * inv_dk2;
                            adj_dens[k] -= av * 2.0 * vf[k] * idk;

                            if (probit_flag) {
                                const double ac = adj_centre[k] - adj_dens[k] * dk * cc[k];
                                if (pk <= 0.5) adj_P[k] += ac * idk;
                                else           adj_U[k] -= ac * idk;
                            } else {
                                adj_P[k] += adj_dens[k] * uk + adj_centre[k] * iPc[k];
                                adj_U[k] += adj_dens[k] * pk - adj_centre[k] * iUc[k];
                            }
                            dot_PU += adj_P[k] * pk + adj_U[k] * uk;
                        }

                        double* suf_P = adj_dens + K;
                        suf_P[K] = 0.0;
                        for (Eigen::Index j = K - 1; j >= 0; --j) suf_P[j] = suf_P[j + 1] + adj_P[j];
                        double pre_U = 0.0;
                        auto& aa = (c == 0) ? arena_a0 : arena_a1;
                        for (Eigen::Index j = 0; j < n_cat; ++j) {
                            if (j > 0) pre_U += adj_U[j - 1];
                            aa.coeffRef(j).adj() += (suf_P[j] + pre_U - dot_PU) * inv_a0 + adj_alpha0;
                        }
                    }
                }
            });
        return Eigen::Matrix<T, -1, 1>(arena_result);
    }
}

}  // namespace pinkney_fast_detail

template <typename T_alpha, typename T_z, typename T_loc_prior, typename T_w, typename T_s,
          typename T_l, typename T_wd, typename T_swd, typename T_ulb, typename T_uwd>
inline Eigen::Matrix<stan::return_type_t<T_alpha, T_z, T_loc_prior>, -1, 1>
nico_pinkney_C_and_log_J_both_from_alpha_precomputed(
        const std::vector<Eigen::Matrix<T_alpha, -1, 1>>& alpha,
        const std::vector<Eigen::Matrix<T_z, -1, 1>>& C_raw,
        const T_loc_prior& baseline_location_prior,
        const std::vector<Eigen::Matrix<T_w, -1, 1>>& Pinkney_w,
        const std::vector<Eigen::Matrix<T_s, -1, 1>>& logit_start,
        const std::vector<Eigen::Matrix<T_l, -1, 1>>& C_raw_lower,
        const std::vector<Eigen::Matrix<T_wd, -1, 1>>& C_raw_width,
        const T_swd& sum_log_width,
        const std::vector<int>& pair_of_each_C,
        const std::vector<int>& is_increment_of_each_C,
        const std::vector<int>& u_k_in,
        const std::vector<int>& u_kp_in,
        const T_ulb& u_lb_in,
        const T_uwd& u_wd_in,
        const int use_probit_link,
        const int softplus,
        std::ostream* /*pstream__*/) {
    (void) C_raw_lower; (void) C_raw_width;
    const auto& blp_ref = stan::math::to_ref(baseline_location_prior);
    const auto& swd_ref = stan::math::to_ref(sum_log_width);
    const Eigen::Matrix<T_alpha, -1, 1>* a_p[2] = {&alpha[0], &alpha[1]};
    const Eigen::Matrix<T_z, -1, 1>* z_p[2] = {&C_raw[0], &C_raw[1]};
    const Eigen::Matrix<T_w, -1, 1>* w_p[2] = {&Pinkney_w[0], &Pinkney_w[1]};
    const Eigen::Matrix<T_s, -1, 1>* s_p[2] = {&logit_start[0], &logit_start[1]};
    const double blp[2] = {stan::math::value_of(blp_ref.coeff(0)), stan::math::value_of(blp_ref.coeff(1))};
    const double swd[2] = {stan::math::value_of(swd_ref.coeff(0)), stan::math::value_of(swd_ref.coeff(1))};
    return pinkney_fast_detail::pinkney_from_alpha_precomputed_impl<2>(
        a_p, z_p, blp, w_p, s_p, swd, pair_of_each_C, is_increment_of_each_C,
        u_k_in, u_kp_in, u_lb_in, u_wd_in, use_probit_link, softplus);
}


// ---- Single-node fused priors + Dirichlet-Se/Sp-Set0 prior block for ordered bivariate random-cutpoint models ----
template <bool propto,
          typename T_phi,
          typename T_bmu,
          typename T_bsd,
          typename T_bl,
          typename T_bz,
          typename T_alpha,
          typename T_lk,
          typename T_kap,
          typename T_tlj,
          typename T_bmu_m,
          typename T_bmu_s,
          typename T_bsd_m,
          typename T_bsd_s,
          typename T_lkj,
          typename T_km,
          typename T_ks,
          typename T_kd,
          typename T_pda,
          typename T_seg,
          typename T_A,
          typename T_beta,
          typename T_lms,
          typename T_blp,
          typename T_v0,
          typename T_bcnd,
          typename T_bcd>
inline stan::return_type_t<T_phi, T_bmu, T_bsd, T_bl, T_bz, T_alpha, T_lk, T_kap, T_tlj>
nico_ord_biv_rc_priors_lpdf(
        const std::vector<Eigen::Matrix<T_phi, -1, 1>>& dirichlet_phi,
        const T_bmu& beta_mu,
        const T_bsd& beta_SD,
        const T_bl& beta_L_Omega,
        const std::vector<Eigen::Matrix<T_bz, 1, -1>>& beta_z,
        const std::vector<Eigen::Matrix<T_alpha, -1, 1>>& alpha,
        const T_lk& log_kappa,
        const T_kap& kappa,
        const T_tlj& total_log_J_of_C,
        const std::vector<int>& n_cov,
        const T_bmu_m& prior_beta_mu_mean,
        const T_bmu_s& prior_beta_mu_SD,
        const T_bsd_m& prior_beta_SD_mean,
        const T_bsd_s& prior_beta_SD_SD,
        const T_lkj& prior_beta_corr_LKJ,
        const T_kd& prior_kappa_dof,
        const T_km& prior_kappa_mean,
        const T_ks& prior_kappa_SD,
        const std::vector<Eigen::Matrix<T_pda, -1, 1>>& prior_dirichlet_alpha,
        const std::vector<Eigen::Matrix<T_seg, -1, 1>>& alpha_of_the_segments,
        const std::vector<std::vector<int>>& chosen_cutoffs,
        const std::vector<int>& n_chosen_cutoffs,
        const int use_Se_Sp_prior,
        const int Se_Sp_prior_type,
        const std::vector<Eigen::Matrix<T_A, -1, 1>>& prior_Se_Sp_dirichlet_A,
        const std::vector<Eigen::Matrix<T_beta, -1, -1>>& prior_beta_td,
        const std::vector<Eigen::Matrix<T_lms, -1, -1>>& prior_link_mean_SD_td,
        const T_blp& baseline_location_prior_td,
        const T_v0& v0_baseline_prior_td,
        const T_bcnd& baseline_case_nd,
        const T_bcd& baseline_case_d,
        const int use_probit_link,
        std::ostream* /*pstream__*/) {
    namespace cmf = nico::math::cutpoint_median_fast_v2;
    namespace s0d = nico::math::set0_marginal_detail;
    namespace prior = nico::math::induced_dirichlet_prior;
    namespace Se_Sp = nico::math::Se_Sp_prior_at_chosen_cutoffs;
    using stan::math::var;
    static constexpr const char* fn_name = "nico_ord_biv_rc_priors_lpdf";
    static constexpr double PI_VAL = 3.1415926535897932384626433832795;

    const int n_cat = static_cast<int>(dirichlet_phi[0].size());
    const int n_thr = n_cat - 1;
    const int n_studies = static_cast<int>(beta_z.size());
    const int n_cov_nd = n_cov[0];
    const int n_cov_d = n_cov[1];

    constexpr bool phi_var = !stan::is_constant_all<T_phi>::value;
    constexpr bool bmu_var = !stan::is_constant_all<T_bmu>::value;
    constexpr bool bsd_var = !stan::is_constant_all<T_bsd>::value;
    constexpr bool bl_var = !stan::is_constant_all<T_bl>::value;
    constexpr bool bz_var = !stan::is_constant_all<T_bz>::value;
    constexpr bool alpha_var = !stan::is_constant_all<T_alpha>::value;
    constexpr bool lk_var = !stan::is_constant_all<T_lk>::value;
    constexpr bool kap_var = !stan::is_constant_all<T_kap>::value;
    constexpr bool tlj_var = !stan::is_constant_all<T_tlj>::value;
    constexpr bool any_var = phi_var || bmu_var || bsd_var || bl_var || bz_var
                          || alpha_var || lk_var || kap_var || tlj_var;

    const auto& bmu_ref = stan::math::to_ref(beta_mu);
    const auto& bsd_ref = stan::math::to_ref(beta_SD);
    const auto& bl_ref = stan::math::to_ref(beta_L_Omega);
    const auto& lk_ref = stan::math::to_ref(log_kappa);
    const auto& kap_ref = stan::math::to_ref(kappa);
    const auto& tlj_ref = stan::math::to_ref(total_log_J_of_C);
    const auto& bmu_m_ref = stan::math::to_ref(prior_beta_mu_mean);
    const auto& bmu_s_ref = stan::math::to_ref(prior_beta_mu_SD);
    const auto& bsd_m_ref = stan::math::to_ref(prior_beta_SD_mean);
    const auto& bsd_s_ref = stan::math::to_ref(prior_beta_SD_SD);
    const auto& kd_ref = stan::math::to_ref(prior_kappa_dof);
    const auto& km_ref = stan::math::to_ref(prior_kappa_mean);
    const auto& ks_ref = stan::math::to_ref(prior_kappa_SD);
    const auto& blp_ref = stan::math::to_ref(baseline_location_prior_td);
    const auto& v0_ref = stan::math::to_ref(v0_baseline_prior_td);
    const auto& bcnd_ref = stan::math::to_ref(baseline_case_nd);
    const auto& bcd_ref = stan::math::to_ref(baseline_case_d);

    double total = 0.0;

    // 1) Allocate compact double scratch on stack (or arena for partials kept for reverse)
    const int n_cov_tot = n_cov_nd + n_cov_d;
    const int mc[2] = {n_chosen_cutoffs[0], n_chosen_cutoffs[1]};
    const int M_tot = (use_Se_Sp_prior == 1) ? (mc[0] + mc[1])
                    : ((use_Se_Sp_prior == 2)
                       ? (((mc[0] > 0) ? n_thr : 0) + ((mc[1] > 0) ? n_thr : 0))
                       : 0);
    const int m_max = std::max(1, std::max(mc[0], mc[1]));
    const int bz_cnt = 2 * n_studies;
    const int max_n = std::max(M_tot, 1);

    // Reverse adjoints buffer allocated cleanly once when any_var
    const std::size_t adj_sz =
        static_cast<std::size_t>(n_cov_tot + 2 + 1 + bz_cnt + 2 + 2 + 2 * n_cat + 2 * n_cat);
    const std::size_t scratch_sz =
        static_cast<std::size_t>(2 * (n_cat + 1) + 8 * std::max(1, M_tot) + 20 * std::max(1, M_tot) + 20
                                 + 6 * std::max(1, n_cat) + s0d::raw_scratch_size(m_max)
                                 + 2 * Se_Sp::buffers_of_the_Se_Sp_prior::n_values(std::max(m_max, n_thr))
                                 + 8 * std::max(n_cat, 1));
    std::vector<double> prim_mem;
    double* rev_mem = nullptr;
    double* work_mem = nullptr;
    if constexpr (!any_var) {
        prim_mem.resize(adj_sz + scratch_sz);
        rev_mem = prim_mem.data();
        work_mem = rev_mem + adj_sz;
    } else {
        auto& arena = stan::math::ChainableStack::instance_->memalloc_;
        rev_mem = arena.alloc_array<double>(adj_sz);
        work_mem = arena.alloc_array<double>(scratch_sz);
        std::fill_n(rev_mem, adj_sz, 0.0);
    }

    double* d_bmu = rev_mem;
    double* d_bsd = d_bmu + n_cov_tot;
    double* d_bl11 = d_bsd + 2;
    double* d_bz = d_bl11 + 1;
    double* d_lk = d_bz + bz_cnt;
    double* d_kap = d_lk + 2;
    double* d_phi = d_kap + 2;       // 2 * n_cat
    double* d_alpha = d_phi + 2 * n_cat; // 2 * n_cat

    // 1a) beta_mu normal priors + Xbeta_baseline[0..1]
    double Xb_base[2] = {0.0, 0.0};
    for (int c = 0; c < 2; ++c) {
        const int nc = (c == 0) ? n_cov_nd : n_cov_d;
        const int off = (c == 0) ? 0 : n_cov_nd;
        for (int j = 0; j < nc; ++j) {
            const double val = stan::math::value_of(bmu_ref.coeff(c, j));
            const double mu = stan::math::value_of(bmu_m_ref.coeff(c, j));
            const double sd = stan::math::value_of(bmu_s_ref.coeff(c, j));
            if (!(sd > 0.0) || !std::isfinite(sd)) {
                stan::math::check_positive_finite(fn_name, "Scale parameter", sd);
            }
            const double inv_sd = 1.0 / sd;
            const double z = (val - mu) * inv_sd;
            total += -0.5 * z * z;
            if constexpr (!propto) total -= std::log(sd) + s0d::HALF_LOG_TWO_PI;
            d_bmu[off + j] = -z * inv_sd;
            const double bc = (c == 0) ? stan::math::value_of(bcnd_ref.coeff(j))
                                       : stan::math::value_of(bcd_ref.coeff(j));
            Xb_base[c] += bc * val;
        }
    }

    // 1b) beta_SD normal priors
    for (int c = 0; c < 2; ++c) {
        const double val = stan::math::value_of(bsd_ref.coeff(c));
        const double mu = stan::math::value_of(bsd_m_ref.coeff(c));
        const double sd = stan::math::value_of(bsd_s_ref.coeff(c));
        if (!(sd > 0.0) || !std::isfinite(sd)) {
            stan::math::check_positive_finite(fn_name, "Scale parameter", sd);
        }
        const double inv_sd = 1.0 / sd;
        const double z = (val - mu) * inv_sd;
        total += -0.5 * z * z;
        if constexpr (!propto) total -= std::log(sd) + s0d::HALF_LOG_TWO_PI;
        d_bsd[c] = -z * inv_sd;
    }

    // 1c) 2x2 LKJ correlation Cholesky
    {
        const double eta = stan::math::value_of(prior_beta_corr_LKJ);
        const double L00 = stan::math::value_of(bl_ref.coeff(0, 0));
        const double L01 = stan::math::value_of(bl_ref.coeff(0, 1));
        const double L11 = stan::math::value_of(bl_ref.coeff(1, 1));
        if (!(eta > 0.0) || !std::isfinite(eta)) {
            stan::math::check_positive_finite(fn_name, "Shape parameter", eta);
        }
        if (!(L00 > 0.0) || !(L11 > 0.0) || (L01 != 0.0)) {
            Eigen::Matrix<double, 2, 2> bad_L;
            bad_L << L00, L01, stan::math::value_of(bl_ref.coeff(1, 0)), L11;
            stan::math::check_cholesky_factor(fn_name, "Random variable", bad_L);
        }
        if constexpr (!propto) total += nico_lkj_constant_from_data(eta, 2U, nullptr);
        if (eta != 1.0) {
            const double coef = 2.0 * eta - 2.0;
            total += coef * std::log(L11);
            *d_bl11 = coef / L11;
        }
    }

    // 1d) beta_z standard Normal prior directly on array[n_studies] row_vector[2] beta_z
    {
        double sq_sum = 0.0;
        for (int s = 0; s < n_studies; ++s) {
            const double z0 = stan::math::value_of(beta_z[s].coeff(0));
            const double z1 = stan::math::value_of(beta_z[s].coeff(1));
            sq_sum += z0 * z0 + z1 * z1;
            d_bz[2 * s] = -z0;
            d_bz[2 * s + 1] = -z1;
        }
        total += -0.5 * sq_sum;
        if constexpr (!propto) total -= static_cast<double>(bz_cnt) * s0d::HALF_LOG_TWO_PI;
    }

    // 1e) log_kappa student_t + (-log_kappa) + total_log_J_of_C
    double lk_val[2], kap_val[2];
    for (int c = 0; c < 2; ++c) {
        const double y = stan::math::value_of(lk_ref.coeff(c));
        const double kp = stan::math::value_of(kap_ref.coeff(c));
        const double nu = stan::math::value_of(kd_ref.coeff(c));
        const double mu = stan::math::value_of(km_ref.coeff(c));
        const double sd = stan::math::value_of(ks_ref.coeff(c));
        if (!(nu > 0.0) || !std::isfinite(nu)) {
            stan::math::check_positive_finite(fn_name, "Degrees of freedom parameter", nu);
        }
        if (!(sd > 0.0) || !std::isfinite(sd)) {
            stan::math::check_positive_finite(fn_name, "Scale parameter", sd);
        }
        if (!std::isfinite(y) || !std::isfinite(mu)) {
            stan::math::check_finite(fn_name, "Location parameter", mu);
        }
        const double z = (y - mu) / sd;
        total += -0.5 * (nu + 1.0) * std::log1p((z * z) / nu) - y;
        if constexpr (!propto) {
            total += std::lgamma(0.5 * (nu + 1.0)) - std::lgamma(0.5 * nu)
                   - 0.5 * std::log(nu * PI_VAL) - std::log(sd);
        }
        d_lk[c] = -(nu + 1.0) * z / (sd * (nu + z * z)) - 1.0;
        lk_val[c] = y;
        kap_val[c] = kp;
        total += stan::math::value_of(tlj_ref.coeff(c));
    }

    // 2) Class Dirichlet / Se_Sp / Set0 priors
    double* cum_alpha[2] = {work_mem, work_mem + (n_cat + 1)};
    double* cur = cum_alpha[1] + (n_cat + 1);
    double alpha_0[2];
    for (int c = 0; c < 2; ++c) {
        cum_alpha[c][0] = 0.0;
        for (int i = 0; i < n_cat; ++i) {
            const double ai = stan::math::value_of(alpha[c].coeff(i));
            if (!(ai > 0.0) || !std::isfinite(ai)) {
                stan::math::check_positive(fn_name, "alpha", ai);
            }
            cum_alpha[c][i + 1] = cum_alpha[c][i] + ai;
        }
        alpha_0[c] = cum_alpha[c][n_cat];
    }

    // Helper: add Dirichlet(phi | prior_alpha) using stack buffer so caller scratch stays untouched
    auto eval_dirichlet_phi = [&](const int c) {
        double stack_pa[512], stack_pk[512], stack_lp[512];
        std::vector<double> heap_dp;
        double *pa = stack_pa, *phi_k = stack_pk, *lphi = stack_lp;
        if (n_cat > 512) {
            heap_dp.resize(3 * static_cast<std::size_t>(n_cat));
            pa = heap_dp.data();
            phi_k = pa + n_cat;
            lphi = phi_k + n_cat;
        }
        double sum_pa = 0.0;
        bool all_unit = true;
        for (int n = 0; n < n_cat; ++n) {
            const double pn = stan::math::value_of(dirichlet_phi[c].coeff(n));
            const double an = stan::math::value_of(prior_dirichlet_alpha[c].coeff(n));
            phi_k[n] = pn;
            pa[n] = an;
            sum_pa += an;
            if (an != 1.0) all_unit = false;
            if constexpr (!propto) total -= std::lgamma(an);
        }
        if constexpr (!propto) total += std::lgamma(sum_pa);
        if (!all_unit) {
            nico::math::copy_n(phi_k, n_cat, lphi);
            nico::math::k_log(lphi, n_cat);
            for (int n = 0; n < n_cat; ++n) {
                const double an = pa[n];
                if (an != 1.0) {
                    total += (an - 1.0) * lphi[n];
                    d_phi[c * n_cat + n] += (an - 1.0) / phi_k[n];
                }
            }
        }
    };

    if (use_Se_Sp_prior == 0 || (mc[0] == 0 && mc[1] == 0)) {
        for (int c = 0; c < 2; ++c) eval_dirichlet_phi(c);
    } else if (use_Se_Sp_prior == 1) {
        // Pack BOTH active classes into ONE medians_and_log_J_packed_raw call (e.g. 4 + 4 = 8 lanes for BAI!)
        int active_cls[2], n_act = 0;
        const double* cum_ptrs[2];
        double a0_ptrs[2], lk_ptrs[2];
        int m_ptrs[2], cuts_storage[2][64];
        const int* cuts_ptrs[2];
        for (int c = 0; c < 2; ++c) {
            if (mc[c] > 0) {
                if (mc[c] > 64) throw std::domain_error(std::string(fn_name) + ": m_c > 64");
                active_cls[n_act] = c;
                cum_ptrs[n_act] = cum_alpha[c];
                a0_ptrs[n_act] = alpha_0[c];
                lk_ptrs[n_act] = lk_val[c];
                m_ptrs[n_act] = mc[c];
                for (int j = 0; j < mc[c]; ++j) {
                    const int cj = chosen_cutoffs[c][j];
                    if (cj < 1 || cj > n_thr || (j > 0 && cj <= chosen_cutoffs[c][j - 1])) {
                        throw std::domain_error(std::string(fn_name) + ": chosen_cutoffs invalid");
                    }
                    cuts_storage[n_act][j] = cj;
                }
                cuts_ptrs[n_act] = cuts_storage[n_act];
                ++n_act;
            } else {
                eval_dirichlet_phi(c);
            }
        }
        double* popmed_all = cur;
        double* logJ_all = popmed_all + M_tot;
        double* dc_da = logJ_all + M_tot;
        double* dc_db = dc_da + M_tot;
        double* dj_da = dc_db + M_tot;
        double* dj_db = dj_da + M_tot;
        double* med_scratch = dj_db + M_tot;
        unsigned char reflect_all[128];

        const double* c_p[2] = {cum_ptrs[0], cum_ptrs[1]};
        const int* k_p[2] = {cuts_ptrs[0], cuts_ptrs[1]};
        cmf::medians_and_log_J_packed_raw<alpha_var>(
            c_p, a0_ptrs, lk_ptrs, k_p, m_ptrs, n_act, M_tot, use_probit_link,
            popmed_all, logJ_all, dc_da, dc_db, dj_da, dj_db, reflect_all, med_scratch);

        double* tmp_after_med = med_scratch;
        int off_m = 0;
        for (int ai = 0; ai < n_act; ++ai) {
            const int c = active_cls[ai];
            const int m = mc[c];
            const double* pop_c = popmed_all + off_m;
            const double* lj_c = logJ_all + off_m;
            double d_pop[64], d_xch[64], x_ch[64], w_seg[65];

            // a) Dirichlet(phi | prior_dirichlet_alpha)
            eval_dirichlet_phi(c);

            // b) Load w_seg & evaluate -Dirichlet(seg_sums(phi) | w_seg)
            bool seg_all_unit = true;
            int start = 0;
            double q_phi[65], lq_phi[65];
            for (int s = 0; s <= m; ++s) {
                const double ws = stan::math::value_of(alpha_of_the_segments[c].coeff(s));
                w_seg[s] = ws;
                if (ws != 1.0) seg_all_unit = false;
                const int end = (s < m) ? cuts_storage[ai][s] : n_cat;
                double q = 0.0;
                for (int n = start; n < end; ++n) q += stan::math::value_of(dirichlet_phi[c].coeff(n));
                q_phi[s] = q;
                lq_phi[s] = q;
                start = end;
            }
            nico::math::k_log(lq_phi, m + 1);
            for (int s = 0; s <= m; ++s) {
                if (w_seg[s] != 1.0) total -= (w_seg[s] - 1.0) * lq_phi[s];
            }
            // Note: include_const for -Dirichlet(q_phi | w_seg) cancels +Dirichlet(q_pop | w_seg) identically!
            if (phi_var) {
                start = 0;
                for (int s = 0; s <= m; ++s) {
                    const int end = (s < m) ? cuts_storage[ai][s] : n_cat;
                    if (w_seg[s] != 1.0) {
                        const double d_seg = (w_seg[s] - 1.0) / q_phi[s];
                        for (int n = start; n < end; ++n) d_phi[c * n_cat + n] -= d_seg;
                    }
                    start = end;
                }
            }

            // c) sum(log_J) + induced Dirichlet prior on popmed_c with A = w_seg (type = 1)
            for (int j = 0; j < m; ++j) {
                total += lj_c[j];
                d_pop[j] = 0.0;
            }
            Se_Sp::buffers_of_the_Se_Sp_prior sbuf;
            sbuf.place(tmp_after_med, m);
            for (int s = 0; s <= m; ++s) {
                sbuf.ones[s] = 1.0;
                sbuf.segments[s] = w_seg[s];
                sbuf.A[s] = w_seg[s];
            }
            if (!Se_Sp::values_of_one_set(pop_c, use_probit_link, sbuf)) {
                throw std::domain_error(std::string(fn_name) + ": C_baseline is not increasing");
            }
            for (int j = 0; j < m; ++j) {
                total += sbuf.log_f[j];
                d_pop[j] += (use_probit_link == 1) ? (-pop_c[j]) : (std::exp(sbuf.lS[j]) - std::exp(sbuf.lF[j]));
            }
            total += Se_Sp::dirichlet_of_the_segments(w_seg, sbuf.log_q, sbuf.log_f, m, false, 1.0, d_pop, nullptr);

            // d) x_chosen = baseline_location_prior_td[c] + pop_c - Xb_base[c]
            const double blp_c = stan::math::value_of(blp_ref.coeff(c));
            const double shift = blp_c - Xb_base[c];
            for (int j = 0; j < m; ++j) {
                x_ch[j] = pop_c[j] + shift;
                d_xch[j] = 0.0;
            }

            // e) -Set0_marginal(x_chosen)
            const double v0_c = stan::math::value_of(v0_ref.coeff(c));
            double* s0_scratch = tmp_after_med + Se_Sp::buffers_of_the_Se_Sp_prior::n_values(m);
            double g_s0[64];
            total -= s0d::eval_set0_marginal_fast_raw<any_var>(
                x_ch, w_seg, m, v0_c, use_probit_link, seg_all_unit, g_s0, s0_scratch);
            for (int j = 0; j < m; ++j) d_xch[j] -= g_s0[j];

            // f) +Se_Sp_prior(x_chosen)
            const int type_of_H = (Se_Sp_prior_type == 1 || Se_Sp_prior_type == 2) ? Se_Sp_prior_type : 4;
            if (type_of_H == 4) {
                for (int j = 0; j < m; ++j) {
                    const double mu = stan::math::value_of(prior_link_mean_SD_td[c].coeff(j, 0));
                    const double sd = stan::math::value_of(prior_link_mean_SD_td[c].coeff(j, 1));
                    const double z = (x_ch[j] - mu) / sd;
                    total += -0.5 * z * z;
                    if constexpr (!propto) total -= std::log(sd) + s0d::HALF_LOG_TWO_PI;
                    d_xch[j] += -z / sd;
                }
            } else {
                const Eigen::Matrix<double, -1, 1> seg_ones = Eigen::Matrix<double, -1, 1>::Ones(m + 1);
                Se_Sp::values_of_the_hyperparameters(
                    fn_name, seg_ones, type_of_H,
                    prior_Se_Sp_dirichlet_A[c].head(m + 1),
                    prior_beta_td[c].topRows(m),
                    prior_link_mean_SD_td[c].topRows(m), sbuf);
                if (type_of_H == 1) {
                    if (!Se_Sp::values_of_one_set(x_ch, use_probit_link, sbuf)) {
                        throw std::domain_error(std::string(fn_name) + ": x_chosen is not increasing");
                    }
                } else {
                    Se_Sp::values_of_one_set_in_any_order(x_ch, use_probit_link, sbuf);
                }
                for (int j = 0; j < m; ++j) {
                    total += sbuf.log_f[j];
                    d_xch[j] += (use_probit_link == 1) ? (-x_ch[j]) : (std::exp(sbuf.lS[j]) - std::exp(sbuf.lF[j]));
                }
                total += Se_Sp::Se_Sp_prior(x_ch, use_probit_link, type_of_H, sbuf, !propto, false, d_xch);
            }

            // g) Backpropagate d_pop and d_xch -> alpha[c], kappa[c], beta_mu[c]
            double sum_dx = 0.0;
            for (int j = 0; j < m; ++j) {
                sum_dx += d_xch[j];
                d_pop[j] += d_xch[j];
            }
            if constexpr (bmu_var) {
                const int nc = (c == 0) ? n_cov_nd : n_cov_d;
                const int off = (c == 0) ? 0 : n_cov_nd;
                const double d_xb = -sum_dx;
                for (int j = 0; j < nc; ++j) {
                    const double bc = (c == 0) ? stan::math::value_of(bcnd_ref.coeff(j))
                                               : stan::math::value_of(bcd_ref.coeff(j));
                    d_bmu[off + j] += d_xb * bc;
                }
            }
            if constexpr (kap_var) {
                d_kap[c] += static_cast<double>(m) / kap_val[c];
            }
            if constexpr (alpha_var) {
                double diff_a[256];
                std::vector<double> diff_heap;
                double* da = (n_cat + 1 <= 256) ? diff_a : (diff_heap.resize(n_cat + 1), diff_heap.data());
                std::fill_n(da, n_cat + 1, 0.0);
                for (int j = 0; j < m; ++j) {
                    const int idx = off_m + j;
                    const double ac = d_pop[j];
                    const double aj = 1.0;
                    double d_araw, d_braw;
                    if (reflect_all[idx] == 0) {
                        d_araw = ac * dc_da[idx] + aj * dj_da[idx];
                        d_braw = ac * dc_db[idx] + aj * dj_db[idx];
                    } else {
                        d_araw = -ac * dc_db[idx] + aj * dj_db[idx];
                        d_braw = -ac * dc_da[idx] + aj * dj_da[idx];
                    }
                    da[0] += d_araw;
                    da[cuts_storage[ai][j]] += d_braw - d_araw;
                }
                double run = 0.0;
                for (int i = 0; i < n_cat; ++i) {
                    run += da[i];
                    d_alpha[c * n_cat + i] += run;
                }
            }
            off_m += m;
        }
    } else {
        // use_Se_Sp_prior == 2
        int active_cls[2], n_act = 0;
        const double* cum_ptrs[2];
        double a0_ptrs[2], lk_ptrs[2];
        int m_ptrs[2];
        std::vector<int> all_cuts(n_thr);
        for (int k = 0; k < n_thr; ++k) all_cuts[k] = k + 1;
        const int* cuts_ptrs[2];
        for (int c = 0; c < 2; ++c) {
            if (mc[c] > 0) {
                active_cls[n_act] = c;
                cum_ptrs[n_act] = cum_alpha[c];
                a0_ptrs[n_act] = alpha_0[c];
                lk_ptrs[n_act] = lk_val[c];
                m_ptrs[n_act] = n_thr;
                cuts_ptrs[n_act] = all_cuts.data();
                ++n_act;
            } else {
                eval_dirichlet_phi(c);
            }
        }
        std::vector<double> p_vec(2 * M_tot + 4 * M_tot);
        std::vector<unsigned char> refl_v(std::max(1, M_tot));
        double* pop_all = p_vec.data();
        double* logJ_all = pop_all + M_tot;
        double* dc_da = logJ_all + M_tot;
        double* dc_db = dc_da + M_tot;
        double* dj_da = dc_db + M_tot;
        double* dj_db = dj_da + M_tot;
        std::vector<double> med_scr(std::max(1, 14 * M_tot + 32));
        cmf::medians_and_log_J_packed_raw<alpha_var>(
            cum_ptrs, a0_ptrs, lk_ptrs, cuts_ptrs, m_ptrs, n_act, M_tot, use_probit_link,
            pop_all, logJ_all, dc_da, dc_db, dj_da, dj_db, refl_v.data(), med_scr.data());

        int off_m = 0;
        for (int ai = 0; ai < n_act; ++ai) {
            const int c = active_cls[ai];
            const int m = mc[c];
            const double* pop_c = pop_all + off_m;
            const double* lj_c = logJ_all + off_m;
            for (int k = 0; k < n_thr; ++k) total += lj_c[k];

            std::vector<double> ind_buf(5 * n_thr + (n_cat) + 2 * m + 64 + Se_Sp::buffers_of_the_Se_Sp_prior::n_values(m));
            double* p_alpha = ind_buf.data();
            double* F_b = p_alpha + n_cat;
            double* p_b = F_b + n_thr;
            double* log_p_b = p_b + n_cat;
            double* J_b = log_p_b + 2 * n_thr + 1;
            double* pop_ch = J_b + n_thr;
            double* d_pop_ch = pop_ch + m;
            double* sp_raw = d_pop_ch + m;
            for (int j = 0; j < n_cat; ++j) p_alpha[j] = stan::math::value_of(prior_dirichlet_alpha[c].coeff(j));
            for (int j = 0; j < m; ++j) {
                const int cj = chosen_cutoffs[c][j];
                if (cj < 1 || cj > n_thr || (j > 0 && cj <= chosen_cutoffs[c][j - 1])) {
                    throw std::domain_error(std::string(fn_name) + ": chosen_cutoffs invalid");
                }
                pop_ch[j] = pop_c[cj - 1];
                d_pop_ch[j] = 0.0;
            }
            double sum_id = 0.0;
            if (!prior::forward_pass<!propto>(pop_c, p_alpha, use_probit_link, n_thr, F_b, p_b, log_p_b, J_b, sum_id)) {
                Eigen::Matrix<double, -1, 1> pv(n_cat);
                for (int j = 0; j < n_cat; ++j) pv[j] = std::exp(log_p_b[j]);
                prior::run_the_checks_of_Stan(pv, Eigen::Map<const Eigen::Matrix<double, -1, 1>>(p_alpha, n_cat));
            }
            total += sum_id;
            if constexpr (!propto) {
                double sum_pa = 0.0, sum_lpa = 0.0;
                for (int j = 0; j < n_cat; ++j) { sum_pa += p_alpha[j]; sum_lpa += std::lgamma(p_alpha[j]); }
                total += std::lgamma(sum_pa) - sum_lpa;
            }
            Se_Sp::buffers_of_the_Se_Sp_prior sbuf;
            sbuf.place(sp_raw, m);
            Se_Sp::values_of_the_hyperparameters(
                fn_name, alpha_of_the_segments[c].head(m + 1), Se_Sp_prior_type,
                prior_Se_Sp_dirichlet_A[c].head(m + 1),
                prior_beta_td[c].topRows(m),
                prior_link_mean_SD_td[c].topRows(m), sbuf);
            if (!Se_Sp::values_of_one_set(pop_ch, use_probit_link, sbuf)) {
                throw std::domain_error(std::string(fn_name) + ": chosen cutpoints not increasing");
            }
            total += Se_Sp::dirichlet_of_the_segments(sbuf.segments, sbuf.log_q, sbuf.log_f, m, !propto, -1.0, d_pop_ch, nullptr);
            total += Se_Sp::Se_Sp_prior(pop_ch, use_probit_link, Se_Sp_prior_type, sbuf, !propto, false, d_pop_ch);

            if constexpr (kap_var) d_kap[c] += static_cast<double>(n_thr) / kap_val[c];
            if constexpr (alpha_var) {
                std::vector<double> d_pop(n_thr), ratios(2 * n_thr), da(n_cat + 1, 0.0);
                prior::reverse_pass(pop_c, p_alpha, use_probit_link, n_thr, F_b, p_b, log_p_b, 1.0,
                                    d_pop.data(), ratios.data(), ratios.data() + n_thr);
                for (int j = 0; j < m; ++j) d_pop[chosen_cutoffs[c][j] - 1] += d_pop_ch[j];
                for (int k = 0; k < n_thr; ++k) {
                    const int idx = off_m + k;
                    const double ac = d_pop[k];
                    double d_araw, d_braw;
                    if (refl_v[idx] == 0) {
                        d_araw = ac * dc_da[idx] + dj_da[idx];
                        d_braw = ac * dc_db[idx] + dj_db[idx];
                    } else {
                        d_araw = -ac * dc_db[idx] + dj_db[idx];
                        d_braw = -ac * dc_da[idx] + dj_da[idx];
                    }
                    da[0] += d_araw;
                    da[k + 1] += d_braw - d_araw;
                }
                double run = 0.0;
                for (int i = 0; i < n_cat; ++i) { run += da[i]; d_alpha[c * n_cat + i] += run; }
            }
            off_m += n_thr;
        }
    }

    if constexpr (!any_var) {
        return total;
    } else {
        stan::arena_t<Eigen::Matrix<stan::scalar_type_t<T_bmu>, -1, -1>> arena_bmu;
        stan::arena_t<Eigen::Matrix<stan::scalar_type_t<T_bsd>, 1, -1>> arena_bsd;
        stan::arena_t<Eigen::Matrix<stan::scalar_type_t<T_bl>, -1, -1>> arena_bl;
        stan::arena_t<Eigen::Matrix<T_bz, -1, -1>> arena_bz(bz_var ? n_studies : 0, bz_var ? 2 : 0);
        stan::arena_t<Eigen::Matrix<stan::scalar_type_t<T_lk>, -1, 1>> arena_lk;
        stan::arena_t<Eigen::Matrix<stan::scalar_type_t<T_kap>, -1, 1>> arena_kap;
        stan::arena_t<Eigen::Matrix<stan::scalar_type_t<T_tlj>, -1, 1>> arena_tlj;
        stan::arena_t<Eigen::Matrix<T_phi, -1, 1>> arena_phi0, arena_phi1;
        stan::arena_t<Eigen::Matrix<T_alpha, -1, 1>> arena_a0, arena_a1;

        if constexpr (bmu_var) arena_bmu = bmu_ref;
        if constexpr (bsd_var) arena_bsd = bsd_ref;
        if constexpr (bl_var) arena_bl = bl_ref;
        if constexpr (bz_var) {
            for (int s = 0; s < n_studies; ++s) {
                arena_bz.coeffRef(s, 0) = beta_z[s].coeff(0);
                arena_bz.coeffRef(s, 1) = beta_z[s].coeff(1);
            }
        }
        if constexpr (lk_var) arena_lk = lk_ref;
        if constexpr (kap_var) arena_kap = kap_ref;
        if constexpr (tlj_var) arena_tlj = tlj_ref;
        if constexpr (phi_var) { arena_phi0 = dirichlet_phi[0]; arena_phi1 = dirichlet_phi[1]; }
        if constexpr (alpha_var) { arena_a0 = alpha[0]; arena_a1 = alpha[1]; }

        return stan::math::make_callback_var(
            total,
            [arena_bmu, arena_bsd, arena_bl, arena_bz, arena_lk, arena_kap, arena_tlj,
             arena_phi0, arena_phi1, arena_a0, arena_a1,
             d_bmu, d_bsd, d_bl11, d_bz, d_lk, d_kap, d_phi, d_alpha,
             n_cov_nd, n_cov_d, n_studies, n_cat](auto& vi) mutable {
                const double adj = vi.adj();
                if constexpr (bmu_var) {
                    for (int j = 0; j < n_cov_nd; ++j) arena_bmu.coeffRef(0, j).adj() += adj * d_bmu[j];
                    for (int j = 0; j < n_cov_d; ++j) arena_bmu.coeffRef(1, j).adj() += adj * d_bmu[n_cov_nd + j];
                }
                if constexpr (bsd_var) {
                    arena_bsd.coeffRef(0).adj() += adj * d_bsd[0];
                    arena_bsd.coeffRef(1).adj() += adj * d_bsd[1];
                }
                if constexpr (bl_var) {
                    arena_bl.coeffRef(1, 1).adj() += adj * (*d_bl11);
                }
                if constexpr (bz_var) {
                    for (int s = 0; s < n_studies; ++s) {
                        arena_bz.coeffRef(s, 0).adj() += adj * d_bz[2 * s];
                        arena_bz.coeffRef(s, 1).adj() += adj * d_bz[2 * s + 1];
                    }
                }
                if constexpr (lk_var) {
                    arena_lk.coeffRef(0).adj() += adj * d_lk[0];
                    arena_lk.coeffRef(1).adj() += adj * d_lk[1];
                }
                if constexpr (kap_var) {
                    arena_kap.coeffRef(0).adj() += adj * d_kap[0];
                    arena_kap.coeffRef(1).adj() += adj * d_kap[1];
                }
                if constexpr (tlj_var) {
                    arena_tlj.coeffRef(0).adj() += adj;
                    arena_tlj.coeffRef(1).adj() += adj;
                }
                if constexpr (phi_var) {
                    for (int i = 0; i < n_cat; ++i) {
                        arena_phi0.coeffRef(i).adj() += adj * d_phi[i];
                        arena_phi1.coeffRef(i).adj() += adj * d_phi[n_cat + i];
                    }
                }
                if constexpr (alpha_var) {
                    for (int i = 0; i < n_cat; ++i) {
                        arena_a0.coeffRef(i).adj() += adj * d_alpha[i];
                        arena_a1.coeffRef(i).adj() += adj * d_alpha[n_cat + i];
                    }
                }
            });
    }
}

}  // namespace nico::math

using nico::math::nico_pinkney_C_and_log_J_both_from_alpha_precomputed;
using nico::math::nico_ord_biv_rc_priors_lpdf;
