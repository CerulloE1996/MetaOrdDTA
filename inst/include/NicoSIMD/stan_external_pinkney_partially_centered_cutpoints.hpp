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
