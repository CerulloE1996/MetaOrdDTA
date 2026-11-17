#pragma once
// ---- Stan math only: a copy of NicoSIMD/stan_external_lkj_corr_cholesky_lpdf.hpp (md5 c907b5a046d63109d9aa50416a6335e3) with its scalar code and
//      its hand-written gradients, the SIMD branches left out, and log, exp, Phi, ... of doubles from
//      std:: and stan::math:: (the buffer helpers k_exp, k_log, ... are in
//      stan_external_functions_on_Stan_math_only.hpp). Where the notes below speak of SIMD or NicoSIMD
//      kernels, the same values come from Stan math here: no -m flags, any CPU.

#include "stan_external_functions_on_Stan_math_only.hpp"

namespace fused_stan_math {

inline double
nico_lkj_constant_from_data(const double eta,
                            const unsigned int K,
                            std::ostream* pstream__) {

    const int Km1 = K - 1;
    double constant;
    if (eta == 1.0) {
        Eigen::Matrix<double, -1, 1> denominator(Km1 / 2);
        for (Eigen::Index k = 1; k <= denominator.rows(); ++k) {
            denominator[k - 1] = stan::math::lgamma(2.0 * k);
        }
        constant = -denominator.sum();
        if (K % 2 == 1) {
            constant -= 0.25 * (K * K - 1) * stan::math::LOG_PI
                        - 0.25 * (Km1 * Km1) * stan::math::LOG_TWO
                        - Km1 * stan::math::lgamma(0.5 * (K + 1));
        } else {
            constant -= 0.25 * K * (K - 2) * stan::math::LOG_PI
                        + 0.25 * (3 * K * K - 4 * K) * stan::math::LOG_TWO
                        + K * stan::math::lgamma(0.5 * K)
                        - Km1 * stan::math::lgamma(static_cast<double>(K));
        }
    } else {
        constant = Km1 * stan::math::lgamma(eta + 0.5 * Km1);
        for (int k = 1; k <= Km1; ++k) {
            constant -= 0.5 * k * stan::math::LOG_PI
                        + stan::math::lgamma(eta + 0.5 * (Km1 - k));
        }
    }
    return constant;

}

template <bool propto__, typename T_L>
inline stan::return_type_t<T_L>
nico_lkj_corr_cholesky_lpdf(const Eigen::Matrix<T_L, Eigen::Dynamic, Eigen::Dynamic>& L,
                           const double eta,
                           std::ostream* pstream__) {

    static_assert(std::is_arithmetic<T_L>::value || stan::is_var<T_L>::value,
                  "The factor must contain arithmetic or reverse-mode scalars.");
    static constexpr const char* function = "lkj_corr_cholesky_lpdf";
    stan::math::check_positive(function, "Shape parameter", eta);
    //// Stan's checks of a Cholesky factor, read on the values where they are; Stan's own check (which copies the
    //// factor) runs only to throw its own error
    bool a_cholesky_factor = L.cols() <= L.rows() && L.cols() > 0;
    for (Eigen::Index n = 1; n < L.cols() && a_cholesky_factor; ++n) {
        for (Eigen::Index m = 0; m < n && m < L.rows(); ++m) {
            if (stan::math::value_of_rec(L.coeff(m, n)) != 0) a_cholesky_factor = false;
        }
    }
    for (Eigen::Index k = 0; k < L.cols() && a_cholesky_factor; ++k) {
        if (!(stan::math::value_of_rec(L.coeff(k, k)) > 0)) a_cholesky_factor = false;
    }
    if (!a_cholesky_factor) stan::math::check_cholesky_factor(function, "Random variable", L);
    const unsigned int K = L.rows();
    if (K == 0) return 0.0;

    stan::return_type_t<T_L> log_prob = 0.0;
    if constexpr (!propto__) {
        log_prob += nico_lkj_constant_from_data(eta, K, pstream__);
    }
    if constexpr (stan::math::include_summand<propto__, T_L>::value) {
        const int Km1 = K - 1;
        const Eigen::Matrix<T_L, Eigen::Dynamic, 1> diagonals = L.diagonal().tail(Km1);
        const auto log_diagonals = stan::math::log(diagonals);
        Eigen::Matrix<stan::return_type_t<T_L>, Eigen::Dynamic, 1> values(Km1);
        for (int k = 0; k < Km1; ++k) {
            values[k] = (Km1 - k - 1) * log_diagonals[k];
        }
        if (eta != 1.0) {
            for (int k = 0; k < Km1; ++k) {
                values[k] += (2.0 * eta - 2.0) * log_diagonals[k];
            }
        }
        log_prob += stan::math::sum(values);
    }
    return log_prob;

}

}  // namespace fused_stan_math

using fused_stan_math::nico_lkj_corr_cholesky_lpdf;
























// ---- Fused stratified bivariate priors across all n_thr thresholds in ONE callback node ----
namespace fused_stan_math {

template <bool propto,
          typename T_mu, typename T_sigma, typename T_L,
          typename T_m0, typename T_m1, typename T_sd0, typename T_sd1, typename T_ssd, typename T_eta,
          stan::require_eigen_col_vector_t<T_m0>* = nullptr,
          stan::require_eigen_col_vector_t<T_m1>* = nullptr,
          stan::require_eigen_col_vector_t<T_sd0>* = nullptr,
          stan::require_eigen_col_vector_t<T_sd1>* = nullptr,
          stan::require_eigen_col_vector_t<T_ssd>* = nullptr,
          stan::require_eigen_col_vector_t<T_eta>* = nullptr>
inline stan::return_type_t<T_mu, T_sigma, T_L>
nico_strat_bivariate_priors_lpdf(
        const std::vector<Eigen::Matrix<T_mu, -1, 1>>& mu,
        const std::vector<Eigen::Matrix<T_sigma, -1, 1>>& sigma,
        const std::vector<Eigen::Matrix<T_L, -1, -1>>& L_Omega,
        const T_m0& prior_mu_mean_0,
        const T_m1& prior_mu_mean_1,
        const T_sd0& prior_mu_SD_0,
        const T_sd1& prior_mu_SD_1,
        const T_ssd& prior_sigma_SD,
        const T_eta& prior_lkj_eta,
        std::ostream* pstream__) {

    static constexpr const char* fn_name = "nico_strat_bivariate_priors_lpdf";
    static constexpr double HALF_LOG_2PI = 0.918938533204672741780329736406;
    using stan::math::var;

    const int n_thr = static_cast<int>(mu.size());
    if (static_cast<int>(sigma.size()) != n_thr || static_cast<int>(L_Omega.size()) != n_thr ||
        prior_mu_mean_0.size() != n_thr || prior_mu_mean_1.size() != n_thr ||
        prior_mu_SD_0.size() != n_thr || prior_mu_SD_1.size() != n_thr ||
        prior_sigma_SD.size() != n_thr || prior_lkj_eta.size() != n_thr) {
        throw std::invalid_argument(std::string(fn_name) + ": inconsistent threshold array/vector sizes");
    }
    if (n_thr == 0) return 0.0;

    const auto& m0_ref  = stan::math::to_ref(prior_mu_mean_0);
    const auto& m1_ref  = stan::math::to_ref(prior_mu_mean_1);
    const auto& sd0_ref = stan::math::to_ref(prior_mu_SD_0);
    const auto& sd1_ref = stan::math::to_ref(prior_mu_SD_1);
    const auto& ssd_ref = stan::math::to_ref(prior_sigma_SD);
    const auto& eta_ref = stan::math::to_ref(prior_lkj_eta);

    constexpr bool is_mu_var    = !stan::is_constant_all<T_mu>::value;
    constexpr bool is_sigma_var = !stan::is_constant_all<T_sigma>::value;
    constexpr bool is_L_var     = !stan::is_constant_all<T_L>::value;
    constexpr bool is_any_var   = is_mu_var || is_sigma_var || is_L_var;

    std::vector<double> d_prim;
    double* d_buf = nullptr;
    if constexpr (!is_any_var) {
        d_prim.resize(static_cast<std::size_t>(1) * n_thr);
        d_buf = d_prim.data();
    } else {
        d_buf = stan::math::ChainableStack::instance_->memalloc_.alloc_array<double>(
            static_cast<std::size_t>(6) * n_thr);
    }

    double* log_l22  = d_buf;
    double* d_mu0    = is_any_var ? (log_l22 + n_thr) : nullptr;
    double* d_mu1    = is_any_var ? (d_mu0 + n_thr) : nullptr;
    double* d_sig0   = is_any_var ? (d_mu1 + n_thr) : nullptr;
    double* d_sig1   = is_any_var ? (d_sig0 + n_thr) : nullptr;
    double* d_l22    = is_any_var ? (d_sig1 + n_thr) : nullptr;

    double sum = 0.0;
    for (int k = 0; k < n_thr; ++k) {
        if (mu[k].size() != 2 || sigma[k].size() != 2 || L_Omega[k].rows() != 2 || L_Omega[k].cols() != 2) {
            throw std::invalid_argument(std::string(fn_name) + ": invalid per-threshold dimensions");
        }
        const double mu0  = stan::math::value_of(mu[k].coeff(0));
        const double mu1  = stan::math::value_of(mu[k].coeff(1));
        const double sig0 = stan::math::value_of(sigma[k].coeff(0));
        const double sig1 = stan::math::value_of(sigma[k].coeff(1));
        const double l00  = stan::math::value_of(L_Omega[k].coeff(0, 0));
        const double l01  = stan::math::value_of(L_Omega[k].coeff(0, 1));
        const double l22  = stan::math::value_of(L_Omega[k].coeff(1, 1));

        const double m0  = stan::math::value_of(m0_ref.coeff(k));
        const double m1  = stan::math::value_of(m1_ref.coeff(k));
        const double sd0 = stan::math::value_of(sd0_ref.coeff(k));
        const double sd1 = stan::math::value_of(sd1_ref.coeff(k));
        const double ssd = stan::math::value_of(ssd_ref.coeff(k));
        const double eta = stan::math::value_of(eta_ref.coeff(k));

        if (!std::isfinite(mu0) || !std::isfinite(mu1) || !std::isfinite(sig0) || !std::isfinite(sig1)) {
            stan::math::check_finite(fn_name, "mu", mu[k]);
            stan::math::check_finite(fn_name, "sigma", sigma[k]);
        }
        if (!(sd0 > 0.0) || !(sd1 > 0.0) || !(ssd > 0.0)) {
            throw std::domain_error(std::string(fn_name) + ": prior SD must be positive");
        }
        if (eta > 0.0) {
            if (!(l00 > 0.0) || !(l22 > 0.0) || l01 != 0.0) {
                stan::math::check_cholesky_factor(fn_name, "L_Omega", L_Omega[k]);
            }
            log_l22[k] = l22;
        } else {
            log_l22[k] = 1.0;
        }

        const double inv_sd0 = 1.0 / sd0;
        const double inv_sd1 = 1.0 / sd1;
        const double inv_ssd = 1.0 / ssd;
        const double zm0 = (mu0 - m0) * inv_sd0;
        const double zm1 = (mu1 - m1) * inv_sd1;
        const double zs0 = sig0 * inv_ssd;
        const double zs1 = sig1 * inv_ssd;

        sum += -0.5 * (zm0 * zm0 + zm1 * zm1 + zs0 * zs0 + zs1 * zs1);
        if constexpr (!propto) {
            sum += -4.0 * HALF_LOG_2PI - std::log(sd0) - std::log(sd1) - 2.0 * std::log(ssd);
            if (eta > 0.0) {
                sum += nico_lkj_constant_from_data(eta, 2, pstream__);
            }
        }
        if constexpr (is_any_var) {
            d_mu0[k]  = -zm0 * inv_sd0;
            d_mu1[k]  = -zm1 * inv_sd1;
            d_sig0[k] = -zs0 * inv_ssd;
            d_sig1[k] = -zs1 * inv_ssd;
            d_l22[k]  = (eta > 0.0 && eta != 1.0) ? ((2.0 * eta - 2.0) / l22) : 0.0;
        }
    }

    k_log(log_l22, n_thr);
    for (int k = 0; k < n_thr; ++k) {
        const double eta = stan::math::value_of(eta_ref.coeff(k));
        if (eta > 0.0 && eta != 1.0) {
            sum += (2.0 * eta - 2.0) * log_l22[k];
        }
    }

    if constexpr (!is_any_var) {
        return propto ? 0.0 : sum;
    } else {
        stan::arena_t<Eigen::Matrix<var, -1, 1>> mu0_vars(is_mu_var ? n_thr : 0);
        stan::arena_t<Eigen::Matrix<var, -1, 1>> mu1_vars(is_mu_var ? n_thr : 0);
        stan::arena_t<Eigen::Matrix<var, -1, 1>> sig0_vars(is_sigma_var ? n_thr : 0);
        stan::arena_t<Eigen::Matrix<var, -1, 1>> sig1_vars(is_sigma_var ? n_thr : 0);
        stan::arena_t<Eigen::Matrix<var, -1, 1>> l22_vars(is_L_var ? n_thr : 0);
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
                l22_vars.coeffRef(k) = L_Omega[k].coeff(1, 1);
            }
        }

        return stan::math::make_callback_var(
            sum,
            [mu0_vars, mu1_vars, sig0_vars, sig1_vars, l22_vars,
             d_mu0, d_mu1, d_sig0, d_sig1, d_l22, n_thr](auto& vi) mutable {
                const double adj = vi.adj();
                for (int k = 0; k < n_thr; ++k) {
                    if constexpr (is_mu_var) {
                        mu0_vars.coeffRef(k).adj() += adj * d_mu0[k];
                        mu1_vars.coeffRef(k).adj() += adj * d_mu1[k];
                    }
                    if constexpr (is_sigma_var) {
                        sig0_vars.coeffRef(k).adj() += adj * d_sig0[k];
                        sig1_vars.coeffRef(k).adj() += adj * d_sig1[k];
                    }
                    if constexpr (is_L_var) {
                        l22_vars.coeffRef(k).adj() += adj * d_l22[k];
                    }
                }
            });
    }
}

}  // namespace fused_stan_math

using fused_stan_math::nico_strat_bivariate_priors_lpdf;
