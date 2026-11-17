#pragma once
// ---- Stan math only: a copy of NicoSIMD/stan_external_Box_Cox_lambda_log_density_of_the_flat_Dirichlet.hpp (md5 9839fb75e2ad434b390ec8431f5cd30a) with its scalar code and
//      its hand-written gradients, the SIMD branches left out, and log, exp, Phi, ... of doubles from
//      std:: and stan::math:: (the buffer helpers k_exp, k_log, ... are in
//      stan_external_functions_on_Stan_math_only.hpp). Where the notes below speak of SIMD or NicoSIMD
//      kernels, the same values come from Stan math here: no -m flags, any CPU.

#include <stan/math/mix.hpp>
#include "stan_external_functions_on_Stan_math_only.hpp"
#include <cmath>
#include <ostream>
#include <vector>

// ---- Log density of the Box-Cox lambda of one class that makes the probabilities of all categories between the
//      chosen cutoffs as close to the flat Dirichlet(1, ..., 1) as the Box-Cox curve allows, the chosen cutoffs held
//      by location and scale:
//        B_k = (x_k^lambda - 1) / lambda,  C_k = B_k or -2 + 4 (B_k - B_min) / (B_max - B_min),
//        z_k = (C_k - location) / scale,  rho_k = (B_k - B_k1) / (B_k2 - B_k1),  z_k = z_k1 + z_gap rho_k;
//      2 chosen cutoffs: log z_gap + 0.5 log sum_i (u_i - u_(i-1))^2, u_k = f(z_k) d rho_k / d lambda (0 at k1, k2);
//      3 or more: log z_gap + log |d rho_k3 / d lambda|.
//      values_and_partials: the value and its partials with respect to (lambda, location, scale), in one pass

namespace fused_stan_math {

static constexpr double BOX_COX_LOG_SQRT_TWO_PI = 0.9189385332046728;

struct Box_Cox_lambda_workspace_struct {
        Eigen::Matrix<double, -1, 1> B;
        Eigen::Matrix<double, -1, 1> dB;
        Eigen::Matrix<double, -1, 1> d2B;
        Eigen::Matrix<double, -1, 1> u;
        Eigen::Matrix<double, -1, 1> du_dlambda;
        Eigen::Matrix<double, -1, 1> du_dlocation;
        Eigen::Matrix<double, -1, 1> du_dscale;
        void resize(const int n_thr) {
                if (B.size() != n_thr) {
                        B.resize(n_thr);
                        dB.resize(n_thr);
                        d2B.resize(n_thr);
                        u.resize(n_thr);
                        du_dlambda.resize(n_thr);
                        du_dlocation.resize(n_thr);
                        du_dscale.resize(n_thr);
                }
        }
};

inline Box_Cox_lambda_workspace_struct& Box_Cox_lambda_workspace() {

        thread_local Box_Cox_lambda_workspace_struct ws;
        return ws;

}

template <bool with_partials, typename T_x>
inline double nico_Box_Cox_lambda_log_density_of_the_flat_Dirichlet_values_and_partials(
                                                const double lambda,
                                                const double location,
                                                const double scale,
                                                const T_x& x,
                                                const std::vector<int>& chosen_cutoffs,
                                                const int use_probit_link,
                                                const int rescale_to_minus_2_to_2,
                                                double* partials) {


        const int n_thr = x.size();
        Box_Cox_lambda_workspace_struct& ws = Box_Cox_lambda_workspace();
        ws.resize(n_thr);
        const int k_1 = chosen_cutoffs[0] - 1;
        const int k_2 = chosen_cutoffs[1] - 1;
        const bool small_lambda = std::abs(lambda) < 1e-4;
        int k_min = 0;
        int k_max = 0;
        //
        // ---- B, dB / d lambda and d2B / d lambda2 at every threshold:
        //
        for (int k = 0; k < n_thr; ++k) {
                const double log_x = stan::math::log(x(k));
                const double E = stan::math::exp(lambda * log_x);
                ws.B(k) = (lambda == 0.0) ? log_x : (E - 1.0) / lambda;
                if (small_lambda) {
                        ws.dB(k) = 0.5 * log_x * log_x + lambda * log_x * log_x * log_x / 3.0;
                        ws.d2B(k) = log_x * log_x * log_x / 3.0;
                } else {
                        ws.dB(k) = (E * (lambda * log_x - 1.0) + 1.0) / (lambda * lambda);
                        ws.d2B(k) = (log_x * log_x * E - 2.0 * ws.dB(k)) / lambda;
                }
                if (x(k) < x(k_min)) k_min = k;
                if (x(k) > x(k_max)) k_max = k;
        }
        //
        // ---- Cutpoints C and dC / d lambda at the chosen cutoffs (rescaled to [-2, 2] or not):
        //
        const double R = ws.B(k_max) - ws.B(k_min);
        const double dR = ws.dB(k_max) - ws.dB(k_min);
        const double multiplier = rescale_to_minus_2_to_2 == 1 ? 4.0 / R : 1.0;
        const double G = ws.B(k_2) - ws.B(k_1);
        const double dG = ws.dB(k_2) - ws.dB(k_1);
        const double d2G = ws.d2B(k_2) - ws.d2B(k_1);
        const double C_gap = multiplier * G;
        const double d_log_C_gap = rescale_to_minus_2_to_2 == 1 ? dG / G - dR / R : dG / G;
        const double inverse_G = 1.0 / G;
        const double log_z_gap = stan::math::log(C_gap) - stan::math::log(scale);
        //
        // ---- 3 or more chosen cutoffs: the Jacobian of the third chosen cutoff:
        //
        if (chosen_cutoffs.size() >= 3) {
                const int k_3 = chosen_cutoffs[2] - 1;
                const double N = ws.B(k_3) - ws.B(k_1);
                const double dN = ws.dB(k_3) - ws.dB(k_1);
                const double d2N = ws.d2B(k_3) - ws.d2B(k_1);
                const double d_rho = (dN * G - N * dG) * inverse_G * inverse_G;
                if constexpr (with_partials) {
                        const double d2_rho = (d2N * G - N * d2G) * inverse_G * inverse_G - 2.0 * dG * d_rho * inverse_G;
                        partials[0] = d_log_C_gap + d2_rho / d_rho;
                        partials[1] = 0.0;
                        partials[2] = -1.0 / scale;
                }
                return log_z_gap + stan::math::log(std::abs(d_rho));
        }
        //
        // ---- 2 chosen cutoffs: u_k = f(z_k) d rho_k / d lambda and its partials:
        //
        const double inverse_scale = 1.0 / scale;
        const double B_lowest = ws.B(k_min);
        const double dB_lowest = ws.dB(k_min);
        for (int k = 0; k < n_thr; ++k) {
                if (k == k_1 || k == k_2) {
                        ws.u(k) = 0.0;
                        if constexpr (with_partials) {
                                ws.du_dlambda(k) = 0.0;
                                ws.du_dlocation(k) = 0.0;
                                ws.du_dscale(k) = 0.0;
                        }
                        continue;
                }
                const double C_k = rescale_to_minus_2_to_2 == 1 ? -2.0 + multiplier * (ws.B(k) - B_lowest) : ws.B(k);
                const double z = (C_k - location) * inverse_scale;
                const double N = ws.B(k) - ws.B(k_1);
                const double dN = ws.dB(k) - ws.dB(k_1);
                const double d_rho = (dN * G - N * dG) * inverse_G * inverse_G;
                double f;
                double df;
                if (use_probit_link == 1) {
                        f = stan::math::exp(-0.5 * z * z - BOX_COX_LOG_SQRT_TWO_PI);
                        df = -z * f;
                } else {
                        const double p = 1.0 / (1.0 + stan::math::exp(-z));
                        f = p * (1.0 - p);
                        df = f * (1.0 - 2.0 * p);
                }
                ws.u(k) = f * d_rho;
                if constexpr (with_partials) {
                        const double d2N = ws.d2B(k) - ws.d2B(k_1);
                        const double d2_rho = (d2N * G - N * d2G) * inverse_G * inverse_G - 2.0 * dG * d_rho * inverse_G;
                        const double dC = rescale_to_minus_2_to_2 == 1 ?
                                          multiplier * ((ws.dB(k) - dB_lowest) - (ws.B(k) - B_lowest) * dR / R) :
                                          ws.dB(k);
                        const double df_d_rho = df * d_rho * inverse_scale;
                        ws.du_dlambda(k) = df_d_rho * dC + f * d2_rho;
                        ws.du_dlocation(k) = -df_d_rho;
                        ws.du_dscale(k) = -df_d_rho * z;
                }
        }
        //
        // ---- S = sum_i (u_i - u_(i-1))^2 and sum_k (2 u_k - u_(k-1) - u_(k+1)) du_k:
        //
        double S = 0.0;
        double A_lambda = 0.0;
        double A_location = 0.0;
        double A_scale = 0.0;
        double u_previous = 0.0;
        for (int k = 0; k < n_thr; ++k) {
                const double u_k = ws.u(k);
                const double D = u_k - u_previous;
                S += D * D;
                if constexpr (with_partials) {
                        const double u_next = (k + 1 < n_thr) ? ws.u(k + 1) : 0.0;
                        const double w = 2.0 * u_k - u_previous - u_next;
                        A_lambda += w * ws.du_dlambda(k);
                        A_location += w * ws.du_dlocation(k);
                        A_scale += w * ws.du_dscale(k);
                }
                u_previous = u_k;
        }
        S += u_previous * u_previous;
        if constexpr (with_partials) {
                const double inverse_S = 1.0 / S;
                partials[0] = d_log_C_gap + A_lambda * inverse_S;
                partials[1] = A_location * inverse_S;
                partials[2] = -inverse_scale + A_scale * inverse_S;
        }

        return log_z_gap + 0.5 * stan::math::log(S);

}

}

template <typename T_lambda,
          typename T_location,
          typename T_scale,
          typename T_x>
inline stan::return_type_t<T_lambda, T_location, T_scale>
nico_Box_Cox_lambda_log_density_of_the_flat_Dirichlet(const T_lambda& lambda,
                                                     const T_location& location,
                                                     const T_scale& scale,
                                                     const T_x& x,
                                                     const std::vector<int>& chosen_cutoffs,
                                                     const int use_probit_link,
                                                     const int rescale_to_minus_2_to_2,
                                                     std::ostream* pstream__) {

        using T = stan::return_type_t<T_lambda, T_location, T_scale>;
        (void) pstream__;
        const auto& x_ref = stan::math::to_ref(stan::math::value_of(x));
        if constexpr (!stan::is_var<T>::value) {
                return fused_stan_math::nico_Box_Cox_lambda_log_density_of_the_flat_Dirichlet_values_and_partials<false>(
                           lambda, location, scale, x_ref, chosen_cutoffs, use_probit_link, rescale_to_minus_2_to_2,
                           nullptr);
        } else {
                double partials[3];
                const double value = fused_stan_math::nico_Box_Cox_lambda_log_density_of_the_flat_Dirichlet_values_and_partials<true>(
                           stan::math::value_of(lambda), stan::math::value_of(location), stan::math::value_of(scale),
                           x_ref, chosen_cutoffs, use_probit_link, rescale_to_minus_2_to_2, partials);
                const double d_lambda = partials[0];
                const double d_location = partials[1];
                const double d_scale = partials[2];
                return stan::math::make_callback_var(value, [lambda, location, scale, d_lambda, d_location, d_scale](
                                                                auto& vi) {

                        if constexpr (stan::is_var<T_lambda>::value) lambda.adj() += vi.adj() * d_lambda;
                        if constexpr (stan::is_var<T_location>::value) location.adj() += vi.adj() * d_location;
                        if constexpr (stan::is_var<T_scale>::value) scale.adj() += vi.adj() * d_scale;

                });
        }

}

namespace fused_stan_math {

template <bool return_raw_and_C, typename T_lambda, typename T_log_x>
inline Eigen::Matrix<stan::return_type_t<T_lambda>, -1, 1>
box_cox_cutpoints_impl(const T_lambda& lambda_arg, const T_log_x& log_x) {

        using T = stan::return_type_t<T_lambda>;
        const int n_thr = static_cast<int>(log_x.size());
        const int out_sz = return_raw_and_C ? (2 * n_thr) : n_thr;
        if (n_thr == 0) {
                return Eigen::Matrix<T, -1, 1>(0);
        }
        const auto& log_x_ref = stan::math::to_ref(stan::math::value_of(log_x));
        const double lam = stan::math::value_of(lambda_arg);

        double stack_B[128];
        double stack_dB[128];
        std::vector<double> heap_B;
        std::vector<double> heap_dB;
        double* B = stack_B;
        double* dB = stack_dB;
        if (n_thr > 128) {
                heap_B.resize(n_thr);
                B = heap_B.data();
                if constexpr (stan::is_var<T>::value) {
                        heap_dB.resize(n_thr);
                        dB = heap_dB.data();
                }
        }

        if (lam == 0.0) {
                for (int k = 0; k < n_thr; ++k) {
                        B[k] = log_x_ref.coeff(k);
                        if constexpr (stan::is_var<T>::value) dB[k] = 0.0;
                }
        } else {
                for (int k = 0; k < n_thr; ++k) {
                        B[k] = lam * log_x_ref.coeff(k);
                }
                fused_stan_math::k_expm1(B, n_thr);
                const double inv_lam = 1.0 / lam;
                const bool small_lam = std::abs(lam) < 1e-5;
                for (int k = 0; k < n_thr; ++k) {
                        const double em1_k = B[k];
                        const double E_k = em1_k + 1.0;
                        B[k] = em1_k * inv_lam;
                        if constexpr (stan::is_var<T>::value) {
                                const double lx_k = log_x_ref.coeff(k);
                                const double lx2 = lx_k * lx_k;
                                const double ll = lam * lx_k;
                                if (small_lam || std::abs(ll) < 1e-4) {
                                        dB[k] = lx2 * (0.5 + ll * ((1.0 / 3.0) + ll * (0.125 + ll * (1.0 / 30.0))));
                                } else {
                                        dB[k] = (E_k * ll - em1_k) * inv_lam * inv_lam;
                                }
                        }
                }
        }

        int k_min = 0;
        int k_max = 0;
        for (int k = 1; k < n_thr; ++k) {
                if (B[k] < B[k_min]) k_min = k;
                if (B[k] > B[k_max]) k_max = k;
        }
        const double B_min = B[k_min];
        const double R = B[k_max] - B_min;

        if constexpr (!stan::is_var<T>::value) {
                Eigen::Matrix<double, -1, 1> out(out_sz);
                if constexpr (return_raw_and_C) {
                        for (int k = 0; k < n_thr; ++k) out.coeffRef(k) = B[k];
                }
                const int c_off = return_raw_and_C ? n_thr : 0;
                if (R > 0.0) {
                        const double mult = 4.0 / R;
                        for (int k = 0; k < n_thr; ++k) {
                                out.coeffRef(c_off + k) = -2.0 + mult * (B[k] - B_min);
                        }
                } else {
                        for (int k = 0; k < n_thr; ++k) out.coeffRef(c_off + k) = 0.0;
                }
                return out;
        } else {
                stan::arena_t<Eigen::Matrix<stan::math::var, -1, 1>> arena_result(out_sz);
                if (R <= 0.0) {
                        if constexpr (return_raw_and_C) {
                                double* arena_dB = stan::math::ChainableStack::instance_->memalloc_.alloc_array<double>(n_thr);
                                for (int k = 0; k < n_thr; ++k) {
                                        arena_result.coeffRef(k) = stan::math::var(new stan::math::vari(B[k], false));
                                        arena_result.coeffRef(n_thr + k) = stan::math::var(new stan::math::vari(0.0, false));
                                        arena_dB[k] = dB[k];
                                }
                                stan::math::var lam_var = lambda_arg;
                                stan::math::reverse_pass_callback(
                                        [lam_var, arena_result, arena_dB, n_thr]() mutable {
                                                double acc = 0.0;
                                                for (int k = 0; k < n_thr; ++k) acc += arena_result.coeff(k).adj() * arena_dB[k];
                                                lam_var.adj() += acc;
                                        });
                        } else {
                                for (int k = 0; k < n_thr; ++k) {
                                        arena_result.coeffRef(k) = stan::math::var(new stan::math::vari(0.0, false));
                                }
                        }
                        return Eigen::Matrix<T, -1, 1>(arena_result);
                }

                const int c_off = return_raw_and_C ? n_thr : 0;
                double* arena_dC = stan::math::ChainableStack::instance_->memalloc_.alloc_array<double>(out_sz);
                double* arena_dB = return_raw_and_C ? (arena_dC + n_thr) : nullptr;
                const double inv_R = 1.0 / R;
                const double mult = 4.0 * inv_R;
                const double dB_min = dB[k_min];
                const double dR = dB[k_max] - dB_min;
                for (int k = 0; k < n_thr; ++k) {
                        const double u_k = (B[k] - B_min) * inv_R;
                        if constexpr (return_raw_and_C) {
                                arena_result.coeffRef(k) = stan::math::var(new stan::math::vari(B[k], false));
                                arena_dB[k] = dB[k];
                        }
                        arena_result.coeffRef(c_off + k) = stan::math::var(new stan::math::vari(-2.0 + 4.0 * u_k, false));
                        arena_dC[k] = mult * ((dB[k] - dB_min) - u_k * dR);
                }
                stan::math::var lam_var = lambda_arg;
                stan::math::reverse_pass_callback(
                        [lam_var, arena_result, arena_dC, arena_dB, n_thr, c_off]() mutable {
                                double acc = 0.0;
                                for (int k = 0; k < n_thr; ++k) {
                                        acc += arena_result.coeff(c_off + k).adj() * arena_dC[k];
                                        if constexpr (return_raw_and_C) {
                                                acc += arena_result.coeff(k).adj() * arena_dB[k];
                                        }
                                }
                                lam_var.adj() += acc;
                        });
                return Eigen::Matrix<T, -1, 1>(arena_result);
        }

}

}  // namespace fused_stan_math

// ---- Fused Box-Cox rescaled cutpoints C in [-2, +2] from precomputed log(cts_thr_shifted):
template <typename T_lambda, typename T_log_x,
          stan::require_stan_scalar_t<T_lambda>* = nullptr,
          stan::require_eigen_col_vector_t<T_log_x>* = nullptr>
inline Eigen::Matrix<stan::return_type_t<T_lambda>, -1, 1>
nico_box_cox_rescaled_cutpoints(const T_lambda& lambda_arg,
                                const T_log_x& log_x,
                                std::ostream* pstream__) {
        (void) pstream__;
        return fused_stan_math::box_cox_cutpoints_impl<false>(lambda_arg, log_x);
}

// ---- Fused Box-Cox returning [C_raw (1..n_thr), C ((n_thr+1)..(2*n_thr))] in one SIMD k_expm1 pass + callback:
template <typename T_lambda, typename T_log_x,
          stan::require_stan_scalar_t<T_lambda>* = nullptr,
          stan::require_eigen_col_vector_t<T_log_x>* = nullptr>
inline Eigen::Matrix<stan::return_type_t<T_lambda>, -1, 1>
nico_box_cox_C_raw_and_C(const T_lambda& lambda_arg,
                         const T_log_x& log_x,
                         std::ostream* pstream__) {
        (void) pstream__;
        return fused_stan_math::box_cox_cutpoints_impl<true>(lambda_arg, log_x);
}
