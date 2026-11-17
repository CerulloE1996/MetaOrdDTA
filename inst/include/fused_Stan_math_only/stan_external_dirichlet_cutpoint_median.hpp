#pragma once
// ---- Stan math only: a copy of NicoSIMD/stan_external_dirichlet_cutpoint_median.hpp (md5 9c41ba767ff028d0c0049d57209948e1) with its scalar code and
//      its hand-written gradients, the SIMD branches left out, and log, exp, Phi, ... of doubles from
//      std:: and stan::math:: (the buffer helpers k_exp, k_log, ... are in
//      stan_external_functions_on_Stan_math_only.hpp). Where the notes below speak of SIMD or NicoSIMD
//      kernels, the same values come from Stan math here: no -m flags, any CPU.

#include <stan/math/mix.hpp>
#include "stan_external_functions_on_Stan_math_only.hpp"
#include <cmath>
#include <limits>
#include <ostream>

namespace fused_stan_math {

static constexpr double CUTPOINT_HALF_LOG_TWO_PI = 0.918938533204672741780329736406;
static constexpr double CUTPOINT_LOG_HALF = -0.693147180559945309417232121458;

// ---- lgamma, digamma and trigamma of x > 0: recurrence up to x >= 10, then the asymptotic series

inline void cutpoint_lgamma_digamma_trigamma(const double x,
                                             double* lgamma_of_x,
                                             double* digamma_of_x,
                                             double* trigamma_of_x) {


        double shifted = x;
        double product = 1.0;
        double digamma_shift = 0.0;
        double trigamma_shift = 0.0;
        while (shifted < 10.0) {
                const double inverse = 1.0 / shifted;
                product *= shifted;
                digamma_shift += inverse;
                trigamma_shift += inverse * inverse;
                shifted += 1.0;
        }
        const double inverse = 1.0 / shifted;
        const double inverse_2 = inverse * inverse;
        const double log_shifted = stan::math::log(shifted);

        *lgamma_of_x = (shifted - 0.5) * log_shifted - shifted + CUTPOINT_HALF_LOG_TWO_PI
                       + inverse * (1.0 / 12.0 - inverse_2 * (1.0 / 360.0 - inverse_2 * (1.0 / 1260.0
                       - inverse_2 * (1.0 / 1680.0 - inverse_2 * (1.0 / 1188.0 - inverse_2 * (691.0 / 360360.0
                       - inverse_2 / 156.0))))))
                       - stan::math::log(product);
        *digamma_of_x = log_shifted - 0.5 * inverse
                        - inverse_2 * (1.0 / 12.0 - inverse_2 * (1.0 / 120.0 - inverse_2 * (1.0 / 252.0
                        - inverse_2 * (1.0 / 240.0 - inverse_2 * (1.0 / 132.0 - inverse_2 * (691.0 / 32760.0
                        - inverse_2 / 12.0))))))
                        - digamma_shift;
        if (trigamma_of_x != nullptr) {
                *trigamma_of_x = inverse + 0.5 * inverse_2
                                 + inverse * inverse_2 * (1.0 / 6.0 - inverse_2 * (1.0 / 30.0 - inverse_2 * (1.0 / 42.0
                                 - inverse_2 * (1.0 / 30.0 - inverse_2 * (5.0 / 66.0 - inverse_2 * (691.0 / 2730.0
                                 - inverse_2 * 7.0 / 6.0))))))
                                 + trigamma_shift;
        }

}

// ---- F = 2F1(1, s; r; x) = sum over n of (s)_n / (r)_n x^n and F_r, and (with_derivatives) the other derivatives the
//      median and its Jacobian need, in one pass over the terms:
//      series = (F, F_s, F_r, F_rr, F_sr, x F_x, x F_rx)

template <bool with_derivatives>
inline void cutpoint_hypergeometric_series(const double s,
                                           const double r,
                                           const double x,
                                           double* series) {


        double term = 1.0;
        double sum_of_inverse_s = 0.0;
        double sum_of_inverse_r = 0.0;
        double sum_of_inverse_r_squared = 0.0;
        double F = 0.0;
        double F_s = 0.0;
        double F_r = 0.0;
        double F_rr = 0.0;
        double F_sr = 0.0;
        double x_F_x = 0.0;
        double x_F_rx = 0.0;
        for (int n = 0; n < 1000000; ++n) {
                F += term;
                F_r -= term * sum_of_inverse_r;
                if constexpr (with_derivatives) {
                        F_s += term * sum_of_inverse_s;
                        F_rr += term * (sum_of_inverse_r * sum_of_inverse_r + sum_of_inverse_r_squared);
                        F_sr -= term * sum_of_inverse_s * sum_of_inverse_r;
                        x_F_x += term * n;
                        x_F_rx -= term * sum_of_inverse_r * n;
                }
                const double s_n = s + n;
                const double r_n = r + n;
                const double weight = with_derivatives
                                      ? 1.0 + n * (1.0 + sum_of_inverse_r) + sum_of_inverse_s * sum_of_inverse_r
                                        + sum_of_inverse_r * sum_of_inverse_r + sum_of_inverse_r_squared
                                      : 1.0 + sum_of_inverse_r;
                if (s_n * x < r_n && term * weight < 1e-17 * F) {
                        break;
                }
                const double inverse_r_n = 1.0 / r_n;
                sum_of_inverse_s += 1.0 / s_n;
                sum_of_inverse_r += inverse_r_n;
                sum_of_inverse_r_squared += inverse_r_n * inverse_r_n;
                term *= s_n * inverse_r_n * x;
        }

        series[0] = F;
        series[2] = F_r;
        if constexpr (with_derivatives) {
                series[1] = F_s;
                series[3] = F_rr;
                series[4] = F_sr;
                series[5] = x_F_x;
                series[6] = x_F_rx;
        }

}

namespace cutpoint_median_fast_v2 {

inline double cutpoint_beta_median_log_halley_scalar(const double a,
                                                     const double b,
                                                     const double log_a,
                                                     const double lbeta);

template <bool with_derivatives>
inline void beta_2F1_hybrid_partials(const double s0,
                                     const double r0,
                                     const double x0,
                                     double* series);

// (the median and log J of one cutpoint by the continued fraction + Halley: used by
//  nico_dirichlet_cutpoint_medians_and_log_J_at_cuts below, defined further down)
template <bool with_partials>
inline void medians_and_log_J_scalar_one(const double a,
                                         const double b,
                                         const double log_a,
                                         const double lbeta,
                                         const double dig_ap1,
                                         const double trig_ap1,
                                         const double dig_b,
                                         const double trig_b,
                                         const double log_kappa,
                                         const double dig_total,
                                         const bool reflect,
                                         const int use_probit_link,
                                         double& popmed_out,
                                         double& logJ_out,
                                         double& d_cut_da_out,
                                         double& d_cut_db_out,
                                         double& d_logJ_da_out,
                                         double& d_logJ_db_out);

}  // namespace cutpoint_median_fast_v2

// ---- log median u of Beta(a, b), a <= b (so u <= log(1 / 2)): delegate to Modified-Lentz CF + Halley

inline double cutpoint_beta_median_log(const double a,
                                       const double b,
                                       const double log_a,
                                       const double lbeta) {

        return cutpoint_median_fast_v2::cutpoint_beta_median_log_halley_scalar(a, b, log_a, lbeta);

}

}

// ---- Median of the Dirichlet cutpoint (the median of Beta(kappa P_lower, kappa P_upper) on the link scale) and
//      the log of its Jacobian, with every partial derivative, in one pass:
//      values = (cutpoint, log J); partials(output, input), inputs (P_lower, P_upper, kappa)

template <bool with_partials>
inline void nico_dirichlet_cutpoint_median_and_log_J_values_and_partials(const double P_lower,
                                                                         const double P_upper,
                                                                         const double kappa,
                                                                         const int use_probit_link,
                                                                         Eigen::Matrix<double, 2, 1>& values,
                                                                         Eigen::Matrix<double, 2, 3, Eigen::DontAlign>& partials) {


        const bool reflect = P_lower > P_upper;
        const double P_a = reflect ? P_upper : P_lower;
        const double P_b = reflect ? P_lower : P_upper;
        const double a = kappa * P_a;
        const double b = kappa * P_b;
        const double log_kappa = stan::math::log(kappa);
        const double log_a = stan::math::log(a);

        double lgamma_a_plus_1, digamma_a_plus_1, trigamma_a_plus_1;
        double lgamma_b, digamma_b, trigamma_b;
        double lgamma_total, digamma_total;
        fused_stan_math::cutpoint_lgamma_digamma_trigamma(a + 1.0, &lgamma_a_plus_1, &digamma_a_plus_1, &trigamma_a_plus_1);
        fused_stan_math::cutpoint_lgamma_digamma_trigamma(b, &lgamma_b, &digamma_b, &trigamma_b);
        fused_stan_math::cutpoint_lgamma_digamma_trigamma(a + b, &lgamma_total, &digamma_total, nullptr);
        const double lbeta = lgamma_a_plus_1 - log_a + lgamma_b - lgamma_total;

        const double u = fused_stan_math::cutpoint_beta_median_log(a, b, log_a, lbeta);
        const double x = stan::math::exp(u);
        const double log_1m_x = stan::math::log1m_exp(u);

        double series[7];
        if ((a + b) * x < a + 1.0) {
                fused_stan_math::cutpoint_median_fast_v2::beta_2F1_hybrid_partials<with_partials>(a + b, a + 1.0, x, series);
        } else {
                fused_stan_math::cutpoint_hypergeometric_series<with_partials>(a + b, a + 1.0, x, series);
        }
        const double F = series[0];
        const double score_difference = log_1m_x - u + digamma_a_plus_1 - digamma_b - series[2] / F;

        double cutpoint;
        double log_link_derivative;
        if (use_probit_link == 1) {
                cutpoint = stan::math::std_normal_log_qf(u);
                if (cutpoint < -fused_stan_math::reference::NORMAL_TAIL_THRESH) {
                        double fraction = -cutpoint;
                        for (int k = fused_stan_math::reference::MILLS_CF_DEPTH; k >= 1; --k) {
                                fraction = -cutpoint + static_cast<double>(k) / fraction;
                        }
                        log_link_derivative = -stan::math::log(fraction);
                } else {
                        log_link_derivative = u + 0.5 * cutpoint * cutpoint + fused_stan_math::CUTPOINT_HALF_LOG_TWO_PI;
                }
        } else {
                cutpoint = u - log_1m_x;
                log_link_derivative = -log_1m_x;
        }

        values[0] = reflect ? -cutpoint : cutpoint;
        values[1] = log_kappa + stan::math::log(score_difference) + fused_stan_math::CUTPOINT_LOG_HALF + lbeta
                    - a * u - (b - 1.0) * log_1m_x + log_link_derivative;

        if constexpr (with_partials) {
                const double F_s = series[1];
                const double F_r = series[2];
                const double F_rr = series[3];
                const double F_sr = series[4];
                const double x_F_x = series[5];
                const double x_F_rx = series[6];
                const double inverse_F_2 = 1.0 / (F * F);
                const double x_over_1m_x = stan::math::exp(u - log_1m_x);

                const double score_a = u - digamma_a_plus_1 + digamma_total + (F_s + F_r) / F;
                const double multiplier = stan::math::exp(fused_stan_math::CUTPOINT_LOG_HALF + lbeta - a * u
                                                               - (b - 1.0) * log_1m_x);
                const double du_da = -score_a * multiplier;
                const double du_db = -(score_a + score_difference) * multiplier;

                const double link_derivative = stan::math::exp(log_link_derivative);
                const double d_log_link_derivative_du = use_probit_link == 1 ? 1.0 + cutpoint * link_derivative
                                                                             : x_over_1m_x;
                const double sign = reflect ? -1.0 : 1.0;

                const double d_difference_da = trigamma_a_plus_1 - ((F_sr + F_rr) * F - F_r * (F_s + F_r)) * inverse_F_2;
                const double d_difference_db = -trigamma_b - (F_sr * F - F_r * F_s) * inverse_F_2;
                const double d_difference_du = -x_over_1m_x - 1.0 - (x_F_rx * F - F_r * x_F_x) * inverse_F_2;

                const double dH_da = d_difference_da / score_difference + digamma_a_plus_1 - 1.0 / a - digamma_total - u;
                const double dH_db = d_difference_db / score_difference + digamma_b - digamma_total - log_1m_x;
                const double dH_du = d_difference_du / score_difference - a + (b - 1.0) * x_over_1m_x
                                     + d_log_link_derivative_du;

                const double d_cutpoint_da = sign * link_derivative * du_da;
                const double d_cutpoint_db = sign * link_derivative * du_db;
                const double d_log_J_da = dH_da + dH_du * du_da;
                const double d_log_J_db = dH_db + dH_du * du_db;

                const int column_a = reflect ? 1 : 0;
                const int column_b = reflect ? 0 : 1;
                partials(0, column_a) = d_cutpoint_da * kappa;
                partials(0, column_b) = d_cutpoint_db * kappa;
                partials(0, 2) = d_cutpoint_da * P_a + d_cutpoint_db * P_b;
                partials(1, column_a) = d_log_J_da * kappa;
                partials(1, column_b) = d_log_J_db * kappa;
                partials(1, 2) = d_log_J_da * P_a + d_log_J_db * P_b + 1.0 / kappa;
        }

}

template <typename T_lower,
          typename T_upper,
          typename T_kappa>
inline Eigen::Matrix<stan::return_type_t<T_lower, T_upper, T_kappa>, -1, 1>
nico_dirichlet_cutpoint_median_and_log_J(const T_lower& P_lower,
                                      const T_upper& P_upper,
                                      const T_kappa& kappa,
                                      const int use_probit_link,
                                      std::ostream* pstream__) {

        using T = stan::return_type_t<T_lower, T_upper, T_kappa>;
        Eigen::Matrix<double, 2, 1> values;
        Eigen::Matrix<double, 2, 3, Eigen::DontAlign> partials;
        (void) pstream__;
        if constexpr (!stan::is_var<T>::value) {
                static_assert(std::is_same<T, double>::value, "double or var arguments only");
                nico_dirichlet_cutpoint_median_and_log_J_values_and_partials<false>(
                    P_lower, P_upper, kappa, use_probit_link, values, partials);
                Eigen::Matrix<double, -1, 1> result(2);
                result[0] = values[0];
                result[1] = values[1];
                return result;
        } else {
                nico_dirichlet_cutpoint_median_and_log_J_values_and_partials<true>(
                    stan::math::value_of(P_lower), stan::math::value_of(P_upper), stan::math::value_of(kappa),
                    use_probit_link, values, partials);
                Eigen::Matrix<T, -1, 1> result(2);
                result[0] = stan::math::var(new stan::math::vari(values[0], false));
                result[1] = stan::math::var(new stan::math::vari(values[1], false));
                const stan::math::var median = result[0];
                const stan::math::var log_J = result[1];
                stan::math::reverse_pass_callback([P_lower, P_upper, kappa, median, log_J, partials]() {

                        if constexpr (stan::is_var<T_lower>::value)
                                P_lower.adj() += median.adj() * partials(0, 0) + log_J.adj() * partials(1, 0);
                        if constexpr (stan::is_var<T_upper>::value)
                                P_upper.adj() += median.adj() * partials(0, 1) + log_J.adj() * partials(1, 1);
                        if constexpr (stan::is_var<T_kappa>::value)
                                kappa.adj() += median.adj() * partials(0, 2) + log_J.adj() * partials(1, 2);

                });
                return result;
        }

}

// ---- The medians of ALL the Dirichlet cutpoints of a class and the sum of the logs of their Jacobians in ONE call:
//      result = (median_1, ..., median_n, sum of log J). The kernel above gives the values and every partial
//      derivative of each cutpoint in one pass; the whole class has one reverse-pass callback.

template <typename T_lower,
          typename T_upper,
          typename T_kappa>
inline Eigen::Matrix<stan::return_type_t<T_lower, T_upper, T_kappa>, -1, 1>
nico_dirichlet_cutpoint_medians_and_sum_of_log_J(const T_lower& P_lower,
                                                 const T_upper& P_upper,
                                                 const T_kappa& kappa,
                                                 const int use_probit_link,
                                                 std::ostream* pstream__) {

        using T = stan::return_type_t<T_lower, T_upper, T_kappa>;
        using T_lower_scalar = stan::value_type_t<T_lower>;
        using T_upper_scalar = stan::value_type_t<T_upper>;
        (void) pstream__;
        const Eigen::Index n = P_lower.size();
        Eigen::Matrix<double, 2, 1> values;
        Eigen::Matrix<double, 2, 3, Eigen::DontAlign> partials;
        if constexpr (!stan::is_var<T>::value) {
                static_assert(std::is_same<T, double>::value, "double or var arguments only");
                const Eigen::Matrix<double, -1, 1> P_lower_values = P_lower;
                const Eigen::Matrix<double, -1, 1> P_upper_values = P_upper;
                Eigen::Matrix<double, -1, 1> result(n + 1);
                double sum_of_log_J = 0.0;
                for (Eigen::Index k = 0; k < n; ++k) {
                        nico_dirichlet_cutpoint_median_and_log_J_values_and_partials<false>(
                            P_lower_values[k], P_upper_values[k], kappa, use_probit_link, values, partials);
                        result[k] = values[0];
                        sum_of_log_J += values[1];
                }
                result[n] = sum_of_log_J;
                return result;
        } else {
                stan::arena_t<Eigen::Matrix<T_lower_scalar, -1, 1>> arena_P_lower = P_lower;
                stan::arena_t<Eigen::Matrix<T_upper_scalar, -1, 1>> arena_P_upper = P_upper;
                stan::arena_t<Eigen::Matrix<double, -1, -1>> arena_partials(6, n);
                stan::arena_t<Eigen::Matrix<stan::math::var, -1, 1>> arena_result(n + 1);
                const double kappa_value = stan::math::value_of(kappa);
                double sum_of_log_J = 0.0;
                for (Eigen::Index k = 0; k < n; ++k) {
                        nico_dirichlet_cutpoint_median_and_log_J_values_and_partials<true>(
                            stan::math::value_of(arena_P_lower.coeff(k)), stan::math::value_of(arena_P_upper.coeff(k)),
                            kappa_value, use_probit_link, values, partials);
                        arena_result.coeffRef(k) = stan::math::var(new stan::math::vari(values[0], false));
                        sum_of_log_J += values[1];
                        for (int j = 0; j < 3; ++j) {
                                arena_partials.coeffRef(j, k) = partials(0, j);
                                arena_partials.coeffRef(3 + j, k) = partials(1, j);
                        }
                }
                arena_result.coeffRef(n) = stan::math::var(new stan::math::vari(sum_of_log_J, false));
                stan::math::reverse_pass_callback([arena_P_lower, arena_P_upper, kappa, arena_partials, arena_result, n]() mutable {

                        const double sum_adj = arena_result.coeff(n).adj();
                        double kappa_adj = 0.0;
                        for (Eigen::Index k = 0; k < n; ++k) {
                                const double median_adj = arena_result.coeff(k).adj();
                                if constexpr (stan::is_var<T_lower_scalar>::value)
                                        arena_P_lower.coeffRef(k).adj() += median_adj * arena_partials.coeff(0, k)
                                                                           + sum_adj * arena_partials.coeff(3, k);
                                if constexpr (stan::is_var<T_upper_scalar>::value)
                                        arena_P_upper.coeffRef(k).adj() += median_adj * arena_partials.coeff(1, k)
                                                                           + sum_adj * arena_partials.coeff(4, k);
                                kappa_adj += median_adj * arena_partials.coeff(2, k) + sum_adj * arena_partials.coeff(5, k);
                        }
                        if constexpr (stan::is_var<T_kappa>::value)
                                kappa.adj() += kappa_adj;

                });
                return Eigen::Matrix<T, -1, 1>(arena_result);
        }

}

// ---- The medians of ALL the Dirichlet cutpoints of a class and the log of the Jacobian of EACH of them in ONE call:
//      result = (median_1, ..., median_n, log J_1, ..., log J_n). The same kernel on every cutpoint; the whole
//      class has one reverse-pass callback.

template <typename T_lower,
          typename T_upper,
          typename T_kappa>
inline Eigen::Matrix<stan::return_type_t<T_lower, T_upper, T_kappa>, -1, 1>
nico_dirichlet_cutpoint_medians_and_log_J_of_each(const T_lower& P_lower,
                                                  const T_upper& P_upper,
                                                  const T_kappa& kappa,
                                                  const int use_probit_link,
                                                  std::ostream* pstream__) {

        using T = stan::return_type_t<T_lower, T_upper, T_kappa>;
        using T_lower_scalar = stan::value_type_t<T_lower>;
        using T_upper_scalar = stan::value_type_t<T_upper>;
        (void) pstream__;
        const Eigen::Index n = P_lower.size();
        Eigen::Matrix<double, 2, 1> values;
        Eigen::Matrix<double, 2, 3, Eigen::DontAlign> partials;
        if constexpr (!stan::is_var<T>::value) {
                static_assert(std::is_same<T, double>::value, "double or var arguments only");
                const Eigen::Matrix<double, -1, 1> P_lower_values = P_lower;
                const Eigen::Matrix<double, -1, 1> P_upper_values = P_upper;
                Eigen::Matrix<double, -1, 1> result(2 * n);
                for (Eigen::Index k = 0; k < n; ++k) {
                        nico_dirichlet_cutpoint_median_and_log_J_values_and_partials<false>(
                            P_lower_values[k], P_upper_values[k], kappa, use_probit_link, values, partials);
                        result[k] = values[0];
                        result[n + k] = values[1];
                }
                return result;
        } else {
                stan::arena_t<Eigen::Matrix<T_lower_scalar, -1, 1>> arena_P_lower = P_lower;
                stan::arena_t<Eigen::Matrix<T_upper_scalar, -1, 1>> arena_P_upper = P_upper;
                stan::arena_t<Eigen::Matrix<double, -1, -1>> arena_partials(6, n);
                stan::arena_t<Eigen::Matrix<stan::math::var, -1, 1>> arena_result(2 * n);
                const double kappa_value = stan::math::value_of(kappa);
                for (Eigen::Index k = 0; k < n; ++k) {
                        nico_dirichlet_cutpoint_median_and_log_J_values_and_partials<true>(
                            stan::math::value_of(arena_P_lower.coeff(k)), stan::math::value_of(arena_P_upper.coeff(k)),
                            kappa_value, use_probit_link, values, partials);
                        arena_result.coeffRef(k) = stan::math::var(new stan::math::vari(values[0], false));
                        arena_result.coeffRef(n + k) = stan::math::var(new stan::math::vari(values[1], false));
                        for (int j = 0; j < 3; ++j) {
                                arena_partials.coeffRef(j, k) = partials(0, j);
                                arena_partials.coeffRef(3 + j, k) = partials(1, j);
                        }
                }
                stan::math::reverse_pass_callback([arena_P_lower, arena_P_upper, kappa, arena_partials, arena_result, n]() mutable {

                        double kappa_adj = 0.0;
                        for (Eigen::Index k = 0; k < n; ++k) {
                                const double median_adj = arena_result.coeff(k).adj();
                                const double log_J_adj = arena_result.coeff(n + k).adj();
                                if constexpr (stan::is_var<T_lower_scalar>::value)
                                        arena_P_lower.coeffRef(k).adj() += median_adj * arena_partials.coeff(0, k)
                                                                           + log_J_adj * arena_partials.coeff(3, k);
                                if constexpr (stan::is_var<T_upper_scalar>::value)
                                        arena_P_upper.coeffRef(k).adj() += median_adj * arena_partials.coeff(1, k)
                                                                           + log_J_adj * arena_partials.coeff(4, k);
                                kappa_adj += median_adj * arena_partials.coeff(2, k) + log_J_adj * arena_partials.coeff(5, k);
                        }
                        if constexpr (stan::is_var<T_kappa>::value)
                                kappa.adj() += kappa_adj;

                });
                return Eigen::Matrix<T, -1, 1>(arena_result);
        }

}


// ---- Fast SIMD (AVX-512 / AVX2 / ARM64 NEON via SIMDe) + Halley Dirichlet cutpoint medians & log-Jacobians directly from (alpha_c, kappa_c, cuts_c) ----

namespace fused_stan_math {
namespace cutpoint_median_fast {

inline double beta_2F1_val_only_scalar(const double s,
                                       const double r,
                                       const double x) {
    double term = 1.0;
    double F = 0.0;
    for (int n = 0; n < 1000000; ++n) {
        F += term;
        const double s_n = s + static_cast<double>(n);
        const double r_n = r + static_cast<double>(n);
        const double s_x = s_n * x;
        if (s_x < r_n && term * r_n < 1e-17 * F * (r_n - s_x)) {
            break;
        }
        term *= s_x / r_n;
    }
    return F;
}

inline double cutpoint_beta_median_log_halley_scalar(const double a,
                                                     const double b,
                                                     const double log_a,
                                                     const double lbeta) {
    double upper_u = CUTPOINT_LOG_HALF;
    double lower_u = -std::numeric_limits<double>::infinity();
    double u = (a >= 1.0) ? std::log((a - 1.0 / 3.0) / (a + b - 2.0 / 3.0))
                          : ((CUTPOINT_LOG_HALF + log_a + lbeta) / a);
    if (!(u < upper_u)) {
        u = upper_u;
    }
    const double s = a + b;
    const double r = a + 1.0;
    for (int iteration = 0; iteration < 200; ++iteration) {
        const double x = std::exp(u);
        const double one_m_x = 1.0 - x;
        const double log_1m_x = std::log1p(-x);
        const double F = beta_2F1_val_only_scalar(s, r, x);
        const double g = a * u + b * log_1m_x - log_a - lbeta + std::log(F) - CUTPOINT_LOG_HALF;
        if (g == 0.0) {
            return u;
        }
        if (g > 0.0) {
            upper_u = u;
        } else {
            lower_u = u;
        }
        const double phi = a / (one_m_x * F);
        const double g_over_gp = g / phi;
        const double dphi_ratio = a - (b - 1.0) * (x / one_m_x) - phi;
        const double halley_corr = 1.0 - 0.5 * g_over_gp * dphi_ratio;
        double next_u = (halley_corr > 0.25) ? (u - g_over_gp / halley_corr) : (u - g_over_gp);
        if (!(next_u < upper_u) || !(next_u > lower_u)) {
            next_u = std::isinf(lower_u) ? (upper_u - 2.0 * (upper_u - u) - 1.0) : (0.5 * (lower_u + upper_u));
        }
        if (std::abs(next_u - u) <= 4e-16 * (1.0 + std::abs(u))) {
            return next_u;
        }
        u = next_u;
    }
    return u;
}


template <bool with_partials>
inline void medians_and_log_J_scalar_one(const double a,
                                         const double b,
                                         const double log_a,
                                         const double lbeta,
                                         const double dig_ap1,
                                         const double trig_ap1,
                                         const double dig_b,
                                         const double trig_b,
                                         const double log_kappa,
                                         const double dig_total,
                                         const bool reflect,
                                         const int use_probit_link,
                                         double& popmed_out,
                                         double& logJ_out,
                                         double& d_cut_da_out,
                                         double& d_cut_db_out,
                                         double& d_logJ_da_out,
                                         double& d_logJ_db_out) {
    const double u = cutpoint_beta_median_log_halley_scalar(a, b, log_a, lbeta);
    const double x = std::exp(u);
    const double log_1m_x = std::log1p(-x);
    double series[7];
    cutpoint_hypergeometric_series<with_partials>(a + b, a + 1.0, x, series);
    const double F = series[0];
    const double score_diff = log_1m_x - u + dig_ap1 - dig_b - series[2] / F;

    double cutpoint;
    double lld;
    if (use_probit_link == 1) {
        cutpoint = stan::math::std_normal_log_qf(u);
        if (cutpoint < -fused_stan_math::reference::NORMAL_TAIL_THRESH) {
            double frac = -cutpoint;
            for (int kk = fused_stan_math::reference::MILLS_CF_DEPTH; kk >= 1; --kk) {
                frac = -cutpoint + static_cast<double>(kk) / frac;
            }
            lld = -std::log(frac);
        } else {
            lld = u + 0.5 * cutpoint * cutpoint + CUTPOINT_HALF_LOG_TWO_PI;
        }
    } else {
        cutpoint = u - log_1m_x;
        lld = -log_1m_x;
    }

    popmed_out = reflect ? -cutpoint : cutpoint;
    logJ_out = log_kappa + std::log(score_diff) + CUTPOINT_LOG_HALF + lbeta - a * u - (b - 1.0) * log_1m_x + lld;
    if constexpr (with_partials) {
        const double F_s = series[1], F_r = series[2], F_rr = series[3], F_sr = series[4];
        const double x_F_x = series[5], x_F_rx = series[6];
        const double inv_F2 = 1.0 / (F * F);
        const double x_over_1mx = std::exp(u - log_1m_x);
        const double score_a = u - dig_ap1 + dig_total + (F_s + F_r) / F;
        const double mult = std::exp(CUTPOINT_LOG_HALF + lbeta - a * u - (b - 1.0) * log_1m_x);
        const double du_da = -score_a * mult;
        const double du_db = -(score_a + score_diff) * mult;
        const double link_d = std::exp(lld);
        const double d_lld_du = (use_probit_link == 1) ? (1.0 + cutpoint * link_d) : x_over_1mx;
        const double d_diff_da = trig_ap1 - ((F_sr + F_rr) * F - F_r * (F_s + F_r)) * inv_F2;
        const double d_diff_db = -trig_b - (F_sr * F - F_r * F_s) * inv_F2;
        const double d_diff_du = -x_over_1mx - 1.0 - (x_F_rx * F - F_r * x_F_x) * inv_F2;
        const double dH_da = d_diff_da / score_diff + dig_ap1 - 1.0 / a - dig_total - u;
        const double dH_db = d_diff_db / score_diff + dig_b - dig_total - log_1m_x;
        const double dH_du = d_diff_du / score_diff - a + (b - 1.0) * x_over_1mx + d_lld_du;
        d_cut_da_out = link_d * du_da;
        d_cut_db_out = link_d * du_db;
        d_logJ_da_out = dH_da + dH_du * du_da;
        d_logJ_db_out = dH_db + dH_du * du_db;
    }
}

}  // namespace cutpoint_median_fast

template <typename T_alpha, typename T_kappa,
          stan::require_eigen_col_vector_t<T_alpha>* = nullptr>
inline Eigen::Matrix<stan::return_type_t<T_alpha, T_kappa>, -1, 1>
nico_dirichlet_cutpoint_medians_and_log_J_at_cuts(
        const T_alpha& alpha_c,
        const T_kappa& kappa_c,
        const std::vector<int>& cuts_c,
        const int use_probit_link,
        std::ostream* /*pstream__*/) {
    // (the median of each cutpoint by cutpoint_median_fast_v2, the continued fraction + Halley: the older
    //  cutpoint_median_fast series of the NicoSIMD file gave NaN log densities and gradients of 1e29 to 1e52
    //  from kappa about 7,100 with the Se / Sp prior on; the same arguments and outputs)
    // namespace cmf = fused_stan_math::cutpoint_median_fast;
    namespace cmf = fused_stan_math::cutpoint_median_fast_v2;
    using T = stan::return_type_t<T_alpha, T_kappa>;
    static constexpr const char* fn_name = "nico_dirichlet_cutpoint_medians_and_log_J_at_cuts";

    const int n_cat = static_cast<int>(alpha_c.size());
    const int n_thr = n_cat - 1;
    const int m = static_cast<int>(cuts_c.size());
    if (n_cat < 2) throw std::domain_error(std::string(fn_name) + ": alpha_c must have length >= 2");
    if (m == 0) return Eigen::Matrix<T, -1, 1>(0);

    const auto& alpha_ref = stan::math::to_ref(alpha_c);
    const double kappa_val = stan::math::value_of(kappa_c);
    if (!(kappa_val > 0.0) || !std::isfinite(kappa_val)) {
        stan::math::check_positive(fn_name, "kappa", kappa_val);
    }

    constexpr bool is_alpha_var = !stan::is_constant_all<T_alpha>::value;
    constexpr bool is_kappa_var = !stan::is_constant_all<T_kappa>::value;
    constexpr bool is_any_var = is_alpha_var || is_kappa_var;

    const std::size_t n_dbl = static_cast<std::size_t>(12) * m + static_cast<std::size_t>(n_cat) + 8;
    std::vector<double> d_prim;
    std::vector<unsigned char> r_prim;
    double* d_buf = nullptr;
    unsigned char* reflect = nullptr;
    int* cuts_buf = nullptr;
    if constexpr (!is_any_var) {
        d_prim.resize(n_dbl);
        r_prim.resize(m);
        d_buf = d_prim.data();
        reflect = r_prim.data();
    } else {
        auto& arena = stan::math::ChainableStack::instance_->memalloc_;
        d_buf = arena.alloc_array<double>(n_dbl);
        reflect = arena.alloc_array<unsigned char>(m);
        cuts_buf = arena.alloc_array<int>(m);
    }

    double* cum_alpha = d_buf;
    double* a_arr = cum_alpha + (n_cat + 1);
    double* b_arr = a_arr + m;
    double* log_a_arr = b_arr + m;
    double* lbeta_arr = log_a_arr + m;
    double* gamma_in = lbeta_arr + m;
    double* gamma_lg = gamma_in + 2 * m + 1;
    double* gamma_dg = gamma_lg + 2 * m + 1;
    double* gamma_tg = gamma_dg + 2 * m + 1;
    cum_alpha[0] = 0.0;
    for (int i = 0; i < n_cat; ++i) {
        const double ai = stan::math::value_of(alpha_ref.coeff(i));
        if (!(ai > 0.0) || !std::isfinite(ai)) {
            stan::math::check_positive(fn_name, "alpha", ai);
        }
        cum_alpha[i + 1] = cum_alpha[i] + ai;
    }
    const double alpha_0 = cum_alpha[n_cat];
    const double log_kappa = std::log(kappa_val);

    for (int j = 0; j < m; ++j) {
        const int cj = cuts_c[j];
        if (cj < 1 || cj > n_thr || (j > 0 && cj <= cuts_c[j - 1])) {
            throw std::domain_error(std::string(fn_name) + ": cuts_c must strictly increase in 1..n_thr");
        }
        if constexpr (is_any_var) cuts_buf[j] = cj;
        const double a_raw = cum_alpha[cj];
        const double b_raw = alpha_0 - a_raw;
        const bool refl = (a_raw > b_raw);
        reflect[j] = refl ? 1u : 0u;
        const double aj = refl ? b_raw : a_raw;
        const double bj = refl ? a_raw : b_raw;
        a_arr[j] = aj;
        b_arr[j] = bj;
        log_a_arr[j] = std::log(aj);
        gamma_in[j] = aj + 1.0;
        gamma_in[m + j] = bj;
    }
    gamma_in[2 * m] = alpha_0;

    for (int j = 0; j < 2 * m + 1; ++j) {
        cutpoint_lgamma_digamma_trigamma(gamma_in[j], gamma_lg + j, gamma_dg + j, gamma_tg + j);
    }
    const double lg_total = gamma_lg[2 * m];
    const double dig_total = gamma_dg[2 * m];
    for (int j = 0; j < m; ++j) {
        lbeta_arr[j] = gamma_lg[j] - log_a_arr[j] + gamma_lg[m + j] - lg_total;
    }

    std::vector<double> out_prim;
    double* out_vals = nullptr;
    double* d_cut_da = nullptr;
    double* d_cut_db = nullptr;
    double* d_logJ_da = nullptr;
    double* d_logJ_db = nullptr;
    if constexpr (!is_any_var) {
        out_prim.resize(2 * m);
        out_vals = out_prim.data();
    } else {
        auto& arena = stan::math::ChainableStack::instance_->memalloc_;
        out_vals = arena.alloc_array<double>(6 * m);
        d_cut_da = out_vals + 2 * m;
        d_cut_db = d_cut_da + m;
        d_logJ_da = d_cut_db + m;
        d_logJ_db = d_logJ_da + m;
    }
    double* popmed_vals = out_vals;
    double* logJ_vals = out_vals + m;

    int j0 = 0;
    for (int j = j0; j < m; ++j) {
        double dummy0 = 0.0, dummy1 = 0.0, dummy2 = 0.0, dummy3 = 0.0;
        cmf::medians_and_log_J_scalar_one<is_any_var>(
            a_arr[j], b_arr[j], log_a_arr[j], lbeta_arr[j],
            gamma_dg[j], gamma_tg[j], gamma_dg[m + j], gamma_tg[m + j],
            log_kappa, dig_total, reflect[j] != 0, use_probit_link,
            popmed_vals[j], logJ_vals[j],
            is_any_var ? d_cut_da[j] : dummy0,
            is_any_var ? d_cut_db[j] : dummy1,
            is_any_var ? d_logJ_da[j] : dummy2,
            is_any_var ? d_logJ_db[j] : dummy3);
    }

    if constexpr (!is_any_var) {
        Eigen::Matrix<double, -1, 1> result(2 * m);
        for (int j = 0; j < 2 * m; ++j) result[j] = out_vals[j];
        return result;
    } else {
        stan::arena_t<Eigen::Matrix<stan::value_type_t<T_alpha>, -1, 1>> arena_alpha;
        if constexpr (is_alpha_var) arena_alpha = alpha_ref;
        stan::arena_t<Eigen::Matrix<stan::math::var, -1, 1>> arena_result(2 * m);
        for (int j = 0; j < 2 * m; ++j) {
            arena_result.coeffRef(j) = stan::math::var(new stan::math::vari(out_vals[j], false));
        }
        stan::math::reverse_pass_callback(
            [arena_alpha, kappa_c, arena_result, cuts_buf, reflect,
             d_cut_da, d_cut_db, d_logJ_da, d_logJ_db,
             kappa_val, m, n_cat]() mutable {
                double adj_k_sum = 0.0;
                double stack_diff[256];
                std::vector<double> heap_diff;
                double* diff = nullptr;
                if constexpr (is_alpha_var) {
                    diff = (n_cat + 1 <= 256) ? stack_diff : (heap_diff.resize(n_cat + 1), heap_diff.data());
                    std::fill_n(diff, n_cat + 1, 0.0);
                }
                for (int j = 0; j < m; ++j) {
                    const double adj_c = arena_result.coeff(j).adj();
                    const double adj_j = arena_result.coeff(m + j).adj();
                    adj_k_sum += adj_j;
                    if constexpr (is_alpha_var) {
                        double da_raw, db_raw;
                        if (reflect[j] == 0) {
                            da_raw = adj_c * d_cut_da[j] + adj_j * d_logJ_da[j];
                            db_raw = adj_c * d_cut_db[j] + adj_j * d_logJ_db[j];
                        } else {
                            da_raw = -adj_c * d_cut_db[j] + adj_j * d_logJ_db[j];
                            db_raw = -adj_c * d_cut_da[j] + adj_j * d_logJ_da[j];
                        }
                        diff[0] += da_raw;
                        diff[cuts_buf[j]] += db_raw - da_raw;
                    }
                }
                if constexpr (is_kappa_var) {
                    kappa_c.adj() += adj_k_sum / kappa_val;
                }
                if constexpr (is_alpha_var) {
                    double running = 0.0;
                    for (int i = 0; i < n_cat; ++i) {
                        running += diff[i];
                        arena_alpha.coeffRef(i).adj() += running;
                    }
                }
            });
        return Eigen::Matrix<T, -1, 1>(arena_result);
    }
}

}  // namespace fused_stan_math

using fused_stan_math::nico_dirichlet_cutpoint_medians_and_log_J_at_cuts;


// ---- Hybrid CF/Series multi-class packed Dirichlet cutpoint medians & log-Jacobians ----
namespace fused_stan_math {
namespace cutpoint_median_fast_v2 {

// Gauss / Modified-Lentz Continued Fraction for F = 2F1(a + b, 1; a + 1; x) = betacf(a, b, x).
// Valid whenever (a + b) * x < a + 1 (always holds at the median for 0 < a <= b).
inline double beta_2F1_cf_scalar(const double a,
                                 const double b,
                                 const double x) {
    const double qab = a + b;
    const double qap = a + 1.0;
    const double qam = a - 1.0;
    double c = 1.0;
    double d = 1.0 - qab * x / qap;
    if (std::abs(d) < 1e-300) d = 1e-300;
    d = 1.0 / d;
    double h = d;
    for (int m = 1; m <= 2000; ++m) {
        const double md = static_cast<double>(m);
        const double m2 = 2.0 * md;
        double aa = md * (b - md) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if (std::abs(d) < 1e-300) d = 1e-300;
        c = 1.0 + aa / c;
        if (std::abs(c) < 1e-300) c = 1e-300;
        d = 1.0 / d;
        h *= d * c;

        aa = -(a + md) * (qab + md) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if (std::abs(d) < 1e-300) d = 1e-300;
        c = 1.0 + aa / c;
        if (std::abs(c) < 1e-300) c = 1e-300;
        d = 1.0 / d;
        const double del = d * c;
        h *= del;
        if (std::abs(del - 1.0) <= 1e-16) {
            break;
        }
    }
    return h;
}

// Evaluate F = 2F1(s, 1; r; x) (series if s * x < 0.7 * r && <= 32 terms; CF if larger)
inline double beta_2F1_val_only_scalar(const double s,
                                       const double r,
                                       const double x) {
    const double a = r - 1.0;
    const double b = s - a;
    if (s * x >= 0.7 * r && s * x < r) {
        return beta_2F1_cf_scalar(a, b, x);
    }
    double term = 1.0;
    double F = 0.0;
    for (int n = 0; n < 32; ++n) {
        F += term;
        const double s_n = s + static_cast<double>(n);
        const double r_n = r + static_cast<double>(n);
        const double s_x = s_n * x;
        if (s_x < r_n && term * r_n < 1e-17 * F * (r_n - s_x)) {
            return F;
        }
        term *= s_x / r_n;
    }
    if (s * x < r) {
        return beta_2F1_cf_scalar(a, b, x);
    }
    for (int n = 32; n < 1000000; ++n) {
        F += term;
        const double s_n = s + static_cast<double>(n);
        const double r_n = r + static_cast<double>(n);
        const double s_x = s_n * x;
        if (s_x < r_n && term * r_n < 1e-17 * F * (r_n - s_x)) {
            break;
        }
        term *= s_x / r_n;
    }
    return F;
}

// 7-component jet for exact first + mixed-r second partials of F = 2F1(s, 1; r; x) via Modified-Lentz CF:
// series = {F, F_s, F_r, F_rr, F_sr, x*F_x, x*F_rx}
struct Jet7 {
    double v, s, r, x, rr, sr, rx;
};
inline Jet7 j_add(const Jet7& a, const Jet7& b) {
    return {a.v + b.v, a.s + b.s, a.r + b.r, a.x + b.x, a.rr + b.rr, a.sr + b.sr, a.rx + b.rx};
}
inline Jet7 j_sub(const Jet7& a, const Jet7& b) {
    return {a.v - b.v, a.s - b.s, a.r - b.r, a.x - b.x, a.rr - b.rr, a.sr - b.sr, a.rx - b.rx};
}
inline Jet7 j_mul(const Jet7& a, const Jet7& b) {
    return {
        a.v * b.v,
        a.s * b.v + a.v * b.s,
        a.r * b.v + a.v * b.r,
        a.x * b.v + a.v * b.x,
        a.rr * b.v + 2.0 * a.r * b.r + a.v * b.rr,
        a.sr * b.v + a.s * b.r + a.r * b.s + a.v * b.sr,
        a.rx * b.v + a.r * b.x + a.x * b.r + a.v * b.rx
    };
}
inline Jet7 j_inv(const Jet7& a) {
    const double iv = 1.0 / a.v;
    const double iv2 = iv * iv;
    const double iv3 = iv2 * iv;
    return {
        iv,
        -a.s * iv2,
        -a.r * iv2,
        -a.x * iv2,
        2.0 * a.r * a.r * iv3 - a.rr * iv2,
        2.0 * a.s * a.r * iv3 - a.sr * iv2,
        2.0 * a.r * a.x * iv3 - a.rx * iv2
    };
}

template <bool with_derivatives>
inline void beta_2F1_hybrid_partials(const double s0,
                                     const double r0,
                                     const double x0,
                                     double* series) {
    if (s0 * x0 < 0.65 * r0) {
        double term = 1.0;
        double sum_s = 0.0, sum_r = 0.0, sum_r2 = 0.0;
        double F = 0.0, Fs = 0.0, Fr = 0.0, Frr = 0.0, Fsr = 0.0, xFx = 0.0, xFrx = 0.0;
        for (int n = 0; n < 48; ++n) {
            F += term;
            Fr -= term * sum_r;
            if constexpr (with_derivatives) {
                Fs += term * sum_s;
                Frr += term * (sum_r * sum_r + sum_r2);
                Fsr -= term * sum_s * sum_r;
                xFx += term * static_cast<double>(n);
                xFrx -= term * sum_r * static_cast<double>(n);
            }
            const double sn = s0 + static_cast<double>(n);
            const double rn = r0 + static_cast<double>(n);
            const double w = with_derivatives
                ? (1.0 + static_cast<double>(n) * (1.0 + sum_r) + sum_s * sum_r + sum_r * sum_r + sum_r2)
                : (1.0 + sum_r);
            if (sn * x0 < rn && term * w < 1e-17 * F) {
                series[0] = F; series[2] = Fr;
                if constexpr (with_derivatives) {
                    series[1] = Fs; series[3] = Frr; series[4] = Fsr;
                    series[5] = xFx; series[6] = xFrx;
                }
                return;
            }
            const double irn = 1.0 / rn;
            if constexpr (with_derivatives) {
                sum_s += 1.0 / sn;
                sum_r2 += irn * irn;
            }
            sum_r += irn;
            term *= (sn * irn) * x0;
        }
    }
    const Jet7 s = {s0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    const Jet7 r = {r0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0};
    const Jet7 x = {x0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0};
    const Jet7 one = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    Jet7 c = one;
    Jet7 d = j_sub(one, j_mul(j_mul(s, x), j_inv(r)));
    if (std::abs(d.v) < 1e-300) d.v = 1e-300;
    d = j_inv(d);
    Jet7 h = d;
    for (int m = 1; m <= 2000; ++m) {
        const double md = static_cast<double>(m);
        const double m2 = 2.0 * md;
        const Jet7 num_e = {
            md * (s0 - r0 + 1.0 - md) * x0,
            md * x0,
            -md * x0,
            md * (s0 - r0 + 1.0 - md),
            0.0, 0.0, -md
        };
        const double rd1 = r0 + m2 - 2.0;
        const double rd2 = r0 + m2 - 1.0;
        const Jet7 den_e = {rd1 * rd2, 0.0, rd1 + rd2, 0.0, 2.0, 0.0, 0.0};
        Jet7 aa = j_mul(num_e, j_inv(den_e));
        d = j_add(one, j_mul(aa, d));
        if (std::abs(d.v) < 1e-300) d.v = 1e-300;
        c = j_add(one, j_mul(aa, j_inv(c)));
        if (std::abs(c.v) < 1e-300) c.v = 1e-300;
        d = j_inv(d);
        h = j_mul(h, j_mul(d, c));

        const double ra = r0 + md - 1.0;
        const double sa = s0 + md;
        const Jet7 num_o = {
            -ra * sa * x0,
            -ra * x0,
            -sa * x0,
            -ra * sa,
            0.0, -x0, -sa
        };
        const double rd3 = r0 + m2 - 1.0;
        const double rd4 = r0 + m2;
        const Jet7 den_o = {rd3 * rd4, 0.0, rd3 + rd4, 0.0, 2.0, 0.0, 0.0};
        aa = j_mul(num_o, j_inv(den_o));
        d = j_add(one, j_mul(aa, d));
        if (std::abs(d.v) < 1e-300) d.v = 1e-300;
        c = j_add(one, j_mul(aa, j_inv(c)));
        if (std::abs(c.v) < 1e-300) c.v = 1e-300;
        d = j_inv(d);
        const Jet7 del = j_mul(d, c);
        h = j_mul(h, del);
        if (std::abs(del.v - 1.0) <= 1e-16 && std::abs(del.r) <= 1e-15 * (1.0 + std::abs(h.r))) {
            break;
        }
    }
    series[0] = h.v;
    series[2] = h.r;
    if constexpr (with_derivatives) {
        series[1] = h.s;
        series[3] = h.rr;
        series[4] = h.sr;
        series[5] = x0 * h.x;
        series[6] = x0 * h.rx;
    }
}

// Bracketed Halley root-finder using CF + symmetry I_x(a, b) = 1 - I_{1-x}(b, a) when x > (a+1)/(a+b)
inline double cutpoint_beta_median_log_halley_scalar(const double a,
                                                     const double b,
                                                     const double log_a,
                                                     const double lbeta) {
    double upper_u = CUTPOINT_LOG_HALF;
    double lower_u = -std::numeric_limits<double>::infinity();
    double u = (a >= 1.0) ? std::log((a - 1.0 / 3.0) / (a + b - 2.0 / 3.0))
                          : ((CUTPOINT_LOG_HALF + log_a + lbeta) / a);
    if (!(u < upper_u)) {
        u = upper_u;
    }
    const double log_b = std::log(b);
    for (int iteration = 0; iteration < 200; ++iteration) {
        const double x = std::exp(u);
        const double one_m_x = 1.0 - x;
        const double log_1m_x = std::log1p(-x);
        double g, phi;
        if ((a + b) * x < a + 1.0) {
            const double F = beta_2F1_val_only_scalar(a + b, a + 1.0, x);
            g = a * u + b * log_1m_x - log_a - lbeta + std::log(F) - CUTPOINT_LOG_HALF;
            phi = a / (one_m_x * F);
        } else {
            const double F_sym = beta_2F1_cf_scalar(b, a, one_m_x);
            const double log_I1mx = b * log_1m_x + a * u - log_b - lbeta + std::log(F_sym);
            const double I1mx = (log_I1mx < 0.0) ? std::exp(log_I1mx) : 1.0;
            const double Ix = std::max(1.0 - I1mx, 1e-300);
            g = std::log1p(-std::min(I1mx, 1.0 - 1e-300)) - CUTPOINT_LOG_HALF;
            phi = std::exp(a * u + (b - 1.0) * log_1m_x - lbeta - std::log(Ix));
        }
        if (g == 0.0) {
            return u;
        }
        if (g > 0.0) {
            upper_u = u;
        } else {
            lower_u = u;
        }
        const double g_over_gp = g / phi;
        const double dphi_ratio = a - (b - 1.0) * (x / one_m_x) - phi;
        const double halley_corr = 1.0 - 0.5 * g_over_gp * dphi_ratio;
        double next_u = (halley_corr > 0.25) ? (u - g_over_gp / halley_corr) : (u - g_over_gp);

        const double tol = 4e-16 * (1.0 + std::abs(u));
        if (std::abs(next_u - u) <= tol) {
            return next_u;
        }
        if (!std::isinf(lower_u) && (upper_u - lower_u) <= 8e-15 * (1.0 + std::abs(u)) && std::abs(g) <= 1e-10) {
            return 0.5 * (lower_u + upper_u);
        }
        if (!(next_u < upper_u) || !(next_u > lower_u)) {
            next_u = std::isinf(lower_u) ? (upper_u - 2.0 * (upper_u - u) - 1.0) : (0.5 * (lower_u + upper_u));
        }
        if (std::abs(next_u - u) <= tol) {
            return next_u;
        }
        u = next_u;
    }
    return u;
}


template <bool with_partials>
inline void medians_and_log_J_scalar_one(const double a,
                                         const double b,
                                         const double log_a,
                                         const double lbeta,
                                         const double dig_ap1,
                                         const double trig_ap1,
                                         const double dig_b,
                                         const double trig_b,
                                         const double log_kappa,
                                         const double dig_total,
                                         const bool reflect,
                                         const int use_probit_link,
                                         double& popmed_out,
                                         double& logJ_out,
                                         double& d_cut_da_out,
                                         double& d_cut_db_out,
                                         double& d_logJ_da_out,
                                         double& d_logJ_db_out) {
    const double u = cutpoint_beta_median_log_halley_scalar(a, b, log_a, lbeta);
    const double x = std::exp(u);
    const double omx = 1.0 - x;
    const double log_1m_x = std::log1p(-x);
    double series[7];
    beta_2F1_hybrid_partials<with_partials>(a + b, a + 1.0, x, series);
    const double F = series[0];
    const double score_diff = log_1m_x - u + dig_ap1 - dig_b - series[2] / F;

    double cutpoint;
    double lld;
    if (use_probit_link == 1) {
        cutpoint = stan::math::std_normal_log_qf(u);
        if (cutpoint < -fused_stan_math::reference::NORMAL_TAIL_THRESH) {
            double frac = -cutpoint;
            for (int kk = fused_stan_math::reference::MILLS_CF_DEPTH; kk >= 1; --kk) {
                frac = -cutpoint + static_cast<double>(kk) / frac;
            }
            lld = -std::log(frac);
        } else {
            lld = u + 0.5 * cutpoint * cutpoint + CUTPOINT_HALF_LOG_TWO_PI;
        }
    } else {
        cutpoint = u - log_1m_x;
        lld = -log_1m_x;
    }

    popmed_out = reflect ? -cutpoint : cutpoint;
    logJ_out = log_kappa + std::log(score_diff) + CUTPOINT_LOG_HALF + lbeta - a * u - (b - 1.0) * log_1m_x + lld;
    if constexpr (with_partials) {
        const double F_s = series[1], F_r = series[2], F_rr = series[3], F_sr = series[4];
        const double x_F_x = series[5], x_F_rx = series[6];
        const double inv_F2 = 1.0 / (F * F);
        const double x_over_1mx = x / omx;
        const double score_a = u - dig_ap1 + dig_total + (F_s + F_r) / F;
        const double mult = std::exp(CUTPOINT_LOG_HALF + lbeta - a * u - (b - 1.0) * log_1m_x);
        const double du_da = -score_a * mult;
        const double du_db = -(score_a + score_diff) * mult;
        const double link_d = (use_probit_link == 1) ? std::exp(lld) : (1.0 / omx);
        const double d_lld_du = (use_probit_link == 1) ? (1.0 + cutpoint * link_d) : x_over_1mx;
        const double d_diff_da = trig_ap1 - ((F_sr + F_rr) * F - F_r * (F_s + F_r)) * inv_F2;
        const double d_diff_db = -trig_b - (F_sr * F - F_r * F_s) * inv_F2;
        const double d_diff_du = -x_over_1mx - 1.0 - (x_F_rx * F - F_r * x_F_x) * inv_F2;
        const double dH_da = d_diff_da / score_diff + dig_ap1 - 1.0 / a - dig_total - u;
        const double dH_db = d_diff_db / score_diff + dig_b - dig_total - log_1m_x;
        const double dH_du = d_diff_du / score_diff - a + (b - 1.0) * x_over_1mx + d_lld_du;
        d_cut_da_out = link_d * du_da;
        d_cut_db_out = link_d * du_db;
        d_logJ_da_out = dH_da + dH_du * du_da;
        d_logJ_db_out = dH_db + dH_du * du_db;
    }
}

template <bool with_partials>
inline void medians_and_log_J_packed_raw(const double* const* cum_alpha_cls,
                                         const double* alpha_0_cls,
                                         const double* log_kappa_cls,
                                         const int* const* cuts_cls,
                                         const int* m_cls,
                                         const int n_cls,
                                         const int M,
                                         const int use_probit_link,
                                         double* popmed_out,
                                         double* logJ_out,
                                         double* d_cut_da_out,
                                         double* d_cut_db_out,
                                         double* d_logJ_da_out,
                                         double* d_logJ_db_out,
                                         unsigned char* reflect_out,
                                         double* scratch) {
    if (M == 0) return;
    double* a_arr = scratch;
    double* b_arr = a_arr + M;
    double* la_arr = b_arr + M;
    double* lbeta_arr = la_arr + M;
    double* lk_arr = lbeta_arr + M;
    double* dgtot_arr = lk_arr + M;
    double* gamma_in = dgtot_arr + M;
    double* gamma_lg = gamma_in + 2 * M + n_cls;
    double* gamma_dg = gamma_lg + 2 * M + n_cls;
    double* gamma_tg = gamma_dg + 2 * M + n_cls;

    int idx = 0;
    for (int c = 0; c < n_cls; ++c) {
        const double* ca = cum_alpha_cls[c];
        const double a0 = alpha_0_cls[c];
        const double lk = log_kappa_cls[c];
        const int mc = m_cls[c];
        const int* cc = cuts_cls[c];
        for (int j = 0; j < mc; ++j) {
            const double ar = ca[cc[j]];
            const double br = a0 - ar;
            const bool refl = (ar > br);
            reflect_out[idx] = refl ? 1u : 0u;
            const double aj = refl ? br : ar;
            const double bj = refl ? ar : br;
            a_arr[idx] = aj;
            b_arr[idx] = bj;
            la_arr[idx] = std::log(aj);
            lk_arr[idx] = lk;
            gamma_in[idx] = aj + 1.0;
            gamma_in[M + idx] = bj;
            ++idx;
        }
        gamma_in[2 * M + c] = a0;
    }

    for (int i = 0; i < 2 * M + n_cls; ++i) {
        cutpoint_lgamma_digamma_trigamma(gamma_in[i], gamma_lg + i, gamma_dg + i, gamma_tg + i);
    }

    idx = 0;
    for (int c = 0; c < n_cls; ++c) {
        const double lg_tot = gamma_lg[2 * M + c];
        const double dg_tot = gamma_dg[2 * M + c];
        const int mc = m_cls[c];
        for (int j = 0; j < mc; ++j) {
            lbeta_arr[idx] = gamma_lg[idx] - la_arr[idx] + gamma_lg[M + idx] - lg_tot;
            dgtot_arr[idx] = dg_tot;
            ++idx;
        }
    }

    int j0 = 0;
    for (int j = j0; j < M; ++j) {
        double dummy0 = 0.0, dummy1 = 0.0, dummy2 = 0.0, dummy3 = 0.0;
        medians_and_log_J_scalar_one<with_partials>(
            a_arr[j], b_arr[j], la_arr[j], lbeta_arr[j],
            gamma_dg[j], gamma_tg[j], gamma_dg[M + j], gamma_tg[M + j],
            lk_arr[j], dgtot_arr[j], reflect_out[j] != 0, use_probit_link,
            popmed_out[j], logJ_out[j],
            with_partials ? d_cut_da_out[j] : dummy0,
            with_partials ? d_cut_db_out[j] : dummy1,
            with_partials ? d_logJ_da_out[j] : dummy2,
            with_partials ? d_logJ_db_out[j] : dummy3);
    }
}

}  // namespace cutpoint_median_fast_v2
}  // namespace fused_stan_math
