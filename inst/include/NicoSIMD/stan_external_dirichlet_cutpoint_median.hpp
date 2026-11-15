#pragma once

#include <stan/math/mix.hpp>
#include "math/nico_double_probability_reference.hpp"
#include "stan_external_scalar_functions.hpp"
#include "stan_external_std_normal_log_qf.hpp"
#include <cmath>
#include <limits>
#include <ostream>

namespace nico::math {

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
        const double log_shifted = nico_log(shifted, nullptr);

        *lgamma_of_x = (shifted - 0.5) * log_shifted - shifted + CUTPOINT_HALF_LOG_TWO_PI
                       + inverse * (1.0 / 12.0 - inverse_2 * (1.0 / 360.0 - inverse_2 * (1.0 / 1260.0
                       - inverse_2 * (1.0 / 1680.0 - inverse_2 * (1.0 / 1188.0 - inverse_2 * (691.0 / 360360.0
                       - inverse_2 / 156.0))))))
                       - nico_log(product, nullptr);
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

// ---- log median u of Beta(a, b), a <= b (so u <= log(1 / 2)): Newton on
//      g(u) = log I_{exp(u)}(a, b) + log(2) = a u + b log(1 - exp(u)) - log(a) - lbeta(a, b) + log(F) + log(2),
//      g'(u) = a / ((1 - exp(u)) F), kept inside the bracket of the root

inline double cutpoint_beta_median_log(const double a,
                                       const double b,
                                       const double log_a,
                                       const double lbeta) {


        double upper_u = CUTPOINT_LOG_HALF;
        double lower_u = -std::numeric_limits<double>::infinity();
        double u = a >= 1.0 ? nico_log((a - 1.0 / 3.0) / (a + b - 2.0 / 3.0), nullptr)
                            : (CUTPOINT_LOG_HALF + log_a + lbeta) / a;
        if (!(u < upper_u)) {
                u = upper_u;
        }
        double series[7];
        for (int iteration = 0; iteration < 200; ++iteration) {
                const double x = nico_exp(u, nullptr);
                const double log_1m_x = nico_log1m_exp(u, nullptr);
                cutpoint_hypergeometric_series<false>(a + b, a + 1.0, x, series);
                const double g = a * u + b * log_1m_x - log_a - lbeta + nico_log(series[0], nullptr)
                                 - CUTPOINT_LOG_HALF;
                if (g == 0.0) {
                        return u;
                }
                if (g > 0.0) {
                        upper_u = u;
                } else {
                        lower_u = u;
                }
                double next_u = u - g * nico_exp(log_1m_x, nullptr) * series[0] / a;
                if (!(next_u < upper_u) || !(next_u > lower_u)) {
                        next_u = std::isinf(lower_u) ? upper_u - 2.0 * (upper_u - u) - 1.0 : 0.5 * (lower_u + upper_u);
                }
                if (std::abs(next_u - u) <= 4e-16 * (1.0 + std::abs(u))) {
                        return next_u;
                }
                u = next_u;
        }

        return u;

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
        const double log_kappa = nico::math::nico_log(kappa, nullptr);
        const double log_a = nico::math::nico_log(a, nullptr);

        double lgamma_a_plus_1, digamma_a_plus_1, trigamma_a_plus_1;
        double lgamma_b, digamma_b, trigamma_b;
        double lgamma_total, digamma_total;
        nico::math::cutpoint_lgamma_digamma_trigamma(a + 1.0, &lgamma_a_plus_1, &digamma_a_plus_1, &trigamma_a_plus_1);
        nico::math::cutpoint_lgamma_digamma_trigamma(b, &lgamma_b, &digamma_b, &trigamma_b);
        nico::math::cutpoint_lgamma_digamma_trigamma(a + b, &lgamma_total, &digamma_total, nullptr);
        const double lbeta = lgamma_a_plus_1 - log_a + lgamma_b - lgamma_total;

        const double u = nico::math::cutpoint_beta_median_log(a, b, log_a, lbeta);
        const double x = nico::math::nico_exp(u, nullptr);
        const double log_1m_x = nico::math::nico_log1m_exp(u, nullptr);

        double series[7];
        nico::math::cutpoint_hypergeometric_series<with_partials>(a + b, a + 1.0, x, series);
        const double F = series[0];
        const double score_difference = log_1m_x - u + digamma_a_plus_1 - digamma_b - series[2] / F;

        double cutpoint;
        double log_link_derivative;
        if (use_probit_link == 1) {
                cutpoint = nico::math::nico_std_normal_log_qf_value(u, nullptr);
                if (cutpoint < -nico::math::reference::NORMAL_TAIL_THRESH) {
                        double fraction = -cutpoint;
                        for (int k = nico::math::reference::MILLS_CF_DEPTH; k >= 1; --k) {
                                fraction = -cutpoint + static_cast<double>(k) / fraction;
                        }
                        log_link_derivative = -nico::math::nico_log(fraction, nullptr);
                } else {
                        log_link_derivative = u + 0.5 * cutpoint * cutpoint + nico::math::CUTPOINT_HALF_LOG_TWO_PI;
                }
        } else {
                cutpoint = u - log_1m_x;
                log_link_derivative = -log_1m_x;
        }

        values[0] = reflect ? -cutpoint : cutpoint;
        values[1] = log_kappa + nico::math::nico_log(score_difference, nullptr) + nico::math::CUTPOINT_LOG_HALF + lbeta
                    - a * u - (b - 1.0) * log_1m_x + log_link_derivative;

        if constexpr (with_partials) {
                const double F_s = series[1];
                const double F_r = series[2];
                const double F_rr = series[3];
                const double F_sr = series[4];
                const double x_F_x = series[5];
                const double x_F_rx = series[6];
                const double inverse_F_2 = 1.0 / (F * F);
                const double x_over_1m_x = nico::math::nico_exp(u - log_1m_x, nullptr);

                const double score_a = u - digamma_a_plus_1 + digamma_total + (F_s + F_r) / F;
                const double multiplier = nico::math::nico_exp(nico::math::CUTPOINT_LOG_HALF + lbeta - a * u
                                                               - (b - 1.0) * log_1m_x, nullptr);
                const double du_da = -score_a * multiplier;
                const double du_db = -(score_a + score_difference) * multiplier;

                const double link_derivative = nico::math::nico_exp(log_link_derivative, nullptr);
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
#if __has_include("stan_external_lgamma_and_digamma_simd.hpp")
#include "stan_external_lgamma_and_digamma_simd.hpp"
#endif

namespace nico::math {
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

#if !defined(NICO_SIMD_SCALAR_KERNELS) && defined(BMVP_LANES) && ((BMVP_LANES == 8) || (BMVP_LANES == 4))

static constexpr int kLanes = BMVP_LANES;
static constexpr unsigned int kFullMask = (1u << kLanes) - 1u;
using vec_t = nico::math::simd_t;

#if BMVP_LANES == 8
using mask_t = __mmask8;
ALWAYS_INLINE vec_t cm_set1(const double a) { return _mm512_set1_pd(a); }
ALWAYS_INLINE vec_t cm_setzero() { return _mm512_setzero_pd(); }
ALWAYS_INLINE vec_t cm_loadu(const double* p) { return _mm512_loadu_pd(p); }
ALWAYS_INLINE vec_t cm_load(const double* p) { return _mm512_load_pd(p); }
ALWAYS_INLINE void  cm_storeu(double* p, const vec_t v) { _mm512_storeu_pd(p, v); }
ALWAYS_INLINE void  cm_store(double* p, const vec_t v) { _mm512_store_pd(p, v); }
ALWAYS_INLINE vec_t cm_add(const vec_t a, const vec_t b) { return _mm512_add_pd(a, b); }
ALWAYS_INLINE vec_t cm_sub(const vec_t a, const vec_t b) { return _mm512_sub_pd(a, b); }
ALWAYS_INLINE vec_t cm_mul(const vec_t a, const vec_t b) { return _mm512_mul_pd(a, b); }
ALWAYS_INLINE vec_t cm_div(const vec_t a, const vec_t b) { return _mm512_div_pd(a, b); }
ALWAYS_INLINE vec_t cm_abs(const vec_t a) { return _mm512_abs_pd(a); }
ALWAYS_INLINE vec_t cm_fmadd(const vec_t a, const vec_t b, const vec_t c) { return _mm512_fmadd_pd(a, b, c); }
ALWAYS_INLINE vec_t cm_fnmadd(const vec_t a, const vec_t b, const vec_t c) { return _mm512_fnmadd_pd(a, b, c); }
ALWAYS_INLINE mask_t cm_cmp_lt(const vec_t a, const vec_t b) { return _mm512_cmp_pd_mask(a, b, _CMP_LT_OQ); }
ALWAYS_INLINE mask_t cm_cmp_le(const vec_t a, const vec_t b) { return _mm512_cmp_pd_mask(a, b, _CMP_LE_OQ); }
ALWAYS_INLINE mask_t cm_cmp_gt(const vec_t a, const vec_t b) { return _mm512_cmp_pd_mask(a, b, _CMP_GT_OQ); }
ALWAYS_INLINE mask_t cm_cmp_eq(const vec_t a, const vec_t b) { return _mm512_cmp_pd_mask(a, b, _CMP_EQ_OQ); }
ALWAYS_INLINE mask_t cm_mask_prefix(const int count) {
    return (count >= 8) ? static_cast<mask_t>(0xFFu) : static_cast<mask_t>((1u << count) - 1u);
}
ALWAYS_INLINE mask_t cm_mask_and(const mask_t a, const mask_t b) { return static_cast<mask_t>(a & b); }
ALWAYS_INLINE mask_t cm_mask_or(const mask_t a, const mask_t b) { return static_cast<mask_t>(a | b); }
ALWAYS_INLINE mask_t cm_mask_andnot(const mask_t keep, const mask_t drop) { return static_cast<mask_t>(keep & ~drop); }
ALWAYS_INLINE bool   cm_mask_any(const mask_t m) { return m != 0; }
ALWAYS_INLINE bool   cm_mask_all(const mask_t m) { return m == 0xFFu; }
ALWAYS_INLINE bool   cm_mask_test(const mask_t m, const int l) { return (m & (1u << l)) != 0; }
ALWAYS_INLINE vec_t cm_mask_blend(const mask_t m, const vec_t if_false, const vec_t if_true) {
    return _mm512_mask_blend_pd(m, if_false, if_true);
}
ALWAYS_INLINE vec_t cm_mask_mov(const vec_t src, const mask_t m, const vec_t v) {
    return _mm512_mask_mov_pd(src, m, v);
}
ALWAYS_INLINE vec_t cm_mask_add(const vec_t src, const mask_t m, const vec_t a, const vec_t b) {
    return _mm512_mask_add_pd(src, m, a, b);
}
ALWAYS_INLINE vec_t cm_mask_sub(const vec_t src, const mask_t m, const vec_t a, const vec_t b) {
    return _mm512_mask_sub_pd(src, m, a, b);
}
ALWAYS_INLINE vec_t cm_mask_mul(const vec_t src, const mask_t m, const vec_t a, const vec_t b) {
    return _mm512_mask_mul_pd(src, m, a, b);
}
ALWAYS_INLINE vec_t cm_log_wo_checks(const vec_t a) { return nico_log_1_double_precision_wo_checks_AVX512(a); }
ALWAYS_INLINE vec_t cm_log(const vec_t a) { return nico_log_1_double_precision_AVX512(a); }
ALWAYS_INLINE vec_t cm_exp(const vec_t a) { return nico_exp_1_AVX512(a); }
ALWAYS_INLINE vec_t cm_log1m_exp(const vec_t a) { return nico_log1m_exp_1_AVX512(a); }
ALWAYS_INLINE vec_t cm_inv_Phi_from_log_p(const vec_t lp, const vec_t lmp) {
    return nico_inv_Phi_from_log_p_wo_checks_AVX512(lp, lmp);
}
#else  // BMVP_LANES == 4: x86-64 AVX2 natively OR ARM64 NEON via SIMDe
struct mask_t {
    __m256d vec;
    unsigned int bits;
};
ALWAYS_INLINE vec_t cm_set1(const double a) { return _mm256_set1_pd(a); }
ALWAYS_INLINE vec_t cm_setzero() { return _mm256_setzero_pd(); }
ALWAYS_INLINE vec_t cm_loadu(const double* p) { return _mm256_loadu_pd(p); }
ALWAYS_INLINE vec_t cm_load(const double* p) { return _mm256_load_pd(p); }
ALWAYS_INLINE void  cm_storeu(double* p, const vec_t v) { _mm256_storeu_pd(p, v); }
ALWAYS_INLINE void  cm_store(double* p, const vec_t v) { _mm256_store_pd(p, v); }
ALWAYS_INLINE vec_t cm_add(const vec_t a, const vec_t b) { return _mm256_add_pd(a, b); }
ALWAYS_INLINE vec_t cm_sub(const vec_t a, const vec_t b) { return _mm256_sub_pd(a, b); }
ALWAYS_INLINE vec_t cm_mul(const vec_t a, const vec_t b) { return _mm256_mul_pd(a, b); }
ALWAYS_INLINE vec_t cm_div(const vec_t a, const vec_t b) { return _mm256_div_pd(a, b); }
ALWAYS_INLINE vec_t cm_abs(const vec_t a) { return _mm256_andnot_pd(_mm256_set1_pd(-0.0), a); }
ALWAYS_INLINE vec_t cm_fmadd(const vec_t a, const vec_t b, const vec_t c) { return _mm256_fmadd_pd(a, b, c); }
ALWAYS_INLINE vec_t cm_fnmadd(const vec_t a, const vec_t b, const vec_t c) { return _mm256_fnmadd_pd(a, b, c); }
ALWAYS_INLINE mask_t cm_make_cmp(const vec_t v) {
    return mask_t{v, static_cast<unsigned int>(_mm256_movemask_pd(v))};
}
ALWAYS_INLINE mask_t cm_cmp_lt(const vec_t a, const vec_t b) { return cm_make_cmp(_mm256_cmp_pd(a, b, _CMP_LT_OQ)); }
ALWAYS_INLINE mask_t cm_cmp_le(const vec_t a, const vec_t b) { return cm_make_cmp(_mm256_cmp_pd(a, b, _CMP_LE_OQ)); }
ALWAYS_INLINE mask_t cm_cmp_gt(const vec_t a, const vec_t b) { return cm_make_cmp(_mm256_cmp_pd(a, b, _CMP_GT_OQ)); }
ALWAYS_INLINE mask_t cm_cmp_eq(const vec_t a, const vec_t b) { return cm_make_cmp(_mm256_cmp_pd(a, b, _CMP_EQ_OQ)); }
ALWAYS_INLINE mask_t cm_mask_prefix(const int count) {
    const vec_t idx = _mm256_setr_pd(0.0, 1.0, 2.0, 3.0);
    return cm_cmp_lt(idx, _mm256_set1_pd(static_cast<double>(count)));
}
ALWAYS_INLINE mask_t cm_mask_and(const mask_t a, const mask_t b) {
    return mask_t{_mm256_and_pd(a.vec, b.vec), a.bits & b.bits};
}
ALWAYS_INLINE mask_t cm_mask_or(const mask_t a, const mask_t b) {
    return mask_t{_mm256_or_pd(a.vec, b.vec), a.bits | b.bits};
}
ALWAYS_INLINE mask_t cm_mask_andnot(const mask_t keep, const mask_t drop) {
    return mask_t{_mm256_andnot_pd(drop.vec, keep.vec), keep.bits & ~drop.bits};
}
ALWAYS_INLINE bool   cm_mask_any(const mask_t m) { return m.bits != 0u; }
ALWAYS_INLINE bool   cm_mask_all(const mask_t m) { return m.bits == 0xFu; }
ALWAYS_INLINE bool   cm_mask_test(const mask_t m, const int l) { return (m.bits & (1u << l)) != 0u; }
ALWAYS_INLINE vec_t cm_mask_blend(const mask_t m, const vec_t if_false, const vec_t if_true) {
    return _mm256_blendv_pd(if_false, if_true, m.vec);
}
ALWAYS_INLINE vec_t cm_mask_mov(const vec_t src, const mask_t m, const vec_t v) {
    return _mm256_blendv_pd(src, v, m.vec);
}
ALWAYS_INLINE vec_t cm_mask_add(const vec_t src, const mask_t m, const vec_t a, const vec_t b) {
    return _mm256_blendv_pd(src, _mm256_add_pd(a, b), m.vec);
}
ALWAYS_INLINE vec_t cm_mask_sub(const vec_t src, const mask_t m, const vec_t a, const vec_t b) {
    return _mm256_blendv_pd(src, _mm256_sub_pd(a, b), m.vec);
}
ALWAYS_INLINE vec_t cm_mask_mul(const vec_t src, const mask_t m, const vec_t a, const vec_t b) {
    return _mm256_blendv_pd(src, _mm256_mul_pd(a, b), m.vec);
}
ALWAYS_INLINE vec_t cm_log_wo_checks(const vec_t a) { return nico_log_1_double_precision_wo_checks_AVX2(a); }
ALWAYS_INLINE vec_t cm_log(const vec_t a) { return nico_log_1_double_precision_AVX2(a); }
ALWAYS_INLINE vec_t cm_exp(const vec_t a) { return nico_exp_1_AVX2(a); }
ALWAYS_INLINE vec_t cm_log1m_exp(const vec_t a) { return nico_log1m_exp_1_AVX2(a); }
ALWAYS_INLINE vec_t cm_inv_Phi_from_log_p(const vec_t lp, const vec_t lmp) {
    return nico_inv_Phi_from_log_p_wo_checks_AVX2(lp, lmp);
}
#endif

ALWAYS_INLINE void lgamma_digamma_trigamma_block(const vec_t x,
                                                 vec_t& lg_out,
                                                 vec_t& dg_out,
                                                 vec_t& tg_out) {
    const vec_t shift_to = cm_set1(10.0);
    const vec_t one = cm_set1(1.0);
    const vec_t zero = cm_setzero();
    vec_t x_s = x;
    vec_t p = one;
    vec_t dp = zero;
    vec_t tg_shift = zero;

    if (cm_mask_any(cm_cmp_lt(x, shift_to))) {
        vec_t factor[10], d_factor[10];
        for (int k = 0; k < 10; ++k) {
            const vec_t xpk = cm_add(x, cm_set1(static_cast<double>(k)));
            const mask_t below = cm_cmp_lt(xpk, shift_to);
            factor[k] = cm_mask_blend(below, one, xpk);
            d_factor[k] = cm_mask_blend(below, zero, one);
            if (cm_mask_any(below)) {
                const vec_t inv2 = cm_div(one, cm_mul(xpk, xpk));
                tg_shift = cm_mask_add(tg_shift, below, tg_shift, inv2);
            }
        }
        vec_t num_shifts = d_factor[0];
        for (int k = 1; k < 10; ++k) num_shifts = cm_add(num_shifts, d_factor[k]);
        x_s = cm_add(x, num_shifts);
        auto join = [&](const int a, const int b) {
            d_factor[a] = cm_fmadd(d_factor[a], factor[b], cm_mul(factor[a], d_factor[b]));
            factor[a] = cm_mul(factor[a], factor[b]);
        };
        join(0, 1); join(2, 3); join(4, 5); join(6, 7); join(8, 9);
        join(0, 2); join(4, 6);
        join(0, 4);
        join(0, 8);
        p = factor[0];
        dp = d_factor[0];
    }

    const vec_t log_x_s = cm_log_wo_checks(x_s);
    const vec_t log_p = cm_log_wo_checks(p);
    const vec_t z = cm_div(one, x_s);
    const vec_t w = cm_mul(z, z);

    vec_t lgs = cm_set1(1.0 / 156.0);
    lgs = cm_fmadd(lgs, w, cm_set1(-691.0 / 360360.0));
    lgs = cm_fmadd(lgs, w, cm_set1(1.0 / 1188.0));
    lgs = cm_fmadd(lgs, w, cm_set1(-1.0 / 1680.0));
    lgs = cm_fmadd(lgs, w, cm_set1(1.0 / 1260.0));
    lgs = cm_fmadd(lgs, w, cm_set1(-1.0 / 360.0));
    lgs = cm_fmadd(lgs, w, cm_set1(1.0 / 12.0));
    lgs = cm_mul(lgs, z);
    lg_out = cm_sub(
        cm_add(cm_sub(cm_mul(cm_sub(x_s, cm_set1(0.5)), log_x_s), x_s),
               cm_add(cm_set1(CUTPOINT_HALF_LOG_TWO_PI), lgs)),
        log_p);

    vec_t dgs = cm_set1(1.0 / 12.0);
    dgs = cm_fmadd(dgs, w, cm_set1(-691.0 / 32760.0));
    dgs = cm_fmadd(dgs, w, cm_set1(1.0 / 132.0));
    dgs = cm_fmadd(dgs, w, cm_set1(-1.0 / 240.0));
    dgs = cm_fmadd(dgs, w, cm_set1(1.0 / 252.0));
    dgs = cm_fmadd(dgs, w, cm_set1(-1.0 / 120.0));
    dgs = cm_fmadd(dgs, w, cm_set1(1.0 / 12.0));
    dgs = cm_mul(dgs, w);
    dg_out = cm_sub(cm_sub(log_x_s, cm_fmadd(cm_set1(0.5), z, dgs)),
                    cm_div(dp, p));

    vec_t tgs = cm_set1(7.0 / 6.0);
    tgs = cm_fmadd(tgs, w, cm_set1(-691.0 / 2730.0));
    tgs = cm_fmadd(tgs, w, cm_set1(5.0 / 66.0));
    tgs = cm_fmadd(tgs, w, cm_set1(-1.0 / 30.0));
    tgs = cm_fmadd(tgs, w, cm_set1(1.0 / 42.0));
    tgs = cm_fmadd(tgs, w, cm_set1(-1.0 / 30.0));
    tgs = cm_fmadd(tgs, w, cm_set1(1.0 / 6.0));
    const vec_t tg_asym = cm_fmadd(cm_mul(tgs, z), w,
                                   cm_fmadd(cm_set1(0.5), w, z));
    tg_out = cm_add(tg_asym, tg_shift);
}

inline void cutpoint_lgamma_digamma_trigamma_simd(const double* x,
                                                  const int n,
                                                  double* lg_out,
                                                  double* dg_out,
                                                  double* tg_out) {
    int first = 0;
    const vec_t smallest = cm_set1(1e-300);
    const vec_t largest = cm_set1(1e300);
    for (; first + kLanes <= n; first += kLanes) {
        const vec_t xb = cm_loadu(x + first);
        const bool all_ok = cm_mask_all(cm_mask_and(cm_cmp_lt(smallest, xb), cm_cmp_lt(xb, largest)));
        if (!all_ok) {
            for (int i = first; i < first + kLanes; ++i) {
                cutpoint_lgamma_digamma_trigamma(x[i], lg_out + i, dg_out + i, tg_out + i);
            }
            continue;
        }
        vec_t lgb, dgb, tgb;
        lgamma_digamma_trigamma_block(xb, lgb, dgb, tgb);
        cm_storeu(lg_out + first, lgb);
        cm_storeu(dg_out + first, dgb);
        cm_storeu(tg_out + first, tgb);
    }
    const int rem = n - first;
    if (rem >= (kLanes / 2)) {
        alignas(64) double xp[kLanes], lgp[kLanes], dgp[kLanes], tgp[kLanes];
        bool all_ok = true;
        for (int l = 0; l < kLanes; ++l) {
            xp[l] = (l < rem) ? x[first + l] : 1.0;
            all_ok &= (xp[l] > 1e-300 && xp[l] < 1e300);
        }
        if (all_ok) {
            vec_t lgb, dgb, tgb;
            lgamma_digamma_trigamma_block(cm_load(xp), lgb, dgb, tgb);
            cm_store(lgp, lgb);
            cm_store(dgp, dgb);
            cm_store(tgp, tgb);
            for (int l = 0; l < rem; ++l) {
                lg_out[first + l] = lgp[l];
                dg_out[first + l] = dgp[l];
                tg_out[first + l] = tgp[l];
            }
            return;
        }
    }
    for (int i = first; i < n; ++i) {
        cutpoint_lgamma_digamma_trigamma(x[i], lg_out + i, dg_out + i, tg_out + i);
    }
}

template <bool with_partials>
inline void medians_and_log_J_block(const double* a_in,
                                    const double* b_in,
                                    const double* log_a_in,
                                    const double* lbeta_in,
                                    const double* dig_ap1_in,
                                    const double* trig_ap1_in,
                                    const double* dig_b_in,
                                    const double* trig_b_in,
                                    const double log_kappa,
                                    const double dig_total,
                                    const unsigned char* reflect_in,
                                    const int count,
                                    const int use_probit_link,
                                    double* popmed_out,
                                    double* logJ_out,
                                    double* d_cut_da_out,
                                    double* d_cut_db_out,
                                    double* d_logJ_da_out,
                                    double* d_logJ_db_out) {
    alignas(64) double ap[kLanes], bp[kLanes], lap[kLanes], lbp[kLanes], u_init[kLanes];
    const mask_t valid_mask = cm_mask_prefix(count);
    for (int l = 0; l < kLanes; ++l) {
        if (l < count) {
            ap[l] = a_in[l];
            bp[l] = b_in[l];
            lap[l] = log_a_in[l];
            lbp[l] = lbeta_in[l];
            const double u0 = (ap[l] >= 1.0)
                ? std::log((ap[l] - 1.0 / 3.0) / (ap[l] + bp[l] - 2.0 / 3.0))
                : ((CUTPOINT_LOG_HALF + lap[l] + lbp[l]) / ap[l]);
            u_init[l] = (u0 < CUTPOINT_LOG_HALF) ? u0 : CUTPOINT_LOG_HALF;
        } else {
            ap[l] = 1.0;
            bp[l] = 1.0;
            lap[l] = 0.0;
            lbp[l] = 0.0;
            u_init[l] = CUTPOINT_LOG_HALF;
        }
    }

    const vec_t a_v = cm_load(ap);
    const vec_t b_v = cm_load(bp);
    const vec_t la_v = cm_load(lap);
    const vec_t lb_v = cm_load(lbp);
    const vec_t s_v = cm_add(a_v, b_v);
    const vec_t r_v = cm_add(a_v, cm_set1(1.0));
    const vec_t bm1_v = cm_sub(b_v, cm_set1(1.0));
    const vec_t one = cm_set1(1.0);
    const vec_t zero = cm_setzero();
    const vec_t log_half = cm_set1(CUTPOINT_LOG_HALF);

    vec_t u_v = cm_load(u_init);
    vec_t lower_u = cm_set1(-std::numeric_limits<double>::infinity());
    vec_t upper_u = log_half;
    mask_t outer_mask = valid_mask;

    for (int iter = 0; iter < 200 && cm_mask_any(outer_mask); ++iter) {
        const vec_t x_v = cm_exp(u_v);
        const vec_t omx_v = cm_sub(one, x_v);
        const vec_t log_1m_x_v = cm_log1m_exp(u_v);

        vec_t F_v = zero;
        vec_t term_v = one;
        mask_t inner_mask = outer_mask;
        for (int n = 0; n < 1000000 && cm_mask_any(inner_mask); ++n) {
            F_v = cm_mask_add(F_v, inner_mask, F_v, term_v);
            const vec_t nv = cm_set1(static_cast<double>(n));
            const vec_t sn = cm_add(s_v, nv);
            const vec_t rn = cm_add(r_v, nv);
            const vec_t lhs = cm_mul(sn, x_v);
            const mask_t r_ok = cm_cmp_lt(lhs, rn);
            const vec_t bound = cm_mul(cm_set1(1e-17), cm_mul(F_v, cm_sub(rn, lhs)));
            const mask_t t_ok = cm_cmp_lt(cm_mul(term_v, rn), bound);
            inner_mask = cm_mask_andnot(inner_mask, cm_mask_and(r_ok, t_ok));
            if (!cm_mask_any(inner_mask)) break;
            term_v = cm_mask_mul(zero, inner_mask, term_v, cm_div(lhs, rn));
        }

        const vec_t log_F = cm_log(F_v);
        const vec_t g_v = cm_sub(
            cm_add(cm_sub(cm_sub(cm_fmadd(a_v, u_v, cm_mul(b_v, log_1m_x_v)),
                                 la_v),
                          lb_v),
                   log_F),
            log_half);

        const mask_t g_zero = cm_cmp_eq(g_v, zero);
        const mask_t g_pos = cm_cmp_gt(g_v, zero);
        const mask_t g_neg = cm_mask_and(outer_mask, cm_cmp_lt(g_v, zero));
        upper_u = cm_mask_mov(upper_u, cm_mask_and(outer_mask, g_pos), u_v);
        lower_u = cm_mask_mov(lower_u, g_neg, u_v);

        const vec_t phi_v = cm_div(a_v, cm_mul(omx_v, F_v));
        const vec_t g_over_gp = cm_div(g_v, phi_v);
        const vec_t dphi_r = cm_sub(cm_sub(a_v, cm_mul(bm1_v, cm_div(x_v, omx_v))),
                                    phi_v);
        const vec_t hcorr = cm_fnmadd(cm_set1(0.5), cm_mul(g_over_gp, dphi_r), one);
        const mask_t use_halley = cm_cmp_gt(hcorr, cm_set1(0.25));
        const vec_t step_h = cm_sub(u_v, cm_div(g_over_gp, hcorr));
        const vec_t step_n = cm_sub(u_v, g_over_gp);
        vec_t next_u = cm_mask_blend(use_halley, step_n, step_h);

        const mask_t in_hi = cm_cmp_lt(next_u, upper_u);
        const mask_t in_lo = cm_cmp_gt(next_u, lower_u);
        const mask_t out_bracket = cm_mask_andnot(outer_mask, cm_mask_and(in_hi, in_lo));
        if (cm_mask_any(out_bracket)) {
            const vec_t inf_v = cm_set1(-std::numeric_limits<double>::infinity());
            const mask_t lo_inf = cm_cmp_eq(lower_u, inf_v);
            const vec_t fallback_inf = cm_sub(
                cm_sub(upper_u, cm_mul(cm_set1(2.0), cm_sub(upper_u, u_v))),
                one);
            const vec_t fallback_mid = cm_mul(cm_set1(0.5), cm_add(lower_u, upper_u));
            const vec_t fb = cm_mask_blend(lo_inf, fallback_mid, fallback_inf);
            next_u = cm_mask_mov(next_u, out_bracket, fb);
        }

        const vec_t tol = cm_mul(cm_set1(4e-16), cm_add(one, cm_abs(u_v)));
        const mask_t conv = cm_cmp_le(cm_abs(cm_sub(next_u, u_v)), tol);
        u_v = cm_mask_mov(u_v, cm_mask_andnot(outer_mask, g_zero), next_u);
        outer_mask = cm_mask_andnot(outer_mask, cm_mask_or(g_zero, conv));
    }

    if (cm_mask_any(outer_mask)) {
        alignas(64) double u_tmp[kLanes];
        cm_store(u_tmp, u_v);
        for (int l = 0; l < count; ++l) {
            if (cm_mask_test(outer_mask, l)) {
                u_tmp[l] = cutpoint_beta_median_log_halley_scalar(ap[l], bp[l], lap[l], lbp[l]);
            }
        }
        u_v = cm_load(u_tmp);
    }

    const vec_t x_v = cm_exp(u_v);
    const vec_t log_1m_x_v = cm_log1m_exp(u_v);
    vec_t F_v = zero, Fs_v = zero, Fr_v = zero, Frr_v = zero, Fsr_v = zero, xFx_v = zero, xFrx_v = zero;
    vec_t sum_inv_s = zero, sum_inv_r = zero, sum_inv_r2 = zero;
    vec_t term_v = one;
    mask_t smask = valid_mask;
    for (int n = 0; n < 1000000 && cm_mask_any(smask); ++n) {
        F_v = cm_mask_add(F_v, smask, F_v, term_v);
        Fr_v = cm_mask_sub(Fr_v, smask, Fr_v, cm_mul(term_v, sum_inv_r));
        if constexpr (with_partials) {
            Fs_v = cm_mask_add(Fs_v, smask, Fs_v, cm_mul(term_v, sum_inv_s));
            const vec_t r_w = cm_fmadd(sum_inv_r, sum_inv_r, sum_inv_r2);
            Frr_v = cm_mask_add(Frr_v, smask, Frr_v, cm_mul(term_v, r_w));
            Fsr_v = cm_mask_sub(Fsr_v, smask, Fsr_v, cm_mul(term_v, cm_mul(sum_inv_s, sum_inv_r)));
            const vec_t nv = cm_set1(static_cast<double>(n));
            xFx_v = cm_mask_add(xFx_v, smask, xFx_v, cm_mul(term_v, nv));
            xFrx_v = cm_mask_sub(xFrx_v, smask, xFrx_v, cm_mul(cm_mul(term_v, sum_inv_r), nv));
        }
        const vec_t nv = cm_set1(static_cast<double>(n));
        const vec_t sn = cm_add(s_v, nv);
        const vec_t rn = cm_add(r_v, nv);
        const vec_t w_chk = with_partials
            ? cm_add(
                  cm_add(one, cm_mul(nv, cm_add(one, sum_inv_r))),
                  cm_add(cm_mul(sum_inv_s, sum_inv_r),
                         cm_fmadd(sum_inv_r, sum_inv_r, sum_inv_r2)))
            : cm_add(one, sum_inv_r);
        const mask_t r_ok = cm_cmp_lt(cm_mul(sn, x_v), rn);
        const mask_t t_ok = cm_cmp_lt(cm_mul(term_v, w_chk), cm_mul(cm_set1(1e-17), F_v));
        smask = cm_mask_andnot(smask, cm_mask_and(r_ok, t_ok));
        if (!cm_mask_any(smask)) break;
        const vec_t inv_rn = cm_div(one, rn);
        if constexpr (with_partials) {
            sum_inv_s = cm_mask_add(sum_inv_s, smask, sum_inv_s, cm_div(one, sn));
            sum_inv_r2 = cm_mask_add(sum_inv_r2, smask, sum_inv_r2, cm_mul(inv_rn, inv_rn));
        }
        sum_inv_r = cm_mask_add(sum_inv_r, smask, sum_inv_r, inv_rn);
        term_v = cm_mask_mul(zero, smask, term_v, cm_mul(cm_mul(sn, inv_rn), x_v));
    }

    vec_t cp_v, lld_v;
    if (use_probit_link == 1) {
        cp_v = cm_inv_Phi_from_log_p(u_v, log_1m_x_v);
        lld_v = cm_add(cm_add(u_v, cm_set1(CUTPOINT_HALF_LOG_TWO_PI)),
                       cm_mul(cm_set1(0.5), cm_mul(cp_v, cp_v)));
        const mask_t deep_tail = cm_mask_and(
            valid_mask, cm_cmp_lt(cp_v, cm_set1(-nico::math::reference::NORMAL_TAIL_THRESH)));
        if (cm_mask_any(deep_tail)) {
            alignas(64) double ctmp[kLanes], ltmp[kLanes];
            cm_store(ctmp, cp_v);
            cm_store(ltmp, lld_v);
            for (int l = 0; l < count; ++l) {
                if (cm_mask_test(deep_tail, l)) {
                    double frac = -ctmp[l];
                    for (int kk = nico::math::reference::MILLS_CF_DEPTH; kk >= 1; --kk) {
                        frac = -ctmp[l] + static_cast<double>(kk) / frac;
                    }
                    ltmp[l] = -std::log(frac);
                }
            }
            lld_v = cm_load(ltmp);
        }
    } else {
        cp_v = cm_sub(u_v, log_1m_x_v);
        lld_v = cm_sub(zero, log_1m_x_v);
    }

    alignas(64) double u_a[kLanes], l1mx_a[kLanes], F_a[kLanes], Fr_a[kLanes], cp_a[kLanes], lld_a[kLanes];
    alignas(64) double Fs_a[kLanes], Frr_a[kLanes], Fsr_a[kLanes], xFx_a[kLanes], xFrx_a[kLanes];
    cm_store(u_a, u_v);
    cm_store(l1mx_a, log_1m_x_v);
    cm_store(F_a, F_v);
    cm_store(Fr_a, Fr_v);
    cm_store(cp_a, cp_v);
    cm_store(lld_a, lld_v);
    if constexpr (with_partials) {
        cm_store(Fs_a, Fs_v);
        cm_store(Frr_a, Frr_v);
        cm_store(Fsr_a, Fsr_v);
        cm_store(xFx_a, xFx_v);
        cm_store(xFrx_a, xFrx_v);
    }

    for (int l = 0; l < count; ++l) {
        const double u = u_a[l];
        const double log_1m_x = l1mx_a[l];
        const double F = F_a[l];
        const double score_diff = log_1m_x - u + dig_ap1_in[l] - dig_b_in[l] - Fr_a[l] / F;
        const double cutpoint = cp_a[l];
        const double lld = lld_a[l];
        popmed_out[l] = reflect_in[l] ? -cutpoint : cutpoint;
        logJ_out[l] = log_kappa + std::log(score_diff) + CUTPOINT_LOG_HALF + lbp[l]
                      - ap[l] * u - (bp[l] - 1.0) * log_1m_x + lld;
        if constexpr (with_partials) {
            const double inv_F2 = 1.0 / (F * F);
            const double x_over_1mx = std::exp(u - log_1m_x);
            const double score_a = u - dig_ap1_in[l] + dig_total + (Fs_a[l] + Fr_a[l]) / F;
            const double mult = std::exp(CUTPOINT_LOG_HALF + lbp[l] - ap[l] * u - (bp[l] - 1.0) * log_1m_x);
            const double du_da = -score_a * mult;
            const double du_db = -(score_a + score_diff) * mult;
            const double link_d = std::exp(lld);
            const double d_lld_du = (use_probit_link == 1) ? (1.0 + cutpoint * link_d) : x_over_1mx;
            const double d_diff_da = trig_ap1_in[l] - ((Fsr_a[l] + Frr_a[l]) * F - Fr_a[l] * (Fs_a[l] + Fr_a[l])) * inv_F2;
            const double d_diff_db = -trig_b_in[l] - (Fsr_a[l] * F - Fr_a[l] * Fs_a[l]) * inv_F2;
            const double d_diff_du = -x_over_1mx - 1.0 - (xFrx_a[l] * F - Fr_a[l] * xFx_a[l]) * inv_F2;
            const double dH_da = d_diff_da / score_diff + dig_ap1_in[l] - 1.0 / ap[l] - dig_total - u;
            const double dH_db = d_diff_db / score_diff + dig_b_in[l] - dig_total - log_1m_x;
            const double dH_du = d_diff_du / score_diff - ap[l] + (bp[l] - 1.0) * x_over_1mx + d_lld_du;
            d_cut_da_out[l] = link_d * du_da;
            d_cut_db_out[l] = link_d * du_db;
            d_logJ_da_out[l] = dH_da + dH_du * du_da;
            d_logJ_db_out[l] = dH_db + dH_du * du_db;
        }
    }
}

#endif  // BMVP_LANES == 8 || BMVP_LANES == 4

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
        cutpoint = nico::math::nico_std_normal_log_qf_value(u, nullptr);
        if (cutpoint < -nico::math::reference::NORMAL_TAIL_THRESH) {
            double frac = -cutpoint;
            for (int kk = nico::math::reference::MILLS_CF_DEPTH; kk >= 1; --kk) {
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
    namespace cmf = nico::math::cutpoint_median_fast;
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

#if !defined(NICO_SIMD_SCALAR_KERNELS) && defined(BMVP_LANES) && ((BMVP_LANES == 8) || (BMVP_LANES == 4))
    cmf::cutpoint_lgamma_digamma_trigamma_simd(gamma_in, 2 * m + 1, gamma_lg, gamma_dg, gamma_tg);
#else
    for (int j = 0; j < 2 * m + 1; ++j) {
        cutpoint_lgamma_digamma_trigamma(gamma_in[j], gamma_lg + j, gamma_dg + j, gamma_tg + j);
    }
#endif
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
#if !defined(NICO_SIMD_SCALAR_KERNELS) && defined(BMVP_LANES) && ((BMVP_LANES == 8) || (BMVP_LANES == 4))
    for (; j0 + cmf::kLanes <= m; j0 += cmf::kLanes) {
        cmf::medians_and_log_J_block<is_any_var>(
            a_arr + j0, b_arr + j0, log_a_arr + j0, lbeta_arr + j0,
            gamma_dg + j0, gamma_tg + j0, gamma_dg + m + j0, gamma_tg + m + j0,
            log_kappa, dig_total, reflect + j0, cmf::kLanes, use_probit_link,
            popmed_vals + j0, logJ_vals + j0,
            is_any_var ? (d_cut_da + j0) : nullptr,
            is_any_var ? (d_cut_db + j0) : nullptr,
            is_any_var ? (d_logJ_da + j0) : nullptr,
            is_any_var ? (d_logJ_db + j0) : nullptr);
    }
    if (m - j0 >= (cmf::kLanes / 2)) {
        const int rem = m - j0;
        cmf::medians_and_log_J_block<is_any_var>(
            a_arr + j0, b_arr + j0, log_a_arr + j0, lbeta_arr + j0,
            gamma_dg + j0, gamma_tg + j0, gamma_dg + m + j0, gamma_tg + m + j0,
            log_kappa, dig_total, reflect + j0, rem, use_probit_link,
            popmed_vals + j0, logJ_vals + j0,
            is_any_var ? (d_cut_da + j0) : nullptr,
            is_any_var ? (d_cut_db + j0) : nullptr,
            is_any_var ? (d_logJ_da + j0) : nullptr,
            is_any_var ? (d_logJ_db + j0) : nullptr);
        j0 = m;
    }
#endif
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

}  // namespace nico::math

using nico::math::nico_dirichlet_cutpoint_medians_and_log_J_at_cuts;


// ---- Hybrid CF/Series multi-class packed Dirichlet cutpoint medians & log-Jacobians ----
namespace nico::math {
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

#if !defined(NICO_SIMD_SCALAR_KERNELS) && defined(BMVP_LANES) && ((BMVP_LANES == 8) || (BMVP_LANES == 4))

static constexpr int kLanes = BMVP_LANES;
static constexpr unsigned int kFullMask = (1u << kLanes) - 1u;
using vec_t = nico::math::simd_t;

#if BMVP_LANES == 8
using mask_t = __mmask8;
ALWAYS_INLINE vec_t cm_set1(const double a) { return _mm512_set1_pd(a); }
ALWAYS_INLINE vec_t cm_setzero() { return _mm512_setzero_pd(); }
ALWAYS_INLINE vec_t cm_loadu(const double* p) { return _mm512_loadu_pd(p); }
ALWAYS_INLINE vec_t cm_load(const double* p) { return _mm512_load_pd(p); }
ALWAYS_INLINE void  cm_storeu(double* p, const vec_t v) { _mm512_storeu_pd(p, v); }
ALWAYS_INLINE void  cm_store(double* p, const vec_t v) { _mm512_store_pd(p, v); }
ALWAYS_INLINE vec_t cm_add(const vec_t a, const vec_t b) { return _mm512_add_pd(a, b); }
ALWAYS_INLINE vec_t cm_sub(const vec_t a, const vec_t b) { return _mm512_sub_pd(a, b); }
ALWAYS_INLINE vec_t cm_mul(const vec_t a, const vec_t b) { return _mm512_mul_pd(a, b); }
ALWAYS_INLINE vec_t cm_div(const vec_t a, const vec_t b) { return _mm512_div_pd(a, b); }
ALWAYS_INLINE vec_t cm_sqrt(const vec_t a) { return _mm512_sqrt_pd(a); }
ALWAYS_INLINE vec_t cm_abs(const vec_t a) { return _mm512_abs_pd(a); }
ALWAYS_INLINE vec_t cm_fmadd(const vec_t a, const vec_t b, const vec_t c) { return _mm512_fmadd_pd(a, b, c); }
ALWAYS_INLINE vec_t cm_fnmadd(const vec_t a, const vec_t b, const vec_t c) { return _mm512_fnmadd_pd(a, b, c); }
ALWAYS_INLINE mask_t cm_cmp_lt(const vec_t a, const vec_t b) { return _mm512_cmp_pd_mask(a, b, _CMP_LT_OQ); }
ALWAYS_INLINE mask_t cm_cmp_le(const vec_t a, const vec_t b) { return _mm512_cmp_pd_mask(a, b, _CMP_LE_OQ); }
ALWAYS_INLINE mask_t cm_cmp_gt(const vec_t a, const vec_t b) { return _mm512_cmp_pd_mask(a, b, _CMP_GT_OQ); }
ALWAYS_INLINE mask_t cm_cmp_eq(const vec_t a, const vec_t b) { return _mm512_cmp_pd_mask(a, b, _CMP_EQ_OQ); }
ALWAYS_INLINE mask_t cm_mask_prefix(const int count) {
    return (count >= 8) ? static_cast<mask_t>(0xFFu) : static_cast<mask_t>((1u << count) - 1u);
}
ALWAYS_INLINE mask_t cm_mask_and(const mask_t a, const mask_t b) { return static_cast<mask_t>(a & b); }
ALWAYS_INLINE mask_t cm_mask_or(const mask_t a, const mask_t b) { return static_cast<mask_t>(a | b); }
ALWAYS_INLINE mask_t cm_mask_andnot(const mask_t keep, const mask_t drop) { return static_cast<mask_t>(keep & ~drop); }
ALWAYS_INLINE bool   cm_mask_any(const mask_t m) { return m != 0; }
ALWAYS_INLINE bool   cm_mask_all(const mask_t m) { return m == 0xFFu; }
ALWAYS_INLINE bool   cm_mask_test(const mask_t m, const int l) { return (m & (1u << l)) != 0; }
ALWAYS_INLINE vec_t cm_mask_blend(const mask_t m, const vec_t if_false, const vec_t if_true) {
    return _mm512_mask_blend_pd(m, if_false, if_true);
}
ALWAYS_INLINE vec_t cm_mask_mov(const vec_t src, const mask_t m, const vec_t v) {
    return _mm512_mask_mov_pd(src, m, v);
}
ALWAYS_INLINE vec_t cm_mask_add(const vec_t src, const mask_t m, const vec_t a, const vec_t b) {
    return _mm512_mask_add_pd(src, m, a, b);
}
ALWAYS_INLINE vec_t cm_mask_sub(const vec_t src, const mask_t m, const vec_t a, const vec_t b) {
    return _mm512_mask_sub_pd(src, m, a, b);
}
ALWAYS_INLINE vec_t cm_mask_mul(const vec_t src, const mask_t m, const vec_t a, const vec_t b) {
    return _mm512_mask_mul_pd(src, m, a, b);
}
ALWAYS_INLINE vec_t cm_log_wo_checks(const vec_t a) { return nico_log_1_double_precision_wo_checks_AVX512(a); }
ALWAYS_INLINE vec_t cm_log(const vec_t a) { return nico_log_1_double_precision_AVX512(a); }
ALWAYS_INLINE vec_t cm_exp(const vec_t a) { return nico_exp_1_AVX512(a); }
ALWAYS_INLINE vec_t cm_log1m_exp(const vec_t a) { return nico_log1m_exp_1_AVX512(a); }
ALWAYS_INLINE vec_t cm_inv_Phi_from_log_p(const vec_t lp, const vec_t lmp) {
    return nico_inv_Phi_from_log_p_wo_checks_AVX512(lp, lmp);
}
#else  // BMVP_LANES == 4
struct mask_t {
    __m256d vec;
    unsigned int bits;
};
ALWAYS_INLINE vec_t cm_set1(const double a) { return _mm256_set1_pd(a); }
ALWAYS_INLINE vec_t cm_setzero() { return _mm256_setzero_pd(); }
ALWAYS_INLINE vec_t cm_loadu(const double* p) { return _mm256_loadu_pd(p); }
ALWAYS_INLINE vec_t cm_load(const double* p) { return _mm256_load_pd(p); }
ALWAYS_INLINE void  cm_storeu(double* p, const vec_t v) { _mm256_storeu_pd(p, v); }
ALWAYS_INLINE void  cm_store(double* p, const vec_t v) { _mm256_store_pd(p, v); }
ALWAYS_INLINE vec_t cm_add(const vec_t a, const vec_t b) { return _mm256_add_pd(a, b); }
ALWAYS_INLINE vec_t cm_sub(const vec_t a, const vec_t b) { return _mm256_sub_pd(a, b); }
ALWAYS_INLINE vec_t cm_mul(const vec_t a, const vec_t b) { return _mm256_mul_pd(a, b); }
ALWAYS_INLINE vec_t cm_div(const vec_t a, const vec_t b) { return _mm256_div_pd(a, b); }
ALWAYS_INLINE vec_t cm_sqrt(const vec_t a) { return _mm256_sqrt_pd(a); }
ALWAYS_INLINE vec_t cm_abs(const vec_t a) { return _mm256_andnot_pd(_mm256_set1_pd(-0.0), a); }
ALWAYS_INLINE vec_t cm_fmadd(const vec_t a, const vec_t b, const vec_t c) { return _mm256_fmadd_pd(a, b, c); }
ALWAYS_INLINE vec_t cm_fnmadd(const vec_t a, const vec_t b, const vec_t c) { return _mm256_fnmadd_pd(a, b, c); }
ALWAYS_INLINE mask_t cm_make_cmp(const vec_t v) {
    return mask_t{v, static_cast<unsigned int>(_mm256_movemask_pd(v))};
}
ALWAYS_INLINE mask_t cm_cmp_lt(const vec_t a, const vec_t b) { return cm_make_cmp(_mm256_cmp_pd(a, b, _CMP_LT_OQ)); }
ALWAYS_INLINE mask_t cm_cmp_le(const vec_t a, const vec_t b) { return cm_make_cmp(_mm256_cmp_pd(a, b, _CMP_LE_OQ)); }
ALWAYS_INLINE mask_t cm_cmp_gt(const vec_t a, const vec_t b) { return cm_make_cmp(_mm256_cmp_pd(a, b, _CMP_GT_OQ)); }
ALWAYS_INLINE mask_t cm_cmp_eq(const vec_t a, const vec_t b) { return cm_make_cmp(_mm256_cmp_pd(a, b, _CMP_EQ_OQ)); }
ALWAYS_INLINE mask_t cm_mask_prefix(const int count) {
    const vec_t idx = _mm256_setr_pd(0.0, 1.0, 2.0, 3.0);
    return cm_cmp_lt(idx, _mm256_set1_pd(static_cast<double>(count)));
}
ALWAYS_INLINE mask_t cm_mask_and(const mask_t a, const mask_t b) {
    return mask_t{_mm256_and_pd(a.vec, b.vec), a.bits & b.bits};
}
ALWAYS_INLINE mask_t cm_mask_or(const mask_t a, const mask_t b) {
    return mask_t{_mm256_or_pd(a.vec, b.vec), a.bits | b.bits};
}
ALWAYS_INLINE mask_t cm_mask_andnot(const mask_t keep, const mask_t drop) {
    return mask_t{_mm256_andnot_pd(drop.vec, keep.vec), keep.bits & ~drop.bits};
}
ALWAYS_INLINE bool   cm_mask_any(const mask_t m) { return m.bits != 0u; }
ALWAYS_INLINE bool   cm_mask_all(const mask_t m) { return m.bits == 0xFu; }
ALWAYS_INLINE bool   cm_mask_test(const mask_t m, const int l) { return (m.bits & (1u << l)) != 0u; }
ALWAYS_INLINE vec_t cm_mask_blend(const mask_t m, const vec_t if_false, const vec_t if_true) {
    return _mm256_blendv_pd(if_false, if_true, m.vec);
}
ALWAYS_INLINE vec_t cm_mask_mov(const vec_t src, const mask_t m, const vec_t v) {
    return _mm256_blendv_pd(src, v, m.vec);
}
ALWAYS_INLINE vec_t cm_mask_add(const vec_t src, const mask_t m, const vec_t a, const vec_t b) {
    return _mm256_blendv_pd(src, _mm256_add_pd(a, b), m.vec);
}
ALWAYS_INLINE vec_t cm_mask_sub(const vec_t src, const mask_t m, const vec_t a, const vec_t b) {
    return _mm256_blendv_pd(src, _mm256_sub_pd(a, b), m.vec);
}
ALWAYS_INLINE vec_t cm_mask_mul(const vec_t src, const mask_t m, const vec_t a, const vec_t b) {
    return _mm256_blendv_pd(src, _mm256_mul_pd(a, b), m.vec);
}
ALWAYS_INLINE vec_t cm_log_wo_checks(const vec_t a) { return nico_log_1_double_precision_wo_checks_AVX2(a); }
ALWAYS_INLINE vec_t cm_log(const vec_t a) { return nico_log_1_double_precision_AVX2(a); }
ALWAYS_INLINE vec_t cm_exp(const vec_t a) { return nico_exp_1_AVX2(a); }
ALWAYS_INLINE vec_t cm_log1m_exp(const vec_t a) { return nico_log1m_exp_1_AVX2(a); }
ALWAYS_INLINE vec_t cm_inv_Phi_from_log_p(const vec_t lp, const vec_t lmp) {
    return nico_inv_Phi_from_log_p_wo_checks_AVX2(lp, lmp);
}
#endif

ALWAYS_INLINE void lgamma_digamma_trigamma_block(const vec_t x,
                                                 vec_t& lg_out,
                                                 vec_t& dg_out,
                                                 vec_t& tg_out) {
    const vec_t shift_to = cm_set1(10.0);
    const vec_t one = cm_set1(1.0);
    const vec_t zero = cm_setzero();
    vec_t x_s = x;
    vec_t p = one;
    vec_t dp = zero;
    vec_t tg_shift = zero;

    if (cm_mask_any(cm_cmp_lt(x, shift_to))) {
        vec_t factor[10], d_factor[10];
        for (int k = 0; k < 10; ++k) {
            const vec_t xpk = cm_add(x, cm_set1(static_cast<double>(k)));
            const mask_t below = cm_cmp_lt(xpk, shift_to);
            factor[k] = cm_mask_blend(below, one, xpk);
            d_factor[k] = cm_mask_blend(below, zero, one);
            if (cm_mask_any(below)) {
                const vec_t inv2 = cm_div(one, cm_mul(xpk, xpk));
                tg_shift = cm_mask_add(tg_shift, below, tg_shift, inv2);
            }
        }
        vec_t num_shifts = d_factor[0];
        for (int k = 1; k < 10; ++k) num_shifts = cm_add(num_shifts, d_factor[k]);
        x_s = cm_add(x, num_shifts);
        auto join = [&](const int a, const int b) {
            d_factor[a] = cm_fmadd(d_factor[a], factor[b], cm_mul(factor[a], d_factor[b]));
            factor[a] = cm_mul(factor[a], factor[b]);
        };
        join(0, 1); join(2, 3); join(4, 5); join(6, 7); join(8, 9);
        join(0, 2); join(4, 6);
        join(0, 4);
        join(0, 8);
        p = factor[0];
        dp = d_factor[0];
    }

    const vec_t log_x_s = cm_log_wo_checks(x_s);
    const vec_t log_p = cm_log_wo_checks(p);
    const vec_t z = cm_div(one, x_s);
    const vec_t w = cm_mul(z, z);

    vec_t lgs = cm_set1(1.0 / 156.0);
    lgs = cm_fmadd(lgs, w, cm_set1(-691.0 / 360360.0));
    lgs = cm_fmadd(lgs, w, cm_set1(1.0 / 1188.0));
    lgs = cm_fmadd(lgs, w, cm_set1(-1.0 / 1680.0));
    lgs = cm_fmadd(lgs, w, cm_set1(1.0 / 1260.0));
    lgs = cm_fmadd(lgs, w, cm_set1(-1.0 / 360.0));
    lgs = cm_fmadd(lgs, w, cm_set1(1.0 / 12.0));
    lgs = cm_mul(lgs, z);
    lg_out = cm_sub(
        cm_add(cm_sub(cm_mul(cm_sub(x_s, cm_set1(0.5)), log_x_s), x_s),
               cm_add(cm_set1(CUTPOINT_HALF_LOG_TWO_PI), lgs)),
        log_p);

    vec_t dgs = cm_set1(1.0 / 12.0);
    dgs = cm_fmadd(dgs, w, cm_set1(-691.0 / 32760.0));
    dgs = cm_fmadd(dgs, w, cm_set1(1.0 / 132.0));
    dgs = cm_fmadd(dgs, w, cm_set1(-1.0 / 240.0));
    dgs = cm_fmadd(dgs, w, cm_set1(1.0 / 252.0));
    dgs = cm_fmadd(dgs, w, cm_set1(-1.0 / 120.0));
    dgs = cm_fmadd(dgs, w, cm_set1(1.0 / 12.0));
    dgs = cm_mul(dgs, w);
    dg_out = cm_sub(cm_sub(log_x_s, cm_fmadd(cm_set1(0.5), z, dgs)),
                    cm_div(dp, p));

    vec_t tgs = cm_set1(7.0 / 6.0);
    tgs = cm_fmadd(tgs, w, cm_set1(-691.0 / 2730.0));
    tgs = cm_fmadd(tgs, w, cm_set1(5.0 / 66.0));
    tgs = cm_fmadd(tgs, w, cm_set1(-1.0 / 30.0));
    tgs = cm_fmadd(tgs, w, cm_set1(1.0 / 42.0));
    tgs = cm_fmadd(tgs, w, cm_set1(-1.0 / 30.0));
    tgs = cm_fmadd(tgs, w, cm_set1(1.0 / 6.0));
    const vec_t tg_asym = cm_fmadd(cm_mul(tgs, z), w,
                                   cm_fmadd(cm_set1(0.5), w, z));
    tg_out = cm_add(tg_asym, tg_shift);
}

inline void cutpoint_lgamma_digamma_trigamma_simd(const double* x,
                                                  const int n,
                                                  double* lg_out,
                                                  double* dg_out,
                                                  double* tg_out) {
    int first = 0;
    const vec_t smallest = cm_set1(1e-300);
    const vec_t largest = cm_set1(1e300);
    for (; first + kLanes <= n; first += kLanes) {
        const vec_t xb = cm_loadu(x + first);
        const bool all_ok = cm_mask_all(cm_mask_and(cm_cmp_lt(smallest, xb), cm_cmp_lt(xb, largest)));
        if (!all_ok) {
            for (int i = first; i < first + kLanes; ++i) {
                cutpoint_lgamma_digamma_trigamma(x[i], lg_out + i, dg_out + i, tg_out + i);
            }
            continue;
        }
        vec_t lgb, dgb, tgb;
        lgamma_digamma_trigamma_block(xb, lgb, dgb, tgb);
        cm_storeu(lg_out + first, lgb);
        cm_storeu(dg_out + first, dgb);
        cm_storeu(tg_out + first, tgb);
    }
    const int rem = n - first;
    if (rem >= (kLanes / 2)) {
        alignas(64) double xp[kLanes], lgp[kLanes], dgp[kLanes], tgp[kLanes];
        bool all_ok = true;
        for (int l = 0; l < kLanes; ++l) {
            xp[l] = (l < rem) ? x[first + l] : 1.0;
            all_ok &= (xp[l] > 1e-300 && xp[l] < 1e300);
        }
        if (all_ok) {
            vec_t lgb, dgb, tgb;
            lgamma_digamma_trigamma_block(cm_load(xp), lgb, dgb, tgb);
            cm_store(lgp, lgb);
            cm_store(dgp, dgb);
            cm_store(tgp, tgb);
            for (int l = 0; l < rem; ++l) {
                lg_out[first + l] = lgp[l];
                dg_out[first + l] = dgp[l];
                tg_out[first + l] = tgp[l];
            }
            return;
        }
    }
    for (int i = first; i < n; ++i) {
        cutpoint_lgamma_digamma_trigamma(x[i], lg_out + i, dg_out + i, tg_out + i);
    }
}

template <bool with_partials>
inline void medians_and_log_J_block(const double* a_in,
                                    const double* b_in,
                                    const double* log_a_in,
                                    const double* lbeta_in,
                                    const double* dig_ap1_in,
                                    const double* trig_ap1_in,
                                    const double* dig_b_in,
                                    const double* trig_b_in,
                                    const double* log_kappa_in,
                                    const double* dig_total_in,
                                    const unsigned char* reflect_in,
                                    const int count,
                                    const int use_probit_link,
                                    double* popmed_out,
                                    double* logJ_out,
                                    double* d_cut_da_out,
                                    double* d_cut_db_out,
                                    double* d_logJ_da_out,
                                    double* d_logJ_db_out) {
    alignas(64) double ap[kLanes], bp[kLanes], lap[kLanes], lbp[kLanes], u_init[kLanes];
    const mask_t valid_mask = cm_mask_prefix(count);
    for (int l = 0; l < kLanes; ++l) {
        if (l < count) {
            ap[l] = a_in[l];
            bp[l] = b_in[l];
            lap[l] = log_a_in[l];
            lbp[l] = lbeta_in[l];
            const double u0 = (ap[l] >= 1.0)
                ? std::log((ap[l] - 1.0 / 3.0) / (ap[l] + bp[l] - 2.0 / 3.0))
                : ((CUTPOINT_LOG_HALF + lap[l] + lbp[l]) / ap[l]);
            u_init[l] = (u0 < CUTPOINT_LOG_HALF) ? u0 : CUTPOINT_LOG_HALF;
        } else {
            ap[l] = 1.0; bp[l] = 1.0; lap[l] = 0.0; lbp[l] = 0.0; u_init[l] = CUTPOINT_LOG_HALF;
        }
    }

    const vec_t a_v = cm_load(ap);
    const vec_t b_v = cm_load(bp);
    const vec_t la_v = cm_load(lap);
    const vec_t lb_v = cm_load(lbp);
    const vec_t s_v = cm_add(a_v, b_v);
    const vec_t r_v = cm_add(a_v, cm_set1(1.0));
    const vec_t bm1_v = cm_sub(b_v, cm_set1(1.0));
    const vec_t one = cm_set1(1.0);
    const vec_t zero = cm_setzero();
    const vec_t log_half = cm_set1(CUTPOINT_LOG_HALF);

    vec_t u_v = cm_load(u_init);
    vec_t lower_u = cm_set1(-std::numeric_limits<double>::infinity());
    vec_t upper_u = log_half;
    mask_t outer_mask = valid_mask;

    // Offload only very skewed / ultra-large kappa lanes ((a+b)*x >= 0.85*(a+1)) to CF scalar!
    {
        const vec_t x_chk = cm_exp(u_v);
        const mask_t slow_lanes = cm_mask_and(outer_mask, cm_cmp_le(cm_mul(cm_set1(0.85), r_v), cm_mul(s_v, x_chk)));
        if (cm_mask_any(slow_lanes)) {
            alignas(64) double u_tmp[kLanes];
            cm_store(u_tmp, u_v);
            for (int l = 0; l < count; ++l) {
                if (cm_mask_test(slow_lanes, l)) {
                    u_tmp[l] = cutpoint_beta_median_log_halley_scalar(ap[l], bp[l], lap[l], lbp[l]);
                }
            }
            u_v = cm_load(u_tmp);
            outer_mask = cm_mask_andnot(outer_mask, slow_lanes);
        }
    }

    for (int iter = 0; iter < 64 && cm_mask_any(outer_mask); ++iter) {
        const vec_t x_v = cm_exp(u_v);
        const vec_t omx_v = cm_mask_blend(valid_mask, one, cm_sub(one, x_v));
        const vec_t log_1m_x_v = cm_log1m_exp(u_v);

        vec_t F_v = zero;
        vec_t term_v = one;
        mask_t inner_mask = outer_mask;
        for (int n = 0; n < 64 && cm_mask_any(inner_mask); ++n) {
            F_v = cm_mask_add(F_v, inner_mask, F_v, term_v);
            const vec_t nv = cm_set1(static_cast<double>(n));
            const vec_t sn = cm_add(s_v, nv);
            const vec_t rn = cm_add(r_v, nv);
            const vec_t lhs = cm_mul(sn, x_v);
            const mask_t r_ok = cm_cmp_lt(lhs, rn);
            const vec_t bound = cm_mul(cm_set1(1e-17), cm_mul(F_v, cm_sub(rn, lhs)));
            const mask_t t_ok = cm_cmp_lt(cm_mul(term_v, rn), bound);
            inner_mask = cm_mask_andnot(inner_mask, cm_mask_and(r_ok, t_ok));
            if (!cm_mask_any(inner_mask)) break;
            term_v = cm_mask_mul(zero, inner_mask, term_v, cm_div(lhs, rn));
        }

        if (cm_mask_any(inner_mask)) {
            alignas(64) double u_tmp[kLanes];
            cm_store(u_tmp, u_v);
            for (int l = 0; l < count; ++l) {
                if (cm_mask_test(inner_mask, l)) {
                    u_tmp[l] = cutpoint_beta_median_log_halley_scalar(ap[l], bp[l], lap[l], lbp[l]);
                }
            }
            u_v = cm_load(u_tmp);
            outer_mask = cm_mask_andnot(outer_mask, inner_mask);
            if (!cm_mask_any(outer_mask)) break;
            continue;
        }

        F_v = cm_mask_blend(outer_mask, one, F_v);
        const vec_t log_F = cm_log(F_v);
        const vec_t g_v = cm_sub(
            cm_add(cm_sub(cm_sub(cm_fmadd(a_v, u_v, cm_mul(b_v, log_1m_x_v)),
                                 la_v),
                          lb_v),
                   log_F),
            log_half);

        const mask_t g_zero = cm_cmp_eq(g_v, zero);
        const mask_t g_pos = cm_cmp_gt(g_v, zero);
        const mask_t g_neg = cm_mask_and(outer_mask, cm_cmp_lt(g_v, zero));
        upper_u = cm_mask_mov(upper_u, cm_mask_and(outer_mask, g_pos), u_v);
        lower_u = cm_mask_mov(lower_u, g_neg, u_v);

        const vec_t phi_v = cm_div(a_v, cm_mul(omx_v, F_v));
        const vec_t g_over_gp = cm_div(g_v, phi_v);
        const vec_t dphi_r = cm_sub(cm_sub(a_v, cm_mul(bm1_v, cm_div(x_v, omx_v))),
                                    phi_v);
        const vec_t hcorr = cm_fnmadd(cm_set1(0.5), cm_mul(g_over_gp, dphi_r), one);
        const mask_t use_halley = cm_cmp_gt(hcorr, cm_set1(0.25));
        const vec_t step_h = cm_sub(u_v, cm_div(g_over_gp, hcorr));
        const vec_t step_n = cm_sub(u_v, g_over_gp);
        vec_t next_u = cm_mask_blend(use_halley, step_n, step_h);

        const vec_t tol = cm_mul(cm_set1(4e-16), cm_add(one, cm_abs(u_v)));
        const mask_t conv_pre = cm_cmp_le(cm_abs(cm_sub(next_u, u_v)), tol);
        const mask_t br_col = cm_mask_and(
            cm_cmp_le(cm_sub(upper_u, lower_u), cm_mul(cm_set1(8e-15), cm_add(one, cm_abs(u_v)))),
            cm_cmp_le(cm_abs(g_v), cm_set1(1e-10)));
        const mask_t done_now = cm_mask_and(outer_mask, cm_mask_or(g_zero, cm_mask_or(conv_pre, br_col)));
        u_v = cm_mask_mov(u_v, cm_mask_andnot(done_now, g_zero), next_u);
        outer_mask = cm_mask_andnot(outer_mask, done_now);
        if (!cm_mask_any(outer_mask)) break;

        const mask_t in_hi = cm_cmp_lt(next_u, upper_u);
        const mask_t in_lo = cm_cmp_gt(next_u, lower_u);
        const mask_t out_bracket = cm_mask_andnot(outer_mask, cm_mask_and(in_hi, in_lo));
        if (cm_mask_any(out_bracket)) {
            const vec_t inf_v = cm_set1(-std::numeric_limits<double>::infinity());
            const mask_t lo_inf = cm_cmp_eq(lower_u, inf_v);
            const vec_t fallback_inf = cm_sub(
                cm_sub(upper_u, cm_mul(cm_set1(2.0), cm_sub(upper_u, u_v))),
                one);
            const vec_t fallback_mid = cm_mul(cm_set1(0.5), cm_add(lower_u, upper_u));
            const vec_t fb = cm_mask_blend(lo_inf, fallback_mid, fallback_inf);
            next_u = cm_mask_mov(next_u, out_bracket, fb);
        }

        const mask_t conv_post = cm_cmp_le(cm_abs(cm_sub(next_u, u_v)), tol);
        u_v = cm_mask_mov(u_v, outer_mask, next_u);
        outer_mask = cm_mask_andnot(outer_mask, conv_post);
    }

    if (cm_mask_any(outer_mask)) {
        alignas(64) double u_tmp[kLanes];
        cm_store(u_tmp, u_v);
        for (int l = 0; l < count; ++l) {
            if (cm_mask_test(outer_mask, l)) {
                u_tmp[l] = cutpoint_beta_median_log_halley_scalar(ap[l], bp[l], lap[l], lbp[l]);
            }
        }
        u_v = cm_load(u_tmp);
    }

    const vec_t x_v = cm_exp(u_v);
    const vec_t log_1m_x_v = cm_log1m_exp(u_v);
    vec_t F_v = zero, Fs_v = zero, Fr_v = zero, Frr_v = zero, Fsr_v = zero, xFx_v = zero, xFrx_v = zero;
    vec_t sum_inv_s = zero, sum_inv_r = zero, sum_inv_r2 = zero;
    vec_t term_v = one;
    const mask_t slow_p = cm_mask_and(valid_mask, cm_cmp_le(cm_mul(cm_set1(0.85), r_v), cm_mul(s_v, x_v)));
    mask_t smask = cm_mask_andnot(valid_mask, slow_p);
    for (int n = 0; n < 64 && cm_mask_any(smask); ++n) {
        F_v = cm_mask_add(F_v, smask, F_v, term_v);
        Fr_v = cm_mask_sub(Fr_v, smask, Fr_v, cm_mul(term_v, sum_inv_r));
        if constexpr (with_partials) {
            Fs_v = cm_mask_add(Fs_v, smask, Fs_v, cm_mul(term_v, sum_inv_s));
            const vec_t r_w = cm_fmadd(sum_inv_r, sum_inv_r, sum_inv_r2);
            Frr_v = cm_mask_add(Frr_v, smask, Frr_v, cm_mul(term_v, r_w));
            Fsr_v = cm_mask_sub(Fsr_v, smask, Fsr_v, cm_mul(term_v, cm_mul(sum_inv_s, sum_inv_r)));
            const vec_t nv = cm_set1(static_cast<double>(n));
            xFx_v = cm_mask_add(xFx_v, smask, xFx_v, cm_mul(term_v, nv));
            xFrx_v = cm_mask_sub(xFrx_v, smask, xFrx_v, cm_mul(cm_mul(term_v, sum_inv_r), nv));
        }
        const vec_t nv = cm_set1(static_cast<double>(n));
        const vec_t sn = cm_add(s_v, nv);
        const vec_t rn = cm_add(r_v, nv);
        const vec_t w_chk = with_partials
            ? cm_add(cm_add(one, cm_mul(nv, cm_add(one, sum_inv_r))),
                     cm_add(cm_mul(sum_inv_s, sum_inv_r), cm_fmadd(sum_inv_r, sum_inv_r, sum_inv_r2)))
            : cm_add(one, sum_inv_r);
        const mask_t r_ok = cm_cmp_lt(cm_mul(sn, x_v), rn);
        const mask_t t_ok = cm_cmp_lt(cm_mul(term_v, w_chk), cm_mul(cm_set1(1e-17), F_v));
        smask = cm_mask_andnot(smask, cm_mask_and(r_ok, t_ok));
        if (!cm_mask_any(smask)) break;
        const vec_t inv_rn = cm_div(one, rn);
        if constexpr (with_partials) {
            sum_inv_s = cm_mask_add(sum_inv_s, smask, sum_inv_s, cm_div(one, sn));
            sum_inv_r2 = cm_mask_add(sum_inv_r2, smask, sum_inv_r2, cm_mul(inv_rn, inv_rn));
        }
        sum_inv_r = cm_mask_add(sum_inv_r, smask, sum_inv_r, inv_rn);
        term_v = cm_mask_mul(zero, smask, term_v, cm_mul(cm_mul(sn, inv_rn), x_v));
    }
    const mask_t need_cf = cm_mask_or(slow_p, smask);

    vec_t cp_v, lld_v;
    if (use_probit_link == 1) {
        cp_v = cm_inv_Phi_from_log_p(u_v, log_1m_x_v);
        lld_v = cm_add(cm_add(u_v, cm_set1(CUTPOINT_HALF_LOG_TWO_PI)),
                       cm_mul(cm_set1(0.5), cm_mul(cp_v, cp_v)));
        const mask_t deep_tail = cm_mask_and(
            valid_mask, cm_cmp_lt(cp_v, cm_set1(-nico::math::reference::NORMAL_TAIL_THRESH)));
        if (cm_mask_any(deep_tail)) {
            alignas(64) double ctmp[kLanes], ltmp[kLanes];
            cm_store(ctmp, cp_v);
            cm_store(ltmp, lld_v);
            for (int l = 0; l < count; ++l) {
                if (cm_mask_test(deep_tail, l)) {
                    double frac = -ctmp[l];
                    for (int kk = nico::math::reference::MILLS_CF_DEPTH; kk >= 1; --kk) {
                        frac = -ctmp[l] + static_cast<double>(kk) / frac;
                    }
                    ltmp[l] = -std::log(frac);
                }
            }
            lld_v = cm_load(ltmp);
        }
    } else {
        cp_v = cm_sub(u_v, log_1m_x_v);
        lld_v = cm_sub(zero, log_1m_x_v);
    }

    alignas(64) double u_a[kLanes], x_a[kLanes], l1mx_a[kLanes], F_a[kLanes], Fr_a[kLanes], cp_a[kLanes], lld_a[kLanes];
    alignas(64) double Fs_a[kLanes], Frr_a[kLanes], Fsr_a[kLanes], xFx_a[kLanes], xFrx_a[kLanes];
    cm_store(u_a, u_v);
    cm_store(x_a, x_v);
    cm_store(l1mx_a, log_1m_x_v);
    cm_store(F_a, F_v);
    cm_store(Fr_a, Fr_v);
    cm_store(cp_a, cp_v);
    cm_store(lld_a, lld_v);
    if constexpr (with_partials) {
        cm_store(Fs_a, Fs_v);
        cm_store(Frr_a, Frr_v);
        cm_store(Fsr_a, Fsr_v);
        cm_store(xFx_a, xFx_v);
        cm_store(xFrx_a, xFrx_v);
    }

    if (cm_mask_any(need_cf)) {
        for (int l = 0; l < count; ++l) {
            if (cm_mask_test(need_cf, l)) {
                double s7[7];
                beta_2F1_hybrid_partials<with_partials>(ap[l] + bp[l], ap[l] + 1.0, x_a[l], s7);
                F_a[l] = s7[0]; Fr_a[l] = s7[2];
                if constexpr (with_partials) {
                    Fs_a[l] = s7[1]; Frr_a[l] = s7[3]; Fsr_a[l] = s7[4];
                    xFx_a[l] = s7[5]; xFrx_a[l] = s7[6];
                }
            }
        }
        F_v = cm_load(F_a);
        Fr_v = cm_load(Fr_a);
    }

    F_v = cm_mask_blend(valid_mask, one, F_v);
    Fr_v = cm_mask_blend(valid_mask, zero, Fr_v);

    alignas(64) double lk_a[kLanes], dt_a[kLanes], dap_a[kLanes], db_a[kLanes];
    for (int l = 0; l < kLanes; ++l) {
        if (l < count) {
            lk_a[l] = log_kappa_in[l];
            dt_a[l] = dig_total_in[l];
            dap_a[l] = dig_ap1_in[l];
            db_a[l] = dig_b_in[l];
        } else {
            lk_a[l] = 0.0; dt_a[l] = 0.0; dap_a[l] = 1.0; db_a[l] = 0.0;
        }
    }
    const vec_t dap_v = cm_load(dap_a);
    const vec_t db_v = cm_load(db_a);
    const vec_t sc_raw = cm_sub(cm_add(cm_sub(log_1m_x_v, u_v), cm_sub(dap_v, db_v)), cm_div(Fr_v, F_v));
    const vec_t sc_diff_v = cm_mask_blend(valid_mask, one, sc_raw);
    const vec_t log_sc_v = cm_log(sc_diff_v);
    const vec_t base_exp_v = cm_sub(cm_add(cm_set1(CUTPOINT_LOG_HALF), lb_v),
                                    cm_fmadd(a_v, u_v, cm_mul(bm1_v, log_1m_x_v)));
    const vec_t logJ_v = cm_add(cm_add(cm_load(lk_a), log_sc_v), cm_add(base_exp_v, lld_v));

    alignas(64) double sd_a[kLanes], lj_a[kLanes], mult_a[kLanes], linkd_a[kLanes];
    cm_store(sd_a, sc_diff_v);
    cm_store(lj_a, logJ_v);
    if constexpr (with_partials) {
        const vec_t mult_v = cm_exp(cm_mask_blend(valid_mask, zero, base_exp_v));
        const vec_t linkd_v = (use_probit_link == 1)
            ? cm_exp(cm_mask_blend(valid_mask, zero, lld_v))
            : cm_div(one, cm_mask_blend(valid_mask, one, cm_sub(one, x_v)));
        cm_store(mult_a, mult_v);
        cm_store(linkd_a, linkd_v);
    }

    for (int l = 0; l < count; ++l) {
        const double cutpoint = cp_a[l];
        popmed_out[l] = reflect_in[l] ? -cutpoint : cutpoint;
        logJ_out[l] = lj_a[l];
        if constexpr (with_partials) {
            const double u = u_a[l];
            const double x = x_a[l];
            const double log_1m_x = l1mx_a[l];
            const double F = F_a[l];
            const double score_diff = sd_a[l];
            const double inv_F2 = 1.0 / (F * F);
            const double x_over_1mx = x / (1.0 - x);
            const double score_a = u - dap_a[l] + dt_a[l] + (Fs_a[l] + Fr_a[l]) / F;
            const double mult = mult_a[l];
            const double du_da = -score_a * mult;
            const double du_db = -(score_a + score_diff) * mult;
            const double link_d = linkd_a[l];
            const double d_lld_du = (use_probit_link == 1) ? (1.0 + cutpoint * link_d) : x_over_1mx;
            const double d_diff_da = trig_ap1_in[l] - ((Fsr_a[l] + Frr_a[l]) * F - Fr_a[l] * (Fs_a[l] + Fr_a[l])) * inv_F2;
            const double d_diff_db = -trig_b_in[l] - (Fsr_a[l] * F - Fr_a[l] * Fs_a[l]) * inv_F2;
            const double d_diff_du = -x_over_1mx - 1.0 - (xFrx_a[l] * F - Fr_a[l] * xFx_a[l]) * inv_F2;
            const double dH_da = d_diff_da / score_diff + dap_a[l] - 1.0 / ap[l] - dt_a[l] - u;
            const double dH_db = d_diff_db / score_diff + db_a[l] - dt_a[l] - log_1m_x;
            const double dH_du = d_diff_du / score_diff - ap[l] + (bp[l] - 1.0) * x_over_1mx + d_lld_du;
            d_cut_da_out[l] = link_d * du_da;
            d_cut_db_out[l] = link_d * du_db;
            d_logJ_da_out[l] = dH_da + dH_du * du_da;
            d_logJ_db_out[l] = dH_db + dH_du * du_db;
        }
    }
}

#endif  // BMVP_LANES == 8 || BMVP_LANES == 4

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
        cutpoint = nico::math::nico_std_normal_log_qf_value(u, nullptr);
        if (cutpoint < -nico::math::reference::NORMAL_TAIL_THRESH) {
            double frac = -cutpoint;
            for (int kk = nico::math::reference::MILLS_CF_DEPTH; kk >= 1; --kk) {
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

#if !defined(NICO_SIMD_SCALAR_KERNELS) && defined(BMVP_LANES) && ((BMVP_LANES == 8) || (BMVP_LANES == 4))
    cutpoint_lgamma_digamma_trigamma_simd(gamma_in, 2 * M + n_cls, gamma_lg, gamma_dg, gamma_tg);
#else
    for (int i = 0; i < 2 * M + n_cls; ++i) {
        cutpoint_lgamma_digamma_trigamma(gamma_in[i], gamma_lg + i, gamma_dg + i, gamma_tg + i);
    }
#endif

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
#if !defined(NICO_SIMD_SCALAR_KERNELS) && defined(BMVP_LANES) && ((BMVP_LANES == 8) || (BMVP_LANES == 4))
    for (; j0 + kLanes <= M; j0 += kLanes) {
        medians_and_log_J_block<with_partials>(
            a_arr + j0, b_arr + j0, la_arr + j0, lbeta_arr + j0,
            gamma_dg + j0, gamma_tg + j0, gamma_dg + M + j0, gamma_tg + M + j0,
            lk_arr + j0, dgtot_arr + j0, reflect_out + j0, kLanes, use_probit_link,
            popmed_out + j0, logJ_out + j0,
            with_partials ? (d_cut_da_out + j0) : nullptr,
            with_partials ? (d_cut_db_out + j0) : nullptr,
            with_partials ? (d_logJ_da_out + j0) : nullptr,
            with_partials ? (d_logJ_db_out + j0) : nullptr);
    }
    if (M - j0 >= (kLanes / 2)) {
        const int rem = M - j0;
        medians_and_log_J_block<with_partials>(
            a_arr + j0, b_arr + j0, la_arr + j0, lbeta_arr + j0,
            gamma_dg + j0, gamma_tg + j0, gamma_dg + M + j0, gamma_tg + M + j0,
            lk_arr + j0, dgtot_arr + j0, reflect_out + j0, rem, use_probit_link,
            popmed_out + j0, logJ_out + j0,
            with_partials ? (d_cut_da_out + j0) : nullptr,
            with_partials ? (d_cut_db_out + j0) : nullptr,
            with_partials ? (d_logJ_da_out + j0) : nullptr,
            with_partials ? (d_logJ_db_out + j0) : nullptr);
        j0 = M;
    }
#endif
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
}  // namespace nico::math
