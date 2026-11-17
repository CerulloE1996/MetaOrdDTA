#pragma once
// ---- Stan math only: a copy of NicoSIMD/stan_external_induced_dirichlet_with_Se_Sp_prior.hpp (md5 5ae32390c4f540d7f80966762f57533c) with its scalar code and
//      its hand-written gradients, the SIMD branches left out, and log, exp, Phi, ... of doubles from
//      std:: and stan::math:: (the buffer helpers k_exp, k_log, ... are in
//      stan_external_functions_on_Stan_math_only.hpp). Where the notes below speak of SIMD or NicoSIMD
//      kernels, the same values come from Stan math here: no -m flags, any CPU.
//
// -| --------- NicoSIMD induced Dirichlet of ordered cutpoints with a Se / Sp prior at chosen cutoffs, and the
// -| --------- Dirichlet of a simplex with the same Se / Sp prior, Stan external functions -------------------
//
//     nico_induced_dirichlet_with_Se_Sp_prior_lpdf(C | alpha, link, location)
//       = nico_induced_dirichlet_prior_of_ordered_cutpoints_lpdf(C | alpha, link, location) (no Se / Sp prior,
//         e.g. the between-study model of random cutpoints),
//     nico_induced_dirichlet_with_Se_Sp_prior_lpdf(C | alpha, link, location, alpha_of_the_segments,
//                                                  chosen_cutoffs, Se_Sp_prior_type, prior_Se_Sp_dirichlet_A,
//                                                  prior_Se_Sp_beta, prior_Se_Sp_link_mean_SD)
//       = for each set c of K cutpoints (a vector C, or each row of a matrix C), at x = c - location: the induced
//         Dirichlet of c, - dirichlet_lpdf(q | alpha_of_the_segments) + H(x[chosen_cutoffs]),
//     nico_dirichlet_with_Se_Sp_prior_lpdf(phi | log_J, C_baseline, alpha, alpha_of_the_segments, chosen_cutoffs,
//                                          Se_Sp_prior_type, prior_Se_Sp_dirichlet_A, prior_Se_Sp_beta,
//                                          prior_Se_Sp_link_mean_SD, link)
//       = dirichlet_lpdf(phi | alpha) - dirichlet_lpdf(q_phi | alpha_of_the_segments) + sum(log_J)
//         + sum(log f(C_baseline)) + H(C_baseline),
//     nico_Se_Sp_prior_at_chosen_cutoffs_lpdf(x | Se_Sp_prior_type, prior_Se_Sp_dirichlet_A, prior_Se_Sp_beta,
//                                                prior_Se_Sp_link_mean_SD, link)
//       = sum(log f(x)) + H(x) of type 1 or 2 when Se_Sp_prior_type is 1 or 2, of type 4 otherwise (a Dirichlet on
//         the segment probabilities, a Beta on F at each value, or a normal on the link; the m values need not be
//         ordered except for type 1),
//
// with m the number of chosen cutoffs, F = Phi (link = 1) or inv_logit, f its density, q the m + 1 segment
// probabilities of F at the chosen cutoffs (F_1, F_2 - F_1, ..., 1 - F_m), q_phi the sums of phi over the same
// segments, log_J and C_baseline given at the chosen cutoffs (m values each), and H the Se / Sp prior on F at the
// chosen cutoffs:
//
//     Se_Sp_prior_type 1: dirichlet_lpdf(q | prior_Se_Sp_dirichlet_A),
//     Se_Sp_prior_type 2: sum_j beta_lpdf(F_j | prior_Se_Sp_beta[j, 1], prior_Se_Sp_beta[j, 2]),
//     Se_Sp_prior_type 4: sum_j normal_lpdf(x_j | prior_Se_Sp_link_mean_SD[j, 1], prior_Se_Sp_link_mean_SD[j, 2])
//                         - sum_j log f(x_j),
//     any other:          sum_j beta_lpdf(q_j / (1 - F_(j-1)) | prior_Se_Sp_beta[j, ]) - sum_j log1m(F_(j-1)),
//
// F_0 = 0. The Se / Sp prior of ordinal models of test accuracy with a prior on Sp or 1 - Se at chosen cutoffs:
// on fixed cutpoints (the first form) or on the population cutpoints of random cutpoints (the second). Each is ONE
// autodiff node, the expensive functions with BayesMVP's SIMD kernels (log Phi, log_inv_logit, log1m_exp, exp)
// on the log scale, as nico_induced_dirichlet_prior_of_ordered_cutpoints_lpdf, whose forward and reverse passes
// it uses. Every argument may be data or parameters except link, location, chosen_cutoffs and Se_Sp_prior_type;
// the constants are kept as Stan's _lpdf, and left out with propto only when their arguments are data.
//
// Declared in the Stan functions block (no body; stanc --allow-undefined) as any of
//
//     real nico_induced_dirichlet_with_Se_Sp_prior_lpdf(vector C, vector alpha, data int link,
//                                                       data real location);
//     real nico_induced_dirichlet_with_Se_Sp_prior_lpdf(matrix C, vector alpha, data int link,
//                                                       data real location);
//     real nico_induced_dirichlet_with_Se_Sp_prior_lpdf(vector C, vector alpha, data int link,
//                                                       data real location, vector alpha_of_the_segments,
//                                                       data array[] int chosen_cutoffs,
//                                                       data int Se_Sp_prior_type,
//                                                       vector prior_Se_Sp_dirichlet_A, matrix prior_Se_Sp_beta,
//                                                       matrix prior_Se_Sp_link_mean_SD);
//     real nico_dirichlet_with_Se_Sp_prior_lpdf(vector phi, vector log_J, vector C_baseline, vector alpha,
//                                               vector alpha_of_the_segments, data array[] int chosen_cutoffs,
//                                               data int Se_Sp_prior_type, vector prior_Se_Sp_dirichlet_A,
//                                               matrix prior_Se_Sp_beta, matrix prior_Se_Sp_link_mean_SD,
//                                               data int link);
//     real nico_Se_Sp_prior_at_chosen_cutoffs_lpdf(vector x, data int Se_Sp_prior_type,
//                                                  vector prior_Se_Sp_dirichlet_A, matrix prior_Se_Sp_beta,
//                                                  matrix prior_Se_Sp_link_mean_SD, data int link);
//
// -----------------------------------------------------------------------------------------------------------

#include "stan_external_induced_dirichlet_prior_of_ordered_cutpoints.hpp"

#include <stan/math/rev.hpp>
#include <stan/math/prim.hpp>

#include <Eigen/Dense>
#include <cmath>
#include <ostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace fused_stan_math {
namespace Se_Sp_prior_at_chosen_cutoffs {

static constexpr double half_log_of_two_pi = 0.91893853320467274178;




// ---- the arena buffers of the Se / Sp prior at m chosen cutoffs: the values of the hyperparameters, the
//      log-scale values of one set and the partials of H (d_x of one set; the hyperparameters' summed over sets)
struct buffers_of_the_Se_Sp_prior {

    int m = 0;
    double* ones = nullptr;           // m + 1
    double* segments = nullptr;       // m + 1: alpha_of_the_segments
    double* A = nullptr;              // m + 1
    double* a = nullptr;              // m
    double* b = nullptr;              // m
    double* mu = nullptr;             // m
    double* sd = nullptr;             // m
    double* lF = nullptr;             // m
    double* lS = nullptr;             // m + 1
    double* log_q = nullptr;          // 2 m + 1
    double* log_f = nullptr;          // m
    double* d_segments = nullptr;     // m + 1
    double* d_A = nullptr;            // m + 1
    double* d_a = nullptr;            // m
    double* d_b = nullptr;            // m
    double* d_mu = nullptr;           // m
    double* d_sd = nullptr;           // m

    inline static std::size_t n_values(const int m) {

        return static_cast<std::size_t>(18) * m + 7;

    }

    inline void place(double* start, const int m_of_the_prior) {

        m = m_of_the_prior;
        double* next = start;
        auto take = [&next](const int n) { double* out = next; next += n; return out; };
        ones = take(m + 1);
        segments = take(m + 1);
        A = take(m + 1);
        a = take(m);
        b = take(m);
        mu = take(m);
        sd = take(m);
        lF = take(m);
        lS = take(m + 1);
        log_q = take(2 * m + 1);
        log_f = take(m);
        d_segments = take(m + 1);
        d_A = take(m + 1);
        d_a = take(m);
        d_b = take(m);
        d_mu = take(m);
        d_sd = take(m);
        for (int j = 0; j <= m; ++j) {
            ones[j] = 1.0;
            d_segments[j] = 0.0;
            d_A[j] = 0.0;
        }
        for (int j = 0; j < m; ++j) {
            d_a[j] = 0.0;
            d_b[j] = 0.0;
            d_mu[j] = 0.0;
            d_sd[j] = 0.0;
        }

    }

};




// ---- sum_i (w_i - 1) log q_i + (lgamma(sum w) - sum lgamma(w)) [include_the_constant] of the m + 1 segment
//      probabilities, its partials d / d x_j (times sign, added to d_x) and d / d w_i (times sign, added to d_w)
inline double dirichlet_of_the_segments(const double* w, const double* log_q, const double* log_f, const int m,
                                        const bool include_the_constant, const double sign, double* d_x,
                                        double* d_w) {

    double total = 0.0;
    double sum_of_w = 0.0;
    for (int i = 0; i <= m; ++i) {
        if (w[i] != 1.0) total += (w[i] - 1.0) * log_q[i];
        sum_of_w += w[i];
    }
    if (include_the_constant) {
        double sum_of_lgamma = 0.0;
        for (int i = 0; i <= m; ++i) sum_of_lgamma += std::lgamma(w[i]);
        total += std::lgamma(sum_of_w) - sum_of_lgamma;
    }
    if (d_x != nullptr) {
        for (int j = 0; j < m; ++j) {
            double d_F = 0.0;
            if (w[j] != 1.0) d_F += (w[j] - 1.0) * std::exp(log_f[j] - log_q[j]);
            if (w[j + 1] != 1.0) d_F -= (w[j + 1] - 1.0) * std::exp(log_f[j] - log_q[j + 1]);
            d_x[j] += sign * d_F;
        }
    }
    if (d_w != nullptr) {
        const double digamma_of_the_sum = stan::math::digamma(sum_of_w);
        for (int i = 0; i <= m; ++i) {
            d_w[i] += sign * (log_q[i] + digamma_of_the_sum - stan::math::digamma(w[i]));
        }
    }
    return sign * total;

}




// ---- H(x) of one set of m link values x (the header), from the buffers' lF, lS, log_q and log_f; its partials d / d x_j added to d_x, those of the hyperparameters to the
//      buffers' d_ arrays when differentiate_the_hyperparameters
inline double Se_Sp_prior(const double* x, const int link, const int Se_Sp_prior_type,
                          buffers_of_the_Se_Sp_prior& buffers, const bool include_the_constants,
                          const bool differentiate_the_hyperparameters, double* d_x) {

    const int m = buffers.m;
    const double* lF = buffers.lF;
    const double* lS = buffers.lS;
    const double* log_q = buffers.log_q;
    const double* log_f = buffers.log_f;
    double total = 0.0;
    if (Se_Sp_prior_type == 1) {
        //// (Dirichlet on the segment probabilities)
        total += dirichlet_of_the_segments(buffers.A, log_q, log_f, m, include_the_constants, 1.0, d_x,
                                           differentiate_the_hyperparameters ? buffers.d_A : nullptr);
    } else if (Se_Sp_prior_type == 2) {
        //// (Beta on each F_j)
        for (int j = 0; j < m; ++j) {
            const double a = buffers.a[j];
            const double b = buffers.b[j];
            total += (a - 1.0) * lF[j] + (b - 1.0) * lS[j];
            if (include_the_constants) total -= stan::math::lbeta(a, b);
            if (d_x != nullptr) {
                d_x[j] += (a - 1.0) * std::exp(log_f[j] - lF[j]) - (b - 1.0) * std::exp(log_f[j] - lS[j]);
            }
            if (differentiate_the_hyperparameters) {
                const double digamma_of_the_sum = stan::math::digamma(a + b);
                buffers.d_a[j] += lF[j] - stan::math::digamma(a) + digamma_of_the_sum;
                buffers.d_b[j] += lS[j] - stan::math::digamma(b) + digamma_of_the_sum;
            }
        }
    } else if (Se_Sp_prior_type == 4) {
        //// (normal on the link values, minus log f: the Jacobian from F to the link)
        for (int j = 0; j < m; ++j) {
            const double mu = buffers.mu[j];
            const double sd = buffers.sd[j];
            const double z = (x[j] - mu) / sd;
            total += -0.5 * z * z - log_f[j];
            if (include_the_constants) total -= std::log(sd) + half_log_of_two_pi;
            const double d_log_f = (link == 1) ? -x[j] : (std::exp(lS[j]) - std::exp(lF[j]));
            if (d_x != nullptr) d_x[j] += -z / sd - d_log_f;
            if (differentiate_the_hyperparameters) {
                buffers.d_mu[j] += z / sd;
                buffers.d_sd[j] += (z * z - 1.0) / sd;
            }
        }
    } else {
        //// (Beta on each (F_j - F_(j-1)) / (1 - F_(j-1)), minus each log(1 - F_(j-1)); F_0 = 0)
        for (int j = 0; j < m; ++j) {
            const double a = buffers.a[j];
            const double b = buffers.b[j];
            const double lS_previous = (j == 0) ? 0.0 : lS[j - 1];
            const double log_u = log_q[j] - lS_previous;
            const double log1m_u = lS[j] - lS_previous;
            total += (a - 1.0) * log_u + (b - 1.0) * log1m_u - lS_previous;
            if (include_the_constants) total -= stan::math::lbeta(a, b);
            if (d_x != nullptr) {
                d_x[j] += (a - 1.0) * std::exp(log_f[j] - log_q[j]) - (b - 1.0) * std::exp(log_f[j] - lS[j]);
                if (j + 1 < m) {
                    const double a_next = buffers.a[j + 1];
                    const double b_next = buffers.b[j + 1];
                    d_x[j] += -(a_next - 1.0) * std::exp(log_f[j] - log_q[j + 1]) +
                              (a_next + b_next - 1.0) * std::exp(log_f[j] - lS[j]);
                }
            }
            if (differentiate_the_hyperparameters) {
                const double digamma_of_the_sum = stan::math::digamma(a + b);
                buffers.d_a[j] += log_u - stan::math::digamma(a) + digamma_of_the_sum;
                buffers.d_b[j] += log1m_u - stan::math::digamma(b) + digamma_of_the_sum;
            }
        }
    }
    return total;

}




// ---- the log-scale values of one set of m increasing link values x into the buffers (lF, lS, log_q, and log_f the
//      log density of F, with its constant); false when x is not increasing or not a number
inline bool values_of_one_set(const double* x, const int link, buffers_of_the_Se_Sp_prior& buffers) {

    double unused = 0.0;
    const bool ok = fused_stan_math::induced_dirichlet_prior::forward_pass<false>(x, buffers.ones, link, buffers.m,
                                                                            buffers.lF, buffers.lS, buffers.log_q,
                                                                            buffers.log_f, unused);
    if (link == 1) {
        for (int j = 0; j < buffers.m; ++j) buffers.log_f[j] -= half_log_of_two_pi;
    }
    return ok;

}




// ---- lF, lS and log_f of m link values in any order (no segment probabilities) into the buffers
inline void values_of_one_set_in_any_order(const double* x, const int link, buffers_of_the_Se_Sp_prior& buffers) {

    const int m = buffers.m;
    copy_n(x, m, buffers.lF);
    for (int j = 0; j < m; ++j) buffers.lS[j] = -x[j];
    if (link == 1) {
        k_log_Phi(buffers.lF, m);
        k_log_Phi(buffers.lS, m);
        for (int j = 0; j < m; ++j) buffers.log_f[j] = -0.5 * x[j] * x[j] - half_log_of_two_pi;
    } else {
        k_log_inv_logit(buffers.lF, m);
        k_log_inv_logit(buffers.lS, m);
        for (int j = 0; j < m; ++j) buffers.log_f[j] = buffers.lF[j] + buffers.lS[j];
    }

}




// ---- the values of the hyperparameters into the buffers, with Stan's checks of their sizes
template <typename T_segments, typename T_A, typename T_beta, typename T_link_mean_SD>
inline void values_of_the_hyperparameters(const char* function_name, const T_segments& alpha_of_the_segments,
                                          const int Se_Sp_prior_type, const T_A& prior_Se_Sp_dirichlet_A,
                                          const T_beta& prior_Se_Sp_beta,
                                          const T_link_mean_SD& prior_Se_Sp_link_mean_SD,
                                          buffers_of_the_Se_Sp_prior& buffers) {

    const int m = buffers.m;
    stan::math::check_size_match(function_name, "Number of segment weights", alpha_of_the_segments.size(),
                                 "Number of chosen cutoffs plus 1", m + 1);
    for (int i = 0; i <= m; ++i) buffers.segments[i] = stan::math::value_of(alpha_of_the_segments.coeff(i));
    if (Se_Sp_prior_type == 1) {
        stan::math::check_size_match(function_name, "Number of Se / Sp Dirichlet weights",
                                     prior_Se_Sp_dirichlet_A.size(), "Number of chosen cutoffs plus 1", m + 1);
        for (int i = 0; i <= m; ++i) buffers.A[i] = stan::math::value_of(prior_Se_Sp_dirichlet_A.coeff(i));
    } else if (Se_Sp_prior_type == 4) {
        stan::math::check_size_match(function_name, "Rows of the link mean and SD", prior_Se_Sp_link_mean_SD.rows(),
                                     "Number of chosen cutoffs", m);
        for (int j = 0; j < m; ++j) {
            buffers.mu[j] = stan::math::value_of(prior_Se_Sp_link_mean_SD.coeff(j, 0));
            buffers.sd[j] = stan::math::value_of(prior_Se_Sp_link_mean_SD.coeff(j, 1));
        }
    } else {
        stan::math::check_size_match(function_name, "Rows of the Beta parameters", prior_Se_Sp_beta.rows(),
                                     "Number of chosen cutoffs", m);
        for (int j = 0; j < m; ++j) {
            buffers.a[j] = stan::math::value_of(prior_Se_Sp_beta.coeff(j, 0));
            buffers.b[j] = stan::math::value_of(prior_Se_Sp_beta.coeff(j, 1));
        }
    }

}




// ---- the hyperparameters' partials (times the adjoint) into their vars, in the reverse pass
template <typename T_segments_vars, typename T_A_vars, typename T_beta_vars, typename T_link_mean_SD_vars>
inline void add_the_partials_of_the_hyperparameters(const int Se_Sp_prior_type,
                                                    const buffers_of_the_Se_Sp_prior& buffers, const double adjoint,
                                                    T_segments_vars& segments_vars, T_A_vars& A_vars,
                                                    T_beta_vars& beta_vars,
                                                    T_link_mean_SD_vars& link_mean_SD_vars) {

    const int m = buffers.m;
    if (segments_vars.size() > 0) {
        for (int i = 0; i <= m; ++i) segments_vars.coeffRef(i).adj() += adjoint * buffers.d_segments[i];
    }
    if ((Se_Sp_prior_type == 1) && (A_vars.size() > 0)) {
        for (int i = 0; i <= m; ++i) A_vars.coeffRef(i).adj() += adjoint * buffers.d_A[i];
    } else if ((Se_Sp_prior_type == 4) && (link_mean_SD_vars.size() > 0)) {
        for (int j = 0; j < m; ++j) {
            link_mean_SD_vars.coeffRef(j, 0).adj() += adjoint * buffers.d_mu[j];
            link_mean_SD_vars.coeffRef(j, 1).adj() += adjoint * buffers.d_sd[j];
        }
    } else if ((Se_Sp_prior_type != 1) && (Se_Sp_prior_type != 4) && (beta_vars.size() > 0)) {
        for (int j = 0; j < m; ++j) {
            beta_vars.coeffRef(j, 0).adj() += adjoint * buffers.d_a[j];
            beta_vars.coeffRef(j, 1).adj() += adjoint * buffers.d_b[j];
        }
    }

}




// ---- the chosen cutoffs (1-based in Stan) checked: increasing, within 1, ..., K
inline void check_the_chosen_cutoffs(const char* function_name, const std::vector<int>& chosen_cutoffs,
                                     const int K) {

    for (std::size_t j = 0; j < chosen_cutoffs.size(); ++j) {
        if ((chosen_cutoffs[j] < 1) || (chosen_cutoffs[j] > K) ||
            ((j > 0) && (chosen_cutoffs[j] <= chosen_cutoffs[j - 1]))) {
            throw std::domain_error(std::string(function_name) + ": chosen_cutoffs must rise within 1, ..., " +
                                    std::to_string(K));
        }
    }

}

}  // namespace Se_Sp_prior_at_chosen_cutoffs
}  // namespace fused_stan_math




// ---- the Stan-facing functions (global scope; stanc calls them with propto first and pstream__ last)
namespace fused_stan_math {

//// (no Se / Sp prior: the induced Dirichlet itself)
template <bool propto, typename T_C, typename T_alpha, stan::require_eigen_t<T_C>* = nullptr,
          stan::require_eigen_col_vector_t<T_alpha>* = nullptr>
inline stan::return_type_t<T_C, T_alpha> nico_induced_dirichlet_with_Se_Sp_prior_lpdf(
        const T_C& C, const T_alpha& alpha, const int link, const double location, std::ostream* pstream__) {

    return nico_induced_dirichlet_prior_of_ordered_cutpoints_lpdf<propto>(C, alpha, link, location, pstream__);

}




//// (with the Se / Sp prior at the chosen cutoffs of every set)
template <bool propto, typename T_C, typename T_alpha, typename T_segments, typename T_A, typename T_beta,
          typename T_link_mean_SD, stan::require_eigen_t<T_C>* = nullptr,
          stan::require_eigen_col_vector_t<T_alpha>* = nullptr>
inline stan::return_type_t<T_C, T_alpha, T_segments, T_A, T_beta, T_link_mean_SD>
nico_induced_dirichlet_with_Se_Sp_prior_lpdf(const T_C& C, const T_alpha& alpha, const int link,
                                             const double location, const T_segments& alpha_of_the_segments,
                                             const std::vector<int>& chosen_cutoffs, const int Se_Sp_prior_type,
                                             const T_A& prior_Se_Sp_dirichlet_A, const T_beta& prior_Se_Sp_beta,
                                             const T_link_mean_SD& prior_Se_Sp_link_mean_SD,
                                             std::ostream* pstream__) {

    namespace prior = fused_stan_math::induced_dirichlet_prior;
    namespace Se_Sp = fused_stan_math::Se_Sp_prior_at_chosen_cutoffs;
    using stan::math::var;
    static constexpr const char* function_name = "nico_induced_dirichlet_with_Se_Sp_prior_lpdf";
    const int m = static_cast<int>(chosen_cutoffs.size());
    if (m == 0) {
        return nico_induced_dirichlet_prior_of_ordered_cutpoints_lpdf<propto>(C, alpha, link, location, pstream__);
    }
    constexpr bool C_is_one_set = stan::is_eigen_vector<T_C>::value;
    constexpr bool C_is_constant = stan::is_constant_all<T_C>::value;
    constexpr bool alpha_is_constant = stan::is_constant_all<T_alpha>::value;
    constexpr bool hyperparameters_are_constant =
            stan::is_constant_all<T_segments, T_A, T_beta, T_link_mean_SD>::value;
    constexpr bool everything_is_constant = C_is_constant && alpha_is_constant && hyperparameters_are_constant;
    constexpr bool include_the_lgamma_terms = !propto || !alpha_is_constant;
    constexpr bool include_the_constant_of_the_normal = !propto;
    constexpr bool include_the_constants_of_the_Se_Sp_prior = !propto || !hyperparameters_are_constant;
    const int R = C_is_one_set ? 1 : static_cast<int>(C.rows());
    const int K = C_is_one_set ? static_cast<int>(C.size()) : static_cast<int>(C.cols());
    stan::math::check_size_match(function_name, "Number of prior weights", alpha.size(),
                                 "Number of cutpoints plus 1", K + 1);
    Se_Sp::check_the_chosen_cutoffs(function_name, chosen_cutoffs, K);
    if (R == 0) return 0.0;
    if constexpr (propto && everything_is_constant) return 0.0;
    const auto& C_ref = stan::math::to_ref(C);
    const auto& alpha_ref = stan::math::to_ref(alpha);
    //// (the buffers, as nico_induced_dirichlet_prior_of_ordered_cutpoints_lpdf's, then those of the Se / Sp prior:
    ////  the chosen link values of every set, the partials of H with respect to them (R m), and one set's work
    ////  space; on the arena when anything is differentiated, in a std::vector otherwise)
    const std::size_t RK = static_cast<std::size_t>(R) * K;
    const std::size_t Rm = static_cast<std::size_t>(R) * m;
    const std::size_t n_buffer_values = 3 * RK + (K + 1) + static_cast<std::size_t>(R) * (K + 1) +
                                        static_cast<std::size_t>(R) * (2 * K + 1) +
                                        2 * static_cast<std::size_t>(K) + 2 * Rm +
                                        Se_Sp::buffers_of_the_Se_Sp_prior::n_values(m);
    std::vector<double> buffer_without_autodiff;
    double* buffers = nullptr;
    if constexpr (everything_is_constant) {
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
    double* x_chosen = ratios + 2 * static_cast<std::size_t>(K);
    double* d_x_chosen = x_chosen + Rm;
    Se_Sp::buffers_of_the_Se_Sp_prior Se_Sp_buffers;
    Se_Sp_buffers.place(d_x_chosen + Rm, m);
    for (int r = 0; r < R; ++r) {
        for (int k = 0; k < K; ++k) {
            const double value_of_C = C_is_one_set ? stan::math::value_of(C_ref.coeff(k))
                                                   : stan::math::value_of(C_ref.coeff(r, k));
            x[static_cast<std::size_t>(r) * K + k] = value_of_C - location;
        }
        for (int j = 0; j < m; ++j) {
            x_chosen[static_cast<std::size_t>(r) * m + j] =
                    x[static_cast<std::size_t>(r) * K + (chosen_cutoffs[j] - 1)];
            d_x_chosen[static_cast<std::size_t>(r) * m + j] = 0.0;
        }
    }
    for (int j = 0; j <= K; ++j) alpha_kept[j] = stan::math::value_of(alpha_ref.coeff(j));
    Se_Sp::values_of_the_hyperparameters(function_name, alpha_of_the_segments, Se_Sp_prior_type,
                                         prior_Se_Sp_dirichlet_A, prior_Se_Sp_beta, prior_Se_Sp_link_mean_SD,
                                         Se_Sp_buffers);
    double sum = 0.0;
    for (int r = 0; r < R; ++r) {
        //// (the induced Dirichlet of the set)
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
        //// (minus the Dirichlet of its segment probabilities, plus the Se / Sp prior, at the chosen cutoffs)
        const double* x_r = x_chosen + static_cast<std::size_t>(r) * m;
        double* d_x_r = d_x_chosen + static_cast<std::size_t>(r) * m;
        if (!Se_Sp::values_of_one_set(x_r, link, Se_Sp_buffers)) {
            throw std::domain_error(std::string(function_name) + ": the cutpoints at the chosen cutoffs are not "
                                    "increasing");
        }
        sum += Se_Sp::dirichlet_of_the_segments(Se_Sp_buffers.segments, Se_Sp_buffers.log_q,
                                                Se_Sp_buffers.log_f, m,
                                                include_the_constants_of_the_Se_Sp_prior, -1.0, d_x_r,
                                                stan::is_constant_all<T_segments>::value
                                                        ? nullptr : Se_Sp_buffers.d_segments);
        sum += Se_Sp::Se_Sp_prior(x_r, link, Se_Sp_prior_type, Se_Sp_buffers,
                                  include_the_constants_of_the_Se_Sp_prior, !hyperparameters_are_constant, d_x_r);
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
    if constexpr (everything_is_constant) {
        return sum;
    } else {
        using C_vars_type = std::conditional_t<C_is_constant, Eigen::Matrix<var, -1, -1>,
                                               stan::plain_type_t<stan::math::promote_scalar_t<var, T_C>>>;
        stan::arena_t<C_vars_type> C_vars;
        stan::arena_t<Eigen::Matrix<var, -1, 1>> alpha_vars;
        stan::arena_t<Eigen::Matrix<var, -1, 1>> segments_vars;
        stan::arena_t<Eigen::Matrix<var, -1, 1>> A_vars;
        stan::arena_t<Eigen::Matrix<var, -1, -1>> beta_vars;
        stan::arena_t<Eigen::Matrix<var, -1, -1>> link_mean_SD_vars;
        if constexpr (!C_is_constant) C_vars = C_ref;
        if constexpr (!alpha_is_constant) alpha_vars = alpha_ref;
        if constexpr (!stan::is_constant_all<T_segments>::value) segments_vars = alpha_of_the_segments;
        if constexpr (!stan::is_constant_all<T_A>::value) A_vars = prior_Se_Sp_dirichlet_A;
        if constexpr (!stan::is_constant_all<T_beta>::value) beta_vars = prior_Se_Sp_beta;
        if constexpr (!stan::is_constant_all<T_link_mean_SD>::value) link_mean_SD_vars = prior_Se_Sp_link_mean_SD;
        int* chosen_kept = stan::math::ChainableStack::instance_->memalloc_.alloc_array<int>(m);
        for (int j = 0; j < m; ++j) chosen_kept[j] = chosen_cutoffs[j] - 1;
        return stan::math::make_callback_var(sum, [C_vars, alpha_vars, segments_vars, A_vars, beta_vars,
                                                   link_mean_SD_vars, x, alpha_kept, link, R, K, m, F, p, log_p,
                                                   J_terms, ratios, d_x_chosen, chosen_kept, Se_Sp_buffers,
                                                   Se_Sp_prior_type](auto& vi) mutable {
            if constexpr (!C_is_constant) {
                for (int r = 0; r < R; ++r) {
                    double* d_x = J_terms + static_cast<std::size_t>(r) * K;
                    prior::reverse_pass(x + static_cast<std::size_t>(r) * K, alpha_kept, link, K,
                                        F + static_cast<std::size_t>(r) * K,
                                        p + static_cast<std::size_t>(r) * (K + 1),
                                        log_p + static_cast<std::size_t>(r) * (2 * K + 1), vi.adj(), d_x,
                                        ratios, ratios + K);
                    for (int j = 0; j < m; ++j) {
                        d_x[chosen_kept[j]] += vi.adj() * d_x_chosen[static_cast<std::size_t>(r) * m + j];
                    }
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
            Se_Sp::add_the_partials_of_the_hyperparameters(Se_Sp_prior_type, Se_Sp_buffers, vi.adj(),
                                                           segments_vars, A_vars, beta_vars, link_mean_SD_vars);
        });
    }

}

}  // namespace fused_stan_math

using fused_stan_math::nico_induced_dirichlet_with_Se_Sp_prior_lpdf;




// ---- the Dirichlet of a simplex phi with the Se / Sp prior at the chosen cutoffs (global scope)
namespace fused_stan_math {

template <bool propto, typename T_phi, typename T_log_J, typename T_x, typename T_alpha, typename T_segments,
          typename T_A, typename T_beta, typename T_link_mean_SD,
          stan::require_eigen_col_vector_t<T_phi>* = nullptr>
inline stan::return_type_t<T_phi, T_log_J, T_x, T_alpha, T_segments, T_A, T_beta, T_link_mean_SD>
nico_dirichlet_with_Se_Sp_prior_lpdf(const T_phi& phi, const T_log_J& log_J, const T_x& C_baseline,
                                     const T_alpha& alpha, const T_segments& alpha_of_the_segments,
                                     const std::vector<int>& chosen_cutoffs, const int Se_Sp_prior_type,
                                     const T_A& prior_Se_Sp_dirichlet_A, const T_beta& prior_Se_Sp_beta,
                                     const T_link_mean_SD& prior_Se_Sp_link_mean_SD, const int link,
                                     std::ostream* /*pstream__*/) {

    namespace Se_Sp = fused_stan_math::Se_Sp_prior_at_chosen_cutoffs;
    using stan::math::var;
    static constexpr const char* function_name = "nico_dirichlet_with_Se_Sp_prior_lpdf";
    const int n_cat = static_cast<int>(phi.size());
    const int m = static_cast<int>(chosen_cutoffs.size());
    constexpr bool phi_is_constant = stan::is_constant_all<T_phi>::value;
    constexpr bool log_J_is_constant = stan::is_constant_all<T_log_J>::value;
    constexpr bool x_is_constant = stan::is_constant_all<T_x>::value;
    constexpr bool alpha_is_constant = stan::is_constant_all<T_alpha>::value;
    constexpr bool segments_are_constant = stan::is_constant_all<T_segments>::value;
    constexpr bool hyperparameters_are_constant =
            stan::is_constant_all<T_segments, T_A, T_beta, T_link_mean_SD>::value;
    constexpr bool everything_is_constant = phi_is_constant && log_J_is_constant && x_is_constant &&
                                            alpha_is_constant && hyperparameters_are_constant;
    constexpr bool include_the_lgamma_terms = !propto || !alpha_is_constant;
    constexpr bool include_the_constants_of_the_Se_Sp_prior = !propto || !hyperparameters_are_constant;
    stan::math::check_size_match(function_name, "Number of prior weights", alpha.size(), "Number of categories",
                                 n_cat);
    stan::math::check_size_match(function_name, "Number of log_J", log_J.size(), "Number of chosen cutoffs", m);
    stan::math::check_size_match(function_name, "Number of C_baseline", C_baseline.size(),
                                 "Number of chosen cutoffs", m);
    Se_Sp::check_the_chosen_cutoffs(function_name, chosen_cutoffs, n_cat - 1);
    if constexpr (propto && everything_is_constant) return 0.0;
    const auto& phi_ref = stan::math::to_ref(phi);
    const auto& log_J_ref = stan::math::to_ref(log_J);
    const auto& x_ref = stan::math::to_ref(C_baseline);
    const auto& alpha_ref = stan::math::to_ref(alpha);
    //// (the buffers: phi, alpha, the log of the m + 1 segment sums of phi, x, d_x (m), then the Se / Sp prior's)
    const std::size_t n_buffer_values = 2 * static_cast<std::size_t>(n_cat) + (m + 1) + 2 * m +
                                        Se_Sp::buffers_of_the_Se_Sp_prior::n_values(m);
    std::vector<double> buffer_without_autodiff;
    double* buffers = nullptr;
    if constexpr (everything_is_constant) {
        buffer_without_autodiff.resize(n_buffer_values);
        buffers = buffer_without_autodiff.data();
    } else {
        buffers = stan::math::ChainableStack::instance_->memalloc_.alloc_array<double>(n_buffer_values);
    }
    double* phi_kept = buffers;
    double* alpha_kept = phi_kept + n_cat;
    double* log_q_phi = alpha_kept + n_cat;
    double* x = log_q_phi + (m + 1);
    double* d_x = x + m;
    Se_Sp::buffers_of_the_Se_Sp_prior Se_Sp_buffers;
    Se_Sp_buffers.place(d_x + m, m);
    for (int n = 0; n < n_cat; ++n) {
        phi_kept[n] = stan::math::value_of(phi_ref.coeff(n));
        alpha_kept[n] = stan::math::value_of(alpha_ref.coeff(n));
    }
    for (int j = 0; j < m; ++j) {
        x[j] = stan::math::value_of(x_ref.coeff(j));
        d_x[j] = 0.0;
    }
    Se_Sp::values_of_the_hyperparameters(function_name, alpha_of_the_segments, Se_Sp_prior_type,
                                         prior_Se_Sp_dirichlet_A, prior_Se_Sp_beta, prior_Se_Sp_link_mean_SD,
                                         Se_Sp_buffers);
    //// (the Dirichlet of phi)
    double sum = 0.0;
    double sum_of_alpha = 0.0;
    for (int n = 0; n < n_cat; ++n) {
        if (alpha_kept[n] != 1.0) sum += (alpha_kept[n] - 1.0) * std::log(phi_kept[n]);
        sum_of_alpha += alpha_kept[n];
    }
    if constexpr (include_the_lgamma_terms) {
        double sum_of_lgamma_of_alpha = 0.0;
        for (int n = 0; n < n_cat; ++n) sum_of_lgamma_of_alpha += std::lgamma(alpha_kept[n]);
        sum += std::lgamma(sum_of_alpha) - sum_of_lgamma_of_alpha;
    }
    //// (minus the Dirichlet of the segment sums of phi: segment s ends at the chosen cutoff s, the last at n_cat)
    {
        int start = 0;
        for (int s = 0; s <= m; ++s) {
            const int end = (s < m) ? chosen_cutoffs[s] : n_cat;
            double q = 0.0;
            for (int n = start; n < end; ++n) q += phi_kept[n];
            log_q_phi[s] = std::log(q);
            start = end;
        }
        sum += Se_Sp::dirichlet_of_the_segments(Se_Sp_buffers.segments, log_q_phi, nullptr, m,
                                                include_the_constants_of_the_Se_Sp_prior, -1.0, nullptr,
                                                segments_are_constant ? nullptr : Se_Sp_buffers.d_segments);
    }
    //// (sum(log_J), then sum(log f(x)) and the Se / Sp prior at x = C_baseline)
    for (int j = 0; j < m; ++j) sum += stan::math::value_of(log_J_ref.coeff(j));
    if (!Se_Sp::values_of_one_set(x, link, Se_Sp_buffers)) {
        throw std::domain_error(std::string(function_name) + ": C_baseline is not increasing");
    }
    for (int j = 0; j < m; ++j) {
        sum += Se_Sp_buffers.log_f[j];
        d_x[j] += (link == 1) ? -x[j] : (std::exp(Se_Sp_buffers.lS[j]) - std::exp(Se_Sp_buffers.lF[j]));
    }
    sum += Se_Sp::Se_Sp_prior(x, link, Se_Sp_prior_type, Se_Sp_buffers, include_the_constants_of_the_Se_Sp_prior,
                              !hyperparameters_are_constant, d_x);
    if constexpr (everything_is_constant) {
        return sum;
    } else {
        stan::arena_t<Eigen::Matrix<var, -1, 1>> phi_vars;
        stan::arena_t<Eigen::Matrix<var, -1, 1>> log_J_vars;
        stan::arena_t<Eigen::Matrix<var, -1, 1>> x_vars;
        stan::arena_t<Eigen::Matrix<var, -1, 1>> alpha_vars;
        stan::arena_t<Eigen::Matrix<var, -1, 1>> segments_vars;
        stan::arena_t<Eigen::Matrix<var, -1, 1>> A_vars;
        stan::arena_t<Eigen::Matrix<var, -1, -1>> beta_vars;
        stan::arena_t<Eigen::Matrix<var, -1, -1>> link_mean_SD_vars;
        if constexpr (!phi_is_constant) phi_vars = phi_ref;
        if constexpr (!log_J_is_constant) log_J_vars = log_J_ref;
        if constexpr (!x_is_constant) x_vars = x_ref;
        if constexpr (!alpha_is_constant) alpha_vars = alpha_ref;
        if constexpr (!segments_are_constant) segments_vars = alpha_of_the_segments;
        if constexpr (!stan::is_constant_all<T_A>::value) A_vars = prior_Se_Sp_dirichlet_A;
        if constexpr (!stan::is_constant_all<T_beta>::value) beta_vars = prior_Se_Sp_beta;
        if constexpr (!stan::is_constant_all<T_link_mean_SD>::value) link_mean_SD_vars = prior_Se_Sp_link_mean_SD;
        int* segment_ends = stan::math::ChainableStack::instance_->memalloc_.alloc_array<int>(m + 1);
        for (int s = 0; s < m; ++s) segment_ends[s] = chosen_cutoffs[s];
        segment_ends[m] = n_cat;
        return stan::math::make_callback_var(sum, [phi_vars, log_J_vars, x_vars, alpha_vars, segments_vars, A_vars,
                                                   beta_vars, link_mean_SD_vars, phi_kept, alpha_kept, log_q_phi,
                                                   d_x, segment_ends, n_cat, m, sum_of_alpha, Se_Sp_buffers,
                                                   Se_Sp_prior_type](auto& vi) mutable {
            const double adjoint = vi.adj();
            if constexpr (!phi_is_constant) {
                int start = 0;
                for (int s = 0; s <= m; ++s) {
                    const double w = Se_Sp_buffers.segments[s];
                    const double d_segment = (w != 1.0) ? (w - 1.0) * std::exp(-log_q_phi[s]) : 0.0;
                    for (int n = start; n < segment_ends[s]; ++n) {
                        const double d_dirichlet = (alpha_kept[n] != 1.0) ? (alpha_kept[n] - 1.0) / phi_kept[n]
                                                                          : 0.0;
                        phi_vars.coeffRef(n).adj() += adjoint * (d_dirichlet - d_segment);
                    }
                    start = segment_ends[s];
                }
            }
            if constexpr (!log_J_is_constant) {
                for (int j = 0; j < m; ++j) log_J_vars.coeffRef(j).adj() += adjoint;
            }
            if constexpr (!x_is_constant) {
                for (int j = 0; j < m; ++j) x_vars.coeffRef(j).adj() += adjoint * d_x[j];
            }
            if constexpr (!alpha_is_constant) {
                const double digamma_of_the_sum = stan::math::digamma(sum_of_alpha);
                for (int n = 0; n < n_cat; ++n) {
                    alpha_vars.coeffRef(n).adj() += adjoint * (std::log(phi_kept[n]) + digamma_of_the_sum -
                                                               stan::math::digamma(alpha_kept[n]));
                }
            }
            Se_Sp::add_the_partials_of_the_hyperparameters(Se_Sp_prior_type, Se_Sp_buffers, adjoint,
                                                           segments_vars, A_vars, beta_vars, link_mean_SD_vars);
        });
    }

}

}  // namespace fused_stan_math

using fused_stan_math::nico_dirichlet_with_Se_Sp_prior_lpdf;




// ---- the Se / Sp prior at m link values in any order (global scope)
namespace fused_stan_math {

template <bool propto, typename T_x, typename T_A, typename T_beta, typename T_link_mean_SD,
          stan::require_eigen_col_vector_t<T_x>* = nullptr>
inline stan::return_type_t<T_x, T_A, T_beta, T_link_mean_SD>
nico_Se_Sp_prior_at_chosen_cutoffs_lpdf(const T_x& x, const int Se_Sp_prior_type, const T_A& prior_Se_Sp_dirichlet_A,
                                        const T_beta& prior_Se_Sp_beta,
                                        const T_link_mean_SD& prior_Se_Sp_link_mean_SD, const int link,
                                        std::ostream* /*pstream__*/) {

    namespace Se_Sp = fused_stan_math::Se_Sp_prior_at_chosen_cutoffs;
    using stan::math::var;
    static constexpr const char* function_name = "nico_Se_Sp_prior_at_chosen_cutoffs_lpdf";
    const int m = static_cast<int>(x.size());
    const int type_of_H = (Se_Sp_prior_type == 1 || Se_Sp_prior_type == 2) ? Se_Sp_prior_type : 4;
    constexpr bool x_is_constant = stan::is_constant_all<T_x>::value;
    constexpr bool hyperparameters_are_constant = stan::is_constant_all<T_A, T_beta, T_link_mean_SD>::value;
    constexpr bool everything_is_constant = x_is_constant && hyperparameters_are_constant;
    constexpr bool include_the_constants = !propto || !hyperparameters_are_constant;
    if (m == 0) return 0.0;
    if constexpr (propto && everything_is_constant) return 0.0;
    const auto& x_ref = stan::math::to_ref(x);
    //// (the buffers: x, d_x (m), then the Se / Sp prior's; the segment weights are not used)
    const std::size_t n_buffer_values = 2 * static_cast<std::size_t>(m) +
                                        Se_Sp::buffers_of_the_Se_Sp_prior::n_values(m);
    std::vector<double> buffer_without_autodiff;
    double* buffers = nullptr;
    if constexpr (everything_is_constant) {
        buffer_without_autodiff.resize(n_buffer_values);
        buffers = buffer_without_autodiff.data();
    } else {
        buffers = stan::math::ChainableStack::instance_->memalloc_.alloc_array<double>(n_buffer_values);
    }
    double* x_kept = buffers;
    double* d_x = x_kept + m;
    Se_Sp::buffers_of_the_Se_Sp_prior Se_Sp_buffers;
    Se_Sp_buffers.place(d_x + m, m);
    for (int j = 0; j < m; ++j) {
        x_kept[j] = stan::math::value_of(x_ref.coeff(j));
        d_x[j] = 0.0;
    }
    const auto& segments_unused = Eigen::Matrix<double, -1, 1>::Ones(m + 1);
    Se_Sp::values_of_the_hyperparameters(function_name, segments_unused, type_of_H, prior_Se_Sp_dirichlet_A,
                                         prior_Se_Sp_beta, prior_Se_Sp_link_mean_SD, Se_Sp_buffers);
    if (type_of_H == 1) {
        if (!Se_Sp::values_of_one_set(x_kept, link, Se_Sp_buffers)) {
            throw std::domain_error(std::string(function_name) + ": x is not increasing");
        }
    } else {
        Se_Sp::values_of_one_set_in_any_order(x_kept, link, Se_Sp_buffers);
    }
    double sum = 0.0;
    for (int j = 0; j < m; ++j) {
        sum += Se_Sp_buffers.log_f[j];
        d_x[j] += (link == 1) ? -x_kept[j] : (std::exp(Se_Sp_buffers.lS[j]) - std::exp(Se_Sp_buffers.lF[j]));
    }
    sum += Se_Sp::Se_Sp_prior(x_kept, link, type_of_H, Se_Sp_buffers, include_the_constants,
                              !hyperparameters_are_constant, d_x);
    if constexpr (everything_is_constant) {
        return sum;
    } else {
        stan::arena_t<Eigen::Matrix<var, -1, 1>> x_vars;
        stan::arena_t<Eigen::Matrix<var, -1, 1>> segments_vars;
        stan::arena_t<Eigen::Matrix<var, -1, 1>> A_vars;
        stan::arena_t<Eigen::Matrix<var, -1, -1>> beta_vars;
        stan::arena_t<Eigen::Matrix<var, -1, -1>> link_mean_SD_vars;
        if constexpr (!x_is_constant) x_vars = x_ref;
        if constexpr (!stan::is_constant_all<T_A>::value) A_vars = prior_Se_Sp_dirichlet_A;
        if constexpr (!stan::is_constant_all<T_beta>::value) beta_vars = prior_Se_Sp_beta;
        if constexpr (!stan::is_constant_all<T_link_mean_SD>::value) link_mean_SD_vars = prior_Se_Sp_link_mean_SD;
        return stan::math::make_callback_var(sum, [x_vars, segments_vars, A_vars, beta_vars, link_mean_SD_vars, d_x,
                                                   m, Se_Sp_buffers, type_of_H](auto& vi) mutable {
            if constexpr (!x_is_constant) {
                for (int j = 0; j < m; ++j) x_vars.coeffRef(j).adj() += vi.adj() * d_x[j];
            }
            Se_Sp::add_the_partials_of_the_hyperparameters(type_of_H, Se_Sp_buffers, vi.adj(), segments_vars, A_vars,
                                                           beta_vars, link_mean_SD_vars);
        });
    }

}

}  // namespace fused_stan_math

using fused_stan_math::nico_Se_Sp_prior_at_chosen_cutoffs_lpdf;


// ---- Fused HSROC between-cutoffs segment + category conditional Dirichlet prior (multi-width k_* kernels) ----
namespace fused_stan_math {

template <bool propto,
          typename T_pop,
          typename T_phi,
          typename T_xb,
          typename T_sc,
          typename T_pa,
          stan::require_eigen_col_vector_t<T_pop>* = nullptr,
          stan::require_eigen_col_vector_t<T_phi>* = nullptr,
          stan::require_eigen_col_vector_t<T_xb>* = nullptr,
          stan::require_eigen_col_vector_t<T_sc>* = nullptr,
          stan::require_eigen_col_vector_t<T_pa>* = nullptr>
inline stan::return_type_t<T_pop, T_phi, T_xb, T_sc>
nico_hsroc_between_cutoffs_prior_lpdf(
        const T_pop& popmed_and_log_J,
        const T_phi& dirichlet_phi,
        const T_xb& Xbeta_of_line,
        const T_sc& scale_of_line,
        const T_pa& prior_dirichlet_alpha,
        const std::vector<int>& chosen_cutoffs_p,
        const int joint_link_prior,
        const double c_shift_segment_ref,
        const int use_probit_link,
        std::ostream* /*pstream__*/) {
    using stan::math::var;
    static constexpr const char* fn_name = "nico_hsroc_between_cutoffs_prior_lpdf";
    static constexpr double HALF_LOG_2PI = 0.918938533204672741780329736406;
    static constexpr double LOG_2        = 0.693147180559945309417232121458;

    const int n_cat = static_cast<int>(dirichlet_phi.size());
    const int n_thr = n_cat - 1;
    const int m_p   = static_cast<int>(chosen_cutoffs_p.size());
    if (n_thr < 1 || m_p >= n_thr) return 0.0;
    if (static_cast<int>(popmed_and_log_J.size()) != 2 * n_thr) {
        throw std::invalid_argument(std::string(fn_name) + ": popmed_and_log_J must have length 2 * n_thr");
    }

    constexpr bool pop_is_const = stan::is_constant_all<T_pop>::value;
    constexpr bool phi_is_const = stan::is_constant_all<T_phi>::value;
    constexpr bool xb_is_const  = stan::is_constant_all<T_xb>::value;
    constexpr bool sc_is_const  = stan::is_constant_all<T_sc>::value;
    const bool xb_active = (!xb_is_const) && (joint_link_prior != 1);
    const bool sc_active = (!sc_is_const) && (joint_link_prior != 1);
    constexpr bool is_any_var   = !pop_is_const || !phi_is_const || !xb_is_const || !sc_is_const;

    const auto& pop_ref = stan::math::to_ref(popmed_and_log_J);
    const auto& phi_ref = stan::math::to_ref(dirichlet_phi);
    const auto& xb_ref  = stan::math::to_ref(Xbeta_of_line);
    const auto& sc_ref  = stan::math::to_ref(scale_of_line);
    const auto& pa_ref  = stan::math::to_ref(prior_dirichlet_alpha);

    const std::size_t dbl_count = static_cast<std::size_t>(13) * n_thr + static_cast<std::size_t>(n_cat) + 8;
    std::vector<double> d_prim;
    double* d_buf = nullptr;
    if constexpr (!is_any_var) {
        d_prim.resize(dbl_count);
        d_buf = d_prim.data();
    } else {
        d_buf = stan::math::ChainableStack::instance_->memalloc_.alloc_array<double>(dbl_count);
    }

    double* P_buf        = d_buf;                // n_thr
    double* Pu_buf       = P_buf + n_thr;        // n_thr
    double* dP_buf       = Pu_buf + n_thr;       // n_thr
    double* dPu_buf      = dP_buf + n_thr;       // n_thr
    double* dlogJP_buf   = dPu_buf + n_thr;      // n_thr
    double* dpop_link    = dlogJP_buf + n_thr;   // n_thr
    double* dxb0_k       = dpop_link + n_thr;    // n_thr
    double* dxb1_k       = dxb0_k + n_thr;       // n_thr
    double* dsc0_k       = dxb1_k + n_thr;       // n_thr
    double* dsc1_k       = dsc0_k + n_thr;       // n_thr
    double* d_J_lsp_buf  = dsc1_k + n_thr;       // n_thr
    double* d_J_lse_buf  = d_J_lsp_buf + n_thr;  // n_thr
    double* w_sp_buf     = d_J_lse_buf + n_thr;  // n_thr
    double* dphi_buf     = w_sp_buf + n_thr;     // n_cat

    const double xb0 = stan::math::value_of(xb_ref.coeff(0));
    const double xb1 = stan::math::value_of(xb_ref.coeff(1));
    const double sc0 = stan::math::value_of(sc_ref.coeff(0));
    const double sc1 = stan::math::value_of(sc_ref.coeff(1));

    double d_xb0_from_logJ = 0.0;
    double d_xb1_from_logJ = 0.0;
    double d_sc0_from_logJ = 0.0;
    double d_sc1_from_logJ = 0.0;
    std::vector<double> work_f(4 * n_thr);
    std::vector<double> work_l(4 * n_thr);

    // Precompute P, P_upper, log_J_P and exact analytical derivatives
    // Re-use dlogJP_buf to hold log_J_P during forward pass, then clear for reverse adjoints!
    double* log_J_P = dlogJP_buf;
    if (joint_link_prior == 1) {
        for (int k = 0; k < n_thr; ++k) {
            const double pop = stan::math::value_of(pop_ref.coeff(k));
            const double lr  = pop - c_shift_segment_ref;
            work_f[k]         = lr;
            work_f[n_thr + k] = -lr;
            work_l[k]         = lr;
            work_l[n_thr + k] = -lr;
        }
        if (use_probit_link == 1) {
            k_Phi(work_f.data(), 2 * n_thr);
            for (int k = 0; k < n_thr; ++k) {
                const double lr = work_l[k];
                const double log_f = -0.5 * lr * lr - HALF_LOG_2PI;
                const double f_k = std::exp(log_f);
                P_buf[k]  = work_f[k];
                Pu_buf[k] = work_f[n_thr + k];
                log_J_P[k] = log_f + stan::math::value_of(pop_ref.coeff(n_thr + k));
                dP_buf[k] = f_k;
                dPu_buf[k] = -f_k;
                dpop_link[k] = -lr;
            }
        } else {
            k_inv_logit(work_f.data(), 2 * n_thr);
            k_log_inv_logit(work_l.data(), 2 * n_thr);
            for (int k = 0; k < n_thr; ++k) {
                const double pk  = work_f[k];
                const double puk = work_f[n_thr + k];
                const double fk  = pk * puk;
                P_buf[k]  = pk;
                Pu_buf[k] = puk;
                log_J_P[k] = work_l[k] + work_l[n_thr + k] + stan::math::value_of(pop_ref.coeff(n_thr + k));
                dP_buf[k] = fk;
                dPu_buf[k] = -fk;
                dpop_link[k] = puk - pk;
            }
        }
    } else {
        const double log_sc0 = std::log(sc0);
        const double log_sc1 = std::log(sc1);
        for (int k = 0; k < n_thr; ++k) {
            const double pop  = stan::math::value_of(pop_ref.coeff(k));
            const double lsp  = (pop + xb0) * sc0;
            const double lse  = (pop + xb1) * sc1;
            work_f[k]             = lsp;
            work_f[n_thr + k]     = -lsp;
            work_f[2 * n_thr + k] = lse;
            work_f[3 * n_thr + k] = -lse;
            work_l[k]             = lsp;
            work_l[n_thr + k]     = -lsp;
            work_l[2 * n_thr + k] = lse;
            work_l[3 * n_thr + k] = -lse;
        }
        std::vector<double> orig_l = work_l;
        if (use_probit_link == 1) {
            k_Phi(work_f.data(), 4 * n_thr);
            for (int k = 0; k < n_thr; ++k) {
                const double lsp = orig_l[k];
                const double lse = orig_l[2 * n_thr + k];
                const double lf_sp = -0.5 * lsp * lsp - HALF_LOG_2PI;
                const double lf_se = -0.5 * lse * lse - HALF_LOG_2PI;
                const double f_sp  = std::exp(lf_sp);
                const double f_se  = std::exp(lf_se);
                P_buf[k]  = 0.5 * (work_f[k] + work_f[2 * n_thr + k]);
                Pu_buf[k] = 0.5 * (work_f[n_thr + k] + work_f[3 * n_thr + k]);
                const double a = lf_sp + log_sc0;
                const double b = lf_se + log_sc1;
                const double hi = std::max(a, b);
                log_J_P[k] = hi + std::log1p(std::exp(-std::abs(a - b))) - LOG_2
                             + stan::math::value_of(pop_ref.coeff(n_thr + k));
                const double w_sp = (a >= b) ? (1.0 / (1.0 + std::exp(b - a))) : (std::exp(a - b) / (1.0 + std::exp(a - b)));
                const double w_se = 1.0 - w_sp;
                const double d_J_lsp = w_sp * (-lsp);
                const double d_J_lse = w_se * (-lse);
                const double pop = stan::math::value_of(pop_ref.coeff(k));
                dP_buf[k]   = 0.5 * f_sp * sc0 + 0.5 * f_se * sc1; // dP/dpop
                dPu_buf[k]  = -dP_buf[k];
                dpop_link[k] = d_J_lsp * sc0 + d_J_lse * sc1;
                dxb0_k[k]   = 0.5 * f_sp * sc0; // dP/dxb0 (dPu = -dxb0_k)
                dxb1_k[k]   = 0.5 * f_se * sc1;
                dsc0_k[k]   = 0.5 * f_sp * (pop + xb0);
                dsc1_k[k]   = 0.5 * f_se * (pop + xb1);
                d_J_lsp_buf[k] = d_J_lsp;
                d_J_lse_buf[k] = d_J_lse;
                w_sp_buf[k]    = w_sp;
            }
        } else {
            k_inv_logit(work_f.data(), 4 * n_thr);
            k_log_inv_logit(work_l.data(), 4 * n_thr);
            for (int k = 0; k < n_thr; ++k) {
                const double s_sp  = work_f[k], sm_sp = work_f[n_thr + k];
                const double s_se  = work_f[2 * n_thr + k], sm_se = work_f[3 * n_thr + k];
                const double f_sp  = s_sp * sm_sp;
                const double f_se  = s_se * sm_se;
                const double lf_sp = work_l[k] + work_l[n_thr + k];
                const double lf_se = work_l[2 * n_thr + k] + work_l[3 * n_thr + k];
                P_buf[k]  = 0.5 * (s_sp + s_se);
                Pu_buf[k] = 0.5 * (sm_sp + sm_se);
                const double a = lf_sp + log_sc0;
                const double b = lf_se + log_sc1;
                const double hi = std::max(a, b);
                log_J_P[k] = hi + std::log1p(std::exp(-std::abs(a - b))) - LOG_2
                             + stan::math::value_of(pop_ref.coeff(n_thr + k));
                const double w_sp = (a >= b) ? (1.0 / (1.0 + std::exp(b - a))) : (std::exp(a - b) / (1.0 + std::exp(a - b)));
                const double w_se = 1.0 - w_sp;
                const double d_J_lsp = w_sp * (sm_sp - s_sp);
                const double d_J_lse = w_se * (sm_se - s_se);
                const double pop = stan::math::value_of(pop_ref.coeff(k));
                dP_buf[k]    = 0.5 * f_sp * sc0 + 0.5 * f_se * sc1;
                dPu_buf[k]   = -dP_buf[k];
                dpop_link[k] = d_J_lsp * sc0 + d_J_lse * sc1;
                dxb0_k[k]    = 0.5 * f_sp * sc0;
                dxb1_k[k]    = 0.5 * f_se * sc1;
                dsc0_k[k]    = 0.5 * f_sp * (pop + xb0);
                dsc1_k[k]    = 0.5 * f_se * (pop + xb1);
                d_J_lsp_buf[k] = d_J_lsp;
                d_J_lse_buf[k] = d_J_lse;
                w_sp_buf[k]    = w_sp;
            }
        }
    }

    // Store per-P, per-Pu, and logJ contributions separately cleanly without aliasing!
    std::vector<double> adj_P(n_thr, 0.0);
    std::vector<double> adj_Pu(n_thr, 0.0);
    std::vector<unsigned char> logJ_in_sum(n_thr, 0);
    std::fill_n(dphi_buf, n_cat, 0.0);

    double sum_val = 0.0;
    int k_lower_end = 0;
    for (int j = 0; j <= m_p; ++j) {
        const int k_upper_end = (j < m_p) ? chosen_cutoffs_p[j] : n_cat;
        if (k_upper_end - k_lower_end > 1) {
            double P_seg = 0.0;
            int seg_p_hi = -1, seg_p_lo = -1, seg_pu_lo = -1, seg_pu_hi = -1;
            if (k_lower_end == 0) {
                P_seg = P_buf[k_upper_end - 1];
                seg_p_hi = k_upper_end - 1;
            } else if (k_upper_end == n_cat) {
                P_seg = Pu_buf[k_lower_end - 1];
                seg_pu_lo = k_lower_end - 1;
            } else if (P_buf[k_lower_end - 1] > 0.5) {
                P_seg = Pu_buf[k_lower_end - 1] - Pu_buf[k_upper_end - 1];
                seg_pu_lo = k_lower_end - 1;
                seg_pu_hi = k_upper_end - 1;
            } else {
                P_seg = P_buf[k_upper_end - 1] - P_buf[k_lower_end - 1];
                seg_p_hi = k_upper_end - 1;
                seg_p_lo = k_lower_end - 1;
            }

            double phi_seg = 0.0;
            double sum_pa = 0.0;
            for (int i = k_lower_end; i < k_upper_end; ++i) {
                phi_seg += stan::math::value_of(phi_ref.coeff(i));
                sum_pa  += stan::math::value_of(pa_ref.coeff(i));
            }
            for (int k = k_lower_end; k <= k_upper_end - 2; ++k) {
                sum_val += log_J_P[k];
                logJ_in_sum[k] = 1u;
            }
            const double sum_pa_m1 = sum_pa - 1.0;
            sum_val -= sum_pa_m1 * (std::log(P_seg) - std::log(phi_seg));
            if constexpr (is_any_var) {
                const double w_P = -sum_pa_m1 / P_seg;
                if (seg_p_hi  >= 0) adj_P[seg_p_hi]   += w_P;
                if (seg_p_lo  >= 0) adj_P[seg_p_lo]   -= w_P;
                if (seg_pu_lo >= 0) adj_Pu[seg_pu_lo] += w_P;
                if (seg_pu_hi >= 0) adj_Pu[seg_pu_hi] -= w_P;
                const double w_phi = sum_pa_m1 / phi_seg;
                for (int i = k_lower_end; i < k_upper_end; ++i) dphi_buf[i] += w_phi;
            }

            for (int i = k_lower_end; i < k_upper_end; ++i) {
                const double pa_i = stan::math::value_of(pa_ref.coeff(i));
                if (pa_i != 1.0) {
                    double P_cat = 0.0;
                    int cp_hi = -1, cp_lo = -1, cpu_lo = -1, cpu_hi = -1;
                    if (i == n_cat - 1) {
                        P_cat = Pu_buf[n_thr - 1];
                        cpu_lo = n_thr - 1;
                    } else if (i == 0) {
                        P_cat = P_buf[0];
                        cp_hi = 0;
                    } else if (P_buf[i - 1] > 0.5) {
                        P_cat = Pu_buf[i - 1] - Pu_buf[i];
                        cpu_lo = i - 1;
                        cpu_hi = i;
                    } else {
                        P_cat = P_buf[i] - P_buf[i - 1];
                        cp_hi = i;
                        cp_lo = i - 1;
                    }
                    const double phi_i = stan::math::value_of(phi_ref.coeff(i));
                    const double pa_m1 = pa_i - 1.0;
                    sum_val += pa_m1 * (std::log(P_cat) - std::log(phi_i));
                    if constexpr (is_any_var) {
                        const double w_cP = pa_m1 / P_cat;
                        if (cp_hi  >= 0) adj_P[cp_hi]   += w_cP;
                        if (cp_lo  >= 0) adj_P[cp_lo]   -= w_cP;
                        if (cpu_lo >= 0) adj_Pu[cpu_lo] += w_cP;
                        if (cpu_hi >= 0) adj_Pu[cpu_hi] -= w_cP;
                        dphi_buf[i] -= pa_m1 / phi_i;
                    }
                }
            }
        }
        k_lower_end = k_upper_end;
    }

    if constexpr (!is_any_var) {
        return sum_val;
    } else {
        // Pre-accumulate non-popmed/xb/sc scalars cleanly on arena or lambda
        double* pop_d  = stan::math::ChainableStack::instance_->memalloc_.alloc_array<double>(n_thr);
        double* logj_d = stan::math::ChainableStack::instance_->memalloc_.alloc_array<double>(n_thr);
        double* phi_d  = stan::math::ChainableStack::instance_->memalloc_.alloc_array<double>(n_cat);
        double dx0 = 0.0, dx1 = 0.0, ds0 = 0.0, ds1 = 0.0;

        for (int k = 0; k < n_thr; ++k) {
            const double in_sum = (logJ_in_sum[k] != 0u) ? 1.0 : 0.0;
            const double aP = adj_P[k], aPu = adj_Pu[k];
            const double diff_P = aP - aPu;
            logj_d[k] = in_sum;
            if (joint_link_prior == 1) {
                // dP = f_k, dPu = -f_k
                pop_d[k] = diff_P * dP_buf[k] + in_sum * dpop_link[k];
            } else {
                const double pop = stan::math::value_of(pop_ref.coeff(k));
                const double d_J_lsp = d_J_lsp_buf[k];
                const double d_J_lse = d_J_lse_buf[k];
                const double w_sp_k  = w_sp_buf[k];
                const double w_se_k  = 1.0 - w_sp_k;
                pop_d[k] = diff_P * dP_buf[k] + in_sum * dpop_link[k];
                dx0 += diff_P * dxb0_k[k] + in_sum * d_J_lsp * sc0;
                dx1 += diff_P * dxb1_k[k] + in_sum * d_J_lse * sc1;
                ds0 += diff_P * dsc0_k[k] + in_sum * (d_J_lsp * (pop + xb0) + w_sp_k / sc0);
                ds1 += diff_P * dsc1_k[k] + in_sum * (d_J_lse * (pop + xb1) + w_se_k / sc1);
            }
        }
        for (int i = 0; i < n_cat; ++i) phi_d[i] = dphi_buf[i];

        stan::arena_t<Eigen::Matrix<var, -1, 1>> pop_vars, phi_vars, xb_vars, sc_vars;
        if constexpr (!pop_is_const) pop_vars = pop_ref;
        if constexpr (!phi_is_const) phi_vars = phi_ref;
        if constexpr (!xb_is_const)  if (xb_active) xb_vars = xb_ref;
        if constexpr (!sc_is_const)  if (sc_active) sc_vars = sc_ref;

        return stan::math::make_callback_var(
            sum_val,
            [pop_vars, phi_vars, xb_vars, sc_vars, pop_d, logj_d, phi_d,
             dx0, dx1, ds0, ds1, n_thr, n_cat, xb_active, sc_active](auto& vi) mutable {
                const double adj = vi.adj();
                if constexpr (!pop_is_const) {
                    for (int k = 0; k < n_thr; ++k) {
                        if (pop_d[k]  != 0.0) pop_vars.coeffRef(k).adj()         += adj * pop_d[k];
                        if (logj_d[k] != 0.0) pop_vars.coeffRef(n_thr + k).adj() += adj * logj_d[k];
                    }
                }
                if constexpr (!phi_is_const) {
                    for (int i = 0; i < n_cat; ++i) {
                        if (phi_d[i] != 0.0) phi_vars.coeffRef(i).adj() += adj * phi_d[i];
                    }
                }
                if constexpr (!xb_is_const) {
                    if (xb_active) {
                        xb_vars.coeffRef(0).adj() += adj * dx0;
                        xb_vars.coeffRef(1).adj() += adj * dx1;
                    }
                }
                if constexpr (!sc_is_const) {
                    if (sc_active) {
                        sc_vars.coeffRef(0).adj() += adj * ds0;
                        sc_vars.coeffRef(1).adj() += adj * ds1;
                    }
                }
            });
    }
}

}  // namespace fused_stan_math

using fused_stan_math::nico_hsroc_between_cutoffs_prior_lpdf;
