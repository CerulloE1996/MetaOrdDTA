#pragma once
//
// -| --------- The fused functions of the MetaOrdDTA models, on Stan math only: the user header ---------------------
//
// The one header to give cmdstanr (or BridgeStan, CmdStan) for the MetaOrdDTA files *_fused_Stan_math_only.stan
// (the MA and NMA models):
//
//     cmdstanr::cmdstan_model(stan_file,
//                             user_header   = <this file>,      (cmdstanr then runs stanc with allow-undefined)
//                             include_paths = <MetaOrdDTA's stan_functions folder>)
//
// Those Stan files call Stan's own built-ins for log, exp, Phi, normal_lpdf, ... and declare only the fused
// functions (no body): each is one autodiff node with its hand-written gradient, and their C++ is the scalar code
// of the fused (NicoSIMD) versions with Stan math in place of the SIMD kernels. Nothing here needs NicoStan,
// NicoSIMD or any -m flag: x86-64 with or without AVX, ARM, any CPU Stan runs on.
//
// The external functions each file gives (the ones the MetaOrdDTA Stan files declare):
//   stan_external_binomial_lpmf.hpp:              fused_binomial_lpmf, fused_strat_bivariate_binomial_lpmf
//   stan_external_binomial_counts_above_successive_cutpoints_of_a_cumulative_link.hpp:
//                                                 fused_binomial_counts_above_successive_cutpoints_of_a_cumulative_link_lpmf
//                                                 (the full name of fused_binomial_counts_above_cuts_cumul_link_lpmf)
//   stan_external_induced_dirichlet_prior_of_ordered_cutpoints.hpp:
//                                                 fused_induced_dirichlet_prior_of_ordered_cutpoints_lpdf (and _ragged_),
//                                                 fused_ord_biv_randthr_ind_dir_and_binomial_lpmf (and _precomputed_),
//                                                 fused_ord_hsroc_randthr_ind_dir_and_binomial_lpmf
//   stan_external_dirichlet_cutpoint_median.hpp:  fused_dirichlet_cutpoint_median_and_log_J,
//                                                 fused_dirichlet_cutpoint_medians_and_sum_of_log_J,
//                                                 fused_dirichlet_cutpoint_medians_and_log_J_of_each,
//                                                 fused_dirichlet_cutpoint_medians_and_log_J_at_cuts
//   stan_external_Box_Cox_lambda_log_density_of_the_flat_Dirichlet.hpp:
//                                                 fused_Box_Cox_lambda_log_density_of_the_flat_Dirichlet,
//                                                 fused_box_cox_C_raw_and_C
//   stan_external_induced_dirichlet_with_Se_Sp_prior.hpp:
//                                                 fused_induced_dirichlet_with_Se_Sp_prior_lpdf,
//                                                 fused_dirichlet_with_Se_Sp_prior_lpdf, fused_Se_Sp_prior_at_chosen_cutoffs_lpdf,
//                                                 fused_hsroc_between_cutoffs_prior_lpdf
//   stan_external_dirichlet_lpdf.hpp, stan_external_beta_lpdf.hpp, stan_external_lkj_corr_cholesky_lpdf.hpp:
//                                                 fused_dirichlet_lpdf, fused_beta_lpdf, fused_lkj_corr_cholesky_lpdf (used by
//                                                 the fused functions), fused_strat_bivariate_priors_lpdf
//   stan_external_log_det_jacobian_of_cutpoints_from_raw_differences.hpp:
//                                                 fused_log_det_jacobian_of_cutpoints_from_raw_differences,
//                                                 fused_raw_C_to_C_log_det_J_lp, fused_construct_C, fused_ord_biv_fixedthr_model_lp
//                                                 (fused_log_det_jac_of_cuts_from_raw_diffs under the full name: NMA)
//   stan_external_factorised_binomial_log_lik_of_ordinal_test_accuracy_data.hpp:
//                                                 fused_compute_log_lik_binomial_fact_data,
//                                                 fused_compute_NMA_log_lik_binomial_fact_data,
//                                                 fused_map_latent_surv_to_fixed_hetero_C (generated quantities, doubles)
//   stan_external_pinkney_partially_centered_cutpoints.hpp:
//                                                 fused_pinkney_C_and_log_J_from_alpha (and _precomputed, _both_ ...),
//                                                 fused_pinkney_C_and_log_J_from_pairs, fused_ord_biv_rc_priors_lpdf
//   stan_external_set0_marginal_at_chosen_cutoffs.hpp:
//                                                 fused_log_Set0_marginal_at_chosen_cutoffs (and _fast)
//
// ------------------------------------------------------------------------------------------------------------------

#include "stan_external_functions_on_Stan_math_only.hpp"
#include "stan_external_binomial_lpmf.hpp"
#include "stan_external_binomial_counts_above_successive_cutpoints_of_a_cumulative_link.hpp"
#include "stan_external_induced_dirichlet_prior_of_ordered_cutpoints.hpp"
#include "stan_external_dirichlet_cutpoint_median.hpp"
#include "stan_external_Box_Cox_lambda_log_density_of_the_flat_Dirichlet.hpp"
#include "stan_external_induced_dirichlet_with_Se_Sp_prior.hpp"
#include "stan_external_dirichlet_lpdf.hpp"
#include "stan_external_beta_lpdf.hpp"
#include "stan_external_lkj_corr_cholesky_lpdf.hpp"
#include "stan_external_log_det_jacobian_of_cutpoints_from_raw_differences.hpp"
#include "stan_external_factorised_binomial_log_lik_of_ordinal_test_accuracy_data.hpp"
#include "stan_external_pinkney_partially_centered_cutpoints.hpp"
#include "stan_external_set0_marginal_at_chosen_cutoffs.hpp"
