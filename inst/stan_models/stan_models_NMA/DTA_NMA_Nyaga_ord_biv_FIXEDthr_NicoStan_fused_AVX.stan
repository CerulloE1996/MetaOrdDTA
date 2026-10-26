//// DTA_NMA_Nyaga_Xu_FIXEDthr_NicoStan_fused_AVX.stan
//// The same model as MetaOrdDTA's DTA_NMA_Nyaga_Xu_FIXEDthr.stan (unchanged; its replaced lines kept here,
//// commented out), for NicoStan with NicoSIMD's AVX functions (MetaOrdDTA's NicoStan option; one file for AVX2
//// and AVX-512): the likelihood of every reported cutpoint of every test, study and class in ONE autodiff node
//// (NicoSIMD's general nico_binomial_counts_above_successive_cutpoints_of_a_cumulative_link_lpmf, the same node
//// as the MA models' and as Paper 1's chunked network meta-analysis's), the induced Dirichlet of each test's
//// cutpoints as NicoSIMD's general nico_induced_dirichlet_prior_of_ordered_cutpoints_lpdf, the normal priors and
//// the standard normal random effects as nico_normal, and the exp() (or softplus) of the cutpoints' raw
//// differences as nico_exp (or nico_log1p_exp). Not chunked and without nuisance blocks (the parameters are the
//// original's, in the same order). The terms of the likelihood are listed once, in transformed data. The terms
//// that are constants of the data are left out of the log density when the sampler asks for it (propto:
//// _lupdf / _lupmf, as the ~ statements do); with propto = FALSE the log density is the model's, every constant
//// included. The generated quantities are the model's, unchanged.
////


functions {
        ////
        //// Include files to compile the necessary custom (user-defined) Stan functions:
        ////
        #include "Stan_fns_basic.stan"
        #include "Stan_fns_corr.stan"
        #include "Stan_fns_ordinal.stan"
        #include "Stan_fns_log_lik.stan"
        #include "Stan_fns_ragged.stan"
        #include "Stan_fns_NMA.stan"
        #include "Stan_fns_model_fit.stan"
        #include "Stan_fns_NMA_bivariate.stan"
        #include "Stan_fns_Jacobian.stan"
        ////
        //// ---- NicoSIMD's AVX functions (declared without a body; stanc --allow-undefined, with NicoSIMD's
        ////      stan_external_functions_of_NicoSIMD.hpp as the user header):
        ////
        real nico_binomial_counts_above_successive_cutpoints_of_a_cumulative_link_lpmf(
              array[] int n, array[] int n_previous, vector latent, data vector is_first, data int link,
              data real bound);
        real nico_induced_dirichlet_prior_of_ordered_cutpoints_lpdf(vector C, vector alpha, data int link,
                                                                    data real location);
        vector nico_exp(vector x);
        vector nico_log1p_exp(vector x);
        real nico_log_det_jacobian_of_cutpoints_from_raw_differences(vector C_raw, data int softplus);
        real nico_log_det_jacobian_of_cutpoints_from_raw_differences(matrix C_raw, data int softplus);
        array[] matrix nico_map_latent_surv_to_fixed_hetero_C(data array[] vector C,
                                                              data matrix locations, data matrix scales,
                                                              data int n_studies, data array[] int n_obs_cutpoints,
                                                              data array[,,] int cutpoint_index);
        array[] matrix nico_map_latent_surv_to_fixed_hetero_C(data array[] vector C,
                                                              data matrix locations,
                                                              data int n_studies, data array[] int n_obs_cutpoints,
                                                              data array[,,] int cutpoint_index);
        array[,] matrix nico_compute_NMA_log_lik_binomial_fact_data(data array[] matrix latent_surv,
                                                                    data int use_probit_link, data int n_thr,
                                                                    data int n_studies, data array[,,] int x_2,
                                                                    data array[,,] int n, data array[,] int N_total,
                                                                    data array[] int n_obs_cutpoints,
                                                                    data array[] int indicator_index_test_t_in_study,
                                                                    data array[] int K_fold_CV_indicator);
        real nico_normal_lpdf(vector y, real mu, real sigma);
        real nico_normal_lpdf(vector y, real mu, vector sigma);
        real nico_normal_lpdf(vector y, vector mu, vector sigma);
        real nico_lkj_corr_cholesky_lpdf(matrix L, data real eta);
        real nico_Phi(real x);
        real nico_inv_logit(real x);
}

 
data {
          int<lower=1> n_studies;
          int<lower=1> n_index_tests;
          ////
          array[n_index_tests] int<lower=1> n_studies_per_test;
          ////
          array[n_index_tests] int<lower=1> n_thr;
          array[n_index_tests] int<lower=2> n_cat;
          int<lower=1> n_thr_max; // corresponding to the test with the most thresholds/cutpoints. 
          int<lower=2> n_cat_max; // corresponding to the test with the most thresholds/cutpoints. 
          ////
          array[n_index_tests, n_studies] int n_obs_cutpoints; //// OBSERVED cutpoints for test t in study s
          ////
          array[n_studies, n_index_tests] int<lower=0, upper=1> indicator_index_test_in_study;   // Binary indicator if test t is in study s
          ////
          //// ---- Data:
          ////
          array[n_index_tests, 2, n_studies, n_thr_max + 1] int x;
          array[n_index_tests, 2, n_studies, n_thr_max + 1] int cutpoint_index;
          // ////
          // //// ---- Covariates (w/ missing cts):
          // ////
          // array[n_index_tests] int n_missing_cov_nd;  // number of missing values per test
          // array[n_index_tests] int n_missing_cov_d;   // number of missing values per test
          // array[n_index_tests, n_studies] int<lower=0, upper=1> has_missing_nd;  // whether a study has ANY missing values
          // array[n_index_tests, n_studies] int<lower=0, upper=1> has_missing_d;   // whether a study has ANY missing values
          // array[n_index_tests, n_studies, max(n_covariates_nd)] int<lower=0, upper=1> is_missing_nd; // which SPECIFIC covariates are missing
          // array[n_index_tests, n_studies, max(n_covariates_d)]  int<lower=0, upper=1> is_missing_d;  // which SPECIFIC covariates are missing
          // int n_covariates_missing_nd;  // total missing values across all tests
          // int n_covariates_missing_d;   // total missing values across all tests
          ////
          array[n_index_tests] int n_covariates_nd; // total INCLUDING ones w/ missing
          array[n_index_tests] int n_covariates_d;  // total INCLUDING ones w/ missing
          array[n_index_tests] matrix[n_studies, max(n_covariates_nd)] X_nd; // "observed" with 999's for missing
          array[n_index_tests] matrix[n_studies, max(n_covariates_d)]  X_d;   // "observed" with 999's for missing
          array[n_index_tests] vector[max(n_covariates_nd)] baseline_case_nd;
          array[n_index_tests] vector[max(n_covariates_d)]  baseline_case_d;
          ////
          //// ---- Priors for beta:
          ////
          array[n_index_tests] matrix[2, max(max(n_covariates_nd), max(n_covariates_d))] prior_beta_mu_mean;
          array[n_index_tests] matrix<lower=0.0>[2, max(max(n_covariates_nd), max(n_covariates_d))] prior_beta_mu_SD;
          matrix<lower=0.0>[n_index_tests, 2] prior_beta_tau_SD;
          vector<lower=0.0>[2] prior_beta_sigma_SD; //// "sigma's (Nyaga et al. notation) are shared between tests.
          ////
          //// Priors (and possible restrictions) for between-study correlations:
          ////
          real<lower=-1.0, upper=+1.0> beta_corr_lb;
          real<lower=beta_corr_lb, upper=+1.0>  beta_corr_ub;
          real<lower=0.0>  prior_beta_corr_LKJ;
          ////
          //// ---- Induced-Dirichlet priors:
          ////
          array[2] matrix<lower=0.0>[n_index_tests, n_thr_max + 1] prior_dirichlet_alpha;
          ////
          //// ---- Other:
          ////
          int<lower=0, upper=1> softplus;
          ////
          int use_probit_link;
          ////
          int hetero_sigma;
          int compound_symmetry;
          ////
          array[n_studies] int<lower=0, upper=1> K_fold_CV_indicator;
}

transformed data {
          ////
          int n_covariates_max = max(max(n_covariates_nd), max(n_covariates_d));
          //
          array[n_index_tests, 2, n_studies, n_thr_max] int x_2;
          array[n_index_tests, 2, n_studies] int N_total;
          array[n_index_tests, 2, n_studies, n_thr_max] int n;
          //// ---- Initialize with -1 (missing indicator):
          for (t in 1:n_index_tests) {
              for (c in 1:2) {
                for (s in 1:n_studies) {
                    N_total[t, c, s] = -1;
                    for (k in 1:n_thr_max) {
                      x_2[t, c, s, k]  = -1;
                      n[t, c, s, k]    = -1;
                    }
                }
              }
          }
          ////
          for (t in 1:n_index_tests) {
                for (s in 1:n_studies) {
                          for (c in 1:2) {
                              N_total[t, c, s] = x[t, c, s, 1];
                              for (k in 1:n_thr[t]) {
                                 x_2[t, c, s, k] = x[t, c, s, k + 1];
                              }
                          }
                          for (c in 1:2) {
                               n[t, c][s, 1] = N_total[t, c, s];
                               for (k in 2:n_obs_cutpoints[t, s]) {
                                          n[t, c, s, k] = x_2[t, c, s, k - 1];
                               }
                          }
                }
          }
          ////
          int n_total_C_per_class_if_fixed = sum(n_thr); // total # of fixed-effect cutpoints
          int n_tests = n_index_tests;
          ////
          //// ---- Calculate indices:
          ////
          array[n_tests] int start_index = calculate_start_indices(n_thr, n_tests);
          array[n_tests] int end_index   = calculate_end_indices(n_thr, n_tests, start_index);
          ////
          int n_test_pairs = n_index_tests * (n_index_tests - 1) / 2;
          ////
          int n_total_study_thresh = 0;
          for (t in 1:n_index_tests) {
              n_total_study_thresh += n_studies_per_test[t] * n_thr[t];
          }
          ////
          int n_total_obs = 0;
          for (t in 1:n_index_tests) {
            for (s in 1:n_studies) {
              if (indicator_index_test_in_study[s, t] == 1) {
                  n_total_obs += 1;
              }
            }
          }
          ////
          //// ---- The terms of the likelihood, one per reported cutpoint of each study (with the test, and
          ////      counted in the cross-validation) and class, in the order of the tests, the classes, the studies
          ////      and their reported cutpoints, made once here: the count above the cutpoint, the count above the
          ////      study's previous reported cutpoint (or its total), whether the cutpoint is the study's first
          ////      reported one, the place of the test, class and study in the locations of all tests (to_vector()
          ////      of each Xbeta[t], one after the other), and the place of the cutpoint in the cutpoints of both
          ////      classes (C_vec[1], then C_vec[2]):
          ////
          int n_terms = 0;
          for (t in 1:n_index_tests) {
                for (c in 1:2) {
                      for (s in 1:n_studies) {
                            if ((indicator_index_test_in_study[s, t] == 1) && (K_fold_CV_indicator[s] == 1)) {
                                  n_terms += n_obs_cutpoints[t, s];
                            }
                      }
                }
          }
          array[n_terms] int count_above_the_cutpoint_of_each_term;
          array[n_terms] int count_above_the_previous_cutpoint_of_each_term;
          vector[n_terms] is_the_first_reported_cutpoint_of_each_term;
          array[n_terms] int index_of_the_location_of_each_term;
          array[n_terms] int index_of_the_cutpoint_of_each_term;
          {
                int term = 0;
                for (t in 1:n_index_tests) {
                      for (c in 1:2) {
                            for (s in 1:n_studies) {
                                  if ((indicator_index_test_in_study[s, t] == 1) &&
                                      (K_fold_CV_indicator[s] == 1)) {
                                        for (cut_i in 1:n_obs_cutpoints[t, s]) {
                                              int k = cutpoint_index[t, c, s, cut_i + 1];
                                              if (k < 1) reject("a reported cutpoint without its index: test ",
                                                                t, ", class ", c, ", study ", s,
                                                                ", reported cutpoint ", cut_i);
                                              term += 1;
                                              count_above_the_cutpoint_of_each_term[term] = x_2[t, c, s, cut_i];
                                              count_above_the_previous_cutpoint_of_each_term[term] =
                                                    n[t, c, s, cut_i];
                                              is_the_first_reported_cutpoint_of_each_term[term] =
                                                    (cut_i == 1) ? 1.0 : 0.0;
                                              index_of_the_location_of_each_term[term] =
                                                    (t - 1) * 2 * n_studies + (c - 1) * n_studies + s;
                                              index_of_the_cutpoint_of_each_term[term] =
                                                    (c - 1) * n_total_C_per_class_if_fixed + start_index[t] +
                                                    k - 1;
                                        }
                                  }
                            }
                      }
                }
          }

}


parameters {
          matrix[2, n_studies] beta_eta_z;
          array[2] vector[n_total_obs] beta_delta_z_vec; //// Standard normal RVs for test-specific effects
          ////
          array[n_index_tests] matrix[2, n_covariates_max] beta_mu; 
          ////
          // array[2] vector<lower=-7.5, upper=2.5>[n_total_C_per_class_if_fixed] C_raw_vec;  //// RAW LOG-DIFFS - Global cutpoints for each test (staggered array/matrix using "n_thr[t]" to index correctly)
          array[2] vector<lower=-15, upper=2.5>[n_total_C_per_class_if_fixed] C_raw_vec;  //// RAW LOG-DIFFS - Global cutpoints for each test (staggered array/matrix using "n_thr[t]" to index correctly)
          ////
          //// ---- "NMA" params:
          ////
          vector<lower=0.0>[2] beta_sigma;  //// Var component - Between-study SD (Nyaga's σ)
          ////
          matrix<lower=0.0>[n_index_tests, 2] beta_tau_raw;
          ////
          real<lower=beta_corr_lb, upper=beta_corr_ub> beta_corr;  //// between-study corr (possibly restricted) 
          // ////
          // //// ---- For missing covariates:
          // ////
          // vector[n_covariates_missing_nd] missing_cts_cov_vals_nd;
          // vector[n_covariates_missing_d]  missing_cts_cov_vals_d;
          // ////
          // vector[n_covariates_max] mu_cov_pop;
          // vector<lower=0>[n_covariates_max] sigma_cov_pop;
}


transformed parameters {
          // array[n_index_tests] matrix[n_studies, max(n_covariates_nd)] X_nd_complete;  // complete
          // array[n_index_tests] matrix[n_studies, max(n_covariates_d)]  X_d_complete;   // complete
          // for (t in 1:n_index_tests) {
          //   X_nd_complete[t] = X_nd[t]; // Initialize with observed
          //   X_d_complete[t]  = X_d[t];   // Initialize with observed
          // }
          // //// ---- Fill in missing cts covariate values:
          // {
          //   int pos_nd = 1;
          //   int pos_d = 1;
          //   for (t in 1:n_index_tests) {
          //     for (s in 1:n_studies) {
          //         if (has_missing_nd[t, s] == 1) {s
          //           for (cov in 1:n_covariates_nd[t]) {
          //             if (is_missing_nd[t, s, cov] == 1) {
          //               X_nd_complete[t][s, cov] = missing_cts_cov_vals_nd[pos_nd];
          //               pos_nd += 1;
          //             }
          //           }
          //         }
          //         if (has_missing_d[t, s] == 1) {
          //           for (cov in 1:n_covariates_d[t]) {
          //             if (is_missing_d[t, s, cov] == 1) {
          //               X_d_complete[t][s, cov] = missing_cts_cov_vals_d[pos_d];
          //               pos_d += 1;
          //             }
          //           }
          //         }
          //     }
          //   }
          // }
          ////
          matrix<lower=0.0>[n_index_tests, 2] beta_tau;
          if (compound_symmetry) {
              for (t in 1:n_index_tests) {
                int test_1_idx = 1;
                for (c in 1:2) {
                  beta_tau[t, c] = beta_tau_raw[test_1_idx, c];
                }
              }
          } else {
              beta_tau = beta_tau_raw;
          }
          ////
          //// ---- Construct (global) cutpoints:
          ////
          array[2] vector[n_total_C_per_class_if_fixed] C_vec;  //// Global cutpoints for each test ("staggered" array/matrix using "n_thr[t]" to index correctly)
          ////
          for (t in 1:n_index_tests) {
              for (c in 1:2) {
                    int n_thr_t = n_thr[t];
                    vector[n_thr_t] C_raw_vec_test_t = get_test_values(C_raw_vec[c], start_index, end_index, t);
                    // vector[n_thr_t] C_vec_test_t =  construct_C(C_raw_vec_test_t, softplus);
                    //// (as construct_C: the first cutpoint, then the cumulative sums of the exp() (or softplus) of
                    ////  the raw differences, with NicoSIMD's nico_exp / nico_log1p_exp)
                    vector[n_thr_t] C_vec_test_t;
                    C_vec_test_t[1] = C_raw_vec_test_t[1];
                    if (n_thr_t > 1) {
                          vector[n_thr_t - 1] increments_of_the_cutpoints =
                                (softplus == 1) ? nico_log1p_exp(C_raw_vec_test_t[2:n_thr_t]) :
                                                  nico_exp(C_raw_vec_test_t[2:n_thr_t]);
                          for (k in 2:n_thr_t) {
                                C_vec_test_t[k] = C_vec_test_t[k - 1] + increments_of_the_cutpoints[k - 1];
                          }
                    }
                    C_vec[c] = update_test_values(C_vec[c], C_vec_test_t, start_index, end_index, t);
              }
          }
          ////
          //// ---- Construct simple 2x2 (bivariate) between-study corr matrices:
          ////
          cholesky_factor_corr[2] beta_L_Omega = make_bivariate_L_Omega(beta_corr);
          cholesky_factor_cov[2]  beta_L_Sigma = diag_pre_multiply(beta_sigma, beta_L_Omega);
          ////
          //// ---- "NMA" params:
          ////
          array[n_index_tests] matrix[n_studies, 2] Xbeta = init_array_of_matrices(n_studies, 2, n_index_tests, 0.0); // NOT local
          ////
          {
              {
                int counter = 0;
                 for (s in 1:n_studies) {
                     vector[2] beta_eta = beta_L_Sigma * beta_eta_z[1:2, s];
                     for (t in 1:n_index_tests) { 
                      if (indicator_index_test_in_study[s, t] == 1) {
                          counter += 1;
                          vector[2] beta_delta;
                          vector[2] beta_random;
                          for (c in 1:2) {
                              beta_delta[c] = beta_tau[t, c] * beta_delta_z_vec[c][counter];
                          }
                          ////
                          beta_random = beta_delta + beta_eta;
                          ////
                          Xbeta[t][s, 1] = X_nd[t][s, 1:n_covariates_nd[t]] * to_vector(beta_mu[t][1, 1:n_covariates_nd[t]]) + 
                                           beta_random[1];
                          Xbeta[t][s, 2] = X_d[t][s, 1:n_covariates_d[t]]   * to_vector(beta_mu[t][2, 1:n_covariates_d[t]])  + 
                                           beta_random[2];
                      }
                  }
                }
              }
          }
}
 
 
model { 
        //// 
        //// ---- Priors:
        ////
        for (t in 1:n_index_tests) {
            // beta_mu[t][1, 1:n_covariates_nd[t]] ~ normal(
            //                                       prior_beta_mu_mean[t][1, 1:n_covariates_nd[t]], 
            //                                       prior_beta_mu_SD[t][1, 1:n_covariates_nd[t]]);  
            // beta_mu[t][2, 1:n_covariates_d[t]]  ~ normal(
            //                                       prior_beta_mu_mean[t][2, 1:n_covariates_d[t]],  
            //                                       prior_beta_mu_SD[t][2, 1:n_covariates_d[t]]);
            to_vector(beta_mu[t][1, 1:n_covariates_nd[t]]) ~ nico_normal(
                                                  to_vector(prior_beta_mu_mean[t][1, 1:n_covariates_nd[t]]),
                                                  to_vector(prior_beta_mu_SD[t][1, 1:n_covariates_nd[t]]));
            to_vector(beta_mu[t][2, 1:n_covariates_d[t]])  ~ nico_normal(
                                                  to_vector(prior_beta_mu_mean[t][2, 1:n_covariates_d[t]]),
                                                  to_vector(prior_beta_mu_SD[t][2, 1:n_covariates_d[t]]));
        }
        ////
        // to_vector(beta_tau_raw) ~ normal(0.0, to_vector(prior_beta_tau_SD));  //// delta_{s, i, t} ~ normal(0, tau_{i, t}):
        //// (delta_{s, i, t} ~ normal(0, tau_{i, t}))
        to_vector(beta_tau_raw) ~ nico_normal(0.0, to_vector(prior_beta_tau_SD));
        beta_L_Omega ~ nico_lkj_corr_cholesky(prior_beta_corr_LKJ);
        ////
        // for (c in 1:2) {
        //     beta_sigma[c] ~ normal(0.0, prior_beta_sigma_SD[c]);
        // }
        //// (both classes at once: the same two normal terms)
        beta_sigma ~ nico_normal(0.0, prior_beta_sigma_SD);
        //// 
        //// ---- Induced-dirichlet ** Prior ** model:
        ////
        for (c in 1:2) {
            // array[n_index_tests] matrix[n_studies, n_thr_max] Ind_Dir_anchor; 
            for (t in 1:n_index_tests) {
                  int n_thr_t = n_thr[t];
                  int n_cat_t = n_thr_t + 1;
                  vector[n_thr_t] C_vec_t = get_test_values(
                                            C_vec[c], start_index, end_index, t);
                  // vector[n_thr_t] Ind_Dir_cumul_prob = (use_probit_link == 1) ? Phi(C_vec_t) : inv_logit(C_vec_t);
                  // vector[n_cat_t] Ind_Dir_ord_prob   = cumul_probs_to_ord_probs(Ind_Dir_cumul_prob);
                  // Ind_Dir_ord_prob ~ induced_dirichlet_given_C(
                  //                    C_vec_t, to_vector(prior_dirichlet_alpha[c][t, 1:n_cat_t]), use_probit_link);
                  //// (NicoSIMD's general node at the prior mean of the intercept: the anchor point)
                  target += nico_induced_dirichlet_prior_of_ordered_cutpoints_lupdf(
                                  C_vec_t | to_vector(prior_dirichlet_alpha[c][t, 1:n_cat_t]), use_probit_link,
                                  prior_beta_mu_mean[t][c, 1]);
            }
        }
        ////
        //// ---- Likelihood / Observational Model:
        ////
        // target += std_normal_lpdf(to_vector(beta_eta_z));
        target += nico_normal_lupdf(to_vector(beta_eta_z) | 0, 1);
        ////
        for (c in 1:2) {
          // target += std_normal_lpdf(to_vector(beta_delta_z_vec[c])); // part of between-study model, NOT prior
          //// (part of between-study model, NOT prior)
          target += nico_normal_lupdf(to_vector(beta_delta_z_vec[c]) | 0, 1);
        }
        ////
        for (t in 1:n_index_tests) {
          for (c in 1:2) {
               int n_thr_t = n_thr[t];
               vector[n_thr_t] C_raw_vec_test_t = get_test_values(C_raw_vec[c], start_index, end_index, t);
               // target += raw_C_to_C_log_det_J_lp(C_raw_vec_test_t, softplus);
               target += nico_log_det_jacobian_of_cutpoints_from_raw_differences(C_raw_vec_test_t, softplus);
          }
        }
        ////
        //// ---- Log-likelihood:
        ////
        // {
        //       for (t in 1:n_index_tests) {
        //             ////
        //             array[2] vector[n_thr[t]] C_vec_t;
        //             for (c in 1:2) {
        //                 vector[n_thr[t]] C_vec_t_given_c = get_test_values(
        //                                                    C_vec[c], start_index, end_index, t);
        //                 C_vec_t[c] = C_vec_t_given_c;
        //             }
        //             ////
        //             array[2] matrix[n_studies, n_thr[t]] latent_surv_t = map_latent_surv_to_fixed_hetero_C(
        //                                                                  C_vec_t, Xbeta[t], n_studies, 
        //                                                                  n_obs_cutpoints[t, ], cutpoint_index[t, ]);
        //             ////
        //             target += compute_NMA_log_lik_binomial_fact_lp(
        //                       latent_surv_t, use_probit_link, n_thr[t], n_studies,
        //                       x_2[t, ], n[t, ], N_total[t, ], n_obs_cutpoints[t, ], 
        //                       indicator_index_test_in_study[, t], 
        //                       K_fold_CV_indicator); 
        //       }
        // }
        {
              ////
              //// ---- every term of every test in ONE autodiff node: NicoSIMD's general
              ////      nico_binomial_counts_above_successive_cutpoints_of_a_cumulative_link_lpmf (the survival
              ////      probability Phi or inv_logit of each term's latent location, its ratio to the previous
              ////      reported cutpoint's, within [1e-10, 1 - 1e-10] as safe_binomial_lpmf, and the binomial log
              ////      probability), the latent location as map_latent_surv_to_fixed_hetero_C: location - C
              ////
              vector[2 * n_studies * n_index_tests] locations_of_all_tests;
              for (t in 1:n_index_tests) {
                    locations_of_all_tests[((t - 1) * 2 * n_studies + 1):(t * 2 * n_studies)] =
                          to_vector(Xbeta[t]);
              }
              vector[2 * n_total_C_per_class_if_fixed] cutpoints_of_both_classes = append_row(C_vec[1], C_vec[2]);
              vector[n_terms] latent_of_each_term =
                    locations_of_all_tests[index_of_the_location_of_each_term] -
                    cutpoints_of_both_classes[index_of_the_cutpoint_of_each_term];
              target += nico_binomial_counts_above_successive_cutpoints_of_a_cumulative_link_lupmf(
                              count_above_the_cutpoint_of_each_term |
                              count_above_the_previous_cutpoint_of_each_term, latent_of_each_term,
                              is_the_first_reported_cutpoint_of_each_term, use_probit_link, 1e-10);
        }
}

 
generated quantities {
          array[2] matrix[n_index_tests, n_thr_max] C_array;
          for (c in 1:2) {
            // array[n_index_tests] matrix[n_studies, n_thr_max] Ind_Dir_anchor; 
            for (t in 1:n_index_tests) {
                  int n_thr_t = n_thr[t]; 
                  vector[n_thr_t] C_vec_t = get_test_values(
                                            C_vec[c], start_index, end_index, t);
                  C_array[c][t, 1:n_thr_t] = to_row_vector(C_vec_t);
            }
          }
          ////
          //// ---- Compute between-study variance-covariance matrices:
          ////
          corr_matrix[2] beta_Omega      = multiply_lower_tri_self_transpose(beta_L_Omega);
          matrix[2, 2] beta_Sigma = multiply_lower_tri_self_transpose(beta_L_Sigma);
          ////
          //// ---- Calculate summary accuracy for each test - "baseline" covariate values:
          ////
          array[n_index_tests] real Xbeta_baseline_nd;
          array[n_index_tests] real Xbeta_baseline_d;
          ////
          for (t in 1:n_index_tests) {
              Xbeta_baseline_nd[t] = dot_product( baseline_case_nd[t][1:n_covariates_nd[t]],
                                                  to_vector(beta_mu[t][1, 1:n_covariates_nd[t]]));
              Xbeta_baseline_d[t]  = dot_product( baseline_case_d[t][1:n_covariates_d[t]], 
                                                  to_vector(beta_mu[t][2, 1:n_covariates_d[t]]));
          }
          ////
          //// ---- Calculate baseline Se/Sp for each test:
          ////
          array[n_index_tests] vector[n_thr_max] Fp_baseline;
          array[n_index_tests] vector[n_thr_max] Sp_baseline;
          array[n_index_tests] vector[n_thr_max] Se_baseline;
          ////
          for (t in 1:n_index_tests) {
              //// ---- Initialize with -1 for unused thresholds
              Fp_baseline[t] = rep_vector(-1.0, n_thr_max);
              Sp_baseline[t] = rep_vector(-1.0, n_thr_max);
              Se_baseline[t] = rep_vector(-1.0, n_thr_max);
              ////
              //// ---- Calculate for actual thresholds
              vector[n_thr[t]] C_vec_t_nd = get_test_values(C_vec[1], start_index, end_index, t);
              vector[n_thr[t]] C_vec_t_d  = get_test_values(C_vec[2], start_index, end_index, t);
              ////
              for (k in 1:n_thr[t]) {
                  Fp_baseline[t][k] = (use_probit_link == 1) ? 
                                      nico_Phi(-(C_vec_t_nd[k] - Xbeta_baseline_nd[t])) : 
                                      nico_inv_logit(-(C_vec_t_nd[k] - Xbeta_baseline_nd[t]));
                  Sp_baseline[t][k] = 1.0 - Fp_baseline[t][k];
                  Se_baseline[t][k] = (use_probit_link == 1) ? 
                                      nico_Phi(-(C_vec_t_d[k] - Xbeta_baseline_d[t])) : 
                                      nico_inv_logit(-(C_vec_t_d[k] - Xbeta_baseline_d[t]));
              }
          }
          ////
          //// ---- Calculate predictive accuracy for each test (NMA structure):
          ////
          array[n_index_tests] vector[n_thr_max] Fp_baseline_pred;
          array[n_index_tests] vector[n_thr_max] Sp_baseline_pred;
          array[n_index_tests] vector[n_thr_max] Se_baseline_pred;
          ////
          vector[2] beta_eta_pred      = multi_normal_cholesky_rng(rep_vector(0.0, 2), beta_L_Sigma);
          for (t in 1:n_index_tests) {
              ////
              //// ---- Initialize with -1 for unused thresholds
              Fp_baseline_pred[t] = rep_vector(-1.0, n_thr_max);
              Sp_baseline_pred[t] = rep_vector(-1.0, n_thr_max);
              Se_baseline_pred[t] = rep_vector(-1.0, n_thr_max);
              //// ---- Draw test-specific deviations
              real beta_delta_pred_nd = normal_rng(0.0, beta_tau[t, 1]);
              real beta_delta_pred_d  = normal_rng(0.0, beta_tau[t, 2]);
              //// ---- Combine to get predictive values
              real Xbeta_baseline_pred_nd = Xbeta_baseline_nd[t] + beta_eta_pred[1] + beta_delta_pred_nd;
              real Xbeta_baseline_pred_d  = Xbeta_baseline_d[t]  + beta_eta_pred[2] + beta_delta_pred_d;
              //// ---- Calculate predictive Se/Sp
              vector[n_thr[t]] C_vec_t_nd = get_test_values(C_vec[1], start_index, end_index, t);
              vector[n_thr[t]] C_vec_t_d  = get_test_values(C_vec[2], start_index, end_index, t);
              ////
              for (k in 1:n_thr[t]) { 
                  Fp_baseline_pred[t][k] = (use_probit_link == 1) ? 
                                           nico_Phi(-(C_vec_t_nd[k] - Xbeta_baseline_pred_nd)) : 
                                           nico_inv_logit(-(C_vec_t_nd[k] - Xbeta_baseline_pred_nd));
                  Sp_baseline_pred[t][k] = 1.0 - Fp_baseline_pred[t][k];
                  Se_baseline_pred[t][k] = (use_probit_link == 1) ? 
                                           nico_Phi(-(C_vec_t_d[k]  - Xbeta_baseline_pred_d))   : 
                                           nico_inv_logit(-(C_vec_t_d[k]  - Xbeta_baseline_pred_d));
              }
          }
          ////
          //// ---- Log-lik + study-specific accuracy computation (FLATTENED):
          ////
          vector[n_total_study_thresh] sp_flat = rep_vector(-1.0, n_total_study_thresh); // global
          vector[n_total_study_thresh] se_flat = rep_vector(-1.0, n_total_study_thresh); // global
          ////
          array[n_index_tests] vector[n_studies] deviance_nd;
          array[n_index_tests] vector[n_studies] deviance_d;
          array[n_index_tests] vector[n_studies] deviance;
          ////
          //// ---- Initialize arrays:
          ////
          for (t in 1:n_index_tests) {
              deviance_nd[t] = rep_vector(0.0, n_studies);
              deviance_d[t]  = rep_vector(0.0, n_studies);
              deviance[t]    = rep_vector(0.0, n_studies);
          }
          ////
          //// ---- Compute accuracy measures and log-likelihood for each test
          ////
          vector[n_studies] log_lik_study = rep_vector(0.0, n_studies);
          {
              int flat_idx = 0;
              for (t in 1:n_index_tests) {
                      array[2] vector[n_thr[t]] C_vec_t_dummy;
                      for (c in 1:2) {
                          vector[n_thr[t]] C_vec_t = get_test_values(C_vec[c], start_index, end_index, t);
                          C_vec_t_dummy[c] = C_vec_t;
                      }
                      ////
                      array[2] matrix[n_studies, n_thr[t]] latent_surv_t = nico_map_latent_surv_to_fixed_hetero_C(
                                                                           C_vec_t_dummy, Xbeta[t], n_studies, 
                                                                           n_obs_cutpoints[t, ], cutpoint_index[t, ]);
                      ////
                      array[n_studies] int dummy_CV_indicator;
                      for (s in 1:n_studies) {
                          dummy_CV_indicator[s] = 1;
                      }
                      array[3, 2] matrix[n_studies, n_thr[t]] outs = nico_compute_NMA_log_lik_binomial_fact_data(
                                                                     latent_surv_t, use_probit_link, n_thr[t], n_studies, 
                                                                     x_2[t, ], n[t, ], N_total[t, ], n_obs_cutpoints[t, ],
                                                                     indicator_index_test_in_study[, t], 
                                                                     dummy_CV_indicator);
                      ////
                      //// ---- Sum log-lik across thresholds and disease status within each study
                      ////
                      for (s in 1:n_studies) {
                        if (indicator_index_test_in_study[s, t] == 1) {
                          for (k in 1:n_obs_cutpoints[t, s]) {
                            log_lik_study[s] += outs[1][1][s, k];  // Non-diseased
                            log_lik_study[s] += outs[1][2][s, k];  // Diseased
                          }
                        }
                      }
                      ////
                      //// ---- Extract results:
                      ////
                      array[2] matrix[n_studies, n_thr[t]] cond_prob = outs[2];
                      array[2] matrix[n_studies, n_thr[t]] surv_prob = outs[3];
                      ////
                      for (s in 1:n_studies) {
                           if (indicator_index_test_in_study[s, t] == 1) {
                                for (k in 1:n_thr[t]) {
                                    flat_idx += 1;
                                    ////
                                    sp_flat[flat_idx] = 1.0 - surv_prob[1][s, k];
                                    se_flat[flat_idx] = surv_prob[2][s, k];
                                }
                          }
                      }
                      ////
                      //// ---- Model fit (deviance):
                      ////
                      if (n_thr[t] <= n_thr_max) {  // Safety check
                          array[4] matrix[n_studies, n_thr[t]] outs_model_fit = compute_deviance(
                                                                                cond_prob, n_thr[t], n_studies,
                                                                                x_2[t, ], n[t, ], n_obs_cutpoints[t, ]);
                          ////
                          matrix[n_studies, n_thr[t]] x_hat_nd = outs_model_fit[1];
                          matrix[n_studies, n_thr[t]] dev_nd   = outs_model_fit[2];
                          matrix[n_studies, n_thr[t]] x_hat_d  = outs_model_fit[3];
                          matrix[n_studies, n_thr[t]] dev_d    = outs_model_fit[4];
                          ////
                          for (s in 1:n_studies) {
                              if (indicator_index_test_in_study[s, t] == 1) {
                                  for (k in 1:n_obs_cutpoints[t, s]) {
                                      deviance_nd[t][s] += dev_nd[s, k];
                                      deviance_d[t][s]  += dev_d[s, k];
                                  }
                                  deviance[t][s] = deviance_nd[t][s] + deviance_d[t][s];
                              }
                          }
                      }
              }
          }
          ////
          //// ---- Between-study heterogeneity metrics for NMA: --------------------------------------------------------
          ////
          //// Total variance for each test (Nyaga et al. formulation)
          array[n_index_tests] vector[2] total_SD_beta;
          ////
          for (t in 1:n_index_tests) {
            for (c in 1:2) {
              total_SD_beta[t][c] = sqrt(square(beta_sigma[c]) + square(beta_tau[t, c]));
            }
          }
          ////
          //// Proportion of variance explained by test-specific component
          array[n_index_tests] vector[2] prop_test_SD_beta;
          ////
          for (t in 1:n_index_tests) {
              prop_test_SD_beta[t][1] = beta_tau[t, 1] / total_SD_beta[t][1];
              prop_test_SD_beta[t][2] = beta_tau[t, 2] / total_SD_beta[t][2];
          }
          array[n_index_tests] vector[2] total_SD_inc_C;
          for (t in 1:n_index_tests) {
            for (c in 1:2) {
               // real median_C_SD = median(estimated_SD_C[c][t, 1:n_thr[t]]);
               total_SD_inc_C[t][c] = total_SD_beta[t][c]; ////  sqrt(square(total_SD_beta[t][c]) + square(median_C_SD));
            }
          }
          ////
          matrix[n_index_tests, n_index_tests] rho12;
          array[2] matrix[n_index_tests, n_index_tests] rho;
          {
             array[3, 2] matrix[n_index_tests, n_index_tests] outs = compute_Nyaga_NMA_summaries(
                                                                     n_index_tests, beta_tau, beta_Sigma, beta_Omega);
             rho12 = outs[3, 1];
             rho   = outs[2, ];
          }
          ////
          //// ---- AUC: -------------------------------------------------------------------------------------------------
          ////
          vector[n_index_tests] AUC;
          ////
          for (t in 1:n_index_tests) {
            real missing_value_marker = -1.0;
            AUC[t] = compute_AUC(
                     Se_baseline[t], Sp_baseline[t], n_thr[t], missing_value_marker);
          }
          ////
          //// ---- Predictive AUC (accounting for uncertainty) ----
          ////
          vector[n_index_tests] AUC_pred;
          ////
          for (t in 1:n_index_tests) {
            real missing_value_marker = -1.0;
            AUC_pred[t] = compute_AUC(
                          Se_baseline_pred[t], Sp_baseline_pred[t], n_thr[t], missing_value_marker);
          }
          ////
          //// ---- Pairwise AUC differences ----
          ////
          //// ---- Flatten these too to avoid waste
          ////
          vector[n_test_pairs] AUC_diff;
          {
            array[n_test_pairs] int AUC_test1_idx;
            array[n_test_pairs] int AUC_test2_idx;
            ////
            int pair_idx = 1;
            for (t1 in 1:(n_index_tests-1)) {
                for (t2 in (t1+1):n_index_tests) {
                    AUC_test1_idx[pair_idx] = t1;
                    AUC_test2_idx[pair_idx] = t2;
                    AUC_diff[pair_idx] = AUC[t1] - AUC[t2];
                    pair_idx += 1;
                }
            }
          }

}





