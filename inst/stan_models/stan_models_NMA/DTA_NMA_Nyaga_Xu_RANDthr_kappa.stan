

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
        #include "Stan_fns_Jacobian.stan"
        #include "Stan_fns_simplex.stan"
        #include "Stan_fns_NMA_bivariate.stan"
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
          array[n_studies, n_index_tests] int<lower=0, upper=1> indicator_index_test_in_study; // Binary ind. if test t is in study s
          ////
          //// ---- Data:
          ////
          array[n_index_tests, 2, n_studies, n_thr_max + 1] int x;
          array[n_index_tests, 2, n_studies, n_thr_max + 1] int cutpoint_index;;
          ////
          //// ---- Covariates:
          ////
          array[n_index_tests] int n_covariates_nd;
          array[n_index_tests] int n_covariates_d;
          array[n_index_tests] matrix[n_studies, max(n_covariates_nd)] X_nd;
          array[n_index_tests] matrix[n_studies, max(n_covariates_d)]  X_d;
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
          array[2] matrix[n_index_tests, n_cat_max] prior_dirichlet_phi;
          array[2] vector[n_index_tests] prior_kappa_mean;
          array[2] vector<lower=0.0>[n_index_tests] prior_kappa_SD;
          array[2] vector[n_index_tests] prior_kappa_df;
          real<lower=0.0> kappa_lb;
          real<lower=kappa_lb>  kappa_ub;
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
          ////
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
          ////
          //// ---- Calculate indices:
          ////
          array[n_index_tests] int start_index = calculate_start_indices(n_thr, n_index_tests);
          array[n_index_tests] int end_index   = calculate_end_indices(n_thr, n_index_tests, start_index);
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
          int n_total_C_if_random = 0;
          for (t in 1:n_index_tests) {
            n_total_C_if_random += n_thr[t]*n_studies_per_test[t];
          }
}


parameters {
          matrix[2, n_studies] beta_eta_z;
          array[2] vector[n_total_obs] beta_delta_z_vec; //// Standard normal RVs for test-specific effects
          ////
          array[n_index_tests] matrix[2, n_covariates_max] beta_mu;
          ////
          //// ---- "NMA" params:
          ////
          vector<lower=0.0>[2] beta_sigma;  //// Var component - Between-study SD (Nyaga's σ)
          matrix<lower=0.0>[n_index_tests, 2] beta_tau_raw;
          ////
          real<lower=beta_corr_lb, upper=beta_corr_ub> beta_corr;  //// between-study corr (possibly restricted)
          ////
          //// ---- Cutpoints and Induced-Dirichlet params:
          //// 
          array[2] vector<lower=-7.5, upper=2.5>[n_total_study_thresh] C_raw_vec;  //// RAW LOG-DIFFS
          array[2] vector<lower=-5.0, upper=5.0>[n_total_C_per_class_if_fixed] raw_dirichlet_phi_vec; // props/probs
          // array[2] vector<lower=0.0, upper=10.0>[n_index_tests] log_kappa; 
          array[2] vector<lower=kappa_lb, upper=kappa_ub>[n_index_tests] kappa;
} 


transformed parameters {
          array[2] vector[n_index_tests] log_kappa;
          for (c in 1:2) {
             log_kappa[c] = log(kappa[c]);
             // kappa[c] = exp(log_kappa[c]);
          }
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
          array[2] matrix[n_index_tests, n_cat_max] dirichlet_phi; // props/probs
          array[2] matrix[n_index_tests, n_cat_max] alpha;
          real log_det_J_simplex = 0.0;
          ////
          for (t in 1:n_index_tests) {
              for (c in 1:2) {
                    vector[n_thr[t]] raw_vec_simplex_t = get_test_values(
                                                         raw_dirichlet_phi_vec[c], start_index, end_index, t);
                    vector[n_cat[t] + 1] out = stickbreaking_logistic_simplex_constrain_inc_J(raw_vec_simplex_t);
                    ////
                    log_det_J_simplex += out[1];
                    ////
                    dirichlet_phi[c][t, 1:n_cat[t]] = to_row_vector(out[2:(n_cat[t] + 1)]);
                    ////
                    alpha[c][t, 1:n_cat[t]] = 0.01 + dirichlet_phi[c][t, 1:n_cat[t]] * kappa[c][t];
              }
          }
          ////
          //// ---- Construct study-specific cutpoints:
          ////
          array[2] vector[n_total_C_if_random] C_vec;
          real log_det_J_cutpoints = 0.0;
          for (c in 1:2) { 
            for (t in 1:n_index_tests) {
                    int n_thr_t = n_thr[t];
                    int study_idx = 0;
                    ////
                    for (s in 1:n_studies) {
                            if (indicator_index_test_in_study[s, t] == 1) {
                                study_idx += 1;
                                vector[n_thr_t] C_raw_vec_test_t = get_test_t_study_s_segment_values(
                                                                   C_raw_vec[c], t, study_idx, n_thr, n_studies_per_test);
                                vector[n_thr_t] C_vec_test_t = construct_C(C_raw_vec_test_t, softplus);
                                C_vec[c] = update_test_t_study_s_segment_values(
                                           C_vec[c], C_vec_test_t, t, study_idx, n_thr, n_studies_per_test);
                                log_det_J_cutpoints += raw_C_to_C_log_det_J_lp(C_raw_vec_test_t, softplus);
                            }
                    }
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
              beta_mu[t][1, 1:n_covariates_nd[t]] ~ normal(
                                                    prior_beta_mu_mean[t][1, 1:n_covariates_nd[t]], 
                                                    prior_beta_mu_SD[t][1, 1:n_covariates_nd[t]]);  
              beta_mu[t][2, 1:n_covariates_d[t]]  ~ normal(
                                                    prior_beta_mu_mean[t][2, 1:n_covariates_d[t]],  
                                                    prior_beta_mu_SD[t][2, 1:n_covariates_d[t]]);
          }
          ////
          to_vector(beta_tau_raw) ~ normal(0.0, to_vector(prior_beta_tau_SD));
          beta_L_Omega ~ lkj_corr_cholesky(prior_beta_corr_LKJ);
          ////
          for (c in 1:2) {
              beta_sigma[c] ~ normal(0.0, prior_beta_sigma_SD[c]);
          }
          ////
          //// ---- Priors for Induced-Dirichlet params for each test:
          ////
          for (t in 1:n_index_tests) {
              for (c in 1:2) {
                 log_kappa[c][t] ~ student_t(prior_kappa_df[c][t], prior_kappa_mean[c][t], prior_kappa_SD[c][t]);
                 target += -log_kappa[c][t];
                 ////
                 // kappa[c][t] ~ normal(prior_kappa_mean[c][t], prior_kappa_SD[c][t]);
                 ////
                 dirichlet_phi[c][t, 1:n_cat[t]] ~ dirichlet(prior_dirichlet_phi[c][t, 1:n_cat[t]]); // ID prior
                 ////
                 //// ---- Induced-dirichlet between study ** model ** (NOT a prior model here but part of the actual 
                 //// ---- likelihood since random-effect cutpoints!):
                 ////
                 vector[n_cat[t]] alpha_vec = to_vector(alpha[c][t, 1:n_cat[t]]);
                 int study_idx = 0;  // Track position among studies with test t
                 ////
                 for (s in 1:n_studies) {
                   if (indicator_index_test_in_study[s, t] == 1) {
                         study_idx += 1;  // Increment sparse index
                         vector[n_thr[t]] C_vec_t = get_test_t_study_s_segment_values(
                                                    C_vec[c], t, study_idx, n_thr, n_studies_per_test);
                         ////
                         vector[n_cat[t]] Ind_Dir_ord_prob = (use_probit_link == 1) ? 
                                                             cumul_probs_to_ord_probs(Phi(C_vec_t)) :
                                                             cumul_probs_to_ord_probs(inv_logit(C_vec_t));
                         ////
                         Ind_Dir_ord_prob ~ induced_dirichlet_given_C(
                                            C_vec_t, alpha_vec, use_probit_link);
                   }
                 }
              }
          }
          ////
          //// ---- Likelihood / Observational Model:
          ////
          target += std_normal_lpdf(to_vector(beta_eta_z));
          ////
          for (c in 1:2) {
            target += std_normal_lpdf(to_vector(beta_delta_z_vec[c])); // part of between-study model, NOT prior
          }
          // ////
          // for (t in 1:n_index_tests) {
          //    for (c in 1:2) {
          //          int study_idx = 0;
          //          for (s in 1:n_studies) {
          //             if (indicator_index_test_in_study[s, t] == 1) {
          //                  study_idx += 1;
          //                  vector[n_thr[t]] C_raw_vec_test_t = get_test_t_study_s_segment_values(
          //                                                      C_raw_vec[c], t, study_idx, n_thr, n_studies_per_test);
          //                  target += raw_C_to_C_log_det_J_lp(C_raw_vec_test_t, softplus);
          //             }
          //          }
          //    }
          // }
          ////
          target += log_det_J_simplex + log_det_J_cutpoints;
          // for (t in 1:n_index_tests) {
          //     for (c in 1:2) {
          //           vector[n_thr[t]] raw_vec_simplex_t = get_test_values(
          //                                                raw_dirichlet_phi_vec[c], start_index, end_index, t);
          //           target += stickbreaking_logistic_simplex_constrain_lp(raw_vec_simplex_t);
          //     }
          // }
          ////
          //// ---- Log-likelihood:
          ////
          for (t in 1:n_index_tests) {
                array[2] matrix[n_studies, n_thr[t]] C_mat;
                for (c in 1:2) {
                    int study_idx = 0;
                    for (s in 1:n_studies) {
                        if (indicator_index_test_in_study[s, t] == 1) {
                             study_idx += 1;
                             vector[n_thr[t]] C_vec_t_study_s = get_test_t_study_s_segment_values(
                                                                C_vec[c], t, study_idx, n_thr, n_studies_per_test);
                             C_mat[c][s, ] = to_row_vector(C_vec_t_study_s);
                        }
                    }
                }
                ////
                array[2] matrix[n_studies, n_thr[t]] latent_surv_t = map_latent_surv_to_random_hetero_C(
                                                                     C_mat, n_thr[t], Xbeta[t], n_studies, 
                                                                     n_obs_cutpoints[t, ], cutpoint_index[t, ]);
                ////
                target += compute_NMA_log_lik_binomial_fact_lp(
                          latent_surv_t, use_probit_link, n_thr[t], n_studies, 
                          x_2[t, ], n[t, ], N_total[t, ], n_obs_cutpoints[t, ], 
                          indicator_index_test_in_study[, t],
                          K_fold_CV_indicator);
          }
}


generated quantities {
          ////
          //// ---- Compute between-study variance-covariance matrices:
          ////
          corr_matrix[2] beta_Omega = multiply_lower_tri_self_transpose(beta_L_Omega);
          matrix[2, 2] beta_Sigma = multiply_lower_tri_self_transpose(beta_L_Sigma);
          ////
          //// ---- SD's for each category probability in a Dirichlet is √(α_k(α_0-α_k)/(α_0^2(α_0+1))) - where α_0 is the sum of all alphas:
          ////
          array[2] matrix[n_index_tests, n_cat_max] dir_cat_SDs_sigma;
          for (t in 1:n_index_tests) {
              for (c in 1:2) {
                  real alpha_0 = sum(alpha[c][t, 1:n_cat[t]]);
                  dir_cat_SDs_sigma[c][t, 1:n_cat[t]] = sqrt((alpha[c][t, 1:n_cat[t]] .* (alpha_0 - alpha[c][t, 1:n_cat[t]]) ) /
                                                                  (square(alpha_0) * (alpha_0 + 1.0)));
              }
          }
          ////
          //// ---- Calculate summary cutpoints:
          ////
          array[2] matrix[n_index_tests, n_thr_max] C_mu;
          for (t in 1:n_index_tests) {
              for (c in 1:2) {
                  for (k in 1:n_thr[t]) {
                      //// ---- Collect cutpoints across studies for this test and threshold
                      vector[n_studies_per_test[t]] C_vals;
                      int study_idx = 0;
                      for (s in 1:n_studies) {
                          if (indicator_index_test_in_study[s, t] == 1) {
                              study_idx += 1;
                              C_vals[study_idx] = get_test_t_study_s_segment_values(
                                                  C_vec[c], t, study_idx, n_thr, n_studies_per_test)[k];
                          }
                      }
                      vector[n_studies_per_test[t]] observed_cutpoints;
                      study_idx = 0;
                      int n_obs = 0;
                      for (s in 1:n_studies) {     
                          if (indicator_index_test_in_study[s, t] == 1) {
                              study_idx += 1;
                              // Search through the cutpoint_index for this study to see if cutpoint k was observed
                              for (j in 2:(n_obs_cutpoints[t, s] + 1)) {  // Start at 2 because index 1 is 0
                                  if (cutpoint_index[t, c, s, j] == k) {
                                      n_obs += 1;
                                      observed_cutpoints[n_obs] = C_vals[study_idx];
                                      break;  // Found it, no need to continue searching
                                  }
                              }
                          }
                      }
                      //// ---- median cutpoint across the n_studies cutpoints
                      if (n_obs > 0) {
                          C_mu[c][t, k] = median(observed_cutpoints[1:n_obs]); // solves "skewness" issue from non-linear transformations
                      }
                  }
              }
          }
          ////
          //// ---- Calculate summary accuracy for each test - "baseline" covariate values:
          ////
          array[n_index_tests] real Xbeta_baseline_nd;
          array[n_index_tests] real Xbeta_baseline_d;
          ////
          for (t in 1:n_index_tests) {
              Xbeta_baseline_nd[t] = dot_product(baseline_case_nd[t][1:n_covariates_nd[t]], to_vector(beta_mu[t][1, 1:n_covariates_nd[t]]));
              Xbeta_baseline_d[t]  = dot_product(baseline_case_d[t][1:n_covariates_d[t]], to_vector(beta_mu[t][2, 1:n_covariates_d[t]]));
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
              //// ---- Calculate for actual thresholds
              for (k in 1:n_thr[t]) {
                  Fp_baseline[t][k] = (use_probit_link == 1) ?
                                      safe_Phi(-(C_mu[1][t, k] - Xbeta_baseline_nd[t])) : 
                                      inv_logit(-(C_mu[1][t, k] - Xbeta_baseline_nd[t]));
                  Sp_baseline[t][k] = 1.0 - Fp_baseline[t][k];
                  Se_baseline[t][k] = (use_probit_link == 1) ?
                                      safe_Phi(-(C_mu[2][t, k] - Xbeta_baseline_d[t])) : 
                                      inv_logit(-(C_mu[2][t, k] - Xbeta_baseline_d[t]));
              }
          }
          ////
          //// ---- Calculate predictive accuracy:
          ////
          array[2] matrix[n_index_tests, n_thr_max] C_pred;
          array[2] matrix[n_index_tests, n_thr_max] estimated_SD_C;
          ////
          for (t in 1:n_index_tests) {
              for (c in 1:2) {
                for (k in 1:n_thr[t]) {
                    real sigma_C_k;
                    real SD_prob_scale = dir_cat_SDs_sigma[c][t, k];
                    if (use_probit_link == 1) {
                        real mean_C_cutpoint_scale = C_mu[c][t, k];
                        sigma_C_k = SD_approx_ID_ord_prob_to_C_probit(mean_C_cutpoint_scale, SD_prob_scale);
                    } else {
                        real mean_prob_scale = dirichlet_phi[c][t, k];
                        sigma_C_k = SD_approx_ID_ord_prob_to_C_logit(mean_prob_scale, SD_prob_scale);
                    }
                    estimated_SD_C[c][t, k] = sigma_C_k;
                    if (!is_nan(sigma_C_k)) {
                       C_pred[c][t, k] = C_mu[c][t, k] + normal_rng(0, sigma_C_k);
                    }
                }
              }
          }
          ////
          //// ---- Calculate predictive accuracy for each test (NMA structure):
          ////
          array[n_index_tests] vector[n_thr_max] Fp_baseline_pred;
          array[n_index_tests] vector[n_thr_max] Sp_baseline_pred;
          array[n_index_tests] vector[n_thr_max] Se_baseline_pred;
          //// ---- Draw shared study-level effects
          vector[2] beta_eta_pred = multi_normal_cholesky_rng(rep_vector(0.0, 2), beta_L_Sigma);
          ////
          for (t in 1:n_index_tests) {
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
              //// ---- Calculate predictive Se/Sp using predictive cutpoints
              for (k in 1:n_thr[t]) {
                  Fp_baseline_pred[t][k] = (use_probit_link == 1) ?
                                           safe_Phi(-(C_pred[1][t, k] - Xbeta_baseline_pred_nd)) :
                                           inv_logit(-(C_pred[1][t, k] - Xbeta_baseline_pred_nd));
                  Sp_baseline_pred[t][k] = 1.0 - Fp_baseline_pred[t][k];
                  Se_baseline_pred[t][k] = (use_probit_link == 1) ?
                                           safe_Phi(-(C_pred[2][t, k] - Xbeta_baseline_pred_d)) :
                                           inv_logit(-(C_pred[2][t, k] - Xbeta_baseline_pred_d));
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
                    array[2] matrix[n_studies, n_thr_max] C_mat;
                    ////
                    for (c in 1:2) {
                        for (s in 1:n_studies) {
                            C_mat[c][s, ] = rep_row_vector(-999.0, n_thr_max);  // Or some other missing indicator
                        }
                    }
                    for (c in 1:2) {
                        int study_idx = 0;
                        for (s in 1:n_studies) {
                            if (indicator_index_test_in_study[s, t] == 1) {
                                 study_idx += 1;
                                 vector[n_thr[t]] C_vec_t_study_s = get_test_t_study_s_segment_values(
                                                                    C_vec[c], t, study_idx, n_thr, n_studies_per_test);
                                 C_mat[c][s, 1:n_thr[t]] = to_row_vector(C_vec_t_study_s);
                            }
                        }
                    }
                    ////
                    array[2] matrix[n_studies, n_thr[t]] latent_surv_t = map_latent_surv_to_random_hetero_C(
                                                                         C_mat, n_thr[t], Xbeta[t], n_studies, 
                                                                         n_obs_cutpoints[t, ], cutpoint_index[t, ]);
                    ////
                    array[n_studies] int dummy_CV_indicator;
                    for (s in 1:n_studies) {
                        dummy_CV_indicator[s] = 1;
                    }
                    array[3, 2] matrix[n_studies, n_thr[t]] outs = compute_NMA_log_lik_binomial_fact_data(
                                                                   latent_surv_t, use_probit_link, n_thr[t], n_studies,
                                                                   x_2[t, ], n[t, ], N_total[t, ], n_obs_cutpoints[t, ],
                                                                   indicator_index_test_in_study[, t],
                                                                   dummy_CV_indicator);
                    ////
                    //// ---- Sum log-lik across thresholds and disease status within each study
                    ////
                    for (s in 1:n_studies) {
                      if (indicator_index_test_in_study[s, t] == 1) {
                        // Sum across all thresholds for this study and test
                        for (k in 1:n_obs_cutpoints[t, s]) {
                          // Both disease states
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
          // for (t in 1:n_index_tests) {
          //   for (c in 1:2) {
          //      real median_C_SD = median(estimated_SD_C[c][t, 1:n_thr[t]]);
          //      total_SD_inc_C[t][c] = sqrt(square(total_SD_beta[t][c]) + square(median_C_SD));
          //   }
          // }
          for (t in 1:n_index_tests) {
              for (c in 1:2) {
                    //// Create filtered vector without NaNs
                    int n_valid = 0;
                    vector[n_thr[t]] filtered_SD_C;
                    ////
                    //// Count valid (non-NaN) values
                    for (k in 1:n_thr[t]) {
                      if (!is_nan(estimated_SD_C[c][t, k])) {
                        n_valid += 1;
                      }
                    }
                    //// Fill filtered vector with valid values
                    if (n_valid > 0) {
                      vector[n_valid] valid_SD_C;
                      int idx = 1;
                      for (k in 1:n_thr[t]) {
                        if (!is_nan(estimated_SD_C[c][t, k])) {
                          valid_SD_C[idx] = estimated_SD_C[c][t, k];
                          idx += 1;
                        }
                      }
                      real median_C_SD = median(valid_SD_C);
                      total_SD_inc_C[t][c] = sqrt(square(total_SD_beta[t][c]) + square(median_C_SD));
                    } else {
                      //// Handle case where all values are NaN
                      total_SD_inc_C[t][c] = total_SD_beta[t][c]; // or set to NaN or some default
                    }
              }
          }
          ////
          matrix[n_index_tests, n_index_tests] rho12;
          array[2] matrix[n_index_tests, n_index_tests] rho;
          {
             array[3, 2] matrix[n_index_tests, n_index_tests] outs = compute_Nyaga_NMA_summaries(
                                                                     n_index_tests, beta_tau, beta_Sigma, beta_Omega);
             rho12 = outs[3, 1];
             rho = outs[2, ];
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






