

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
          array[n_index_tests] vector<lower=0.0>[n_thr_max + 1] prior_dirichlet_alpha;
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
}


parameters {
          matrix[2, n_studies] beta_eta_z;
          array[2] vector[n_total_obs] beta_delta_z_vec; //// Standard normal RVs for test-specific effects
          ////
          array[n_index_tests] matrix[2, n_covariates_max] beta_mu; 
          ////
          array[2] vector<lower=-7.5, upper=2.5>[n_total_C_per_class_if_fixed] C_raw_vec;  //// RAW LOG-DIFFS - Global cutpoints for each test (staggered array/matrix using "n_thr[t]" to index correctly)
          ////
          //// ---- "NMA" params:
          ////
          vector<lower=0.0>[2] beta_sigma;  //// Var component - Between-study SD (Nyaga's σ)
          ////
          matrix<lower=0.0>[n_index_tests, 2] beta_tau_raw;
          ////
          real<lower=beta_corr_lb, upper=beta_corr_ub> beta_corr;  //// between-study corr (possibly restricted) 
}

 
transformed parameters {
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
                    vector[n_thr_t] C_vec_test_t =  construct_C(C_raw_vec_test_t, softplus);
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
            beta_mu[t][1, 1:n_covariates_nd[t]] ~ normal(
                                                  prior_beta_mu_mean[t][1, 1:n_covariates_nd[t]], 
                                                  prior_beta_mu_SD[t][1, 1:n_covariates_nd[t]]);  
            beta_mu[t][2, 1:n_covariates_d[t]]  ~ normal(
                                                  prior_beta_mu_mean[t][2, 1:n_covariates_d[t]],  
                                                  prior_beta_mu_SD[t][2, 1:n_covariates_d[t]]);
        }
        ////
        to_vector(beta_tau_raw) ~ normal(0.0, to_vector(prior_beta_tau_SD));  //// delta_{s, i, t} ~ normal(0, tau_{i, t}):
        beta_L_Omega ~ lkj_corr_cholesky(prior_beta_corr_LKJ);
        ////
        for (c in 1:2) {
            beta_sigma[c] ~ normal(0.0, prior_beta_sigma_SD[c]);
        }
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
                  vector[n_thr_t] Ind_Dir_cumul_prob = (use_probit_link == 1) ? Phi(C_vec_t) : inv_logit(C_vec_t);
                  vector[n_cat_t] Ind_Dir_ord_prob   = cumul_probs_to_ord_probs(Ind_Dir_cumul_prob);
                  Ind_Dir_ord_prob ~ induced_dirichlet_given_C(
                                     C_vec_t, prior_dirichlet_alpha[t][1:n_cat_t], use_probit_link);
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
        ////
        for (t in 1:n_index_tests) {
          for (c in 1:2) {
               int n_thr_t = n_thr[t];
               vector[n_thr_t] C_raw_vec_test_t = get_test_values(C_raw_vec[c], start_index, end_index, t);
               target += raw_C_to_C_log_det_J_lp(C_raw_vec_test_t, softplus);
          }
        }
        ////
        //// ---- Log-likelihood:
        ////
        {
              for (t in 1:n_index_tests) {
                    ////
                    array[2] vector[n_thr[t]] C_vec_t;
                    for (c in 1:2) {
                        vector[n_thr[t]] C_vec_t_given_c = get_test_values(
                                                           C_vec[c], start_index, end_index, t);
                        C_vec_t[c] = C_vec_t_given_c;
                    }
                    ////
                    array[2] matrix[n_studies, n_thr[t]] latent_surv_t = map_latent_surv_to_fixed_hetero_C(
                                                                         C_vec_t, Xbeta[t], n_studies, 
                                                                         n_obs_cutpoints[t, ], cutpoint_index[t, ]);
                    ////
                    target += compute_NMA_log_lik_binomial_fact_lp(
                              latent_surv_t, use_probit_link, n_thr[t], n_studies,
                              x_2[t, ], n[t, ], N_total[t, ], n_obs_cutpoints[t, ], 
                              indicator_index_test_in_study[, t], 
                              K_fold_CV_indicator); 
              }
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
                                      Phi(-(C_vec_t_nd[k] - Xbeta_baseline_nd[t])) : 
                                      inv_logit(-(C_vec_t_nd[k] - Xbeta_baseline_nd[t]));
                  Sp_baseline[t][k] = 1.0 - Fp_baseline[t][k];
                  Se_baseline[t][k] = (use_probit_link == 1) ? 
                                      Phi(-(C_vec_t_d[k] - Xbeta_baseline_d[t])) : 
                                      inv_logit(-(C_vec_t_d[k] - Xbeta_baseline_d[t]));
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
                                           Phi(-(C_vec_t_nd[k] - Xbeta_baseline_pred_nd)) : 
                                           inv_logit(-(C_vec_t_nd[k] - Xbeta_baseline_pred_nd));
                  Sp_baseline_pred[t][k] = 1.0 - Fp_baseline_pred[t][k];
                  Se_baseline_pred[t][k] = (use_probit_link == 1) ? 
                                           Phi(-(C_vec_t_d[k]  - Xbeta_baseline_pred_d))   : 
                                           inv_logit(-(C_vec_t_d[k]  - Xbeta_baseline_pred_d));
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
                      array[2] matrix[n_studies, n_thr[t]] latent_surv_t = map_latent_surv_to_fixed_hetero_C(
                                                                           C_vec_t_dummy, Xbeta[t], n_studies, 
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





