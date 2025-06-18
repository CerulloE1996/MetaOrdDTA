

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
          array[n_index_tests, 2] matrix[n_studies, n_thr_max + 1] x;
          array[n_index_tests, 2] matrix[n_studies, n_thr_max + 1] cutpoint_index;
          ////
          //// ---- Covariates:
          ////
          array[n_index_tests] int n_covariates_nd;
          array[n_index_tests] int n_covariates_d;
          int n_covariates_max; // this is the max across ALL index tests
          array[n_index_tests] matrix[n_studies, max(n_covariates_nd)] X_nd;
          array[n_index_tests] matrix[n_studies, max(n_covariates_d)]  X_d;
          array[n_index_tests] vector[max(n_covariates_nd)] baseline_case_nd;  // must be user-inputted - e.g. could be {0, 1, 45.3} for 2 binary covariates and 1 cts one (e.g. age)
          array[n_index_tests] vector[max(n_covariates_d)]  baseline_case_d;   // must be user-inputted - e.g. could be {0, 1, 45.3} for 2 binary covariates and 1 cts one (e.g. age)
          ////
          //// ---- Priors for beta:
          ////
          array[n_index_tests] matrix[2, n_covariates_max] prior_beta_mu_mean;
          array[n_index_tests] matrix<lower=0.0>[2, n_covariates_max] prior_beta_mu_SD;
          ////
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
}

transformed data {
          int use_probit_link = 1;
          ////
          array[n_index_tests, 2] matrix[n_studies, n_thr_max] x_2;
          array[n_index_tests, 2] vector[n_studies] N_total;
          array[n_index_tests, 2] matrix[n_studies, n_thr_max] n;
          //// ---- Initialize with -1 (missing indicator):
          for (t in 1:n_index_tests) {
              for (c in 1:2) {
                  x_2[t, c]     = rep_matrix(-1, n_studies, n_thr_max);
                  n[t, c]       = rep_matrix(-1, n_studies, n_thr_max);
                  N_total[t, c] = rep_vector(-1, n_studies);
              }
          }
          ////
          for (t in 1:n_index_tests) {
                for (s in 1:n_studies) {
                          for (c in 1:2) {
                              N_total[t, c][s] = x[t, c][s, 1];
                              for (k in 1:n_thr[t]) {
                                 x_2[t, c][s, k] = x[t, c][s, k + 1];
                              }
                          }
                          for (c in 1:2) {
                               n[t, c][s, 1] = N_total[t, c][s];
                               for (k in 2:n_obs_cutpoints[t, s]) {
                                          n[t, c][s, k] = x_2[t, c][s, k - 1];
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
}


parameters {
          array[n_index_tests] matrix[2, n_covariates_max] beta_mu; 
          ////
          array[2] vector<lower=-15.0, upper=15.0>[n_total_C_per_class_if_fixed] C_raw_vec;  //// RAW LOG-DIFFS - Global cutpoints for each test (staggered array/matrix using "n_thr[t]" to index correctly)
          ////
          //// ---- "NMA" params:
          ////
          vector<lower=0.0>[2] beta_sigma;              //// Variance component - Between-study SD (Nyaga's σ)
          matrix<lower=0.0>[n_index_tests, 2] beta_tau; //// Variance component - Test-specific SD (Nyaga's τ) - delta_{s, c, t} ~ normal(0, tau_{c, t}).
          matrix[n_studies, 2] beta_eta_z;              //// Standard normal RVs for study-level effects - eta[s, 1:2] ~ multi_normal({0, 0}, Sigma).
          array[n_index_tests] matrix[n_studies, 2] beta_delta_z; //// Standard normal RVs for test-specific effects
          ////
          real<lower=beta_corr_lb, upper=beta_corr_ub> beta_corr;  //// between-study corr (possibly restricted) 
}

 
transformed parameters {
          ////
          //// ---- Construct simple 2x2 (bivariate) between-study corr matrices:
          ////
          cholesky_factor_corr[2] beta_L_Omega = make_bivariate_L_Omega(beta_corr);
          cholesky_factor_cov[2]  beta_L_Sigma = diag_pre_multiply(beta_sigma, beta_L_Omega);
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
          //// ---- "NMA" params:
          ////
          matrix[n_studies, 2] beta_eta = rep_matrix(0.0, n_studies, 2);
          array[n_index_tests] matrix[n_studies, 2] beta_delta = init_array_of_matrices(n_studies, 2, n_index_tests, 0.0);
          ////
          //// ---- Compute Study-level random effects (eta in Nyaga notation):
          ////
          for (s in 1:n_studies) {
              beta_eta[s, 1:2] = to_row_vector( beta_L_Sigma * to_vector(beta_eta_z[s, 1:2]) );  //// beta_eta[s, 1:2] ~ normal({0.0, 0.0}, Sigma);
          }
          for (c in 1:2) {
              //// Compute test-specific deviations ("delta" in Nyaga notation):
              for (t in 1:n_index_tests) { 
                 beta_delta[t][, c] = beta_tau[t, c] * beta_delta_z[t][, c]; //// beta_delta[t][s, c] ~ normal(0.0, beta_tau[t, c]);
              }
          }
          ////
          //// ---- Study-level random effects (after Cholesky decomposition):
          ////
          array[n_index_tests] matrix[n_studies, 2] beta_random;
          for (t in 1:n_index_tests) {
                for (s in 1:n_studies) {
                     beta_random[t][s, 1:2] =  beta_delta[t][s, 1:2] + beta_eta[s, 1:2];
                }
          }
          ////
          //// ---- Linear predictors for each disease statu + apply covariates for non-diseased and diseased groups:
          ////
          array[n_index_tests] matrix[n_studies, 2] Xbeta; // NOT local
          ////
          for (t in 1:n_index_tests) {
              Xbeta[t][, 1]       = X_nd[t][1:n_studies, 1:n_covariates_nd[t]]  * to_vector(beta_mu[t][1, 1:n_covariates_nd[t]]) + beta_random[t][, 1];
              Xbeta[t][, 2]       = X_d[t][1:n_studies,  1:n_covariates_d[t]]   * to_vector(beta_mu[t][2, 1:n_covariates_d[t]])  + beta_random[t][, 2];
          }
}
 
 
model { 
        //// 
        //// ---- Priors:
        ////
        for (t in 1:n_index_tests) {
            beta_mu[t][1, 1:n_covariates_nd[t]] ~ normal(
                                                  prior_beta_mu_mean[t][1, 1:n_covariates_nd[t]], prior_beta_mu_SD[t][1, 1:n_covariates_nd[t]]);  
            beta_mu[t][2, 1:n_covariates_d[t]]  ~ normal(
                                                  prior_beta_mu_mean[t][2, 1:n_covariates_d[t]],  prior_beta_mu_SD[t][2, 1:n_covariates_d[t]]);
        }
        ////
        to_vector(beta_tau) ~ normal(0.0, to_vector(prior_beta_tau_SD));  //// delta_{s, i, t} ~ normal(0, tau_{i, t}):
        beta_L_Omega ~ lkj_corr_cholesky(prior_beta_corr_LKJ);
        ////
        for (c in 1:2) {
          beta_sigma[c] ~ normal(0.0, prior_beta_sigma_SD[c]);      //// eta[s, i] ~ normal(0, sigma[i]):
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
        //// (part of between-test / between-study model, NOT prior) - eta[s, i] ~ normal(0, sigma[i]):
        target += std_normal_lpdf(to_vector(beta_eta_z)); 
        //// (part of between-test / between-study model, NOT prior)  -   delta_{s, i, t} ~ normal(0, tau_{i, t}):
        for (t in 1:n_index_tests) {
          target += std_normal_lpdf(to_vector(beta_delta_z[t])); 
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
                    array[2] matrix[n_studies, n_thr[t]] latent_surv_t = map_latent_surv_prob_to_fixed_hetero_C(
                                                                         C_vec_t, Xbeta[t], 1.0, n_studies, n_obs_cutpoints[t, ], 
                                                                         cutpoint_index[t, ]);
                    ////
                    target += compute_log_lik_binomial_fact_lp(
                              latent_surv_t, use_probit_link, n_thr[t], x_2[t, ], n[t, ], 
                              N_total[t, ], n_obs_cutpoints[t, ], indicator_index_test_in_study[, t]); 
              }
        }
}

 
generated quantities {
          ////
          //// ---- Compute between-study variance-covariance matrices:
          ////
          corr_matrix[2] beta_Omega      = multiply_lower_tri_self_transpose(beta_L_Omega);
          matrix[2, 2]   beta_Sigma      = multiply_lower_tri_self_transpose(beta_L_Sigma);
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
              
              //// ---- Calculate for actual thresholds
              vector[n_thr[t]] C_vec_t_nd = get_test_values(C_vec[1], start_index, end_index, t);
              vector[n_thr[t]] C_vec_t_d  = get_test_values(C_vec[2], start_index, end_index, t);
              
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
          //// ---- Draw shared study-level effects
          vector[2] beta_eta_pred      = multi_normal_cholesky_rng(rep_vector(0.0, 2), beta_L_Sigma);
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
          //// ---- Log-lik + study-specific accuracy computation:
          ////
          array[n_index_tests, 2] matrix[n_studies, n_thr_max] log_lik;
          array[n_index_tests] matrix[n_studies, n_thr_max] fp;
          array[n_index_tests] matrix[n_studies, n_thr_max] sp;
          array[n_index_tests] matrix[n_studies, n_thr_max] se;
          ////
          array[n_index_tests] vector[n_studies] deviance_nd;
          array[n_index_tests] vector[n_studies] deviance_d;
          array[n_index_tests] vector[n_studies] deviance;
          ////
          //// ---- Initialize arrays:
          ////
          for (t in 1:n_index_tests) {
              log_lik[t] = init_array_of_matrices(n_studies, n_thr_max, 2, 0.0);
              fp[t] = rep_matrix(-1.0, n_studies, n_thr_max);
              sp[t] = rep_matrix(-1.0, n_studies, n_thr_max);
              se[t] = rep_matrix(-1.0, n_studies, n_thr_max);
              deviance_nd[t] = rep_vector(0.0, n_studies);
              deviance_d[t]  = rep_vector(0.0, n_studies);
              deviance[t]    = rep_vector(0.0, n_studies);
          }
          ////
          //// ---- Compute accuracy measures and log-likelihood for each test
          ////
          {
              for (t in 1:n_index_tests) {
                  array[2] vector[n_thr[t]] C_vec_t_dummy;
                  for (c in 1:2) {
                      vector[n_thr[t]] C_vec_t = get_test_values(C_vec[c], start_index, end_index, t);
                      C_vec_t_dummy[c] = C_vec_t;
                  }
                  ////
                  array[2] matrix[n_studies, n_thr[t]] latent_surv_t = map_latent_surv_prob_to_fixed_hetero_C(
                                                                       C_vec_t_dummy, Xbeta[t], 1.0, n_studies, n_obs_cutpoints[t, ], 
                                                                       cutpoint_index[t, ]);
                  ////
                  array[3, 2] matrix[n_studies, n_thr[t]] outs = compute_log_lik_binomial_fact_data(
                                                                 latent_surv_t, use_probit_link, n_thr[t], x_2[t, ], n[t, ], 
                                                                 N_total[t, ], n_obs_cutpoints[t, ]);
                  ////
                  //// ---- Extract results:
                  ////
                  for (c in 1:2) {
                      for (s in 1:n_studies) {
                          for (k in 1:n_thr[t]) {
                              log_lik[t, c][s, k] = outs[1][c][s, k];
                          }
                      }
                  }
                  ////
                  array[2] matrix[n_studies, n_thr[t]] cond_prob = outs[2];
                  array[2] matrix[n_studies, n_thr[t]] surv_prob = outs[3];
                  ////
                  //// ---- Store accuracy measures:
                  ////
                  for (s in 1:n_studies) {
                      for (k in 1:n_thr[t]) {
                          fp[t][s, k] = surv_prob[1][s, k];
                          sp[t][s, k] = 1.0 - fp[t][s, k];
                          se[t][s, k] = surv_prob[2][s, k];
                      }
                  }
                  ////
                  //// ---- Model fit (deviance):
                  ////
                  if (n_thr[t] <= n_thr_max) {  // Safety check
                      array[4] matrix[n_studies, n_thr[t]] outs_model_fit = compute_deviance(
                                                                            cond_prob, n_thr[t], x_2[t, ], n[t, ], n_obs_cutpoints[t, ]);
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
          //// ---- NMA-specific outputs: Comparative measures between tests:
          ////
          array[n_index_tests, n_index_tests] matrix[n_thr_max, n_thr_max] diff_Se;
          array[n_index_tests, n_index_tests] matrix[n_thr_max, n_thr_max] diff_Sp;
          array[n_index_tests, n_index_tests] matrix[n_thr_max, n_thr_max] ratio_Se;
          array[n_index_tests, n_index_tests] matrix[n_thr_max, n_thr_max] ratio_Sp;
          ////
          //// ---- Initialize comparison matrices:
          ////
          for (t1 in 1:n_index_tests) {
              for (t2 in 1:n_index_tests) {
                  diff_Se[t1, t2]  = rep_matrix(-999.0, n_thr_max, n_thr_max);
                  diff_Sp[t1, t2]  = rep_matrix(-999.0, n_thr_max, n_thr_max);
                  ratio_Se[t1, t2] = rep_matrix(-999.0, n_thr_max, n_thr_max);
                  ratio_Sp[t1, t2] = rep_matrix(-999.0, n_thr_max, n_thr_max);
              }
          }
          ////
          //// ---- Compute pairwise comparisons:
          ////
          for (t1 in 1:(n_index_tests-1)) {
              for (t2 in (t1+1):n_index_tests) {
                  for (k1 in 1:n_thr[t1]) {
                      for (k2 in 1:n_thr[t2]) {
                          // Differences
                          diff_Se[t1, t2][k1, k2] = Se_baseline[t1][k1] - Se_baseline[t2][k2];
                          diff_Sp[t1, t2][k1, k2] = Sp_baseline[t1][k1] - Sp_baseline[t2][k2];
                          // Ratios
                          if (Se_baseline[t2][k2] > 0) {
                              ratio_Se[t1, t2][k1, k2] = Se_baseline[t1][k1] / Se_baseline[t2][k2];
                          }
                          if (Sp_baseline[t2][k2] > 0) {
                              ratio_Sp[t1, t2][k1, k2] = Sp_baseline[t1][k1] / Sp_baseline[t2][k2];
                          }
                      }
                  }
              }
          }
          ////
          //// ---- Diagnostic performance metrics for each test:
          ////
          array[n_index_tests] vector[n_thr_max] DOR;
          array[n_index_tests] vector[n_thr_max] LRpos;
          array[n_index_tests] vector[n_thr_max] LRneg;
          array[n_index_tests] vector[n_thr_max] Youden;
          ////
          for (t in 1:n_index_tests) {
              DOR[t]    = rep_vector(-1.0, n_thr_max);
              LRpos[t]  = rep_vector(-1.0, n_thr_max);
              LRneg[t]  = rep_vector(-1.0, n_thr_max);
              Youden[t] = rep_vector(-1.0, n_thr_max);
              
              for (k in 1:n_thr[t]) {
                  if (Se_baseline[t][k] > 0 && Sp_baseline[t][k] > 0) {
                      //// ---- Diagnostic odds ratio
                      DOR[t][k] = (Se_baseline[t][k] * Sp_baseline[t][k]) / ((1 - Se_baseline[t][k]) * (1 - Sp_baseline[t][k]));
                      //// ---- Likelihood ratios
                      LRpos[t][k] = Se_baseline[t][k] / (1 - Sp_baseline[t][k]);
                      LRneg[t][k] = (1 - Se_baseline[t][k]) / Sp_baseline[t][k];
                      //// ---- Youden index
                      Youden[t][k] = Se_baseline[t][k] + Sp_baseline[t][k] - 1;
                  }
              }
          }
          ////
          //// ---- Between-study heterogeneity metrics for NMA:
          ////
          //// Total variance for each test (Nyaga et al. formulation)
          array[n_index_tests] vector[2] total_var_beta;
          array[n_index_tests] vector[2] total_var_raw_scale;
          ////
          for (t in 1:n_index_tests) {
              total_var_beta[t][1] = square(beta_sigma[1]) + square(beta_tau[t, 1]);
              total_var_beta[t][2] = square(beta_sigma[2]) + square(beta_tau[t, 2]);
          }
          ////
          //// Proportion of variance explained by test-specific component
          array[n_index_tests] vector[2] prop_test_var_beta;
          array[n_index_tests] vector[2] prop_test_var_raw_scale;
          ////
          for (t in 1:n_index_tests) {
              prop_test_var_beta[t][1] = square(beta_tau[t, 1]) / total_var_beta[t][1];
              prop_test_var_beta[t][2] = square(beta_tau[t, 2]) / total_var_beta[t][2];
          }
}





