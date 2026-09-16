

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
          array[n_index_tests, 2, n_studies, n_thr_max + 1] int x;
          array[n_index_tests, 2, n_studies, n_thr_max + 1] int cutpoint_index;
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
          ////
          int use_probit_link;
}

parameters {
    // Empty - parameters come from fitted model
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
          //// ---- Log-lik + study-specific accuracy computation (FLATTENED):
          ////
          matrix[2, n_total_study_thresh] log_lik_flat = rep_matrix(0.0, 2, n_total_study_thresh); // global
          ////
          vector[n_total_study_thresh] fp_flat = rep_vector(-1.0, n_total_study_thresh); // global
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
          {
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
                      array[3, 2] matrix[n_studies, n_thr[t]] outs = compute_NMA_log_lik_binomial_fact_data(
                                                                     latent_surv_t, use_probit_link, n_thr[t], n_studies, 
                                                                     x_2[t, ], n[t, ], N_total[t, ], n_obs_cutpoints[t, ],
                                                                     indicator_index_test_in_study[, t]);
                      ////
                      //// ---- Extract results:
                      ////
                      array[2] matrix[n_studies, n_thr[t]] cond_prob = outs[2];
                      array[2] matrix[n_studies, n_thr[t]] surv_prob = outs[3];
                      ////
                      int flat_idx = 1;
                          for (s in 1:n_studies) {
                              for (k in 1:n_thr[t]) {
                                  ////
                                  log_lik_flat[1, flat_idx] = outs[1][1][s, k];
                                  log_lik_flat[2, flat_idx] = outs[1][2][s, k];
                                  ////
                                  fp_flat[flat_idx] = surv_prob[1][s, k];
                                  sp_flat[flat_idx] = 1.0 - fp_flat[flat_idx];
                                  se_flat[flat_idx] = surv_prob[2][s, k];
                                  ////
                                  flat_idx += 1;
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
          //// ---- NMA-specific outputs: Comparative measures between tests:
          ////
          vector[n_total_comparisons] diff_Se_flat;
          vector[n_total_comparisons] diff_Sp_flat;
          vector[n_total_comparisons] ratio_Se_flat;
          vector[n_total_comparisons] ratio_Sp_flat;
          ////
          {;
              ////
              //// ---- Fill the flattened arrays:
              ////
              int comp_idx = 1;
              for (t1 in 1:(n_index_tests-1)) {
                  for (t2 in (t1+1):n_index_tests) {
                      for (k1 in 1:n_thr[t1]) {
                          for (k2 in 1:n_thr[t2]) {
                              ////
                              diff_Se_flat[comp_idx] = Se_baseline[t1][k1] - Se_baseline[t2][k2];
                              diff_Sp_flat[comp_idx] = Sp_baseline[t1][k1] - Sp_baseline[t2][k2];
                              ////
                              if (Se_baseline[t2][k2] > 0) {
                                  ratio_Se_flat[comp_idx] = Se_baseline[t1][k1] / Se_baseline[t2][k2];
                              } else {
                                  ratio_Se_flat[comp_idx] = -999.0; // Missing indicator
                              }
                              ////
                              if (Sp_baseline[t2][k2] > 0) {
                                  ratio_Sp_flat[comp_idx] = Sp_baseline[t1][k1] / Sp_baseline[t2][k2];
                              } else {
                                  ratio_Sp_flat[comp_idx] = -999.0;
                              }
                              ////
                              comp_idx += 1;
                          }
                      }
                  }
              }
          }
          ////
          //// ---- Diagnostic performance metrics for each test (FLATTENED):
          ////
          vector[sum(n_thr)] DOR_flat;
          vector[sum(n_thr)] LRpos_flat;
          vector[sum(n_thr)] LRneg_flat;
          vector[sum(n_thr)] Youden_flat;
          {
              int perf_idx = 1;
              for (t in 1:n_index_tests) {
                  for (k in 1:n_thr[t]) {
                      ////
                      if (Se_baseline[t][k] > 0 && Sp_baseline[t][k] > 0) {
                          DOR_flat[perf_idx] = (Se_baseline[t][k] * Sp_baseline[t][k]) / 
                                               ((1 - Se_baseline[t][k]) * (1 - Sp_baseline[t][k]));
                          LRpos_flat[perf_idx] = Se_baseline[t][k] / (1 - Sp_baseline[t][k]);
                          LRneg_flat[perf_idx] = (1 - Se_baseline[t][k]) / Sp_baseline[t][k];
                          Youden_flat[perf_idx] = Se_baseline[t][k] + Sp_baseline[t][k] - 1;
                      } else {
                          DOR_flat[perf_idx] = -999.0;
                          LRpos_flat[perf_idx] = -999.0;
                          LRneg_flat[perf_idx] = -999.0;
                          Youden_flat[perf_idx] = -999.0;
                      }
                      ////
                      perf_idx += 1;
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
          ////
          //// ---- AUC computation for each test ----
          ////
          vector[n_index_tests] AUC;
          ////
          for (t in 1:n_index_tests) {
              //// Initialize AUC
              AUC[t] = 0.0;
              ////
              //// We need to add points (0,0) and (1,1) to complete the ROC curve
              //// Create extended arrays with these points
              ////
              int n_points = n_thr[t] + 2;
              vector[n_points] FPR; // False Positive Rate = 1 - Specificity
              vector[n_points] TPR; // True Positive Rate = Sensitivity
              ////
              //// Add (0,0) point
              ////
              FPR[1] = 0.0;
              TPR[1] = 0.0;
              ////
              //// Add actual threshold points
              ////
              for (k in 1:n_thr[t]) {
                  FPR[k+1] = 1.0 - Sp_baseline[t][k];  // FPR = 1 - Specificity
                  TPR[k+1] = Se_baseline[t][k];         // TPR = Sensitivity
              }
              ////
              //// Add (1,1) point
              ////
              FPR[n_points] = 1.0;
              TPR[n_points] = 1.0;
              ////
              //// Sort by FPR (in case thresholds aren't ordered by FPR)
              //// Simple bubble sort since n_points is typically small
              ////
              for (i in 1:(n_points-1)) {
                  for (j in 1:(n_points-i)) {
                      if (FPR[j] > FPR[j+1]) {
                          // Swap FPR
                          real temp_fpr = FPR[j];
                          FPR[j] = FPR[j+1];
                          FPR[j+1] = temp_fpr;
                          // Swap TPR
                          real temp_tpr = TPR[j];
                          TPR[j] = TPR[j+1];
                          TPR[j+1] = temp_tpr;
                      }
                  }
              }
              ////
              //// Calculate AUC using trapezoidal rule
              ////
              for (i in 1:(n_points-1)) {
                  AUC[t] += 0.5 * (TPR[i] + TPR[i+1]) * (FPR[i+1] - FPR[i]);
              }
          }
          ////
          //// ---- Predictive AUC (accounting for uncertainty) ----
          ////
          vector[n_index_tests] AUC_pred;
          ////
          for (t in 1:n_index_tests) {
              AUC_pred[t] = 0.0;
              ////
              int n_points = n_thr[t] + 2;
              vector[n_points] FPR_pred;
              vector[n_points] TPR_pred;
              ////
              FPR_pred[1] = 0.0;
              TPR_pred[1] = 0.0;
              ////
              for (k in 1:n_thr[t]) {
                  FPR_pred[k+1] = 1.0 - Sp_baseline_pred[t][k];
                  TPR_pred[k+1] = Se_baseline_pred[t][k];
              }
              ////
              FPR_pred[n_points] = 1.0;
              TPR_pred[n_points] = 1.0;
              ////
              //// Sort by FPR
              ////
              for (i in 1:(n_points-1)) {
                  for (j in 1:(n_points-i)) {
                      if (FPR_pred[j] > FPR_pred[j+1]) {
                          real temp_fpr = FPR_pred[j];
                          FPR_pred[j] = FPR_pred[j+1];
                          FPR_pred[j+1] = temp_fpr;
                          real temp_tpr = TPR_pred[j];
                          TPR_pred[j] = TPR_pred[j+1];
                          TPR_pred[j+1] = temp_tpr;
                      }
                  }
              }
              ////
              //// Calculate AUC
              ////
              for (i in 1:(n_points-1)) {
                  AUC_pred[t] += 0.5 * (TPR_pred[i] + TPR_pred[i+1]) * (FPR_pred[i+1] - FPR_pred[i]);
              }
          }
          ////
          //// ---- Pairwise AUC differences ----
          ////
          //// ---- Flatten these too to avoid waste
          ////
          vector[n_test_pairs] AUC_diff;
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






