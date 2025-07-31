

functions {
        ////
        //// ---- Include files to compile the necessary custom (user-defined) Stan functions:
        ////
        #include "Stan_fns_basic.stan"
        #include "Stan_fns_Box_Cox.stan"
        #include "Stan_fns_corr.stan"
        #include "Stan_fns_ordinal.stan"
        #include "Stan_fns_log_lik.stan"
        #include "Stan_fns_ragged.stan"
        #include "Stan_fns_NMA.stan"
        #include "Stan_fns_model_fit.stan"
        #include "Stan_fns_NMA_bivariate.stan"
}

 
data {
          int<lower=1> n_studies;
          int<lower=1> n_index_tests;
          ////
          array[n_index_tests] int<lower=1> n_thr;
          int<lower=1> n_thr_max; // corresponding to the test with the most thresholds/cutpoints. 
          int<lower=1> n_cat_max; // corresponding to the test with the most thresholds/cutpoints. 
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
          array[n_index_tests] vector[max(n_covariates_nd)] baseline_case_nd;  // must be user-inputted - e.g. could be {0, 1, 45.3} for 2 binary covariates and 1 cts one (e.g. age)
          array[n_index_tests] vector[max(n_covariates_d)]  baseline_case_d;   // must be user-inputted - e.g. could be {0, 1, 45.3} for 2 binary covariates and 1 cts one (e.g. age)
          ////
          //// ---- Priors for beta:
          ////
          array[n_index_tests] matrix[2, max(max(n_covariates_nd), max(n_covariates_d))] prior_beta_mu_mean;
          array[n_index_tests] matrix<lower=0.0>[2, max(max(n_covariates_nd), max(n_covariates_d))] prior_beta_mu_SD;
          ////
          matrix<lower=0.0>[n_index_tests, 2] prior_beta_tau_SD;
          vector<lower=0.0>[2] prior_beta_sigma_SD; //// "sigma's (Nyaga et al. notation) are shared between tests.
          ////
          //// ---- Priors for raw_scale:
          ////
          array[n_index_tests] matrix[2, max(max(n_covariates_nd), max(n_covariates_d))] prior_raw_scale_mu_mean;
          array[n_index_tests] matrix<lower=0.0>[2, max(max(n_covariates_nd), max(n_covariates_d))] prior_raw_scale_mu_SD;
          ////
          matrix<lower=0.0>[n_index_tests, 2] prior_raw_scale_tau_SD; 
          vector<lower=0.0>[2] prior_raw_scale_sigma_SD; //// "sigma's (Nyaga et al. notation) are shared between tests.
          ////
          //// ---- Priors for box-cox:
          ////
          vector[n_index_tests] prior_boxcox_lambda_mean;
          vector<lower=0.0>[n_index_tests] prior_boxcox_lambda_SD;
          ////
          //// ---- Priors (and possible restrictions) for between-study correlations:
          ////
          real<lower=-1.0, upper=+1.0> beta_corr_lb;
          real<lower=beta_corr_lb, upper=+1.0>  beta_corr_ub;
          real<lower=0.0>  prior_beta_corr_LKJ;
          ////
          real<lower=-1.0, upper=+1.0> raw_scale_corr_lb;
          real<lower=raw_scale_corr_lb, upper=+1.0>  raw_scale_corr_ub;
          real<lower=0.0>  prior_raw_scale_corr_LKJ;
          ////
          //// ---- Other:
          ////
          int<lower=0, upper=1> softplus;
          int<lower=0, upper=1> box_cox;
          array[n_index_tests] vector[n_thr_max] cts_thr_values;
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
          int n_total_C_if_fixed = sum(n_thr); // total # of fixed-effect cutpoints
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
              n_total_study_thresh += n_studies * n_thr[t];
          }
          ////
          int n_total_comparisons = 0;
          for (t1 in 1:(n_index_tests - 1)) {
              for (t2 in (t1 + 1):n_index_tests) {
                  n_total_comparisons += n_thr[t1] * n_thr[t2]; 
              }
          }
}


parameters {
          array[n_index_tests] matrix[2, n_covariates_max] beta_mu;   
          array[n_index_tests] matrix[2, n_covariates_max] raw_scale_mu; 
          ////
          array[2] vector<lower=-5.0, upper=5.0>[n_index_tests] lambda;
          ////
          //// ---- "NMA" params:
          ////
          // vector<lower=0.0>[2] beta_sigma; 
          vector[2] log_beta_sigma_MU;
          vector<lower=0.0>[2] log_beta_sigma_SD;
          matrix[n_index_tests, 2] log_beta_sigma_z;
          ////
          matrix<lower=0.0>[n_index_tests, 2] beta_tau_raw;
          ////
          array[n_index_tests] matrix[n_studies, 2] beta_eta_z;  
          array[n_index_tests] matrix[n_studies, 2] beta_delta_z; 
          ////
          // vector<lower=0.0>[2] raw_scale_sigma;
          vector[2] log_raw_scale_sigma_MU;
          vector<lower=0.0>[2] log_raw_scale_sigma_SD;
          matrix[n_index_tests, 2] log_raw_scale_sigma_z;
          ////
          matrix<lower=0.0>[n_index_tests, 2] raw_scale_tau_raw;
          ////
          array[n_index_tests] matrix[n_studies, 2] raw_scale_eta_z;
          array[n_index_tests] matrix[n_studies, 2] raw_scale_delta_z;
          ////
          //// ---- Between-study corr's:
          ////
          real<lower=beta_corr_lb, upper=beta_corr_ub> beta_corr;       //// between-study corr (possibly restricted)
          real<lower=raw_scale_corr_lb, upper=raw_scale_corr_ub> raw_scale_corr;  //// between-study corr (possibly restricted)
}


transformed parameters {
          matrix<lower=0.0>[n_index_tests, 2] beta_tau;
          matrix<lower=0.0>[n_index_tests, 2] raw_scale_tau;
          if (compound_symmetry) {
              for (t in 1:n_index_tests) {
                for (c in 1:2) {
                   beta_tau[t, c] = beta_tau_raw[1, c];
                   raw_scale_tau[t, c] = raw_scale_tau_raw[1, c];
                }
              }
          } else {
              beta_tau = beta_tau_raw;
              raw_scale_tau = raw_scale_tau_raw;
          }
          ////
          array[n_index_tests] vector[2] log_beta_sigma;
          array[n_index_tests] vector[2] log_raw_scale_sigma;
          array[n_index_tests] vector[2] beta_sigma;
          array[n_index_tests] vector[2] raw_scale_sigma;
          ////
          for (t in 1:n_index_tests) {
             for (c in 1:2) {
               log_beta_sigma[t][c]      = log_beta_sigma_MU[c] + log_beta_sigma_SD[c] * log_beta_sigma_z[t, c];
               log_raw_scale_sigma[t][c] = log_raw_scale_sigma_MU[c] + log_raw_scale_sigma_SD[c] * log_raw_scale_sigma_z[t, c];
               if (hetero_sigma) {
                    beta_sigma[t][c] = exp(log_beta_sigma[t][c]);
                    raw_scale_sigma[t][c] = exp(log_raw_scale_sigma[t][c]);
               } else {
                    beta_sigma[t][c] = exp(log_beta_sigma[1][c]);
                    raw_scale_sigma[t][c] = exp(log_raw_scale_sigma[1][c]);
               }
             }
          }
          ////
          //// ---- Construct simple 2x2 (bivariate) between-study corr matrices:
          ////
          cholesky_factor_corr[2] beta_L_Omega = make_bivariate_L_Omega(beta_corr);
          cholesky_factor_corr[2] raw_scale_L_Omega = make_bivariate_L_Omega(raw_scale_corr);
          array[n_index_tests] cholesky_factor_cov[2]  beta_L_Sigma;
          array[n_index_tests] cholesky_factor_cov[2]  raw_scale_L_Sigma;
          ////
          for (t in 1:n_index_tests) {
               beta_L_Sigma[t] = diag_pre_multiply(beta_sigma[t], beta_L_Omega);
               raw_scale_L_Sigma[t] = diag_pre_multiply(raw_scale_sigma[t], raw_scale_L_Omega);
          }
          ////
          //// ---- Cutpoint params:
          ////
          array[2] vector[n_total_C_if_fixed] C_vec;  //// Global cutpoints for each test ("staggered" array/matrix using "n_thr[t]" to index correctly)
          ////
          //// ---- Construct cutpoints for each test:
          ////
          for (t in 1:n_index_tests) {
            for (c in 1:2) {
                  int n_thr_t = n_thr[t];
                  vector[n_thr_t] C_vec_test_t = ((box_cox == 0) ? 
                                                 log(cts_thr_values[t][1:n_thr[t]]) : 
                                                 fn_Stan_box_cox(cts_thr_values[t][1:n_thr[t]], lambda[c][t]));
                  C_vec[c] = update_test_values(C_vec[c], C_vec_test_t, start_index, end_index, t);
            }
          }
          ////
          //// ---- "NMA" params:
          //// ---- Declare Study-level random effects (eta in Nyaga notation) - eta[s, 1:2] ~ multi_normal({0, 0}, Sigma):
          ////
          array[n_index_tests] matrix[n_studies, 2] beta_eta        = init_array_of_matrices(n_studies, 2, n_index_tests, 0.0);
          array[n_index_tests] matrix[n_studies, 2] raw_scale_eta   = init_array_of_matrices(n_studies, 2, n_index_tests, 0.0);
          array[n_index_tests] matrix[n_studies, 2] beta_delta      = init_array_of_matrices(n_studies, 2, n_index_tests, 0.0);
          array[n_index_tests] matrix[n_studies, 2] raw_scale_delta = init_array_of_matrices(n_studies, 2, n_index_tests, 0.0);
          ////
          //// ---- Compute Study-level random effects (eta in Nyaga notation):
          ////
          for (t in 1:n_index_tests) { 
            for (s in 1:n_studies) {
              if (indicator_index_test_in_study[s, t] == 1) {
                 if (hetero_sigma) {
                      beta_eta[t][s, 1:2] = to_row_vector( beta_L_Sigma[t] * to_vector(beta_eta_z[t][s, 1:2]) );
                      raw_scale_eta[t][s, 1:2] = to_row_vector( raw_scale_L_Sigma[t] * to_vector(raw_scale_eta_z[t][s, 1:2]) );
                 } else  {
                      beta_eta[t][s, 1:2] = to_row_vector( beta_L_Sigma[1] * to_vector(beta_eta_z[1][s, 1:2]) );
                      raw_scale_eta[t][s, 1:2] = to_row_vector( raw_scale_L_Sigma[1] * to_vector(raw_scale_eta_z[1][s, 1:2]) );
                 }
              }
            }
          }
          ////
          for (t in 1:n_index_tests) {
              for (c in 1:2) {
                  beta_delta[t][1:n_studies, c]      = beta_tau[t, c]      * beta_delta_z[t][1:n_studies, c];       
                  raw_scale_delta[t][1:n_studies, c] = raw_scale_tau[t, c] * raw_scale_delta_z[t][1:n_studies, c];
              }
          }
          ////
          //// ---- Study-level random effects (after Cholesky decomposition):
          ////
          array[n_index_tests] matrix[n_studies, 2] beta_random;
          array[n_index_tests] matrix[n_studies, 2] raw_scale_random;
          for (t in 1:n_index_tests) {
                for (s in 1:n_studies) {
                     beta_random[t][s, 1:2]      =  beta_delta[t][s, 1:2]      + beta_eta[t][s, 1:2];
                     raw_scale_random[t][s, 1:2] =  raw_scale_delta[t][s, 1:2] + raw_scale_eta[t][s, 1:2];
                }
          }
          ////
          //// ---- Linear predictors for each disease statu + apply covariates for non-diseased and diseased groups:
          ////
          array[n_index_tests] matrix[n_studies, 2] Xbeta; // NOT local
          array[n_index_tests] matrix[n_studies, 2] Xraw_scale; // NOT local
          array[n_index_tests] matrix[n_studies, 2] Xscale; // NOT local
          ////
          for (t in 1:n_index_tests) {
                for (s in 1:n_studies) {
                    if (indicator_index_test_in_study[s, t] == 1) {
                        Xbeta[t][s, 1] = X_nd[t][s, 1:n_covariates_nd[t]] * to_vector(beta_mu[t][1, 1:n_covariates_nd[t]]) + 
                                         beta_random[t][s, 1];
                        Xbeta[t][s, 2] = X_d[t][s, 1:n_covariates_d[t]]   * to_vector(beta_mu[t][2, 1:n_covariates_d[t]])  + 
                                         beta_random[t][s, 2];
                        ////
                        Xraw_scale[t][s, 1] = X_nd[t][s, 1:n_covariates_nd[t]] * to_vector(raw_scale_mu[t][1, 1:n_covariates_nd[t]]) + 
                                              raw_scale_random[t][s, 1];
                        Xraw_scale[t][s, 2] = X_d[t][s, 1:n_covariates_d[t]]   * to_vector(raw_scale_mu[t][2, 1:n_covariates_d[t]])  + 
                                              raw_scale_random[t][s, 2];
                        ////
                        Xscale[t][s, ] = ((softplus == 1) ? softplus_scaled(Xraw_scale[t][s, ]) : exp(Xraw_scale[t][s, ]));  // local
                    }
                }
          }
}


model {
        ////
        //// ---- Priors:
        ////
        for (c in 1:2) {
           lambda[c] ~ normal(prior_boxcox_lambda_mean, prior_boxcox_lambda_SD);
        }
        ////
        for (t in 1:n_index_tests) {
            beta_mu[t][1, 1:n_covariates_nd[t]] ~ normal(
                                                  prior_beta_mu_mean[t][1, 1:n_covariates_nd[t]], 
                                                  prior_beta_mu_SD[t][1, 1:n_covariates_nd[t]]);  
            beta_mu[t][2, 1:n_covariates_d[t]]  ~ normal(
                                                  prior_beta_mu_mean[t][2, 1:n_covariates_d[t]],  
                                                  prior_beta_mu_SD[t][2, 1:n_covariates_d[t]]);
        }
        // to_vector(beta_mu)  ~ normal(to_vector(prior_beta_mu_mean), to_vector(prior_beta_mu_SD));
        ////
        to_vector(beta_tau) ~ normal(0.0, to_vector(prior_beta_tau_SD));  //// delta_{s, i, t} ~ normal(0, tau_{i, t}):
        beta_L_Omega ~ lkj_corr_cholesky(prior_beta_corr_LKJ);
        ////
        for (c in 1:2) {
            log_beta_sigma_MU[c] ~ normal(-1, 0.50);
            log_beta_sigma_SD[c] ~ normal(0, 0.25);
            // for (t in 1:n_index_tests) {
            //     target += log_beta_sigma[t][c]; //// ??? - J not needed as putting prior on heirarchical components.
            // }
            // beta_sigma_raw[c] ~ normal(0.0, prior_beta_sigma_SD[c]);
        }
        ////
        // to_vector(raw_scale_mu)  ~ normal(to_vector(prior_raw_scale_mu_mean), to_vector(prior_raw_scale_mu_SD));
        for (t in 1:n_index_tests) {
            raw_scale_mu[t][1, 1:n_covariates_nd[t]] ~ normal(
                                                       prior_raw_scale_mu_mean[t][1, 1:n_covariates_nd[t]], 
                                                       prior_raw_scale_mu_SD[t][1, 1:n_covariates_nd[t]]);  
            raw_scale_mu[t][2, 1:n_covariates_d[t]]  ~ normal(
                                                       prior_raw_scale_mu_mean[t][2, 1:n_covariates_d[t]],  
                                                       prior_raw_scale_mu_SD[t][2, 1:n_covariates_d[t]]);
        }
        ////
        to_vector(raw_scale_tau) ~ normal(0.0, to_vector(prior_beta_tau_SD));  //// delta_{s, i, t} ~ normal(0, tau_{i, t}):
        raw_scale_L_Omega ~ lkj_corr_cholesky(prior_raw_scale_corr_LKJ);
        ////
        for (c in 1:2) {
            log_raw_scale_sigma_MU[c] ~ normal(-1, 0.50);
            log_raw_scale_sigma_SD[c] ~ normal(0, 0.25);
            // for (t in 1:n_index_tests) {
            //     target += log_beta_sigma[t][c]; //// ??? - J not needed as putting prior on heirarchical components.
            // }
            // beta_sigma_raw[c] ~ normal(0.0, prior_beta_sigma_SD[c]);
        }
        target += std_normal_lpdf(to_vector(log_beta_sigma_z));
        target += std_normal_lpdf(to_vector(log_raw_scale_sigma_z));
        ////
        for (t in 1:n_index_tests) {
            target += std_normal_lpdf(to_vector(beta_eta_z[t]));   
            target += std_normal_lpdf(to_vector(raw_scale_eta_z[t]));   
            ////
            target += std_normal_lpdf(to_vector(beta_delta_z[t]));   
            target += std_normal_lpdf(to_vector(raw_scale_delta_z[t])); 
        }
        ////
        //// ---- Log-likelihood:
        ////
        {
              for (t in 1:n_index_tests) { 
                    ////
                    array[2] vector[n_thr[t]] C_vec_t;
                    for (c in 1:2) {
                        vector[n_thr[t]] C_vec_t_given_c = get_test_values(C_vec[c], start_index, end_index, t);
                        C_vec_t[c] = C_vec_t_given_c;
                    }
                    ////
                    array[2] matrix[n_studies, n_thr[t]] latent_surv_t = map_latent_surv_to_fixed_hetero_C(
                                                                         C_vec_t, Xbeta[t], Xscale[t], n_studies, 
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
          ////
          //// ---- Compute between-study variance-covariance matrices:
          ////
          corr_matrix[2] beta_Omega      = multiply_lower_tri_self_transpose(beta_L_Omega);
          corr_matrix[2] raw_scale_Omega = multiply_lower_tri_self_transpose(raw_scale_L_Omega);
          array[n_index_tests] matrix[2, 2] beta_Sigma;
          array[n_index_tests] matrix[2, 2] raw_scale_Sigma;
          for (t in 1:n_index_tests) {
              beta_Sigma[t]      = multiply_lower_tri_self_transpose(beta_L_Sigma[t]);
              raw_scale_Sigma[t] = multiply_lower_tri_self_transpose(raw_scale_L_Sigma[t]);
          }
          ////
          //// ---- Calculate summary accuracy for each test - "baseline" covariate values:
          ////
          array[n_index_tests] real Xbeta_baseline_nd;
          array[n_index_tests] real Xbeta_baseline_d;
          array[n_index_tests] real Xraw_scale_baseline_nd;
          array[n_index_tests] real Xraw_scale_baseline_d;
          array[n_index_tests] real scale_nd_baseline;
          array[n_index_tests] real scale_d_baseline;
          ////
          for (t in 1:n_index_tests) {
              Xbeta_baseline_nd[t] = dot_product( baseline_case_nd[t][1:n_covariates_nd[t]], 
                                                  to_vector(beta_mu[t][1, 1:n_covariates_nd[t]]));
              Xbeta_baseline_d[t] = dot_product( baseline_case_d[t][1:n_covariates_d[t]], 
                                                 to_vector(beta_mu[t][2, 1:n_covariates_d[t]]));
              ////
              Xraw_scale_baseline_nd[t] = dot_product( baseline_case_nd[t][1:n_covariates_nd[t]], 
                                                       to_vector(raw_scale_mu[t][1, 1:n_covariates_nd[t]]));
              Xraw_scale_baseline_d[t] = dot_product( baseline_case_d[t][1:n_covariates_d[t]], 
                                                      to_vector(raw_scale_mu[t][2, 1:n_covariates_d[t]]));
              ////
              scale_nd_baseline[t] = ((softplus == 1) ? 
                                     softplus_scaled(Xraw_scale_baseline_nd[t]) : exp(Xraw_scale_baseline_nd[t]));
              scale_d_baseline[t] = ((softplus == 1) ? 
                                    softplus_scaled(Xraw_scale_baseline_d[t]) : exp(Xraw_scale_baseline_d[t]));
          }
          ////
          //// ---- Calculate baseline Se/Sp for each test:
          ////
          array[n_index_tests] vector[n_thr_max] Fp_baseline;
          array[n_index_tests] vector[n_thr_max] Sp_baseline;
          array[n_index_tests] vector[n_thr_max] Se_baseline;
          ////
          for (t in 1:n_index_tests) {
              //// Initialize with -1 for unused thresholds
              Fp_baseline[t] = rep_vector(-1.0, n_thr_max);
              Sp_baseline[t] = rep_vector(-1.0, n_thr_max);
              Se_baseline[t] = rep_vector(-1.0, n_thr_max);
              
              //// Calculate for actual thresholds
              vector[n_thr[t]] C_vec_t_nd = get_test_values(C_vec[1], start_index, end_index, t);
              vector[n_thr[t]] C_vec_t_d  = get_test_values(C_vec[2], start_index, end_index, t);
              
              for (k in 1:n_thr[t]) {
                  Fp_baseline[t][k] = (use_probit_link == 1) ? 
                                      Phi(-(C_vec_t_nd[k] - Xbeta_baseline_nd[t])/scale_nd_baseline[t]) : 
                                      inv_logit(-(C_vec_t_nd[k] - Xbeta_baseline_nd[t])/scale_nd_baseline[t]);
                  Sp_baseline[t][k] = 1.0 - Fp_baseline[t][k];
                  Se_baseline[t][k] = (use_probit_link == 1) ? 
                                      Phi(-(C_vec_t_d[k] - Xbeta_baseline_d[t])/scale_d_baseline[t]) : 
                                      inv_logit(-(C_vec_t_d[k] - Xbeta_baseline_d[t])/scale_d_baseline[t]);
              }
          }
          ////
          //// ---- Calculate predictive accuracy for each test (NMA structure):
          ////
          array[n_index_tests] vector[n_thr_max] Fp_baseline_pred;
          array[n_index_tests] vector[n_thr_max] Sp_baseline_pred;
          array[n_index_tests] vector[n_thr_max] Se_baseline_pred;
          ////
          for (t in 1:n_index_tests) {
              //// ---- Draw shared study-level effects
              vector[2] beta_eta_pred      = multi_normal_cholesky_rng(rep_vector(0.0, 2), beta_L_Sigma[t]);
              vector[2] raw_scale_eta_pred = multi_normal_cholesky_rng(rep_vector(0.0, 2), raw_scale_L_Sigma[t]);
              //// ---- Initialize with -1 for unused thresholds
              Fp_baseline_pred[t] = rep_vector(-1.0, n_thr_max);
              Sp_baseline_pred[t] = rep_vector(-1.0, n_thr_max);
              Se_baseline_pred[t] = rep_vector(-1.0, n_thr_max);
              //// ---- Draw test-specific deviations
              real beta_delta_pred_nd = normal_rng(0.0, beta_tau[t, 1]);
              real beta_delta_pred_d  = normal_rng(0.0, beta_tau[t, 2]);
              real raw_scale_delta_pred_nd = normal_rng(0.0, raw_scale_tau[t, 1]);
              real raw_scale_delta_pred_d  = normal_rng(0.0, raw_scale_tau[t, 2]);
              //// ---- Combine to get predictive values
              real Xbeta_baseline_pred_nd = Xbeta_baseline_nd[t] + beta_eta_pred[1] + beta_delta_pred_nd;
              real Xbeta_baseline_pred_d  = Xbeta_baseline_d[t]  + beta_eta_pred[2] + beta_delta_pred_d;
              ////
              real Xraw_scale_baseline_pred_nd = Xraw_scale_baseline_nd[t] + raw_scale_eta_pred[1] + raw_scale_delta_pred_nd;
              real Xraw_scale_baseline_pred_d  = Xraw_scale_baseline_d[t]  + raw_scale_eta_pred[2] + raw_scale_delta_pred_d;
              ////
              real scale_nd_baseline_pred = ((softplus == 1) ? 
                                            softplus_scaled(Xraw_scale_baseline_pred_nd) : exp(Xraw_scale_baseline_pred_nd));
              real scale_d_baseline_pred  = ((softplus == 1) ? 
                                            softplus_scaled(Xraw_scale_baseline_pred_d)  : exp(Xraw_scale_baseline_pred_d));
              //// ---- Calculate predictive Se/Sp
              vector[n_thr[t]] C_vec_t_nd = get_test_values(C_vec[1], start_index, end_index, t);
              vector[n_thr[t]] C_vec_t_d  = get_test_values(C_vec[2], start_index, end_index, t);
              ////
              for (k in 1:n_thr[t]) {
                  Fp_baseline_pred[t][k] = (use_probit_link == 1) ? 
                                           Phi(-(C_vec_t_nd[k] - Xbeta_baseline_pred_nd)/scale_nd_baseline_pred) : 
                                           inv_logit(-(C_vec_t_nd[k] - Xbeta_baseline_pred_nd)/scale_nd_baseline_pred);
                  Sp_baseline_pred[t][k] = 1.0 - Fp_baseline_pred[t][k];
                  Se_baseline_pred[t][k] = (use_probit_link == 1) ? 
                                           Phi(-(C_vec_t_d[k]  - Xbeta_baseline_pred_d)/scale_d_baseline_pred)   : 
                                           inv_logit(-(C_vec_t_d[k]  - Xbeta_baseline_pred_d)/scale_d_baseline_pred);
              }
          }
          ////
          //// ---- Log-lik + study-specific accuracy computation:
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
          vector[n_studies] log_lik_study = rep_vector(0.0, n_studies);
          {
              for (t in 1:n_index_tests) {
                  array[2] vector[n_thr[t]] C_vec_t_dummy;
                  for (c in 1:2) {
                      vector[n_thr[t]] C_vec_t = get_test_values(C_vec[c], start_index, end_index, t);
                      C_vec_t_dummy[c] = C_vec_t;
                  }
                  ////
                  array[2] matrix[n_studies, n_thr[t]] latent_surv_t = map_latent_surv_to_fixed_hetero_C(
                                                                       C_vec_t_dummy, Xbeta[t], Xscale[t], n_studies,
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
                  int flat_idx = 1;
                  for (s in 1:n_studies) {
                      for (k in 1:n_thr[t]) {
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
          //// ---- Between-study heterogeneity metrics for NMA: --------------------------------------------------------
          ////
          //// Total variance for each test (Nyaga et al. formulation)
          array[n_index_tests] vector[2] total_SD_beta;
          array[n_index_tests] vector[2] total_SD_raw_scale;
          ////
          for (t in 1:n_index_tests) {
            for (c in 1:2) {
              total_SD_beta[t][c]      = sqrt(square(beta_sigma[t][c]) + square(beta_tau[t, c]));
              total_SD_raw_scale[t][c] = sqrt(square(raw_scale_sigma[t][c]) + square(raw_scale_tau[t, c]));
            }
          }
          ////
          //// Proportion of variance explained by test-specific component
          array[n_index_tests] vector[2] prop_test_SD_beta;
          array[n_index_tests] vector[2] prop_test_SD_raw_scale;
          ////
          for (t in 1:n_index_tests) {
            for (c in 1:2) {
               prop_test_SD_beta[t][c]      = beta_tau[t, c] / total_SD_beta[t][c];
               prop_test_SD_raw_scale[t][c] = raw_scale_tau[t, c] / total_SD_raw_scale[t][c];
            }
          }
          //// Estimate of "total_SD_inc_C":
          array[n_index_tests] vector[2] total_SD_inc_C;
          ////
          for (t in 1:n_index_tests) {
            for (c in 1:2) {
              ////
              //// Baseline scale
              real scale_baseline = (c == 1) ? scale_nd_baseline[t] : scale_d_baseline[t];
              ////
              //// SD on latent scale from location parameters
              real latent_SD_from_location = total_SD_beta[t][c] / scale_baseline;
              ////
              //// Approximate CV of scale (SD on log scale ≈ CV on original scale for small SD)
              real scale_relative_SD = total_SD_raw_scale[t][c];
              ////
              //// First-order approximation: multiply by (1 + CV²)^0.5
              total_SD_inc_C[t][c] = latent_SD_from_location * sqrt(1 + square(scale_relative_SD));
            }
          }
          ////
          matrix[n_index_tests, n_index_tests] rho12;
          array[2] matrix[n_index_tests, n_index_tests] rho;
          {
             array[3, 2] matrix[n_index_tests, n_index_tests] outs = compute_Nyaga_NMA_summaries_hetero_SDs(
                                                                     n_index_tests, beta_tau, beta_Sigma, beta_Omega);
             rho12 = outs[3, 1];
             rho   = outs[2, ];
          }
          ////
          //// ---- AUC computation for each test: ----------------------------------------------------------------------
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






