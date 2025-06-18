

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
          array[n_index_tests] int n_covariates; //// For R+G - cannot have SEPERATE covariates for D+ and D-, so only one n_covariates!
          int n_covariates_max; // this is the max across ALL index tests
          array[n_index_tests] matrix[n_studies, max(n_covariates)] X;
          array[n_index_tests] vector[max(n_covariates)] baseline_case; // must be user-inputted - e.g. could be {0, 1, 45.3} for 2 binary covariates and 1 cts one (e.g. age)s
          ////
          //// ---- Priors for SHARED "beta":
          ////
          array[n_index_tests] vector[n_covariates_max] prior_beta_mu_mean;
          array[n_index_tests] vector<lower=0.0>[n_covariates_max] prior_beta_mu_SD;
          vector<lower=0.0>[n_index_tests] prior_beta_tau_SD;
          real<lower=0.0> prior_beta_sigma_SD;
          ////
          //// ---- Priors for SHARED "raw_scale":
          //// 
          array[n_index_tests] vector[n_covariates_max] prior_raw_scale_mu_mean;
          array[n_index_tests] vector<lower=0.0>[n_covariates_max] prior_raw_scale_mu_SD;
          vector<lower=0.0>[n_index_tests] prior_raw_scale_tau_SD; 
          real<lower=0.0> prior_raw_scale_sigma_SD; //// "sigma's (Nyaga et al. notation) are shared between tests.
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
          int n_total_C_if_fixed = sum(n_thr); // total # of fixed-effect cutpoints
          int n_tests = n_index_tests;
          ////
          //// ---- Calculate indices:
          ////
          array[n_tests] int start_index = calculate_start_indices(n_thr, n_tests);
          array[n_tests] int end_index   = calculate_end_indices(n_thr, n_tests, start_index);
          ////
          //// ---- HSROC stuff:
          ////
          real mult_nd = (1)*(-1.0);
          real mult_d  = (1)*(+1.0); //// hence: mult_d > mult_nd
}


parameters {
          array[n_index_tests] vector[n_index_tests] beta_mu;
          array[n_index_tests] vector[n_index_tests] raw_scale_mu;
          ////
          vector[n_total_C_if_fixed] C_raw_vec;  //// RAW LOG-DIFFERENCES - Global cutpoints for each test (staggered array/matrix using "n_thr[t]" to index correctly)
          ////
          //// ---- "NMA" params for SHARED "beta":
          ////
          real<lower=0.0> beta_sigma;              //// Variance component - Between-study SD (Nyaga's σ)
          vector<lower=0.0>[n_index_tests] beta_tau; //// Variance component - Test-specific SD (Nyaga's τ) - delta_{s, c, t} ~ normal(0, tau_{c, t}).
          vector[n_studies] beta_eta_z;              //// Standard normal RVs for study-level effects - eta[s, 1:2] ~ multi_normal({0, 0}, Sigma).
          array[n_index_tests] vector[n_studies] beta_delta_z; //// Standard normal RVs for test-specific effects
          ////
          //// ---- "NMA" params for SHARED "raw_scale":
          ////
          real<lower=0.0> raw_scale_sigma;              //// Variance component - Between-study SD (Nyaga's σ)
          vector<lower=0.0>[n_index_tests] raw_scale_tau; //// Variance component - Test-specific SD (Nyaga's τ) - delta_{s, c, t} ~ normal(0, tau_{c, t}).
          vector[n_studies] raw_scale_eta_z;              //// Standard normal RVs for study-level effects - eta[s, 1:2] ~ multi_normal({0, 0}, Sigma).
          array[n_index_tests] vector[n_studies] raw_scale_delta_z; //// Standard normal RVs for test-specific effects
}


transformed parameters {
          ////
          //// ---- Construct (global) cutpoints:
          ////
          vector[n_total_C_if_fixed] C_vec;  //// Global cutpoints for each test ("staggered" array/matrix using "n_thr[t]" to index correctly)
          {
            int counter = 1;
            for (t in 1:n_index_tests) {
                  int n_thr_t = n_thr[t];
                  vector[n_thr_t] C_raw_vec_test_t = get_test_values(C_raw_vec, start_index, end_index, t);
                  vector[n_thr_t] C_vec_test_t = construct_C(C_raw_vec_test_t, softplus);
                  C_vec = update_test_values(C_vec, C_vec_test_t, start_index, end_index, t);
            }
          }
          ////
          //// ---- "NMA" params:
          ////
          vector[n_studies] beta_eta      = rep_vector(0.0, n_studies);
          vector[n_studies] raw_scale_eta = rep_vector(0.0, n_studies);
          array[n_index_tests] vector[n_studies] beta_delta;
          array[n_index_tests] vector[n_studies] raw_scale_delta;
          ////
          //// ---- Compute Study-level random effects (eta in Nyaga notation):
          ////
          for (s in 1:n_studies) {
              beta_eta[s]      = beta_sigma      * beta_eta_z[s];       //// beta_eta[s] ~ normal({0.0}, sigma);
              raw_scale_eta[s] = raw_scale_sigma * raw_scale_eta_z[s];  //// beta_eta[s] ~ normal({0.0}, sigma);
          }
          for (c in 1:2) {
              //// Compute test-specific deviations ("delta" in Nyaga notation):
              for (t in 1:n_index_tests) { 
                 beta_delta[t]      = beta_tau[t]      * beta_delta_z[t]; //// beta_delta[t][s] ~ normal(0.0, beta_tau[t]);
                 raw_scale_delta[t] = raw_scale_tau[t] * raw_scale_delta_z[t];
              }
          }
          ////
          //// ---- Study-level random effects (after Cholesky decomposition):
          ////
          array[n_index_tests] vector[n_studies] beta_random;
          array[n_index_tests] vector[n_studies] raw_scale_random;
          for (t in 1:n_index_tests) {
                for (s in 1:n_studies) {
                     beta_random[t][s]      =  beta_delta[t][s]      + beta_eta[s];
                     raw_scale_random[t][s] =  raw_scale_delta[t][s] + raw_scale_eta[s];
                }
          }
          ////
          //// ---- Linear predictors for each disease statu + apply covariates for non-diseased and diseased groups:
          ////
          array[n_index_tests] matrix[n_studies, 2] Xlocations; // NOT local
          array[n_index_tests] matrix[n_studies, 2] Xraw_scale; // NOT local
          array[n_index_tests] matrix[n_studies, 2] Xscale;     // NOT local
          ////
          for (t in 1:n_index_tests) {
                for (s in 1:n_studies) {
                    real Xbeta_baseline  = sum(to_vector(X[t][s, 1:n_covariates[t]]) .* beta_mu[t][1:n_covariates[t]]) + beta_random[t][s];
                    Xlocations[t][s, 1] = mult_nd*Xbeta_baseline; 
                    Xlocations[t][s, 2] = mult_d*Xbeta_baseline; 
                    ////
                    real Xraw_scale_baseline = sum(to_vector(X[t][s, 1:n_covariates[t]]) .* raw_scale_mu[t][1:n_covariates[t]]) + raw_scale_random[t][s]; // "Xgamma"
                    Xscale[t][s, 1] =  ((softplus == 1) ? softplus_scaled( mult_nd*Xraw_scale_baseline) : exp( mult_nd*Xraw_scale_baseline)); 
                    Xscale[t][s, 2] =  ((softplus == 1) ? softplus_scaled( mult_d*Xraw_scale_baseline)  : exp( mult_d*Xraw_scale_baseline));
                }
          }
}


model {
          ////
          //// ---- Priors:
          ////
          for (t in 1:n_index_tests) {
              beta_mu[t][1:n_covariates[t]] ~ normal(
                                              prior_beta_mu_mean[t][1:n_covariates[t]], 
                                              prior_beta_mu_SD[t][1:n_covariates[t]]);
          }
          ////
          to_vector(beta_tau) ~ normal(0.0, to_vector(prior_beta_tau_SD));  //// delta_{s, i, t} ~ normal(0, tau_{i, t}):
          beta_sigma ~ normal(0.0, prior_beta_sigma_SD);      //// eta[s, i] ~ normal(0, sigma[i]):
          ////
          for (t in 1:n_index_tests) {
              raw_scale_mu[t][1:n_covariates[t]] ~ normal(
                                                   prior_raw_scale_mu_mean[t][1:n_covariates[t]], 
                                                   prior_raw_scale_mu_SD[t][1:n_covariates[t]]);
          }
          to_vector(raw_scale_tau) ~ normal(0.0, to_vector(prior_raw_scale_tau_SD));  //// delta_{s, i, t} ~ normal(0, tau_{i, t}):
          raw_scale_sigma ~ normal(0.0, prior_raw_scale_sigma_SD);      //// eta[s, i] ~ normal(0, sigma[i]):
          ////
          //// ---- Induced-dirichlet ** Prior ** model:
          ////
          {
              // array[n_index_tests] matrix[n_studies, n_thr_max] Ind_Dir_anchor;
              for (t in 1:n_index_tests) {
                    int n_thr_t = n_thr[t];
                    int n_cat_t = n_thr_t + 1;
                    vector[n_thr_t] C_vec_t = get_test_values(C_vec, start_index, end_index, t);
                    ////
                    vector[n_cat_t] Ind_Dir_ord_prob = (use_probit_link == 1) ?
                                                       cumul_probs_to_ord_probs(Phi(C_vec_t)) : 
                                                       cumul_probs_to_ord_probs(inv_logit(C_vec_t));
                    ////
                    Ind_Dir_ord_prob ~ induced_dirichlet_given_C(
                                       C_vec_t, prior_dirichlet_alpha[t][1:n_cat_t], use_probit_link);
              }
          }
          ////
          //// ---- Likelihood / Observational Model:
          ////
          //// (part of between-test / between-study model, NOT prior) - eta[s, i] ~ normal(0, sigma[i]):
          target += std_normal_lpdf(to_vector(beta_eta_z)); 
          target += std_normal_lpdf(to_vector(raw_scale_eta_z)); 
          ////
          //// (part of between-test / between-study model, NOT prior)  -   delta_{s, i, t} ~ normal(0, tau_{i, t}):
          for (t in 1:n_index_tests) {
              target += std_normal_lpdf(to_vector(beta_delta_z[t]));
              target += std_normal_lpdf(to_vector(raw_scale_delta_z[t])); 
              ////
              int n_thr_t = n_thr[t];
              vector[n_thr_t] C_raw_vec_test_t = get_test_values(C_raw_vec, start_index, end_index, t);
              target += raw_C_to_C_log_det_J_lp(C_raw_vec_test_t, softplus);
          }
          ////
          //// ---- Log-likelihood:
          ////
          {
                for (t in 1:n_index_tests) {
                      ////
                      vector[n_thr[t]] C_vec_t;
                      {
                          vector[n_thr[t]] C_vec_t_given_c = get_test_values(C_vec, start_index, end_index, t);
                          C_vec_t = C_vec_t_given_c;
                      }
                      ////
                      array[2] matrix[n_studies, n_thr[t]] latent_surv_t = map_latent_surv_prob_to_fixed_C(
                                                                           C_vec_t, Xlocations[t], Xscale[t], n_studies, n_obs_cutpoints[t, ], 
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
      //// ---- Calculate summary accuracy (using mean parameters):
      ////
      array[n_index_tests] vector[n_covariates_max] location_nd_mu;
      array[n_index_tests] vector[n_covariates_max] location_d_mu;
      array[n_index_tests] vector[n_covariates_max] scale_nd_mu;
      array[n_index_tests] vector[n_covariates_max] scale_d_mu;
      ////
      for (t in 1:n_index_tests) {
          location_nd_mu[t][1:n_covariates[t]] = mult_nd * beta_mu[t][1:n_covariates[t]];
          location_d_mu[t][1:n_covariates[t]]  = mult_d * beta_mu[t][1:n_covariates[t]];
          ////
          scale_nd_mu[t][1:n_covariates[t]] = ((softplus == 1) ? 
                                              softplus_scaled(mult_nd * raw_scale_mu[t][1:n_covariates[t]]) : 
                                              exp(mult_nd * raw_scale_mu[t][1:n_covariates[t]]));
          scale_d_mu[t][1:n_covariates[t]]  = ((softplus == 1) ? 
                                              softplus_scaled(mult_d * raw_scale_mu[t][1:n_covariates[t]])  : 
                                              exp(mult_d * raw_scale_mu[t][1:n_covariates[t]]));
      }
      ////
      //// ---- Calculate summary accuracy for each test - "baseline" covariate values:
      ////
      array[n_index_tests] real Xbeta_baseline;
      array[n_index_tests] real Xbeta_baseline_nd;
      array[n_index_tests] real Xbeta_baseline_d;
      ////
      array[n_index_tests] real Xraw_scale_baseline;
      array[n_index_tests] real scale_nd_baseline;
      array[n_index_tests] real scale_d_baseline;
      ////
      for (t in 1:n_index_tests) {
          Xbeta_baseline[t] = dot_product(baseline_case[t][1:n_covariates[t]], beta_mu[t][1:n_covariates[t]]);
          Xbeta_baseline_nd[t] = mult_nd * Xbeta_baseline[t];
          Xbeta_baseline_d[t]  = mult_d * Xbeta_baseline[t];
          ////
          Xraw_scale_baseline[t] = dot_product(baseline_case[t][1:n_covariates[t]], raw_scale_mu[t][1:n_covariates[t]]);
          scale_nd_baseline[t] = ((softplus == 1) ? softplus_scaled(mult_nd * Xraw_scale_baseline[t]) : exp(mult_nd * Xraw_scale_baseline[t]));
          scale_d_baseline[t]  = ((softplus == 1) ? softplus_scaled(mult_d * Xraw_scale_baseline[t])  : exp(mult_d * Xraw_scale_baseline[t]));
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
          //// ---- Get cutpoints for this test
          vector[n_thr[t]] C_vec_t = get_test_values(C_vec, start_index, end_index, t);
          ////
          //// ---- Calculate for actual thresholds
          for (k in 1:n_thr[t]) {
              Fp_baseline[t][k] = (use_probit_link == 1) ? 
                                  Phi(-(C_vec_t[k] - Xbeta_baseline_nd[t])/scale_nd_baseline[t]) : 
                                  inv_logit(-(C_vec_t[k] - Xbeta_baseline_nd[t])/scale_nd_baseline[t]);
              Sp_baseline[t][k] = 1.0 - Fp_baseline[t][k];
              Se_baseline[t][k] = (use_probit_link == 1) ? 
                                  Phi(-(C_vec_t[k] - Xbeta_baseline_d[t])/scale_d_baseline[t]) : 
                                  inv_logit(-(C_vec_t[k] - Xbeta_baseline_d[t])/scale_d_baseline[t]);
          }
      }
      ////
      //// ---- Calculate predictive accuracy for each test:
      ////
      array[n_index_tests] vector[n_thr_max] Fp_baseline_pred;
      array[n_index_tests] vector[n_thr_max] Sp_baseline_pred;
      array[n_index_tests] vector[n_thr_max] Se_baseline_pred;
      ////
      //// ---- Draw shared study-level effects
      real beta_eta_pred      = normal_rng(0.0, beta_sigma);
      real raw_scale_eta_pred = normal_rng(0.0, raw_scale_sigma);
      ////
      for (t in 1:n_index_tests) {
          //// ---- Initialize with -1 for unused thresholds
          Fp_baseline_pred[t] = rep_vector(-1.0, n_thr_max);
          Sp_baseline_pred[t] = rep_vector(-1.0, n_thr_max);
          Se_baseline_pred[t] = rep_vector(-1.0, n_thr_max);
          ////
          //// ---- Draw test-specific deviations
          real beta_delta_pred      = normal_rng(0.0, beta_tau[t]);
          real raw_scale_delta_pred = normal_rng(0.0, raw_scale_tau[t]);
          ////
          //// ---- Combine to get predictive values
          real Xbeta_baseline_pred = dot_product(baseline_case[t][1:n_covariates[t]], beta_mu[t][1:n_covariates[t]]) + 
                                     beta_eta_pred + beta_delta_pred;
          real Xbeta_baseline_pred_nd = mult_nd * Xbeta_baseline_pred;
          real Xbeta_baseline_pred_d  = mult_d * Xbeta_baseline_pred;
          ////
          real Xraw_scale_baseline_pred = dot_product(baseline_case[t][1:n_covariates[t]], raw_scale_mu[t][1:n_covariates[t]]) + 
                                          raw_scale_eta_pred + raw_scale_delta_pred;
          real scale_baseline_pred_nd = ((softplus == 1) ? 
                                        softplus_scaled(mult_nd * Xraw_scale_baseline_pred) : exp(mult_nd * Xraw_scale_baseline_pred));
          real scale_baseline_pred_d  = ((softplus == 1) ? 
                                        softplus_scaled(mult_d * Xraw_scale_baseline_pred)  : exp(mult_d * Xraw_scale_baseline_pred));
          ////
          //// ---- Get cutpoints for this test
          vector[n_thr[t]] C_vec_t = get_test_values(C_vec, start_index, end_index, t);
          ////
          //// ---- Calculate predictive Se/Sp
          for (k in 1:n_thr[t]) {
              Fp_baseline_pred[t][k] = (use_probit_link == 1) ? 
                                       Phi(-(C_vec_t[k] - Xbeta_baseline_pred_nd)/scale_baseline_pred_nd) : 
                                       inv_logit(-(C_vec_t[k] - Xbeta_baseline_pred_nd)/scale_baseline_pred_nd);
              Sp_baseline_pred[t][k] = 1.0 - Fp_baseline_pred[t][k];
              Se_baseline_pred[t][k] = (use_probit_link == 1) ? 
                                       Phi(-(C_vec_t[k] - Xbeta_baseline_pred_d)/scale_baseline_pred_d) : 
                                       inv_logit(-(C_vec_t[k] - Xbeta_baseline_pred_d)/scale_baseline_pred_d);
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
              vector[n_thr[t]] C_vec_t = get_test_values(C_vec, start_index, end_index, t);
              ////
              array[2] matrix[n_studies, n_thr[t]] latent_surv_t = map_latent_surv_prob_to_fixed_C(
                                                                   C_vec_t, Xlocations[t], Xscale[t], n_studies, n_obs_cutpoints[t, ], 
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
      //// ---- "Ordinal-bivariate equivalent" params:
      ////
      array[n_index_tests, 2] matrix[n_thr_max, n_covariates_max] biv_equiv_C;
      for (t in 1:n_index_tests) {
          vector[n_thr[t]] C_vec_t = get_test_values(C_vec, start_index, end_index, t);
          for (x_i in 1:n_covariates[t]) {
              for (k in 1:n_thr[t]) {
                  biv_equiv_C[t, 1][k, x_i] = C_vec_t[k] / scale_nd_mu[t][x_i];
                  biv_equiv_C[t, 2][k, x_i] = C_vec_t[k] / scale_d_mu[t][x_i];
              }
          }
      }
}






