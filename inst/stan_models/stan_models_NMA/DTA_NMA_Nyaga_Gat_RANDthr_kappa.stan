

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
          //// ---- Priors for beta:
          ////
          matrix[n_index_tests, 2] prior_beta_mu_mean;
          matrix<lower=0.0>[n_index_tests, 2] prior_beta_mu_SD;
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
          // array[n_index_tests] vector<lower=0.0>[n_cat_max] prior_dirichlet_alpha;
          array[n_index_tests] vector<lower=0.0>[n_cat_max] prior_alpha_mean;
          array[n_index_tests] vector<lower=0.0>[n_cat_max] prior_alpha_SD;
          ////
          //// ---- Other:
          ////
          int<lower=0, upper=1> softplus;
          real<lower=0.0> alpha_lb;
          int n_total_C_if_random;
          int n_total_pooled_cat;
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
          int n_tests = n_index_tests;
          ////
          //// ---- HSROC stuff:
          ////
          real mult_nd = (1)*(-0.5);
          real mult_d  = (1)*(+0.5); //// hence: mult_d > mult_nd
}


parameters {
          array[n_index_tests] vector[n_index_tests] beta_mu;
          array[n_index_tests] vector[n_index_tests] raw_scale_mu;
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
          ////
          //// ---- Induced-Dirichlet:
          ////
          vector[n_total_C_if_random] C_raw_vec;  //// RAW LOG-DIFFERENCES - Global cutpoints for each test (staggered array/matrix using "n_thr[t]" to index correctly)
          vector<lower=alpha_lb, upper=1000.0>[n_total_pooled_cat] alpha_vec;
          // vector[n_total_pooled_thr] phi_raw_vec; //// for "kappa" or "SD" parameterisations (not "alpha" param)
}


transformed parameters {
          ////
          //// ---- Construct study-specific cutpoints:
          ////
          vector[n_total_C_if_random] C_vec;  //// Global cutpoints for each test ("staggered" array/matrix using "n_thr[t]" to index correctly)
          for (t in 1:n_index_tests) {
                  int n_thr_t = n_thr[t];
                  // int n_thr_t_all_studies = n_thr_t * n_studies;
                  // int start_index_test_t  = 1 + n_thr_t_all_studies*(t - 1); // ????? check
                  // //// get all threshold params for test t for all studies:
                  // vector[n_thr_t_all_studies] C_raw_vec_test_t_all_studies = segment(C_raw_vec, start_index_test_t, n_thr_t_all_studies); 
                  //// empty containers to fill in inner loop:
                  vector[n_thr_t] C_raw_vec_test_t;
                  vector[n_thr_t] C_vec_test_t;
                  ////
                  for (s in 1:n_studies) {
                            C_raw_vec_test_t = get_test_t_study_s_segment_values(C_raw_vec, t, s, n_thr, n_studies);
                            // int start_index_test_t_study_s =  1 + (s - 1)*n_thr_t;
                            // C_raw_vec_test_t = segment(C_raw_vec_test_t_all_studies, start_index_test_t_study_s, n_thr_t);
                            C_vec_test_t = construct_C(C_raw_vec_test_t, softplus);
                            C_vec = update_test_t_study_s_segment_values(C_vec, C_vec_test_t, t, s, n_thr, n_studies);
                  }
          }
          ////
          //// ---- "NMA" params:
          ////
          vector[n_studies] beta_eta = rep_vector(0.0, n_studies);
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
              beta_mu[t][1:n_covariates[t]] ~ normal(prior_beta_mu_mean[t][1:n_covariates[t]], prior_beta_mu_SD[t][1:n_covariates[t]]);
          }
          ////
          to_vector(beta_tau) ~ normal(0.0, to_vector(prior_beta_tau_SD));  //// delta_{s, i, t} ~ normal(0, tau_{i, t}):
          ////
          beta_sigma ~ normal(0.0, prior_beta_sigma_SD);      //// eta[s, i] ~ normal(0, sigma[i]):
          ////
          //// ---- Induced-dirichlet ** Prior ** model:
          ////
          // array[n_index_tests] matrix[n_studies, n_thr_max] Ind_Dir_anchor;
          for (t in 1:n_index_tests) {
                
                    int n_thr_t = n_thr[t];
                    int n_cat_t = n_thr_t + 1;
                    ////
                    //// ---- Prior for Induced-Dirichlet pooled "alpha" parameters:
                    ////
                    int start_alpha_index = ((t == 1) ? (1 + sum(n_cat[1:(t - 1)])) : 1);
                    vector[n_cat_t] alpha_t = segment(alpha_vec, start_alpha_index, n_cat_t);
                    alpha_t ~ normal(prior_alpha_mean[t][1:n_cat_t], prior_alpha_SD[t][1:n_cat_t]);
                    //// Containers to fill:
                    vector[n_thr_t] Ind_Dir_cumul;
                    vector[n_thr_t] Ind_Dir_cumul_prob;
                    vector[n_cat_t] Ind_Dir_ord_prob;
                    vector[n_thr_t] C_vec_t_study_s;
                    for (s in 1:n_studies) {
                            C_vec_t_study_s = get_test_t_study_s_segment_values(C_vec, t, s, n_thr, n_studies);
                            Ind_Dir_cumul      = C_vec_t_study_s; // - Ind_Dir_anchor;
                            Ind_Dir_cumul_prob = Phi_approx(Ind_Dir_cumul);
                            Ind_Dir_ord_prob   = cumul_probs_to_ord_probs(Ind_Dir_cumul_prob);
                            ////
                            vector[n_thr_t] rho =  std_normal_approx_pdf(Ind_Dir_cumul, Ind_Dir_cumul_prob);
                            target += induced_dirichlet_given_rho_lpdf(Ind_Dir_ord_prob | rho, alpha_t);
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
          }
          ////
          //// ---- Log-likelihood:
          ////
          matrix[n_studies, n_thr_max] C_mat;
          ////
          for (t in 1:n_index_tests) {
                ////
                for (s in 1:n_studies) {
                         vector[n_thr[t]] C_vec_t_study_s = get_test_t_study_s_segment_values(C_vec, t, s, n_thr, n_studies);
                         C_mat[s, 1:n_thr[t]] = to_row_vector(C_vec_t_study_s);
                }
                ////
                array[2] matrix[n_studies, n_thr[t]] latent_surv_t = map_latent_surv_prob_to_random_C(
                                                                     C_mat, n_thr[t], Xlocations[t], Xscale[t], n_studies, n_obs_cutpoints[t, ], cutpoint_index[t, ]);
                ////
                target += compute_log_lik_binomial_fact_lp(
                          latent_surv_t, use_probit_link, n_thr[t], x_2[t, ], n[t, ], N_total[t, ], n_obs_cutpoints[t, ], indicator_index_test_in_study[, t]);
          }
         
}


generated quantities {
      ////
      //// ---- Calculate summary cutpoints:
      ////
      array[n_index_tests] vector[n_thr_max] C_mu;
      for (t in 1:n_index_tests) {
          for (k in 1:n_thr[t]) { 
              // Collect cutpoints across studies for this test and threshold
              vector[n_studies] C_vals;
              for (s in 1:n_studies) {
                  C_vals[s] = get_test_t_study_s_segment_values(C_vec, t, s, n_thr, n_studies)[k];
              }
              C_mu[t][k] = median(C_vals); // median cutpoint across the n_studies cutpoints
          }
      }
      ////
      //// ---- Compute Induced-Dirichlet "SD" params:
      ////
      array[n_index_tests] vector[n_cat_max] dirichlet_cat_SDs_sigma;
      array[n_index_tests] real alpha_0;
      for (t in 1:n_index_tests) {
          int n_cat_t = n_cat[t];
          int start_alpha_index = ((t == 1) ? 1 : (1 + sum(n_cat[1:(t - 1)])));
          vector[n_cat_t] alpha_t = segment(alpha_vec, start_alpha_index, n_cat_t);
          ////
          alpha_0[t] = sum(alpha_t);
          dirichlet_cat_SDs_sigma[t][1:n_cat_t] = sqrt((alpha_t .* (alpha_0[t] - alpha_t) ) ./ (square(alpha_0[t]) * (alpha_0[t] + 1.0)));
      }
      ////
      //// ---- Calculate summary accuracy (using mean parameters):
      ////
      array[n_index_tests] vector[n_covariates_max] location_nd_mu;
      array[n_index_tests] vector[n_covariates_max] location_d_mu;
      array[n_index_tests] vector[n_covariates_max] scale_nd_mu;
      array[n_index_tests] vector[n_covariates_max] scale_d_mu;
      ////
      for (t in 1:n_index_tests) {
          location_nd_mu[t] = mult_nd * beta_mu[t];
          location_d_mu[t]  = mult_d * beta_mu[t];
          ////
          scale_nd_mu[t] = ((softplus == 1) ? softplus_scaled(mult_nd * raw_scale_mu[t]) : exp(mult_nd * raw_scale_mu[t]));
          scale_d_mu[t]  = ((softplus == 1) ? softplus_scaled(mult_d * raw_scale_mu[t])  : exp(mult_d * raw_scale_mu[t]));
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
          //// ---- Calculate for actual thresholds
          for (k in 1:n_thr[t]) {
              Fp_baseline[t][k] = (use_probit_link == 1) ? 
                                  Phi(-(C_mu[t][k] - Xbeta_baseline_nd[t])/scale_nd_baseline[t]) : 
                                  inv_logit(-(C_mu[t][k] - Xbeta_baseline_nd[t])/scale_nd_baseline[t]);
              Sp_baseline[t][k] = 1.0 - Fp_baseline[t][k];
              Se_baseline[t][k] = (use_probit_link == 1) ? 
                                  Phi(-(C_mu[t][k] - Xbeta_baseline_d[t])/scale_d_baseline[t]) : 
                                  inv_logit(-(C_mu[t][k] - Xbeta_baseline_d[t])/scale_d_baseline[t]);
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
          //// ---- Draw predictive cutpoints for this test
          int n_cat_t = n_cat[t];
          int start_alpha_index = ((t == 1) ? 1 : (1 + sum(n_cat[1:(t - 1)])));
          vector[n_cat_t] alpha_t = segment(alpha_vec, start_alpha_index, n_cat_t);
          ////
          vector[n_cat_t] prob_ord_pred   = dirichlet_rng(alpha_t); // Simulate from Dirichlet
          vector[n_thr[t]] prob_cumul_pred = ord_probs_to_cumul_probs(prob_ord_pred); // Compute PREDICTED cumulative probabilities
          vector[n_thr[t]] C_pred_t = ID_cumul_probs_to_C(prob_cumul_pred, use_probit_link); // Compute PREDICTED cutpoints
          ////
          //// ---- Combine to get predictive values
          real Xbeta_baseline_pred = dot_product(baseline_case[t][1:n_covariates[t]], beta_mu[t][1:n_covariates[t]]) + beta_eta_pred + beta_delta_pred;
          real Xbeta_baseline_pred_nd = mult_nd * Xbeta_baseline_pred;
          real Xbeta_baseline_pred_d  = mult_d * Xbeta_baseline_pred;
          ////
          real Xraw_scale_baseline_pred = dot_product(baseline_case[t][1:n_covariates[t]], raw_scale_mu[t][1:n_covariates[t]]) + raw_scale_eta_pred + raw_scale_delta_pred;
          real scale_baseline_pred_nd = ((softplus == 1) ? softplus_scaled(mult_nd * Xraw_scale_baseline_pred) : exp(mult_nd * Xraw_scale_baseline_pred));
          real scale_baseline_pred_d  = ((softplus == 1) ? softplus_scaled(mult_d * Xraw_scale_baseline_pred)  : exp(mult_d * Xraw_scale_baseline_pred));
          ////
          //// ---- Calculate predictive Se/Sp
          for (k in 1:n_thr[t]) {
              Fp_baseline_pred[t][k] = (use_probit_link == 1) ? 
                                       Phi(-(C_pred_t[k] - Xbeta_baseline_pred_nd)/scale_baseline_pred_nd) : 
                                       inv_logit(-(C_pred_t[k] - Xbeta_baseline_pred_nd)/scale_baseline_pred_nd);
              Sp_baseline_pred[t][k] = 1.0 - Fp_baseline_pred[t][k];
              Se_baseline_pred[t][k] = (use_probit_link == 1) ? 
                                       Phi(-(C_pred_t[k] - Xbeta_baseline_pred_d)/scale_baseline_pred_d) : 
                                       inv_logit(-(C_pred_t[k] - Xbeta_baseline_pred_d)/scale_baseline_pred_d);
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
              //// ---- Extract cutpoints for all studies for this test
              array[n_studies] vector[n_thr[t]] C_vec_t_all_studies;
              for (s in 1:n_studies) {
                  C_vec_t_all_studies[s] = get_test_t_study_s_segment_values(C_vec, t, s, n_thr, n_studies);
              }
              ////
              //// ---- Convert to matrix format expected by map_latent_surv_prob_to_random_C
              matrix[n_studies, n_thr[t]] C_mat_t;
              for (s in 1:n_studies) {
                  for (k in 1:n_thr[t]) { 
                      C_mat_t[s, k] = C_vec_t_all_studies[s][k];
                  }
              }
              ////
              array[2] matrix[n_studies, n_thr[t]] latent_surv_t = map_latent_surv_prob_to_random_C(
                                                                     C_mat_t, n_thr[t], Xlocations[t], Xscale[t], 
                                                                     n_studies, n_obs_cutpoints[t, ], cutpoint_index[t, ]);
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
                                                                         cond_prob, n_thr[t], x_2[t, ], n[t, ], 
                                                                         n_obs_cutpoints[t, ]);
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
}



// generated quantities {
// 
//           array[n_index_tests] vector[n_thr_max] Se;
//           array[n_index_tests] vector[n_thr_max] Sp;
//           array[n_index_tests] vector[n_thr_max] Fp;
//           ////
//           array[n_index_tests] vector[n_thr_max] Se_pred;
//           array[n_index_tests] vector[n_thr_max] Sp_pred;
//           array[n_index_tests] vector[n_thr_max] Fp_pred;
//           ////
//           array[n_index_tests] matrix[n_studies, n_thr_max] se;
//           array[n_index_tests] matrix[n_studies, n_thr_max] sp;
//           array[n_index_tests] matrix[n_studies, n_thr_max] fp;
//           ////
//           array[n_index_tests] matrix[n_studies, n_thr_max] x_hat_nd = init_array_of_matrices(n_studies, n_thr_max, n_index_tests, -1.0);
//           array[n_index_tests] matrix[n_studies, n_thr_max] x_hat_d  = init_array_of_matrices(n_studies, n_thr_max, n_index_tests, -1.0);
//           array[n_index_tests] matrix[n_studies, n_thr_max] dev_nd   = init_array_of_matrices(n_studies, n_thr_max, n_index_tests, -1.0);
//           array[n_index_tests] matrix[n_studies, n_thr_max] dev_d    = init_array_of_matrices(n_studies, n_thr_max, n_index_tests, -1.0);
//           ////
//           //// ---- Compute mean / pooled summary-level cutpoints (from Induced-Dirichlet between-study model):
//           ////
//           array[n_index_tests] vector[n_cat_max] prob_ord_mu;
//           array[n_index_tests] vector[n_thr_max] prob_cumul_mu;
//           array[n_index_tests] vector[n_thr_max] C_mu;
//           ////
//           for (t in 1:n_index_tests) {
//                 int n_thr_t = n_thr[t];
//                 int n_cat_t = n_cat[t];
//                 int start_alpha_index = ((t == 1) ? (1 + sum(n_cat[1:(t - 1)])) : 1);
//                 vector[n_cat_t] alpha_t = segment(alpha_vec, start_alpha_index, n_cat_t);
//                 ////
//                 prob_ord_mu[t][1:n_cat_t]   = alpha_t / sum(alpha_t);  //// dirichlet_cat_means_phi[c];
//                 prob_cumul_mu[t][1:n_thr_t] = ord_probs_to_cumul_probs(prob_ord_mu[t][1:n_cat_t]);
//                 C_mu[t][1:n_thr_t]          = cumul_probs_to_C(prob_cumul_mu[t][1:n_thr_t], 0.0, 1.0);
//           }
//           ////
//           //// ---- Compute Induced-Dirichlet "SD" params:
//           ////
//           array[n_index_tests] vector[n_cat_max] dirichlet_cat_SDs_sigma;
//           vector[n_index_tests] alpha_0;
//           for (t in 1:n_index_tests) {
//                 int n_cat_t = n_cat[t];
//                 int start_alpha_index = ((t == 1) ? (1 + sum(n_cat[1:(t - 1)])) : 1);
//                 vector[n_cat_t] alpha_t = segment(alpha_vec, start_alpha_index, n_cat_t);
//                 ////
//                 // vector[n_cat_t] alpha_t = get_test_values(alpha, start_index_pooled, end_index_pooled, t);
//                 alpha_0[t] = sum(alpha_t);
//                 dirichlet_cat_SDs_sigma[t][1:n_cat_t] = sqrt((alpha_t .* (alpha_0[t] - alpha_t) ) ./ (square(alpha_0[t]) * (alpha_0[t] + 1.0)));
//           }
//           ////
//           //// ---- Calculate study-specific accuracy:
//           ////
//           for (t in 1:n_index_tests) {
//                   int n_thr_t = n_thr[t];
//                   fp[t][1:n_studies, 1:n_thr_t] = 1.0 - cumul_prob[t, 1][1:n_studies, 1:n_thr_t];
//                   sp[t][1:n_studies, 1:n_thr_t] = 1.0 - fp[t][1:n_studies, 1:n_thr_t];
//                   se[t][1:n_studies, 1:n_thr_t] = 1.0 - cumul_prob[t, 2][1:n_studies, 1:n_thr_t];
//           }
//           ////
//           //// ---- Calculate summary accuracy (using mean parameters):
//           ////
//           for (t in 1:n_index_tests) {
//                   int n_thr_t = n_thr[t];
//                   vector[n_thr_t] C_mu_t = C_mu[t][1:n_thr_t];
//                   // vector[n_thr_t] C_vec_t = get_test_values(C_vec, start_index_study_s, end_index_study_s, t);
//                   Fp[t][1:n_thr_t] = 1.0 - Phi_approx((C_mu_t - beta_mu[t, 1]));
//                   Sp[t][1:n_thr_t] = 1.0 - Fp[t][1:n_thr_t];
//                   Se[t][1:n_thr_t] = 1.0 - Phi_approx((C_mu_t - beta_mu[t, 2]));
//           }
//           ////
//           //// ---- Calculate predictive accuracy:
//           ////
//           {
//             vector[2] beta_eta_pred      = to_vector(normal_rng(rep_vector(0.0, 2), beta_sigma[1:2]));  //// shared between tests
//             // ////
//             // array[n_index_tests] vector[n_cat_max] prob_ord_pred   = dirichlet_rng(alpha[t][1:n_thr_t]); //// Simulate from Dirichlet by using the summary "alpha" parameters.
//             // array[n_index_tests] vector[n_thr_max] prob_cumul_pred = ord_probs_to_cumul_probs(prob_ord_pred);  //// Compute PREDICTED cumulative probabilities.
//             // array[n_index_tests] vector[n_thr_max] C_pred = cumul_probs_to_C(prob_cumul_pred, 0.0, 1.0);  //// Compute PREDICTED cutpoints.
//             ////
//             for (t in 1:n_index_tests) {
//               
//                   int n_thr_t = n_thr[t];
//                   int n_cat_t = n_cat[t];
//                   int start_alpha_index = ((t == 1) ? (1 + sum(n_cat[1:(t - 1)])) : 1);
//                   vector[n_cat_t] alpha_t = segment(alpha_vec, start_alpha_index, n_cat_t);
//                   
//                   vector[n_cat_t] prob_ord_pred_t   = dirichlet_rng(alpha_t); //// Simulate from Dirichlet by using the summary "alpha" parameters.
//                   vector[n_thr_t] prob_cumul_pred_t = ord_probs_to_cumul_probs(prob_ord_pred_t);  //// Compute PREDICTED cumulative probabilities.
//                   vector[n_thr_t] C_pred_t          = cumul_probs_to_C(prob_cumul_pred_t, 0.0, 1.0);  //// Compute PREDICTED cutpoints.
//                   ////
//                   vector[2] beta_delta_t_pred = to_vector(normal_rng(rep_vector(0.0, 2), beta_tau[t, 1:2]));
//                   vector[2] beta_t_pred       = to_vector(beta_mu[t, 1:2]) + beta_eta_pred[1:2] + beta_delta_t_pred[1:2];
//                   ////
//                   Fp_pred[t][1:n_thr_t] = 1.0 - Phi_approx(C_pred_t - beta_t_pred[1]);
//                   Sp_pred[t][1:n_thr_t] = 1.0 - Fp_pred[t][1:n_thr_t];
//                   Se_pred[t][1:n_thr_t] = 1.0 - Phi_approx(C_pred_t - beta_t_pred[2]);
//             }
//           }
//           ////
//           //// ---- Model-predicted ("re-constructed") data:
//           ////
//           {
//                   array[n_index_tests, 2] matrix[n_studies, n_thr_max] x_hat = init_nested_array_of_matrices(n_studies, n_thr_max, n_index_tests, 2, -1.0);
//                   array[n_index_tests, 2] matrix[n_studies, n_thr_max] dev   = init_nested_array_of_matrices(n_studies, n_thr_max, n_index_tests, 2, -1.0);
//                   ////
//                   for (t in 1:n_index_tests) {
//                       for (s in 1:n_studies) {
//                             if (indicator_index_test_in_study[s, t] == 1) {
//                                       for (c in 1:2) {
//                                          for (cut_i in 1:to_int(n_obs_cutpoints[s, t])) {
//         
//                                                   //// Model-estimated data:
//                                                   x_hat[t, c][s, cut_i] = cond_prob[t, c][s, cut_i] * n[t, c][s, cut_i];  	 //// Fitted values
//         
//                                                   //// Compute residual deviance contribution:
//                                                   real n_i =  (n[t, c][s, cut_i]);
//                                                   real x_i =  (x[t, c][s, cut_i]);
//                                                   real x_hat_i =  (x_hat[t, c][s, cut_i]);
//                                                   real log_x_minus_log_x_hat = log(x_i) - log(x_hat_i);
//                                                   real log_diff_n_minus_x = log(n_i - x_i);
//                                                   real log_diff_n_minus_x_hat = log(abs(n_i - x_hat_i));
//         
//                                                   // array[n_index_tests, 2] matrix[n_studies, n_thr_max] n;
//         
//                                                   dev[t, c][s, cut_i] = 2.0 * ( x_i * log_x_minus_log_x_hat + (n_i - x_i) * (log_diff_n_minus_x - log_diff_n_minus_x_hat) );
//         
//                                          }
//                                       }
//                             }
//                       }
//                   }
//                   ////
//                   //// ---- Store deviance and x_hat split by disease status:
//                   ////
//                   for (t in 1:n_index_tests) {
//                        x_hat_nd[t] = x_hat[t, 1];
//                        dev_nd[t]   = dev[t, 1];
//                        x_hat_d[t]  = x_hat[t, 2];
//                        dev_d[t]    = dev[t, 2];
//                   }
//           }
//           ////
//           //// ---- NMA: Compute between-study heterogeneity + correlations for "beta":
//           ////
//           matrix[n_index_tests, 2] beta_tau_sq = square(beta_tau);
//           array[2] matrix[n_index_tests, n_index_tests] beta_sigma_sq;
//           array[2] matrix[n_index_tests, n_index_tests] beta_rho;
//           array[2] matrix[n_index_tests, n_index_tests] beta_rho12;
//           cov_matrix[2] beta_Sigma = multiply_lower_tri_self_transpose(beta_L_Sigma);
//           {
//             array[3, 2] matrix[n_tests, n_tests] NMA_Nyaga_outs_beta = compute_Nyaga_NMA_summaries(n_tests, beta_tau_sq, beta_Sigma, beta_Omega);
//             beta_sigma_sq = NMA_Nyaga_outs_beta[1];
//             beta_rho      = NMA_Nyaga_outs_beta[2];
//             beta_rho12    = NMA_Nyaga_outs_beta[3];
//           }
//           ////
//           //// ---- NMA: Compute "rank statistics" + other summary (non-comparative) estimates:
//           ////
//           array[n_index_tests] vector[n_thr_max] LRpos;  
//           array[n_index_tests] vector[n_thr_max] LRneg; 
//           array[n_index_tests] vector[n_thr_max] DOR;  
//           array[n_index_tests] vector[n_thr_max] Youden_index; 
//           array[n_index_tests] vector[n_thr_max] Youden_index_weighted; 
//           for (t in 1:n_index_tests) {
//             
//                     int n_thr_t = n_thr[t];
//                     array[n_thr_t] int n_thr_t_index = linspaced_int_array(n_thr_t, 1, n_thr_t);
//                     vector[n_thr_t] Se_test_t = Se[t][n_thr_t_index];
//                     vector[n_thr_t] Sp_test_t = Sp[t][n_thr_t_index];
//                     vector[n_thr_t] Fp_test_t = Fp[t][n_thr_t_index];
//                     ////
//                     //// ---- Youden index:
//                     ////
//                     Youden_index[t][n_thr_t_index] = Se_test_t + Sp_test_t - 1.0;
//                     ////
//                     //// ---- Weighted Youden index (weight = 0.5 -> same as youden, weight > 0.5 -> more importance on Se):
//                     ////
//                     real youden_weight = 0.50;
//                     Youden_index_weighted[t][n_thr_t_index] = 2.0 * (youden_weight .* Se_test_t + (1.0 - youden_weight) .* Sp_test_t);
//                     ////
//                     //// ---- DOR:
//                     ////
//                   	DOR[t][n_thr_t_index] = (Se_test_t .* Sp_test_t) ./ ((1.0 - Se_test_t) .* (1.0 - Sp_test_t));
//                     ////
//                     //// ---- "Likelihood Ratio's" (LR's):
//                     ////
//                    	LRpos[t][n_thr_t_index] = Se_test_t ./ Sp_test_t;
//     	              LRneg[t][n_thr_t_index] = (1.0 - Se_test_t) ./ Sp_test_t;
//     	              // ////
//     	              // //// PPV and NPV (using user-inputted prev data vector)
//     	              // ////
//     	              // PPV[t][n_thr_t_index] = Se[t] * prev[Test[n]] / ( Se[t]*prev[Test[n]] + (1-Sp[t])*(1.0-prev[Test[n]]) );
//     	              // NPV[t][n_thr_t_index] = Sp[t] * prev[Test[n]] / ( (1.0-Se[t])*prev[Test[n]] + Sp[t]*(1.0-prev[Test[n]]) );
//     	              
//           }
//           ////
//           //// ---- NMA: Compute pairwise accuracy differences + ratios:
//           ////
//           array[n_index_tests, n_index_tests] matrix[n_thr_max, n_thr_max] diff_Se;
//           array[n_index_tests, n_index_tests] matrix[n_thr_max, n_thr_max] diff_Sp;
//           array[n_index_tests, n_index_tests] matrix[n_thr_max, n_thr_max] ratio_Se;
//           array[n_index_tests, n_index_tests] matrix[n_thr_max, n_thr_max] ratio_Sp;
//           for (t1 in 1:(n_tests - 1)) {
//               
//                     int n_thr_t1 = n_thr[t1];
//                     array[n_thr_t1] int n_thr_t1_index = linspaced_int_array(n_thr_t1, 1, n_thr_t1);
//                           
//                     for (t2 in (t1 + 1):n_tests) {
//                       
//                             int n_thr_t2 = n_thr[t2];
//                             array[n_thr_t2] int n_thr_t2_index = linspaced_int_array(n_thr_t2, 1, n_thr_t2);
//                             ////
//                             //// Get vectors of accuracies for all thresholds:
//                             ////
//                             vector[n_thr_t1] Se_test_t1 = Se[t1][n_thr_t1_index];
//                             vector[n_thr_t2] Se_test_t2 = Se[t2][n_thr_t2_index];
//                             vector[n_thr_t1] Sp_test_t1 = Sp[t1][n_thr_t1_index];
//                             vector[n_thr_t2] Sp_test_t2 = Sp[t2][n_thr_t2_index];
//                             ////
//                             //// Compute all pairwise differences:
//                             ////
//                             diff_Se[t1, t2][n_thr_t1_index, n_thr_t2_index] = compute_between_test_diffs( Se_test_t1, Se_test_t2, n_thr_t1, n_thr_t2);
//                             diff_Sp[t1, t2][n_thr_t1_index, n_thr_t2_index] = compute_between_test_diffs( Sp_test_t1, Sp_test_t2, n_thr_t1, n_thr_t2);
//                             ////
//                             //// Compute all pairwise ratios:
//                             ////
//                             ratio_Se[t1, t2][n_thr_t1_index, n_thr_t2_index] = compute_between_test_ratios(Se_test_t1, Se_test_t2, n_thr_t1, n_thr_t2);
//                             ratio_Sp[t1, t2][n_thr_t1_index, n_thr_t2_index] = compute_between_test_ratios(Sp_test_t1, Sp_test_t2, n_thr_t1, n_thr_t2);
//                       
//                     }
//                   
//           }
// 
// 
// }
// 
