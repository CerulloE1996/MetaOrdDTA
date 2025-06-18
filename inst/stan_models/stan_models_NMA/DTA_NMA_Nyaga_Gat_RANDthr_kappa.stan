

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
          //// ---- baseline_case_nd/baseline_case_d must be user-inputted - 
          //// ---- e.g. could be {0, 1, 45.3} for 2 binary covs and 1 cts one (e.g. age)
          ////
          array[n_index_tests] int n_covariates; //// For R+G - cannot have SEPERATE covariates for D+ and D-, so only one n_covariates!
          int n_covariates_max; // this is the max across ALL index tests
          array[n_index_tests] matrix[n_studies, max(n_covariates)] X;
          array[n_index_tests] vector[max(n_covariates)] baseline_case; 
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
          array[n_index_tests] vector<lower=0.0>[n_thr_max + 1] prior_dirichlet_phi;
          array[n_index_tests] real prior_kappa_mean;
          array[n_index_tests] real<lower=0.0> prior_kappa_SD;
          ////
          //// ---- Other:
          ////
          real<lower=0.0> kappa_lb;
          int<lower=0, upper=1> softplus;
          int n_total_C_if_random;
          ////
          int use_probit_link;
}

transformed data {
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
          };
          ////
          //// ---- HSROC stuff:
          ////
          real mult_nd = -1.0;
          real mult_d  = +1.0; //// hence: mult_d > mult_nd
          ////
          int n_total_C_per_class_if_fixed = sum(n_thr); // total # of fixed-effect cutpoints
          ////
          //// ---- Calculate indices:
          ////
          array[n_index_tests] int start_index = calculate_start_indices(n_thr, n_index_tests);
          array[n_index_tests] int end_index   = calculate_end_indices(n_thr, n_index_tests, start_index);
}


parameters {
          array[n_index_tests] vector[n_covariates_max] beta_mu;
          array[n_index_tests] vector[n_covariates_max] raw_scale_mu;
          ////
          //// ---- "NMA" params for SHARED "beta":
          ////
          real<lower=0.0> beta_sigma;              //// Variance component - Between-study SD (Nyaga's σ)
          vector<lower=0.0>[n_index_tests] beta_tau; //// Variance component - Test-specific SD (Nyaga's τ) - delta_{s,c,t} ~ normal(0, tau_{c, t}).
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
          //// ---- Induced-Dirichlet params:
          ////
          vector<lower=-5.0, upper=3.0>[n_total_C_if_random] C_raw_vec;  //// RAW LOG-DIFFERENCES -
          vector<lower=-5.0, upper=5.0>[n_total_C_per_class_if_fixed] raw_dirichlet_phi_vec; // props/probs
          array[n_index_tests] real<lower=1.0> kappa;
          
}


transformed parameters {
          matrix[n_index_tests, n_cat_max] dirichlet_phi; // props/probs
          for (t in 1:n_index_tests) {
                  vector[n_thr[t]] raw_vec_simplex_t = get_test_values(
                                                       raw_dirichlet_phi_vec, start_index, end_index, t);
                  vector[n_cat[t]] vec_simplex_t = stickbreaking_logistic_simplex_constrain_jacobian(raw_vec_simplex_t);
                  dirichlet_phi[t, 1:n_cat[t]] = to_row_vector(vec_simplex_t);
          }
          ////
          matrix[n_index_tests, n_cat_max] alpha;
          matrix[n_index_tests, n_cat_max] dir_cat_SDs_sigma;
          ////
          //// ---- Calculate alpha parameters from target SDs for each test:
          ////
          for (t in 1:n_index_tests) {
              ////
              alpha[t, 1:n_cat[t]] = 0.01 + dirichlet_phi[t, 1:n_cat[t]] * kappa[t]; //// Set alpha using the approximated kappa:
              ////
              real alpha_0 = sum(alpha[t, 1:n_cat[t]]);
              dir_cat_SDs_sigma[t, 1:n_cat[t]] = sqrt((alpha[t, 1:n_cat[t]] .* (alpha_0 - alpha[t, 1:n_cat[t]])) / 
                                                 (square(alpha_0) * (alpha_0 + 1.0)));
          }
          ////
          //// ---- Construct study-specific cutpoints:
          ////
          vector[n_total_C_if_random] C_vec;  //// Global cutpoints for each test ("staggered" array/matrix using "n_thr[t]" to index correctly)
          for (t in 1:n_index_tests) {
                int n_thr_t = n_thr[t];
                ////
                for (s in 1:n_studies) {
                      vector[n_thr_t] C_raw_vec_test_t = get_test_t_study_s_segment_values(
                                                         C_raw_vec, t, s, n_thr, n_studies);
                      vector[n_thr_t] C_vec_test_t = construct_C(C_raw_vec_test_t, softplus);
                      C_vec = update_test_t_study_s_segment_values(
                              C_vec, C_vec_test_t, t, s, n_thr, n_studies);
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
              beta_eta[s]      = beta_sigma      * beta_eta_z[s];       //// beta_eta[s] ~ normal(0.0, sigma);
              raw_scale_eta[s] = raw_scale_sigma * raw_scale_eta_z[s];  //// raw_scale_eta[s] ~ normal(0.0, sigma);
          }
          //// Compute test-specific deviations ("delta" in Nyaga notation):
          for (t in 1:n_index_tests) { 
             beta_delta[t]      = beta_tau[t]      * beta_delta_z[t]; //// beta_delta[t][s] ~ normal(0.0, beta_tau[t]);
             raw_scale_delta[t] = raw_scale_tau[t] * raw_scale_delta_z[t];
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
          //// ---- Linear predictors for each disease status + apply covariates for non-diseased and diseased groups:
          ////
          array[n_index_tests] matrix[n_studies, 2] Xlocations; // NOT local
          array[n_index_tests] matrix[n_studies, 2] Xraw_scale; // NOT local
          array[n_index_tests] matrix[n_studies, 2] Xscale;     // NOT local
          ////
          for (t in 1:n_index_tests) {
                // for (s in 1:n_studies) {
                    vector[n_studies] Xbeta_baseline  = X[t][, 1:n_covariates[t]] * to_vector(beta_mu[t][1:n_covariates[t]]) + 
                                                        beta_random[t];
                    Xlocations[t][, 1] = mult_nd*Xbeta_baseline; 
                    Xlocations[t][, 2] = mult_d*Xbeta_baseline; 
                    ////
                    vector[n_studies] Xraw_scale_baseline = X[t][, 1:n_covariates[t]] * to_vector(raw_scale_mu[t][1:n_covariates[t]]) + 
                                                            raw_scale_random[t]; // "Xgamma"
                    Xscale[t][, 1] =  ((softplus == 1) ? 
                                       softplus_scaled( mult_nd*Xraw_scale_baseline) : 
                                       exp( mult_nd*Xraw_scale_baseline)); 
                    Xscale[t][, 2] =  ((softplus == 1) ? 
                                       softplus_scaled( mult_d*Xraw_scale_baseline)  :
                                       exp( mult_d*Xraw_scale_baseline));
                // }
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
          //// ---- Priors for Induced-Dirichlet params for each test:
          ////
          for (t in 1:n_index_tests) {
              kappa[t] ~ lognormal(prior_kappa_mean[t], prior_kappa_SD[t]);
              dirichlet_phi[t, 1:n_cat[t]] ~ dirichlet(prior_dirichlet_phi[t, 1:n_cat[t]]); // ID prior
          }
          ////
          //// ---- Induced-dirichlet between study ** model ** (NOT a prior model here but part of the actual
          //// ---- likelihood since random-effect cutpoints!):
          ////
          for (t in 1:n_index_tests) {
                for (s in 1:n_studies) {
                        vector[n_thr[t]] C_vec_t = get_test_t_study_s_segment_values(C_vec, t, s, n_thr, n_studies);
                        ////
                        vector[n_cat[t]] Ind_Dir_ord_prob = (use_probit_link == 1) ?
                                                           cumul_probs_to_ord_probs(Phi(C_vec_t)) : 
                                                           cumul_probs_to_ord_probs(inv_logit(C_vec_t));
                        ////
                        Ind_Dir_ord_prob ~ induced_dirichlet_given_C(
                                           C_vec_t, to_vector(alpha[t][1:n_cat[t]]), use_probit_link);
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
          for (t in 1:n_index_tests) {
                for (s in 1:n_studies) {
                    vector[n_thr[t]] C_raw_vec_test_t = get_test_t_study_s_segment_values(C_raw_vec, t, s, n_thr, n_studies);
                    target += raw_C_to_C_log_det_J_lp(C_raw_vec_test_t, softplus);
                }
          }
          ////
          //// ---- Log-likelihood:
          ////
          matrix[n_studies, n_thr_max] C_mat;
          for (t in 1:n_index_tests) {
               for (c in 1:2) {
                    for (s in 1:n_studies) {
                             vector[n_thr[t]] C_vec_t_study_s = get_test_t_study_s_segment_values(C_vec, t, s, n_thr, n_studies);
                             C_mat[s, 1:n_thr[t]] = to_row_vector(C_vec_t_study_s);
                    }
                }
                // ////
                // vector[n_thr[t]] C_vec_t;
                // {
                //     vector[n_thr[t]] C_vec_t_given_c = get_test_values(C_vec, start_index, end_index, t);
                //     C_vec_t = C_vec_t_given_c;
                // }
                // ////
                array[2] matrix[n_studies, n_thr[t]] latent_surv_t = map_latent_surv_prob_to_random_C(
                                                                     C_mat, n_thr[t], Xlocations[t], Xscale[t], n_studies, 
                                                                     n_obs_cutpoints[t, ], cutpoint_index[t, ]);
                ////
                target += compute_log_lik_binomial_fact_lp(
                          latent_surv_t, use_probit_link, n_thr[t], x_2[t, ], n[t, ], 
                          N_total[t, ], n_obs_cutpoints[t, ], indicator_index_test_in_study[, t]);
          }
}


generated quantities {
      ////
      //// ---- SD's for each category probability in a Dirichlet is √(α_k(α_0-α_k)/(α_0^2(α_0+1))) - where α_0 is the sum of all alphas:
      ////
      array[n_index_tests] vector[n_cat_max] dirichlet_cat_SDs_sigma_full;
      for (t in 1:n_index_tests) {
          real alpha_0 = sum(alpha[t][1:n_cat[t]]);
          dirichlet_cat_SDs_sigma_full[t][1:n_cat[t]] = to_vector( sqrt((alpha[t][1:n_cat[t]] .* (alpha_0 - alpha[t][1:n_cat[t]]) ) /
                                                                   (square(alpha_0) * (alpha_0 + 1.0))));
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
      //// ---- Calculate summary cutpoints:
      ////
      array[n_index_tests] vector[n_thr_max] C_mu;
      for (t in 1:n_index_tests) {
          for (k in 1:n_thr[t]) {
              //// ---- Collect cutpoints across studies for this test and threshold
              vector[n_studies] C_vals;
              for (s in 1:n_studies) {
                  C_vals[s] = get_test_t_study_s_segment_values(C_vec, t, s, n_thr, n_studies)[k];
              }
              //// ---- median cutpoint across the n_studies cutpoints
              C_mu[t][k] = median(C_vals); // solves "skewness" issue from non-linear transformations
          }
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
          scale_nd_baseline[t] = ((softplus == 1) ? 
                                 softplus_scaled(mult_nd * Xraw_scale_baseline[t]) : 
                                 exp(mult_nd * Xraw_scale_baseline[t]));
          scale_d_baseline[t]  = ((softplus == 1) ? 
                                 softplus_scaled(mult_d * Xraw_scale_baseline[t])  : 
                                 exp(mult_d * Xraw_scale_baseline[t]));
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
      //// ---- Calculate predictive accuracy:
      ////
      array[n_index_tests] vector[n_thr_max] C_pred;
      ////
      for (t in 1:n_index_tests) {
        for (k in 1:n_thr[t]) {
            real sigma_C_k;
            real SD_prob_scale = dirichlet_cat_SDs_sigma_full[t][k];
            if (use_probit_link == 1) {
                real mean_C_cutpoint_scale = C_mu[t][k];
                sigma_C_k = SD_approx_ID_ord_prob_to_C_probit(mean_C_cutpoint_scale, SD_prob_scale);
            } else { 
                real mean_prob_scale = dirichlet_phi[t][k];
                sigma_C_k = SD_approx_ID_ord_prob_to_C_logit(mean_prob_scale, SD_prob_scale);
            }
            C_pred[t][k] = C_mu[t][k] + normal_rng(0, sigma_C_k);
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
                                        softplus_scaled(mult_nd * Xraw_scale_baseline_pred) : 
                                        exp(mult_nd * Xraw_scale_baseline_pred));
          real scale_baseline_pred_d  = ((softplus == 1) ? 
                                        softplus_scaled(mult_d * Xraw_scale_baseline_pred)  : 
                                        exp(mult_d * Xraw_scale_baseline_pred));
          ////
          //// ---- Calculate predictive Se/Sp
          for (k in 1:n_thr[t]) {
              Fp_baseline_pred[t][k] = (use_probit_link == 1) ? 
                                       Phi(-(C_pred[t][k] - Xbeta_baseline_pred_nd)/scale_baseline_pred_nd) : 
                                       inv_logit(-(C_pred[t][k] - Xbeta_baseline_pred_nd)/scale_baseline_pred_nd);
              Sp_baseline_pred[t][k] = 1.0 - Fp_baseline_pred[t][k];
              Se_baseline_pred[t][k] = (use_probit_link == 1) ? 
                                       Phi(-(C_pred[t][k] - Xbeta_baseline_pred_d)/scale_baseline_pred_d) : 
                                       inv_logit(-(C_pred[t][k] - Xbeta_baseline_pred_d)/scale_baseline_pred_d);
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
          matrix[n_studies, n_thr_max] C_mat;
          for (t in 1:n_index_tests) {
              for (s in 1:n_studies) {
                   vector[n_thr[t]] C_vec_t_study_s = get_test_t_study_s_segment_values(C_vec, t, s, n_thr, n_studies);
                   C_mat[s, 1:n_thr[t]] = to_row_vector(C_vec_t_study_s);
              }
              vector[n_thr[t]] C_vec_t = get_test_values(C_vec, start_index, end_index, t);
              ////
              array[2] matrix[n_studies, n_thr[t]] latent_surv_t = map_latent_surv_prob_to_random_C(
                                                                   C_mat, n_thr[t], Xlocations[t], Xscale[t], n_studies, 
                                                                   n_obs_cutpoints[t, ], cutpoint_index[t, ]);
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
                              deviance_nd[t][s] += is_nan(dev_nd[s, k]) ? 0.0 : dev_nd[s, k];
                              deviance_d[t][s]  += is_nan(dev_d[s, k]) ? 0.0 : dev_d[s, k];
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
          for (x_i in 1:n_covariates[t]) {
              for (k in 1:n_thr[t]) {
                  biv_equiv_C[t, 1][k, x_i] = C_mu[t][k] / scale_nd_mu[t][x_i];
                  biv_equiv_C[t, 2][k, x_i] = C_mu[t][k] / scale_d_mu[t][x_i];
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
      array[n_index_tests] real total_var_beta;
      array[n_index_tests] real total_var_raw_scale;
      ////
      for (t in 1:n_index_tests) {
          total_var_beta[t] = square(beta_sigma) + square(beta_tau[t]);
          total_var_raw_scale[t] = square(raw_scale_sigma) + square(raw_scale_tau[t]);
      }
      ////
      //// Proportion of variance explained by test-specific component
      array[n_index_tests] real prop_test_var_beta;
      array[n_index_tests] real prop_test_var_raw_scale;
      ////
      for (t in 1:n_index_tests) {
          prop_test_var_beta[t] = square(beta_tau[t]) / total_var_beta[t];
          prop_test_var_raw_scale[t] = square(raw_scale_tau[t]) / total_var_raw_scale[t];
      }
}

