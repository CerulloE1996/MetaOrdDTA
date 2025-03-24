

functions {
        ////
        //// Include files to compile the necessary custom (user-defined) Stan functions:
        ////
        #include "Stan_fns_basic.stan"
        #include "Stan_fns_corr.stan"
        #include "Stan_fns_ordinal_and_cutpoints.stan"
        #include "Stan_fns_log_lik.stan"
        #include "Stan_fns_ragged.stan"
        #include "Stan_fns_NMA.stan"
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
          array[n_studies, n_index_tests] int n_obs_cutpoints; //// OBSERVED cutpoints for test t in study s
          ////
          array[n_studies, n_index_tests] int<lower=0, upper=1> indicator_index_test_in_study;   // Binary indicator if test t is in study s
          ////
          //// ---- Data:
          ////
          array[n_index_tests, 2] matrix[n_studies, n_thr_max] x;
          array[n_index_tests, 2] matrix[n_studies, n_thr_max] n;
          array[n_index_tests, 2] matrix[n_studies, n_thr_max] cutpoint_index;
          ////
          // matrix[n_studies, n_index_tests] obs_indicator; // BINARY indicator for whether test t is observed in study s
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
 
             int n_tests = n_index_tests;
             // int n_total_pooled_cat = sum(n_cat);
//           ////
//           array[n_index_tests] int<lower=1> n_thr_random;
//           // array[n_index_tests] int<lower=2> n_cat_random;
//           for (t in 1:n_index_tests) {
//              n_thr_random[t] = n_studies * n_thr[t];
//              // n_cat_random[t] = n_studies * n_cat[t];
//           }
//           int n_total_C_if_random   = sum(n_thr_random);
//           // // int n_total_cat_if_random = sum(n_cat_random);
//           // ////
//           // //// ---- Calculate indices for thresholds:
//           // ////
//           // array[n_tests] int start_thr_i = calculate_start_indices(n_thr, n_tests);
//           // array[n_tests] int end_thr_i   = calculate_end_indices(n_thr, n_tests, start_thr_i);
//           // ////
//           // // array[n_tests] int start_rand_thr_i = calculate_start_indices_study(n_studies, n_thr, n_tests);
//           // // array[n_tests] int end_thr_index   = calculate_end_indices_study(n_studies, n_thr, n_tests, start_thr_index);
//           // ////
//           // //// ---- Calculate indices for categories:
//           // ////
//           // array[n_tests] int start_cat_i = calculate_start_indices(n_cat, n_tests);
//           // array[n_tests] int end_cat_i   = calculate_end_indices(n_cat, n_tests, start_cat_i);
//           // ////
//           // // array[n_tests] int start_cat_index = calculate_start_indices_study(n_studies, n_cat, n_tests);
//           // // array[n_tests] int end_cat_index   = calculate_end_indices_study(n_studies, n_cat, n_tests, start_cat_index);
      
}


parameters {
  
          matrix[n_index_tests, 2] beta_mu;
          ////
          //// ---- "NMA" params:
          ////
          vector<lower=0.0>[2] beta_sigma;              //// Variance component - Between-study SD (Nyaga's σ)
          matrix<lower=0.0>[n_index_tests, 2] beta_tau; //// Variance component - Test-specific SD (Nyaga's τ) - delta_{s, c, t} ~ normal(0, tau_{c, t}).
          ////
          matrix[n_studies, 2] beta_eta_z;              //// Standard normal RVs for study-level effects - eta[s, 1:2] ~ multi_normal({0, 0}, Sigma).
          array[n_index_tests] matrix[n_studies, 2] beta_delta_z; //// Standard normal RVs for test-specific effects
          ////
          real beta_corr;  //// between-study corr (possibly restricted)
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
                            C_vec_test_t = ((softplus == 1) ? construct_C_using_SP_jacobian(C_raw_vec_test_t) : construct_C_using_exp_jacobian(C_raw_vec_test_t));
                            C_vec = update_test_t_study_s_segment_values(C_vec, C_vec_test_t, t, s, n_thr, n_studies);
                  }
          }
          ////
          //// ---- Construct simple 2x2 (bivariate) between-study corr matrices:
          ////
          cholesky_factor_corr[2] beta_L_Omega = make_restricted_bivariate_L_Omega_jacobian(beta_corr, beta_corr_lb, beta_corr_ub);
          corr_matrix[2] beta_Omega      = multiply_lower_tri_self_transpose(beta_L_Omega);
          cholesky_factor_cov[2]  beta_L_Sigma      = diag_pre_multiply(beta_sigma, beta_L_Omega);
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
          //// ---- Log-likelihood stuff:
          ////
          array[n_index_tests, 2] matrix[n_studies, n_thr_max] cumul_prob = init_nested_array_of_matrices(n_studies, n_thr_max, n_index_tests, 2, 1.0);
          array[n_index_tests, 2] matrix[n_studies, n_thr_max] cond_prob  = init_nested_array_of_matrices(n_studies, n_thr_max, n_index_tests, 2, 1.0);
          array[n_index_tests, 2] matrix[n_studies, n_thr_max] log_lik    = init_nested_array_of_matrices(n_studies, n_thr_max, n_index_tests, 2, 0.0);
          {
                    array[n_index_tests, 2] matrix[n_studies, n_thr_max] latent_cumul_prob = init_nested_array_of_matrices(n_studies, n_thr_max, n_index_tests, 2, positive_infinity());
                    ////
                    array[n_index_tests] matrix[n_studies, 2] locations = init_array_of_matrices(n_studies, 2, n_index_tests, 0.0);
                    ////
                    for (t in 1:n_index_tests) {
                          ////
                          int n_thr_t = n_thr[t];
                          // int start_index_study_s = 1 + start_index*(s - 1);
                          // int end_index_study_s   = end_index*s;
                          // vector[n_thr_t] C_vec_t = get_test_values(C_vec, start_index_study_s, end_index_study_s, t);
                          ////
                          for (s in 1:n_studies) {
                                vector[n_thr_t] C_vec_t = get_test_t_study_s_segment_values(C_vec, t, s, n_thr, n_studies);
                                ////
                                if (indicator_index_test_in_study[s, t] == 1) {
                                    for (cut_i in 1:n_obs_cutpoints[s, t]) { //// only loop through the OBSERVED cutpoints for test t in study s
                                        ////
                                        int k = to_int(cutpoint_index[t, 1][s, cut_i]);
                                        if (k < n_thr[t]) {
                                            ////
                                            for (c in 1:2) {
                                                //// pi_{s, c, t} = Phi(mu_{c, t} + eta_{s, c} + delta_{s, c, t}):
                                                locations[t][s, c] =  beta_mu[t, c] + beta_eta[s, c] + beta_delta[t][s, c];
                                                ////
                                                latent_cumul_prob[t, c][s, cut_i] = (C_vec_t[k] - locations[t][s, c]);
                                            }
                                        }
                                    }
                                }
                         }
                    }
                    ////
                    //// ---- Calculate CUMULATIVE probabilities:
                    ////
                    for (t in 1:n_index_tests) {
                          int n_thr_t = n_thr[t];
                          for (s in 1:n_studies) {
                             if (indicator_index_test_in_study[s, t] == 1) {
                                 for (c in 1:2) {
                                    cumul_prob[t, c][s, 1:n_thr_t] = Phi_approx(latent_cumul_prob[t, c][s, 1:n_thr_t]); //// INCREASING sequence (as C_k > C_{k - 1})
                                 }
                             }
                          }
                    }
                    ////
                    //// ---- Multinomial (factorised binomial likelihood)
                    ////
                    for (t in 1:n_index_tests) {
                        array[2, 2] matrix[n_studies, n_thr_max] log_lik_outs;
                        log_lik_outs  = compute_log_lik_binomial_fact_NMA( cumul_prob[t, ], 
                                                                          x[t, ],
                                                                          n[t, ],
                                                                          n_obs_cutpoints[, t], 
                                                                          indicator_index_test_in_study[, t]);
                        log_lik[t, ]   = log_lik_outs[1];
                        cond_prob[t, ] = log_lik_outs[2];
                    }
          }
          
                    
 
}


model {
          
          ////
          //// ---- Priors:
          ////
          to_vector(beta_mu)  ~ normal(to_vector(prior_beta_mu_mean), to_vector(prior_beta_mu_SD));
          to_vector(beta_tau) ~ normal(0.0, to_vector(prior_beta_tau_SD));  //// delta_{s, i, t} ~ normal(0, tau_{i, t}):
          beta_Omega ~ lkj_corr(prior_beta_corr_LKJ); //// possibly truncated
          ////
          for (c in 1:2) {
            beta_sigma[c] ~ normal(0.0, prior_beta_sigma_SD[c]);      //// eta[s, i] ~ normal(0, sigma[i]):
          }
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
                    // ////
                    // vector[n_thr_t] rho = std_normal_approx_pdf(Ind_Dir_cumul, Ind_Dir_cumul_prob);
                    // Ind_Dir_ord_prob ~ induced_dirichlet_given_rho(rho, prior_dirichlet_alpha[t][1:n_cat_t]);
          }
          ////
          //// ---- Likelihood / Observational Model:
          ////
          //// (part of between-test / between-study model, NOT prior) - eta[s, i] ~ normal(0, sigma[i]):
          target += std_normal_lpdf(to_vector(beta_eta_z)); 
          ////
          //// (part of between-test / between-study model, NOT prior)  -   delta_{s, i, t} ~ normal(0, tau_{i, t}):
          for (t in 1:n_index_tests) {
             target += std_normal_lpdf(to_vector(beta_delta_z[t])); 
          }
          ////
          ////
          //// ---- Increment the log-likelihood:
          ////
          for (t in 1:n_index_tests) {
               for (s in 1:n_studies) {
                   int observed = to_int(indicator_index_test_in_study[s, t]);
                   if (observed == 1) {
                      for (c in 1:2) {
                         target += sum(log_lik[t, c][s, 1:n_thr[t]]);
                      }
                   }
               }
          }
          
}

 
generated quantities {

          array[n_index_tests] vector[n_thr_max] Se;
          array[n_index_tests] vector[n_thr_max] Sp;
          array[n_index_tests] vector[n_thr_max] Fp;
          ////
          array[n_index_tests] vector[n_thr_max] Se_pred;
          array[n_index_tests] vector[n_thr_max] Sp_pred;
          array[n_index_tests] vector[n_thr_max] Fp_pred;
          ////
          array[n_index_tests] matrix[n_studies, n_thr_max] se;
          array[n_index_tests] matrix[n_studies, n_thr_max] sp;
          array[n_index_tests] matrix[n_studies, n_thr_max] fp;
          ////
          array[n_index_tests] matrix[n_studies, n_thr_max] x_hat_nd = init_array_of_matrices(n_studies, n_thr_max, n_index_tests, -1.0);
          array[n_index_tests] matrix[n_studies, n_thr_max] x_hat_d  = init_array_of_matrices(n_studies, n_thr_max, n_index_tests, -1.0);
          array[n_index_tests] matrix[n_studies, n_thr_max] dev_nd   = init_array_of_matrices(n_studies, n_thr_max, n_index_tests, -1.0);
          array[n_index_tests] matrix[n_studies, n_thr_max] dev_d    = init_array_of_matrices(n_studies, n_thr_max, n_index_tests, -1.0);
          ////
          //// ---- Compute mean / pooled summary-level cutpoints (from Induced-Dirichlet between-study model):
          ////
          array[n_index_tests] vector[n_cat_max] prob_ord_mu;
          array[n_index_tests] vector[n_thr_max] prob_cumul_mu;
          array[n_index_tests] vector[n_thr_max] C_mu;
          ////
          for (t in 1:n_index_tests) {
                int n_thr_t = n_thr[t];
                int n_cat_t = n_cat[t];
                int start_alpha_index = ((t == 1) ? (1 + sum(n_cat[1:(t - 1)])) : 1);
                vector[n_cat_t] alpha_t = segment(alpha_vec, start_alpha_index, n_cat_t);
                ////
                prob_ord_mu[t][1:n_cat_t]   = alpha_t / sum(alpha_t);  //// dirichlet_cat_means_phi[c];
                prob_cumul_mu[t][1:n_thr_t] = ord_probs_to_cumul_probs(prob_ord_mu[t][1:n_cat_t]);
                C_mu[t][1:n_thr_t]          = cumul_probs_to_C(prob_cumul_mu[t][1:n_thr_t], 0.0, 1.0);
          }
          ////
          //// ---- Compute Induced-Dirichlet "SD" params:
          ////
          array[n_index_tests] vector<lower=0.0>[n_cat_max] dirichlet_cat_SDs_sigma;
          vector[n_index_tests] alpha_0;
          for (t in 1:n_index_tests) {
                int n_cat_t = n_cat[t];
                int start_alpha_index = ((t == 1) ? (1 + sum(n_cat[1:(t - 1)])) : 1);
                vector[n_cat_t] alpha_t = segment(alpha_vec, start_alpha_index, n_cat_t);
                ////
                // vector[n_cat_t] alpha_t = get_test_values(alpha, start_index_pooled, end_index_pooled, t);
                alpha_0[t] = sum(alpha_t);
                dirichlet_cat_SDs_sigma[t][1:n_cat_t] = sqrt((alpha_t .* (alpha_0[t] - alpha_t) ) ./ (square(alpha_0[t]) * (alpha_0[t] + 1.0)));
          }
          ////
          //// ---- Calculate study-specific accuracy:
          ////
          for (t in 1:n_index_tests) {
                  int n_thr_t = n_thr[t];
                  fp[t][1:n_studies, 1:n_thr_t] = 1.0 - cumul_prob[t, 1][1:n_studies, 1:n_thr_t];
                  sp[t][1:n_studies, 1:n_thr_t] = 1.0 - fp[t][1:n_studies, 1:n_thr_t];
                  se[t][1:n_studies, 1:n_thr_t] = 1.0 - cumul_prob[t, 2][1:n_studies, 1:n_thr_t];
          }
          ////
          //// ---- Calculate summary accuracy (using mean parameters):
          ////
          for (t in 1:n_index_tests) {
                  int n_thr_t = n_thr[t];
                  vector[n_thr_t] C_mu_t = C_mu[t][1:n_thr_t];
                  // vector[n_thr_t] C_vec_t = get_test_values(C_vec, start_index_study_s, end_index_study_s, t);
                  Fp[t][1:n_thr_t] = 1.0 - Phi_approx((C_mu_t - beta_mu[t, 1]));
                  Sp[t][1:n_thr_t] = 1.0 - Fp[t][1:n_thr_t];
                  Se[t][1:n_thr_t] = 1.0 - Phi_approx((C_mu_t - beta_mu[t, 2]));
          }
          ////
          //// ---- Calculate predictive accuracy:
          ////
          {
            vector[2] beta_eta_pred      = to_vector(normal_rng(rep_vector(0.0, 2), beta_sigma[1:2]));  //// shared between tests
            // ////
            // array[n_index_tests] vector[n_cat_max] prob_ord_pred   = dirichlet_rng(alpha[t][1:n_thr_t]); //// Simulate from Dirichlet by using the summary "alpha" parameters.
            // array[n_index_tests] vector[n_thr_max] prob_cumul_pred = ord_probs_to_cumul_probs(prob_ord_pred);  //// Compute PREDICTED cumulative probabilities.
            // array[n_index_tests] vector[n_thr_max] C_pred = cumul_probs_to_C(prob_cumul_pred, 0.0, 1.0);  //// Compute PREDICTED cutpoints.
            ////
            for (t in 1:n_index_tests) {
              
                  int n_thr_t = n_thr[t];
                  int n_cat_t = n_cat[t];
                  int start_alpha_index = ((t == 1) ? (1 + sum(n_cat[1:(t - 1)])) : 1);
                  vector[n_cat_t] alpha_t = segment(alpha_vec, start_alpha_index, n_cat_t);
                  
                  vector[n_cat_t] prob_ord_pred_t   = dirichlet_rng(alpha_t); //// Simulate from Dirichlet by using the summary "alpha" parameters.
                  vector[n_thr_t] prob_cumul_pred_t = ord_probs_to_cumul_probs(prob_ord_pred_t);  //// Compute PREDICTED cumulative probabilities.
                  vector[n_thr_t] C_pred_t          = cumul_probs_to_C(prob_cumul_pred_t, 0.0, 1.0);  //// Compute PREDICTED cutpoints.
                  ////
                  vector[2] beta_delta_t_pred = to_vector(normal_rng(rep_vector(0.0, 2), beta_tau[t, 1:2]));
                  vector[2] beta_t_pred       = to_vector(beta_mu[t, 1:2]) + beta_eta_pred[1:2] + beta_delta_t_pred[1:2];
                  ////
                  Fp_pred[t][1:n_thr_t] = 1.0 - Phi_approx(C_pred_t - beta_t_pred[1]);
                  Sp_pred[t][1:n_thr_t] = 1.0 - Fp_pred[t][1:n_thr_t];
                  Se_pred[t][1:n_thr_t] = 1.0 - Phi_approx(C_pred_t - beta_t_pred[2]);
            }
          }
          ////
          //// ---- Model-predicted ("re-constructed") data:
          ////
          {
                  array[n_index_tests, 2] matrix[n_studies, n_thr_max] x_hat = init_nested_array_of_matrices(n_studies, n_thr_max, n_index_tests, 2, -1.0);
                  array[n_index_tests, 2] matrix[n_studies, n_thr_max] dev   = init_nested_array_of_matrices(n_studies, n_thr_max, n_index_tests, 2, -1.0);
                  ////
                  for (t in 1:n_index_tests) {
                      for (s in 1:n_studies) {
                            if (indicator_index_test_in_study[s, t] == 1) {
                                      for (c in 1:2) {
                                         for (cut_i in 1:to_int(n_obs_cutpoints[s, t])) {
        
                                                  //// Model-estimated data:
                                                  x_hat[t, c][s, cut_i] = cond_prob[t, c][s, cut_i] * n[t, c][s, cut_i];  	 //// Fitted values
        
                                                  //// Compute residual deviance contribution:
                                                  real n_i =  (n[t, c][s, cut_i]);
                                                  real x_i =  (x[t, c][s, cut_i]);
                                                  real x_hat_i =  (x_hat[t, c][s, cut_i]);
                                                  real log_x_minus_log_x_hat = log(x_i) - log(x_hat_i);
                                                  real log_diff_n_minus_x = log(n_i - x_i);
                                                  real log_diff_n_minus_x_hat = log(abs(n_i - x_hat_i));
        
                                                  // array[n_index_tests, 2] matrix[n_studies, n_thr_max] n;
        
                                                  dev[t, c][s, cut_i] = 2.0 * ( x_i * log_x_minus_log_x_hat + (n_i - x_i) * (log_diff_n_minus_x - log_diff_n_minus_x_hat) );
        
                                         }
                                      }
                            }
                      }
                  }
                  ////
                  //// ---- Store deviance and x_hat split by disease status:
                  ////
                  for (t in 1:n_index_tests) {
                       x_hat_nd[t] = x_hat[t, 1];
                       dev_nd[t]   = dev[t, 1];
                       x_hat_d[t]  = x_hat[t, 2];
                       dev_d[t]    = dev[t, 2];
                  }
          }
          ////
          //// ---- NMA: Compute between-study heterogeneity + correlations for "beta":
          ////
          matrix[n_index_tests, 2] beta_tau_sq = square(beta_tau);
          array[2] matrix[n_index_tests, n_index_tests] beta_sigma_sq;
          array[2] matrix[n_index_tests, n_index_tests] beta_rho;
          array[2] matrix[n_index_tests, n_index_tests] beta_rho12;
          cov_matrix[2] beta_Sigma = multiply_lower_tri_self_transpose(beta_L_Sigma);
          {
            array[3, 2] matrix[n_tests, n_tests] NMA_Nyaga_outs_beta = compute_Nyaga_NMA_summaries(n_tests, beta_tau_sq, beta_Sigma, beta_Omega);
            beta_sigma_sq = NMA_Nyaga_outs_beta[1];
            beta_rho      = NMA_Nyaga_outs_beta[2];
            beta_rho12    = NMA_Nyaga_outs_beta[3];
          }
          ////
          //// ---- NMA: Compute "rank statistics" + other summary (non-comparative) estimates:
          ////
          array[n_index_tests] vector[n_thr_max] LRpos;  
          array[n_index_tests] vector[n_thr_max] LRneg; 
          array[n_index_tests] vector[n_thr_max] DOR;  
          array[n_index_tests] vector[n_thr_max] Youden_index; 
          array[n_index_tests] vector[n_thr_max] Youden_index_weighted; 
          for (t in 1:n_index_tests) {
            
                    int n_thr_t = n_thr[t];
                    array[n_thr_t] int n_thr_t_index = linspaced_int_array(n_thr_t, 1, n_thr_t);
                    vector[n_thr_t] Se_test_t = Se[t][n_thr_t_index];
                    vector[n_thr_t] Sp_test_t = Sp[t][n_thr_t_index];
                    vector[n_thr_t] Fp_test_t = Fp[t][n_thr_t_index];
                    ////
                    //// ---- Youden index:
                    ////
                    Youden_index[t][n_thr_t_index] = Se_test_t + Sp_test_t - 1.0;
                    ////
                    //// ---- Weighted Youden index (weight = 0.5 -> same as youden, weight > 0.5 -> more importance on Se):
                    ////
                    real youden_weight = 0.50;
                    Youden_index_weighted[t][n_thr_t_index] = 2.0 * (youden_weight .* Se_test_t + (1.0 - youden_weight) .* Sp_test_t);
                    ////
                    //// ---- DOR:
                    ////
                  	DOR[t][n_thr_t_index] = (Se_test_t .* Sp_test_t) ./ ((1.0 - Se_test_t) .* (1.0 - Sp_test_t));
                    ////
                    //// ---- "Likelihood Ratio's" (LR's):
                    ////
                   	LRpos[t][n_thr_t_index] = Se_test_t ./ Sp_test_t;
    	              LRneg[t][n_thr_t_index] = (1.0 - Se_test_t) ./ Sp_test_t;
    	              // ////
    	              // //// PPV and NPV (using user-inputted prev data vector)
    	              // ////
    	              // PPV[t][n_thr_t_index] = Se[t] * prev[Test[n]] / ( Se[t]*prev[Test[n]] + (1-Sp[t])*(1.0-prev[Test[n]]) );
    	              // NPV[t][n_thr_t_index] = Sp[t] * prev[Test[n]] / ( (1.0-Se[t])*prev[Test[n]] + Sp[t]*(1.0-prev[Test[n]]) );
    	              
          }
          ////
          //// ---- NMA: Compute pairwise accuracy differences + ratios:
          ////
          array[n_index_tests, n_index_tests] matrix[n_thr_max, n_thr_max] diff_Se;
          array[n_index_tests, n_index_tests] matrix[n_thr_max, n_thr_max] diff_Sp;
          array[n_index_tests, n_index_tests] matrix[n_thr_max, n_thr_max] ratio_Se;
          array[n_index_tests, n_index_tests] matrix[n_thr_max, n_thr_max] ratio_Sp;
          for (t1 in 1:(n_tests - 1)) {
              
                    int n_thr_t1 = n_thr[t1];
                    array[n_thr_t1] int n_thr_t1_index = linspaced_int_array(n_thr_t1, 1, n_thr_t1);
                          
                    for (t2 in (t1 + 1):n_tests) {
                      
                            int n_thr_t2 = n_thr[t2];
                            array[n_thr_t2] int n_thr_t2_index = linspaced_int_array(n_thr_t2, 1, n_thr_t2);
                            ////
                            //// Get vectors of accuracies for all thresholds:
                            ////
                            vector[n_thr_t1] Se_test_t1 = Se[t1][n_thr_t1_index];
                            vector[n_thr_t2] Se_test_t2 = Se[t2][n_thr_t2_index];
                            vector[n_thr_t1] Sp_test_t1 = Sp[t1][n_thr_t1_index];
                            vector[n_thr_t2] Sp_test_t2 = Sp[t2][n_thr_t2_index];
                            ////
                            //// Compute all pairwise differences:
                            ////
                            diff_Se[t1, t2][n_thr_t1_index, n_thr_t2_index] = compute_between_test_diffs( Se_test_t1, Se_test_t2, n_thr_t1, n_thr_t2);
                            diff_Sp[t1, t2][n_thr_t1_index, n_thr_t2_index] = compute_between_test_diffs( Sp_test_t1, Sp_test_t2, n_thr_t1, n_thr_t2);
                            ////
                            //// Compute all pairwise ratios:
                            ////
                            ratio_Se[t1, t2][n_thr_t1_index, n_thr_t2_index] = compute_between_test_ratios(Se_test_t1, Se_test_t2, n_thr_t1, n_thr_t2);
                            ratio_Sp[t1, t2][n_thr_t1_index, n_thr_t2_index] = compute_between_test_ratios(Sp_test_t1, Sp_test_t2, n_thr_t1, n_thr_t2);
                      
                    }
                  
          }


}

