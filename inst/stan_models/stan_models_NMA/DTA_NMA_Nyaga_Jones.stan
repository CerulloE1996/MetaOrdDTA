

functions {
        ////
        //// Include files to compile the necessary custom (user-defined) Stan functions:
        ////
        #include "Stan_fns_basic.stan"
        #include "Stan_fns_Box_Cox.stan"
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
          int<lower=1> n_cat_max; // corresponding to the test with the most thresholds/cutpoints. 
          ////
          array[n_studies, n_index_tests] int n_obs_cutpoints; //// OBSERVED cutpoints for test t in study s
          ////
          // array[n_studies] int<lower=0, upper=n_index_tests> n_index_tests_per_study;  // Tests per study
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
          //// ---- Priors for raw_scale:
          ////
          matrix[n_index_tests, 2] prior_raw_scale_mu_mean;
          matrix<lower=0.0>[n_index_tests, 2] prior_raw_scale_mu_SD;
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
  
}

transformed data {
  
          int n_total_C_if_fixed = sum(n_thr); // total # of fixed-effect cutpoints
          int n_tests = n_index_tests;
          ////
          //// ---- Calculate indices:
          ////
          array[n_tests] int start_index = calculate_start_indices(n_thr, n_tests);
          array[n_tests] int end_index   = calculate_end_indices(n_thr, n_tests, start_index);
          
}


parameters {

          matrix[n_index_tests, 2] beta_mu;
          matrix[n_index_tests, 2] raw_scale_mu;
          ////
          vector<lower=-5.0, upper=5.0>[n_index_tests] lambda;
          ////
          //// ---- "NMA" params:
          ////
          vector<lower=0.0>[2] beta_sigma;              //// Variance component - Between-study SD (Nyaga's σ)
          matrix<lower=0.0>[n_index_tests, 2] beta_tau; //// Variance component - Test-specific SD (Nyaga's τ) - delta_{s, c, t} ~ normal(0, tau_{c, t}).
          matrix[n_studies, 2] beta_eta_z;              //// Standard normal RVs for study-level effects - eta[s, 1:2] ~ multi_normal({0, 0}, Sigma).
          array[n_index_tests] matrix[n_studies, 2] beta_delta_z; //// Standard normal RVs for test-specific effects
          ////
          vector<lower=0.0>[2] raw_scale_sigma;
          matrix<lower=0.0>[n_index_tests, 2] raw_scale_tau;
          matrix[n_studies, 2] raw_scale_eta_z;
          array[n_index_tests] matrix[n_studies, 2] raw_scale_delta_z;
          ////
          //// ---- Between-study corr's:
          ////
          real beta_corr;       //// between-study corr (possibly restricted)
          real raw_scale_corr;  //// between-study corr (possibly restricted)
        
}


transformed parameters {
    
          ////
          //// ---- Construct simple 2x2 (bivariate) between-study corr matrices:
          ////
          cholesky_factor_corr[2] beta_L_Omega      = make_restricted_bivariate_L_Omega_jacobian(beta_corr, beta_corr_lb, beta_corr_ub);
          cholesky_factor_cov[2]  beta_L_Sigma      = diag_pre_multiply(beta_sigma, beta_L_Omega);
          corr_matrix[2] beta_Omega      = multiply_lower_tri_self_transpose(beta_L_Omega);
          ////
          cholesky_factor_corr[2] raw_scale_L_Omega = make_restricted_bivariate_L_Omega_jacobian(raw_scale_corr, raw_scale_corr_lb, raw_scale_corr_ub);
          cholesky_factor_cov[2]  raw_scale_L_Sigma = diag_pre_multiply(raw_scale_sigma, raw_scale_L_Omega);
          corr_matrix[2] raw_scale_Omega = multiply_lower_tri_self_transpose(raw_scale_L_Omega);
          ////
          //// ---- Cutpoint params:
          ////
          vector[n_total_C_if_fixed] C_vec;  //// Global cutpoints for each test ("staggered" array/matrix using "n_thr[t]" to index correctly)
          ////
          //// ---- Construct cutpoints for each test:
          ////
          for (t in 1:n_index_tests) {
                int n_thr_t = n_thr[t];
                vector[n_thr_t] C_vec_test_t = ((box_cox == 0) ? log(cts_thr_values[t][1:n_thr[t]]) : fn_Stan_box_cox(cts_thr_values[t][1:n_thr[t]], lambda[t]));
                C_vec = update_test_values(C_vec, C_vec_test_t, start_index, end_index, t);
                //// C[t][1:n_thr[t]] = ((box_cox == 0) ? log(cts_thr_values[t][1:n_thr[t]]) : fn_Stan_box_cox(cts_thr_values[t][1:n_thr[t]], lambda[t]));
          }
          ////
          //// ---- "NMA" params:
          //// ---- Declare Study-level random effects (eta in Nyaga notation) - eta[s, 1:2] ~ multi_normal({0, 0}, Sigma):
          ////
          matrix[n_studies, 2] beta_eta      = rep_matrix(0.0, n_studies, 2);
          matrix[n_studies, 2] raw_scale_eta = rep_matrix(0.0, n_studies, 2);
          ////
          //// ---- Compute Study-level random effects (eta in Nyaga notation):
          ////
          for (s in 1:n_studies) {
              beta_eta[s, 1:2]      = to_row_vector( beta_L_Sigma      * to_vector(beta_eta_z[s, 1:2]) );  //// beta_eta[s, 1:2] ~ normal({0.0, 0.0}, Sigma);
              raw_scale_eta[s, 1:2] = to_row_vector( raw_scale_L_Sigma * to_vector(raw_scale_eta_z[s, 1:2]) );  //// beta_eta[s, 1:2] ~ normal({0.0, 0.0}, Sigma);
          }
          ////
          //// ---- Declare test-specific deviations ("delta" in Nyaga notation) - delta_{s, c, t} ~ normal(0, tau_{c, t}):
          ////
          array[n_index_tests] matrix[n_studies, 2] beta_delta      = init_array_of_matrices(n_studies, 2, n_index_tests, 0.0);
          array[n_index_tests] matrix[n_studies, 2] raw_scale_delta = init_array_of_matrices(n_studies, 2, n_index_tests, 0.0);
          ////
          //// ---- Compute test-specific deviations ("delta" in Nyaga notation) - delta_{s, c, t} ~ normal(0, tau_{c, t}):
          ////
          for (t in 1:n_index_tests) {
                // for (s in 1:n_studies) {
                      for (c in 1:2) {
                          //// delta's ("beta_delta") correspond to shared (between tests) component of "beta"   -   delta_{s, c, t} ~ normal(0, tau_{c, t}):
                          beta_delta[t][1:n_studies, c]      = beta_tau[t, c] * beta_delta_z[t][1:n_studies, c];           //// beta_delta[t][s, c]      ~ normal(0.0, beta_tau[t, c]);
                          raw_scale_delta[t][1:n_studies, c] = raw_scale_tau[t, c] * raw_scale_delta_z[t][1:n_studies, c]; //// raw_scale_delta[t][s, c] ~ normal(0.0, raw_scale_tau[t, c]);
                      }
                // }
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
                    array[n_index_tests] matrix[n_studies, 2] scales    = init_array_of_matrices(n_studies, 2, n_index_tests, 1.0);
                    ////
                    for (t in 1:n_index_tests) {
                          vector[n_thr[t]] C_vec_t = get_test_values(C_vec, start_index, end_index, t);
                          for (s in 1:n_studies) {
                                if (indicator_index_test_in_study[s, t] == 1) {
          
                                        for (cut_i in 1:n_obs_cutpoints[s, t]) { //// only loop through the OBSERVED cutpoints for test t in study s
                                                ////
                                                int k = to_int(cutpoint_index[t, 1][s, cut_i]);
                                                if (k < n_thr[t]) {
                                                    ////
                                                    for (c in 1:2) {
                                                        //// pi_{s, c, t} = Phi(mu_{c, t} + eta_{s, c} + delta_{s, c, t}).
                                                        locations[t][s, c]     =  beta_mu[t, c]      + beta_eta[s, c]      + beta_delta[t][s, c]; 
                                                        real raw_scale_given_c =  raw_scale_mu[t, c] + raw_scale_eta[s, c] + raw_scale_delta[t][s, c];
                                                        scales[t][s, c] = ((softplus == 1) ? softplus_scaled_jacobian(raw_scale_given_c) : exp_jacobian(raw_scale_given_c));
                                                    }
                                                    ////
                                                    //// latent_cumul_prob[c, t, s, cut_i] = (C[t][k] - location)/scale;
                                                    for (c in 1:2) {
                                                        latent_cumul_prob[t, c][s, cut_i] = (C_vec_t[k] - locations[t][s, c])/scales[t][s, c];
                                                    }
                                                }
                                        }
                                }
                         }
                    }
                    // if (abs(sum(raw_scale_SD)) != 0.0) target += log(abs(sum(raw_scale_SD)));      // double-checked the log-derivative of this by hand (correct)
                    // if (abs(sum(raw_scale_z)) != 0.0)  target += log(abs(sum(raw_scale_z)));  // double-checked the log-derivative of this by hand (correct)
                    ////
                    //// ---- Calculate CUMULATIVE probabilities (vectorised):
                    ////
                    for (t in 1:n_index_tests) {
                          for (s in 1:n_studies) {
                             if (indicator_index_test_in_study[s, t] == 1) {
                                for (c in 1:2) {
                                    cumul_prob[t, c][s, 1:n_thr[t]] = Phi_approx(latent_cumul_prob[t, c][s, 1:n_thr[t]]); //// INCREASING sequence (as C_k > C_{k - 1})
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
          lambda ~ normal(prior_boxcox_lambda_mean, prior_boxcox_lambda_SD);
          ////
          to_vector(beta_mu)  ~ normal(to_vector(prior_beta_mu_mean), to_vector(prior_beta_mu_SD));
          to_vector(beta_tau) ~ normal(0.0, to_vector(prior_beta_tau_SD));  //// delta_{s, i, t} ~ normal(0, tau_{i, t}):
          beta_Omega ~ lkj_corr(prior_beta_corr_LKJ);
          ////
          to_vector(raw_scale_mu)  ~ normal(to_vector(prior_raw_scale_mu_mean), to_vector(prior_raw_scale_mu_SD));
          to_vector(raw_scale_tau) ~ normal(0.0, to_vector(prior_beta_tau_SD));  //// delta_{s, i, t} ~ normal(0, tau_{i, t}):
          raw_scale_Omega ~ lkj_corr(prior_raw_scale_corr_LKJ);
          ////
          for (c in 1:2) {
            beta_sigma[c]      ~ normal(0.0, prior_beta_sigma_SD[c]);      //// eta[s, i] ~ normal(0, sigma[i]):
            raw_scale_sigma[c] ~ normal(0.0, prior_raw_scale_sigma_SD[c]); //// eta[s, i] ~ normal(0, sigma[i]):
          }
          ////
          target += std_normal_lpdf(to_vector(beta_eta_z));      //// (part of between-test / between-study model, NOT prior) - eta[s, i] ~ normal(0, sigma[i]):
          target += std_normal_lpdf(to_vector(raw_scale_eta_z)); //// (part of between-test / between-study model, NOT prior) - eta[s, i] ~ normal(0, sigma[i]):
          ////
          for (t in 1:n_index_tests) {
             target += std_normal_lpdf(to_vector(beta_delta_z[t])); ////    (part of between-test / between-study model, NOT prior)  -   delta_{s, i, t} ~ normal(0, tau_{i, t}):
             target += std_normal_lpdf(to_vector(raw_scale_delta_z[t])); ////   (part of between-test / between-study model, NOT prior)  -   delta_{s, i, t} ~ normal(0, tau_{i, t}):
          }
          ////
          //// ---- Increment the log-likelihood:
          ////
          for (t in 1:n_index_tests) {
               for (s in 1:n_studies) {
                   int observed = to_int(indicator_index_test_in_study[s, t]);
                   if (observed == 1) {
                     for (c in 1:2) {
                       target +=  sum(log_lik[t, c][s, 1:n_thr[t]]);
                     }
                   }
               }
          }
    
}


generated quantities {
  
          ////
          //// ---- Summary accuracy parameters (at each threshold):
          ////
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
          //// ---- Calculate study-specific accuracy:
          ////
          for (t in 1:n_index_tests) {
                  int n_thr_t = n_thr[t];
                  fp[t][1:n_studies, 1:n_thr_t] =   1.0 - cumul_prob[t, 1][1:n_studies, 1:n_thr_t];
                  sp[t][1:n_studies, 1:n_thr_t] =   1.0 - fp[t][1:n_studies, 1:n_thr_t];
                  se[t][1:n_studies, 1:n_thr_t] =   1.0 - cumul_prob[t, 2][1:n_studies, 1:n_thr_t];
          }
          ////
          //// ---- Calculate summary accuracy (using mean parameters):
          ////
          for (t in 1:n_index_tests) {
                  int n_thr_t = n_thr[t];
                  vector[n_thr_t] C_vec_t = get_test_values(C_vec, start_index, end_index, t);
                  ////
                  real scale_mu_nd_t = ((softplus == 1) ? softplus_scaled(raw_scale_mu[t, 1]) : exp(raw_scale_mu[t, 1]));
                  real scale_mu_d_t  = ((softplus == 1) ? softplus_scaled(raw_scale_mu[t, 2]) : exp(raw_scale_mu[t, 2]));
                  ////
                  Fp[t][1:n_thr_t] =   1.0 - Phi_approx((C_vec_t - beta_mu[t, 1])/scale_mu_nd_t);
                  Sp[t][1:n_thr_t] =   1.0 - Fp[t][1:n_thr_t];
                  Se[t][1:n_thr_t] =   1.0 - Phi_approx((C_vec_t - beta_mu[t, 2])/scale_mu_d_t);
          }
          ////
          //// ---- Calculate predictive accuracy:
          ////
          vector[2] beta_eta_pred      = to_vector(normal_rng(rep_vector(0.0, 2), beta_sigma[1:2]));  //// shared between tests
          vector[2] raw_scale_eta_pred = to_vector(normal_rng(rep_vector(0.0, 2), raw_scale_sigma[1:2]));  //// shared between tests
          ////
          for (t in 1:n_index_tests) {
                int n_thr_t = n_thr[t];
                vector[n_thr_t] C_vec_t = get_test_values(C_vec, start_index, end_index, t);
                ////
                real beta_nd_delta_t_pred      = normal_rng(0.0, beta_tau[t, 1]);
                real beta_d_delta_t_pred       = normal_rng(0.0, beta_tau[t, 2]);
                real raw_scale_nd_delta_t_pred = normal_rng(0.0, raw_scale_tau[t, 1]);
                real raw_scale_d_delta_t_pred  = normal_rng(0.0, raw_scale_tau[t, 2]);
                ////
                real beta_nd_pred_t      = beta_mu[t, 1] + beta_eta_pred[1] + beta_nd_delta_t_pred;
                real beta_d_pred_t       = beta_mu[t, 2] + beta_eta_pred[2] + beta_d_delta_t_pred;
                ////
                real raw_scale_nd_pred_t = raw_scale_mu[t, 1] + raw_scale_eta_pred[1]  + raw_scale_nd_delta_t_pred;
                real raw_scale_d_pred_t  = raw_scale_mu[t, 2] + raw_scale_eta_pred[2]  + raw_scale_d_delta_t_pred;
                real scale_pred_nd_t = ((softplus == 1) ? softplus_scaled(raw_scale_nd_pred_t) : exp(raw_scale_nd_pred_t));
                real scale_pred_d_t  = ((softplus == 1) ? softplus_scaled(raw_scale_d_pred_t)  : exp(raw_scale_d_pred_t));
                ////
                {
                      Fp_pred[t][1:n_thr_t] =   1.0 - Phi_approx((C_vec_t - beta_nd_pred_t)/scale_pred_nd_t);
                      Sp_pred[t][1:n_thr_t] =   1.0 - Fp_pred[t][1:n_thr_t];
                      Se_pred[t][1:n_thr_t] =   1.0 - Phi_approx((C_vec_t - beta_d_pred_t)/scale_pred_d_t);
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
          //// ---- NMA: Compute between-study heterogeneity + correlations for "raw_scale":
          ////
          matrix[n_index_tests, 2] raw_scale_tau_sq = square(raw_scale_tau);
          array[2] matrix[n_index_tests, n_index_tests] raw_scale_sigma_sq;
          array[2] matrix[n_index_tests, n_index_tests] raw_scale_rho;
          array[2] matrix[n_index_tests, n_index_tests] raw_scale_rho12;
          ////
          cov_matrix[2] raw_scale_Sigma = multiply_lower_tri_self_transpose(raw_scale_L_Sigma);
          {
            array[3, 2] matrix[n_tests, n_tests] NMA_Nyaga_outs_beta = compute_Nyaga_NMA_summaries(n_tests, raw_scale_tau_sq, raw_scale_Sigma, raw_scale_Omega);
            raw_scale_sigma_sq = NMA_Nyaga_outs_beta[1];
            raw_scale_rho      = NMA_Nyaga_outs_beta[2];
            raw_scale_rho12    = NMA_Nyaga_outs_beta[3];
          }
          ////
          //// ---- NMA: Compute "rank statistics" + other summary (non-comparative) estimates:
          ////
          array[n_index_tests] vector[n_thr_max] LRpos;  
          array[n_index_tests] vector[n_thr_max] LRneg; 
          array[n_index_tests] vector[n_thr_max] DOR; 
          array[n_index_tests] vector[n_thr_max] Youden_index; 
          array[n_index_tests] vector[n_thr_max] Youden_index_weighted; 
          //// Initialise:
          for (t in 1:n_index_tests) {
                    LRpos[t] = rep_vector(-1.0, n_thr_max);
                    LRneg[t] = rep_vector(-1.0, n_thr_max);
                    DOR[t]   = rep_vector(-1.0, n_thr_max);
                    Youden_index[t]          = rep_vector(-1.0, n_thr_max);
                    Youden_index_weighted[t] = rep_vector(-1.0, n_thr_max);
          }
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
          array[n_index_tests, n_index_tests] matrix[n_thr_max, n_thr_max] diff_Se  = init_nested_array_of_matrices(n_thr_max, n_thr_max, n_index_tests, n_index_tests, -1.0);
          array[n_index_tests, n_index_tests] matrix[n_thr_max, n_thr_max] diff_Sp  = init_nested_array_of_matrices(n_thr_max, n_thr_max, n_index_tests, n_index_tests, -1.0);
          array[n_index_tests, n_index_tests] matrix[n_thr_max, n_thr_max] ratio_Se = init_nested_array_of_matrices(n_thr_max, n_thr_max, n_index_tests, n_index_tests, -1.0);
          array[n_index_tests, n_index_tests] matrix[n_thr_max, n_thr_max] ratio_Sp = init_nested_array_of_matrices(n_thr_max, n_thr_max, n_index_tests, n_index_tests, -1.0);
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
















