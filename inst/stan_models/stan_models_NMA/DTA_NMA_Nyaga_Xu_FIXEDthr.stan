

functions {
        ////
        //// Include files to compile the necessary custom (user-defined) Stan functions:
        ////
        #include "Stan_fns_basic.stan"
        #include "Stan_fns_corr.stan"
        #include "Stan_fns_ordinal_and_cutpoints.stan"
        #include "Stan_fns_log_lik.stan"
        #include "Stan_fns_ragged.stan"
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
          array[n_index_tests] vector<lower=0.0>[n_thr_max + 1] prior_dirichlet_alpha;
          ////
          //// ---- Other:
          ////
          int<lower=0, upper=1> softplus;
  
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
          ////
          vector[n_total_C_if_fixed] C_raw_vec;  //// RAW LOG-DIFFERENCES - Global cutpoints for each test (staggered array/matrix using "n_thr[t]" to index correctly)
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
        
}


transformed parameters {

          ////
          //// ---- Construct simple 2x2 (bivariate) between-study corr matrices:
          ////
          cholesky_factor_corr[2] beta_L_Omega = make_restricted_bivariate_L_Omega_jacobian(beta_corr, beta_corr_lb, beta_corr_ub);
          corr_matrix[2] beta_Omega      = multiply_lower_tri_self_transpose(beta_L_Omega);
          ////
          //// ---- Construct (global) cutpoints:
          ////
          vector[n_total_C_if_fixed] C_vec;  //// Global cutpoints for each test ("staggered" array/matrix using "n_thr[t]" to index correctly)
          {
            int counter = 1;
            for (t in 1:n_index_tests) {
                  int n_thr_t = n_thr[t];
                  vector[n_thr_t] C_raw_vec_test_t = get_test_values(C_raw_vec, start_index, end_index, t);
                  vector[n_thr_t] C_vec_test_t = ((softplus == 1) ? construct_C_using_SP_jacobian(C_raw_vec_test_t) : construct_C_using_exp_jacobian(C_raw_vec_test_t));
                  C_vec = update_test_values(C_vec, C_vec_test_t, start_index, end_index, t);
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
              beta_eta[s, 1:2] = to_row_vector( diag_pre_multiply(beta_sigma, beta_L_Omega) * to_vector(beta_eta_z[s, 1:2]) );  //// beta_eta[s, 1:2] ~ normal({0.0, 0.0}, Sigma);
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
                          vector[n_thr[t]] C_vec_t = get_test_values(C_vec, start_index, end_index, t);
                          for (s in 1:n_studies) {
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
          {
              // array[n_index_tests] matrix[n_studies, n_thr_max] Ind_Dir_anchor;
              for (t in 1:n_index_tests) {
                    int n_thr_t = n_thr[t];
                    int n_cat_t = n_thr_t + 1;
                    vector[n_thr_t] C_vec_t = get_test_values(C_vec, start_index, end_index, t);
                    ////
                    vector[n_thr_t] Ind_Dir_cumul      = C_vec_t; // - Ind_Dir_anchor;
                    vector[n_thr_t] Ind_Dir_cumul_prob = Phi_approx(Ind_Dir_cumul);
                    vector[n_cat_t] Ind_Dir_ord_prob   = cumul_probs_to_ord_probs(Ind_Dir_cumul_prob);
                    ////
                    vector[n_thr_t] rho = std_normal_approx_pdf(Ind_Dir_cumul, Ind_Dir_cumul_prob);
                    Ind_Dir_ord_prob ~ induced_dirichlet_given_rho(rho, prior_dirichlet_alpha[t][1:n_cat_t]);
              }
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

          //// Summary accuracy parameters (at each threshold):
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
                  Fp[t][1:n_thr_t] =   1.0 - Phi_approx((C_vec_t - beta_mu[t, 1]));
                  Sp[t][1:n_thr_t] =   1.0 - Fp[t][1:n_thr_t];
                  Se[t][1:n_thr_t] =   1.0 - Phi_approx((C_vec_t - beta_mu[t, 2]));
          }
          ////
          //// ---- Calculate predictive accuracy:
          ////
          vector[2] beta_eta_pred      = to_vector(normal_rng(rep_vector(0.0, 2), beta_sigma[1:2]));  //// shared between tests
          ////
          for (t in 1:n_index_tests) {
                int n_thr_t = n_thr[t];
                vector[n_thr_t] C_vec_t = get_test_values(C_vec, start_index, end_index, t);
                ////
                vector[2] beta_delta_t_pred      = to_vector(normal_rng(rep_vector(0.0, 2), beta_tau[t, 1:2]));
                ////
                vector[2] beta_pred = to_vector(beta_mu[t, 1:2]) + beta_eta_pred[1:2] + beta_delta_t_pred[1:2];
                ////
                Fp_pred[t][1:n_thr_t] = 1.0 - Phi_approx(C_vec_t - beta_pred[1]);
                Sp_pred[t][1:n_thr_t] = 1.0 - Fp_pred[t][1:n_thr_t];
                Se_pred[t][1:n_thr_t] = 1.0 - Phi_approx(C_vec_t - beta_pred[2]);
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


}

