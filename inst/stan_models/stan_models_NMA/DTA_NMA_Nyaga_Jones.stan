

functions {
        ////
        //// Include files to compile the necessary custom (user-defined) Stan functions:
        ////
        #include "Stan_fns_basic.stan"
        #include "Stan_fns_Box_Cox.stan"
        #include "Stan_fns_corr.stan"
        #include "Stan_fns_ordinal_and_cutpoints.stan"
        #include "Stan_fns_log_lik.stan"
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
          //// ---- Other:
          ////
          int<lower=0, upper=1> box_cox;
          int<lower=0, upper=1> softplus;
          array[n_index_tests] vector[n_thr_max] cts_thr_values;
          // int n_total_cutpoints;
          // int n_total_raw_simplex_elements;
  
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
        
}


transformed parameters {
    
          // ////
          // //// ---- Cutpoint params:
          // //// Global cutpoints for each test ("staggered" array/matrix using "n_thr[t]" to index correctly) - "staggered" array
          // ////
          // array[n_index_tests] vector[n_thr_max] C;
          // ////
          // //// ---- Construct cutpoints for each test:
          // ////
          // for (t in 1:n_index_tests) {
          //       C[t][1:n_thr[t]] = ((box_cox == 0) ? log(cts_thr_values[t][1:n_thr[t]]) : fn_Stan_box_cox(cts_thr_values[t][1:n_thr[t]], lambda[t])); 
          // }
          // ////
          // array[n_index_tests, 2] matrix[n_studies, n_thr_max] latent_cumul_prob;
          // array[n_index_tests, 2] matrix[n_studies, n_thr_max] cumul_prob;
          // array[n_index_tests, 2] matrix[n_studies, n_thr_max] cond_prob;
          // array[n_index_tests, 2] matrix[n_studies, n_thr_max] log_lik;
          // ////
          // //// ---- Initialise matrices:
          // ////
          // for (t in 1:n_index_tests) { 
          //     for (c in 1:2) {
          //              latent_cumul_prob[t, c]  = rep_matrix(positive_infinity(), n_studies, n_thr_max);
          //              cumul_prob[t, c]         = rep_matrix(1.0, n_studies, n_thr_max);
          //              cond_prob[t, c]          = rep_matrix(1.0, n_studies, n_thr_max);
          //              log_lik[t, c]            = rep_matrix(0.0, n_studies, n_thr_max);
          //     }
          // }
          // ////
          // //// ---- "NMA" params:
          // ////
          // //// ---- Declare Study-level random effects (eta in Nyaga notation) -  eta[s, 1:2] ~ multi_normal({0, 0}, Sigma):
          // ////
          // matrix[n_studies, 2] beta_eta      = rep_matrix(0.0, n_studies, 2);
          // matrix[n_studies, 2] raw_scale_eta = rep_matrix(0.0, n_studies, 2);
          // ////
          // //// ---- Compute Study-level random effects (eta in Nyaga notation):
          // ////
          // for (s in 1:n_studies) {
          //   for (c in 1:2) {
          //     //// eta's ("beta_eta") correspond to shared (between tests) component of "beta" - eta_{s, i} ~ normal(0, sigma_{i}).
          //     beta_eta[s, c]      = 0.0 + beta_sigma[c]      * beta_eta_z[s, c];      
          //     raw_scale_eta[s, c] = 0.0 + raw_scale_sigma[c] * raw_scale_eta_z[s, c];
          //   }
          // }
          // ////
          // //// ---- Declare test-specific deviations ("delta" in Nyaga notation) - delta_{s, c, t} ~ normal(0, tau_{c, t}):
          // ////
          // array[n_index_tests] matrix[n_studies, 2] beta_delta;
          // array[n_index_tests] matrix[n_studies, 2] raw_scale_delta;
          // for (t in 1:n_index_tests) {
          //    beta_delta[t]      = rep_matrix(0.0, n_studies, 2);
          //    raw_scale_delta[t] = rep_matrix(0.0, n_studies, 2);
          // }
          // ////
          // //// ---- Compute test-specific deviations ("delta" in Nyaga notation) - delta_{s, c, t} ~ normal(0, tau_{c, t}):
          // ////
          // for (t in 1:n_index_tests) {
          //     for (s in 1:n_studies) {
          //         if (indicator_index_test_in_study[s, t] == 1) {
          //           for(c in 1:2) {
          //               //// delta's ("beta_delta") correspond to shared (between tests) component of "beta"   -   delta_{s, c, t} ~ normal(0, tau_{c, t}):
          //               beta_delta[t][s, c]      = 0.0 + beta_delta_z[t][s, c]      * beta_tau[t, c];     
          //               raw_scale_delta[t][s, c] = 0.0 + raw_scale_delta_z[t][s, c] * raw_scale_tau[t, c];  
          //           }
          //         }
          //     }
          // }
          // ////  
          // //// Compute probit CUMULATIVE probabilities:
          // ////
          // //// real Jacobian_raw_scale_to_scale = 0.0;
          // array[n_index_tests] matrix[n_studies, 2] locations;
          // array[n_index_tests] matrix[n_studies, 2] scales;
          // ////
          // for (t in 1:n_index_tests) { 
          //   locations[t] = rep_matrix(0.0, 2, n_studies);
          //   scales[t]    = rep_matrix(0.0, 2, n_studies);
          // }
          // ////
          // for (t in 1:n_index_tests) {
          //       for (s in 1:n_studies) {
          //             if (indicator_index_test_in_study[s, t] == 1) {
          //               
          //                     for (cut_i in 1:n_obs_cutpoints[s, t]) { //// only loop through the OBSERVED cutpoints for test t in study s
          //                           
          //                           int k = to_int(cutpoint_index[t, 1][s, cut_i]);
          //                           ////
          //                           locations[t][s, 1] =  beta_mu[t, 1] + beta_eta[s, 1] + beta_delta[t][s, 1]; //// pi_{s, c, t} = Phi(mu_{c, t} + eta_{s, c} + delta_{s, c, t}).
          //                           locations[t][s, 2] =  beta_mu[t, 2] + beta_eta[s, 2] + beta_delta[t][s, 2]; //// pi_{s, c, t} = Phi(mu_{c, t} + eta_{s, c} + delta_{s, c, t}).
          //                           ////
          //                           real raw_scale_nd_baseline = raw_scale_mu[t, 1] + raw_scale_eta[s, 1] + raw_scale_delta[t][s, 1];
          //                           real raw_scale_d_baseline  = raw_scale_mu[t, 2] + raw_scale_eta[s, 2] + raw_scale_delta[t][s, 2]; 
          //                           ////
          //                           scales[t][s, 1] = ((softplus == 1) ? softplus_scaled_jacobian(raw_scale_nd_baseline) : exp_jacobian(raw_scale_nd_baseline)); 
          //                           scales[t][s, 2] = ((softplus == 1) ? softplus_scaled_jacobian(raw_scale_d_baseline)  : exp_jacobian(raw_scale_d_baseline)); 
          //                           // //// Jacobian for raw_scale -> scale:
          //                           // jacobian += +log(2) - 0.5*raw_scale_baseline; //// deriv of exp(-0.5*raw_scale_baseline) = -0.5*exp(-0.5*raw_scale_baseline) -> log_abs_deriv = +log(2) -0.5*raw_scale_baseline; 
          //                           // jacobian += +log(2) + 0.5*raw_scale_baseline; //// deriv of exp(+0.5*raw_scale_baseline) = +0.5*exp(+0.5*raw_scale_baseline) -> log_abs_deriv = -log(2) +0.5*raw_scale_baseline;
          //                           // 
          //                           //// latent_cumul_prob[c, t, s, cut_i] = (C[t][k] - location)/scale;
          //                           for (c in 1:2) {
          //                               latent_cumul_prob[t, c][s, cut_i] = (C[t][k] - locations[t][s, c])/scales[t][s, c];
          //                           }
          //                       
          //                     }
          //             }
          //             
          //      } 
          // }
          // 
          //   // if (abs(sum(raw_scale_SD)) != 0.0) target += log(abs(sum(raw_scale_SD)));      // double-checked the log-derivative of this by hand (correct)
          //   // if (abs(sum(raw_scale_z)) != 0.0)  target += log(abs(sum(raw_scale_z)));  // double-checked the log-derivative of this by hand (correct)
          //   
          //   
          // ////
          // //// ---- Calculate CUMULATIVE probabilities (vectorised):
          // ////
          // for (t in 1:n_index_tests) {
          //     for (c in 1:2) {
          //         cumul_prob[t, c] = Phi_approx(latent_cumul_prob[t, c]); //// INCREASING sequence (as C_k > C_{k - 1})
          //     }
          // }
          // ////
          // //// ---- Multinomial (factorised binomial likelihood)
          // ////
          // for (t in 1:n_index_tests) {
          //     log_lik[t, ] = compute_log_lik_binomial_fact(cumul_prob[t, ], x[t,,,], n[t,,,], n_obs_cutpoints[, t]);
          // }
      
}


model {
          
          ////
          //// ---- Priors:
          ////
          for (t in 1:n_index_tests) {
              lambda[t] ~ normal(prior_boxcox_lambda_mean[t], prior_boxcox_lambda_SD[t]);
              beta_mu[t, 1:2]      ~ normal(prior_beta_mu_mean[t, 1:2], prior_beta_mu_SD[t, 1:2]);
              raw_scale_mu[t, 1:2] ~ normal(prior_raw_scale_mu_mean[t, 1:2], prior_raw_scale_mu_SD[t, 1:2]);
          }
          ////
          for (t in 1:n_index_tests) {
               beta_tau[t, 1:2]        ~ normal(0.0, prior_beta_tau_SD[t, 1:2]);       //// delta_{s, i, t} ~ normal(0, tau_{i, t}):
               raw_scale_tau[t, 1:2]   ~ normal(0.0, prior_raw_scale_tau_SD[t, 1:2]);  //// delta_{s, i, t} ~ normal(0, tau_{i, t}):
          }
          ////
          beta_sigma[1:2]      ~ normal(0.0, prior_beta_sigma_SD[1:2]);      //// eta[s, i] ~ normal(0, sigma[i]):
          raw_scale_sigma[1:2] ~ normal(0.0, prior_raw_scale_sigma_SD[1:2]); //// eta[s, i] ~ normal(0, sigma[i]):
          ////
          //// ---- Likelihood / Model:
          ////
          {
              target += std_normal_lpdf(to_vector(beta_eta_z)); ////    // (part of between-test / between-study model, NOT prior) - eta[s, i] ~ normal(0, sigma[i]):
              target += std_normal_lpdf(to_vector(raw_scale_eta_z)); //  // (part of between-test / between-study model, NOT prior) - eta[s, i] ~ normal(0, sigma[i]):
          }
          ////
          for (t in 1:n_index_tests) {
             target += std_normal_lpdf(to_vector(beta_delta_z[t])); ////     ~ std_normal(); // (part of between-test / between-study model, NOT prior)  -   delta_{s, i, t} ~ normal(0, tau_{i, t}):
             target += std_normal_lpdf(to_vector(raw_scale_delta_z[t])); //// ~ std_normal(); // (part of between-test / between-study model, NOT prior)  -   delta_{s, i, t} ~ normal(0, tau_{i, t}):
          }
          // //
          // // ---- Increment the log-likelihood:
          // ////
          // for (t in 1:n_index_tests) {
          //      for (s in 1:n_studies) {
          //          int observed = to_int(indicator_index_test_in_study[s, t]);
          //          if (observed == 1) {
          //            for (c in 1:2) {
          //              target +=  sum(log_lik[t, c][s, 1:n_thr[t]]);
          //            }
          //          }
          //      }
          // }
    
}


// generated quantities {
// 
//           //// Summary accuracy parameters (at each threshold):
//           array[n_index_tests] vector[n_thr_max] Se;
//           array[n_index_tests] vector[n_thr_max] Sp;
//           array[n_index_tests] vector[n_thr_max] Fp;
//           ////
//           array[n_index_tests] vector[n_thr_max] Se_pred;
//           array[n_index_tests] vector[n_thr_max] Sp_pred;
//           array[n_index_tests] vector[n_thr_max] Fp_pred;
//           array[n_index_tests] vector[n_thr_max] C_mu_pred;
//           ////
//           array[n_index_tests] matrix[n_studies, n_thr_max] se;
//           array[n_index_tests] matrix[n_studies, n_thr_max] sp;
//           array[n_index_tests] matrix[n_studies, n_thr_max] fp;
//           ////
//           array[n_index_tests] matrix[n_studies, n_thr_max] x_hat_nd;// = rep_matrix(-1, n_studies, n_thr_max);
//           array[n_index_tests] matrix[n_studies, n_thr_max] x_hat_d;//  = rep_matrix(-1, n_studies, n_thr_max);
//           array[n_index_tests] matrix[n_studies, n_thr_max] dev_nd;//    = rep_matrix(-1, n_studies, n_thr_max);
//           array[n_index_tests] matrix[n_studies, n_thr_max] dev_d;//     = rep_matrix(-1, n_studies, n_thr_max);
//           array[2, n_index_tests] matrix[n_studies, n_thr_max] x_hat;
//           array[2, n_index_tests] matrix[n_studies, n_thr_max] dev;
//           ////
//           //// ---- Initialise containers:
//           ////
//           for (t in 1:n_index_tests) {
//             for (c in 1:2) {
//                 x_hat[t, c]  = rep_matrix(-1, n_studies, n_thr_max);
//                 dev[t, c]    = rep_matrix(-1, n_studies, n_thr_max);
//             }
//           }
//           ////
//           //// ---- Calculate study-specific accuracy:
//           ////
//           for (s in 1:n_studies) {
//               for (t in 1:n_index_tests) {
//                   for (k in 1:n_thr[t]) {
//                       fp[t][s, k] =   1.0 - cumul_prob[t, 1][s, k];
//                       sp[t][s, k] =   1.0 - fp[t][s, k];
//                       se[t][s, k] =   1.0 - cumul_prob[t, 2][s, k];
//                   }
//              }
//           }
//           ////
//           //// ---- Calculate summary accuracy (using mean parameters):
//           ////
//           for (t in 1:n_index_tests) {
//             for (k in 1:n_thr[t]) {
//                   Fp[t][k] =   1.0 - Phi_approx((C[t][k] - beta_mu[t, 1])/softplus_scaled(raw_scale_mu[t, 1]));
//                   Sp[t][k] =   1.0 - Fp[t][k];
//                   Se[t][k] =   1.0 - Phi_approx((C[t][k] - beta_mu[t, 2])/softplus_scaled(raw_scale_mu[t, 2]));
//             }
//           }
//           ////
//           //// ---- Calculate predictive accuracy:
//           ////
//           real beta_nd_eta_pred      = normal_rng(0.0, beta_sigma[1]);      //// shared between tests
//           real beta_d_eta_pred       = normal_rng(0.0, beta_sigma[2]);      //// shared between tests
//           real raw_scale_nd_eta_pred = normal_rng(0.0, raw_scale_sigma[1]); //// shared between tests
//           real raw_scale_d_eta_pred  = normal_rng(0.0, raw_scale_sigma[2]); //// shared between tests
//           ////
//           for (t in 1:n_index_tests) {
//                 ////
//                 real beta_nd_delta_t_pred      = normal_rng(0.0, beta_tau[t, 1]);
//                 real beta_d_delta_t_pred       = normal_rng(0.0, beta_tau[t, 2]);
//                 real raw_scale_nd_delta_t_pred = normal_rng(0.0, raw_scale_tau[t, 1]);
//                 real raw_scale_d_delta_t_pred  = normal_rng(0.0, raw_scale_tau[t, 2]);
//                 ////
//                 real beta_nd_pred      = beta_nd_eta_pred       + beta_nd_delta_t_pred;
//                 real beta_d_pred       = beta_d_eta_pred        + beta_d_delta_t_pred;
//                 real raw_scale_nd_pred = raw_scale_nd_eta_pred  + raw_scale_nd_delta_t_pred;
//                 real raw_scale_d_pred  = raw_scale_d_eta_pred   + raw_scale_d_delta_t_pred;
//                 ////
//                 for (k in 1:n_thr[t]) {
//                       Fp_pred[t][k] =   1.0 - Phi_approx((C_mu_pred[t][k] - beta_nd_pred)/softplus_scaled(raw_scale_nd_pred));
//                       Sp_pred[t][k] =   1.0 - Fp_pred[t][k];
//                       Se_pred[t][k] =   1.0 - Phi_approx((C_mu_pred[t][k] - beta_d_pred)/softplus_scaled(raw_scale_d_pred));
//                 }
//           }
//           ////
//           //// ---- Model-predicted ("re-constructed") data:
//           ////
//           for (t in 1:n_index_tests) {
//               for (s in 1:n_studies) {
//                     if (indicator_index_test_in_study[s, t] == 1) {
// 
//                               for (c in 1:2) {
//                                  for (cut_i in 1:to_int(n_obs_cutpoints[s, t])) {
// 
//                                           //// Model-estimated data:
//                                           x_hat[c][t][s, cut_i] = cond_prob[c][t][s, cut_i] * n[c][t][s, cut_i];  	 //// Fitted values
// 
//                                           //// Compute residual deviance contribution:
//                                           real n_i =  (n[t, c, s, cut_i]);
//                                           real x_i =  (x[t, c, s, cut_i]);
//                                           real x_hat_i =  (x_hat[c][t][s, cut_i]);
//                                           real log_x_minus_log_x_hat = log(x_i) - log(x_hat_i);
//                                           real log_diff_n_minus_x = log(n_i - x_i);
//                                           real log_diff_n_minus_x_hat = log(abs(n_i - x_hat_i));
// 
//                                           dev[c][t][s, cut_i] = 2.0 * ( x_i * log_x_minus_log_x_hat + (n_i - x_i) * (log_diff_n_minus_x - log_diff_n_minus_x_hat) );
// 
//                                  }
//                               }
// 
//                     }
//               }
//           }
//           ////
//           //// ---- Store deviance and x_hat split by disease status:
//           ////
//           for (t in 1:n_index_tests) {
//                x_hat_nd[t] = x_hat[t][1];
//                dev_nd[t]   = dev[t][1];
//                x_hat_d[t]  = x_hat[t][2];
//                dev_d[t]    = dev[t][2];
//           }
// 
// }













// 
