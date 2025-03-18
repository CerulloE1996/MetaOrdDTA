


functions {
        ////
        //// Include files to compile the necessary custom (user-defined) Stan functions:
        ////
        #include "Stan_fns_basic.stan"
        #include "Stan_fns_Pinkney_corr.stan"
        #include "Stan_fns_ordinal_and_cutpoints.stan"
}


data {
    
        ////
        //// Data:
        ////
        int<lower=1> n_studies;
        int<lower=1> n_thr;
        array[n_studies] int n_obs_cutpoints;  //// OBSERVED cutpoints for test in study s
        ////
        // array[2] matrix[n_studies, n_thr] x_with_missings;
        array[2] matrix[n_studies, n_thr] n;
        array[2] matrix[n_studies, n_thr] x;
        array[2] matrix[n_studies, n_thr] cutpoint_index;
        ////
        //// Priors for location:
        ////
        real prior_beta_mu_mean;
        real prior_beta_mu_SD;
        real prior_beta_SD_mean;
        real prior_beta_SD_SD;
        ////
        //// Priors for raw_scale:
        ////
        real prior_raw_scale_mu_mean;
        real prior_raw_scale_mu_SD;
        real prior_raw_scale_SD_mean;
        real prior_raw_scale_SD_SD;
        ////
        //// Priors for cutpoints (using "Induced-Dirichlet" ordinal probs):
        ////
        vector[n_thr + 1] prior_dirichlet_cat_means_alpha;
        //// Dirichlet priors for the transformed Dirichlet "SD" params:
        vector[n_thr + 1] prior_dirichlet_cat_SDs_mean;
        vector<lower=0.0>[n_thr + 1] prior_dirichlet_cat_SDs_SD;
        ////
        //// Other:
        ////
        real kappa_lb;
        int<lower=0, upper=1> softplus;
    
  
}

transformed data { 
    
        int n_cat = n_thr + 1; //// Number of ordinal categories for index test
        
}
  


parameters {
  
        real beta_mu;    
        real<lower=0.0> beta_SD;   
        vector[n_studies] beta_z;
        ////
        real raw_scale_mu;    
        real<lower=0.0> raw_scale_SD;   
        vector[n_studies] raw_scale_z;
        ////
        array[n_studies] ordered[n_thr] C_array; // study-specific cutpoints
        simplex[n_cat] dirichlet_cat_means_phi;
        real<lower=kappa_lb> kappa;  
        
}


transformed parameters { 
  
        matrix[2, n_studies] locations;
        matrix[2, n_studies] scales;
        ////
        array[2] matrix[n_studies, n_thr] latent_cumul_prob; // Ordinal probs for the likelihood
        array[2] matrix[n_studies, n_thr] cumul_prob;  // Ordinal probs for the likelihood
        array[2] matrix[n_studies, n_thr] cond_prob;   // Ordinal probs for the likelihood
        array[2] matrix[n_studies, n_thr] log_lik;     // log_lik storage
        ////
        matrix[n_studies, n_thr] C = convert_C_array_to_mat(C_array);
        ////
        //// Induced-Dirichlet ** model ** stuff:
        ////
        matrix[n_studies, n_thr] Ind_Dir_anchor = rep_matrix(0.0, n_studies, n_thr);
        matrix[n_studies, n_thr] Ind_Dir_cumul  = C - Ind_Dir_anchor;
        matrix[n_studies, n_thr] Ind_Dir_cumul_prob = Phi(Ind_Dir_cumul);
        matrix[n_studies, n_cat] Ind_Dir_ord_prob = cumul_probs_to_ord_probs(Ind_Dir_cumul_prob);
        ////
        real log_kappa = log(kappa);
        vector[n_cat] log_dirichlet_cat_means_phi = log(dirichlet_cat_means_phi);
        vector[n_cat] log_alpha = log_kappa + log_dirichlet_cat_means_phi;
        vector[n_cat] alpha = exp(log_alpha);
        real alpha_0 = sum(alpha);
        ////
        ////  Perform "alpha" -> "dirichlet_cat_SDs_sigma" transformation:
        ////
        vector[n_cat] dirichlet_cat_SDs_sigma = convert_Dirichlet_alpha_to_SDs_jacobian(alpha);
        ////
        //// Between-study model for location and scale:
        ////
        {
            for (s in 1:n_studies) {
                
                  //// locations:
                  real raw_beta_baseline  = beta_mu + beta_SD * beta_z[s];
                  locations[1, s] = (-1)*(-0.5)*raw_beta_baseline;
                  locations[2, s] = (-1)*(+0.5)*raw_beta_baseline;
                  //// scales:
                  real raw_scale_baseline;
                  // if (fix_scale_between_studies == 1) {
                  //      raw_scale_baseline = raw_scale_mu;
                  // } else { 
                       raw_scale_baseline = raw_scale_mu + raw_scale_SD * raw_scale_z[s];
                  // }
                  ////
                  // if (fix_scale_between_studies == 1) {
                  //   
                  // } else { 
                       scales[1, s] = ((softplus == 1) ? softplus_scaled(-0.5*raw_scale_baseline) : exp_jacobian(-0.5*raw_scale_baseline));
                       scales[2, s] = ((softplus == 1) ? softplus_scaled(+0.5*raw_scale_baseline) : exp_jacobian(+0.5*raw_scale_baseline));
                  // }
            }
        }
        ////
        //// Likelihood using binomial factorization:
        ////
        for (s in 1:n_studies) {
                  for (c in 1:2) {
                      for (cut_i in 1:to_int(n_obs_cutpoints[s])) {
                             int k = to_int(cutpoint_index[c][s, cut_i]);
                             latent_cumul_prob[c][s, cut_i] = (C[s, k] - locations[c, s])/scales[c, s];
                      }
                  }
        }
        ////
        //// Calculate CUMULATIVE probabilities (vectorised):
        ////
        for (c in 1:2) {
            cumul_prob[c] = Phi(latent_cumul_prob[c]); //// INCREASING sequence (as C_k > C_{k - 1})
        }
        ////
        //// ------- Binomial likelihood:
        ////
        for (s in 1:n_studies) {
                for (c in 1:2) {
                      for (cut_i in 1:n_obs_cutpoints[s]) {
                              // Current and next cumulative counts
                              int x_current = to_int(x[c][s, cut_i]);
                              int x_next    = to_int(n[c][s, cut_i]);
                              
                              // Skip if the current count is zero (no observations to classify)
                              if (x_current != 0)  {
                              
                                    // Conditional probability of being at or below the current cutpoint - given being at or below the next cutpoint
                                    if (cut_i == n_obs_cutpoints[s]) { 
                                             cond_prob[c][s, cut_i]  = cumul_prob[c][s, cut_i] / 1.0;
                                    } else {
                                          if (x_next > 0) { 
                                             cond_prob[c][s, cut_i] = cumul_prob[c][s, cut_i] / cumul_prob[c][s, cut_i + 1];
                                          } else { 
                                             cond_prob[c][s, cut_i] = 1.0;
                                          }
                                    }
                                    
                                    // Binomial for observations at or below current cutpoint out of those at or below next cutpoint
                                    log_lik[c][s, cut_i] = binomial_lpmf(x_current | x_next, cond_prob[c][s, cut_i]);
                                    
                              }
                      }
                }
        }
      
}


model {
      
        ////
        //// Priors:
        ////
        {
            //// locations:
            beta_mu ~ normal(prior_beta_mu_mean, prior_beta_mu_SD);
            beta_SD ~ normal(prior_beta_SD_mean, prior_beta_SD_SD);
            //// scales:
            raw_scale_mu ~ normal(prior_raw_scale_mu_mean, prior_raw_scale_mu_SD);
            raw_scale_SD ~ normal(prior_raw_scale_SD_mean, prior_raw_scale_SD_SD);
        }
        ////
        //// Induced-dirichlet between study ** model ** (NOT a prior model here but part of the actual likelihood since random-effect cutpoints!):
        ////
        {
            //// target += normal_lpdf(kappa | prior_kappa_mean, prior_kappa_SD);
            target += dirichlet_lpdf( dirichlet_cat_means_phi | prior_dirichlet_cat_means_alpha ); // "flat" prior on the simplex dirichlet_cat_means_phi. 
            target += normal_lpdf(dirichlet_cat_SDs_sigma | prior_dirichlet_cat_SDs_mean, prior_dirichlet_cat_SDs_SD );
        }
        ////
        {
            for (s in 1:n_studies) {
                row_vector[n_thr] rho =  normal_pdf(Ind_Dir_cumul[s, 1:n_thr], 0.0, 1.0);
                target += induced_dirichlet_v2_lpdf(Ind_Dir_ord_prob[s, 1:n_cat] | rho, to_row_vector(alpha[1:n_cat]));
            }
            for (k in 1:n_cat) {
                target += log_kappa;
                target += log_dirichlet_cat_means_phi[k];
            }
        }
        ////
        //// Likelihood / Model:
        ////
        // target += Jacobian_raw_scale_to_scale;
        {
            to_vector(beta_z) ~ std_normal();   // (part of between-study model, NOT prior)
            to_vector(raw_scale_z) ~ std_normal();  // (part of between-study model, NOT prior)
        }
        ////
        //// Increment the log-likelihood:
        ////
        for (c in 1:2) {
          target +=  sum(log_lik[c]);
        }
  
}


generated quantities {

          ////
          //// Calculate summary accuracy (using mean parameters):
          ////
          real location_nd_mu = (-1)*(-0.5)*beta_mu;
          real location_d_mu  = (-1)*(+0.5)*beta_mu;
          real scale_nd_mu =  ((softplus == 1) ? softplus_scaled(-0.5*raw_scale_mu) : exp(-0.5*raw_scale_mu));
          real scale_d_mu  =  ((softplus == 1) ? softplus_scaled(+0.5*raw_scale_mu) : exp(+0.5*raw_scale_mu));// This corresponds to the MEDIAN of log-normal (could also use mean which is e.g.,: "scale_nd = exp(-0.5*(raw_scale_mu + 0.5*raw_scale_SD"))"
          ////
          vector[n_cat] prob_ord_mu   = dirichlet_cat_means_phi;
          vector[n_thr] prob_cumul_mu = ord_probs_to_cumul_probs(prob_ord_mu);
          vector[n_thr] C_mu = cumul_probs_to_C(prob_cumul_mu, 0.0, 1.0);   
          ////
          vector[n_thr] Fp = 1.0 - Phi((C_mu - location_nd_mu)/scale_nd_mu);
          vector[n_thr] Sp = 1.0 - Fp;
          vector[n_thr] Se = 1.0 - Phi((C_mu - location_d_mu)/scale_d_mu);
          ////
          //// Calculate predictive accuracy:
          ////
          vector[n_cat] prob_ord_pred   = dirichlet_rng(alpha); //// Simulate from Dirichlet by using the summary "alpha" parameters.
          vector[n_thr] prob_cumul_pred = ord_probs_to_cumul_probs(prob_ord_pred);  //// Compute PREDICTED cumulative probabilities.
          vector[n_thr] C_pred = cumul_probs_to_C(prob_cumul_pred, 0.0, 1.0);  //// Compute PREDICTED cutpoints.
          ////
          real beta_pred      = normal_rng(beta_mu, beta_SD);
          real raw_scale_pred = normal_rng(raw_scale_mu, raw_scale_SD); //// raw_scale_mu;
          real location_nd_pred = (-1)*(-0.5)*beta_pred;
          real location_d_pred  = (-1)*(+0.5)*beta_pred;
          real scale_nd_pred = ((softplus == 1) ? softplus_scaled(-0.5*raw_scale_pred) : exp(-0.5*raw_scale_pred));
          real scale_d_pred  = ((softplus == 1) ? softplus_scaled(+0.5*raw_scale_pred) : exp(+0.5*raw_scale_pred)); // if using exp(), this = MEDIAN of log-normal (could also use mean which is e.g.,: "exp(-0.5*(raw_scale_pred + 0.5*raw_scale_SD"))"
          ////
          vector[n_thr] Fp_pred =   1.0 - Phi((C_pred - location_nd_pred)/scale_nd_pred);
          vector[n_thr] Sp_pred =   1.0 - Fp_pred;
          vector[n_thr] Se_pred =   1.0 - Phi((C_pred - location_d_pred)/scale_d_pred);
          ////
          //// Calculate study-specific accuracy:
          ////
          matrix[n_studies, n_thr] se = 1.0 - cumul_prob[1];
          matrix[n_studies, n_thr] fp = 1.0 - fp;
          matrix[n_studies, n_thr] sp = 1.0 - cumul_prob[2];
          ////
          matrix[n_studies, n_thr] x_hat_nd = rep_matrix(-1, n_studies, n_thr);
          matrix[n_studies, n_thr] x_hat_d  = rep_matrix(-1, n_studies, n_thr);
          matrix[n_studies, n_thr] dev_nd   = rep_matrix(-1, n_studies, n_thr);
          matrix[n_studies, n_thr] dev_d    = rep_matrix(-1, n_studies, n_thr);
          array[2] matrix[n_studies, n_thr] x_hat = init_array_of_matrices(n_studies, n_thr, 2, -1);
          array[2] matrix[n_studies, n_thr] dev   = init_array_of_matrices(n_studies, n_thr, 2, -1);
          ////
          vector[n_thr] C_MU_empirical = rowMedians(C);  //// Empirical-mean cutpoints:
          ////
          //// Model-predicted ("re-constructed") data:
          ////
          for (s in 1:n_studies) {
              for (c in 1:2) {
                 for (cut_i in 1:to_int(n_obs_cutpoints[s])) {
                  
                      //// Model-estimated data:
                      x_hat[c][s, cut_i] = cond_prob[c][s, cut_i] * n[c][s, cut_i];  	 // Fitted values
                    
                      //// Compute residual deviance contribution:
                      real n_i =  (n[c][s, cut_i]);
                      real x_i =  (x[c][s, cut_i]);
                      real x_hat_i =  (x_hat[c][s, cut_i]);
                      real log_x_minus_log_x_hat = log(x_i) - log(x_hat_i);
                      real log_diff_n_minus_x = log(n_i - x_i);
                      real log_diff_n_minus_x_hat = log(abs(n_i - x_hat_i));
                       
                      dev[c][s, cut_i] = 2.0 * ( x_i * log_x_minus_log_x_hat + (n_i - x_i) * (log_diff_n_minus_x - log_diff_n_minus_x_hat) ); 
                      
                 }
              }
          }
          
         x_hat_nd = x_hat[1];
         dev_nd = dev[1];
         x_hat_d = x_hat[2];
         dev_d = dev[2];
       
}














