

functions {
        ////
        //// Include files to compile the necessary custom (user-defined) Stan functions:
        ////
        #include "Stan_fns_basic.stan"
        #include "Stan_fns_corr.stan"
        #include "Stan_fns_ordinal_and_cutpoints.stan"
}


data {

        ////
        //// Data:
        ////
        int<lower=1> n_studies;
        int<lower=1> n_thr;
        array[n_studies] int n_obs_cutpoints;
        ////
        // array[2] matrix[n_studies, n_thr] x_with_missings;
        array[2] matrix[n_studies, n_thr] n;
        array[2] matrix[n_studies, n_thr] x;
        array[2] matrix[n_studies, n_thr] cutpoint_index;
        ////
        //// Priors for locations:
        ////
        vector[2] prior_beta_mu_mean;
        vector[2] prior_beta_mu_SD;
        vector[2] prior_beta_SD_mean;
        vector[2] prior_beta_SD_SD;
        ////
        //// Priors for between-study correlation matrix:
        ////
        real<lower=-1.0> beta_corr_lb;
        real<lower=beta_corr_lb, upper=1.0>  beta_corr_ub;
        real<lower=0.0>  prior_beta_corr_LKJ;
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
  
        vector[2] beta_mu;    
        vector<lower=0.0>[2] beta_SD;   
        matrix[2, n_studies] beta_z;    //// Study-specific random effects (off-centered parameterisation)
        ////
        real beta_corr;  //// between-study corr (possibly restricted)
        ////
        array[n_studies] ordered[n_thr] C_array; //// study-specific cutpoints
        simplex[n_cat] dirichlet_cat_means_phi;
        real<lower=kappa_lb> kappa;  
        
}


transformed parameters { 
    
        matrix[2, n_studies] beta;
        ////
        //// Construct simple 2x2 (bivariate) between-study corr matrices:
        ////
        cholesky_factor_corr[2] beta_L_Omega = make_restricted_bivariate_L_Omega_jacobian(beta_corr, beta_corr_lb, beta_corr_ub);
        cholesky_factor_cov[2]  beta_L_Sigma = diag_pre_multiply(beta_SD, beta_L_Omega);
        ////
        //// Initialise matrices:
        ////
        array[2] matrix[n_studies, n_thr] logit_cumul_prob = init_array_of_matrices(n_studies, n_studies, 2, positive_infinity());
        array[2] matrix[n_studies, n_thr] cumul_prob       = init_array_of_matrices(n_studies, n_studies, 2, 1.0);
        array[2] matrix[n_studies, n_thr] cond_prob        = init_array_of_matrices(n_studies, n_studies, 2, 1.0);
        array[2] matrix[n_studies, n_thr] log_lik          = init_array_of_matrices(n_studies, n_studies, 2, 0.0);
        ////
        //// Cutpoints + Induced-Dirichlet ** model ** stuff:
        ////
        matrix[n_studies, n_thr] C = convert_C_array_to_mat(C_array);
        matrix[n_studies, n_thr] Ind_Dir_anchor = rep_matrix(0.0, n_studies, n_thr);
        matrix[n_studies, n_thr] Ind_Dir_cumul  = C - Ind_Dir_anchor;
        matrix[n_studies, n_thr] Ind_Dir_cumul_prob = Phi(Ind_Dir_cumul);
        matrix[n_studies, n_cat] Ind_Dir_ord_prob = cumul_probs_to_ord_probs(Ind_Dir_cumul_prob);
        ////
        real log_kappa = log(kappa);
        vector[n_cat] log_dirichlet_cat_means_phi = log(dirichlet_cat_means_phi);
        vector[n_cat] log_alpha = log_kappa + log_dirichlet_cat_means_phi;
        vector[n_cat] alpha = exp(log_alpha);
        ////
        //// Perform "alpha" -> "dirichlet_cat_SDs_sigma" transformation:
        ////
        vector[n_cat] dirichlet_cat_SDs_sigma = convert_Dirichlet_alpha_to_SDs_jacobian(alpha);
        ////
        //// Between-study model for the location parameters ("beta") - models between-study correlation:
        ////
        for (s in 1:n_studies) {
             beta[, s] = beta_mu + beta_L_Sigma * beta_z[, s];
        }
        ////
        //// Likelihood using binomial factorization:
        ////
        for (s in 1:n_studies) {
                  for (c in 1:2) {
                      for (cut_i in 1:to_int(n_obs_cutpoints[s])) {
                             int k = to_int(cutpoint_index[c][s, cut_i]);
                             logit_cumul_prob[c][s, cut_i] = (C[s, k] - beta[c, s]);
                      }
                  }
        }
        ////
        //// Calculate CUMULATIVE probabilities (vectorised):
        ////
        for (c in 1:2) {
            cumul_prob[c] = Phi(logit_cumul_prob[c]); //// INCREASING sequence (as C_k > C_{k - 1})
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
        //// Priors for locations:
        ////
        beta_mu ~ normal(prior_beta_mu_mean, prior_beta_mu_SD);
        beta_SD ~ normal(prior_beta_SD_mean, prior_beta_SD_SD);
        beta_Omega ~ lkj_corr(prior_beta_corr_LKJ);
       
        //// Induced-dirichlet between study ** model ** (NOT a prior model here but part of the actual likelihood since random-effect cutpoints!):
        {
            //// target += normal_lpdf(kappa | prior_kappa_mean, prior_kappa_SD);
            target += dirichlet_lpdf( dirichlet_cat_means_phi | prior_dirichlet_cat_means_alpha ); // "flat" prior on the simplex dirichlet_cat_means_phi. 
            target += normal_lpdf( dirichlet_cat_SDs_sigma | prior_dirichlet_cat_SDs_mean, prior_dirichlet_cat_SDs_SD );
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
        {
            to_vector(beta_z) ~ std_normal();   // (part of between-study model, NOT prior)
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
          vector[n_cat] prob_ord_mu   = dirichlet_cat_means_phi;
          vector[n_thr] prob_cumul_mu = ord_probs_to_cumul_probs(prob_ord_mu);
          vector[n_thr] C_mu = cumul_probs_to_C(prob_cumul_mu, 0.0, 1.0);
          //// 
          vector[n_thr] Fp = 1.0 - Phi(C_mu - beta_mu[1]);
          vector[n_thr] Sp = 1.0 - Fp;
          vector[n_thr] Se = 1.0 - Phi(C_mu - beta_mu[2]);
          ////
          //// Calculate predictive accuracy:
          ////
          vector[n_cat] prob_ord_pred   = dirichlet_rng(alpha); //// Simulate from Dirichlet by using the summary "alpha" parameters.
          vector[n_thr] prob_cumul_pred = ord_probs_to_cumul_probs(prob_ord_pred);  //// Compute PREDICTED cumulative probabilities.
          vector[n_thr] C_pred = cumul_probs_to_C(prob_cumul_pred, 0.0, 1.0);  //// Compute PREDICTED cutpoints.
          ////
          vector[2] beta_pred =  multi_normal_cholesky_rng(beta_mu, beta_L_Sigma);
          ////
          vector[n_thr] Fp_pred = 1.0 - Phi(C_pred - beta_pred[1]);
          vector[n_thr] Sp_pred = 1.0 - Fp_pred;
          vector[n_thr] Se_pred = 1.0 - Phi(C_pred - beta_pred[2]);
          ////
          //// Calculate study-specific accuracy:
          ////
          array[3] matrix[n_studies, n_thr] out = compute_ss_accuracy_from_cumul_prob(cumul_prob, n_studies, n_thr);
          matrix[n_studies, n_thr] fp = out[1];
          matrix[n_studies, n_thr] sp = out[2];
          matrix[n_studies, n_thr] se = out[3];
          ////
          matrix[n_studies, n_thr] x_hat_nd = rep_matrix(-1, n_studies, n_thr);
          matrix[n_studies, n_thr] x_hat_d  = rep_matrix(-1, n_studies, n_thr);
          matrix[n_studies, n_thr] dev_nd   = rep_matrix(-1, n_studies, n_thr);
          matrix[n_studies, n_thr] dev_d    = rep_matrix(-1, n_studies, n_thr);
          ////
          array[2] matrix[n_studies, n_thr] x_hat = init_array_of_matrices(n_studies, n_thr, 2, -1);
          array[2] matrix[n_studies, n_thr] dev   = init_array_of_matrices(n_studies, n_thr, 2, -1);
          ////
          corr_matrix[2] beta_Omega = multiply_lower_tri_self_transpose(beta_L_Omega); //// Compute between-study correlation matrix for location parameters:
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














