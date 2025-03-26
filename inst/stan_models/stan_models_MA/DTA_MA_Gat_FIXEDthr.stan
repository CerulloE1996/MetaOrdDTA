


functions {
        ////
        //// ---- Include files to compile the necessary custom (user-defined) Stan functions:
        ////
        #include "Stan_fns_basic.stan"
        #include "Stan_fns_corr.stan"
        #include "Stan_fns_ordinal_and_cutpoints.stan"
        #include "Stan_fns_log_lik.stan"
}


data {
        ////
        //// ---- Data:
        ////
        int<lower=1> n_studies;
        int<lower=1> n_thr;
        array[n_studies] int n_obs_cutpoints;  //// OBSERVED cutpoints for test in study s
        ////
        // array[2] matrix[n_studies, n_thr] x_with_missings;
        array[2] matrix[n_studies, n_thr] x;
        array[2] matrix[n_studies, n_thr] n;
        array[2] matrix[n_studies, n_thr] cutpoint_index;
        ////
        //// ---- Priors for location ("beta"):
        ////
        real prior_beta_mu_mean;
        real prior_beta_mu_SD;
        real prior_beta_SD_mean;
        real prior_beta_SD_SD;
        ////
        //// ---- Priors for raw_scale ("gamma"):
        ////
        real prior_raw_scale_mu_mean;
        real prior_raw_scale_mu_SD;
        real prior_raw_scale_SD_mean;
        real prior_raw_scale_SD_SD;
        ////
        //// ---- Induced-Dirichlet priors:
        ////
        vector<lower=0.0>[n_thr + 1] prior_dirichlet_alpha; // Induced-Dirichlet prior vector 
        ////
        //// ---- Other:
        ////
        int<lower=0, upper=1> softplus;
  
}


transformed data { 
    
        int n_cat = n_thr + 1; //// Number of ordinal categories for index test
        vector[n_thr] Ind_Dir_anchor = rep_vector(0.0, n_thr);
        real mult_nd = (1)*(-1.0);
        real mult_d  = (1)*(+1.0); //// hence: mult_d > mult_nd
        
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
        // vector[n_thr] C_raw_vec;   //// Global cutpoints ("raw" / unconstrained)
        ordered[n_thr] C;   //// Global cutpoints
        
}
 

transformed parameters { 
  
        ////
        //// ---- Construct (global) cutpoints:
        ////
        // vector[n_thr] C = ((softplus == 1) ? construct_C_using_SP_jacobian(C_raw_vec) : construct_C_using_exp_jacobian(C_raw_vec));
        ////
        //// ---- Likelihood stuff:
        ////
        array[2] matrix[n_studies, n_thr] cumul_prob = init_array_of_matrices(n_studies, n_thr, 2, 1.0);
        array[2] matrix[n_studies, n_thr] cond_prob  = init_array_of_matrices(n_studies, n_thr, 2, 0.0);
        array[2] matrix[n_studies, n_thr] log_lik    = init_array_of_matrices(n_studies, n_thr, 2, 0.0);
        {
                matrix[2, n_studies] locations; // local
                matrix[2, n_studies] scales;    // local
                ////
                array[2] matrix[n_studies, n_thr] latent_cumul_prob = init_array_of_matrices(n_studies, n_thr, 2, positive_infinity()); // local
                {
                    ////
                    //// ---- Between-study model for location and scale:
                    ////
                    for (s in 1:n_studies) {
                        ////
                        real raw_beta_baseline  = beta_mu + beta_SD * beta_z[s];
                        real raw_beta_baseline_nd = mult_nd*raw_beta_baseline;
                        real raw_beta_baseline_d  = mult_d*raw_beta_baseline;
                        locations[1, s] = raw_beta_baseline_nd;
                        locations[2, s] = raw_beta_baseline_d;
                        ////
                        real raw_scale_baseline = raw_scale_mu + raw_scale_SD * raw_scale_z[s]; // "gamma"
                        real raw_scale_baseline_nd = mult_nd*raw_scale_baseline; // "gamma_nd"
                        real raw_scale_baseline_d  = mult_d*raw_scale_baseline; // "gamma_d"
                        scales[1, s] =  ((softplus == 1) ? softplus_scaled(raw_scale_baseline_nd) : exp_jacobian(raw_scale_baseline_nd));
                        scales[2, s] =  ((softplus == 1) ? softplus_scaled(raw_scale_baseline_d)  : exp_jacobian(raw_scale_baseline_d));
                        ////
                    }
                }
                ////
                //// ---- Get the cutpoint index (k) to map "latent_cumul_prob[c][s, cut_i]" to correct cutpoint "C[k]":
                ////
                latent_cumul_prob = map_latent_cumul_prob_to_fixed_C(C, locations, scales, n_studies, n_obs_cutpoints, cutpoint_index);
                ////
                //// ---- Calculate CUMULATIVE probabilities (vectorised):
                ////
                for (c in 1:2) {
                    cumul_prob[c] = Phi_approx(latent_cumul_prob[c]); //// INCREASING sequence (as C_k > C_{k - 1})
                }
                // ////
                // //// ---- Multinomial (factorised binomial likelihood)
                // ////
                log_lik = compute_log_lik_binomial_fact(cumul_prob, x, n, n_obs_cutpoints);
                // array[2, 2] matrix[n_studies, n_thr] log_lik_outs;
                // log_lik_outs = compute_log_lik_binomial_fact(cumul_prob, x, n, n_obs_cutpoints);
                // log_lik   = log_lik_outs[1];
                // cond_prob = log_lik_outs[2];
        }
      
}


model {
      
        //// Priors:
        {
            //// locations:
            beta_mu ~ normal(prior_beta_mu_mean, prior_beta_mu_SD);
            beta_SD ~ normal(prior_beta_SD_mean, prior_beta_SD_SD);
            //// scales:
            raw_scale_mu ~ normal(prior_raw_scale_mu_mean, prior_raw_scale_mu_SD);
            raw_scale_SD ~ normal(prior_raw_scale_SD_mean, prior_raw_scale_SD_SD);
        }
        ////
        //// Induced-dirichlet ** Prior ** model:
        ////
        {
           vector[n_thr] Ind_Dir_cumul  = C; //- Ind_Dir_anchor;
           vector[n_thr] Ind_Dir_cumul_prob = Phi_approx(Ind_Dir_cumul);
           vector[n_cat] Ind_Dir_ord_prob = cumul_probs_to_ord_probs(Ind_Dir_cumul_prob);
           ////
           vector[n_thr] rho = std_normal_approx_pdf(Ind_Dir_cumul[1:n_thr], Ind_Dir_cumul_prob[1:n_thr]);
           Ind_Dir_ord_prob[1:n_cat] ~ induced_dirichlet_given_rho(rho, prior_dirichlet_alpha);
        }
        ////
        //// Likelihood / Model:
        ////
        target += std_normal_lpdf(to_vector(beta_z)); // part of between-study model, NOT prior
        target += std_normal_lpdf(to_vector(raw_scale_z)); // part of between-study model, NOT prior
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
          real location_nd_mu = mult_nd*beta_mu;
          real location_d_mu  = mult_d*beta_mu;
          ////
          real raw_scale_nd_mu = mult_nd*raw_scale_mu;
          real raw_scale_d_mu  = mult_d*raw_scale_mu;
          real scale_nd_mu =  ((softplus == 1) ? softplus_scaled(raw_scale_nd_mu) : exp(raw_scale_nd_mu));
          real scale_d_mu  =  ((softplus == 1) ? softplus_scaled(raw_scale_d_mu)  : exp(raw_scale_d_mu));
          ////
          vector[n_thr] Fp = 1.0 - Phi_approx((C - location_nd_mu)/scale_nd_mu);
          vector[n_thr] Sp = 1.0 - Fp;
          vector[n_thr] Se = 1.0 - Phi_approx((C - location_d_mu)/scale_d_mu);
          ////
          //// Calculate predictive accuracy:
          ////
          real beta_pred      = normal_rng(beta_mu, beta_SD);
          real location_nd_pred = mult_nd*beta_pred;
          real location_d_pred  = mult_d*beta_pred;
          ////
          real raw_scale_pred = normal_rng(raw_scale_mu, raw_scale_SD); //// raw_scale_mu;
          real raw_scale_nd_pred = mult_nd*beta_pred;
          real raw_scale_d_pred  = mult_d*beta_pred;
          real scale_nd_pred = ((softplus == 1) ? softplus_scaled(raw_scale_nd_pred) : exp(raw_scale_nd_pred));
          real scale_d_pred  = ((softplus == 1) ? softplus_scaled(raw_scale_d_pred)  : exp(raw_scale_d_pred));
          ////
          vector[n_thr] Fp_pred =   1.0 - Phi_approx((C - location_nd_pred)/scale_nd_pred);
          vector[n_thr] Sp_pred =   1.0 - Fp_pred;
          vector[n_thr] Se_pred =   1.0 - Phi_approx((C - location_d_pred)/scale_d_pred);
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
          {
              array[2] matrix[n_studies, n_thr] x_hat = init_array_of_matrices(n_studies, n_thr, 2, -1);
              array[2] matrix[n_studies, n_thr] dev   = init_array_of_matrices(n_studies, n_thr, 2, -1);
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
       
}


















