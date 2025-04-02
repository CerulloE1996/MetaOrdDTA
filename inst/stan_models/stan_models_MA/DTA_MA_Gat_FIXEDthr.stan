


functions {
        ////
        //// ---- Include files to compile the necessary custom (user-defined) Stan functions:
        ////
        #include "Stan_fns_basic.stan"
        #include "Stan_fns_corr.stan"
        #include "Stan_fns_ordinal_and_cutpoints.stan"
        #include "Stan_fns_log_lik.stan"
        #include "Stan_fns_Jacobian.stan"
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
        vector[n_thr] C_raw_vec;   //// Global cutpoints ("raw" / unconstrained)
}
 

transformed parameters { 
        real mult_nd = -1.0;
        real mult_d  = +1.0; //// hence: mult_d > mult_nd
        ////
        //// ---- Construct (global) cutpoints:
        ////
        vector[n_thr] C = construct_C(C_raw_vec, softplus);
        array[2] matrix[n_studies, n_thr] cumul_prob = init_array_of_matrices(n_studies, n_thr, 2, 1.0);
        array[2] matrix[n_studies, n_thr] cond_prob  = init_array_of_matrices(n_studies, n_thr, 2, 1.0);
        array[2] matrix[n_studies, n_thr] log_lik    = init_array_of_matrices(n_studies, n_thr, 2, 0.0);
                    ////
                    //// ---- Between-study model for location and scale:
                    ////
                    {
                        matrix[n_studies, 2] locations = rep_matrix(0.0, n_studies, 2); // local
                        matrix[n_studies, 2] scales    = rep_matrix(1.0, n_studies, 2); // local
                        for (s in 1:n_studies) {
                            ////
                            real raw_beta_baseline  = beta_mu + beta_SD * beta_z[s];
                            locations[s, 1] = mult_nd*raw_beta_baseline; 
                            locations[s, 2] = mult_d*raw_beta_baseline; 
                            ////
                            real raw_scale_baseline = raw_scale_mu + raw_scale_SD * raw_scale_z[s]; // "gamma"
                            scales[s, 1] =  ((softplus == 1) ? softplus_scaled( mult_nd*raw_scale_baseline) : exp( mult_nd*raw_scale_baseline)); 
                            scales[s, 2] =  ((softplus == 1) ? softplus_scaled( mult_d*raw_scale_baseline)  : exp( mult_d*raw_scale_baseline));
                            ////
                        }
                        //// 
                        //// ---- Get the cutpoint index (k) to map "latent_cumul_prob[c][s, cut_i]" to correct cutpoint "C[k]":
                        ////
                        array[2] matrix[n_studies, n_thr] latent_cumul = map_latent_cumul_prob_to_fixed_C(C, locations, scales, n_studies, n_obs_cutpoints, cutpoint_index);
                        ////
                        //// ---- Multinomial (factorised binomial likelihood):
                        ////
                        int use_probit_link = 0;
                        array[2, 3] matrix[n_studies, n_thr] log_lik_outs = compute_log_lik_binomial_probit_fact(latent_cumul, use_probit_link, x, n, n_obs_cutpoints);
                        log_lik    = log_lik_outs[, 1];
                        cond_prob  = log_lik_outs[, 2];
                        cumul_prob = log_lik_outs[, 3];
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
           // vector[n_thr] Ind_Dir_cumul  = C; //- Ind_Dir_anchor;
           vector[n_thr] Ind_Dir_cumul_prob = Phi(C);
           vector[n_cat] Ind_Dir_ord_prob = cumul_probs_to_ord_probs(Ind_Dir_cumul_prob);
           ////
           vector[n_thr] rho = std_normal_pdf(C) //// , Ind_Dir_cumul_prob);
           Ind_Dir_ord_prob ~ induced_dirichlet_given_rho(rho, prior_dirichlet_alpha);
        }
        ////
        //// Likelihood / Model:
        ////
        target += std_normal_lpdf(to_vector(beta_z)); // part of between-study model, NOT prior
        target += std_normal_lpdf(to_vector(raw_scale_z)); // part of between-study model, NOT prior
        // target += raw_scale_to_scale_log_det_J_lp(raw_scale, softplus); // no Jacobian needed for R&G models as it cancels out (since + then -).
        ////
        //// Increment the log-likelihood:
        ////
        target += raw_C_to_C_log_det_J_lp(C_raw_vec, softplus);
        for (c in 1:2) {
          target +=  sum(log_lik[c]);
        }
}
 

generated quantities {
          ////
          //// ---- Calculate summary accuracy (using mean parameters):
          ////
          real location_nd_mu = mult_nd*beta_mu;
          real location_d_mu  = mult_d*beta_mu;
          ////
          real scale_nd_mu =  ((softplus == 1) ? softplus_scaled(mult_nd*raw_scale_mu) : exp(mult_nd*raw_scale_mu));
          real scale_d_mu  =  ((softplus == 1) ? softplus_scaled(mult_d*raw_scale_mu)  : exp(mult_d*raw_scale_mu));
          ////
          vector[n_thr] Fp = 1.0 - Phi((C - location_nd_mu)/scale_nd_mu);
          vector[n_thr] Sp = 1.0 - Fp;
          vector[n_thr] Se = 1.0 - Phi((C - location_d_mu)/scale_d_mu);
          ////
          //// ---- Calculate predictive accuracy:
          ////
          real beta_pred      = normal_rng(beta_mu, beta_SD);
          real location_nd_pred = mult_nd*beta_pred;
          real location_d_pred  = mult_d*beta_pred;
          ////
          real raw_scale_pred = normal_rng(raw_scale_mu, raw_scale_SD); //// raw_scale_mu;
          real scale_nd_pred = ((softplus == 1) ? softplus_scaled( mult_nd*beta_pred) : exp( mult_nd*beta_pred));
          real scale_d_pred  = ((softplus == 1) ? softplus_scaled( mult_d*beta_pred)  : exp( mult_d*beta_pred));
          ////
          vector[n_thr] Fp_pred =   1.0 - Phi((C - location_nd_pred)/scale_nd_pred);
          vector[n_thr] Sp_pred =   1.0 - Fp_pred;
          vector[n_thr] Se_pred =   1.0 - Phi((C - location_d_pred)/scale_d_pred);
          ////
          //// ---- Calculate study-specific accuracy:
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
              //// ---- Model-predicted ("re-constructed") data:
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
          
          ////
          //// ---- "Ordinal-bivariate" params:
          ////
          array[2] vector[n_thr] biv_equiv_C;
          biv_equiv_C[1] = C / scale_nd_mu;
          biv_equiv_C[2] = C / scale_d_mu;
       
}


















