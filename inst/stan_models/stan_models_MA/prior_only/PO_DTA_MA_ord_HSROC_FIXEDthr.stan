


functions {
        ////
        //// Include files to compile the necessary custom (user-defined) Stan functions:
        ////
        #include "Stan_fns_basic.stan"
        #include "Stan_fns_Pinkney_corr.stan"
        #include "Stan_fns_ordinal_and_cutpoints.stan"
}


data {
    
        int<lower=1> n_studies;
        int<lower=1> n_thr;
        array[n_studies] int n_obs_cutpoints;  //// OBSERVED cutpoints for test in study s
        // ////
        array[2] matrix[n_studies, n_thr] cutpoint_index;
        //// Priors for beta:
        real prior_beta_mu_mean;
        real prior_beta_mu_SD;
        real prior_beta_SD_mean;
        real prior_beta_SD_SD;
        //// Priors for raw_scale:
        real prior_raw_scale_mu_mean;
        real prior_raw_scale_mu_SD;
        real prior_raw_scale_SD_mean;
        real prior_raw_scale_SD_SD;
        //// Induced-Dirichlet priors:
        vector<lower=0.0>[n_thr + 1] prior_alpha; // Induced-Dirichlet prior vector 
        //// Other:
        int use_softplus_for_scales;
  
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
        ordered[n_thr] C;  // Global cutpoints
        
}


transformed parameters { 
  
        matrix[2, n_studies] locations;
        matrix[2, n_studies] scales;
        ////
        //// Induced-Dirichlet ** prior model ** stuff:
        ////
        vector[n_thr] Ind_Dir_anchor = rep_vector(0.0, n_thr);
        vector[n_thr] Ind_Dir_cumul  = C - Ind_Dir_anchor;
        vector[n_thr] Ind_Dir_cumul_prob = Phi(Ind_Dir_cumul);
        vector[n_cat] Ind_Dir_ord_prob = cumul_probs_to_ord_probs(Ind_Dir_cumul_prob);
        ////
        array[2] matrix[n_studies, n_thr] latent_cumul_prob = init_array_of_matrices(n_studies, n_thr, 2, positive_infinity());
        array[2] matrix[n_studies, n_thr] cumul_prob        = init_array_of_matrices(n_studies, n_thr, 2, 1.0);
        array[2] matrix[n_studies, n_thr] cond_prob         = init_array_of_matrices(n_studies, n_thr, 2, 0.0);
        array[2] matrix[n_studies, n_thr] log_lik           = init_array_of_matrices(n_studies, n_thr, 2, 0.0);
        ////
        //// Between-study model for location and scale:
        ////
        {
            ////
            for (s in 1:n_studies) {
                real raw_beta_baseline  = beta_mu      + beta_SD      * beta_z[s];
                real raw_scale_baseline = raw_scale_mu + raw_scale_SD * raw_scale_z[s];
                ////
                locations[1, s] = (-1)*(-0.5)*raw_beta_baseline;
                locations[2, s] = (-1)*(+0.5)*raw_beta_baseline;
                ////
                scales[1, s] =  ((use_softplus_for_scales == 1) ? softplus_scaled(-0.5*raw_scale_baseline) : exp(-0.5*raw_scale_baseline));
                scales[2, s] =  ((use_softplus_for_scales == 1) ? softplus_scaled(-0.5*raw_scale_baseline) : exp(-0.5*raw_scale_baseline));
                //// Jacobian for raw_scale -> scale:
                //// No Jacobian needed as due to the +/- signs they cancel each other out!
                //// Jacobian_raw_scale_to_scale += +log(2) - 0.5*raw_scale_baseline; //// deriv of exp(-0.5*raw_scale_baseline) = -0.5*exp(-0.5*raw_scale_baseline) -> log_abs_deriv = +log(2) -0.5*raw_scale_baseline; //// log_inv_logit(raw_scale_nd);
                //// Jacobian_raw_scale_to_scale += +log(2) + 0.5*raw_scale_baseline; //// deriv of exp(+0.5*raw_scale_baseline) = +0.5*exp(+0.5*raw_scale_baseline) -> log_abs_deriv = -log(2) +0.5*raw_scale_baseline; //// log_inv_logit(raw_scale_nd);
            }
        }
        ////
        //// Likelihood using binomial factorization:
        ////
        for (s in 1:n_studies) {
                  for (c in 1:2) {
                      for (cut_i in 1:to_int(n_obs_cutpoints[s])) {
                             int k = to_int(cutpoint_index[c][s, cut_i]);
                             latent_cumul_prob[c][s, cut_i] = (C[k] - locations[c, s])/scales[c, s];
                      }
                  }
        }
        ////
        //// Calculate CUMULATIVE probabilities (vectorised):
        ////
        for (c in 1:2) {
            cumul_prob[c] = Phi(latent_cumul_prob[c]); //// INCREASING sequence (as C_k > C_{k - 1})
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
               vector[n_thr] rho =  normal_pdf(Ind_Dir_cumul[1:n_thr], 0.0, 1.0);   //  p_cumul[k - 1] * (1.0 - p_cumul[k - 1]); // original
               target += induced_dirichlet_v2_lpdf( Ind_Dir_ord_prob[1:n_cat] | rho, prior_alpha);
        }
        ////
        //// Likelihood / Model:
        ////
        {
            to_vector(beta_z) ~ std_normal();       // (part of between-study model, NOT prior)
            to_vector(raw_scale_z) ~ std_normal();  // (part of between-study model, NOT prior)
        }
  
}
 

generated quantities {

          ////
          //// Calculate summary accuracy (using mean parameters):
          ////
          real location_nd_mu = (-1)*(-0.5)*beta_mu;
          real location_d_mu  = (-1)*(+0.5)*beta_mu;
          real scale_nd_mu = softplus_scaled(-0.5*raw_scale_mu);
          real scale_d_mu  = softplus_scaled(+0.5*raw_scale_mu);
          ////
          vector[n_thr] Fp = 1.0 - Phi((C - location_nd_mu)/scale_nd_mu);
          vector[n_thr] Sp = 1.0 - Fp;
          vector[n_thr] Se = 1.0 - Phi((C - location_d_mu)/scale_d_mu);
          ////
          //// Calculate predictive accuracy:
          ////
          real beta_pred      = normal_rng(beta_mu, beta_SD);
          real raw_scale_pred = normal_rng(raw_scale_mu, raw_scale_SD); //// raw_scale_mu;
          real location_nd_pred = (-1)*(-0.5)*beta_pred;
          real location_d_pred  = (-1)*(+0.5)*beta_pred;
          real scale_nd_pred = softplus_scaled(-0.5*raw_scale_pred);
          real scale_d_pred  = softplus_scaled(+0.5*raw_scale_pred); // if using exp(), this = MEDIAN of log-normal (could also use mean which is e.g.,: "exp(-0.5*(raw_scale_pred + 0.5*raw_scale_SD"))"
          ////
          vector[n_thr] Fp_pred =   1.0 - Phi((C - location_nd_pred)/scale_nd_pred);
          vector[n_thr] Sp_pred =   1.0 - Fp_pred;
          vector[n_thr] Se_pred =   1.0 - Phi((C - location_d_pred)/scale_d_pred);
          ////
          //// Calculate study-specific accuracy:
          ////
          array[3] matrix[n_studies, n_thr] out = compute_ss_accuracy_from_cumul_prob(cumul_prob, n_studies, n_thr);
          matrix[n_studies, n_thr] fp = out[1];
          matrix[n_studies, n_thr] sp = out[2];
          matrix[n_studies, n_thr] se = out[3];
       
}


















