


functions {
        ////
        //// ---- Include files to compile the necessary custom (user-defined) Stan functions:
        ////
        #include "Stan_fns_basic.stan"
        #include "Stan_fns_corr.stan"
        #include "Stan_fns_ordinal.stan"
        #include "Stan_fns_log_lik.stan"
        #include "Stan_fns_Jacobian.stan"
        #include "Stan_fns_model_fit.stan"
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
        array[2] matrix[n_studies, n_thr + 1] x;
        array[2] matrix[n_studies, n_thr + 1] cutpoint_index; 
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
        int use_probit_link = 1;
        ////
        real mult_nd = -1.0;
        real mult_d  = +1.0; //// hence: mult_d > mult_nd
        ////
        array[2] matrix[n_studies, n_thr] x_2;
        array[2] vector[n_studies] N_total;
        array[2] matrix[n_studies, n_thr] n;
        ////
        for (s in 1:n_studies) {
          
                  for (c in 1:2) {
                      N_total[c][s] = x[c][s, 1];
                      for (k in 1:n_thr) {
                         x_2[c][s, k] = x[c][s, k + 1];
                      }
                  }
                  
                  for (c in 1:2) {
                       n[c][s, 1] = N_total[c][s];
                       ////
                       for (k in 2:n_obs_cutpoints[s]) {
                                  n[c][s, k] = x_2[c][s, k - 1];
                       }
                  }
                    
       }
        
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
        ////
        //// ---- Construct (global) cutpoints:
        ////
        vector[n_thr] C = construct_C(C_raw_vec, softplus);
        ////
        //// ---- Between-study model for location and scale:
        ////
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
      
}


model {
        //// ---- Priors:
        beta_mu ~ normal(prior_beta_mu_mean, prior_beta_mu_SD);
        beta_SD ~ normal(prior_beta_SD_mean, prior_beta_SD_SD);
        ////
        raw_scale_mu ~ normal(prior_raw_scale_mu_mean, prior_raw_scale_mu_SD);
        raw_scale_SD ~ normal(prior_raw_scale_SD_mean, prior_raw_scale_SD_SD);
        ////
        //// ---- Induced-dirichlet ** Prior ** model:
        ////
        {
           vector[n_cat] Ind_Dir_ord_prob; //// = cumul_probs_to_ord_probs(Phi(C));
           if (use_probit_link == 1) Ind_Dir_ord_prob = cumul_probs_to_ord_probs(Phi(C));
           else                      Ind_Dir_ord_prob = cumul_probs_to_ord_probs(inv_logit(C));
           ////
           Ind_Dir_ord_prob ~ induced_dirichlet_given_C(C, prior_dirichlet_alpha, use_probit_link);
        }
        ////
        //// ---- Likelihood / Model:
        ////
        target += std_normal_lpdf(to_vector(beta_z)); // part of between-study model, NOT prior
        target += std_normal_lpdf(to_vector(raw_scale_z)); // part of between-study model, NOT prior
        // target += raw_scale_to_scale_log_det_J_lp(raw_scale, softplus); // no Jacobian needed for R&G models as it cancels out (since + then -).
        target += raw_C_to_C_log_det_J_lp(C_raw_vec, softplus);
        ////
        //// ---- Increment the log-likelihood:
        ////
        {
            ////
            //// ---- Get the cutpoint index (k) to map "latent_surv[c][s, cut_i]" to correct cutpoint "C[k]":
            ////
            real scale = 1.0; // since using "Xu-like"" param.
            array[2] matrix[n_studies, n_thr] latent_surv = map_latent_surv_prob_to_fixed_C(C, locations, scales, n_studies, n_obs_cutpoints, cutpoint_index);
            target += compute_log_lik_binomial_fact_lp(latent_surv, use_probit_link, x_2, n, N_total, n_obs_cutpoints);
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
          vector[n_thr] Fp = (use_probit_link == 1) ? Phi(-(C - location_nd_mu)/scale_nd_mu) : inv_logit(-(C - location_nd_mu)/scale_nd_mu);
          vector[n_thr] Sp = 1.0 - Fp;
          vector[n_thr] Se = (use_probit_link == 1) ? Phi(-(C - location_d_mu)/scale_d_mu)   : inv_logit(-(C - location_d_mu)/scale_d_mu);
          ////
          //// ---- Calculate predictive accuracy:
          ////
          real beta_pred        = normal_rng(beta_mu, beta_SD);
          real location_nd_pred = mult_nd*beta_pred;
          real location_d_pred  = mult_d*beta_pred;
          ////
          real raw_scale_pred = normal_rng(raw_scale_mu, raw_scale_SD); //// raw_scale_mu;
          real scale_nd_pred  = ((softplus == 1) ? softplus_scaled( mult_nd*beta_pred) : exp( mult_nd*beta_pred));
          real scale_d_pred   = ((softplus == 1) ? softplus_scaled( mult_d*beta_pred)  : exp( mult_d*beta_pred));
          ////
          // vector[n_thr] Fp_pred =   Phi(-(C - location_nd_pred)/scale_nd_pred);
          // vector[n_thr] Sp_pred =   1.0 - Fp_pred;
          // vector[n_thr] Se_pred =   Phi(-(C - location_d_pred)/scale_d_pred);
          vector[n_thr] Fp_pred = (use_probit_link == 1) ? Phi(-(C - location_nd_pred)/scale_nd_pred) : inv_logit(-(C - location_nd_pred)/scale_nd_pred);
          vector[n_thr] Sp_pred = 1.0 - Fp_pred;
          vector[n_thr] Se_pred = (use_probit_link == 1) ? Phi(-(C - location_d_pred)/scale_d_pred)   : inv_logit(-(C - location_d_pred)/scale_d_pred);
          ////
          //// ---- Log-lik + study-specific accuracy computation (using "data" / double-precision fn for efficiency):
          ////
          array[2] matrix[n_studies, n_thr] log_lik; // global (for e.g. LOO)
          ////
          // matrix[n_studies, n_thr] fp; // global
          // matrix[n_studies, n_thr] sp; // global
          // matrix[n_studies, n_thr] se; // global
          ////
          vector[n_studies] deviance_nd; // global
          vector[n_studies] deviance_d;  // global
          vector[n_studies] deviance;    // global
          {
              real scale = 1.0; // since using "Xu-like"" param.
              array[2] matrix[n_studies, n_thr] latent_surv = map_latent_surv_prob_to_fixed_C(C, locations, scales, n_studies, n_obs_cutpoints, cutpoint_index);
              ////
              array[3, 2] matrix[n_studies, n_thr] outs = compute_log_lik_binomial_fact_data(latent_surv, use_probit_link, x_2, n, N_total, n_obs_cutpoints);
              ////
              log_lik   = outs[1];
              array[2] matrix[n_studies, n_thr] cond_prob = outs[2];
              array[2] matrix[n_studies, n_thr] surv_prob = outs[3];
              ////
              // fp = surv_prob[1];
              // sp = 1.0 - fp;
              // se = surv_prob[2];
              ////
              //// ---- Model fit (deviance):
              ////
              array[4] matrix[n_studies, n_thr] outs_model_fit = compute_deviance(cond_prob, x_2, n, n_obs_cutpoints);
              ////
              matrix[n_studies, n_thr] x_hat_nd = outs_model_fit[1];
              matrix[n_studies, n_thr] dev_nd   = outs_model_fit[2];
              matrix[n_studies, n_thr] x_hat_d  = outs_model_fit[3];
              matrix[n_studies, n_thr] dev_d    = outs_model_fit[4];
              ////
              for (s in 1:n_studies) {
                  for (k in 1:to_int(n_obs_cutpoints[s])) {
                     deviance_nd[s] += dev_nd[s, k];
                     deviance_d[s]  += dev_d[s, k];
                  }
                  deviance[s] = deviance_nd[s] + deviance_d[s];
              }
          }
          ////
          //// ---- "Ordinal-bivariate equivelent" params:
          ////
          array[2] vector[n_thr] biv_equiv_C;
          biv_equiv_C[1] = C / scale_nd_mu;
          biv_equiv_C[2] = C / scale_d_mu;
       
}


















