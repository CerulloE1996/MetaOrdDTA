


functions {
        ////
        //// ---- Include files to compile the necessary custom (user-defined) Stan functions:
        ////
        #include "Stan_fns_basic.stan"
        #include "Stan_fns_corr.stan"
        #include "Stan_fns_ordinal.stan"
        #include "Stan_fns_log_lik.stan"
        #include "Stan_fns_simplex.stan"
        #include "Stan_fns_Jacobian.stan"
        #include "Stan_fns_model_fit.stan"
}


data {
        ////
        //// ---- Data:
        ////
        int<lower=1> n_studies;
        int<lower=1> n_thr;
        array[n_studies] int n_obs_cutpoints;
        ////
        //// ---- Covariates:
        ////
        int n_covariates; //// For R+G - cannot have SEPERATE covariates for D+ and D-, so only one n_covariates!
        matrix[n_studies, n_covariates] X;
        vector[n_covariates] baseline_case; // must be user-inputted - e.g. could be {0, 1, 45.3} for 2 binary covariates and 1 cts one (e.g. age)
        ////
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
        vector[n_thr + 1] prior_dirichlet_phi;
        real<lower=0.0> prior_dirichlet_prob_SDs;
        ////
        //// ---- Other:
        ////
        real<lower=0.0> kappa_lb;
        int<lower=0, upper=1> softplus;
}


transformed data {  
        real mult_nd = -1.0;
        real mult_d  = +1.0; //// hence: mult_d > mult_nd
        ////
        int n_cat = n_thr + 1; //// Number of ordinal categories for index test
        ////
        int use_probit_link = 1;
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
                       for (k in 2:n_obs_cutpoints[s]) {
                                  n[c][s, k] = x_2[c][s, k - 1];
                       }
                  }
       }
       ////
       int summary_method = 2;
       ////
       array[n_studies] int dummy_ind_test_in_study;
       for (s in 1:n_studies) {
          dummy_ind_test_in_study[s] = 1;
       }
}
  

parameters {
        vector[n_covariates] beta_mu;  
        real<lower=0.0> beta_SD;   
        vector[n_studies] beta_z;
        ////
        vector[n_covariates] raw_scale_mu;    
        real<lower=0.0> raw_scale_SD;   
        vector[n_studies] raw_scale_z;
        ////
        matrix[n_studies, n_thr] C_raw;
        ////
        // real<lower=(kappa_lb)> kappa;
        simplex[n_cat] dirichlet_phi;
        ////
        vector<lower=0>[n_cat] target_sds;
}


transformed parameters { 
        vector<lower=0.0>[n_cat] alpha; 
        real<lower=(kappa_lb)> kappa;
        vector<lower=0.0>[n_cat] dir_cat_SDs_sigma;
        ////
        real avg_sd;
        {
            ////
            //// ---- Use a simple approximation based on the average SD:
            ////
            avg_sd = mean(target_sds);
            kappa = max_reals(kappa_lb, (1.0/square(avg_sd)) - 1.0 ); ////  For a symmetric Dirichlet, SD ≈ sqrt(1/(α_0+1)):
            alpha = dirichlet_phi * kappa; //// Set alpha using the approximated kappa:
            ////
            //// ---- Calculate the actual SDs:
            ////
            real alpha_0 = sum(alpha);
            dir_cat_SDs_sigma = sqrt((alpha .* (alpha_0 - alpha)) / (square(alpha_0) * (alpha_0 + 1.0)));
        }
        ////
        //// ---- Construct study-specific cutpoints:
        ////
        matrix[n_studies, n_thr] C; ////= convert_C_array_to_mat(C_array);
        for (s in 1:n_studies) {
                     C[s, ] =  construct_C(C_raw[s, 1:n_thr], softplus);
        }
        ////
        //// ---- Between-study model for location and scale:
        ////
        vector[n_studies] beta_random;
        vector[n_studies] raw_scale_random;
        for (s in 1:n_studies) {
            beta_random[s]      = beta_SD      * beta_z[s];
            raw_scale_random[s] = raw_scale_SD * raw_scale_z[s];
        }
        ////
        matrix[n_studies, 2] Xlocations; // NOT local
        matrix[n_studies, 2] Xraw_scale; // NOT local
        matrix[n_studies, 2] Xscale; // NOT local
        ////
        for (s in 1:n_studies) {
            real Xbeta_baseline  = sum(to_vector(X[s, 1:n_covariates]) .* beta_mu[1:n_covariates]) + beta_random[s];
            Xlocations[s, 1] = mult_nd*Xbeta_baseline; 
            Xlocations[s, 2] = mult_d*Xbeta_baseline; 
            ////
            real Xraw_scale_baseline = sum(to_vector(X[s, 1:n_covariates]) .* raw_scale_mu[1:n_covariates]) + raw_scale_random[s]; // "Xgamma"
            Xscale[s, 1] =  ((softplus == 1) ? softplus_scaled( mult_nd*Xraw_scale_baseline) : exp( mult_nd*Xraw_scale_baseline)); 
            Xscale[s, 2] =  ((softplus == 1) ? softplus_scaled( mult_d*Xraw_scale_baseline)  : exp( mult_d*Xraw_scale_baseline));
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
        {
              // kappa[c] ~ normal(prior_kappa_mean, prior_kappa_SD);  // ID prior
              dirichlet_phi ~ dirichlet(prior_dirichlet_phi); // ID prior
              // if (prior_only == 0) target +=  n_cat * log(kappa[c]);
              ////
              target_sds ~ normal(0.0, prior_dirichlet_prob_SDs);
              target += -sum(square(dir_cat_SDs_sigma - target_sds));
              ////
              //// Jacobian adjustment for the transformation from target_sds to kappa:
              ////
              //// Avoid including the Jacobian adjustment if we hit the max_reals constraint
              if (kappa > kappa_lb) {
                   target += log(2.0) - 3.0 * log(avg_sd);
                   target += - log(n_cat);
              }
        }
        ////
        //// ---- Induced-dirichlet between study ** model ** (NOT a prior model here but part of the actual likelihood since random-effect cutpoints!):
        ////
        for (s in 1:n_studies) {
             vector[n_cat] Ind_Dir_ord_prob;
             if (use_probit_link == 1) Ind_Dir_ord_prob = to_vector(cumul_probs_to_ord_probs(Phi(C[s, ])));
             else                      Ind_Dir_ord_prob = to_vector(cumul_probs_to_ord_probs(inv_logit(C[s, ])));
             ////
             Ind_Dir_ord_prob ~ induced_dirichlet_given_C(to_vector(C[s, ]), alpha, use_probit_link);
        }
        ////
        //// ---- Likelihood / Model:
        ////
        target += std_normal_lpdf(to_vector(beta_z)); // part of between-study model, NOT prior
        target += std_normal_lpdf(to_vector(raw_scale_z)); // part of between-study model, NOT prior
        // target += raw_scale_to_scale_log_det_J_lp(raw_scale, softplus); // no Jacobian needed for R&G models as it cancels out (since + then -).
        ////
        //// ---- Increment the log-likelihood:
        ////
        target += raw_C_to_C_log_det_J_lp(C_raw, softplus);
        ////
        //// ---- Log-likelihood:
        ////
        {
            ////
            //// ---- Get the cutpoint index (k) to map "latent_surv[c][s, cut_i]" to correct cutpoint "C[k]":
            ////
            array[2] matrix[n_studies, n_thr] latent_surv = map_latent_surv_prob_to_random_C(C, n_thr, Xlocations, Xscale, n_studies, n_obs_cutpoints, cutpoint_index);
            target += compute_log_lik_binomial_fact_lp(latent_surv, use_probit_link, n_thr, x_2, n, N_total, n_obs_cutpoints, dummy_ind_test_in_study);
        }
}


generated quantities {
      ////
      //// SD's for each category probability in a Dirichlet is √(α_k(α_0-α_k)/(α_0^2(α_0+1))) - where α_0 is the sum of all alphas:
      ////
      vector<lower=0.0>[n_cat] dirichlet_cat_SDs_sigma;
      {
          real alpha_0 = sum(alpha);
          dirichlet_cat_SDs_sigma = sqrt((alpha .* (alpha_0 - alpha) ) / (square(alpha_0) * (alpha_0 + 1.0)));
      }
      ////
      //// ---- Calculate summary accuracy (using mean parameters):
      ////
      vector[n_thr] C_mu;
      {
         if (summary_method == 1) { // this method produces SKEWED MEAN CUTPOINTS - DONT USE
           
            vector[n_cat] prob_ord_mu   = dirichlet_phi; // alpha[c] / sum(alpha[c]);  //// dirichlet_cat_means_phi[c];
            vector[n_thr] prob_cumul_mu = ord_probs_to_cumul_probs(prob_ord_mu);
            C_mu = ID_cumul_probs_to_C(prob_cumul_mu, use_probit_link);
            
         } else if (summary_method == 2) {
           
            for (k in 1:n_thr) {
                C_mu[k] = median(C[1:n_studies, k]); // solves "skewness" issue from non-linear transformations
            }
            
         } else if (summary_method == 3) { // this method produces SKEWED MEAN CUTPOINTS - DONT USE
           
            for (k in 1:n_thr) {
                C_mu[k] = mean(C[1:n_studies, k]); // ID_cumul_probs_to_C(prob_cumul_mu, 0.0, 1.0); this method produces SKEWED MEAN CUTPOINTS - DONT USE
            }
           
         }
      }
      ////
      //// ---- Calculate summary accuracy (using mean parameters):
      ////
      vector[n_covariates] location_nd_mu = mult_nd*beta_mu;
      vector[n_covariates] location_d_mu  = mult_d*beta_mu;
      ////
      vector[n_covariates] scale_nd_mu =  ((softplus == 1) ? softplus_scaled(mult_nd*raw_scale_mu) : exp(mult_nd*raw_scale_mu));
      vector[n_covariates] scale_d_mu  =  ((softplus == 1) ? softplus_scaled(mult_d*raw_scale_mu)  : exp(mult_d*raw_scale_mu));
      ////
      real Xbeta_baseline = dot_product(baseline_case, beta_mu);
      real Xbeta_baseline_nd = mult_nd * Xbeta_baseline;
      real Xbeta_baseline_d  = mult_d * Xbeta_baseline;
      ////
      real Xraw_scale_baseline = dot_product(baseline_case, raw_scale_mu);
      real scale_nd_baseline = ((softplus == 1) ? softplus_scaled(mult_nd * Xraw_scale_baseline) : exp(mult_nd * Xraw_scale_baseline));
      real scale_d_baseline  = ((softplus == 1) ? softplus_scaled(mult_d  * Xraw_scale_baseline) : exp(mult_d  * Xraw_scale_baseline));
      ////
      //// ---- Calculate baseline Se/Sp:
      ////
      vector[n_thr] Fp_baseline = (use_probit_link == 1) ? Phi(-(C_mu - Xbeta_baseline_nd)/scale_nd_baseline)  : inv_logit(-(C_mu - Xbeta_baseline_nd)/scale_nd_baseline);
      vector[n_thr] Sp_baseline = 1.0 - Fp_baseline;
      vector[n_thr] Se_baseline = (use_probit_link == 1) ? Phi(-(C_mu - Xbeta_baseline_d)/scale_d_baseline)    : inv_logit(-(C_mu - Xbeta_baseline_d)/scale_d_baseline);
      ////
      //// ---- Calculate predictive accuracy:
      ////
      vector[n_thr] C_pred;
      {
        vector[n_cat] prob_ord_pred   = dirichlet_rng(alpha); //// Simulate from Dirichlet by using the summary "alpha" parameters.
        vector[n_thr] prob_cumul_pred = ord_probs_to_cumul_probs(prob_ord_pred);  //// Compute PREDICTED cumulative probabilities.
        C_pred = ID_cumul_probs_to_C(prob_cumul_pred, use_probit_link);  //// Compute PREDICTED cutpoints.
      }
      ////
      real beta_random_pred      = normal_rng(0.0, beta_SD);
      real raw_scale_random_pred = normal_rng(0.0, raw_scale_SD);
      ////
      real Xbeta_baseline_pred = dot_product(baseline_case, beta_mu) + beta_random_pred;
      real Xbeta_baseline_pred_nd = mult_nd * Xbeta_baseline_pred;
      real Xbeta_baseline_pred_d  = mult_d  * Xbeta_baseline_pred;
      ////
      real Xraw_scale_baseline_pred = dot_product(baseline_case, raw_scale_mu) + raw_scale_random_pred;
      real scale_baseline_pred_nd = ((softplus == 1) ? softplus_scaled(mult_nd * Xraw_scale_baseline_pred) : exp(mult_nd * Xraw_scale_baseline_pred));
      real scale_baseline_pred_d  = ((softplus == 1) ? softplus_scaled(mult_d  * Xraw_scale_baseline_pred) : exp(mult_d  * Xraw_scale_baseline_pred));
      ////
      vector[n_thr] Fp_baseline_pred = (use_probit_link == 1) ? Phi(-(C_pred - Xbeta_baseline_pred_nd)/scale_baseline_pred_nd) : inv_logit(-(C_pred - Xbeta_baseline_pred_nd)/scale_baseline_pred_nd);
      vector[n_thr] Sp_baseline_pred = 1.0 - Fp_baseline_pred;
      vector[n_thr] Se_baseline_pred = (use_probit_link == 1) ? Phi(-(C_pred - Xbeta_baseline_pred_d)/scale_baseline_pred_d)   : inv_logit(-(C_pred - Xbeta_baseline_pred_d)/scale_baseline_pred_d);
      ////
      //// ---- Log-lik + study-specific accuracy computation (using "data" / double-precision fn for efficiency):
      ////
      array[2] matrix[n_studies, n_thr] log_lik; // global (for e.g. LOO)
      ////
      matrix[n_studies, n_thr] fp; // global
      matrix[n_studies, n_thr] sp; // global
      matrix[n_studies, n_thr] se; // global
      ////
      vector[n_studies] deviance_nd; // global
      vector[n_studies] deviance_d;  // global
      vector[n_studies] deviance;    // global
      {
          real scale = 1.0; // since using "Xu-like"" param.
          array[2] matrix[n_studies, n_thr] latent_surv = map_latent_surv_prob_to_random_C(C, n_thr, Xlocations, Xscale, n_studies, n_obs_cutpoints, cutpoint_index);
          ////
          array[3, 2] matrix[n_studies, n_thr] outs = compute_log_lik_binomial_fact_data(latent_surv, use_probit_link, n_thr, x_2, n, N_total, n_obs_cutpoints);
          ////
          log_lik   = outs[1];
          array[2] matrix[n_studies, n_thr] cond_prob = outs[2];
          array[2] matrix[n_studies, n_thr] surv_prob = outs[3];
          ////
          fp = surv_prob[1];
          sp = 1.0 - fp;
          se = surv_prob[2];
          ////
          //// ---- Model fit (deviance):
          ////
          array[4] matrix[n_studies, n_thr] outs_model_fit = compute_deviance(cond_prob, n_thr, x_2, n, n_obs_cutpoints);
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

}














