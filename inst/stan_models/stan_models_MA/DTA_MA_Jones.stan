


functions {
        ////
        //// ---- Include files to compile the necessary custom (user-defined) Stan functions:
        ////
        #include "Stan_fns_basic.stan"
        #include "Stan_fns_Box_Cox.stan"
        #include "Stan_fns_corr.stan"
        #include "Stan_fns_ordinal.stan"
        #include "Stan_fns_log_lik.stan"
        #include "Stan_fns_simplex.stan"
        #include "Stan_fns_Jacobian.stan"
        #include "Stan_fns_model_fit.stan"
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
        // array[2] matrix[n_studies, n_thr] n;
        array[2] matrix[n_studies, n_thr + 1] x;
        array[2] matrix[n_studies, n_thr] cutpoint_index;
        ////
        vector[n_thr] cts_thr_values;
        ////
        //// Priors for locations ("beta"):
        ////
        vector[2] prior_beta_mu_mean;
        vector[2] prior_beta_mu_SD;
        vector[2] prior_beta_SD_mean;
        vector[2] prior_beta_SD_SD;
        ////
        //// Priors for raw scales:
        ////
        vector[2] prior_raw_scale_mu_mean;
        vector[2] prior_raw_scale_mu_SD;
        vector[2] prior_raw_scale_SD_mean;
        vector[2] prior_raw_scale_SD_SD;
        ////
        //// Priors for box-cox ("lambda"):
        ////
        real prior_boxcox_lambda_mean;
        real prior_boxcox_lambda_SD;
        ////
        //// Priors (and possible restrictions) for between-study correlations:
        ////
        real<lower=-1.0, upper=+1.0> beta_corr_lb;
        real<lower=beta_corr_lb, upper=+1.0>  beta_corr_ub;
        real<lower=0.0>  prior_beta_corr_LKJ;
        ////
        real<lower=-1.0, upper=+1.0> raw_scale_corr_lb;
        real<lower=raw_scale_corr_lb, upper=+1.0>  raw_scale_corr_ub;
        real<lower=0.0>  prior_raw_scale_corr_LKJ;
        ////
        //// Other:
        ////
        int<lower=0, upper=1> box_cox;
        int<lower=0, upper=1> softplus;
  
}


parameters {
    
        row_vector[2] beta_mu;    
        row_vector<lower=0.0>[2] beta_SD;   
        matrix[n_studies, 2] beta_z; //// Study-specific random effects (off-centered parameterisation)
        ////
        row_vector[2] raw_scale_mu;    
        row_vector<lower=0.0>[2] raw_scale_SD;
        matrix[n_studies, 2] raw_scale_z; //// Study-specific random effects (off-centered parameterisation)
        ////
        real<lower=beta_corr_lb, upper=beta_corr_ub> beta_corr;  //// between-study corr (possibly restricted)
        real<lower=raw_scale_corr_lb, upper=raw_scale_corr_ub> raw_scale_corr;  //// between-study corr (possibly restricted)
        ////
        real<lower=-5.0, upper=5.0> lambda; //// box-cox params
   
}


transformed parameters {
  
        ////
        //// ---- Construct simple 2x2 (bivariate) between-study corr matrices:
        ////
        cholesky_factor_corr[2] beta_L_Omega      = make_bivariate_L_Omega(beta_corr);
        cholesky_factor_corr[2] raw_scale_L_Omega = make_bivariate_L_Omega(raw_scale_corr);
        cholesky_factor_cov[2]  beta_L_Sigma       = diag_pre_multiply(beta_SD, beta_L_Omega);
        cholesky_factor_cov[2]  raw_scale_L_Sigma  = diag_pre_multiply(raw_scale_SD, raw_scale_L_Omega);
        ////
        //// ---- Compute (global) cutpoints using either log or box-cox transform:
        //// 
        vector[n_thr] C = ((box_cox == 0) ? log(cts_thr_values) : fn_Stan_box_cox(cts_thr_values, lambda));  
        ////
        //// ---- Likelihood stuff:
        ////
        array[2] matrix[n_studies, n_thr] surv_prob  = init_array_of_matrices(n_studies, n_thr, 2, 1.0); 
        array[2] matrix[n_studies, n_thr] cond_prob  = init_array_of_matrices(n_studies, n_thr, 2, 0.0);
        array[2] matrix[n_studies, n_thr] log_lik    = init_array_of_matrices(n_studies, n_thr, 2, 0.0); 
        matrix[n_studies, 2] raw_scale; // NOT local - needed for log_get_J later 
        {
                ////
                //// ---- Between-study model for location and scale:
                ////
                matrix[n_studies, 2] beta;  // local
                for (s in 1:n_studies) {
                     beta[s, ]      = beta_mu      + beta_z[s, ] * beta_L_Sigma;
                     raw_scale[s, ] = raw_scale_mu + raw_scale_z[s, ] * raw_scale_L_Sigma;
                }
                //// 
                //// ---- Compute scales and Jacobian adjustment for raw_scale -> scale transformation (w/ Jacobian for raw_scale -> scale)
                ////
                matrix[n_studies, 2] scale = ((softplus == 1) ? softplus_scaled(raw_scale) : exp(raw_scale));  // local
                ////
                //// ---- Get the cutpoint index (k) to map "latent_surv[c][s, cut_i]" to correct cutpoint "C[k]":
                ////
                array[2] matrix[n_studies, n_thr] latent_surv = map_latent_surv_prob_to_fixed_C(C, beta, scale, n_studies, n_obs_cutpoints, cutpoint_index);  // local
                ////
                //// ---- Multinomial (factorised binomial likelihood)
                ////
                int use_probit_link = 1;
                array[2, 3] matrix[n_studies, n_thr] log_lik_outs = compute_log_lik_binomial_probit_fact(latent_surv, use_probit_link, x, n_obs_cutpoints);
                log_lik    = log_lik_outs[, 1];
                cond_prob  = log_lik_outs[, 2];
                surv_prob  = log_lik_outs[, 3];
        }
      
}


model {
       
        ////
        //// ---- Priors for locations:
        ////
        beta_mu    ~ normal(prior_beta_mu_mean, prior_beta_mu_SD);
        beta_SD    ~ normal(prior_beta_SD_mean, prior_beta_SD_SD); 
        beta_L_Omega ~ lkj_corr_cholesky(prior_beta_corr_LKJ);
        ////
        //// ---- Priors for scales:
        ////
        raw_scale_mu    ~ normal(prior_raw_scale_mu_mean, prior_raw_scale_mu_SD);
        raw_scale_SD    ~ normal(prior_raw_scale_SD_mean, prior_raw_scale_SD_SD);
        raw_scale_L_Omega ~ lkj_corr_cholesky(prior_raw_scale_corr_LKJ);
        ////
        //// ---- Prior for box-cox parameters:
        ////
        lambda ~ normal(prior_boxcox_lambda_mean, prior_boxcox_lambda_SD);
        ////
        //// ---- Likelihood / Model:
        ////
        target += std_normal_lpdf(to_vector(beta_z)); // part of between-study model, NOT prior
        target += std_normal_lpdf(to_vector(raw_scale_z)); // part of between-study model, NOT prior
        target += raw_scale_to_scale_log_det_J_lp(raw_scale, softplus);
        ////
        //// ---- Increment the log-likelihood:
        //// 
        for (c in 1:2) {
          target +=  sum(log_lik[c]);
        } 
  
}





generated quantities { 
          ////
          //// ---- Compute between-study variance-covariance matrix for location parameters:
          ////
          corr_matrix[2] beta_Omega      = multiply_lower_tri_self_transpose(beta_L_Omega);
          corr_matrix[2] raw_scale_Omega = multiply_lower_tri_self_transpose(raw_scale_L_Omega);
          ////
          //// ---- Calculate summary accuracy (using mean parameters):
          ////
          row_vector[2] scale_mu = (softplus == 1) ? softplus_scaled(raw_scale_mu) : exp(raw_scale_mu);
          ////
          vector[n_thr] Fp = Phi(-(C - beta_mu[1])/scale_mu[1]);
          vector[n_thr] Sp = 1.0 - Fp;
          vector[n_thr] Se = Phi(-(C - beta_mu[2])/scale_mu[2]);
          ////
          //// ---- Calculate predictive accuracy:
          ////
          matrix[2, 2] beta_Sigma      = multiply_lower_tri_self_transpose(beta_L_Sigma);
          matrix[2, 2] raw_scale_Sigma = multiply_lower_tri_self_transpose(raw_scale_L_Sigma);
          ////
          row_vector[2] beta_pred      = to_row_vector(multi_normal_cholesky_rng(beta_mu,      raw_scale_L_Sigma));
          row_vector[2] raw_scale_pred = to_row_vector(multi_normal_cholesky_rng(raw_scale_mu, beta_L_Sigma));
          row_vector[2] scale_pred = (softplus == 1) ? softplus_scaled(raw_scale_pred) : exp(raw_scale_pred);
          ////
          vector[n_thr] Fp_pred = Phi(-(C - beta_pred[1])/scale_pred[1]);
          vector[n_thr] Sp_pred = 1.0 - Fp_pred;
          vector[n_thr] Se_pred = Phi(-(C - beta_pred[2])/scale_pred[2]);
          ////
          //// ---- Calculate study-specific accuracy:
          ////
          matrix[n_studies, n_thr] fp = surv_prob[1];
          matrix[n_studies, n_thr] sp = 1.0 - fp;
          matrix[n_studies, n_thr] se = surv_prob[2];
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
                          x_hat[c][s, cut_i] = cond_prob[c][s, cut_i] * x[c][s, cut_i];  	 // Fitted values
                        
                          //// Compute residual deviance contribution:
                          real n_i =  (x[c][s, cut_i]);
                          real x_i =  (x[c][s, cut_i + 1]);
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























