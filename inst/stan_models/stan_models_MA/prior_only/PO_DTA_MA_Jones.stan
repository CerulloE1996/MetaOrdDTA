


functions {
     ////
     //// Include files to compile the necessary custom (user-defined) Stan functions:
     ////
     #include "Stan_fns_basic.stan"
     #include "Stan_fns_Box_Cox.stan"
     #include "Stan_fns_Pinkney_corr.stan"
}


data {
        ////
        //// Data:
        ////
        int<lower=1> n_studies;
        int<lower=1> n_thr;
        array[n_studies] int n_obs_cutpoints;
        ////
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
        int use_box_cox;
        int use_softplus_for_scales;
  
}


parameters {
    
        vector[2] beta_mu;    
        vector<lower=0.0>[2] beta_SD;   
        matrix[2, n_studies] beta_z; //// Study-specific random effects (off-centered parameterisation)
        ////
        vector[2] raw_scale_mu;    
        vector<lower=0.0>[2] raw_scale_SD;
        matrix[2, n_studies] raw_scale_z; //// Study-specific random effects (off-centered parameterisation)
        ////
        real beta_corr;  //// between-study corr (possibly restricted)
        real raw_scale_corr;  //// between-study corr (possibly restricted)
        ////
        real<lower=-5.0, upper=5.0> lambda; //// box-cox params
  
}


transformed parameters { 
  
        ////
        //// Construct the study-specific random effects (off-centered param.):
        ////
        matrix[2, n_studies] beta;
        matrix[2, n_studies] raw_scale;
        matrix[2, n_studies] scale;
        ////
        //// Construct simple 2x2 (bivariate) between-study corr matrices:
        ////
        cholesky_factor_corr[2] beta_L_Omega      = make_restricted_bivariate_L_Omega_jacobian(beta_corr, beta_corr_lb, beta_corr_ub);
        cholesky_factor_corr[2] raw_scale_L_Omega = make_restricted_bivariate_L_Omega_jacobian(raw_scale_corr, raw_scale_corr_lb, raw_scale_corr_ub);
        cholesky_factor_cov[2] beta_L_Sigma      = diag_pre_multiply(beta_SD, beta_L_Omega);
        cholesky_factor_cov[2] raw_scale_L_Sigma = diag_pre_multiply(raw_scale_SD, raw_scale_L_Omega);
        corr_matrix[2] beta_Omega      = multiply_lower_tri_self_transpose(beta_L_Omega);
        corr_matrix[2] raw_scale_Omega = multiply_lower_tri_self_transpose(raw_scale_L_Omega);
        ////
        //// Compute cutpoints using either log or box-cox transform:
        ////
        vector[n_thr] C = ((use_box_cox == 0) ? log(cts_thr_values) : box_cox(cts_thr_values, lambda)); 
        ////
        //// Between-study model for the location parameters ("beta") and the raw_scale parameters - models between-study correlation:
        ////
        for (s in 1:n_studies) {
             beta[, s]      =  beta_mu      + diag_pre_multiply(beta_SD, beta_L_Omega) * beta_z[, s];
             raw_scale[, s] =  raw_scale_mu + diag_pre_multiply(raw_scale_SD, raw_scale_L_Omega) * raw_scale_z[, s];
        }
        //// 
        //// Compute scales and Jacobian adjustment for raw_scale -> scale transformation (w/ Jacobian for raw_scale -> scale)
        ////
        scale  =  ((use_softplus_for_scales == 1) ? softplus_scaled_jacobian(raw_scale) : exp_jacobian(raw_scale));
        // //// Also need the Jacobian to go from (raw_scale_z, raw_scale_mu, raw_scale_SD) -> raw_scale! (using chain rule):
        // if (abs(sum(raw_scale_SD)) != 0.0) jacobian += log(abs(sum(raw_scale_SD)));      // double-checked the log-derivative of this by hand (correct)
        // if (abs(sum(raw_scale_z)) != 0.0)  jacobian += log(abs(sum(raw_scale_z)));       // double-checked the log-derivative of this by hand (correct)
        ////
        //// Likelihood using binomial factorization:
        ////
        for (s in 1:n_studies) {
                  for (c in 1:2) {
                      for (cut_i in 1:to_int(n_obs_cutpoints[s])) {
                             int k = to_int(cutpoint_index[c][s, cut_i]);
                             latent_cumul_prob[c][s, cut_i] = (C[k] - beta[c, s])/scale[c, s];
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
       
        ////
        //// Priors for locations:
        ////
        beta_mu    ~ normal(prior_beta_mu_mean, prior_beta_mu_SD);
        beta_SD    ~ normal(prior_beta_SD_mean, prior_beta_SD_SD);
        beta_Omega ~ lkj_corr(prior_beta_corr_LKJ); //// possibly truncated
        ////
        //// Priors for scales:
        ////
        raw_scale_mu    ~ normal(prior_raw_scale_mu_mean, prior_raw_scale_mu_SD);
        raw_scale_SD    ~ normal(prior_raw_scale_SD_mean, prior_raw_scale_SD_SD);
        raw_scale_Omega ~ lkj_corr(prior_raw_scale_corr_LKJ); //// possibly truncated
        ////
        //// Prior for box-cox parameters:
        ////
        lambda ~ normal(prior_boxcox_lambda_mean, prior_boxcox_lambda_SD);
        ////
        //// Likelihood / Model:
        ////
        to_vector(beta_z)      ~ std_normal();
        to_vector(raw_scale_z) ~ std_normal();
  
}





generated quantities {
  
          ////
          //// Calculate summary accuracy (using mean parameters):
          ////
          vector[2] scale_mu =  (use_softplus_for_scales == 1) ? softplus_scaled(raw_scale_mu) : exp(raw_scale_mu);
          ////
          vector[n_thr] Fp = 1.0 - Phi((C - beta_mu[1])/scale_mu[1]);
          vector[n_thr] Sp = 1.0 - Fp;
          vector[n_thr] Se = 1.0 - Phi((C - beta_mu[2])/scale_mu[2]);
          ////
          //// Calculate predictive accuracy:
          ////
          matrix[2, 2] beta_Sigma      = multiply_lower_tri_self_transpose(beta_L_Sigma);
          matrix[2, 2] raw_scale_Sigma = multiply_lower_tri_self_transpose(raw_scale_L_Sigma);
          ////
          vector[2] beta_pred      = multi_normal_cholesky_rng(beta_mu,      raw_scale_L_Sigma);
          vector[2] raw_scale_pred = multi_normal_cholesky_rng(raw_scale_mu, beta_L_Sigma);
          vector[2] scale_pred = (use_softplus_for_scales == 1) ? softplus_scaled(raw_scale_pred) : exp(raw_scale_pred);
          ////
          vector[n_thr] Fp_pred = 1.0 - Phi((C - beta_pred[1])/scale_pred[1]);
          vector[n_thr] Sp_pred = 1.0 - Fp_pred;
          vector[n_thr] Se_pred = 1.0 - Phi((C - beta_pred[2])/scale_pred[2]);
          ////
          //// Calculate study-specific accuracy:
          ////
          matrix[n_studies, n_thr] fp = 1.0 - cumul_prob[1];
          matrix[n_studies, n_thr] sp = 1.0 - fp;
          matrix[n_studies, n_thr] se = 1.0 - cumul_prob[2];
          ////
          matrix[n_studies, n_thr] x_hat_nd = rep_matrix(-1, n_studies, n_thr);
          matrix[n_studies, n_thr] x_hat_d  = rep_matrix(-1, n_studies, n_thr);
          matrix[n_studies, n_thr] dev_nd   = rep_matrix(-1, n_studies, n_thr);
          matrix[n_studies, n_thr] dev_d    = rep_matrix(-1, n_studies, n_thr);
          array[2] matrix[n_studies, n_thr] x_hat = init_array_of_matrices(n_studies, n_thr, 2, -1);
          array[2] matrix[n_studies, n_thr] dev   = init_array_of_matrices(n_studies, n_thr, 2, -1);
    
}























