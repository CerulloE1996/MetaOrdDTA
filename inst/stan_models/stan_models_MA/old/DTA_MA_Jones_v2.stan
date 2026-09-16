


functions {
        ////
        //// ---- Include files to compile the necessary custom (user-defined) Stan functions:
        ////
        #include "Stan_fns_basic.stan"
        #include "Stan_fns_Box_Cox.stan"
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
        array[n_studies] int n_obs_cutpoints;
        ////
        array[2, n_studies, n_thr + 1] int x;
        array[2, n_studies, n_thr + 1] int cutpoint_index;
        ////
        vector[n_thr] cts_thr_values; 
        ////
        //// ---- Covariates:
        ////
        int n_covariates_nd;
        int n_covariates_d;
        int n_covariates_max;
        matrix[n_studies, n_covariates_nd] X_nd;
        matrix[n_studies, n_covariates_d]  X_d;
        vector[n_covariates_nd] baseline_case_nd;  // must be user-inputted - e.g. could be {0, 1, 45.3} for 2 binary covariates and 1 cts one (e.g. age)
        vector[n_covariates_d]  baseline_case_d;   // must be user-inputted - e.g. could be {0, 1, 45.3} for 2 binary covariates and 1 cts one (e.g. age)
        ////
        //// ---- Priors for locations ("beta"):
        ////
        matrix[2, n_covariates_max] prior_beta_mu_mean;
        matrix[2, n_covariates_max] prior_beta_mu_SD;  
        vector[2] prior_beta_SD_mean;   
        vector[2] prior_beta_SD_SD;   
        ////
        //// ---- Priors for raw scales: 
        ////
        matrix[2, n_covariates_max] prior_raw_scale_mu_mean;
        matrix[2, n_covariates_max] prior_raw_scale_mu_SD;  
        vector[2] prior_raw_scale_SD_mean; 
        vector[2] prior_raw_scale_SD_SD; 
        //// 
        //// ---- Priors for box-cox ("lambda"):
        ////
        real prior_boxcox_lambda_mean;
        real prior_boxcox_lambda_SD;
        ////
        //// ---- Priors (and possible restrictions) for between-study correlations:
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
        ////
        int use_probit_link;
  
}

transformed data { 
        ////
        array[2, n_studies, n_thr] int x_2;
        array[2, n_studies] int N_total;
        array[2, n_studies, n_thr] int n;
        ////
        for (s in 1:n_studies) {
              for (c in 1:2) {
                  N_total[c, s] = x[c, s, 1];
                  for (k in 1:n_thr) { 
                     x_2[c, s, k] = x[c, s, k + 1];
                  }
              }
              for (c in 1:2) {
                   n[c, s, 1] = N_total[c, s];
                   for (k in 2:n_obs_cutpoints[s]) {
                              n[c, s, k] = x_2[c, s, k - 1];
                   }
              }
       }
}


parameters {
        matrix[2, n_covariates_max] beta_mu;  
        row_vector<lower=0.0>[2] beta_SD;   
        matrix[n_studies, 2] beta_z; //// Study-specific random effects (off-centered parameterisation)
        ////
        matrix[2, n_covariates_max] raw_scale_mu;
        row_vector<lower=0.0>[2] raw_scale_SD;
        matrix[n_studies, 2] raw_scale_z; //// Study-specific random effects (off-centered parameterisation)
        ////
        real<lower=beta_corr_lb, upper=beta_corr_ub> beta_corr;  //// between-study corr (possibly restricted)
        real<lower=raw_scale_corr_lb, upper=raw_scale_corr_ub> raw_scale_corr;  //// between-study corr (possibly restricted)
        ////
        vector<lower=-5.0, upper=5.0>[2] lambda; //// box-cox params
}


transformed parameters {
        ////
        //// ---- Construct simple 2x2 (bivariate) between-study corr matrices:
        ////
        cholesky_factor_corr[2] beta_L_Omega       = make_bivariate_L_Omega(beta_corr);
        cholesky_factor_corr[2] raw_scale_L_Omega  = make_bivariate_L_Omega(raw_scale_corr);
        cholesky_factor_cov[2]  beta_L_Sigma       = diag_pre_multiply(beta_SD, beta_L_Omega);
        cholesky_factor_cov[2]  raw_scale_L_Sigma  = diag_pre_multiply(raw_scale_SD, raw_scale_L_Omega);
        ////
        //// ---- Compute (global) cutpoints using either log or box-cox transform:
        //// 
        array[2] vector[n_thr] C;
        C[1] = ((box_cox == 0) ? log(cts_thr_values) : fn_Stan_box_cox(cts_thr_values, lambda[1]));  
        C[2] = ((box_cox == 0) ? log(cts_thr_values) : fn_Stan_box_cox(cts_thr_values, lambda[2]));
        ////
        //// ---- Study-level random effects (after Cholesky decomposition):
        ////
        matrix[n_studies, 2] beta_random;
        matrix[n_studies, 2] raw_scale_random;
        for (s in 1:n_studies) {
             beta_random[s, ]      =  beta_z[s, ]      * beta_L_Sigma;
             raw_scale_random[s, ] =  raw_scale_z[s, ] * raw_scale_L_Sigma;
        }
        ////
        //// ---- Linear predictors for each disease statu + apply covariates for non-diseased and diseased groups:
        ////
        matrix[n_studies, 2] Xbeta; // NOT local
        matrix[n_studies, 2] Xraw_scale; // NOT local
        matrix[n_studies, 2] Xscale; // NOT local
        ////
        Xbeta[, 1]       = X_nd[1:n_studies, 1:n_covariates_nd]  * to_vector(beta_mu[1, 1:n_covariates_nd]) + beta_random[, 1];
        Xbeta[, 2]       = X_d[1:n_studies,  1:n_covariates_d]   * to_vector(beta_mu[2, 1:n_covariates_d])  + beta_random[, 2];
        ////
        Xraw_scale[, 1]  = X_nd[1:n_studies, 1:n_covariates_nd]  * to_vector(raw_scale_mu[1, 1:n_covariates_nd]) + raw_scale_random[, 1];
        Xraw_scale[, 2]  = X_d[1:n_studies,  1:n_covariates_d]   * to_vector(raw_scale_mu[2, 1:n_covariates_d])  + raw_scale_random[, 2];
        ////
        Xscale = ((softplus == 1) ? softplus_scaled(Xraw_scale) : exp(Xraw_scale));  // local
}


model {
        ////
        //// ---- Priors for locations:
        ////
        beta_mu[1, 1:n_covariates_nd] ~ normal(prior_beta_mu_mean[1, 1:n_covariates_nd], prior_beta_mu_SD[1, 1:n_covariates_nd]);  
        beta_mu[2, 1:n_covariates_d]  ~ normal(prior_beta_mu_mean[2, 1:n_covariates_d],  prior_beta_mu_SD[2, 1:n_covariates_d]);
        ////
        beta_SD    ~ normal(prior_beta_SD_mean, prior_beta_SD_SD); 
        beta_L_Omega ~ lkj_corr_cholesky(prior_beta_corr_LKJ);
        ////
        //// ---- Priors for scales:
        ////
        raw_scale_mu[1, 1:n_covariates_nd] ~ normal(prior_raw_scale_mu_mean[1, 1:n_covariates_nd], prior_raw_scale_mu_SD[1, 1:n_covariates_nd]);  
        raw_scale_mu[2, 1:n_covariates_d]  ~ normal(prior_raw_scale_mu_mean[2, 1:n_covariates_d],  prior_raw_scale_mu_SD[2, 1:n_covariates_d]);
        ////
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
        ////
        target += raw_scale_to_scale_log_det_J_lp(Xraw_scale, softplus);
        ////
        //// ---- Log-likelihood:
        ////
        {
            ////
            //// ---- Get the cutpoint index (k) to map "latent_surv[c][s, cut_i]" to correct cutpoint "C[k]":
            ////
            array[2] matrix[n_studies, n_thr] latent_surv = map_latent_surv_prob_to_fixed_hetero_C(
                                                            C, Xbeta, Xscale, n_studies, 
                                                            n_obs_cutpoints, cutpoint_index);
            target += compute_log_lik_binomial_fact_lp(
                      latent_surv, use_probit_link, n_thr, n_studies, 
                      x_2, n, N_total, n_obs_cutpoints);
        }
}


generated quantities { 
          ////
          //// ---- Compute between-study variance-covariance matrix for location parameters:
          ////
          corr_matrix[2] beta_Omega      = multiply_lower_tri_self_transpose(beta_L_Omega);
          corr_matrix[2] raw_scale_Omega = multiply_lower_tri_self_transpose(raw_scale_L_Omega);
          ////
          matrix[2, 2] beta_Sigma      = multiply_lower_tri_self_transpose(beta_L_Sigma);
          matrix[2, 2] raw_scale_Sigma = multiply_lower_tri_self_transpose(raw_scale_L_Sigma);
          ////
          //// ---- Calculate summary accuracy - "baseline" covariate values:
          ////  
          real Xbeta_baseline_nd = dot_product(baseline_case_nd, to_vector(beta_mu[1, 1:n_covariates_nd]));
          real Xbeta_baseline_d  = dot_product(baseline_case_d,  to_vector(beta_mu[2, 1:n_covariates_d]));
          ////
          real Xraw_scale_baseline_nd = dot_product(baseline_case_nd, raw_scale_mu[1, 1:n_covariates_nd]);
          real Xraw_scale_baseline_d  = dot_product(baseline_case_d,  raw_scale_mu[2, 1:n_covariates_d]);
          real scale_nd_baseline = ((softplus == 1) ? 
                                   softplus_scaled(Xraw_scale_baseline_nd) : 
                                   exp(Xraw_scale_baseline_nd - 0.5*square(beta_SD[1])));
          real scale_d_baseline  = ((softplus == 1) ? 
                                   softplus_scaled(Xraw_scale_baseline_d)  : 
                                   exp(Xraw_scale_baseline_d  - 0.5*square(beta_SD[2])));
          ////
          //// ---- Calculate baseline Se/Sp:
          ////
          vector[n_thr] Fp_baseline = (use_probit_link == 1) ? 
                                      Phi(-(C[1] - Xbeta_baseline_nd)/scale_nd_baseline) : 
                                      inv_logit(-(C[1] - Xbeta_baseline_nd)/scale_nd_baseline);
          vector[n_thr] Sp_baseline = 1.0 - Fp_baseline;
          vector[n_thr] Se_baseline = (use_probit_link == 1) ?
                                      Phi(-(C[2] - Xbeta_baseline_d)/scale_d_baseline) : 
                                      inv_logit(-(C[2] - Xbeta_baseline_d)/scale_d_baseline);
          ////
          //// ---- Calculate predictive accuracy:
          ////
          vector[2] beta_random_pred      = multi_normal_cholesky_rng(rep_vector(0.0, 2), beta_L_Sigma);
          vector[2] raw_scale_random_pred = multi_normal_cholesky_rng(rep_vector(0.0, 2), raw_scale_L_Sigma);
          ////
          real Xbeta_baseline_pred_nd = dot_product(baseline_case_nd, beta_mu[1, 1:n_covariates_nd]) + beta_random_pred[1];
          real Xbeta_baseline_pred_d  = dot_product(baseline_case_d,  beta_mu[2, 1:n_covariates_d])  + beta_random_pred[2];
          ////
          real Xraw_scale_baseline_pred_nd = dot_product(baseline_case_nd, raw_scale_mu[1, 1:n_covariates_nd]) + raw_scale_random_pred[1];
          real Xraw_scale_baseline_pred_d  = dot_product(baseline_case_d,  raw_scale_mu[2, 1:n_covariates_d]) + raw_scale_random_pred[2];
          real scale_nd_baseline_pred = ((softplus == 1) ? 
                                        softplus_scaled(Xraw_scale_baseline_pred_nd) : exp(Xraw_scale_baseline_pred_nd));
          real scale_d_baseline_pred  = ((softplus == 1) ? 
                                        softplus_scaled(Xraw_scale_baseline_pred_d)  : exp(Xraw_scale_baseline_pred_d));
          ////
          vector[n_thr] Fp_baseline_pred = (use_probit_link == 1) ? 
                                           Phi(-(C[1] - Xbeta_baseline_pred_nd)/scale_nd_baseline_pred)  : 
                                           inv_logit(-(C[1] - Xbeta_baseline_pred_nd)/scale_nd_baseline_pred);
          vector[n_thr] Sp_baseline_pred = 1.0 - Fp_baseline_pred;
          vector[n_thr] Se_baseline_pred = (use_probit_link == 1) ? 
                                           Phi(-(C[2] - Xbeta_baseline_pred_d)/scale_d_baseline_pred) : 
                                           inv_logit(-(C[2] - Xbeta_baseline_pred_d)/scale_d_baseline_pred);
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
              array[2] matrix[n_studies, n_thr] latent_surv = map_latent_surv_prob_to_fixed_hetero_C(
                                                              C, Xbeta, Xscale, n_studies, 
                                                              n_obs_cutpoints, cutpoint_index);
              ////
              array[3, 2] matrix[n_studies, n_thr] outs = compute_log_lik_binomial_fact_data(
                                                          latent_surv, use_probit_link, n_thr, n_studies,
                                                          x_2, n, N_total, n_obs_cutpoints);
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
              array[4] matrix[n_studies, n_thr] outs_model_fit = compute_deviance(
                                                                 cond_prob, n_thr, n_studies, 
                                                                 x_2, n, n_obs_cutpoints);
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















