


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
        array[2] matrix[n_studies, n_thr + 1] x;
        array[2] matrix[n_studies, n_thr + 1] cutpoint_index;
        ////
        //// ---- Covariates:
        ////
        int n_covariates_nd;
        int n_covariates_d;
        int n_covariates_max;
        matrix[n_studies, n_covariates_nd] X_nd; // must be user-inoutted - study-level covariates for D-
        matrix[n_studies, n_covariates_d]  X_d;  // must be user-inoutted - study-level covariates for D+
        vector[n_covariates_nd] baseline_case_nd; // must be user-inputted - e.g. could be {0, 1, 45.3} for 2 binary covariates and 1 cts one (e.g. age)
        vector[n_covariates_d] baseline_case_d;   // must be user-inputted - e.g. could be {0, 1, 45.3} for 2 binary covariates and 1 cts one (e.g. age)
        ////
        //// ---- Priors for locations ("beta"):
        ////
        matrix[2, n_covariates_max] prior_beta_mu_mean;
        matrix[2, n_covariates_max] prior_beta_mu_SD;
        vector[2] prior_beta_SD_mean;
        vector[2] prior_beta_SD_SD;
        ////
        //// ---- Priors for between-study correlation matrix: 
        ////
        real<lower=-1.0> beta_corr_lb;
        real<lower=beta_corr_lb, upper=1.0>  beta_corr_ub;
        real<lower=0.0>  prior_beta_corr_LKJ;
        ////
        //// ---- Induced-Dirichlet priors:
        ////
        array[2] vector[n_thr + 1] prior_dirichlet_phi;
        real<lower=0.0> prior_dirichlet_prob_SDs;
        ////
        //// ---- Other:
        ////
        real<lower=0.0> kappa_lb;
        int<lower=0, upper=1> softplus;
}


transformed data {  
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
       int prior_only = 0;
       ////
       array[n_studies] int dummy_ind_test_in_study;
       for (s in 1:n_studies) {
          dummy_ind_test_in_study[s] = 1;
       }
}
  

parameters {
        matrix[2, n_covariates_max] beta_mu;      
        row_vector<lower=0.0>[2] beta_SD;   
        matrix[n_studies, 2] beta_z;    //// Study-specific random effects (off-centered parameterisation)
        ////
        real<lower=beta_corr_lb, upper=beta_corr_ub> beta_corr;  //// between-study corr (possibly restricted)
        ////
        array[2] matrix[n_studies, n_thr] C_raw;
        ////w;
        // vector<lower=(kappa_lb)>[2] kappa;
        array[2] simplex[n_cat] dirichlet_phi; // props/probs
        ////
        array[2] vector<lower=0>[n_cat] target_sds;
}

 
transformed parameters { 
        array[2] vector<lower=0.0>[n_cat] alpha; 
        vector<lower=(kappa_lb)>[2] kappa;
        array[2] vector<lower=0.0>[n_cat] dir_cat_SDs_sigma;
        ////
        vector[2] avg_sd;
        for (c in 1:2) {
            ////
            //// ---- Use a simple approximation based on the average SD:
            ////
            avg_sd[c] = mean(target_sds[c]);
            kappa[c] = max_reals(kappa_lb, (1.0/square(avg_sd[c])) - 1.0 ); ////  For a symmetric Dirichlet, SD ≈ sqrt(1/(α_0+1)):
            alpha[c] = dirichlet_phi[c] * kappa[c]; //// Set alpha using the approximated kappa:
            ////
            //// ---- Calculate the actual SDs:
            ////
            real alpha_0 = sum(alpha[c]);
            dir_cat_SDs_sigma[c] = sqrt((alpha[c] .* (alpha_0 - alpha[c])) / (square(alpha_0) * (alpha_0 + 1.0)));
        }
        ////
        //// ---- Construct study-specific cutpoints:
        //// 
        array[2] matrix[n_studies, n_thr] C; ////= convert_C_array_to_mat(C_array);
        for (c in 1:2) {
          for (s in 1:n_studies) {
                     C[c][s, ] =  construct_C(C_raw[c][s, 1:n_thr], softplus);
          }
        }
        ////
        //// ---- Construct simple 2x2 (bivariate) between-study corr matrices:
        ////
        cholesky_factor_corr[2] beta_L_Omega = make_bivariate_L_Omega(beta_corr);
        cholesky_factor_cov[2]  beta_L_Sigma = diag_pre_multiply(beta_SD, beta_L_Omega);
        ////
        //// ---- Study-level random effects (after Cholesky decomposition):
        ////
        matrix[n_studies, 2] beta_random;
        for (s in 1:n_studies) {
             beta_random[s, ] =  beta_z[s, ] * beta_L_Sigma;
        }
        ////
        //// ---- Linear predictors for each disease statu + apply covariates for non-diseased and diseased groups:
        ////
        matrix[n_studies, 2] Xbeta;
        Xbeta[, 1]  = X_nd[1:n_studies, 1:n_covariates_nd] * to_vector(beta_mu[1, 1:n_covariates_nd]) + beta_random[, 1];
        Xbeta[, 2]  = X_d[1:n_studies, 1:n_covariates_d]   * to_vector(beta_mu[2, 1:n_covariates_d])  + beta_random[, 2];
        // ////
        // //// ---- Between-study model for the location parameters ("beta") - models between-study correlation:
        // ////
        // matrix[n_studies, 2] beta;
        // for (s in 1:n_studies) {
        //      beta[s, ] = beta_mu +  beta_z[s, ] * beta_L_Sigma; 
        // }
}


model {
        ////
        //// ---- Priors: 
        ////
        beta_mu[1, 1:n_covariates_nd] ~ normal(prior_beta_mu_mean[1, 1:n_covariates_nd], prior_beta_mu_SD[1, 1:n_covariates_nd]);  
        beta_mu[2, 1:n_covariates_d]  ~ normal(prior_beta_mu_mean[2, 1:n_covariates_d],  prior_beta_mu_SD[2, 1:n_covariates_d]);
        ////
        beta_SD ~ normal(prior_beta_SD_mean, prior_beta_SD_SD); 
        beta_L_Omega ~ lkj_corr_cholesky(prior_beta_corr_LKJ); 
        ////
        for (c in 1:2) {
              // kappa[c] ~ normal(prior_kappa_mean, prior_kappa_SD);  // ID prior
              dirichlet_phi[c] ~ dirichlet(prior_dirichlet_phi[c]); // ID prior
              // if (prior_only == 0) target +=  n_cat * log(kappa[c]);
              ////
              target_sds[c] ~ normal(0.0, prior_dirichlet_prob_SDs);
              target += -sum(square(dir_cat_SDs_sigma[c] - target_sds[c]));
              ////
              //// Jacobian adjustment for the transformation from target_sds to kappa:
              ////
              //// Avoid including the Jacobian adjustment if we hit the max_reals constraint
              if (kappa[c] > kappa_lb) {
                   target += log(2.0) - 3.0 * log(avg_sd[c]);
                   target += - log(n_cat);
              }
        }
        ////
        //// ---- Induced-dirichlet between study ** model ** (NOT a prior model here but part of the actual likelihood since random-effect cutpoints!):
        ////
        if (prior_only == 0) {
          for (c in 1:2) {
               for (s in 1:n_studies) {
                   vector[n_cat] Ind_Dir_ord_prob;
                   if (use_probit_link == 1) Ind_Dir_ord_prob = to_vector(cumul_probs_to_ord_probs(Phi(C[c][s, ])));
                   else                      Ind_Dir_ord_prob = to_vector(cumul_probs_to_ord_probs(inv_logit(C[c][s, ])));
                   ////
                   Ind_Dir_ord_prob ~ induced_dirichlet_given_C(to_vector(C[c][s, ]), alpha[c], use_probit_link);
                   // to_vector(C[c][s,]) ~ induced_dirichlet(alpha[c], 0.0);
               }
          }
        }
        ////
        //// ---- Likelihood / Model:
        ////
        if (prior_only == 0) {
          target += std_normal_lpdf(to_vector(beta_z)); // part of between-study model, NOT prior // beta_{c, s} ~ normal(beta_mu_{c}, beta_SD_{c})
          ////
          for (c in 1:2) {
             target += raw_C_to_C_log_det_J_lp(C_raw[c], softplus);
          }
        }
        ////
        //// ---- Log-likelihood:
        ////
        if (prior_only == 0) {
            ////
            //// ---- Get the cutpoint index (k) to map "latent_surv[c][s, cut_i]" to correct cutpoint "C[k]":
            ////
            real scale = 1.0; // since using "Xu-like"" param.
            array[2] matrix[n_studies, n_thr] latent_surv = map_latent_surv_to_random_hetero_C(C, n_thr, Xbeta, scale, n_studies, n_obs_cutpoints, cutpoint_index);
            target += compute_log_lik_binomial_fact_lp(latent_surv, use_probit_link, n_thr, x_2, n, N_total, n_obs_cutpoints, dummy_ind_test_in_study);
        }

}


generated quantities {
      ////
      //// ---- Compute between-study variance-covariance matrix for location parameters:
      ////
      cov_matrix[2]  beta_Sigma = multiply_lower_tri_self_transpose(beta_L_Sigma);
      corr_matrix[2] beta_Omega = multiply_lower_tri_self_transpose(beta_L_Omega);
      ////
      //// ---- SD's for each category probability in a Dirichlet is √(α_k(α_0-α_k)/(α_0^2(α_0+1))) - where α_0 is the sum of all alphas:
      ////
      array[2] vector<lower=0.0>[n_cat] dirichlet_cat_SDs_sigma;
      for (c in 1:2) {
          real alpha_0 = sum(alpha[c]);
          dirichlet_cat_SDs_sigma[c] = sqrt((alpha[c] .* (alpha_0 - alpha[c]) ) / (square(alpha_0) * (alpha_0 + 1.0)));
      }
      ////
      //// ---- Calculate summary cutpoints:
      ////
      array[2] vector[n_thr] C_mu;
      for (c in 1:2) {
            for (k in 1:n_thr) {
                C_mu[c][k] = median(C[c][1:n_studies, k]); // solves "skewness" issue from non-linear transformations
            }
      }
      ////
      //// ---- Calculate summary accuracy - "baseline" covariate values:
      ////  
      real Xbeta_baseline_nd = dot_product(baseline_case_nd, to_vector(beta_mu[1, 1:n_covariates_nd])); 
      real Xbeta_baseline_d  = dot_product(baseline_case_d,  to_vector(beta_mu[2, 1:n_covariates_d]));
      ////
      //// ---- Calculate baseline Se/Sp:
      ////
      vector[n_thr] Fp_baseline = (use_probit_link == 1) ? Phi(-(C_mu[1] - Xbeta_baseline_nd))  : inv_logit(-(C_mu[1] - Xbeta_baseline_nd));
      vector[n_thr] Sp_baseline = 1.0 - Fp_baseline;
      vector[n_thr] Se_baseline = (use_probit_link == 1) ? Phi(-(C_mu[2] - Xbeta_baseline_d))   : inv_logit(-(C_mu[2] - Xbeta_baseline_d));
      ////
      //// ---- Calculate predictive accuracy:
      ////
      array[2] vector[n_thr] C_pred;
      ////
      for (c in 1:2) { 
        for (k in 1:n_thr) {
            real sigma_C_k;
            real SD_prob_scale = dirichlet_cat_SDs_sigma[c][k];
            if (use_probit_link == 1) {
                real mean_C_cutpoint_scale = C_mu[c][k];
                sigma_C_k = SD_approx_ID_ord_prob_to_C_probit(mean_C_cutpoint_scale, SD_prob_scale);
            } else { 
                real mean_prob_scale = dirichlet_phi[c][k];
                sigma_C_k = SD_approx_ID_ord_prob_to_C_probit(mean_prob_scale, SD_prob_scale);
            }
            C_pred[c][k] = C_mu[c][k] + normal_rng(0, sigma_C_k);
        }
      } 
      ////
      vector[2] beta_random_pred = multi_normal_cholesky_rng(rep_vector(0.0, 2), beta_L_Sigma); 
      ////
      //// ---- Use baseline covariates with predicted random effects
      ////
      real Xbeta_baseline_pred_nd = dot_product(baseline_case_nd, beta_mu[1, 1:n_covariates_nd]) + beta_random_pred[1];
      real Xbeta_baseline_pred_d  = dot_product(baseline_case_d,  beta_mu[2, 1:n_covariates_d])  + beta_random_pred[2];
      ////
      vector[n_thr] Fp_baseline_pred = (use_probit_link == 1) ? Phi(-(C_pred[1] - Xbeta_baseline_pred_nd)) : inv_logit(-(C_pred[1] - Xbeta_baseline_pred_nd));
      vector[n_thr] Sp_baseline_pred = 1.0 - Fp_baseline_pred;
      vector[n_thr] Se_baseline_pred = (use_probit_link == 1) ? Phi(-(C_pred[2] - Xbeta_baseline_pred_d))  : inv_logit(-(C_pred[2] - Xbeta_baseline_pred_d));
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
          array[2] matrix[n_studies, n_thr] latent_surv = map_latent_surv_to_random_hetero_C(C, n_thr, Xbeta, scale, n_studies, n_obs_cutpoints, cutpoint_index);
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
                 deviance_nd[s] += is_nan(dev_nd[s, k]) ? 0.0 : dev_nd[s, k];
                 deviance_d[s]  += is_nan(dev_d[s, k]) ? 0.0 : dev_d[s, k];
              }
              deviance[s] = deviance_nd[s] + deviance_d[s];
          }
      }

}








