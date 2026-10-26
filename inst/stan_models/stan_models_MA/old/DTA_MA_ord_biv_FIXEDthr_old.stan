


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
        int n_covariates_max_nd;
        int n_covariates_max_d;
        int n_covariates_max;
        matrix[n_studies, n_covariates_max_nd] X_nd; //// covariate array (can have  DIFFERENT NUMBERS of covariates for each  outcome - fill rest of array with 999999 if they vary between outcomes)
        matrix[n_studies, n_covariates_max_d]  X_d;  //// covariate array (can have  DIFFERENT NUMBERS of covariates for each  outcome - fill rest of array with 999999 if they vary between outcomes)
        array[2] int n_covs_per_outcome;
        ////
        //// ---- Priors for locations:
        ////
        matrix[2, n_covariates_max] prior_beta_mu_mean;
        matrix[2, n_covariates_max] prior_beta_mu_SD;  
        vector[2] prior_beta_SD_mean;
        vector[2] prior_beta_SD_SD;
        ////
        //// ---- Priors for between-study correlation:
        ////
        real<lower=-1.0> beta_corr_lb;
        real<lower=beta_corr_lb, upper=1.0>  beta_corr_ub;
        real<lower=0.0>  prior_beta_corr_LKJ;
        //// 
        //// ---- Priors for cutpoints (using "induced-Dirichlet"):
        //// 
        array[2] vector<lower=0.0>[n_thr + 1] prior_dirichlet_alpha; 
        ////
        //// ---- Other:
        //// 
        int<lower=0, upper=1> softplus; 
}

 
transformed data {   
        int n_cat = n_thr + 1; //// Number of ordinal categories for index test 
        // vector[n_thr] Ind_Dir_anchor = rep_vector(0.0, n_thr);} 
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
          
        // int n_total_non_missing_obs = 0;
        // 
        // for (c in 1:2) {
        //   for (s in 1:n_studies) {
        //           for (i in 1:(n_obs_cutpoints[s])) {
        //              
        //                   //// ---- Current and PREVIOUS NON-MISSING counts //// next
        //                   int x_current  = to_int(x[c][s, i + 1]);
        //                   int x_previous = to_int(x[c][s, i + 0]); //// to_int(n[c][s, i]); // x_previous > x_next (counts are DECREASING)
        //                   ////
        //                   int missing_indicator = 0;
        //                   if (x_current == -1)  missing_indicator = 1;
        //                   if (x_previous == -1) missing_indicator = 1;
        //                   
        //                   if ((x_current != 0)  && (missing_indicator == 0)) { 
        //                       n_total_non_missing_obs += 1;
        //                   }
        //           }
        //   }
        // }
        // // vector[2] w = rep_vector(0.25, 2); 
}

 
parameters { 
        matrix[2, n_covariates_max] beta_mu;       
        row_vector<lower=0.0>[2] beta_SD;    
        matrix[n_studies, 2] beta_z;    //// Study-specific random-effects (raw params)
        //// 
        real<lower=beta_corr_lb, upper=beta_corr_ub> beta_corr;  //// between-study corr (possibly restricted)  
        //// 
        array[2] vector[n_thr] C_raw_vec;   //// Global cutpoints ("raw" / unconstrained) 
}  
 
 
transformed parameters { 
            ////   
            //// ---- Construct (global) cutpoints:
            ////
            array[2] vector[n_thr] C; 
            for (c in 1:2) { 
                C[c] = construct_C(C_raw_vec[c], softplus);  
            }
            ////
            //// ---- Construct simple 2x2 (bivariate) between-study corr matrices for between-study model:
            ////  
            cholesky_factor_corr[2] beta_L_Omega      = make_bivariate_L_Omega(beta_corr);
            cholesky_factor_cov[2]  beta_L_Sigma      = diag_pre_multiply(beta_SD, beta_L_Omega);
            ////
            // row_vector[2] w = beta_SD ./ (1.0 + beta_SD);
            ////
            //// Partially centered parameterization for beta:
            matrix[n_studies, 2] beta; 
            for (s in 1:n_studies) {
                // beta[s, ] = (rep_row_vector(1.0, 2) - w) .* beta_mu + beta_z[s, ] * beta_L_Sigma; // Matrix multiplication version (using Cholesky decomposition)
                 beta[s, ] = beta_mu + beta_z[s, ] * beta_L_Sigma; // Matrix multiplication version (using Cholesky decomposition)
            }
            // for (c in 1:2) {
            //    beta[, c] = (1.0 - w[c]) * beta_mu[c] + beta_z[, c] * beta_SD[c];
            // }
}
 
 
model { 
        ////
        //// ---- Priors: 
        ////
        beta_mu ~ normal(prior_beta_mu_mean, prior_beta_mu_SD);  
        beta_SD ~ normal(prior_beta_SD_mean, prior_beta_SD_SD); 
        beta_L_Omega ~ lkj_corr_cholesky(prior_beta_corr_LKJ); 
        //// 
        //// ---- Induced-dirichlet ** Prior ** model:
        //// 
        for (c in 1:2) {
             // vector[n_thr] Ind_Dir_cumul  = C[c];// - Ind_Dir_anchor;
             // vector[n_thr] Ind_Dir_cumul_prob = Phi(C[c]); 
              vector[n_cat] Ind_Dir_ord_prob; //// = cumul_probs_to_ord_probs(Phi(C));
             if (use_probit_link == 1) Ind_Dir_ord_prob = cumul_probs_to_ord_probs(Phi(C[c]));
             else                      Ind_Dir_ord_prob = cumul_probs_to_ord_probs(inv_logit(C[c]));
             ////
             Ind_Dir_ord_prob ~ induced_dirichlet_given_C(C[c], prior_dirichlet_alpha[c], use_probit_link); // more efficient than Betancourt et al. and seems to work fine. 
        }
        //// 
        //// ---- Likelihood / Model:
        ////
        for (c in 1:2) {
               // beta_z[, c] ~ normal(w[c] * beta_mu[c], 1.0);
                beta_z[, c] ~ normal(0.0, 1.0);
                target += raw_C_to_C_log_det_J_lp(C_raw_vec[c], softplus);
        }
        ////
        //// ---- Log-likelihood:
        ////
        {
            ////
            //// ---- Get the cutpoint index (k) to map "latent_surv[c][s, cut_i]" to correct cutpoint "C[k]":
            ////
            real scale = 1.0; // since using "Xu-like"" param.
            array[2] matrix[n_studies, n_thr] latent_surv = map_latent_surv_prob_to_fixed_hetero_C(C, beta, scale, n_studies, n_obs_cutpoints, cutpoint_index);
            target += compute_log_lik_binomial_fact_lp(latent_surv, use_probit_link, x_2, n, N_total, n_obs_cutpoints);
        }
}


generated quantities {
          // vector[2] pcp_weights = to_vector(w); // smaller w = more NON-centered (as w -> 1 we approach centered param.)
          ////
          //// ---- Compute between-study variance-covariance matrix for location parameters: 
          //// 
          cov_matrix[2]  beta_Sigma = multiply_lower_tri_self_transpose(beta_L_Sigma);
          corr_matrix[2] beta_Omega = multiply_lower_tri_self_transpose(beta_L_Omega);
          ////
          //// ---- Calculate summary accuracy (using mean parameters):
          ////  
          vector[n_thr] Fp = (use_probit_link == 1) ? Phi(-(C[1] - beta_mu[1])) : inv_logit(-(C[1] - beta_mu[1]));
          vector[n_thr] Sp = 1.0 - Fp;
          vector[n_thr] Se = (use_probit_link == 1) ? Phi(-(C[2] - beta_mu[2])) : inv_logit(-(C[2] - beta_mu[2]));
          //// 
          //// ---- Calculate predictive accuracy:
          ////
          vector[2] beta_pred =  multi_normal_cholesky_rng(beta_mu, beta_L_Sigma); 
          //// 
          vector[n_thr] Fp_pred = (use_probit_link == 1) ? Phi(-(C[1] - beta_pred[1])) : inv_logit(-(C[1] - beta_pred[1]));
          vector[n_thr] Sp_pred = 1.0 - Fp_pred;
          vector[n_thr] Se_pred = (use_probit_link == 1) ? Phi(-(C[2] - beta_pred[2])) : inv_logit(-(C[2] - beta_pred[2]));
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
              array[2] matrix[n_studies, n_thr] latent_surv = map_latent_surv_prob_to_fixed_hetero_C(C, beta, scale, n_studies, n_obs_cutpoints, cutpoint_index);
              ////
              array[3, 2] matrix[n_studies, n_thr] outs = compute_log_lik_binomial_fact_data(latent_surv, use_probit_link, x_2, n, N_total, n_obs_cutpoints);
              ////
              log_lik   = outs[1];
              array[2] matrix[n_studies, n_thr] cond_prob = outs[2];
              array[2] matrix[n_studies, n_thr] surv_prob = outs[3];
              // ////
              // fp = surv_prob[1];
              // sp = 1.0 - fp;
              // se = surv_prob[2];
              // ////
              // Fp = rowMedians(fp);
              // Sp = rowMedians(sp);
              // Se = rowMedians(se);
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

}






