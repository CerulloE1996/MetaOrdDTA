


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
        // array[2] matrix[n_studies, n_thr] x_with_missings;
        array[2] matrix[n_studies, n_thr + 1] x;  
        array[2] matrix[n_studies, n_thr + 1] cutpoint_index;
        ////
        //// ---- Priors for locations:
        ////
        vector[2] prior_beta_mu_mean;
        vector[2] prior_beta_mu_SD;  
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
        ////
        int compute_sim_study_metrics; // for "internal"" use only.
        int n_inner_thr;
        int n_outer_thr;
        array[n_inner_thr] int vec_index_inner_thr;
        array[n_outer_thr] int vec_index_outer_thr;
}

 
transformed data {  
        int n_cat = n_thr + 1; //// Number of ordinal categories for index test
        // vector[n_thr] Ind_Dir_anchor = rep_vector(0.0, n_thr);} 
        int use_probit_link = 1;
}

 
parameters { 
        row_vector[2] beta_mu;       
        row_vector<lower=0.0>[2] beta_SD;    
        matrix[n_studies, 2] beta_z;    //// Study-specific random-effects (off-centered parameterisation) 
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
            array[2] matrix[n_studies, n_thr] surv_prob   = init_array_of_matrices(n_studies, n_thr, 2, 1.0); 
            array[2] matrix[n_studies, n_thr] cond_prob   = init_array_of_matrices(n_studies, n_thr, 2, 0.0);
            array[2] matrix[n_studies, n_thr] log_lik     = init_array_of_matrices(n_studies, n_thr, 2, 0.0);
            ////
            {
                  matrix[n_studies, 2] beta; // local
                  for (s in 1:n_studies) {  
                     beta[s, ] = beta_mu + beta_z[s, ] * beta_L_Sigma;
                  }
                  ////
                  //// ---- Get the cutpoint index (k) to map "latent_surv[c][s, cut_i]" to correct cutpoint "C[k]":
                  ////
                  real scale = 1.0; // since using "Xu-like"" param.
                  array[2] matrix[n_studies, n_thr] latent_surv = map_latent_surv_prob_to_fixed_hetero_C(C, beta, scale, n_studies, n_obs_cutpoints, cutpoint_index);
                  ////
                  //// ---- Multinomial (factorised binomial likelihood)
                  ////
                  array[2, 3] matrix[n_studies, n_thr] log_lik_outs = compute_log_lik_binomial_probit_fact(latent_surv, use_probit_link, x, n_obs_cutpoints);
                  log_lik    = log_lik_outs[, 1];
                  cond_prob  = log_lik_outs[, 2];
                  surv_prob  = log_lik_outs[, 3];
            }
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
           // vector[n_thr] log_J_p_to_C; // Jacobian matrix for p_ord_ID -> C transformation is diagonal
           //    for (k in 1:n_thr) {
           //       log_J_p_to_C[k] = - std_normal_lpdf(C[c][k]);
           //    }
           //    target += sum(log_J_p_to_C);  // Jacobian matrix for p_ord_ID -> C transformation is diagonal
           vector[n_thr] Ind_Dir_cumul  = C[c];// - Ind_Dir_anchor; 
           vector[n_thr] Ind_Dir_cumul_prob = Phi(C[c]);
           vector[n_cat] Ind_Dir_ord_prob = cumul_probs_to_ord_probs(Ind_Dir_cumul_prob);
           ////
           Ind_Dir_ord_prob ~ induced_dirichlet_given_C(C[c], prior_dirichlet_alpha[c], use_probit_link);
         
        } 
        //// 
        //// ---- Likelihood / Model:
        ////
        target += std_normal_lpdf(to_vector(beta_z)); // part of between-study model, NOT prior
        ////
        //// ---- Increment the log-likelihood:
        ////
        for (c in 1:2) {
          target += raw_C_to_C_log_det_J_lp(C_raw_vec[c], softplus);
          target += sum(log_lik[c]);
        }
}


generated quantities {
          ////
          //// ---- Compute between-study variance-covariance matrix for location parameters:
          ////
          cov_matrix[2]  beta_Sigma = multiply_lower_tri_self_transpose(beta_L_Sigma);
          corr_matrix[2] beta_Omega = multiply_lower_tri_self_transpose(beta_L_Omega);
          ////
          //// ---- Calculate summary accuracy (using mean parameters):
          ////
          vector[n_thr] Fp = Phi(-(C[1] - beta_mu[1]));
          vector[n_thr] Sp = 1.0 - Fp;
          vector[n_thr] Se = Phi(-(C[2] - beta_mu[2]));
          // for (k in 1:(n_thr + 1)) {
          //     Fp[k] = mean(surv_prob[1][, k]);
          //     Sp[k] = 1.0 - Fp[k];
          //     Se[k] = mean(surv_prob[2][, k]);
          // }
          ////
          //// ---- Calculate predictive accuracy:
          ////
          vector[2] beta_pred =  multi_normal_cholesky_rng(beta_mu, beta_L_Sigma);
          ////
          vector[n_thr] Fp_pred = Phi(-(C[1] - beta_pred[1]));
          vector[n_thr] Sp_pred = 1.0 - Fp_pred;
          vector[n_thr] Se_pred = Phi(-(C[2] - beta_pred[2]));
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
                     for (i in 1:to_int(n_obs_cutpoints[s])) {
    
                          //// Model-estimated data:
                          x_hat[c][s, i] = cond_prob[c][s, i] * x[c][s, i];  	 // Fitted values
    
                          //// Compute residual deviance contribution:
                          real n_i =  (x[c][s, i]);
                          real x_i =  (x[c][s, i + 1]);
                          real x_hat_i =  (x_hat[c][s, i]);
                          real log_x_minus_log_x_hat = log(x_i) - log(x_hat_i);
                          real log_diff_n_minus_x = log(n_i - x_i);
                          real log_diff_n_minus_x_hat = log(abs(n_i - x_hat_i));
    
                          dev[c][s, i] = 2.0 * ( x_i * log_x_minus_log_x_hat + (n_i - x_i) * (log_diff_n_minus_x - log_diff_n_minus_x_hat) );
    
                     }
                  }
              }
    
              x_hat_nd = x_hat[1];
              dev_nd = dev[1];
              x_hat_d = x_hat[2];
              dev_d = dev[2];
          }
         ////
         //// ---- Measures for simulation study ("internal" use only!):
         ////
         ////
         real Se_mean_inner_thr = 0.0;
         real Sp_mean_inner_thr = 0.0;
         real Se_mean_outer_thr = 0.0;
         real Sp_mean_outer_thr = 0.0;
         real Se_mean_overall_measure = 0.0;
         real Sp_mean_overall_measure = 0.0;
         ////
         if (compute_sim_study_metrics == 1) {
               ////
               Se_mean_inner_thr = mean(Se[vec_index_inner_thr]);
               Sp_mean_inner_thr = mean(Sp[vec_index_inner_thr]);
               ////
               Se_mean_outer_thr = mean(Se[vec_index_outer_thr]);
               Sp_mean_outer_thr = mean(Sp[vec_index_outer_thr]);
               ////
               Se_mean_overall_measure = 0.75*Se_mean_inner_thr + 0.25*Se_mean_outer_thr;
               Sp_mean_overall_measure = 0.75*Sp_mean_inner_thr + 0.25*Sp_mean_outer_thr;
         }

}






