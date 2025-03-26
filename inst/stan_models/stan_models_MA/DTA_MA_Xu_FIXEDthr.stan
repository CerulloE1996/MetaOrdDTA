


functions {
        ////
        //// ---- Include files to compile the necessary custom (user-defined) Stan functions:
        ////
        #include "Stan_fns_basic.stan"
        #include "Stan_fns_corr.stan"
        #include "Stan_fns_ordinal_and_cutpoints.stan"
        #include "Stan_fns_log_lik.stan"
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
        array[2] matrix[n_studies, n_thr] n;
        array[2] matrix[n_studies, n_thr] x;
        array[2] matrix[n_studies, n_thr] cutpoint_index;
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
        vector<lower=0.0>[n_thr + 1] prior_dirichlet_alpha;
        ////
        //// ---- Other:
        ////
        int<lower=0, upper=1> softplus;
    
}


transformed data { 
    
        int n_cat = n_thr + 1; //// Number of ordinal categories for index test
        // vector[n_thr] Ind_Dir_anchor = rep_vector(0.0, n_thr);
}
  


parameters {
  
        vector[2] beta_mu;    
        vector<lower=0.0>[2] beta_SD;   
        matrix[2, n_studies] beta_z;    //// Study-specific random-effects (off-centered parameterisation)
        ////
        real<lower=-1.0, upper=1.0> beta_corr;  //// between-study corr (possibly restricted)
        ////
        array[2] vector[n_thr] C_raw_vec;   //// Global cutpoints ("raw" / unconstrained)
        // array[2] ordered[n_thr] C;   //// Global cutpoints
      
}


transformed parameters { 
  
        ////
        //// ---- Construct (global) cutpoints:
        ////
        array[2] vector[n_thr] C;
        for (c in 1:2) {
            if (softplus == 1) C[c] = construct_C_using_SP_jacobian( C_raw_vec[c]);
            else               C[c] = construct_C_using_exp_jacobian(C_raw_vec[c]);
        }
        ////
        //// ---- Construct simple 2x2 (bivariate) between-study corr matrices for between-study model:
        ////
        cholesky_factor_corr[2] beta_L_Omega      = make_restricted_bivariate_L_Omega_jacobian(beta_corr, beta_corr_lb, beta_corr_ub);
        cholesky_factor_cov[2]  beta_L_Sigma      = diag_pre_multiply(beta_SD, beta_L_Omega);
        corr_matrix[2] beta_Omega = multiply_lower_tri_self_transpose(beta_L_Omega);
        ////
        //// ---- Likelihood stuff:
        ////
        array[2] matrix[n_studies, n_thr] cumul_prob = init_array_of_matrices(n_studies, n_thr, 2, 1.0); // global
        array[2] matrix[n_studies, n_thr] cond_prob  = init_array_of_matrices(n_studies, n_thr, 2, 0.0); // global
        array[2] matrix[n_studies, n_thr] log_lik    = init_array_of_matrices(n_studies, n_thr, 2, 0.0); // global
        {
            matrix[2, n_studies] beta; // local
            for (s in 1:n_studies) {
               beta[1:2, s] = beta_mu[1:2] + beta_L_Sigma * beta_z[1:2, s];
            }
            ////
            array[2] matrix[n_studies, n_thr] latent_cumul_prob = init_array_of_matrices(n_studies, n_thr, 2, positive_infinity()); // local
            ////
            //// ---- Get the cutpoint index (k) to map "latent_cumul_prob[c][s, cut_i]" to correct cutpoint "C[k]":
            ////
            real scale = 1.0; // since using Xu param.
            latent_cumul_prob = map_latent_cumul_prob_to_fixed_hetero_C(C, beta, scale, n_studies, n_obs_cutpoints, cutpoint_index);
            ////
            //// ---- Calculate CUMULATIVE probabilities (vectorised):
            ////
            for (c in 1:2) {
                cumul_prob[c] = Phi_approx(latent_cumul_prob[c]);
            }
            ////
            //// ---- Multinomial (factorised binomial likelihood)
            ////
            array[2, 2] matrix[n_studies, n_thr] log_lik_outs;
            log_lik_outs = compute_log_lik_binomial_fact(cumul_prob, x, n, n_obs_cutpoints);
            log_lik = log_lik_outs[1];
            cond_prob = log_lik_outs[2];
        }
      
}


model {

        ////
        //// ---- Priors:
        ////
        beta_mu ~ normal(prior_beta_mu_mean, prior_beta_mu_SD);
        beta_SD ~ normal(prior_beta_SD_mean, prior_beta_SD_SD);
        beta_Omega ~ lkj_corr(prior_beta_corr_LKJ);
        ////
        //// ---- Induced-dirichlet ** Prior ** model:
        ////
        for (c in 1:2) {
           vector[n_thr] Ind_Dir_cumul  = C[c];// - Ind_Dir_anchor;
           vector[n_thr] Ind_Dir_cumul_prob = Phi_approx(Ind_Dir_cumul);
           vector[n_cat] Ind_Dir_ord_prob = cumul_probs_to_ord_probs(Ind_Dir_cumul_prob);
           ////
           vector[n_thr] rho = std_normal_approx_pdf(Ind_Dir_cumul, Ind_Dir_cumul_prob);
           Ind_Dir_ord_prob ~ induced_dirichlet_given_rho(rho, prior_dirichlet_alpha);
        }
        ////
        //// ---- Likelihood / Model:
        ////
        target += std_normal_lpdf(to_vector(beta_z)); // part of between-study model, NOT prior
        ////
        //// ---- Increment the log-likelihood:
        ////
        for (c in 1:2) {
          target +=  sum(log_lik[c]);
        }
  
}


generated quantities {

          ////
          //// Calculate summary accuracy (using mean parameters):
          ////
          vector[n_thr] Fp = 1.0 - Phi_approx(C[1] - beta_mu[1]);
          vector[n_thr] Sp = 1.0 - Fp;
          vector[n_thr] Se = 1.0 - Phi_approx(C[2] - beta_mu[2]);
          ////
          //// Calculate predictive accuracy:
          ////
          vector[2] beta_pred =  multi_normal_cholesky_rng(beta_mu, beta_L_Sigma);
          ////
          vector[n_thr] Fp_pred = 1.0 - Phi_approx(C[1] - beta_pred[1]);
          vector[n_thr] Sp_pred = 1.0 - Fp_pred;
          vector[n_thr] Se_pred = 1.0 - Phi_approx(C[2] - beta_pred[2]);
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
          {
              array[2] matrix[n_studies, n_thr] x_hat = init_array_of_matrices(n_studies, n_thr, 2, -1);
              array[2] matrix[n_studies, n_thr] dev   = init_array_of_matrices(n_studies, n_thr, 2, -1);
              ////
              // corr_matrix[2] beta_Omega = multiply_lower_tri_self_transpose(beta_L_Omega);
              ////
              //// Model-predicted ("re-constructed") data:
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

}






