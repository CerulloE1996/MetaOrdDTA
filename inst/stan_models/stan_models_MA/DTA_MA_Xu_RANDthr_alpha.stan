


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
        //// ---- Priors for locations ("beta"):
        ////
        vector[2] prior_beta_mu_mean;
        vector[2] prior_beta_mu_SD;
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
        // vector[n_thr + 1] prior_dirichlet_cat_means_alpha;
        array[2] vector<lower=0.0>[n_thr + 1] prior_alpha_mean;
        array[2] vector<lower=0.0>[n_thr + 1] prior_alpha_SD; 
        ////
        //// ---- Other:
        ////
        real<lower=0.0> alpha_lb;
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
}
  


parameters {
        row_vector[2] beta_mu;    
        row_vector<lower=0.0>[2] beta_SD;   
        matrix[n_studies, 2] beta_z;    //// Study-specific random effects (off-centered parameterisation)
        ////
        real<lower=beta_corr_lb, upper=beta_corr_ub> beta_corr;  //// between-study corr (possibly restricted)
        ////
        array[2] matrix[n_studies, n_thr] C_raw;
        ////
        array[2] vector<lower=alpha_lb>[n_cat] alpha;
}

 
transformed parameters { 
         // array[2] vector[n_cat] alpha;
         // for (c in 1:2) {
         //    alpha[c] = log1p_exp(soft_alpha[c]);
         // }
        // array[2] vector[n_cat] log_dirichlet_phi;
        //  for (c in 1:2) {
        // // //  dirichlet_phi[c] = stickbreaking_power_logistic_simplex_constrain_jacobian(dirichlet_phi_raw[c]); // "power logistic"
        //    dirichlet_phi[c] = stickbricking_angular_simplex_constrain_jacobian(dirichlet_phi_raw[c]); // angular
        // //     //  dirichlet_phi[c] = stickbreaking_logistic_simplex_constrain_jacobian(dirichlet_phi_raw[c]); // logistic
        // //    //   dirichlet_phi[c] = stickbreaking_power_normal_simplex_constrain_jacobian(dirichlet_phi_raw[c]); // "power normal"
        // //      dirichlet_phi[c] = stickbreaking_normal_simplex_constrain_jacobian(dirichlet_phi_raw[c]); // normal
        // //      ////
        //       log_dirichlet_phi[c] = log(dirichlet_phi[c]);
        // }
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
        //// ---- Between-study model for the location parameters ("beta") - models between-study correlation:
        ////
        matrix[n_studies, 2] beta;
        for (s in 1:n_studies) {
             beta[s, ] = beta_mu +  beta_z[s, ] * beta_L_Sigma; 
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
        for (c in 1:2) {
              alpha[c] ~ normal(prior_alpha_mean[c], prior_alpha_SD[c]);  // ID prior
        }
        ////
        //// ---- Induced-dirichlet between study ** model ** (NOT a prior model here but part of the actual likelihood since random-effect cutpoints!):
        ////
        for (c in 1:2) {
             for (s in 1:n_studies) {
                 // vector[n_cat] Ind_Dir_ord_prob;
                 // if (use_probit_link == 1) Ind_Dir_ord_prob = to_vector(cumul_probs_to_ord_probs(Phi(C[c][s, ])));
                 // else                      Ind_Dir_ord_prob = to_vector(cumul_probs_to_ord_probs(inv_logit(C[c][s, ])));
                 // ////
                 // //// Ind_Dir_ord_prob ~ induced_dirichlet_given_C(to_vector(C[c][s, ]), alpha[c], use_probit_link);
                 // // Ind_Dir_ord_prob ~ induced_dirichlet_given_C(to_vector(C[c][s, ]), alpha[c], use_probit_link);
                 // vector[n_thr] rho = to_vector(std_normal_pdf(C[c][s, ]));
                 // Ind_Dir_ord_prob ~ induced_dirichlet_given_rho(rho, alpha[c]);
                 ////
                 to_vector(C[c][s,]) ~ induced_dirichlet(alpha[c], 0.0);
             }
        }
        ////
        //// ---- Likelihood / Model:
        ////
        target += std_normal_lpdf(to_vector(beta_z)); // part of between-study model, NOT prior // beta_{c, s} ~ normal(beta_mu_{c}, beta_SD_{c})
        ////
        for (c in 1:2) {
          target += raw_C_to_C_log_det_J_lp(C_raw[c], softplus);
        }
        ////
        //// ---- Log-likelihood:
        ////
        {
            ////
            //// ---- Get the cutpoint index (k) to map "latent_surv[c][s, cut_i]" to correct cutpoint "C[k]":
            ////
            real scale = 1.0; // since using "Xu-like"" param.
            array[2] matrix[n_studies, n_thr] latent_surv = map_latent_surv_to_random_hetero_C(C, n_thr, beta, scale, n_studies, n_obs_cutpoints, cutpoint_index);
            target += compute_log_lik_binomial_fact_lp(latent_surv, use_probit_link, x_2, n, N_total, n_obs_cutpoints);
        }
}


generated quantities {
          ////
          //// ---- Compute between-study variance-covariance matrix for location parameters:
          ////
          cov_matrix[2]  beta_Sigma = multiply_lower_tri_self_transpose(beta_L_Sigma);
          corr_matrix[2] beta_Omega = multiply_lower_tri_self_transpose(beta_L_Omega);
          ////
          //// SD's for each category probability in a Dirichlet is √(α_k(α_0-α_k)/(α_0^2(α_0+1))) - where α_0 is the sum of all alphas:
          ////
          array[2] vector<lower=0.0>[n_cat] dirichlet_cat_SDs_sigma;
          for (c in 1:2) {
              real alpha_0 = sum(alpha[c]);
              dirichlet_cat_SDs_sigma[c] = sqrt((alpha[c] .* (alpha_0 - alpha[c]) ) / (square(alpha_0) * (alpha_0 + 1.0)));
          }
          ////
          //// ---- Calculate summary accuracy (using mean parameters):
          ////
          array[2] vector[n_thr] C_mu;
          // for (c in 1:2) {
          //   vector[n_cat] prob_ord_mu   = alpha[c]/sum(alpha[c]);  //// dirichlet_cat_means_phi[c];
          //   vector[n_thr] prob_cumul_mu = ord_probs_to_cumul_probs(prob_ord_mu);
          //   C_mu[c] = ID_cumul_probs_to_C(prob_cumul_mu, use_probit_link);
          //   // for (k in 1:n_thr) {
          //   //     C_mu[c][k] = mean(C[c][1:n_studies, k]); // ID_cumul_probs_to_C(prob_cumul_mu, 0.0, 1.0);
          //   // }
          // }
          ////
          for (c in 1:2) {
               matrix[1000, n_thr + 1] p_dm_sim;
               vector[n_thr + 1] p_dm;
               ////
               for (i in 1:1000) 
                   p_dm_sim[i,]  =  to_row_vector(dirichlet_rng(alpha[c]));
                  
               for (i in 1:(n_thr+1)) {
                   p_dm[i] = mean(p_dm_sim[,i]); 
               }
          
                 C_mu[c][1] =   inv_Phi(p_dm_sim[1, 1]); 
                for (i in 2:n_thr) {
                  C_mu[c][i] =    inv_Phi(p_dm_sim[1, i] + Phi(C_mu[c][i - 1]));
                }
          }
          ////
          vector[n_thr] Fp = (use_probit_link == 1) ? Phi(-(C_mu[1] - beta_mu[1])) : inv_logit(-(C_mu[1] - beta_mu[1]));
          vector[n_thr] Sp = 1.0 - Fp;
          vector[n_thr] Se = (use_probit_link == 1) ? Phi(-(C_mu[2] - beta_mu[2])) : inv_logit(-(C_mu[2] - beta_mu[2]));
          ////
          //// ---- Calculate predictive accuracy:
          ////
          array[2] vector[n_thr] C_pred;
          for (c in 1:2) {
            vector[n_cat] prob_ord_pred   = dirichlet_rng(alpha[c]); //// Simulate from Dirichlet by using the summary "alpha" parameters.
            vector[n_thr] prob_cumul_pred = ord_probs_to_cumul_probs(prob_ord_pred);  //// Compute PREDICTED cumulative probabilities.
            C_pred[c] = ID_cumul_probs_to_C(prob_cumul_pred, use_probit_link);  //// Compute PREDICTED cutpoints.
          }
          ////
          vector[2] beta_pred =  multi_normal_cholesky_rng(beta_mu, beta_L_Sigma);
          ////
          vector[n_thr] Fp_pred = (use_probit_link == 1) ? Phi(-(C_pred[1] - beta_pred[1])) : inv_logit(-(C_pred[1] - beta_pred[1]));
          vector[n_thr] Sp_pred = 1.0 - Fp_pred;
          vector[n_thr] Se_pred = (use_probit_link == 1) ? Phi(-(C_pred[2] - beta_pred[2])) : inv_logit(-(C_pred[2] - beta_pred[2]));
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
              array[2] matrix[n_studies, n_thr] latent_surv = map_latent_surv_to_random_hetero_C(C, n_thr, beta, scale, n_studies, n_obs_cutpoints, cutpoint_index);
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

}






