


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
        array[2] matrix[n_studies, n_thr] cutpoint_index;
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
        array[2] vector[n_thr + 1] prior_dirichlet_phi;
        real<lower=0.0> prior_kappa_mean;
        real<lower=0.0> prior_kappa_SD;
        // array[2] vector<lower=0.0>[n_thr + 1] prior_alpha_mean;
        // array[2] vector<lower=0.0>[n_thr + 1] prior_alpha_SD; 
        ////
        //// ---- Other:
        ////
        real<lower=0.0> kappa_lb;
        int<lower=0, upper=1> softplus;
}

transformed data { 
        int n_cat = n_thr + 1; //// Number of ordinal categories for index test
        //// matrix[n_studies, n_thr] Ind_Dir_anchor = rep_matrix(0.0, n_studies, n_thr);
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
        // array[2] vector[n_cat - 1] dirichlet_phi_raw;
        vector<lower=(kappa_lb)>[2] kappa;
        // array[2] vector<lower=alpha_lb>[n_cat] alpha;
        array[2] simplex[n_cat] dirichlet_phi;
}

 
transformed parameters { 
        vector[2] log_kappa = log(kappa);
        // array[2] simplex[n_cat] dirichlet_phi;
        array[2] vector[n_cat] log_dirichlet_phi;
         for (c in 1:2) {
        // //  dirichlet_phi[c] = stickbreaking_power_logistic_simplex_constrain_jacobian(dirichlet_phi_raw[c]); // "power logistic"
        // // dirichlet_phi[c] = stickbricking_angular_simplex_constrain_jacobian(dirichlet_phi_raw[c]); // angular
        //     //  dirichlet_phi[c] = stickbreaking_logistic_simplex_constrain_jacobian(dirichlet_phi_raw[c]); // logistic
        //    //   dirichlet_phi[c] = stickbreaking_power_normal_simplex_constrain_jacobian(dirichlet_phi_raw[c]); // "power normal"
        //      dirichlet_phi[c] = stickbreaking_normal_simplex_constrain_jacobian(dirichlet_phi_raw[c]); // normal
        //      ////
              log_dirichlet_phi[c] = log(dirichlet_phi[c]);
        }
        array[2] vector[n_cat] alpha;
        array[2] vector[n_cat] log_alpha;
        ////
        for (c in 1:2) {
            // alpha[c]  = kappa[c] * dirichlet_phi[c];
            log_alpha[c] = log(kappa[c]) + log_dirichlet_phi[c];
            //// log_alpha[c] = log_sum_exp(log_kappa[c], log_dirichlet_phi[c]);
            alpha[c] = exp(log_alpha[c]);
            // // jacobian += sum(log_alpha[c]);
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
        //// ---- ID stuff:
        ////
        array[2] matrix[n_studies, n_thr] Ind_Dir_cumul;
        array[2] matrix[n_studies, n_thr] Ind_Dir_cumul_prob;
        array[2] matrix[n_studies, n_cat] Ind_Dir_ord_prob;
        for (c in 1:2) {
           Ind_Dir_cumul[c] = C[c];
           Ind_Dir_cumul_prob[c] = Phi(Ind_Dir_cumul[c]);
           for (s in 1:n_studies) {
               Ind_Dir_ord_prob[c][s, 1:n_cat] = cumul_probs_to_ord_probs(Ind_Dir_cumul_prob[c][s, 1:n_thr]);
           }
        }
        ////
        //// ---- Likelihood stuff:
        ////
        array[2] matrix[n_studies, n_thr] surv_prob  = init_array_of_matrices(n_studies, n_thr, 2, 1.0); // global
        array[2] matrix[n_studies, n_thr] cond_prob  = init_array_of_matrices(n_studies, n_thr, 2, 0.0); // global 
        array[2] matrix[n_studies, n_thr] log_lik    = init_array_of_matrices(n_studies, n_thr, 2, 0.0); // global
        ////
        //// Construct simple 2x2 (bivariate) between-study corr matrices:
        ////
        cholesky_factor_corr[2] beta_L_Omega = make_bivariate_L_Omega(beta_corr);
        cholesky_factor_cov[2]  beta_L_Sigma = diag_pre_multiply(beta_SD, beta_L_Omega);
        ////
        {
            ////
            //// Between-study model for the location parameters ("beta") - models between-study correlation:
            ////
            matrix[n_studies, 2] beta;
            for (s in 1:n_studies) {
                 beta[s, ] = beta_mu +  beta_z[s, ] * beta_L_Sigma; 
            }
            ////
            //// ---- Map the log-lik contributions to the correct cutpoint indexes: 
            ////
            real scale = 1.0;
            array[2] matrix[n_studies, n_thr] latent_surv = map_latent_surv_to_random_hetero_C(C, n_thr, beta, scale, n_studies, n_obs_cutpoints, cutpoint_index);
            ////
            //// ---- Multinomial (factorised binomial likelihood):
            ////
            int use_probit_link = 1;
            array[2, 3] matrix[n_studies, n_thr] log_lik_outs = compute_log_lik_binomial_probit_fact(latent_surv, use_probit_link, x, n_obs_cutpoints);
            log_lik    = log_lik_outs[, 1];
            cond_prob  = log_lik_outs[, 2];
            surv_prob  = log_lik_outs[, 3];
        }
}


model {
        beta_mu ~ normal(prior_beta_mu_mean, prior_beta_mu_SD);
        beta_SD ~ normal(prior_beta_SD_mean, prior_beta_SD_SD); 
        beta_L_Omega ~ lkj_corr_cholesky(prior_beta_corr_LKJ); 
        ////
        //// ---- Induced-dirichlet between study ** model ** (NOT a prior model here but part of the actual likelihood since random-effect cutpoints!):
        ////
        for (c in 1:2) {
              
                //// can put a "flat" prior on the simplex dirichlet_cat_means_phi.
                // alpha[c] ~ normal(prior_alpha_mean[c], prior_alpha_SD[c]);
                kappa[c] ~ normal(prior_kappa_mean, prior_kappa_SD);
                dirichlet_phi[c] ~ dirichlet(prior_dirichlet_phi[c]);
                ////
                for (s in 1:n_studies) {
                      // vector[n_thr] Ind_Dir_cumul = to_vector(C[c][s, 1:n_thr]); ////- Ind_Dir_anchor;
                      // vector[n_thr] Ind_Dir_cumul_prob = to_vector(Phi(C[c][s, 1:n_thr]));
                      // vector[n_cat] Ind_Dir_ord_prob = cumul_probs_to_ord_probs(Ind_Dir_cumul_prob);
                      for (k in 1:n_thr) {
                          target += normal_lpdf(C[c][s, k] | 0.0, 1.0);
                      }
                      target += dirichlet_lpdf(to_vector(Ind_Dir_ord_prob[c][s, 1:n_cat]) | alpha[c][1:n_cat]);
                      // // row_vector[n_thr] rho = std_normal_pdf(C[c][s, ]); //// , Ind_Dir_cumul_prob[s, 1:n_thr]);
                      // target += induced_dirichlet_given_C_lpdf(Ind_Dir_ord_prob[s, 1:n_cat] | C[c][s, 1:n_thr], to_row_vector(alpha[c][1:n_cat]));
                }
        }
        ////
        //// ---- Likelihood / Model:
        ////
        target += std_normal_lpdf(to_vector(beta_z)); // part of between-study model, NOT prior // beta_{c, s} ~ normal(beta_mu_{c}, beta_SD_{c})
        ////
        //// ---- Increment the log-likelihood:
        ////
        for (c in 1:2) {
          target += raw_C_to_C_log_det_J_lp(C_raw[c], softplus);
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
          array[2] vector[n_thr] C_mu; //// = ID_cumul_probs_to_C_probit(prob_cumul_mu, 0.0, 1.0);
          for (c in 1:2) {
            vector[n_cat] prob_ord_mu   = dirichlet_phi[c]; // alpha[c] / sum(alpha[c]);  //// dirichlet_cat_means_phi[c];
            vector[n_thr] prob_cumul_mu = ord_probs_to_cumul_probs(prob_ord_mu);
            C_mu[c] = ID_cumul_probs_to_C_probit(prob_cumul_mu);
            // for (k in 1:n_thr) {
            //     C_mu[c][k] = mean(C[c][1:n_studies, k]); // ID_cumul_probs_to_C_probit(prob_cumul_mu, 0.0, 1.0);
            // }
          }
          ////
          vector[n_thr] Fp = Phi(-(C_mu[1] - beta_mu[1]));
          vector[n_thr] Sp = 1.0 - Fp;
          vector[n_thr] Se = Phi(-(C_mu[2] - beta_mu[2]));
          ////
          //// ---- Calculate predictive accuracy:
          ////
          array[2] vector[n_thr] C_pred;
          for (c in 1:2) {
            vector[n_cat] prob_ord_pred   = dirichlet_rng(alpha[c]); //// Simulate from Dirichlet by using the summary "alpha" parameters.
            vector[n_thr] prob_cumul_pred = ord_probs_to_cumul_probs(prob_ord_pred);  //// Compute PREDICTED cumulative probabilities.
            C_pred[c] = ID_cumul_probs_to_C_probit(prob_cumul_pred);  //// Compute PREDICTED cutpoints.
          }
          ////
          vector[2] beta_pred =  multi_normal_cholesky_rng(beta_mu, beta_L_Sigma);
          ////
          vector[n_thr] Fp_pred = Phi(-(C_pred[1] - beta_pred[1]));
          vector[n_thr] Sp_pred = 1.0 - Fp_pred;
          vector[n_thr] Se_pred = Phi(-(C_pred[2] - beta_pred[2]));
          ////
          //// ---- Calculate study-specific accuracy:
          ////
          matrix[n_studies, n_thr] fp = surv_prob[1];
          matrix[n_studies, n_thr] sp = 1.0 - fp;
          matrix[n_studies, n_thr] se = surv_prob[2];
          ////
          //// vector[n_thr] C_MU_empirical; //// = rowMedians(C);  //// Empirical-mean cutpoints:
          ////
          matrix[n_studies, n_thr] x_hat_nd = rep_matrix(-1, n_studies, n_thr);
          matrix[n_studies, n_thr] x_hat_d  = rep_matrix(-1, n_studies, n_thr);
          matrix[n_studies, n_thr] dev_nd   = rep_matrix(-1, n_studies, n_thr);
          matrix[n_studies, n_thr] dev_d    = rep_matrix(-1, n_studies, n_thr);
          {
                array[2] matrix[n_studies, n_thr] x_hat = init_array_of_matrices(n_studies, n_thr, 2, -1);
                array[2] matrix[n_studies, n_thr] dev   = init_array_of_matrices(n_studies, n_thr, 2, -1);
                ////
                //// Model-predicted ("re-constructed") data:
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














