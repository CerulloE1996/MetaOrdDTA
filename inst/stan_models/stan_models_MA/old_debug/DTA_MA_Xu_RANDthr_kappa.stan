


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
        //// Data:
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
        //// Priors for locations:
        ////
        vector[2] prior_beta_mu_mean;
        vector[2] prior_beta_mu_SD;
        vector[2] prior_beta_SD_mean;
        vector[2] prior_beta_SD_SD;
        ////
        //// Priors for between-study correlation matrix:
        ////
        real<lower=-1.0> beta_corr_lb;
        real<lower=beta_corr_lb, upper=1.0>  beta_corr_ub;
        real<lower=0.0>  prior_beta_corr_LKJ;
        ////
        //// Priors for cutpoints (using "Induced-Dirichlet" ordinal probs):
        ////
        // array[2] vector[n_thr + 1] prior_dirichlet_cat_means_alpha;
        array[2] real<lower=0.0> prior_kappa_mean;
        array[2] real<lower=0.0> prior_kappa_SD;
        // vector[n_thr + 1] prior_dirichlet_cat_SDs_mean;
        // vector<lower=0.0>[n_thr + 1] prior_dirichlet_cat_SDs_SD;
        ////
        //// Other:
        ////
        real<lower=0.0> kappa_lb;
        int<lower=0, upper=1> softplus;
}

transformed data { 
    
        int n_cat = n_thr + 1; //// Number of ordinal categories for index test
        matrix[n_studies, n_thr] Ind_Dir_anchor = rep_matrix(0.0, n_studies, n_thr);
        int n_sets_of_C = 1;
        
}
  


parameters {
  
        vector[2] beta_mu;    
        vector<lower=0.0>[2] beta_SD;   
        matrix[2, n_studies] beta_z;    //// Study-specific random effects (off-centered parameterisation)
        ////
        real beta_corr;  //// between-study corr (possibly restricted)
        ////
        array[n_sets_of_C] matrix[n_studies, n_thr] C_raw;
        ////
        // array[n_sets_of_C] simplex[n_cat] dirichlet_cat_means_phi;
        vector[n_cat - 1] phi_raw;
        array[n_sets_of_C] real<lower=kappa_lb> kappa;
        // array[n_sets_of_C] vector<lower=kappa_lb, upper=1000.0>[n_cat] alpha;
        
}


transformed parameters { 
  
        ////
        //// ---- Construct (global) cutpoints + simplex:
        ////
        array[n_sets_of_C] matrix[n_studies, n_thr] C; ////= convert_C_array_to_mat(C_array);
        for (c in 1:n_sets_of_C) {
           for (s in 1:n_studies) {
                 vector[n_thr] C_raw_vec_temp = to_vector(C_raw[c][s, 1:n_thr]);
                 vector[n_thr] C_vec_temp = ((softplus == 1) ? construct_C_using_SP_jacobian(C_raw_vec_temp) : construct_C_using_exp_jacobian(C_raw_vec_temp));
                 C[c][s, ] = to_row_vector(C_vec_temp);
           }
        }
         array[n_sets_of_C] simplex[n_cat] dirichlet_cat_means_phi = stickbreaking_logistic_simplex_constrain(phi_raw);
        ////
        //// ---- Likelihood stuff:
        ////
        array[2] matrix[n_studies, n_thr] cumul_prob = init_array_of_matrices(n_studies, n_thr, 2, 1.0); // global
        array[2] matrix[n_studies, n_thr] cond_prob  = init_array_of_matrices(n_studies, n_thr, 2, 0.0); // global
        array[2] matrix[n_studies, n_thr] log_lik    = init_array_of_matrices(n_studies, n_thr, 2, 0.0); // global
        ////
        //// Construct simple 2x2 (bivariate) between-study corr matrices:
        ////
        corr_matrix[2] beta_Omega = make_restricted_bivariate_Omega_jacobian(beta_corr, beta_corr_lb, beta_corr_ub);
        cholesky_factor_corr[2] beta_L_Omega = cholesky_decompose(beta_Omega);
        cholesky_factor_cov[2]  beta_L_Sigma = diag_pre_multiply(beta_SD, beta_L_Omega);
        ////
        array[n_sets_of_C] vector[n_cat] alpha; //// = kappa .* dirichlet_cat_means_phi;
        // // array[2] real log_kappa; //// = log(kappa);
        for (c in 1:n_sets_of_C) {
           // log_kappa[c] = log(kappa[c]);
           for (k in 1:n_cat) {
               alpha[c][k] = kappa[c] * dirichlet_cat_means_phi[c][k];
           }
        }
        // vector[n_cat] log_dirichlet_cat_means_phi = log(dirichlet_cat_means_phi);
        // vector[n_cat] log_alpha = log_kappa + log_dirichlet_cat_means_phi;
        // vector[n_cat] alpha = exp(log_alpha);
        ////
        //// SD's for each category probability in a Dirichlet is √(α_k(α_0-α_k)/(α_0^2(α_0+1))) - where α_0 is the sum of all alphas:
        ////
        // array[n_sets_of_C] vector<lower=0.0>[n_cat] dirichlet_cat_SDs_sigma;
        // for (c in 1:n_sets_of_C) {
        //   real alpha_0 = sum(alpha[c]);
        //   dirichlet_cat_SDs_sigma[c] = convert_Dirichlet_alpha_to_SDs_jacobian(alpha[c]);
        //   //// dirichlet_cat_SDs_sigma[c] = sqrt((alpha[c] .* (alpha_0 - alpha[c]) ) ./ (square(alpha_0) * (alpha_0 + 1.0)));
        // }
        // // ////
        {
                array[2] matrix[n_studies, n_thr] latent_cumul_prob = init_array_of_matrices(n_studies, n_thr, 2, positive_infinity());
                ////
                //// Between-study model for the location parameters ("beta") - models between-study correlation:
                ////
                matrix[2, n_studies] beta;
                for (s in 1:n_studies) {
                     beta[, s] = beta_mu + beta_L_Sigma * beta_z[, s];
                }
                ////
                //// ---- Map the log-lik contributions to the correct cutpoint indexes:
                ////
                real scale = 1.0;
                latent_cumul_prob = map_latent_cumul_prob_to_random_C(C[1], beta, scale, n_studies, n_thr, n_obs_cutpoints, cutpoint_index);
                ////
                //// ---- Calculate CUMULATIVE probabilities (vectorised):
                ////
                for (c in 1:2) {
                    cumul_prob[c] = Phi_approx(latent_cumul_prob[c]); //// INCREASING sequence (as C_k > C_{k - 1})
                }
                ////
                //// ---- Multinomial (factorised Binomial) likelihood:
                ////
                array[2, 2] matrix[n_studies, n_thr] log_lik_outs;
                log_lik_outs = compute_log_lik_binomial_fact(cumul_prob, x, n, n_obs_cutpoints);
                log_lik = log_lik_outs[1];
                cond_prob = log_lik_outs[2];
        }
        
}


model {
      
        {
            //// locations:
            beta_mu ~ normal(prior_beta_mu_mean, prior_beta_mu_SD);
            beta_SD ~ normal(prior_beta_SD_mean, prior_beta_SD_SD);
            beta_Omega ~ lkj_corr(prior_beta_corr_LKJ);
        }
        //// Induced-dirichlet between study ** model ** (NOT a prior model here but part of the actual likelihood since random-effect cutpoints!):
        for (c in 1:n_sets_of_C) {
          
                matrix[n_studies, n_thr] Ind_Dir_cumul = C[c] - Ind_Dir_anchor;
                matrix[n_studies, n_thr] Ind_Dir_cumul_prob = Phi_approx(Ind_Dir_cumul);
                matrix[n_studies, n_cat] Ind_Dir_ord_prob = cumul_probs_to_ord_probs(Ind_Dir_cumul_prob);
                ////
                //// "flat" prior on the simplex dirichlet_cat_means_phi.
                //// //// dirichlet_cat_means_phi ~ dirichlet(prior_dirichlet_cat_means_alpha); ///// target += dirichlet_lpdf( dirichlet_cat_means_phi | prior_dirichlet_cat_means_alpha ); 
                kappa[c] ~ normal(0.0, prior_kappa_SD[1]);
                //// dirichlet_cat_SDs_sigma[c] ~ normal(prior_dirichlet_cat_SDs_mean, prior_dirichlet_cat_SDs_SD);
                //// alpha[c] ~ normal(0.0, prior_kappa_SD[1]);
                for (k in 1:n_cat) {
                  target += log(abs(dirichlet_cat_means_phi[c][k]));
                  target += log(abs(kappa[c]));
                }
                ////
                for (s in 1:n_studies) {
                    vector[n_thr] rho =  std_normal_approx_pdf(to_vector(Ind_Dir_cumul[s, 1:n_thr]), to_vector(Ind_Dir_cumul_prob[s, 1:n_thr]));
                    target += induced_dirichlet_given_rho_lpdf(to_vector(Ind_Dir_ord_prob[s, 1:n_cat]) | rho, alpha[c][1:n_cat]);
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
          target +=  sum(log_lik[c]);
        }

}


generated quantities {

          ////
          //// Calculate summary accuracy (using mean parameters):
          ////
          vector[n_thr] Fp; //// = 1.0 - Phi_approx(C_mu - beta_mu[1]);
          vector[n_thr] Sp; //// = 1.0 - Fp;
          vector[n_thr] Se; //// = 1.0 - Phi_approx(C_mu - beta_mu[2]);
          ////
          for (c in 1:n_sets_of_C) {
              vector[n_cat] prob_ord_mu   = alpha[c] / sum(alpha[c]);  //// dirichlet_cat_means_phi[c];
              vector[n_thr] prob_cumul_mu = ord_probs_to_cumul_probs(prob_ord_mu);
              vector[n_thr] C_mu = cumul_probs_to_C(prob_cumul_mu, 0.0, 1.0);
              ////
              Fp = 1.0 - Phi_approx(C_mu - beta_mu[1]);
              Sp = 1.0 - Fp;
              Se = 1.0 - Phi_approx(C_mu - beta_mu[2]);
          }
          ////
          //// Calculate predictive accuracy:
          ////
          vector[n_thr] Fp_pred; //// = 1.0 - Phi_approx(C_mu - beta_mu[1]);
          vector[n_thr] Sp_pred; //// = 1.0 - Fp;
          vector[n_thr] Se_pred; //// = 1.0 - Phi_approx(C_mu - beta_mu[2]);
          ////
          for (c in 1:n_sets_of_C) {
              vector[n_cat] prob_ord_pred   = dirichlet_rng(alpha[c]); //// Simulate from Dirichlet by using the summary "alpha" parameters.
              vector[n_thr] prob_cumul_pred = ord_probs_to_cumul_probs(prob_ord_pred);  //// Compute PREDICTED cumulative probabilities.
              vector[n_thr] C_pred = cumul_probs_to_C(prob_cumul_pred, 0.0, 1.0);  //// Compute PREDICTED cutpoints.
              ////
              vector[2] beta_pred =  multi_normal_cholesky_rng(beta_mu, beta_L_Sigma);
              ////
              Fp_pred = 1.0 - Phi_approx(C_pred - beta_pred[1]);
              Sp_pred = 1.0 - Fp_pred;
              Se_pred = 1.0 - Phi_approx(C_pred - beta_pred[2]);
          }
          ////
          //// Calculate study-specific accuracy:
          ////
          matrix[n_studies, n_thr] fp = 1.0 - cumul_prob[1];
          matrix[n_studies, n_thr] sp = 1.0 - fp;
          matrix[n_studies, n_thr] se = 1.0 - cumul_prob[2];
          // ////
          // //// Compute between-study correlation matrix for location parameters:
          // ////
          // corr_matrix[2] beta_Omega = multiply_lower_tri_self_transpose(beta_L_Omega);
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














