


functions {
        ////
        //// Include files to compile the necessary custom (user-defined) Stan functions:
        ////
        #include "/Stan_fns_basic.stan"
        #include "Stan_fns_Pinkney_corr.stan"
        #include "Stan_fns_ordinal_and_cutpoints.stan"
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
        real<lower=1.0>  beta_corr_ub;
        real<lower=0.0>  prior_beta_corr_LKJ;
        ////
        //// Priors for cutpoints (using "Induced-Dirichlet" ordinal probs):
        ////
        vector[n_thr + 1] prior_dirichlet_cat_means_alpha;
        real<lower=0.0> prior_kappa_mean;
        real<lower=0.0> prior_kappa_SD;
        ////
        //// Other:
        ////
        real kappa_lb;
        
}


transformed data { 
    
        int n_cat = n_thr + 1; //// Number of ordinal categories for index test
        
}
  

parameters {
  
        vector[2] beta_mu;    
        vector<lower=0.0>[2] beta_SD;   
        matrix[2, n_studies] beta_z;    // Study-specific random effects for beta (off-centered parameterisation)
        ////
        real<lower=beta_corr_lb, upper=beta_corr_ub> beta_corr;  //// between-study corr (possibly restricted)
        ////
        array[n_studies] ordered[n_thr] C_array;  // study-specific cutpoints
        simplex[n_cat] dirichlet_cat_means_phi;
        real<lower=kappa_lb> kappa;  
        
}


transformed parameters { 
    
        matrix[2, n_studies] beta;
        ////
        //// Construct simple 2x2 (bivariate) between-study corr matrices:
        ////
        cholesky_factor_corr[2] beta_L_Omega = make_restricted_bivariate_L_Omega_jacobian(beta_corr, beta_corr_lb, beta_corr_ub);
        cholesky_factor_cov[2]  beta_L_Sigma = diag_pre_multiply(beta_SD, beta_L_Omega);
        ////
        matrix[n_studies, n_thr] C = convert_C_array_to_mat(C_array);
        ////
        array[2] matrix[n_studies, n_thr] logit_cumul_prob = init_array_of_matrices(n_studies, n_studies, 2, positive_infinity());
        array[2] matrix[n_studies, n_thr] cumul_prob       = init_array_of_matrices(n_studies, n_studies, 2, 1.0);
        array[2] matrix[n_studies, n_thr] cond_prob        = init_array_of_matrices(n_studies, n_studies, 2, 1.0);
        ////
        //// Induced-Dirichlet ** model ** stuff:
        ////
        matrix[n_studies, n_thr] Ind_Dir_anchor = rep_matrix(0.0, n_studies, n_thr);
        matrix[n_studies, n_thr] Ind_Dir_cumul  = C - Ind_Dir_anchor;
        matrix[n_studies, n_thr] Ind_Dir_cumul_prob = Phi(Ind_Dir_cumul);
        matrix[n_studies, n_cat] Ind_Dir_ord_prob = cumul_probs_to_ord_probs(Ind_Dir_cumul_prob);
        ////
        real log_kappa = log(kappa);
        vector[n_cat] log_dirichlet_cat_means_phi = log(dirichlet_cat_means_phi);
        ////
        //// Compute alpha:
        ////
        vector[n_cat] log_alpha = log_kappa + log_dirichlet_cat_means_phi;
        vector[n_cat] alpha = exp(log_alpha);
        real alpha_0 = sum(alpha);
        ////
        //// NOTE: SD for each category probability in a Dirichlet is √(α_k(α_0-α_k)/(α_0^2(α_0+1))) - where α_0 is the sum of all alphas:
        ////
        real alpha_0_sq = square(alpha_0);
        vector[n_cat] dirichlet_cat_SDs_sigma = sqrt((alpha .* (alpha_0 - alpha) ) ./ (alpha_0_sq * (alpha_0 + 1.0)));
        ////
        //// Between-study model for the location parameters ("beta") - models between-study correlation:
        ////
        for (s in 1:n_studies) {
             beta[, s] = beta_mu + beta_L_Sigma * beta_z[, s];
        }
        ////
        //// Cutpoints:
        ////
        for (s in 1:n_studies) {
               for (k in 1:n_thr) {
                    C[s, k] = C_array[s][k];
               }
        }
        ////
        //// Induced-Dirichlet model cumulative probs:
        ////
        {
                Ind_Dir_cumul_prob = Phi(C - Ind_Dir_anchor);
                for (s in 1:n_studies) {
                        //// Induced-Dirichlet ordinal probs:
                        Ind_Dir_ord_prob[s, 1] = Ind_Dir_cumul_prob[s, 1] - 0.0;
                        for (k in 2:n_thr) {
                           Ind_Dir_ord_prob[s, k] = Ind_Dir_cumul_prob[s, k] - Ind_Dir_cumul_prob[s, k - 1]; // since cutpoints are increasing with k
                        }
                        Ind_Dir_ord_prob[s, n_cat] = 1.0 - Ind_Dir_cumul_prob[s, n_cat - 1];
                }
        }
        ////
        //// Likelihood using binomial factorization:
        ////
        for (s in 1:n_studies) {
                  for (c in 1:2) {
                      for (cut_i in 1:to_int(n_obs_cutpoints[s])) {
                             int k = to_int(cutpoint_index[c][s, cut_i]);
                             logit_cumul_prob[c][s, cut_i] = (C[s, k] - beta[c, s]);
                      }
                  }
        }
        ////
        //// Calculate CUMULATIVE probabilities (vectorised):
        ////
        for (c in 1:2) {
            cumul_prob[c] = Phi(logit_cumul_prob[c]); //// INCREASING sequence (as C_k > C_{k - 1})
        }
      
}


model {
      
      
        ////
        //// Priors for locations:
        ////
        beta_mu ~ normal(prior_beta_mu_mean, prior_beta_mu_SD);
        beta_SD ~ normal(prior_beta_SD_mean, prior_beta_SD_SD);
        beta_Omega ~ lkj_corr(prior_beta_corr_LKJ);
        ////
        //// Priors for cutpoints / induced-Dirichlet:
        ////
        {
            target += dirichlet_lpdf(dirichlet_cat_means_phi | prior_dirichlet_cat_means_alpha ); // "flat" prior on the simplex dirichlet_cat_means_phi. 
            target += normal_lpdf(kappa | prior_kappa_mean, prior_kappa_SD);
        }
        ////
        //// Induced-dirichlet between study ** model ** (NOT a prior model here but part of the actual likelihood since random-effect cutpoints!):
        ////
        {
            for (s in 1:n_studies) {
                vector[n_thr] rho =  normal_pdf(to_vector(Ind_Dir_cumul_prob[s, 1:n_thr]), 0.0, 1.0);
                target += induced_dirichlet_v2_lpdf(to_vector(Ind_Dir_ord_prob[s, 1:n_cat]) | rho, alpha[1:n_cat]);
            }
            // if (prior_only == 0) {
                  for (k in 1:n_cat) {
                      target += log_kappa;
                      target += log_dirichlet_cat_means_phi[k];
                  }
            // }
        }
  
}

 
generated quantities {

          vector[n_thr] Se;
          vector[n_thr] Sp;
          vector[n_thr] Fp;
          vector[n_thr] Se_pred;
          vector[n_thr] Sp_pred;
          vector[n_thr] Fp_pred;
          matrix[n_studies, n_thr] se;
          matrix[n_studies, n_thr] sp;
          matrix[n_studies, n_thr] fp;
          matrix[n_studies, n_thr] x_hat_nd = rep_matrix(-1, n_studies, n_thr);
          matrix[n_studies, n_thr] x_hat_d  = rep_matrix(-1, n_studies, n_thr);
          matrix[n_studies, n_thr] dev_nd   = rep_matrix(-1, n_studies, n_thr);
          matrix[n_studies, n_thr] dev_d    = rep_matrix(-1, n_studies, n_thr);
          array[2] matrix[n_studies, n_thr] x_hat = init_array_of_matrices(n_studies, n_thr, 2, -1);
          array[2] matrix[n_studies, n_thr] dev = = init_array_of_matrices(n_studies, n_thr, 2, -1);
          vector[n_thr] C_MU;       
          vector[n_thr] C_MU_empirical; 
          ////
          //// Empirical-mean cutpoints:
          ////
          for (k in 1:n_thr) {
                    C_MU_empirical[k] = median(C[, k]);
          }
          ////
          //// Calculate study-specific accuracy:
          ////
          for (s in 1:n_studies) {
                for (k in 1:n_thr) { 
                    fp[s, k] =   1.0 - cumul_prob[1][s, k];
                    sp[s, k] =   1.0 - fp[s, k];
                    se[s, k] =   1.0 - cumul_prob[2][s, k];
                }
          }
          ////
          //// Calculate summary accuracy (using mean parameters):
          ////
          real location_nd_mu = (-1)*(-0.5)*beta_mu;
          real location_d_mu  = (-1)*(+0.5)*beta_mu;
          real scale_nd_mu = softplus_scaled(-0.5*raw_scale_mu);
          real scale_d_mu  = softplus_scaled(+0.5*raw_scale_mu); // This corresponds to the MEDIAN of log-normal (could also use mean which is e.g.,: "scale_nd = exp(-0.5*(raw_scale_mu + 0.5*raw_scale_SD"))"
          ////
          for (k in 1:n_thr) {
                Fp[k] =   1.0 - Phi((C_MU[k] - location_nd_mu)/scale_nd_mu);
                Sp[k] =   1.0 - Fp[k];
                Se[k] =   1.0 - Phi((C_MU[k] - location_d_mu)/scale_d_mu);
          }
          ////
          //// Calculate predictive accuracy:
          ////
          real beta_pred      = normal_rng(beta_mu, beta_SD);
          real raw_scale_pred = normal_rng(raw_scale_mu, raw_scale_SD); //// raw_scale_mu;
          ////
          real location_nd_pred = (-1)*(-0.5)*beta_pred;
          real location_d_pred  = (-1)*(+0.5)*beta_pred;
          real scale_nd_pred = softplus_scaled(-0.5*raw_scale_pred);
          real scale_d_pred  = softplus_scaled(+0.5*raw_scale_pred); // if using exp(), this corresponds to the MEDIAN of log-normal (could also use mean which is e.g.,: "exp(-0.5*(raw_scale_pred + 0.5*raw_scale_SD"))"
          ////
          //// Generate a prediction for each cutpoint:
          ////
          vector[n_cat] prob_ord_mu_pred  =  dirichlet_rng(alpha);   //// Simulate from Dirichlet by using the summary "alpha" parameters:
          vector[n_thr] prob_cumul_mu_pred = ord_probs_to_cumul_probs(prob_ord_mu_pred);  //// Compute PREDICTED cumulative probabilities:
          vector[n_thr] C_mu_pred = cumul_probs_to_C(prob_cumul_mu_pred, 0.0, 1.0);   //// Transform to PREDICTED cutpoints:
          ////
          for (k in 1:n_thr) {
                  Fp_pred[k] = 1.0 - Phi((C_mu_pred[k] - location_nd_pred)/scale_nd_pred);
                  Sp_pred[k] = 1.0 - Fp_pred[k];
                  Se_pred[k] = 1.0 - Phi((C_mu_pred[k] - location_d_pred)/scale_d_pred);
          }
       
}













