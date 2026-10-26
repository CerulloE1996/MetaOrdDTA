//// DTA_MA_ord_biv_FIXEDthr_NicoStan_fused_AVX.stan
//// The model of MetaOrdDTA's DTA_MA_ord_biv_FIXEDthr.stan with NicoSIMD's AVX functions:
//// the likelihood and the induced Dirichlet of the cutpoints each in one autodiff node.
////



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
        ////
        //// ---- NicoSIMD's AVX functions (their bodies in NicoSIMD's user header):
        ////
        real nico_binomial_counts_above_successive_cutpoints_of_a_cumulative_link_lpmf(
              array[] int n, array[] int n_previous, vector latent, data vector is_first, data int link,
              data real bound);
        real nico_induced_dirichlet_prior_of_ordered_cutpoints_lpdf(vector C, vector alpha, data int link,
                                                                    data real location);
        real nico_induced_dirichlet_prior_of_ordered_cutpoints_lpdf(matrix C, vector alpha, data int link,
                                                                    data real location);
        vector nico_exp(vector x);
        vector nico_log1p_exp(vector x);
        vector nico_log_inv_logit(vector x);
        vector nico_log1m_inv_logit(vector x);
        real nico_normal_lpdf(vector y, real mu, real sigma);
        real nico_normal_lpdf(vector y, vector mu, vector sigma);
        
        // // Induced Dirichlet log-density given cutpoints C:
        // real induced_dirichlet_lpdf( vector ord_probs, 
        //                              vector C_vec, 
        //                              vector alpha, 
        //                              int use_probit,
        //                              real anchor) {
        //                                
        //           int K = num_elements(ord_probs);
        //           int n_thr = K - 1;
        //           real lp = 0.0;
        //           // Dirichlet on category probs:
        //           lp += dirichlet_lpdf(ord_probs | alpha);
        //           // Jacobian: d(cumul_prob)/d(C) = pdf(C)
        //           for (k in 1:n_thr) {
        //               if (use_probit == 1) {
        //                   lp += std_normal_lpdf(C_vec[k] - anchor);
        //               } else {
        //                   // logistic density
        //                   lp += log_inv_logit(C_vec[k] - anchor) + log1m_inv_logit(C_vec[k] - anchor);
        //               }
        //           }
        //           return lp;
        //       
        // }
        
        // real induced_dirichlet_lpdf( vector c, 
        //                              vector alpha, 
        //                              real anchor) {
        //                                
        //           int K = num_elements(c) + 1;
        //           vector[K - 1] sigma = Phi(c - anchor);
        //           vector[K] p;
        //           matrix[K, K] J = rep_matrix(0, K, K);
        //           
        //           // Induced ordinal probabilities
        //           p[1] = 1 - sigma[1];
        //           for (k in 2:(K - 1))
        //             p[k] = sigma[k - 1] - sigma[k];
        //           p[K] = sigma[K - 1];
        //           
        //           // Baseline column of Jacobian
        //           for (k in 1:K) J[k, 1] = 1;
        //           
        //           // Diagonal entries of Jacobian
        //           for (k in 2:K) {
        //             real rho = std_normal_lpdf(c - anchor[k - 1]); // sigma[k - 1] * (1 - sigma[k - 1]);
        //             J[k, k] = - rho;
        //             J[k - 1, k] = rho;
        //           }
        //           
        //           return   dirichlet_lpdf(p | alpha) + log_determinant(J);
        //           
        // }

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
        //// ---- Covariates:
        ////
        int n_covariates_nd;
        int n_covariates_d;
        int n_covariates_max;
        matrix[n_studies, n_covariates_nd] X_nd; // must be user-inoutted - study-level covariates for D-
        matrix[n_studies, n_covariates_d]  X_d;  // must be user-inoutted - study-level covariates for D+
        vector[n_covariates_nd] baseline_case_nd;  // must be user-inputted - e.g. could be {0, 1, 45.3} for 2 binary covariates and 1 cts one (e.g. age)
        vector[n_covariates_d]  baseline_case_d;   // must be user-inputted - e.g. could be {0, 1, 45.3} for 2 binary covariates and 1 cts one (e.g. age)
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
        ////
        int use_probit_link;
        ////
        array[n_studies] int<lower=0, upper=1> K_fold_CV_indicator;
}

 
transformed data {   
        int n_cat = n_thr + 1; //// Number of ordinal categories for index test 
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
        ////
        //// ---- The terms of the likelihood, one per reported cutpoint of each study and class (its counts and
        ////      its places in the locations and the cutpoints):
        ////
        int n_terms = 0;
        for (c in 1:2) {
              for (s in 1:n_studies) {
                    if (K_fold_CV_indicator[s] == 1) n_terms += n_obs_cutpoints[s];
              }
        }
        array[n_terms] int count_above_the_cutpoint_of_each_term;
        array[n_terms] int count_above_the_previous_cutpoint_of_each_term;
        vector[n_terms] is_the_first_reported_cutpoint_of_each_term;
        array[n_terms] int index_of_the_location_of_each_term;
        array[n_terms] int index_of_the_cutpoint_of_each_term;
        {
              int t = 0;
              for (c in 1:2) {
                    for (s in 1:n_studies) {
                          if (K_fold_CV_indicator[s] == 1) {
                                for (cut_i in 1:n_obs_cutpoints[s]) {
                                      int k = cutpoint_index[c, s, cut_i + 1];
                                      if (k < 1) reject("a reported cutpoint without its index: class ", c,
                                                        ", study ", s, ", reported cutpoint ", cut_i);
                                      t += 1;
                                      count_above_the_cutpoint_of_each_term[t] = x_2[c, s, cut_i];
                                      count_above_the_previous_cutpoint_of_each_term[t] = n[c, s, cut_i];
                                      is_the_first_reported_cutpoint_of_each_term[t] = (cut_i == 1) ? 1.0 : 0.0;
                                      index_of_the_location_of_each_term[t] = (c - 1) * n_studies + s;
                                      index_of_the_cutpoint_of_each_term[t] = (c - 1) * n_thr + k;
                                }
                          }
                    }
              }
        }

}


parameters { 
        matrix[2, n_covariates_max] beta_mu;       
        row_vector<lower=0.0>[2] beta_SD;    
        matrix[n_studies, 2] beta_z;    //// Study-specific random-effects (raw params)
        //// 
        real<lower=beta_corr_lb, upper=beta_corr_ub> beta_corr;  //// between-study corr (possibly restricted)  
        //// 
        array[2] vector<lower=-7.5, upper=2.5>[n_thr] C_raw_vec;   //// Global cutpoints ("raw" / unconstrained) 
        /// array[2] vector[n_thr] C_raw_vec; 
        ////
        // matrix[n_studies, 2] beta_random;
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
        cholesky_factor_corr[2] beta_L_Omega = make_bivariate_L_Omega(beta_corr);
        cholesky_factor_cov[2]  beta_L_Sigma = diag_pre_multiply(beta_SD, beta_L_Omega);
        ////
        //// ---- Study-level random effects (after Cholesky decomposition):
        ////
        // matrix[n_studies, 2] beta_random;
        // for (s in 1:n_studies) {
        //      beta_random[s, ] =  beta_z[s, ] * beta_L_Sigma;
        // }

        //// (every study's row times beta_L_Sigma at once)
        matrix[n_studies, 2] beta_random = beta_z * beta_L_Sigma;
        ////
        //// ---- Linear predictors for each disease statu + apply covariates for non-diseased and diseased groups:
        ////
        matrix[n_studies, 2] Xbeta;
        Xbeta[, 1]  = X_nd[, 1:n_covariates_nd] * to_vector(beta_mu[1, 1:n_covariates_nd]) + beta_random[, 1];
        Xbeta[, 2]  = X_d[, 1:n_covariates_d]   * to_vector(beta_mu[2, 1:n_covariates_d])  + beta_random[, 2];
} 
 
 
model { 
        ////
        //// ---- Priors: 
        ////
        // beta_mu[1, 1:n_covariates_nd] ~ normal(
        //                                 prior_beta_mu_mean[1, 1:n_covariates_nd],
        //                                 prior_beta_mu_SD[1, 1:n_covariates_nd]);
        // beta_mu[2, 1:n_covariates_d]  ~ normal(
        //                                 prior_beta_mu_mean[2, 1:n_covariates_d],
        //                                 prior_beta_mu_SD[2, 1:n_covariates_d]);
        // ////
        // beta_SD ~ normal(prior_beta_SD_mean, prior_beta_SD_SD);

        to_vector(beta_mu[1, 1:n_covariates_nd]) ~ nico_normal(
                                                    to_vector(prior_beta_mu_mean[1, 1:n_covariates_nd]),
                                                    to_vector(prior_beta_mu_SD[1, 1:n_covariates_nd]));
        to_vector(beta_mu[2, 1:n_covariates_d])  ~ nico_normal(
                                                    to_vector(prior_beta_mu_mean[2, 1:n_covariates_d]),
                                                    to_vector(prior_beta_mu_SD[2, 1:n_covariates_d]));
        ////
        to_vector(beta_SD) ~ nico_normal(prior_beta_SD_mean, prior_beta_SD_SD);
        beta_L_Omega ~ lkj_corr_cholesky(prior_beta_corr_LKJ); 
        //// 
        //// ---- Induced-dirichlet ** Prior ** model:
        //// 
        for (c in 1:2) {
             // vector[n_cat] Ind_Dir_ord_prob; //// = cumul_probs_to_ord_probs(Phi(C));
             // if (use_probit_link == 1) Ind_Dir_ord_prob = cumul_probs_to_ord_probs(Phi(C[c]));
             // else                      Ind_Dir_ord_prob = cumul_probs_to_ord_probs(inv_logit(C[c]));
             ////
             real anchor_ct = prior_beta_mu_mean[c, 1]; // prior intercepts as anchor point
             // target += induced_dirichlet_lpdf(Ind_Dir_ord_prob | C[c], prior_dirichlet_alpha[c], 1, anchor_ct);
             // target += induced_dirichlet_lpdf(C[c] | prior_dirichlet_alpha[c], anchor_ct);

             //// (as induced_dirichlet_lpdf, at the anchor prior_beta_mu_mean[c, 1], with Phi whatever the link)
             target += nico_induced_dirichlet_prior_of_ordered_cutpoints_lupdf(C[c] | prior_dirichlet_alpha[c], 1,
                                                                               prior_beta_mu_mean[c, 1]);
             target += log(n_cat);
             // Ind_Dir_ord_prob ~ induced_dirichlet_given_C(
             //                    C[c], prior_dirichlet_alpha[c], use_probit_link); // more efficient than Betancourt et al. and seems to work fine. 
        }
        //// 
        //// ---- Likelihood / Model:
        ////
        // for (s in 1:n_studies) {
        //    beta_random[s, 1:2] ~ multi_normal_cholesky(rep_vector(0.0, 2), beta_L_Sigma);
        // }
        ////
        // target += std_normal_lpdf(to_vector(beta_z));

        target += nico_normal_lupdf(to_vector(beta_z) | 0, 1);
        for (c in 1:2) {
                target += raw_C_to_C_log_det_J_lp(C_raw_vec[c], softplus);
        }
        ////
        //// ---- Log-likelihood:
        ////
        // {
        //     ////
        //     //// ---- Get the cutpoint index (k) to map "latent_surv[c][s, cut_i]" to correct cutpoint "C[k]":
        //     ////
        //     // real scale = 1.0; // since using "Xu-like"" param.
        //     array[2] matrix[n_studies, n_thr] latent_surv = map_latent_surv_to_fixed_hetero_C(
        //                                                     C, Xbeta, n_studies,
        //                                                     n_obs_cutpoints, cutpoint_index);
        //     target += compute_log_lik_binomial_fact_lp(
        //               latent_surv, use_probit_link, n_thr, n_studies,
        //               x_2, n, N_total, n_obs_cutpoints,
        //               K_fold_CV_indicator);
        // }

        {
            ////
            //// ---- every term in one autodiff node:
            ////      nico_binomial_counts_above_successive_cutpoints_of_a_cumulative_link_lpmf
            ////
            vector[2 * n_thr] cutpoints_of_both_classes = append_row(C[1], C[2]);
            vector[2 * n_studies] locations_of_both_classes = to_vector(Xbeta);
            //// (as map_latent_surv_to_fixed_hetero_C: location - C)
            vector[n_terms] latent_of_each_term =
                  locations_of_both_classes[index_of_the_location_of_each_term] -
                  cutpoints_of_both_classes[index_of_the_cutpoint_of_each_term];
            target += nico_binomial_counts_above_successive_cutpoints_of_a_cumulative_link_lupmf(
                            count_above_the_cutpoint_of_each_term |
                            count_above_the_previous_cutpoint_of_each_term, latent_of_each_term,
                            is_the_first_reported_cutpoint_of_each_term, use_probit_link, 1e-10);
        }
}


generated quantities {
      ////
      //// ---- Compute between-study variance-covariance matrix for location parameters: 
      //// 
      cov_matrix[2]  beta_Sigma = multiply_lower_tri_self_transpose(beta_L_Sigma);
      corr_matrix[2] beta_Omega = multiply_lower_tri_self_transpose(beta_L_Omega);
      ////
      //// ---- Calculate summary accuracy - "baseline" covariate values:
      ////  
      real Xbeta_baseline_nd = dot_product(baseline_case_nd, to_vector(beta_mu[1, 1:n_covariates_nd]));
      real Xbeta_baseline_d  = dot_product(baseline_case_d,  to_vector(beta_mu[2, 1:n_covariates_d]));
      ////
      //// ---- Calculate baseline Se/Sp:
      ////
      vector[n_thr] Fp_baseline = (use_probit_link == 1) ? 
                                  Phi(-(C[1]  - Xbeta_baseline_nd))  : 
                                  inv_logit(-(C[1] - Xbeta_baseline_nd));
      vector[n_thr] Sp_baseline = 1.0 - Fp_baseline;
      vector[n_thr] Se_baseline = (use_probit_link == 1) ? 
                                  Phi(-(C[2] - Xbeta_baseline_d))   : 
                                  inv_logit(-(C[2] - Xbeta_baseline_d));
      //// 
      //// ---- Calculate predictive accuracy:
      ////
      vector[2] beta_random_pred = multi_normal_cholesky_rng(rep_vector(0.0, 2), beta_L_Sigma); 
      ////
      //// ---- Use baseline covariates with predicted random effects
      ////
      real Xbeta_baseline_pred_nd = dot_product(baseline_case_nd, beta_mu[1, 1:n_covariates_nd]) + beta_random_pred[1];
      real Xbeta_baseline_pred_d  = dot_product(baseline_case_d,  beta_mu[2, 1:n_covariates_d])  + beta_random_pred[2];
      ////
      vector[n_thr] Fp_baseline_pred = (use_probit_link == 1) ? 
                                       Phi(-(C[1]  - Xbeta_baseline_pred_nd)) : 
                                       inv_logit(-(C[1] - Xbeta_baseline_pred_nd));
      vector[n_thr] Sp_baseline_pred = 1.0 - Fp_baseline_pred;
      vector[n_thr] Se_baseline_pred = (use_probit_link == 1) ? 
                                       Phi(-(C[2]  - Xbeta_baseline_pred_d))  : 
                                       inv_logit(-(C[2] - Xbeta_baseline_pred_d));
      ////
      //// ---- Log-lik + study-specific accuracy computation (using "data" / double-precision fn for efficiency):
      ////
      vector[n_studies] log_lik_study = rep_vector(0.0, n_studies);
      ////
      // matrix[n_studies, n_thr] fp; // global
      // matrix[n_studies, n_thr] sp; // global
      // matrix[n_studies, n_thr] se; // global
      // ////
      // vector[n_studies] deviance_nd; // global
      // vector[n_studies] deviance_d;  // global
      // vector[n_studies] deviance;    // global
      {
          // real scale = 1.0; // since using "Xu-like"" param.
          array[2] matrix[n_studies, n_thr] latent_surv = map_latent_surv_to_fixed_hetero_C(
                                                          C, Xbeta, n_studies, 
                                                          n_obs_cutpoints, cutpoint_index);
          ////
          array[n_studies] int dummy_CV_indicator;
          for (s in 1:n_studies) {
              dummy_CV_indicator[s] = 1;
          }
          array[3, 2] matrix[n_studies, n_thr] outs = compute_log_lik_binomial_fact_data(
                                                      latent_surv, use_probit_link, n_thr, n_studies,
                                                      x_2, n, N_total, n_obs_cutpoints,
                                                      dummy_CV_indicator);
          ////
          //// ---- Sum log-lik across thresholds and disease status within each study
          ////
          for (s in 1:n_studies) {
              for (k in 1:n_obs_cutpoints[s]) {
                log_lik_study[s] += outs[1][1][s, k];  // Non-diseased
                log_lik_study[s] += outs[1][2][s, k];  // Diseased
              }
          }
          ////
          array[2] matrix[n_studies, n_thr] cond_prob = outs[2];
          array[2] matrix[n_studies, n_thr] surv_prob = outs[3];
          // ////
          // fp = surv_prob[1];
          // sp = 1.0 - fp;
          // se = surv_prob[2];
          // ////
          // //// ---- Model fit (deviance):
          // ////
          // array[4] matrix[n_studies, n_thr] outs_model_fit = compute_deviance(
          //                                                    cond_prob, n_thr, n_studies, 
          //                                                    x_2, n, n_obs_cutpoints);
          // ////
          // matrix[n_studies, n_thr] x_hat_nd = outs_model_fit[1];
          // matrix[n_studies, n_thr] dev_nd   = outs_model_fit[2];
          // matrix[n_studies, n_thr] x_hat_d  = outs_model_fit[3];
          // matrix[n_studies, n_thr] dev_d    = outs_model_fit[4];
          // ////
          // // for (s in 1:n_studies) {
          // //     for (k in 1:to_int(n_obs_cutpoints[s])) {
          // //        deviance_nd[s] += dev_nd[s, k];
          // //        deviance_d[s]  += dev_d[s, k];
          // //     }
          // //     deviance[s] = deviance_nd[s] + deviance_d[s];
          // // }
      }

}














