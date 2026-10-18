//// DTA_MA_Gat_RANDthr_kappa_NicoStan_fused_AVX.stan
//// The same model as MetaOrdDTA's DTA_MA_Gat_RANDthr_kappa.stan (unchanged; its replaced lines kept here, commented out), for
//// NicoStan with NicoSIMD's AVX functions (Enzo's simulation study, arXiv 2505.23393): the likelihood of every
//// reported cutpoint of every study and class in ONE autodiff node (NicoSIMD's general
//// nico_binomial_counts_above_successive_cutpoints_of_a_cumulative_link_lpmf, the same node as the network
//// meta-analysis's), the induced Dirichlet of the cutpoints in one node (NicoSIMD's general
//// nico_induced_dirichlet_prior_of_ordered_cutpoints_lpdf, every set of cutpoints at once), the normal priors
//// and the standard normal random effects as nico_normal, and the exp() (or softplus) of the cutpoints' raw
//// differences as nico_exp (or nico_log1p_exp). The terms of the likelihood are listed once, in transformed
//// data. The terms that are constants of the data are left out of the log density when the sampler asks for it
//// (propto: _lupdf / _lupmf, as the ~ statements do); with propto = FALSE the log density is the model's, every
//// constant included. The generated quantities are the model's, unchanged.
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
        //// ---- NicoSIMD's AVX functions (declared without a body; stanc --allow-undefined, with NicoSIMD's
        ////      stan_external_functions_of_NicoSIMD.hpp as the user header):
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
        real nico_log_det_jacobian_of_cutpoints_from_raw_differences(vector C_raw, data int softplus);
        real nico_log_det_jacobian_of_cutpoints_from_raw_differences(matrix C_raw, data int softplus);
        array[] matrix nico_map_latent_surv_prob_to_random_C_HSROC(data matrix C, data int n_thr,
                                                                   data matrix locations, data matrix scales,
                                                                   data int n_studies,
                                                                   data array[] int n_obs_cutpoints,
                                                                   data array[,,] int cutpoint_index);
        array[,] matrix nico_compute_log_lik_binomial_fact_data(data array[] matrix latent_surv,
                                                                data int use_probit_link, data int n_thr,
                                                                data int n_studies, data array[,,] int x_2,
                                                                data array[,,] int n, data array[,] int N_total,
                                                                data array[] int n_obs_cutpoints,
                                                                data array[] int K_fold_CV_indicator);
        vector nico_Phi(vector x);
        vector nico_inv_logit(vector x);
        real nico_normal_lpdf(vector y, real mu, real sigma);
        real nico_normal_lpdf(vector y, vector mu, vector sigma);
}


data {
        ////
        //// ---- Data:
        ////
        int<lower=1> n_studies;
        int<lower=1> n_thr;
        array[n_studies] int n_obs_cutpoints;
        ////
        //// ---- Covariates:
        ////
        int n_covariates; //// For R+G - cannot have SEPERATE covariates for D+ and D-, so only one n_covariates!
        matrix[n_studies, n_covariates] X;
        vector[n_covariates] baseline_case; // must be user-inputted - e.g. could be {0, 1, 45.3} for 2 binary covariates and 1 cts one (e.g. age)
        ////
        array[2, n_studies, n_thr + 1] int x;
        array[2, n_studies, n_thr + 1] int cutpoint_index;
        ////
        //// ---- Priors for location ("beta"):
        ////
        vector[n_covariates] prior_beta_mu_mean;
        vector[n_covariates] prior_beta_mu_SD;
        real prior_beta_SD_mean;
        real prior_beta_SD_SD;
        ////
        //// ---- Priors for raw_scale ("gamma"):
        ////
        vector[n_covariates] prior_raw_scale_mu_mean;
        vector[n_covariates] prior_raw_scale_mu_SD;
        real prior_raw_scale_SD_mean;
        real prior_raw_scale_SD_SD;
        ////
        //// ---- Induced-Dirichlet priors:
        ////
        vector[n_thr + 1] prior_dirichlet_alpha;
        real prior_kappa_mean;
        real<lower=0.0> prior_kappa_SD;
        real<lower=0.0> prior_kappa_dof;
        ////
        //// ---- Other:
        ////
        real<lower=0.0> kappa_lb;
        int<lower=0, upper=1> softplus;
        ////
        int use_probit_link;
        ////
        array[n_studies] int<lower=0, upper=1> K_fold_CV_indicator;
}


transformed data {  
        int n_cat = n_thr + 1; //// Number of ordinal categories for index test
        ////
        real mult_nd = -1.0;
        real mult_d  = +1.0; //// hence: mult_d > mult_nd
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
        //// ---- The terms of the likelihood, one per reported cutpoint of each study (counted in the
        ////      cross-validation) and class, in the order of the classes, the studies and their reported
        ////      cutpoints, made once here: the count above the cutpoint, the count above the study's previous
        ////      reported cutpoint (or its total), whether the cutpoint is the study's first reported one, the
        ////      place of the study and class in to_vector() of the locations (and scales), and the place of the
        ////      cutpoint in the cutpoints of the model:
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
                                      index_of_the_cutpoint_of_each_term[t] = (k - 1) * n_studies + s;
                                }
                          }
                    }
              }
        }

}
  

parameters {
        vector[n_covariates] beta_mu;  
        real<lower=0.0> beta_SD;   
        vector[n_studies] beta_z;
        ////
        vector[n_covariates] raw_scale_mu;    
        real<lower=0.0> raw_scale_SD;   
        vector[n_studies] raw_scale_z;
        ////
        matrix<lower=-7.5, upper=2.5>[n_studies, n_thr] C_raw; //// restricting solves divs/numerical stab. issues w/ lots of sparse categories!
        real<lower=1.0> kappa;
        simplex[n_cat] dirichlet_phi; // props/probs
}


transformed parameters { 
        real log_kappa = log(kappa);
        vector<lower=0.0>[n_cat] alpha = 0.01 + dirichlet_phi * kappa;
        vector<lower=0.0>[n_cat] dirichlet_cat_SDs_sigma;
        {
           real alpha_0 = sum(alpha);
           dirichlet_cat_SDs_sigma = sqrt((alpha .* (alpha_0 - alpha) ) / 
                                     (square(alpha_0) * (alpha_0 + 1.0)));
        }
        ////
        //// ---- Construct study-specific cutpoints:
        ////
        // matrix[n_studies, n_thr] C;
        // for (s in 1:n_studies) {
        //              C[s, ] =  construct_C(C_raw[s, 1:n_thr], softplus);
        // }

        //// (as construct_C on each study's row: the first cutpoint, then the cumulative sums of the exp() (or
        ////  softplus) of the raw differences, for all studies at once)
        matrix[n_studies, n_thr] C;
        C[, 1] = C_raw[, 1];
        if (n_thr > 1) {
              vector[n_studies * (n_thr - 1)] increments_of_the_cutpoints =
                    (softplus == 1) ? nico_log1p_exp(to_vector(C_raw[, 2:n_thr])) :
                                      nico_exp(to_vector(C_raw[, 2:n_thr]));
              for (k in 2:n_thr) {
                    C[, k] = C[, k - 1] +
                              increments_of_the_cutpoints[((k - 2) * n_studies + 1):((k - 1) * n_studies)];
              }
        }
        ////
        //// ---- Between-study model for location and scale:
        ////
        vector[n_studies] beta_random = beta_SD * beta_z;
        vector[n_studies] raw_scale_random = raw_scale_SD * raw_scale_z;
        ////
        matrix[n_studies, 2] Xlocations; // NOT local
        matrix[n_studies, 2] Xscale;     // NOT local
        ////
        {
            vector[n_studies] Xbeta_baseline = X[, 1:n_covariates] * to_vector(beta_mu[1:n_covariates]) + beta_random;
            Xlocations[, 1] = mult_nd*Xbeta_baseline; 
            Xlocations[, 2] = mult_d*Xbeta_baseline; 
            ////
            vector[n_studies] Xraw_scale_baseline = X[, 1:n_covariates] * to_vector(raw_scale_mu[1:n_covariates]) + raw_scale_random;
            // Xscale[, 1] =  ((softplus == 1) ? softplus_scaled( -mult_nd*Xraw_scale_baseline) : exp( -mult_nd*Xraw_scale_baseline));
            Xscale[, 1] =  ((softplus == 1) ? softplus_scaled( -mult_nd*Xraw_scale_baseline) :
                                                nico_exp( -mult_nd*Xraw_scale_baseline)); 
            // Xscale[, 2] =  ((softplus == 1) ? softplus_scaled( -mult_d*Xraw_scale_baseline)  : exp( -mult_d*Xraw_scale_baseline));
            Xscale[, 2] =  ((softplus == 1) ? softplus_scaled( -mult_d*Xraw_scale_baseline)  :
                                                nico_exp( -mult_d*Xraw_scale_baseline));
        }
}


model {
        ////
        //// ---- Priors:
        ////
        // beta_mu[1:n_covariates] ~ normal(
        //                           prior_beta_mu_mean[1:n_covariates], 
        //                           prior_beta_mu_SD[1:n_covariates]);
        beta_mu[1:n_covariates] ~ nico_normal(prior_beta_mu_mean[1:n_covariates],
                                              prior_beta_mu_SD[1:n_covariates]);
        beta_SD ~ normal(
                  prior_beta_SD_mean, 
                  prior_beta_SD_SD);
        ////
        // raw_scale_mu[1:n_covariates] ~ normal(
        //                                prior_raw_scale_mu_mean[1:n_covariates], 
        //                                prior_raw_scale_mu_SD[1:n_covariates]);
        raw_scale_mu[1:n_covariates] ~ nico_normal(prior_raw_scale_mu_mean[1:n_covariates],
                                                   prior_raw_scale_mu_SD[1:n_covariates]);
        raw_scale_SD ~ normal(
                       prior_raw_scale_SD_mean, 
                       prior_raw_scale_SD_SD);
        ////
        {
           log_kappa ~ student_t(prior_kappa_dof, prior_kappa_mean, prior_kappa_SD);
           target += -log_kappa;
           dirichlet_phi ~ dirichlet(prior_dirichlet_alpha);
        }
        ////
        //// ---- Induced-dirichlet between study ** model ** (NOT a prior model here but part of the actual likelihood since random-effect cutpoints!):
        ////
        // for (s in 1:n_studies) {
        //      vector[n_cat] Ind_Dir_ord_prob;
        //      if (use_probit_link == 1) Ind_Dir_ord_prob = to_vector(cumul_probs_to_ord_probs(Phi(C[s, ])));
        //      else                      Ind_Dir_ord_prob = to_vector(cumul_probs_to_ord_probs(inv_logit(C[s, ])));
        //      ////
        //      Ind_Dir_ord_prob ~ induced_dirichlet_given_C(
        //                         to_vector(C[s, ]), alpha, use_probit_link);
        // }

        //// (every study's cutpoints in ONE node: NicoSIMD's general node at the location 0, with the parameters
        ////  alpha, as induced_dirichlet_given_C on each study's row)
        target += nico_induced_dirichlet_prior_of_ordered_cutpoints_lupdf(C | alpha, use_probit_link, 0.0);
        ////
        //// ---- Likelihood / Model:
        ////
        // target += std_normal_lpdf(to_vector(beta_z)); // part of between-study model, NOT prior

        target += nico_normal_lupdf(to_vector(beta_z) | 0, 1); // part of between-study model, NOT prior
        // target += std_normal_lpdf(to_vector(raw_scale_z)); // part of between-study model, NOT prior

        target += nico_normal_lupdf(to_vector(raw_scale_z) | 0, 1); // part of between-study model, NOT prior
        // target += raw_scale_to_scale_log_det_J_lp(raw_scale, softplus); // no Jacobian needed for R&G models as it cancels out (since + then -).
        // target += raw_C_to_C_log_det_J_lp(C_raw, softplus);
        target += nico_log_det_jacobian_of_cutpoints_from_raw_differences(C_raw, softplus);
        ////
        //// ---- Log-likelihood:
        ////
        // {
        //     ////
        //     //// ---- Get the cutpoint index (k) to map "latent_surv[c][s, cut_i]" to correct cutpoint "C[k]":
        //     ////
        //     array[2] matrix[n_studies, n_thr] latent_surv = map_latent_surv_prob_to_random_C_HSROC(
        //                                                     C, n_thr, Xlocations, Xscale, n_studies, 
        //                                                     n_obs_cutpoints, cutpoint_index);
        //     target += compute_log_lik_binomial_fact_lp(
        //               latent_surv, use_probit_link, n_thr, n_studies, 
        //               x_2, n, N_total, n_obs_cutpoints,
        //               K_fold_CV_indicator);
        // }

        {
            ////
            //// ---- every term in ONE autodiff node: NicoSIMD's general
            ////      nico_binomial_counts_above_successive_cutpoints_of_a_cumulative_link_lpmf (the survival
            ////      probability Phi or inv_logit of each term's latent location, its ratio to the previous
            ////      reported cutpoint's, within [1e-10, 1 - 1e-10] as safe_binomial_lpmf, and the binomial log
            ////      probability):
            ////
            vector[n_studies * n_thr] cutpoints_of_every_study = to_vector(C);
            vector[2 * n_studies] locations_of_both_classes = to_vector(Xlocations);
            vector[2 * n_studies] scales_of_both_classes = to_vector(Xscale);
            //// (as map_latent_surv_prob_to_random_C_HSROC: -(C + location) * scale)
            vector[n_terms] latent_of_each_term =
                  -(cutpoints_of_every_study[index_of_the_cutpoint_of_each_term] +
                    locations_of_both_classes[index_of_the_location_of_each_term]) .*
                  scales_of_both_classes[index_of_the_location_of_each_term];
            target += nico_binomial_counts_above_successive_cutpoints_of_a_cumulative_link_lupmf(
                            count_above_the_cutpoint_of_each_term |
                            count_above_the_previous_cutpoint_of_each_term, latent_of_each_term,
                            is_the_first_reported_cutpoint_of_each_term, use_probit_link, 1e-10);
        }
}


generated quantities {
      ////
      //// ---- SD's for each category probability in a Dirichlet is √(α_k(α_0-α_k)/(α_0^2(α_0+1))) - where α_0 is the sum of all alphas:
      ////
      vector<lower=0.0>[n_cat] dirichlet_cat_SDs_sigma_full;
      {
          real alpha_0 = sum(alpha);
          dirichlet_cat_SDs_sigma_full = sqrt((alpha .* (alpha_0 - alpha) ) / (square(alpha_0) * (alpha_0 + 1.0)));
      }
      ////
      //// ---- Calculate summary cutpoints:
      ////
      vector[n_thr] C_mu;
      ////
      // {
      //     for (k in 1:n_thr) {
      //         vector[n_studies] observed_cutpoints;
      //         int n_obs = 0;
      //         
      //         for (s in 1:n_studies) {
      //             // Search through the cutpoint_index for this study to see if cutpoint k was observed
      //             for (j in 2:(n_obs_cutpoints[s] + 1)) {  // Start at 2 because index 1 is 0
      //                 if (cutpoint_index[1, s, j] == k) {
      //                     n_obs += 1;
      //                     observed_cutpoints[n_obs] = C[s, k];
      //                     break;  // Found it, no need to continue searching
      //                 }
      //             }
      //         }
      //         
      //         if (n_obs > 0) {
      //             C_mu[k] = median(observed_cutpoints[1:n_obs]);
      //         }
      //     }
      // }
      for (k in 1:n_thr) {
        C_mu[k] = median(C[, k]);
      }
      ////
      //// ---- Calculate summary accuracy (using mean parameters):
      ////
      // vector[n_covariates] location_nd_mu = mult_nd*beta_mu;
      // vector[n_covariates] location_d_mu  = mult_d*beta_mu; 
      ////
      // vector[n_covariates] scale_nd_mu =  ((softplus == 1) ? softplus_scaled(-mult_nd*raw_scale_mu) : exp(-mult_nd*raw_scale_mu));
      // vector[n_covariates] scale_d_mu  =  ((softplus == 1) ? softplus_scaled(-mult_d*raw_scale_mu)  : exp(-mult_d*raw_scale_mu));
      ////
      real Xbeta_baseline = dot_product(baseline_case, beta_mu);
      real Xbeta_baseline_nd = mult_nd * Xbeta_baseline;
      real Xbeta_baseline_d  = mult_d  * Xbeta_baseline;
      ////
      real Xraw_scale_baseline = dot_product(baseline_case, raw_scale_mu);
      real scale_nd_baseline = ((softplus == 1) ? 
                               softplus_scaled(-mult_nd * Xraw_scale_baseline) : exp(-mult_nd * Xraw_scale_baseline));
      real scale_d_baseline  = ((softplus == 1) ? 
                               softplus_scaled(-mult_d  * Xraw_scale_baseline) : exp(-mult_d  * Xraw_scale_baseline));
      ////
      //// ---- Calculate baseline Se/Sp:
      ////
      vector[n_thr] Fp_baseline = (use_probit_link == 1) ? 
                                  nico_Phi(-(C_mu + Xbeta_baseline_nd)*scale_nd_baseline)  : 
                                  nico_inv_logit(-(C_mu + Xbeta_baseline_nd)*scale_nd_baseline);
      vector[n_thr] Sp_baseline = 1.0 - Fp_baseline;
      vector[n_thr] Se_baseline = (use_probit_link == 1) ? 
                                  nico_Phi(-(C_mu + Xbeta_baseline_d)*scale_d_baseline)    : 
                                  nico_inv_logit(-(C_mu + Xbeta_baseline_d)*scale_d_baseline);
      ////
      //// ---- Calculate predictive accuracy:
      ////
      vector[n_thr] C_pred;
      {
        for (k in 1:n_thr) {
            real sigma_C_k;
            real SD_prob_scale = dirichlet_cat_SDs_sigma_full[k];
            if (use_probit_link == 1) {
                real mean_C_cutpoint_scale = C_mu[k];
                sigma_C_k = SD_approx_ID_ord_prob_to_C_probit(mean_C_cutpoint_scale, SD_prob_scale);
            } else { 
                real mean_prob_scale = dirichlet_phi[k];
                sigma_C_k = SD_approx_ID_ord_prob_to_C_logit(mean_prob_scale, SD_prob_scale);
            }
            if (!is_nan(sigma_C_k)) {
               C_pred[k] = C_mu[k] + normal_rng(0, sigma_C_k);
            }
        }
      }
      ////
      real beta_random_pred      = normal_rng(0.0, beta_SD);
      real raw_scale_random_pred = normal_rng(0.0, raw_scale_SD);
      ////
      real Xbeta_baseline_pred = dot_product(baseline_case, beta_mu) + beta_random_pred;
      real Xbeta_baseline_pred_nd = mult_nd * Xbeta_baseline_pred;
      real Xbeta_baseline_pred_d  = mult_d  * Xbeta_baseline_pred;
      ////
      real Xraw_scale_baseline_pred = dot_product(baseline_case, raw_scale_mu) + raw_scale_random_pred;
      real scale_baseline_pred_nd = ((softplus == 1) ? 
                                    softplus_scaled(-mult_nd * Xraw_scale_baseline_pred) : exp(-mult_nd * Xraw_scale_baseline_pred));
      real scale_baseline_pred_d  = ((softplus == 1) ? 
                                    softplus_scaled(-mult_d  * Xraw_scale_baseline_pred) : exp(-mult_d  * Xraw_scale_baseline_pred));
      ////
      vector[n_thr] Fp_baseline_pred = (use_probit_link == 1) ? 
                                       nico_Phi(-(C_pred + Xbeta_baseline_pred_nd)*scale_baseline_pred_nd) : 
                                       nico_inv_logit(-(C_pred + Xbeta_baseline_pred_nd)*scale_baseline_pred_nd);
      vector[n_thr] Sp_baseline_pred = 1.0 - Fp_baseline_pred;
      vector[n_thr] Se_baseline_pred = (use_probit_link == 1) ? 
                                       nico_Phi(-(C_pred + Xbeta_baseline_pred_d)*scale_baseline_pred_d)   : 
                                       nico_inv_logit(-(C_pred + Xbeta_baseline_pred_d)*scale_baseline_pred_d);
      ////
      //// ---- Log-lik + study-specific accuracy computation (using "data" / double-precision fn for efficiency):
      ////
      vector[n_studies] log_lik_study = rep_vector(0.0, n_studies);
      ////
      // matrix[n_studies, n_thr] fp; // global
      // matrix[n_studies, n_thr] sp; // global
      // matrix[n_studies, n_thr] se; // global
      ////
      // vector[n_studies] deviance_nd; // global
      // vector[n_studies] deviance_d;  // global
      // vector[n_studies] deviance;    // global
      {
          array[2] matrix[n_studies, n_thr] latent_surv = nico_map_latent_surv_prob_to_random_C_HSROC(
                                                          C, n_thr, Xlocations, Xscale, n_studies, 
                                                          n_obs_cutpoints, cutpoint_index);
          ////
          array[n_studies] int dummy_CV_indicator;
          for (s in 1:n_studies) {
              dummy_CV_indicator[s] = 1;
          }
          array[3, 2] matrix[n_studies, n_thr] outs = nico_compute_log_lik_binomial_fact_data(
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
          ////
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
          // for (s in 1:n_studies) {
          //     for (k in 1:to_int(n_obs_cutpoints[s])) {
          //        deviance_nd[s] += is_nan(dev_nd[s, k]) ? 0.0 : dev_nd[s, k];
          //        deviance_d[s]  += is_nan(dev_d[s, k]) ? 0.0 : dev_d[s, k];
          //     }
          //     deviance[s] = deviance_nd[s] + deviance_d[s];
          // }
      }
      // ////
      // //// ---- "Ordinal-bivariate equivalent" params:
      // ////
      // array[2] matrix[n_thr, n_covariates] biv_equiv_C;
      // for (x_i in 1:n_covariates) {
      //     biv_equiv_C[1][, x_i] = to_vector(C_mu ./ scale_nd_mu[x_i]);
      //     biv_equiv_C[2][, x_i] = to_vector(C_mu ./ scale_d_mu[x_i]);
      // }
}

