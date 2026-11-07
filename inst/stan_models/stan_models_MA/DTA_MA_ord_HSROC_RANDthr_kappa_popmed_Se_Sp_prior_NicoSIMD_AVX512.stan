// NicoSIMD_AVX512 version of DTA_MA_ord_HSROC_RANDthr_kappa_popmed_Se_Sp_prior.stan written by NicoStan (source model md5 f2afc4ea88fa3b9aeb0ce9f790315f60; NicoSIMD declarations md5 e3c567fd115d2ccf49ec504688009e6b; writer implementation md5 aef1406a74d6694a1d260c4cd31c7a35)
// NicoSIMD functions used: exp, log, log1m, log1p_exp, log1m_exp, inv_logit, log_inv_logit, Phi, inv_Phi, tanh, logit, log1m_inv_logit, normal_rng, normal, lgamma, beta, dirichlet, log_diff_exp, std_normal_log_qf, student_t
//// DTA_MA_ord_HSROC_RANDthr_kappa_popmed_Se_Sp_prior_Stan_math_only.stan
//// The same data, parameters and target as the fused AVX file of this model with the Dirichlet prior of dirichlet_phi, on Stan::math only.
////

functions {
    real bmvp_simd_lanes();
    real nico_exp(real x);
    real nico_log(real x);
    real nico_log1m(real x);
    real nico_log1p_exp(real x);
    real nico_log1m_exp(real x);
    real nico_inv_logit(real x);
    real nico_log_inv_logit(real x);
    real nico_Phi(real x);
    real nico_inv_Phi(real x);
    real nico_tanh(real x);
    real nico_logit(real x);
    real nico_log1m_inv_logit(real x);
    vector nico_exp(vector x);
    vector nico_log(vector x);
    vector nico_log1m(vector x);
    vector nico_log1p_exp(vector x);
    vector nico_log1m_exp(vector x);
    vector nico_inv_logit(vector x);
    vector nico_log_inv_logit(vector x);
    vector nico_Phi(vector x);
    vector nico_inv_Phi(vector x);
    vector nico_tanh(vector x);
    vector nico_logit(vector x);
    vector nico_log1m_inv_logit(vector x);
    row_vector nico_exp(row_vector x);
    row_vector nico_log(row_vector x);
    row_vector nico_log1m(row_vector x);
    row_vector nico_log1p_exp(row_vector x);
    row_vector nico_log1m_exp(row_vector x);
    row_vector nico_inv_logit(row_vector x);
    row_vector nico_log_inv_logit(row_vector x);
    row_vector nico_Phi(row_vector x);
    row_vector nico_inv_Phi(row_vector x);
    row_vector nico_tanh(row_vector x);
    row_vector nico_logit(row_vector x);
    row_vector nico_log1m_inv_logit(row_vector x);
    matrix nico_exp(matrix x);
    matrix nico_log(matrix x);
    matrix nico_log1m(matrix x);
    matrix nico_log1p_exp(matrix x);
    matrix nico_log1m_exp(matrix x);
    matrix nico_inv_logit(matrix x);
    matrix nico_log_inv_logit(matrix x);
    matrix nico_Phi(matrix x);
    matrix nico_inv_Phi(matrix x);
    matrix nico_tanh(matrix x);
    matrix nico_logit(matrix x);
    matrix nico_log1m_inv_logit(matrix x);
    real nico_normal_rng(real mu, real sigma);
    array[] real nico_normal_rng(vector mu, real sigma);
    array[] real nico_normal_rng(real mu, vector sigma);
    array[] real nico_normal_rng(vector mu, vector sigma);
    array[] real nico_normal_rng(row_vector mu, real sigma);
    array[] real nico_normal_rng(array[] real mu, real sigma);
    real nico_normal_lpdf(vector y, vector mu, vector sigma);
    real nico_normal_lpdf(vector y, vector mu, real sigma);
    real nico_normal_lpdf(vector y, real mu, vector sigma);
    real nico_normal_lpdf(vector y, real mu, real sigma);
    real nico_normal_lpdf(row_vector y, row_vector mu, row_vector sigma);
    real nico_normal_lpdf(array[] real y, vector mu, vector sigma);
    real nico_normal_lpdf(array[] real y, vector mu, real sigma);
    real nico_normal_lpdf(array[] real y, real mu, vector sigma);
    real nico_normal_lpdf(array[] real y, real mu, real sigma);
    real nico_normal_lpdf(real y, real mu, real sigma);
    real nico_lgamma(data real x);
    vector nico_lgamma(vector x);
    real nico_beta_lpdf(real y, data real alpha, data real beta);
    real nico_beta_lpdf(vector y, vector alpha, vector beta);
    real nico_dirichlet_lpdf(vector theta, vector alpha);
    real nico_log_diff_exp(real x, real y);
    vector nico_log_diff_exp(vector x, vector y);
    real nico_std_normal_log_qf(real x);
    vector nico_std_normal_log_qf(vector x);
    real nico_student_t_lpdf(vector y, real nu, vector mu, vector sigma);
    real nico_student_t_lpdf(vector y, real nu, vector mu, real sigma);
    real nico_student_t_lpdf(vector y, real nu, real mu, vector sigma);
    real nico_student_t_lpdf(vector y, real nu, real mu, real sigma);
    real nico_student_t_lpdf(array[] real y, real nu, vector mu, vector sigma);
    real nico_student_t_lpdf(array[] real y, real nu, vector mu, real sigma);
    real nico_student_t_lpdf(array[] real y, real nu, real mu, vector sigma);
    real nico_student_t_lpdf(array[] real y, real nu, real mu, real sigma);
    real nico_student_t_lpdf(real y, real nu, real mu, real sigma);
        ////
        //// ---- Include files to compile the necessary custom (user-defined) Stan functions:
        ////
        #include "Stan_fns_basic_NicoSIMD_AVX.stan"
        #include "Stan_fns_corr_NicoSIMD_AVX.stan"
        #include "Stan_fns_ordinal_NicoSIMD_AVX.stan"
        #include "Stan_fns_log_lik_NicoSIMD_AVX.stan"
        #include "Stan_fns_simplex_NicoSIMD_AVX.stan"
        #include "Stan_fns_Jacobian_NicoSIMD_AVX.stan"
        #include "Stan_fns_model_fit_NicoSIMD_AVX.stan"
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
        //// ---- Prior on baseline Sp and Se at chosen cutoffs (swapped for the Dirichlet prior of dirichlet_phi):
        ////
        int<lower=0, upper=1> use_Se_Sp_prior;     //// 0 = original dirichlet_phi prior
        int<lower=1, upper=4> Se_Sp_prior_type;    //// 1 = Dirichlet, 2-3 = Betas, 4 = link-scale Normal
        int<lower=0, upper=n_thr> n_chosen_cutoffs_max;
        array[2] int<lower=0, upper=n_chosen_cutoffs_max> n_chosen_cutoffs;
        array[2, n_chosen_cutoffs_max] int<lower=1, upper=n_thr> chosen_cutoffs;
        array[2] vector<lower=0.0>[n_chosen_cutoffs_max + 1] prior_Se_Sp_dirichlet_A;
        array[2] matrix<lower=0.0>[n_chosen_cutoffs_max, 2] prior_Se_Sp_beta;
        array[2] matrix[n_chosen_cutoffs_max, 2] prior_Se_Sp_link_mean_SD;  //// Mean and SD of link(Sp), link(Se)
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
        if (use_Se_Sp_prior == 1) reject("this version of the model has no Se / Sp prior");
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
        //// ---- Sums of prior_dirichlet_alpha over the segments between chosen cutoffs:
        ////
        array[2] vector[n_chosen_cutoffs_max + 1] alpha_of_the_segments;
        for (c in 1:2) {
              alpha_of_the_segments[c] = rep_vector(0.0, n_chosen_cutoffs_max + 1);
              int k_previous = 0;
              for (j in 1:n_chosen_cutoffs[c]) {
                   int k = chosen_cutoffs[c, j];
                   if (k <= k_previous) reject("chosen_cutoffs must rise: c = ", c);
                   if (use_Se_Sp_prior == 1 && Se_Sp_prior_type == 4 &&
                       prior_Se_Sp_link_mean_SD[c][j, 2] <= 0.0) {
                        reject("The link-scale prior SD must be positive: c = ", c, ", j = ", j);
                   }
                   alpha_of_the_segments[c][j] = sum(prior_dirichlet_alpha[(k_previous + 1):k]);
                   k_previous = k;
              }
              alpha_of_the_segments[c][n_chosen_cutoffs[c] + 1] =
                   sum(prior_dirichlet_alpha[(k_previous + 1):n_cat]);
        }
        ////
        //// ---- Class on the cutpoints (p_class) and class on the baseline location / raw scale (q_class):
        ////
        int p_class = (n_chosen_cutoffs[2] > n_chosen_cutoffs[1]) ? 2 : 1;
        int q_class = 3 - p_class;
        int n_q = (use_Se_Sp_prior == 1) ? n_chosen_cutoffs[q_class] : 0;
        int n_q_exact = min(n_q, 2);
        for (j in 1:n_q_exact) {
             int found = 0;
             for (i in 1:n_chosen_cutoffs[p_class]) {
                  if (chosen_cutoffs[p_class, i] == chosen_cutoffs[q_class, j]) found = 1;
             }
             if (found == 0) reject("chosen_cutoffs of c = ", q_class, " must be among those of c = ", p_class);
        }
        if (n_q > 0 && baseline_case[1] == 0.0) reject("baseline_case[1] must not be 0");
        int first_beta_mu = (n_q >= 1) ? 2 : 1;
        int first_raw_scale_mu = (n_q >= 2) ? 2 : 1;

}


parameters {
        vector[n_studies] beta_z;  //// the nuisance block, declared first
        vector[n_covariates] beta_mu;  
        real<lower=0.0> beta_SD;   
        // vector[n_studies] beta_z;
        ////
        vector[n_covariates] raw_scale_mu;    
        real<lower=0.0> raw_scale_SD;   
        vector[n_studies] raw_scale_z;
        ////
        matrix<lower=-15, upper=2.5>[n_studies, n_thr] C_raw; //// restricting solves divs/numerical stab. issues w/ lots of sparse categories!
        real<lower=1.0> kappa;
        simplex[n_cat] dirichlet_phi; // props/probs
}


transformed parameters { 
        real log_kappa = nico_log(kappa);
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

        //// (as construct_C, for all studies at once)
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
            Xscale[, 1] =  ((softplus == 1) ? 1.4426950408889634 * nico_log1p_exp( -mult_nd*Xraw_scale_baseline) :
                                                nico_exp( -mult_nd*Xraw_scale_baseline)); 
            // Xscale[, 2] =  ((softplus == 1) ? softplus_scaled( -mult_d*Xraw_scale_baseline)  : exp( -mult_d*Xraw_scale_baseline));
            Xscale[, 2] =  ((softplus == 1) ? 1.4426950408889634 * nico_log1p_exp( -mult_d*Xraw_scale_baseline)  :
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
        if (first_beta_mu <= n_covariates) {
             beta_mu[first_beta_mu:n_covariates] ~ nico_normal( prior_beta_mu_mean[first_beta_mu:n_covariates],
                                                           prior_beta_mu_SD[first_beta_mu:n_covariates]);
        }
        beta_SD ~ nico_normal( prior_beta_SD_mean, prior_beta_SD_SD);
        ////
        // raw_scale_mu[1:n_covariates] ~ normal(
        //                                prior_raw_scale_mu_mean[1:n_covariates], 
        //                                prior_raw_scale_mu_SD[1:n_covariates]);
        if (first_raw_scale_mu <= n_covariates) {
             raw_scale_mu[first_raw_scale_mu:n_covariates] ~ nico_normal( prior_raw_scale_mu_mean[first_raw_scale_mu:n_covariates],
                                                                     prior_raw_scale_mu_SD[first_raw_scale_mu:n_covariates]);
        }
        raw_scale_SD ~ nico_normal( prior_raw_scale_SD_mean, prior_raw_scale_SD_SD);
        ////
        {
           log_kappa ~ nico_student_t(prior_kappa_dof, prior_kappa_mean, prior_kappa_SD);
           target += -log_kappa;
           ////
           dirichlet_phi ~ nico_dirichlet(prior_dirichlet_alpha);
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

        //// ---- Induced-dirichlet between study ** model ** (NOT a prior: part of the likelihood since
        ////      random-effect cutpoints), every study's cutpoints in one node, as
        ////      induced_dirichlet_given_C on each row:
        //// (anchored at the prior mean of the mid-location)
        for (s in 1:n_studies) {
             target += induced_dirichlet_given_C_on_the_log_scale_lupdf(C[s, ] -
                                                                        0.5 * (mult_nd + mult_d) *
                                                                        prior_beta_mu_mean[1] |
                                                                        alpha, use_probit_link);
        }
        ////
        //// ---- Likelihood / Model:
        ////
        target += nico_normal_lupdf(to_vector(beta_z) | 0, 1); // part of between-study model, NOT prior
        target += nico_normal_lupdf(to_vector(raw_scale_z) | 0, 1); // part of between-study model, NOT prior
        ////
        target += raw_C_to_C_log_det_J_lp(C_raw, softplus);
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
            //// ---- Get the cutpoint index (k) to map "latent_surv[c][s, cut_i]" to correct cutpoint "C[k]":
            ////
            array[2] matrix[n_studies, n_thr] latent_surv = map_latent_surv_prob_to_random_C_HSROC(
                                                            C, n_thr, Xlocations, Xscale, n_studies,
                                                            n_obs_cutpoints, cutpoint_index);
            target += compute_log_lik_binomial_fact_lp(
                      latent_surv, use_probit_link, n_thr, n_studies,
                      x_2, n, N_total, n_obs_cutpoints,
                      K_fold_CV_indicator);
        }
}


generated quantities {
    real simd_lanes = bmvp_simd_lanes();
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
      // for (k in 1:n_thr) {
      //   C_mu[k] = median(C[, k]);
      // }
      {
         //// (the shapes of the Beta of each cutpoint from alpha itself: sum(alpha) = kappa + 0.01 * n_cat; the
         ////  study cutpoints are anchored at 0.5 * (mult_nd + mult_d) * prior_beta_mu_mean[1] = 0)
         real alpha_0 = sum(alpha);
         vector[n_thr] P = cumulative_sum(alpha[1:n_thr]) / alpha_0;
         vector[n_thr] P_upper = reverse(cumulative_sum(reverse(alpha[2:n_cat]))) / alpha_0;
         for (k in 1:n_thr) {
            C_mu[k] = dirichlet_cutpoint_median(P[k], P_upper[k], alpha_0, use_probit_link);
         }
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
                               1.4426950408889634 * nico_log1p_exp(-mult_nd * Xraw_scale_baseline) :
                                                  nico_exp(-mult_nd * Xraw_scale_baseline));
      real scale_d_baseline  = ((softplus == 1) ? 
                               1.4426950408889634 * nico_log1p_exp(-mult_d  * Xraw_scale_baseline) :
                                                  nico_exp(-mult_d  * Xraw_scale_baseline));
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
               C_pred[k] = C_mu[k] + nico_normal_rng(0, sigma_C_k);
            }
        }
      }
      ////
      real beta_random_pred      = nico_normal_rng(0.0, beta_SD);
      real raw_scale_random_pred = nico_normal_rng(0.0, raw_scale_SD);
      ////
      real Xbeta_baseline_pred = dot_product(baseline_case, beta_mu) + beta_random_pred;
      real Xbeta_baseline_pred_nd = mult_nd * Xbeta_baseline_pred;
      real Xbeta_baseline_pred_d  = mult_d  * Xbeta_baseline_pred;
      ////
      real Xraw_scale_baseline_pred = dot_product(baseline_case, raw_scale_mu) + raw_scale_random_pred;
      real scale_baseline_pred_nd = ((softplus == 1) ? 
                                    1.4426950408889634 * nico_log1p_exp(-mult_nd * Xraw_scale_baseline_pred) :
                                                       nico_exp(-mult_nd * Xraw_scale_baseline_pred));
      real scale_baseline_pred_d  = ((softplus == 1) ? 
                                    1.4426950408889634 * nico_log1p_exp(-mult_d  * Xraw_scale_baseline_pred) :
                                                       nico_exp(-mult_d  * Xraw_scale_baseline_pred));
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
          array[2] matrix[n_studies, n_thr] latent_surv = map_latent_surv_prob_to_random_C_HSROC(
                                                          C, n_thr, Xlocations, Xscale, n_studies, 
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
