//// DTA_MA_ord_biv_FIXEDthr_Se_Sp_prior_NicoStan_fused_AVX.stan
//// The model of MetaOrdDTA's DTA_MA_ord_biv_FIXEDthr_Se_Sp_prior.stan with NicoSIMD's AVX functions:
//// the likelihood and the induced Dirichlet of the cutpoints each in one autodiff node.
////



functions {
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
        ////
        //// ---- NicoSIMD's AVX functions (their bodies in NicoSIMD's user header):
        ////
        real nico_binomial_counts_above_successive_cutpoints_of_a_cumulative_link_lpmf(
              array[] int n, array[] int n_previous, vector latent, data vector is_first, data int link,
              data real bound);
        real nico_induced_dirichlet_with_Se_Sp_prior_lpdf(vector C, vector alpha, data int link,
                                                          data real location);
        real nico_induced_dirichlet_with_Se_Sp_prior_lpdf(vector C, vector alpha, data int link,
                                                          data real location, vector alpha_of_the_segments,
                                                          data array[] int chosen_cutoffs,
                                                          data int Se_Sp_prior_type,
                                                          vector prior_Se_Sp_dirichlet_A, matrix prior_Se_Sp_beta,
                                                          matrix prior_Se_Sp_link_mean_SD);
        real nico_log_Set0_marginal_at_chosen_cutoffs(vector x_chosen, vector alpha_of_the_segments,
                                                      real v0, data int use_probit_link);
        real nico_Se_Sp_prior_at_chosen_cutoffs_lpdf(vector x, data int Se_Sp_prior_type,
                                                     vector prior_Se_Sp_dirichlet_A, matrix prior_Se_Sp_beta,
                                                     matrix prior_Se_Sp_link_mean_SD, data int link);
        vector nico_exp(vector x);
        vector nico_log1p_exp(vector x);
        vector nico_log_inv_logit(vector x);
        vector nico_log1m_inv_logit(vector x);
        real nico_log_det_jacobian_of_cutpoints_from_raw_differences(vector C_raw, data int softplus);
        real nico_log_det_jacobian_of_cutpoints_from_raw_differences(matrix C_raw, data int softplus);
        real nico_raw_C_to_C_log_det_J_lp(vector C_raw, data int softplus);
        real nico_raw_C_to_C_log_det_J_lp(matrix C_raw, data int softplus);
        array[] matrix nico_map_latent_surv_to_fixed_hetero_C(data array[] vector C,
                                                              data matrix locations, data matrix scales,
                                                              data int n_studies, data array[] int n_obs_cutpoints,
                                                              data array[,,] int cutpoint_index);
        array[] matrix nico_map_latent_surv_to_fixed_hetero_C(data array[] vector C,
                                                              data matrix locations,
                                                              data int n_studies, data array[] int n_obs_cutpoints,
                                                              data array[,,] int cutpoint_index);
        array[,] matrix nico_compute_log_lik_binomial_fact_data(data array[] matrix latent_surv,
                                                                data int use_probit_link, data int n_thr,
                                                                data int n_studies, data array[,,] int x_2,
                                                                data array[,,] int n, data array[,] int N_total,
                                                                data array[] int n_obs_cutpoints,
                                                                data array[] int K_fold_CV_indicator);
        real nico_induced_dirichlet_prior_of_ordered_cutpoints_lpdf( vector C, 
                                                                     vector alpha,
                                                                     data int link,
                                                                     data real location);
        vector nico_Phi(vector x);
        real nico_lkj_corr_cholesky_lpdf(matrix L, data real eta);
        real nico_dirichlet_lpdf(vector theta, vector alpha);
        real nico_beta_lpdf(vector y, vector alpha, vector beta);
        vector nico_log1m(vector x);
        vector nico_inv_logit(vector x);
        real nico_normal_lpdf(vector y, real mu, real sigma);
        real nico_normal_lpdf(vector y, vector mu, vector sigma);
        //// (the nico_ functions of the NicoSIMD versions of the included files:)
        real nico_std_normal_log_qf(real x);
        vector nico_std_normal_log_qf(vector x);
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
        vector nico_log(vector x);
        vector nico_log1m_exp(vector x);
        vector nico_inv_Phi(vector x);
        vector nico_tanh(vector x);
        vector nico_logit(vector x);
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
        real nico_normal_lpdf(vector y, vector mu, real sigma);
        real nico_normal_lpdf(vector y, real mu, vector sigma);
        real nico_normal_lpdf(row_vector y, row_vector mu, row_vector sigma);
        real nico_normal_lpdf(array[] real y, vector mu, vector sigma);
        real nico_normal_lpdf(array[] real y, vector mu, real sigma);
        real nico_normal_lpdf(array[] real y, real mu, vector sigma);
        real nico_normal_lpdf(array[] real y, real mu, real sigma);
        real nico_normal_lpdf(real y, real mu, real sigma);
        real nico_lgamma(data real x);
        vector nico_lgamma(vector x);
        real nico_beta_lpdf(real y, data real alpha, data real beta);
        real nico_log_diff_exp(real x, real y);
        vector nico_log_diff_exp(vector x, vector y);
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
        //// ---- Priors directly on Sp (c = 1) and Se (c = 2) at chosen cutoffs:
        //// ---- Types 1 and 3 use P = Sp and 1 - Se; types 2 and 4 use Sp and Se:
        ////
        //// ---- alpha = (1, ..., 1) through this prior (the same density, up to a constant): use_Se_Sp_prior = 1,
        //// ---- Se_Sp_prior_type = 4, every cutoff chosen in both classes (n_chosen_cutoffs = (n_thr, n_thr)), and
        //// ---- prior_Se_Sp_link_mean_SD = (0, 1) at every cutoff, i.e. C_k ~ normal(baseline_location_prior, 1) for every k, kept
        //// ---- in rising order (baseline_location_prior = prior_beta_mu_mean[c, 1]; a single chosen cutoff is NOT alpha = 1):
        ////
        int<lower=0, upper=2> use_Se_Sp_prior;     //// 0 = original cutpoint prior; 1 = Se/Sp prior; 2 = Set 0 via Se/Sp prior
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
        int<lower=0, upper=1> softplus; 
        int<lower=0, upper=1> log_lik_indicator; //// 1 = compute log_lik_study; 0 = skip it
        real C_raw_lower;
        real C_raw_upper;
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
                   alpha_of_the_segments[c][j] = sum(prior_dirichlet_alpha[c][(k_previous + 1):k]);
                   k_previous = k;
              }
              alpha_of_the_segments[c][n_chosen_cutoffs[c] + 1] =
                   sum(prior_dirichlet_alpha[c][(k_previous + 1):n_cat]);
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

        if (use_Se_Sp_prior == 1) {
              if (use_probit_link != 1) reject("use_Se_Sp_prior == 1 requires use_probit_link == 1");
              for (c in 1:2) {
                    for (j in 1:n_cat) {
                          if (abs(prior_dirichlet_alpha[c][j] - 1.0) > 1e-12) {
                                reject("use_Se_Sp_prior == 1 requires prior_dirichlet_alpha == 1");
                          }
                    }
              }
        }
        vector[2] baseline_location_prior;  //// the baseline case at the prior means (the same under every prior set)
        baseline_location_prior[1] = dot_product(baseline_case_nd, to_vector(prior_beta_mu_mean[1, 1:n_covariates_nd]));
        baseline_location_prior[2] = dot_product(baseline_case_d,  to_vector(prior_beta_mu_mean[2, 1:n_covariates_d]));
        vector[2] v0_baseline_prior;
        v0_baseline_prior[1] = dot_self(baseline_case_nd .* to_vector(prior_beta_mu_SD[1, 1:n_covariates_nd]));
        v0_baseline_prior[2] = dot_self(baseline_case_d  .* to_vector(prior_beta_mu_SD[2, 1:n_covariates_d]));
}


parameters { 
        matrix[2, n_covariates_max] beta_mu;       
        row_vector<lower=0.0>[2] beta_SD;    
        matrix[n_studies, 2] beta_z;    //// Study-specific random-effects (raw params)
        //// 
        real<lower=beta_corr_lb, upper=beta_corr_ub> beta_corr;  //// between-study corr (possibly restricted)  
        //// 
        array[2] vector<lower=C_raw_lower, upper=C_raw_upper>[n_thr] C_raw_vec; //// Global raw-cutpoints (unconstrained - nuisance block #2) 
}  


transformed parameters { 
        ////   
        //// ---- Construct (global) cutpoints:
        ////
        array[2] vector[n_thr] C; 
        for (c in 1:2) { 
            // C[c] = construct_C(C_raw_vec[c], softplus);
            //// (as construct_C: the first cutpoint, then the cumulative sums of the exp() (or softplus) of
            ////  the raw differences, with NicoSIMD's nico_exp / nico_log1p_exp)
            C[c][1] = C_raw_vec[c][1];
            if (n_thr > 1) {
                  vector[n_thr - 1] increments_of_the_cutpoints = (softplus == 1) ? 
                                                                  nico_log1p_exp(C_raw_vec[c][2:n_thr]) :
                                                                  nico_exp(C_raw_vec[c][2:n_thr]);
                  for (k in 2:n_thr) {
                        C[c][k] = C[c][k - 1] + increments_of_the_cutpoints[k - 1];
                  }
            }
        }
        ////
        //// ---- Construct simple 2x2 (bivariate) between-study corr matrices for between-study model:
        ////  
        cholesky_factor_corr[2] beta_L_Omega = make_bivariate_L_Omega(beta_corr);
        cholesky_factor_cov[2]  beta_L_Sigma = diag_pre_multiply(beta_SD, beta_L_Omega);
        ////
        //// ---- Study-level random effects (after Cholesky decomposition):
        ////
        matrix[n_studies, 2] beta_random;
        for (s in 1:n_studies) {
             beta_random[s, ] =  beta_z[s, ] * beta_L_Sigma;
        }
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
        beta_mu[1, 1:n_covariates_nd] ~ nico_normal(
                                          prior_beta_mu_mean[1, 1:n_covariates_nd],
                                          prior_beta_mu_SD[1, 1:n_covariates_nd]);
        beta_mu[2, 1:n_covariates_d]  ~ nico_normal(
                                          prior_beta_mu_mean[2, 1:n_covariates_d],
                                          prior_beta_mu_SD[2, 1:n_covariates_d]);
        ////
        to_vector(beta_SD) ~ nico_normal(prior_beta_SD_mean, prior_beta_SD_SD);
        beta_L_Omega ~ nico_lkj_corr_cholesky(prior_beta_corr_LKJ); 
        //// 
        //// ---- Induced-dirichlet ** Prior ** model:
        //// 
        // for (c in 1:2) {
        //      vector[n_cat] Ind_Dir_ord_prob; //// = cumul_probs_to_ord_probs(Phi(C));
        //      if (use_probit_link == 1) Ind_Dir_ord_prob = cumul_probs_to_ord_probs(Phi(C[c]));
        //      else                      Ind_Dir_ord_prob = cumul_probs_to_ord_probs(inv_logit(C[c]));
        //      ////
        //      Ind_Dir_ord_prob ~ induced_dirichlet_given_C(
        //                         C[c], prior_dirichlet_alpha[c], use_probit_link); // more efficient than Betancourt et al. and seems to work fine.
        // }
        for (c in 1:2) {
             if (use_Se_Sp_prior >= 1 && n_chosen_cutoffs[c] > 0) {
             
                    matrix[n_chosen_cutoffs[c], 2] prior_beta_of_this_class = prior_Se_Sp_beta[c][1:n_chosen_cutoffs[c], ];
                    matrix[n_chosen_cutoffs[c], 2] prior_link_mean_SD_of_this_class = prior_Se_Sp_link_mean_SD[c][1:n_chosen_cutoffs[c], ];
                    ////
                    if (c == 2) {
                         prior_link_mean_SD_of_this_class[, 1] = - prior_link_mean_SD_of_this_class[, 1];
                         if (Se_Sp_prior_type == 2) {
                              prior_beta_of_this_class = append_col( prior_beta_of_this_class[, 2],
                                                                     prior_beta_of_this_class[, 1]);
                         }
                    }
                    if (use_Se_Sp_prior == 2) {
                         vector[n_thr] C_mu_minus_baseline_location_prior = C[c] - baseline_location_prior[c];
                         target += nico_induced_dirichlet_with_Se_Sp_prior_lupdf(
                                   C_mu_minus_baseline_location_prior | prior_dirichlet_alpha[c],
                                   use_probit_link, 
                                   0.0,
                                   alpha_of_the_segments[c][1:(n_chosen_cutoffs[c] + 1)],
                                   chosen_cutoffs[c, 1:n_chosen_cutoffs[c]], 
                                   Se_Sp_prior_type,
                                   prior_Se_Sp_dirichlet_A[c][1:(n_chosen_cutoffs[c] + 1)],
                                   prior_beta_of_this_class,
                                   prior_link_mean_SD_of_this_class);
                    } else {
                    real Xbeta_baseline = (c == 1) ?
                                          dot_product(baseline_case_nd, to_vector(beta_mu[1, 1:n_covariates_nd])) :
                                          dot_product(baseline_case_d, to_vector(beta_mu[2, 1:n_covariates_d]));
                    ////
                    //// ---- C_mu - Xbeta_baseline of class c: the quantity that Se and Sp are made of (generated
                    ////      quantities: Phi(-(C - Xbeta_baseline))). The cutpoints are fixed, so C is itself the summary
                    ////      cutpoint C_mu. The random-cutpoint model (DTA_MA_ord_biv_RANDthr_kappa_popmed_Se_Sp_prior.stan)
                    ////      has this same prior on its own C_mu - Xbeta_baseline:
                         vector[n_thr] C_mu_minus_baseline_location_prior = C[c] - baseline_location_prior[c];
                         vector[n_chosen_cutoffs[c]] x_chosen = C[c][chosen_cutoffs[c, 1:n_chosen_cutoffs[c]]] - Xbeta_baseline;
                         target += nico_induced_dirichlet_prior_of_ordered_cutpoints_lupdf(
                                   C_mu_minus_baseline_location_prior | prior_dirichlet_alpha[c], use_probit_link, 0.0);
                         target += - nico_log_Set0_marginal_at_chosen_cutoffs(
                                   x_chosen, alpha_of_the_segments[c][1:(n_chosen_cutoffs[c] + 1)],
                                   v0_baseline_prior[c], use_probit_link);
                         target += nico_Se_Sp_prior_at_chosen_cutoffs_lupdf(
                                   x_chosen | Se_Sp_prior_type,
                                   prior_Se_Sp_dirichlet_A[c][1:(n_chosen_cutoffs[c] + 1)],
                                   prior_beta_of_this_class,
                                   prior_link_mean_SD_of_this_class,
                                   use_probit_link);
                    }
             } else {
                    ////
                    //// (the cutpoints measured from baseline_location_prior, the baseline case at the prior means)
                    vector[n_thr] C_mu_minus_baseline_location_prior = C[c] - baseline_location_prior[c];
                    target += nico_induced_dirichlet_prior_of_ordered_cutpoints_lupdf(
                              C_mu_minus_baseline_location_prior | prior_dirichlet_alpha[c], use_probit_link, 0.0);
             }
        }
        //// 
        //// ---- Likelihood / Model:
        ////
        target += nico_normal_lupdf(to_vector(to_matrix(beta_z)) | 0.0, 1.0);
        ////
        for (c in 1:2) {
           target += nico_raw_C_to_C_log_det_J_lp(C_raw_vec[c], softplus);
        }
        ////
        //// ---- Log-likelihood:
        ////
        {
            ////
            //// ---- every term in one autodiff node: nico_binomial_counts_above_successive_cutpoints_of_a_cumulative_link_lpmf
            ////
            vector[2 * n_thr] cutpoints_of_both_classes = append_row(C[1], C[2]);
            vector[2 * n_studies] locations_of_both_classes = to_vector(Xbeta);
            //// (as map_latent_surv_to_fixed_hetero_C: location - C)
            vector[n_terms] latent_of_each_term = locations_of_both_classes[index_of_the_location_of_each_term] -
                                                  cutpoints_of_both_classes[index_of_the_cutpoint_of_each_term];
            ////
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
                                  nico_Phi(-(C[1]  - Xbeta_baseline_nd))  : 
                                  nico_inv_logit(-(C[1] - Xbeta_baseline_nd));
      vector[n_thr] Sp_baseline = 1.0 - Fp_baseline;
      vector[n_thr] Se_baseline = (use_probit_link == 1) ? 
                                  nico_Phi(-(C[2] - Xbeta_baseline_d))   : 
                                  nico_inv_logit(-(C[2] - Xbeta_baseline_d));
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
                                       nico_Phi(-(C[1]  - Xbeta_baseline_pred_nd)) : 
                                       nico_inv_logit(-(C[1] - Xbeta_baseline_pred_nd));
      vector[n_thr] Sp_baseline_pred = 1.0 - Fp_baseline_pred;
      vector[n_thr] Se_baseline_pred = (use_probit_link == 1) ? 
                                       nico_Phi(-(C[2]  - Xbeta_baseline_pred_d))  : 
                                       nico_inv_logit(-(C[2] - Xbeta_baseline_pred_d));
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
      if (log_lik_indicator == 1) {
          // real scale = 1.0; // since using "ord_bivariate-like"" param.
          array[2] matrix[n_studies, n_thr] latent_surv = nico_map_latent_surv_to_fixed_hetero_C(
                                                          C, Xbeta, n_studies, 
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














