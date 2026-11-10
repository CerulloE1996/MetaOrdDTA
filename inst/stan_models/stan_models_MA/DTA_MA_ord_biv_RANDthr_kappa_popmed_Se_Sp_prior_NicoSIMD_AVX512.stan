// NicoSIMD_AVX512 version of DTA_MA_ord_biv_RANDthr_kappa_popmed_Se_Sp_prior.stan written by NicoStan (source model md5 56ffdd058ec61ce9de06dbe22d01db55; NicoSIMD declarations md5 7f272e2ac00d6303123b13fb8f59653e; writer implementation md5 f6843a34264d79c6ed96effa2d04b0c4)
// NicoSIMD functions used: exp, log, log1m, log1p_exp, log1m_exp, inv_logit, log_inv_logit, Phi, inv_Phi, tanh, logit, log1m_inv_logit, expm1, normal_rng, normal, lgamma, beta, dirichlet, lkj_corr_cholesky, log_diff_exp, std_normal_log_qf, student_t
//// DTA_MA_ord_biv_RANDthr_kappa_popmed_Se_Sp_prior_Stan_math_only.stan
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
    real nico_expm1(real x);
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
    vector nico_expm1(vector x);
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
    row_vector nico_expm1(row_vector x);
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
    matrix nico_expm1(matrix x);
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
    real nico_lkj_corr_cholesky_lpdf(matrix L, data real eta);
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
        vector[n_covariates_nd] baseline_case_nd; // must be user-inputted - e.g. could be {0, 1, 45.3} for 2 binary covariates and 1 cts one (e.g. age)
        vector[n_covariates_d] baseline_case_d;   // must be user-inputted - e.g. could be {0, 1, 45.3} for 2 binary covariates and 1 cts one (e.g. age)
        ////
        //// ---- Priors for locations ("beta"):
        ////
        matrix[2, n_covariates_max] prior_beta_mu_mean;
        matrix[2, n_covariates_max] prior_beta_mu_SD;
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
        array[2] vector[n_thr + 1] prior_dirichlet_alpha;
        vector[2] prior_kappa_mean;
        vector<lower=0.0>[2] prior_kappa_SD;
        vector<lower=0.0>[2] prior_kappa_dof;
        ////
        //// ---- Priors directly on Sp (c = 1) and Se (c = 2) at chosen cutoffs:
        //// ---- Types 1 and 3 use P = Sp and 1 - Se; types 2 and 4 use Sp and Se:
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
        // if (use_Se_Sp_prior == 1) reject("this version of the model has no Se / Sp prior");
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
        //// ---- The reported cutpoints of each study, study after study (C_raw and C hold only these):
        ////
        array[n_studies] int n_reported_cutpoints;
        for (s in 1:n_studies) {
              n_reported_cutpoints[s] = 0;
              for (cut_i in 1:n_obs_cutpoints[s]) {
                   if (cutpoint_index[1, s, cut_i + 1] != cutpoint_index[2, s, cut_i + 1]) {
                        reject("cutpoint_index must be the same for both classes: s = ", s);
                   }
                   if (cutpoint_index[1, s, cut_i + 1] != -1) n_reported_cutpoints[s] += 1;
              }
        }
        int n_C = sum(n_reported_cutpoints);  //// the reported cutpoints of all studies
        int n_studies_with_a_cutpoint = 0;
        for (s in 1:n_studies) {
              if (n_reported_cutpoints[s] > 0) n_studies_with_a_cutpoint += 1;
        }
        array[n_studies] int first_C_of_each_study;
        array[n_studies, n_thr] int position_of_each_C;  //// 0 when cut_i is not reported
        array[n_C - n_studies_with_a_cutpoint] int position_of_each_increment;
        array[2] vector[n_C] C_raw_lower_of_each_C;
        array[2] vector[n_C] C_raw_upper_of_each_C;
        array[n_C + n_studies] int upper_index_of_each_segment;  //// a segment = the categories k_(j - 1) + 1 to k_j
        array[n_C + n_studies] int lower_index_of_each_segment;
        {
              //// (each one-cutpoint step lies between these, as each column of the full C_raw did; a gap of g
              ////  cutpoints is g such steps)
              real step_lower = (softplus == 1) ? nico_log1p_exp(C_raw_lower) : nico_exp(C_raw_lower);
              real step_upper = (softplus == 1) ? nico_log1p_exp(C_raw_upper) : nico_exp(C_raw_upper);
              int p = 0;
              int i = 0;
              for (s in 1:n_studies) {
                    first_C_of_each_study[s] = p + 1;
                    int k_previous = 0;
                    for (cut_i in 1:n_thr) position_of_each_C[s, cut_i] = 0;
                    for (cut_i in 1:n_obs_cutpoints[s]) {
                         int k = cutpoint_index[1, s, cut_i + 1];
                         if (k != -1) {
                              if (k <= k_previous) reject("cutpoint_index must rise within each study: s = ", s);
                              p += 1;
                              position_of_each_C[s, cut_i] = p;
                              upper_index_of_each_segment[p + s - 1] = k + 1;
                              lower_index_of_each_segment[p + s - 1] = k_previous + 1;
                              if (k_previous == 0) {
                                   C_raw_lower_of_each_C[1][p] = C_raw_lower + (k - 1) * step_lower;
                                   C_raw_upper_of_each_C[1][p] = C_raw_upper + (k - 1) * step_upper;
                              } else {
                                   i += 1;
                                   position_of_each_increment[i] = p;
                                   C_raw_lower_of_each_C[1][p] = (softplus == 1) ? nico_log(nico_expm1((k - k_previous) * step_lower)) :
                                                                                 C_raw_lower + nico_log(k - k_previous);
                                   C_raw_upper_of_each_C[1][p] = (softplus == 1) ? nico_log(nico_expm1((k - k_previous) * step_upper)) :
                                                                                 C_raw_upper + nico_log(k - k_previous);
                              }
                              k_previous = k;
                         }
                    }
                    upper_index_of_each_segment[p + s] = n_cat + 1;
                    lower_index_of_each_segment[p + s] = k_previous + 1;
              }
              C_raw_lower_of_each_C[2] = C_raw_lower_of_each_C[1];
              C_raw_upper_of_each_C[2] = C_raw_upper_of_each_C[1];
        }

}
  

parameters {
        array[n_studies] row_vector[2] beta_z;  //// the nuisance block, declared first
        matrix[2, n_covariates_max] beta_mu;      
        row_vector<lower=0.0>[2] beta_SD;   
        // matrix[n_studies, 2] beta_z;    //// Study-specific random effects (off-centered parameterisation)
        ////
        real<lower=beta_corr_lb, upper=beta_corr_ub> beta_corr;  //// between-study corr (possibly restricted)
        ////
        array[2] vector<lower=C_raw_lower_of_each_C, upper=C_raw_upper_of_each_C>[n_C] C_raw; //// the reported cutpoints of each study only (restricting solves divs/numerical stab. issues w/ lots of sparse categories!)
        vector<lower=1.0>[2] kappa;
        array[2] simplex[n_cat] dirichlet_phi;
        // array[2] simplex[n_thr] raw_dirichlet_phi;
}

 
transformed parameters { 
        // array[2] simplex[n_cat] dirichlet_phi;
        // for (c in 1:2) {
        //   dirichlet_phi[c] = stickbreaking_logistic_simplex_constrain(raw_dirichlet_phi[c]); // THIS IS WRONG!!!!!!
        // }
        ////
        vector[2] log_kappa = nico_log(kappa);
        ////
        array[2] vector<lower=0.0>[n_cat] alpha;
        for (c in 1:2) {
          alpha[c] = 0.01 + dirichlet_phi[c] * kappa[c];
        }
        ////
        array[2] vector<lower=0.0>[n_cat] dirichlet_cat_SDs_sigma;
        for (c in 1:2) {
           real alpha_0 = sum(alpha[c]);
           dirichlet_cat_SDs_sigma[c] = sqrt((alpha[c] .* (alpha_0 - alpha[c]) ) / 
                                       (square(alpha_0) * (alpha_0 + 1.0)));
        }
        ////
        //// ---- Construct study-specific cutpoints:
        //// 
        // array[2] matrix[n_studies, n_thr] C;
        // for (c in 1:2) {
        //   for (s in 1:n_studies) {
        //        C[c][s, ] =  construct_C(C_raw[c][s, 1:n_thr], softplus);
        //   }
        // }

        //// (the reported cutpoints of each study only, study after study, as construct_C on each study's own)
        array[2] vector[n_C] C;
        for (c in 1:2) {
              vector[n_C - n_studies_with_a_cutpoint] increments_of_the_cutpoints =
                    (softplus == 1) ? nico_log1p_exp(C_raw[c][position_of_each_increment]) :
                                      nico_exp(C_raw[c][position_of_each_increment]);
              C[c] = C_raw[c];
              for (i in 1:(n_C - n_studies_with_a_cutpoint)) {
                    int p = position_of_each_increment[i];
                    C[c][p] = C[c][p - 1] + increments_of_the_cutpoints[i];
              }
        }
        ////
        //// ---- Under the Se / Sp prior: cutpoints relative to the intercept of the baseline case (the Se / Sp
        ////      prior and the likelihood hold C - Xbeta only, so the raw cutpoints are then free of beta_mu):
        ////
        if (use_Se_Sp_prior == 1) {
             C[1] += dot_product(baseline_case_nd, to_vector(beta_mu[1, 1:n_covariates_nd]));
             C[2] += dot_product(baseline_case_d,  to_vector(beta_mu[2, 1:n_covariates_d]));
        }
        ////
        //// ---- Construct simple 2x2 (bivariate) between-study corr matrices:
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
        // matrix[n_studies, 2] beta_random = beta_z * beta_L_Sigma;
        matrix[n_studies, 2] beta_random = to_matrix(beta_z) * beta_L_Sigma;
        ////
        //// ---- Linear predictors for each disease statu + apply covariates for non-diseased and diseased groups:
        ////
        matrix[n_studies, 2] Xbeta;
        Xbeta[, 1]  = X_nd[1:n_studies, 1:n_covariates_nd] * to_vector(beta_mu[1, 1:n_covariates_nd]) + beta_random[, 1];
        Xbeta[, 2]  = X_d[1:n_studies, 1:n_covariates_d]   * to_vector(beta_mu[2, 1:n_covariates_d])  + beta_random[, 2];
}


model { 
        // for (c in 1:2) {
        //   target += stickbreaking_logistic_simplex_constrain_lp(raw_dirichlet_phi[c]);
        // }
        ////
        //// ---- Priors: 
        ////
        beta_mu[1, 1:n_covariates_nd] ~ nico_normal(
                                          prior_beta_mu_mean[1, 1:n_covariates_nd],
                                          prior_beta_mu_SD[1, 1:n_covariates_nd]);
        beta_mu[2, 1:n_covariates_d]  ~ nico_normal(
                                          prior_beta_mu_mean[2, 1:n_covariates_d],
                                          prior_beta_mu_SD[2, 1:n_covariates_d]);
        // ////
        // beta_SD ~ normal(prior_beta_SD_mean, prior_beta_SD_SD);

        // to_vector(beta_mu[1, 1:n_covariates_nd]) ~ normal(
        //                                             to_vector(prior_beta_mu_mean[1, 1:n_covariates_nd]),
        //                                             to_vector(prior_beta_mu_SD[1, 1:n_covariates_nd]));
        // to_vector(beta_mu[2, 1:n_covariates_d])  ~ normal(
        //                                             to_vector(prior_beta_mu_mean[2, 1:n_covariates_d]),
        //                                             to_vector(prior_beta_mu_SD[2, 1:n_covariates_d]));
        ////
        to_vector(beta_SD) ~ nico_normal(prior_beta_SD_mean, prior_beta_SD_SD);
        beta_L_Omega ~ nico_lkj_corr_cholesky(prior_beta_corr_LKJ);
        ////
        vector[2] baseline_location;  //// the prior mean of the intercept (prior intercepts as anchor point)
        baseline_location[1] = prior_beta_mu_mean[1, 1];
        baseline_location[2] = prior_beta_mu_mean[2, 1];
        ////
        real Xbeta_baseline_nd = dot_product(baseline_case_nd, to_vector(beta_mu[1, 1:n_covariates_nd]));
        real Xbeta_baseline_d  = dot_product(baseline_case_d,  to_vector(beta_mu[2, 1:n_covariates_d]));
        ////
        for (c in 1:2) {
           log_kappa[c] ~ nico_student_t(prior_kappa_dof[c], prior_kappa_mean[c], prior_kappa_SD[c]);
           // kappa[c] ~ normal(0, 100);
           target += -log_kappa[c];
           ////
           if (use_Se_Sp_prior == 1 && n_chosen_cutoffs[c] > 0) {

                    matrix[n_chosen_cutoffs[c], 2] prior_beta_of_this_class = prior_Se_Sp_beta[c][1:n_chosen_cutoffs[c], ];
                    matrix[n_chosen_cutoffs[c], 2] prior_link_mean_SD_of_this_class = prior_Se_Sp_link_mean_SD[c][1:n_chosen_cutoffs[c], ];
                    ////
                    if (c == 2) {
                         prior_link_mean_SD_of_this_class[, 1] = - prior_link_mean_SD_of_this_class[, 1];
                         if (Se_Sp_prior_type == 2) {
                              prior_beta_of_this_class = append_col(prior_beta_of_this_class[, 2],
                                                                    prior_beta_of_this_class[, 1]);
                         }
                    }
                    real Xbeta_baseline = (c == 1) ? Xbeta_baseline_nd : Xbeta_baseline_d;
                    ////
                    //// (the shapes of the Beta of each cutpoint from alpha itself: sum(alpha) = kappa + 0.01 * n_cat)
                    ////
                    real alpha_0 = sum(alpha[c]);
                    vector[n_thr] P_of_alpha = cumulative_sum(alpha[c][1:n_thr]) / alpha_0;
                    vector[n_thr] P_upper_of_alpha = reverse(cumulative_sum(reverse(alpha[c][2:n_cat]))) / alpha_0;
                    vector[n_chosen_cutoffs[c]] C_mu;
                    vector[n_chosen_cutoffs[c]] log_J;
                    ////
                    for (j in 1:n_chosen_cutoffs[c]) {
                         int k = chosen_cutoffs[c, j];
                         vector[2] median_and_log_J = dirichlet_cutpoint_median_and_log_J( P_of_alpha[k],
                                                                                           P_upper_of_alpha[k],
                                                                                           alpha_0,
                                                                                           use_probit_link);
                         C_mu[j]  = Xbeta_baseline + median_and_log_J[1]; // Xbeta_baseline = anchor of the Se / Sp prior route.
                         log_J[j] = median_and_log_J[2] + nico_log(kappa[c] / alpha_0);
                    }
                    ////
                    target += dirichlet_with_Se_Sp_prior_lupdf(
                              dirichlet_phi[c] | log_J,
                              C_mu - Xbeta_baseline,
                              prior_dirichlet_alpha[c],
                              alpha_of_the_segments[c][1:(n_chosen_cutoffs[c] + 1)],
                              chosen_cutoffs[c, 1:n_chosen_cutoffs[c]], Se_Sp_prior_type,
                              prior_Se_Sp_dirichlet_A[c][1:(n_chosen_cutoffs[c] + 1)], prior_beta_of_this_class,
                              prior_link_mean_SD_of_this_class,
                              use_probit_link);

           } else {

                     dirichlet_phi[c] ~ nico_dirichlet(prior_dirichlet_alpha[c]);

           }
        }
        ////
        //// ---- Induced-dirichlet between study ** model ** (NOT a prior model here but part of the actual likelihood since random-effect cutpoints!):
        ////
        // for (c in 1:2) {
        //      for (s in 1:n_studies) {
        //          vector[n_cat] Ind_Dir_ord_prob;
        //          if (use_probit_link == 1) Ind_Dir_ord_prob = to_vector(cumul_probs_to_ord_probs(Phi(C[c][s, ])));
        //          else                      Ind_Dir_ord_prob = to_vector(cumul_probs_to_ord_probs(inv_logit(C[c][s, ])));
        //          ////
        //          Ind_Dir_ord_prob ~ induced_dirichlet_given_C(
        //                             to_vector(C[c][s, ]), alpha[c], use_probit_link);
        //      }
        // }

        //// ---- Induced-dirichlet between study ** model ** (NOT a prior: part of the likelihood since
        ////      random-effect cutpoints), every study's cutpoints in one node per class, as
        ////      induced_dirichlet_given_C on each row:
        for (c in 1:2) {
             //// (anchored at the prior mean of the intercept; at Xbeta_baseline under the Se / Sp prior route, as
             ////  the fixed-cutpoint model's Se / Sp route is)
             real C_anchor = (use_Se_Sp_prior == 1) ? ((c == 1) ? Xbeta_baseline_nd : Xbeta_baseline_d) : baseline_location[c];
             //// (the reported cutpoints of each study only: the alpha of each segment between them is the sum of
             ////  alpha[c] over its categories, the same distribution with the other cutpoints integrated out)
             vector[n_cat + 1] cumulative_alpha = append_row(0.0, cumulative_sum(alpha[c]));
             vector[n_C + n_studies] alpha_of_each_segment = cumulative_alpha[upper_index_of_each_segment] -
                                                             cumulative_alpha[lower_index_of_each_segment];
             for (s in 1:n_studies) {
                  if (n_reported_cutpoints[s] > 0) {
                       int p_1 = first_C_of_each_study[s];
                       int p_2 = p_1 + n_reported_cutpoints[s] - 1;
                       target += induced_dirichlet_given_C_on_the_log_scale_lupdf(C[c][p_1:p_2] - C_anchor |
                                                                                  alpha_of_each_segment[(p_1 + s - 1):(p_2 + s)],
                                                                                  use_probit_link);
                  }
             }
        }
        ////
        //// ---- Likelihood / Model:
        ////
        // target += std_normal_lpdf(to_vector(beta_z)); // part of between-study model, NOT prior // beta_{c, s} ~ normal(beta_mu_{c}, beta_SD_{c})

        // target += normal_lupdf(to_vector(beta_z) | 0, 1); // part of between-study model, NOT prior
        target += nico_normal_lupdf(to_vector(to_matrix(beta_z)) | 0, 1); // part of between-study model, NOT prior
        ////
        for (c in 1:2) {
           //// (a study's first reported cutpoint is C_raw itself; each later one adds an exp() or softplus step)
           target += (softplus == 1) ? sum(nico_log_inv_logit(C_raw[c][position_of_each_increment])) :
                                       sum(C_raw[c][position_of_each_increment]);
        }
        ////
        //// ---- Log-likelihood:
        ////
        // {
        //     ////
        //     //// ---- Get the cutpoint index (k) to map "latent_surv[c][s, cut_i]" to correct cutpoint "C[k]":
        //     ////
        //     // real scale = 1.0; // since using "ord_bivariate-like"" param.
        //     array[2] matrix[n_studies, n_thr] latent_surv = map_latent_surv_to_random_hetero_C(
        //                                                     C, n_thr, Xbeta, n_studies,
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
            array[2] matrix[n_studies, n_thr] latent_surv = map_latent_surv_to_the_reported_C(
                                                            C, n_thr, Xbeta, n_studies,
                                                            n_obs_cutpoints, position_of_each_C);
            target += compute_log_lik_binomial_fact_lp(
                      latent_surv, use_probit_link, n_thr, n_studies,
                      x_2, n, N_total, n_obs_cutpoints,
                      K_fold_CV_indicator);
        }

}

generated quantities {
    real simd_lanes = bmvp_simd_lanes();
      ////
      //// ---- Compute between-study variance-covariance matrix for location parameters:
      ////
      cov_matrix[2]  beta_Sigma = multiply_lower_tri_self_transpose(beta_L_Sigma);
      corr_matrix[2] beta_Omega = multiply_lower_tri_self_transpose(beta_L_Omega);
      // ////
      // //// ---- SD's for each category probability in a Dirichlet is √(α_k(α_0-α_k)/(α_0^2(α_0+1))) - where α_0 is the sum of all alphas:
      // ////
      // array[2] vector<lower=0.0>[n_cat] dirichlet_cat_SDs_sigma;
      // for (c in 1:2) {
      //     real alpha_0 = sum(alpha[c]);
      //     dirichlet_cat_SDs_sigma[c] = sqrt((alpha[c] .* (alpha_0 - alpha[c]) ) / 
      //                                  (square(alpha_0) * (alpha_0 + 1.0)));
      // }
      ////
      // //// ---- Calculate summary cutpoints:
      // ////
      array[2] vector[n_thr] C_mu;
      // ////
      // for (c in 1:2) {
      //     for (k in 1:n_thr) {
      //         vector[n_studies] observed_cutpoints;
      //         int n_obs = 0;
      //         ////
      //         for (s in 1:n_studies) {
      //             // Search through the cutpoint_index for this study to see if cutpoint k was observed
      //             for (j in 2:(n_obs_cutpoints[s] + 1)) {  // Start at 2 because index 1 is 0
      //                 if (cutpoint_index[c, s, j] == k) {
      //                     n_obs += 1;
      //                     observed_cutpoints[n_obs] = C[c][s, k];
      //                     break;  // Found it, no need to continue searching
      //                 }
      //             }
      //         }
      //         if (n_obs > 0) {
      //             C_mu[c][k] = median(observed_cutpoints[1:n_obs]);
      //         }
      //     }
      // }
      for (c in 1:2) {
         //// (the shapes of the Beta of each cutpoint from alpha itself: sum(alpha) = kappa + 0.01 * n_cat)
         real alpha_0 = sum(alpha[c]);
         vector[n_thr] P = cumulative_sum(alpha[c][1:n_thr]) / alpha_0;
         vector[n_thr] P_upper = reverse(cumulative_sum(reverse(alpha[c][2:n_cat]))) / alpha_0;
         //// (anchored at the prior mean of the intercept; at Xbeta_baseline under the Se / Sp prior route)
         real C_anchor = (use_Se_Sp_prior == 1) ?
              ((c == 1) ? dot_product(baseline_case_nd, to_vector(beta_mu[1, 1:n_covariates_nd])) :
                          dot_product(baseline_case_d, to_vector(beta_mu[2, 1:n_covariates_d]))) :
              prior_beta_mu_mean[c, 1];
         for (k in 1:n_thr) {
            C_mu[c][k] = C_anchor +
                         dirichlet_cutpoint_median(P[k], P_upper[k], alpha_0, use_probit_link);
         }
      }
      //
      // for (c in 1:2) {
      //   vector[n_cat - 1] cumul_probs = ord_probs_to_cumul_probs(dirichlet_cat_SDs_sigma[c]);
      //   C_mu[c] = ID_cumul_probs_to_C(cumul_probs, use_probit_link);
      // }
      ////
      //// ---- Calculate summary accuracy - "baseline" covariate values:
      ////  
      real Xbeta_baseline_nd = dot_product(baseline_case_nd, to_vector(beta_mu[1, 1:n_covariates_nd])); 
      real Xbeta_baseline_d  = dot_product(baseline_case_d,  to_vector(beta_mu[2, 1:n_covariates_d]));
      ////
      //// ---- Calculate baseline Se/Sp:
      ////
      vector[n_thr] Fp_baseline = (use_probit_link == 1) ? 
                                  nico_Phi( -(C_mu[1] - Xbeta_baseline_nd))  : 
                                  nico_inv_logit(-(C_mu[1] - Xbeta_baseline_nd));
      vector[n_thr] Sp_baseline = 1.0 - Fp_baseline;
      vector[n_thr] Se_baseline = (use_probit_link == 1) ? 
                                  nico_Phi( -(C_mu[2] - Xbeta_baseline_d))   : 
                                  nico_inv_logit(-(C_mu[2] - Xbeta_baseline_d));
      ////
      //// ---- Calculate predictive accuracy:
      ////
      array[2] vector[n_thr] C_pred;
      ////
      for (c in 1:2) {
        for (k in 1:n_thr) {
            real sigma_C_k;
            real SD_prob_scale = dirichlet_cat_SDs_sigma[c][k];
            if (use_probit_link == 1) {
                real mean_C_cutpoint_scale = C_mu[c][k];
                sigma_C_k = SD_approx_ID_ord_prob_to_C_probit(mean_C_cutpoint_scale, SD_prob_scale);
            } else {
                real mean_prob_scale = dirichlet_phi[c][k];
                sigma_C_k = SD_approx_ID_ord_prob_to_C_probit(mean_prob_scale, SD_prob_scale);
            }
            if (!is_nan(sigma_C_k)) {
              C_pred[c][k] = C_mu[c][k] + nico_normal_rng(0.0, sigma_C_k);
            }
        }
      }
      ////
      vector[2] beta_random_pred = multi_normal_cholesky_rng(rep_vector(0.0, 2), beta_L_Sigma);
      ////
      //// ---- Use baseline covariates with predicted random effects
      ////
      real Xbeta_baseline_pred_nd = dot_product(baseline_case_nd, beta_mu[1, 1:n_covariates_nd]) + beta_random_pred[1];
      real Xbeta_baseline_pred_d  = dot_product(baseline_case_d,  beta_mu[2, 1:n_covariates_d])  + beta_random_pred[2];
      ////
      vector[n_thr] Fp_baseline_pred = (use_probit_link == 1) ?
                                       nico_Phi(-(C_pred[1] - Xbeta_baseline_pred_nd)) :
                                       nico_inv_logit(-(C_pred[1] - Xbeta_baseline_pred_nd));
      vector[n_thr] Sp_baseline_pred = 1.0 - Fp_baseline_pred;
      vector[n_thr] Se_baseline_pred = (use_probit_link == 1) ?
                                       nico_Phi(-(C_pred[2] - Xbeta_baseline_pred_d))  :
                                       nico_inv_logit(-(C_pred[2] - Xbeta_baseline_pred_d));
      ////
      //// ---- Log-lik + study-specific accuracy computation (using "data" / double-precision fn for efficiency):
      ////
      vector[n_studies] log_lik_study = rep_vector(0.0, n_studies);
      // ////
      // matrix[n_studies, n_thr] fp; // global
      // matrix[n_studies, n_thr] sp; // global
      // matrix[n_studies, n_thr] se; // global
      // ////
      // vector[n_studies] deviance_nd; // global
      // vector[n_studies] deviance_d;  // global
      // vector[n_studies] deviance;    // global
      if (log_lik_indicator == 1) {
          real scale = 1.0; // since using "ord_bivariate-like"" param.
          array[2] matrix[n_studies, n_thr] latent_surv = map_latent_surv_to_the_reported_C(
                                                          C, n_thr, Xbeta, n_studies,
                                                          n_obs_cutpoints, position_of_each_C);
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
          ////
          // for (s in 1:n_studies) {
          //     for (k in 1:to_int(n_obs_cutpoints[s])) {
          //        deviance_nd[s] += is_nan(dev_nd[s, k]) ? 0.0 : dev_nd[s, k];
          //        deviance_d[s]  += is_nan(dev_d[s, k]) ? 0.0 : dev_d[s, k];
          //     }
          //     deviance[s] = deviance_nd[s] + deviance_d[s];
          // }
      }

}







