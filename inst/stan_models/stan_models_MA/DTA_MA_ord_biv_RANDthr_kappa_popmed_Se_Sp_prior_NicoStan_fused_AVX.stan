// NicoSIMD_AVX512 version of DTA_MA_ord_biv_RANDthr_kappa_popmed_Se_Sp_prior_fused_base.stan written by NicoStan (source model md5 a62a5611d9cfc416b997934fd8574159; NicoSIMD declarations md5 7f272e2ac00d6303123b13fb8f59653e; writer implementation md5 f6843a34264d79c6ed96effa2d04b0c4)
// NicoSIMD functions used: exp, log, log1m, log1p_exp, log1m_exp, inv_logit, log_inv_logit, Phi, inv_Phi, tanh, logit, log1m_inv_logit, expm1, normal_rng, normal, lgamma, beta, dirichlet, lkj_corr_cholesky, log_diff_exp, std_normal_log_qf, student_t, binomial_counts_above_successive_cutpoints_of_a_cumulative_link
//// DTA_MA_ord_biv_RANDthr_kappa_popmed_Se_Sp_prior_fused_base.stan
//// Base model for O-biv-RC fused AVX generation via NicoStan::R_fn_write_BayesMVP_AVX_version_of_Stan_model_and_user_header.
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
    real nico_binomial_counts_above_successive_cutpoints_of_a_cumulative_link_lpmf(
          array[] int n, array[] int n_previous, vector latent, data vector is_first, data int link,
          data real bound);
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
        //// ---- Fused NicoSIMD C++ functions (bodies in NicoSIMD):
        ////
        vector nico_pinkney_C_and_log_J_from_alpha(vector alpha,
                                                   vector C_raw_c,
                                                   real baseline_location_prior,
                                                   data vector Pinkney_w_c,
                                                   data vector logit_start,
                                                   data vector C_raw_lower_c,
                                                   data vector C_raw_upper_c,
                                                   data array[] int k_of_each_C,
                                                   data array[] int k_previous_of_each_C,
                                                   data int use_probit_link,
                                                   data int softplus);
        real nico_induced_dirichlet_prior_of_ordered_cutpoints_ragged_lpdf(
                vector C,
                vector alpha,
                data array[] int n_cutpoints_per_group,
                data array[] int cutpoint_k_of_each_C,
                data int link,
                real location);
        real nico_binomial_counts_above_successive_cutpoints_of_a_cumulative_link_lpmf(
                array[,,] int n,
                array[,,] int n_previous,
                matrix locations,
                array[] vector C,
                data array[,] int position_of_each_C,
                data array[] int n_obs_cutpoints,
                data array[] int K_fold_CV_indicator,
                data int link,
                data real bound);
        real nico_induced_dirichlet_with_Se_Sp_prior_lpdf(vector C, vector alpha, data int link,
                                                          data real location, vector alpha_of_the_segments,
                                                          data array[] int chosen_cutoffs,
                                                          data int Se_Sp_prior_type,
                                                          vector prior_Se_Sp_dirichlet_A, matrix prior_Se_Sp_beta,
                                                          matrix prior_Se_Sp_link_mean_SD);
        vector nico_dirichlet_cutpoint_medians_and_sum_of_log_J(vector P_lower, vector P_upper, real kappa,
                                                                data int use_probit_link);
        vector nico_dirichlet_cutpoint_medians_and_log_J_of_each(vector P_lower, vector P_upper, real kappa,
                                                                 data int use_probit_link);
        real nico_log_Set0_marginal_at_chosen_cutoffs(vector x_chosen, vector alpha_of_the_segments,
                                                      real v0, int use_probit_link);
        real nico_Se_Sp_prior_at_chosen_cutoffs_lpdf(vector x, data int Se_Sp_prior_type,
                                                     vector prior_Se_Sp_dirichlet_A, matrix prior_Se_Sp_beta,
                                                     matrix prior_Se_Sp_link_mean_SD, data int link);
        real nico_dirichlet_with_Se_Sp_prior_lpdf(vector phi,
                                                  vector log_J,
                                                  vector C_baseline,
                                                  vector alpha,
                                                  vector alpha_of_the_segments,
                                                  data array[] int chosen_cutoffs,
                                                  data int Se_Sp_prior_type,
                                                  vector prior_Se_Sp_dirichlet_A,
                                                  matrix prior_Se_Sp_beta,
                                                  matrix prior_Se_Sp_link_mean_SD,
                                                  data int link);
        array[,] matrix nico_compute_log_lik_binomial_fact_data(data array[] matrix latent_surv,
                                                                data int use_probit_link, data int n_thr,
                                                                data int n_studies, data array[,,] int x_2,
                                                                data array[,,] int n, data array[,] int N_total,
                                                                data array[] int n_obs_cutpoints,
                                                                data array[] int K_fold_CV_indicator);
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
        int<lower=0, upper=2> use_Se_Sp_prior;     //// 0 = original dirichlet_phi prior; 1 = Se / Sp prior; 2 = prior set 0 through the Se / Sp prior
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
        //// ---- Under the Se / Sp prior the set 0 prior of (summary cutpoints, intercept) is kept and what it implies
        ////      for C_mu - Xbeta_baseline is divided out (model block). Both are written for the probit link and for
        ////      prior_dirichlet_alpha = (1, ..., 1) only:
        ////
        if (use_Se_Sp_prior == 1) {
              if (use_probit_link != 1) reject("the Se / Sp prior of this model needs the probit link");
              for (c in 1:2) {
                    for (j in 1:n_cat) {
                          if (prior_dirichlet_alpha[c][j] != 1.0) {
                                reject("the Se / Sp prior of this model needs prior_dirichlet_alpha = (1, ..., 1): c = ", c);
                          }
                    }
              }
        }
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
                   if (use_Se_Sp_prior >= 1 && Se_Sp_prior_type == 4 &&
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
        //// (the cutpoint index k of each reported cutpoint and that of the reported cutpoint before it in its study, 0 for
        ////  a study's first: for the centre and the spread of each element of C_raw)
        array[n_C] int k_of_each_C;
        array[n_C] int k_previous_of_each_C;
        //// (the study of each reported cutpoint)
        array[n_C] int study_of_each_C;
        //// (where a free value of 0 puts C_raw when the weight is 0: the first cutpoint at start_of_C_raw and every
        ////  one-cutpoint step at exp(start_of_C_raw), the start the package gives C_raw itself)
        real start_of_C_raw = -2.0;
        vector[n_C] logit_of_the_place_of_the_start;
        if (softplus == 1) reject("this copy writes the centre and the spread of C_raw for the exp() steps only");
        array[n_C + n_studies] int upper_index_of_each_segment;  //// a segment = the categories k_(j - 1) + 1 to k_j
        array[n_C + n_studies] int lower_index_of_each_segment;
        {
              //// (each one-cutpoint step lies between these, as each column of the full C_raw did; a gap of g
              ////  cutpoints is g such steps)
              real step_lower = (softplus == 1) ? log1p_exp(C_raw_lower) : exp(C_raw_lower);
              real step_upper = (softplus == 1) ? log1p_exp(C_raw_upper) : exp(C_raw_upper);
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
                              k_of_each_C[p] = k;
                              k_previous_of_each_C[p] = k_previous;
                              study_of_each_C[p] = s;
                              upper_index_of_each_segment[p + s - 1] = k + 1;
                              lower_index_of_each_segment[p + s - 1] = k_previous + 1;
                              if (k_previous == 0) {
                                   C_raw_lower_of_each_C[1][p] = C_raw_lower + (k - 1) * step_lower;
                                   C_raw_upper_of_each_C[1][p] = C_raw_upper + (k - 1) * step_upper;
                              } else {
                                   i += 1;
                                   position_of_each_increment[i] = p;
                                   C_raw_lower_of_each_C[1][p] = (softplus == 1) ? log(expm1((k - k_previous) * step_lower)) :
                                                                                 C_raw_lower + log(k - k_previous);
                                   C_raw_upper_of_each_C[1][p] = (softplus == 1) ? log(expm1((k - k_previous) * step_upper)) :
                                                                                 C_raw_upper + log(k - k_previous);
                              }
                              k_previous = k;
                         }
                    }
                    upper_index_of_each_segment[p + s] = n_cat + 1;
                    lower_index_of_each_segment[p + s] = k_previous + 1;
              }
              C_raw_lower_of_each_C[2] = C_raw_lower_of_each_C[1];
              C_raw_upper_of_each_C[2] = C_raw_upper_of_each_C[1];
              for (q in 1:n_C) {
                    real start_of_this_C_raw = (k_previous_of_each_C[q] == 0) ?
                                               start_of_C_raw + (k_of_each_C[q] - 1) * exp(start_of_C_raw) :
                                               start_of_C_raw + log(k_of_each_C[q] - k_previous_of_each_C[q]);
                    real place_of_the_start = (start_of_this_C_raw - C_raw_lower_of_each_C[1][q]) /
                                              (C_raw_upper_of_each_C[1][q] - C_raw_lower_of_each_C[1][q]);
                    logit_of_the_place_of_the_start[q] = logit(fmin(fmax(place_of_the_start, 0.001), 0.999));
              }
        }
        //// ---- Data-only partial-centering weights Pinkney_w in [0, 1] and empirical origins Pinkney_eta0
        ////      (Sean Pinkney, partially_centered_example_Sean_Pinkney.stan & spinkney/hierarchical_model_tutorial
        ////       stan/hierarchical_pcp.stan, R/1_hier_code.R: w = 0 is centered, w = 1 is non-centered;
        ////       dimensionless precision-variance ratio w = 1 / (1 + J_i * tau2_pos) computed once from x and
        ////       K_fold_CV_indicator only, never from parameters, so dw / dkappa = 0 everywhere):
        array[2] vector<lower=0.0, upper=1.0>[n_C] Pinkney_w;
        array[2] vector[n_C] Pinkney_eta0;
        {
              int max_dl_pairs = n_thr * (n_thr + 1);
              for (c in 1:2) {
                    vector[n_C] J_of_each_C = rep_vector(0.0, n_C);
                    array[n_C] int active_of_each_C;
                    matrix[n_thr, n_thr + 1] dl_pair_lookup = rep_matrix(0.0, n_thr, n_thr + 1);
                    int n_dl_pairs = 0;
                    array[max_dl_pairs] int n_studies_of_pair = rep_array(0, max_dl_pairs);
                    vector[max_dl_pairs] W1_of_pair  = rep_vector(0.0, max_dl_pairs);
                    vector[max_dl_pairs] W2_of_pair  = rep_vector(0.0, max_dl_pairs);
                    vector[max_dl_pairs] WY_of_pair  = rep_vector(0.0, max_dl_pairs);
                    vector[max_dl_pairs] WY2_of_pair = rep_vector(0.0, max_dl_pairs);
                    int p_idx = 0;
                    for (s in 1:n_studies) {
                          int N_cs = N_total[c, s];
                          int k_prev = 0;
                          real p_prev = 0.0;
                          real z_prev = 0.0;
                          real f_prev = 1.0;
                          int x_ab_prev = N_cs;
                          for (cut_i in 1:n_obs_cutpoints[s]) {
                                int k = cutpoint_index[c, s, cut_i + 1];
                                if (k != -1) {
                                      p_idx += 1;
                                      real lb_p = C_raw_lower_of_each_C[c][p_idx];
                                      real wd_p = C_raw_upper_of_each_C[c][p_idx] - lb_p;
                                      if (K_fold_CV_indicator[s] == 0 || N_cs <= 0) {
                                            active_of_each_C[p_idx] = 0;
                                            Pinkney_eta0[c][p_idx] = logit_of_the_place_of_the_start[p_idx];
                                      } else {
                                            active_of_each_C[p_idx] = 1;
                                            int x_ab = x_2[c, s, cut_i];
                                            real p_raw = (N_cs - x_ab + 0.5) / (N_cs + 1.0);
                                            real p_below = (k_prev > 0 && p_raw <= p_prev) ?
                                                           fmin(0.999999, p_prev + 0.5 / ((N_cs + 1.0) * (n_thr + 1.0))) :
                                                           p_raw;
                                            real z_hat = (use_probit_link == 1) ? inv_Phi(p_below) : logit(p_below);
                                            real fz = (use_probit_link == 1) ?
                                                      (exp(-0.5 * square(z_hat)) / sqrt(2.0 * pi())) :
                                                      (p_below * (1.0 - p_below));
                                            real vf_curr = p_below * (1.0 - p_below) / square(fz);
                                            int u = to_int(dl_pair_lookup[k, k_prev + 1]);
                                            if (u == 0) {
                                                  n_dl_pairs += 1;
                                                  u = n_dl_pairs;
                                                  dl_pair_lookup[k, k_prev + 1] = u;
                                            }
                                            if (k_prev == 0) {
                                                  real q_hat = fmin(fmax((z_hat - lb_p) / wd_p, 0.005), 0.995);
                                                  Pinkney_eta0[c][p_idx] = logit(q_hat);
                                                  J_of_each_C[p_idx] = N_cs / fmax(vf_curr, 1e-12);
                                            } else {
                                                  real vf_prev = p_prev * (1.0 - p_prev) / square(f_prev);
                                                  real cov_term = p_prev * (1.0 - p_below) / (f_prev * fz);
                                                  real var_gap = fmax((vf_curr + vf_prev - 2.0 * cov_term) / N_cs, 1e-10);
                                                  real gap = fmax(z_hat - z_prev, 1e-3);
                                                  real J_gap = 1.0 / var_gap;
                                                  real r_hat = log(gap);
                                                  real q_hat = fmin(fmax((r_hat - lb_p) / wd_p, 0.005), 0.995);
                                                  Pinkney_eta0[c][p_idx] = logit(q_hat);
                                                  J_of_each_C[p_idx] = J_gap;
                                                  if ((x_ab_prev - x_ab) >= 1) {
                                                        n_studies_of_pair[u] += 1;
                                                        W1_of_pair[u]  += J_gap;
                                                        W2_of_pair[u]  += square(J_gap);
                                                        WY_of_pair[u]  += J_gap * gap;
                                                        WY2_of_pair[u] += J_gap * square(gap);
                                                  }
                                            }
                                            p_prev = p_below;
                                            z_prev = z_hat;
                                            f_prev = fz;
                                            x_ab_prev = x_ab;
                                      }
                                      k_prev = k;
                                }
                          }
                    }
                    real Q_tot = 0.0;
                    real df_tot = 0.0;
                    real C_tot = 0.0;
                    for (u in 1:n_dl_pairs) {
                          if (n_studies_of_pair[u] >= 2 && W1_of_pair[u] > 0.0) {
                                real Q_u = fmax(0.0, WY2_of_pair[u] - square(WY_of_pair[u]) / W1_of_pair[u]);
                                real C_u = W1_of_pair[u] - W2_of_pair[u] / W1_of_pair[u];
                                if (C_u > 0.0) {
                                      Q_tot  += Q_u;
                                      df_tot += n_studies_of_pair[u] - 1.0;
                                      C_tot  += C_u;
                                }
                          }
                    }
                    real tau2_pos = 0.0;
                    if (C_tot > 0.0 && df_tot > 0.0) {
                          real tau2_raw = (Q_tot - df_tot) / C_tot;
                          real se_tau2  = sqrt(2.0 * df_tot) / C_tot;
                          tau2_pos = 0.5 * (tau2_raw + sqrt(square(tau2_raw) + square(se_tau2)));
                    }
                    for (p in 1:n_C) {
                          Pinkney_w[c][p] = (active_of_each_C[p] == 0) ? 0.0 : (1.0 / (1.0 + J_of_each_C[p] * tau2_pos));
                    }
              }
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
        //// (one free value for each reported cutpoint of each study: C_raw, between the same bounds as before, is built
        ////  from it in transformed parameters as centre + spread x free value, with the Jacobian in the model block)
        array[2] vector[n_C] C_raw;
        vector<lower=1.0>[2] kappa;
        array[2] simplex[n_cat] dirichlet_phi;
        // array[2] simplex[n_thr] raw_dirichlet_phi;
}


transformed parameters {
        vector[2] log_kappa = nico_log(kappa);
        array[2] vector<lower=0.0>[n_cat] alpha;
        for (c in 1:2) {
              alpha[c] = 0.01 + dirichlet_phi[c] * kappa[c];
        }
        //// ---- Fused Pinkney partially-centered study cutpoints and Jacobian directly from alpha:
        array[2] vector[n_C] C;
        vector[2] total_log_J_of_C = rep_vector(0.0, 2);
        {
              vector[2] baseline_location_prior_of_C_raw;
              baseline_location_prior_of_C_raw[1] = dot_product(baseline_case_nd, to_vector(prior_beta_mu_mean[1, 1:n_covariates_nd]));
              baseline_location_prior_of_C_raw[2] = dot_product(baseline_case_d,  to_vector(prior_beta_mu_mean[2, 1:n_covariates_d]));
              for (c in 1:2) {
                    vector[n_C + 1] C_and_log_J = nico_pinkney_C_and_log_J_from_alpha(
                          alpha[c],
                          C_raw[c],
                          baseline_location_prior_of_C_raw[c],
                          Pinkney_w[c],
                          Pinkney_eta0[c],
                          C_raw_lower_of_each_C[c],
                          C_raw_upper_of_each_C[c],
                          k_of_each_C,
                          k_previous_of_each_C,
                          use_probit_link,
                          softplus);
                    C[c] = C_and_log_J[1:n_C];
                    total_log_J_of_C[c] = C_and_log_J[n_C + 1];
              }
        }
        ////
        cholesky_factor_corr[2] beta_L_Omega = make_bivariate_L_Omega(beta_corr);
        cholesky_factor_cov[2]  beta_L_Sigma = diag_pre_multiply(beta_SD, beta_L_Omega);
        matrix[n_studies, 2] beta_random = to_matrix(beta_z) * beta_L_Sigma;
        matrix[n_studies, 2] Xbeta;
        Xbeta[, 1] = X_nd[1:n_studies, 1:n_covariates_nd] * to_vector(beta_mu[1, 1:n_covariates_nd]) + beta_random[, 1];
        Xbeta[, 2] = X_d[1:n_studies, 1:n_covariates_d]   * to_vector(beta_mu[2, 1:n_covariates_d])  + beta_random[, 2];
}


model {
        //// ---- Priors:
        beta_mu[1, 1:n_covariates_nd] ~ nico_normal(
                                          prior_beta_mu_mean[1, 1:n_covariates_nd],
                                          prior_beta_mu_SD[1, 1:n_covariates_nd]);
        beta_mu[2, 1:n_covariates_d]  ~ nico_normal(
                                          prior_beta_mu_mean[2, 1:n_covariates_d],
                                          prior_beta_mu_SD[2, 1:n_covariates_d]);
        to_vector(beta_SD) ~ nico_normal(prior_beta_SD_mean, prior_beta_SD_SD);
        beta_L_Omega ~ nico_lkj_corr_cholesky(prior_beta_corr_LKJ);
        ////
        vector[2] baseline_location_prior;
        baseline_location_prior[1] = dot_product(baseline_case_nd, to_vector(prior_beta_mu_mean[1, 1:n_covariates_nd]));
        baseline_location_prior[2] = dot_product(baseline_case_d,  to_vector(prior_beta_mu_mean[2, 1:n_covariates_d]));
        ////
        for (c in 1:2) {
           log_kappa[c] ~ nico_student_t(prior_kappa_dof[c], prior_kappa_mean[c], prior_kappa_SD[c]);
           target += -log_kappa[c];
           ////
           if (use_Se_Sp_prior >= 1 && n_chosen_cutoffs[c] > 0) {
                    matrix[n_chosen_cutoffs[c], 2] prior_beta_of_this_class = prior_Se_Sp_beta[c][1:n_chosen_cutoffs[c], ];
                    matrix[n_chosen_cutoffs[c], 2] prior_link_mean_SD_of_this_class = prior_Se_Sp_link_mean_SD[c][1:n_chosen_cutoffs[c], ];
                    if (c == 2) {
                         prior_link_mean_SD_of_this_class[, 1] = - prior_link_mean_SD_of_this_class[, 1];
                         if (Se_Sp_prior_type == 2) {
                              prior_beta_of_this_class = append_col(prior_beta_of_this_class[, 2],
                                                                    prior_beta_of_this_class[, 1]);
                         }
                    }
                    real Xbeta_baseline = (c == 1) ?
                        dot_product(baseline_case_nd, to_vector(beta_mu[1, 1:n_covariates_nd])) :
                        dot_product(baseline_case_d,  to_vector(beta_mu[2, 1:n_covariates_d]));
                    real alpha_0 = sum(alpha[c]);
                    if (use_Se_Sp_prior == 2) {
                         vector[n_thr] P_of_alpha = cumulative_sum(alpha[c][1:n_thr]) / alpha_0;
                         vector[n_thr] P_upper_of_alpha = reverse(cumulative_sum(reverse(alpha[c][2:n_cat]))) / alpha_0;
                         vector[n_thr + 1] popmed_of_every_cutpoint_and_sum_of_log_J =
                               nico_dirichlet_cutpoint_medians_and_sum_of_log_J(P_of_alpha, P_upper_of_alpha, alpha_0, use_probit_link);
                         real log_J_of_every_median = popmed_of_every_cutpoint_and_sum_of_log_J[n_thr + 1] +
                                                      n_thr * nico_log(kappa[c] / alpha_0);
                         target += log_J_of_every_median;
                         target += nico_induced_dirichlet_with_Se_Sp_prior_lupdf(
                                   popmed_of_every_cutpoint_and_sum_of_log_J[1:n_thr] | prior_dirichlet_alpha[c],
                                   use_probit_link,
                                   0.0,
                                   alpha_of_the_segments[c][1:(n_chosen_cutoffs[c] + 1)],
                                   chosen_cutoffs[c, 1:n_chosen_cutoffs[c]],
                                   Se_Sp_prior_type,
                                   prior_Se_Sp_dirichlet_A[c][1:(n_chosen_cutoffs[c] + 1)],
                                   prior_beta_of_this_class,
                                   prior_link_mean_SD_of_this_class);
                    } else {
                         int m_c = n_chosen_cutoffs[c];
                         array[m_c] int cuts_c = chosen_cutoffs[c, 1:m_c];
                         vector[m_c + 1] alpha_seg_c = alpha_of_the_segments[c][1:(m_c + 1)];
                         vector[m_c + 1] seg_alpha_var;
                         {
                              int k_prev_c = 0;
                              for (j in 1:m_c) {
                                   int k_c = cuts_c[j];
                                   seg_alpha_var[j] = sum(alpha[c][(k_prev_c + 1):k_c]);
                                   k_prev_c = k_c;
                              }
                              seg_alpha_var[m_c + 1] = sum(alpha[c][(k_prev_c + 1):n_cat]);
                         }
                         vector[m_c] P_chosen = cumulative_sum(seg_alpha_var[1:m_c]) / alpha_0;
                         vector[m_c] P_upper_chosen;
                         {
                              real rem = seg_alpha_var[m_c + 1];
                              for (j_rev in 1:m_c) {
                                   int j = m_c + 1 - j_rev;
                                   P_upper_chosen[j] = rem / alpha_0;
                                   rem += seg_alpha_var[j];
                              }
                         }
                         vector[2 * m_c] popmed_and_J_chosen =
                               nico_dirichlet_cutpoint_medians_and_log_J_of_each(P_chosen, P_upper_chosen,
                                                                                 alpha_0, use_probit_link);
                         vector[m_c] popmed_chosen = popmed_and_J_chosen[1:m_c];
                         vector[m_c] log_J_chosen  = popmed_and_J_chosen[(m_c + 1):(2 * m_c)] + nico_log(kappa[c] / alpha_0);
                         vector[m_c] x_chosen = (baseline_location_prior[c] + popmed_chosen) - Xbeta_baseline;
                         real v0 = (c == 1) ?
                                   dot_self(baseline_case_nd .* to_vector(prior_beta_mu_SD[1, 1:n_covariates_nd])) :
                                   dot_self(baseline_case_d  .* to_vector(prior_beta_mu_SD[2, 1:n_covariates_d]));
                         target += nico_dirichlet_with_Se_Sp_prior_lupdf(
                                   dirichlet_phi[c] | log_J_chosen, popmed_chosen,
                                   prior_dirichlet_alpha[c], alpha_seg_c, cuts_c, 1,
                                   alpha_seg_c, prior_beta_of_this_class, prior_link_mean_SD_of_this_class, use_probit_link);
                         target += - nico_log_Set0_marginal_at_chosen_cutoffs(x_chosen, alpha_seg_c, v0, use_probit_link);
                         target += nico_Se_Sp_prior_at_chosen_cutoffs_lupdf(
                                   x_chosen | Se_Sp_prior_type,
                                   prior_Se_Sp_dirichlet_A[c][1:(m_c + 1)],
                                   prior_beta_of_this_class, prior_link_mean_SD_of_this_class, use_probit_link);
                    }
           } else {
                     dirichlet_phi[c] ~ nico_dirichlet(prior_dirichlet_alpha[c]);
           }
        }
        ////
        //// ---- Induced-dirichlet between study model (all studies per class in ONE ragged node):
        ////
        for (c in 1:2) {
             target += nico_induced_dirichlet_prior_of_ordered_cutpoints_ragged_lupdf(
                             C[c] | alpha[c], n_reported_cutpoints, k_of_each_C,
                             use_probit_link, baseline_location_prior[c]);
        }
        ////
        target += nico_normal_lupdf(to_vector(to_matrix(beta_z)) | 0.0, 1.0);
        target += sum(total_log_J_of_C);
        ////
        //// ---- Log-likelihood (all classes, studies & reported cutpoints directly from Xbeta and C in ONE node):
        ////
        target += nico_binomial_counts_above_successive_cutpoints_of_a_cumulative_link_lupmf(
                        x_2 | n, Xbeta, C, position_of_each_C, n_obs_cutpoints,
                        K_fold_CV_indicator, use_probit_link, 1e-10);
}


generated quantities {
    real simd_lanes = bmvp_simd_lanes();
      cov_matrix[2]  beta_Sigma = multiply_lower_tri_self_transpose(beta_L_Sigma);
      corr_matrix[2] beta_Omega = multiply_lower_tri_self_transpose(beta_L_Omega);
      ////
      array[2] vector<lower=0.0>[n_cat] dirichlet_cat_SDs_sigma;
      for (c in 1:2) {
          real alpha_0 = sum(alpha[c]);
          dirichlet_cat_SDs_sigma[c] = sqrt((alpha[c] .* (alpha_0 - alpha[c]) ) /
                                       (square(alpha_0) * (alpha_0 + 1.0)));
      }
      ////
      array[2] vector[n_thr] C_mu;
      for (c in 1:2) {
         real alpha_0 = sum(alpha[c]);
         vector[n_thr] P = cumulative_sum(alpha[c][1:n_thr]) / alpha_0;
         vector[n_thr] P_upper = reverse(cumulative_sum(reverse(alpha[c][2:n_cat]))) / alpha_0;
         real baseline_location_prior_of_class_c =
                    (c == 1) ? dot_product(baseline_case_nd, to_vector(prior_beta_mu_mean[1, 1:n_covariates_nd])) :
                               dot_product(baseline_case_d,  to_vector(prior_beta_mu_mean[2, 1:n_covariates_d]));
         C_mu[c] = baseline_location_prior_of_class_c +
                   nico_dirichlet_cutpoint_medians_and_sum_of_log_J(P, P_upper, alpha_0, use_probit_link)[1:n_thr];
      }
      ////
      real Xbeta_baseline_nd = dot_product(baseline_case_nd, to_vector(beta_mu[1, 1:n_covariates_nd]));
      real Xbeta_baseline_d  = dot_product(baseline_case_d,  to_vector(beta_mu[2, 1:n_covariates_d]));
      ////
      vector[n_thr] Fp_baseline = (use_probit_link == 1) ?
                                  nico_Phi( -(C_mu[1] - Xbeta_baseline_nd))  :
                                  nico_inv_logit(-(C_mu[1] - Xbeta_baseline_nd));
      vector[n_thr] Sp_baseline = 1.0 - Fp_baseline;
      vector[n_thr] Se_baseline = (use_probit_link == 1) ?
                                  nico_Phi( -(C_mu[2] - Xbeta_baseline_d))   :
                                  nico_inv_logit(-(C_mu[2] - Xbeta_baseline_d));
      ////
      array[2] vector[n_thr] C_pred;
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
      vector[n_studies] log_lik_study = rep_vector(0.0, n_studies);
      if (log_lik_indicator == 1) {
          array[2] matrix[n_studies, n_thr] latent_surv = map_latent_surv_to_the_reported_C(
                                                          C, n_thr, Xbeta, n_studies,
                                                          n_obs_cutpoints, position_of_each_C);
          array[n_studies] int dummy_CV_indicator;
          for (s in 1:n_studies) {
              dummy_CV_indicator[s] = 1;
          }
          array[3, 2] matrix[n_studies, n_thr] outs = nico_compute_log_lik_binomial_fact_data(
                                                      latent_surv, use_probit_link, n_thr, n_studies,
                                                      x_2, n, N_total, n_obs_cutpoints,
                                                      dummy_CV_indicator);
          for (s in 1:n_studies) {
              for (k in 1:n_obs_cutpoints[s]) {
                log_lik_study[s] += outs[1][1][s, k];
                log_lik_study[s] += outs[1][2][s, k];
              }
          }
      }
}
