//// DTA_MA_ord_biv_RANDthr_kappa_popmed_Se_Sp_prior_fused_Stan_math_only.stan
//// The fused model of MetaOrdDTA on Stan math only: the same model as
//// DTA_MA_ord_biv_RANDthr_kappa_popmed_Se_Sp_prior_NicoStan_fused_AVX.stan,
//// with Stan's built-in functions in place of its elementary functions (log, exp, Phi, normal_lpdf, ...)
//// and the plain Stan math include files. Only the fused functions stay external (declared below with no
//// body): their C++ is in inst/include/fused_Stan_math_only (Stan math only: no SIMD, any CPU), given to
//// cmdstanr as its user header with stanc's allow-undefined option.



functions {
    real nico_binomial_counts_above_successive_cutpoints_of_a_cumulative_link_lpmf(
          array[] int n, array[] int n_previous, vector latent, data vector is_first, data int link,
          data real bound);
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
        //// ---- The fused functions (their C++ in inst/include/fused_Stan_math_only, Stan math only):
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
        vector nico_pinkney_C_and_log_J_from_alpha_precomputed(
                vector alpha,
                vector C_raw_c,
                real baseline_location_prior,
                data vector Pinkney_w_c,
                data vector logit_start,
                data vector C_raw_lower_c,
                data vector C_raw_width_c,
                data real sum_log_width_c,
                data array[] int pair_of_each_C,
                data array[] int is_increment_of_each_C,
                data array[] int u_k,
                data array[] int u_kp,
                data vector u_lb,
                data vector u_wd,
                data int use_probit_link,
                data int softplus);
        vector nico_dirichlet_cutpoint_medians_and_log_J_at_cuts(
                vector alpha_c,
                real kappa_c,
                data array[] int cuts_c,
                data int use_probit_link);
        real nico_log_Set0_marginal_at_chosen_cutoffs_fast(
                vector x_chosen,
                data vector alpha_of_the_segments,
                data real v0,
                data int use_probit_link);
        real nico_ord_biv_randthr_ind_dir_and_binomial_lpmf(
                array[,,] int n,
                array[,,] int n_previous,
                matrix locations,
                array[] vector C,
                array[] vector alpha,
                data array[,] int position_of_each_C,
                data array[] int n_obs_cutpoints,
                data array[] int n_reported_cutpoints,
                data array[] int id_u_seg_lo,
                data array[] int id_u_seg_hi,
                data array[] int id_seg_uid,
                data array[] int id_seg_own_of_C,
                data array[] int id_seg_next_of_C,
                data array[] int K_fold_CV_indicator,
                data int link,
                data vector baseline_location_prior,
                data real bound);
        vector nico_pinkney_C_and_log_J_both_from_alpha_precomputed(
                array[] vector alpha,
                array[] vector C_raw,
                data vector baseline_location_prior,
                data array[] vector Pinkney_w,
                data array[] vector logit_start,
                data array[] vector C_raw_lower,
                data array[] vector C_raw_width,
                data vector sum_log_width,
                data array[] int pair_of_each_C,
                data array[] int is_increment_of_each_C,
                data array[] int u_k,
                data array[] int u_kp,
                data vector u_lb,
                data vector u_wd,
                data int use_probit_link,
                data int softplus);
        real nico_ord_biv_rc_priors_lpdf(
                array[] vector dirichlet_phi,
                matrix beta_mu,
                row_vector beta_SD,
                matrix beta_L_Omega,
                array[] row_vector beta_z,
                array[] vector alpha,
                vector log_kappa,
                vector kappa,
                vector total_log_J_of_C,
                data array[] int n_cov,
                data matrix prior_beta_mu_mean,
                data matrix prior_beta_mu_SD,
                data vector prior_beta_SD_mean,
                data vector prior_beta_SD_SD,
                data real prior_beta_corr_LKJ,
                data vector prior_kappa_dof,
                data vector prior_kappa_mean,
                data vector prior_kappa_SD,
                data array[] vector prior_dirichlet_alpha,
                data array[] vector alpha_of_the_segments,
                data array[,] int chosen_cutoffs,
                data array[] int n_chosen_cutoffs,
                data int use_Se_Sp_prior,
                data int Se_Sp_prior_type,
                data array[] vector prior_Se_Sp_dirichlet_A,
                data array[] matrix prior_beta_td,
                data array[] matrix prior_link_mean_SD_td,
                data vector baseline_location_prior_td,
                data vector v0_baseline_prior_td,
                data vector baseline_case_nd,
                data vector baseline_case_d,
                data int use_probit_link);
        real nico_ord_biv_randthr_ind_dir_and_binomial_precomputed_lpmf(
                data array[] int bin_n,
                data array[] int bin_n_prev,
                data array[] int bin_sc,
                data array[] int bin_Cpos,
                data vector bin_is_first,
                data int K_bin,
                data int bin_is_seq,
                matrix Xbeta,
                array[] vector C,
                array[] vector alpha,
                data array[] int id_u_seg_lo,
                data array[] int id_u_seg_hi,
                data array[] int id_u_seg_cnt,
                data array[] int id_seg_uid,
                data array[] int id_seg_own_of_C,
                data array[] int id_seg_next_of_C,
                data array[] int id_int_p,
                data array[] int id_seg_p,
                data array[] int id_seg_w,
                data int id_n_interior,
                data int R_nonempty,
                data int link,
                data vector baseline_location_prior,
                data real bound);
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
        //// ---- Precomputed width, Pinkney pair tables, Induced-Dirichlet unique segments, and prior constants:
        array[2] vector[n_C] C_raw_width_of_each_C;
        vector[2] sum_log_width_of_each_C;
        for (c in 1:2) {
              C_raw_width_of_each_C[c] = C_raw_upper_of_each_C[c] - C_raw_lower_of_each_C[c];
              sum_log_width_of_each_C[c] = sum(log(C_raw_width_of_each_C[c]));
        }
        int n_unique_C_pairs = 0;
        int max_C_pairs = n_thr * (n_thr + 1);
        array[n_C] int pair_of_each_C;
        array[n_C] int is_increment_of_each_C;
        array[max_C_pairs] int raw_u_k = rep_array(0, max_C_pairs);
        array[max_C_pairs] int raw_u_kp = rep_array(0, max_C_pairs);
        vector[max_C_pairs] raw_u_lb = rep_vector(0.0, max_C_pairs);
        vector[max_C_pairs] raw_u_wd = rep_vector(1.0, max_C_pairs);
        {
              matrix[n_thr, n_thr + 1] pair_lookup = rep_matrix(0.0, n_thr, n_thr + 1);
              for (p in 1:n_C) {
                    int k = k_of_each_C[p];
                    int kp = k_previous_of_each_C[p];
                    is_increment_of_each_C[p] = (kp > 0) ? 1 : 0;
                    int uid = to_int(pair_lookup[k, kp + 1]);
                    if (uid == 0) {
                          n_unique_C_pairs += 1;
                          uid = n_unique_C_pairs;
                          pair_lookup[k, kp + 1] = uid;
                          raw_u_k[uid] = k;
                          raw_u_kp[uid] = kp;
                          raw_u_lb[uid] = C_raw_lower_of_each_C[1][p];
                          raw_u_wd[uid] = C_raw_width_of_each_C[1][p];
                    }
                    pair_of_each_C[p] = uid;
              }
        }
        int n_unique_C_pairs_kept = (n_unique_C_pairs > 0) ? n_unique_C_pairs : 1;
        array[n_unique_C_pairs_kept] int u_k_of_each_pair = raw_u_k[1:n_unique_C_pairs_kept];
        array[n_unique_C_pairs_kept] int u_kp_of_each_pair = raw_u_kp[1:n_unique_C_pairs_kept];
        vector[n_unique_C_pairs_kept] u_lb_of_each_pair = raw_u_lb[1:n_unique_C_pairs_kept];
        vector[n_unique_C_pairs_kept] u_wd_of_each_pair = raw_u_wd[1:n_unique_C_pairs_kept];

        int n_id_total_segs = n_C + n_studies_with_a_cutpoint;
        int n_id_total_segs_kept = (n_id_total_segs > 0) ? n_id_total_segs : 1;
        int n_id_u_segs = 0;
        array[n_id_total_segs_kept] int raw_id_u_lo = rep_array(0, n_id_total_segs_kept);
        array[n_id_total_segs_kept] int raw_id_u_hi = rep_array(1, n_id_total_segs_kept);
        array[n_id_total_segs_kept] int id_seg_uid = rep_array(0, n_id_total_segs_kept);
        array[n_C] int id_seg_own_of_C;
        array[n_C] int id_seg_next_of_C;
        {
              matrix[n_cat + 1, n_cat + 1] seg_lookup = rep_matrix(0.0, n_cat + 1, n_cat + 1);
              int off = 0;
              int s_idx = 0;
              for (r in 1:n_studies) {
                    int Kr = n_reported_cutpoints[r];
                    if (Kr > 0) {
                          for (s in 0:Kr) {
                                int lo = (s == 0) ? 0 : k_of_each_C[off + s];
                                int hi = (s == Kr) ? n_cat : k_of_each_C[off + s + 1];
                                if (s < Kr) {
                                      id_seg_own_of_C[off + s + 1] = s_idx;
                                      id_seg_next_of_C[off + s + 1] = s_idx + 1;
                                }
                                int uid = to_int(seg_lookup[lo + 1, hi + 1]);
                                if (uid == 0) {
                                      n_id_u_segs += 1;
                                      uid = n_id_u_segs;
                                      seg_lookup[lo + 1, hi + 1] = uid;
                                      raw_id_u_lo[uid] = lo;
                                      raw_id_u_hi[uid] = hi;
                                }
                                s_idx += 1;
                                id_seg_uid[s_idx] = uid - 1;
                          }
                          off += Kr;
                    }
              }
        }
        int n_id_u_segs_kept = (n_id_u_segs > 0) ? n_id_u_segs : 1;
        array[n_id_u_segs_kept] int id_u_seg_lo = raw_id_u_lo[1:n_id_u_segs_kept];
        array[n_id_u_segs_kept] int id_u_seg_hi = raw_id_u_hi[1:n_id_u_segs_kept];

        array[2] int n_cov_td = {n_covariates_nd, n_covariates_d};
        vector[2] baseline_location_prior_td;
        baseline_location_prior_td[1] = dot_product(baseline_case_nd, to_vector(prior_beta_mu_mean[1, 1:n_covariates_nd]));
        baseline_location_prior_td[2] = dot_product(baseline_case_d,  to_vector(prior_beta_mu_mean[2, 1:n_covariates_d]));
        //// ---- Precomputed flat Induced-Dirichlet and Binomial indices:
        array[n_id_u_segs_kept] int id_u_seg_cnt = rep_array(0, n_id_u_segs_kept);
        for (si in 1:n_id_total_segs) {
              id_u_seg_cnt[id_seg_uid[si] + 1] += 1;
        }
        int id_n_interior_td = (n_C > n_studies_with_a_cutpoint) ? (n_C - n_studies_with_a_cutpoint) : 0;
        int id_n_interior_kept = (id_n_interior_td > 0) ? id_n_interior_td : 1;
        array[id_n_interior_kept] int id_int_p_td = rep_array(0, id_n_interior_kept);
        array[n_id_total_segs_kept] int id_seg_p_td = rep_array(0, n_id_total_segs_kept);
        array[n_id_total_segs_kept] int id_seg_w_td = rep_array(0, n_id_total_segs_kept);
        {
              int off = 0;
              int s_idx = 0;
              int int_idx = 0;
              for (r in 1:n_studies) {
                    int Kr = n_reported_cutpoints[r];
                    if (Kr > 0) {
                          for (m in 1:(Kr - 1)) {
                                int_idx += 1;
                                id_int_p_td[int_idx] = off + m;
                          }
                          for (s in 0:Kr) {
                                s_idx += 1;
                                if (s == 0) {
                                      id_seg_p_td[s_idx] = off;
                                      id_seg_w_td[s_idx] = -1;
                                } else if (s == Kr) {
                                      id_seg_p_td[s_idx] = off + Kr - 1;
                                      id_seg_w_td[s_idx] = -2;
                                } else {
                                      id_seg_p_td[s_idx] = off + s;
                                      id_seg_w_td[s_idx] = int_idx - (Kr - 1) + s - 1;
                                }
                          }
                          off += Kr;
                    }
              }
        }
        int K_bin_td = 0;
        for (s in 1:n_studies) {
              if (K_fold_CV_indicator[s] == 1) K_bin_td += 2 * n_obs_cutpoints[s];
        }
        int K_bin_kept = (K_bin_td > 0) ? K_bin_td : 1;
        array[K_bin_kept] int bin_n_td = rep_array(0, K_bin_kept);
        array[K_bin_kept] int bin_n_prev_td = rep_array(0, K_bin_kept);
        array[K_bin_kept] int bin_sc_td = rep_array(0, K_bin_kept);
        array[K_bin_kept] int bin_Cpos_td = rep_array(0, K_bin_kept);
        vector[K_bin_kept] bin_is_first_td = rep_vector(0.0, K_bin_kept);
        int bin_is_seq_td = (K_bin_td == 2 * n_C) ? 1 : 0;
        {
              int t = 0;
              for (c in 1:2) {
                    for (s in 1:n_studies) {
                          if (K_fold_CV_indicator[s] == 1) {
                                for (cut_i in 1:n_obs_cutpoints[s]) {
                                      int p0 = position_of_each_C[s, cut_i] - 1;
                                      if (p0 < 0 || p0 >= n_C) reject("invalid position_of_each_C");
                                      t += 1;
                                      bin_n_td[t] = x_2[c, s, cut_i];
                                      bin_n_prev_td[t] = n[c, s, cut_i];
                                      bin_sc_td[t] = (c - 1) * n_studies + (s - 1);
                                      bin_Cpos_td[t] = (c - 1) * n_C + p0;
                                      bin_is_first_td[t] = (cut_i == 1) ? 1.0 : 0.0;
                                      if (bin_is_seq_td == 1 && bin_Cpos_td[t] != (t - 1)) bin_is_seq_td = 0;
                                }
                          }
                    }
              }
        }
        vector[2] v0_baseline_prior_td;
        v0_baseline_prior_td[1] = dot_self(baseline_case_nd .* to_vector(prior_beta_mu_SD[1, 1:n_covariates_nd]));
        v0_baseline_prior_td[2] = dot_self(baseline_case_d  .* to_vector(prior_beta_mu_SD[2, 1:n_covariates_d]));
        array[2] matrix[n_chosen_cutoffs_max, 2] prior_beta_td = prior_Se_Sp_beta;
        array[2] matrix[n_chosen_cutoffs_max, 2] prior_link_mean_SD_td = prior_Se_Sp_link_mean_SD;
        if (n_chosen_cutoffs[2] > 0) {
              prior_link_mean_SD_td[2][1:n_chosen_cutoffs[2], 1] =
                    -prior_link_mean_SD_td[2][1:n_chosen_cutoffs[2], 1];
              if (Se_Sp_prior_type == 2) {
                    prior_beta_td[2][1:n_chosen_cutoffs[2], ] =
                          append_col(prior_Se_Sp_beta[2][1:n_chosen_cutoffs[2], 2],
                                     prior_Se_Sp_beta[2][1:n_chosen_cutoffs[2], 1]);
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
        vector[2] log_kappa = log(kappa);
        array[2] vector<lower=0.0>[n_cat] alpha;
        for (c in 1:2) {
              alpha[c] = 0.01 + dirichlet_phi[c] * kappa[c];
        }
        //// ---- Fused dual-class Pinkney partially-centered study cutpoints and Jacobian directly from alpha:
        array[2] vector[n_C] C;
        vector[2] total_log_J_of_C;
        {
              vector[2 * n_C + 2] C_and_J_both = nico_pinkney_C_and_log_J_both_from_alpha_precomputed(
                    alpha,
                    C_raw,
                    baseline_location_prior_td,
                    Pinkney_w,
                    Pinkney_eta0,
                    C_raw_lower_of_each_C,
                    C_raw_width_of_each_C,
                    sum_log_width_of_each_C,
                    pair_of_each_C,
                    is_increment_of_each_C,
                    u_k_of_each_pair,
                    u_kp_of_each_pair,
                    u_lb_of_each_pair,
                    u_wd_of_each_pair,
                    use_probit_link,
                    softplus);
              C[1] = C_and_J_both[1:n_C];
              C[2] = C_and_J_both[(n_C + 1):(2 * n_C)];
              total_log_J_of_C[1] = C_and_J_both[2 * n_C + 1];
              total_log_J_of_C[2] = C_and_J_both[2 * n_C + 2];
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
        //// ---- Fused model priors + Dirichlet / Se-Sp / Set-0 marginal / cutpoint-median priors + Pinkney log-det-J:
        target += nico_ord_biv_rc_priors_lupdf(
                        dirichlet_phi | beta_mu, beta_SD, beta_L_Omega, beta_z,
                        alpha, log_kappa, kappa, total_log_J_of_C,
                        n_cov_td, prior_beta_mu_mean, prior_beta_mu_SD,
                        prior_beta_SD_mean, prior_beta_SD_SD, prior_beta_corr_LKJ,
                        prior_kappa_dof, prior_kappa_mean, prior_kappa_SD,
                        prior_dirichlet_alpha, alpha_of_the_segments,
                        chosen_cutoffs, n_chosen_cutoffs,
                        use_Se_Sp_prior, Se_Sp_prior_type,
                        prior_Se_Sp_dirichlet_A, prior_beta_td, prior_link_mean_SD_td,
                        baseline_location_prior_td, v0_baseline_prior_td,
                        baseline_case_nd, baseline_case_d, use_probit_link);
        ////
        //// ---- Fused Induced-Dirichlet between-study model + Binomial cumulative-link likelihood (all studies/classes in ONE node):
        ////
        target += nico_ord_biv_randthr_ind_dir_and_binomial_precomputed_lupmf(
                        bin_n_td | bin_n_prev_td, bin_sc_td, bin_Cpos_td, bin_is_first_td,
                        K_bin_td, bin_is_seq_td, Xbeta, C, alpha,
                        id_u_seg_lo, id_u_seg_hi, id_u_seg_cnt, id_seg_uid,
                        id_seg_own_of_C, id_seg_next_of_C, id_int_p_td, id_seg_p_td, id_seg_w_td,
                        id_n_interior_td, n_studies_with_a_cutpoint, use_probit_link,
                        baseline_location_prior_td, 1e-10);
}


generated quantities {
    // real simd_lanes = bmvp_simd_lanes();  //// (the number of SIMD lanes of the fused version: none here)
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
                                  Phi( -(C_mu[1] - Xbeta_baseline_nd))  :
                                  inv_logit(-(C_mu[1] - Xbeta_baseline_nd));
      vector[n_thr] Sp_baseline = 1.0 - Fp_baseline;
      vector[n_thr] Se_baseline = (use_probit_link == 1) ?
                                  Phi( -(C_mu[2] - Xbeta_baseline_d))   :
                                  inv_logit(-(C_mu[2] - Xbeta_baseline_d));
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
              C_pred[c][k] = C_mu[c][k] + normal_rng(0.0, sigma_C_k);
            }
        }
      }
      ////
      vector[2] beta_random_pred = multi_normal_cholesky_rng(rep_vector(0.0, 2), beta_L_Sigma);
      real Xbeta_baseline_pred_nd = dot_product(baseline_case_nd, beta_mu[1, 1:n_covariates_nd]) + beta_random_pred[1];
      real Xbeta_baseline_pred_d  = dot_product(baseline_case_d,  beta_mu[2, 1:n_covariates_d])  + beta_random_pred[2];
      ////
      vector[n_thr] Fp_baseline_pred = (use_probit_link == 1) ?
                                       Phi(-(C_pred[1] - Xbeta_baseline_pred_nd)) :
                                       inv_logit(-(C_pred[1] - Xbeta_baseline_pred_nd));
      vector[n_thr] Sp_baseline_pred = 1.0 - Fp_baseline_pred;
      vector[n_thr] Se_baseline_pred = (use_probit_link == 1) ?
                                       Phi(-(C_pred[2] - Xbeta_baseline_pred_d))  :
                                       inv_logit(-(C_pred[2] - Xbeta_baseline_pred_d));
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
