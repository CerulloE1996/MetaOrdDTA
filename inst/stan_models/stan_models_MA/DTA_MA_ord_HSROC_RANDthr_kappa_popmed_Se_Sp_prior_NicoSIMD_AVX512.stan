// NicoSIMD_AVX512 version of DTA_MA_ord_HSROC_RANDthr_kappa_popmed_Se_Sp_prior.stan written by NicoStan (source md5 2a7720b0d5d69ef8916972d7d403dc0a)
// NicoSIMD functions used: exp, log, log1p_exp, inv_logit, log_inv_logit, Phi, inv_Phi, logit, log1m_inv_logit, expm1, normal_rng, normal, dirichlet, student_t
//// DTA_MA_ord_HSROC_RANDthr_kappa_popmed_Se_Sp_prior_Stan_math_only.stan
//// The same data, parameters and target as the fused AVX file of this model with the Dirichlet prior of dirichlet_phi, on Stan::math only.
////

functions {
    real bmvp_simd_lanes();
    real nico_exp(real x);
    real nico_log(real x);
    real nico_log1p_exp(real x);
    real nico_inv_logit(real x);
    real nico_log_inv_logit(real x);
    real nico_Phi(real x);
    real nico_inv_Phi(real x);
    real nico_logit(real x);
    real nico_log1m_inv_logit(real x);
    real nico_expm1(real x);
    vector nico_exp(vector x);
    vector nico_log(vector x);
    vector nico_log1p_exp(vector x);
    vector nico_inv_logit(vector x);
    vector nico_log_inv_logit(vector x);
    vector nico_Phi(vector x);
    vector nico_inv_Phi(vector x);
    vector nico_logit(vector x);
    vector nico_log1m_inv_logit(vector x);
    vector nico_expm1(vector x);
    row_vector nico_exp(row_vector x);
    row_vector nico_log(row_vector x);
    row_vector nico_log1p_exp(row_vector x);
    row_vector nico_inv_logit(row_vector x);
    row_vector nico_log_inv_logit(row_vector x);
    row_vector nico_Phi(row_vector x);
    row_vector nico_inv_Phi(row_vector x);
    row_vector nico_logit(row_vector x);
    row_vector nico_log1m_inv_logit(row_vector x);
    row_vector nico_expm1(row_vector x);
    matrix nico_exp(matrix x);
    matrix nico_log(matrix x);
    matrix nico_log1p_exp(matrix x);
    matrix nico_inv_logit(matrix x);
    matrix nico_log_inv_logit(matrix x);
    matrix nico_Phi(matrix x);
    matrix nico_inv_Phi(matrix x);
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
    real nico_dirichlet_lpdf(vector theta, vector alpha);
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
        vector[n_C] C_raw_lower_of_each_C;
        vector[n_C] C_raw_upper_of_each_C;
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
                                   C_raw_lower_of_each_C[p] = C_raw_lower + (k - 1) * step_lower;
                                   C_raw_upper_of_each_C[p] = C_raw_upper + (k - 1) * step_upper;
                              } else {
                                   i += 1;
                                   position_of_each_increment[i] = p;
                                   C_raw_lower_of_each_C[p] = (softplus == 1) ? nico_log(nico_expm1((k - k_previous) * step_lower)) :
                                                                                 C_raw_lower + nico_log(k - k_previous);
                                   C_raw_upper_of_each_C[p] = (softplus == 1) ? nico_log(nico_expm1((k - k_previous) * step_upper)) :
                                                                                 C_raw_upper + nico_log(k - k_previous);
                              }
                              k_previous = k;
                         }
                    }
                    upper_index_of_each_segment[p + s] = n_cat + 1;
                    lower_index_of_each_segment[p + s] = k_previous + 1;
              }
        }
        vector[n_C] C_raw_width_of_each_C = C_raw_upper_of_each_C - C_raw_lower_of_each_C;
        real sum_log_width_of_each_C = sum(nico_log(C_raw_width_of_each_C));
        array[n_C] int k_of_each_C;
        array[n_C] int k_previous_of_each_C;
        array[n_C] int is_increment_of_each_C = rep_array(0, n_C);
        for (i in 1:(n_C - n_studies_with_a_cutpoint)) {
              is_increment_of_each_C[position_of_each_increment[i]] = 1;
        }
        int n_unique_C_pairs = 0;
        int n_first_C_pairs = 0;
        int n_step_C_pairs = 0;
        int max_C_pairs = n_thr * (n_thr + 1);
        vector[max_C_pairs] lower_of_each_unique_pair = rep_vector(0.0, max_C_pairs);
        vector[max_C_pairs] width_of_each_unique_pair = rep_vector(1.0, max_C_pairs);
        array[max_C_pairs] int pos_first_pair = rep_array(0, max_C_pairs);
        array[max_C_pairs] int k_of_first_pair = rep_array(0, max_C_pairs);
        array[max_C_pairs] int pos_step_pair = rep_array(0, max_C_pairs);
        array[max_C_pairs] int k_of_step_pair = rep_array(0, max_C_pairs);
        array[max_C_pairs] int kprev_of_step_pair = rep_array(0, max_C_pairs);
        array[n_C] int pair_of_each_C;
        real start_of_C_raw = -2.0;
        vector[n_C] logit_of_the_place_of_the_start;
        {
              matrix[n_thr, n_thr + 1] pair_lookup = rep_matrix(0.0, n_thr, n_thr + 1);
              int p_idx = 0;
              for (s in 1:n_studies) {
                    int k_prev = 0;
                    for (cut_i in 1:n_obs_cutpoints[s]) {
                          int k = cutpoint_index[1, s, cut_i + 1];
                          if (k != -1) {
                                p_idx += 1;
                                k_of_each_C[p_idx] = k;
                                k_previous_of_each_C[p_idx] = k_prev;
                                real step_start = (softplus == 1) ? nico_log1p_exp(start_of_C_raw) : nico_exp(start_of_C_raw);
                                real start_this = (k_prev == 0) ?
                                                  (start_of_C_raw + (k - 1) * step_start) :
                                                  ((softplus == 1) ? nico_log(nico_expm1((k - k_prev) * step_start)) : (start_of_C_raw + nico_log(k - k_prev)));
                                real place_start = (start_this - C_raw_lower_of_each_C[p_idx]) / C_raw_width_of_each_C[p_idx];
                                logit_of_the_place_of_the_start[p_idx] = nico_logit(fmin(fmax(place_start, 0.001), 0.999));
                                int u = to_int(pair_lookup[k, k_prev + 1]);
                                if (u == 0) {
                                      n_unique_C_pairs += 1;
                                      u = n_unique_C_pairs;
                                      pair_lookup[k, k_prev + 1] = u;
                                      lower_of_each_unique_pair[u] = C_raw_lower_of_each_C[p_idx];
                                      width_of_each_unique_pair[u] = C_raw_width_of_each_C[p_idx];
                                      if (k_prev == 0) {
                                            n_first_C_pairs += 1;
                                            pos_first_pair[n_first_C_pairs] = u;
                                            k_of_first_pair[n_first_C_pairs] = k;
                                      } else {
                                            n_step_C_pairs += 1;
                                            pos_step_pair[n_step_C_pairs] = u;
                                            k_of_step_pair[n_step_C_pairs] = k;
                                            kprev_of_step_pair[n_step_C_pairs] = k_prev;
                                      }
                                }
                                pair_of_each_C[p_idx] = u;
                                k_prev = k;
                          }
                    }
              }
        }
        //// ---- Data-only partial-centering weights Pinkney_w in [0, 1] and empirical origins Pinkney_eta0
        ////      (Sean Pinkney, partially_centered_example_Sean_Pinkney.stan & spinkney/hierarchical_model_tutorial
        ////       stan/hierarchical_pcp.stan, R/1_hier_code.R: shared cutpoints across c = 1, 2 combine Fisher
        ////       precisions J_1 + J_2 from x and K_fold_CV_indicator only, never from parameters):
        vector<lower=0.0, upper=1.0>[n_C] Pinkney_w;
        vector[n_C] Pinkney_eta0;
        {
              vector[n_C] J_tot_of_each_C = rep_vector(0.0, n_C);
              vector[n_C] J_eta_num_of_each_C = rep_vector(0.0, n_C);
              array[n_C] int active_of_each_C = rep_array(0, n_C);
              array[max_C_pairs] int n_studies_of_pair = rep_array(0, max_C_pairs);
              vector[max_C_pairs] W1_of_pair  = rep_vector(0.0, max_C_pairs);
              vector[max_C_pairs] W2_of_pair  = rep_vector(0.0, max_C_pairs);
              vector[max_C_pairs] WY_of_pair  = rep_vector(0.0, max_C_pairs);
              vector[max_C_pairs] WY2_of_pair = rep_vector(0.0, max_C_pairs);
              int p_idx = 0;
              for (s in 1:n_studies) {
                    int N1 = N_total[1, s];
                    int N2 = N_total[2, s];
                    int k_prev = 0;
                    array[2] real p_prev = {0.0, 0.0};
                    array[2] real z_prev = {0.0, 0.0};
                    array[2] real f_prev = {1.0, 1.0};
                    array[2] int x_ab_prev = {N1, N2};
                    for (cut_i in 1:n_obs_cutpoints[s]) {
                          int k = cutpoint_index[1, s, cut_i + 1];
                          if (k != -1) {
                                p_idx += 1;
                                real lb_p = C_raw_lower_of_each_C[p_idx];
                                real wd_p = C_raw_width_of_each_C[p_idx];
                                int u = pair_of_each_C[p_idx];
                                if (K_fold_CV_indicator[s] == 0 || (N1 <= 0 && N2 <= 0)) {
                                      active_of_each_C[p_idx] = 0;
                                      Pinkney_eta0[p_idx] = logit_of_the_place_of_the_start[p_idx];
                                } else {
                                      real gap_num = 0.0;
                                      real gap_den = 0.0;
                                      int has_pos_cell = 0;
                                      for (c in 1:2) {
                                            int N_cs = N_total[c, s];
                                            if (N_cs > 0) {
                                                  active_of_each_C[p_idx] = 1;
                                                  int x_ab = x_2[c, s, cut_i];
                                                  real p_raw = (N_cs - x_ab + 0.5) / (N_cs + 1.0);
                                                  real p_below = (k_prev > 0 && p_raw <= p_prev[c]) ?
                                                                 fmin(0.999999, p_prev[c] + 0.5 / ((N_cs + 1.0) * (n_thr + 1.0))) :
                                                                 p_raw;
                                                  real z_hat = (use_probit_link == 1) ? nico_inv_Phi(p_below) : nico_logit(p_below);
                                                  real fz = (use_probit_link == 1) ?
                                                            (nico_exp(-0.5 * square(z_hat)) / sqrt(2.0 * pi())) :
                                                            (p_below * (1.0 - p_below));
                                                  real vf_curr = p_below * (1.0 - p_below) / square(fz);
                                                  if (k_prev == 0) {
                                                        real q_hat = fmin(fmax((z_hat - lb_p) / wd_p, 0.005), 0.995);
                                                        real J_c = N_cs / fmax(vf_curr, 1e-12);
                                                        J_tot_of_each_C[p_idx] += J_c;
                                                        J_eta_num_of_each_C[p_idx] += J_c * nico_logit(q_hat);
                                                  } else {
                                                        real vf_prev = p_prev[c] * (1.0 - p_prev[c]) / square(f_prev[c]);
                                                        real cov_term = p_prev[c] * (1.0 - p_below) / (f_prev[c] * fz);
                                                        real var_gap = fmax((vf_curr + vf_prev - 2.0 * cov_term) / N_cs, 1e-10);
                                                        real gap = fmax(z_hat - z_prev[c], 1e-3);
                                                        real J_gap = 1.0 / var_gap;
                                                        real r_hat = (softplus == 1) ? nico_log(fmax(nico_expm1(gap), 1e-3)) : nico_log(gap);
                                                        real q_hat = fmin(fmax((r_hat - lb_p) / wd_p, 0.005), 0.995);
                                                        J_tot_of_each_C[p_idx] += J_gap;
                                                        J_eta_num_of_each_C[p_idx] += J_gap * nico_logit(q_hat);
                                                        if ((x_ab_prev[c] - x_ab) >= 1) {
                                                              has_pos_cell = 1;
                                                              gap_num += J_gap * gap;
                                                              gap_den += J_gap;
                                                        }
                                                  }
                                                  p_prev[c] = p_below;
                                                  z_prev[c] = z_hat;
                                                  f_prev[c] = fz;
                                                  x_ab_prev[c] = x_ab;
                                            }
                                      }
                                      Pinkney_eta0[p_idx] = (J_tot_of_each_C[p_idx] > 0.0) ?
                                                            (J_eta_num_of_each_C[p_idx] / J_tot_of_each_C[p_idx]) :
                                                            logit_of_the_place_of_the_start[p_idx];
                                      if (k_prev > 0 && has_pos_cell == 1 && gap_den > 0.0) {
                                            real gap_bar = gap_num / gap_den;
                                            n_studies_of_pair[u] += 1;
                                            W1_of_pair[u]  += gap_den;
                                            W2_of_pair[u]  += square(gap_den);
                                            WY_of_pair[u]  += gap_den * gap_bar;
                                            WY2_of_pair[u] += gap_den * square(gap_bar);
                                      }
                                }
                                k_prev = k;
                          }
                    }
              }
              real Q_tot = 0.0;
              real df_tot = 0.0;
              real C_tot = 0.0;
              for (u in 1:n_unique_C_pairs) {
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
                    Pinkney_w[p] = (active_of_each_C[p] == 0) ? 0.0 : (1.0 / (1.0 + J_tot_of_each_C[p] * tau2_pos));
              }
        }
        ////
        //// ---- Class on the cutpoints (p_class):
        ////
        int p_class = (n_chosen_cutoffs[2] > n_chosen_cutoffs[1]) ? 2 : 1;
        ////
        //// ---- Unit probability priors at every cutoff, with exponential scales and the probit link:
        ////
        int joint_unit_probability_prior = use_Se_Sp_prior == 1 &&
                                             (Se_Sp_prior_type == 1 || Se_Sp_prior_type == 2) &&
                                             use_probit_link == 1 && softplus == 0 && n_thr > 1 &&
                                             n_chosen_cutoffs[1] == n_thr && n_chosen_cutoffs[2] == n_thr;
        if (joint_unit_probability_prior == 1) {
              for (c in 1:2) {
                    for (j in 1:n_thr) {
                          if (chosen_cutoffs[c, j] != j) joint_unit_probability_prior = 0;
                          if (Se_Sp_prior_type == 2 && (prior_Se_Sp_beta[c][j, 1] != 1.0 ||
                              prior_Se_Sp_beta[c][j, 2] != 1.0)) joint_unit_probability_prior = 0;
                    }
                    if (Se_Sp_prior_type == 1) {
                          for (j in 1:(n_thr + 1)) {
                                if (prior_Se_Sp_dirichlet_A[c][j] != 1.0)
                                      joint_unit_probability_prior = 0;
                          }
                    }
              }
        }
        ////
        //// ---- Joint link-scale priors at shared cutoffs:
        ////
        int joint_link_prior = use_Se_Sp_prior == 1 && Se_Sp_prior_type == 4 &&
                               n_chosen_cutoffs[1] >= 2 &&
                               n_chosen_cutoffs[1] == n_chosen_cutoffs[2];
        if (joint_link_prior == 1) {
              for (j in 1:n_chosen_cutoffs[1]) {
                    if (chosen_cutoffs[1, j] != chosen_cutoffs[2, j]) joint_link_prior = 0;
              }
              if (baseline_case[1] == 0.0) joint_link_prior = 0;
        }

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
        vector[n_C] C_raw; //// one free value per reported study cutpoint (Sean Pinkney partially-centered parameterization)
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

        //// (the reported cutpoints of each study only, study after study, as construct_C on each study's own)
        vector[n_C] C_raw_between_its_bounds;
        real log_J_of_the_free_values = 0.0;
        {
              real baseline_location_prior_of_C_raw = 0.5 * (mult_nd + mult_d) * dot_product(baseline_case, prior_beta_mu_mean);
              real alpha_0 = sum(alpha);
              vector[n_thr] P = cumulative_sum(alpha[1:n_thr]) / alpha_0;
              vector[n_thr] P_upper = reverse(cumulative_sum(reverse(alpha[2:n_cat]))) / alpha_0;
              vector[n_thr] centre_of_each_cutpoint;
              vector[n_thr] density_at_each_cutpoint;
              if (use_probit_link == 1) {
                    for (k in 1:n_thr) {
                          centre_of_each_cutpoint[k] = (P[k] <= 0.5) ? nico_inv_Phi(P[k]) : -nico_inv_Phi(P_upper[k]);
                    }
                    density_at_each_cutpoint = nico_exp(-0.5 * square(centre_of_each_cutpoint)) / sqrt(2.0 * pi());
              } else {
                    centre_of_each_cutpoint = nico_log(P) - nico_log(P_upper);
                    density_at_each_cutpoint = P .* P_upper;
              }
              for (p in 1:n_C) {
                    int k = k_of_each_C[p];
                    int k_previous = k_previous_of_each_C[p];
                    real centre;
                    real spread;
                    if (k_previous == 0) {
                          centre = centre_of_each_cutpoint[k] + baseline_location_prior_of_C_raw;
                          spread = sqrt(P[k] * P_upper[k] / (alpha_0 + 1.0)) / density_at_each_cutpoint[k];
                    } else {
                          real step_between_the_centres = centre_of_each_cutpoint[k] - centre_of_each_cutpoint[k_previous];
                          real variance_of_the_step =
                                (P[k] * P_upper[k] / square(density_at_each_cutpoint[k]) +
                                 P[k_previous] * P_upper[k_previous] / square(density_at_each_cutpoint[k_previous]) -
                                 2.0 * P[k_previous] * P_upper[k] / (density_at_each_cutpoint[k] * density_at_each_cutpoint[k_previous])) /
                                (alpha_0 + 1.0);
                          if (softplus == 1) {
                                real expm1_step = fmax(nico_expm1(step_between_the_centres), 1e-12);
                                centre = nico_log(expm1_step);
                                spread = sqrt(fmax(variance_of_the_step, 1e-14)) * nico_exp(step_between_the_centres) / expm1_step;
                          } else {
                                centre = nico_log(step_between_the_centres);
                                spread = sqrt(fmax(variance_of_the_step, 1e-14)) / step_between_the_centres;
                          }
                    }
                    real lb_this = C_raw_lower_of_each_C[p];
                    real width = C_raw_width_of_each_C[p];
                    real place_as_it_is = (centre - lb_this) / width;
                    real place_of_the_centre = 0.001 + (nico_log1p_exp(500.0 * (place_as_it_is - 0.001)) -
                                                        nico_log1p_exp(500.0 * (place_as_it_is - 0.999))) / 500.0;
                    real tau = spread / sqrt(square(width * place_of_the_centre * (1.0 - place_of_the_centre)) + square(spread));
                    real w = Pinkney_w[p];
                    real tau_wt = tau - w * (tau - 1.0);
                    real eta_shifted = (1.0 - w) * Pinkney_eta0[p] + C_raw[p];
                    real logit_of_the_place = (w * nico_logit(place_of_the_centre) + tau * eta_shifted) / tau_wt;
                    C_raw_between_its_bounds[p] = lb_this + width * nico_inv_logit(logit_of_the_place);
                    log_J_of_the_free_values += nico_log(width) + nico_log_inv_logit(logit_of_the_place) + nico_log1m_inv_logit(logit_of_the_place) +
                                                nico_log(tau) - nico_log(tau_wt);
              }
        }
        vector[n_C] C = C_raw_between_its_bounds;
        {
              vector[n_C - n_studies_with_a_cutpoint] increments_of_the_cutpoints =
                    (softplus == 1) ? nico_log1p_exp(C_raw_between_its_bounds[position_of_each_increment]) :
                                      nico_exp(C_raw_between_its_bounds[position_of_each_increment]);
              for (i in 1:(n_C - n_studies_with_a_cutpoint)) {
                    int p = position_of_each_increment[i];
                    C[p] = C[p - 1] + increments_of_the_cutpoints[i];
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
        beta_SD ~ nico_normal( prior_beta_SD_mean, prior_beta_SD_SD);
        raw_scale_SD ~ nico_normal( prior_raw_scale_SD_mean, prior_raw_scale_SD_SD);
        ////
        {
           log_kappa ~ nico_student_t(prior_kappa_dof, prior_kappa_mean, prior_kappa_SD);
           target += -log_kappa;
           ////
           if (use_Se_Sp_prior == 1 && n_chosen_cutoffs[p_class] > 0) {
             
                  ////
                  //// ---- Population medians and their derivatives, from the Beta of each cutpoint under
                  //// ---- the alpha of the study cutpoints (0.01 + phi * kappa), at the chosen cutoffs of p_class
                  //// ---- (C_mu is read nowhere else):
                  ////
                  real alpha_0 = sum(alpha);
                  vector[n_thr] P_of_alpha = cumulative_sum(alpha[1:n_thr]) / alpha_0;
                  vector[n_thr] P_upper_of_alpha = reverse(cumulative_sum(reverse(alpha[2:n_cat]))) / alpha_0;
                  vector[n_chosen_cutoffs[p_class]] C_mu;
                  vector[n_chosen_cutoffs[p_class]] log_J;
                  ////
                  for (j in 1:n_chosen_cutoffs[p_class]) {
                       int k = chosen_cutoffs[p_class, j];
                       vector[2] median_and_log_J = dirichlet_cutpoint_median_and_log_J(  P_of_alpha[k],
                                                                                               P_upper_of_alpha[k],
                                                                                               alpha_0,
                                                                                               use_probit_link);
                       C_mu[j] = median_and_log_J[1];
                       //// (chain rule: d C_mu / d P = d C_mu / d P_of_alpha * d P_of_alpha / d P = d C_mu / d P_of_alpha * kappa / alpha_0)
                       log_J[j] = median_and_log_J[2] + nico_log(kappa / alpha_0);
                  }
                  ////
                  //// ---- Baseline lines: probit(Sp) (c = 1) and probit(1 - Se) (c = 2) = (C_mu + Xbeta) * scale,
                  //// ---- at the prior means of beta_mu and raw_scale_mu (the constant anchor):
                  ////
                  real Xbeta_baseline_c      = dot_product(baseline_case, beta_mu);
                  real Xraw_scale_baseline_c = dot_product(baseline_case, raw_scale_mu);
                  vector[2] Xbeta_of_line = [mult_nd, mult_d]' * Xbeta_baseline_c;
                  vector[2] scale_of_line;
                  scale_of_line[1] = ((softplus == 1) ? 1.4426950408889634 * nico_log1p_exp(-mult_nd * Xraw_scale_baseline_c) :
                                                        nico_exp(-mult_nd * Xraw_scale_baseline_c));
                  scale_of_line[2] = ((softplus == 1) ? 1.4426950408889634 * nico_log1p_exp(-mult_d  * Xraw_scale_baseline_c) :
                                                        nico_exp(-mult_d  * Xraw_scale_baseline_c));
                  ////
                  array[2] vector[n_chosen_cutoffs[p_class]] C_baseline;
                  for (c in 1:2) {
                       C_baseline[c] = (C_mu + Xbeta_of_line[c]) * scale_of_line[c];
                  }
                  ////
                  //// ---- Between the chosen cutoffs: the Dirichlet(prior_dirichlet_alpha) of each segment on the increments of
                  //// ---- P_of_C_mu = (Sp + 1 - Se) / 2 at C_mu, both classes' own probabilities at the population medians, in place of the
                  //// ---- increments of the cumulative dirichlet_phi (the prior at the chosen cutoffs is unchanged):
                  ////
                  if (n_chosen_cutoffs[p_class] < n_thr) {
                       vector[n_thr] P_of_C_mu;
                       vector[n_thr] P_upper_of_C_mu;
                       vector[n_thr] log_J_of_P_of_C_mu;
                       real log_scale_nd = nico_log(scale_of_line[1]);
                       real log_scale_d = nico_log(scale_of_line[2]);
                       real log_kappa_over_alpha_0 = nico_log(kappa / alpha_0);
                       for (k in 1:n_thr) {
                            vector[2] median_and_log_J_of_k = dirichlet_cutpoint_median_and_log_J(P_of_alpha[k],
                                                                                                        P_upper_of_alpha[k],
                                                                                                        alpha_0,
                                                                                                        use_probit_link);
                            real link_Sp = (median_and_log_J_of_k[1] + Xbeta_of_line[1]) * scale_of_line[1];
                            real link_1m_Se = (median_and_log_J_of_k[1] + Xbeta_of_line[2]) * scale_of_line[2];
                            real log_f_Sp;
                            real log_f_1m_Se;
                            if (use_probit_link == 1) {
                                 P_of_C_mu[k] = 0.5 * (nico_Phi(link_Sp) + nico_Phi(link_1m_Se));
                                 P_upper_of_C_mu[k] = 0.5 * (nico_Phi(-link_Sp) + nico_Phi(-link_1m_Se));
                                 log_f_Sp = -0.5 * square(link_Sp) - 0.9189385332046728;
                                 log_f_1m_Se = -0.5 * square(link_1m_Se) - 0.9189385332046728;
                            } else {
                                 P_of_C_mu[k] = 0.5 * (nico_inv_logit(link_Sp) + nico_inv_logit(link_1m_Se));
                                 P_upper_of_C_mu[k] = 0.5 * (nico_inv_logit(-link_Sp) + nico_inv_logit(-link_1m_Se));
                                 log_f_Sp = nico_log_inv_logit(link_Sp) + nico_log1m_inv_logit(link_Sp);
                                 log_f_1m_Se = nico_log_inv_logit(link_1m_Se) + nico_log1m_inv_logit(link_1m_Se);
                            }
                            //// (log d P_of_C_mu / d cumulative phi = log((f(link_Sp) * scale_nd + f(link_1m_Se) * scale_d) / 2) + log dC_mu / d P)
                            real log_term_Sp = log_f_Sp + log_scale_nd;
                            real log_term_1m_Se = log_f_1m_Se + log_scale_d;
                            log_J_of_P_of_C_mu[k] = fmax(log_term_Sp, log_term_1m_Se) + nico_log1p_exp(-abs(log_term_Sp - log_term_1m_Se)) -
                                             0.6931471805599453 + median_and_log_J_of_k[2] + log_kappa_over_alpha_0;
                       }
                       ////
                       //// ---- Each segment with cutoffs inside it: log(new conditional) - log(Dirichlet of dirichlet_phi given the segment sum):
                       ////
                       int k_lower_end = 0;
                       for (j in 1:(n_chosen_cutoffs[p_class] + 1)) {
                            int k_upper_end = (j <= n_chosen_cutoffs[p_class]) ? chosen_cutoffs[p_class, j] : n_cat;
                            if (k_upper_end - k_lower_end > 1) {
                                 real P_of_C_mu_of_the_segment = (k_lower_end == 0) ? P_of_C_mu[k_upper_end] :
                                                         (k_upper_end == n_cat) ? P_upper_of_C_mu[k_lower_end] :
                                                         (P_of_C_mu[k_lower_end] > 0.5) ? (P_upper_of_C_mu[k_lower_end] - P_upper_of_C_mu[k_upper_end]) :
                                                                                        (P_of_C_mu[k_upper_end] - P_of_C_mu[k_lower_end]);
                                 real phi_of_the_segment = sum(dirichlet_phi[(k_lower_end + 1):k_upper_end]);
                                 target += sum(log_J_of_P_of_C_mu[(k_lower_end + 1):(k_upper_end - 1)]) -
                                           (sum(prior_dirichlet_alpha[(k_lower_end + 1):k_upper_end]) - 1.0) *
                                           (nico_log(P_of_C_mu_of_the_segment) - nico_log(phi_of_the_segment));
                                 for (i in (k_lower_end + 1):k_upper_end) {
                                      if (prior_dirichlet_alpha[i] != 1.0) {
                                           real P_of_C_mu_of_the_category = (i == n_cat) ? P_upper_of_C_mu[n_thr] :
                                                                    (i == 1) ? P_of_C_mu[1] :
                                                                    (P_of_C_mu[i - 1] > 0.5) ? (P_upper_of_C_mu[i - 1] - P_upper_of_C_mu[i]) :
                                                                                             (P_of_C_mu[i] - P_of_C_mu[i - 1]);
                                           target += (prior_dirichlet_alpha[i] - 1.0) *
                                                     (nico_log(P_of_C_mu_of_the_category) - nico_log(dirichlet_phi[i]));
                                      }
                                 }
                            }
                            k_lower_end = k_upper_end;
                       }
                  }
                  ////
                  //// ---- Class p_class: the Dirichlet of dirichlet_phi with its prior at the chosen cutoffs:
                  ////
                  {
                     matrix[n_chosen_cutoffs[p_class], 2] prior_beta_of_this_class = prior_Se_Sp_beta[p_class][1:n_chosen_cutoffs[p_class], ];
                     matrix[n_chosen_cutoffs[p_class], 2] prior_link_mean_SD_of_this_class = prior_Se_Sp_link_mean_SD[p_class][1:n_chosen_cutoffs[p_class], ];
                     if (p_class == 2) {
                          prior_link_mean_SD_of_this_class[, 1] = - prior_link_mean_SD_of_this_class[, 1];
                          if (Se_Sp_prior_type == 2) {
                               prior_beta_of_this_class = append_col(prior_beta_of_this_class[, 2],
                                                                     prior_beta_of_this_class[, 1]);
                          }
                     }
                     target += dirichlet_with_Se_Sp_prior_lupdf(
                               dirichlet_phi | log_J, C_baseline[p_class], 
                               prior_dirichlet_alpha,
                               alpha_of_the_segments[p_class][1:(n_chosen_cutoffs[p_class] + 1)],
                               chosen_cutoffs[p_class, 1:n_chosen_cutoffs[p_class]], Se_Sp_prior_type,
                               prior_Se_Sp_dirichlet_A[p_class][1:(n_chosen_cutoffs[p_class] + 1)], prior_beta_of_this_class,
                               prior_link_mean_SD_of_this_class,
                               use_probit_link);
                     ////
                     //// ---- Jacobian from C_mu to the scaled line at the chosen cutpoints:
                     ////
                     if (joint_unit_probability_prior == 0 && joint_link_prior == 0) {
                          target += n_chosen_cutoffs[p_class] * nico_log(scale_of_line[p_class]);
                     }
                  }
                  ////
                  //// ---- The other class: its Se/Sp prior at the same cutoff, through Xbeta_baseline_c (beta_mu[1]);
                  //// ---- the Jacobian of (P, b) -> (C_baseline[1], C_baseline[2]) is 2 * scale_nd * scale_d * |dm/dP|,
                  //// ---- of which the terms of p_class above hold scale_of_line[p_class] * |dm/dP|:
                  ////
                  {
                     int q_class = 3 - p_class;
                     matrix[n_chosen_cutoffs[q_class], 2] prior_beta_of_this_class = prior_Se_Sp_beta[q_class][1:n_chosen_cutoffs[q_class], ];
                     matrix[n_chosen_cutoffs[q_class], 2] prior_link_mean_SD_of_this_class = prior_Se_Sp_link_mean_SD[q_class][1:n_chosen_cutoffs[q_class], ];
                     if (q_class == 2) {
                          prior_link_mean_SD_of_this_class[, 1] = - prior_link_mean_SD_of_this_class[, 1];
                          if (Se_Sp_prior_type == 2) {
                               prior_beta_of_this_class = append_col(prior_beta_of_this_class[, 2],
                                                                     prior_beta_of_this_class[, 1]);
                          }
                     }
                     target += Se_Sp_prior_at_chosen_cutoffs_lupdf(
                               C_baseline[q_class] | Se_Sp_prior_type,
                               prior_Se_Sp_dirichlet_A[q_class][1:(n_chosen_cutoffs[q_class] + 1)], prior_beta_of_this_class,
                               prior_link_mean_SD_of_this_class,
                               use_probit_link);
                     if (joint_unit_probability_prior == 1) {
                          //// Both ordered Normal lines share their centred direction, not their means or radii.
                          target += nico_log(4.0) + 0.5 * (n_thr - 1) *
                                    nico_log(dot_self(C_mu - mean(C_mu)));
                     } else if (joint_link_prior == 1) {
                          //// Link-space volume of both lines with respect to C_mu, location and raw scale.
                          real a = scale_of_line[1];
                          real d = scale_of_line[2];
                          real da = a;
                          real dd = -d;
                          if (softplus == 1) {
                               da = 1.4426950408889634 * nico_inv_logit(Xraw_scale_baseline_c);
                               dd = -1.4426950408889634 * nico_inv_logit(-Xraw_scale_baseline_c);
                          }
                          target += nico_log(2.0 * a * d * (da * d - a * dd)) +
                                    0.5 * nico_log(n_chosen_cutoffs[p_class] *
                                              dot_self(C_mu - mean(C_mu))) +
                                    0.5 * (n_chosen_cutoffs[p_class] - 2) *
                                    nico_log(square(a) + square(d));
                     } else {
                          target += nico_log(2.0 * scale_of_line[q_class]);
                     }
                  }
                  if (n_covariates > 1) {
                       beta_mu[2:n_covariates] ~ nico_normal( prior_beta_mu_mean[2:n_covariates],
                                                         prior_beta_mu_SD[2:n_covariates]);
                  }
                  if (joint_unit_probability_prior == 1 || joint_link_prior == 1) {
                       if (n_covariates > 1) {
                            raw_scale_mu[2:n_covariates] ~ nico_normal(
                                  prior_raw_scale_mu_mean[2:n_covariates],
                                  prior_raw_scale_mu_SD[2:n_covariates]);
                       }
                  } else {
                       raw_scale_mu[1:n_covariates] ~ nico_normal(
                             prior_raw_scale_mu_mean[1:n_covariates],
                             prior_raw_scale_mu_SD[1:n_covariates]);
                  }
                  
           } else {
             
                  beta_mu[1:n_covariates] ~ nico_normal( prior_beta_mu_mean[1:n_covariates],
                                                         prior_beta_mu_SD[1:n_covariates]);
                  raw_scale_mu[1:n_covariates] ~ nico_normal( prior_raw_scale_mu_mean[1:n_covariates],
                                                              prior_raw_scale_mu_SD[1:n_covariates]);
                  ////
                  dirichlet_phi ~ nico_dirichlet( prior_dirichlet_alpha);
                  
           }
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
        //// (anchored at the mid-location of the baseline case at the prior means, under every prior set)
        {
             //// (the reported cutpoints of each study only: the alpha of each segment between them is the sum of
             ////  alpha over its categories, the same distribution with the other cutpoints integrated out)
             vector[n_cat + 1] cumulative_alpha = append_row(0.0, cumulative_sum(alpha));
             vector[n_C + n_studies] alpha_of_each_segment = cumulative_alpha[upper_index_of_each_segment] -
                                                             cumulative_alpha[lower_index_of_each_segment];
             for (s in 1:n_studies) {
                  if (n_reported_cutpoints[s] > 0) {
                       int p_1 = first_C_of_each_study[s];
                       int p_2 = p_1 + n_reported_cutpoints[s] - 1;
                       target += induced_dirichlet_given_C_on_the_log_scale_lupdf(C[p_1:p_2] -
                                                                                  0.5 * (mult_nd + mult_d) *
                                                                                  dot_product(baseline_case, prior_beta_mu_mean) |
                                                                                  alpha_of_each_segment[(p_1 + s - 1):(p_2 + s)],
                                                                                  use_probit_link);
                  }
             }
        }
        ////
        //// ---- Likelihood / Model:
        ////
        target += nico_normal_lupdf(to_vector(beta_z) | 0, 1); // part of between-study model, NOT prior
        target += nico_normal_lupdf(to_vector(raw_scale_z) | 0, 1); // part of between-study model, NOT prior
        ////
        //// (a study's first reported cutpoint is C_raw itself; each later one adds an exp() or softplus step)
        target += (softplus == 1) ? sum(nico_log_inv_logit(C_raw_between_its_bounds[position_of_each_increment])) :
                                    sum(C_raw_between_its_bounds[position_of_each_increment]);
        target += log_J_of_the_free_values;
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
            array[2] matrix[n_studies, n_thr] latent_surv = map_latent_surv_prob_to_the_reported_C_HSROC(
                                                            C, n_thr, Xlocations, Xscale, n_studies,
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
         ////  study cutpoints are anchored at 0.5 * (mult_nd + mult_d) * dot_product(baseline_case, prior_beta_mu_mean) = 0)
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
      if (log_lik_indicator == 1) {
          array[2] matrix[n_studies, n_thr] latent_surv = map_latent_surv_prob_to_the_reported_C_HSROC(
                                                          C, n_thr, Xlocations, Xscale, n_studies,
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

