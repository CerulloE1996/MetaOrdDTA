//// DTA_MA_ord_biv_RANDthr_kappa_popmed_Se_Sp_prior_Stan_math_only.stan
//// The same data, parameters and target as the fused AVX file of this model with the Dirichlet prior of dirichlet_phi, on Stan::math only.
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

        //// ---- Exact normalized Set-0 marginal log-density at chosen cutoffs via 24-node Laplace-centered Gauss-Hermite:
        real log_Set0_marginal_at_chosen_cutoffs(vector x_chosen,
                                                 vector alpha_of_the_segments,
                                                 real v0,
                                                 int use_probit_link) {
              int m = num_elements(x_chosen);
              if (m == 0) return 0.0;
              vector[24] gh_t = [
                  -8.50780351919527433e+00, -7.43789066602167104e+00, -6.54167500509862876e+00, -5.73274717525120181e+00,
                  -4.97804137463911367e+00, -4.26038360501990354e+00, -3.56930676407355918e+00, -2.89772864322331181e+00,
                  -2.24046785169175244e+00, -1.59348042981641846e+00, -9.53421922932108146e-01, -3.17370096629455034e-01,
                   3.17370096629453258e-01,  9.53421922932111698e-01,  1.59348042981641846e+00,  2.24046785169175156e+00,
                   2.89772864322331269e+00,  3.56930676407356007e+00,  4.26038360501990709e+00,  4.97804137463911989e+00,
                   5.73274717525120003e+00,  6.54167500509863320e+00,  7.43789066602166216e+00,  8.50780351919525657e+00 ]';
              vector[24] gh_lw = [
                  -3.69042806607304854e+01, -2.86212344883108756e+01, -2.24843031451895179e+01, -1.76020197436324572e+01,
                  -1.36185797828409356e+01, -1.03468101270751411e+01, -7.67407824125704963e+00, -5.52732585071174043e+00,
                  -3.85723446710271611e+00, -2.63012623460007200e+00, -1.82350086207048112e+00, -1.42349743041213506e+00,
                  -1.42349743041215793e+00, -1.82350086207048467e+00, -2.63012623460007866e+00, -3.85723446710268059e+00,
                  -5.52732585071170757e+00, -7.67407824125703453e+00, -1.03468101270751269e+01, -1.36185797828409267e+01,
                  -1.76020197436324182e+01, -2.24843031451895001e+01, -2.86212344883107441e+01, -3.69042806607303291e+01 ]';
              real sum_seg = sum(alpha_of_the_segments);
              real lg_const = lgamma(sum_seg) - sum(lgamma(alpha_of_the_segments));
              real inv_v0 = 1.0 / v0;
              real b_star = 0.0;
              {
                    real cum_w = 0.0;
                    for (j in 1:m) {
                          cum_w += alpha_of_the_segments[j];
                          real p_m = cum_w / sum_seg;
                          real q_m = (use_probit_link == 1) ? inv_Phi(p_m) : logit(p_m);
                          b_star += (q_m - x_chosen[j]);
                    }
                    b_star /= m;
              }
              real sig_star = 1.0;
              for (it in 1:11) {
                    vector[m] z = x_chosen + b_star;
                    vector[m] lF;
                    vector[m] lS;
                    vector[m] lf;
                    if (use_probit_link == 1) {
                          for (j in 1:m) {
                                lF[j] = std_normal_lcdf(z[j] |);
                                lS[j] = std_normal_lccdf(z[j] |);
                                lf[j] = -0.5 * square(z[j]) - 0.91893853320467274178;
                          }
                    } else {
                          lF = log_inv_logit(z);
                          lS = log1m_inv_logit(z);
                          lf = lF + lS;
                    }
                    vector[m + 1] log_q;
                    log_q[1] = lF[1];
                    log_q[m + 1] = lS[m];
                    for (j in 2:m) {
                          log_q[j] = (z[j - 1] >= 0.0) ? log_diff_exp(lS[j - 1], lS[j]) : log_diff_exp(lF[j], lF[j - 1]);
                    }
                    vector[m] r_own = exp(lf - log_q[1:m]);
                    vector[m] r_nxt = exp(lf - log_q[2:(m + 1)]);
                    real sum_dg = 0.0;
                    real sum_d2 = -inv_v0;
                    for (j in 1:m) {
                          real df_j = (use_probit_link == 1) ? (-z[j]) : (exp(lS[j]) - exp(lF[j]));
                          sum_dg += (alpha_of_the_segments[j] - 1.0) * r_own[j] - (alpha_of_the_segments[j + 1] - 1.0) * r_nxt[j] + df_j;
                          sum_d2 -= (use_probit_link == 1) ? 1.0 : (2.0 * exp(lf[j]));
                    }
                    if (alpha_of_the_segments[1] != 1.0) {
                          real df1 = (use_probit_link == 1) ? (-z[1]) : (exp(lS[1]) - exp(lF[1]));
                          sum_d2 += (alpha_of_the_segments[1] - 1.0) * (df1 * r_own[1] - square(r_own[1]));
                    }
                    if (alpha_of_the_segments[m + 1] != 1.0) {
                          real dfm = (use_probit_link == 1) ? (-z[m]) : (exp(lS[m]) - exp(lF[m]));
                          sum_d2 += (alpha_of_the_segments[m + 1] - 1.0) * (-dfm * r_nxt[m] - square(r_nxt[m]));
                    }
                    for (i in 2:m) {
                          if (alpha_of_the_segments[i] != 1.0) {
                                real df_hi = (use_probit_link == 1) ? (-z[i]) : (exp(lS[i]) - exp(lF[i]));
                                real df_lo = (use_probit_link == 1) ? (-z[i - 1]) : (exp(lS[i - 1]) - exp(lF[i - 1]));
                                real ui = r_own[i] - r_nxt[i - 1];
                                sum_d2 += (alpha_of_the_segments[i] - 1.0) * (df_hi * r_own[i] - df_lo * r_nxt[i - 1] - square(ui));
                          }
                    }
                    if (it == 11) {
                          sig_star = inv_sqrt(-sum_d2);
                    } else {
                          real dg_db = -b_star * inv_v0 + sum_dg;
                          real step = -dg_db / sum_d2;
                          if (step >  2.0) step =  2.0;
                          if (step < -2.0) step = -2.0;
                          b_star += step;
                    }
              }
              vector[24] lw_terms;
              for (i in 1:24) {
                    real b_i = b_star + sig_star * gh_t[i];
                    vector[m] z = x_chosen + b_i;
                    vector[m] lF;
                    vector[m] lS;
                    vector[m] lf;
                    if (use_probit_link == 1) {
                          for (j in 1:m) {
                                lF[j] = std_normal_lcdf(z[j] |);
                                lS[j] = std_normal_lccdf(z[j] |);
                                lf[j] = -0.5 * square(z[j]) - 0.91893853320467274178;
                          }
                    } else {
                          lF = log_inv_logit(z);
                          lS = log1m_inv_logit(z);
                          lf = lF + lS;
                    }
                    vector[m + 1] log_q;
                    log_q[1] = lF[1];
                    log_q[m + 1] = lS[m];
                    for (j in 2:m) {
                          log_q[j] = (z[j - 1] >= 0.0) ? log_diff_exp(lS[j - 1], lS[j]) : log_diff_exp(lF[j], lF[j - 1]);
                    }
                    real val_i = -0.5 * square(b_i) * inv_v0 + lg_const + sum(lf);
                    for (s_i in 1:(m + 1)) {
                          if (alpha_of_the_segments[s_i] != 1.0) {
                                val_i += (alpha_of_the_segments[s_i] - 1.0) * log_q[s_i];
                          }
                    }
                    lw_terms[i] = gh_lw[i] + val_i + 0.5 * square(gh_t[i]);
              }
              return -0.5 * log(v0) + log(sig_star) + log_sum_exp(lw_terms);
        }
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
        // array[2] simplex[n_cat] dirichlet_phi;
        // for (c in 1:2) {
        //   dirichlet_phi[c] = stickbreaking_logistic_simplex_constrain(raw_dirichlet_phi[c]); // THIS IS WRONG!!!!!!
        // }
        ////
        vector[2] log_kappa = log(kappa);
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
        //// ---- C_raw (a study's first reported cutpoint; the log of each later step), between the same bounds as before,
        ////      from its free value: on the logit scale of its place between its bounds,
        ////           logit((C_raw - lower) / (upper - lower)) = logit((centre - lower) / (upper - lower)) + spread x free value.
        ////      The centre and the spread are those of the Dirichlet(alpha) of the study cutpoints when alpha is large:
        ////      the cumulative probability of cutpoint k has mean P_k = sum(alpha[1:k]) / alpha_0 and variance
        ////      P_k (1 - P_k) / (alpha_0 + 1); two cutpoints j < k have covariance P_j (1 - P_k) / (alpha_0 + 1); each is
        ////      carried to the link scale by the delta method. Any centre and any spread give the same model: this is
        ////      a change of variables, with its Jacobian (log_J_of_the_free_values) added to the target in the model block.
        ////
        array[2] vector[n_C] C_raw_between_its_bounds;
        vector[2] log_J_of_the_free_values = rep_vector(0.0, 2);
        {
              vector[2] baseline_location_prior_of_C_raw;
              baseline_location_prior_of_C_raw[1] = dot_product(baseline_case_nd, to_vector(prior_beta_mu_mean[1, 1:n_covariates_nd]));
              baseline_location_prior_of_C_raw[2] = dot_product(baseline_case_d,  to_vector(prior_beta_mu_mean[2, 1:n_covariates_d]));
              for (c in 1:2) {
                    real alpha_0 = sum(alpha[c]);
                    vector[n_thr] P = cumulative_sum(alpha[c][1:n_thr]) / alpha_0;
                    //// (1 - P from the sum of the upper categories, never as 1 - P: a cutpoint in the upper half is taken
                    ////  from its upper tail, so a P that rounds to 1 gives no infinite centre)
                    vector[n_thr] P_upper = reverse(cumulative_sum(reverse(alpha[c][2:n_cat]))) / alpha_0;
                    vector[n_thr] centre_of_each_cutpoint;
                    vector[n_thr] density_at_each_cutpoint;
                    if (use_probit_link == 1) {
                          for (k in 1:n_thr) {
                                centre_of_each_cutpoint[k] = (P[k] <= 0.5) ? inv_Phi(P[k]) : -inv_Phi(P_upper[k]);
                          }
                          density_at_each_cutpoint = exp(-0.5 * square(centre_of_each_cutpoint)) / sqrt(2.0 * pi());
                    } else {
                          centre_of_each_cutpoint = log(P) - log(P_upper);
                          density_at_each_cutpoint = P .* P_upper;
                    }
                    for (p in 1:n_C) {
                          int k = k_of_each_C[p];
                          int k_previous = k_previous_of_each_C[p];
                          real centre;
                          real spread;
                          if (k_previous == 0) {
                                centre = centre_of_each_cutpoint[k] + baseline_location_prior_of_C_raw[c];
                                spread = sqrt(P[k] * P_upper[k] / (alpha_0 + 1.0)) / density_at_each_cutpoint[k];
                          } else {
                                real step_between_the_centres = centre_of_each_cutpoint[k] - centre_of_each_cutpoint[k_previous];
                                real variance_of_the_step =
                                      (P[k] * P_upper[k] / square(density_at_each_cutpoint[k]) +
                                       P[k_previous] * P_upper[k_previous] / square(density_at_each_cutpoint[k_previous]) -
                                       2.0 * P[k_previous] * P_upper[k] / (density_at_each_cutpoint[k] * density_at_each_cutpoint[k_previous])) /
                                      (alpha_0 + 1.0);
                                centre = log(step_between_the_centres);
                                spread = sqrt(fmax(variance_of_the_step, 1e-14)) / step_between_the_centres;
                          }
                          real lower_bound_of_this_C_raw = C_raw_lower_of_each_C[c][p];
                          real width = C_raw_upper_of_each_C[c][p] - lower_bound_of_this_C_raw;
                          //// (the place of the centre between the bounds, held inside (0.001, 0.999) by a smooth function that
                          ////  leaves a place well inside as it is; the spread on the logit scale never above 4, smoothly: a
                          ////  centre at a bound, as the log of a step between two cutpoints that almost coincide, then gives a
                          ////  free value that spans the whole width instead of an overflow. Any centre and any spread give
                          ////  the same model, so neither changes the target; both depend on kappa and dirichlet_phi only.)
                          real place_as_it_is = (centre - lower_bound_of_this_C_raw) / width;
                          real place_of_the_centre = 0.001 + (log1p_exp(500.0 * (place_as_it_is - 0.001)) -
                                                              log1p_exp(500.0 * (place_as_it_is - 0.999))) / 500.0;
                          //// (the spread of the free-value form on the logit scale, never above 1, the scale of C_raw itself)
                          real spread_of_the_free_value_form = spread / sqrt(square(width * place_of_the_centre * (1.0 - place_of_the_centre)) +
                                                                             square(spread));
                          //// (Sean Pinkney partially-centered parameterization, partially_centered_example_Sean_Pinkney.stan
                          ////  & spinkney/hierarchical_model_tutorial stan/hierarchical_pcp.stan:
                          ////  tau_wt = tau - w * (tau - 1); theta = (w * mu + tau * eta) / tau_wt;
                          ////  w in [0, 1] is fixed data: w = 0 gives centered, w = 1 gives non-centered,
                          ////  with exact log-Jacobian log(tau) - log(tau_wt))
                          real w = Pinkney_w[c][p];
                          real tau = spread_of_the_free_value_form;
                          real tau_wt = tau - w * (tau - 1.0);
                          real eta_shifted = (1.0 - w) * Pinkney_eta0[c][p] + C_raw[c][p];
                          real logit_of_the_place = (w * logit(place_of_the_centre) + tau * eta_shifted) / tau_wt;
                          C_raw_between_its_bounds[c][p] = lower_bound_of_this_C_raw + width * inv_logit(logit_of_the_place);
                          log_J_of_the_free_values[c] += log(width) + log_inv_logit(logit_of_the_place) + log1m_inv_logit(logit_of_the_place) +
                                               log(tau) - log(tau_wt);
                    }
              }
        }
        ////
        //// ---- Construct study-specific cutpoints:
        ////
        // array[2] matrix[n_studies, n_thr] C;
        // for (c in 1:2) {
        //   for (s in 1:n_studies) {
        //        C[c][s, ] =  construct_C(C_raw_between_its_bounds[c][s, 1:n_thr], softplus);
        //   }
        // }
        //// (the reported cutpoints of each study only, study after study, as construct_C on each study's own)
        array[2] vector[n_C] C;
        for (c in 1:2) {
              vector[n_C - n_studies_with_a_cutpoint] increments_of_the_cutpoints =
                    (softplus == 1) ? log1p_exp(C_raw_between_its_bounds[c][position_of_each_increment]) :
                                      exp(C_raw_between_its_bounds[c][position_of_each_increment]);
              C[c] = C_raw_between_its_bounds[c];
              for (i in 1:(n_C - n_studies_with_a_cutpoint)) {
                    int p = position_of_each_increment[i];
                    C[c][p] = C[c][p - 1] + increments_of_the_cutpoints[i];
              }
        }
        ////
        //// ---- (C holds the study cutpoints themselves under every prior set: they are never moved by beta_mu)
        ////
        //// ---- Construct simple 2x2 (bivariate) between-study corr matrices:
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
        beta_mu[1, 1:n_covariates_nd] ~ normal(
                                          prior_beta_mu_mean[1, 1:n_covariates_nd],
                                          prior_beta_mu_SD[1, 1:n_covariates_nd]);
        beta_mu[2, 1:n_covariates_d]  ~ normal(
                                          prior_beta_mu_mean[2, 1:n_covariates_d],
                                          prior_beta_mu_SD[2, 1:n_covariates_d]);
        ////
        to_vector(beta_SD) ~ normal(prior_beta_SD_mean, prior_beta_SD_SD);
        beta_L_Omega ~ lkj_corr_cholesky(prior_beta_corr_LKJ);
        ////
        vector[2] baseline_location_prior;  //// the baseline case at the prior means (the same under every prior set)
        baseline_location_prior[1] = dot_product(baseline_case_nd, to_vector(prior_beta_mu_mean[1, 1:n_covariates_nd]));
        baseline_location_prior[2] = dot_product(baseline_case_d,  to_vector(prior_beta_mu_mean[2, 1:n_covariates_d]));
        ////
        real Xbeta_baseline_nd = dot_product(baseline_case_nd, to_vector(beta_mu[1, 1:n_covariates_nd]));
        real Xbeta_baseline_d  = dot_product(baseline_case_d,  to_vector(beta_mu[2, 1:n_covariates_d]));
        ////
        for (c in 1:2) {
           log_kappa[c] ~ student_t(prior_kappa_dof[c], prior_kappa_mean[c], prior_kappa_SD[c]); // kappa[c] ~ normal(0, 100);
           target += -log_kappa[c];
           ////
           if (use_Se_Sp_prior >= 1 && n_chosen_cutoffs[c] > 0) {

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
                    //// ---- C_mu - Xbeta_baseline of class c at EVERY threshold: the quantity that Se and Sp are made of
                    ////      (generated quantities: Phi(-(C_mu - Xbeta_baseline))). C_mu is baseline_location_prior plus the
                    ////      population median of the cutpoint, as in generated quantities.
                    ////      Its prior is the fixed-cutpoint model's own prior of C - Xbeta_baseline (the Se / Sp prior at
                    ////      the chosen cutoffs, the Dirichlet of each segment between them), with the Jacobian of every median:
                    if (use_Se_Sp_prior == 2) {
                         vector[n_thr] C_mu_minus_baseline_location_prior;
                         real log_J_of_every_median = 0.0;
                         for (k in 1:n_thr) {
                              vector[2] median_and_log_J = dirichlet_cutpoint_median_and_log_J( P_of_alpha[k],
                                                                                                P_upper_of_alpha[k],
                                                                                                alpha_0,
                                                                                                use_probit_link);
                              C_mu_minus_baseline_location_prior[k] = median_and_log_J[1];
                              log_J_of_every_median += median_and_log_J[2] + log(kappa[c] / alpha_0);
                         }
                         target += log_J_of_every_median;
                         target += induced_dirichlet_with_Se_Sp_prior_lupdf(
                                   C_mu_minus_baseline_location_prior | prior_dirichlet_alpha[c],
                                   alpha_of_the_segments[c][1:(n_chosen_cutoffs[c] + 1)],
                                   chosen_cutoffs[c, 1:n_chosen_cutoffs[c]], Se_Sp_prior_type,
                                   prior_Se_Sp_dirichlet_A[c][1:(n_chosen_cutoffs[c] + 1)],
                                   prior_beta_of_this_class,
                                   prior_link_mean_SD_of_this_class,
                                   use_probit_link);
                    } else {
                         int m_c = n_chosen_cutoffs[c];
                         array[m_c] int cuts_c = chosen_cutoffs[c, 1:m_c];
                         vector[m_c + 1] alpha_seg_c = alpha_of_the_segments[c][1:(m_c + 1)];
                         vector[m_c] popmed_chosen;
                         vector[m_c] log_J_chosen;
                         real log_kappa_over_alpha0 = log(kappa[c] / alpha_0);
                         for (j in 1:m_c) {
                              int k = cuts_c[j];
                              vector[2] median_and_log_J = dirichlet_cutpoint_median_and_log_J( P_of_alpha[k],
                                                                                                P_upper_of_alpha[k],
                                                                                                alpha_0,
                                                                                                use_probit_link);
                              popmed_chosen[j] = median_and_log_J[1];
                              log_J_chosen[j]  = median_and_log_J[2] + log_kappa_over_alpha0;
                         }
                         vector[m_c] x_chosen = (baseline_location_prior[c] + popmed_chosen) - Xbeta_baseline;
                         real v0 = (c == 1) ?
                                   dot_self(baseline_case_nd .* to_vector(prior_beta_mu_SD[1, 1:n_covariates_nd])) :
                                   dot_self(baseline_case_d  .* to_vector(prior_beta_mu_SD[2, 1:n_covariates_d]));
                         real norm_const_dir = (lgamma(sum(prior_dirichlet_alpha[c])) - sum(lgamma(prior_dirichlet_alpha[c])))
                                               - ((use_probit_link == 1) ? (0.91893853320467274178 * m_c) : 0.0);
                         target += dirichlet_with_Se_Sp_prior_lupdf(
                                   dirichlet_phi[c] | log_J_chosen, popmed_chosen,
                                   prior_dirichlet_alpha[c], alpha_seg_c, cuts_c, 1,
                                   alpha_seg_c, prior_beta_of_this_class, prior_link_mean_SD_of_this_class, use_probit_link)
                                   - norm_const_dir;
                         target += - log_Set0_marginal_at_chosen_cutoffs(x_chosen, alpha_seg_c, v0, use_probit_link);
                         target += Se_Sp_prior_at_chosen_cutoffs_lupdf(
                                   x_chosen | Se_Sp_prior_type,
                                   prior_Se_Sp_dirichlet_A[c][1:(m_c + 1)],
                                   prior_beta_of_this_class, prior_link_mean_SD_of_this_class, use_probit_link);
                    }

           } else {

                     dirichlet_phi[c] ~ dirichlet(prior_dirichlet_alpha[c]);

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
             //// (the study cutpoints are measured from baseline_location_prior, the baseline case at the prior means,
             ////  under every prior set)
             //// (the reported cutpoints of each study only: the alpha of each segment between them is the sum of
             ////  alpha[c] over its categories, the same distribution with the other cutpoints integrated out)
             vector[n_cat + 1] cumulative_alpha = append_row(0.0, cumulative_sum(alpha[c]));
             vector[n_C + n_studies] alpha_of_each_segment = cumulative_alpha[upper_index_of_each_segment] -
                                                             cumulative_alpha[lower_index_of_each_segment];
             for (s in 1:n_studies) {
                  if (n_reported_cutpoints[s] > 0) {
                       int p_1 = first_C_of_each_study[s];
                       int p_2 = p_1 + n_reported_cutpoints[s] - 1;
                       target += induced_dirichlet_given_C_on_the_log_scale_lupdf(C[c][p_1:p_2] - baseline_location_prior[c] |
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
        target += normal_lupdf(to_vector(to_matrix(beta_z)) | 0, 1); // part of between-study model, NOT prior
        ////
        for (c in 1:2) {
           //// (a study's first reported cutpoint is C_raw itself; each later one adds an exp() or softplus step)
           target += (softplus == 1) ? sum(log_inv_logit(C_raw_between_its_bounds[c][position_of_each_increment])) :
                                       sum(C_raw_between_its_bounds[c][position_of_each_increment]);
        }
        //// (the Jacobian of C_raw as a function of its free value; before, Stan added that of the bounds of C_raw itself)
        target += sum(log_J_of_the_free_values);
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
         //// (the baseline case at the prior means, as baseline_location_prior in the model block, under every prior set)
         real baseline_location_prior_of_class_c =
                    (c == 1) ? dot_product(baseline_case_nd, to_vector(prior_beta_mu_mean[1, 1:n_covariates_nd])) :
                               dot_product(baseline_case_d,  to_vector(prior_beta_mu_mean[2, 1:n_covariates_d]));
         for (k in 1:n_thr) {
            C_mu[c][k] = baseline_location_prior_of_class_c +
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
                                  Phi( -(C_mu[1] - Xbeta_baseline_nd))  : 
                                  inv_logit(-(C_mu[1] - Xbeta_baseline_nd));
      vector[n_thr] Sp_baseline = 1.0 - Fp_baseline;
      vector[n_thr] Se_baseline = (use_probit_link == 1) ? 
                                  Phi( -(C_mu[2] - Xbeta_baseline_d))   : 
                                  inv_logit(-(C_mu[2] - Xbeta_baseline_d));
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
              C_pred[c][k] = C_mu[c][k] + normal_rng(0.0, sigma_C_k);
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
                                       Phi(-(C_pred[1] - Xbeta_baseline_pred_nd)) :
                                       inv_logit(-(C_pred[1] - Xbeta_baseline_pred_nd));
      vector[n_thr] Sp_baseline_pred = 1.0 - Fp_baseline_pred;
      vector[n_thr] Se_baseline_pred = (use_probit_link == 1) ?
                                       Phi(-(C_pred[2] - Xbeta_baseline_pred_d))  :
                                       inv_logit(-(C_pred[2] - Xbeta_baseline_pred_d));
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









