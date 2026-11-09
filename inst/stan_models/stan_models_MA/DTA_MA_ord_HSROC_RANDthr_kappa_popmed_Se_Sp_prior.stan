//// DTA_MA_ord_HSROC_RANDthr_kappa_popmed_Se_Sp_prior_Stan_math_only.stan
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
                              upper_index_of_each_segment[p + s - 1] = k + 1;
                              lower_index_of_each_segment[p + s - 1] = k_previous + 1;
                              if (k_previous == 0) {
                                   C_raw_lower_of_each_C[p] = C_raw_lower + (k - 1) * step_lower;
                                   C_raw_upper_of_each_C[p] = C_raw_upper + (k - 1) * step_upper;
                              } else {
                                   i += 1;
                                   position_of_each_increment[i] = p;
                                   C_raw_lower_of_each_C[p] = (softplus == 1) ? log(expm1((k - k_previous) * step_lower)) :
                                                                                 C_raw_lower + log(k - k_previous);
                                   C_raw_upper_of_each_C[p] = (softplus == 1) ? log(expm1((k - k_previous) * step_upper)) :
                                                                                 C_raw_upper + log(k - k_previous);
                              }
                              k_previous = k;
                         }
                    }
                    upper_index_of_each_segment[p + s] = n_cat + 1;
                    lower_index_of_each_segment[p + s] = k_previous + 1;
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
        vector<lower=C_raw_lower_of_each_C, upper=C_raw_upper_of_each_C>[n_C] C_raw; //// the reported cutpoints of each study only (restricting solves divs/numerical stab. issues w/ lots of sparse categories!)
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

        //// (the reported cutpoints of each study only, study after study, as construct_C on each study's own)
        vector[n_C] C = C_raw;
        {
              vector[n_C - n_studies_with_a_cutpoint] increments_of_the_cutpoints =
                    (softplus == 1) ? log1p_exp(C_raw[position_of_each_increment]) :
                                      exp(C_raw[position_of_each_increment]);
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
            Xscale[, 1] =  ((softplus == 1) ? 1.4426950408889634 * log1p_exp( -mult_nd*Xraw_scale_baseline) :
                                                exp( -mult_nd*Xraw_scale_baseline)); 
            // Xscale[, 2] =  ((softplus == 1) ? softplus_scaled( -mult_d*Xraw_scale_baseline)  : exp( -mult_d*Xraw_scale_baseline));
            Xscale[, 2] =  ((softplus == 1) ? 1.4426950408889634 * log1p_exp( -mult_d*Xraw_scale_baseline)  :
                                                exp( -mult_d*Xraw_scale_baseline));
        }
}


model {
        ////
        //// ---- Priors:
        ////
        beta_SD ~ normal( prior_beta_SD_mean, prior_beta_SD_SD);
        raw_scale_SD ~ normal( prior_raw_scale_SD_mean, prior_raw_scale_SD_SD);
        ////
        {
           log_kappa ~ student_t(prior_kappa_dof, prior_kappa_mean, prior_kappa_SD);
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
                       log_J[j] = median_and_log_J[2] + log(kappa / alpha_0);
                  }
                  ////
                  //// ---- Baseline lines: probit(Sp) (c = 1) and probit(1 - Se) (c = 2) = (C_mu + Xbeta) * scale,
                  //// ---- at the prior means of beta_mu and raw_scale_mu (the constant anchor):
                  ////
                  real Xbeta_baseline_c      = dot_product(baseline_case, beta_mu);
                  real Xraw_scale_baseline_c = dot_product(baseline_case, raw_scale_mu);
                  vector[2] Xbeta_of_line = [mult_nd, mult_d]' * Xbeta_baseline_c;
                  vector[2] scale_of_line;
                  scale_of_line[1] = ((softplus == 1) ? 1.4426950408889634 * log1p_exp(-mult_nd * Xraw_scale_baseline_c) :
                                                        exp(-mult_nd * Xraw_scale_baseline_c));
                  scale_of_line[2] = ((softplus == 1) ? 1.4426950408889634 * log1p_exp(-mult_d  * Xraw_scale_baseline_c) :
                                                        exp(-mult_d  * Xraw_scale_baseline_c));
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
                       real log_scale_nd = log(scale_of_line[1]);
                       real log_scale_d = log(scale_of_line[2]);
                       real log_kappa_over_alpha_0 = log(kappa / alpha_0);
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
                                 P_of_C_mu[k] = 0.5 * (Phi(link_Sp) + Phi(link_1m_Se));
                                 P_upper_of_C_mu[k] = 0.5 * (Phi(-link_Sp) + Phi(-link_1m_Se));
                                 log_f_Sp = -0.5 * square(link_Sp) - 0.9189385332046728;
                                 log_f_1m_Se = -0.5 * square(link_1m_Se) - 0.9189385332046728;
                            } else {
                                 P_of_C_mu[k] = 0.5 * (inv_logit(link_Sp) + inv_logit(link_1m_Se));
                                 P_upper_of_C_mu[k] = 0.5 * (inv_logit(-link_Sp) + inv_logit(-link_1m_Se));
                                 log_f_Sp = log_inv_logit(link_Sp) + log1m_inv_logit(link_Sp);
                                 log_f_1m_Se = log_inv_logit(link_1m_Se) + log1m_inv_logit(link_1m_Se);
                            }
                            //// (log d P_of_C_mu / d cumulative phi = log((f(link_Sp) * scale_nd + f(link_1m_Se) * scale_d) / 2) + log dC_mu / d P)
                            real log_term_Sp = log_f_Sp + log_scale_nd;
                            real log_term_1m_Se = log_f_1m_Se + log_scale_d;
                            log_J_of_P_of_C_mu[k] = fmax(log_term_Sp, log_term_1m_Se) + log1p_exp(-abs(log_term_Sp - log_term_1m_Se)) -
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
                                           (log(P_of_C_mu_of_the_segment) - log(phi_of_the_segment));
                                 for (i in (k_lower_end + 1):k_upper_end) {
                                      if (prior_dirichlet_alpha[i] != 1.0) {
                                           real P_of_C_mu_of_the_category = (i == n_cat) ? P_upper_of_C_mu[n_thr] :
                                                                    (i == 1) ? P_of_C_mu[1] :
                                                                    (P_of_C_mu[i - 1] > 0.5) ? (P_upper_of_C_mu[i - 1] - P_upper_of_C_mu[i]) :
                                                                                             (P_of_C_mu[i] - P_of_C_mu[i - 1]);
                                           target += (prior_dirichlet_alpha[i] - 1.0) *
                                                     (log(P_of_C_mu_of_the_category) - log(dirichlet_phi[i]));
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
                          target += n_chosen_cutoffs[p_class] * log(scale_of_line[p_class]);
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
                          target += log(4.0) + 0.5 * (n_thr - 1) *
                                    log(dot_self(C_mu - mean(C_mu)));
                     } else if (joint_link_prior == 1) {
                          //// Link-space volume of both lines with respect to C_mu, location and raw scale.
                          real a = scale_of_line[1];
                          real d = scale_of_line[2];
                          real da = a;
                          real dd = -d;
                          if (softplus == 1) {
                               da = 1.4426950408889634 * inv_logit(Xraw_scale_baseline_c);
                               dd = -1.4426950408889634 * inv_logit(-Xraw_scale_baseline_c);
                          }
                          target += log(2.0 * a * d * (da * d - a * dd)) +
                                    0.5 * log(n_chosen_cutoffs[p_class] *
                                              dot_self(C_mu - mean(C_mu))) +
                                    0.5 * (n_chosen_cutoffs[p_class] - 2) *
                                    log(square(a) + square(d));
                     } else {
                          target += log(2.0 * scale_of_line[q_class]);
                     }
                  }
                  if (n_covariates > 1) {
                       beta_mu[2:n_covariates] ~ normal( prior_beta_mu_mean[2:n_covariates],
                                                         prior_beta_mu_SD[2:n_covariates]);
                  }
                  if (joint_unit_probability_prior == 1 || joint_link_prior == 1) {
                       if (n_covariates > 1) {
                            raw_scale_mu[2:n_covariates] ~ normal(
                                  prior_raw_scale_mu_mean[2:n_covariates],
                                  prior_raw_scale_mu_SD[2:n_covariates]);
                       }
                  } else {
                       raw_scale_mu[1:n_covariates] ~ normal(
                             prior_raw_scale_mu_mean[1:n_covariates],
                             prior_raw_scale_mu_SD[1:n_covariates]);
                  }
                  
           } else {
             
                  beta_mu[1:n_covariates] ~ normal( prior_beta_mu_mean[1:n_covariates],
                                                         prior_beta_mu_SD[1:n_covariates]);
                  raw_scale_mu[1:n_covariates] ~ normal( prior_raw_scale_mu_mean[1:n_covariates],
                                                              prior_raw_scale_mu_SD[1:n_covariates]);
                  ////
                  dirichlet_phi ~ dirichlet( prior_dirichlet_alpha);
                  
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
        //// (anchored at the prior mean of the mid-location)
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
                                                                                  prior_beta_mu_mean[1] |
                                                                                  alpha_of_each_segment[(p_1 + s - 1):(p_2 + s)],
                                                                                  use_probit_link);
                  }
             }
        }
        ////
        //// ---- Likelihood / Model:
        ////
        target += normal_lupdf(to_vector(beta_z) | 0, 1); // part of between-study model, NOT prior
        target += normal_lupdf(to_vector(raw_scale_z) | 0, 1); // part of between-study model, NOT prior
        ////
        //// (a study's first reported cutpoint is C_raw itself; each later one adds an exp() or softplus step)
        target += (softplus == 1) ? sum(log_inv_logit(C_raw[position_of_each_increment])) :
                                    sum(C_raw[position_of_each_increment]);
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
                               1.4426950408889634 * log1p_exp(-mult_nd * Xraw_scale_baseline) :
                                                  exp(-mult_nd * Xraw_scale_baseline));
      real scale_d_baseline  = ((softplus == 1) ? 
                               1.4426950408889634 * log1p_exp(-mult_d  * Xraw_scale_baseline) :
                                                  exp(-mult_d  * Xraw_scale_baseline));
      ////
      //// ---- Calculate baseline Se/Sp:
      ////
      vector[n_thr] Fp_baseline = (use_probit_link == 1) ? 
                                  Phi(-(C_mu + Xbeta_baseline_nd)*scale_nd_baseline)  : 
                                  inv_logit(-(C_mu + Xbeta_baseline_nd)*scale_nd_baseline);
      vector[n_thr] Sp_baseline = 1.0 - Fp_baseline;
      vector[n_thr] Se_baseline = (use_probit_link == 1) ? 
                                  Phi(-(C_mu + Xbeta_baseline_d)*scale_d_baseline)    : 
                                  inv_logit(-(C_mu + Xbeta_baseline_d)*scale_d_baseline);
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
                                    1.4426950408889634 * log1p_exp(-mult_nd * Xraw_scale_baseline_pred) :
                                                       exp(-mult_nd * Xraw_scale_baseline_pred));
      real scale_baseline_pred_d  = ((softplus == 1) ? 
                                    1.4426950408889634 * log1p_exp(-mult_d  * Xraw_scale_baseline_pred) :
                                                       exp(-mult_d  * Xraw_scale_baseline_pred));
      ////
      vector[n_thr] Fp_baseline_pred = (use_probit_link == 1) ? 
                                       Phi(-(C_pred + Xbeta_baseline_pred_nd)*scale_baseline_pred_nd) : 
                                       inv_logit(-(C_pred + Xbeta_baseline_pred_nd)*scale_baseline_pred_nd);
      vector[n_thr] Sp_baseline_pred = 1.0 - Fp_baseline_pred;
      vector[n_thr] Se_baseline_pred = (use_probit_link == 1) ? 
                                       Phi(-(C_pred + Xbeta_baseline_pred_d)*scale_baseline_pred_d)   : 
                                       inv_logit(-(C_pred + Xbeta_baseline_pred_d)*scale_baseline_pred_d);
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

