//// DTA_MA_ord_HSROC_FIXEDthr_NicoStan_fused_AVX.stan
//// The model of MetaOrdDTA's DTA_MA_ord_HSROC_FIXEDthr.stan with NicoSIMD's AVX functions:
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
        real nico_log1p_exp(real x);
        vector nico_log_inv_logit(vector x);
        vector nico_log1m_inv_logit(vector x);
        real nico_log_det_jacobian_of_cutpoints_from_raw_differences(vector C_raw, data int softplus);
        real nico_log_det_jacobian_of_cutpoints_from_raw_differences(matrix C_raw, data int softplus);
        array[] matrix nico_map_latent_surv_prob_to_fixed_C_HSROC(data vector C,
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
        real nico_normal_lpdf(real y, real mu, real sigma);
        real nico_log(real x);
        real nico_exp(real x);
        real nico_inv_logit(real x);
        real nico_log_inv_logit(real x);
        real nico_log1m_inv_logit(real x);
        vector nico_log1m(vector x);
        real nico_dirichlet_lpdf(vector theta, vector alpha);
        real nico_beta_lpdf(vector y, vector alpha, vector beta);
        
        ////
        //// ---- Induced Dirichlet with the Dirichlet of the segment sums (between chosen cutoffs) swapped for
        //// ---- a prior on baseline Sp (c = 1) or 1 - Se (c = 2) at the chosen cutoffs:
        ////
        real induced_dirichlet_with_Se_Sp_prior_lpdf(vector C_shifted,
                                                     vector C_baseline,
                                                     vector alpha,
                                                     vector alpha_seg,
                                                     array[] int chosen,
                                                     int prior_type,
                                                     vector A_user,
                                                     matrix beta_ab,
                                                     matrix link_mean_SD,
                                                     int use_probit_link) {

                int n_thr = num_elements(C_shifted);
                int m = num_elements(chosen);
                vector[n_thr] P = (use_probit_link == 1) ? nico_Phi(C_shifted) : nico_inv_logit(C_shifted);
                vector[m] P_chosen = (use_probit_link == 1) ? nico_Phi(C_baseline[chosen]) :
                                                            nico_inv_logit(C_baseline[chosen]);
                vector[m + 1] q = cumul_probs_to_ord_probs(P_chosen);
                vector[m + 1] q_anchor = cumul_probs_to_ord_probs(P[chosen]);
                ////
                //// ---- induced Dirichlet minus the Dirichlet of the segment sums:
                ////
                real log_prob = nico_induced_dirichlet_prior_of_ordered_cutpoints_lpdf(
                                C_shifted | alpha, use_probit_link, 0.0);
                log_prob += - nico_dirichlet_lpdf(q_anchor | alpha_seg);
                ////
                //// ---- Jacobian between selected anchor and baseline probabilities:
                ////
                for (j in 1:m) {
                     int k = chosen[j];
                     if (use_probit_link == 1) {
                          log_prob += - nico_normal_lpdf(C_shifted[k] | 0, 1);
                          if (prior_type != 4) {
                               log_prob += nico_normal_lpdf(C_baseline[k] | 0, 1);
                          }
                     } else {
                          log_prob += - nico_log_inv_logit(C_shifted[k]) - nico_log1m_inv_logit(C_shifted[k]);
                          if (prior_type != 4) {
                               log_prob += nico_log_inv_logit(C_baseline[k]) +
                                           nico_log1m_inv_logit(C_baseline[k]);
                          }
                     }
                }
                ////
                //// ---- prior on P at the chosen cutoffs:
                ////
                if (prior_type == 1) {
                     //// Dirichlet on the segment sums:
                     log_prob += nico_dirichlet_lpdf(q | A_user);
                } else if (prior_type == 2) {
                     //// Beta on each P, restricted to rising P (normaliser constant):
                     log_prob += nico_beta_lpdf(P_chosen | beta_ab[, 1], beta_ab[, 2]);
                } else if (prior_type == 4) {
                     //// Normal on the link scale, restricted to rising P:
                     log_prob += nico_normal_lpdf(C_baseline[chosen] | link_mean_SD[, 1], link_mean_SD[, 2]);
                } else {
                     //// Beta on each (P_j - P_{j-1}) / (1 - P_{j-1}):
                     vector[m] P_prev = append_row(0.0, P_chosen[1:(m - 1)]);
                     log_prob += nico_beta_lpdf(q[1:m] ./ (1.0 - P_prev) | beta_ab[, 1], beta_ab[, 2]);
                     log_prob += - sum(nico_log1m(P_prev));
                }
                return log_prob;

        }
        ////
        //// ---- Prior on baseline Sp (c = 1) or 1 - Se (c = 2) at the chosen cutoffs, on its own (no cutpoint
        //// ---- law swapped; used for the class whose prior goes on the baseline location and raw scale):
        ////
        real Se_Sp_prior_at_chosen_cutoffs_lpdf(vector C_baseline_chosen,
                                                int prior_type,
                                                vector A_user,
                                                matrix beta_ab,
                                                matrix link_mean_SD,
                                                int use_probit_link) {

                int m = num_elements(C_baseline_chosen);
                vector[m] P_chosen = (use_probit_link == 1) ? nico_Phi(C_baseline_chosen) :
                                                            nico_inv_logit(C_baseline_chosen);
                vector[m + 1] q = cumul_probs_to_ord_probs(P_chosen);
                real log_prob = 0.0;
                ////
                //// ---- Jacobian between the probabilities and the baseline line:
                ////
                if (prior_type != 4) {
                     if (use_probit_link == 1) {
                          log_prob += nico_normal_lpdf(C_baseline_chosen | 0, 1);
                     } else {
                          log_prob += sum(nico_log_inv_logit(C_baseline_chosen) +
                                          nico_log1m_inv_logit(C_baseline_chosen));
                     }
                }
                ////
                //// ---- prior on P at the chosen cutoffs:
                ////
                if (prior_type == 1) {
                     //// Dirichlet on the segment sums:
                     log_prob += nico_dirichlet_lpdf(q | A_user);
                } else if (prior_type == 2) {
                     //// Beta on each P, restricted to rising P (normaliser constant):
                     log_prob += nico_beta_lpdf(P_chosen | beta_ab[, 1], beta_ab[, 2]);
                } else if (prior_type == 4) {
                     //// Normal on the link scale, restricted to rising P:
                     log_prob += nico_normal_lpdf(C_baseline_chosen | link_mean_SD[, 1], link_mean_SD[, 2]);
                } else {
                     //// Beta on each (P_j - P_{j-1}) / (1 - P_{j-1}):
                     vector[m] P_prev = append_row(0.0, P_chosen[1:(m - 1)]);
                     log_prob += nico_beta_lpdf(q[1:m] ./ (1.0 - P_prev) | beta_ab[, 1], beta_ab[, 2]);
                     log_prob += - sum(nico_log1m(P_prev));
                }
                return log_prob;

        }
} 


data {
        ////
        //// ---- Data:
        //// 
        int<lower=1> n_studies;
        int<lower=1> n_thr;
        array[n_studies] int n_obs_cutpoints;  //// OBSERVED cutpoints for test in study s
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
        vector<lower=0.0>[n_thr + 1] prior_dirichlet_alpha; // Induced-Dirichlet prior vector 
        ////
        //// ---- Priors directly on Sp (c = 1) and Se (c = 2) at chosen cutoffs:
        //// ---- Types 1 and 3 use P = Sp and 1 - Se; types 2 and 4 use Sp and Se:
        //// ---- The class with more chosen cutoffs (Sp on a tie) goes on the cutpoints; the other class at <= 2
        //// ---- of the same cutoffs replaces the prior on the baseline location (1 cutoff) and raw scale (2):
        ////
        int<lower=0, upper=1> use_Se_Sp_prior;     //// 0 = original cutpoint prior
        int<lower=0, upper=1> anchor_at_location;  //// 1 = cutpoints shifted and scaled at the baseline line
        int<lower=1, upper=4> Se_Sp_prior_type;    //// 1 = Dirichlet, 2-3 = Betas, 4 = link-scale Normal
        int<lower=0, upper=n_thr> n_chosen_max;
        array[2] int<lower=0, upper=n_chosen_max> n_chosen;
        array[2, n_chosen_max] int<lower=1, upper=n_thr> chosen_cutoffs;
        array[2] vector<lower=0.0>[n_chosen_max + 1] prior_Se_Sp_dirichlet_A;
        array[2] matrix<lower=0.0>[n_chosen_max, 2] prior_Se_Sp_beta;
        array[2] matrix[n_chosen_max, 2] prior_Se_Sp_link_mean_SD;  //// Mean and SD of link(Sp), link(Se)
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
                                      index_of_the_cutpoint_of_each_term[t] = k;
                                }
                          }
                    }
              }
        }
        ////
        //// ---- Sums of prior_dirichlet_alpha over the segments between chosen cutoffs:
        ////
        array[2] vector[n_chosen_max + 1] alpha_seg;
        for (c in 1:2) {
              alpha_seg[c] = rep_vector(0.0, n_chosen_max + 1);
              int k_prev = 0;
              for (j in 1:n_chosen[c]) {
                   int k_j = chosen_cutoffs[c, j];
                   if (k_j <= k_prev) reject("chosen_cutoffs must rise: c = ", c);
                   if (use_Se_Sp_prior == 1 && Se_Sp_prior_type == 4 &&
                       prior_Se_Sp_link_mean_SD[c][j, 2] <= 0.0) {
                        reject("The link-scale prior SD must be positive: c = ", c, ", j = ", j);
                   }
                   alpha_seg[c][j] = sum(prior_dirichlet_alpha[(k_prev + 1):k_j]);
                   k_prev = k_j;
              }
              alpha_seg[c][n_chosen[c] + 1] = sum(prior_dirichlet_alpha[(k_prev + 1):n_cat]);
        }
        ////
        //// ---- Class on the cutpoints (p_class) and class on the baseline location / raw scale (q_class):
        ////
        int p_class = (n_chosen[2] > n_chosen[1]) ? 2 : 1;
        int q_class = 3 - p_class;
        int n_q = (use_Se_Sp_prior == 1) ? n_chosen[q_class] : 0;
        if (n_q > 2) reject("At most 2 chosen cutoffs for the class with fewer: c = ", q_class);
        for (j in 1:n_q) {
             int found = 0;
             for (i in 1:n_chosen[p_class]) {
                  if (chosen_cutoffs[p_class, i] == chosen_cutoffs[q_class, j]) found = 1;
             }
             if (found == 0) reject("chosen_cutoffs of c = ", q_class, " must be among those of c = ", p_class);
        }
        if (n_q > 0 && baseline_case[1] == 0.0) reject("baseline_case[1] must not be 0");
        int first_beta_mu = (n_q >= 1) ? 2 : 1;
        int first_raw_scale_mu = (n_q == 2) ? 2 : 1;
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
        vector<lower=-15, upper=2.5>[n_thr] C_raw_vec; 
}
 

transformed parameters { 
        ////
        //// ---- Construct (global) cutpoints:
        ////
        // vector[n_thr] C = construct_C(C_raw_vec, softplus);
        //// (as construct_C: the first cutpoint, then the cumulative sums of the exp() (or softplus) of
        ////  the raw differences, with NicoSIMD's nico_exp / nico_log1p_exp)
        vector[n_thr] C;
        C[1] = C_raw_vec[1];
        if (n_thr > 1) {
              vector[n_thr - 1] increments_of_the_cutpoints =
                    (softplus == 1) ? nico_log1p_exp(C_raw_vec[2:n_thr]) :
                                      nico_exp(C_raw_vec[2:n_thr]);
              for (k in 2:n_thr) {
                    C[k] = C[k - 1] + increments_of_the_cutpoints[k - 1];
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
            
            vector[n_studies] Xraw_scale_baseline = X[, 1:n_covariates] * to_vector(raw_scale_mu[1:n_covariates]) + raw_scale_random;
            // Xscale[, 1] =  ((softplus == 1) ? softplus_scaled(-mult_nd*Xraw_scale_baseline) : exp(-mult_nd*Xraw_scale_baseline));
            Xscale[, 1] =  ((softplus == 1) ? 1.4426950408889634 * nico_log1p_exp(-mult_nd*Xraw_scale_baseline) :
                                                nico_exp(-mult_nd*Xraw_scale_baseline)); 
            // Xscale[, 2] =  ((softplus == 1) ? softplus_scaled(-mult_d*Xraw_scale_baseline)  : exp(-mult_d*Xraw_scale_baseline));
            Xscale[, 2] =  ((softplus == 1) ? 1.4426950408889634 * nico_log1p_exp(-mult_d*Xraw_scale_baseline)  :
                                                nico_exp(-mult_d*Xraw_scale_baseline));
        }
}


model {
        ////
        //// ---- Priors for locations:
        ////
        // beta_mu[1:n_covariates] ~ normal(
        //                           prior_beta_mu_mean[1:n_covariates], 
        //                           prior_beta_mu_SD[1:n_covariates]);  
        if (first_beta_mu <= n_covariates) {
             beta_mu[first_beta_mu:n_covariates] ~ nico_normal(prior_beta_mu_mean[first_beta_mu:n_covariates],
                                                                prior_beta_mu_SD[first_beta_mu:n_covariates]);
        }
        beta_SD ~ nico_normal(
                  prior_beta_SD_mean, 
                  prior_beta_SD_SD);
        ////
        // raw_scale_mu[1:n_covariates] ~ normal(
        //                                prior_raw_scale_mu_mean[1:n_covariates], 
        //                                prior_raw_scale_mu_SD[1:n_covariates]); 
        if (first_raw_scale_mu <= n_covariates) {
             raw_scale_mu[first_raw_scale_mu:n_covariates] ~ nico_normal(
                                  prior_raw_scale_mu_mean[first_raw_scale_mu:n_covariates],
                                  prior_raw_scale_mu_SD[first_raw_scale_mu:n_covariates]);
        }
        raw_scale_SD ~ nico_normal(
                       prior_raw_scale_SD_mean, 
                       prior_raw_scale_SD_SD);
        ////
        //// ---- Induced-dirichlet ** Prior ** model:
        ////
        // {
        //    vector[n_cat] Ind_Dir_ord_prob; //// = cumul_probs_to_ord_probs(Phi(C));
        //    if (use_probit_link == 1) Ind_Dir_ord_prob = cumul_probs_to_ord_probs(Phi(C));
        //    else                      Ind_Dir_ord_prob = cumul_probs_to_ord_probs(inv_logit(C));
        //    ////
        //    Ind_Dir_ord_prob ~ induced_dirichlet_given_C(
        //                       C, prior_dirichlet_alpha, use_probit_link);
        // }

        if (use_Se_Sp_prior == 1 && n_chosen[p_class] > 0) {
           ////
           //// ---- Baseline lines: probit(Sp) (c = 1) and probit(1 - Se) (c = 2) = (C + Xbeta) * scale:
           ////
           real Xbeta_baseline_c = dot_product(baseline_case, beta_mu);
           real Xraw_scale_baseline_c = dot_product(baseline_case, raw_scale_mu);
           vector[2] Xbeta_of_line = [mult_nd, mult_d]' * Xbeta_baseline_c;
           vector[2] scale_of_line;
           scale_of_line[1] = ((softplus == 1) ? 1.4426950408889634 *
                                                 nico_log1p_exp(-mult_nd * Xraw_scale_baseline_c) :
                                                 nico_exp(-mult_nd * Xraw_scale_baseline_c));
           scale_of_line[2] = ((softplus == 1) ? 1.4426950408889634 *
                                                 nico_log1p_exp(-mult_d  * Xraw_scale_baseline_c) :
                                                 nico_exp(-mult_d  * Xraw_scale_baseline_c));
           array[2] vector[n_thr] C_baseline;
           for (c in 1:2) {
                C_baseline[c] = (C + Xbeta_of_line[c]) * scale_of_line[c];
           }
           vector[n_thr] C_shifted = (anchor_at_location == 1) ? C_baseline[p_class] : C;
           ////
           //// ---- Class p_class: induced Dirichlet with its prior at the chosen cutoffs:
           ////
           {
              int m = n_chosen[p_class];
              matrix[m, 2] beta_ab = prior_Se_Sp_beta[p_class][1:m, ];
              matrix[m, 2] link_mean_SD = prior_Se_Sp_link_mean_SD[p_class][1:m, ];
              if (p_class == 2) {
                   link_mean_SD[, 1] = - link_mean_SD[, 1];
                   if (Se_Sp_prior_type == 2) {
                        beta_ab = append_col(beta_ab[, 2], beta_ab[, 1]);
                   }
              }
              target += induced_dirichlet_with_Se_Sp_prior_lpdf(
                        C_shifted | C_baseline[p_class], prior_dirichlet_alpha,
                        alpha_seg[p_class][1:(m + 1)],
                        chosen_cutoffs[p_class, 1:m], Se_Sp_prior_type,
                        prior_Se_Sp_dirichlet_A[p_class][1:(m + 1)], beta_ab,
                        link_mean_SD,
                        use_probit_link);
              ////
              //// ---- Jacobian from C to the scaled line (all cutpoints when shifted there, else the chosen):
              ////
              target += ((anchor_at_location == 1) ? n_thr : m) * nico_log(scale_of_line[p_class]);
           }
           ////
           //// ---- Class q_class at <= 2 of the same cutoffs, in place of the priors on the baseline location
           //// ---- (1 cutoff) and raw scale (2 cutoffs):
           ////
           if (n_q > 0) {
              matrix[n_q, 2] beta_ab = prior_Se_Sp_beta[q_class][1:n_q, ];
              matrix[n_q, 2] link_mean_SD = prior_Se_Sp_link_mean_SD[q_class][1:n_q, ];
              if (q_class == 2) {
                   link_mean_SD[, 1] = - link_mean_SD[, 1];
                   if (Se_Sp_prior_type == 2) {
                        beta_ab = append_col(beta_ab[, 2], beta_ab[, 1]);
                   }
              }
              target += Se_Sp_prior_at_chosen_cutoffs_lpdf(
                        C_baseline[q_class][chosen_cutoffs[q_class, 1:n_q]] | Se_Sp_prior_type,
                        prior_Se_Sp_dirichlet_A[q_class][1:(n_q + 1)], beta_ab,
                        link_mean_SD,
                        use_probit_link);
              ////
              //// ---- Jacobian from beta_mu[1] (and raw_scale_mu[1]) to line q_class, line p_class held:
              ////
              target += nico_log(2.0) + nico_log(scale_of_line[q_class]) + nico_log(abs(baseline_case[1]));
              if (n_q == 2) {
                   int k_a = chosen_cutoffs[q_class, 1];
                   int k_b = chosen_cutoffs[q_class, 2];
                   real r = Xraw_scale_baseline_c;
                   real d_of_scales = 2.0;
                   if (softplus == 1) {
                        d_of_scales = nico_inv_logit(r) * 1.4426950408889634 * nico_log1p_exp(-r) +
                                      1.4426950408889634 * nico_log1p_exp(r) * nico_inv_logit(-r);
                        d_of_scales = d_of_scales / nico_log(2.0);
                   }
                   target += nico_log(C[k_b] - C[k_a]) + nico_log(d_of_scales) - nico_log(scale_of_line[p_class]);
                   target += nico_log(abs(baseline_case[1]));
              }
           }
        } else {
           //// (as induced_dirichlet_given_C)
           target += nico_induced_dirichlet_prior_of_ordered_cutpoints_lupdf(C | prior_dirichlet_alpha,
                                                                             use_probit_link, 0.0);
        }
        ////
        //// ---- Likelihood / Model:
        ////
        // target += std_normal_lpdf(to_vector(beta_z)); // part of between-study model, NOT prior

        target += nico_normal_lupdf(to_vector(beta_z) | 0, 1); // part of between-study model, NOT prior
        // target += std_normal_lpdf(to_vector(raw_scale_z)); // part of between-study model, NOT prior

        target += nico_normal_lupdf(to_vector(raw_scale_z) | 0, 1); // part of between-study model, NOT prior
        // target += raw_scale_to_scale_log_det_J_lp(raw_scale, softplus); // no Jacobian needed for R&G models as it cancels out (since + then -).
        // target += raw_C_to_C_log_det_J_lp(C_raw_vec, softplus);
        target += nico_log_det_jacobian_of_cutpoints_from_raw_differences(C_raw_vec, softplus);
        ////
        //// ---- Increment the log-likelihood:
        ////
        // {
        //     ////
        //     //// ---- Get the cutpoint index (k) to map "latent_surv[c][s, cut_i]" to correct cutpoint "C[k]":
        //     ////
        //     array[2] matrix[n_studies, n_thr] latent_surv = map_latent_surv_prob_to_fixed_C_HSROC(
        //                                                     C, Xlocations, Xscale, n_studies, 
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
            vector[2 * n_studies] locations_of_both_classes = to_vector(Xlocations);
            vector[2 * n_studies] scales_of_both_classes = to_vector(Xscale);
            //// (as map_latent_surv_prob_to_fixed_C_HSROC: -(C + location) * scale)
            vector[n_terms] latent_of_each_term =
                  -(C[index_of_the_cutpoint_of_each_term] +
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
      //// ---- Calculate summary accuracy (using mean parameters):
      ////
      // vector[n_covariates] location_nd_mu = mult_nd*beta_mu;
      // vector[n_covariates] location_d_mu  = mult_d*beta_mu;
      // ////
      // vector[n_covariates] scale_nd_mu =  ((softplus == 1) ? softplus_scaled(mult_nd*raw_scale_mu) : exp(mult_nd*raw_scale_mu));
      // vector[n_covariates] scale_d_mu  =  ((softplus == 1) ? softplus_scaled(mult_d*raw_scale_mu)  : exp(mult_d*raw_scale_mu));
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
                                  nico_Phi(-(C + Xbeta_baseline_nd)*scale_nd_baseline) :
                                  nico_inv_logit(-(C + Xbeta_baseline_nd)*scale_nd_baseline);
      vector[n_thr] Sp_baseline = 1.0 - Fp_baseline;
      vector[n_thr] Se_baseline = (use_probit_link == 1) ? 
                                  nico_Phi(-(C + Xbeta_baseline_d)*scale_d_baseline) : 
                                  nico_inv_logit(-(C + Xbeta_baseline_d)*scale_d_baseline);
      ////
      //// ---- Calculate predictive accuracy:
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
                                    1.4426950408889634 * nico_log1p_exp(-mult_nd * Xraw_scale_baseline_pred) :
                                    nico_exp(-mult_nd * Xraw_scale_baseline_pred));
      real scale_baseline_pred_d  = ((softplus == 1) ? 
                                    1.4426950408889634 * nico_log1p_exp(-mult_d  * Xraw_scale_baseline_pred) :
                                    nico_exp(-mult_d  * Xraw_scale_baseline_pred));
      ////
      vector[n_thr] Fp_baseline_pred = (use_probit_link == 1) ? 
                                       nico_Phi(-(C + Xbeta_baseline_pred_nd)*scale_baseline_pred_nd) : 
                                       nico_inv_logit(-(C + Xbeta_baseline_pred_nd)*scale_baseline_pred_nd);
      vector[n_thr] Sp_baseline_pred = 1.0 - Fp_baseline_pred;
      vector[n_thr] Se_baseline_pred = (use_probit_link == 1) ? 
                                       nico_Phi(-(C + Xbeta_baseline_pred_d)*scale_baseline_pred_d)   : 
                                       nico_inv_logit(-(C + Xbeta_baseline_pred_d)*scale_baseline_pred_d);
      ////
      //// ---- Log-lik + study-specific accuracy computation (using "data" / double-precision fn for efficiency):
      ////
      vector[n_studies] log_lik_study = rep_vector(0.0, n_studies);
      ////
      matrix[n_studies, n_thr] fp; // global
      matrix[n_studies, n_thr] sp; // global
      matrix[n_studies, n_thr] se; // global
      ////
      vector[n_studies] deviance_nd; // global
      vector[n_studies] deviance_d;  // global
      vector[n_studies] deviance;    // global
      {
          real scale = 1.0; // since using "ord_bivariate-like"" param.
          array[2] matrix[n_studies, n_thr] latent_surv = nico_map_latent_surv_prob_to_fixed_C_HSROC(
                                                          C, Xlocations, Xscale, n_studies, 
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
          fp = surv_prob[1];
          sp = 1.0 - fp;
          se = surv_prob[2];
          ////
          //// ---- Model fit (deviance):
          ////
          array[4] matrix[n_studies, n_thr] outs_model_fit = compute_deviance(
                                                             cond_prob, n_thr, n_studies, 
                                                             x_2, n, n_obs_cutpoints);
          ////
          matrix[n_studies, n_thr] x_hat_nd = outs_model_fit[1];
          matrix[n_studies, n_thr] dev_nd   = outs_model_fit[2];
          matrix[n_studies, n_thr] x_hat_d  = outs_model_fit[3];
          matrix[n_studies, n_thr] dev_d    = outs_model_fit[4];
          ////
          for (s in 1:n_studies) {
              for (k in 1:to_int(n_obs_cutpoints[s])) {
                 deviance_nd[s] += dev_nd[s, k];
                 deviance_d[s]  += dev_d[s, k];
              }
              deviance[s] = deviance_nd[s] + deviance_d[s];
          }
      }
      // ////
      // //// ---- "Ordinal-bivariate equivelent" params:
      // ////
      // array[2] matrix[n_thr, n_covariates] biv_equiv_C;
      // for (x_i in 1:n_covariates) {
      //     biv_equiv_C[1][, x_i] = to_vector(C ./ scale_nd_mu[x_i]);
      //     biv_equiv_C[2][, x_i] = to_vector(C ./ scale_d_mu[x_i]);
      // }
}













