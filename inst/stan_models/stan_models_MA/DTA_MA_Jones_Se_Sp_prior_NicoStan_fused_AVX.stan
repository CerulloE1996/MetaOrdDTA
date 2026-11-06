//// DTA_MA_Jones_Se_Sp_prior_NicoStan_fused_AVX.stan
//// The model of MetaOrdDTA's DTA_MA_Jones_Se_Sp_prior.stan with NicoSIMD's AVX functions:
//// the likelihood in one autodiff node, the Box-Cox cutpoints and the scales on nico_exp / nico_log1p_exp.
////



functions {
        ////
        //// ---- Include files to compile the necessary custom (user-defined) Stan functions:
        ////
        #include "Stan_fns_basic.stan"
        #include "Stan_fns_Box_Cox.stan"
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
        real nico_exp(real x);
        vector nico_exp(vector x);
        matrix nico_exp(matrix x);
        real nico_log1p_exp(real x);
        matrix nico_log1p_exp(matrix x);
        real nico_log(real x);
        vector nico_log(vector x);
        real nico_log_inv_logit(real x);
        vector nico_log_inv_logit(vector x);
        vector nico_log1m_inv_logit(vector x);
        vector nico_Phi(vector x);
        vector nico_inv_logit(vector x);
        array[] matrix nico_map_latent_surv_to_fixed_hetero_C(data array[] vector C,
                                                              data matrix locations, data matrix scales,
                                                              data int n_studies, data array[] int n_obs_cutpoints,
                                                              data array[,,] int cutpoint_index);
        array[,] matrix nico_compute_log_lik_binomial_fact_data(data array[] matrix latent_surv,
                                                                data int use_probit_link, data int n_thr,
                                                                data int n_studies, data array[,,] int x_2,
                                                                data array[,,] int n, data array[,] int N_total,
                                                                data array[] int n_obs_cutpoints,
                                                                data array[] int K_fold_CV_indicator);
        real nico_lkj_corr_cholesky_lpdf(matrix L, data real eta);
        real nico_beta_lpdf(vector y, vector alpha, vector beta);
        real nico_normal_lpdf(vector y, real mu, real sigma);
        real nico_normal_lpdf(vector y, vector mu, vector sigma);
        real nico_normal_lpdf(real y, real mu, real sigma);
        ////
        //// ---- Vectorized Box-Cox transform function (fn_Stan_box_cox of Stan_fns_Box_Cox.stan, with nico_log):
        ////
        vector fn_nico_box_cox( vector x, 
                                real lambda) {
              
              int N = num_elements(x);
              vector[N] result;
              
              for (n in 1:N) {
                if (lambda == 0.0) {
                    if (x[n] != 0.0) {
                      result[n] = nico_log(x[n]);
                    } else {
                      result[n] = -700.0;
                    }
                } else {
                    result[n] = (pow(x[n], lambda) - 1.0) / lambda;
                }
              }
              
              return result;
              
        }
        ////
        //// ---- Probability-scale Beta or link-scale Normal at chosen cutoffs:
        ////
        real Se_Sp_prior_at_chosen_cutoffs_lpdf(vector C_baseline_chosen,
                                                int prior_type,
                                                matrix beta_ab,
                                                matrix link_mean_SD,
                                                int use_probit_link) {

                real log_prob = 0.0;
                //// Selected baseline values follow the ordered curve.
                if (prior_type == 2) {
                     vector[num_elements(C_baseline_chosen)] P_chosen = (use_probit_link == 1) ?
                           nico_Phi(C_baseline_chosen) : nico_inv_logit(C_baseline_chosen);
                     log_prob += nico_beta_lpdf(P_chosen | beta_ab[, 1], beta_ab[, 2]);
                     if (use_probit_link == 1) {
                          log_prob += nico_normal_lpdf(C_baseline_chosen | 0, 1);
                     } else {
                          log_prob += sum(nico_log_inv_logit(C_baseline_chosen) +
                                          nico_log1m_inv_logit(C_baseline_chosen));
                     }
                } else {
                     log_prob += nico_normal_lpdf(C_baseline_chosen | link_mean_SD[, 1], link_mean_SD[, 2]);
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
        array[n_studies] int n_obs_cutpoints;
        ////
        array[2, n_studies, n_thr + 1] int x;
        array[2, n_studies, n_thr + 1] int cutpoint_index;
        ////
        vector[n_thr] cts_thr_values;
        ////
        //// ---- Covariates:
        ////
        int n_covariates_nd;
        int n_covariates_d;
        int n_covariates_max;
        matrix[n_studies, n_covariates_nd] X_nd;
        matrix[n_studies, n_covariates_d]  X_d;
        vector[n_covariates_nd] baseline_case_nd;  // must be user-inputted - e.g. could be {0, 1, 45.3} for 2 binary covariates and 1 cts one (e.g. age)
        vector[n_covariates_d]  baseline_case_d;   // must be user-inputted - e.g. could be {0, 1, 45.3} for 2 binary covariates and 1 cts one (e.g. age)
        ////
        //// ---- Priors for locations ("beta"):
        ////
        matrix[2, n_covariates_max] prior_beta_mu_mean;
        matrix[2, n_covariates_max] prior_beta_mu_SD;
        vector[2] prior_beta_SD_mean;
        vector[2] prior_beta_SD_SD;
        ////
        //// ---- Priors for raw scales:
        ////
        matrix[2, n_covariates_max] prior_raw_scale_mu_mean;
        matrix[2, n_covariates_max] prior_raw_scale_mu_SD;
        vector[2] prior_raw_scale_SD_mean;
        vector[2] prior_raw_scale_SD_SD;
        ////
        //// ---- Priors for box-cox ("lambda"):
        ////
        vector[2] prior_boxcox_lambda_mean;
        vector[2] prior_boxcox_lambda_SD;
        ////
        //// ---- Priors (and possible restrictions) for between-study correlations:
        ////
        real<lower=-1.0, upper=+1.0> beta_corr_lb;
        real<lower=beta_corr_lb, upper=+1.0>  beta_corr_ub;
        real<lower=0.0>  prior_beta_corr_LKJ;
        ////
        real<lower=-1.0, upper=+1.0> raw_scale_corr_lb;
        real<lower=raw_scale_corr_lb, upper=+1.0>  raw_scale_corr_ub;
        real<lower=0.0>  prior_raw_scale_corr_LKJ;
        ////
        ////
        //// ---- Direct baseline Sp (class 1) and Se (class 2) priors:
        ////
        int<lower=0, upper=1> use_Se_Sp_prior;
        int<lower=1, upper=4> Se_Sp_prior_type;
        int<lower=0, upper=n_thr> n_chosen_max;
        array[2] int<lower=0, upper=n_chosen_max> n_chosen;
        array[2, n_chosen_max] int<lower=1, upper=n_thr> chosen_cutoffs;
        array[2] matrix<lower=0.0>[n_chosen_max, 2] prior_Se_Sp_beta;
        array[2] matrix[n_chosen_max, 2] prior_Se_Sp_link_mean_SD;
        ////
        //// ---- Other:
        ////
        int<lower=0, upper=1> box_cox;
        int<lower=0, upper=1> softplus;
        ////
        int use_probit_link;
        ////
        array[n_studies] int<lower=0, upper=1> K_fold_CV_indicator;

}

transformed data {
        ////
        //// ---- One cutoff replaces location; multiple cutoffs also replace scale:
        ////
        if (use_Se_Sp_prior == 1) {
             if (Se_Sp_prior_type != 2 && Se_Sp_prior_type != 4) {
                  reject("Jones direct priors support type 2 (Beta) or type 4 (link-scale Normal)");
             }
             if (use_probit_link != 0 && use_probit_link != 1) reject("use_probit_link must be 0 or 1");
             for (c in 1:2) {
                  if (n_chosen[c] > 0) {
                       if (c == 1 && baseline_case_nd[1] == 0.0) reject("baseline_case_nd[1] must not be 0");
                       if (c == 2 && baseline_case_d[1] == 0.0) reject("baseline_case_d[1] must not be 0");
                       for (j in 1:n_chosen[c]) {
                            if (j > 1) {
                                 if (chosen_cutoffs[c, j] <= chosen_cutoffs[c, j - 1]) {
                                      reject("Chosen cutoffs must rise");
                                 }
                                 if (cts_thr_values[chosen_cutoffs[c, j]] <=
                                     cts_thr_values[chosen_cutoffs[c, j - 1]]) {
                                      reject("Chosen threshold values must rise");
                                 }
                            }
                            if (Se_Sp_prior_type == 2 &&
                                min(prior_Se_Sp_beta[c][j, ]) <= 0.0) {
                                 reject("Beta shapes must be positive");
                            }
                            if (Se_Sp_prior_type == 4 &&
                                prior_Se_Sp_link_mean_SD[c][j, 2] <= 0.0) {
                                 reject("The link-scale prior SD must be positive");
                            }
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
        //// ---- The terms of the likelihood, one per reported cutpoint of each study and class (its counts and
        ////      its places in the locations, the scales and the cutpoints):
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
}


parameters {
        matrix[2, n_covariates_max] beta_mu;
        row_vector<lower=0.0>[2] beta_SD;
        matrix[n_studies, 2] beta_z; //// Study-specific random effects (off-centered parameterisation)
        ////
        matrix[2, n_covariates_max] raw_scale_mu;
        row_vector<lower=0.0>[2] raw_scale_SD;
        matrix[n_studies, 2] raw_scale_z; //// Study-specific random effects (off-centered parameterisation)
        ////
        real<lower=beta_corr_lb, upper=beta_corr_ub> beta_corr;  //// between-study corr (possibly restricted)
        real<lower=raw_scale_corr_lb, upper=raw_scale_corr_ub> raw_scale_corr;  //// between-study corr (possibly restricted)
        ////
        vector<lower=-5.0, upper=5.0>[2] lambda; //// box-cox params
}


transformed parameters {
        ////
        //// ---- Construct simple 2x2 (bivariate) between-study corr matrices:
        ////
        cholesky_factor_corr[2] beta_L_Omega       = make_bivariate_L_Omega(beta_corr);
        cholesky_factor_corr[2] raw_scale_L_Omega  = make_bivariate_L_Omega(raw_scale_corr);
        cholesky_factor_cov[2]  beta_L_Sigma       = diag_pre_multiply(beta_SD, beta_L_Omega);
        cholesky_factor_cov[2]  raw_scale_L_Sigma  = diag_pre_multiply(raw_scale_SD, raw_scale_L_Omega);
        ////
        //// ---- Compute (global) cutpoints using either log or box-cox transform:
        ////
        array[2] vector[n_thr] C_raw;  // Unconstrained Box-Cox output
        array[2] vector[n_thr] C;      // Scaled to [-2, +2]
        C_raw[1] = ((box_cox == 0) ? nico_log(cts_thr_values) : fn_nico_box_cox(cts_thr_values, lambda[1]));
        C_raw[2] = ((box_cox == 0) ? nico_log(cts_thr_values) : fn_nico_box_cox(cts_thr_values, lambda[2]));
        ////
        // Step 1: Shift to positive if needed
        real min_thr = min(cts_thr_values);
        real shift = min_thr <= 0 ? abs(min_thr) + 1.0 : 0.0;
        vector[n_thr] cts_thr_shifted = cts_thr_values + shift;

        // Step 2: Apply Box-Cox
        C_raw[1] = fn_nico_box_cox(cts_thr_shifted, lambda[1]);
        C_raw[2] = fn_nico_box_cox(cts_thr_shifted, lambda[2]);

        // Step 3: Scale to [-2, +2] range
        for (c in 1:2) {
            real C_min = min(C_raw[c]);
            real C_max = max(C_raw[c]);
            real C_range = C_max - C_min;

            if (C_range > 0) {
                // Linear rescaling: map [C_min, C_max] → [-2, +2]
                C[c] = -2.0 + 4.0 * (C_raw[c] - C_min) / C_range;
            } else {
                C[c] = rep_vector(0.0, n_thr);  // All same value → center at 0
            }
        }
        ////
        //// ---- Study-level random effects (after Cholesky decomposition):
        ////
        matrix[n_studies, 2] beta_random;
        matrix[n_studies, 2] raw_scale_random;
        for (s in 1:n_studies) {
             beta_random[s, ]      =  beta_z[s, ]      * beta_L_Sigma;
             raw_scale_random[s, ] =  raw_scale_z[s, ] * raw_scale_L_Sigma;
        };
        ////
        //// ---- Linear predictors for each disease statu + apply covariates for non-diseased and diseased groups:
        ////
        matrix[n_studies, 2] Xbeta; // NOT local
        matrix[n_studies, 2] Xraw_scale; // NOT local
        matrix[n_studies, 2] Xscale; // NOT local
        ////
        Xbeta[, 1]       = X_nd[1:n_studies, 1:n_covariates_nd]  * to_vector(beta_mu[1, 1:n_covariates_nd]) + beta_random[, 1];
        Xbeta[, 2]       = X_d[1:n_studies,  1:n_covariates_d]   * to_vector(beta_mu[2, 1:n_covariates_d])  + beta_random[, 2];
        ////
        Xraw_scale[, 1]  = X_nd[1:n_studies, 1:n_covariates_nd]  * to_vector(raw_scale_mu[1, 1:n_covariates_nd]) + raw_scale_random[, 1];
        Xraw_scale[, 2]  = X_d[1:n_studies,  1:n_covariates_d]   * to_vector(raw_scale_mu[2, 1:n_covariates_d])  + raw_scale_random[, 2];
        ////
        Xscale = ((softplus == 1) ? 1.4426950408889634 * nico_log1p_exp(Xraw_scale) : nico_exp(Xraw_scale));  // local
}


model {
        ////
        //// ---- Priors for locations:
        ////
        to_vector(beta_mu[1, 1:n_covariates_nd]) ~ nico_normal(
                                                    to_vector(prior_beta_mu_mean[1, 1:n_covariates_nd]),
                                                    to_vector(prior_beta_mu_SD[1, 1:n_covariates_nd]));
        to_vector(beta_mu[2, 1:n_covariates_d])  ~ nico_normal(
                                                    to_vector(prior_beta_mu_mean[2, 1:n_covariates_d]),
                                                    to_vector(prior_beta_mu_SD[2, 1:n_covariates_d]));
        ////
        to_vector(beta_SD) ~ nico_normal(prior_beta_SD_mean, prior_beta_SD_SD);
        beta_L_Omega ~ nico_lkj_corr_cholesky(prior_beta_corr_LKJ);
        ////
        //// ---- Priors for scales:
        ////
        to_vector(raw_scale_mu[1, 1:n_covariates_nd]) ~ nico_normal(
                                                         to_vector(prior_raw_scale_mu_mean[1, 1:n_covariates_nd]),
                                                         to_vector(prior_raw_scale_mu_SD[1, 1:n_covariates_nd]));
        to_vector(raw_scale_mu[2, 1:n_covariates_d])  ~ nico_normal(
                                                         to_vector(prior_raw_scale_mu_mean[2, 1:n_covariates_d]),
                                                         to_vector(prior_raw_scale_mu_SD[2, 1:n_covariates_d]));
        ////
        to_vector(raw_scale_SD) ~ nico_normal(prior_raw_scale_SD_mean, prior_raw_scale_SD_SD);
        raw_scale_L_Omega ~ nico_lkj_corr_cholesky(prior_raw_scale_corr_LKJ);
        ////
        ////
        //// ---- Replace selected baseline location and scale priors:
        ////
        if (use_Se_Sp_prior == 1) {
             for (c in 1:2) {
                  int m = n_chosen[c];
                  if (m > 0) {
                       int n_covariates_c = (c == 1) ? n_covariates_nd : n_covariates_d;
                       vector[n_covariates_c] baseline = (c == 1) ? baseline_case_nd : baseline_case_d;
                       real location = dot_product(baseline, to_vector(beta_mu[c, 1:n_covariates_c]));
                       real raw_scale = dot_product(baseline, to_vector(raw_scale_mu[c, 1:n_covariates_c]));
                       // real scale = (softplus == 1) ? softplus_scaled(raw_scale) : exp(raw_scale);
                       real scale = (softplus == 1) ? 1.4426950408889634 * nico_log1p_exp(raw_scale) :
                                                      nico_exp(raw_scale);
                       vector[m] C_baseline = (C[c][chosen_cutoffs[c, 1:m]] - location) / scale;
                       matrix[m, 2] beta_ab = prior_Se_Sp_beta[c][1:m, ];
                       matrix[m, 2] link_mean_SD = prior_Se_Sp_link_mean_SD[c][1:m, ];
                       if (c == 2) {
                            link_mean_SD[, 1] = -link_mean_SD[, 1];
                            if (Se_Sp_prior_type == 2) {
                                 beta_ab = append_col(beta_ab[, 2], beta_ab[, 1]);
                            }
                       }
                       target += -nico_normal_lpdf(beta_mu[c, 1] | prior_beta_mu_mean[c, 1], prior_beta_mu_SD[c, 1]);
                       target += Se_Sp_prior_at_chosen_cutoffs_lpdf(
                                 C_baseline | Se_Sp_prior_type, beta_ab, link_mean_SD, use_probit_link);
                       target += nico_log(abs(baseline[1])) - nico_log(scale);
                       if (m >= 2) {
                            real log_scale_derivative = (softplus == 1) ?
                                  nico_log_inv_logit(raw_scale) - nico_log(nico_log(2.0)) : raw_scale;
                            real C_gap = C[c][chosen_cutoffs[c, 2]] - C[c][chosen_cutoffs[c, 1]];
                            target += -nico_normal_lpdf(raw_scale_mu[c, 1] |
                                                   prior_raw_scale_mu_mean[c, 1],
                                                   prior_raw_scale_mu_SD[c, 1]);
                            target += nico_log(abs(baseline[1])) + nico_log(C_gap) +
                                      log_scale_derivative - 2.0 * nico_log(scale);
                       }
                  }
             }
        }
        ////
        //// ---- Prior for box-cox parameters:
        ////
        lambda ~ nico_normal(prior_boxcox_lambda_mean, prior_boxcox_lambda_SD);
        ////
        //// ---- Likelihood / Model:
        ////
        target += nico_normal_lpdf(to_vector(beta_z) | 0, 1); // part of between-study model, NOT prior
        target += nico_normal_lpdf(to_vector(raw_scale_z) | 0, 1); // part of between-study model, NOT prior
        ////
        //// target += raw_scale_to_scale_log_det_J_lp(Xraw_scale, softplus); // dont actually need this as not putting a prior on "Xscale"
        ////
        //// ---- Log-likelihood:
        ////
        // {
        //     ////
        //     //// ---- Get the cutpoint index (k) to map "latent_surv[c][s, cut_i]" to correct cutpoint "C[k]":
        //     ////
        //     array[2] matrix[n_studies, n_thr] latent_surv = map_latent_surv_to_fixed_hetero_C(
        //                                                     C, Xbeta, Xscale, n_studies,
        //                                                     n_obs_cutpoints, cutpoint_index);
        //     target += compute_log_lik_binomial_fact_lp(
        //               latent_surv, use_probit_link, n_thr, n_studies,
        //               x_2, n, N_total, n_obs_cutpoints,
        //               K_fold_CV_indicator);
        // }
        {
            ////
            //// ---- every term in one autodiff node: nico_binomial_counts_above_successive_cutpoints_of_a_cumulative_link_lpmf
            ////
            vector[2 * n_thr] cutpoints_of_both_classes = append_row(C[1], C[2]);
            vector[2 * n_studies] locations_of_both_classes = to_vector(Xbeta);
            vector[2 * n_studies] scales_of_both_classes = to_vector(Xscale);
            //// (as map_latent_surv_to_fixed_hetero_C: -(C - location) / scale)
            vector[n_terms] latent_of_each_term =
                  (locations_of_both_classes[index_of_the_location_of_each_term] -
                   cutpoints_of_both_classes[index_of_the_cutpoint_of_each_term]) ./
                  scales_of_both_classes[index_of_the_location_of_each_term];
            ////
            target += nico_binomial_counts_above_successive_cutpoints_of_a_cumulative_link_lpmf(
                            count_above_the_cutpoint_of_each_term |
                            count_above_the_previous_cutpoint_of_each_term, latent_of_each_term,
                            is_the_first_reported_cutpoint_of_each_term, use_probit_link, 1e-10);
        }
}


generated quantities {
          ////
          //// ---- Compute between-study variance-covariance matrix for location parameters:
          ////
          corr_matrix[2] beta_Omega      = multiply_lower_tri_self_transpose(beta_L_Omega);
          corr_matrix[2] raw_scale_Omega = multiply_lower_tri_self_transpose(raw_scale_L_Omega);
          ////
          matrix[2, 2] beta_Sigma      = multiply_lower_tri_self_transpose(beta_L_Sigma);
          matrix[2, 2] raw_scale_Sigma = multiply_lower_tri_self_transpose(raw_scale_L_Sigma);
          ////
          //// ---- Calculate summary accuracy - "baseline" covariate values:
          ////
          real Xbeta_baseline_nd = dot_product(baseline_case_nd, to_vector(beta_mu[1, 1:n_covariates_nd]));
          real Xbeta_baseline_d  = dot_product(baseline_case_d,  to_vector(beta_mu[2, 1:n_covariates_d]));
          ////
          real Xraw_scale_baseline_nd = dot_product(baseline_case_nd, raw_scale_mu[1, 1:n_covariates_nd]);
          real Xraw_scale_baseline_d  = dot_product(baseline_case_d,  raw_scale_mu[2, 1:n_covariates_d]);
          real scale_nd_baseline = ((softplus == 1) ?
                                    1.4426950408889634 * nico_log1p_exp(Xraw_scale_baseline_nd) :
                                    nico_exp(Xraw_scale_baseline_nd));
          real scale_d_baseline  = ((softplus == 1) ?
                                    1.4426950408889634 * nico_log1p_exp(Xraw_scale_baseline_d)  :
                                    nico_exp(Xraw_scale_baseline_d));
          ////
          //// ---- Calculate baseline Se/Sp:
          ////
          vector[n_thr] Fp_baseline = (use_probit_link == 1) ?
                                      nico_Phi(-(C[1] - Xbeta_baseline_nd)/scale_nd_baseline) :
                                      nico_inv_logit(-(C[1] - Xbeta_baseline_nd)/scale_nd_baseline);
          vector[n_thr] Sp_baseline = 1.0 - Fp_baseline;
          vector[n_thr] Se_baseline = (use_probit_link == 1) ?
                                      nico_Phi(-(C[2] - Xbeta_baseline_d)/scale_d_baseline) :
                                      nico_inv_logit(-(C[2] - Xbeta_baseline_d)/scale_d_baseline);
          ////
          //// ---- Calculate predictive accuracy:
          ////
          vector[2] beta_random_pred      = multi_normal_cholesky_rng(rep_vector(0.0, 2), beta_L_Sigma);
          vector[2] raw_scale_random_pred = multi_normal_cholesky_rng(rep_vector(0.0, 2), raw_scale_L_Sigma);
          ////
          real Xbeta_baseline_pred_nd = dot_product(baseline_case_nd, beta_mu[1, 1:n_covariates_nd]) + beta_random_pred[1];
          real Xbeta_baseline_pred_d  = dot_product(baseline_case_d,  beta_mu[2, 1:n_covariates_d])  + beta_random_pred[2];
          ////
          real Xraw_scale_baseline_pred_nd = dot_product(baseline_case_nd, raw_scale_mu[1, 1:n_covariates_nd]) + raw_scale_random_pred[1];
          real Xraw_scale_baseline_pred_d  = dot_product(baseline_case_d,  raw_scale_mu[2, 1:n_covariates_d]) + raw_scale_random_pred[2];
          real scale_nd_baseline_pred = ((softplus == 1) ?
                                        1.4426950408889634 * nico_log1p_exp(Xraw_scale_baseline_pred_nd) :
                                        nico_exp(Xraw_scale_baseline_pred_nd));
          real scale_d_baseline_pred  = ((softplus == 1) ?
                                        1.4426950408889634 * nico_log1p_exp(Xraw_scale_baseline_pred_d)  :
                                        nico_exp(Xraw_scale_baseline_pred_d));
          ////
          vector[n_thr] Fp_baseline_pred = (use_probit_link == 1) ?
                                           nico_Phi(-(C[1] - Xbeta_baseline_pred_nd)/scale_nd_baseline_pred)  :
                                           nico_inv_logit(-(C[1] - Xbeta_baseline_pred_nd)/scale_nd_baseline_pred);
          vector[n_thr] Sp_baseline_pred = 1.0 - Fp_baseline_pred;
          vector[n_thr] Se_baseline_pred = (use_probit_link == 1) ?
                                           nico_Phi(-(C[2] - Xbeta_baseline_pred_d)/scale_d_baseline_pred) :
                                           nico_inv_logit(-(C[2] - Xbeta_baseline_pred_d)/scale_d_baseline_pred);
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
          {
              array[2] matrix[n_studies, n_thr] latent_surv = nico_map_latent_surv_to_fixed_hetero_C(
                                                              C, Xbeta, Xscale, n_studies,
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
              //        deviance_nd[s] += dev_nd[s, k];
              //        deviance_d[s]  += dev_d[s, k];
              //     }
              //     deviance[s] = deviance_nd[s] + deviance_d[s];
              // }
          }
}























