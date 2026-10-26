//// DTA_MA_strat_bivariate_Se_Sp_prior_NicoStan_fused_AVX.stan
//// The model of MetaOrdDTA's DTA_MA_strat_bivariate_Se_Sp_prior.stan with NicoSIMD's AVX functions:
//// the binomial likelihood of Se and Sp at every threshold in one autodiff node.
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
        vector nico_Phi(vector x);
        vector nico_inv_logit(vector x);
}


data {
        int<lower=1> n_studies;
        int<lower=1> n_thr;
        array[n_studies] int n_obs_cutpoints;
        ////
        // Ordinal data format
        array[2, n_studies, n_thr + 1] int x;  // x[1,s,k] = non-diseased, x[2,s,k] = diseased
        array[2, n_studies, n_thr + 1] int cutpoint_index;
        ////
        // Priors of the mean position: [1] = probit Se, [2] = probit Sp, each per threshold
        array[2] vector[n_thr] prior_mu_mean;
        array[2] vector<lower=0>[n_thr] prior_mu_SD;
        ////
        //// ---- Priors directly on Se (c = 1) and Sp (c = 2) at chosen thresholds:
        //// ---- (hyperparameters on the scale P = Se for c = 1 and P = Sp for c = 2)
        ////
        int<lower=0, upper=1> use_Se_Sp_prior;     //// 0 = the location prior at every threshold
        int<lower=0, upper=n_thr> n_chosen_max;
        array[2] int<lower=0, upper=n_chosen_max> n_chosen;
        array[2, n_chosen_max] int<lower=1, upper=n_thr> chosen_cutoffs;
        array[2] matrix<lower=0.0>[n_chosen_max, 2] prior_Se_Sp_beta;
        int<lower=2, upper=4> Se_Sp_prior_type;  //// 2 = Betas, 4 = link-scale Normal
        array[2] matrix[n_chosen_max, 2] prior_Se_Sp_link_mean_SD;  //// Mean and SD of link(Se), link(Sp)
        ////
        int use_probit_link;  // for consistency with your other models
}


transformed data {
    // Arrays to store binary data for each threshold
    array[n_thr] int n_studies_per_thr;
    array[n_thr, n_studies] int TP;
    array[n_thr, n_studies] int FN;
    array[n_thr, n_studies] int TN;
    array[n_thr, n_studies] int FP;
    array[n_thr, n_studies] int study_has_thr;
    ////
    // The terms of the fused likelihood: one per (threshold, study, Se or Sp) that is reported
    int n_terms;
    array[n_thr * n_studies * 2] int term_counts;
    array[n_thr * n_studies * 2] int term_totals;
    array[n_thr * n_studies * 2] int term_k;
    array[n_thr * n_studies * 2] int term_s;
    array[n_thr * n_studies * 2] int term_c;
    vector[n_thr * n_studies * 2] term_is_first;

    // Initialize
    n_terms = 0;
    for (k in 1:n_thr) {
        n_studies_per_thr[k] = 0;
        for (s in 1:n_studies) {
            study_has_thr[k, s] = 0;
            TP[k, s] = 0;
            FN[k, s] = 0;
            TN[k, s] = 0;
            FP[k, s] = 0;
        }
    }

    // Process each study
    for (s in 1:n_studies) {
        for (k in 1:n_thr) {
            int pos_nd = -1;
            int pos_d = -1;

            // cutpoint_index values are 1, 2, ..., n_thr (not 0-based!)
            for (j in 1:n_obs_cutpoints[s]) {
                if (cutpoint_index[1, s, j + 1] == k) {
                    pos_nd = j + 1;
                }
                if (cutpoint_index[2, s, j + 1] == k) {
                    pos_d = j + 1;
                }
            }

            // Only proceed if both diseased and non-diseased have this threshold
            if (pos_nd > 0 && pos_d > 0 && x[1, s, pos_nd] >= 0 && x[2, s, pos_d] >= 0) {
                study_has_thr[k, s] = 1;
                n_studies_per_thr[k] += 1;

                // The data is cumulative: x[c,s,j] is number testing positive at threshold k or higher
                TP[k, s] = x[2, s, pos_d];
                FN[k, s] = x[2, s, 1] - TP[k, s];

                FP[k, s] = x[1, s, pos_nd];
                TN[k, s] = x[1, s, 1] - FP[k, s];

                // Sanity checks
                if (TP[k, s] < 0 || FN[k, s] < 0 || FP[k, s] < 0 || TN[k, s] < 0) {
                    reject("Negative count detected for study ", s, " threshold ", k);
                }
                if (TP[k, s] > x[2, s, 1]) {
                    reject("TP > total diseased for study ", s, " threshold ", k);
                }
                if (FP[k, s] > x[1, s, 1]) {
                    reject("FP > total non-diseased for study ", s, " threshold ", k);
                }
            }
        }
    }

    //// ---- The terms of the fused likelihood, in the order the original model adds them:
    for (k in 1:n_thr) {
        for (s in 1:n_studies) {
            if (study_has_thr[k, s] == 1) {
                //// ---- the diseased counts (Se), only when the study has diseased subjects there:
                if (TP[k, s] + FN[k, s] > 0) {
                    n_terms += 1;
                    term_counts[n_terms] = TP[k, s];
                    term_totals[n_terms] = TP[k, s] + FN[k, s];
                    term_k[n_terms] = k;
                    term_s[n_terms] = s;
                    term_c[n_terms] = 1;
                    term_is_first[n_terms] = 1.0;
                }
                //// ---- the non-diseased counts (Sp):
                if (TN[k, s] + FP[k, s] > 0) {
                    n_terms += 1;
                    term_counts[n_terms] = TN[k, s];
                    term_totals[n_terms] = TN[k, s] + FP[k, s];
                    term_k[n_terms] = k;
                    term_s[n_terms] = s;
                    term_c[n_terms] = 2;
                    term_is_first[n_terms] = 1.0;
                }
            }
        }
    }

    ////
    //// ---- Thresholds that carry the direct Se/Sp prior instead of the location prior:
    ////
    array[2, n_thr] int is_chosen;
    for (c in 1:2) {
        for (k in 1:n_thr) {
            is_chosen[c, k] = 0;
        }
        int k_prev = 0;
        for (j in 1:n_chosen[c]) {
            int k_j = chosen_cutoffs[c, j];
            if (k_j <= k_prev) reject("chosen_cutoffs must rise: c = ", c);
            is_chosen[c, k_j] = 1;
            k_prev = k_j;
        }
    }
    ////
    //// ---- Only the two prior types that act on mu[k, c] apply here: this model has no cutpoint vector
    //// ---- and no segment between cutoffs, so types 1 and 3 have nothing to act on; and the SD of
    //// ---- type 4 must be positive:
    ////
    if (use_Se_Sp_prior == 1 && Se_Sp_prior_type != 2 && Se_Sp_prior_type != 4) {
        reject("Se_Sp_prior_type must be 2 or 4 here, got ", Se_Sp_prior_type);
    }
    if (use_Se_Sp_prior == 1 && Se_Sp_prior_type == 4) {
        for (c in 1:2) {
            for (j in 1:n_chosen[c]) {
                if (prior_Se_Sp_link_mean_SD[c][j, 2] <= 0.0) {
                    reject("prior_Se_Sp_link_mean_SD must be positive: c = ", c, ", j = ", j);
                }
            }
        }
    }
}


parameters {
    //// Separate parameters for each threshold
    array[n_thr] vector[2] mu;  // [threshold][1=sens, 2=spec]
    array[n_thr] cholesky_factor_corr[2] L_Omega;
    array[n_thr] vector<lower=0.0>[2] sigma;
    array[n_thr, n_studies] vector[2] z;  // random effects
}


model {
    ////
    //// ---- Priors for each threshold:
    ////
    for (k in 1:n_thr) {

          //// ---- The location prior of every threshold that does not carry the Se/Sp prior:
          if (use_Se_Sp_prior == 0) {
                mu[k, 1] ~ normal(prior_mu_mean[1][k], prior_mu_SD[1][k]);
                mu[k, 2] ~ normal(prior_mu_mean[2][k], prior_mu_SD[2][k]);
          } else {
                if (is_chosen[1, k] == 0) {
                      mu[k, 1] ~ normal(prior_mu_mean[1][k], prior_mu_SD[1][k]);
                }
                if (is_chosen[2, k] == 0) {
                      mu[k, 2] ~ normal(prior_mu_mean[2][k], prior_mu_SD[2][k]);
                }
          }
          if (n_studies_per_thr[k] > 2) {

                sigma[k, 1:2] ~ normal(0, 1);
                L_Omega[k] ~ lkj_corr_cholesky(2);

          } else if (n_studies_per_thr[k] == 2) {

                sigma[k, 1:2] ~ normal(0, 0.25);
                L_Omega[k] ~ lkj_corr_cholesky(24);

          } else if (n_studies_per_thr[k] == 1) {

                sigma[k, 1:2] ~ normal(0, 0.01);
                L_Omega[k] ~ lkj_corr_cholesky(50);

          } else {
                sigma[k, 1:2] ~ normal(0.0, 1.0);
          }
          ////
          //// ---- Random effects of every study:
          ////
          for (s in 1:n_studies) {
                  z[k, s] ~ std_normal();
          }

    }
    ////
    ////
    //// ---- The direct Se/Sp prior at the chosen thresholds. Type 4 is a normal straight on mu[k, c],
    //// ---- which IS the link scale of Se (c = 1) and Sp (c = 2), so it carries no Jacobian. Type 2 is
    //// ---- a Beta on the scale P = Se (c = 1) or P = Sp (c = 2), with the Jacobian of the link:
    ////
    if (use_Se_Sp_prior == 1) {
          for (c in 1:2) {
                for (j in 1:n_chosen[c]) {
                      int k = chosen_cutoffs[c, j];
                      if (Se_Sp_prior_type == 4) {
                            target += normal_lpdf(mu[k, c] |
                                                  prior_Se_Sp_link_mean_SD[c][j, 1],
                                                  prior_Se_Sp_link_mean_SD[c][j, 2]);
                      } else {
                            real P = (use_probit_link == 1) ? Phi(mu[k, c]) : inv_logit(mu[k, c]);
                            target += beta_lpdf(P | prior_Se_Sp_beta[c, j, 1], prior_Se_Sp_beta[c, j, 2]);
                            if (use_probit_link == 1) {
                                  target += std_normal_lpdf(mu[k, c]);
                            } else {
                                  target += log_inv_logit(mu[k, c]) + log1m_inv_logit(mu[k, c]);
                            }
                      }
                }
          }
    }
    ////
    //// ---- Likelihood: every reported count of Se and Sp in one node. The probit of each term is its element of
    ////      mu[k] + diag_pre_multiply(sigma[k], L_Omega[k]) * z[k, s], since L_Omega[k, 1] = [1, 0];
    ////      is_first is 1 for every term, so its probability is Phi(probit) exactly:
    ////
    {
          vector[n_terms] latent_of_the_terms;
          for (i in 1:n_terms) {
                int k = term_k[i];
                int s = term_s[i];
                latent_of_the_terms[i] = (term_c[i] == 1) ?
                                         mu[k, 1] + sigma[k, 1] * z[k, s, 1] :
                                         mu[k, 2] + sigma[k, 2] * (L_Omega[k, 2, 1] * z[k, s, 1] +
                                                                   L_Omega[k, 2, 2] * z[k, s, 2]);
          }
          target += nico_binomial_counts_above_successive_cutpoints_of_a_cumulative_link_lupmf(
                          term_counts[1:n_terms] | term_totals[1:n_terms], latent_of_the_terms,
                          term_is_first[1:n_terms], 1, 0.0);
    }
}


generated quantities {
        // Summary statistics for each threshold (matching Jones model output)
        vector[n_thr] Se_baseline;
        vector[n_thr] Sp_baseline;
        vector[n_thr] Fp_baseline;

        for (k in 1:n_thr) {
            // Summary accuracy
            if (use_probit_link == 1) {
                Se_baseline[k] = safe_Phi(mu[k, 1]);
                Sp_baseline[k] = safe_Phi(mu[k, 2]);
            } else {
                Se_baseline[k] = inv_logit(mu[k, 1]);
                Sp_baseline[k] = inv_logit(mu[k, 2]);
            }
            Fp_baseline[k] = 1 - Sp_baseline[k];
        }
}






















