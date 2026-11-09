// NicoSIMD_AVX512 version of DTA_MA_strat_bivariate.stan written by NicoStan (source model md5 06e80e705a257f0f41820bdcedec2962; NicoSIMD declarations md5 e3c567fd115d2ccf49ec504688009e6b; writer implementation md5 d387df684f8005ca53ed8808961860dc)
// NicoSIMD functions used: inv_logit, Phi, normal, lkj_corr_cholesky
// Fused by hand from that file: the binomial likelihood of Se and Sp of every threshold in one node.



functions {
    real bmvp_simd_lanes();
    real nico_inv_logit(real x);
    real nico_Phi(real x);
    vector nico_inv_logit(vector x);
    vector nico_Phi(vector x);
    row_vector nico_inv_logit(row_vector x);
    row_vector nico_Phi(row_vector x);
    matrix nico_inv_logit(matrix x);
    matrix nico_Phi(matrix x);
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
    real nico_lkj_corr_cholesky_lpdf(matrix L, data real eta);
    real nico_binomial_counts_above_successive_cutpoints_of_a_cumulative_link_lpmf(
          array[] int n, array[] int n_previous, vector latent, data vector is_first, data int link,
          data real bound);
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
        // // Priors (similar to original bivariate model)
        // real MA_prior_mean_sens_mu;
        // real<lower=0> MA_prior_mean_sens_sd;
        // real MA_prior_mean_spec_mu;
        // real<lower=0> MA_prior_mean_spec_sd;
        // real<lower=0> MA_prior_SD_sens_sd;
        // real<lower=0> MA_prior_SD_spec_sd;
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
    // The terms of the fused likelihood: one per (threshold, study, Se or Sp) that is reported
    int n_terms = 0;
    // The (threshold, study) pairs that are reported: z has 2 values per pair
    int n_pairs = 0;
    array[n_thr, n_studies] int pair_of_thr_and_study;
    array[n_thr * n_studies * 2] int term_counts;
    array[n_thr * n_studies * 2] int term_totals;
    array[n_thr * n_studies * 2] int term_k;
    array[n_thr * n_studies * 2] int term_s;
    array[n_thr * n_studies * 2] int term_c;
    vector[n_thr * n_studies * 2] term_is_first;

    // Initialize
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
    
    for (k in 1:n_thr) {
        n_studies_per_thr[k] = 0;
        for (s in 1:n_studies) {
            if (study_has_thr[k, s] == 1) {
                n_studies_per_thr[k] += 1;
            }
        }
    }
    
    // Process each study
    for (s in 1:n_studies) {
        // For each threshold k (1 to n_thr)
        for (k in 1:n_thr) {
            // Find the position in the data corresponding to threshold k
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
                
                // Debug info for first few
                if (s <= 3 && k <= 2) {
                    print("Study ", s, ", Threshold ", k, ":");
                    print("  Position in data: nd=", pos_nd, ", d=", pos_d);
                    print("  cutpoint_index[1,s,pos_nd]=", cutpoint_index[1, s, pos_nd]);
                    print("  Total diseased: ", x[2, s, 1]);
                    print("  Diseased positive at thr ", k, " or higher: ", x[2, s, pos_d]);
                    print("  TP: ", TP[k, s], ", FN: ", FN[k, s]);
                    print("  Total non-diseased: ", x[1, s, 1]);
                    print("  Non-diseased positive at thr ", k, " or higher: ", x[1, s, pos_nd]);
                    print("  FP: ", FP[k, s], ", TN: ", TN[k, s]);
                }
                
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
    
    print("Studies per threshold: ", n_studies_per_thr);
    // The terms of the fused likelihood, in the order the likelihood adds them:
    for (k in 1:n_thr) {
        for (s in 1:n_studies) {
            pair_of_thr_and_study[k, s] = 0;
            if (study_has_thr[k, s] == 1) {
                n_pairs += 1;
                pair_of_thr_and_study[k, s] = n_pairs;
                if (TP[k, s] + FN[k, s] > 0) {
                    n_terms += 1;
                    term_counts[n_terms] = TP[k, s];
                    term_totals[n_terms] = TP[k, s] + FN[k, s];
                    term_k[n_terms] = k;
                    term_s[n_terms] = s;
                    term_c[n_terms] = 1;
                    term_is_first[n_terms] = 1.0;
                }
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
}




parameters {
    //// Separate parameters for each threshold
    array[n_thr] vector[2] mu;  // [threshold][1=sens, 2=spec]
    array[n_thr] cholesky_factor_corr[2] L_Omega;
    array[n_thr] vector<lower=0.0>[2] sigma;
    vector[2 * n_pairs] z;  // random effects: 2 per (threshold, study) pair that is reported
}

model {
    //// ---- Priors for each threshold:
    for (k in 1:n_thr) {
      
          mu[k, 1:2] ~ nico_normal(0, 1);
          if (n_studies_per_thr[k] > 2) { 
            
                sigma[k, 1:2] ~ nico_normal(0, 0.5); 
                L_Omega[k] ~ nico_lkj_corr_cholesky(2);
              
          } else if (n_studies_per_thr[k] == 2) {
               
                sigma[k, 1:2] ~ nico_normal(0, 0.25); 
                L_Omega[k] ~ nico_lkj_corr_cholesky(12);
                
          } else if (n_studies_per_thr[k] == 1) {
            
                sigma[k, 1:2] ~ nico_normal(0, 0.01); 
                L_Omega[k] ~ nico_lkj_corr_cholesky(25);
                
          } else { 
                sigma[k, 1:2] ~ nico_normal(0.0, 1.0);
          }
            
    }
    // Random effects of the studies that report each threshold
    z ~ nico_normal(0, 1);
    //// ---- Likelihood: every TP ~ binomial(n_diseased, Phi(probit_pi[k, s, 1])) and
    ////      TN ~ binomial(n_nondiseased, Phi(probit_pi[k, s, 2])) term in one node (probit link, is_first = 1).
    ////      probit_pi[k, s] = mu[k] + diag_pre_multiply(sigma[k], L_Omega[k]) * z[k, s], element by element
    ////      (L_Omega[k, 1] = [1, 0]):
    {
          vector[n_terms] latent_of_the_terms;
          for (i in 1:n_terms) {
                int k = term_k[i];
                int p = pair_of_thr_and_study[k, term_s[i]];
                latent_of_the_terms[i] = (term_c[i] == 1) ?
                                         mu[k, 1] + sigma[k, 1] * z[2 * p - 1] :
                                         mu[k, 2] + sigma[k, 2] * (L_Omega[k, 2, 1] * z[2 * p - 1] +
                                                                   L_Omega[k, 2, 2] * z[2 * p]);
          }
          target += nico_binomial_counts_above_successive_cutpoints_of_a_cumulative_link_lupmf(
                          term_counts[1:n_terms] | term_totals[1:n_terms], latent_of_the_terms,
                          term_is_first[1:n_terms], 1, 0.0);
    }
}

generated quantities {
    real simd_lanes = bmvp_simd_lanes();
        // Summary statistics for each threshold (matching Jones model output)
        vector[n_thr] Se_baseline;
        vector[n_thr] Sp_baseline;
        vector[n_thr] Fp_baseline;
        ////
        // // Predictive values
        // vector[n_thr] Se_baseline_pred;
        // vector[n_thr] Sp_baseline_pred;
        // vector[n_thr] Fp_baseline_pred;
        ////
        // // Between-study correlations and variances
        // array[n_thr] corr_matrix[2] Omega;
        // array[n_thr] matrix[2, 2] Sigma;
        ////
        // Log-likelihood for model comparison
        // array[n_thr, n_studies, 2] real log_lik;
        
        for (k in 1:n_thr) {
            // Summary accuracy
            if (use_probit_link == 1) {
                Se_baseline[k] = safe_Phi(mu[k, 1]);
                Sp_baseline[k] = safe_Phi(mu[k, 2]);
            } else {
                Se_baseline[k] = nico_inv_logit(mu[k, 1]);
                Sp_baseline[k] = nico_inv_logit(mu[k, 2]);
            }
            Fp_baseline[k] = 1 - Sp_baseline[k];
        }
            
        // for (k in 1:n_thr) {
        //     if (n_studies_per_thr[k] > 1) {
        //         // Between-study correlation and covariance
        //         Omega[k] = multiply_lower_tri_self_transpose(L_Omega[k]);
        //         Sigma[k] = quad_form_diag(Omega[k], sigma[k]);
        //         
        //         // Predictive accuracy
        //         vector[2] pred = multi_normal_rng(mu[k], Sigma[k]);
        //         if (use_probit_link == 1) {
        //             Se_baseline_pred[k] = safe_Phi(pred[1]);
        //             Sp_baseline_pred[k] = safe_Phi(pred[2]);
        //         } else {
        //             Se_baseline_pred[k] = inv_logit(pred[1]);
        //             Sp_baseline_pred[k] = inv_logit(pred[2]);
        //         }
        //         Fp_baseline_pred[k] = 1 - Sp_baseline_pred[k];
        //     }
        //     
        //     // // Diagnostic measures
        //     // DOR[k] = (Se_baseline[k] * Sp_baseline[k]) / ((1 - Se_baseline[k]) * (1 - Sp_baseline[k]));
        //     // LRp[k] = Se_baseline[k] / (1 - Sp_baseline[k]);
        //     // LRn[k] = (1 - Se_baseline[k]) / Sp_baseline[k];
        //     
        //     // // Log-likelihood for each study at this threshold
        //     // for (i in 1:n_studies_per_thr[k]) {
        //     //     int s = study_index[k, i];
        //     //     int pos = TP[k, i] + FN[k, i];
        //     //     int neg = TN[k, i] + FP[k, i];
        //     //     
        //     //     log_lik[k, s, 1] = binomial_logit_lpmf(TP[k, i] | pos, logit_pi[k, s, 1]);
        //     //     log_lik[k, s, 2] = binomial_logit_lpmf(TN[k, i] | neg, logit_pi[k, s, 2]);
        //     // }
        // }
}










