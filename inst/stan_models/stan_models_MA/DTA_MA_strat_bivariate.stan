


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
}




parameters {
        // Separate parameters for each threshold
        array[n_thr] vector[2] mu;  // [threshold][1=sens, 2=spec]
        array[n_thr] cholesky_factor_corr[2] L_Omega;
        array[n_thr] vector<lower=0>[2] sigma;
        array[n_thr, n_studies] vector[2] z;  // random effects
}

transformed parameters {
    // Study-specific logit probabilities
    array[n_thr, n_studies] vector[2] logit_pi;
    
    for (k in 1:n_thr) {
        for (s in 1:n_studies) {
            if (study_has_thr[k, s] == 1) {
                logit_pi[k, s] = mu[k] + diag_pre_multiply(sigma[k], L_Omega[k]) * z[k, s];
            }
        }
    }
}

model {
    // Priors for each threshold
    for (k in 1:n_thr) {
        mu[k, 1] ~ normal(0, 1);
        mu[k, 2] ~ normal(0, 1);
        
        sigma[k, 1] ~ normal(0, 1);
        sigma[k, 2] ~ normal(0, 1);
        
        L_Omega[k] ~ lkj_corr_cholesky(2);
        
        // Random effects for studies with this threshold
        for (s in 1:n_studies) {
            // if (study_has_thr[k, s] == 1) {
                z[k, s] ~ std_normal();
            // }
        }
    }
    
    // Likelihood
    for (k in 1:n_thr) {
        for (s in 1:n_studies) {
            if (study_has_thr[k, s] == 1) {
                int n_diseased = TP[k, s] + FN[k, s];
                int n_nondiseased = TN[k, s] + FP[k, s];
                
                if (n_diseased > 0) {
                    TP[k, s] ~ binomial(n_diseased, Phi(logit_pi[k, s, 1]));
                }
                if (n_nondiseased > 0) {
                    TN[k, s] ~ binomial(n_nondiseased, Phi(logit_pi[k, s, 2]));
                }
            }
        }
    }
}

generated quantities {
        // Summary statistics for each threshold (matching Jones model output)
        vector[n_thr] Se_baseline;
        vector[n_thr] Sp_baseline;
        vector[n_thr] Fp_baseline;
        ////
        // Predictive values
        vector[n_thr] Se_baseline_pred;
        vector[n_thr] Sp_baseline_pred;
        vector[n_thr] Fp_baseline_pred;
        ////
        // Between-study correlations and variances
        array[n_thr] corr_matrix[2] Omega;
        array[n_thr] matrix[2, 2] Sigma;
        ////
        // Log-likelihood for model comparison
        // array[n_thr, n_studies, 2] real log_lik;
        
        for (k in 1:n_thr) {
            // Summary accuracy
            if (use_probit_link == 1) {
                Se_baseline[k] = Phi(mu[k, 1]);
                Sp_baseline[k] = Phi(mu[k, 2]);
            } else {
                Se_baseline[k] = inv_logit(mu[k, 1]);
                Sp_baseline[k] = inv_logit(mu[k, 2]);
            }
            Fp_baseline[k] = 1 - Sp_baseline[k];
            
            // Between-study correlation and covariance
            Omega[k] = multiply_lower_tri_self_transpose(L_Omega[k]);
            Sigma[k] = quad_form_diag(Omega[k], sigma[k]);
            
            // Predictive accuracy
            vector[2] pred = multi_normal_rng(mu[k], Sigma[k]);
            if (use_probit_link == 1) {
                Se_baseline_pred[k] = Phi(pred[1]);
                Sp_baseline_pred[k] = Phi(pred[2]);
            } else {
                Se_baseline_pred[k] = inv_logit(pred[1]);
                Sp_baseline_pred[k] = inv_logit(pred[2]);
            }
            Fp_baseline_pred[k] = 1 - Sp_baseline_pred[k];
            
            // // Diagnostic measures
            // DOR[k] = (Se_baseline[k] * Sp_baseline[k]) / ((1 - Se_baseline[k]) * (1 - Sp_baseline[k]));
            // LRp[k] = Se_baseline[k] / (1 - Sp_baseline[k]);
            // LRn[k] = (1 - Se_baseline[k]) / Sp_baseline[k];
            
            // // Log-likelihood for each study at this threshold
            // for (i in 1:n_studies_per_thr[k]) {
            //     int s = study_index[k, i];
            //     int pos = TP[k, i] + FN[k, i];
            //     int neg = TN[k, i] + FP[k, i];
            //     
            //     log_lik[k, s, 1] = binomial_logit_lpmf(TP[k, i] | pos, logit_pi[k, s, 1]);
            //     log_lik[k, s, 2] = binomial_logit_lpmf(TN[k, i] | neg, logit_pi[k, s, 2]);
            // }
        }
}












