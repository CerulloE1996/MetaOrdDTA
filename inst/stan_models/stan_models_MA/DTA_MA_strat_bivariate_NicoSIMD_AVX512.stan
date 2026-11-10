// NicoSIMD_AVX512 version of DTA_MA_strat_bivariate.stan written by NicoStan (source model md5 d5fbe211aa8b666a70c27bd0eaffbad1; NicoSIMD declarations md5 e3c567fd115d2ccf49ec504688009e6b; writer implementation md5 aef1406a74d6694a1d260c4cd31c7a35)
// NicoSIMD functions used: exp, log, log1m, log1p_exp, log1m_exp, inv_logit, log_inv_logit, Phi, inv_Phi, tanh, logit, log1m_inv_logit, normal, lgamma, beta, dirichlet, lkj_corr_cholesky, log_diff_exp, std_normal_log_qf



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
        ////
        //// ---- Include files to compile the necessary custom (user-defined) Stan functions:
        ////
        #include "Stan_fns_basic_NicoSIMD_AVX.stan"
        #include "Stan_fns_Box_Cox_NicoSIMD_AVX.stan"
        #include "Stan_fns_corr_NicoSIMD_AVX.stan"
        #include "Stan_fns_ordinal_NicoSIMD_AVX.stan"
        #include "Stan_fns_log_lik_NicoSIMD_AVX.stan"
        #include "Stan_fns_Jacobian_NicoSIMD_AVX.stan"
        #include "Stan_fns_model_fit_NicoSIMD_AVX.stan"
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
        ////
        // Prior of mu at each threshold: [1=sens, 2=spec][threshold]
        array[2] vector[n_thr] prior_mu_mean;
        array[2] vector<lower=0>[n_thr] prior_mu_SD;
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
}




parameters {
    //// Separate parameters for each threshold
    array[n_thr] vector[2] mu;  // [threshold][1=sens, 2=spec]
    array[n_thr] cholesky_factor_corr[2] L_Omega;
    array[n_thr] vector<lower=0.0>[2] sigma;
    array[n_thr, n_studies] vector[2] z;  // random effects
}

transformed parameters {
    //// Study-specific logit probabilities
    array[n_thr, n_studies] vector[2] probit_pi;
    
    for (k in 1:n_thr) {
        for (s in 1:n_studies) {
            if (study_has_thr[k, s] == 1) {
                probit_pi[k, s] = mu[k] + diag_pre_multiply(sigma[k], L_Omega[k]) * z[k, s];
            }
        }
    }
}

model {
    //// ---- Priors for each threshold:
    for (k in 1:n_thr) {
      
          mu[k, 1:2] ~ nico_normal([prior_mu_mean[1][k], prior_mu_mean[2][k]]',
                                   [prior_mu_SD[1][k], prior_mu_SD[2][k]]');
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
            
          // Random effects for studies with this threshold
          for (s in 1:n_studies) {
               // if (study_has_thr[k, s] == 1) {
                  z[k, s] ~ nico_normal(0, 1);
               // }
          }
          
    }
    //// ---- Likelihood:
    for (k in 1:n_thr) {
        for (s in 1:n_studies) {
            if (study_has_thr[k, s] == 1) {
              
                int n_diseased = TP[k, s] + FN[k, s];
                int n_nondiseased = TN[k, s] + FP[k, s];
                
                if (n_diseased > 0) {
                    TP[k, s] ~ binomial(n_diseased, nico_Phi(probit_pi[k, s, 1]));
                }
                if (n_nondiseased > 0) {
                    TN[k, s] ~ binomial(n_nondiseased, nico_Phi(probit_pi[k, s, 2]));
                }
                
            }
        }
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








