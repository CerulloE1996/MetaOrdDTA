

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
        int<lower=1> n_thr;
        vector[n_thr] cts_thr_values;
        ////
        //// ---- The trace of the reference prior to copy (the prior-only ordinal-bivariate fit): probit(Sp) and probit(Se)
        //// ---- at each threshold, one row per draw:
        ////
        int<lower=1> N;
        matrix[N, n_thr] probit_Sp_draws_to_copy;
        matrix[N, n_thr] probit_Se_draws_to_copy;
        ////
        //// ---- Other:
        ////
        int<lower=0, upper=1> box_cox;
        int<lower=0, upper=1> softplus;
}


transformed data {
        ////
        //// ---- The shifted thresholds (as DTA_MA_Jones.stan):
        ////
        real min_thr = min(cts_thr_values);
        real shift = min_thr <= 0 ? abs(min_thr) + 1.0 : 0.0;
        vector[n_thr] cts_thr_shifted = cts_thr_values + shift;
}


parameters {
        ////
        //// ---- The Jones priors to find (the SDs of raw_scale_mu and lambda stay at their defaults):
        ////
        vector[2] prior_beta_mu_mean;
        vector<lower=0.0>[2] prior_beta_mu_SD;
        vector[2] prior_raw_scale_mu_mean;
        real<lower=-5.0, upper=5.0> prior_boxcox_lambda_mean;
}


model {
        ////
        //// ---- Hyperpriors:
        ////
        prior_beta_mu_mean       ~ normal(0, 5);
        prior_beta_mu_SD         ~ normal(0, 5);
        ////
        prior_raw_scale_mu_mean  ~ normal(0, 1);
        ////
        prior_boxcox_lambda_mean ~ normal(0, 1);
        ////
        //// ---- Cutpoints (as DTA_MA_Jones.stan), with lambda of both classes at its prior mean:
        ////
        vector[n_thr] C_raw = ((box_cox == 0) ? log(cts_thr_shifted)
                                              : fn_Stan_box_cox(cts_thr_shifted, prior_boxcox_lambda_mean));
        vector[n_thr] C;      // Scaled to [-2, +2]
        real C_min = min(C_raw);
        real C_max = max(C_raw);
        real C_range = C_max - C_min;
        if (C_range > 0) {
            // Linear rescaling: map [C_min, C_max] → [-2, +2]
            C = -2.0 + 4.0 * (C_raw - C_min) / C_range;
        } else {
            C = rep_vector(0.0, n_thr);  // All same value → center at 0
        }
        ////
        //// ---- Scales (as DTA_MA_Jones.stan), with raw_scale_mu at its prior mean:
        ////
        real scale_nd_baseline = ((softplus == 1) ? softplus_scaled(prior_raw_scale_mu_mean[1])
                                                  : exp(prior_raw_scale_mu_mean[1]));
        real scale_d_baseline  = ((softplus == 1) ? softplus_scaled(prior_raw_scale_mu_mean[2])
                                                  : exp(prior_raw_scale_mu_mean[2]));
        ////
        //// ---- beta_mu[c] ~ normal(prior_beta_mu_mean[c], prior_beta_mu_SD[c]), so probit(Sp) and probit(Se)
        //// ---- at each threshold are normal with these means and SDs:
        ////
        for (i in 1:N) {
            target += normal_lpdf( to_vector(probit_Sp_draws_to_copy[i, ]) |
                                   (C - prior_beta_mu_mean[1]) / scale_nd_baseline,
                                   prior_beta_mu_SD[1] / scale_nd_baseline);
            target += normal_lpdf( to_vector(probit_Se_draws_to_copy[i, ]) |
                                   -(C - prior_beta_mu_mean[2]) / scale_d_baseline,
                                   prior_beta_mu_SD[2] / scale_d_baseline);
        }

}
