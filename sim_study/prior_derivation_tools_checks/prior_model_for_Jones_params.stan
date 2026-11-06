

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
        //// ---- The Jones priors (the hyperparameters to find): the prior of the location of each class, of the log-scale
        //// ---- of each class and of the Box-Cox lambda (one prior for both classes, as DTA_MA_Jones.stan):
        ////
        vector[2] prior_beta_mu_mean;
        vector<lower=0.0>[2] prior_beta_mu_SD;
        vector[2] prior_raw_scale_mu_mean;
        vector<lower=0.0>[2] prior_raw_scale_mu_SD;
        real prior_boxcox_lambda_mean;
        real<lower=0.0> prior_boxcox_lambda_SD;
        ////
        //// ---- The Jones parameters of each draw of the reference prior (as DTA_MA_Jones.stan):
        ////
        matrix[N, 2] beta_mu;
        matrix[N, 2] raw_scale_mu;
        matrix<lower=-5.0, upper=5.0>[N, 2] lambda; //// box-cox params
        ////
        //// ---- The spread of the reference draws around the Jones curve of each class:
        ////
        vector<lower=0.0>[2] sigma;
}


model {
        ////
        //// ---- Hyperpriors:
        ////
        prior_beta_mu_mean ~ normal(0, 5);
        prior_beta_mu_SD ~ normal(0, 5);
        prior_raw_scale_mu_mean ~ normal(0, 5);
        prior_raw_scale_mu_SD ~ normal(0, 5);
        prior_boxcox_lambda_mean ~ normal(0, 5);
        prior_boxcox_lambda_SD ~ normal(0, 5);
        sigma ~ normal(0, 1);
        ////
        //// ---- The Jones priors on the parameters of each draw (the priors of DTA_MA_Jones.stan):
        ////
        for (c in 1:2) {
            beta_mu[, c] ~ normal(prior_beta_mu_mean[c], prior_beta_mu_SD[c]);
            raw_scale_mu[, c] ~ normal(prior_raw_scale_mu_mean[c], prior_raw_scale_mu_SD[c]);
            lambda[, c] ~ normal(prior_boxcox_lambda_mean, prior_boxcox_lambda_SD);
        }
        ////
        //// ---- Each draw of the reference prior: probit(Sp) = (C[1] - beta_mu[1]) / scale_nd and probit(Se) =
        //// ---- -(C[2] - beta_mu[2]) / scale_d (the baseline Se / Sp of DTA_MA_Jones.stan), with the cutpoints and the
        //// ---- scales computed from that draw's Jones parameters:
        ////
        for (i in 1:N) {
            ////
            //// ---- Compute (global) cutpoints using either log or box-cox transform (as DTA_MA_Jones.stan):
            ////
            array[2] vector[n_thr] C_raw;  // Unconstrained Box-Cox output
            array[2] vector[n_thr] C;      // Scaled to [-2, +2]
            C_raw[1] = ((box_cox == 0) ? log(cts_thr_shifted) : fn_Stan_box_cox(cts_thr_shifted, lambda[i, 1]));
            C_raw[2] = ((box_cox == 0) ? log(cts_thr_shifted) : fn_Stan_box_cox(cts_thr_shifted, lambda[i, 2]));
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
            //// ---- The scales at the baseline (as DTA_MA_Jones.stan):
            ////
            real scale_nd_baseline = ((softplus == 1) ? softplus_scaled(raw_scale_mu[i, 1]) : exp(raw_scale_mu[i, 1]));
            real scale_d_baseline  = ((softplus == 1) ? softplus_scaled(raw_scale_mu[i, 2]) : exp(raw_scale_mu[i, 2]));
            ////
            target += normal_lpdf( to_vector(probit_Sp_draws_to_copy[i, ]) | (C[1] - beta_mu[i, 1]) / scale_nd_baseline,
                                                                             sigma[1]);
            target += normal_lpdf( to_vector(probit_Se_draws_to_copy[i, ]) | -(C[2] - beta_mu[i, 2]) / scale_d_baseline,
                                                                             sigma[2]);
        }

}

