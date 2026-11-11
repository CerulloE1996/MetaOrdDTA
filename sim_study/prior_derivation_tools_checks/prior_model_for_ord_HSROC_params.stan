

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
}


data {
        int<lower=1> n_thr;
        ////
        //// ---- The trace of the reference prior to copy (the prior-only ordinal-bivariate fit): probit(Sp) and probit(Se)
        //// ---- at each threshold, one row per draw:
        ////
        int<lower=1> N;
        matrix[N, n_thr] probit_Sp_draws_to_copy;
        matrix[N, n_thr] probit_Se_draws_to_copy;
        ////
        //// ---- The prior on the cutpoints of each draw (as DTA_MA_ord_HSROC_RANDthr_kappa.stan, alpha fixed):
        ////
        vector<lower=0.0>[n_thr + 1] prior_dirichlet_alpha;
        ////
        //// ---- Other:
        ////
        int<lower=0, upper=1> softplus;
        int<lower=0, upper=1> use_probit_link;
}


transformed data {
        int n_cat = n_thr + 1;
        real mult_nd = -1.0;
        real mult_d  = +1.0; //// hence: mult_d > mult_nd
        ////
        //// ---- The cutpoints at the mean of their induced-Dirichlet prior (ordinal probabilities
        //// ---- prior_dirichlet_alpha / sum(prior_dirichlet_alpha)):
        ////
        vector[n_thr] cumul_prob_at_the_prior_mean = cumulative_sum(prior_dirichlet_alpha[1:n_thr]) /
                                                     sum(prior_dirichlet_alpha);
        vector[n_thr] C_at_the_prior_mean = ((use_probit_link == 1) ? inv_Phi(cumul_prob_at_the_prior_mean)
                                                                    : logit(cumul_prob_at_the_prior_mean));
}


parameters {
        ////
        //// ---- The ordinal-HSROC priors to find (the SD of raw_scale_mu stays at its default):
        ////
        real prior_beta_mu_mean;
        real<lower=0.0> prior_beta_mu_SD;
        real prior_raw_scale_mu_mean;
}


model {
        ////
        //// ---- Hyperpriors:
        ////
        prior_beta_mu_mean ~ normal(0, 5);
        prior_beta_mu_SD ~ normal(0, 5);
        prior_raw_scale_mu_mean ~ normal(0, 5);
        ////
        //// ---- Cutpoints at their prior mean, anchored at the prior mean of the mid-location (as
        //// ---- DTA_MA_ord_HSROC_FIXEDthr_Se_Sp_prior.stan):
        ////
        real anchor_ct = 0.5 * (mult_nd + mult_d) * prior_beta_mu_mean; // prior mean of the mid-location
        vector[n_thr] C = C_at_the_prior_mean + anchor_ct;
        ////
        //// ---- Scales (as DTA_MA_ord_HSROC_FIXEDthr_Se_Sp_prior.stan), with raw_scale_mu at its prior mean:
        ////
        real scale_nd_baseline = ((softplus == 1) ? softplus_scaled(-mult_nd * prior_raw_scale_mu_mean)
                                                  : exp(-mult_nd * prior_raw_scale_mu_mean));
        real scale_d_baseline  = ((softplus == 1) ? softplus_scaled(-mult_d  * prior_raw_scale_mu_mean)
                                                  : exp(-mult_d  * prior_raw_scale_mu_mean));
        ////
        //// ---- beta_mu ~ normal(prior_beta_mu_mean, prior_beta_mu_SD), so probit(Sp) and probit(Se) at each
        //// ---- threshold are normal with these means and SDs:
        ////
        for (i in 1:N) {
            target += normal_lpdf( to_vector(probit_Sp_draws_to_copy[i, ]) |
                                   (C + mult_nd * prior_beta_mu_mean) * scale_nd_baseline,
                                   prior_beta_mu_SD * scale_nd_baseline);
            target += normal_lpdf( to_vector(probit_Se_draws_to_copy[i, ]) |
                                   -(C + mult_d * prior_beta_mu_mean) * scale_d_baseline,
                                   prior_beta_mu_SD * scale_d_baseline);
        }

}
