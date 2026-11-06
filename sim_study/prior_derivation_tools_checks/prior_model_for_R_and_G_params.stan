

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
}


parameters {
        ////
        //// ---- The ordinal-HSROC priors (the hyperparameters to find): the prior of the location and of the log-scale:
        ////
        real prior_beta_mu_mean;
        real<lower=0.0> prior_beta_mu_SD;
        real prior_raw_scale_mu_mean;
        real<lower=0.0> prior_raw_scale_mu_SD;
        ////
        //// ---- The ordinal-HSROC parameters of each draw of the reference prior (as DTA_MA_ord_HSROC_RANDthr_kappa.stan):
        ////
        vector[N] beta_mu;
        vector[N] raw_scale_mu;
        matrix<lower=-15, upper=2.5>[N, n_thr] C_raw; //// restricting solves divs/numerical stab. issues w/ lots of sparse categories!
        ////
        //// ---- The spread of the reference draws around the ordinal-HSROC curve, one for both classes (a spread per
        //// ---- class lets the free cutpoints copy one class exactly and that class's spread go to zero):
        ////
        real<lower=0.0> sigma;
}


transformed parameters {
        ////
        //// ---- Construct the cutpoints of each draw:
        ////
        matrix[N, n_thr] C;
        for (i in 1:N) {
                     C[i, ] =  construct_C(C_raw[i, 1:n_thr], softplus);
        }
}


model {
        ////
        //// ---- Hyperpriors:
        ////
        prior_beta_mu_mean ~ normal(0, 5);
        prior_beta_mu_SD ~ normal(0, 5);
        prior_raw_scale_mu_mean ~ normal(0, 5);
        prior_raw_scale_mu_SD ~ normal(0, 5);
        sigma ~ normal(0, 1);
        ////
        //// ---- The ordinal-HSROC priors on the parameters of each draw (the priors of DTA_MA_ord_HSROC_RANDthr_kappa.stan):
        ////
        beta_mu ~ normal(prior_beta_mu_mean, prior_beta_mu_SD);
        raw_scale_mu ~ normal(prior_raw_scale_mu_mean, prior_raw_scale_mu_SD);
        ////
        //// ---- Induced-dirichlet prior on the cutpoints of each draw (as DTA_MA_ord_HSROC_RANDthr_kappa.stan):
        ////
        for (i in 1:N) {
             real anchor_ct = 0.5 * (mult_nd + mult_d) * prior_beta_mu_mean; // prior mean of the mid-location
             vector[n_thr] C_shifted = to_vector(C[i, ]) - anchor_ct;
             vector[n_cat] Ind_Dir_ord_prob;
             if (use_probit_link == 1) Ind_Dir_ord_prob = cumul_probs_to_ord_probs(Phi(C_shifted));
             else                      Ind_Dir_ord_prob = cumul_probs_to_ord_probs(inv_logit(C_shifted));
             ////
             Ind_Dir_ord_prob ~ induced_dirichlet_given_C(
                                C_shifted, prior_dirichlet_alpha, use_probit_link);
        }
        target += raw_C_to_C_log_det_J_lp(C_raw, softplus);
        ////
        //// ---- Each draw of the reference prior: probit(Sp) = (C + mult_nd * beta_mu) * scale_nd and probit(Se) =
        //// ---- -(C + mult_d * beta_mu) * scale_d (the baseline Se / Sp of DTA_MA_ord_HSROC_RANDthr_kappa.stan), with
        //// ---- the scales computed from that draw's log-scale:
        ////
        for (i in 1:N) {
            real scale_nd_baseline = ((softplus == 1) ?
                                     softplus_scaled(-mult_nd * raw_scale_mu[i]) : exp(-mult_nd * raw_scale_mu[i]));
            real scale_d_baseline  = ((softplus == 1) ?
                                     softplus_scaled(-mult_d  * raw_scale_mu[i]) : exp(-mult_d  * raw_scale_mu[i]));
            ////
            target += normal_lpdf( to_vector(probit_Sp_draws_to_copy[i, ]) |
                                   (to_vector(C[i, ]) + mult_nd * beta_mu[i]) * scale_nd_baseline, sigma);
            target += normal_lpdf( to_vector(probit_Se_draws_to_copy[i, ]) |
                                   -(to_vector(C[i, ]) + mult_d * beta_mu[i]) * scale_d_baseline, sigma);
        }

}

