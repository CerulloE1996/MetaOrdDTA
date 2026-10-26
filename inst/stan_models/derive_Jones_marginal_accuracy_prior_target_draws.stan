////
//// ---- Marginal accuracy likelihood under the native baseline Jones prior:
////
//// The integration bank contains paired cutpoint, lower and upper lambda, and raw baseline-scale draws from
//// the native prior-only model. Its proposal inputs are supplied explicitly. The normal baseline location
//// is integrated analytically; the remaining integral uses self-normalised importance sampling.
//// Summing these marginal log densities is a composite likelihood, not a joint curve likelihood.
//// Normalising constants of the bounded lambda laws cancel in self-normalised weights.
//// Endpoint probabilities are not included; the probit branch covers the native Phi interior only.
//// The returned normal inputs describe baseline projections, not arbitrary regression coefficients.
//// HMC includes the parameter-transform Jacobian; the native integration bank remains fixed.
//// The free baseline location inputs use their flat constrained-space measure; the free raw-scale and Box-Cox
//// inputs carry the normal and half-normal hyperpriors supplied as data.
//// The sum over the targets of each status is one node with its gradient by hand, in the user header
//// fn_Jones_log_density_of_the_accuracy_targets_of_one_status.hpp (compiled with this model).

functions {

    vector fn_Jones_log_accuracy_components(vector cutpoint_bank,
                                            vector scale_bank,
                                            vector log_scale_bank,
                                            vector log_importance_weight,
                                            real target_z,
                                            real target_log_link_density,
                                            real accuracy_sign,
                                            real beta_mean,
                                            real beta_SD) {

        vector[num_elements(cutpoint_bank)] beta_value = cutpoint_bank
                                                        + accuracy_sign * target_z * scale_bank;

        return log_importance_weight + log_scale_bank - target_log_link_density - log(beta_SD)
               - 0.5 * log(2.0 * pi()) - 0.5 * square((beta_value - beta_mean) / beta_SD);

    }

    real fn_Jones_log_density_of_the_accuracy_targets_of_one_status(vector log_weight_plus_log_scale,
                                                                    data matrix cutpoint_bank_of_the_status,
                                                                    data vector scale_bank_of_the_status,
                                                                    data array[] int target_cutoff_index,
                                                                    data vector target_z,
                                                                    array[] int targets_of_this_status,
                                                                    data real accuracy_sign,
                                                                    real beta_mean,
                                                                    real beta_SD);

    real fn_Jones_accuracy_partial_sum(array[] int target_index_slice,
                                       int start,
                                       int end,
                                       data array[] int target_status,
                                       data array[] int target_cutoff_index,
                                       data array[] matrix cutpoint_bank,
                                       data matrix scale_bank,
                                       data matrix log_scale_bank,
                                       data vector target_z,
                                       data vector target_log_link_density,
                                       matrix log_importance_weight,
                                       vector log_weight_sum,
                                       vector prior_beta_baseline_mean,
                                       vector prior_beta_baseline_SD) {

        int n_slice = size(target_index_slice);
        real log_density = 0.0;

        for (c in 1:2) {

            int n_of_this_status = 0;
            for (s in 1:n_slice) {
                n_of_this_status += target_status[target_index_slice[s]] == c;
            }

            if (n_of_this_status > 0) {

                array[n_of_this_status] int targets_of_this_status;
                real sum_of_log_link_densities = 0.0;
                int m = 0;

                for (s in 1:n_slice) {
                    int i = target_index_slice[s];
                    if (target_status[i] == c) {
                        m += 1;
                        targets_of_this_status[m] = i;
                        sum_of_log_link_densities += target_log_link_density[i];
                    }
                }

                log_density += fn_Jones_log_density_of_the_accuracy_targets_of_one_status(
                                   log_importance_weight[, c] + log_scale_bank[, c],
                                   cutpoint_bank[c],
                                   scale_bank[, c],
                                   target_cutoff_index,
                                   target_z,
                                   targets_of_this_status,
                                   c == 1 ? -1.0 : 1.0,
                                   prior_beta_baseline_mean[c],
                                   prior_beta_baseline_SD[c])
                               - sum_of_log_link_densities
                               - n_of_this_status * (log(prior_beta_baseline_SD[c]) + 0.5 * log(2.0 * pi())
                                                     + log_weight_sum[c]);

            }

        }

        return log_density;

    }

}

data {

    int<lower=1> n_thr;
    vector[n_thr] cts_thr_values;
    int<lower=2> n_bank;
    array[2] matrix[n_bank, n_thr] cutpoint_bank;
    matrix<lower=-5.0, upper=5.0>[n_bank, 2] lambda_bank;
    matrix<lower=-5.0, upper=5.0>[n_bank, 2] lambda_upper_bank;
    matrix[n_bank, 2] raw_scale_bank;
    vector[2] bank_raw_scale_mean;
    vector<lower=0.0>[2] bank_raw_scale_SD;
    vector[2] bank_boxcox_mean;
    vector<lower=0.0>[2] bank_boxcox_SD;

    int<lower=1> n_target;
    int<lower=1> target_grain_size;
    array[n_target] int<lower=1, upper=2> target_status;
    array[n_target] int<lower=1, upper=n_thr> target_cutoff_index;
    vector<lower=0.0, upper=1.0>[n_target] target_probability;

    int<lower=0, upper=1> softplus;
    int<lower=0, upper=1> use_probit_link;
    int<lower=1, upper=1> marginal_composite_likelihood;
    array[2] int<lower=0, upper=1> estimate_beta_baseline_mean;
    array[2] int<lower=0, upper=1> estimate_beta_baseline_SD;
    array[2] int<lower=0, upper=1> estimate_raw_scale_baseline_mean;
    array[2] int<lower=0, upper=1> estimate_raw_scale_baseline_SD;
    vector[2] fixed_beta_baseline_mean;
    vector<lower=0.0>[2] fixed_beta_baseline_SD;
    vector[2] fixed_raw_scale_baseline_mean;
    vector<lower=0.0>[2] fixed_raw_scale_baseline_SD;
    array[2] int<lower=0, upper=1> estimate_boxcox_prior;
    vector[2] fixed_boxcox_mean;
    vector<lower=0.0>[2] fixed_boxcox_SD;


    //// ---- the hyperpriors of the free raw-scale and Box-Cox prior inputs (a normal on each mean, a half-normal
    ////      on each SD through its lower bound of zero); unused where the input is fixed:
    real hyperprior_raw_scale_baseline_mean_mean;
    real<lower=0.0> hyperprior_raw_scale_baseline_mean_SD;
    real<lower=0.0> hyperprior_raw_scale_baseline_SD_scale;
    real hyperprior_boxcox_mean_mean;
    real<lower=0.0> hyperprior_boxcox_mean_SD;
    real<lower=0.0> hyperprior_boxcox_SD_scale;

}

transformed data {

    matrix[n_bank, 2] scale_bank;
    matrix[n_bank, 2] log_scale_bank;
    vector[n_target] target_z;
    vector[n_target] target_log_link_density;
    array[n_target] int target_indices;
    array[2] int n_targets_of_each_status = rep_array(0, 2);
    array[2] int has_interior_cutoff_target = rep_array(0, 2);

    if (min(bank_boxcox_SD) <= 0.0 || min(fixed_boxcox_SD) <= 0.0
        || min(bank_raw_scale_SD) <= 0.0 || min(fixed_beta_baseline_SD) <= 0.0
        || min(fixed_raw_scale_baseline_SD) <= 0.0) {
        reject("The supplied proposal and fixed standard deviations must be positive.");
    }

    for (k in 1:n_thr) {
        if (is_nan(cts_thr_values[k]) || is_inf(cts_thr_values[k])) {
            reject("The full native threshold grid must be finite.");
        }
        if (k > 1 && cts_thr_values[k] <= cts_thr_values[k - 1]) {
            reject("The full native threshold grid must be strictly increasing.");
        }
    }

    for (c in 1:2) {
        for (j in 1:n_bank) {
            if (is_nan(lambda_bank[j, c]) || is_inf(lambda_bank[j, c])
                || is_nan(lambda_upper_bank[j, c]) || is_inf(lambda_upper_bank[j, c])
                || is_nan(raw_scale_bank[j, c]) || is_inf(raw_scale_bank[j, c])) {
                reject("The native integration bank must contain finite lambdas and raw scales.");
            }

            for (k in 1:n_thr) {
                if (is_nan(cutpoint_bank[c][j, k]) || is_inf(cutpoint_bank[c][j, k])) {
                    reject("The paired native cutpoint bank must be finite.");
                }
                if (k > 1 && cutpoint_bank[c][j, k] < cutpoint_bank[c][j, k - 1]) {
                    reject("The paired native cutpoint bank must retain its original cutoff order.");
                }
            }

            scale_bank[j, c] = softplus == 1
                               ? 1.4426950408889634 * log1p_exp(raw_scale_bank[j, c])
                               : exp(raw_scale_bank[j, c]);

            if (scale_bank[j, c] <= 0.0 || is_nan(scale_bank[j, c]) || is_inf(scale_bank[j, c])) {
                reject("Every bank draw must produce a finite positive native scale.");
            }

            log_scale_bank[j, c] = log(scale_bank[j, c]);
        }
    }

    for (i in 1:n_target) {
        int c = target_status[i];
        real native_link_argument;

        target_indices[i] = i;

        if (target_probability[i] <= 0.0 || target_probability[i] >= 1.0
            || is_nan(target_probability[i])) {
            reject("Target probabilities must be strictly between zero and one; no clipping is used.");
        }

        target_z[i] = use_probit_link == 1 ? inv_Phi(target_probability[i]) : logit(target_probability[i]);
        native_link_argument = c == 1 ? -target_z[i] : target_z[i];

        if (use_probit_link == 1 && (native_link_argument < -37.5 || native_link_argument > 8.25)) {
            reject("A target lies outside the native interior Phi branch; endpoint masses are not fitted.");
        }

        target_log_link_density[i] = use_probit_link == 1
                                     ? std_normal_lpdf(target_z[i])
                                     : logistic_lpdf(target_z[i] | 0.0, 1.0);
        n_targets_of_each_status[c] += 1;

        if (target_cutoff_index[i] > 1 && target_cutoff_index[i] < n_thr) {
            has_interior_cutoff_target[c] = 1;
        }
    }

    if (min(n_targets_of_each_status) < 1) {
        reject("Supply specificity targets for status one and sensitivity targets for status two.");
    }

    for (c in 1:2) {
        if (estimate_boxcox_prior[c] == 1 && has_interior_cutoff_target[c] == 0) {
            reject("An estimated class lambda prior needs an interior cutoff target.");
        }
    }

}

parameters {

    vector[sum(estimate_beta_baseline_mean)] estimated_beta_baseline_mean;
    vector<lower=0.0>[sum(estimate_beta_baseline_SD)] estimated_beta_baseline_SD;
    vector[sum(estimate_raw_scale_baseline_mean)] estimated_raw_scale_baseline_mean;
    vector<lower=0.0>[sum(estimate_raw_scale_baseline_SD)] estimated_raw_scale_baseline_SD;
    vector[sum(estimate_boxcox_prior)] estimated_boxcox_mean;
    vector<lower=0.0>[sum(estimate_boxcox_prior)] estimated_boxcox_SD;

}

transformed parameters {

    vector[2] prior_beta_baseline_mean = fixed_beta_baseline_mean;
    vector<lower=0.0>[2] prior_beta_baseline_SD = fixed_beta_baseline_SD;
    vector[2] prior_raw_scale_baseline_mean = fixed_raw_scale_baseline_mean;
    vector<lower=0.0>[2] prior_raw_scale_baseline_SD = fixed_raw_scale_baseline_SD;
    vector[2] prior_boxcox_lambda_mean = fixed_boxcox_mean;
    vector<lower=0.0>[2] prior_boxcox_lambda_SD = fixed_boxcox_SD;
    matrix[n_bank, 2] log_importance_weight;
    vector[2] log_weight_sum;

    {
        int beta_mean_index = 0;
        int beta_SD_index = 0;
        int raw_scale_mean_index = 0;
        int raw_scale_SD_index = 0;
        int boxcox_index = 0;

        for (c in 1:2) {
            if (estimate_beta_baseline_mean[c] == 1) {
                beta_mean_index += 1;
                prior_beta_baseline_mean[c] = estimated_beta_baseline_mean[beta_mean_index];
            }
            if (estimate_beta_baseline_SD[c] == 1) {
                beta_SD_index += 1;
                prior_beta_baseline_SD[c] = estimated_beta_baseline_SD[beta_SD_index];
            }
            if (estimate_raw_scale_baseline_mean[c] == 1) {
                raw_scale_mean_index += 1;
                prior_raw_scale_baseline_mean[c] = estimated_raw_scale_baseline_mean[raw_scale_mean_index];
            }
            if (estimate_raw_scale_baseline_SD[c] == 1) {
                raw_scale_SD_index += 1;
                prior_raw_scale_baseline_SD[c] = estimated_raw_scale_baseline_SD[raw_scale_SD_index];
            }
            if (estimate_boxcox_prior[c] == 1) {
                boxcox_index += 1;
                prior_boxcox_lambda_mean[c] = estimated_boxcox_mean[boxcox_index];
                prior_boxcox_lambda_SD[c] = estimated_boxcox_SD[boxcox_index];
            }
        }
    }

    for (c in 1:2) {
        log_importance_weight[, c] = log(bank_raw_scale_SD[c]) - log(prior_raw_scale_baseline_SD[c])
            - 0.5 * square((raw_scale_bank[, c] - prior_raw_scale_baseline_mean[c])
                           / prior_raw_scale_baseline_SD[c])
            + 0.5 * square((raw_scale_bank[, c] - bank_raw_scale_mean[c]) / bank_raw_scale_SD[c])
            + log(bank_boxcox_SD[c]) - log(prior_boxcox_lambda_SD[c])
            - 0.5 * square((lambda_bank[, c] - prior_boxcox_lambda_mean[c])
                           / prior_boxcox_lambda_SD[c])
            + 0.5 * square((lambda_bank[, c] - bank_boxcox_mean[c]) / bank_boxcox_SD[c])
            + log(bank_boxcox_SD[c]) - log(prior_boxcox_lambda_SD[c])
            - 0.5 * square((lambda_upper_bank[, c] - prior_boxcox_lambda_mean[c])
                           / prior_boxcox_lambda_SD[c])
            + 0.5 * square((lambda_upper_bank[, c] - bank_boxcox_mean[c])
                           / bank_boxcox_SD[c]);
        log_weight_sum[c] = log_sum_exp(log_importance_weight[, c]);
    }

}

model {

    target += reduce_sum(fn_Jones_accuracy_partial_sum,
                         target_indices,
                         target_grain_size,
                         target_status,
                         target_cutoff_index,
                         cutpoint_bank,
                         scale_bank,
                         log_scale_bank,
                         target_z,
                         target_log_link_density,
                         log_importance_weight,
                         log_weight_sum,
                         prior_beta_baseline_mean,
                         prior_beta_baseline_SD);
    ////
    //// ---- the hyperpriors of the free raw-scale and Box-Cox prior inputs (empty vectors add nothing):
    ////
    target += normal_lpdf(estimated_raw_scale_baseline_mean | hyperprior_raw_scale_baseline_mean_mean,
                                                             hyperprior_raw_scale_baseline_mean_SD);
    target += normal_lpdf(estimated_raw_scale_baseline_SD | 0.0, hyperprior_raw_scale_baseline_SD_scale);
    target += normal_lpdf(estimated_boxcox_mean | hyperprior_boxcox_mean_mean, hyperprior_boxcox_mean_SD);
    target += normal_lpdf(estimated_boxcox_SD | 0.0, hyperprior_boxcox_SD_scale);

}

generated quantities {

    matrix[n_bank, 2] normalised_importance_weight;
    vector[2] importance_weight_ESS;
    vector[2] maximum_normalised_importance_weight;
    vector[n_target] marginal_log_density;
    vector[n_target] conditional_integration_ESS;

    for (c in 1:2) {
        normalised_importance_weight[, c] = softmax(log_importance_weight[, c]);
        importance_weight_ESS[c] = 1.0 / dot_self(normalised_importance_weight[, c]);
        maximum_normalised_importance_weight[c] = max(normalised_importance_weight[, c]);
    }

    for (i in 1:n_target) {
        int c = target_status[i];
        real accuracy_sign = c == 1 ? -1.0 : 1.0;
        vector[n_bank] log_components = fn_Jones_log_accuracy_components(
                                           cutpoint_bank[c][, target_cutoff_index[i]],
                                           scale_bank[, c],
                                           log_scale_bank[, c],
                                           log_importance_weight[, c],
                                           target_z[i],
                                           target_log_link_density[i],
                                           accuracy_sign,
                                           prior_beta_baseline_mean[c],
                                           prior_beta_baseline_SD[c]);
        vector[n_bank] conditional_weight = softmax(log_components);

        marginal_log_density[i] = log_sum_exp(log_components) - log_weight_sum[c];
        conditional_integration_ESS[i] = 1.0 / dot_self(conditional_weight);
    }

}






















