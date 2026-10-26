////
//// ---- Marginal accuracy likelihood under the native ord_bivariate fixed-cutpoint prior:
////
//// The integration bank holds paired cutpoint draws of each class from the native likelihood-off model. The
//// location of each class is integrated analytically; the Dirichlet concentrations of the cutpoints are
//// changed through self-normalised importance weights on the bank's category probabilities. The sum over the
//// targets of each class is one node with its gradient by hand (fn_Jones_log_density_of_the_accuracy_targets_
//// of_one_status, the native scale set to one). The free location SDs carry a flat measure; the free
//// concentrations carry a half-normal with the scale given as data.
////

functions {

    real fn_Jones_log_density_of_the_accuracy_targets_of_one_status(vector log_weight_plus_log_scale,
                                                                    data matrix cutpoint_bank_of_the_status,
                                                                    data vector scale_bank_of_the_status,
                                                                    data array[] int target_cutoff_index,
                                                                    data vector target_z,
                                                                    array[] int targets_of_this_status,
                                                                    data real accuracy_sign,
                                                                    real beta_mean,
                                                                    real beta_SD);

    real fn_ord_biv_fixed_accuracy_partial_sum(array[] int target_index_slice,
                                          int start,
                                          int end,
                                          data array[] int target_status,
                                          data array[] int target_cutoff_index,
                                          data array[] matrix cutpoint_bank,
                                          data vector unit_scale_bank,
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
                                   log_importance_weight[, c],
                                   cutpoint_bank[c],
                                   unit_scale_bank,
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
    int<lower=2> n_bank;
    array[2] matrix[n_bank, n_thr] cutpoint_bank; // from PRIOR-ONLY model Ordinal-bivariate model
    vector[2] anchor;
    array[2] vector<lower=0.0>[n_thr + 1] bank_dirichlet_alpha;

    int<lower=1> n_target;
    int<lower=1> target_grain_size;
    array[n_target] int<lower=1, upper=2> target_status;
    array[n_target] int<lower=1, upper=n_thr> target_cutoff_index;
    vector<lower=0.0, upper=1.0>[n_target] target_probability;

    array[2] int<lower=0, upper=1> estimate_beta_baseline_SD;
    vector<lower=0.0>[2] fixed_beta_baseline_SD;
    array[2] int<lower=0, upper=1> estimate_dirichlet_alpha;
    array[2] int<lower=1, upper=n_thr + 1> n_alpha_groups;
    array[2, n_thr + 1] int<lower=1, upper=n_thr + 1> alpha_group;
    real<lower=0.0> hyperprior_alpha_scale;

}

transformed data {

    int n_cat = n_thr + 1;
    array[2] matrix[n_bank, n_cat] log_phi_bank;
    vector[n_bank] unit_scale_bank = rep_vector(1.0, n_bank);
    vector[n_target] target_z;
    vector[n_target] target_log_link_density;
    array[n_target] int target_indices;
    array[2] int n_targets_of_each_status = rep_array(0, 2);
    int n_free_alpha_groups = 0;

    for (c in 1:2) {
        if (is_nan(anchor[c]) || is_inf(anchor[c])) {
            reject("The anchor of each class must be finite.");
        }
        if (min(bank_dirichlet_alpha[c]) <= 0.0) {
            reject("Every bank Dirichlet concentration must be positive.");
        }
        for (k in 1:n_cat) {
            if (alpha_group[c, k] > n_alpha_groups[c]) {
                reject("Every category must belong to one of the alpha groups of its class.");
            }
        }
        if (estimate_dirichlet_alpha[c] == 1) n_free_alpha_groups += n_alpha_groups[c];
        for (j in 1:n_bank) {
            vector[n_thr] C = cutpoint_bank[c][j]' - anchor[c];
            vector[n_thr] log_F;
            vector[n_thr] log_S;
            for (k in 1:n_thr) {
                if (is_nan(C[k]) || is_inf(C[k])) {
                    reject("The paired native cutpoint bank must be finite.");
                }
                if (k > 1 && C[k] < C[k - 1]) {
                    reject("The paired native cutpoint bank must retain its original cutoff order.");
                }
                log_F[k] = std_normal_lcdf(C[k] |);
                log_S[k] = std_normal_lccdf(C[k] |);
            }
            log_phi_bank[c][j, 1] = log_F[1];
            log_phi_bank[c][j, n_cat] = log_S[n_thr];
            if (n_thr > 1) {
                for (k in 2:n_thr) {
                    log_phi_bank[c][j, k] = (C[k - 1] >= 0)
                                            ? log_diff_exp(log_S[k - 1], log_S[k])
                                            : log_diff_exp(log_F[k], log_F[k - 1]);
                }
            }
            for (k in 1:n_cat) {
                if (is_nan(log_phi_bank[c][j, k]) || is_inf(log_phi_bank[c][j, k])) {
                    reject("A bank category log probability is not finite: class ", c, ", draw ", j,
                           ", category ", k);
                }
            }
        }
    }

    for (i in 1:n_target) {
        target_indices[i] = i;
        if (target_probability[i] <= 0.0 || target_probability[i] >= 1.0 || is_nan(target_probability[i])) {
            reject("Target probabilities must be strictly between zero and one; no clipping is used.");
        }
        target_z[i] = inv_Phi(target_probability[i]);
        target_log_link_density[i] = std_normal_lpdf(target_z[i]);
        n_targets_of_each_status[target_status[i]] += 1;
    }

    if (min(n_targets_of_each_status) < 1) {
        reject("Supply specificity targets for status one and sensitivity targets for status two.");
    }

}

parameters {

    vector<lower=0.0>[sum(estimate_beta_baseline_SD)] estimated_beta_baseline_SD;
    vector<lower=0.0>[n_free_alpha_groups] estimated_alpha_group_concentration;

}

transformed parameters {

    vector<lower=0.0>[2] prior_beta_baseline_SD = fixed_beta_baseline_SD;
    array[2] vector<lower=0.0>[n_cat] prior_dirichlet_alpha = bank_dirichlet_alpha;
    matrix[n_bank, 2] log_importance_weight;
    vector[2] log_weight_sum;

    {
        int SD_index = 0;
        int group_offset = 0;
        for (c in 1:2) {
            if (estimate_beta_baseline_SD[c] == 1) {
                SD_index += 1;
                prior_beta_baseline_SD[c] = estimated_beta_baseline_SD[SD_index];
            }
            if (estimate_dirichlet_alpha[c] == 1) {
                for (k in 1:n_cat) {
                    prior_dirichlet_alpha[c][k] = estimated_alpha_group_concentration[group_offset
                                                                                       + alpha_group[c, k]];
                }
                group_offset += n_alpha_groups[c];
            }
        }
    }

    for (c in 1:2) {
        log_importance_weight[, c] = log_phi_bank[c] * (prior_dirichlet_alpha[c] - bank_dirichlet_alpha[c]);
        log_weight_sum[c] = log_sum_exp(log_importance_weight[, c]);
    }

}

model {

    target += reduce_sum(fn_ord_biv_fixed_accuracy_partial_sum,
                         target_indices,
                         target_grain_size,
                         target_status,
                         target_cutoff_index,
                         cutpoint_bank,
                         unit_scale_bank,
                         target_z,
                         target_log_link_density,
                         log_importance_weight,
                         log_weight_sum,
                         anchor,
                         prior_beta_baseline_SD);
    ////
    //// ---- the half-normal hyperprior of the free concentrations (an empty vector adds nothing):
    ////
    target += normal_lpdf(estimated_alpha_group_concentration | 0.0, hyperprior_alpha_scale);

}

generated quantities {

    vector[2] importance_weight_ESS;
    vector[2] maximum_normalised_importance_weight;

    for (c in 1:2) {
        vector[n_bank] normalised_importance_weight = softmax(log_importance_weight[, c]);
        importance_weight_ESS[c] = 1.0 / dot_self(normalised_importance_weight);
        maximum_normalised_importance_weight[c] = max(normalised_importance_weight);
    }

}
