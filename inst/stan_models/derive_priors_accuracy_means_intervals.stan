////
//// ---- Derive selected prior inputs from accuracy means and interval endpoints:
////
//// The bank is sampled from the actual likelihood-off model. Every row remains paired.
//// Gaussian and Dirichlet density ratios reweight that bank; unchanged factors and constraints
//// remain those of the original model. All Gaussian variables and category probabilities must
//// be the variables to which the receiving model applies the corresponding prior factors.
//// The raw cutpoint anchor, model structure and every omitted prior factor are held fixed.
//// This is a finite-bank summary discrepancy, not an observable-data likelihood or a posterior.
//// The caller supplies each discrepancy weight and every fixed/free choice explicitly.
//// No full accuracy distribution is inferred from its mean and interval endpoints.


data {

    int<lower=2> n_bank;
    int<lower=1> n_accuracy;
    matrix<lower=0.0, upper=1.0>[n_bank, n_accuracy] accuracy_bank;
    vector<lower=0.0, upper=1.0>[n_accuracy] target_mean;
    vector<lower=0.0, upper=1.0>[n_accuracy] target_lower;
    vector<lower=0.0, upper=1.0>[n_accuracy] target_upper;
    vector<lower=0.0, upper=1.0>[n_accuracy] interval_level;
    matrix<lower=0.0>[n_accuracy, 3] discrepancy_weight;

    int<lower=0> n_normal;
    matrix[n_bank, n_normal] normal_draws;
    vector[n_normal] reference_normal_mean;
    vector<lower=0.0>[n_normal] reference_normal_SD;
    vector[n_normal] fixed_normal_mean;
    vector<lower=0.0>[n_normal] fixed_normal_SD;
    array[n_normal] int<lower=0, upper=1> estimate_normal_mean;
    array[n_normal] int<lower=0, upper=1> estimate_normal_SD;
    int<lower=0, upper=n_normal> n_free_normal_mean;
    int<lower=0, upper=n_normal> n_free_normal_SD;

    int<lower=0> n_alpha;
    matrix[n_bank, n_alpha] log_category_probability;
    vector<lower=0.0>[n_alpha] reference_alpha;
    vector<lower=0.0>[n_alpha] fixed_alpha;
    array[n_alpha] int<lower=0, upper=1> estimate_alpha;
    int<lower=0, upper=n_alpha> n_free_alpha;

}


transformed data {

    matrix[n_bank, n_accuracy] below_lower;
    matrix[n_bank, n_accuracy] below_upper;
    vector[n_accuracy] target_lower_probability = (1.0 - interval_level) / 2.0;
    vector[n_accuracy] target_upper_probability = 1.0 - target_lower_probability;

    if (sum(estimate_normal_mean) != n_free_normal_mean
        || sum(estimate_normal_SD) != n_free_normal_SD || sum(estimate_alpha) != n_free_alpha) {
        reject("The free-parameter counts must equal the supplied entrywise selection masks.");
    }

    if (sum(to_vector(discrepancy_weight)) <= 0.0) {
        reject("Supply a positive discrepancy weight for at least one requested summary.");
    }

    if (n_normal > 0) {
        if (min(reference_normal_SD) <= 0.0 || min(fixed_normal_SD) <= 0.0) {
            reject("Every supplied Gaussian prior standard deviation must be positive.");
        }
    }

    if (n_alpha > 0) {
        if (min(reference_alpha) <= 0.0 || min(fixed_alpha) <= 0.0) {
            reject("Every supplied Dirichlet concentration must be positive.");
        }
    }

    for (k in 1:n_accuracy) {

        if (target_lower[k] >= target_upper[k]
            || interval_level[k] <= 0.0 || interval_level[k] >= 1.0) {
            reject("Interval endpoints must be increasing and their probability level strictly interior.");
        }

        for (j in 1:n_bank) {
            below_lower[j, k] = accuracy_bank[j, k] <= target_lower[k];
            below_upper[j, k] = accuracy_bank[j, k] <= target_upper[k];
        }

    }

}


parameters {

    vector[n_free_normal_mean] free_normal_mean;
    vector<lower=0.0>[n_free_normal_SD] free_normal_SD;
    vector<lower=0.0>[n_free_alpha] free_alpha;

}


transformed parameters {

    vector[n_normal] derived_normal_mean = fixed_normal_mean;
    vector[n_normal] derived_normal_SD = fixed_normal_SD;
    vector[n_alpha] derived_alpha = fixed_alpha;
    vector[n_bank] log_weight = rep_vector(0.0, n_bank);
    vector[n_bank] bank_weight;
    vector[n_accuracy] achieved_mean;
    vector[n_accuracy] achieved_lower_probability;
    vector[n_accuracy] achieved_upper_probability;

    {

        int mean_i = 0;
        int SD_i = 0;
        int alpha_i = 0;

        if (n_normal > 0) {
            for (k in 1:n_normal) {

                if (estimate_normal_mean[k] == 1) {
                    mean_i += 1;
                    derived_normal_mean[k] = free_normal_mean[mean_i];
                }

                if (estimate_normal_SD[k] == 1) {
                    SD_i += 1;
                    derived_normal_SD[k] = free_normal_SD[SD_i];
                }

                log_weight += -0.5 * (
                      square((normal_draws[, k] - derived_normal_mean[k]) / derived_normal_SD[k])
                    - square((normal_draws[, k] - reference_normal_mean[k]) / reference_normal_SD[k]));

            }
        }

        if (n_alpha > 0) {
            for (k in 1:n_alpha) {
                if (estimate_alpha[k] == 1) {
                    alpha_i += 1;
                    derived_alpha[k] = free_alpha[alpha_i];
                }
            }

            log_weight += log_category_probability * (derived_alpha - reference_alpha);
        }

    }

    bank_weight = softmax(log_weight);
    achieved_mean = accuracy_bank' * bank_weight;
    achieved_lower_probability = below_lower' * bank_weight;
    achieved_upper_probability = below_upper' * bank_weight;

}


model {

    target += -dot_product(discrepancy_weight[, 1], square(achieved_mean - target_mean));
    target += -dot_product(discrepancy_weight[, 2],
                          square(achieved_lower_probability - target_lower_probability));
    target += -dot_product(discrepancy_weight[, 3],
                          square(achieved_upper_probability - target_upper_probability));

}


generated quantities {

    real importance_ESS = inv(dot_self(bank_weight));
    real largest_importance_weight = max(bank_weight);
    int n_free_parameters = n_free_normal_mean + n_free_normal_SD + n_free_alpha;

}






















