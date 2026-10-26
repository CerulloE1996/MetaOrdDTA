##
## ---- Derive prior inputs by sampling the target-distribution likelihood:
##
## Source this file to load the functions and editable HADS-A belief table. Sourcing does not fit a model.
## Call fn_derive_HADS_prior_inputs with your actual sampled native model, selected target constructor,
## fixed/free input lists, compiled derivation model and explicit sampling/check controls.
##


prior_function_folder <- "/home/enzocerullo/Documents/Work/PhD_work/R_packages/MetaOrdDTA/R"
prior_function_files <- c("R_fns_derive_priors_from_accuracy_means_and_intervals.R",
                           "R_fn_accuracy_prior_bank_from_MetaOrdDTA.R",
                           "R_fn_accuracy_prior_bank_from_a_MetaOrdDTA_model.R",
                           "R_fn_data_for_deriving_the_Jones_marginal_accuracy_prior.R",
                           "R_fn_data_for_deriving_ord_HSROC_marginal_accuracy_prior.R",
                           "R_fns_target_distributions_from_accuracy_beliefs.R",
                           "R_fn_marginal_prior_request_from_accuracy_targets.R",
                           "R_fns_derive_prior_inputs_from_target_draws.R",
                           "R_fn_MetaOrdDTA_inputs_from_sampled_baseline_priors.R",
                           "R_fn_run_prior_predictive_check_of_Se_and_Sp.R",
                           "R_fn_check_derived_priors_against_accuracy_targets.R")

for (prior_function_file in prior_function_files) {
    source(file = file.path(prior_function_folder, prior_function_file), local = FALSE)
}


##
## ---- Edit your probability means and interval endpoints here:
##

targets <- data.frame(target_id = c("Se_at_8", "Sp_at_8"),
                      test = c("HADS-A", "HADS-A"),
                      quantity = c("Se", "Sp"),
                      raw_cutoff = c(8, 8),
                      estimand = c("baseline", "baseline"),
                      mean = c(0.80, 0.80),
                      lower = c(0.65, 0.65),
                      upper = c(0.92, 0.92),
                      interval_level = c(0.95, 0.95))


##
## ---- Run the selected derivation and inspect its actual receiving-model prior:
##
## estimate and prior_values use beta_baseline_mean, beta_baseline_SD, raw_scale_baseline_mean,
## raw_scale_baseline_SD, then boxcox_mean/boxcox_SD for Jones or dirichlet_alpha for ord_HSROC.
## Jones normal entries have two statuses; ord_HSROC shares one normal entry between statuses.
## A FALSE mask retains the exact prior_values entry supplied by you, even if the proposal differs.
## The current flat-measure finite-bank models refuse free importance-weight parameters.
## Native checking is separate from target construction and reports the achieved means and endpoints.
##


fn_derive_HADS_prior_inputs <- function(model_samples_obj,
                                        targets,
                                        family_of_each_target,
                                        constructor_arguments,
                                        constructors,
                                        n_target_draws,
                                        seed_of_target_draws,
                                        objective,
                                        target_grain_size,
                                        estimate,
                                        prior_values,
                                        compiled_model,
                                        initial_values_per_chain,
                                        sampling_arguments,
                                        posterior_probs,
                                        native_check_arguments,
                                        probability_tolerance) {

    raw_cutoff_grid <- seq_len(21)
    verified <- fn_accuracy_prior_bank_from_a_MetaOrdDTA_model(model_samples_obj = model_samples_obj,
                                                               raw_cutoff_grid = raw_cutoff_grid,
                                                               targets = targets)
    model_code <- verified$provenance$model_code
    model <- switch(model_code, DTA_MA_Jones = "Jones", DTA_MA_ord_HSROC_RANDthr_kappa = "R_and_G")

    if (is.null(model)) stop("This example requires the actual Jones or ord_HSROC random-cutpoint native fit.")

    native_data <- verified$provenance$model_data
    current_priors <- model_samples_obj$priors
    recorded_prior_fields <- intersect(names(current_priors), names(native_data))
    named_native_priors <- grep(pattern = "^prior_", x = names(native_data), value = TRUE)

    if (length(recorded_prior_fields) == 0 || !all(named_native_priors %in% names(current_priors))) {
        stop("The sampled object's recorded prior inputs are unavailable or incomplete.")
    }

    for (field in recorded_prior_fields) {

        unchanged <- identical(dim(current_priors[[field]]), dim(native_data[[field]]))
        unchanged <- unchanged && length(current_priors[[field]]) == length(native_data[[field]])
        unchanged <- unchanged && isTRUE(all.equal(current_priors[[field]], native_data[[field]],
                                                    tolerance = 0, check.attributes = FALSE))
        if (!unchanged) {
            stop(paste0("The sampled object's current prior input differs from its recorded sampling data: ",
                        field, ". Supply the original unmodified sampled object."))
        }

    }

    native_fit <- model_samples_obj$internal_obj$outs_stan_sampling$stan_mod_samples
    native_draws <- native_fit$draws(variables = NULL, inc_warmup = FALSE, format = "matrix")
    native_bank <- if (model == "Jones") {
        fn_Jones_native_prior_integration_bank(native_prior_draws = native_draws, native_stan_data = native_data)
    } else {
        fn_ord_HSROC_native_prior_integration_bank(native_prior_draws = native_draws,
                                                   native_stan_data = native_data,
                                                   raw_cutoff_grid = raw_cutoff_grid)
    }

    target_distributions <- fn_target_distributions_from_accuracy_beliefs(
                                targets = targets,
                                family_of_each_target = family_of_each_target,
                                constructor_arguments = constructor_arguments,
                                constructors = constructors,
                                n_target_draws = n_target_draws,
                                seed = seed_of_target_draws)
    request <- fn_marginal_prior_request_from_accuracy_targets(target_distributions = target_distributions,
                                                               native_prior_bank = native_bank,
                                                               model = model,
                                                               objective = objective,
                                                               target_grain_size = target_grain_size,
                                                               estimate = estimate,
                                                               prior_values = prior_values)
    derive <- if (model == "Jones") {
        fn_derive_the_Jones_marginal_accuracy_prior_by_sampling
    } else fn_derive_the_R_and_G_marginal_accuracy_prior_by_sampling

    result <- derive(request = request,
                      initial_values_per_chain = initial_values_per_chain,
                      compiled_model = compiled_model,
                      sampling_arguments = sampling_arguments,
                      probs = posterior_probs)
    inputs <- fn_MetaOrdDTA_inputs_from_sampled_baseline_priors(result = result,
                                                               request = request,
                                                               native_stan_data = native_data,
                                                               model = model)
    receiving_priors <- current_priors
    for (field in names(inputs$prior_inputs)) receiving_priors[[field]] <- inputs$prior_inputs[[field]]

    checked <- fn_check_derived_priors_against_accuracy_targets(model_prep_obj = model_samples_obj,
                                                               priors = receiving_priors,
                                                               targets = targets,
                                                               check_arguments = native_check_arguments,
                                                               probability_tolerance = probability_tolerance)
    print(target_distributions$comparison)
    print(result$posterior_summary)
    print(result$posterior_quantiles)
    print(inputs$changes)
    print(checked$comparison)

    list(target_distributions = target_distributions,
         request = request,
         derivation = result,
         native_inputs = inputs,
         native_check = checked)

}






















