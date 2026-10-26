##
## ---- Accuracy targets and entrywise fixed or estimated prior inputs:
##
## The integration bank comes from the actual likelihood-off model. Its rows are never resampled,
## sorted independently or replaced by draws from a nominal Dirichlet approximation.
## Gaussian variables and category probabilities must be those in the model's own prior factors.
## Model-specific anchors, constraints and every unrepresented prior factor remain fixed.
##


fn_accuracy_prior_targets <- function(targets) {

    required <- c("target_id", "test", "quantity", "raw_cutoff", "estimand", "mean", "lower", "upper",
                  "interval_level")

    if (!is.data.frame(targets) || nrow(targets) < 1 || !all(required %in% names(targets))) {
        stop(paste0("targets must contain: ", paste(required, collapse = ", "), "."))
    }

    targets <- targets[, required, drop = FALSE]
    rownames(targets) <- NULL

    for (nm in c("target_id", "test", "quantity", "estimand")) {
        if (!is.character(targets[[nm]]) || anyNA(targets[[nm]]) || any(!nzchar(targets[[nm]]))) {
            stop(paste0(nm, " must contain explicit nonempty character values."))
        }
    }

    if (anyDuplicated(targets$target_id) || !all(targets$quantity %in% c("Se", "Sp"))) {
        stop("Each target_id must be unique; quantity must be Se or Sp.")
    }

    for (nm in c("raw_cutoff", "mean", "lower", "upper", "interval_level")) {
        if (!is.numeric(targets[[nm]]) || any(!is.finite(targets[[nm]]))) {
            stop(paste0(nm, " must contain finite numeric values."))
        }
    }

    if (any(targets$mean < 0 | targets$mean > 1 | targets$lower < 0 | targets$upper > 1)
        || any(targets$lower >= targets$upper)
        || any(targets$interval_level <= 0 | targets$interval_level >= 1)) {
        stop("Means must be probabilities; endpoints must increase within [0, 1]; levels must be interior.")
    }

    key <- targets[, c("test", "quantity", "raw_cutoff", "estimand"), drop = FALSE]

    if (anyDuplicated(key)) {
        stop("Do not supply the same test, quantity, raw cutoff and estimand more than once.")
    }

    targets

}


fn_accuracy_prior_integration_bank <- function(accuracy_draws,
                                               normal_draws,
                                               log_category_probability,
                                               reference_priors,
                                               provenance) {

    matrices <- list(accuracy_draws = accuracy_draws,
                     normal_draws = normal_draws,
                     log_category_probability = log_category_probability)

    if (!is.matrix(accuracy_draws) || nrow(accuracy_draws) < 2 || ncol(accuracy_draws) < 1) {
        stop("accuracy_draws must be a named matrix with at least two paired native prior draws.")
    }

    for (nm in names(matrices)) {
        x <- matrices[[nm]]

        if (!is.matrix(x) || !is.numeric(x) || nrow(x) != nrow(accuracy_draws) || any(!is.finite(x))) {
            stop(paste0(nm, " must be a finite numeric matrix with the same paired native rows."))
        }

        if (ncol(x) > 0 && (is.null(colnames(x)) || anyNA(colnames(x))
                           || any(!nzchar(colnames(x))) || anyDuplicated(colnames(x)))) {
            stop(paste0(nm, " must have unique, explicit column names."))
        }

    }

    if (any(accuracy_draws < 0 | accuracy_draws > 1) || any(log_category_probability > 0)) {
        stop("The bank must contain actual accuracy probabilities and log category probabilities.")
    }

    blocks <- c("normal_mean", "normal_SD", "alpha")

    if (!is.list(reference_priors) || !identical(names(reference_priors), blocks)) {
        stop("reference_priors must explicitly contain normal_mean, normal_SD and alpha, in that order.")
    }

    expected_names <- list(normal_mean = colnames(normal_draws),
                           normal_SD = colnames(normal_draws),
                           alpha = colnames(log_category_probability))

    for (nm in blocks) {
        x <- reference_priors[[nm]]
        expected <- expected_names[[nm]]

        if (!is.numeric(x) || !is.null(dim(x)) || any(!is.finite(x))
            || length(x) != length(expected) || (length(x) > 0 && !identical(names(x), expected))) {
            stop(paste0("The reference ", nm, " must name the bank columns exactly, in their supplied order."))
        }

        if (nm != "normal_mean" && any(x <= 0)) {
            stop(paste0("Every reference ", nm, " must be strictly positive."))
        }

    }

    provenance_fields <- c("model_file", "model_data", "draw_source", "fixed_inputs", "target_map", "can_change")

    if (!is.list(provenance) || !all(provenance_fields %in% names(provenance))) {
        stop(paste0("Supply provenance: ", paste(provenance_fields, collapse = ", "), "."))
    }

    if (!is.list(provenance$can_change) || !identical(names(provenance$can_change), blocks)) {
        stop("The native model adapter must explicitly identify which prior entries can change.")
    }

    for (nm in blocks) {
        allowed <- provenance$can_change[[nm]]
        reference <- reference_priors[[nm]]

        if (!is.logical(allowed) || anyNA(allowed) || length(allowed) != length(reference)
            || (length(allowed) > 0 && !identical(names(allowed), names(reference)))) {
            stop(paste0("can_change$", nm, " must explicitly name every corresponding prior entry."))
        }

    }

    list(accuracy_draws = accuracy_draws,
         normal_draws = normal_draws,
         log_category_probability = log_category_probability,
         reference_priors = reference_priors,
         provenance = provenance)

}


fn_select_prior_inputs_for_derivation <- function(prior_values,
                                                  estimate,
                                                  prior_bank) {

    blocks <- c("normal_mean", "normal_SD", "alpha")

    if (!is.list(prior_values) || !is.list(estimate)
        || !identical(names(prior_values), blocks) || !identical(names(estimate), blocks)) {
        stop("prior_values and estimate must each explicitly contain normal_mean, normal_SD and alpha.")
    }

    for (nm in blocks) {

        value <- prior_values[[nm]]
        selected <- estimate[[nm]]
        reference <- prior_bank$reference_priors[[nm]]

        if (!is.numeric(value) || !is.null(dim(value)) || any(!is.finite(value))
            || length(value) != length(reference)
            || (length(value) > 0 && !identical(names(value), names(reference)))) {
            stop(paste0(nm, " values must name every corresponding bank column explicitly, in order."))
        }

        if (nm != "normal_mean" && any(value <= 0)) {
            stop(paste0("Every supplied ", nm, " must be positive, including fixed entries."))
        }

        if (!is.logical(selected) || !is.null(dim(selected)) || anyNA(selected)
            || length(selected) != length(value)
            || (length(selected) > 0 && !identical(names(selected), names(value)))) {
            stop(paste0("estimate$", nm, " must explicitly mark every named entry TRUE or FALSE."))
        }

        allowed <- prior_bank$provenance$can_change[[nm]]

        if (any(selected & !allowed) || any(value[!allowed] != reference[!allowed])) {
            stop(paste0(nm, ": an entry changes an input held fixed by the native bank, such as its anchor."))
        }

    }

    initial <- list(free_normal_mean = unname(prior_values$normal_mean[estimate$normal_mean]),
                    free_normal_SD = unname(prior_values$normal_SD[estimate$normal_SD]),
                    free_alpha = unname(prior_values$alpha[estimate$alpha]))

    list(prior_values = prior_values,
         estimate = estimate,
         initial = initial,
         n_free = sum(lengths(initial)))

}


fn_data_for_deriving_priors_from_Se_Sp <- function(targets,
                                                  prior_bank,
                                                  prior_values,
                                                  estimate,
                                                  discrepancy_weights) {

    targets <- fn_accuracy_prior_targets(targets = targets)

    prior_bank <- fn_accuracy_prior_integration_bank(
          accuracy_draws = prior_bank$accuracy_draws,
          normal_draws = prior_bank$normal_draws,
          log_category_probability = prior_bank$log_category_probability,
          reference_priors = prior_bank$reference_priors,
          provenance = prior_bank$provenance)

    if (!identical(colnames(prior_bank$accuracy_draws), targets$target_id)) {
        stop("The native accuracy bank columns must match target_id exactly, in the supplied target order.")
    }

    mapping_columns <- c("target_id", "test", "quantity", "raw_cutoff", "estimand")
    target_map <- prior_bank$provenance$target_map

    if (!is.data.frame(target_map) || !all(mapping_columns %in% names(target_map))
        || !identical(target_map[, mapping_columns, drop = FALSE], targets[, mapping_columns, drop = FALSE])) {
        stop("The native bank target map must match the requested tests, quantities, raw cutoffs and estimands.")
    }

    selection <- fn_select_prior_inputs_for_derivation(prior_values = prior_values,
                                                       estimate = estimate,
                                                       prior_bank = prior_bank)

    if (!is.matrix(discrepancy_weights) || !is.numeric(discrepancy_weights)
        || !identical(dim(discrepancy_weights), c(nrow(targets), 3L))
        || !identical(rownames(discrepancy_weights), targets$target_id)
        || !identical(colnames(discrepancy_weights), c("mean", "lower_CDF", "upper_CDF"))
        || any(!is.finite(discrepancy_weights)) || any(discrepancy_weights < 0)
        || sum(discrepancy_weights) <= 0) {
        stop("Supply nonnegative weights named by target_id and mean, lower_CDF, upper_CDF; not all zero.")
    }

    stan_data <- list(n_bank = nrow(prior_bank$accuracy_draws),
                      n_accuracy = nrow(targets),
                      accuracy_bank = unname(prior_bank$accuracy_draws),
                      target_mean = targets$mean,
                      target_lower = targets$lower,
                      target_upper = targets$upper,
                      interval_level = targets$interval_level,
                      discrepancy_weight = unname(discrepancy_weights),
                      n_normal = ncol(prior_bank$normal_draws),
                      normal_draws = unname(prior_bank$normal_draws),
                      reference_normal_mean = unname(prior_bank$reference_priors$normal_mean),
                      reference_normal_SD = unname(prior_bank$reference_priors$normal_SD),
                      fixed_normal_mean = unname(prior_values$normal_mean),
                      fixed_normal_SD = unname(prior_values$normal_SD),
                      estimate_normal_mean = as.integer(estimate$normal_mean),
                      estimate_normal_SD = as.integer(estimate$normal_SD),
                      n_free_normal_mean = sum(estimate$normal_mean),
                      n_free_normal_SD = sum(estimate$normal_SD),
                      n_alpha = ncol(prior_bank$log_category_probability),
                      log_category_probability = unname(prior_bank$log_category_probability),
                      reference_alpha = unname(prior_bank$reference_priors$alpha),
                      fixed_alpha = unname(prior_values$alpha),
                      estimate_alpha = as.integer(estimate$alpha),
                      n_free_alpha = sum(estimate$alpha))

    list(stan_data = stan_data,
         targets = targets,
         prior_bank = prior_bank,
         prior_values = prior_values,
         estimate = estimate,
         initial = selection$initial,
         n_free = selection$n_free,
         discrepancy_weights = discrepancy_weights,
         scope = "Finite native-bank mean and endpoint-CDF matching; not a full-distribution match.")

}


fn_evaluate_accuracy_prior_targets <- function(request,
                                               prior_values,
                                               probability_tolerance,
                                               CDF_tolerance) {

    fn_select_prior_inputs_for_derivation(prior_values = prior_values,
                                          estimate = request$estimate,
                                          prior_bank = request$prior_bank)

    for (nm in names(prior_values)) {
        held <- !request$estimate[[nm]]

        if (!identical(prior_values[[nm]][held], request$prior_values[[nm]][held])) {
            stop(paste0("The evaluation changed a fixed ", nm, " entry."))
        }

    }

    tolerances <- list(probability_tolerance = probability_tolerance,
                       CDF_tolerance = CDF_tolerance)

    for (nm in names(tolerances)) {
        x <- tolerances[[nm]]

        if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x < 0) {
            stop(paste0(nm, " must be an explicit finite nonnegative probability difference."))
        }

    }

    bank <- request$prior_bank
    reference <- bank$reference_priors
    log_weight <- numeric(nrow(bank$accuracy_draws))

    for (k in seq_len(ncol(bank$normal_draws))) {

        value <- bank$normal_draws[, k]
        log_weight <- log_weight - 0.5 * (
              ((value - prior_values$normal_mean[k]) / prior_values$normal_SD[k])^2
            - ((value - reference$normal_mean[k]) / reference$normal_SD[k])^2)

    }

    if (length(reference$alpha) > 0) {
        log_weight <- log_weight + as.vector(bank$log_category_probability %*%
                                             (prior_values$alpha - reference$alpha))
    }

    if (any(!is.finite(log_weight))) {
        stop("The proposed prior gives nonfinite native-bank log weights; no weights are clipped.")
    }

    bank_weight <- exp(log_weight - max(log_weight))
    bank_weight <- bank_weight / sum(bank_weight)
    results <- request$targets
    achieved_mean <- as.vector(crossprod(bank_weight, bank$accuracy_draws))
    lower_CDF <- upper_CDF <- achieved_lower <- achieved_upper <- numeric(nrow(results))
    lower_CDF_fixed_by_bank <- upper_CDF_fixed_by_bank <- logical(nrow(results))
    lower_probability <- (1 - results$interval_level) / 2
    upper_probability <- 1 - lower_probability

    for (k in seq_len(nrow(results))) {

        value <- bank$accuracy_draws[, k]
        lower_indicator <- value <= results$lower[k]
        upper_indicator <- value <= results$upper[k]
        lower_CDF_fixed_by_bank[k] <- all(lower_indicator) || !any(lower_indicator)
        upper_CDF_fixed_by_bank[k] <- all(upper_indicator) || !any(upper_indicator)
        lower_CDF[k] <- sum(bank_weight[value <= results$lower[k]])
        upper_CDF[k] <- sum(bank_weight[value <= results$upper[k]])
        order_i <- order(value)
        cumulative <- cumsum(bank_weight[order_i])
        lower_i <- which(cumulative >= lower_probability[k])
        upper_i <- which(cumulative >= upper_probability[k])

        if (length(lower_i) == 0 || length(upper_i) == 0) {
            stop("The weighted empirical CDF did not reach a requested quantile; no endpoint is replaced.")
        }

        achieved_lower[k] <- value[order_i[lower_i[1]]]
        achieved_upper[k] <- value[order_i[upper_i[1]]]

    }

    results$achieved_mean <- achieved_mean
    results$achieved_lower <- achieved_lower
    results$achieved_upper <- achieved_upper
    results$achieved_lower_CDF <- lower_CDF
    results$achieved_upper_CDF <- upper_CDF
    results$mean_error <- achieved_mean - results$mean
    results$lower_error <- achieved_lower - results$lower
    results$upper_error <- achieved_upper - results$upper
    results$lower_CDF_error <- lower_CDF - lower_probability
    results$upper_CDF_error <- upper_CDF - upper_probability
    results$summaries_match <- abs(results$mean_error) <= probability_tolerance &
                              abs(results$lower_error) <= probability_tolerance &
                              abs(results$upper_error) <= probability_tolerance
    results$CDFs_match <- abs(results$lower_CDF_error) <= CDF_tolerance &
                         abs(results$upper_CDF_error) <= CDF_tolerance
    results$bank_minimum <- apply(X = bank$accuracy_draws, MARGIN = 2, FUN = min)
    results$bank_maximum <- apply(X = bank$accuracy_draws, MARGIN = 2, FUN = max)
    results$target_outside_bank_support <- results$mean < results$bank_minimum |
                                          results$mean > results$bank_maximum |
                                          results$lower < results$bank_minimum |
                                          results$upper > results$bank_maximum
    results$lower_CDF_fixed_by_bank <- lower_CDF_fixed_by_bank
    results$upper_CDF_fixed_by_bank <- upper_CDF_fixed_by_bank

    message("Requested and achieved summaries; interval endpoints use the inverse weighted empirical CDF.")
    print(results, row.names = FALSE)
    message(paste0("Importance ESS: ", format(1 / sum(bank_weight^2)),
                   " of ", length(bank_weight), "; largest weight: ", format(max(bank_weight)), "."))
    message("A finite-bank result requires an independent check in the actual likelihood-off model.")

    if (any(results$target_outside_bank_support
            | results$lower_CDF_fixed_by_bank | results$upper_CDF_fixed_by_bank)) {
        message(paste0("The bank does not cover some requested summaries, or an endpoint indicator is constant. ",
                       "Changing weights cannot create native draws outside the observed bank support."))
    }

    n_free <- request$n_free
    n_constraints <- sum(request$discrepancy_weights > 0)

    if (n_free > n_constraints) {
        message(paste0("Not identified by these summaries: ", n_free, " free entries but only ",
                       n_constraints, " weighted summary constraints. No unique prior is established."))
    }

    list(comparison = results,
         prior_values = prior_values,
         importance_ESS = 1 / sum(bank_weight^2),
         largest_weight = max(bank_weight),
         bank_weight = bank_weight,
         probability_tolerance = probability_tolerance,
         CDF_tolerance = CDF_tolerance,
         n_free_parameters = n_free,
         n_summary_constraints = n_constraints,
         unique_solution_established = if (n_free > n_constraints) FALSE else NA,
         native_model_check_passed = NA,
         scope = request$scope)

}


fn_derive_priors_from_Se_Sp <- function(request,
                                       compiled_model,
                                       optimisation_arguments,
                                       probability_tolerance,
                                       CDF_tolerance) {

    if (request$n_free == 0) {

        message("Every prior input is fixed: evaluating the requested summaries without optimisation.")

        return(fn_evaluate_accuracy_prior_targets(request = request,
                                                  prior_values = request$prior_values,
                                                  probability_tolerance = probability_tolerance,
                                                  CDF_tolerance = CDF_tolerance))

    }

    if (is.null(compiled_model) || !is.function(compiled_model$optimize)
        || !is.function(compiled_model$stan_file)) {
        stop("Supply the compiled accuracy-summary derivation model, not a fitted native data model.")
    }

    if (basename(compiled_model$stan_file()) != "derive_priors_accuracy_means_intervals.stan") {
        stop("The compiled model must be derive_priors_accuracy_means_intervals.stan.")
    }

    if (!is.list(optimisation_arguments) || is.null(names(optimisation_arguments))
        || anyNA(names(optimisation_arguments)) || any(!nzchar(names(optimisation_arguments)))
        || anyDuplicated(names(optimisation_arguments))
        || any(c("data", "init") %in% names(optimisation_arguments))) {
        stop("Supply uniquely named optimiser arguments; data and init come from the explicit request.")
    }

    required_arguments <- c("seed", "algorithm", "iter", "jacobian", "refresh", "show_messages")

    if (!all(required_arguments %in% names(optimisation_arguments))
        || !identical(optimisation_arguments$jacobian, FALSE)
        || !identical(optimisation_arguments$show_messages, TRUE)
        || !is.numeric(optimisation_arguments$refresh)
        || length(optimisation_arguments$refresh) != 1
        || !is.finite(optimisation_arguments$refresh) || optimisation_arguments$refresh <= 0) {
        stop("Explicit optimiser controls must include seed, algorithm, iter, jacobian FALSE and live output.")
    }

    message(paste0("Native bank rows: ", request$stan_data$n_bank,
                   "; target accuracies: ", request$stan_data$n_accuracy,
                   "; free prior inputs: ", request$n_free, "."))
    message("Entrywise selection: TRUE is derived; FALSE remains exactly as supplied.")
    print(request$estimate)
    message("Supplied prior values; selected entries also supply the starting values.")
    print(request$prior_values)
    message("User-supplied discrepancy weights and optimiser arguments.")
    print(request$discrepancy_weights)
    print(optimisation_arguments)

    fit <- do.call(what = compiled_model$optimize,
                   args = c(list(data = request$stan_data, init = list(request$initial)),
                            optimisation_arguments))

    if (any(fit$return_codes() != 0)) {
        stop("The optimisation did not finish successfully; no prior values are cleared for use.")
    }

    derived <- request$prior_values
    variables <- c(normal_mean = "free_normal_mean", normal_SD = "free_normal_SD", alpha = "free_alpha")

    for (nm in names(derived)) {
        selected <- request$estimate[[nm]]

        if (any(selected)) {

            expected <- paste0(variables[[nm]], "[", seq_len(sum(selected)), "]")
            tab <- fit$summary(variables = unname(variables[[nm]]))

            if (anyDuplicated(tab$variable) || !all(expected %in% tab$variable)) {
                stop(paste0("The optimisation output does not identify every ", nm, " entry exactly."))
            }

            candidate <- tab$estimate[match(expected, tab$variable)]

            if (any(!is.finite(candidate))) {
                stop(paste0("The output contains nonfinite estimated ", nm, " entries."))
            }

            derived[[nm]][selected] <- candidate

        }
    }

    report <- fn_evaluate_accuracy_prior_targets(request = request,
                                                 prior_values = derived,
                                                 probability_tolerance = probability_tolerance,
                                                 CDF_tolerance = CDF_tolerance)

    report$optimisation_fit <- fit
    report$optimisation_converged <- NA
    report$estimate <- request$estimate
    report$objective <- "User-weighted squared differences of native-bank means and endpoint CDFs."
    message(paste0("The optimiser process exited successfully. Its termination output must also confirm ",
                   "convergence; a successful exit alone does not exclude the iteration limit."))
    report

}






















