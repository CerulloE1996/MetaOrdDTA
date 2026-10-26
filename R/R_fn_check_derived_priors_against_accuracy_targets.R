##
## ---- Compare user-supplied prior inputs with the requested observable accuracy beliefs:
##
## The existing native likelihood-off checker supplies the realised prior. Agreement of its
## summaries with the targets is reported separately from sampling and model qualification.
##


fn_check_derived_priors_against_accuracy_targets <- function(model_prep_obj,
                                                             priors,
                                                             targets,
                                                             check_arguments,
                                                             probability_tolerance) {

    targets <- fn_accuracy_prior_targets(targets = targets)

    if (any(targets$estimand != "baseline")) {
        stop("The native checker reports baseline accuracy; predictive accuracy needs its own native check.")
    }

    if (!is.list(priors) || length(priors) == 0 || is.null(names(priors))
        || anyNA(names(priors)) || any(!nzchar(names(priors))) || anyDuplicated(names(priors))) {
        stop("Supply the complete receiving model's prior list, with your selected derived entries applied.")
    }

    required <- setdiff(names(formals(run_prior_predictive_check_of_Se_and_Sp)), c("model_prep_obj", "priors"))

    if (!is.list(check_arguments) || is.null(names(check_arguments)) || anyNA(names(check_arguments))
        || any(!nzchar(names(check_arguments))) || anyDuplicated(names(check_arguments))
        || !setequal(names(check_arguments), required)) {
        stop("Supply every native prior-check argument explicitly, except model_prep_obj and priors.")
    }

    if (!is.numeric(probability_tolerance) || length(probability_tolerance) != 1
        || !is.finite(probability_tolerance) || probability_tolerance < 0 || probability_tolerance > 1) {
        stop("Supply the absolute probability tolerance explicitly.")
    }

    probs <- check_arguments$probs
    expected_lower <- (1 - targets$interval_level) / 2
    expected_upper <- (1 + targets$interval_level) / 2

    ## The machine epsilon allowance concerns decimal representation, not accuracy agreement.
    if (!is.numeric(probs) || length(probs) != 2 || any(!is.finite(probs))
        || any(abs(expected_lower - probs[1]) > .Machine$double.eps)
        || any(abs(expected_upper - probs[2]) > .Machine$double.eps)) {
        stop("Use the same explicitly supplied equal-tail interval probabilities for the targets and check.")
    }

    message("Receiving-model prior inputs supplied by the caller:")
    print(priors)
    message("Explicit native likelihood-off check arguments:")
    print(check_arguments)
    message(paste0("Absolute mean and endpoint tolerance: ", probability_tolerance, "."))
    progress <- utils::txtProgressBar(min = 0, max = 1, initial = 0, style = 3)
    on.exit(close(progress), add = TRUE)

    checked <- do.call(what = run_prior_predictive_check_of_Se_and_Sp,
                       args = c(list(model_prep_obj = model_prep_obj, priors = priors), check_arguments))
    native_summary <- checked$summary
    columns <- c("test", "quantity", "score_at_or_above", "prior_mean", "prior_lower", "prior_upper",
                 "ESS_bulk", "ESS_tail", "ESS_SD")

    if (!is.data.frame(native_summary) || !all(columns %in% names(native_summary))) {
        stop("The native checker did not return its mean, interval and all three ESS diagnostics.")
    }

    comparison <- targets
    comparison$achieved_mean <- NA_real_
    comparison$achieved_lower <- NA_real_
    comparison$achieved_upper <- NA_real_
    comparison$ESS_bulk <- NA_real_
    comparison$ESS_tail <- NA_real_
    comparison$ESS_SD <- NA_real_

    for (i in seq_len(nrow(targets))) {

        match_row <- which(as.character(native_summary$test) == targets$test[i]
                           & native_summary$quantity == targets$quantity[i]
                           & native_summary$score_at_or_above == targets$raw_cutoff[i])

        if (length(match_row) != 1) {
            stop(paste0("The native result must identify exactly one row for ", targets$target_id[i], "."))
        }

        row <- native_summary[match_row, , drop = FALSE]
        comparison$achieved_mean[i] <- row$prior_mean
        comparison$achieved_lower[i] <- row$prior_lower
        comparison$achieved_upper[i] <- row$prior_upper
        comparison$ESS_bulk[i] <- row$ESS_bulk
        comparison$ESS_tail[i] <- row$ESS_tail
        comparison$ESS_SD[i] <- row$ESS_SD

    }

    achieved <- as.matrix(comparison[, c("achieved_mean", "achieved_lower", "achieved_upper"), drop = FALSE])
    if (any(!is.finite(achieved))) stop("The native accuracy summaries contain nonfinite entries.")

    comparison$mean_error <- comparison$achieved_mean - comparison$mean
    comparison$lower_error <- comparison$achieved_lower - comparison$lower
    comparison$upper_error <- comparison$achieved_upper - comparison$upper
    comparison$mean_matches <- abs(comparison$mean_error) <= probability_tolerance
    comparison$lower_matches <- abs(comparison$lower_error) <= probability_tolerance
    comparison$upper_matches <- abs(comparison$upper_error) <= probability_tolerance
    comparison$all_three_match <- comparison$mean_matches & comparison$lower_matches & comparison$upper_matches

    dimensions <- dim(checked$draws_array)
    divergent <- checked$n_divergent_transitions

    if (length(dimensions) != 3 || any(dimensions < 1)
        || !is.numeric(divergent) || length(divergent) != 1 || !is.finite(divergent)
        || divergent < 0 || divergent != as.integer(divergent)) {
        stop("The actual native draw dimensions and divergent-transition count must be available.")
    }

    retained_transitions <- dimensions[1] * dimensions[2]
    if (divergent > retained_transitions) stop("Divergent transitions exceed the retained sampling transitions.")
    divergent_percent <- 100 * divergent / retained_transitions

    utils::setTxtProgressBar(pb = progress, value = 1)
    print(comparison)
    message(paste0("Native divergent transitions: ", divergent, " / ", retained_transitions,
                   " (", divergent_percent, "%)."))
    message("Summary agreement alone does not qualify gradients, posterior convergence or integration accuracy.")

    list(comparison = comparison,
         matches_requested_summaries = all(comparison$all_three_match),
         probability_tolerance = probability_tolerance,
         divergent_percent = divergent_percent,
         native_check = checked,
         fully_qualified = NA)

}






















