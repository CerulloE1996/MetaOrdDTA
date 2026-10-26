##
## ---- Use elicited target marginals without inventing joint accuracy curves:
##
## The native integration bank remains paired. The target distributions enter separately through
## the selected marginal composite likelihood, so target rows are not sorted across cutoffs.
##


fn_marginal_prior_request_from_accuracy_targets <- function(target_distributions,
                                                            native_prior_bank,
                                                            model,
                                                            objective,
                                                            target_grain_size,
                                                            estimate,
                                                            prior_values) {

    if (!is.character(model) || length(model) != 1 || is.na(model)
        || !(model %in% c("Jones", "R_and_G")) || !identical(objective, "marginal_composite")) {
        stop("Select Jones or R_and_G and the marginal_composite objective explicitly.")
    }

    if (!is.list(target_distributions) || !is.data.frame(target_distributions$target_draws)
        || !is.data.frame(target_distributions$targets)) {
        stop("Prepare the explicit marginal target distributions and retain their complete target table.")
    }

    targets <- fn_accuracy_prior_targets(targets = target_distributions$targets)
    draws <- target_distributions$target_draws
    required_draws <- c("target_id", "draw_row", "test", "quantity", "raw_cutoff", "estimand", "probability")

    if (!all(required_draws %in% names(draws)) || nrow(draws) < 2
        || anyNA(draws[, required_draws, drop = FALSE])
        || !setequal(unique(draws$target_id), targets$target_id)
        || anyDuplicated(draws[, c("target_id", "draw_row"), drop = FALSE])) {
        stop("Every target must have uniquely identified complete marginal draws.")
    }

    if (length(unique(targets$test)) != 1 || any(targets$estimand != "baseline")
        || !setequal(unique(targets$quantity), c("Se", "Sp"))) {
        stop("These models require baseline Se and Sp targets for one test per request.")
    }

    if (!is.numeric(draws$probability) || any(!is.finite(draws$probability))
        || any(draws$probability <= 0 | draws$probability >= 1)) {
        stop("Every target draw must be strictly inside zero and one; no rows are clipped or discarded.")
    }

    for (i in seq_len(nrow(targets))) {

        selected <- draws[draws$target_id == targets$target_id[i], , drop = FALSE]
        if (nrow(selected) < 2 || any(selected$test != targets$test[i])
            || any(selected$quantity != targets$quantity[i]) || any(selected$raw_cutoff != targets$raw_cutoff[i])
            || any(selected$estimand != targets$estimand[i])) {
            stop("Each target identifier must retain its actual test, quantity, cutoff and estimand.")
        }

    }

    if (!is.numeric(target_grain_size) || length(target_grain_size) != 1 || !is.finite(target_grain_size)
        || target_grain_size < 1 || target_grain_size != as.integer(target_grain_size)) {
        stop("Supply target_grain_size explicitly in target accuracy points per minimum task.")
    }

    is_Jones <- model == "Jones"
    common_bank <- c("n_thr", "n_bank", "cutpoint_bank", "raw_scale_bank", "bank_raw_scale_mean",
                     "bank_raw_scale_SD", "softplus", "use_probit_link", "baseline_case", "n_covariates", "scope")
    additional_bank <- if (is_Jones) {
        c("cts_thr_values", "lambda_bank", "lambda_upper_bank", "bank_boxcox_mean", "bank_boxcox_SD")
    } else c("raw_cutoff_grid", "phi_bank", "bank_dirichlet_alpha", "n_studies")

    if (!is.list(native_prior_bank) || !all(c(common_bank, additional_bank) %in% names(native_prior_bank))) {
        stop("Prepare the complete integration bank using the selected model's native bank adapter.")
    }

    raw_grid <- native_prior_bank[[if (is_Jones) "cts_thr_values" else "raw_cutoff_grid"]]
    cutoff_index <- match(draws$raw_cutoff, raw_grid)

    if (anyNA(cutoff_index)) stop("Every target cutoff must occur exactly in the full native raw-score grid.")

    fields <- c("beta_baseline_mean", "beta_baseline_SD", "raw_scale_baseline_mean", "raw_scale_baseline_SD",
                if (is_Jones) c("boxcox_mean", "boxcox_SD") else "dirichlet_alpha")

    if (!is.list(estimate) || !identical(names(estimate), fields)
        || !is.list(prior_values) || !identical(names(prior_values), fields)) {
        stop(paste0("Name every fixed/free mask and prior input in this order: ",
                    paste0(fields, collapse = ", "), "."))
    }

    initial <- stats::setNames(vector("list", length(fields)), paste0("estimated_", fields))
    masks <- list()
    fixed <- list()

    for (field in fields) {

        size <- if (field == "dirichlet_alpha") {
            native_prior_bank$n_thr + 1
        } else if (is_Jones) 2 else 1
        mask <- estimate[[field]]
        value <- prior_values[[field]]
        positive <- endsWith(field, "_SD") || field == "dirichlet_alpha"

        if (!is.logical(mask) || !is.null(dim(mask)) || length(mask) != size || anyNA(mask)
            || !is.numeric(value) || !is.null(dim(value)) || length(value) != size
            || any(!is.finite(value)) || (positive && any(value <= 0))) {
            stop(paste0("Supply every explicit mask and finite prior input for ", field,
                        "; SDs and concentrations must be positive."))
        }

        masks[[paste0("estimate_", field)]] <- as.integer(mask)
        fixed[[paste0("fixed_", field)]] <- unname(value)
        initial[[paste0("estimated_", field)]] <- unname(value[mask])

    }

    if (is_Jones) {

        if (!identical(estimate$boxcox_mean, estimate$boxcox_SD)) {
            stop("This Jones model estimates each Box-Cox mean and SD together or keeps both fixed.")
        }
        for (status in seq_len(2)) {
            quantity <- if (status == 1) "Sp" else "Se"
            cutoff_by_status <- cutoff_index[draws$quantity == quantity]
            if (estimate$boxcox_mean[status] &&
                !any(cutoff_by_status > 1 & cutoff_by_status < native_prior_bank$n_thr)) {
                stop("Endpoint cutoffs do not identify the Box-Cox prior; keep those inputs fixed.")
            }
        }
        masks$estimate_boxcox_mean <- NULL
        masks$estimate_boxcox_SD <- NULL
        masks$estimate_boxcox_prior <- as.integer(estimate$boxcox_mean)

    }

    metadata <- c("baseline_case", "n_covariates", "n_studies", "scope")
    data <- native_prior_bank[setdiff(c(common_bank, additional_bank), metadata)]
    data$n_target <- nrow(draws)
    data$target_status <- ifelse(draws$quantity == "Sp", 1, 2)
    data$target_cutoff_index <- cutoff_index
    data$target_probability <- draws$probability
    data$target_grain_size <- as.integer(target_grain_size)
    data$marginal_composite_likelihood <- 1
    data <- c(data, masks, fixed)
    target_map <- draws[, setdiff(required_draws, "probability"), drop = FALSE]
    target_map$native_cutoff_index <- cutoff_index

    message(paste0(model, ": ", nrow(draws), " marginal target observations at ", length(unique(cutoff_index)),
                   " raw cutoffs; ", sum(lengths(initial)), " free prior inputs."))
    message("The target marginals are separate; their row indices do not specify a joint accuracy curve.")

    list(stan_data = data,
         initial = initial,
         n_free = sum(lengths(initial)),
         target_map = target_map,
         targets = targets,
         target_distribution_comparison = target_distributions$comparison,
         baseline_case = native_prior_bank$baseline_case,
         n_covariates = native_prior_bank$n_covariates,
         n_studies = native_prior_bank$n_studies,
         native_stan_data = native_prior_bank$native_stan_data,
         scope = native_prior_bank$scope,
         objective = objective,
         approximation = "self-normalised importance integration over the paired native prior bank")

}






















