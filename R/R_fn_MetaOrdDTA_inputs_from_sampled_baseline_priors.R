##
## ---- Map sampled baseline-prior inputs to an intercept-only receiving model:
##
## Fixed entries take their exact user-supplied request values; free entries take sampled medians.
## Other native inputs remain exactly as supplied. Proposal priors need not equal the requested priors.
## A multi-coefficient baseline projection has no unique inverse and is not allocated automatically.
##


fn_MetaOrdDTA_inputs_from_sampled_baseline_priors <- function(result,
                                                             request,
                                                             native_stan_data,
                                                             model) {

    if (!is.character(model) || length(model) != 1 || is.na(model)
        || !(model %in% c("Jones", "R_and_G")) || !is.list(result$prior_values)
        || !is.list(request$stan_data) || !is.list(native_stan_data)) {
        stop("Supply the sampled derivation result, its request, actual native data and explicit model.")
    }

    if (!is.list(request$native_stan_data) || !identical(native_stan_data, request$native_stan_data)) {
        stop("The receiving inputs must be the actual complete native-data snapshot retained by this request.")
    }

    if (!is.list(result$stan_data) || !identical(result$stan_data, request$stan_data)) {
        stop("The sampled result must retain this exact derivation request, including its target and fixed inputs.")
    }

    is_Jones <- model == "Jones"
    n_covariates <- if (is_Jones) {
        c(native_stan_data$n_covariates_nd, native_stan_data$n_covariates_d)
    } else native_stan_data$n_covariates
    baseline <- if (is_Jones) {
        list(native_stan_data$baseline_case_nd, native_stan_data$baseline_case_d)
    } else list(native_stan_data$baseline_case)
    n_status <- if (is_Jones) 2 else 1

    if (!is.numeric(n_covariates) || length(n_covariates) != n_status || anyNA(n_covariates)
        || any(n_covariates != 1) || length(baseline) != n_status
        || any(!vapply(X = baseline,
                       FUN = function(value) {

                           is.numeric(value) && is.null(dim(value)) && length(value) == 1 &&
                               is.finite(value) && value == 1

                       }, FUN.VALUE = logical(1)))) {
        stop("Automatic export requires one intercept per status and actual baseline value one.")
    }

    fields <- c("beta_baseline_mean", "beta_baseline_SD", "raw_scale_baseline_mean", "raw_scale_baseline_SD",
                if (is_Jones) c("boxcox_mean", "boxcox_SD") else "dirichlet_alpha")
    destination <- c(beta_baseline_mean = "prior_beta_mu_mean",
                     beta_baseline_SD = "prior_beta_mu_SD",
                     raw_scale_baseline_mean = "prior_raw_scale_mu_mean",
                     raw_scale_baseline_SD = "prior_raw_scale_mu_SD",
                     boxcox_mean = "prior_boxcox_lambda_mean",
                     boxcox_SD = "prior_boxcox_lambda_SD",
                     dirichlet_alpha = "prior_dirichlet_alpha")
    output <- list()
    changes <- list()

    for (field in fields) {

        is_boxcox <- startsWith(field, "boxcox_")
        size <- if (field == "dirichlet_alpha") request$stan_data$n_thr + 1 else n_status
        mask <- request$stan_data[[if (is_boxcox) "estimate_boxcox_prior" else paste0("estimate_", field)]]
        derived <- result$prior_values[[paste0("prior_", field)]]
        fixed <- request$stan_data[[paste0("fixed_", field)]]

        if (!is.numeric(mask) || length(mask) != size || anyNA(mask) || any(!mask %in% c(0, 1))
            || !is.numeric(derived) || length(derived) != size || any(!is.finite(derived))
            || !is.numeric(fixed) || length(fixed) != size || any(!is.finite(fixed))
            || any(derived[mask == 0] != fixed[mask == 0])
            || ((endsWith(field, "_SD") || field == "dirichlet_alpha") && any(derived <= 0))) {
            stop(paste0("The result must preserve every fixed input and valid support for ", field, "."))
        }

        selected <- seq_len(size)
        native_name <- unname(destination[field])
        original <- native_stan_data[[native_name]]

        if (is_Jones && !is_boxcox) {

            if (!is.matrix(original) || !is.numeric(original) || nrow(original) != 2 || ncol(original) < 1) {
                stop(paste0("The actual native coefficient-prior matrix is unavailable: ", native_name, "."))
            }
            old_values <- original[selected, 1]
            replacement <- original
            replacement[selected, 1] <- derived[selected]
            input_names <- paste0(native_name, "[", selected, ",1]")

        } else {

            if (!is.numeric(original) || !is.null(dim(original)) || length(original) != size) {
                stop(paste0("The actual native prior-input vector is unavailable: ", native_name, "."))
            }
            old_values <- original[selected]
            replacement <- original
            replacement[selected] <- derived[selected]
            input_names <- paste0(native_name, "[", selected, "]")

        }

        output[[native_name]] <- replacement
        changes[[field]] <- data.frame(native_input = input_names,
                                       original = unname(old_values),
                                       derived = unname(derived[selected]),
                                       estimated = mask[selected] == 1,
                                       changed = unname(old_values != derived[selected]))

    }

    table <- if (length(changes) > 0) {
        do.call(what = rbind, args = unname(changes))
    } else data.frame(native_input = character(0), original = numeric(0), derived = numeric(0),
                       estimated = logical(0), changed = logical(0))

    print(table, row.names = FALSE)
    message("Apply only these entries to the complete receiving model's prior list, then run its native check.")

    list(prior_inputs = output,
         changes = table,
         model = model,
         scope = "intercept-only baseline; requested fixed inputs exact and other native input fields unchanged",
         native_model_check_passed = NA)

}






















