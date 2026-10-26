##
## ---- Retain paired cutpoints, lambdas and scales from the native prior-only Jones model:
##


fn_Jones_native_prior_integration_bank <- function(native_prior_draws,
                                                   native_stan_data) {

      if (!is.matrix(native_prior_draws) || !is.numeric(native_prior_draws)
            || nrow(native_prior_draws) < 2 || is.null(colnames(native_prior_draws))) {
            stop("native_prior_draws must be a named numeric matrix containing at least two native prior draws.");
      }

      if (!is.list(native_stan_data)) {
            stop("native_stan_data must be the actual data list used by the prior-only Jones model.");
      }

      required_data <- c("n_thr", "cts_thr_values", "n_covariates_nd", "n_covariates_d",
                         "baseline_case_nd", "baseline_case_d", "prior_raw_scale_mu_mean",
                         "prior_raw_scale_mu_SD", "prior_boxcox_lambda_mean", "prior_boxcox_lambda_SD",
                         "softplus", "use_probit_link", "K_fold_CV_indicator")

      if (!all(required_data %in% names(native_stan_data))) {
            stop(paste0("The native data are missing: ",
                        paste0(setdiff(required_data, names(native_stan_data)), collapse = ", "), "."));
      }

      if (length(native_stan_data$K_fold_CV_indicator) < 1
            || anyNA(native_stan_data$K_fold_CV_indicator)
            || any(native_stan_data$K_fold_CV_indicator != 0)) {
            stop("The supplied native model data must have every likelihood indicator set to zero.");
      }

      n_thr <- native_stan_data$n_thr
      raw_grid <- native_stan_data$cts_thr_values

      if (!is.numeric(n_thr) || length(n_thr) != 1 || !is.finite(n_thr)
            || n_thr < 1 || n_thr != as.integer(n_thr)
            || !is.numeric(raw_grid) || !is.null(dim(raw_grid)) || length(raw_grid) != n_thr
            || any(!is.finite(raw_grid)) || any(diff(raw_grid) <= 0)) {
            stop("The full native threshold grid must agree with n_thr and be finite and strictly increasing.");
      }

      native_flags <- list(native_stan_data$softplus, native_stan_data$use_probit_link)

      for (flag in native_flags) {

            if (!(is.numeric(flag) || is.logical(flag)) || length(flag) != 1
                  || !is.finite(flag) || !(flag %in% c(0, 1))) {
                  stop("The native softplus and use_probit_link inputs must each be exactly zero or one.");
            }

      }

      required_draws <- c("lambda[1]", "lambda[2]", "lambda_upper[1]", "lambda_upper[2]",
                         "Xraw_scale_baseline_nd", "Xraw_scale_baseline_d")

      if (!all(required_draws %in% colnames(native_prior_draws))
            || anyDuplicated(colnames(native_prior_draws))) {
            stop(paste0("Supply full native draws with unique columns including ",
                        paste0(required_draws, collapse = ", "), ". No columns are inferred."));
      }

      lambda_bank <- unname(native_prior_draws[, required_draws[1:2], drop = FALSE])
      lambda_upper_bank <- unname(native_prior_draws[, required_draws[3:4], drop = FALSE])
      raw_scale_bank <- unname(native_prior_draws[, required_draws[5:6], drop = FALSE])

      if (any(!is.finite(lambda_bank)) || any(lambda_bank < -5 | lambda_bank > 5)
            || any(!is.finite(lambda_upper_bank)) || any(lambda_upper_bank < -5 | lambda_upper_bank > 5)
            || any(!is.finite(raw_scale_bank))) {
            stop("The paired native lambda and raw-scale draws must be finite and within the native support.");
      }

      cutpoint_bank <- array(data = NA_real_, dim = c(2, nrow(native_prior_draws), n_thr))

      for (status in seq_len(2)) {

            cutpoint_columns <- paste0("C[", status, ",", seq_len(n_thr), "]")

            if (!all(cutpoint_columns %in% colnames(native_prior_draws))) {
                  stop(paste0("The native prior draws are missing their actual C columns for status ", status,
                              ". Cutpoints are not reconstructed from lambda."));
            }

            cutpoints <- native_prior_draws[, cutpoint_columns, drop = FALSE]

            if (any(!is.finite(cutpoints))) {
                  stop("The actual native cutpoint draws must be finite; no rows are removed or replaced.");
            }

            cutpoint_bank[status, , ] <- cutpoints

      }

      ##
      ## ---- Project the independent normal coefficient priors at the native baseline covariates:
      ##

      n_covariates <- c(native_stan_data$n_covariates_nd, native_stan_data$n_covariates_d)
      baseline_case <- list(native_stan_data$baseline_case_nd, native_stan_data$baseline_case_d)
      prior_mean <- native_stan_data$prior_raw_scale_mu_mean
      prior_SD <- native_stan_data$prior_raw_scale_mu_SD

      if (!is.numeric(n_covariates) || length(n_covariates) != 2
            || any(!is.finite(n_covariates)) || any(n_covariates < 1)
            || any(n_covariates != as.integer(n_covariates))
            || !is.matrix(prior_mean) || !is.numeric(prior_mean)
            || !is.matrix(prior_SD) || !is.numeric(prior_SD)
            || !identical(dim(prior_mean), dim(prior_SD)) || nrow(prior_mean) != 2
            || ncol(prior_mean) < max(n_covariates)) {
            stop("The native coefficient-prior matrices and per-status covariate counts do not agree.");
      }

      bank_raw_scale_mean <- numeric(2)
      bank_raw_scale_SD <- numeric(2)

      for (status in seq_len(2)) {

            columns <- seq_len(n_covariates[status])
            baseline <- baseline_case[[status]]
            means <- prior_mean[status, columns]
            SDs <- prior_SD[status, columns]

            if (!is.numeric(baseline) || !is.null(dim(baseline))
                  || length(baseline) != n_covariates[status] || any(!is.finite(baseline))
                  || any(!is.finite(means)) || any(!is.finite(SDs)) || any(SDs <= 0)) {
                  stop("Each native baseline vector and its active normal coefficient-prior inputs must agree.");
            }

            bank_raw_scale_mean[status] <- sum(baseline * means)
            bank_raw_scale_SD[status] <- sqrt(sum((baseline * SDs)^2))

      }

      if (any(!is.finite(bank_raw_scale_mean)) || any(!is.finite(bank_raw_scale_SD))
            || any(bank_raw_scale_SD <= 0)) {
            stop("Each native raw-scale baseline projection must have a finite positive standard deviation.");
      }

      bank_boxcox_mean <- native_stan_data$prior_boxcox_lambda_mean
      bank_boxcox_SD <- native_stan_data$prior_boxcox_lambda_SD

      if (!is.numeric(bank_boxcox_mean) || !is.null(dim(bank_boxcox_mean)) || length(bank_boxcox_mean) != 2
            || any(!is.finite(bank_boxcox_mean)) || !is.numeric(bank_boxcox_SD)
            || !is.null(dim(bank_boxcox_SD)) || length(bank_boxcox_SD) != 2
            || any(!is.finite(bank_boxcox_SD)) || any(bank_boxcox_SD <= 0)) {
            stop("Supply both native class-specific Box-Cox prior means and positive standard deviations.");
      }

      message(paste0("Native Jones integration bank: ", nrow(lambda_bank),
                     " paired rows; full raw threshold grid: ", paste0(raw_grid, collapse = ", "), "."))

      list(n_thr = as.integer(n_thr),
           cts_thr_values = unname(raw_grid),
           n_bank = nrow(lambda_bank),
           cutpoint_bank = cutpoint_bank,
           lambda_bank = lambda_bank,
           lambda_upper_bank = lambda_upper_bank,
           raw_scale_bank = raw_scale_bank,
           bank_raw_scale_mean = bank_raw_scale_mean,
           bank_raw_scale_SD = bank_raw_scale_SD,
           bank_boxcox_mean = bank_boxcox_mean,
           bank_boxcox_SD = bank_boxcox_SD,
           softplus = as.integer(native_stan_data$softplus),
           use_probit_link = as.integer(native_stan_data$use_probit_link),
           baseline_case = baseline_case,
           n_covariates = n_covariates,
           native_stan_data = native_stan_data,
           scope = "native baseline marginal accuracy; not predictive accuracy or joint curve matching")

}


##
## ---- Prepare the explicitly selected marginal composite likelihood of the reference accuracy draws:
##


fn_data_for_deriving_the_Jones_marginal_accuracy_prior <- function(reference_Se_draws,
                                                                   reference_Sp_draws,
                                                                   raw_cutoffs,
                                                                   native_prior_bank,
                                                                   objective,
                                                                   target_grain_size,
                                                                   estimate_beta_baseline_mean,
                                                                   estimate_beta_baseline_SD,
                                                                   estimate_raw_scale_baseline_mean,
                                                                   estimate_raw_scale_baseline_SD,
                                                                   fixed_beta_baseline_mean,
                                                                   fixed_beta_baseline_SD,
                                                                   fixed_raw_scale_baseline_mean,
                                                                   fixed_raw_scale_baseline_SD,
                                                                   estimate_boxcox_prior,
                                                                   fixed_boxcox_mean,
                                                                   fixed_boxcox_SD) {

      if (!identical(objective, "marginal_composite")) {
            stop("This tool implements only the explicitly requested marginal_composite objective.");
      }

      if (!is.numeric(target_grain_size) || length(target_grain_size) != 1
            || !is.finite(target_grain_size) || target_grain_size < 1
            || target_grain_size != as.integer(target_grain_size)) {
            stop("Supply target_grain_size explicitly, in target accuracy points per minimum task.");
      }

      masks <- list(estimate_beta_baseline_mean = estimate_beta_baseline_mean,
                    estimate_beta_baseline_SD = estimate_beta_baseline_SD,
                    estimate_raw_scale_baseline_mean = estimate_raw_scale_baseline_mean,
                    estimate_raw_scale_baseline_SD = estimate_raw_scale_baseline_SD)
      fixed <- list(fixed_beta_baseline_mean = fixed_beta_baseline_mean,
                    fixed_beta_baseline_SD = fixed_beta_baseline_SD,
                    fixed_raw_scale_baseline_mean = fixed_raw_scale_baseline_mean,
                    fixed_raw_scale_baseline_SD = fixed_raw_scale_baseline_SD)

      for (nm in names(masks)) {
            value <- masks[[nm]]

            if (!is.logical(value) || !is.null(dim(value)) || length(value) != 2 || anyNA(value)) {
                  stop(paste0(nm, " must explicitly select status one and status two with TRUE or FALSE."));
            }

            masks[[nm]] <- as.integer(value)
      }

      for (nm in names(fixed)) {
            value <- fixed[[nm]]

            if (!is.numeric(value) || !is.null(dim(value)) || length(value) != 2
                  || any(!is.finite(value)) || (endsWith(nm, "_SD") && any(value <= 0))) {
                  stop(paste0(nm, " must supply two finite values, with strictly positive SDs."));
            }
      }

      if (!is.matrix(reference_Se_draws) || !is.numeric(reference_Se_draws)
            || !is.matrix(reference_Sp_draws) || !is.numeric(reference_Sp_draws)
            || !identical(dim(reference_Se_draws), dim(reference_Sp_draws))
            || nrow(reference_Se_draws) < 2 || ncol(reference_Se_draws) < 1) {
            stop(paste0("Reference Se and Sp must be paired numeric matrices: ",
                        "complete draws in rows, cutoffs in columns."));
      }

      target_probability <- c(as.vector(reference_Sp_draws), as.vector(reference_Se_draws))

      if (any(!is.finite(target_probability)) || any(target_probability <= 0 | target_probability >= 1)) {
            stop("Every target probability must be finite and strictly inside zero and one; none are clipped.");
      }

      if (!is.numeric(raw_cutoffs) || !is.null(dim(raw_cutoffs))
            || length(raw_cutoffs) != ncol(reference_Se_draws)
            || any(!is.finite(raw_cutoffs)) || any(diff(raw_cutoffs) <= 0)) {
            stop("raw_cutoffs must give the strictly increasing raw scores of the reference matrix columns.");
      }

      required_bank <- c("n_thr", "cts_thr_values", "n_bank", "cutpoint_bank", "lambda_bank", "lambda_upper_bank",
                         "raw_scale_bank",
                         "bank_raw_scale_mean", "bank_raw_scale_SD", "bank_boxcox_mean", "bank_boxcox_SD",
                         "softplus", "use_probit_link", "baseline_case", "n_covariates", "scope")

      if (!is.list(native_prior_bank) || !all(required_bank %in% names(native_prior_bank))) {
            stop("Prepare native_prior_bank with fn_Jones_native_prior_integration_bank first.");
      }

      cutoff_index <- match(raw_cutoffs, native_prior_bank$cts_thr_values)

      if (anyNA(cutoff_index)) {
            stop("Every reference cutoff must occur exactly in the full native Jones threshold grid.");
      }

      if (ncol(reference_Se_draws) > 1) {
            following <- seq.int(2, ncol(reference_Se_draws))
            preceding <- following - 1

            if (any(reference_Se_draws[, following, drop = FALSE]
                        >= reference_Se_draws[, preceding, drop = FALSE])
                  || any(reference_Sp_draws[, following, drop = FALSE]
                        <= reference_Sp_draws[, preceding, drop = FALSE])) {
                  stop(paste0("Each paired reference curve must have decreasing Se and increasing Sp ",
                              "across raw cutoffs."));
            }
      }

      if (!is.logical(estimate_boxcox_prior) || !is.null(dim(estimate_boxcox_prior))
            || length(estimate_boxcox_prior) != 2 || anyNA(estimate_boxcox_prior)) {
            stop("estimate_boxcox_prior must select both classes explicitly with TRUE or FALSE.");
      }

      if (!is.numeric(fixed_boxcox_mean) || !is.null(dim(fixed_boxcox_mean)) || length(fixed_boxcox_mean) != 2
            || any(!is.finite(fixed_boxcox_mean)) || !is.numeric(fixed_boxcox_SD)
            || !is.null(dim(fixed_boxcox_SD)) || length(fixed_boxcox_SD) != 2
            || any(!is.finite(fixed_boxcox_SD)) || any(fixed_boxcox_SD <= 0)) {
            stop("Supply both fixed Box-Cox means and positive standard deviations explicitly.");
      }

      if (any(estimate_boxcox_prior)
            && !any(cutoff_index > 1 & cutoff_index < native_prior_bank$n_thr)) {
            stop("Endpoint cutoffs do not identify the Box-Cox prior; supply its fixed prior inputs.");
      }

      n_reference <- nrow(reference_Se_draws)
      n_selected <- ncol(reference_Se_draws)
      n_per_status <- n_reference * n_selected
      target_status <- rep(c(1, 2), each = n_per_status)
      target_cutoff_index <- rep(rep(cutoff_index, each = n_reference), times = 2)
      target_draw_row <- rep(seq_len(n_reference), times = 2 * n_selected)
      stan_data <- native_prior_bank[setdiff(required_bank, c("baseline_case", "n_covariates", "scope"))]

      stan_data$n_target <- length(target_probability)
      stan_data$target_status <- target_status
      stan_data$target_cutoff_index <- target_cutoff_index
      stan_data$target_probability <- target_probability
      stan_data$marginal_composite_likelihood <- 1
      stan_data$estimate_boxcox_prior <- as.integer(estimate_boxcox_prior)
      stan_data$fixed_boxcox_mean <- fixed_boxcox_mean
      stan_data$fixed_boxcox_SD <- fixed_boxcox_SD
      stan_data$target_grain_size <- as.integer(target_grain_size)
      stan_data <- c(stan_data, masks, fixed)

      target_map <- data.frame(draw_row = target_draw_row,
                               quantity = ifelse(target_status == 1, "Sp", "Se"),
                               raw_cutoff = native_prior_bank$cts_thr_values[target_cutoff_index],
                               native_cutoff_index = target_cutoff_index)

      initial <- list(
            estimated_beta_baseline_mean = unname(fixed_beta_baseline_mean[estimate_beta_baseline_mean]),
            estimated_beta_baseline_SD = unname(fixed_beta_baseline_SD[estimate_beta_baseline_SD]),
            estimated_raw_scale_baseline_mean = unname(
                  fixed_raw_scale_baseline_mean[estimate_raw_scale_baseline_mean]),
            estimated_raw_scale_baseline_SD = unname(fixed_raw_scale_baseline_SD[estimate_raw_scale_baseline_SD]),
            estimated_boxcox_mean = unname(fixed_boxcox_mean[estimate_boxcox_prior]),
            estimated_boxcox_SD = unname(fixed_boxcox_SD[estimate_boxcox_prior]))

      message(paste0("Jones target: ", n_reference, " paired reference rows at ", n_selected,
                     " raw cutoffs; ", length(target_probability), " marginal observations. ",
                     "This is a finite-bank marginal composite likelihood, not a joint curve likelihood."))

      list(stan_data = stan_data,
           initial = initial,
           n_free = sum(lengths(initial)),
           target_map = target_map,
           baseline_case = native_prior_bank$baseline_case,
           n_covariates = native_prior_bank$n_covariates,
           native_stan_data = native_prior_bank$native_stan_data,
           scope = native_prior_bank$scope,
           objective = objective,
           approximation = "self-normalised importance integration over the supplied native prior bank")

}






















