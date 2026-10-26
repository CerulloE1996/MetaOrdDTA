##
## ---- Connect native meta-analysis prior draws to the common accuracy-prior derivation functions:
##
## Keep the full native prior draw rows paired. Raw cutoff values are supplied explicitly, never
## inferred from a category index. The fixed-cutpoint ordinal-bivariate
## anchor and every omitted prior remain fixed.
## Jones has no Dirichlet alpha; its Box-Cox prior is unchanged by this adapter.
##


fn_accuracy_prior_bank_from_MetaOrdDTA <- function(native_prior_draws,
                                                   native_stan_data,
                                                   model_code,
                                                   raw_cutoff_grid,
                                                   targets,
                                                   model_file,
                                                   draw_source) {

    supported <- c("DTA_MA_ord_biv_FIXEDthr", "DTA_MA_ord_HSROC_RANDthr_kappa", "DTA_MA_Jones")

    if (!is.character(model_code) || length(model_code) != 1 || !(model_code %in% supported)) {
        stop(paste0("This adapter covers ordinal-bivariate fixed-cutpoint, ",
                    "ordinal-HSROC random-cutpoint and Jones meta-analysis models only."))
    }

    if (!is.character(model_file) || length(model_file) != 1 || !file.exists(model_file)
        || basename(model_file) != paste0(model_code, ".stan")) {
        stop("Supply the actual original native Stan file corresponding to model_code.")
    }

    if (!is.character(draw_source) || length(draw_source) != 1 || is.na(draw_source) || !nzchar(draw_source)) {
        stop("Identify the actual native prior-fit draw source explicitly.")
    }

    if (!is.matrix(native_prior_draws) || !is.numeric(native_prior_draws)
        || nrow(native_prior_draws) < 2 || is.null(colnames(native_prior_draws))
        || anyDuplicated(colnames(native_prior_draws))) {
        stop("Supply the full named numeric matrix of paired native prior-only draws.")
    }

    if (!is.list(native_stan_data)
        || !all(c("K_fold_CV_indicator", "n_studies", "n_thr") %in% names(native_stan_data))
        || !is.numeric(native_stan_data$n_studies) || length(native_stan_data$n_studies) != 1
        || !is.finite(native_stan_data$n_studies) || native_stan_data$n_studies < 1
        || is.null(native_stan_data$K_fold_CV_indicator)
        || anyNA(native_stan_data$K_fold_CV_indicator)
        || length(native_stan_data$K_fold_CV_indicator) != native_stan_data$n_studies
        || any(native_stan_data$K_fold_CV_indicator != 0)) {
        stop("Supply the matching native Stan data that marks every study as prior-only.")
    }

    n_thr <- native_stan_data$n_thr

    if (!is.numeric(n_thr) || length(n_thr) != 1 || !is.finite(n_thr)
        || n_thr < 1 || n_thr != as.integer(n_thr)
        || !is.numeric(raw_cutoff_grid) || !is.null(dim(raw_cutoff_grid))
        || length(raw_cutoff_grid) != n_thr || any(!is.finite(raw_cutoff_grid))
        || any(diff(raw_cutoff_grid) <= 0)) {
        stop("Supply the actual full increasing raw cutoff grid, one value per native threshold.")
    }

    if (model_code == "DTA_MA_Jones"
        && !identical(as.numeric(raw_cutoff_grid), as.numeric(native_stan_data$cts_thr_values))) {
        stop("The Jones raw cutoff grid must equal the actual cts_thr_values in its Stan data.")
    }

    targets <- fn_accuracy_prior_targets(targets = targets)

    if (length(unique(targets$test)) != 1
        || !all(targets$estimand %in% c("baseline", "baseline_predictive"))) {
        stop("Supply one test and explicitly choose baseline or baseline_predictive accuracy for every row.")
    }

    cutoff_index <- match(targets$raw_cutoff, raw_cutoff_grid)

    if (anyNA(cutoff_index)) {
        stop("Every requested raw cutoff must occur exactly in the supplied native grid.")
    }

    fn_columns <- function(column_names) {

        if (!all(column_names %in% colnames(native_prior_draws))) {
            stop(paste0("The native prior bank lacks: ",
                        paste(setdiff(column_names, colnames(native_prior_draws)), collapse = ", "), "."))
        }

        value <- native_prior_draws[, column_names, drop = FALSE]

        if (any(!is.finite(value))) {
            stop("The requested native prior-draw columns contain nonfinite values.")
        }

        value

    }

    quantity_names <- paste0(targets$quantity,
                             ifelse(targets$estimand == "baseline", "_baseline[", "_baseline_pred["),
                             cutoff_index, "]")
    accuracy_draws <- fn_columns(column_names = quantity_names)
    colnames(accuracy_draws) <- targets$target_id
    normal_names <- mean_names <- SD_names <- character(0)
    reference_mean <- reference_SD <- numeric(0)
    mean_can_change <- logical(0)

    normal_factors <- if (model_code == "DTA_MA_ord_biv_FIXEDthr") {
        "beta_mu"
    } else {
        c("beta_mu", "raw_scale_mu")
    }

    for (factor in normal_factors) {

        mean_field <- paste0("prior_", factor, "_mean")
        SD_field <- paste0("prior_", factor, "_SD")

        if (model_code == "DTA_MA_ord_HSROC_RANDthr_kappa") {

            n_cov <- native_stan_data$n_covariates

            if (!is.numeric(n_cov) || length(n_cov) != 1 || !is.finite(n_cov)
                || n_cov < 1 || n_cov != as.integer(n_cov)) {
                stop("The native number of covariates must be a positive integer.")
            }

            cols <- paste0(factor, "[", seq_len(n_cov), "]")
            mu <- native_stan_data[[mean_field]]
            sigma <- native_stan_data[[SD_field]]

            if (length(mu) != n_cov || length(sigma) != n_cov) {
                stop("The native Gaussian coefficient priors do not match the active covariates.")
            }

            normal_names <- c(normal_names, cols)
            mean_names <- c(mean_names, paste0(mean_field, "[", seq_len(n_cov), "]"))
            SD_names <- c(SD_names, paste0(SD_field, "[", seq_len(n_cov), "]"))
            reference_mean <- c(reference_mean, mu)
            reference_SD <- c(reference_SD, sigma)
            mean_can_change <- c(mean_can_change, rep(TRUE, n_cov))

        } else {

            for (status in seq_len(2)) {

                n_cov <- if (status == 1) native_stan_data$n_covariates_nd else native_stan_data$n_covariates_d

                if (!is.numeric(n_cov) || length(n_cov) != 1 || !is.finite(n_cov)
                    || n_cov < 1 || n_cov != as.integer(n_cov)) {
                    stop("Each native number of covariates must be a positive integer.")
                }

                cols <- paste0(factor, "[", status, ",", seq_len(n_cov), "]")
                mu_matrix <- native_stan_data[[mean_field]]
                SD_matrix <- native_stan_data[[SD_field]]

                if (!is.matrix(mu_matrix) || !is.matrix(SD_matrix)
                    || nrow(mu_matrix) != 2 || nrow(SD_matrix) != 2
                    || ncol(mu_matrix) < n_cov || ncol(SD_matrix) < n_cov) {
                    stop("The native Gaussian coefficient prior matrices do not match the active covariates.")
                }

                normal_names <- c(normal_names, cols)
                mean_names <- c(mean_names, paste0(mean_field, "[", status, ",", seq_len(n_cov), "]"))
                SD_names <- c(SD_names, paste0(SD_field, "[", status, ",", seq_len(n_cov), "]"))
                reference_mean <- c(reference_mean, mu_matrix[status, seq_len(n_cov)])
                reference_SD <- c(reference_SD, SD_matrix[status, seq_len(n_cov)])
                allowed <- rep(TRUE, n_cov)

                if (model_code == "DTA_MA_ord_biv_FIXEDthr" && factor == "beta_mu") {
                    allowed[1] <- FALSE
                }

                mean_can_change <- c(mean_can_change, allowed)

            }

        }

    }

    normal_draws <- fn_columns(column_names = normal_names)
    reference_mean <- stats::setNames(reference_mean, normal_names)
    reference_SD <- stats::setNames(reference_SD, normal_names)
    mean_can_change <- stats::setNames(mean_can_change, normal_names)
    alpha_names <- character(0)
    reference_alpha <- numeric(0)
    log_category_probability <- matrix(numeric(0), nrow = nrow(native_prior_draws), ncol = 0)

    if (model_code == "DTA_MA_ord_HSROC_RANDthr_kappa") {

        phi_names <- paste0("dirichlet_phi[", seq_len(n_thr + 1), "]")
        phi <- fn_columns(column_names = phi_names)

        if (any(phi <= 0 | phi > 1)) {
            stop("The native dirichlet_phi bank must contain strictly positive simplex entries.")
        }

        alpha_names <- paste0("prior_dirichlet_alpha[", seq_len(n_thr + 1), "]")
        reference_alpha <- native_stan_data$prior_dirichlet_alpha
        log_category_probability <- log(phi)
        colnames(log_category_probability) <- alpha_names

    } else if (model_code == "DTA_MA_ord_biv_FIXEDthr") {

        log_blocks <- vector(mode = "list", length = 2)

        for (status in seq_len(2)) {

            C_names <- paste0("C[", status, ",", seq_len(n_thr), "]")
            C <- fn_columns(column_names = C_names) - native_stan_data$prior_beta_mu_mean[status, 1]
            log_F <- matrix(stats::pnorm(q = as.vector(C), lower.tail = TRUE, log.p = TRUE),
                             nrow = nrow(C), ncol = ncol(C))
            log_S <- matrix(stats::pnorm(q = as.vector(C), lower.tail = FALSE, log.p = TRUE),
                             nrow = nrow(C), ncol = ncol(C))
            log_p <- matrix(NA_real_, nrow = nrow(C), ncol = n_thr + 1)
            log_p[, 1] <- log_F[, 1]
            log_p[, n_thr + 1] <- log_S[, n_thr]

            if (n_thr > 1) {
                for (k in seq.int(from = 2, to = n_thr)) {

                    upper_tail <- C[, k - 1] >= 0
                    large <- ifelse(upper_tail, log_S[, k - 1], log_F[, k])
                    small <- ifelse(upper_tail, log_S[, k], log_F[, k - 1])
                    log_p[, k] <- large + log(-expm1(small - large))

                }
            }

            alpha_names <- c(alpha_names,
                             paste0("prior_dirichlet_alpha[", status, ",", seq_len(n_thr + 1), "]"))
            alpha_data <- native_stan_data$prior_dirichlet_alpha
            alpha_values <- if (is.matrix(alpha_data)) alpha_data[status, ] else alpha_data[[status]]

            if (!is.numeric(alpha_values) || length(alpha_values) != n_thr + 1) {
                stop("Each native ordinal-bivariate Dirichlet prior must contain one alpha per full category.")
            }

            reference_alpha <- c(reference_alpha, alpha_values)
            log_blocks[[status]] <- log_p

        }

        log_category_probability <- do.call(what = cbind, args = log_blocks)
        colnames(log_category_probability) <- alpha_names

    }

    reference_alpha <- stats::setNames(reference_alpha, alpha_names)
    can_change <- list(normal_mean = mean_can_change,
                       normal_SD = stats::setNames(rep(TRUE, length(normal_names)), normal_names),
                       alpha = stats::setNames(rep(TRUE, length(alpha_names)), alpha_names))
    reference_priors <- list(normal_mean = reference_mean, normal_SD = reference_SD, alpha = reference_alpha)
    parameter_map <- list(normal_mean = stats::setNames(mean_names, normal_names),
                          normal_SD = stats::setNames(SD_names, normal_names),
                          alpha = stats::setNames(alpha_names, alpha_names))
    mapping_columns <- c("target_id", "test", "quantity", "raw_cutoff", "estimand")
    provenance <- list(model_file = model_file,
                       model_data = native_stan_data,
                       draw_source = draw_source,
                       fixed_inputs = native_stan_data,
                       target_map = targets[, mapping_columns, drop = FALSE],
                       can_change = can_change,
                       parameter_map = parameter_map,
                       model_code = model_code,
                       raw_cutoff_grid = raw_cutoff_grid)

    fn_accuracy_prior_integration_bank(accuracy_draws = accuracy_draws,
                                       normal_draws = normal_draws,
                                       log_category_probability = log_category_probability,
                                       reference_priors = reference_priors,
                                       provenance = provenance)

}


fn_MetaOrdDTA_prior_inputs_from_accuracy_derivation <- function(derived_prior_values,
                                                               prior_bank) {

    source_data <- prior_bank$provenance$model_data
    map <- prior_bank$provenance$parameter_map
    output <- list()
    changes <- list()
    change_i <- 0

    if (!is.list(map) || !identical(names(map), c("normal_mean", "normal_SD", "alpha"))) {
        stop("The bank lacks the native prior-input mapping.")
    }

    for (block in names(map)) {

        value <- derived_prior_values[[block]]

        if (!is.numeric(value) || !identical(names(value), names(prior_bank$reference_priors[[block]]))
            || any(!is.finite(value)) || (block != "normal_mean" && any(value <= 0))) {
            stop(paste0("Supply every named derived ", block, " entry exactly as returned."))
        }

        allowed <- prior_bank$provenance$can_change[[block]]
        reference <- prior_bank$reference_priors[[block]]

        if (any(value[!allowed] != reference[!allowed])) {
            stop(paste0(block, ": the exported inputs change an entry held fixed by the native bank."))
        }

        for (k in seq_along(value)) {

            destination <- unname(map[[block]][k])
            field <- sub(pattern = "\\[.*$", replacement = "", x = destination)
            inside <- sub(pattern = "^.*\\[([0-9,]+)\\]$", replacement = "\\1", x = destination)

            if (!grepl(pattern = "^[0-9,]+$", x = inside) || is.null(source_data[[field]])) {
                stop("The native prior destination is not an explicitly indexed Stan-data entry.")
            }

            index <- as.integer(strsplit(x = inside, split = ",", fixed = TRUE)[[1]])

            if (is.null(output[[field]])) {
                output[[field]] <- source_data[[field]]
            }

            if (length(index) == 1 && is.numeric(output[[field]]) && is.null(dim(output[[field]]))) {
                old <- output[[field]][index[1]]
                output[[field]][index[1]] <- value[k]
            } else if (length(index) == 2 && is.matrix(output[[field]])) {
                old <- output[[field]][index[1], index[2]]
                output[[field]][index[1], index[2]] <- value[k]
            } else if (length(index) == 2 && is.list(output[[field]])) {
                old <- output[[field]][[index[1]]][index[2]]
                output[[field]][[index[1]]][index[2]] <- value[k]
            } else {
                stop("The native prior destination does not match its actual scalar, vector or array layout.")
            }

            change_i <- change_i + 1
            changes[[change_i]] <- data.frame(native_input = destination,
                                              original = old,
                                              derived = unname(value[k]),
                                              changed = old != value[k])

        }
    }

    table <- do.call(what = rbind, args = changes)
    print(table, row.names = FALSE)

    list(prior_inputs = output,
         changes = table,
         model_code = prior_bank$provenance$model_code,
         native_model_check_passed = NA)

}






















