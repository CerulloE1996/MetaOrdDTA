##
## ---- Independent evaluation of the marginal Jones accuracy likelihood and its gradient:
##
## The supplied parameter point is evaluated without optimisation or a transform Jacobian.
## Positive free standard deviations use the logarithmic unconstrained parameterisation.
## The bank and fixed/free masks are the data of the corresponding Stan model.


fn_log_density_and_gradient_of_the_Jones_marginal_accuracy_prior <- function(theta,
                                                                           stan_data) {

    if (!is.list(theta) || !is.list(stan_data)) {
        stop("Supply a named parameter-point list and the corresponding Stan-data list.")
    }

    fn_numeric_vector <- function(value,
                                  size,
                                  label) {

        if (!is.numeric(value) || !is.null(dim(value)) || length(value) != size
            || any(!is.finite(value))) {
            stop(paste0(label, " must be a finite numeric vector of length ", size, "."))
        }

        return(as.numeric(value))

    }

    fn_integer <- function(value,
                           lower,
                           label) {

        value <- fn_numeric_vector(value, 1, label)
        if (value != floor(value) || value < lower) {
            stop(paste0(label, " must be an integer at least ", lower, "."))
        }

        return(value)

    }

    fn_mask <- function(value,
                        size,
                        label) {

        if (is.logical(value) && !anyNA(value)) value <- as.numeric(value)
        value <- fn_numeric_vector(value, size, label)
        if (any(!value %in% c(0, 1))) stop(paste0(label, " must contain only zero or one."))

        return(value)

    }

    fn_log_sum_exp <- function(value) {

        largest <- max(value)
        if (!is.finite(largest)) stop("The supplied point has no finite integration mass.")

        return(largest + log(sum(exp(value - largest))))

    }


    ##
    ## ---- Validate the native bank and the target map:
    ##

    n_bank <- fn_integer(stan_data$n_bank, 2, "n_bank")
    n_thr <- fn_integer(stan_data$n_thr, 1, "n_thr")
    n_target <- fn_integer(stan_data$n_target, 1, "n_target")
    target_grain_size <- fn_integer(stan_data$target_grain_size, 1, "target_grain_size")
    softplus <- fn_mask(stan_data$softplus, 1, "softplus")
    use_probit_link <- fn_mask(stan_data$use_probit_link, 1, "use_probit_link")

    if (!identical(as.numeric(stan_data$marginal_composite_likelihood), 1)) {
        stop("The objective must be explicitly identified as a marginal composite likelihood.")
    }

    cts_thr_values <- fn_numeric_vector(stan_data$cts_thr_values, n_thr, "cts_thr_values")
    if (any(diff(cts_thr_values) <= 0)) stop("The full native threshold grid must be increasing.")

    cutpoint_bank <- stan_data$cutpoint_bank
    if (!is.numeric(cutpoint_bank) || !identical(dim(cutpoint_bank), as.integer(c(2, n_bank, n_thr)))
        || any(!is.finite(cutpoint_bank))) {
        stop("cutpoint_bank must have dimensions status by bank draw by native cutoff.")
    }

    raw_scale_bank <- stan_data$raw_scale_bank
    lambda_bank <- stan_data$lambda_bank
    lambda_upper_bank <- stan_data$lambda_upper_bank
    for (bank in list(raw_scale_bank, lambda_bank, lambda_upper_bank)) {
        if (!is.numeric(bank) || !identical(dim(bank), as.integer(c(n_bank, 2)))
            || any(!is.finite(bank))) {
            stop("The raw-scale and lambda banks must be finite n_bank by two matrices.")
        }
    }
    if (any(lambda_bank < -5 | lambda_bank > 5)
        || any(lambda_upper_bank < -5 | lambda_upper_bank > 5)) {
        stop("The native lambda banks lie outside their support.")
    }
    for (status in seq_len(2)) {
        for (j in seq_len(n_bank)) {
            if (any(diff(cutpoint_bank[status, j, ]) < 0)) {
                stop("The native cutpoint bank must retain its original cutoff order.")
            }
        }
    }

    bank_raw_scale_mean <- fn_numeric_vector(stan_data$bank_raw_scale_mean, 2, "bank_raw_scale_mean")
    bank_raw_scale_SD <- fn_numeric_vector(stan_data$bank_raw_scale_SD, 2, "bank_raw_scale_SD")
    bank_boxcox_mean <- fn_numeric_vector(stan_data$bank_boxcox_mean, 2, "bank_boxcox_mean")
    bank_boxcox_SD <- fn_numeric_vector(stan_data$bank_boxcox_SD, 2, "bank_boxcox_SD")
    if (any(bank_raw_scale_SD <= 0) || any(bank_boxcox_SD <= 0)) {
        stop("Every proposal standard deviation must be positive.")
    }

    target_status <- fn_numeric_vector(stan_data$target_status, n_target, "target_status")
    target_cutoff_index <- fn_numeric_vector(stan_data$target_cutoff_index, n_target, "target_cutoff_index")
    target_probability <- fn_numeric_vector(stan_data$target_probability, n_target, "target_probability")
    if (any(!target_status %in% c(1, 2)) || !all(c(1, 2) %in% target_status)) {
        stop("Supply specificity targets for status one and sensitivity targets for status two.")
    }
    if (any(target_cutoff_index != floor(target_cutoff_index))
        || any(target_cutoff_index < 1 | target_cutoff_index > n_thr)) {
        stop("Every target cutoff index must refer to the full native threshold grid.")
    }
    if (any(target_probability <= 0 | target_probability >= 1)) {
        stop("Target probabilities must be strictly between zero and one; no clipping is used.")
    }


    ##
    ## ---- Merge fixed and free inputs in the declared Stan parameter order:
    ##

    blocks <- c("beta_baseline_mean", "beta_baseline_SD", "raw_scale_baseline_mean",
                "raw_scale_baseline_SD", "boxcox_mean", "boxcox_SD")
    positive_block <- c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE)
    estimate_boxcox_prior <- fn_mask(stan_data$estimate_boxcox_prior, 2, "estimate_boxcox_prior")
    masks <- list(fn_mask(stan_data$estimate_beta_baseline_mean, 2, "estimate_beta_baseline_mean"),
                  fn_mask(stan_data$estimate_beta_baseline_SD, 2, "estimate_beta_baseline_SD"),
                  fn_mask(stan_data$estimate_raw_scale_baseline_mean, 2, "estimate_raw_scale_baseline_mean"),
                  fn_mask(stan_data$estimate_raw_scale_baseline_SD, 2, "estimate_raw_scale_baseline_SD"),
                  estimate_boxcox_prior,
                  estimate_boxcox_prior)
    prior <- vector("list", length(blocks))
    names(prior) <- blocks

    for (b in seq_along(blocks)) {
        fixed_name <- paste0("fixed_", blocks[b])
        parameter_name <- paste0("estimated_", blocks[b])
        prior[[b]] <- fn_numeric_vector(stan_data[[fixed_name]], length(masks[[b]]), fixed_name)
        if (positive_block[b] && any(prior[[b]] <= 0)) stop(paste0(fixed_name, " must be positive."))
        n_free <- sum(masks[[b]])
        free_value <- if (n_free == 0 && is.null(theta[[parameter_name]])) numeric(0) else theta[[parameter_name]]
        free_value <- fn_numeric_vector(free_value, n_free, parameter_name)
        if (positive_block[b] && any(free_value <= 0)) stop(paste0(parameter_name, " must be positive."))
        prior[[b]][masks[[b]] == 1] <- free_value
    }

    for (status in seq_len(2)) {
        has_interior_target <- any(target_status == status & target_cutoff_index > 1
                                   & target_cutoff_index < n_thr)
        if (estimate_boxcox_prior[status] == 1 && !has_interior_target) {
            stop("Endpoint cutoffs do not identify lambda; supply its fixed prior inputs explicitly.")
        }
    }


    ##
    ## ---- Independent normal-density ratios and analytic normal scores:
    ##

    scale_bank <- if (softplus == 1) {
        1.4426950408889634 * (pmax(raw_scale_bank, 0) + log1p(exp(-abs(raw_scale_bank))))
    } else {
        exp(raw_scale_bank)
    }
    if (any(!is.finite(scale_bank)) || any(scale_bank <= 0)) {
        stop("Every bank draw must produce a finite positive native scale.")
    }

    target_z <- if (use_probit_link == 1) {
        stats::qnorm(target_probability, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
    } else {
        stats::qlogis(target_probability, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE)
    }
    native_link_argument <- ifelse(target_status == 1, -target_z, target_z)
    if (use_probit_link == 1 && any(native_link_argument < -37.5 | native_link_argument > 8.25)) {
        stop("A target lies outside the native interior Phi branch; endpoint masses are not fitted.")
    }
    target_log_link_density <- if (use_probit_link == 1) {
        stats::dnorm(target_z, mean = 0, sd = 1, log = TRUE)
    } else {
        stats::dlogis(target_z, location = 0, scale = 1, log = TRUE)
    }

    log_weight <- matrix(0, nrow = n_bank, ncol = 2)
    bank_weight <- matrix(0, nrow = n_bank, ncol = 2)
    raw_mean_score <- raw_SD_score <- lambda_mean_score <- lambda_SD_score <- log_weight
    log_weight_sum <- numeric(2)
    score <- lapply(prior, function(value) {

        return(rep(0, length(value)))

    })

    for (status in seq_len(2)) {
        raw <- raw_scale_bank[, status]
        lambda <- lambda_bank[, status]
        lambda_upper <- lambda_upper_bank[, status]
        raw_mean <- prior$raw_scale_baseline_mean[status]
        raw_SD <- prior$raw_scale_baseline_SD[status]
        lambda_mean <- prior$boxcox_mean[status]
        lambda_SD <- prior$boxcox_SD[status]

        log_weight[, status] <- stats::dnorm(raw, mean = raw_mean, sd = raw_SD, log = TRUE) -
            stats::dnorm(raw, mean = bank_raw_scale_mean[status], sd = bank_raw_scale_SD[status], log = TRUE) +
            stats::dnorm(lambda, mean = lambda_mean, sd = lambda_SD, log = TRUE) +
            stats::dnorm(lambda_upper, mean = lambda_mean, sd = lambda_SD, log = TRUE) -
            stats::dnorm(lambda, mean = bank_boxcox_mean[status], sd = bank_boxcox_SD[status], log = TRUE) -
            stats::dnorm(lambda_upper, mean = bank_boxcox_mean[status],
                         sd = bank_boxcox_SD[status], log = TRUE)
        log_weight_sum[status] <- fn_log_sum_exp(log_weight[, status])
        bank_weight[, status] <- exp(log_weight[, status] - log_weight_sum[status])

        raw_mean_score[, status] <- (raw - raw_mean) / raw_SD^2
        raw_SD_score[, status] <- (raw - raw_mean)^2 / raw_SD^3 - 1 / raw_SD
        lambda_mean_score[, status] <- ((lambda - lambda_mean) + (lambda_upper - lambda_mean)) /
            lambda_SD^2
        lambda_SD_score[, status] <- ((lambda - lambda_mean)^2 + (lambda_upper - lambda_mean)^2) /
            lambda_SD^3 - 2 / lambda_SD
    }


    ##
    ## ---- Conditional location integration and the composite score:
    ##

    marginal_log_density <- numeric(n_target)
    conditional_integration_ESS <- numeric(n_target)

    for (i in seq_len(n_target)) {
        status <- target_status[i]
        accuracy_sign <- if (status == 1) -1 else 1
        beta <- cutpoint_bank[status, , target_cutoff_index[i]] +
            accuracy_sign * scale_bank[, status] * target_z[i]
        beta_mean <- prior$beta_baseline_mean[status]
        beta_SD <- prior$beta_baseline_SD[status]
        log_components <- log_weight[, status] + log(scale_bank[, status]) - target_log_link_density[i] +
            stats::dnorm(beta, mean = beta_mean, sd = beta_SD, log = TRUE)
        component_log_sum <- fn_log_sum_exp(log_components)
        conditional_weight <- exp(log_components - component_log_sum)
        difference_weight <- conditional_weight - bank_weight[, status]

        marginal_log_density[i] <- component_log_sum - log_weight_sum[status]
        conditional_integration_ESS[i] <- 1 / sum(conditional_weight^2)
        score$beta_baseline_mean[status] <- score$beta_baseline_mean[status] +
            sum(conditional_weight * (beta - beta_mean) / beta_SD^2)
        score$beta_baseline_SD[status] <- score$beta_baseline_SD[status] +
            sum(conditional_weight * ((beta - beta_mean)^2 / beta_SD^3 - 1 / beta_SD))
        score$raw_scale_baseline_mean[status] <- score$raw_scale_baseline_mean[status] +
            sum(difference_weight * raw_mean_score[, status])
        score$raw_scale_baseline_SD[status] <- score$raw_scale_baseline_SD[status] +
            sum(difference_weight * raw_SD_score[, status])
        score$boxcox_mean[status] <- score$boxcox_mean[status] +
            sum(difference_weight * lambda_mean_score[, status])
        score$boxcox_SD[status] <- score$boxcox_SD[status] +
            sum(difference_weight * lambda_SD_score[, status])
    }


    ##
    ## ---- Map the score to the unconstrained Stan coordinate order without a Jacobian term:
    ##

    gradient_unconstrained <- numeric(0)
    for (b in seq_along(blocks)) {
        selected <- which(masks[[b]] == 1)
        if (length(selected) == 0) next
        component <- score[[b]][selected]
        if (positive_block[b]) component <- component * prior[[b]][selected]
        names(component) <- paste0("estimated_", blocks[b], "[", seq_along(selected), "]")
        gradient_unconstrained <- c(gradient_unconstrained, component)
    }

    if (any(!is.finite(marginal_log_density)) || any(!is.finite(gradient_unconstrained))) {
        stop("The supplied parameter point produces a non-finite log density or analytic gradient.")
    }

    return(list(log_density = sum(marginal_log_density),
                gradient_unconstrained = gradient_unconstrained,
                gradient_constrained = score,
                prior_inputs = prior,
                marginal_log_density = marginal_log_density,
                normalised_importance_weight = bank_weight,
                importance_weight_ESS = 1 / colSums(bank_weight^2),
                conditional_integration_ESS = conditional_integration_ESS))

}






















