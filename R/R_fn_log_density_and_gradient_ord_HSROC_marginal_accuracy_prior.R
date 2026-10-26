##
## ---- Independent evaluation of the native ord_HSROC marginal composite likelihood and gradient:
##
## The supplied point is evaluated without optimisation or a transform Jacobian. Positive free
## concentrations and standard deviations use logarithmic unconstrained coordinates. The location
## is integrated analytically; the paired native cutpoint, simplex and raw-scale bank is held fixed.


fn_log_density_and_gradient_ord_HSROC_marginal_accuracy_prior <- function(theta,
                                                                           stan_data) {

    if (!is.list(theta) || !is.list(stan_data)) stop("Supply the parameter-point and Stan-data lists.")

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
        if (value < lower || value != floor(value)) stop(paste0(label, " is not a valid integer."))
        return(value)

    }

    fn_mask <- function(value,
                        size,
                        label) {

        if (is.logical(value) && !anyNA(value)) value <- as.numeric(value)
        value <- fn_numeric_vector(value, size, label)
        if (any(!value %in% c(0, 1))) stop(paste0(label, " must contain only zero and one."))
        return(value)

    }

    fn_log_sum_exp <- function(value) {

        largest <- max(value)
        if (!is.finite(largest)) stop("The supplied point has no finite integration mass.")
        return(largest + log(sum(exp(value - largest))))

    }

    n_bank <- fn_integer(stan_data$n_bank, 2, "n_bank")
    n_thr <- fn_integer(stan_data$n_thr, 1, "n_thr")
    n_target <- fn_integer(stan_data$n_target, 1, "n_target")
    target_grain_size <- fn_integer(stan_data$target_grain_size, 1, "target_grain_size")
    softplus <- fn_mask(stan_data$softplus, 1, "softplus")
    use_probit_link <- fn_mask(stan_data$use_probit_link, 1, "use_probit_link")
    if (!identical(as.numeric(stan_data$marginal_composite_likelihood), 1)) {
        stop("Explicitly identify the objective as a marginal composite likelihood.")
    }

    grid <- fn_numeric_vector(stan_data$raw_cutoff_grid, n_thr, "raw_cutoff_grid")
    if (any(diff(grid) <= 0)) stop("The complete raw-cutoff grid must be strictly increasing.")
    cutpoints <- stan_data$cutpoint_bank
    phi <- stan_data$phi_bank
    if (!is.numeric(cutpoints) || !identical(dim(cutpoints), as.integer(c(n_bank, n_thr)))
        || any(!is.finite(cutpoints)) || !is.numeric(phi)
        || !identical(dim(phi), as.integer(c(n_bank, n_thr + 1)))
        || any(!is.finite(phi)) || any(phi <= 0 | phi > 1)) {
        stop("Retain finite paired native C_mu and strictly positive dirichlet_phi bank matrices.")
    }
    for (j in seq_len(n_bank)) {
        if (any(diff(cutpoints[j, ]) < 0)) stop("The native cutpoint bank must retain its cutoff order.")
    }
    raw <- fn_numeric_vector(stan_data$raw_scale_bank, n_bank, "raw_scale_bank")
    proposal_mean <- fn_numeric_vector(stan_data$bank_raw_scale_mean, 1, "bank_raw_scale_mean")
    proposal_SD <- fn_numeric_vector(stan_data$bank_raw_scale_SD, 1, "bank_raw_scale_SD")
    proposal_alpha <- fn_numeric_vector(stan_data$bank_dirichlet_alpha, n_thr + 1, "bank_dirichlet_alpha")
    if (proposal_SD <= 0 || any(proposal_alpha <= 0)) {
        stop("The native proposal SD and concentrations must be strictly positive.")
    }

    status <- fn_numeric_vector(stan_data$target_status, n_target, "target_status")
    cutoff <- fn_numeric_vector(stan_data$target_cutoff_index, n_target, "target_cutoff_index")
    probability <- fn_numeric_vector(stan_data$target_probability, n_target, "target_probability")
    if (any(!status %in% c(1, 2)) || !all(c(1, 2) %in% status)
        || any(cutoff != floor(cutoff)) || any(cutoff < 1 | cutoff > n_thr)) {
        stop("Supply both statuses and valid indices into the complete native cutoff grid.")
    }
    if (any(probability <= 0 | probability >= 1)) {
        stop("Every target must be strictly inside zero and one; no clipping is used.")
    }

    ##
    ## ---- Retain fixed values and insert free entries in the Stan parameter order:
    ##

    blocks <- c("beta_baseline_mean", "beta_baseline_SD", "raw_scale_baseline_mean",
                "raw_scale_baseline_SD", "dirichlet_alpha")
    positive <- c(FALSE, TRUE, FALSE, TRUE, TRUE)
    masks <- prior <- vector(mode = "list", length = length(blocks))
    names(prior) <- names(masks) <- blocks
    for (b in seq_along(blocks)) {
        size <- if (b == length(blocks)) n_thr + 1 else 1
        masks[[b]] <- fn_mask(stan_data[[paste0("estimate_", blocks[b])]], size, blocks[b])
        prior[[b]] <- fn_numeric_vector(stan_data[[paste0("fixed_", blocks[b])]], size, blocks[b])
        parameter <- paste0("estimated_", blocks[b])
        n_free <- sum(masks[[b]])
        free <- if (n_free == 0 && is.null(theta[[parameter]])) numeric(0) else theta[[parameter]]
        free <- fn_numeric_vector(free, n_free, parameter)
        if (positive[b] && (any(prior[[b]] <= 0) || any(free <= 0))) {
            stop(paste0(blocks[b], " must remain strictly positive."))
        }
        prior[[b]][masks[[b]] == 1] <- free
    }

    ##
    ## ---- Native reciprocal-scale transformations and self-normalised density ratios:
    ##

    signed_raw <- cbind(raw, -raw)
    scale <- if (softplus == 1) {
        1.4426950408889634 * (pmax(signed_raw, 0) + log1p(exp(-abs(signed_raw))))
    } else exp(signed_raw)
    if (any(!is.finite(scale)) || any(scale <= 0) || any(!is.finite(1 / scale))) {
        stop("Every native bank row must give finite positive scales and finite inverse scales.")
    }
    z <- if (use_probit_link == 1) {
        stats::qnorm(probability, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
    } else {
        stats::qlogis(probability, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE)
    }
    native_argument <- ifelse(status == 1, -z, z)
    if (use_probit_link == 1 && any(native_argument < -37.5 | native_argument > 8.25)) {
        stop("A target is outside the native interior Phi branch; endpoint masses are not fitted.")
    }
    log_link_density <- if (use_probit_link == 1) {
        stats::dnorm(z, mean = 0, sd = 1, log = TRUE)
    } else {
        stats::dlogis(z, location = 0, scale = 1, log = TRUE)
    }
    log_phi <- log(phi)
    log_weight <- stats::dnorm(raw, mean = prior$raw_scale_baseline_mean,
                              sd = prior$raw_scale_baseline_SD, log = TRUE) -
        stats::dnorm(raw, mean = proposal_mean, sd = proposal_SD, log = TRUE) +
        as.vector(log_phi %*% (prior$dirichlet_alpha - proposal_alpha))
    log_weight_sum <- fn_log_sum_exp(log_weight)
    bank_weight <- exp(log_weight - log_weight_sum)
    raw_mean_score <- (raw - prior$raw_scale_baseline_mean) / prior$raw_scale_baseline_SD^2
    raw_SD_score <- (raw - prior$raw_scale_baseline_mean)^2 / prior$raw_scale_baseline_SD^3 -
        1 / prior$raw_scale_baseline_SD
    score <- lapply(X = prior,
                     FUN = function(value) {

                         return(rep(0, length(value)))

                     })
    marginal_log_density <- conditional_integration_ESS <- numeric(n_target)

    for (i in seq_len(n_target)) {
        c <- status[i]
        beta <- (if (c == 1) cutpoints[, cutoff[i]] else -cutpoints[, cutoff[i]]) - z[i] / scale[, c]
        residual <- beta - prior$beta_baseline_mean
        beta_SD <- prior$beta_baseline_SD
        log_components <- log_weight - log(scale[, c]) - log_link_density[i] +
            stats::dnorm(beta, mean = prior$beta_baseline_mean, sd = beta_SD, log = TRUE)
        component_sum <- fn_log_sum_exp(log_components)
        conditional_weight <- exp(log_components - component_sum)
        difference <- conditional_weight - bank_weight
        marginal_log_density[i] <- component_sum - log_weight_sum
        conditional_integration_ESS[i] <- 1 / sum(conditional_weight^2)
        score$beta_baseline_mean <- score$beta_baseline_mean + sum(conditional_weight * residual / beta_SD^2)
        score$beta_baseline_SD <- score$beta_baseline_SD +
            sum(conditional_weight * (residual^2 / beta_SD^3 - 1 / beta_SD))
        score$raw_scale_baseline_mean <- score$raw_scale_baseline_mean + sum(difference * raw_mean_score)
        score$raw_scale_baseline_SD <- score$raw_scale_baseline_SD + sum(difference * raw_SD_score)
        score$dirichlet_alpha <- score$dirichlet_alpha + as.vector(crossprod(difference, log_phi))
    }

    gradient <- numeric(0)
    for (b in seq_along(blocks)) {
        selected <- which(masks[[b]] == 1)
        if (length(selected) == 0) next
        component <- score[[b]][selected]
        if (positive[b]) component <- component * prior[[b]][selected]
        names(component) <- paste0("estimated_", blocks[b], "[", seq_along(selected), "]")
        gradient <- c(gradient, component)
    }
    if (any(!is.finite(marginal_log_density)) || any(!is.finite(gradient))) {
        stop("The supplied point produces a nonfinite likelihood or analytic gradient.")
    }

    list(log_density = sum(marginal_log_density),
         gradient_unconstrained = gradient,
         gradient_constrained = score,
         prior_inputs = prior,
         marginal_log_density = marginal_log_density,
         normalised_importance_weight = bank_weight,
         importance_weight_ESS = 1 / sum(bank_weight^2),
         conditional_integration_ESS = conditional_integration_ESS)

}






















