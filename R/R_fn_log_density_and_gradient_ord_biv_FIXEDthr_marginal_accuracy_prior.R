##
## ---- Independent evaluation of the ord_biv fixed-cutpoint accuracy likelihood and gradient:
##
## The supplied parameter point is evaluated without optimisation or a transform Jacobian. Positive free
## parameters use the logarithmic unconstrained parameterisation. The bank, masks and groups are the data of
## derive_ord_biv_FIXEDthr_marginal_accuracy_prior_target_draws.stan.
##
fn_log_density_and_gradient_ord_biv_FIXEDthr_marginal_accuracy_prior <- function(theta,
                                                                                 stan_data) {

    if (!is.list(theta) || !is.list(stan_data)) stop("Supply the parameter-point and Stan-data lists.")
    fn_log_sum_exp <- function(value) {

        largest <- max(value)
        if (!is.finite(largest)) stop("The supplied point has no finite integration mass.")
        return(largest + log(sum(exp(value - largest))))

    }
    n_thr <- as.integer(stan_data$n_thr)
    n_cat <- n_thr + 1L
    n_bank <- as.integer(stan_data$n_bank)
    n_target <- as.integer(stan_data$n_target)
    anchor <- as.numeric(stan_data$anchor)
    cutpoint_bank <- stan_data$cutpoint_bank
    bank_alpha <- stan_data$bank_dirichlet_alpha
    estimate_SD <- as.numeric(stan_data$estimate_beta_baseline_SD)
    fixed_SD <- as.numeric(stan_data$fixed_beta_baseline_SD)
    estimate_alpha <- as.numeric(stan_data$estimate_dirichlet_alpha)
    n_groups <- as.integer(stan_data$n_alpha_groups)
    alpha_group <- stan_data$alpha_group
    scale_of_the_hyperprior <- as.numeric(stan_data$hyperprior_alpha_scale)
    target_status <- as.integer(stan_data$target_status)
    target_cutoff_index <- as.integer(stan_data$target_cutoff_index)
    target_probability <- as.numeric(stan_data$target_probability)
    ##
    ## ---- the free parameters in the declared Stan order:
    ##
    SD <- fixed_SD
    free_SD <- theta$estimated_beta_baseline_SD
    if (sum(estimate_SD) == 0 && is.null(free_SD)) free_SD <- numeric(0)
    if (length(free_SD) != sum(estimate_SD) || any(free_SD <= 0)) stop("Supply positive free location SDs.")
    SD[estimate_SD == 1] <- free_SD
    gamma <- theta$estimated_alpha_group_concentration
    n_free_groups <- sum(n_groups[estimate_alpha == 1])
    if (n_free_groups == 0 && is.null(gamma)) gamma <- numeric(0)
    if (length(gamma) != n_free_groups || any(gamma <= 0)) stop("Supply positive free concentrations.")
    alpha <- vector("list", 2)
    group_offset <- 0
    for (c in 1:2) {
        alpha[[c]] <- as.numeric(bank_alpha[c, ])
        if (estimate_alpha[c] == 1) {
            alpha[[c]] <- gamma[group_offset + as.integer(alpha_group[c, ])]
            group_offset <- group_offset + n_groups[c]
        }
    }
    ##
    ## ---- the bank's category log probabilities, the importance weights and the composite likelihood:
    ##
    log_density <- 0
    gradient_SD <- numeric(2)
    gradient_alpha <- vector("list", 2)
    target_z <- stats::qnorm(target_probability)
    target_log_link_density <- stats::dnorm(target_z, log = TRUE)
    for (c in 1:2) {
        C <- cutpoint_bank[c, , , drop = TRUE] - anchor[c]
        if (!is.matrix(C)) C <- matrix(C, nrow = n_bank)
        log_F <- stats::pnorm(C, log.p = TRUE)
        log_S <- stats::pnorm(C, lower.tail = FALSE, log.p = TRUE)
        log_phi <- matrix(NA_real_, n_bank, n_cat)
        log_phi[, 1] <- log_F[, 1]
        log_phi[, n_cat] <- log_S[, n_thr]
        if (n_thr > 1) {
            for (k in 2:n_thr) {
                use_S <- C[, k - 1] >= 0
                log_phi[use_S, k] <- log_S[use_S, k - 1] + log1p(-exp(log_S[use_S, k] - log_S[use_S, k - 1]))
                log_phi[!use_S, k] <- log_F[!use_S, k] + log1p(-exp(log_F[!use_S, k - 1] - log_F[!use_S, k]))
            }
        }
        if (any(!is.finite(log_phi))) stop("A bank category log probability is not finite.")
        log_w <- as.numeric(log_phi %*% (alpha[[c]] - as.numeric(bank_alpha[c, ])))
        log_w_sum <- fn_log_sum_exp(log_w)
        bank_weight <- exp(log_w - log_w_sum)
        weighted_log_phi <- as.numeric(crossprod(bank_weight, log_phi))
        gradient_alpha[[c]] <- numeric(n_cat)
        sign <- if (c == 1) -1 else 1
        for (i in which(target_status == c)) {
            k <- target_cutoff_index[i]
            beta <- C[, k] + anchor[c] + sign * target_z[i]
            standardised <- (beta - anchor[c]) / SD[c]
            log_components <- log_w - target_log_link_density[i] - log(SD[c]) - 0.5 * log(2 * pi) -
                  0.5 * standardised^2
            log_components_sum <- fn_log_sum_exp(log_components)
            log_density <- log_density + log_components_sum - log_w_sum
            p <- exp(log_components - log_components_sum)
            gradient_SD[c] <- gradient_SD[c] + sum(p * (standardised^2 / SD[c] - 1 / SD[c]))
            gradient_alpha[[c]] <- gradient_alpha[[c]] + as.numeric(crossprod(p, log_phi)) - weighted_log_phi
        }
    }
    ##
    ## ---- the half-normal hyperprior of the free concentrations, then the unconstrained-scale gradient:
    ##
    log_density <- log_density + sum(stats::dnorm(gamma, mean = 0, sd = scale_of_the_hyperprior, log = TRUE))
    gradient_gamma <- numeric(n_free_groups)
    group_offset <- 0
    for (c in 1:2) {
        if (estimate_alpha[c] == 1) {
            for (g in seq_len(n_groups[c])) {
                in_group <- as.integer(alpha_group[c, ]) == g
                gradient_gamma[group_offset + g] <- sum(gradient_alpha[[c]][in_group])
            }
            group_offset <- group_offset + n_groups[c]
        }
    }
    gradient_gamma <- gradient_gamma - gamma / scale_of_the_hyperprior^2
    gradient <- c(gradient_SD[estimate_SD == 1] * free_SD, gradient_gamma * gamma)
    names(gradient) <- c(paste0("estimated_beta_baseline_SD[", seq_len(sum(estimate_SD)), "]"),
                         paste0("estimated_alpha_group_concentration[", seq_len(n_free_groups), "]"))
    return(list(log_density = log_density,
                gradient_unconstrained = gradient))

}






















