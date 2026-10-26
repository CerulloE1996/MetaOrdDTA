##
## ---- Independent analytic derivative of the native-bank accuracy-summary discrepancy:
##
## Differentiate the normalised weights through their Gaussian and Dirichlet scores.
## Positive free parameters use their log-coordinate derivatives, without a transform Jacobian.
## No finite-difference step, target distribution or prior values are chosen here.
##


fn_log_density_and_gradient_of_accuracy_summary_derivation <- function(point,
                                                                       request) {

    free_blocks <- c("free_normal_mean", "free_normal_SD", "free_alpha")
    prior_blocks <- c("normal_mean", "normal_SD", "alpha")

    if (!is.list(point) || !identical(names(point), free_blocks)) {
        stop("Supply the three free-parameter vectors in the Stan declaration order, including empty vectors.")
    }

    prior <- request$prior_values
    free_names <- character(0)
    positive_multiplier <- numeric(0)

    for (block_i in seq_along(free_blocks)) {

        free_name <- free_blocks[block_i]
        prior_name <- prior_blocks[block_i]
        value <- point[[free_name]]
        selected <- request$estimate[[prior_name]]

        if (!is.numeric(value) || !is.null(dim(value)) || length(value) != sum(selected)
            || any(!is.finite(value)) || (block_i > 1 && any(value <= 0))) {
            stop(paste0(free_name, " must contain exactly the finite selected entries, in native order."))
        }

        prior[[prior_name]][selected] <- value

        if (length(value) > 0) {
            free_names <- c(free_names, paste0(free_name, "[", seq_along(value), "]"))
            positive_multiplier <- c(positive_multiplier, if (block_i == 1) rep(1, length(value)) else value)
        }

    }

    fn_select_prior_inputs_for_derivation(prior_values = prior,
                                          estimate = request$estimate,
                                          prior_bank = request$prior_bank)

    bank <- request$prior_bank
    reference <- bank$reference_priors
    n_bank <- nrow(bank$accuracy_draws)
    log_weight <- numeric(n_bank)

    for (k in seq_len(ncol(bank$normal_draws))) {

        x <- bank$normal_draws[, k]
        log_weight <- log_weight - 0.5 * (
              ((x - prior$normal_mean[k]) / prior$normal_SD[k])^2
            - ((x - reference$normal_mean[k]) / reference$normal_SD[k])^2)

    }

    if (length(prior$alpha) > 0) {
        log_weight <- log_weight + as.vector(bank$log_category_probability %*%
                                             (prior$alpha - reference$alpha))
    }

    if (any(!is.finite(log_weight))) {
        stop("The supplied point produces nonfinite log weights; no clipping or replacement is used.")
    }

    weight <- exp(log_weight - max(log_weight))
    weight <- weight / sum(weight)
    score <- matrix(numeric(0), nrow = n_bank, ncol = 0)

    for (k in which(request$estimate$normal_mean)) {
        derivative <- (bank$normal_draws[, k] - prior$normal_mean[k]) / prior$normal_SD[k]^2
        score <- cbind(score, derivative)
    }

    for (k in which(request$estimate$normal_SD)) {
        derivative <- (bank$normal_draws[, k] - prior$normal_mean[k])^2 / prior$normal_SD[k]^3
        score <- cbind(score, derivative)
    }

    if (any(request$estimate$alpha)) {
        score <- cbind(score, bank$log_category_probability[, request$estimate$alpha, drop = FALSE])
    }

    lower <- sweep(x = bank$accuracy_draws, MARGIN = 2, STATS = request$targets$lower, FUN = "<=")
    upper <- sweep(x = bank$accuracy_draws, MARGIN = 2, STATS = request$targets$upper, FUN = "<=")
    features <- cbind(bank$accuracy_draws, lower, upper)
    achieved <- as.vector(crossprod(weight, features))
    tail <- (1 - request$targets$interval_level) / 2
    desired <- c(request$targets$mean, tail, 1 - tail)
    residual <- achieved - desired
    discrepancy_weight <- as.vector(request$discrepancy_weights)
    centred <- sweep(x = features, MARGIN = 2, STATS = achieved, FUN = "-")
    summary_gradient <- crossprod(centred, score * weight)
    constrained_gradient <- as.vector(-2 * crossprod(discrepancy_weight * residual, summary_gradient))
    unconstrained_gradient <- constrained_gradient * positive_multiplier

    list(log_density = -sum(discrepancy_weight * residual^2),
         gradient_unconstrained = stats::setNames(unconstrained_gradient, free_names),
         gradient_constrained = stats::setNames(constrained_gradient, free_names),
         importance_ESS = 1 / sum(weight^2),
         largest_importance_weight = max(weight))

}






















