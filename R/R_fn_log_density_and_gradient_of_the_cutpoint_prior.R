##
## ---- Independent density and derivative of the joint induced-Dirichlet cutpoint target:
##
## The target draws, anchor and link are held fixed. The concentration derivative is analytic;
## the positive-parameter log transform contributes its chain rule but no transform Jacobian.
## No prior values are estimated and no finite-difference step or probability replacement is used.
##


fn_log_density_and_gradient_of_the_cutpoint_prior <- function(point,
                                                              stan_data) {

      if (!is.list(point) || !identical(names(point), "alpha") || !is.list(stan_data)) {
            stop("Supply the named alpha parameter block and its matching Stan-data list.")
      }

      for (field in c("n_draws", "n_thr")) {

            value <- stan_data[[field]]

            if (!is.numeric(value) || length(value) != 1 || !is.null(dim(value))
                  || !is.finite(value) || value != floor(value)) {
                  stop(paste0(field, " must be a finite integer scalar from the Stan data."))
            }

      }

      n_draws <- stan_data$n_draws
      n_thr <- stan_data$n_thr

      if (n_draws < 2 || n_thr < 1) {
            stop("At least two complete draws and one cutoff are required, as in the Stan model.")
      }

      n_cat <- n_thr + 1
      alpha <- point$alpha
      cutpoints <- stan_data$cutpoint_draws
      anchor <- stan_data$anchor
      use_probit_link <- stan_data$use_probit_link

      if (!is.numeric(alpha) || !is.null(dim(alpha)) || length(alpha) != n_cat
            || any(!is.finite(alpha)) || any(alpha <= 0)) {
            stop("alpha must contain one finite positive concentration per target category.")
      }

      if (!is.matrix(cutpoints) || !is.numeric(cutpoints)
            || !identical(dim(cutpoints), as.integer(c(n_draws, n_thr)))
            || any(!is.finite(cutpoints))) {
            stop("cutpoint_draws must be the finite n_draws by n_thr matrix supplied to Stan.")
      }

      if (n_thr > 1 && any(cutpoints[, 2:n_thr, drop = FALSE]
                           <= cutpoints[, 1:(n_thr - 1), drop = FALSE])) {
            stop("Each complete target cutpoint vector must be strictly increasing.")
      }

      if (!is.numeric(anchor) || length(anchor) != 1 || !is.null(dim(anchor)) || !is.finite(anchor)) {
            stop("The cutpoint anchor must be one finite scalar from the Stan data.")
      }

      if (!is.numeric(use_probit_link) || length(use_probit_link) != 1 || !is.null(dim(use_probit_link))
            || !is.finite(use_probit_link) || !(use_probit_link %in% c(0, 1))) {
            stop("use_probit_link must be exactly zero or one, as supplied to Stan.")
      }

      anchored_cutpoints <- cutpoints - anchor

      if (any(!is.finite(anchored_cutpoints))) {
            stop("The anchored target cutpoints are not finite; no values are replaced.")
      }

      if (use_probit_link == 1) {

            log_F <- stats::pnorm(q = anchored_cutpoints,
                                   mean = 0,
                                   sd = 1,
                                   lower.tail = TRUE,
                                   log.p = TRUE)
            log_S <- stats::pnorm(q = anchored_cutpoints,
                                   mean = 0,
                                   sd = 1,
                                   lower.tail = FALSE,
                                   log.p = TRUE)
            log_density_at_cutpoints <- stats::dnorm(x = anchored_cutpoints,
                                                       mean = 0,
                                                       sd = 1,
                                                       log = TRUE)

      } else {

            log_F <- stats::plogis(q = anchored_cutpoints,
                                    location = 0,
                                    scale = 1,
                                    lower.tail = TRUE,
                                    log.p = TRUE)
            log_S <- stats::plogis(q = anchored_cutpoints,
                                    location = 0,
                                    scale = 1,
                                    lower.tail = FALSE,
                                    log.p = TRUE)
            log_density_at_cutpoints <- log_F + log_S

      }

      dim(log_F) <- c(n_draws, n_thr)
      dim(log_S) <- c(n_draws, n_thr)
      dim(log_density_at_cutpoints) <- c(n_draws, n_thr)

      log_probability <- matrix(data = 0, nrow = n_draws, ncol = n_cat)
      log_probability[, 1] <- log_F[, 1]
      log_probability[, n_cat] <- log_S[, n_thr]

      if (n_thr > 1) {
            for (k in 2:n_thr) {

                  use_survival <- anchored_cutpoints[, k - 1] >= 0
                  log_larger <- ifelse(use_survival, log_S[, k - 1], log_F[, k])
                  log_smaller <- ifelse(use_survival, log_S[, k], log_F[, k - 1])

                  if (any(!is.finite(log_larger)) || any(!is.finite(log_smaller))
                        || any(log_smaller >= log_larger)) {
                        stop("A target category has no finite positive mass; no probabilities are replaced.")
                  }

                  log_probability[, k] <- log_larger + log(-expm1(log_smaller - log_larger))

            }
      }

      log_cutpoint_Jacobian <- rowSums(log_density_at_cutpoints)

      if (any(!is.finite(log_probability)) || any(!is.finite(log_cutpoint_Jacobian))) {
            stop("The target category probabilities or cutpoint Jacobian are not finite.")
      }

      first_draw <- matrix(data = log_probability[1, ],
                             nrow = n_draws,
                             ncol = n_cat,
                             byrow = TRUE,
                             dimnames = NULL)

      if (all(log_probability == first_draw)) {
            stop("Identical target draws have no finite maximum-likelihood Dirichlet concentration.")
      }

      log_normalising_constant <- lgamma(sum(alpha)) - sum(lgamma(alpha))
      log_density <- n_draws * log_normalising_constant +
            sum(sweep(x = log_probability, MARGIN = 2, STATS = alpha - 1, FUN = "*")) +
            sum(log_cutpoint_Jacobian)
      gradient_constrained <- n_draws * (digamma(sum(alpha)) - digamma(alpha)) +
            colSums(log_probability)
      gradient_unconstrained <- alpha * gradient_constrained

      if (!is.finite(log_density) || any(!is.finite(gradient_unconstrained))) {
            stop("The supplied saved point has a non-finite density or gradient; no replacement is used.")
      }

      parameter_names <- paste0("alpha[", seq_len(n_cat), "]")

      list(log_density = log_density,
           gradient_unconstrained = stats::setNames(gradient_unconstrained, parameter_names),
           gradient_constrained = stats::setNames(gradient_constrained, parameter_names))

}






















