##
## ---- Prepare a joint cutpoint target from sensitivity or specificity draws:
##


fn_data_for_deriving_the_cutpoint_prior <- function(target_accuracy_draws,
                                                    raw_cutoffs,
                                                    accuracy_measure,
                                                    anchor,
                                                    use_probit_link,
                                                    target_is_conditional_on_location_and_scale) {

      if (!is.logical(target_is_conditional_on_location_and_scale)
            || length(target_is_conditional_on_location_and_scale) != 1
            || is.na(target_is_conditional_on_location_and_scale)
            || !target_is_conditional_on_location_and_scale) {
            stop(paste0("This component needs accuracy draws conditional on the supplied location ",
                        "and unit scale. ",
                        "Marginal model-prior draws need the model-specific translation."));
      }

      if (!is.matrix(target_accuracy_draws) || !is.numeric(target_accuracy_draws)) {
            stop("target_accuracy_draws must be a numeric matrix: draws in rows, cutoffs in columns.");
      }

      n_draws <- nrow(target_accuracy_draws)
      n_thr <- ncol(target_accuracy_draws)

      if (n_draws < 2 || n_thr < 1) {
            stop("At least two complete target draws and one cutoff are required.");
      }

      if (any(!is.finite(target_accuracy_draws))
            || any(target_accuracy_draws <= 0 | target_accuracy_draws >= 1)) {
            stop("Target draws must be finite and strictly between 0 and 1; no draws are clipped.");
      }

      if (!is.numeric(raw_cutoffs) || !is.null(dim(raw_cutoffs)) || length(raw_cutoffs) != n_thr
            || any(!is.finite(raw_cutoffs)) || any(diff(raw_cutoffs) <= 0)) {
            stop("raw_cutoffs must give the strictly increasing raw scores of the matrix columns.");
      }

      if (!is.character(accuracy_measure) || length(accuracy_measure) != 1
            || is.na(accuracy_measure) || !(accuracy_measure %in% c("Se", "Sp"))) {
            stop("accuracy_measure must be exactly Se or Sp, for one test and one disease status.");
      }

      if (!is.numeric(anchor) || length(anchor) != 1 || !is.finite(anchor)) {
            stop("anchor must be one explicitly supplied finite number.");
      }

      if (!is.numeric(use_probit_link) || length(use_probit_link) != 1
            || !is.finite(use_probit_link) || !(use_probit_link %in% c(0, 1))) {
            stop("use_probit_link must be supplied explicitly: 1 for probit, 0 for logit.");
      }

      ##
      ## ---- Preserve each complete target curve; do not sort or pair draws independently:
      ##

      if (n_thr > 1) {
            differences <- target_accuracy_draws[, 2:n_thr, drop = FALSE] -
                  target_accuracy_draws[, 1:(n_thr - 1), drop = FALSE]

            target_is_not_ordered <- if (accuracy_measure == "Se") {
                  any(differences >= 0)
            } else {
                  any(differences <= 0)
            }

            if (target_is_not_ordered) {
                  stop("Every target curve must be strictly ordered across the supplied raw cutoffs.");
            }
      }

      ##
      ## ---- Invert the unit-scale link at the explicitly supplied cutpoint anchor:
      ##

      target_cutpoints <- if (use_probit_link == 1) {
            anchor + stats::qnorm(p = target_accuracy_draws,
                                  mean = 0,
                                  sd = 1,
                                  lower.tail = accuracy_measure == "Sp",
                                  log.p = FALSE)
      } else {
            anchor + stats::qlogis(p = target_accuracy_draws,
                                   location = 0,
                                   scale = 1,
                                   lower.tail = accuracy_measure == "Sp",
                                   log.p = FALSE)
      }

      dim(target_cutpoints) <- c(n_draws, n_thr)

      if (any(!is.finite(target_cutpoints))) {
            stop("The transformed target cutpoints must be finite; no values are replaced.");
      }

      if (all(target_cutpoints == matrix(data = target_cutpoints[1, ],
                                         nrow = n_draws,
                                         ncol = n_thr,
                                         byrow = TRUE,
                                         dimnames = NULL))) {
            stop("Identical target draws have no finite maximum-likelihood Dirichlet concentration.");
      }

      ##
      ## ---- This target is the cutpoint law, not a marginal over additional random locations or scales:
      ##

      stan_data <- list(n_draws = n_draws,
                        n_thr = n_thr,
                        cutpoint_draws = unname(target_cutpoints),
                        anchor = anchor,
                        use_probit_link = as.integer(use_probit_link))

      message(paste0("Cutpoint target: ", n_draws, " complete ", accuracy_measure,
                     " draws at raw cutoffs ", paste0(raw_cutoffs, collapse = ", "),
                     "; anchor = ", anchor, "; use_probit_link = ", use_probit_link, "."))

      list(stan_data = stan_data,
           raw_cutoffs = raw_cutoffs,
           accuracy_measure = accuracy_measure,
           alpha_scope = "intervals between supplied cutoffs; no concentrations for omitted categories",
           scope = "joint cutpoint distribution under the supplied anchor and unit-scale link")

}






















