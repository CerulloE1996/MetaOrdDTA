


##
## ---- Draws of one parameter family (from the main, transformed or generated-quantities trace) --------------------------------------
##
#' R_fn_draws_of_one_parameter_family
#' @keywords internal
R_fn_draws_of_one_parameter_family <- function(model_summary_and_trace_obj,
                                               family_name) {

        traces_as_arrays <- model_summary_and_trace_obj$internal_obj$traces$traces_as_arrays
        for (trace in list(traces_as_arrays$trace_main,
                           traces_as_arrays$trace_tp_wo_log_lik,
                           traces_as_arrays$trace_gq)) {

              if (is.null(trace)) next
              if (!any(startsWith(dimnames(trace)[[3]], paste0(family_name, "[")))) next
              trace_of_family <- MetaOrdDTA:::extract_params_from_array_batch(debugging = FALSE,
                                                                              array_3d = trace,
                                                                              param_strings_vec = c(family_name),
                                                                              condition = "exact_match")
              n_draws <- dim(trace_of_family)[1] * dim(trace_of_family)[2]
              draws <- matrix(trace_of_family, nrow = n_draws)
              colnames(draws) <- dimnames(trace_of_family)[[3]]
              return(draws)

        }
        stop(paste0("no draws of ", family_name, " in the traces of the fit"))

}



##
## ---- Summary Se / Sp of a meta-regression, weighted over covariate rows (MA and NMA) -----------------------------------------------
##
#' R_fn_weighted_summary_accuracy
#' @description
#' One summary Se / Sp per test and threshold from a meta-regression (MA or NMA, ordinal bivariate models):
#' Se / Sp are computed at each covariate row with the formula of the Stan models' Se_baseline / Sp_baseline
#' (random effects at 0),
#' then averaged over the rows with the weights, draw by draw, on the probability scale.
#' weighting = "equal": every row the same weight (e.g. rows = the studies: each study counts once);
#' "study_size": rows weighted by their study's size (weights = study sizes);
#' "user": the user's rows and weights (a target population); one row with weight 1 = Se_baseline / Sp_baseline at
#' that row (post hoc, without refitting).
#' Giving every test the same rows gives every test the same population, so that the tests can be compared.
#' @param model_summary_and_trace_obj The fitted MetaOrdDTA model (with its traces).
#' @param network TRUE for NMA (X_nd, X_d, weights, n_covariates_nd, n_covariates_d: lists / vectors over the
#' tests),
#' FALSE for MA (one matrix / vector / number each).
#' @param random_thresholds TRUE for the random-cutpoint models (summary cutpoints C_mu), FALSE for the fixed ones
#' (C_array for NMA, C for MA).
#' @param use_probit_link TRUE (Phi) or FALSE (inv_logit), as in the fit.
#' @param n_thr Number of thresholds (NMA: one per test).
#' @param n_covariates_nd Number of covariates of the non-diseased class, intercept included (NMA: one per test).
#' @param n_covariates_d Number of covariates of the diseased class, intercept included (NMA: one per test).
#' @param X_nd Covariate rows of the non-diseased class (rows = studies or target profiles; no missing values).
#' @param X_d Covariate rows of the diseased class.
#' @param weighting "equal", "study_size" or "user".
#' @param weights Study sizes (weighting "study_size") or the user's weights (weighting "user"); NULL for "equal".
#' @return A tibble: test, threshold, Se and Sp mean, median and 95% interval, weighting.
#' @export
R_fn_weighted_summary_accuracy <- function(model_summary_and_trace_obj,
                                           network,
                                           random_thresholds,
                                           use_probit_link,
                                           n_thr,
                                           n_covariates_nd,
                                           n_covariates_d,
                                           X_nd,
                                           X_d,
                                           weighting,
                                           weights) {

        stopifnot(weighting %in% c("equal", "study_size", "user"))
        stopifnot(weighting == "equal" || !is.null(weights))
        ##
        ## ---- MA as one test --------------------------------------------------------------------------------------------------------
        ##
        if (network == FALSE) {
          X_nd <- list(X_nd)
          X_d <- list(X_d)
          if (!is.null(weights)) weights <- list(weights)
        }
        n_index_tests <- length(X_nd)
        ##
        ## ---- Draws of beta_mu and of the summary cutpoints -------------------------------------------------------------------------
        ##
        beta_mu_draws <- R_fn_draws_of_one_parameter_family(
                           model_summary_and_trace_obj = model_summary_and_trace_obj,
                           family_name = "beta_mu")
        C_name <- if (random_thresholds) "C_mu" else if (network) "C_array" else "C"
        C_draws <- R_fn_draws_of_one_parameter_family(
                     model_summary_and_trace_obj = model_summary_and_trace_obj,
                     family_name = C_name)
        link_fn <- if (use_probit_link) stats::pnorm else stats::plogis
        ##
        rows_of_output <- list()
        for (t in 1:n_index_tests) {

              ##
              ## ---- Rows and weights of test t --------------------------------------------------------------------------------------
              ##
              n_cov_nd <- n_covariates_nd[t]
              n_cov_d <- n_covariates_d[t]
              X_nd_t <- as.matrix(X_nd[[t]])[, 1:n_cov_nd, drop = FALSE]
              X_d_t <- as.matrix(X_d[[t]])[, 1:n_cov_d, drop = FALSE]
              if (any(X_nd_t == 999) || any(X_d_t == 999)) stop("missing covariate values (999) in X_nd / X_d")
              n_rows <- nrow(X_d_t)
              w <- if (weighting == "equal") rep(1, n_rows) else weights[[t]]
              stopifnot(length(w) == n_rows, all(w >= 0), sum(w) > 0)
              w <- w / sum(w)
              ##
              ## ---- Xbeta of each row and draw (draws x rows), as Xbeta_baseline in the Stan models ---------------------------------
              ##
              if (network) {
                beta_names_nd <- paste0("beta_mu[", t, ",1,", 1:n_cov_nd, "]")
                beta_names_d <- paste0("beta_mu[", t, ",2,", 1:n_cov_d, "]")
              } else {
                beta_names_nd <- paste0("beta_mu[1,", 1:n_cov_nd, "]")
                beta_names_d <- paste0("beta_mu[2,", 1:n_cov_d, "]")
              }
              Xbeta_nd <- beta_mu_draws[, beta_names_nd, drop = FALSE] %*% t(X_nd_t)
              Xbeta_d <- beta_mu_draws[, beta_names_d, drop = FALSE] %*% t(X_d_t)
              ##
              n_thr_t <- n_thr[t]
              for (k in 1:n_thr_t) {

                    if (network) {
                      C_nd <- C_draws[, paste0(C_name, "[1,", t, ",", k, "]")]
                      C_d <- C_draws[, paste0(C_name, "[2,", t, ",", k, "]")]
                    } else {
                      C_nd <- C_draws[, paste0(C_name, "[1,", k, "]")]
                      C_d <- C_draws[, paste0(C_name, "[2,", k, "]")]
                    }
                    ##
                    ## ---- Fp = link(-(C_nd - Xbeta_nd)), Sp = 1 - Fp, Se = link(-(C_d - Xbeta_d)), weighted -------------------------
                    ##
                    Fp_of_rows <- link_fn(-(C_nd - Xbeta_nd))
                    Se_of_rows <- link_fn(-(C_d - Xbeta_d))
                    Sp_weighted <- 1.0 - as.vector(Fp_of_rows %*% w)
                    Se_weighted <- as.vector(Se_of_rows %*% w)
                    q_Se <- stats::quantile(Se_weighted, probs = c(0.025, 0.5, 0.975), names = FALSE)
                    q_Sp <- stats::quantile(Sp_weighted, probs = c(0.025, 0.5, 0.975), names = FALSE)
                    rows_of_output[[length(rows_of_output) + 1]] <- tibble::tibble(test = t,
                                                                                   threshold = k,
                                                                                   Se_mean = mean(Se_weighted),
                                                                                   Se_median = q_Se[2],
                                                                                   Se_lower = q_Se[1],
                                                                                   Se_upper = q_Se[3],
                                                                                   Sp_mean = mean(Sp_weighted),
                                                                                   Sp_median = q_Sp[2],
                                                                                   Sp_lower = q_Sp[1],
                                                                                   Sp_upper = q_Sp[3],
                                                                                   weighting = weighting)

              }

        }
        return(dplyr::bind_rows(rows_of_output))

}
