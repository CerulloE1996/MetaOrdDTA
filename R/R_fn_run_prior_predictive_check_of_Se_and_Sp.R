##
## ==== R_fn_run_prior_predictive_check_of_Se_and_Sp.R ===========================================================
##
## ---- The prior predictive check of a MetaOrdDTA model (MA or NMA): the model's OWN Stan file sampled
##      prior-only (every study's K_fold_CV_indicator = 0, so that no study contributes a likelihood
##      term), so that the draws of the summary sensitivity and specificity (Se_baseline, Sp_baseline) of every
##      test at every threshold are their prior, with every transformation, bound, Jacobian and constraint of the
##      model. The prior is never simulated outside the model (e.g. Dirichlet draws pushed through the cutpoint
##      map in R), because only the model itself gives the real induced prior.
##
##      Optionally, beside it, what the induced-Dirichlet alphas alone say at each threshold (fixed-cutpoint
##      models): Beta(sum of alpha above k, sum of alpha at or below k) for Se and the reverse for Sp, by the
##      Dirichlet's aggregation property, at the anchor and without the location prior; labelled as such, it is
##      never the prior itself.
##
#' run_prior_predictive_check_of_Se_and_Sp
#'
#' @description
#' Prior predictive check of the summary sensitivity and specificity of a MetaOrdDTA model: the model's own Stan
#' file is sampled prior-only (every study's \code{K_fold_CV_indicator} set to 0), so the
#' draws of \code{Se_baseline} and \code{Sp_baseline} are their prior, for every test and every threshold. Every
#' transformation, bound, Jacobian and constraint of the model is included, so this is the real induced prior.
#'
#' @param model_prep_obj A \code{MetaOrd_model} object (from \code{MetaOrdDTA::MetaOrd_model$new()}). It is
#'   cloned, so the object given is left as it is.
#' @param priors The priors to check (e.g. \code{model_prep_obj$priors} with some of them changed); NULL means the
#'   priors of \code{model_prep_obj}.
#' @param init_lists_per_chain The initial values of each chain; NULL means those of \code{model_prep_obj},
#'   resized to \code{n_chains}.
#' @param test_names The name of every test, in the order of the data (one name for a meta-analysis of one test,
#'   e.g. \code{"HADS-A"}; one per test for a network meta-analysis, e.g. \code{c("GAD-7", "HADS-A")}).
#' @param cutoffs_of_interest_of_each_test A list with one numeric vector per test, in the order of
#'   \code{test_names} (or named by them): the cutoffs whose prior is shown, each a score at or above which a test
#'   counts as positive, e.g. \code{list("HADS-A" = c(8, 11))}.
#' @param lowest_possible_score_of_each_test The lowest possible score of each test's scale, in the order of
#'   \code{test_names} (0 for HADS-A, GAD-2 and the BAI; 20 for the STAI). The model's threshold k of a test is a
#'   score of its lowest possible score + k or more.
#' @param n_chains The number of chains (2 is enough for a prior predictive check).
#' @param n_burnin The number of burn-in iterations.
#' @param n_iter The number of sampling iterations of each chain.
#' @param seed The seed of the sampler.
#' @param algorithm "NicoStan" (MetaOrdDTA's NicoStan option) or "Stan" (cmdstanr).
#' @param NicoStan_args The arguments of NicoStan's \code{$sample()} for \code{algorithm = "NicoStan"}, e.g.
#'   \code{list(adapt_delta = 0.80)}; NULL for \code{algorithm = "Stan"}. The nuisance blocks are chosen by
#'   MetaOrdDTA itself (fn_nuisance_params_of_a_MetaOrdDTA_model), never given here.
#' @param adapt_delta The target acceptance rate for \code{algorithm = "Stan"}; NULL for "NicoStan".
#' @param max_treedepth The maximum tree depth for \code{algorithm = "Stan"}; NULL for "NicoStan".
#' @param probs The lower and upper probabilities of the prior interval shown, e.g. \code{c(0.025, 0.975)}.
#' @param show_what_the_alphas_alone_say Whether to add, for fixed-cutpoint models, what the induced-Dirichlet
#'   alphas alone say at each threshold (at the anchor, without the location prior).
#'
#' @return A list:
#'   \code{at_the_cutoffs_of_interest} (a tibble with one row per test and cutoff of interest: the prior median
#'   and interval of Se and of Sp, also printed),
#'   \code{plot} (Se and Sp against the score, one panel per test, the cutoffs of interest marked with their
#'   median and interval),
#'   \code{summary} (a tibble with one row per test, quantity and threshold: the score, the prior median, mean,
#'   interval and bulk / tail / SD ESS of Se or Sp),
#'   \code{min_ESS} (the smallest of the bulk, tail and SD ESS over every Se and Sp), \code{min_ESS_of_each_type},
#'   \code{n_divergent_transitions}, \code{seconds_of_the_fit}, \code{draws_of_Se_and_Sp}, \code{draws_array}
#'   (every draw of the fit, e.g. for the cutpoints) and \code{priors} (the priors checked prior-only).
#'
#' @examples
#' \dontrun{
#' model <-  MetaOrdDTA::MetaOrd_model$new( x = x, ... )
#' prior_check <-  model$run_prior_predictive_check_of_Se_and_Sp(
#'                       priors = model$priors,
#'                       init_lists_per_chain = NULL,
#'                       test_names = c("HADS-A"),
#'                       cutoffs_of_interest_of_each_test = list("HADS-A" = c(8, 11)),
#'                       lowest_possible_score_of_each_test = c(0),
#'                       n_chains = 2,
#'                       n_burnin = 500,
#'                       n_iter = 2500,
#'                       seed = 123,
#'                       algorithm = "NicoStan",
#'                       NicoStan_args = list(adapt_delta = 0.80),
#'                       adapt_delta = NULL,
#'                       max_treedepth = NULL,
#'                       probs = c(0.025, 0.975),
#'                       show_what_the_alphas_alone_say = TRUE)
#' prior_check$at_the_cutoffs_of_interest
#' prior_check$plot
#' }
#' @export
# run_prior_predictive_check_of_Se_and_Sp  <-  function(  model_prep_obj,
#                                                         priors = NULL,
#                                                         init_lists_per_chain = NULL,
#                                                         n_chains = 2,
#                                                         n_burnin = 500,
#                                                         n_iter = 2500,
#                                                         seed = 123,
#                                                         algorithm = "NicoStan",
#                                                         NicoStan_args = list(),
#                                                         adapt_delta = 0.80,
#                                                         max_treedepth = 10,
#                                                         test_names = NULL,
#                                                         probs = c(0.025, 0.975),
#                                                         show_what_the_alphas_alone_say = TRUE) {
run_prior_predictive_check_of_Se_and_Sp  <-  function(  model_prep_obj,
                                                        priors,
                                                        init_lists_per_chain,
                                                        test_names,
                                                        cutoffs_of_interest_of_each_test,
                                                        lowest_possible_score_of_each_test,
                                                        n_chains,
                                                        n_burnin,
                                                        n_iter,
                                                        seed,
                                                        algorithm,
                                                        NicoStan_args,
                                                        adapt_delta,
                                                        max_treedepth,
                                                        probs,
                                                        show_what_the_alphas_alone_say) {

        ##
        stopifnot(algorithm %in% c("NicoStan", "Stan"), length(probs) == 2, probs[1] < probs[2])
        stopifnot( is.character(test_names),
                   length(test_names) >= 1,
                   !anyDuplicated(test_names),
                   is.list(cutoffs_of_interest_of_each_test),
                   length(cutoffs_of_interest_of_each_test) == length(test_names),
                   all(vapply(cutoffs_of_interest_of_each_test, is.numeric, logical(1))),
                   is.numeric(lowest_possible_score_of_each_test),
                   length(lowest_possible_score_of_each_test) == length(test_names))
        ## (the cutoffs of interest in the order of test_names; by name when the list is named)
        if (!is.null(names(cutoffs_of_interest_of_each_test))) {
              if (!setequal(names(cutoffs_of_interest_of_each_test), test_names)) {
                    stop(paste0("run_prior_predictive_check_of_Se_and_Sp: the names of ",
                                "cutoffs_of_interest_of_each_test (",
                                paste(names(cutoffs_of_interest_of_each_test), collapse = ", "),
                                ") are not the test names (", paste(test_names, collapse = ", "), ")."))
              }
              cutoffs_of_interest_of_each_test <-  cutoffs_of_interest_of_each_test[test_names]
        }
        ## (the bulk, tail and SD ESS come from NicoStan's C++ diagnostics, for both algorithms)
        if (!requireNamespace("NicoStan", quietly = TRUE)) {
              stop("run_prior_predictive_check_of_Se_and_Sp needs NicoStan (its C++ diagnostics give the ESS).")
        }
        ##
        ## ---- The model clone, its prior-only inputs and initial values:
        ##
        model_obj <-  model_prep_obj$clone(deep = TRUE)
        priors_checked <-  if (is.null(priors)) model_obj$priors else priors
        n_studies <-  length(priors_checked$K_fold_CV_indicator)
        if (n_studies == 0) {
              stop(paste0("run_prior_predictive_check_of_Se_and_Sp: the priors have no K_fold_CV_indicator, ",
                          "so a prior-only check cannot be prepared."))
        }
        ## (no study contributes a likelihood term: Stan_fns_log_lik.stan and the fused models only add the terms
        ##  of the studies whose K_fold_CV_indicator is 1)
        priors_checked$K_fold_CV_indicator <-  rep(0, n_studies)
        ##
        if (is.null(init_lists_per_chain)) {
              init_lists_per_chain <-  MetaOrdDTA:::resize_init_list(
                                            init_lists_per_chain = model_obj$init_lists_per_chain,
                                            n_chains_new = n_chains)
        }
        ##
        message(paste0("\033[36m", "Prior-only check (", n_studies,
                       " studies): ", algorithm, ", ", n_chains, " chains, ", n_burnin, " burn-in and ", n_iter,
                       " sampling iterations, seed ", seed, "\033[0m"))
        ##
        ## ---- The fit:
        ##
        start_of_the_fit <-  Sys.time()
        if (algorithm == "NicoStan") {

              outs_data <-  model_obj$internal_obj$outs_data
              native_data <-  outs_data$stan_data_list
              native_data$box_cox <-  model_obj$advanced_model_options$box_cox
              native_data$softplus <-  model_obj$advanced_model_options$softplus
              native_data$n_studies <-  outs_data$n_studies
              native_data$n_thr <-  outs_data$n_thr
              for (name in names(priors_checked)) {

                    native_data[[name]] <-  priors_checked[[name]]

              }
              validated_data <-  model_obj$internal_obj$validated_prior_check_data
              if (!is.null(validated_data) && !identical(native_data, validated_data)) {
                    stop("The prior-only Stan inputs differ from the validated inputs.")
              }
              native_file <-  model_obj$internal_obj$outs_stan_model_name$stan_model_file_name
              model_obj$internal_obj$outs_stan_sampling <-  MetaOrdDTA:::R_fn_sample_model_with_NicoStan(
                    stan_data_list = native_data,
                    stan_model_file_name = native_file,
                    n_chains = n_chains,
                    init_lists_per_chain = init_lists_per_chain,
                    seed = seed,
                    n_burnin = n_burnin,
                    n_iter = n_iter,
                    NicoStan_args = NicoStan_args)

        } else {
              model_obj$sample( seed = seed,
                                n_burnin = n_burnin,
                                n_iter = n_iter,
                                n_chains = n_chains,
                                priors = priors_checked,
                                init_lists_per_chain = init_lists_per_chain,
                                algorithm = "Stan",
                                adapt_delta = adapt_delta,
                                max_treedepth = max_treedepth)
        }
        seconds_of_the_fit <-  as.numeric(difftime(Sys.time(), start_of_the_fit, units = "secs"))
        outs <-  model_obj$internal_obj$outs_stan_sampling
        ##
        ## ---- The draws of Se_baseline and Sp_baseline (iterations x chains x quantities), and the divergences:
        ##
        if (algorithm == "NicoStan") {
              draws_array <-  outs$draws_array
              if (is.null(outs$divergences$n_divs)) {
                    stop(paste0("run_prior_predictive_check_of_Se_and_Sp: the NicoStan fit has no count of its ",
                                "divergent transitions (outs$divergences$n_divs)."))
              }
              n_divergent_transitions <-  sum(unlist(outs$divergences$n_divs))
        } else {
              draws_array <-  outs$stan_mod_samples$draws(variables = c("Se_baseline", "Sp_baseline"),
                                                          format = "draws_array")
              n_divergent_transitions <-  sum(outs$stan_mod_samples$diagnostic_summary()$num_divergent)
        }
        names_of_the_draws <-  dimnames(draws_array)[[3]]
        names_of_Se_and_Sp <-  grep("^(Se|Sp)_baseline\\[", names_of_the_draws, value = TRUE)
        ##
        ## ---- One row for each test, quantity and threshold (an NMA test's thresholds above its own number of
        ##      thresholds are padded with -1 by the model, and left out):
        ##
        rows <-  list()
        draws_of_each <-  list()
        for (name in names_of_Se_and_Sp) {
              draws <-  as.numeric(draws_array[, , name])
              if (all(draws == -1)) next
              indices <-  as.numeric(strsplit(sub("^.*\\[(.*)\\]$", "\\1", name), ",")[[1]])
              test_index <-  if (length(indices) == 2) indices[1] else 1
              threshold <-  indices[length(indices)]
              quantity <-  sub("_baseline\\[.*$", "", name)
              ## (the bulk and tail ESS - the rank-normalised split ESS - and the SD ESS - the split ESS of the
              ##  centred squares - with NicoStan's C++ diagnostics, as the Paper 3 runners compute them)
              # draws_as_a_list_of_one_matrix <-  list(matrix(draws, nrow = dim(draws_array)[1]))
              # ESS <-  BayesMVP:::Rcpp_compute_MCMC_diagnostics( draws_as_a_list_of_one_matrix,
              #                                                   diagnostic = "split_ESS",
              #                                                   n_threads = 1)$diagnostics[1, 1]
              draws_matrix <-  matrix(draws, nrow = dim(draws_array)[1])
              bulk_and_tail_ESS <-  NicoStan:::Rcpp_compute_MCMC_diagnostics( list(draws_matrix),
                                                                              diagnostic = "split_ESS_rank",
                                                                              n_threads = 1)$diagnostics
              SD_ESS <-  NicoStan:::Rcpp_compute_MCMC_diagnostics( list((draws_matrix - mean(draws_matrix))^2),
                                                                   diagnostic = "split_ESS",
                                                                   n_threads = 1)$diagnostics[1, 1]
              ESS <-  min(bulk_and_tail_ESS[1, 1], bulk_and_tail_ESS[1, 2], SD_ESS)
              rows[[length(rows) + 1]] <-  data.frame( test_index = test_index,
                                                      quantity = quantity,
                                                      threshold = threshold,
                                                      prior_median = stats::median(draws),
                                                      prior_mean = mean(draws),
                                                      prior_lower = unname(stats::quantile(draws, probs[1])),
                                                      prior_upper = unname(stats::quantile(draws, probs[2])),
                                                      ESS_bulk = bulk_and_tail_ESS[1, 1],
                                                      ESS_tail = bulk_and_tail_ESS[1, 2],
                                                      ESS_SD = SD_ESS,
                                                      ESS = ESS,
                                                      stringsAsFactors = FALSE)
              draws_of_each[[name]] <-  draws
        }
        summary_of_the_prior <-  do.call(rbind, rows)
        n_tests <-  max(summary_of_the_prior$test_index)
        # if (is.null(test_names)) test_names <-  paste("test", seq_len(n_tests))
        if (n_tests != length(test_names)) {
              stop(paste0("run_prior_predictive_check_of_Se_and_Sp: the model has ", n_tests, " test(s), but ",
                          length(test_names), " test name(s) were given (", paste(test_names, collapse = ", "),
                          ")."))
        }
        summary_of_the_prior$test <-  factor(test_names[summary_of_the_prior$test_index], levels = test_names)
        ## (the score of each threshold: the model's threshold k of a test is a score of its lowest possible score
        ##  + k or more)
        lowest_score_of_each_row <-  lowest_possible_score_of_each_test[summary_of_the_prior$test_index]
        summary_of_the_prior$score_at_or_above <-  lowest_score_of_each_row + summary_of_the_prior$threshold
        ##
        ## ---- What the alphas alone say (fixed-cutpoint models; MA: a list of 2 vectors, NMA: a list of 2
        ##      matrices with one row per test):
        ##
        summary_of_the_prior$alphas_alone_median <-  NA_real_
        summary_of_the_prior$alphas_alone_lower  <-  NA_real_
        summary_of_the_prior$alphas_alone_upper  <-  NA_real_
        alphas <-  priors_checked$prior_dirichlet_alpha
        model_has_fixed_cutpoints <-  !isTRUE(model_obj$advanced_model_options$random_thresholds)
        if (show_what_the_alphas_alone_say && model_has_fixed_cutpoints &&
            is.list(alphas) && length(alphas) == 2) {
              for (row_index in seq_len(nrow(summary_of_the_prior))) {
                    class_index <-  if (summary_of_the_prior$quantity[row_index] == "Se") 2 else 1
                    alpha_of_the_class <-  if (is.matrix(alphas[[class_index]])) {
                          alphas[[class_index]][summary_of_the_prior$test_index[row_index], ]
                    } else {
                          alphas[[class_index]]
                    }
                    k <-  summary_of_the_prior$threshold[row_index]
                    n_categories <-  sum(alpha_of_the_class > 0)
                    alpha_above <-  sum(alpha_of_the_class[(k + 1):n_categories])
                    alpha_at_or_below <-  sum(alpha_of_the_class[1:k])
                    shape_1 <-  if (class_index == 2) alpha_above else alpha_at_or_below
                    shape_2 <-  if (class_index == 2) alpha_at_or_below else alpha_above
                    summary_of_the_prior$alphas_alone_median[row_index] <-  stats::qbeta(0.5, shape_1, shape_2)
                    summary_of_the_prior$alphas_alone_lower[row_index]  <-  stats::qbeta(probs[1], shape_1,
                                                                                         shape_2)
                    summary_of_the_prior$alphas_alone_upper[row_index]  <-  stats::qbeta(probs[2], shape_1,
                                                                                         shape_2)
              }
        }
        ##
        ## ---- The prior at the cutoffs of interest of each test (one row per test and cutoff):
        ##
        rows_at_the_cutoffs <-  list()
        for (test_index in seq_along(test_names)) {
              rows_of_the_test <-  summary_of_the_prior[summary_of_the_prior$test_index == test_index, ]
              for (cutoff in cutoffs_of_interest_of_each_test[[test_index]]) {
                    row_Se <-  rows_of_the_test[rows_of_the_test$quantity == "Se" &
                                                rows_of_the_test$score_at_or_above == cutoff, ]
                    row_Sp <-  rows_of_the_test[rows_of_the_test$quantity == "Sp" &
                                                rows_of_the_test$score_at_or_above == cutoff, ]
                    if (nrow(row_Se) != 1 || nrow(row_Sp) != 1) {
                          stop(paste0("run_prior_predictive_check_of_Se_and_Sp: the cutoff ", cutoff, " of ",
                                      test_names[test_index], " is not one of its thresholds (scores ",
                                      min(rows_of_the_test$score_at_or_above), " to ",
                                      max(rows_of_the_test$score_at_or_above), " or more)."))
                    }
                    rows_at_the_cutoffs[[length(rows_at_the_cutoffs) + 1]] <-  data.frame(
                          test = test_names[test_index],
                          score_at_or_above = cutoff,
                          Se_median = row_Se$prior_median,
                          Se_lower = row_Se$prior_lower,
                          Se_upper = row_Se$prior_upper,
                          Sp_median = row_Sp$prior_median,
                          Sp_lower = row_Sp$prior_lower,
                          Sp_upper = row_Sp$prior_upper,
                          stringsAsFactors = FALSE)
              }
        }
        at_the_cutoffs_of_interest <-  do.call(rbind, rows_at_the_cutoffs)
        at_the_cutoffs_of_interest$test <-  factor(at_the_cutoffs_of_interest$test, levels = test_names)
        fn_median_and_interval <-  function(median_value, lower_value, upper_value) {
              return(paste0(formatC(median_value, format = "f", digits = 2), " [",
                            formatC(lower_value, format = "f", digits = 2), ", ",
                            formatC(upper_value, format = "f", digits = 2), "]"))
        }
        ##
        ## ---- The plot: Se (red) and Sp (green) against the score, the prior median and interval, one panel per
        ##      test; what the alphas alone say as dashed lines; each cutoff of interest as a dotted line, with
        ##      its median and interval written beside its points:
        ##
        alphas_shown <-  any(!is.na(summary_of_the_prior$alphas_alone_median))
        points_at_the_cutoffs <-  rbind(
              data.frame( test = at_the_cutoffs_of_interest$test,
                          quantity = "Se",
                          score_at_or_above = at_the_cutoffs_of_interest$score_at_or_above,
                          prior_median = at_the_cutoffs_of_interest$Se_median,
                          label = paste0("Se ", fn_median_and_interval(at_the_cutoffs_of_interest$Se_median,
                                                                       at_the_cutoffs_of_interest$Se_lower,
                                                                       at_the_cutoffs_of_interest$Se_upper)),
                          vertical_justification = -0.8,
                          stringsAsFactors = FALSE),
              data.frame( test = at_the_cutoffs_of_interest$test,
                          quantity = "Sp",
                          score_at_or_above = at_the_cutoffs_of_interest$score_at_or_above,
                          prior_median = at_the_cutoffs_of_interest$Sp_median,
                          label = paste0("Sp ", fn_median_and_interval(at_the_cutoffs_of_interest$Sp_median,
                                                                       at_the_cutoffs_of_interest$Sp_lower,
                                                                       at_the_cutoffs_of_interest$Sp_upper)),
                          vertical_justification = 1.8,
                          stringsAsFactors = FALSE))
        plot_of_the_prior <-  ggplot2::ggplot( summary_of_the_prior,
                                               ggplot2::aes( x = score_at_or_above,
                                                             colour = quantity,
                                                             fill = quantity)) +
              ggplot2::geom_ribbon(ggplot2::aes(ymin = prior_lower, ymax = prior_upper), alpha = 0.2,
                                   colour = NA) +
              ggplot2::geom_line(ggplot2::aes(y = prior_median), linewidth = 0.9) +
              ggplot2::scale_colour_manual(values = c(Se = "red", Sp = "darkgreen")) +
              ggplot2::scale_fill_manual(values = c(Se = "red", Sp = "darkgreen")) +
              ggplot2::facet_wrap(~ test, scales = "free_x") +
              ggplot2::geom_vline( data = at_the_cutoffs_of_interest,
                                   mapping = ggplot2::aes(xintercept = score_at_or_above),
                                   linetype = "dotted",
                                   colour = "grey40") +
              ggplot2::geom_point( data = points_at_the_cutoffs,
                                   mapping = ggplot2::aes( x = score_at_or_above,
                                                           y = prior_median,
                                                           colour = quantity),
                                   size = 2.5,
                                   inherit.aes = FALSE) +
              ggplot2::geom_text( data = points_at_the_cutoffs,
                                  mapping = ggplot2::aes( x = score_at_or_above,
                                                          y = prior_median,
                                                          label = label,
                                                          colour = quantity,
                                                          vjust = vertical_justification),
                                  hjust = -0.05,
                                  size = 3,
                                  show.legend = FALSE,
                                  inherit.aes = FALSE) +
              ggplot2::labs(x = "Score (a score at or above it counts as positive)", y = "Prior of Se / Sp",
                            title = "Prior predictive check: the model's own prior-only Stan fit",
                            subtitle = paste0("Median and ", 100 * (probs[2] - probs[1]), "% interval",
                                              if (alphas_shown) {
                                                    "; dashed: what the alphas alone say (at the anchor)"
                                              } else {
                                                    ""
                                              },
                                              "; dotted: the cutoffs of interest")) +
              ggplot2::coord_cartesian(ylim = c(0, 1)) +
              ggplot2::theme_bw()
        if (alphas_shown) {
              plot_of_the_prior <-  plot_of_the_prior +
                    ggplot2::geom_line(ggplot2::aes(y = alphas_alone_median), linetype = "dashed",
                                       linewidth = 0.6)
        }
        ##
        min_ESS <-  min(summary_of_the_prior$ESS, na.rm = TRUE)
        name_of_each_row <-  paste0(summary_of_the_prior$quantity, " at threshold ",
                                    summary_of_the_prior$threshold, " of ", summary_of_the_prior$test)
        row_of_the_min_bulk_ESS <-  which.min(summary_of_the_prior$ESS_bulk)
        row_of_the_min_tail_ESS <-  which.min(summary_of_the_prior$ESS_tail)
        row_of_the_min_SD_ESS   <-  which.min(summary_of_the_prior$ESS_SD)
        min_ESS_of_each_type <-  list( min_bulk_ESS = summary_of_the_prior$ESS_bulk[row_of_the_min_bulk_ESS],
                                       min_tail_ESS = summary_of_the_prior$ESS_tail[row_of_the_min_tail_ESS],
                                       min_SD_ESS = summary_of_the_prior$ESS_SD[row_of_the_min_SD_ESS],
                                       quantity_of_the_min_bulk_ESS = name_of_each_row[row_of_the_min_bulk_ESS],
                                       quantity_of_the_min_tail_ESS = name_of_each_row[row_of_the_min_tail_ESS],
                                       quantity_of_the_min_SD_ESS = name_of_each_row[row_of_the_min_SD_ESS])
        message(paste0("\033[32m", "Prior predictive check done in ",
                       formatC(seconds_of_the_fit, format = "f", digits = 1),
                       " s: min ESS over every Se / Sp ", round(min_ESS), " (bulk ",
                       round(min_ESS_of_each_type$min_bulk_ESS), ", tail ",
                       round(min_ESS_of_each_type$min_tail_ESS), ", SD ",
                       round(min_ESS_of_each_type$min_SD_ESS), "), divergent transitions ",
                       n_divergent_transitions, "\033[0m"))
        message(paste0("\033[36m", "The prior of Se and Sp at the cutoffs of interest, median [",
                       100 * probs[1], "%, ", 100 * probs[2], "%]:", "\033[0m"))
        for (row_index in seq_len(nrow(at_the_cutoffs_of_interest))) {
              row_shown <-  at_the_cutoffs_of_interest[row_index, ]
              message(paste0("  ", row_shown$test, ", a score of ", row_shown$score_at_or_above, " or more:  Se ",
                             fn_median_and_interval(row_shown$Se_median, row_shown$Se_lower, row_shown$Se_upper),
                             ";  Sp ",
                             fn_median_and_interval(row_shown$Sp_median, row_shown$Sp_lower, row_shown$Sp_upper)))
        }
        ##
        NicoStan_provenance <-  NULL
        if (algorithm == "NicoStan") {
              NicoStan_provenance <-  list(
                    Stan_model_file_path = outs$NicoStan_model_obj$Stan_model_file_path,
                    Stan_data_list = outs$NicoStan_model_obj$Stan_data_list,
                    model_so_file = outs$NicoStan_model_obj$init_object$model_so_file,
                    efficiency = outs$NicoStan_summary_obj$get_efficiency_metrics())
        }
        ##
        return(list( at_the_cutoffs_of_interest = tibble::as_tibble(at_the_cutoffs_of_interest),
                     summary                  = tibble::as_tibble(summary_of_the_prior),
                     plot                     = plot_of_the_prior,
                     min_ESS                  = min_ESS,
                     min_ESS_of_each_type     = min_ESS_of_each_type,
                     n_divergent_transitions  = n_divergent_transitions,
                     draws_of_Se_and_Sp       = draws_of_each,
                     draws_array              = draws_array,
                     seconds_of_the_fit       = seconds_of_the_fit,
                     NicoStan_provenance      = NicoStan_provenance,
                     priors                   = priors_checked))
##

}
























