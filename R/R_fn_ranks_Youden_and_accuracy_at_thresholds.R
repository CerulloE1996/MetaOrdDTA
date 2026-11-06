


##
## ---- Quantiles (median, 2.5%, 97.5%) of a vector of draws --------------------------------------------------------------------------
##
#' R_fn_median_and_95_interval
#' @keywords internal
R_fn_median_and_95_interval <- function(draws) {

        return(stats::quantile(draws, probs = c(0.5, 0.025, 0.975), names = FALSE))

}



##
## ---- Ranks of the columns of a draws matrix (rank 1 = best), per draw --------------------------------------------------------------
##
#' R_fn_rank_summary
#' @description
#' Ranks of the items (columns of draws) in each draw, rank 1 = the best (the highest value if higher_is_better),
#' and their median, 95% interval, mean and probability of rank 1 over the draws. With the draws of AUC this
#' gives the ranks with 95% CrI; with the draws of AUC_pred, the ranks with 95% PrI.
#' @param draws Matrix of draws (draws x items).
#' @param names_of_the_items Names of the items (columns of draws), e.g. the test names.
#' @param higher_is_better TRUE: the highest value is rank 1.
#' @keywords internal
R_fn_rank_summary <- function(draws,
                              names_of_the_items,
                              higher_is_better) {

        stopifnot(ncol(draws) == length(names_of_the_items))
        sign_of_the_values <- if (higher_is_better) -1 else 1
        ranks <- t(apply(sign_of_the_values * draws, 1, rank, ties.method = "random"))
        rows <- list()
        for (i in seq_along(names_of_the_items)) {

              q <- stats::quantile(ranks[, i], probs = c(0.5, 0.025, 0.975), names = FALSE, type = 1)
              rows[[i]] <- tibble::tibble(item = names_of_the_items[i],
                                          rank_median = q[1],
                                          rank_lower = q[2],
                                          rank_upper = q[3],
                                          rank_mean = mean(ranks[, i]),
                                          prob_rank_1 = mean(ranks[, i] == 1))

        }
        return(dplyr::bind_rows(rows))

}



##
## ---- Se / Sp (median, 95% CrI, 95% PrI) of each test at the chosen thresholds ------------------------------------------------------
##
#' R_fn_accuracy_at_thresholds
#' @description
#' Summary Se and Sp of each test at the chosen thresholds, with 95% CrI (from the draws of Se_baseline /
#' Sp_baseline) and 95% PrI (from the draws of Se_baseline_pred / Sp_baseline_pred), for the NMA models
#' (columns named "<family>[t,k]").
#' @param thresholds_list Named list: for each test name, the thresholds k (Se = P(score >= k)).
#' @keywords internal
R_fn_accuracy_at_thresholds <- function(Se_draws,
                                        Sp_draws,
                                        Se_pred_draws,
                                        Sp_pred_draws,
                                        thresholds_list,
                                        test_names) {

        rows <- list()
        for (test_name in names(thresholds_list)) {

              t <- match(test_name, test_names)
              stopifnot(!is.na(t))
              for (k in thresholds_list[[test_name]]) {

                    name_of_the_column <- paste0("[", t, ",", k, "]")
                    q_Se <- R_fn_median_and_95_interval(Se_draws[, paste0("Se_baseline", name_of_the_column)])
                    q_Sp <- R_fn_median_and_95_interval(Sp_draws[, paste0("Sp_baseline", name_of_the_column)])
                    name_Se_pred <- paste0("Se_baseline_pred", name_of_the_column)
                    name_Sp_pred <- paste0("Sp_baseline_pred", name_of_the_column)
                    q_Se_pred <- R_fn_median_and_95_interval(Se_pred_draws[, name_Se_pred])
                    q_Sp_pred <- R_fn_median_and_95_interval(Sp_pred_draws[, name_Sp_pred])
                    rows[[length(rows) + 1]] <- tibble::tibble(test = test_name,
                                                               threshold = k,
                                                               Se_median = q_Se[1],
                                                               Se_lower = q_Se[2],
                                                               Se_upper = q_Se[3],
                                                               Se_pred_lower = q_Se_pred[2],
                                                               Se_pred_upper = q_Se_pred[3],
                                                               Sp_median = q_Sp[1],
                                                               Sp_lower = q_Sp[2],
                                                               Sp_upper = q_Sp[3],
                                                               Sp_pred_lower = q_Sp_pred[2],
                                                               Sp_pred_upper = q_Sp_pred[3])

              }

        }
        return(dplyr::bind_rows(rows))

}



##
## ---- Threshold of each test that maximises the Youden index (Se + Sp - 1) ----------------------------------------------------------
##
#' R_fn_Youden_optimal_thresholds
#' @description
#' For each test, the threshold k with the highest posterior mean of the Youden index J = Se + Sp - 1 (from the
#' draws of Se_baseline / Sp_baseline), the probability that k has the highest J (share of draws), and J at k.
#' @keywords internal
R_fn_Youden_optimal_thresholds <- function(Se_draws,
                                           Sp_draws,
                                           n_thr,
                                           test_names) {

        rows <- list()
        for (t in seq_along(test_names)) {

              names_Se <- paste0("Se_baseline[", t, ",", 1:n_thr[t], "]")
              names_Sp <- paste0("Sp_baseline[", t, ",", 1:n_thr[t], "]")
              J <- Se_draws[, names_Se, drop = FALSE] + Sp_draws[, names_Sp, drop = FALSE] - 1.0
              k_best <- which.max(colMeans(J))
              q_J <- R_fn_median_and_95_interval(J[, k_best])
              rows[[t]] <- tibble::tibble(test = test_names[t],
                                          threshold = k_best,
                                          prob_highest_J = mean(apply(J, 1, which.max) == k_best),
                                          J_median = q_J[1],
                                          J_lower = q_J[2],
                                          J_upper = q_J[3])

        }
        return(dplyr::bind_rows(rows))

}



##
## ---- Differences in Se / Sp between each pair of tests, each at its chosen threshold -----------------------------------------------
##
#' R_fn_pairwise_accuracy_differences
#' @description
#' Se(test 1) - Se(test 2) and Sp(test 1) - Sp(test 2), each test at its threshold in thresholds_of_each_test,
#' median and 95% CrI over the draws of Se_baseline / Sp_baseline, for every pair of tests.
#' @param thresholds_of_each_test Named vector: one threshold k per test name.
#' @keywords internal
R_fn_pairwise_accuracy_differences <- function(Se_draws,
                                               Sp_draws,
                                               thresholds_of_each_test,
                                               test_names) {

        names_of_the_tests <- names(thresholds_of_each_test)
        rows <- list()
        for (i in 1:(length(names_of_the_tests) - 1)) {

              for (j in (i + 1):length(names_of_the_tests)) {

                    t1 <- match(names_of_the_tests[i], test_names)
                    t2 <- match(names_of_the_tests[j], test_names)
                    k1 <- thresholds_of_each_test[[i]]
                    k2 <- thresholds_of_each_test[[j]]
                    diff_Se <- Se_draws[, paste0("Se_baseline[", t1, ",", k1, "]")] -
                               Se_draws[, paste0("Se_baseline[", t2, ",", k2, "]")]
                    diff_Sp <- Sp_draws[, paste0("Sp_baseline[", t1, ",", k1, "]")] -
                               Sp_draws[, paste0("Sp_baseline[", t2, ",", k2, "]")]
                    q_Se <- R_fn_median_and_95_interval(diff_Se)
                    q_Sp <- R_fn_median_and_95_interval(diff_Sp)
                    rows[[length(rows) + 1]] <- tibble::tibble(test1 = names_of_the_tests[i],
                                                               threshold1 = k1,
                                                               test2 = names_of_the_tests[j],
                                                               threshold2 = k2,
                                                               diff_Se_median = q_Se[1],
                                                               diff_Se_lower = q_Se[2],
                                                               diff_Se_upper = q_Se[3],
                                                               diff_Sp_median = q_Sp[1],
                                                               diff_Sp_lower = q_Sp[2],
                                                               diff_Sp_upper = q_Sp[3])

              }

        }
        return(dplyr::bind_rows(rows))

}



##
## ---- LaTeX: a value with its 95% interval(s) ---------------------------------------------------------------------------------------
##
#' R_fn_latex_value_and_intervals
#' @keywords internal
R_fn_latex_value_and_intervals <- function(median,
                                           lower,
                                           upper,
                                           pred_lower,
                                           pred_upper,
                                           digits) {

        format_one <- function(x) formatC(round(x, digits) + 0, format = "f", digits = digits)
        text <- paste0(format_one(median), " (", format_one(lower), ", ", format_one(upper), ")")
        if (!is.null(pred_lower)) {
          text <- paste0(text, " [", format_one(pred_lower), ", ", format_one(pred_upper), "]")
        }
        return(text)

}



##
## ---- LaTeX: AUC and the rank of each test (95% CrI and 95% PrI) --------------------------------------------------------------------
##
#' create_AUC_rank_latex_table
#' @param auc_outs Output of R_fn_extract_NMA_AUC_summary (auc, auc_pred).
#' @param ranks_CrI Output of R_fn_rank_summary on the draws of AUC.
#' @param ranks_PrI Output of R_fn_rank_summary on the draws of AUC_pred.
#' @keywords internal
create_AUC_rank_latex_table <- function(auc_outs,
                                        ranks_CrI,
                                        ranks_PrI,
                                        caption,
                                        label) {

        order_of_the_rows <- order(ranks_CrI$rank_mean)
        lines <- c("\\begin{table}[H]",
                   "\\centering",
                   paste0("\\caption{", caption, "}\\label{", label, "}"),
                   "\\small",
                   "\\begin{tabular}{lcccc}",
                   "\\toprule",
                   paste0("\\textbf{Test} & \\textbf{AUC (\\%) (95\\% CrI) [95\\% PrI]} & ",
                          "\\textbf{Rank (95\\% CrI)} & \\textbf{Rank (95\\% PrI)} & \\textbf{P(rank 1)} \\\\"),
                   "\\midrule")
        for (i in order_of_the_rows) {

              t <- match(ranks_CrI$item[i], auc_outs$auc$test_name)
              lines <- c(lines,
                         sprintf("%s & %.1f (%.1f, %.1f) [%.1f, %.1f] & %g (%g, %g) & %g (%g, %g) & %.2f \\\\",
                                 gsub("_", "\\\\_", ranks_CrI$item[i]),
                                 auc_outs$auc$AUC_median[t], auc_outs$auc$AUC_lower[t], auc_outs$auc$AUC_upper[t],
                                 auc_outs$auc_pred$AUC_lower[t], auc_outs$auc_pred$AUC_upper[t],
                                 ranks_CrI$rank_median[i], ranks_CrI$rank_lower[i], ranks_CrI$rank_upper[i],
                                 ranks_PrI$rank_median[i], ranks_PrI$rank_lower[i], ranks_PrI$rank_upper[i],
                                 ranks_CrI$prob_rank_1[i]))

        }
        lines <- c(lines,
                   "\\bottomrule",
                   "\\end{tabular}",
                   "\\begin{tablenotes}",
                   "\\footnotesize",
                   paste0("\\item Rank 1 = the highest AUC. Rank CrI from the summary AUCs; ",
                          "rank PrI from the predicted AUCs of a new study."),
                   "\\item P(rank 1) = posterior probability that the test has the highest summary AUC.",
                   "\\end{tablenotes}",
                   "\\end{table}")
        return(lines)

}



##
## ---- LaTeX: Se / Sp with 95% CrI and 95% PrI at the chosen thresholds --------------------------------------------------------------
##
#' create_accuracy_CrI_PrI_latex_table
#' @param accuracy Output of R_fn_accuracy_at_thresholds (optionally with a column "note" per row).
#' @keywords internal
create_accuracy_CrI_PrI_latex_table <- function(accuracy,
                                                caption,
                                                label,
                                                digits) {

        lines <- c("\\begin{table}[H]",
                   "\\centering",
                   paste0("\\caption{", caption, "}\\label{", label, "}"),
                   "\\small",
                   "\\begin{tabular}{llcc}",
                   "\\toprule",
                   paste0("\\textbf{Test} & \\textbf{Threshold} & ",
                          "\\textbf{Sensitivity (95\\% CrI) [95\\% PrI]} & ",
                          "\\textbf{Specificity (95\\% CrI) [95\\% PrI]} \\\\"),
                   "\\midrule")
        for (i in seq_len(nrow(accuracy))) {

              note <- if ("note" %in% names(accuracy)) accuracy$note[i] else ""
              row <- accuracy[i, ]
              text_Se <- R_fn_latex_value_and_intervals(row$Se_median, row$Se_lower, row$Se_upper,
                                                        row$Se_pred_lower, row$Se_pred_upper, digits)
              text_Sp <- R_fn_latex_value_and_intervals(row$Sp_median, row$Sp_lower, row$Sp_upper,
                                                        row$Sp_pred_lower, row$Sp_pred_upper, digits)
              lines <- c(lines,
                         paste0(gsub("_", "\\\\_", row$test), " & $\\geq ", row$threshold, "$", note, " & ",
                                text_Se, " & ", text_Sp, " \\\\"))

        }
        lines <- c(lines,
                   "\\bottomrule",
                   "\\end{tabular}",
                   "\\begin{tablenotes}",
                   "\\footnotesize",
                   "\\item CrI = credible interval (summary estimate); PrI = prediction interval (a new study).",
                   "\\end{tablenotes}",
                   "\\end{table}")
        return(lines)

}



##
## ---- LaTeX: ranks of Se and of Sp across the tests, each test at its chosen threshold ----------------------------------------------
##
#' create_Se_Sp_rank_latex_table
#' @param ranks_Se,ranks_Sp Outputs of R_fn_rank_summary on the draws of Se / Sp (one column per test).
#' @param thresholds_of_each_test Named vector: the threshold of each test.
#' @keywords internal
create_Se_Sp_rank_latex_table <- function(ranks_Se,
                                          ranks_Sp,
                                          thresholds_of_each_test,
                                          caption,
                                          label) {

        lines <- c("\\begin{table}[H]",
                   "\\centering",
                   paste0("\\caption{", caption, "}\\label{", label, "}"),
                   "\\small",
                   "\\begin{tabular}{llcccc}",
                   "\\toprule",
                   paste0("\\textbf{Test} & \\textbf{Threshold} & ",
                          "\\textbf{Se rank (95\\% CrI)} & \\textbf{P(Se rank 1)} & ",
                          "\\textbf{Sp rank (95\\% CrI)} & \\textbf{P(Sp rank 1)} \\\\"),
                   "\\midrule")
        for (i in seq_len(nrow(ranks_Se))) {

              lines <- c(lines,
                         sprintf("%s & $\\geq %g$ & %g (%g, %g) & %.2f & %g (%g, %g) & %.2f \\\\",
                                 gsub("_", "\\\\_", ranks_Se$item[i]),
                                 thresholds_of_each_test[[ranks_Se$item[i]]],
                                 ranks_Se$rank_median[i], ranks_Se$rank_lower[i], ranks_Se$rank_upper[i],
                                 ranks_Se$prob_rank_1[i],
                                 ranks_Sp$rank_median[i], ranks_Sp$rank_lower[i], ranks_Sp$rank_upper[i],
                                 ranks_Sp$prob_rank_1[i]))

        }
        lines <- c(lines,
                   "\\bottomrule",
                   "\\end{tabular}",
                   "\\begin{tablenotes}",
                   "\\footnotesize",
                   paste0("\\item Rank 1 = the highest summary sensitivity (specificity) among the tests, ",
                          "each at its threshold."),
                   "\\end{tablenotes}",
                   "\\end{table}")
        return(lines)

}



##
## ---- LaTeX: pairwise differences in Se and Sp --------------------------------------------------------------------------------------
##
#' create_pairwise_accuracy_differences_latex_table
#' @param differences Output of R_fn_pairwise_accuracy_differences.
#' @keywords internal
create_pairwise_accuracy_differences_latex_table <- function(differences,
                                                             caption,
                                                             label,
                                                             digits) {

        lines <- c("\\begin{table}[H]",
                   "\\centering",
                   paste0("\\caption{", caption, "}\\label{", label, "}"),
                   "\\small",
                   "\\begin{tabular}{lcc}",
                   "\\toprule",
                   paste0("\\textbf{Comparison} & \\textbf{Difference in Se (95\\% CrI)} & ",
                          "\\textbf{Difference in Sp (95\\% CrI)} \\\\"),
                   "\\midrule")
        for (i in seq_len(nrow(differences))) {

              star_Se <- if (differences$diff_Se_lower[i] > 0 || differences$diff_Se_upper[i] < 0) "*" else ""
              star_Sp <- if (differences$diff_Sp_lower[i] > 0 || differences$diff_Sp_upper[i] < 0) "*" else ""
              row <- differences[i, ]
              text_Se <- R_fn_latex_value_and_intervals(row$diff_Se_median, row$diff_Se_lower, row$diff_Se_upper,
                                                        NULL, NULL, digits)
              text_Sp <- R_fn_latex_value_and_intervals(row$diff_Sp_median, row$diff_Sp_lower, row$diff_Sp_upper,
                                                        NULL, NULL, digits)
              lines <- c(lines,
                         paste0(gsub("_", "\\\\_", row$test1), " ($\\geq ", row$threshold1, "$) vs ",
                                gsub("_", "\\\\_", row$test2), " ($\\geq ", row$threshold2, "$) & ",
                                text_Se, star_Se, " & ", text_Sp, star_Sp, " \\\\"))

        }
        lines <- c(lines,
                   "\\bottomrule",
                   "\\end{tabular}",
                   "\\begin{tablenotes}",
                   "\\footnotesize",
                   "\\item Difference = test 1 minus test 2; * = 95\\% CrI excludes zero.",
                   "\\end{tablenotes}",
                   "\\end{table}")
        return(lines)

}



##
## ---- Between-study SDs of meta-regressions and their % reduction from the model without covariates (CS models) ---------------------
##
#' R_fn_heterogeneity_reduction
#' @description
#' Posterior medians of the between-study SDs (beta_sigma, and beta_tau and total_SD_inc_C of test 1, the same
#' for every test with compound symmetry) of the model without covariates and of each meta-regression, and the
#' percentage reduction 100 * (1 - SD of the meta-regression / SD without covariates).
#' @param Intercept_only_model_summary_and_trace_obj The model without covariates.
#' @param list_of_MR_model_summary_and_trace_objs Named list of the meta-regressions (names = their covariates).
#' @keywords internal
R_fn_heterogeneity_reduction <- function(Intercept_only_model_summary_and_trace_obj,
                                         list_of_MR_model_summary_and_trace_objs) {

        names_of_the_SDs <- c("beta_sigma[1]", "beta_sigma[2]", "beta_tau[1,1]", "beta_tau[1,2]",
                              "total_SD_inc_C[1,1]", "total_SD_inc_C[1,2]")
        fn_medians <- function(model_summary_and_trace_obj) {

              families <- c("beta_sigma", "beta_tau", "total_SD_inc_C")
              summary_of_the_SDs <- dplyr::bind_rows(lapply(families, function(family) {

                    model_summary_and_trace_obj$extract_params(params = c(family))

              }))
              return(summary_of_the_SDs$`50%`[match(names_of_the_SDs, summary_of_the_SDs$parameter)])

        }
        SD_without_covariates <- fn_medians(Intercept_only_model_summary_and_trace_obj)
        output <- tibble::tibble(parameter = names_of_the_SDs,
                                 without_covariates = SD_without_covariates)
        for (name in names(list_of_MR_model_summary_and_trace_objs)) {

              SD_of_the_MR <- fn_medians(list_of_MR_model_summary_and_trace_objs[[name]])
              output[[paste0(name, ": SD")]] <- SD_of_the_MR
              output[[paste0(name, ": % reduction")]] <- 100 * (1 - SD_of_the_MR / SD_without_covariates)

        }
        return(output)

}



##
## ---- LaTeX: between-study SDs and their % reduction by each meta-regression --------------------------------------------------------
##
#' create_heterogeneity_reduction_latex_table
#' @param heterogeneity_reduction Output of R_fn_heterogeneity_reduction.
#' @keywords internal
create_heterogeneity_reduction_latex_table <- function(heterogeneity_reduction,
                                                       caption,
                                                       label) {

        names_of_the_MRs <- sub(": SD$", "", grep(": SD$", names(heterogeneity_reduction), value = TRUE))
        labels_of_the_SDs <- c("$\\sigma_{\\beta}^{(d-)}$", "$\\sigma_{\\beta}^{(d+)}$",
                               "$\\tau^{(d-)}$", "$\\tau^{(d+)}$",
                               "$\\sigma_{total}^{(d-)}$", "$\\sigma_{total}^{(d+)}$")
        lines <- c("\\begin{table}[H]",
                   "\\centering",
                   paste0("\\caption{", caption, "}\\label{", label, "}"),
                   "\\small",
                   "\\resizebox{\\textwidth}{!}{%",
                   paste0("\\begin{tabular}{l", strrep("c", 1 + length(names_of_the_MRs)), "}"),
                   "\\toprule",
                   paste0("\\textbf{SD} & \\textbf{No covariates} & ",
                          paste0("\\textbf{", names_of_the_MRs, "}", collapse = " & "), " \\\\"),
                   "\\midrule")
        for (i in seq_len(nrow(heterogeneity_reduction))) {

              cells <- sapply(names_of_the_MRs, function(name) {

                    sprintf("%.3f (%+.0f\\%%)", heterogeneity_reduction[[paste0(name, ": SD")]][i],
                            round(-heterogeneity_reduction[[paste0(name, ": % reduction")]][i]) + 0)

              })
              lines <- c(lines,
                         paste0(labels_of_the_SDs[i], " & ",
                                sprintf("%.3f", heterogeneity_reduction$without_covariates[i]),
                                " & ", paste(cells, collapse = " & "), " \\\\"))

        }
        lines <- c(lines,
                   "\\bottomrule",
                   "\\end{tabular}",
                   "}",
                   "\\begin{tablenotes}",
                   "\\footnotesize",
                   paste0("\\item Posterior medians; in brackets the change from the model without covariates ",
                          "(negative = reduction)."),
                   "\\item Compound symmetry: $\\tau$ and $\\sigma_{total}$ are the same for every test.",
                   "\\end{tablenotes}",
                   "\\end{table}")
        return(lines)

}



##
## ---- AUC and Se / Sp of two fits of the same NMA model (e.g. all studies and a sensitivity analysis), side by side -----------------
##
#' R_fn_compare_accuracy_of_two_fits
#' @description
#' Median and 95% CrI of AUC and of Se / Sp at the chosen thresholds, from the draws of two fits (lists with the
#' draws matrices AUC, Se_baseline and Sp_baseline, as from R_fn_draws_of_one_parameter_family).
#' @param thresholds_list Named list: for each test name, the thresholds k (Se = P(score >= k)).
#' @keywords internal
R_fn_compare_accuracy_of_two_fits <- function(draws_1,
                                              draws_2,
                                              thresholds_list,
                                              test_names) {

        fn_q <- function(draws, name) R_fn_median_and_95_interval(draws[, name])
        rows <- list()
        for (test_name in names(thresholds_list)) {

              t <- match(test_name, test_names)
              stopifnot(!is.na(t))
              for (k in thresholds_list[[test_name]]) {

                    name_AUC <- paste0("AUC[", t, "]")
                    name_Se <- paste0("Se_baseline[", t, ",", k, "]")
                    name_Sp <- paste0("Sp_baseline[", t, ",", k, "]")
                    rows[[length(rows) + 1]] <- tibble::tibble(test = test_name,
                                                               threshold = k,
                                                               AUC_1 = list(fn_q(draws_1$AUC, name_AUC)),
                                                               AUC_2 = list(fn_q(draws_2$AUC, name_AUC)),
                                                               Se_1 = list(fn_q(draws_1$Se_baseline, name_Se)),
                                                               Se_2 = list(fn_q(draws_2$Se_baseline, name_Se)),
                                                               Sp_1 = list(fn_q(draws_1$Sp_baseline, name_Sp)),
                                                               Sp_2 = list(fn_q(draws_2$Sp_baseline, name_Sp)))

              }

        }
        return(dplyr::bind_rows(rows))

}



##
## ---- LaTeX: AUC and Se / Sp of two fits side by side -------------------------------------------------------------------------------
##
#' create_two_fit_comparison_latex_table
#' @param comparison Output of R_fn_compare_accuracy_of_two_fits.
#' @param name_1,name_2 Short names of the two fits for the column headers.
#' @keywords internal
create_two_fit_comparison_latex_table <- function(comparison,
                                                  name_1,
                                                  name_2,
                                                  caption,
                                                  label,
                                                  digits) {

        fn_text <- function(q) R_fn_latex_value_and_intervals(q[1], q[2], q[3], NULL, NULL, digits)
        fn_text_AUC <- function(q) sprintf("%.1f (%.1f, %.1f)", 100 * q[1], 100 * q[2], 100 * q[3])
        lines <- c("\\begin{table}[H]",
                   "\\centering",
                   paste0("\\caption{", caption, "}\\label{", label, "}"),
                   "\\small",
                   "\\resizebox{\\textwidth}{!}{%",
                   "\\begin{tabular}{llcccccc}",
                   "\\toprule",
                   paste0(" & & \\multicolumn{2}{c}{\\textbf{AUC (\\%)}} & ",
                          "\\multicolumn{2}{c}{\\textbf{Sensitivity}} & ",
                          "\\multicolumn{2}{c}{\\textbf{Specificity}} \\\\"),
                   paste0("\\textbf{Test} & \\textbf{Threshold} & ", paste(rep(paste0(name_1, " & ", name_2), 3),
                                                                             collapse = " & "), " \\\\"),
                   "\\midrule")
        for (i in seq_len(nrow(comparison))) {

              row <- comparison[i, ]
              first_row_of_the_test <- i == 1 || comparison$test[i - 1] != row$test
              text_AUC <- if (first_row_of_the_test) {
                paste0(fn_text_AUC(row$AUC_1[[1]]), " & ", fn_text_AUC(row$AUC_2[[1]]))
              } else {
                " & "
              }
              lines <- c(lines,
                         paste0(if (first_row_of_the_test) gsub("_", "\\\\_", row$test) else "", " & $\\geq ",
                                row$threshold, "$ & ", text_AUC, " & ",
                                fn_text(row$Se_1[[1]]), " & ", fn_text(row$Se_2[[1]]), " & ",
                                fn_text(row$Sp_1[[1]]), " & ", fn_text(row$Sp_2[[1]]), " \\\\"))

        }
        lines <- c(lines,
                   "\\bottomrule",
                   "\\end{tabular}",
                   "}",
                   "\\begin{tablenotes}",
                   "\\footnotesize",
                   "\\item Median (95\\% CrI).",
                   "\\end{tablenotes}",
                   "\\end{table}")
        return(lines)

}
