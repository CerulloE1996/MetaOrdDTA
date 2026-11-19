##
## ---- QC functions: every results table of the simulation study, made again straight from the scored outputs ---------------------
##
## What this file does:
##   One function per results table of the simulation study (the overall ranking; the ranking by share of thresholds missing;
##   the ranking by test; the RMSE in the nine sets; coverage and width; RMSE and bias by prior set; the models by DGM; the
##   appendix rankings of prior sets 0, 1 and 2; the detailed appendix tables of each test and share missing). Each function takes
##   the folder of the scored outputs (the folder holding table_of_every_cell_of_the_nine_sets.rds) and returns a data frame of
##   the numbers and symbols of the table, one row per cell (table, row, column, number, text, marks). It is written from the
##   definitions in the notes of the tables of the paper, NOT from the scripts that made the tables (it is an independent check).
##   The second part reads the tables of a LaTeX file (found by their \label, or the whole file), takes out their numbers and
##   symbols, and compares them cell by cell with the data frames; QC_compare_every_table_with_the_paper() writes a report.
##
## The data it reads: only the scored outputs of the simulation study (table_of_every_cell_of_the_nine_sets.rds,
## table_of_every_cell_of_the_redo.rds, fits_kept_and_lost_of_every_cell_of_the_redo.csv) and the LaTeX of the tables.
##
## How to run it:
##   Rscript QC_functions_that_make_every_table_from_the_outputs_19_Nov.R <outputs folder> <paper .tex file> <report .md file> [old|new]
##   (the last argument is the grouping of the DGMs of the .tex being checked: "new" (the default) = two families, ordinal and
##   continuous, the grouping of the paper since 19 Nov 2026 23:21; "old" = the Jones / ordinal / competing DGMs of the earlier versions)
##   or, in R: source() this file, then e.g. QC_table_overall_ranking("<outputs folder>").
##
## Needs base R only (utils, stats). Written 19 Nov 2026.
##
## Definitions (from the notes of the paper):
##   RMSE of a cell = RMSE(Se) + RMSE(Sp), bias of a cell = |bias(Se)| + |bias(Sp)| (percentage points);
##   MCSE of a cell = the larger of the MCSEs of Se and Sp (for RMSE and for bias, each);
##   4-group classification of a cell within its scenario, relative to the best (lowest) model of the scenario:
##     statistically different = difference >= 4 x MCSE (taken here as 4 x the larger of the MCSEs of the cell and of the best cell:
##     the paper does not say which MCSE; see the report), practically different = at least 10% higher;
##     Group 1 = neither, 2 = statistically only, 3 = practically only, 4 = both;
##   main-text means: MCSE of a mean = sqrt(sum of MCSE_i^2) / n over its n scenarios; star = Group 4, dagger = Group 2, vs the
##     best-ranked model of the column, over the scenarios in which every model of the column has a result (like for like);
##   appendix rankings (notes of those tables): star = a difference of at least 4 x MCSE, the MCSE of a mean = the mean MCSE
##     across the scenarios, no practical criterion; double dagger = a model with no result in some scenarios, ranked and compared
##     with the best model against the means of the other models across the same scenarios.
##   Prior sets 1 and 2 are for the five Bayesian models only: the frequentist models have no priors (their rows are left out by
##     default; frequentist_models = "included_with_their_set_0_cells" reproduces what the appendix tables of 19 Nov did).
##
QC_relative_difference_of_practical_significance <-  10
QC_paper_name_of_the_model <-  c("strat. bivariate" = "Strat-biv", "Jones (matching priors)" = "Jones (matching priors)", "Jones (own priors)" = "Jones (own priors)",
                                 "O-biv-FC" = "O-biv-FC", "O-biv-RC" = "O-biv-RC", "ord-HSROC" = "O-HSROC-RC",
                                 "Hoyer 2018" = "Hoyer", "logit GLMM" = "Logit-GLMM", "diagmeta" = "Steinhauser (DIDS)", "diagmeta (CICS)" = "Steinhauser (CICS)")
QC_frequentist_models <-  c("Hoyer 2018", "logit GLMM", "diagmeta", "diagmeta (CICS)")
QC_paper_name_of_the_DGM <-  c("Jones" = "Jones", "ord_biv_fixed_C" = "O-biv-FC", "ord_biv_random_C" = "O-biv-RC", "R&G_random_C" = "O-HSROC-RC",
                               "Hoyer_2018" = "Hoyer", "logit_GLMM" = "Logit-GLMM", "diagmeta" = "Steinhauser")
## (the model that generated the data of each DGM: its row is in bold in the detailed tables; both Jones models under the Jones DGM)
QC_model_of_the_DGM <-  c("Jones" = "Jones (matching priors)", "ord_biv_fixed_C" = "O-biv-FC", "ord_biv_random_C" = "O-biv-RC", "R&G_random_C" = "ord-HSROC",
                          "Hoyer_2018" = "Hoyer 2018", "logit_GLMM" = "logit GLMM", "diagmeta" = "diagmeta")
QC_paper_name_of_the_test <-  c(GAD_2 = "GAD-2", HADS = "HADS", BAI = "BAI")
QC_tests <-  c("GAD_2", "HADS", "BAI")
QC_shares_missing <-  c(15, 40, 70)
QC_DGMs <-  c("Jones", "ord_biv_fixed_C", "ord_biv_random_C", "R&G_random_C", "Hoyer_2018", "logit_GLMM", "diagmeta")
##
## ---- The family of each DGM. "old" = the grouping of the paper of 19 Nov 2026 (the Jones DGM; the three ordinal DGMs; the three
## ---- competing DGMs); "new" (19 Nov 2026) = two families: ordinal (O-biv-FC, O-biv-RC, O-HSROC-RC) and continuous
## ---- (Jones, Hoyer, Logit-GLMM, Steinhauser).
##
QC_family_of_each_DGM <- function(DGM,
                                  family_grouping = "new") {

      ordinal <-  DGM %in% c("ord_biv_fixed_C", "ord_biv_random_C", "R&G_random_C")
      if (family_grouping == "old") return(ifelse(DGM == "Jones", "Jones", ifelse(ordinal, "Ordinal", "Competing")))
      if (family_grouping == "new") return(ifelse(ordinal, "Ordinal", "Continuous"))
      stop("family_grouping must be \"old\" or \"new\"")

}
QC_families <- function(family_grouping = "new") {
      if (family_grouping == "old") c("Jones", "Ordinal", "Competing") else c("Ordinal", "Continuous")
}
##
## ---- The scored cells of the paper (the nine sets, the 14 failed cells of diagmeta DIDS left out), with the measures of the paper:
##
QC_read_the_cells_of_the_paper <- function(outputs_folder,
                                           family_grouping = "new") {

      cells <-  readRDS(file.path(outputs_folder, "table_of_every_cell_of_the_nine_sets.rds"))
      cells$scenario <-  paste(cells$test, cells$pct, cells$DGM, cells$n_studies)
      cells$RMSE <-  cells$avg_RMSE_Se + cells$avg_RMSE_Sp
      cells$bias <-  cells$avg_abs_bias_Se + cells$avg_abs_bias_Sp
      cells$MCSE_RMSE <-  pmax(cells$avg_MCSE_RMSE_Se, cells$avg_MCSE_RMSE_Sp)
      cells$MCSE_bias <-  pmax(cells$avg_MCSE_bias_Se, cells$avg_MCSE_bias_Sp)
      cells$paper_model <-  unname(QC_paper_name_of_the_model[cells$model])
      cells$frequentist <-  cells$model %in% QC_frequentist_models
      cells$family <-  QC_family_of_each_DGM(cells$DGM, family_grouping)
      stopifnot(!any(is.na(cells$paper_model)), !any(is.na(cells$RMSE)), !any(is.na(cells$MCSE_RMSE)), !any(is.na(cells$MCSE_bias)))
      cells

}
##
## ---- The cells left out of the paper: every cell of the redo that is not in the table of the nine sets, with its share of fits not
## ---- kept, and whether the table of the nine sets is exactly the table of the redo less them.
##
QC_cells_left_out_of_the_paper <- function(outputs_folder) {

      every_cell <-  readRDS(file.path(outputs_folder, "table_of_every_cell_of_the_redo.rds"))
      cells_of_the_paper <-  readRDS(file.path(outputs_folder, "table_of_every_cell_of_the_nine_sets.rds"))
      key <- function(x) paste(x$model, x$prior_set, x$prior_set_of_the_run, x$test, x$pct, x$DGM, x$n_studies)
      left_out <-  every_cell[!(key(every_cell) %in% key(cells_of_the_paper)), c("model", "prior_set", "test", "pct", "DGM", "n_studies", "kept_fits", "fits_tried", "percent_of_fits_not_kept")]
      kept <-  every_cell[key(every_cell) %in% key(cells_of_the_paper), ]
      attr(left_out, "rest_identical_to_the_redo") <-  isTRUE(all.equal(kept, cells_of_the_paper, check.attributes = FALSE))
      attr(left_out, "cells_with_over_50_pct_not_kept") <-  sum(every_cell$percent_of_fits_not_kept > 50, na.rm = TRUE)
      left_out

}
##
## ---- The cells of one prior set. Prior set 0 = prior_set "set_0" (the Jones model with its own priors, run as set_4, is a model of
## ---- set 0 only). Prior sets 1 and 2: the five Bayesian models; the frequentist models are left out unless asked for.
##
QC_cells_of_one_prior_set <- function(cells,
                                      prior_set,
                                      frequentist_models = "left_out") {

      if (prior_set == "set_0") return(cells[cells$prior_set == "set_0", ])
      chosen <-  cells[cells$prior_set == prior_set, ]
      if (frequentist_models == "included_with_their_set_0_cells") chosen <-  rbind(chosen, cells[cells$frequentist & cells$prior_set == "set_0", ])
      chosen

}
##
## ---- The 4-group classification of every cell within its scenario (1 = Best, 2 = Stat only, 3 = Practical only, 4 = Worse):
##
QC_group_of_each_cell_in_its_scenario <- function(value,
                                                  MCSE,
                                                  scenario) {

      group <-  rep(NA_integer_, length(value))
      for (one_scenario in unique(scenario)) {
            i <-  which(scenario == one_scenario)
            best <-  i[which.min(value[i])]
            difference <-  value[i] - value[best]
            statistically <-  difference >= 4 * pmax(MCSE[i], MCSE[best])
            practically <-  100 * difference / value[best] >= QC_relative_difference_of_practical_significance
            group[i] <-  ifelse(!statistically & !practically, 1L, ifelse(statistically & !practically, 2L, ifelse(!statistically & practically, 3L, 4L)))
      }
      group

}
QC_add_the_groups <- function(cells) {
      cells$group_RMSE <-  QC_group_of_each_cell_in_its_scenario(cells$RMSE, cells$MCSE_RMSE, cells$scenario)
      cells$group_bias <-  QC_group_of_each_cell_in_its_scenario(cells$bias, cells$MCSE_bias, cells$scenario)
      cells
}
##
## ---- The MCSE of a mean across scenarios:
##
QC_MCSE_of_a_mean <- function(MCSE_of_each_scenario,
                              rule_of_the_MCSE_of_a_mean) {
      if (rule_of_the_MCSE_of_a_mean == "root_of_the_sum_of_squares_over_n") return(sqrt(sum(MCSE_of_each_scenario^2)) / length(MCSE_of_each_scenario))
      if (rule_of_the_MCSE_of_a_mean == "mean_of_the_MCSEs") return(mean(MCSE_of_each_scenario))
      stop("unknown rule of the MCSE of a mean")
}
##
## ---- The ranking of the models of one column by the mean of one measure ("RMSE" or "bias"). The cells must carry their groups.
## ---- like_for_like = TRUE: only the scenarios in which every model of the column has a result (main text). FALSE: a model with
## ---- fewer scenarios gets a double dagger and is placed and judged against the means of the other models over its scenarios.
##
QC_rank_the_means_of_one_column <- function(cells,
                                            measure,
                                            like_for_like,
                                            rule_of_the_MCSE_of_a_mean,
                                            practical_criterion) {

      other_measure <-  if (measure == "RMSE") "bias" else "RMSE"
      if (like_for_like) {
            common_scenarios <-  Reduce(intersect, split(cells$scenario, cells$model))
            cells <-  cells[cells$scenario %in% common_scenarios, ]
      }
      scenarios_of_the_column <-  length(unique(cells$scenario))
      models <-  unique(cells$model)
      ranking <-  data.frame(model = models, stringsAsFactors = FALSE)
      ranking$n_scenarios <-  sapply(models, function(m) sum(cells$model == m))
      ranking$mean <-  sapply(models, function(m) mean(cells[[measure]][cells$model == m]))
      ranking$MCSE_of_the_mean <-  sapply(models, function(m) QC_MCSE_of_a_mean(cells[[paste0("MCSE_", measure)]][cells$model == m], rule_of_the_MCSE_of_a_mean))
      ranking$mean_of_the_other_measure <-  sapply(models, function(m) mean(cells[[other_measure]][cells$model == m]))
      ranking$cells_in_group_1 <-  sapply(models, function(m) sum(cells[[paste0("group_", measure)]][cells$model == m] == 1L))
      ranking$fewer_scenarios <-  ranking$n_scenarios < scenarios_of_the_column
      complete <-  ranking[!ranking$fewer_scenarios, ]
      complete <-  complete[order(complete$mean), ]
      complete$position <-  seq_len(nrow(complete))
      complete$mean_of_the_best <-  complete$mean[1]
      complete$MCSE_of_the_best <-  complete$MCSE_of_the_mean[1]
      fewer <-  ranking[ranking$fewer_scenarios, ]
      if (nrow(fewer) > 0) {
            fewer$position <-  NA; fewer$mean_of_the_best <-  NA; fewer$MCSE_of_the_best <-  NA
            for (k in seq_len(nrow(fewer))) {
                  same_scenarios <-  cells$scenario %in% cells$scenario[cells$model == fewer$model[k]]
                  means_over_the_same <-  sapply(complete$model, function(m) mean(cells[[measure]][cells$model == m & same_scenarios]))
                  best_over_the_same <-  complete$model[which.min(means_over_the_same)]
                  fewer$position[k] <-  sum(means_over_the_same < fewer$mean[k]) + 0.5
                  fewer$mean_of_the_best[k] <-  min(means_over_the_same)
                  fewer$MCSE_of_the_best[k] <-  QC_MCSE_of_a_mean(cells[[paste0("MCSE_", measure)]][cells$model == best_over_the_same & same_scenarios], rule_of_the_MCSE_of_a_mean)
            }
      }
      ranking <-  rbind(complete, fewer)
      ranking <-  ranking[order(ranking$position), ]
      ranking$rank <-  seq_len(nrow(ranking))
      difference <-  ranking$mean - ranking$mean_of_the_best
      ranking$statistically_worse <-  difference > 0 & difference >= 4 * pmax(ranking$MCSE_of_the_mean, ranking$MCSE_of_the_best)
      ranking$practically_worse <-  difference > 0 & 100 * difference / ranking$mean_of_the_best >= QC_relative_difference_of_practical_significance
      if (practical_criterion) {
            ranking$mark <-  ifelse(ranking$statistically_worse & ranking$practically_worse, "*",
                                    ifelse(ranking$statistically_worse, "dagger", ifelse(ranking$practically_worse, "practically_only_no_symbol", "")))
      } else {
            ranking$mark <-  ifelse(ranking$statistically_worse, "*", "")
      }
      ranking$position <-  NULL
      rownames(ranking) <-  NULL
      ranking

}
##
## ---- One cell of a table, as a row of a data frame:
##
QC_cell <- function(table,
                    row,
                    column,
                    number = NA_real_,
                    text = NA_character_,
                    marks = "") {
      data.frame(table = table, row = row, column = column, number = number, text = text, marks = marks, stringsAsFactors = FALSE)
}
QC_join_the_marks <- function(...) {
      marks <-  unlist(list(...))
      marks <-  sort(unique(marks[!is.na(marks) & marks != ""]))
      paste(marks, collapse = ",")
}
QC_cells_of_one_ranking <- function(ranking,
                                    table,
                                    group,
                                    name_of_the_value,
                                    name_of_the_best_count,
                                    name_of_the_other_measure,
                                    mark_the_other_measure = NULL) {

      out <-  NULL
      for (k in seq_len(nrow(ranking))) {
            r <-  ranking[k, ]
            m <-  QC_paper_name_of_the_model[r$model]
            out <-  rbind(out,
                          QC_cell(table, m, paste0(group, ": Rank"), number = r$rank, marks = if (r$fewer_scenarios) "ddagger" else ""),
                          QC_cell(table, m, paste0(group, ": ", name_of_the_value), number = r$mean, marks = r$mark),
                          QC_cell(table, m, paste0(group, ": ", name_of_the_best_count), text = sprintf("%d/%d", r$cells_in_group_1, r$n_scenarios)))
            if (!is.null(name_of_the_other_measure)) {
                  other_marks <-  if (is.null(mark_the_other_measure)) "" else mark_the_other_measure$mark[mark_the_other_measure$model == r$model]
                  out <-  rbind(out, QC_cell(table, m, paste0(group, ": ", name_of_the_other_measure), number = r$mean_of_the_other_measure, marks = other_marks))
            }
      }
      out

}
##
## ==== The tables of the main text ===================================================================================================
##
## ---- Table: overall ranking (prior set 0), with the mean RMSE of each family of DGMs:
##
QC_table_overall_ranking <- function(outputs_folder,
                                     family_grouping = "new") {

      table <-  "Table:ranking_of_the_models_overall"
      cells <-  QC_add_the_groups(QC_cells_of_one_prior_set(QC_read_the_cells_of_the_paper(outputs_folder, family_grouping), "set_0"))
      ranking_by_RMSE <-  QC_rank_the_means_of_one_column(cells, "RMSE", TRUE, "root_of_the_sum_of_squares_over_n", TRUE)
      ranking_by_bias <-  QC_rank_the_means_of_one_column(cells, "bias", TRUE, "root_of_the_sum_of_squares_over_n", TRUE)
      out <-  QC_cells_of_one_ranking(ranking_by_RMSE, table, "overall", "RMSE", "Best RMSE", NULL)
      common_scenarios <-  Reduce(intersect, split(cells$scenario, cells$model))
      for (k in seq_len(nrow(ranking_by_RMSE))) {
            model <-  ranking_by_RMSE$model[k]
            b <-  ranking_by_bias[ranking_by_bias$model == model, ]
            out <-  rbind(out,
                          QC_cell(table, QC_paper_name_of_the_model[model], "overall: Bias", number = b$mean, marks = b$mark),
                          QC_cell(table, QC_paper_name_of_the_model[model], "overall: Best bias", text = sprintf("%d/%d", b$cells_in_group_1, b$n_scenarios)))
            for (f in QC_families(family_grouping)) {
                  rows <-  cells$model == model & cells$family == f & cells$scenario %in% common_scenarios
                  out <-  rbind(out, QC_cell(table, QC_paper_name_of_the_model[model], paste0("overall: ", f), number = mean(cells$RMSE[rows])))
            }
      }
      ## (the numbers of the notes: scenarios like for like, and of each family)
      out <-  rbind(out, QC_cell(table, "note", "scenarios in which every model has a result", number = length(common_scenarios)))
      for (f in QC_families(family_grouping)) out <-  rbind(out, QC_cell(table, "note", paste("scenarios of the family", f), number = length(intersect(unique(cells$scenario[cells$family == f]), common_scenarios))))
      out

}
##
## ---- Tables: ranking by share of thresholds missing, and by test (main text, by RMSE with the mean bias beside it):
##
QC_table_ranking_by_share_missing <- function(outputs_folder) {

      table <-  "Table:ranking_of_the_models_by_share_missing_RMSE"
      cells <-  QC_add_the_groups(QC_cells_of_one_prior_set(QC_read_the_cells_of_the_paper(outputs_folder), "set_0"))
      out <-  NULL
      for (p in QC_shares_missing) {
            ranking <-  QC_rank_the_means_of_one_column(cells[cells$pct == p, ], "RMSE", TRUE, "root_of_the_sum_of_squares_over_n", TRUE)
            out <-  rbind(out, QC_cells_of_one_ranking(ranking, table, paste0(p, "% missing"), "RMSE", "Best", "Bias"))
      }
      out

}
QC_table_ranking_by_test <- function(outputs_folder) {

      table <-  "Table:ranking_of_the_models_by_test_RMSE"
      cells <-  QC_add_the_groups(QC_cells_of_one_prior_set(QC_read_the_cells_of_the_paper(outputs_folder), "set_0"))
      out <-  NULL
      for (t in QC_tests) {
            ranking <-  QC_rank_the_means_of_one_column(cells[cells$test == t, ], "RMSE", TRUE, "root_of_the_sum_of_squares_over_n", TRUE)
            out <-  rbind(out, QC_cells_of_one_ranking(ranking, table, QC_paper_name_of_the_test[t], "RMSE", "Best", "Bias"))
      }
      out

}
##
## ---- Table: mean RMSE of each model in each of the nine sets (bold = best of the column; "--" = no result):
##
QC_table_RMSE_in_the_nine_sets <- function(outputs_folder) {

      table <-  "Table:RMSE_of_the_models_in_the_nine_sets"
      cells <-  QC_add_the_groups(QC_cells_of_one_prior_set(QC_read_the_cells_of_the_paper(outputs_folder), "set_0"))
      out <-  NULL
      for (t in QC_tests) for (p in QC_shares_missing) {
            column <-  paste(QC_paper_name_of_the_test[t], paste0(p, "%"))
            of_the_set <-  cells[cells$test == t & cells$pct == p, ]
            ranking <-  QC_rank_the_means_of_one_column(of_the_set, "RMSE", TRUE, "root_of_the_sum_of_squares_over_n", TRUE)
            for (model in names(QC_paper_name_of_the_model)) {
                  i <-  which(ranking$model == model)
                  if (!length(i)) { out <-  rbind(out, QC_cell(table, QC_paper_name_of_the_model[model], column, text = "--")); next }
                  out <-  rbind(out, QC_cell(table, QC_paper_name_of_the_model[model], column, number = ranking$mean[i],
                                             marks = QC_join_the_marks(ranking$mark[i], if (ranking$rank[i] == 1) "bold", if (ranking$fewer_scenarios[i]) "ddagger")))
            }
      }
      out

}
##
## ---- Table: mean coverage (%) and mean width of the 95% intervals, for each share missing (like for like):
##
QC_table_coverage_and_width <- function(outputs_folder) {

      table <-  "Table:coverage_and_width_by_share_missing"
      cells <-  QC_cells_of_one_prior_set(QC_read_the_cells_of_the_paper(outputs_folder), "set_0")
      out <-  NULL
      for (p in QC_shares_missing) {
            of_the_share <-  cells[cells$pct == p, ]
            common_scenarios <-  Reduce(intersect, split(of_the_share$scenario, of_the_share$model))
            for (model in unique(of_the_share$model)) {
                  rows <-  of_the_share[of_the_share$model == model & of_the_share$scenario %in% common_scenarios, ]
                  m <-  QC_paper_name_of_the_model[model]
                  out <-  rbind(out,
                                QC_cell(table, m, paste0(p, "%: Cvg. Se"), number = 100 * mean(rows$avg_coverage_Se)),
                                QC_cell(table, m, paste0(p, "%: Cvg. Sp"), number = 100 * mean(rows$avg_coverage_Sp)),
                                QC_cell(table, m, paste0(p, "%: Width Se"), number = mean(rows$avg_width_Se)),
                                QC_cell(table, m, paste0(p, "%: Width Sp"), number = mean(rows$avg_width_Sp)))
            }
      }
      out

}
##
## ---- Table: mean RMSE and mean bias of the five Bayesian models under prior sets 0, 1 and 2, with the fits kept and lost:
##
QC_table_RMSE_and_bias_by_prior_set <- function(outputs_folder) {

      table <-  "Table:RMSE_and_bias_by_prior_set"
      every_cell <-  QC_read_the_cells_of_the_paper(outputs_folder)
      ## (the five Bayesian models fitted under prior sets 0, 1 and 2: not the Jones model with its own priors, run as set_4)
      bayesian <-  every_cell[!every_cell$frequentist & every_cell$prior_set_of_the_run %in% c("set_0", "set_1", "set_2"), ]
      out <-  NULL
      for (k in c("set_0", "set_1", "set_2")) {
            of_the_set <-  QC_add_the_groups(bayesian[bayesian$prior_set == k, ])
            for (measure in c("RMSE", "bias")) {
                  ranking <-  QC_rank_the_means_of_one_column(of_the_set, measure, TRUE, "root_of_the_sum_of_squares_over_n", TRUE)
                  for (i in seq_len(nrow(ranking))) {
                        out <-  rbind(out, QC_cell(table, QC_paper_name_of_the_model[ranking$model[i]], paste(if (measure == "RMSE") "RMSE" else "Bias", sub("set_", "Set ", k)),
                                                   number = ranking$mean[i], marks = QC_join_the_marks(ranking$mark[i], if (ranking$fewer_scenarios[i]) "ddagger")))
                  }
            }
      }
      ## (fits kept and lost to divergent transitions, from the count of every cell; all scenarios and prior sets 0, 1 and 2)
      fits <-  utils::read.csv(file.path(outputs_folder, "fits_kept_and_lost_of_every_cell_of_the_redo.csv"), stringsAsFactors = FALSE, check.names = FALSE)
      fits <-  fits[fits$prior_set %in% c("set_0", "set_1", "set_2"), ]
      for (m in unique(fits$model)) {
            rows <-  fits[fits$model == m, ]
            out <-  rbind(out,
                          QC_cell(table, m, "Kept", number = 100 * sum(rows$kept) / sum(rows$fits)),
                          QC_cell(table, m, "Div.", number = 100 * sum(rows$lost_to_divergences) / sum(rows$fits)))
      }
      ## (the numbers of the notes)
      percent_lost_in_each_cell <-  100 * fits$lost_to_divergences / fits$fits
      largest <-  which.max(percent_lost_in_each_cell)
      over_the_target <-  every_cell[every_cell$MCSE_val > 0.125 + 0.0005, ]
      out <-  rbind(out,
                    QC_cell(table, "note", "fits of the five models", number = sum(fits$fits)),
                    QC_cell(table, "note", "% discarded due to divergent transitions", number = 100 * sum(fits$lost_to_divergences) / sum(fits$fits)),
                    QC_cell(table, "note", "% under the minimum ESS", number = 100 * sum(fits$under_the_min_ESS) / sum(fits$fits)),
                    QC_cell(table, "note", "% over the Rhat limit", number = 100 * sum(fits$over_the_Rhat_limit) / sum(fits$fits)),
                    QC_cell(table, "note", "% kept", number = 100 * sum(fits$kept) / sum(fits$fits)),
                    QC_cell(table, "note", "largest % discarded due to divergent transitions in one cell", number = percent_lost_in_each_cell[largest],
                            text = paste(fits$model[largest], fits$DGM[largest], fits$test[largest], fits$pct[largest], fits$n_studies[largest], fits$prior_set[largest])),
                    QC_cell(table, "note", "cells over the MCSE target", number = nrow(over_the_target), text = paste(names(table(over_the_target$paper_model)), table(over_the_target$paper_model), collapse = "; ")),
                    QC_cell(table, "note", "cells of the ten models", number = nrow(every_cell)),
                    QC_cell(table, "note", "scenarios", number = length(unique(bayesian$scenario))))
      out

}
##
## ---- Table: mean RMSE (and bias) of the models by DGM (panel A) and by family of DGMs and share missing (panel B); prior set 0.
## ---- Bold = the lowest of the column among the models with a result in every scenario of the column; double dagger = fewer.
##
QC_table_models_by_DGM <- function(outputs_folder,
                                   family_grouping = "new") {

      table <-  "Table:RMSE_of_the_models_by_DGM"
      cells <-  QC_add_the_groups(QC_cells_of_one_prior_set(QC_read_the_cells_of_the_paper(outputs_folder, family_grouping), "set_0"))
      out <-  NULL
      fn_column <- function(rows, measure, column) {
            scenarios_of_the_column <-  length(unique(cells$scenario[rows]))
            models <-  unique(cells$model)
            means <-  sapply(models, function(m) mean(cells[[measure]][rows & cells$model == m]))
            complete <-  sapply(models, function(m) sum(rows & cells$model == m)) == scenarios_of_the_column
            best <-  models[complete][which.min(means[complete])]
            do.call(rbind, lapply(models, function(m) QC_cell(table, QC_paper_name_of_the_model[m], column, number = means[m],
                                                               marks = QC_join_the_marks(if (m == best) "bold", if (!complete[m]) "ddagger"))))
      }
      for (g in QC_DGMs) out <-  rbind(out, fn_column(cells$DGM == g, "RMSE", paste0("A: ", QC_paper_name_of_the_DGM[g])))
      for (f in QC_families(family_grouping)) {
            for (p in QC_shares_missing) out <-  rbind(out, fn_column(cells$family == f & cells$pct == p, "RMSE", paste0("B: ", f, " RMSE ", p, "%")))
            out <-  rbind(out, fn_column(cells$family == f, "RMSE", paste0("B: ", f, " RMSE All")), fn_column(cells$family == f, "bias", paste0("B: ", f, " Bias All")))
            for (m in unique(cells$model)) {
                  rows <-  cells$family == f & cells$model == m
                  out <-  rbind(out, QC_cell(table, QC_paper_name_of_the_model[m], paste0("B: ", f, " Best"), text = sprintf("%d/%d", sum(cells$group_RMSE[rows] == 1L), sum(rows))))
            }
      }
      out

}
##
## ==== The tables of the appendix ====================================================================================================
##
## ---- Appendix rankings (prior set 0, 1 or 2), overall, by share missing or by test, by RMSE or by bias, with the rule of the notes
## ---- of those tables (mean MCSE, statistical criterion only) unless rule = "as_in_the_main_text".
##
QC_table_appendix_ranking <- function(outputs_folder,
                                      prior_set,
                                      measure,
                                      by,
                                      frequentist_models = "left_out",
                                      rule = "as_in_the_notes_of_the_appendix") {

      name_of_the_measure <-  if (measure == "RMSE") "RMSE" else "Bias"
      name_of_the_other <-  if (measure == "RMSE") "Bias" else "RMSE"
      by_in_the_label <-  c(overall = "overall", share = "by_share_missing", test = "by_test")[by]
      table <-  sprintf("Table:ranking_of_the_models_%s_%s%s", by_in_the_label, measure, if (by == "overall" && prior_set != "set_0") paste0("_", prior_set) else "")
      cells <-  QC_add_the_groups(QC_cells_of_one_prior_set(QC_read_the_cells_of_the_paper(outputs_folder), prior_set, frequentist_models))
      like_for_like <-  rule == "as_in_the_main_text"
      rule_of_the_MCSE <-  if (rule == "as_in_the_main_text") "root_of_the_sum_of_squares_over_n" else "mean_of_the_MCSEs"
      practical <-  rule == "as_in_the_main_text"
      groups <-  switch(by, overall = list(overall = rep(TRUE, nrow(cells))),
                            share = stats::setNames(lapply(QC_shares_missing, function(p) cells$pct == p), paste0(QC_shares_missing, "% missing")),
                            test = stats::setNames(lapply(QC_tests, function(t) cells$test == t), QC_paper_name_of_the_test[QC_tests]))
      out <-  NULL
      for (g in names(groups)) {
            ranking <-  QC_rank_the_means_of_one_column(cells[groups[[g]], ], measure, like_for_like, rule_of_the_MCSE, practical)
            out <-  rbind(out, QC_cells_of_one_ranking(ranking, table, g, name_of_the_measure, "Best", if (by == "overall") name_of_the_other else NULL))
      }
      out

}
##
## ---- Detailed appendix table of one test and share missing (prior set 0): per DGM, number of studies and model, RMSE and |bias| of
## ---- Se and Sp with their MCSE, coverage and width, and the group of the 4-group classification by RMSE and by bias. Bold = the
## ---- model that generated the data; a dagger on the model = MCSE over the target (0.125, with an allowance of 0.0005).
##
QC_table_detailed <- function(outputs_folder,
                              test,
                              pct,
                              MCSE_target = 0.125,
                              allowance = 0.0005) {

      table <-  sprintf("table:sim_results_%s_%d_pct", sub("GAD_2", "GAD2", test), pct)
      cells <-  QC_add_the_groups(QC_cells_of_one_prior_set(QC_read_the_cells_of_the_paper(outputs_folder), "set_0"))
      cells <-  cells[cells$test == test & cells$pct == pct, ]
      out <-  vector("list", nrow(cells))
      for (i in seq_len(nrow(cells))) {
            r <-  cells[i, ]
            row <-  paste(QC_paper_name_of_the_DGM[r$DGM], r$n_studies, r$paper_model, sep = " / ")
            bold <-  if (r$model == QC_model_of_the_DGM[r$DGM] || (r$DGM == "Jones" && r$model == "Jones (own priors)")) "bold" else ""
            out[[i]] <-  rbind(
                          QC_cell(table, row, "Model", text = r$paper_model, marks = if (r$MCSE_val > MCSE_target + allowance) "dagger" else ""),
                          QC_cell(table, row, "RMSE Se", number = r$avg_RMSE_Se, marks = bold), QC_cell(table, row, "RMSE Se (MCSE)", number = r$avg_MCSE_RMSE_Se),
                          QC_cell(table, row, "RMSE Sp", number = r$avg_RMSE_Sp, marks = bold), QC_cell(table, row, "RMSE Sp (MCSE)", number = r$avg_MCSE_RMSE_Sp),
                          QC_cell(table, row, "Bias Se", number = r$avg_abs_bias_Se, marks = bold), QC_cell(table, row, "Bias Se (MCSE)", number = r$avg_MCSE_bias_Se),
                          QC_cell(table, row, "Bias Sp", number = r$avg_abs_bias_Sp, marks = bold), QC_cell(table, row, "Bias Sp (MCSE)", number = r$avg_MCSE_bias_Sp),
                          QC_cell(table, row, "Cvg. Se", number = 100 * r$avg_coverage_Se, marks = bold), QC_cell(table, row, "Cvg. Sp", number = 100 * r$avg_coverage_Sp, marks = bold),
                          QC_cell(table, row, "Width Se", number = r$avg_width_Se, marks = bold), QC_cell(table, row, "Width Sp", number = r$avg_width_Sp, marks = bold),
                          QC_cell(table, row, "Group RMSE", number = r$group_RMSE, marks = bold), QC_cell(table, row, "Group bias", number = r$group_bias, marks = bold))
      }
      do.call(rbind, out)

}
##
## ---- Every table, as a named list of data frames (family_grouping applies to the tables that group the DGMs):
##
QC_make_every_table <- function(outputs_folder,
                                family_grouping = "new",
                                frequentist_models_in_sets_1_and_2 = "left_out") {

      tables <-  list()
      tables[["Table:ranking_of_the_models_overall"]] <-  QC_table_overall_ranking(outputs_folder, family_grouping)
      tables[["Table:ranking_of_the_models_by_share_missing_RMSE"]] <-  QC_table_ranking_by_share_missing(outputs_folder)
      tables[["Table:ranking_of_the_models_by_test_RMSE"]] <-  QC_table_ranking_by_test(outputs_folder)
      tables[["Table:RMSE_of_the_models_in_the_nine_sets"]] <-  QC_table_RMSE_in_the_nine_sets(outputs_folder)
      tables[["Table:coverage_and_width_by_share_missing"]] <-  QC_table_coverage_and_width(outputs_folder)
      tables[["Table:RMSE_and_bias_by_prior_set"]] <-  QC_table_RMSE_and_bias_by_prior_set(outputs_folder)
      tables[["Table:RMSE_of_the_models_by_DGM"]] <-  QC_table_models_by_DGM(outputs_folder, family_grouping)
      for (k in c("set_0", "set_1", "set_2")) for (measure in c("RMSE", "bias")) for (by in c("overall", "share", "test")) {
            one <-  QC_table_appendix_ranking(outputs_folder, k, measure, by, if (k == "set_0") "left_out" else frequentist_models_in_sets_1_and_2)
            tables[[paste(k, one$table[1])]] <-  one
      }
      for (t in QC_tests) for (p in QC_shares_missing) {
            one <-  QC_table_detailed(outputs_folder, t, p)
            tables[[one$table[1]]] <-  one
      }
      tables

}
##
## ---- A table of the list, laid out wide (rows x columns) as text, for reading:
##
QC_show_one_table_wide <- function(cells_of_the_table) {

      shown <-  ifelse(is.na(cells_of_the_table$number), cells_of_the_table$text, sprintf("%.2f", cells_of_the_table$number))
      shown <-  ifelse(cells_of_the_table$marks != "", paste0(shown, " [", cells_of_the_table$marks, "]"), shown)
      wide <-  tapply(shown, list(factor(cells_of_the_table$row, unique(cells_of_the_table$row)), factor(cells_of_the_table$column, unique(cells_of_the_table$column))), function(x) x[1])
      as.data.frame(wide, stringsAsFactors = FALSE)

}
##
## ==== Reading the tables of a LaTeX file ===========================================================================================
##
QC_active_lines_of_a_tex_file <- function(tex_file) {
      lines <-  readLines(tex_file, warn = FALSE)
      lines[!grepl("^\\s*%", lines)]
}
## (the lines of the table with this \label; with label = NULL, of the first table of the lines)
QC_lines_of_one_table <- function(lines,
                                  label = NULL) {

      at <-  if (is.null(label)) grep("\\\\begin\\{(table|longtable)\\}", lines)[1] else grep(paste0("\\label{", label, "}"), lines, fixed = TRUE)[1]
      if (is.na(at)) return(NULL)
      begin <-  max(grep("\\\\begin\\{(table|longtable)\\}", lines[seq_len(at)]))
      environment <-  if (grepl("\\\\begin\\{longtable\\}", lines[begin])) "longtable" else "table"
      end <-  begin - 1 + grep(sprintf("\\end{%s}", environment), lines[begin:length(lines)], fixed = TRUE)[1]
      lines[begin:end]

}
QC_split_a_row <- function(line) {
      line <-  sub("\\\\\\\\\\s*$", "", line)
      line <-  gsub("\\&", "<AMPERSAND>", line, fixed = TRUE)
      cells <-  trimws(strsplit(line, "&", fixed = TRUE)[[1]])
      gsub("<AMPERSAND>", "\\&", cells, fixed = TRUE)
}
## (the body rows of each tabular, or of the longtable, as vectors of cells)
QC_rows_of_each_block <- function(table_lines) {

      if (any(grepl("\\endlastfoot", table_lines, fixed = TRUE))) {
            from <-  grep("\\endlastfoot", table_lines, fixed = TRUE)[1] + 1
            to <-  grep("\\end{longtable}", table_lines, fixed = TRUE)[1] - 1
            body <-  table_lines[from:to]
            return(list(lapply(body[grepl("&", body) & grepl("\\\\\\\\\\s*$", body)], QC_split_a_row)))
      }
      starts <-  grep("\\begin{tabular}", table_lines, fixed = TRUE)
      blocks <-  list()
      for (s in starts) {
            mid <-  s - 1 + grep("\\midrule", table_lines[s:length(table_lines)], fixed = TRUE)[1]
            bottom <-  s - 1 + grep("\\bottomrule", table_lines[s:length(table_lines)], fixed = TRUE)[1]
            body <-  table_lines[(mid + 1):(bottom - 1)]
            blocks[[length(blocks) + 1]] <-  lapply(body[grepl("&", body)], QC_split_a_row)
      }
      blocks

}
## (the number, its decimals, the text, the marks and the MCSE in brackets of one cell)
QC_parse_a_cell <- function(x) {

      marks <-  c(if (grepl("^{*}", x, fixed = TRUE)) "*", if (grepl("\\ddagger", x, fixed = TRUE)) "ddagger",
                  if (grepl("\\dagger", gsub("\\ddagger", "", x, fixed = TRUE), fixed = TRUE)) "dagger",
                  if (grepl("\\\\mathbf|\\\\textbf", x)) "bold")
      MCSE <-  NA_real_; MCSE_decimals <-  NA_integer_
      in_brackets <-  regmatches(x, regexec("\\\\scriptsize\\(([0-9.]+)\\)", x))[[1]]
      if (length(in_brackets)) {
            MCSE <-  as.numeric(in_brackets[2]); MCSE_decimals <-  nchar(sub("^[^.]*\\.?", "", in_brackets[2]))
            x <-  sub("\\{?\\\\scriptsize\\([0-9.]+\\)\\}?", "", x)
      }
      bare <-  gsub("\\^\\{[^}]*\\}", "", x)
      bare <-  gsub("\\\\mathbf|\\\\textbf", "", bare)
      bare <-  trimws(gsub("[${}]", "", bare))
      number <-  NA_real_; decimals <-  NA_integer_; text <-  bare
      if (grepl("^[0-9]+/[0-9]+$", bare)) {
            text <-  bare
      } else if (grepl("^-?[0-9]+(\\.[0-9]+)?$", bare)) {
            number <-  as.numeric(bare); decimals <-  nchar(sub("^[^.]*\\.?", "", bare)); text <-  NA_character_
      }
      list(number = number, decimals = decimals, text = text, marks = QC_join_the_marks(marks), MCSE = MCSE, MCSE_decimals = MCSE_decimals)

}
QC_cell_of_the_paper <- function(table,
                                 row,
                                 column,
                                 parsed,
                                 only_the_number = FALSE) {
      data.frame(table = table, row = row, column = column, number = parsed$number, decimals = parsed$decimals,
                 text = if (only_the_number) NA_character_ else parsed$text, marks = parsed$marks, stringsAsFactors = FALSE)
}
## (a ranking table: the names of the columns in order, the group of each column; the model cell gives the row and the double dagger)
QC_paper_cells_of_a_ranking_table <- function(rows,
                                              table,
                                              names_of_the_columns,
                                              group_of_each_column,
                                              columns_left_out = c("Reference", "R package")) {

      out <-  NULL
      for (cells in rows) {
            if (length(cells) != length(names_of_the_columns)) stop(sprintf("%s: a row of %d cells, %d expected", table, length(cells), length(names_of_the_columns)))
            rank <-  QC_parse_a_cell(cells[names_of_the_columns == "Rank"])$number
            for (g in unique(group_of_each_column[!is.na(group_of_each_column)])) {
                  of_g <-  which(group_of_each_column == g)
                  model_cell <-  QC_parse_a_cell(cells[of_g[names_of_the_columns[of_g] == "Model"]])
                  model <-  trimws(gsub("\\s+", " ", model_cell$text))
                  if (is.na(model) || model == "") next
                  out <-  rbind(out, data.frame(table = table, row = model, column = paste0(g, ": Rank"), number = rank, decimals = 0L, text = NA_character_, marks = model_cell$marks, stringsAsFactors = FALSE))
                  for (j in of_g[!(names_of_the_columns[of_g] %in% c("Model", "Rank", columns_left_out))]) {
                        out <-  rbind(out, QC_cell_of_the_paper(table, model, paste0(g, ": ", names_of_the_columns[j]), QC_parse_a_cell(cells[j])))
                  }
            }
      }
      out

}
## (a table with one model a row: the first cell is the model)
QC_paper_cells_of_a_model_table <- function(rows,
                                            table,
                                            names_of_the_columns) {

      out <-  NULL
      for (cells in rows) {
            if (length(cells) != length(names_of_the_columns) + 1) stop(sprintf("%s: a row of %d cells, %d expected", table, length(cells), length(names_of_the_columns) + 1))
            model <-  trimws(QC_parse_a_cell(cells[1])$text)
            for (j in seq_along(names_of_the_columns)) out <-  rbind(out, QC_cell_of_the_paper(table, model, names_of_the_columns[j], QC_parse_a_cell(cells[j + 1])))
      }
      out

}
## (a detailed table: DGM and n are carried down from the rows above)
QC_paper_cells_of_a_detailed_table <- function(rows,
                                               table) {

      names_of_the_values <-  c("RMSE Se", "RMSE Sp", "Bias Se", "Bias Sp", "Cvg. Se", "Cvg. Sp", "Width Se", "Width Sp", "Group RMSE", "Group bias")
      out <-  list()
      DGM <-  NA; n <-  NA
      for (cells in rows) {
            if (length(cells) != 13) stop(sprintf("%s: a row of %d cells, 13 expected", table, length(cells)))
            if (cells[1] != "") DGM <-  cells[1]
            if (cells[2] != "") n <-  cells[2]
            model_cell <-  QC_parse_a_cell(cells[3])
            row <-  paste(DGM, n, model_cell$text, sep = " / ")
            out[[length(out) + 1]] <-  data.frame(table = table, row = row, column = "Model", number = NA_real_, decimals = NA_integer_, text = model_cell$text, marks = model_cell$marks, stringsAsFactors = FALSE)
            for (j in 4:13) {
                  parsed <-  QC_parse_a_cell(cells[j])
                  out[[length(out) + 1]] <-  QC_cell_of_the_paper(table, row, names_of_the_values[j - 3], parsed)
                  if (!is.na(parsed$MCSE)) out[[length(out) + 1]] <-  (data.frame(table = table, row = row, column = paste0(names_of_the_values[j - 3], " (MCSE)"), number = parsed$MCSE,
                                                                          decimals = parsed$MCSE_decimals, text = NA_character_, marks = "", stringsAsFactors = FALSE))
            }
      }
      do.call(rbind, out)

}
## (numbers of the notes of a table, by patterns; a pattern not found gives no cell)
QC_paper_cells_of_the_notes <- function(table_lines,
                                        table,
                                        patterns) {

      text <-  paste(table_lines, collapse = " ")
      out <-  NULL
      for (name in names(patterns)) {
            found <-  regmatches(text, regexec(patterns[[name]], text, perl = TRUE))[[1]]
            if (!length(found)) next
            value <-  gsub("\\{,\\}", "", found[2])
            out <-  rbind(out, data.frame(table = table, row = "note", column = name, number = as.numeric(value),
                                          decimals = nchar(sub("^[^.]*\\.?", "", value)), text = NA_character_, marks = "", stringsAsFactors = FALSE))
      }
      out

}
##
## ---- The cells of one table of a .tex file, by the kind of the table. label = the \label to find; NULL = the first table of the file.
##
QC_read_one_table_of_a_tex_file <- function(lines,
                                            label,
                                            kind,
                                            table = label,
                                            family_grouping = "old") {

      table_lines <-  QC_lines_of_one_table(lines, label)
      if (is.null(table_lines)) return(NULL)
      blocks <-  QC_rows_of_each_block(table_lines)
      rows <-  blocks[[1]]
      g3 <- function(heads, fields) rep(heads, each = length(fields))
      if (kind == "overall_main") {
            families <-  QC_families(family_grouping)
            out <-  QC_paper_cells_of_a_ranking_table(rows, table, c("Rank", "Model", "Reference", "R package", "RMSE", "Best RMSE", "Bias", "Best bias", families),
                                                     c(NA, rep("overall", 7 + length(families))))
            patterns <-  list("scenarios in which every model has a result" = "across the \\$(\\d+)\\$ scenarios in which every model has a result")
            if (family_grouping == "old") patterns <-  c(patterns, list("scenarios of the family Jones" = "across the \\$(\\d+)\\$ scenarios of the Jones DGM",
                                                                        "scenarios of the family Ordinal" = "the \\$(\\d+)\\$ scenarios of the three ordinal DGMs",
                                                                        "scenarios of the family Competing" = "the \\$(\\d+)\\$ scenarios of the three competing DGMs"))
            return(rbind(out, QC_paper_cells_of_the_notes(table_lines, table, patterns)))
      }
      if (kind %in% c("share_main", "test_main")) {
            heads <-  if (kind == "share_main") paste0(QC_shares_missing, "% missing") else unname(QC_paper_name_of_the_test[QC_tests])
            fields <-  c("Model", "RMSE", "Best", "Bias")
            return(QC_paper_cells_of_a_ranking_table(rows, table, c("Rank", rep(fields, 3)), c(NA, g3(heads, fields))))
      }
      if (kind %in% c("share_appendix_RMSE", "share_appendix_bias", "test_appendix_RMSE", "test_appendix_bias")) {
            heads <-  if (grepl("^share", kind)) paste0(QC_shares_missing, "% missing") else unname(QC_paper_name_of_the_test[QC_tests])
            fields <-  c("Model", if (grepl("RMSE$", kind)) "RMSE" else "Bias", "Best")
            return(QC_paper_cells_of_a_ranking_table(rows, table, c("Rank", rep(fields, 3)), c(NA, g3(heads, fields))))
      }
      if (kind %in% c("overall_appendix_RMSE", "overall_appendix_bias")) {
            fields <-  if (kind == "overall_appendix_RMSE") c("RMSE", "Best", "Bias") else c("Bias", "Best", "RMSE")
            return(QC_paper_cells_of_a_ranking_table(rows, table, c("Rank", "Model", "Reference", "R package", fields), c(NA, rep("overall", 6))))
      }
      if (kind == "nine_sets") {
            return(QC_paper_cells_of_a_model_table(rows, table, as.vector(t(outer(QC_paper_name_of_the_test[QC_tests], paste0(QC_shares_missing, "%"), paste)))))
      }
      if (kind == "coverage") {
            return(QC_paper_cells_of_a_model_table(rows, table, as.vector(sapply(QC_shares_missing, function(p) paste0(p, "%: ", c("Cvg. Se", "Cvg. Sp", "Width Se", "Width Sp"))))))
      }
      if (kind == "prior_sets") {
            out <-  QC_paper_cells_of_a_model_table(rows, table, c("RMSE Set 0", "RMSE Set 1", "RMSE Set 2", "Bias Set 0", "Bias Set 1", "Bias Set 2", "Kept", "Div."))
            patterns <-  list("fits of the five models" = "all five models: \\$([0-9{},]+)\\$ fits",
                              "% discarded due to divergent transitions" = "fits, \\$([0-9.]+)\\\\%\\$ discarded due to divergent",
                              "% under the minimum ESS" = "\\$([0-9.]+)\\\\%\\$ because the minimum",
                              "% over the Rhat limit" = "\\$([0-9.]+)\\\\%\\$ because the largest",
                              "% kept" = "and \\$([0-9.]+)\\\\%\\$ kept",
                              "largest % discarded due to divergent transitions in one cell" = "in one cell: \\$([0-9.]+)\\\\%\\$",
                              "cells over the MCSE target" = "stopped: \\$(\\d+)\\$ of the",
                              "cells of the ten models" = "of the \\$([0-9{},]+)\\$ cells of the",
                              "scenarios" = "across the \\$(\\d+)\\$ scenarios, of RMSE")
            return(rbind(out, QC_paper_cells_of_the_notes(table_lines, table, patterns)))
      }
      if (kind == "models_by_DGM") {
            ## (panel A: one column a DGM, in the order of the families: "old" = Jones first, as QC_DGMs; "new" = the ordinal DGMs, then the
            ## continuous ones; the header row of the table must name the DGMs in that order, or the cells would be read into the wrong DGM)
            DGMs_in_the_order_of_the_columns <-  if (family_grouping == "new") {
                  c(QC_DGMs[QC_family_of_each_DGM(QC_DGMs, "new") == "Ordinal"], QC_DGMs[QC_family_of_each_DGM(QC_DGMs, "new") == "Continuous"])
            } else QC_DGMs
            header <-  grep("\\\\textbf\\{Model\\}", table_lines, value = TRUE)[1]
            if (!is.na(header)) {
                  named <-  regmatches(header, gregexpr("\\\\textbf\\{[^}]*\\}", header))[[1]]
                  named <-  gsub("^\\\\textbf\\{|\\}$", "", named)[-1]
                  expected <-  unname(QC_paper_name_of_the_DGM[DGMs_in_the_order_of_the_columns])
                  if (!identical(named, expected)) stop(sprintf("%s: the header of panel A names the DGMs %s, the %s grouping expects %s", table,
                                                                paste(named, collapse = ", "), family_grouping, paste(expected, collapse = ", ")))
            }
            out <-  QC_paper_cells_of_a_model_table(blocks[[1]], table, paste0("A: ", QC_paper_name_of_the_DGM[DGMs_in_the_order_of_the_columns]))
            families <-  QC_families(family_grouping)
            names_B <-  unlist(lapply(families, function(f) c(paste0("B: ", f, " RMSE ", QC_shares_missing, "%"), paste0("B: ", f, " RMSE All"), paste0("B: ", f, " Bias All"))))
            ## (the old table has one Best column, after the last family)
            names_B <-  c(names_B, paste0("B: ", tail(families, 1), " Best"))
            return(rbind(out, QC_paper_cells_of_a_model_table(blocks[[2]], table, names_B)))
      }
      if (kind == "detailed") return(QC_paper_cells_of_a_detailed_table(rows, table))
      stop("unknown kind of table: ", kind)

}
##
## ---- Comparison of the cells of a table of a .tex file with the cells made here. A number matches when it is within the rounding
## ---- of the decimals shown (half a unit of the last decimal); a text ("x/y", a model, "--") matches when it is the same; the marks
## ---- (star, dagger, double dagger, bold) must be the same.
##
QC_compare_one_table <- function(cells_made_here,
                                 cells_of_the_tex) {

      if (is.null(cells_of_the_tex) || !nrow(cells_of_the_tex)) return(NULL)
      p <-  cells_of_the_tex
      j <-  match(paste(p$row, p$column, sep = " || "), paste(cells_made_here$row, cells_made_here$column, sep = " || "))
      h <-  cells_made_here[j, ]
      found <-  !is.na(j)
      difference <-  ifelse(found & !is.na(p$number), h$number - p$number, NA_real_)
      tolerance <-  0.5 * 10^(-ifelse(is.na(p$decimals), 0, p$decimals)) + 1e-9
      number_matches <-  ifelse(!is.na(p$number), !is.na(difference) & abs(difference) <= tolerance,
                                found & !is.na(h$text) & trimws(as.character(p$text)) == trimws(as.character(h$text)))
      number_matches[is.na(number_matches)] <-  FALSE
      marks_match <-  found & p$marks == ifelse(is.na(h$marks), "", h$marks)
      shown <-  mapply(QC_shown, p$number, p$decimals, p$text, p$marks)
      made_here <-  ifelse(found, mapply(QC_shown, h$number, ifelse(is.na(p$decimals), 4L, p$decimals + 2L), h$text, h$marks), "(no such cell made here)")
      data.frame(table = p$table, row = p$row, column = p$column, shown = unname(shown), made_here = unname(made_here),
                 number_matches = number_matches, marks_match = marks_match, difference = difference, stringsAsFactors = FALSE)

}
QC_shown <- function(number,
                     decimals,
                     text,
                     marks) {
      value <-  if (!is.na(number)) sprintf(paste0("%.", if (is.na(decimals)) 2 else decimals, "f"), number) else as.character(text)
      if (!is.na(marks) && marks != "") value <-  paste0(value, " [", marks, "]")
      value
}
##
## ---- Every table of the paper, with its label, kind and the cells made here; and the .tex file the tables were generated in.
##
QC_tables_of_the_paper <- function() {

      tables <-  data.frame(label = c("Table:ranking_of_the_models_overall", "Table:ranking_of_the_models_by_share_missing_RMSE", "Table:ranking_of_the_models_by_test_RMSE",
                                      "Table:RMSE_of_the_models_in_the_nine_sets", "Table:coverage_and_width_by_share_missing", "Table:RMSE_and_bias_by_prior_set",
                                      "Table:ranking_of_the_models_by_share_missing_bias", "Table:ranking_of_the_models_by_test_bias",
                                      "Table:ranking_of_the_models_overall_set_1", "Table:ranking_of_the_models_overall_set_2"),
                            kind = c("overall_main", "share_main", "test_main", "nine_sets", "coverage", "prior_sets", "share_appendix_bias", "test_appendix_bias",
                                     "overall_appendix_RMSE", "overall_appendix_RMSE"),
                            key_made_here = c("Table:ranking_of_the_models_overall", "Table:ranking_of_the_models_by_share_missing_RMSE", "Table:ranking_of_the_models_by_test_RMSE",
                                              "Table:RMSE_of_the_models_in_the_nine_sets", "Table:coverage_and_width_by_share_missing", "Table:RMSE_and_bias_by_prior_set",
                                              "set_0 Table:ranking_of_the_models_by_share_missing_bias", "set_0 Table:ranking_of_the_models_by_test_bias",
                                              "set_1 Table:ranking_of_the_models_overall_RMSE_set_1", "set_2 Table:ranking_of_the_models_overall_RMSE_set_2"),
                            generated_file = c("main_text_tables_four_groups/main_table_overall_ranking.tex", "main_text_tables_four_groups/main_table_ranking_by_share_missing.tex",
                                               "main_text_tables_four_groups/main_table_ranking_by_test.tex", "main_text_tables_four_groups/main_table_RMSE_in_the_nine_sets.tex",
                                               "main_text_tables_four_groups/main_table_coverage_and_width.tex", "main_text_tables_four_groups/main_table_prior_sets.tex",
                                               "ranking_tables_set_0/ranking_of_the_models_by_share_missing_bias.tex", "ranking_tables_set_0/ranking_of_the_models_by_test_bias.tex",
                                               "ranking_tables_set_1/ranking_of_the_models_overall_RMSE.tex", "ranking_tables_set_2/ranking_of_the_models_overall_RMSE.tex"),
                            stringsAsFactors = FALSE)
      for (t in QC_tests) for (p in QC_shares_missing) {
            label <-  sprintf("table:sim_results_%s_%d_pct", sub("GAD_2", "GAD2", t), p)
            tables <-  rbind(tables, data.frame(label = label, kind = "detailed", key_made_here = label, generated_file = sprintf("detailed_tables_of_the_appendix/detailed_table_%s_%d_pct.tex", t, p), stringsAsFactors = FALSE))
      }
      tables

}
## (the generated tables that are not in the paper: checked against their own .tex files only)
QC_generated_tables_not_in_the_paper <- function() {

      out <-  data.frame(kind = "models_by_DGM", key_made_here = "Table:RMSE_of_the_models_by_DGM", generated_file = "ranking_tables_set_0/table_of_the_models_by_DGM.tex", stringsAsFactors = FALSE)
      for (k in c("set_0", "set_1", "set_2")) for (measure in c("RMSE", "bias")) for (by in c("overall", "share", "test")) {
            if (k == "set_0" && measure == "bias" && by != "overall") next
            if (k != "set_0" && measure == "RMSE" && by == "overall") next
            by_in_the_label <-  c(overall = "overall", share = "by_share_missing", test = "by_test")[by]
            key <-  sprintf("%s Table:ranking_of_the_models_%s_%s%s", k, by_in_the_label, measure, if (by == "overall" && k != "set_0") paste0("_", k) else "")
            kind <-  if (by == "overall") paste0("overall_appendix_", measure) else paste0(by, "_appendix_", measure)
            out <-  rbind(out, data.frame(kind = kind, key_made_here = key, generated_file = sprintf("ranking_tables_%s/ranking_of_the_models_%s_%s.tex", k, by_in_the_label, measure), stringsAsFactors = FALSE))
      }
      out

}
##
## ---- The whole check: every table of the paper against the cells made here (and the generated .tex files of the outputs folder too),
## ---- and a report in Markdown. family_grouping_of_the_tex = the grouping of the DGMs of the .tex being checked.
##
QC_compare_every_table_with_the_paper <- function(outputs_folder,
                                                  paper_tex_file,
                                                  report_file,
                                                  family_grouping_of_the_tex = "old") {

      made_here <-  QC_make_every_table(outputs_folder, family_grouping_of_the_tex, "left_out")
      made_here_with_the_frequentist_cells <-  QC_make_every_table(outputs_folder, family_grouping_of_the_tex, "included_with_their_set_0_cells")
      made_here_with_the_main_text_rule <-  list()
      for (k in c("set_0", "set_1", "set_2")) for (measure in c("RMSE", "bias")) for (by in c("overall", "share", "test")) {
            one <-  QC_table_appendix_ranking(outputs_folder, k, measure, by, if (k == "set_0") "left_out" else "included_with_their_set_0_cells", "as_in_the_main_text")
            made_here_with_the_main_text_rule[[paste(k, one$table[1])]] <-  one
      }
      paper_lines <-  QC_active_lines_of_a_tex_file(paper_tex_file)
      tables <-  QC_tables_of_the_paper()
      every_comparison <-  list()
      for (i in seq_len(nrow(tables))) {
            label <-  tables$label[i]
            mine <-  made_here[[tables$key_made_here[i]]]
            mine$table <-  label
            paper <-  QC_read_one_table_of_a_tex_file(paper_lines, label, tables$kind[i], label, family_grouping_of_the_tex)
            every_comparison[[paste("paper:", label)]] <-  QC_compare_one_table(mine, paper)
            generated_file <-  file.path(outputs_folder, tables$generated_file[i])
            if (file.exists(generated_file)) {
                  generated <-  QC_read_one_table_of_a_tex_file(QC_active_lines_of_a_tex_file(generated_file), NULL, tables$kind[i], label, family_grouping_of_the_tex)
                  every_comparison[[paste("generated .tex:", label)]] <-  QC_compare_one_table(mine, generated)
            }
            if (label %in% c("Table:ranking_of_the_models_overall_set_1", "Table:ranking_of_the_models_overall_set_2")) {
                  with_frequentist <-  made_here_with_the_frequentist_cells[[tables$key_made_here[i]]]; with_frequentist$table <-  label
                  every_comparison[[paste("paper, against the frequentist cells of set 0 put in as the 19 Nov script did:", label)]] <-  QC_compare_one_table(with_frequentist, paper)
            }
            if (tables$kind[i] %in% c("share_appendix_bias", "test_appendix_bias", "overall_appendix_RMSE")) {
                  main_rule <-  made_here_with_the_main_text_rule[[tables$key_made_here[i]]]; main_rule$table <-  label
                  every_comparison[[paste("paper, against the rule of the main text (like for like, sqrt(sum MCSE^2)/n, star and dagger):", label)]] <-  QC_compare_one_table(main_rule, paper)
            }
      }
      others <-  QC_generated_tables_not_in_the_paper()
      for (i in seq_len(nrow(others))) {
            generated_file <-  file.path(outputs_folder, others$generated_file[i])
            if (!file.exists(generated_file)) next
            mine <-  made_here[[others$key_made_here[i]]]
            generated <-  QC_read_one_table_of_a_tex_file(QC_active_lines_of_a_tex_file(generated_file), NULL, others$kind[i], mine$table[1], family_grouping_of_the_tex)
            every_comparison[[paste("generated .tex only (not in the paper):", others$generated_file[i])]] <-  QC_compare_one_table(mine, generated)
            if (grepl("^set_[12]", others$key_made_here[i])) {
                  with_frequentist <-  made_here_with_the_frequentist_cells[[others$key_made_here[i]]]
                  every_comparison[[paste("generated .tex only (not in the paper), against the frequentist cells of set 0 put in as the 19 Nov script did:", others$generated_file[i])]] <-  QC_compare_one_table(with_frequentist, generated)
            }
      }
      left_out <-  QC_cells_left_out_of_the_paper(outputs_folder)
      QC_write_the_report(every_comparison, left_out, report_file, paper_tex_file, outputs_folder)
      invisible(every_comparison)

}
QC_write_the_report <- function(every_comparison,
                                left_out,
                                report_file,
                                paper_tex_file,
                                outputs_folder) {

      lines <-  c("# QC of every results table of the simulation study against the scored outputs", "",
                  sprintf("Made by QC_functions_that_make_every_table_from_the_outputs_19_Nov.R on %s.", format(Sys.time(), "%d %b %Y %H:%M")),
                  sprintf("Outputs: %s", outputs_folder), sprintf("Paper: %s", paper_tex_file), "",
                  "A number matches when it is within the rounding of the decimals shown; a text (x/y, model, --) when it is the same; the marks (*, dagger, ddagger, bold) must be the same.", "",
                  "## Counts", "", "| comparison | cells checked | numbers or texts not matching | marks not matching |", "|---|---|---|---|")
      for (name in names(every_comparison)) {
            x <-  every_comparison[[name]]
            if (is.null(x)) { lines <-  c(lines, sprintf("| %s | table not found | | |", name)); next }
            lines <-  c(lines, sprintf("| %s | %d | %d | %d |", name, nrow(x), sum(!x$number_matches), sum(!x$marks_match)))
      }
      lines <-  c(lines, "", "## Every cell not matching", "")
      for (name in names(every_comparison)) {
            x <-  every_comparison[[name]]
            if (is.null(x)) next
            bad <-  x[!x$number_matches | !x$marks_match, ]
            if (!nrow(bad)) next
            lines <-  c(lines, sprintf("### %s (%d cells)", name, nrow(bad)), "", "| row | column | in the .tex | made here |", "|---|---|---|---|",
                        sprintf("| %s | %s | %s | %s |", bad$row, bad$column, bad$shown, bad$made_here), "")
      }
      lines <-  c(lines, "## Cells of the redo left out of the paper", "",
                  sprintf("Rest of the table of the nine sets identical to the table of the redo: %s; cells of the redo with over 50%% of fits not kept: %d.",
                          attr(left_out, "rest_identical_to_the_redo"), attr(left_out, "cells_with_over_50_pct_not_kept")), "",
                  "| model | test | pct | DGM | n | kept fits | fits tried | % not kept |", "|---|---|---|---|---|---|---|---|",
                  sprintf("| %s | %s | %s | %s | %s | %s | %s | %s |", left_out$model, left_out$test, left_out$pct, left_out$DGM, left_out$n_studies, left_out$kept_fits, left_out$fits_tried, left_out$percent_of_fits_not_kept), "")
      writeLines(lines, report_file)

}
##
## ---- Run from the command line: Rscript <this file> <outputs folder> <paper .tex file> <report .md file> [old|new]
##
if (sys.nframe() == 0L) {
      arguments <-  commandArgs(trailingOnly = TRUE)
      if (length(arguments) >= 3) {
            family_grouping_of_the_tex <-  if (length(arguments) >= 4) arguments[4] else "new"
            comparisons <-  QC_compare_every_table_with_the_paper(arguments[1], arguments[2], arguments[3], family_grouping_of_the_tex)
            for (name in names(comparisons)) {
                  x <-  comparisons[[name]]
                  cat(sprintf("%-150s checked %5s  numbers/texts off %4s  marks off %4s\n", substr(name, 1, 150), if (is.null(x)) "-" else nrow(x),
                              if (is.null(x)) "-" else sum(!x$number_matches), if (is.null(x)) "-" else sum(!x$marks_match)))
            }
      }
}
