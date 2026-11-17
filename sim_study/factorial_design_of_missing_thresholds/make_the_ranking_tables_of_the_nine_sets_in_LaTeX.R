## Makes the ranking tables of the paper (LaTeX) from the table of every cell of the nine sets of the factorial design
## (3 tests x 3 shares of thresholds missing): ranking by mean RMSE and by mean bias (1) for each share missing, (2) for each test,
## (3) overall; each for one prior set. Also writes every number of the tables as one .csv file.
##
## The star: a model is marked when its mean is worse than the mean of the best-ranked model of the same column by at least
## 4 x MCSE, i.e. the rule for statistical significance of the 4-group tables of the paper, applied to the means:
##   MCSE of one cell  =  max(MCSE of Se, MCSE of Sp)      (as in the 4-group tables)
##   MCSE of a model   =  mean of that over the scenarios
##   star if (mean of the model - mean of the best) >= 4 * max(MCSE of the model, MCSE of the best).
## "Best" = the number of scenarios in which the model is in "Group 1 - Best" of the 4-group table of that scenario.
## Every scored cell counts (a cell over the MCSE target too). A model's mean is over the scenarios it has a result in; a model
## with fewer scenarios than the column has carries a double dagger. The frequentist models have no priors: their cells of prior
## set 0 are used at every prior set.
##
## Usage: Rscript <this file> <table of every cell of the nine sets, .rds> <folder to write in> <prior set, e.g. set_0>
##                <relative difference of practical significance in percent, e.g. 10>
##
arguments <-  commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 4)
table_of_the_cells <-  readRDS(arguments[1])
folder <-  arguments[2]
prior_set <-  arguments[3]
relative_difference_threshold <-  as.numeric(arguments[4])
dir.create(folder, recursive = TRUE, showWarnings = FALSE)
##
frequentist <-  table_of_the_cells$model %in% c("Hoyer 2018", "logit GLMM", "diagmeta")
d <-  table_of_the_cells[table_of_the_cells$prior_set == prior_set | (frequentist & table_of_the_cells$prior_set == "set_0"), ]
d$scenario <-  paste(d$test, d$pct, d$DGM, d$n_studies)
cat("cells:", nrow(d), " scenarios:", length(unique(d$scenario)), " models:", length(unique(d$model)), "\n")
print(table(d$model))
##
## ---- Names, as the paper has them:
##
name_of_the_model <-  c("strat. bivariate" = "Strat-biv", "Jones" = "Jones", "O-biv-FC" = "O-biv-FC", "O-biv-RC" = "O-biv-RC", "ord-HSROC" = "O-HSROC-RC",
                        "Hoyer 2018" = "Hoyer", "logit GLMM" = "Logit-GLMM", "diagmeta" = "Steinhauser")
reference_of_the_model <-  c("strat. bivariate" = "Reitsma et al, 2005\\supercite{Reitsma2005}", "Jones" = "Jones et al, 2019\\supercite{Jones2019}",
                             "O-biv-FC" = "Proposed", "O-biv-RC" = "Proposed", "ord-HSROC" = "Proposed",
                             "Hoyer 2018" = "Hoyer et al, 2018\\supercite{Hoyer_et_al_2018_mult_thr}",
                             "logit GLMM" = "Hoyer \\& Kuss, 2018\\supercite{Hoyer_Kuss_2018_GLMM}",
                             "diagmeta" = "Steinhauser et al, 2016\\supercite{Steinhauser2016}")
package_of_the_model <-  c("strat. bivariate" = "MetaOrdDTA", "Jones" = "MetaOrdDTA", "O-biv-FC" = "MetaOrdDTA", "O-biv-RC" = "MetaOrdDTA", "ord-HSROC" = "MetaOrdDTA",
                           "Hoyer 2018" = "dtametaTMB", "logit GLMM" = "lme4", "diagmeta" = "diagmeta")
stopifnot(all(d$model %in% names(name_of_the_model)))
tests <-  c("GAD_2", "HADS", "BAI")
shares <-  c(15, 55, 85)
name_of_the_test <-  c(GAD_2 = "GAD-2", HADS = "HADS", BAI = "BAI")
name_of_the_set <-  c(set_0 = "prior set 0", set_1 = "prior set 1 (well-specified)", set_2 = "prior set 2 (misspecified)")
fn_2 <- function(x) sprintf("%.2f", x)
##
## ---- The value and the MCSE of every cell for one measure, and whether the cell is in "Group 1 - Best" (the rule of the 4-group tables):
##
fn_cells_of_one_measure <- function(d,
                                    measure) {

      if (measure == "RMSE") {
            d$value <-  d$avg_RMSE
            d$MCSE_actual <-  pmax(d$avg_MCSE_RMSE_Se, d$avg_MCSE_RMSE_Sp, na.rm = TRUE)
            d$other <-  d$avg_abs_bias
      } else {
            d$value <-  d$avg_abs_bias
            d$MCSE_actual <-  pmax(d$avg_MCSE_bias_Se, d$avg_MCSE_bias_Sp, na.rm = TRUE)
            d$other <-  d$avg_RMSE
      }
      d$in_the_best_group <-  NA
      for (s in unique(d$scenario)) {
            i <-  which(d$scenario == s)
            best <-  i[which.min(d$value[i])]
            statistically <-  (d$value[i] - d$value[best]) >= 4 * pmax(d$MCSE_actual[i], d$MCSE_actual[best])
            practically <-  100 * (d$value[i] - d$value[best]) / d$value[best] >= relative_difference_threshold
            d$in_the_best_group[i] <-  !statistically & !practically
      }
      d

}
##
## ---- The order of the models of one column. A model with a result in every scenario is ranked by its mean, with a star when it
## ---- is worse than the best such model by at least 4 x max(MCSE of the model, MCSE of the best). A model with no result in some
## ---- scenarios (e.g. one that cannot be fitted in them) is placed, and given its star, like for like: against the means of the
## ---- other models over the scenarios it has a result in.
##
fn_rank_like_for_like <- function(r,
                                  model,
                                  scenario,
                                  value,
                                  MCSE) {

      complete <-  r[!r$fewer_scenarios, ]
      complete <-  complete[order(complete$mean), ]
      complete$place <-  seq_len(nrow(complete))
      complete$mean_of_the_best <-  complete$mean[1]
      complete$critical_value <-  4 * pmax(complete$mean_MCSE, complete$mean_MCSE[1])
      others <-  r[r$fewer_scenarios, ]
      if (nrow(others) > 0) {
            others$place <-  NA
            others$mean_of_the_best <-  NA
            others$critical_value <-  NA
            for (k in seq_len(nrow(others))) {
                  same <-  scenario %in% scenario[model == others$model[k]]
                  means <-  sapply(complete$model, function(m) mean(value[model == m & same]))
                  MCSEs <-  sapply(complete$model, function(m) mean(MCSE[model == m & same]))
                  others$place[k] <-  sum(means < others$mean[k]) + 0.5
                  others$mean_of_the_best[k] <-  min(means)
                  others$critical_value[k] <-  4 * max(others$mean_MCSE[k], MCSEs[which.min(means)])
            }
      }
      r <-  rbind(complete, others)
      r <-  r[order(r$place, r$mean), ]
      r$rank <-  seq_len(nrow(r))
      r$difference_to_the_best <-  r$mean - r$mean_of_the_best
      r$star <-  r$difference_to_the_best > 0 & r$difference_to_the_best >= r$critical_value
      r$place <-  NULL
      r

}
##
## ---- The ranking of one measure within one set of scenarios (one test, or all):
##
fn_ranking <- function(cells,
                       measure,
                       scope) {

      models <-  unique(cells$model)
      r <-  data.frame(measure = measure, scope = scope, model = models, stringsAsFactors = FALSE)
      r$n_scenarios <-  sapply(models, function(m) sum(cells$model == m))
      r$mean <-  sapply(models, function(m) mean(cells$value[cells$model == m]))
      r$mean_of_the_other_measure <-  sapply(models, function(m) mean(cells$other[cells$model == m]))
      r$mean_MCSE <-  sapply(models, function(m) mean(cells$MCSE_actual[cells$model == m]))
      r$best_group <-  sapply(models, function(m) sum(cells$in_the_best_group[cells$model == m]))
      r$scenarios_of_the_scope <-  length(unique(cells$scenario))
      r$fewer_scenarios <-  r$n_scenarios < r$scenarios_of_the_scope
      fn_rank_like_for_like(r, cells$model, cells$scenario, cells$value, cells$MCSE_actual)

}
##
## ---- LaTeX of the table by test and of the overall table:
##
fn_number_with_star <- function(x,
                                star) {
      ifelse(star, sprintf("$%s^{*}$", fn_2(x)), sprintf("$%s$", fn_2(x)))
}
double_dagger <-  "$^{\\ddagger}$"
fn_table_by_test <- function(rankings,
                             measure,
                             by = "test") {

      name <-  if (measure == "RMSE") "RMSE" else "Bias"
      n <-  rankings[[1]]$scenarios_of_the_scope[1]
      columns <-  names(rankings)
      heads <-  if (by == "test") name_of_the_test[columns] else sprintf("$\\sim %s\\%%$ missing", columns)
      lines <-  c("\\begin{table}[H]",
                  "\\centering",
                  "\\caption{",
                  sprintf("\\textcolor{orange}{Ranking of the eight models by mean %s for each %s (%s)}", if (measure == "RMSE") "RMSE" else "bias",
                          if (by == "test") "test" else "proportion of missing thresholds", name_of_the_set[prior_set]),
                  "}",
                  "\\footnotesize",
                  "\\setlength{\\tabcolsep}{4pt}",
                  "\\begin{tabular}{c|lcc|lcc|lcc}",
                  "\\toprule",
                  paste0(" & ", paste(sprintf("\\multicolumn{3}{c%s}{\\textbf{%s}}", c("|", "|", ""), heads), collapse = " & "), " \\\\"),
                  paste0("\\textbf{Rank} & ", paste(rep(sprintf("\\textbf{Model} & \\textbf{%s} & \\textbf{Best}", name), 3), collapse = " & "), " \\\\"),
                  "\\midrule")
      any_fewer <-  any(sapply(rankings, function(r) any(r$fewer_scenarios)))
      for (k in seq_len(max(sapply(rankings, nrow)))) {
            cells <-  sapply(columns, function(t) { if (k > nrow(rankings[[t]])) return(" & & ")
                                                    r <-  rankings[[t]][k, ]
                                                    sprintf("%s%s & %s & $%d/%d$", name_of_the_model[r$model], if (r$fewer_scenarios) double_dagger else "",
                                                            fn_number_with_star(r$mean, r$star), r$best_group, r$n_scenarios) })
            lines <-  c(lines, paste0(sprintf("$%d$ & ", k), paste(cells, collapse = " & "), " \\\\"))
      }
      lines <-  c(lines,
                  "\\bottomrule",
                  "\\end{tabular}",
                  "\\begin{tablenotes}",
                  "\\footnotesize",
                  sprintf("\\item %s = mean, across the $%d$ scenarios of each %s ($3$ %s $\\times$ $7$ DGMs $\\times$ $N_{studies} \\in \\{10, 50\\}$), of %s (percentage points).",
                          name, n, if (by == "test") "test" else "proportion of missing thresholds", if (by == "test") "proportions of missing thresholds" else "tests",
                          if (measure == "RMSE") "RMSE(Se) + RMSE(Sp)" else "$|$bias(Se)$|$ + $|$bias(Sp)$|$"),
                  sprintf("\\item Best = number of scenarios in which the model was in \"Group 1 - Best\" of the 4-group tables (ordered by %s).", if (measure == "RMSE") "RMSE" else "bias"),
                  sprintf("\\item $^{*}$ = significantly worse than the best-ranked model of the same column, i.e., a difference in mean %s of at least $4 \\times$ MCSE", name),
                  "(the criterion for statistical significance of the 4-group tables, using the mean MCSE across the scenarios).",
                  if (any_fewer) sprintf("\\item %s = no result in some scenarios (the denominator of \"Best\" gives the number of scenarios with a result): the mean is across those scenarios, and the model is ranked - and compared with the best model - against the means of the other models across the same scenarios.", double_dagger),
                  "\\end{tablenotes}",
                  sprintf("\\label{Table:ranking_of_the_models_by_%s_%s}", if (by == "test") "test" else "share_missing", measure),
                  "\\end{table}")
      lines

}
fn_table_overall <- function(r,
                             measure) {

      name <-  if (measure == "RMSE") "RMSE" else "Bias"
      other_name <-  if (measure == "RMSE") "Bias" else "RMSE"
      lines <-  c("\\begin{table}[H]",
                  "\\centering",
                  "\\caption{",
                  sprintf("\\textcolor{orange}{Overall ranking of the eight models by mean %s across all three tests and all three proportions of missing thresholds (%s)}", if (measure == "RMSE") "RMSE" else "bias", name_of_the_set[prior_set]),
                  "}",
                  "\\footnotesize",
                  "\\setlength{\\tabcolsep}{4pt}",
                  "\\begin{tabular}{clllccc}",
                  "\\toprule",
                  sprintf("\\textbf{Rank} & \\textbf{Model} & \\textbf{Reference} & \\textbf{R package} & \\textbf{%s} & \\textbf{Best} & \\textbf{%s} \\\\", name, other_name),
                  "\\midrule")
      for (k in seq_len(nrow(r))) {
            lines <-  c(lines, sprintf("$%d$ & %s%s & %s & %s & %s & $%d/%d$ & $%s$ \\\\", k, name_of_the_model[r$model[k]], if (r$fewer_scenarios[k]) double_dagger else "",
                                       reference_of_the_model[r$model[k]], package_of_the_model[r$model[k]],
                                       fn_number_with_star(r$mean[k], r$star[k]), r$best_group[k], r$n_scenarios[k], fn_2(r$mean_of_the_other_measure[k])))
      }
      lines <-  c(lines,
                  "\\bottomrule",
                  "\\end{tabular}",
                  "\\begin{tablenotes}",
                  "\\footnotesize",
                  sprintf("\\item RMSE and Bias = mean, across the $%d$ scenarios ($3$ tests $\\times$ $3$ proportions of missing thresholds $\\times$ $7$ DGMs $\\times$ $N_{studies} \\in \\{10, 50\\}$), of RMSE(Se) + RMSE(Sp) and of $|$bias(Se)$|$ + $|$bias(Sp)$|$ (percentage points).",
                          r$scenarios_of_the_scope[1]),
                  sprintf("\\item Best = number of scenarios in which the model was in \"Group 1 - Best\" of the 4-group tables (ordered by %s).", if (measure == "RMSE") "RMSE" else "bias"),
                  sprintf("\\item $^{*}$ = significantly worse than the best-ranked model (%s), i.e., a difference in mean %s of at least $4 \\times$ MCSE", name_of_the_model[r$model[1]], name),
                  "(the criterion for statistical significance of the 4-group tables, using the mean MCSE across the scenarios).",
                  if (any(r$fewer_scenarios)) sprintf("\\item %s = no result in some scenarios (the denominator of \"Best\" gives the number of scenarios with a result): the mean is across those scenarios, and the model is ranked - and compared with the best model - against the means of the other models across the same scenarios.", double_dagger),
                  "\\end{tablenotes}",
                  sprintf("\\label{Table:ranking_of_the_models_overall_%s}", measure),
                  "\\end{table}")
      lines

}
##
## ---- Make and write:
##
every_number <-  NULL
for (measure in c("RMSE", "bias")) {
      cells <-  fn_cells_of_one_measure(d, measure)
      rankings <-  lapply(tests, function(t) fn_ranking(cells[cells$test == t, ], measure, t)); names(rankings) <-  tests
      by_share <-  lapply(shares, function(p) fn_ranking(cells[cells$pct == p, ], measure, paste0(p, "% missing"))); names(by_share) <-  as.character(shares)
      by_set <-  lapply(tests, function(t) lapply(shares, function(p) fn_ranking(cells[cells$test == t & cells$pct == p, ], measure, paste0(t, " ", p, "% missing"))))
      overall <-  fn_ranking(cells, measure, "overall")
      writeLines(fn_table_by_test(by_share, measure, by = "share"), file.path(folder, sprintf("ranking_of_the_models_by_share_missing_%s.tex", measure)))
      writeLines(fn_table_by_test(rankings, measure), file.path(folder, sprintf("ranking_of_the_models_by_test_%s.tex", measure)))
      writeLines(fn_table_overall(overall, measure), file.path(folder, sprintf("ranking_of_the_models_overall_%s.tex", measure)))
      every_number <-  rbind(every_number, do.call(rbind, by_share), do.call(rbind, rankings), do.call(rbind, unlist(by_set, recursive = FALSE)), overall)
}
rownames(every_number) <-  NULL
every_number$model <-  name_of_the_model[every_number$model]
utils::write.csv(every_number, file.path(folder, "every_number_of_the_ranking_tables.csv"), row.names = FALSE)
shown <-  every_number[, c("measure", "scope", "rank", "model", "mean", "star", "difference_to_the_best", "critical_value", "best_group", "n_scenarios", "mean_of_the_other_measure")]
for (column in c("mean", "difference_to_the_best", "critical_value", "mean_of_the_other_measure")) shown[[column]] <-  round(shown[[column]], 3)
print(shown, row.names = FALSE)
