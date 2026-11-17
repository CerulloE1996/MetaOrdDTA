##
## ---- The results tables of the main manuscript (LaTeX) from the table of every cell of the nine sets of the factorial design -------
##
##   (1) overall ranking of the eight models by mean RMSE, with the mean bias beside it (all scenarios, prior set 0);
##   (2) ranking by mean RMSE for each proportion of thresholds missing;   (3) ranking by mean RMSE for each test;
##   (4) mean RMSE of each model in each of the nine sets (3 tests x 3 proportions missing);
##   (5) mean coverage and mean interval width of Se and Sp of each model, for each proportion missing;
##   (6) mean RMSE and mean bias of the five Bayesian models under prior sets 0, 1 and 2.
## Rules, as in the paper's ranking tables: RMSE = RMSE(Se) + RMSE(Sp) and bias = |bias(Se)| + |bias(Sp)| over the thresholds of
## interest, averaged over the scenarios; MCSE of a cell = max(MCSE of Se, MCSE of Sp), of a model = its mean over the scenarios;
## a star when the mean is worse than the best of the same column by at least 4 x max(MCSE of the model, MCSE of the best);
## "Best" = the number of scenarios in which the model is in "Group 1 - Best" of the 4-group rule (relative difference under 10%
## and not statistically different). Every scored cell counts. A model's mean is over the scenarios it has a result in; a double
## dagger marks a model with fewer scenarios than the column has.
## Written: one .tex file per table and every number as .csv.
##
## Usage: Rscript <this file> <table of every cell of the nine sets, .rds> <folder to write in> <fits kept and lost of the nine sets, .json>
##
arguments <-  commandArgs(trailingOnly = TRUE)
d <-  readRDS(arguments[1])
folder <-  arguments[2]
fits <-  jsonlite::fromJSON(arguments[3])
dir.create(folder, recursive = TRUE, showWarnings = FALSE)
relative_difference_threshold <-  10
d$scenario <-  paste(d$test, d$pct, d$DGM, d$n_studies)
d$RMSE <-  d$avg_RMSE
d$bias <-  d$avg_abs_bias
d$MCSE_RMSE <-  pmax(d$avg_MCSE_RMSE_Se, d$avg_MCSE_RMSE_Sp, na.rm = TRUE)
d$MCSE_bias <-  pmax(d$avg_MCSE_bias_Se, d$avg_MCSE_bias_Sp, na.rm = TRUE)
frequentist <-  d$model %in% c("Hoyer 2018", "logit GLMM", "diagmeta")
set_0 <-  d[d$prior_set == "set_0", ]
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
tests <-  c("GAD_2", "HADS", "BAI")
name_of_the_test <-  c(GAD_2 = "GAD-2", HADS = "HADS", BAI = "BAI")
shares <-  c(15, 55, 85)
fn_2 <- function(x) sprintf("%.2f", x)
fn_1 <- function(x) sprintf("%.1f", x)
double_dagger <-  "$^{\\ddagger}$"
fn_number <- function(x, star) ifelse(star, sprintf("$%s^{*}$", fn_2(x)), sprintf("$%s$", fn_2(x)))
##
## ---- Whether each cell is in "Group 1 - Best" of its scenario, for one measure:
##
fn_best_group <- function(cells,
                          measure) {

      value <-  cells[[measure]]
      MCSE <-  cells[[paste0("MCSE_", measure)]]
      in_the_best_group <-  rep(NA, nrow(cells))
      for (s in unique(cells$scenario)) {
            i <-  which(cells$scenario == s)
            best <-  i[which.min(value[i])]
            statistically <-  (value[i] - value[best]) >= 4 * pmax(MCSE[i], MCSE[best])
            practically <-  100 * (value[i] - value[best]) / value[best] >= relative_difference_threshold
            in_the_best_group[i] <-  !statistically & !practically
      }
      in_the_best_group

}
set_0$best_RMSE <-  fn_best_group(set_0, "RMSE")
set_0$best_bias <-  fn_best_group(set_0, "bias")
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
## ---- The ranking of one measure within one set of scenarios:
##
fn_ranking <- function(cells,
                       measure,
                       scope) {

      models <-  unique(cells$model)
      r <-  data.frame(measure = measure, scope = scope, model = models, stringsAsFactors = FALSE)
      r$n_scenarios <-  sapply(models, function(m) sum(cells$model == m))
      r$mean <-  sapply(models, function(m) mean(cells[[measure]][cells$model == m]))
      r$mean_MCSE <-  sapply(models, function(m) mean(cells[[paste0("MCSE_", measure)]][cells$model == m]))
      r$best_group <-  sapply(models, function(m) sum(cells[[paste0("best_", measure)]][cells$model == m]))
      r$mean_bias <-  sapply(models, function(m) mean(cells$bias[cells$model == m]))
      r$scenarios_of_the_scope <-  length(unique(cells$scenario))
      r$fewer_scenarios <-  r$n_scenarios < r$scenarios_of_the_scope
      fn_rank_like_for_like(r, cells$model, cells$scenario, cells[[measure]], cells[[paste0("MCSE_", measure)]])

}
every_number <-  NULL
note_on_the_star <-  c("\\item $^{*}$ = significantly worse than the best-ranked model of the same column, i.e., a difference in means of at least $4 \\times$ MCSE",
                       "(the mean MCSE across the scenarios; the criterion for statistical significance of the 4-group classification).")
note_on_fewer <-  sprintf("\\item %s = no result in some scenarios (the denominator of \"Best\" gives the number of scenarios with a result): the mean is across those scenarios, and the model is ranked - and compared with the best model - against the means of the other models across the same scenarios.", double_dagger)
##
## ---- (1) the overall ranking:
##
set_0$family <-  ifelse(set_0$DGM == "Jones", "Jones", ifelse(set_0$DGM %in% c("ord_biv_fixed_C", "ord_biv_random_C", "R&G_random_C"), "ordinal", "competing"))
families <-  c("Jones", "ordinal", "competing")
## (mean RMSE of a model under a family of DGMs; a double dagger when it has no result in some scenarios of the family)
fn_family <- function(m, f) {
      rows <-  set_0$family == f
      sprintf("$%s%s$", fn_2(mean(set_0$RMSE[rows & set_0$model == m])), if (sum(rows & set_0$model == m) < length(unique(set_0$scenario[rows]))) "^{\\ddagger}" else "")
}
overall_RMSE <-  fn_ranking(set_0, "RMSE", "overall")
overall_bias <-  fn_ranking(set_0, "bias", "overall")
every_number <-  rbind(every_number, overall_RMSE, overall_bias)
lines <-  c("\\begin{table}[H]", "\\centering", "\\caption{",
            "\\textcolor{orange}{Overall ranking of the eight models by mean RMSE - with the mean bias - across all three tests and all three proportions of missing thresholds (prior set 0)}",
            "}", "\\footnotesize", "\\setlength{\\tabcolsep}{3pt}", "\\begin{tabular}{clllcccc|ccc}", "\\toprule",
            " & & & & & & & & \\multicolumn{3}{c}{\\textbf{RMSE by DGM}} \\\\",
            "\\textbf{Rank} & \\textbf{Model} & \\textbf{Reference} & \\textbf{R package} & \\textbf{RMSE} & \\textbf{Best} & \\textbf{Bias} & \\textbf{Best} & \\textbf{Jones} & \\textbf{Ordinal} & \\textbf{Competing} \\\\", "\\midrule")
for (k in seq_len(nrow(overall_RMSE))) {
      r <-  overall_RMSE[k, ]
      b <-  overall_bias[overall_bias$model == r$model, ]
      lines <-  c(lines, sprintf("$%d$ & %s%s & %s & %s & %s & $%d/%d$ & %s & $%d/%d$ & %s \\\\", k, name_of_the_model[r$model], if (r$fewer_scenarios) double_dagger else "",
                                 reference_of_the_model[r$model], package_of_the_model[r$model], fn_number(r$mean, r$star), r$best_group, r$n_scenarios,
                                 fn_number(b$mean, b$star), b$best_group, b$n_scenarios, paste(sapply(families, function(f) fn_family(r$model, f)), collapse = " & ")))
}
lines <-  c(lines, "\\bottomrule", "\\end{tabular}", "\\begin{tablenotes}", "\\footnotesize",
            sprintf(paste0("\\item RMSE and Bias = mean, across the $%d$ scenarios ($3$ tests $\\times$ $3$ proportions of missing thresholds $\\times$ $7$ DGMs $\\times$ $N_{studies} \\in \\{10, 50\\}$), ",
                           "of RMSE(Se) + RMSE(Sp) and of $|$bias(Se)$|$ + $|$bias(Sp)$|$ (percentage points)."), overall_RMSE$scenarios_of_the_scope[1]),
            "\\item Best = number of scenarios in which the model was in \"Group 1 - Best\" of the 4-group classification (for RMSE and for bias, respectively).",
            sprintf("\\item RMSE by DGM = mean RMSE across the $%d$ scenarios of the Jones DGM, the $%d$ scenarios of the three ordinal DGMs (O-biv-FC, O-biv-RC and O-HSROC-RC) and the $%d$ scenarios of the three competing DGMs (Hoyer, Logit-GLMM and Steinhauser).",
                    length(unique(set_0$scenario[set_0$family == "Jones"])), length(unique(set_0$scenario[set_0$family == "ordinal"])), length(unique(set_0$scenario[set_0$family == "competing"]))),
            note_on_the_star, if (any(overall_RMSE$fewer_scenarios)) note_on_fewer, "\\end{tablenotes}", "\\label{Table:ranking_of_the_models_overall}", "\\end{table}")
writeLines(lines, file.path(folder, "main_table_overall_ranking.tex"))
##
## ---- (2), (3) the ranking by proportion missing and by test (RMSE):
##
fn_table_by <- function(rankings,
                        heads,
                        what,
                        other,
                        label) {

      n <-  rankings[[1]]$scenarios_of_the_scope[1]
      any_fewer <-  any(sapply(rankings, function(r) any(r$fewer_scenarios)))
      lines <-  c("\\begin{table}[H]", "\\centering", "\\caption{",
                  sprintf("\\textcolor{orange}{Ranking of the eight models by mean RMSE for each %s (prior set 0)}", what),
                  "}", "\\footnotesize", "\\setlength{\\tabcolsep}{3pt}", "\\begin{tabular}{c|lccc|lccc|lccc}", "\\toprule",
                  paste0(" & ", paste(sprintf("\\multicolumn{4}{c%s}{\\textbf{%s}}", c("|", "|", ""), heads), collapse = " & "), " \\\\"),
                  paste0("\\textbf{Rank} & ", paste(rep("\\textbf{Model} & \\textbf{RMSE} & \\textbf{Best} & \\textbf{Bias}", 3), collapse = " & "), " \\\\"), "\\midrule")
      for (k in seq_len(max(sapply(rankings, nrow)))) {
            cells <-  sapply(rankings, function(r) { if (k > nrow(r)) return(" & & & ")
                                                     sprintf("%s%s & %s & $%d/%d$ & $%s$", name_of_the_model[r$model[k]], if (r$fewer_scenarios[k]) double_dagger else "",
                                                             fn_number(r$mean[k], r$star[k]), r$best_group[k], r$n_scenarios[k], fn_2(r$mean_bias[k])) })
            lines <-  c(lines, paste0(sprintf("$%d$ & ", k), paste(cells, collapse = " & "), " \\\\"))
      }
      c(lines, "\\bottomrule", "\\end{tabular}", "\\begin{tablenotes}", "\\footnotesize",
        sprintf("\\item RMSE = mean, across the $%d$ scenarios of each %s ($3$ %s $\\times$ $7$ DGMs $\\times$ $N_{studies} \\in \\{10, 50\\}$), of RMSE(Se) + RMSE(Sp) (percentage points).", n, what, other),
        "\\item Best = number of scenarios in which the model was in \"Group 1 - Best\" of the 4-group classification (ordered by RMSE).",
        "\\item Bias = mean, across the same scenarios, of $|$bias(Se)$|$ + $|$bias(Sp)$|$ (percentage points); the rankings by bias are in appendix \\ref{appendix_E_rankings_by_bias}.",
        note_on_the_star, if (any_fewer) note_on_fewer, "\\end{tablenotes}", sprintf("\\label{%s}", label), "\\end{table}")

}
by_share <-  lapply(shares, function(p) fn_ranking(set_0[set_0$pct == p, ], "RMSE", paste0(p, "% missing")))
by_test <-  lapply(tests, function(t) fn_ranking(set_0[set_0$test == t, ], "RMSE", t))
every_number <-  rbind(every_number, do.call(rbind, by_share), do.call(rbind, by_test))
writeLines(fn_table_by(by_share, sprintf("$\\sim %d\\%%$ missing", shares), "proportion of missing thresholds", "tests", "Table:ranking_of_the_models_by_share_missing_RMSE"),
           file.path(folder, "main_table_ranking_by_share_missing.tex"))
writeLines(fn_table_by(by_test, name_of_the_test[tests], "test", "proportions of missing thresholds", "Table:ranking_of_the_models_by_test_RMSE"),
           file.path(folder, "main_table_ranking_by_test.tex"))
##
## ---- (4) the mean RMSE of each model in each of the nine sets:
##
by_set <-  list()
for (t in tests) for (p in shares) by_set[[paste(t, p)]] <-  fn_ranking(set_0[set_0$test == t & set_0$pct == p, ], "RMSE", paste0(t, " ", p, "% missing"))
every_number <-  rbind(every_number, do.call(rbind, by_set))
order_of_the_rows <-  overall_RMSE$model
lines <-  c("\\begin{table}[H]", "\\centering", "\\caption{",
            "\\textcolor{orange}{Mean RMSE of the eight models for each test and each proportion of missing thresholds (prior set 0)}",
            "}", "\\footnotesize", "\\setlength{\\tabcolsep}{4pt}", "\\begin{tabular}{l|ccc|ccc|ccc}", "\\toprule",
            paste0(" & ", paste(sprintf("\\multicolumn{3}{c%s}{\\textbf{%s}}", c("|", "|", ""), name_of_the_test[tests]), collapse = " & "), " \\\\"),
            paste0("\\textbf{Model} & ", paste(rep(paste(sprintf("$\\sim %d\\%%$", shares), collapse = " & "), 3), collapse = " & "), " \\\\"), "\\midrule")
any_fewer <-  FALSE
any_dash <-  FALSE
for (m in order_of_the_rows) {
      cells <-  sapply(by_set, function(r) { i <-  which(r$model == m)
                                             if (!length(i)) return("--")
                                             text <-  fn_2(r$mean[i])
                                             if (r$rank[i] == 1) text <-  sprintf("\\mathbf{%s}", text)
                                             sprintf("$%s%s$%s", text, if (r$star[i]) "^{*}" else "", if (r$fewer_scenarios[i]) double_dagger else "") })
      any_fewer <-  any_fewer || any(grepl("ddagger", cells))
      any_dash <-  any_dash || any(cells == "--")
      lines <-  c(lines, paste0(name_of_the_model[m], " & ", paste(cells, collapse = " & "), " \\\\"))
}
lines <-  c(lines, "\\bottomrule", "\\end{tabular}", "\\begin{tablenotes}", "\\footnotesize",
            "\\item Each value = mean, across the $14$ scenarios of the column ($7$ DGMs $\\times$ $N_{studies} \\in \\{10, 50\\}$), of RMSE(Se) + RMSE(Sp) (percentage points); bold = the best model of the column.",
            "\\item Column heads = proportion of missing thresholds; models ordered as in table \\ref{Table:ranking_of_the_models_overall}.",
            note_on_the_star,
            if (any_fewer) sprintf("\\item %s = no result in some scenarios of the column (mean across the scenarios with a result).", double_dagger),
            if (any_dash) "\\item -- = the model could not be fitted (fewer than two thresholds reported per study on average; table \\ref{table:Sim_study_the_three_tests}).",
            "\\end{tablenotes}", "\\label{Table:RMSE_of_the_models_in_the_nine_sets}", "\\end{table}")
writeLines(lines, file.path(folder, "main_table_RMSE_in_the_nine_sets.tex"))
##
## ---- (5) coverage and interval width, for each proportion missing:
##
coverage <-  NULL
lines <-  c("\\begin{table}[H]", "\\centering", "\\caption{",
            "\\textcolor{orange}{Mean coverage ($\\%$) and mean width of the $95\\%$ intervals of the eight models for each proportion of missing thresholds (prior set 0)}",
            "}", "\\footnotesize", "\\setlength{\\tabcolsep}{3pt}", "\\begin{tabular}{l|cccc|cccc|cccc}", "\\toprule",
            paste0(" & ", paste(sprintf("\\multicolumn{4}{c%s}{\\textbf{$\\sim %d\\%%$ missing}}", c("|", "|", ""), shares), collapse = " & "), " \\\\"),
            paste0(" & ", paste(rep("\\multicolumn{2}{c}{\\textbf{Cvg.}} & \\multicolumn{2}{c}{\\textbf{Width}}", 3), collapse = " & "), " \\\\"),
            paste0("\\textbf{Model} & ", paste(rep("\\textbf{Se} & \\textbf{Sp} & \\textbf{Se} & \\textbf{Sp}", 3), collapse = " & "), " \\\\"), "\\midrule")
for (m in order_of_the_rows) {
      cells <-  character(0)
      for (p in shares) {
            rows <-  set_0[set_0$model == m & set_0$pct == p, ]
            values <-  c(100 * mean(rows$avg_coverage_Se), 100 * mean(rows$avg_coverage_Sp), mean(rows$avg_width_Se), mean(rows$avg_width_Sp))
            coverage <-  rbind(coverage, data.frame(model = name_of_the_model[m], pct = p, n_scenarios = nrow(rows), coverage_Se = values[1], coverage_Sp = values[2], width_Se = values[3], width_Sp = values[4]))
            cells <-  c(cells, sprintf("$%s$", c(fn_1(values[1]), fn_1(values[2]), fn_2(values[3]), fn_2(values[4]))))
      }
      lines <-  c(lines, paste0(name_of_the_model[m], " & ", paste(cells, collapse = " & "), " \\\\"))
}
lines <-  c(lines, "\\bottomrule", "\\end{tabular}", "\\begin{tablenotes}", "\\footnotesize",
            "\\item Cvg. = mean, across the scenarios of each proportion of missing thresholds ($3$ tests $\\times$ $7$ DGMs $\\times$ $N_{studies} \\in \\{10, 50\\}$), of the coverage of the $95\\%$ intervals (nominal: $95\\%$).",
            "\\item Width = mean width of the $95\\%$ intervals (percentage points); models ordered as in table \\ref{Table:ranking_of_the_models_overall}.",
            "\\end{tablenotes}", "\\label{Table:coverage_and_width_by_share_missing}", "\\end{table}")
writeLines(lines, file.path(folder, "main_table_coverage_and_width.tex"))
##
## ---- (6) the five Bayesian models under the three prior sets:
##
bayesian <-  d[!frequentist, ]
by_prior_set <-  NULL
for (k in c("set_0", "set_1", "set_2")) for (measure in c("RMSE", "bias")) {
      cells <-  bayesian[bayesian$prior_set == k, ]
      cells$best_RMSE <-  NA; cells$best_bias <-  NA
      cells[[paste0("best_", measure)]] <-  fn_best_group(cells, measure)
      r <-  fn_ranking(cells, measure, paste("Bayesian models,", k))
      r$prior_set <-  k
      by_prior_set <-  rbind(by_prior_set, r)
}
rows_of_the_table <-  by_prior_set$model[by_prior_set$prior_set == "set_0" & by_prior_set$measure == "RMSE"]
lines <-  c("\\begin{table}[H]", "\\centering", "\\caption{",
            "\\textcolor{orange}{Mean RMSE and mean bias of the five Bayesian models under the three prior sets (all tests and proportions of missing thresholds)}",
            "}", "\\footnotesize", "\\setlength{\\tabcolsep}{5pt}", "\\begin{tabular}{l|ccc|ccc|cc}", "\\toprule",
            " & \\multicolumn{3}{c|}{\\textbf{RMSE}} & \\multicolumn{3}{c|}{\\textbf{Bias}} & \\multicolumn{2}{c}{\\textbf{Fits (\\%)}} \\\\",
            "\\textbf{Model} & \\textbf{Set 0} & \\textbf{Set 1} & \\textbf{Set 2} & \\textbf{Set 0} & \\textbf{Set 1} & \\textbf{Set 2} & \\textbf{Kept} & \\textbf{Div.} \\\\", "\\midrule")
any_fewer <-  FALSE
for (m in rows_of_the_table) {
      cells <-  character(0)
      for (measure in c("RMSE", "bias")) for (k in c("set_0", "set_1", "set_2")) {
            r <-  by_prior_set[by_prior_set$model == m & by_prior_set$measure == measure & by_prior_set$prior_set == k, ]
            cells <-  c(cells, paste0(fn_number(r$mean, r$star), if (r$fewer_scenarios) double_dagger else ""))
            any_fewer <-  any_fewer || r$fewer_scenarios
      }
      cells <-  c(cells, sprintf("$%s$", fn_1(fits[[paste0("pct_kept_", name_of_the_model[m])]])), sprintf("$%s$", fn_1(fits[[paste0("pct_lost_to_divergences_", name_of_the_model[m])]])))
      lines <-  c(lines, paste0(name_of_the_model[m], " & ", paste(cells, collapse = " & "), " \\\\"))
}
over <-  d[d$MCSE_val > 0.125 + 0.0005, ]
cell_of_the_largest <-  regmatches(fits$cell_with_the_largest_pct_lost, regexec("^(\\S+), DGM (\\S+), (\\S+) at (\\d+)% missing, (\\d+) studies, set_(\\d)$", fits$cell_with_the_largest_pct_lost))[[1]]
name_of_the_DGM <-  c("Jones" = "Jones", "ord_biv_fixed_C" = "O-biv-FC", "ord_biv_random_C" = "O-biv-RC", "ord_HSROC_random_C" = "O-HSROC-RC", "R&G_random_C" = "O-HSROC-RC",
                      "Hoyer_2018" = "Hoyer", "logit_GLMM" = "Logit-GLMM", "diagmeta" = "Steinhauser")
lines <-  c(lines, "\\bottomrule", "\\end{tabular}", "\\begin{tablenotes}", "\\footnotesize",
            sprintf(paste0("\\item Fits: Kept and Div. = percentage of the fits of the model (all scenarios and prior sets) which were kept, and which were discarded due to divergent transitions (more than $1\\%%$ of the transitions after warm-up); ",
                           "all five models: $%s$ fits, $%s\\%%$ discarded due to divergent transitions, $%s\\%%$ because the minimum $ESS$ was below the minimum required, and $%s\\%%$ kept; ",
                           "largest proportion discarded due to divergent transitions in one cell: $%s\\%%$ (%s model, %s DGM, %s, $\\sim %s\\%%$ missing thresholds, $N_{studies} = %s$, prior set $%s$)."),
                    format(fits$fits, big.mark = "{,}"), fn_1(fits$pct_lost_to_divergences), fn_1(fits$pct_under_the_min_ESS), fn_1(fits$pct_kept), fn_1(fits$largest_pct_lost_to_divergences_in_a_cell),
                    cell_of_the_largest[2], name_of_the_DGM[cell_of_the_largest[3]], sub("GAD_2", "GAD-2", cell_of_the_largest[4]), cell_of_the_largest[5], cell_of_the_largest[6], cell_of_the_largest[7]),
            sprintf("\\item MCSE target not reached when the simulation study was stopped: $%d$ of the $%s$ cells of the eight models (%s), all with $N_{studies} = 10$.",
                    nrow(over), format(nrow(d), big.mark = "{,}"), paste(sprintf("$%d$ %s", as.integer(table(over$model)), name_of_the_model[names(table(over$model))]), collapse = " and ")),
            sprintf("\\item RMSE and Bias = mean, across the $%d$ scenarios, of RMSE(Se) + RMSE(Sp) and of $|$bias(Se)$|$ + $|$bias(Sp)$|$ (percentage points); the prior sets are given in table \\ref{table:Sim_study_Se_Sp_prior_of_each_model_by_prior_set}.",
                    max(by_prior_set$scenarios_of_the_scope)),
            "\\item $^{*}$ = significantly worse than the best of the five models under the same prior set (a difference in means of at least $4 \\times$ MCSE).",
            if (any(over$n_studies != 10)) stop("a cell over the MCSE target with 50 studies: the note says all have 10"),
            if (any_fewer) sprintf("\\item %s = no result in some scenarios (mean across the scenarios with a result).", double_dagger),
            "\\end{tablenotes}", "\\label{Table:RMSE_and_bias_by_prior_set}", "\\end{table}")
writeLines(lines, file.path(folder, "main_table_prior_sets.tex"))
every_number <-  rbind(cbind(every_number, prior_set = "set_0"), by_prior_set[, c(names(every_number), "prior_set")])
every_number$model <-  name_of_the_model[every_number$model]
utils::write.csv(every_number, file.path(folder, "every_number_of_the_main_text_tables.csv"), row.names = FALSE)
utils::write.csv(coverage, file.path(folder, "every_number_of_the_coverage_table.csv"), row.names = FALSE)
cat("tables written:", length(list.files(folder, pattern = "^main_table_.*[.]tex$")), "\n")
print(overall_RMSE[, c("rank", "model", "mean", "star", "best_group", "n_scenarios")], row.names = FALSE, digits = 4)
