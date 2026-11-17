##
## ---- The detailed results tables of the appendix (LaTeX), one for each of the nine sets of the factorial design (prior set 0) -------
##
## One long table a set (test x proportion of thresholds missing): for each DGM, number of studies and model, the RMSE and |bias| of
## Se and Sp (each with its MCSE), the coverage and the width of the 95% intervals, and the group of the model in the 4-group
## classification of that scenario by RMSE and by bias (the paper's rule: 1 = best; 2 = statistically different only; 3 = practically
## different only, i.e. at least 10% worse; 4 = both). The row of the model that generated the data is in bold; a dagger marks a cell
## whose MCSE is over the target. Needs the longtable package.
##
## Usage: Rscript <this file> <table of every cell of the nine sets, .rds> <folder to write in> <MCSE target> <allowance over it>
##
arguments <-  commandArgs(trailingOnly = TRUE)
d <-  readRDS(arguments[1])
folder <-  arguments[2]
MCSE_target <-  as.numeric(arguments[3])
allowance <-  as.numeric(arguments[4])
dir.create(folder, recursive = TRUE, showWarnings = FALSE)
relative_difference_threshold <-  10
d <-  d[d$prior_set == "set_0", ]
d$over_the_target <-  d$MCSE_val > MCSE_target + allowance
order_of_the_models <-  c("strat. bivariate", "Jones", "O-biv-FC", "O-biv-RC", "ord-HSROC", "Hoyer 2018", "logit GLMM", "diagmeta")
name_of_the_model <-  c("strat. bivariate" = "Strat-biv", "Jones" = "Jones", "O-biv-FC" = "O-biv-FC", "O-biv-RC" = "O-biv-RC", "ord-HSROC" = "O-HSROC-RC",
                        "Hoyer 2018" = "Hoyer", "logit GLMM" = "Logit-GLMM", "diagmeta" = "Steinhauser")
DGMs <-  c("Jones", "ord_biv_fixed_C", "ord_biv_random_C", "R&G_random_C", "Hoyer_2018", "logit_GLMM", "diagmeta")
name_of_the_DGM <-  c("Jones" = "Jones", "ord_biv_fixed_C" = "O-biv-FC", "ord_biv_random_C" = "O-biv-RC", "R&G_random_C" = "O-HSROC-RC",
                      "Hoyer_2018" = "Hoyer", "logit_GLMM" = "Logit-GLMM", "diagmeta" = "Steinhauser")
model_of_the_DGM <-  c("Jones" = "Jones", "ord_biv_fixed_C" = "O-biv-FC", "ord_biv_random_C" = "O-biv-RC", "R&G_random_C" = "ord-HSROC",
                       "Hoyer_2018" = "Hoyer 2018", "logit_GLMM" = "logit GLMM", "diagmeta" = "diagmeta")
name_of_the_test <-  c(GAD_2 = "GAD-2", HADS = "HADS", BAI = "BAI")
fn_2 <- function(x) sprintf("%.2f", x)
fn_1 <- function(x) sprintf("%.1f", x)
fn_3_figures <- function(x) if (x >= 10) sprintf("%.1f", x) else sprintf("%.2f", x)
dagger <-  "$^{\\dagger}$"
fn_groups <- function(value,
                      MCSE) {

      best <-  which.min(value)
      statistically <-  (value - value[best]) >= 4 * pmax(MCSE, MCSE[best])
      practically <-  100 * (value - value[best]) / value[best] >= relative_difference_threshold
      ifelse(!statistically & !practically, 1L, ifelse(statistically & !practically, 2L, ifelse(!statistically & practically, 3L, 4L)))

}
head_of_the_columns <-  c("\\toprule",
                          "\\multirow{2}{*}{\\textbf{DGM}} & \\multirow{2}{*}{\\textbf{$n$}} & \\multirow{2}{*}{\\textbf{Model}} & \\multicolumn{2}{c|}{\\textbf{RMSE (\\%)}} & \\multicolumn{2}{c|}{\\textbf{$|$Bias$|$ (\\%)}} & \\multicolumn{2}{c|}{\\textbf{Cvg. (\\%)}} & \\multicolumn{2}{c|}{\\textbf{Width}} & \\multicolumn{2}{c}{\\textbf{Group}} \\\\",
                          "\\cmidrule(lr){4-5} \\cmidrule(lr){6-7} \\cmidrule(lr){8-9} \\cmidrule(lr){10-11} \\cmidrule(lr){12-13}",
                          "& & & \\textbf{Se} & \\textbf{Sp} & \\textbf{Se} & \\textbf{Sp} & \\textbf{Se} & \\textbf{Sp} & \\textbf{Se} & \\textbf{Sp} & \\textbf{RMSE} & \\textbf{Bias} \\\\",
                          "\\midrule")
everything <-  character(0)
for (test in c("GAD_2", "HADS", "BAI")) for (pct in c(15, 55, 85)) {
      rows_of_the_set <-  d[d$test == test & d$pct == pct, ]
      body <-  character(0)
      for (DGM in DGMs) {
            for (n_studies in c(10, 50)) {
                  rows <-  rows_of_the_set[rows_of_the_set$DGM == DGM & rows_of_the_set$n_studies == n_studies, ]
                  rows <-  rows[order(match(rows$model, order_of_the_models)), ]
                  if (nrow(rows) == 0) next
                  group_RMSE <-  fn_groups(rows$avg_RMSE, pmax(rows$avg_MCSE_RMSE_Se, rows$avg_MCSE_RMSE_Sp, na.rm = TRUE))
                  group_bias <-  fn_groups(rows$avg_abs_bias, pmax(rows$avg_MCSE_bias_Se, rows$avg_MCSE_bias_Sp, na.rm = TRUE))
                  for (i in seq_len(nrow(rows))) {
                        r <-  rows[i, ]
                        values <-  c(fn_2(r$avg_RMSE_Se), fn_2(r$avg_RMSE_Sp), fn_2(r$avg_abs_bias_Se), fn_2(r$avg_abs_bias_Sp),
                                     fn_1(100 * r$avg_coverage_Se), fn_1(100 * r$avg_coverage_Sp), fn_3_figures(r$avg_width_Se), fn_3_figures(r$avg_width_Sp),
                                     as.character(group_RMSE[i]), as.character(group_bias[i]))
                        if (r$model == model_of_the_DGM[DGM]) values <-  sprintf("\\textbf{%s}", values)
                        MCSEs <-  c(fn_2(r$avg_MCSE_RMSE_Se), fn_2(r$avg_MCSE_RMSE_Sp), fn_2(r$avg_MCSE_bias_Se), fn_2(r$avg_MCSE_bias_Sp))
                        values[1:4] <-  sprintf("%s {\\scriptsize(%s)}", values[1:4], MCSEs)
                        body <-  c(body, sprintf("%s & %s & %s%s & %s \\\\", if (i == 1 && n_studies == 10) name_of_the_DGM[DGM] else "", if (i == 1) as.character(n_studies) else "",
                                                 name_of_the_model[r$model], if (r$over_the_target) dagger else "", paste(values, collapse = " & ")))
                  }
                  body <-  c(body, if (n_studies == 10) "\\cmidrule(lr){2-13}" else if (DGM != tail(DGMs, 1)) "\\midrule")
            }
      }
      notes <-  c("\\item $n$ = number of studies; RMSE and $|$Bias$|$ are given with their MCSE in brackets; Cvg. = coverage of the $95\\%$ intervals; Width = mean width of the $95\\%$ intervals; the row of the model that generated the data is in bold.",
                  "\\item Group = group of the model in the 4-group classification of the scenario (by RMSE and by bias): $1$ = best; $2$ = statistically different from the best only; $3$ = practically different only (at least $10\\%$ worse); $4$ = both.",
                  if (any(rows_of_the_set$over_the_target)) sprintf("\\item %s = MCSE of the RMSE above the target of $%s$ when the simulation study was stopped.", dagger, format(MCSE_target)),
                  if (test == "GAD_2" && pct == 85) "\\item The Steinhauser model could not be fitted at this proportion of missing thresholds (fewer than two thresholds reported per study on average, i.e., more random effects than observations in its most general variant).")
      label <-  sprintf("table:sim_results_%s_%d_pct", sub("GAD_2", "GAD2", test), pct)
      caption <-  sprintf("%s simulation results, $\\sim %d\\%%$ of thresholds missing, prior set 0: performance and 4-group classification of each model, by DGM and number of studies",
                          name_of_the_test[test], pct)
      table <-  c("{\\footnotesize", "\\setlength{\\tabcolsep}{2pt}", "\\begin{longtable}{lcl|cc|cc|cc|cc|cc}",
                  sprintf("\\caption{%s}", caption), sprintf("\\label{%s} \\\\", label), head_of_the_columns, "\\endfirsthead",
                  sprintf("\\caption[]{%s (continued)} \\\\", caption), head_of_the_columns, "\\endhead",
                  "\\midrule", "\\multicolumn{13}{r}{\\textit{continued on the next page}} \\\\", "\\endfoot", "\\bottomrule", "\\endlastfoot",
                  body, "\\end{longtable}", "\\begin{tablenotes}", "\\footnotesize", notes, "\\end{tablenotes}", "}")
      writeLines(table, file.path(folder, sprintf("detailed_table_%s_%d_pct.tex", test, pct)))
      everything <-  c(everything, "%%%%", sprintf("%%%% ==== %s, %d%% missing", test, pct), "%%%%", table, "")
      cat(test, pct, ": rows", sum(grepl("&", body)), "\n")
}
writeLines(everything, file.path(folder, "every_detailed_table_in_the_order_of_the_paper.tex"))
