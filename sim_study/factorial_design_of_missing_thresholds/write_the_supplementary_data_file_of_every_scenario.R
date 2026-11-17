##
## ---- The supplementary data file of the simulation study: one row for each model, prior set and scenario of the nine sets ----------
##
## Columns: prior set, test, proportion of thresholds missing (%), DGM, number of studies, model, simulations kept, RMSE, |bias|,
## coverage and interval width of Se and Sp (percentage points) with the MCSE of the RMSE and of the bias, and whether the MCSE
## target was met. Names of the models and DGMs as in the paper.
##
## Usage: Rscript <this file> <table of every cell of the nine sets, .rds> <csv to write>
##
arguments <-  commandArgs(trailingOnly = TRUE)
d <-  readRDS(arguments[1])
name_of_the_model <-  c("strat. bivariate" = "Strat-biv", "Jones" = "Jones", "O-biv-FC" = "O-biv-FC", "O-biv-RC" = "O-biv-RC", "ord-HSROC" = "O-HSROC-RC",
                        "Hoyer 2018" = "Hoyer", "logit GLMM" = "Logit-GLMM", "diagmeta" = "Steinhauser")
name_of_the_DGM <-  c("Jones" = "Jones", "ord_biv_fixed_C" = "O-biv-FC", "ord_biv_random_C" = "O-biv-RC", "R&G_random_C" = "O-HSROC-RC",
                      "Hoyer_2018" = "Hoyer", "logit_GLMM" = "Logit-GLMM", "diagmeta" = "Steinhauser")
out <-  data.frame(prior_set = sub("set_", "", d$prior_set), test = sub("GAD_2", "GAD-2", d$test), pct_thresholds_missing = d$pct, DGM = name_of_the_DGM[d$DGM],
                   n_studies = d$n_studies, model = name_of_the_model[d$model], n_sims = d$n_sims,
                   RMSE_Se = d$avg_RMSE_Se, RMSE_Sp = d$avg_RMSE_Sp, MCSE_RMSE_Se = d$avg_MCSE_RMSE_Se, MCSE_RMSE_Sp = d$avg_MCSE_RMSE_Sp,
                   abs_bias_Se = d$avg_abs_bias_Se, abs_bias_Sp = d$avg_abs_bias_Sp, MCSE_bias_Se = d$avg_MCSE_bias_Se, MCSE_bias_Sp = d$avg_MCSE_bias_Sp,
                   coverage_Se = 100 * d$avg_coverage_Se, coverage_Sp = 100 * d$avg_coverage_Sp, width_Se = d$avg_width_Se, width_Sp = d$avg_width_Sp,
                   MCSE_target_met = d$MCSE_val <= 0.125 + 0.0005, stringsAsFactors = FALSE)
stopifnot(!anyNA(out$model), !anyNA(out$DGM))
out <-  out[order(out$prior_set, match(out$test, c("GAD-2", "HADS", "BAI")), out$pct_thresholds_missing, out$DGM, out$n_studies, out$model), ]
numeric_columns <-  sapply(out, is.numeric) & !names(out) %in% c("pct_thresholds_missing", "n_studies", "n_sims")
out[numeric_columns] <-  lapply(out[numeric_columns], function(x) round(x, 4))
utils::write.csv(out, arguments[2], row.names = FALSE)
cat("rows written:", nrow(out), "; MCSE target not met:", sum(!out$MCSE_target_met), "\n")
