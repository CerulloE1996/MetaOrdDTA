##
## ---- Every number that the text of the paper quotes from the simulation study, as one list of key -> text (.json) --------------
##
## Means of RMSE (RMSE(Se) + RMSE(Sp)), bias (|bias(Se)| + |bias(Sp)|), coverage and interval width by model within groups of
## scenarios (all; each proportion of thresholds missing; each test; each of the nine sets; each family of DGM; each number of
## studies), under prior set 0, with the means of the five Bayesian models under prior sets 1 and 2, the order of the models
## (ranked like for like, as in the ranking tables) and the counts the text gives. Keys: <measure>.<group>.<level>.<model>.
##
## Usage: Rscript <this file> <table of every cell of the nine sets, .rds> <csv of every number of the main text tables>
##        <json of the fits kept and lost> <json to write>
##
arguments <-  commandArgs(trailingOnly = TRUE)
d <-  readRDS(arguments[1])
rankings <-  utils::read.csv(arguments[2], stringsAsFactors = FALSE)
fits <-  jsonlite::fromJSON(arguments[3])
out <-  arguments[4]
id_of_the_model <-  c("strat. bivariate" = "Strat", "Jones" = "Jones", "O-biv-FC" = "FC", "O-biv-RC" = "RC", "ord-HSROC" = "HSROC",
                      "Hoyer 2018" = "Hoyer", "logit GLMM" = "GLMM", "diagmeta" = "Stein")
name_in_the_text <-  c(Strat = "stratified-bivariate", Jones = "Jones", FC = "O-biv-FC", RC = "O-biv-RC", HSROC = "O-HSROC-RC", Hoyer = "Hoyer", GLMM = "logit GLMM", Stein = "Steinhauser")
id_from_the_table_name <-  c("Strat-biv" = "Strat", "Jones" = "Jones", "O-biv-FC" = "FC", "O-biv-RC" = "RC", "O-HSROC-RC" = "HSROC", "Hoyer" = "Hoyer", "Logit-GLMM" = "GLMM", "Steinhauser" = "Stein")
d$id <-  id_of_the_model[d$model]
d$family <-  ifelse(d$DGM == "Jones", "Jones", ifelse(d$DGM %in% c("ord_biv_fixed_C", "ord_biv_random_C", "R&G_random_C"), "ordinal", "competing"))
d$scenario <-  paste(d$test, d$pct, d$DGM, d$n_studies)
d$MCSE_RMSE <-  pmax(d$avg_MCSE_RMSE_Se, d$avg_MCSE_RMSE_Sp, na.rm = TRUE)
d$MCSE_bias <-  pmax(d$avg_MCSE_bias_Se, d$avg_MCSE_bias_Sp, na.rm = TRUE)
s0 <-  d[d$prior_set == "set_0", ]
numbers <-  list()
fn_2 <- function(x) sprintf("%.2f", x)
fn_1 <- function(x) sprintf("%.1f", x)
fn_0 <- function(x) sprintf("%.0f", x)
fn_list <- function(x) if (length(x) == 1) x else paste0(paste(x[-length(x)], collapse = ", "), " and ", x[length(x)])
measures <-  list(RMSE = list(s0$avg_RMSE, fn_2), bias = list(s0$avg_abs_bias, fn_2), cvgSe = list(100 * s0$avg_coverage_Se, fn_1), cvgSp = list(100 * s0$avg_coverage_Sp, fn_1),
                  widthSe = list(s0$avg_width_Se, fn_2), widthSp = list(s0$avg_width_Sp, fn_2))
groups <-  list(all = rep("all", nrow(s0)), pct = s0$pct, test = s0$test, set = paste0(s0$test, "_", s0$pct), family = s0$family, familypct = paste0(s0$family, "_", s0$pct),
                n = s0$n_studies, npct = paste0(s0$n_studies, "_", s0$pct), DGM = s0$DGM, testn = paste0(s0$test, "_", s0$n_studies))
means <-  list()
for (measure in names(measures)) for (group in names(groups)) {
      t <-  tapply(measures[[measure]][[1]], list(s0$id, groups[[group]]), mean)
      means[[paste(measure, group)]] <-  t
      for (m in rownames(t)) for (level in colnames(t)) if (!is.na(t[m, level])) numbers[[paste(measure, group, level, m, sep = ".")]] <-  measures[[measure]][[2]](t[m, level])
}
##
## ---- The order of the models, from the rankings of the main text tables (like for like), and which are starred:
##
rankings$id <-  id_from_the_table_name[rankings$model]
fn_order <- function(measure, scope, prior_set = "set_0") {
      r <-  rankings[rankings$measure == measure & rankings$scope == scope & rankings$prior_set == prior_set, ]
      r[order(r$rank), ]
}
scopes <-  c(all = "overall", pct15 = "15% missing", pct55 = "55% missing", pct85 = "85% missing", GAD_2 = "GAD_2", HADS = "HADS", BAI = "BAI")
for (k in names(scopes)) {
      r <-  fn_order("RMSE", scopes[k])
      numbers[[paste0("order.RMSE.", k)]] <-  fn_list(name_in_the_text[r$id])
      numbers[[paste0("ids.RMSE.", k)]] <-  paste(r$id, collapse = " ")
      numbers[[paste0("starred.RMSE.", k)]] <-  paste(r$id[r$star], collapse = " ")
      numbers[[paste0("n_starred.RMSE.", k)]] <-  as.character(sum(r$star))
      for (i in seq_len(nrow(r))) numbers[[paste0("rank.RMSE.", k, ".", r$id[i])]] <-  as.character(r$rank[i])
      for (i in seq_len(nrow(r))) numbers[[paste0("best.RMSE.", k, ".", r$id[i])]] <-  as.character(r$best_group[i])
      for (i in seq_len(nrow(r))) numbers[[paste0("n.RMSE.", k, ".", r$id[i])]] <-  as.character(r$n_scenarios[i])
}
r <-  fn_order("bias", "overall")
numbers[["order.bias.all"]] <-  fn_list(name_in_the_text[r$id])
numbers[["ids.bias.all"]] <-  paste(r$id, collapse = " ")
numbers[["starred.bias.all"]] <-  paste(r$id[r$star], collapse = " ")
numbers[["n_starred.bias.all"]] <-  as.character(sum(r$star))
for (i in seq_len(nrow(r))) numbers[[paste0("best.bias.all.", r$id[i])]] <-  as.character(r$best_group[i])
numbers[["n_scenarios"]] <-  as.character(length(unique(s0$scenario)))
numbers[["n_scenarios_Stein"]] <-  as.character(sum(s0$id == "Stein"))
##
## ---- The two O-bivariate models against one another, scenario by scenario (not different = under 4 x the larger MCSE):
##
fc <-  s0[s0$id == "FC", ]; rc <-  s0[s0$id == "RC", ]
both <-  merge(fc[, c("scenario", "avg_RMSE", "avg_abs_bias", "MCSE_RMSE", "MCSE_bias")], rc[, c("scenario", "avg_RMSE", "avg_abs_bias", "MCSE_RMSE", "MCSE_bias")], by = "scenario")
numbers[["FC_vs_RC.n"]] <-  as.character(nrow(both))
numbers[["FC_vs_RC.not_different.RMSE"]] <-  as.character(sum(abs(both$avg_RMSE.x - both$avg_RMSE.y) < 4 * pmax(both$MCSE_RMSE.x, both$MCSE_RMSE.y)))
numbers[["FC_vs_RC.not_different.bias"]] <-  as.character(sum(abs(both$avg_abs_bias.x - both$avg_abs_bias.y) < 4 * pmax(both$MCSE_bias.x, both$MCSE_bias.y)))
##
## ---- The stratified-bivariate model against the best model, for each proportion missing (ratio of the mean RMSEs):
##
RMSE_pct <-  means[["RMSE pct"]]
for (p in colnames(RMSE_pct)) {
      complete <-  RMSE_pct[rownames(RMSE_pct) != "Stein", p]
      numbers[[paste0("Strat_over_best.pct.", p)]] <-  fn_1(RMSE_pct["Strat", p] / min(complete))
      numbers[[paste0("Strat_pct_worse.pct.", p)]] <-  fn_0(100 * (RMSE_pct["Strat", p] / min(complete) - 1))
      numbers[[paste0("best_id.pct.", p)]] <-  names(complete)[which.min(complete)]
}
RMSE_npct <-  means[["RMSE npct"]]
for (k in colnames(RMSE_npct)) numbers[[paste0("Strat_pct_worse.npct.", k)]] <-  fn_0(100 * (RMSE_npct["Strat", k] / min(RMSE_npct[rownames(RMSE_npct) != "Stein", k]) - 1))
##
## ---- In how many of the scenarios of a family of DGMs each model is in "Group 1 - Best" by RMSE:
##
fn_best_group <- function(cells) {
      in_the_best_group <-  rep(NA, nrow(cells))
      for (s in unique(cells$scenario)) {
            i <-  which(cells$scenario == s)
            best <-  i[which.min(cells$avg_RMSE[i])]
            in_the_best_group[i] <-  !((cells$avg_RMSE[i] - cells$avg_RMSE[best]) >= 4 * pmax(cells$MCSE_RMSE[i], cells$MCSE_RMSE[best])) &
                                     !(100 * (cells$avg_RMSE[i] - cells$avg_RMSE[best]) / cells$avg_RMSE[best] >= 10)
      }
      in_the_best_group
}
s0$best <-  fn_best_group(s0)
t <-  tapply(s0$best, list(s0$id, s0$family), sum)
for (m in rownames(t)) for (f in colnames(t)) numbers[[paste("best.family", f, m, sep = ".")]] <-  as.character(t[m, f])
for (f in unique(s0$family)) numbers[[paste0("n_scenarios.family.", f)]] <-  as.character(length(unique(s0$scenario[s0$family == f])))
t <-  tapply(s0$best, list(s0$id, s0$test), sum)
for (m in rownames(t)) for (f in colnames(t)) numbers[[paste("best.test", f, m, sep = ".")]] <-  as.character(t[m, f])
##
## ---- Pairwise: is model b significantly worse than model a (a difference in means of at least 4 x the larger of their mean MCSEs),
## ---- over the scenarios of a group in which both have a result:
##
fn_pair <- function(rows, a, b, value, MCSE) {
      common <-  intersect(s0$scenario[rows & s0$id == a], s0$scenario[rows & s0$id == b])
      x <-  rows & s0$id == a & s0$scenario %in% common
      y <-  rows & s0$id == b & s0$scenario %in% common
      as.character(mean(value[y]) - mean(value[x]) >= 4 * max(mean(MCSE[x]), mean(MCSE[y])))
}
for (m in c("RC", "FC")) {
      numbers[[paste0("sig.RMSE.test.BAI.", m, "_vs_Jones")]] <-  fn_pair(s0$test == "BAI", m, "Jones", s0$avg_RMSE, s0$MCSE_RMSE)
      for (o in c("Jones", "Stein", "GLMM", "Hoyer")) numbers[[paste0("sig.bias.all.", m, "_vs_", o)]] <-  fn_pair(rep(TRUE, nrow(s0)), m, o, s0$avg_abs_bias, s0$MCSE_bias)
      for (p in c(15, 55, 85)) numbers[[paste0("sig.RMSE.set.BAI_", p, ".", m, "_vs_Jones")]] <-  fn_pair(s0$test == "BAI" & s0$pct == p, m, "Jones", s0$avg_RMSE, s0$MCSE_RMSE)
}
##
## ---- The lower bias of the O-bivariate models against the Jones model for the BAI (number of scenarios):
##
bai <-  s0[s0$test == "BAI", ]
jones <-  bai[bai$id == "Jones", c("scenario", "avg_abs_bias")]
for (m in c("FC", "RC")) {
      both <-  merge(bai[bai$id == m, c("scenario", "avg_abs_bias")], jones, by = "scenario")
      numbers[[paste0("BAI_lower_bias_than_Jones.", m)]] <-  as.character(sum(both$avg_abs_bias.x < both$avg_abs_bias.y))
      numbers[[paste0("BAI_lower_bias_than_Jones.n.", m)]] <-  as.character(nrow(both))
}
##
## ---- Ranges of the mean coverage over the three proportions missing, by model:
##
for (measure in c("cvgSe", "cvgSp")) {
      t <-  means[[paste(measure, "pct")]]
      for (m in rownames(t)) {
            numbers[[paste0(measure, ".min.", m)]] <-  fn_1(min(t[m, ], na.rm = TRUE))
            numbers[[paste0(measure, ".max.", m)]] <-  fn_1(max(t[m, ], na.rm = TRUE))
      }
      frequentist <-  t[c("Stein", "GLMM", "Hoyer"), ]
      numbers[[paste0(measure, ".min.frequentist")]] <-  fn_1(min(frequentist, na.rm = TRUE))
      numbers[[paste0(measure, ".max.frequentist")]] <-  fn_1(max(frequentist, na.rm = TRUE))
}
##
## ---- The five Bayesian models under prior sets 1 and 2 (means over all scenarios):
##
for (k in c("set_1", "set_2")) {
      cells <-  d[d$prior_set == k, ]
      for (m in unique(cells$id)) {
            numbers[[paste("RMSE", k, m, sep = ".")]] <-  fn_2(mean(cells$avg_RMSE[cells$id == m]))
            numbers[[paste("bias", k, m, sep = ".")]] <-  fn_2(mean(cells$avg_abs_bias[cells$id == m]))
      }
      r <-  rankings[rankings$measure == "RMSE" & rankings$prior_set == k, ]
      r <-  r[order(r$rank), ]
      numbers[[paste0("order.RMSE.", k)]] <-  fn_list(name_in_the_text[r$id])
      numbers[[paste0("ids.RMSE.", k)]] <-  paste(r$id, collapse = " ")
      numbers[[paste0("n_scenarios.", k)]] <-  as.character(length(unique(cells$scenario)))
}
##
## ---- The cells: how many, how many over the MCSE target, and of which models:
##
numbers[["cells"]] <-  format(nrow(d), big.mark = "{,}")
## (over the target = MCSE of the RMSE above 0.125 by more than 0.0005, the rule of the daggers of the detailed tables)
over <-  d[d$MCSE_val > 0.125 + 0.0005, ]
numbers[["cells_over_the_target"]] <-  as.character(nrow(over))
numbers[["cells_over_the_target.Strat"]] <-  as.character(sum(over$id == "Strat"))
numbers[["cells_over_the_target.largest_MCSE"]] <-  fn_2(max(c(over$MCSE_val, 0)))
numbers[["cells_over_the_target.by_model"]] <-  paste(sprintf("%s %d", names(table(over$id)), as.integer(table(over$id))), collapse = "; ")
numbers[["cells_over_the_target.by_set"]] <-  paste(sprintf("%s %d", names(table(paste(over$test, over$pct, over$n_studies))), as.integer(table(paste(over$test, over$pct, over$n_studies)))), collapse = "; ")
##
## ---- The fits kept and lost:
##
numbers[["fits"]] <-  format(fits$fits, big.mark = "{,}")
for (k in names(fits)) if (grepl("^pct_|^largest_pct|^median_seconds", k)) numbers[[paste0("fits.", k)]] <-  fn_1(fits[[k]])
numbers[["fits.cell_with_the_largest_pct_lost"]] <-  fits$cell_with_the_largest_pct_lost
numbers[["fits.cells_over_10_pct_lost"]] <-  as.character(fits$cells_over_10_pct_lost)
writeLines(jsonlite::toJSON(numbers, auto_unbox = TRUE, pretty = TRUE), out)
cat("numbers written:", length(numbers), "\n")
for (k in c("order.RMSE.all", "starred.RMSE.all", "order.bias.all", "starred.bias.all", "ids.RMSE.pct15", "ids.RMSE.pct55", "ids.RMSE.pct85", "starred.RMSE.pct85", "ids.RMSE.GAD_2", "ids.RMSE.HADS", "ids.RMSE.BAI",
            "FC_vs_RC.n", "FC_vs_RC.not_different.RMSE", "FC_vs_RC.not_different.bias", "Strat_over_best.pct.15", "Strat_over_best.pct.55", "Strat_over_best.pct.85",
            "cells", "cells_over_the_target", "cells_over_the_target.by_model", "cells_over_the_target.by_set", "order.RMSE.set_1", "order.RMSE.set_2"))
      cat(sprintf("%-34s %s\n", k, numbers[[k]]))
