##
## ---- One table of every cell of the nine sets of the factorial design (3 tests x 3 shares of thresholds missing) -----------------
##
## GAD-2 at 15% and BAI at 55% are the cells of the earlier sim (each test at the share missing of its real data: 16.1% and 54.4%),
## scored with the same rule (a threshold counts when at least 1 study reports it); the other seven sets are the factorial sim.
## The earlier sim's HADS set (41% missing) is outside the design and is left out.
## The ordinal-HSROC model and its DGM are named as the paper's scripts name them ("R&G", "R&G_random_C").
## Written: the joined table, and one table for each share missing (with the three tests), for the paper's table and plot scripts.
##
## Usage: Rscript <this file> <table of the earlier sim, .rds> <table of the factorial sim, .rds> <folder to write in>
##
arguments <-  commandArgs(trailingOnly = TRUE)
earlier <-  readRDS(arguments[1])
factorial <-  readRDS(arguments[2])
folder <-  arguments[3]
earlier <-  earlier[earlier$test %in% c("GAD_2", "BAI"), ]
earlier$pct <-  ifelse(earlier$test == "GAD_2", 15, 55)
earlier$sim <-  "earlier"
factorial$sim <-  "factorial"
factorial$DGM[factorial$DGM == "ord_HSROC_random_C"] <-  "R&G_random_C"
factorial$model_parameterisation[factorial$model_parameterisation == "ord_HSROC"] <-  "R&G"
stopifnot(setequal(names(earlier), names(factorial)))
joined <-  rbind(earlier[, names(factorial)], factorial)
stopifnot(!anyDuplicated(joined[, c("model", "prior_set", "test", "pct", "n_studies", "DGM")]))
saveRDS(joined, file.path(folder, "table_of_every_cell_of_the_nine_sets.rds"))
for (p in c(15, 55, 85)) saveRDS(joined[joined$pct == p, ], file.path(folder, sprintf("table_of_every_cell_%d_pct_missing.rds", p)))
cat("cells:", nrow(joined), "\n")
print(with(joined, table(paste(test, pct), prior_set)))
print(with(joined, table(model, prior_set)))
cat("cells over the MCSE target:", sum(!joined$mcse_met), "\n")
