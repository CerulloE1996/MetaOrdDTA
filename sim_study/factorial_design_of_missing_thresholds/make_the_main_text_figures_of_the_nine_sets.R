##
## ---- The three results figures of the main manuscript, from the table of every cell of the nine sets of the factorial design -------
##
## RMSE, bias and coverage against the proportion of thresholds missing: one line a model, one panel a test (columns) and a number
## of studies (rows); each point = mean over the 7 DGMs of the scenario's value at prior set 0 (RMSE = RMSE(Se) + RMSE(Sp), bias =
## |bias(Se)| + |bias(Sp)|, as in the tables; coverage of Se and of Sp apart). The bars of the RMSE and bias figures are
## +/- 1.96 x (mean MCSE of the 7 cells / sqrt(7)). Colours: those of the paper's other simulation figures.
##
## Usage: Rscript <this file> <table of every cell of the nine sets, .rds> <folder to save the figures in>
##
arguments <-  commandArgs(trailingOnly = TRUE)
d <-  readRDS(arguments[1])
folder <-  arguments[2]
suppressPackageStartupMessages({ library(ggplot2) })
dir.create(folder, recursive = TRUE, showWarnings = FALSE)
d <-  d[d$prior_set == "set_0", ]
models <-  c("strat. bivariate" = "Strat-biv", "Jones" = "Jones", "O-biv-FC" = "O-biv-FC", "O-biv-RC" = "O-biv-RC", "ord-HSROC" = "O-HSROC-RC",
             "Hoyer 2018" = "Hoyer", "logit GLMM" = "Logit-GLMM", "diagmeta" = "Steinhauser")
colours <-  stats::setNames(RColorBrewer::brewer.pal(9, "Set1")[c(1:5, 7, 8, 9)], unname(models))
d$Model <-  factor(models[d$model], levels = unname(models))
d$Test <-  factor(c(GAD_2 = "GAD-2 (7 categories)", HADS = "HADS (22 categories)", BAI = "BAI (64 categories)")[d$test],
                  levels = c("GAD-2 (7 categories)", "HADS (22 categories)", "BAI (64 categories)"))
d$Studies <-  factor(paste(d$n_studies, "studies"), levels = c("10 studies", "50 studies"))
fn_means <- function(value,
                     MCSE = NULL) {

      groups <-  split(seq_len(nrow(d)), list(d$Model, d$Test, d$Studies, d$pct), drop = TRUE)
      do.call(rbind, lapply(groups, function(i) data.frame(Model = d$Model[i[1]], Test = d$Test[i[1]], Studies = d$Studies[i[1]], pct = d$pct[i[1]], n = length(i),
                                                           mean = mean(value[i]), half_width = if (is.null(MCSE)) NA else 1.96 * mean(MCSE[i]) / sqrt(length(i)))))

}
fn_plot <- function(rows,
                    y_name,
                    bars = TRUE) {

      plot <-  ggplot(rows, aes(x = pct, y = mean, colour = Model, group = Model)) +
               geom_line(linewidth = 0.7) +
               geom_point(size = 2.2) +
               scale_x_continuous(breaks = c(15, 55, 85), labels = function(x) paste0("~", x, "%"), limits = c(8, 92)) +
               scale_colour_manual(values = colours) +
               labs(x = "Proportion of thresholds missing", y = y_name, colour = NULL) +
               theme_bw(base_size = 13) +
               theme(legend.position = "bottom", panel.grid.minor = element_blank(), strip.text = element_text(face = "bold")) +
               guides(colour = guide_legend(nrow = 1))
      if (bars) plot <-  plot + geom_errorbar(aes(ymin = mean - half_width, ymax = mean + half_width), width = 3, linewidth = 0.4)
      plot

}
MCSE_RMSE <-  pmax(d$avg_MCSE_RMSE_Se, d$avg_MCSE_RMSE_Sp, na.rm = TRUE)
MCSE_bias <-  pmax(d$avg_MCSE_bias_Se, d$avg_MCSE_bias_Sp, na.rm = TRUE)
rows_RMSE <-  fn_means(d$avg_RMSE, MCSE_RMSE)
rows_bias <-  fn_means(d$avg_abs_bias, MCSE_bias)
plot_RMSE <-  fn_plot(rows_RMSE, "Mean RMSE(Se) + RMSE(Sp), percentage points") + facet_grid(Studies ~ Test, scales = "free_y")
plot_bias <-  fn_plot(rows_bias, "Mean |bias(Se)| + |bias(Sp)|, percentage points") + facet_grid(Studies ~ Test, scales = "free_y")
rows_coverage <-  rbind(cbind(fn_means(100 * d$avg_coverage_Se), Measure = "Sensitivity"), cbind(fn_means(100 * d$avg_coverage_Sp), Measure = "Specificity"))
rows_coverage$Row <-  factor(paste0(rows_coverage$Measure, ", ", rows_coverage$Studies),
                             levels = c("Sensitivity, 10 studies", "Sensitivity, 50 studies", "Specificity, 10 studies", "Specificity, 50 studies"))
plot_coverage <-  fn_plot(rows_coverage, "Mean coverage of the 95% intervals (%)", bars = FALSE) + geom_hline(yintercept = 95, linetype = "dashed", colour = "grey35") +
                  facet_grid(Row ~ Test)
ggsave(file.path(folder, "Sim_study_RMSE_by_proportion_of_thresholds_missing.png"), plot_RMSE, width = 11.5, height = 7, dpi = 200)
ggsave(file.path(folder, "Sim_study_Bias_by_proportion_of_thresholds_missing.png"), plot_bias, width = 11.5, height = 7, dpi = 200)
ggsave(file.path(folder, "Sim_study_Coverage_by_proportion_of_thresholds_missing.png"), plot_coverage, width = 11.5, height = 11, dpi = 200)
utils::write.csv(rbind(cbind(what = "RMSE", rows_RMSE, Measure = NA, Row = NA), cbind(what = "bias", rows_bias, Measure = NA, Row = NA), cbind(what = "coverage", rows_coverage)),
                 file.path(folder, "every_number_of_the_main_text_figures.csv"), row.names = FALSE)
cat("figures written: 3; points with fewer than 7 DGMs:", sum(rows_RMSE$n < 7), "of", nrow(rows_RMSE), "\n")
