##
## ---- The figures of the paper's appendices B, C, D (RMSE, bias, coverage; GAD-2, HADS, BAI) from the table of every cell of the
##      whole sim, with EVERY model of the sim on them (the five MetaOrdDTA models and the three competing models) -----------------
##
## The plots are made by the paper's own functions (create_RMSE_plot, create_bias_plot, create_coverage_plot of
## thr_sim_study_helper_fns.R), as make_the_plots_of_the_paper_from_the_table_of_the_cells_of_the_whole_sim.R does. Those functions
## name five models only, so three things are set here on top of them, and nothing in the package is changed:
##   - the models on the x axis and their colours (the paper's five colours, and three more of the same palette): the models that have rows at the prior set given, in the order below (the three
##     competing models are fitted as their authors give them, so they have rows at prior set 0 only);
##   - the names of the DGMs in the strips, on two lines, in a size that fits six rows of panels;
##   - the name of each model (the paper's function gives no name to the three competing models).
## The dashed line of each panel (the upper limit of the best model) is drawn by the paper's function from the rows it is given,
## which are exactly the models shown.
## Usage: Rscript <this file> <folder of the MetaOrdDTA package> <table of every cell, .rds> <prior set, e.g. set_0>
##        <TRUE: only the cells with the MCSE target met; FALSE: every cell of the table> <folder to save the figures in>
##        <text added to the end of each file name before .png, or the word none>
##
arguments <-  commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 6)
local_pkg_dir <-  arguments[1]
table_of_the_cells <-  readRDS(arguments[2])
prior_set <-  arguments[3]
only_the_cells_with_the_MCSE_target_met <-  as.logical(arguments[4])
folder_of_the_figures <-  arguments[5]
end_of_the_file_names <-  if (arguments[6] == "none") "" else arguments[6]
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(tidyr); library(purrr); library(patchwork) })
source(file.path(local_pkg_dir, "sim_study", "thr_sim_study_helper_fns.R"))
dir.create(folder_of_the_figures, recursive = TRUE, showWarnings = FALSE)
rows <-  table_of_the_cells[table_of_the_cells$prior_set == prior_set, ]
if (only_the_cells_with_the_MCSE_target_met) rows <-  rows[rows$mcse_met, ]
stopifnot(nrow(rows) > 0)
##
## ---- The names of the models (the first five are the paper's own) and of the DGMs:
##
names_of_the_models <-  c("Bivariate" = "Bivariate", "Jones" = "Jones", "ord_bivariate FALSE" = "Ord-Biv-FC", "ord_bivariate TRUE" = "Ord-Biv-RC",
                          "R&G" = "Ord-HSROC-RC", "Hoyer_2018" = "Hoyer-2018", "logit_GLMM" = "Logit-GLMM", "diagmeta" = "Steinhauser-2016")
names_of_the_DGMs <-  c("Jones" = "DGM:\nJones", "ord_biv_fixed_C" = "DGM:\nOrd-Biv-FC", "ord_biv_random_C" = "DGM:\nOrd-Biv-RC", "R&G_random_C" = "DGM:\nOrd-HSROC-RC",
                        "Hoyer_2018" = "DGM:\nHoyer 2018", "logit_GLMM" = "DGM:\nLogit GLMM", "diagmeta" = "DGM:\nSteinhauser\n2016")
fn_rows_for_the_plots <- function(rows_of_the_test) {

      key <-  ifelse(rows_of_the_test$model_parameterisation == "ord_bivariate",
                     paste(rows_of_the_test$model_parameterisation, rows_of_the_test$random_thresholds), rows_of_the_test$model_parameterisation)
      stopifnot(all(key %in% names(names_of_the_models)), all(rows_of_the_test$DGM %in% names(names_of_the_DGMs)))
      rows_of_the_test$Model <-  factor(unname(names_of_the_models[key]), levels = unname(names_of_the_models))
      rows_of_the_test$DGM <-  factor(rows_of_the_test$DGM, levels = names(names_of_the_DGMs))
      rows_of_the_test

}
## (the five colours of the paper's figures for its five models; brown, pink and grey of the same palette for the three competing
##  models: the palette's sixth colour, yellow, does not show on white)
colours_of_the_models <-  stats::setNames(RColorBrewer::brewer.pal(9, "Set1")[c(1:5, 7, 8, 9)], unname(names_of_the_models))
base_size <-  14
dpi <-  200
mult <-  1.1
for (test in c("GAD_2", "HADS", "BAI")) {

      rows_of_the_test <-  dplyr::as_tibble(rows[rows$test == test, ])
      if (nrow(rows_of_the_test) == 0) next
      overall_sim_data <-  fn_rows_for_the_plots(rows_of_the_test)
      models_shown <-  levels(droplevels(overall_sim_data$Model))
      ## (the height of a figure follows its rows of panels: 16 for the six DGMs the figures were first made with)
      height <-  16 * length(unique(overall_sim_data$DGM)) / 6
      fn_with_every_model <- function(plot) {

            suppressMessages(plot +
                             scale_x_discrete(limits = models_shown) +
                             scale_color_manual(values = colours_of_the_models[models_shown], limits = models_shown) +
                             facet_grid(DGM ~ n_studies,
                                        scales = "free_y",
                                        labeller = labeller(n_studies = function(x) paste(x, "studies"),
                                                            DGM = function(x) unname(names_of_the_DGMs[as.character(x)]))) +
                             theme(strip.text.y = element_text(size = base_size * 0.55, face = "bold", lineheight = 0.9)))

      }
      plot_RMSE <-  fn_with_every_model(create_RMSE_plot(data = overall_sim_data, metric_type = "Se", base_size = base_size)) +
                    fn_with_every_model(create_RMSE_plot(data = overall_sim_data, metric_type = "Sp", base_size = base_size)) +
                    plot_layout(ncol = 1)
      plot_Bias <-  fn_with_every_model(create_bias_plot(data = overall_sim_data, metric_type = "Se", base_size = base_size)) +
                    fn_with_every_model(create_bias_plot(data = overall_sim_data, metric_type = "Sp", base_size = base_size)) +
                    plot_layout(ncol = 1)
      plot_Cvg <-   fn_with_every_model(create_coverage_plot(data = overall_sim_data, metric_type = "Se", base_size = base_size)) +
                    fn_with_every_model(create_coverage_plot(data = overall_sim_data, metric_type = "Sp", base_size = base_size)) +
                    plot_layout(ncol = 1)
      ggsave(file.path(folder_of_the_figures, paste0("Sim_study_RMSE_", test, end_of_the_file_names, ".png")), plot_RMSE,
             width = 11 * mult, height = height * mult, dpi = dpi)
      ggsave(file.path(folder_of_the_figures, paste0("Sim_study_Bias_", test, end_of_the_file_names, ".png")), plot_Bias,
             width = 11 * mult, height = height * mult, dpi = dpi)
      ggsave(file.path(folder_of_the_figures, paste0("Sim_study_Coverage_", test, end_of_the_file_names, ".png")), plot_Cvg,
             width = 11 * mult, height = height * mult, dpi = dpi)
      cat(test, prior_set, ": 3 figures from", nrow(rows_of_the_test), "cells; models:", paste(models_shown, collapse = ", "), "\n")

}
