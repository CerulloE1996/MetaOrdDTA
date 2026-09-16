



rm(list = ls())
.rs.restartR()  # In RStudio


require(ggplot2)


##
## ---- For REAL data (NOTE: THIS MAY CHANGE WHEN INCLUDE ALL 8 TESTS!!! AS MORE STUDIES TOO!!!)  -------------------------------------------------
##
# "While study setting is often considered in diagnostic meta-analyses, our k-fold cross-validation 
# revealed it actually degraded predictive performance (ΔELPD = -973), suggesting these broad categorizations
# may not capture meaningful heterogeneity."
##
# Of four pre-specified covariates, only reference test type improved predictive performance (ΔELPD = +416). 
# Disease prevalence had negligible effect (+23), while study setting (-973) and risk of bias (-288) substantially
# degraded predictions, suggesting these categorical variables introduce noise rather than explain
# meaningful heterogeneity
##
## ---- How to explain sim study reasoning in paper:
##
# "In the simulated example (generated from a model with all covariates having non-zero effects), k-fold CV correctly
# identified that disease prevalence improved predictive performance (ΔELPD = 3557), unlike in real data where it provided 
# negligible benefit. This demonstrates our methodology's ability to detect dataset-specific covariate importance."
##
##
# You're exactly right about the real data - stop there and use intercept + reference test as your main model!
# For the dummy data, I'd recommend this strategy:
#   Skip RoB (only 128 ELPD) and work with the 3 useful covariates:
#   Option 1: Incremental approach (I recommend this)
# 
# You already have the univariate results
# Fit one model: intercept + prev_GAD + reference_test + study_setting
# See if it beats prev_GAD alone (which was your best univariate)
# 
# My recommendation: Go with Option 1 because:
#   
#   K-fold already told you which covariates help individually
# You're demonstrating the methodology, not doing exhaustive model selection
# It mirrors what you'd do in practice: "these 3 helped, let's try them together"
# 
# Then for your final models:
#   Real data: intercept + reference_test
# 
# SROC curves by reference test type
# Se/Sp tables/plots
# AUC differences between tests
# 
# Dummy data: intercept + prev_GAD + reference_test + study_setting (if it wins)
# 
# Or just prev_GAD if the combined model doesn't improve
# Same outputs as above
# 
# This gives you a clean story: "We used k-fold to identify useful covariates, combined those that helped,
## and present results for the best model." Perfect for a thesis demonstration! 🎯
# 

 


# # In your paper:
# "All plots were created using the MetaOrdDTA R package. Base plots were 
# generated using the following R6 class methods:
# 
# - sROC curves: `model$plot$sroc()`
# - Se/Sp vs threshold: `model$plot$threshold()`  
# - Forest plots: `model$plot$forest()`
# - Pairwise comparisons: `model$plot$comparisons()`
# 
# Custom layouts were created using patchwork. Full code for all figures 
# is available at: https://github.com/yourname/MetaOrdDTA-paper"


# 
# ## If WTE reduction kicks in 1st September 2025:
# 4*(0.3/0.2) ## = 6, so get extra 2 months (so ends 28th feb rather than 31st december)
# ##
# ## If WTE reduction kicks in 1st October 2025:
# 3*(0.3/0.2) ## = 4.5 ## so get extra 14th feb
# 
# 

# 
 


 

##
## ----------- Plot sROC curve(s) (NMA - for best model WITHOUT covariates): ------------------------------------------
##
relevant_thresholds = list( "GAD-2" = c(1:6), 
                            "GAD-7" =  c(3:18), 
                            "HADS" =  c(3:18),
                            "BAI" = c(3:38))
# stan_model_file_name = model_summary_and_trace_obj_DUMMY$internal_obj$outs_stan_model_name$stan_model_file_name
# stan_mod_samples     = model_summary_and_trace_obj_DUMMY$internal_obj$outs_stan_sampling$stan_mod_samples
# ##
# test_names = test_names
# # covariates = covariates,
# # new_cov_data = new_cov_data,
# ##
# # df_true = df_true
# ##
# conf_region_colour = "blue"
# pred_region_colour = "green"
# ##
# base_size = 24
# point_size = 4.0
# linewidth = 1.0
# ##
# n_index_tests = n_index_tests
# n_thr = n_thr
#
#
##
# plots_REAL <- model_summary_and_trace_obj_REAL$plot_sROC(test_names = test_names)
##
# plots_REAL <- R_fn_sROC_plot_NMA(     stan_model_file_name = model_summary_and_trace_obj_REAL$internal_obj$outs_stan_model_name$stan_model_file_name,
#                                       stan_mod_samples     =  model_summary_and_trace_obj_REAL$internal_obj$outs_stan_sampling$stan_mod_samples,
#                                       ##
#                                       test_names = test_names,
#                                       ##
#                                       relevant_thresholds = relevant_thresholds,
#                                       ##
#                                       # covariates = covariates,
#                                       # new_cov_data = new_cov_data,
#                                       ##
#                                       df_true = NULL,
#                                       ##
#                                       conf_region_colour = "blue", 
#                                       pred_region_colour = "green",
#                                       ##
#                                       conf_region_alpha = 0.30,
#                                       pred_region_alpha = 0.15,
#                                       ##
#                                       base_size = 24,
#                                       point_size = 4.0,
#                                       linewidth = 1.0,
#                                       ##
#                                       n_index_tests = n_index_tests,
#                                       n_thr = n_thr)
# ##
# plots_REAL$plot_list[[1]]
# #### plots_REAL$plot_list[[2]]
# #### plots_REAL$plot_list[[3]]
# #### plots_REAL$plot_list[[4]]
# plots_REAL$plot_list[[5]]
# #### plots_REAL$plot_list[[6]]
# plots_REAL$plot_list[[7]]
# plots_REAL$plot_list[[8]]
# plots_REAL$plot_list[[9]]
# ##
# ggsave("REAL_Application_baseline_sROC.png", plots_REAL$plot_list[[1]], width = 16, height = 9, dpi = 400)
# ggsave("REAL_Application_baseline_sROC_panel_w_CrI_PrI.png", plots_REAL$plot_list[[6]], width = 16, height = 9, dpi = 400)
# ##
# ggsave("REAL_Application_baseline_Se_Sp_vs_thr.png", plots_REAL$plot_list[[7]], width = 16, height = 16, dpi = 400)
# ggsave("REAL_Application_baseline_Se_Sp_vs_thr_w_CrI.png", plots_REAL$plot_list[[8]], width = 16, height = 16, dpi = 400)
# ggsave("REAL_Application_baseline_Se_Sp_vs_thr_panel_w_CrI_PrI.png", plots_REAL$plot_list[[9]], width = 16, height = 16, dpi = 400)
# ##
##
##
# plots_DUMMY <- model_summary_and_trace_obj_DUMMY$plot_sROC(test_names = test_names)
##
n_thr <- c(6, 21, 21, 63)
##
plots_DUMMY <- R_fn_sROC_plot_NMA(     stan_model_file_name = model_summary_and_trace_obj_DUMMY$internal_obj$outs_stan_model_name$stan_model_file_name,
                                       stan_mod_samples     = model_summary_and_trace_obj_DUMMY$internal_obj$outs_stan_sampling$stan_mod_samples,
                                       ##
                                       test_names = test_names,
                                       ##
                                       relevant_thresholds = relevant_thresholds,
                                       ##
                                       # covariates = covariates,
                                       # new_cov_data = new_cov_data,
                                       ##
                                       df_true = NULL,
                                       ##
                                       conf_region_colour = "blue", 
                                       pred_region_colour = "green",
                                       ##
                                       conf_region_alpha = 0.20,
                                       pred_region_alpha = 0.10,
                                       ##
                                       base_size = 24,
                                       point_size = 4.0,
                                       linewidth = 1.0,
                                       ##
                                       n_index_tests = n_index_tests,
                                       n_thr = n_thr)
##
plots_DUMMY$plot_list[[1]]
#### plots_DUMMY$plot_list[[2]]
#### plots_DUMMY$plot_list[[3]]
#### plots_DUMMY$plot_list[[4]]
plots_DUMMY$plot_list[[5]]
#### plots_DUMMY$plot_list[[6]]
plots_DUMMY$plot_list[[7]]
plots_DUMMY$plot_list[[8]]
plots_DUMMY$plot_list[[9]]
##
ggsave("DUMMY_Application_baseline_sROC.png", plots_DUMMY$plot_list[[1]], width = 16, height = 9, dpi = 400)
ggsave("DUMMY_Application_baseline_sROC_panel_w_CrI_PrI.png", plots_DUMMY$plot_list[[5]], width = 16, height = 16, dpi = 400)
##
ggsave("DUMMY_Application_baseline_Se_Sp_vs_thr.png", plots_DUMMY$plot_list[[7]], width = 16, height = 16, dpi = 400)
ggsave("DUMMY_Application_baseline_Se_Sp_vs_thr_w_CrI.png", plots_DUMMY$plot_list[[8]], width = 16, height = 16, dpi = 400)
ggsave("DUMMY_Application_baseline_Se_Sp_vs_thr_panel_w_CrI_PrI.png", plots_DUMMY$plot_list[[9]], width = 16, height = 16, dpi = 400)
# ##
# ## ---- Compare real to simulated/dummy data (for personal interest NOT FOR PAPER/THESIS):
# ##
# require(patchwork)
# plots_DUMMY$plot_list[[1]] + plots_REAL$plot_list[[1]]




##
## ---- Extract NMA Se/Sp comparison metrics for tables + plots  ------------------------------------------------------------------------
##




draws_array <- model_summary_and_trace_obj_DUMMY$get_all_traces_list()$traces_as_arrays$draws_array
# trace_gq <-  model_summary_and_trace_obj_DUMMY$get_all_traces_list()$traces_as_arrays$trace_gq
# str(trace_gq)
##
outs_NMA_comp_DUMMY <-  R_fn_using_Rcpp_compute_NMA_comparisons_posthoc(  trace_gq = draws_array,
                                                                          n_thr = n_thr,
                                                                          test_names = test_names)
##






# Create separate plots for each clinical context

# Define clinical threshold groups
get_threshold_groups <- function() {
  
  list(
  
    screening = list(
      "GAD-2" = c(2, 3),         # 3 is SCREENING standard for GAD
      "GAD-7" = c(9, 10),    # 10 is SCREENING standard  
      # "HADS" = c(7, 8, 9),       # 8 for SCREENING ("possible") (11 is "probable")
      "HADS" = c(9, 10, 11),       # 8 for SCREENING ("possible") (11 is "probable")
      # "BAI" = c(7, 8, 9)         # 8 for mild (screening??)
        "BAI" = c(17, 18, 19)      # 16 is standard
    ),
    
    ## More thresholds here, use this for more detailed (e.g. appendix for most/all
    ## of them) pairwise plots (choose(4, 2) = 6 plots):
    screening_2 = list(
      "GAD-2" = c(2, 3, 4),      # 3 is SCREENING standard for GAD
      "GAD-7" = c(8, 9, 10, 11),    # 10 is SCREENING standard for GAD
      "HADS" = c(8, 9, 10, 11, 12),    # 11 is standard
      "BAI" = c(16, 17, 18, 19, 20)      # 16 is standard
    )
 
  )
  
}

##
threshold_groups <- get_threshold_groups()
threshold_groups
##
threshold_groups$screening
##
screening_comps <- filter_clinically_relevant_comparisons( tibble_comp = outs_NMA_comp_DUMMY,
                                                           relevant_thresholds = threshold_groups$screening) ## %>%
                   # identify_important_comparisons(n_per_type = 100)
screening_comps %>% print(n = 1000)



##
## ---- Create forest plots for Se + Sp - screening thresholds for GAD:
##
plot_dpi <- 100
##
plot_screening_Se <- plot_se_differences( filtered_tibble = screening_comps,
                                          show_top_n = 1000,
                                          order_by = "default",
                                          plot_title = "Sensitivity (Se) Differences: Screening Thresholds")
##
plot_screening_Sp <- plot_sp_differences( filtered_tibble = screening_comps,
                                          show_top_n = 1000,
                                          order_by = "default",
                                          plot_title = "Specificity (Sp) Differences: Screening Thresholds")
## 
plot_screening_Se_and_Sp <- plot_screening_Se + plot_screening_Sp
plot_screening_Se_and_Sp
##
ggsave("DUMMY_Application_Se_and_Sp_pairwise_diffs_screening_GAD.png", 
       plot_screening_Se_and_Sp, 
       width = 16, height = 16, dpi = plot_dpi)


##
## ---- Se/Sp pairwise diff's subsetted by test-pair combos's (choose(4, 2) = 6 plots in total):
##
test_names
##
screening_comps_2 <- filter_clinically_relevant_comparisons( tibble_comp = outs_NMA_comp_DUMMY,
                                                           relevant_thresholds = threshold_groups$screening_2) ## %>%




##
## ---- GAD-2 vs. GAD-7:
##
test_pair <- c("GAD-2", "GAD-7")
screening_comps_GAD_2_vs_GAD_7 <- dplyr::filter(screening_comps_2, (test1_name %in% test_pair) & (test2_name %in% test_pair))
##
plot_screening_Se <- plot_se_differences( filtered_tibble = screening_comps_GAD_2_vs_GAD_7,
                                          show_top_n = 1000,
                                          order_by = "default",
                                          plot_title = "Sensitivity (Se) Differences: Screening Thresholds")
##
plot_screening_Sp <- plot_sp_differences( filtered_tibble = screening_comps_GAD_2_vs_GAD_7,
                                          show_top_n = 1000,
                                          order_by = "default",
                                          plot_title = "Specificity (Sp) Differences: Screening Thresholds")
## 
plot_screening_Se_and_Sp_GAD_2_vs_GAD_7 <- plot_screening_Se + plot_screening_Sp
plot_screening_Se_and_Sp_GAD_2_vs_GAD_7
##
ggsave("DUMMY_Application_Se_and_Sp_pairwise_diffs_screening_GAD_GAD_2_vs_GAD_7.png", 
       plot_screening_Se_and_Sp_GAD_2_vs_GAD_7, 
       width = 16, height = 16, dpi = plot_dpi)



##
## ---- GAD-2 vs. HADS
##
test_pair <- c("GAD-2", "HADS")
screening_comps_GAD_2_vs_HADS <- dplyr::filter(screening_comps_2, (test1_name %in% test_pair) & (test2_name %in% test_pair))
##
plot_screening_Se <- plot_se_differences( filtered_tibble = screening_comps_GAD_2_vs_HADS,
                                          show_top_n = 1000,
                                          order_by = "default",
                                          plot_title = "Sensitivity (Se) Differences: Screening Thresholds")
##
plot_screening_Sp <- plot_sp_differences( filtered_tibble = screening_comps_GAD_2_vs_HADS,
                                          show_top_n = 1000,
                                          order_by = "default",
                                          plot_title = "Specificity (Sp) Differences: Screening Thresholds")
## 
plot_screening_Se_and_Sp_GAD_2_vs_HADS <- plot_screening_Se + plot_screening_Sp
plot_screening_Se_and_Sp_GAD_2_vs_HADS
##
ggsave("DUMMY_Application_Se_and_Sp_pairwise_diffs_screening_GAD_GAD_2_vs_HADS.png", 
       plot_screening_Se_and_Sp_GAD_2_vs_HADS, 
       width = 16, height = 16, dpi = plot_dpi)



##
## ---- GAD-7 vs. HADS:
##
test_pair <- c("GAD-7", "HADS")
screening_comps_GAD_7_vs_HADS <- dplyr::filter(screening_comps_2, (test1_name %in% test_pair) & (test2_name %in% test_pair))
##
plot_screening_Se <- plot_se_differences( filtered_tibble = screening_comps_GAD_7_vs_HADS,
                                          show_top_n = 1000,
                                          order_by = "default",
                                          plot_title = "Sensitivity (Se) Differences: Screening Thresholds")
##
plot_screening_Sp <- plot_sp_differences( filtered_tibble = screening_comps_GAD_7_vs_HADS,
                                          show_top_n = 1000,
                                          order_by = "default",
                                          plot_title = "Specificity (Sp) Differences: Screening Thresholds")
## 
plot_screening_Se_and_Sp_GAD_7_vs_HADS <- plot_screening_Se + plot_screening_Sp
plot_screening_Se_and_Sp_GAD_7_vs_HADS
##
ggsave("DUMMY_Application_Se_and_Sp_pairwise_diffs_screening_GAD_GAD_7_vs_HADS.png", 
       plot_screening_Se_and_Sp_GAD_7_vs_HADS, 
       width = 16, height = 16, dpi = plot_dpi)



##
## ---- GAD-2 vs. BAI:
##
test_pair <- c("GAD-2", "BAI")
screening_comps_GAD_2_vs_BAI <- dplyr::filter(screening_comps_2, (test1_name %in% test_pair) & (test2_name %in% test_pair))
##
plot_screening_Se <- plot_se_differences( filtered_tibble = screening_comps_GAD_2_vs_BAI,
                                          show_top_n = 1000,
                                          order_by = "default",
                                          plot_title = "Sensitivity (Se) Differences: Screening Thresholds")
##
plot_screening_Sp <- plot_sp_differences( filtered_tibble = screening_comps_GAD_2_vs_BAI,
                                          show_top_n = 1000,
                                          order_by = "default",
                                          plot_title = "Specificity (Sp) Differences: Screening Thresholds")
## 
plot_screening_Se_and_Sp_GAD_2_vs_BAI <- plot_screening_Se + plot_screening_Sp
plot_screening_Se_and_Sp_GAD_2_vs_BAI
##
ggsave("DUMMY_Application_Se_and_Sp_pairwise_diffs_screening_GAD_GAD_2_vs_BAI.png", 
       plot_screening_Se_and_Sp_GAD_2_vs_BAI, 
       width = 16, height = 16, dpi = plot_dpi)



##
## ---- GAD-7 vs. BAI:
##
test_pair <- c("GAD-7", "BAI")
screening_comps_GAD_7_vs_BAI <- dplyr::filter(screening_comps_2, (test1_name %in% test_pair) & (test2_name %in% test_pair))
##
plot_screening_Se <- plot_se_differences( filtered_tibble = screening_comps_GAD_7_vs_BAI,
                                          show_top_n = 1000,
                                          order_by = "default",
                                          plot_title = "Sensitivity (Se) Differences: Screening Thresholds")
##
plot_screening_Sp <- plot_sp_differences( filtered_tibble = screening_comps_GAD_7_vs_BAI,
                                          show_top_n = 1000,
                                          order_by = "default",
                                          plot_title = "Specificity (Sp) Differences: Screening Thresholds")
## 
plot_screening_Se_and_Sp_GAD_7_vs_BAI <- plot_screening_Se + plot_screening_Sp
plot_screening_Se_and_Sp_GAD_7_vs_BAI
##
ggsave("DUMMY_Application_Se_and_Sp_pairwise_diffs_screening_GAD_GAD_7_vs_BAI.png", 
       plot_screening_Se_and_Sp_GAD_7_vs_BAI, 
       width = 16, height = 16, dpi = plot_dpi)



##
## ---- HADS vs. BAI:
##
test_pair <- c("HADS", "BAI")
screening_comps_HADS_vs_BAI <- dplyr::filter(screening_comps_2, (test1_name %in% test_pair) & (test2_name %in% test_pair))
##
plot_screening_Se <- plot_se_differences( filtered_tibble = screening_comps_HADS_vs_BAI,
                                          show_top_n = 1000,
                                          order_by = "default",
                                          plot_title = "Sensitivity (Se) Differences: Screening Thresholds")
##
plot_screening_Sp <- plot_sp_differences( filtered_tibble = screening_comps_HADS_vs_BAI,
                                          show_top_n = 1000,
                                          order_by = "default",
                                          plot_title = "Specificity (Sp) Differences: Screening Thresholds")
## 
plot_screening_Se_and_Sp_HADS_vs_BAI <- plot_screening_Se + plot_screening_Sp
plot_screening_Se_and_Sp_HADS_vs_BAI
##
ggsave("DUMMY_Application_Se_and_Sp_pairwise_diffs_screening_GAD_HADS_vs_BAI.png", 
       plot_screening_Se_and_Sp_HADS_vs_BAI, 
       width = 16, height = 16, dpi = plot_dpi)






# # Check what comparisons exist in the original data
# outs_NMA_comp_DUMMY %>%
#   filter(
#     test1_name %in% c("GAD-2", "GAD-7", "HADS", "BAI") &
#       test2_name %in% c("GAD-2", "GAD-7", "HADS", "BAI")
#   ) %>%
#   count(test1_name, test2_name) %>%
#   print(n = 20)
# 
# # See which specific threshold pairs exist
# screening_comps %>%
#   count(test1_name, threshold1, test2_name, threshold2) %>%
#   print(n = 100)
# 


# ##
# ##
# ## ---- Extract NMA performance metrics  -----------------------------------------------------------------------------
# ## NOTE: (not doing these for dummy analysis - only do for real if Klaus et al insist):
# ##
# # tibble_NMA_performance_metrics <- model_summary_and_trace_obj_REAL$extract_NMA_performance()
# # tibble_NMA_performance_metrics
# # ##
# # tibble_NMA_performance_metrics <- model_summary_and_trace_obj_DUMMY$extract_NMA_performance()
# # tibble_NMA_performance_metrics
# ##
# outs <- R_fn_using_Rcpp_compute_NMA_performance_posthoc(  trace_gq = model_summary_and_trace_obj_REAL$get_trace_generated_quantities(),
#                                                           test_names = test_names,
#                                                           n_thr = n_thr)
# outs
# ##
# ##
# outs <- R_fn_using_Rcpp_compute_NMA_performance_posthoc(  trace_gq = model_summary_and_trace_obj_DUMMY$get_trace_generated_quantities(),
#                                                           test_names = test_names,
#                                                           n_thr = n_thr)
# outs
















##
## ----------- Meta-regression results tables + plots: --------------------------------------------------------------------------------------------------------------------------
##

 

# # 
# #  
# ##
# all_MR_scenarios <- list(
#     ##
#     ## ---- Ref test = MINI:
#     ##
#     "Low GAD prev. (~6%), MINI"   = list( logit_prev_GAD = -1.0, 
#                                                              Ref_test_clean_SCID = 0, 
#                                                              Ref_test_clean_Structured = 0),
#     ##
#     "Median GAD prev. (~15%), MINI" = list( logit_prev_GAD =  0.0, 
#                                                              Ref_test_clean_SCID = 0, 
#                                                              Ref_test_clean_Structured = 0),
#     ##
#     "High GAD prev. (~33%), MINI"   = list( logit_prev_GAD = +1.0,
#                                                              Ref_test_clean_SCID = 0, 
#                                                              Ref_test_clean_Structured = 0),
#     ##
#     ## ---- Ref test = SCID:
#     ##
#     "Low GAD prev. (~6%), SCID"   = list( logit_prev_GAD = -1.0, 
#                                                              Ref_test_clean_SCID = 1, 
#                                                              Ref_test_clean_Structured = 0),
#     ##
#     "Median GAD prev. (~15%), SCID" = list( logit_prev_GAD =  0.0, 
#                                                              Ref_test_clean_SCID = 1, 
#                                                              Ref_test_clean_Structured = 0),
#     ##
#     "High GAD prev. (~33%), SCID"   = list( logit_prev_GAD = +1.0,
#                                                              Ref_test_clean_SCID = 1, 
#                                                              Ref_test_clean_Structured = 0),
#     ##
#     ## ---- Ref test = Structured:
#     ##
#     "Low GAD prev. (~6%), Structured"   = list( logit_prev_GAD = -1.0, 
#                                                                    Ref_test_clean_SCID = 0, 
#                                                                    Ref_test_clean_Structured = 1),
#     ##
#     "Median GAD prev. (~15%), Structured" = list( logit_prev_GAD =  0.0, 
#                                                                    Ref_test_clean_SCID = 0, 
#                                                                    Ref_test_clean_Structured = 1),
#     ##
#     "High GAD prev. (~33%), Structured"   = list( logit_prev_GAD = +1.0,
#                                                                    Ref_test_clean_SCID = 0, 
#                                                                    Ref_test_clean_Structured = 1)
# )
# ##
# baseline_case_nd_list <- baseline_case_d_list <- list()
# ##
# for (k in 1:length(all_MR_scenarios)) {
#   
#   baseline_case_d <- baseline_case_nd <- list()
#   ##
#   baseline_case_d[["GAD-2"]] <- c( 1, 
#                                    all_MR_scenarios[[k]]$logit_prev_GAD,  ## mean prev_GAD
#                                    all_MR_scenarios[[k]]$Ref_test_clean_SCID,
#                                    all_MR_scenarios[[k]]$Ref_test_clean_Structured) 
#   baseline_case_nd[["GAD-2"]] <- baseline_case_d[["GAD-2"]]
#   ##
#   baseline_case_d[["GAD-7"]] <- c( 1, 
#                                    all_MR_scenarios[[k]]$logit_prev_GAD,  ## mean prev_GAD
#                                    all_MR_scenarios[[k]]$Ref_test_clean_SCID,
#                                    all_MR_scenarios[[k]]$Ref_test_clean_Structured) 
#   baseline_case_nd[["GAD-7"]] <- baseline_case_d[["GAD-7"]]
#   ##
#   baseline_case_d[["HADS"]] <- c( 1,
#                                   all_MR_scenarios[[k]]$logit_prev_GAD,  ## mean prev_GAD
#                                   all_MR_scenarios[[k]]$Ref_test_clean_SCID,
#                                   all_MR_scenarios[[k]]$Ref_test_clean_Structured) 
#   baseline_case_nd[["HADS"]] <- baseline_case_d[["HADS"]]
#   ##
#   baseline_case_d[["BAI"]] <- c( 1,
#                                  all_MR_scenarios[[k]]$logit_prev_GAD,  ## mean prev_GAD
#                                  all_MR_scenarios[[k]]$Ref_test_clean_SCID,
#                                  all_MR_scenarios[[k]]$Ref_test_clean_Structured) 
#   baseline_case_nd[["BAI"]] <- baseline_case_d[["BAI"]]
#   ##
#   
#   baseline_case_nd_list[[k]] <-  baseline_case_nd
#   baseline_case_d_list[[k]]  <- baseline_case_d
# }
# 
# str(baseline_case_d_list)
# 



####
#### ---- MR - sROC plots ---------------------------------------------------------------------------------------
####

 


##
## ---- First compute results for all scenarios
##
all_scenario_results <- R_fn_compute_MR_multiple_scenarios(
  model_prep_obj = MR_Model_dummy$model_prep_obj,
  model_summary_and_trace_obj = MR_Model_dummy$model_summary_and_trace_obj,
  baseline_case_nd_list = baseline_case_nd_list,
  baseline_case_d_list = baseline_case_d_list,
  use_probit_link = TRUE,
  scenario_names = names(all_MR_scenarios),
  test_names = test_names,
  compute_comparisons = TRUE,
  save_results = TRUE,
  output_dir = MetaOrdDTA::metaorddta_output_path("MR_results"),
  parallel = FALSE
)

 

# Get the data
combined_data_sROC <- combine_sroc_data_scenarios(all_scenario_results)

  
# Check what scenario names you actually have
unique(combined_data_sROC$scenario)

# Or to see them sorted
sort(unique(combined_data_sROC$scenario))

combined_data_sROC <- combined_data_sROC %>%
  mutate(
    # Order scenarios the way YOU want
    scenario = factor(scenario, levels = c(
      "Low GAD prev. (~6%), MINI",
      "Low GAD prev. (~6%), SCID", 
      "Low GAD prev. (~6%), Structured",
      "Median GAD prev. (~15%), MINI",
      "Median GAD prev. (~15%), SCID",
      "Median GAD prev. (~15%), Structured", 
      "High GAD prev. (~33%), MINI",
      "High GAD prev. (~33%), SCID",
      "High GAD prev. (~33%), Structured"
    ))
  )

names(combined_data_sROC)

if ("prev_level" %in% names(combined_data_sROC)) {
  combined_data_sROC <- combined_data_sROC %>%
    mutate(
      prev_level = factor(prev_level, levels = c("Low (6%)", "Medium (15%)", "High (33%)")),
      ref_test = factor(ref_test, levels = c("MINI", "SCID", "Structured"))
    )
}

 

##
## ---- Create basic sROC plots:
##
separate_sROC_settings <- list( base_size = 14,
                                point_size = 2,
                                line_size = 0.8,
                                rows_in_legend = 3,
                                legend_position = "none",
                                add_regions = TRUE,
                                show_conf = TRUE,
                                show_pred = TRUE,
                                conf_region_alpha = 0.20,
                                pred_region_alpha = 0.10,
                                conf_region_colour = "blue", 
                                pred_region_colour = "green")
##
##
##
by_test_sROC_settings <- list(base_size = 14,
                              point_size = 2,
                              line_size = 0.8,
                              rows_in_legend = 3,
                              legend_position = "bottom",
                              add_regions = TRUE,
                              show_conf = TRUE,
                              show_pred = FALSE,
                              conf_region_alpha = 0.20,
                              pred_region_alpha = 0.10,
                              conf_region_colour = NULL, 
                              pred_region_colour = NULL)
##
by_scenario_sROC_settings <- list(base_size = 14,
                                  point_size = 2,
                                  line_size = 0.8,
                                  rows_in_legend = 2,
                                  legend_position = "bottom",
                                  add_regions = TRUE,
                                  show_conf = TRUE,
                                  show_pred = FALSE,
                                  conf_region_alpha = 0.20,
                                  pred_region_alpha = 0.10,
                                  conf_region_colour = NULL, 
                                  pred_region_colour = NULL)
##
##
##
by_test_grid_sROC_settings <- list(base_size = 14,
                                   point_size = 2,
                                   line_size = 0.8,
                                   rows_in_legend = 3,
                                   legend_position = "bottom",
                                   add_regions = FALSE,
                                   show_conf = FALSE,
                                   show_pred = FALSE,
                                   conf_region_alpha = 0.20,
                                   pred_region_alpha = 0.10,
                                   conf_region_colour = NULL, 
                                   pred_region_colour = NULL,
                                   ncol = 2)
##
by_scenario_grid_sROC_settings <- list(base_size = 12,
                                   point_size = 2,
                                   line_size = 0.8,
                                   rows_in_legend = 1,
                                   legend_position = "bottom",
                                   add_regions = FALSE,
                                   show_conf = FALSE,
                                   show_pred = FALSE,
                                   conf_region_alpha = 0.20,
                                   pred_region_alpha = 0.10,
                                   conf_region_colour = NULL, 
                                   pred_region_colour = NULL,
                                   ncol = 3)
##
##
##
by_test_grid_separate_sROC_settings <- list(base_size = 14,
                                   point_size = 2,
                                   line_size = 0.8,
                                   rows_in_legend = 3,
                                   legend_position = "bottom",
                                   add_regions = TRUE,
                                   show_conf = TRUE,
                                   show_pred = TRUE,
                                   conf_region_alpha = 0.20,
                                   pred_region_alpha = 0.10,
                                   conf_region_colour = "blue", 
                                   pred_region_colour = "green",
                                   ncol = 3)
##
by_scenario_grid_separate_sROC_settings <- list(base_size = 14,
                                       point_size = 2,
                                       line_size = 0.8,
                                       rows_in_legend = 1,
                                       legend_position = "bottom",
                                       add_regions = TRUE,
                                       show_conf = TRUE,
                                       show_pred = TRUE,
                                       conf_region_alpha = 0.20,
                                       pred_region_alpha = 0.10,
                                       conf_region_colour = "blue", 
                                       pred_region_colour = "green",
                                       ncol = 2)

##
sROC_plots <- create_sROC_plots_flexible(
              combined_data = combined_data_sROC,
              test_names = test_names,
              layout = c( "separate",
                         "by_test", "by_scenario",
                         "grid", "grid_separate"),
              ##
              by_test_sROC_settings = by_test_sROC_settings,
              by_scenario_sROC_settings = by_scenario_sROC_settings,
              ##
              separate_sROC_settings = separate_sROC_settings,
              ##
              by_test_grid_sROC_settings = by_test_grid_sROC_settings,
              by_scenario_grid_sROC_settings = by_scenario_grid_sROC_settings)


sROC_plots$separate$`Low GAD prev. (~6%), MINI_GAD-2`  ## this is OK - BUT DONT INCLUDE IN PAPER
##
sROC_plots$by_test$`by_test_GAD-2` ## regions not working - - BUT DONT INCLUDE IN PAPER
##
plot <- sROC_plots$grid$grid_by_test     ## this is OK + very useful plot - INCLUDE IN PAPER (but nicer version below!)
plot
ggsave("DUMMY_Application_MR_sROC_2x2_ALL_TESTS_ALL_SCENARIOS.png", plot, width = 16, height = 9, dpi = 100)
##
sROC_plots$by_scenario$`by_scenario_Low GAD prev. (~6%), MINI` ## this is OK - BUT DONT INCLUDE IN PAPER
##
plot <- sROC_plots$grid$grid_by_scenario ## this is OK + very useful plot - INCLUDE IN PAPER
plot
ggsave("DUMMY_Application_MR_sROC_3x3_ALL_TESTS_ALL_SCENARIOS.png", plot, width = 16, height = 9, dpi = 100)
##
plot <- sROC_plots$grid_separate$`grid_separate_test_GAD-2` ## pthis is OK + THIS IS A VERY USEFUL PLOT - INCLUDE IN PAPER
plot
ggsave("DUMMY_Application_MR_sROC_3x3_GAD_2_w_CrI_PrI.png", plot, width = 16, height = 9, dpi = 100)
##
plot <- sROC_plots$grid_separate$`grid_separate_test_GAD-7` ## pthis is OK + THIS IS A VERY USEFUL PLOT - INCLUDE IN PAPER
plot
ggsave("DUMMY_Application_MR_sROC_3x3_GAD_7_w_CrI_PrI.png", plot, width = 16, height = 9, dpi = 100)
##
plot <- sROC_plots$grid_separate$`grid_separate_test_HADS` ## pthis is OK + THIS IS A VERY USEFUL PLOT - INCLUDE IN PAPER
plot
ggsave("DUMMY_Application_MR_sROC_3x3_HADS_w_CrI_PrI.png", plot, width = 16, height = 9, dpi = 100)
##
plot <- sROC_plots$grid_separate$`grid_separate_test_BAI` ## pthis is OK + THIS IS A VERY USEFUL PLOT - INCLUDE IN PAPER
plot
ggsave("DUMMY_Application_MR_sROC_3x3_BAI_w_CrI_PrI.png", plot, width = 16, height = 9, dpi = 100)
##
# sROC_plots$`grid_separate_scenario_Low GAD prev. (~6%), Ref = MINI` ## this is OK - BUT DONT INCLUDE IN PAPER
 




##
## ---- ENHANCED MR PLOT - sROC 2x2 panel plot
##

# Define aesthetic mapping
my_aesthetic_mapping <- list(
  # Parsing function to create new variables
  parse_fn = function(data) {
    data %>%
      mutate(
        prev_level = case_when(
          grepl("Low.*prev", scenario) ~ "Low (6%)",
          grepl("Median.*prev", scenario) ~ "Medium (15%)",
          grepl("High.*prev", scenario) ~ "High (33%)",
          TRUE ~ "Unknown"
        ),
        ref_test = case_when(
          grepl("MINI", scenario) ~ "MINI",
          grepl("SCID", scenario) ~ "SCID",
          grepl("Structured", scenario) ~ "Structured",
          TRUE ~ "Unknown"
        )
      )
  },
  # Variable mappings
  color_var = "prev_level",
  linetype_var = "ref_test",
  shape_var = "ref_test",
  # Custom scales
  color_scale = scale_color_manual(
    values = c("Low (6%)" = "#1f78b4",
               "Medium (15%)" = "#33a02c", 
               "High (33%)" = "#e31a1c"),
    name = "GAD Prevalence"
  ),
  linetype_scale = scale_linetype_manual(
    values = c("MINI" = "solid", 
               "SCID" = "dashed", 
               "Structured" = "dotted"),
    name = "Reference Test"
  ),
  shape_scale = scale_shape_manual(
    values = c("MINI" = 16, 
               "SCID" = 17, 
               "Structured" = 15),
    name = "Reference Test"
  )
)




 

by_test_sROC_settings$show_conf  <- FALSE
##
sROC_plots_enhanced <- create_sROC_plots_flexible(
  combined_data = combined_data,
  test_names = test_names,
  layout = "by_test",
  aesthetic_mapping = my_aesthetic_mapping,
  ##
  by_test_sROC_settings = by_test_sROC_settings,
  by_scenario_sROC_settings = by_scenario_sROC_settings,
  ##
  separate_sROC_settings = separate_sROC_settings,
  ##
  by_test_grid_sROC_settings = by_test_grid_sROC_settings,
  by_scenario_grid_sROC_settings = by_scenario_grid_sROC_settings)
##
sROC_plots_enhanced
##


# Create grid with the enhanced plots
sROC_custom_grid <- (sROC_plots_enhanced$by_test[["by_test_GAD-2"]] + sROC_plots_enhanced$by_test[["by_test_GAD-7"]]) /
               (sROC_plots_enhanced$by_test[["by_test_HADS"]] + sROC_plots_enhanced$by_test[["by_test_BAI"]]) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Meta-regression sROC curves across all scenarios",
    subtitle = "Colored by prev, line type by reference test"
  ) & 
  theme_bw() + 
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        # nrow = 1,
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 14)) 

##
sROC_custom_grid
##
ggsave("DUMMY_Application_MR_sROC_2x2_ALL_TESTS_ALL_SCENARIOS_custom.png", 
       sROC_custom_grid, width = 16, height = 9, dpi = 100)



## Compare clarity:
sROC_custom_grid
sROC_plots$grid$grid_by_test






























 