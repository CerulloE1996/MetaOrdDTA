



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
##
{
    os <- .Platform$OS.type
    
    
    if (os == "unix") {
      user_root_dir <- Sys.getenv("PWD")
    } else if (os == "windows") {
      user_root_dir <- Sys.getenv("USERPROFILE")
    }
    local_pkg_dir <- file.path(user_root_dir, "Documents/Work/PhD_work/R_packages/MetaOrdDTA")
    
    
    setwd(local_pkg_dir)
}
##
# Model_RAND <- readRDS(
#   file.path("application_results_dummy_data", 
#             "seed_1_dummy_data_1intercept_only_1N_covs_1param_Xurand_thr_TRUECS_0het_sigma_0_applied_results.RDS"))
# 
# model_summary_and_trace_obj <- Model_RAND$model_summary_and_trace_obj
# model_prep_obj <- Model_RAND$model_prep_obj

##
## ---- Extract results + plots for FINAL INTERCEPT-ONLY chosen model (i.e. FIXED-cutpoints + CS) --------------------------
##
Baseline_Model_best <- readRDS(
  file.path("application_results_dummy_data", 
            "seed_1_dummy_data_1intercept_only_1N_covs_1param_Xurand_thr_FALSECS_1het_sigma_0_applied_results.RDS"))
##
# model_summary_and_trace_obj_DUMMY <- Model_dummy$model_summary_and_trace_obj
##
model_summary_and_trace_obj <- Baseline_Model_best$model_summary_and_trace_obj
model_prep_obj <- Baseline_Model_best$model_prep_obj
##
test_names <- c("GAD-2",
                "GAD-7", 
                "HADS", 
                "BAI")
##
n_thr <- c(6, 21, 21, 63)
##
# # # 
# stan_mod_samples <- model_summary_and_trace_obj_DUMMY$internal_obj$outs_stan_sampling$stan_mod_samples
# fitted_model <- stan_mod_samples
# stan_data_list <- model_summary_and_trace_obj_DUMMY$internal_obj$outs_data$stan_data_list
# 



# Model_real <- readRDS(
#   file.path("application_results_real_data", 
#             "seed_1_dummy_data_0intercept_only_1N_covs_1param_Xurand_thr_FALSECS_1het_sigma_0_applied_results.RDS"))
# 
# model_summary_and_trace_obj_REAL <- Model_real$model_summary_and_trace_obj
# 
# 




 


##
## ---- Extract NMA AUC metrics: ----------------------------------------------------------------------------------------
##
# tibble_NMA_AUC_metrics <- model_summary_and_trace_obj_REAL$extract_AUC(new_cov_data = NULL,
#                                                                        test_names = test_names,
#                                                                        use_BayesMVP_for_faster_summaries = TRUE)
# tibble_NMA_AUC_metrics
# ##
# tibble_NMA_AUC_metrics <- model_summary_and_trace_obj_DUMMY$extract_AUC(new_cov_data = NULL,
#                                                                         test_names = test_names,
#                                                                         use_BayesMVP_for_faster_summaries = TRUE)
# tibble_NMA_AUC_metrics
# # # tibble_NMA_AUC_metrics
# # ##
# source(file.path("R", "R_fn_NMA_computations.R"))
# ##
# AUC_outs_REAL <- R_fn_extract_NMA_AUC_summary( tibble_gq = model_summary_and_trace_obj$get_summary_generated_quantities(),
#                                                n_thr = n_thr,
#                                                test_names = test_names,
#                                                network = TRUE)
# ##
# AUC_outs_REAL$auc
# AUC_outs_REAL$auc_diff
# AUC_outs_REAL$auc_pred
# # # A tibble: 4 × 7
# test test_name AUC_mean AUC_median AUC_lower AUC_upper AUC_sd
# <dbl> <chr>        <dbl>      <dbl>     <dbl>     <dbl>  <dbl>
# 1     1 GAD-2        0.834      0.834     0.800     0.865 0.0165
# 2     2 GAD-7        0.858      0.858     0.824     0.887 0.0159
# 3     3 HADS         0.846      0.846     0.813     0.876 0.0162
# 4     4 BAI          0.836      0.839     0.757     0.895 0.0352
# > AUC_outs_REAL$auc_diff
# # A tibble: 6 × 10
# test1 test2 test1_name test2_name AUC_diff_mean AUC_diff_median AUC_diff_lower AUC_diff_upper AUC_diff_sd prob_test1_better
# <int> <int> <chr>      <chr>              <dbl>           <dbl>          <dbl>          <dbl>       <dbl>             <dbl>
# 1     1     2 GAD-2      GAD-7           -0.0241         -0.0239         -0.0613         0.0126      0.0187                NA
# 2     1     3 GAD-2      HADS            -0.0119         -0.0121         -0.0548         0.0314      0.0218                NA
# 3     2     3 GAD-7      HADS            -0.00196        -0.00451        -0.0693         0.0801      0.0382                NA
# 4     1     4 GAD-2      BAI              0.0122          0.0122         -0.0305         0.0538      0.0214                NA
# 5     2     4 GAD-7      BAI              0.0221          0.0199         -0.0450         0.104       0.0380                NA
# 6     3     4 HADS       BAI              0.00995         0.00697        -0.0553         0.0912      0.0376                NA
# AUC_outs_REAL$auc_pred
# # A tibble: 4 × 7
# test test_name AUC_mean AUC_median AUC_lower AUC_upper AUC_sd
# <dbl> <chr>        <dbl>      <dbl>     <dbl>     <dbl>  <dbl>
# 1     1 GAD-2        0.808      0.823     0.588     0.946 0.0929
# 2     2 GAD-7        0.830      0.848     0.586     0.966 0.0989
# 3     3 HADS         0.826      0.840     0.607     0.954 0.0892
# 4     4 BAI          0.817      0.835     0.574     0.955 0.0986
##
# ## ---- Make LaTeX table:
# ##
# AUC_latex_table_REAL <- create_AUC_latex_table( auc_outs = AUC_outs_REAL, 
#                                                 digits = 3)
# cat(AUC_latex_table_REAL)


rm(MetaOrd_model)
 
MetaOrd_model <- NULL
source(file.path("R", "R_R6_class_aaa_init_and_sample.R"))
# Get the object's environment
updated_methods <- c("extract_AUC", 
                     "extract_NMA_performance", 
                     "extract_NMA_comparisons",
                     "plot_sROC",
                     "run_k_fold_CV")
##
## Get the object's environment
obj_env <- environment(Baseline_Model_best$model_summary_and_trace_obj$initialize)
obj_as_env <- as.environment(Baseline_Model_best$model_summary_and_trace_obj)
##
for (method in updated_methods) {
  
      # Get the new method
      new_method <- MetaOrd_model$public_methods[[method]]
      
      # Bind it to the object's environment
      environment(new_method) <- obj_env
      
      # Unlock and assign using the object as environment
      unlockBinding(method, obj_as_env)
      assign(method, new_method, envir = obj_as_env)
  
}
##
# str(model_summary_and_trace_obj$get_trace_generated_quantities())
##
## ---- For dummy applied data:
##
outs_AUC <- Baseline_Model_best$model_summary_and_trace_obj$extract_AUC( test_names = test_names,
                                                                         make_latex_table = TRUE)
outs_AUC$outs_AUC
##
## ---- Make LaTeX table:
##
# AUC_latex_table_DUMMY <- create_AUC_latex_table( auc_outs = AUC_outs_DUMMY, 
#                                                  digits = 3)
# cat(AUC_latex_table_DUMMY)
###
cat(outs_AUC$outs_AUC_latex_table)
##
##
# \begin{table}[!h]
# \centering
# \caption{AUC values, pairwise differences, and predictive intervals}
# \small
# \begin{tabular}{lcc}
# \toprule
# \multicolumn{3}{l}{\textbf{AUC estimates (95\% CI):}} \\
# \midrule
# GAD-2 & 84.5 & (80.8, 87.7) \\
# GAD-7 & 89.2 & (86.8, 91.4) \\
# HADS & 85.5 & (82.4, 88.2) \\
# BAI & 86.7 & (82.1, 90.6) \\
# \midrule
# \multicolumn{3}{l}{\textbf{Pairwise AUC differences:}} \\
# \midrule
# GAD-2 vs GAD-7 & -4.8 & (-8.4, -1.3)* \\
# GAD-2 vs HADS & -1.1 & (-5.2, 3.1) \\
# GAD-2 vs BAI & -2.2 & (-7.6, 3.3) \\
# GAD-7 vs HADS & 3.7 & (0.5, 7.0)* \\
# GAD-7 vs BAI & 2.6 & (-2.0, 7.4) \\
# HADS vs BAI & -1.2 & (-6.0, 4.0) \\
# \midrule
# \multicolumn{3}{l}{\textbf{AUC predictive intervals:}} \\
# \midrule
# GAD-2 & 82.6 & (58.3, 96.0) \\
# GAD-7 & 86.8 & (64.6, 97.6) \\
# HADS & 83.3 & (58.3, 96.8) \\
# BAI & 84.4 & (59.7, 97.3) \\
# \bottomrule
# \end{tabular}
# \begin{tablenotes}
# \footnotesize
# \item * indicates 95\% CI excludes zero (significant difference)
# \item Predictive intervals incorporate between-study heterogeneity
# \end{tablenotes}
# \end{table}
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
##
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
# #
# ##
# # plots_REAL <- model_summary_and_trace_obj_REAL$plot_sROC(test_names = test_names)
##
source(file.path("R", "R_fn_test_accuracy_plots_NMA.R"))
##
# 
# # relevant_thresholds <- NULL
# tibbles_provided = FALSE
# tibbles = NULL
# ##
# base_size = 16
# point_size = 3
# linewidth = 0.5
# ##
# ##
# conf_region_colour = "blue"
# pred_region_colour = "green"
# ##
# conf_region_alpha = 0.30
# pred_region_alpha = 0.15
# 
# stan_model_file_name = model_summary_and_trace_obj$internal_obj$outs_stan_model_name$stan_model_file_name
# stan_mod_samples     =  model_summary_and_trace_obj$internal_obj$outs_stan_sampling$stan_mod_samples
# 
# plots_REAL <- R_fn_sROC_plot_NMA(     stan_model_file_name =  model_summary_and_trace_obj$internal_obj$outs_stan_model_name$stan_model_file_name,
#                                       stan_mod_samples     =  model_summary_and_trace_obj$internal_obj$outs_stan_sampling$stan_mod_samples,
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
#                                       conf_region_alpha = 0.15,
#                                       pred_region_alpha = 0.30,
#                                       ##
#                                       base_size = 24,
#                                       point_size = 4.0,
#                                       linewidth = 1.0,
#                                       ##
#                                       n_index_tests = n_index_tests,
#                                       n_thr = n_thr)
# # ##
# # plots_REAL$plot_list[[1]]
# # #### plots_REAL$plot_list[[2]]
# # #### plots_REAL$plot_list[[3]]
# # #### plots_REAL$plot_list[[4]]
# # plots_REAL$plot_list[[5]]
# # #### plots_REAL$plot_list[[6]]
# # plots_REAL$plot_list[[7]]
# # plots_REAL$plot_list[[8]]
# # plots_REAL$plot_list[[9]]
# # ##
# # ggsave("REAL_Application_baseline_sROC.png", plots_REAL$plot_list[[1]], width = 16, height = 9, dpi = 400)
# # ggsave("REAL_Application_baseline_sROC_panel_w_CrI_PrI.png", plots_REAL$plot_list[[6]], width = 16, height = 9, dpi = 400)
# # ##
# # ggsave("REAL_Application_baseline_Se_Sp_vs_thr.png", plots_REAL$plot_list[[7]], width = 16, height = 16, dpi = 400)
# # ggsave("REAL_Application_baseline_Se_Sp_vs_thr_w_CrI.png", plots_REAL$plot_list[[8]], width = 16, height = 16, dpi = 400)
# # ggsave("REAL_Application_baseline_Se_Sp_vs_thr_panel_w_CrI_PrI.png", plots_REAL$plot_list[[9]], width = 16, height = 16, dpi = 400)
# # ##
# ##
# ##
# # plots_DUMMY <- model_summary_and_trace_obj_DUMMY$plot_sROC(test_names = test_names)
# ##
# n_thr <- c(6, 21, 21, 63)
# ##
# plots_DUMMY <- R_fn_sROC_plot_NMA(     stan_model_file_name = model_summary_and_trace_obj_DUMMY$internal_obj$outs_stan_model_name$stan_model_file_name,
#                                        stan_mod_samples     = model_summary_and_trace_obj_DUMMY$internal_obj$outs_stan_sampling$stan_mod_samples,
#                                        ##
#                                        test_names = test_names,
#                                        ##
#                                        relevant_thresholds = relevant_thresholds,
#                                        ##
#                                        # covariates = covariates,
#                                        # new_cov_data = new_cov_data,
#                                        ##
#                                        df_true = NULL,
#                                        ##
#                                        conf_region_colour = "blue", 
#                                        pred_region_colour = "green",
#                                        ##
#                                        conf_region_alpha = 0.20,
#                                        pred_region_alpha = 0.10,
#                                        ##
#                                        base_size = 24,
#                                        point_size = 4.0,
#                                        linewidth = 1.0,
#                                        ##
#                                        n_index_tests = n_index_tests,
#                                        n_thr = n_thr)
# ##
# plots_DUMMY$plot_list[[1]]
# #### plots_DUMMY$plot_list[[2]]
# #### plots_DUMMY$plot_list[[3]]
# #### plots_DUMMY$plot_list[[4]]
# plots_DUMMY$plot_list[[5]]
# #### plots_DUMMY$plot_list[[6]]
# plots_DUMMY$plot_list[[7]]
# plots_DUMMY$plot_list[[8]]
# plots_DUMMY$plot_list[[9]]
# ##
# ggsave("DUMMY_Application_baseline_sROC.png", plots_DUMMY$plot_list[[1]], width = 16, height = 9, dpi = 400)
# ggsave("DUMMY_Application_baseline_sROC_panel_w_CrI_PrI.png", plots_DUMMY$plot_list[[5]], width = 16, height = 16, dpi = 400)
# ##
# ggsave("DUMMY_Application_baseline_Se_Sp_vs_thr.png", plots_DUMMY$plot_list[[7]], width = 16, height = 16, dpi = 400)
# ggsave("DUMMY_Application_baseline_Se_Sp_vs_thr_w_CrI.png", plots_DUMMY$plot_list[[8]], width = 16, height = 16, dpi = 400)
# ggsave("DUMMY_Application_baseline_Se_Sp_vs_thr_panel_w_CrI_PrI.png", plots_DUMMY$plot_list[[9]], width = 16, height = 16, dpi = 400)
# # ##
# # ## ---- Compare real to simulated/dummy data (for personal interest NOT FOR PAPER/THESIS):
# # ##
# # require(patchwork)
# # plots_DUMMY$plot_list[[1]] + plots_REAL$plot_list[[1]]
# 
# 
# 
# 
# 
# 










##
## ----------- Meta-regression results tables + plots: --------------------------------------------------------------------------------------------------------------------------
##

 

# # 
# #  
# ##
all_MR_scenarios <- list(
    ##
    ## ---- Ref test = MINI:
    ##
    "Low GAD prev. (~6%), MINI"   = list( logit_prev_GAD = -1.0,
                                                             Ref_test_clean_SCID = 0,
                                                             Ref_test_clean_Structured = 0),
    ##
    "Median GAD prev. (~15%), MINI" = list( logit_prev_GAD =  0.0,
                                                             Ref_test_clean_SCID = 0,
                                                             Ref_test_clean_Structured = 0),
    ##
    "High GAD prev. (~33%), MINI"   = list( logit_prev_GAD = +1.0,
                                                             Ref_test_clean_SCID = 0,
                                                             Ref_test_clean_Structured = 0),
    ##
    ## ---- Ref test = SCID:
    ##
    "Low GAD prev. (~6%), SCID"   = list( logit_prev_GAD = -1.0,
                                                             Ref_test_clean_SCID = 1,
                                                             Ref_test_clean_Structured = 0),
    ##
    "Median GAD prev. (~15%), SCID" = list( logit_prev_GAD =  0.0,
                                                             Ref_test_clean_SCID = 1,
                                                             Ref_test_clean_Structured = 0),
    ##
    "High GAD prev. (~33%), SCID"   = list( logit_prev_GAD = +1.0,
                                                             Ref_test_clean_SCID = 1,
                                                             Ref_test_clean_Structured = 0),
    ##
    ## ---- Ref test = Structured:
    ##
    "Low GAD prev. (~6%), Structured"   = list( logit_prev_GAD = -1.0,
                                                                   Ref_test_clean_SCID = 0,
                                                                   Ref_test_clean_Structured = 1),
    ##
    "Median GAD prev. (~15%), Structured" = list( logit_prev_GAD =  0.0,
                                                                   Ref_test_clean_SCID = 0,
                                                                   Ref_test_clean_Structured = 1),
    ##
    "High GAD prev. (~33%), Structured"   = list( logit_prev_GAD = +1.0,
                                                                   Ref_test_clean_SCID = 0,
                                                                   Ref_test_clean_Structured = 1)
)
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
# # 



####
#### ---- MR - sROC plots ---------------------------------------------------------------------------------------
####


##
## ----------- Plot sROC curve(s) (NMA - for best model WITH covariates): ---------------------------------------------
##
# MR_Model_dummy <- readRDS(
#   file.path("application_results_dummy_data", 
#             "seed_1_dummy_data_1intercept_only_0N_covs_4param_Xurand_thr_FALSECS_TRUEhet_sigma_0_MR_model_6_applied_results.RDS"))
# ##
# model_summary_and_trace_obj  <- MR_Model_dummy$model_summary_and_trace_obj

# model_summary_and_trace_obj_META_REG_DUMMY <- MR_Model_best$model_summary_and_trace_obj
# 
# MR_Model_dummy_X <- subset_X_covariates( X = sim_results$X_mat,
#                                          covariates_to_keep = c("intercept",
#                                                                 "logit_prev_GAD",
#                                                                 "Ref_test_clean_SCID", "Ref_test_clean_Structured"))
 


# ##
# ## ---- First compute results for all scenarios
# ##
# all_scenario_results <- R_fn_compute_MR_multiple_scenarios(
#   model_prep_obj = MR_Model_best$model_prep_obj,
#   model_summary_and_trace_obj = MR_Model_best$model_summary_and_trace_obj,
#   baseline_case_nd_list = baseline_case_nd_list,
#   baseline_case_d_list = baseline_case_d_list,
#   use_probit_link = TRUE,
#   scenario_names = names(all_MR_scenarios),
#   test_names = test_names,
#   compute_comparisons = TRUE,
#   save_results = TRUE,
#   output_dir = "MR_results",
#   parallel = FALSE
# )
# 
#  
# 
# # Get the data
# combined_data <- combine_sroc_data_scenarios(all_scenario_results)
# 
#   
# # Check what scenario names you actually have
# unique(combined_data$scenario)
# 
# # Or to see them sorted
# sort(unique(combined_data$scenario))
# 
# combined_data <- combined_data %>%
#   mutate(
#     # Order scenarios the way YOU want
#     scenario = factor(scenario, levels = c(
#       "Low GAD prev. (~6%), MINI",
#       "Low GAD prev. (~6%), SCID", 
#       "Low GAD prev. (~6%), Structured",
#       "Median GAD prev. (~15%), MINI",
#       "Median GAD prev. (~15%), SCID",
#       "Median GAD prev. (~15%), Structured", 
#       "High GAD prev. (~33%), MINI",
#       "High GAD prev. (~33%), SCID",
#       "High GAD prev. (~33%), Structured"
#     ))
#   )
# 
# names(combined_data)
# 
# if ("prev_level" %in% names(combined_data)) {
#   combined_data <- combined_data %>%
#     mutate(
#       prev_level = factor(prev_level, levels = c("Low (6%)", "Medium (15%)", "High (33%)")),
#       ref_test = factor(ref_test, levels = c("MINI", "SCID", "Structured"))
#     )
# }
# 
#  

##
## ---- Create basic sROC plots:  --------------------------------------------------------------
##
# ##
# plots <- create_sROC_plots_flexible(
#               combined_data = combined_data,
#               test_names = test_names,
#               layout = c( "separate",
#                          "by_test", "by_scenario",
#                          "grid", "grid_separate"),
#               ##
#               by_test_sROC_settings = NULL,
#               by_scenario_sROC_settings = NULL,
#               ##
#               separate_sROC_settings = NULL,
#               ##
#               by_test_grid_sROC_settings = NULL,
#               by_scenario_grid_sROC_settings = NULL)
# 
# 
# 
# 

MR_Model_dummy <- readRDS(
  file.path("application_results_dummy_data", 
            "seed_1_dummy_data_1intercept_only_0N_covs_4param_Xurand_thr_FALSECS_TRUEhet_sigma_0_MR_model_6_applied_results.RDS"))
##
model_summary_and_trace_obj  <- MR_Model_dummy$model_summary_and_trace_obj



##
## ---- Extract results + plots for FINAL INTERCEPT-ONLY chosen model (i.e. FIXED-cutpoints + CS) -----------------------------------
##
Baseline_Model_best <- readRDS(
  file.path("application_results_dummy_data", 
            "seed_1_dummy_data_1intercept_only_1N_covs_1param_Xurand_thr_FALSECS_1het_sigma_0_applied_results.RDS"))
##
# model_summary_and_trace_obj_DUMMY <- Model_dummy$model_summary_and_trace_obj
##
model_summary_and_trace_obj <- Baseline_Model_best$model_summary_and_trace_obj
model_prep_obj <- Baseline_Model_best$model_prep_obj
##
test_names <- c("GAD-2", "GAD-7", "HADS", "BAI")
n_thr <- c(6, 21, 21, 63)


 
##
rm(create_sROC_plots_flexible)
rm(MetaOrd_model)
##
MetaOrd_model <- NULL
source(file.path("R", "R_R6_class_aaa_init_and_sample.R"))
##
updated_methods <- c("extract_AUC", 
                     "extract_NMA_performance", 
                     "extract_NMA_comparisons",
                     "plot_sROC",
                     "run_k_fold_CV")
##
obj_env <- environment(MR_Model_dummy$model_summary_and_trace_obj$initialize)
obj_as_env <- as.environment(MR_Model_dummy$model_summary_and_trace_obj)
##
for (method in updated_methods) {
  
      # Get the new method
      new_method <- MetaOrd_model$public_methods[[method]]
      
      # Bind it to the object's environment
      environment(new_method) <- obj_env
      
      # Unlock and assign using the object as environment
      unlockBinding(method, obj_as_env)
      assign(method, new_method, envir = obj_as_env)
  
}
##
## ---- Non-meta regression sROC plots first:  -------------------------------------------------------------------------------------
##
relevant_thresholds = list( "GAD-2" = c(1:6), 
                            "GAD-7" =  c(3:18), 
                            "HADS" =  c(3:18),
                            "BAI" = c(3:38))
##
by_scenario_sROC_settings <- list()
by_scenario_sROC_settings$base_size <- 20
by_scenario_sROC_settings$point_size <- 4
##
# by_scenario_sROC_settings$add_regions <- TRUE ## not recommended (too muh overlapping unless 2-3 tests maybe...)
# by_scenario_sROC_settings$show_conf <- TRUE ## not recommended (too muh overlapping unless 2-3 tests maybe...)
# by_scenario_sROC_settings$conf_region_alpha <- 0.5 ## not recommended (too muh overlapping unless 2-3 tests maybe...)
# by_scenario_sROC_settings$conf_region_colour <- NULL ## not recommended (too muh overlapping unless 2-3 tests maybe...)
##
by_scenario_grid_separate_sROC_settings <- list()
by_scenario_grid_separate_sROC_settings$base_size <- 20
by_scenario_grid_separate_sROC_settings$point_size <- 4
# by_scenario_sROC_settings$conf_region_colour <- NULL
##
outs_plots_sROC <- Baseline_Model_best$model_summary_and_trace_obj$plot_sROC(
                   test_names = test_names, 
                   relevant_thresholds = relevant_thresholds,
                   by_scenario_sROC_settings = by_scenario_sROC_settings,
                   by_scenario_grid_separate_sROC_settings = by_scenario_grid_separate_sROC_settings)
##
outs_plots_sROC$all_tests
outs_plots_sROC$panels_with_regions
##
## save plots:
##
ggsave("DUMMY_Application_baseline_sROC.png",
       outs_plots_sROC$all_tests,
       width = 16*0.75, height = 9*0.75, dpi = 200)
##
ggsave("DUMMY_Application_baseline_sROC_panel_w_CrI_PrI.png", 
       outs_plots_sROC$panels_with_regions,
       width = 16*0.75, height = 16*0.75, dpi = 200)






##
## ---- Meta-regression sROC plots second:  -------------------------------------------------------------------------------------
##
relevant_thresholds = list( "GAD-2" = c(1:6), 
                            "GAD-7" =  c(3:18), 
                            "HADS" =  c(3:18),
                            "BAI" = c(3:38))
##
by_scenario_sROC_settings <- list()
by_scenario_sROC_settings$base_size <- 20
by_scenario_sROC_settings$point_size <- 4
# by_scenario_sROC_settings$add_regions <- TRUE ## not recommended (too muh overlapping unless 2-3 tests maybe...)
# by_scenario_sROC_settings$show_conf <- TRUE ## not recommended (too muh overlapping unless 2-3 tests maybe...)
# by_scenario_sROC_settings$conf_region_alpha <- 0.5 ## not recommended (too muh overlapping unless 2-3 tests maybe...)
# by_scenario_sROC_settings$conf_region_colour <- NULL ## not recommended (too muh overlapping unless 2-3 tests maybe...)
##
by_scenario_grid_sROC_settings <- list()
by_scenario_grid_sROC_settings$ncol <- 3 ## change to 3 (default = 2)
##
by_scenario_grid_separate_sROC_settings <- list()
by_scenario_grid_separate_sROC_settings$base_size <- 20
by_scenario_grid_separate_sROC_settings$point_size <- 4
by_scenario_grid_separate_sROC_settings$ncol <- 3
# by_scenario_sROC_settings$conf_region_colour <- NULL
##

## ---- make sure recompute_baseline() is run first:
new_MR_baseline_outs <-  MR_Model_dummy$model_summary_and_trace_obj$recompute_baseline(
  test_names = test_names,
  MR_scenarios = all_MR_scenarios)
##
outs_plots_sROC <- MR_Model_dummy$model_summary_and_trace_obj$plot_sROC( test_names = test_names, 
                                                                         relevant_thresholds = relevant_thresholds,
                                                                         by_scenario_sROC_settings = NULL,
                                                                         by_scenario_grid_sROC_settings = by_scenario_grid_sROC_settings,
                                                                         by_scenario_grid_separate_sROC_settings = NULL)
##
# outs_plots_sROC$separate
# outs_plots_sROC$panels_with_regions
##
## save plots:
##
# 
# 
# outs_plots_sROC$separate$`Low GAD prev. (~6%), MINI_GAD-2`  ## this is OK - BUT DONT INCLUDE IN PAPER
# ##
# outs_plots_sROC$by_test$`by_test_GAD-2` ## regions not working - - BUT DONT INCLUDE IN PAPER
##
outs_plots_sROC$grid$grid_by_test     ## this is OK + very useful plot - INCLUDE IN PAPER (but nicer version below!)
ggsave("DUMMY_Application_MR_sROC_2x2_ALL_TESTS_ALL_SCENARIOS.png",
       outs_plots_sROC$grid$grid_by_test,
       width = 16, height = 9, dpi = 200)
##
# outs_plots_sROC$by_scenario$`by_scenario_Low GAD prev. (~6%), MINI` ## this is OK - BUT DONT INCLUDE IN PAPER
##
outs_plots_sROC$grid$grid_by_scenario ## this is OK + very useful plot - INCLUDE IN PAPER
ggsave("DUMMY_Application_MR_sROC_3x3_ALL_TESTS_ALL_SCENARIOS.png",
       outs_plots_sROC$grid$grid_by_scenario,
       width = 16, height = 9, dpi = 200)
##
##
##
# outs_plots_sROC$grid_separate$`grid_separate_test_GAD-2` ## pthis is OK + THIS IS A VERY USEFUL PLOT - INCLUDE IN PAPER
# ggsave("DUMMY_Application_MR_sROC_3x3_GAD_2_w_CrI_PrI.png", outs_plots_sROC$`grid_separate_test_GAD-2`, width = 16, height = 9, dpi = 100)
##
outs_plots_sROC$grid_separate$`grid_separate_test_GAD-7` ## pthis is OK + THIS IS A VERY USEFUL PLOT - INCLUDE IN PAPER
ggsave("DUMMY_Application_MR_sROC_3x3_GAD_7_w_CrI_PrI.png", 
       outs_plots_sROC$grid_separate$`grid_separate_test_GAD-7`, 
       width = 16, height = 9, dpi = 200)
# ##
# outs_plots_sROC$grid_separate$`grid_separate_test_HADS` ## pthis is OK + THIS IS A VERY USEFUL PLOT - INCLUDE IN PAPER
# ggsave("DUMMY_Application_MR_sROC_3x3_HADS_w_CrI_PrI.png", outs_plots_sROC$`grid_separate_test_HADS`, width = 16, height = 9, dpi = 100)
# ##
# outs_plots_sROC$grid_separate$`grid_separate_test_BAI` ## pthis is OK + THIS IS A VERY USEFUL PLOT - INCLUDE IN PAPER
# ggsave("DUMMY_Application_MR_sROC_3x3_BAI_w_CrI_PrI.png", outs_plots_sROC$`grid_separate_test_BAI`, width = 16, height = 9, dpi = 100)
# ##
# # outs_plots_sROC$`grid_separate_scenario_Low GAD prev. (~6%), Ref = MINI` ## this is OK - BUT DONT INCLUDE IN PAPER
#  


 

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
  color_var = "prev_level", linetype_var = "ref_test", shape_var = "ref_test",
  # Custom scales
  color_scale = scale_color_manual(
    values = c("Low (6%)" = "#1f78b4", "Medium (15%)" = "#33a02c", "High (33%)" = "#e31a1c"),
    name = "GAD Prevalence"
  ),
  linetype_scale = scale_linetype_manual(
    values = c("MINI" = "solid", "SCID" = "dashed", "Structured" = "dotted"),
    name = "Reference Test"
  ),
  shape_scale = scale_shape_manual(
    values = c("MINI" = 16,  "SCID" = 17, "Structured" = 15),
    name = "Reference Test"
  )
)

# new_MR_baseline_outs <-  MR_Model_dummy$model_summary_and_trace_obj$recompute_baseline(
#   test_names = test_names,
#   MR_scenarios = all_MR_scenarios)
##
outs_plots_sROC <- MR_Model_dummy$model_summary_and_trace_obj$plot_sROC(  test_names = test_names, 
                                                                         layout = "by_test",
                                                                         relevant_thresholds = relevant_thresholds,
                                                                         aesthetic_mapping = my_aesthetic_mapping)
##
## ---- ENHANCED MR PLOT - sROC 2x2 panel plot
##
# by_test_sROC_settings$show_conf  <- FALSE
##
scenario_results <- new_MR_baseline_outs
combined_data <- combine_sroc_data_scenarios(scenario_results)
##
outs_plots_sROC <- create_sROC_plots_flexible(
  combined_data = combined_data,
  # scenario_results = scenario_results,
  ##
  test_names = test_names,
  relevant_thresholds = relevant_thresholds,
  ##
  layout = "by_test",
  aesthetic_mapping = my_aesthetic_mapping)
  ##
  # by_test_sROC_settings = by_test_sROC_settings,
  # by_scenario_sROC_settings = by_scenario_sROC_settings,
  # ##
  # separate_sROC_settings = separate_sROC_settings,
  ##
  # by_test_grid_sROC_settings = by_test_grid_sROC_settings,
  # by_scenario_grid_sROC_settings = by_scenario_grid_sROC_settings,
  ##
  # conf_region_alpha = 0.15,
  # pred_region_alpha = 0.30,
  # ##
  # conf_region_colour = NULL, 
  # pred_region_colour = NULL)
##
outs_plots_sROC
##
plots_enhanced_by_test <- outs_plots_sROC$by_test


# Create grid with the enhanced plots:
custom_grid <- (outs_plots_sROC$by_test[["by_test_GAD-2"]] + outs_plots_sROC$by_test[["by_test_GAD-7"]]) /
               (outs_plots_sROC$by_test[["by_test_HADS"]]  + outs_plots_sROC$by_test[["by_test_BAI"]]) +
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
custom_grid
##
ggsave("DUMMY_Application_MR_sROC_2x2_ALL_TESTS_ALL_SCENARIOS_custom.png", 
       custom_grid, width = 16, height = 9,
       dpi = 200)








# 
# # First, add the prevalence and ref test variables to your combined data
# combined_data_enhanced <- combined_data %>%
#   mutate(
#     prev_level = case_when(
#       grepl("Low.*prev.", scenario) ~ "Low (6%)",
#       grepl("Median.*prev", scenario) ~ "Medium (15%)",
#       grepl("High.*prev", scenario) ~ "High (33%)",
#       TRUE ~ "Unknown"
#     ),
#     ref_test = case_when(
#       grepl("MINI", scenario) ~ "MINI",
#       grepl("SCID", scenario) ~ "SCID",
#       grepl("Structured", scenario) ~ "Structured",
#       TRUE ~ "Unknown"
#     )
#   )
# 
# 
# # Now recreate the plots with color AND linetype aesthetics
# plots_enhanced <- list()
# 
# for (t in unique(combined_data_enhanced$test)) {
# 
#   plot_data <- combined_data_enhanced %>%
#     filter(test == t)
# 
#   test_name <- ifelse(!is.null(test_names), test_names[t], paste0("Test ", t))
# 
#   p <- ggplot(plot_data,
#               aes(x = Fp_median, y = Se_median,
#                   color = prev_level,      # Color by prevalence
#                   linetype = ref_test,     # Linetype by ref test
#                   group = scenario)) +     # Still group by scenario
#     geom_line(size = 0.8) +
#     geom_point(aes(shape = ref_test), size = 2) +  # Also use shape for ref test
#     scale_color_manual(values = c("Low (6%)" = "#1f78b4",
#                                   "Medium (15%)" = "#33a02c",
#                                   "High (33%)" = "#e31a1c"),
#                        name = "GAD Prevalence") +
#     scale_linetype_manual(values = c("MINI" = "solid",
#                                      "SCID" = "dashed",
#                                      "Structured" = "dotted"),
#                           name = "Reference Test") +
#     scale_shape_manual(values = c("MINI" = 16,
#                                   "SCID" = 17,
#                                   "Structured" = 15),
#                        name = "Reference Test") +
#     scale_x_continuous(labels = scales::percent_format(accuracy = 1),
#                        limits = c(0, 1)) +
#     scale_y_continuous(labels = scales::percent_format(accuracy = 1),
#                        limits = c(0, 1)) +
#     labs(x = "False positive rate",
#          y = "Sensitivity",
#          title = test_name) +
#     theme_bw(base_size = 14)
# 
#   plots_enhanced[[paste0("by_test_", test_name)]] <- p
# 
# }
# plots_enhanced











##
## ---- sROC plots by scenario (3x3 panel if 9 scenarios)
##
# Create your own custom grid
library(patchwork)

# Customize the shared legend
custom_grid <- (
  ( plots$`by_scenario_Low GAD prev. (~6%), MINI` + 
     plots$`by_scenario_Median GAD prev. (~15%), MINI` +  
     plots$`by_scenario_High GAD prev. (~33%), MINI` ) / 
    (plots$`by_scenario_Low GAD prev. (~6%), SCID` +
       plots$`by_scenario_Median GAD prev. (~15%), SCID` +
       plots$`by_scenario_High GAD prev. (~33%), SCID`) /
    (plots$`by_scenario_Low GAD prev. (~6%), Structured` +
       plots$`by_scenario_Median GAD prev. (~15%), Structured` +
       plots$`by_scenario_High GAD prev. (~33%), Structured`)
) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Meta-regression sROC curves across all scenarios",
    subtitle = "Each panel shows one scenario with all 4 tests"
  ) & 
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10)) &
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

custom_grid

# Save any plot you want
ggsave("DUMMY_Application_MR_sROC_3x3_ALL_TESTS_ALL_SCENARIOS_custom.png",
       custom_grid, width = 12, height = 10)




####
#### ---- MR - Se/Sp diff forest plots ---------------------------------------------------------------------------------------
####



# Example usage for your specific case:
# First, add the group variables to your data if they don't exist
all_comparisons_enhanced <- all_comparisons %>%
  mutate(
    # Parse scenario names - this part is specific to your naming convention
    prev_level = case_when(
      grepl("Low.*prev.", scenario) ~ "Low (6%)",
      grepl("Median.*prev.", scenario) ~ "Medium (15%)",
      grepl("High.*prev.", scenario) ~ "High (33%)",
      TRUE ~ "Unknown"
    ),
    ref_test = case_when(
      grepl("MINI", scenario) ~ "MINI",
      grepl("SCID", scenario) ~ "SCID",
      grepl("Structured", scenario) ~ "Structured",
      TRUE ~ "Unknown"
    )
  )

# Now use the general function
forest_plot <- R_fn_plot_forest_comparisons_general(
  all_comparisons_data = all_comparisons_enhanced,
  comparison_type = "Se",
  scenario_var = "scenario",
  group_var = "prev_level",    # Color by prev.
  shape_var = "ref_test",       # Shape by reference test
  n_auto_select = 6
)

# Customize colors and shapes after the fact
forest_plot_custom <- forest_plot +
  scale_color_manual(values = c("Low (6%)" = "#1f78b4",
                                "Medium (15%)" = "#33a02c", 
                                "High (33%)" = "#e31a1c"),
                     name = "GAD Prevalence") +
  scale_shape_manual(values = c("MINI" = 16, 
                                "SCID" = 17, 
                                "Structured" = 15),
                     name = "Reference Test")


forest_plot_custom












####
#### ---- MR - Se/Sp vs. threshold plots---------------------------------------------------------------------------------------
####




# # Get the combined data for customization
combined_data <- extract_threshold_data_scenarios(all_scenario_results)


# Check what scenario names you actually have
unique(combined_data$scenario)

# Or to see them sorted
sort(unique(combined_data$scenario))

combined_data <- combined_data %>%
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

names(combined_data)

if ("prev_level" %in% names(combined_data)) {
  combined_data <- combined_data %>%
    mutate(
      prev_level = factor(prev_level, levels = c("Low (6%)", "Medium (15%)", "High (33%)")),
      ref_test = factor(ref_test, levels = c("MINI", "SCID", "Structured"))
    )
}


# 

# Add your specific variables AFTER
combined_data_enhanced <- combined_data %>%
  mutate(
    prev_level = case_when(
      grepl("Low.*prev.", scenario) ~ "Low (6%)",
      grepl("Median.*prev.", scenario) ~ "Medium (15%)",
      grepl("High.*prev.", scenario) ~ "High (33%)"
    ),
    ref_test = case_when(
      grepl("MINI", scenario) ~ "MINI",
      grepl("SCID", scenario) ~ "SCID",  
      grepl("Structured", scenario) ~ "Structured"
    ),
    test_name = test_names[test]
  )



# 
 

##
## ---- Create Se/Sp vs. threshold plots:
##
separate_threshold_settings <- list( base_size = 14,
                                point_size = 2,
                                line_size = 0.8,
                                rows_in_legend = 3,
                                legend_position = "none",
                                ##
                                add_regions = TRUE,
                                show_conf = TRUE,
                                show_pred = TRUE,
                                ##
                                conf_region_alpha = 0.20,
                                pred_region_alpha = 0.10,
                                ##
                                conf_region_colour_Se = "red", 
                                conf_region_colour_Sp = "green", 
                                ##
                                pred_region_colour_Se = "red",
                                pred_region_colour_Sp = "green")
                         
##
##
##
by_test_threshold_settings <- list(base_size = 14,
                              point_size = 2,
                              line_size = 0.8,
                              rows_in_legend = 3,
                              legend_position = "bottom",
                              ##
                              add_regions = TRUE,
                              show_conf = TRUE,
                              show_pred = FALSE,
                              ##
                              conf_region_alpha = 0.20,
                              pred_region_alpha = 0.10,
                              ##
                              conf_region_colour_Se = "red", 
                              conf_region_colour_Sp = "green", 
                              ##
                              pred_region_colour_Se = "red",
                              pred_region_colour_Sp = "green")
##
by_scenario_threshold_settings <- list(base_size = 14,
                                  point_size = 2,
                                  line_size = 0.8,
                                  rows_in_legend = 2,
                                  legend_position = "bottom",
                                  add_regions = TRUE,
                                  show_conf = TRUE,
                                  show_pred = FALSE,
                                  conf_region_alpha = 0.20,
                                  pred_region_alpha = 0.10,
                                  ##
                                  conf_region_colour_Se = "red", 
                                  conf_region_colour_Sp = "green", 
                                  ##
                                  pred_region_colour_Se = "red",
                                  pred_region_colour_Sp = "green")
##
##
##
by_test_grid_threshold_settings <- list(base_size = 14,
                                        point_size = 2,
                                        line_size = 0.8,
                                        rows_in_legend = 3,
                                        legend_position = "bottom",
                                        add_regions = FALSE,
                                        show_conf = FALSE,
                                        show_pred = FALSE,
                                        conf_region_alpha = NULL,
                                        pred_region_alpha = NULL,
                                        ##
                                        conf_region_colour_Se = "red", 
                                        conf_region_colour_Sp = "green", 
                                        ##
                                        pred_region_colour_Se = "red",
                                        pred_region_colour_Sp = "green")
##
by_scenario_grid_threshold_settings <- list(base_size = 14,
                                            point_size = 2,
                                            line_size = 0.8,
                                            rows_in_legend = 1,
                                            legend_position = "bottom",
                                            add_regions = FALSE,
                                            show_conf = FALSE,
                                            show_pred = FALSE,
                                            ##
                                            conf_region_alpha = NULL,
                                            pred_region_alpha = NULL,
                                            ##
                                            conf_region_colour_Se = "red", 
                                            conf_region_colour_Sp = "green", 
                                            ##
                                            pred_region_colour_Se = "red",
                                            pred_region_colour_Sp = "green")
##
##
##
by_test_grid_separate_threshold_settings <- list( base_size = 14,
                                                  point_size = 2,
                                                  line_size = 0.8,
                                                  rows_in_legend = 3,
                                                  legend_position = "none",
                                                  add_bands = TRUE,
                                                  show_conf = TRUE,
                                                  show_pred = TRUE,
                                                  ##
                                                  conf_region_alpha = 0.20,
                                                  pred_region_alpha = 0.10,
                                                  ##
                                                  conf_region_colour_Se = "red", 
                                                  conf_region_colour_Sp = "green", 
                                                  ##
                                                  pred_region_colour_Se = "red",
                                                  pred_region_colour_Sp = "green")
##
by_scenario_grid_separate_threshold_settings <- list( base_size = 14,
                                                      point_size = 2,
                                                      line_size = 0.8,
                                                      rows_in_legend = 1,
                                                      legend_position = "none",
                                                      add_bands = TRUE,
                                                      show_conf = TRUE,
                                                      show_pred = FALSE,
                                                      ##
                                                      conf_region_alpha = 0.20,
                                                      pred_region_alpha = 0.10,
                                                      ##
                                                      conf_region_colour_Se = "red", 
                                                      conf_region_colour_Sp = "green", 
                                                      ##
                                                      pred_region_colour_Se = "red",
                                                      pred_region_colour_Sp = "green")

##
require(grid)
##
plots <- create_threshold_plots_flexible(
            combined_data = combined_data,
            test_names = test_names,
            plot_type = c("Se", "Sp", "both"),
            layout = c( "separate",
                        "by_test", "by_scenario",
                        "grid", "grid_separate"),
            ##
            separate_threshold_settings = separate_threshold_settings,
            ##
            by_test_threshold_settings = by_test_threshold_settings,
            by_scenario_threshold_settings = by_scenario_threshold_settings,
            ##
            by_test_grid_threshold_settings = by_test_grid_threshold_settings,
            by_scenario_grid_threshold_settings = by_scenario_grid_threshold_settings,
            ##
            by_test_grid_separate_threshold_settings = by_test_grid_separate_threshold_settings,
            by_scenario_grid_separate_threshold_settings = by_scenario_grid_separate_threshold_settings,
            ##
            ## Aesthetic mapping options (same as sROC)
            aesthetic_mapping = NULL,
            ##
            show_actual_thresholds = FALSE)


plots$separate$`both_Low GAD prev. (~6%), MINI_GAD-2`  ## this is OK - BUT DONT INCLUDE IN PAPER
##
plots$by_test$`both_by_test_GAD-2` 
##
# plots$both_grid_by_test     ## this is OK + very useful plot - INCLUDE IN PAPER (but nicer version below!)
# plots$Se_grid_by_test
# plots$Sp_grid_by_test
##
both_grid_by_test <- plots$grid[["Se_grid_by_test"]] + theme(
  axis.text.x = element_blank(),
  axis.title.x = element_blank(),
  axis.ticks.x = element_blank()
) + ggtitle("Se/Sp vs. threshold by test") + 
                     plots$grid[["Sp_grid_by_test"]] + ggtitle(" ") + 
                     plot_layout(ncol = 1, guides = "collect") &
                     theme(legend.position = "bottom")
both_grid_by_test
##
ggsave("DUMMY_Application_MR_threshold_2x2_ALL_TESTS_ALL_SCENARIOS.png", both_grid_by_test, width = 16, height = 18, dpi = 100)
##
# plots$`both_by_scenario_Low GAD prev. (~6%), MINI` ## this is OK - BUT DONT INCLUDE IN PAPER
##
# plots$grid$Se_grid_by_scenario ## this is OK + very useful plot - INCLUDE IN PAPER
# plots$grid$Sp_grid_by_scenario ## this is OK + very useful plot - INCLUDE IN PAPER
plot <- plots$grid$Se_grid_by_scenario
plot
ggsave("DUMMY_Application_MR_threshold_3x3_ALL_TESTS_ALL_SCENARIOS_Se.png", plot, width = 16, height = 9, dpi = 100)
##
plot <- plots$grid$Sp_grid_by_scenario
plot
ggsave("DUMMY_Application_MR_threshold_3x3_ALL_TESTS_ALL_SCENARIOS_Sp.png", plot, width = 16, height = 9, dpi = 100)
##
plot <- plots$grid$both_joined_grid_by_scenario
plot
ggsave("DUMMY_Application_MR_threshold_3x3_ALL_TESTS_ALL_SCENARIOS_Se_and_Sp.png", plot, width = 16, height = 9, dpi = 100)
##
## ---- 3x3 grid's by test:
##
plot <- plots$grid_separate$`both_joined_grid_separate_test_GAD-2`
plot
ggsave("DUMMY_Application_MR_threshold_3x3_GAD_2_w_CrI_PrI.png", plot, width = 16, height = 9, dpi = 100)
##
plot <- plots$grid_separate$`both_joined_grid_separate_test_GAD-7`
plot
ggsave("DUMMY_Application_MR_threshold_3x3_GAD_7_w_CrI_PrI.png", plot, width = 16, height = 9, dpi = 100)
##
plot <- plots$grid_separate$`both_joined_grid_separate_test_HADS`
plot
ggsave("DUMMY_Application_MR_threshold_3x3_HADS_w_CrI_PrI.png", plot, width = 16, height = 9, dpi = 100)
##
plot <- plots$grid_separate$`both_joined_grid_separate_test_BAI`
plot
ggsave("DUMMY_Application_MR_threshold_3x3_BAI_w_CrI_PrI.png", plot, width = 16, height = 9, dpi = 100)
##
# plots$`grid_separate_scenario_Low GAD prev. (~6%), Ref = MINI` ## this is OK - BUT DONT INCLUDE IN PAPER




plots$grid_separate$`both_joined_grid_separate_test_GAD-2`
plots$grid$Se_grid_by_scenario



is.factor(combined_data$scenario)




separate
separate # r


plots$grid_separate$`Sp_grid_separate_scenario_High GAD prev. (~33%), MINI`




# # # Now customize a specific plot
# custom_se_gad2 <- ggplot(combined_data_enhanced %>% filter(test == 1),
#                          aes(x = threshold_norm, y = Se_median,
#                              color = prev_level, 
#                              linetype = ref_test,
#                              group = scenario)) +
#   geom_line(size = 0.8) +
#   geom_point(aes(shape = ref_test), size = 2) +
#   scale_color_manual(values = c("Low (5%)" = "#1f78b4",
#                                 "Medium (10%)" = "#33a02c", 
#                                 "High (22%)" = "#e31a1c"),
#                      name = "GAD Prevalence") +
#   scale_linetype_manual(values = c("MINI" = "solid", 
#                                    "SCID" = "dashed", 
#                                    "Structured" = "dotted"),
#                         name = "Reference Test") +
#   scale_shape_manual(values = c("MINI" = 16, "SCID" = 17, "Structured" = 15),
#                      name = "Reference Test") +
#   scale_y_continuous(labels = scales::percent_format(accuracy = 1),
#                      limits = c(0, 1)) +
#   labs(x = "Normalized threshold", 
#        y = "Sensitivity",
#        title = "GAD-2 Sensitivity across Meta-regression Scenarios") +
#   theme_bw(base_size = 14) +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         legend.position = "bottom")
# 
# 
# custom_se_gad2





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






















































# # Create heatmap showing all effects
# heatmap_se <- R_fn_plot_MR_comparison_heatmap(
#   all_comparisons_data = all_comparisons,
#   comparison_type = "Se",
#   n_comparisons = 30,
#   save_plot = TRUE
# )
# 
# heatmap_se





# 
 
# 
# ## If WTE reduction kicks in 1st September 2025:
# 4*(0.3/0.2) ## = 6, so get extra 2 months (so ends 28th feb rather than 31st december)
# ##
# ## If WTE reduction kicks in 1st October 2025:
# 3*(0.3/0.2) ## = 4.5 ## so get extra 14th feb
# 
# 

# 
# 
# MR_Model_dummy$model_prep_obj$internal_obj$outs_stan_compile$stan_model_file_path
# 
# 
# str(stan_mod_samples)
# 
# 
# # Method 3: Check summary
# summary_df <- stan_mod_samples$summary()
# print(head(summary_df))
# print(unique(summary_df$variable)[1:20])  # First 20 unique variable names
# 
# # Method 4: Try to see the actual error from failed chain
# # This is important - it will tell us exactly what went wrong
# gq_model$generate_quantities(
#   fitted_params = fitted_draws[1:10, ],  # Just 10 iterations
#   data = gq_data,
#   seed = 123
#   # chains = 1,
#   # show_messages = TRUE  # This might show more info
# )
# 
# # Then immediately check the output
# gq_model$output(1)
# 
# 
# 
main_trace =trace_main_as_list
nuisance_trace = trace_nuisance_as_list
str(nuisance_trace)
n_nuisance <- 1
include_nuisance

summarise <- model_summary_and_trace_obj_META_REG_DUMMY$summary(  compute_main_params = TRUE,
                                                     compute_transformed_parameters = FALSE, 
                                                     compute_generated_quantities = TRUE,
                                                     ##
                                                     save_log_lik_trace = FALSE,
                                                     ##
                                                     use_BayesMVP_for_faster_summaries = TRUE,
                                                     ##
                                                     compute_nested_rhat = FALSE)
##
summarise

##
## ----------- Plot sROC curve(s) (NMA - for best model WITH covariates): ---------------------------------------------
##
MR_Model_dummy <- readRDS(
  file.path("application_results_dummy_data", 
            "seed_1_dummy_data_1intercept_only_0N_covs_4param_Xurand_thr_FALSECS_1het_sigma_0_MR_model_6_applied_results.RDS"))

model_summary_and_trace_obj_META_REG_DUMMY <- MR_Model_dummy$model_summary_and_trace_obj
 
MR_Model_dummy_X <- subset_X_covariates( X = sim_results$X_mat,
                          covariates_to_keep = c("intercept",
                                                 "logit_prev_GAD",
                                                 "Ref_test_clean_SCID", "Ref_test_clean_Structured"))

baseline_case <- c(1, 
                   0.0,
                   0, 0)


 

stan_mod_samples     = model_summary_and_trace_obj_META_REG_DUMMY$internal_obj$outs_stan_sampling$stan_mod_samples


draws_array <- model_summary_and_trace_obj_META_REG_DUMMY$get_all_traces_list()$traces_as_arrays$draws_array
# trace_gq <-  model_summary_and_trace_obj_META_REG_DUMMY$get_all_traces_list()$traces_as_arrays$trace_gq
# str(trace_gq)
##
outs_NMA_comp_DUMMY <-  R_fn_using_Rcpp_compute_NMA_comparisons_posthoc(  trace_gq = draws_array,
                                                                          n_thr = n_thr,
                                                                          test_names = test_names)
##
str(outs_NMA_comp_DUMMY)
outs_NMA_comp_DUMMY

 

n_params_full <- n_par_inc_tp_and_gq
all_param_outs_trace <-    (BayesMVPtest:::fn_compute_param_constrain_from_trace_parallel(     unc_params_trace_input_main = trace_main_as_list,
                                                                                               unc_params_trace_input_nuisance = trace_nuisance_as_list,
                                                                                               pars_indicies_to_track = pars_indicies_to_track,
                                                                                               n_params_full = n_params_full,
                                                                                               n_nuisance = 0,
                                                                                               n_params_main = n_params_main,
                                                                                               include_nuisance = FALSE,
                                                                                               model_so_file = model_so_file,
                                                                                               json_file_path = json_file_path))




str(all_param_outs_trace)


## R_fn_using_Rcpp_compute_NMA_comparisons_posthoc
 

 
# 
# # 

# 
# # bs_model$param_constrain(theta_unc = )
# 
# prev_val <- -1.0
# 
# baseline_case_d <- baseline_case_nd <- list()
# ##
# baseline_case_d[["GAD-2"]] <- c( 1, prev_val,  ## mean prev_GAD
#                                  0, 0) ## MINI
# baseline_case_nd[["GAD-2"]] <- baseline_case_d[["GAD-2"]]
# ##
# baseline_case_d[["GAD-7"]] <- c( 1, prev_val,  ## mean prev_GAD
#                                  0, 0) ## MINI
# baseline_case_nd[["GAD-7"]] <- baseline_case_d[["GAD-7"]]
# ##
# baseline_case_d[["HADS"]] <- c( 1, prev_val,  ## mean prev_GAD
#                                 0, 0) ## MINI
# baseline_case_nd[["HADS"]] <- baseline_case_d[["HADS"]]
# ##
# baseline_case_d[["BAI"]] <- c( 1, prev_val,  ## mean prev_GAD
#                                0, 0) ## MINI
# baseline_case_nd[["BAI"]] <- baseline_case_d[["BAI"]]
# ##
# ## ---- check:
# ##
# baseline_case_d
# baseline_case_nd


model_prep_obj <- MR_Model_dummy$model_prep_obj
model_summary_and_trace_obj <- MR_Model_dummy$model_summary_and_trace_obj




##
stan_mod_samples_to_use <- R_fn_GQ_covariates_NMA(  baseline_case_nd = baseline_case_nd,
                                                    baseline_case_d = baseline_case_d,
                                                    model_summary_and_trace_obj = model_summary_and_trace_obj,
                                                    model_prep_obj = model_prep_obj)



##
## ------------- Compute new Se/Sp baselines + PI's + AUC's: -------------------------------------------------------
##
new_baseline_outs <- R_fn_using_Rcpp_compute_MR_Se_Sp_AUC_baseline(  
                                      model_prep_obj = MR_Model_dummy$model_prep_obj,
                                      model_summary_and_trace_obj = MR_Model_dummy$model_summary_and_trace_obj,
                                      baseline_case_nd = baseline_case_nd,
                                      baseline_case_d = baseline_case_d,
                                      use_probit_link = TRUE,
                                      test_names = test_names)

##
## ---- New AUC outs for new MR covariate baseline:
##
new_baseline_outs$AUC_summary
new_baseline_outs$AUC_pred_summary
new_baseline_outs$AUC_diff_summary
##
new_baseline_outs$Se_summaries %>% print(n=25)
new_baseline_outs$Sp_summaries %>% print(n=25)
##
new_baseline_outs$Se_pred_summaries %>% print(n=25)
new_baseline_outs$Sp_pred_summaries %>% print(n=25)
# ##
# ## sanity check to make sure matches stan -- ONLY WILL MATCH IF USING SAME MR-COVARIATE BASELINE!!!!
# ##
# options(scipen = 999)
# ##
# outs$Se_summaries %>% print(n=25)
# outs_Stan <- stan_mod_samples$summary(c("Se_baseline"),
#                                       mean, 
#                                       quantiles = ~ quantile(., na.rm = TRUE, probs = c(0.025, 0.50, 0.975))) %>%  print(n=25)
# ##
# outs$Se_summaries$`50%` - outs_Stan$`50%`
# ##
# outs$Se_pred_summaries %>%  print(n=25)
# outs_Stan_pred <- stan_mod_samples$summary(c("Se_baseline_pred"), 
#                                            mean, 
#                                            quantiles = ~ quantile(., na.rm = TRUE, probs = c(0.025, 0.50, 0.975))) %>%  print(n=25)
# # Compare predictive values
# outs$Se_pred_summaries$`50%` - outs_Stan_pred$`50%`

##
require(tictoc)
{
  tic()
  new_baseline_outs_orig <- R_fn_compute_MR_complete_new_baseline_analysis_v1(
                                    debugging = FALSE,
                                    model_prep_obj = MR_Model_dummy$model_prep_obj,
                                    model_summary_and_trace_obj = MR_Model_dummy$model_summary_and_trace_obj,
                                    baseline_case_nd = baseline_case_nd,
                                    baseline_case_d = baseline_case_d,
                                    use_probit_link = TRUE,
                                    test_names = test_names,
                                    compute_comparisons = TRUE)
  toc()
}

new_baseline_outs_orig$comparison_results
##
{
  tic()
  new_baseline_outs <- R_fn_compute_MR_complete_new_baseline_analysis_v2(
                                    debugging = FALSE,
                                    model_prep_obj = MR_Model_dummy$model_prep_obj,
                                    model_summary_and_trace_obj = MR_Model_dummy$model_summary_and_trace_obj,
                                    baseline_case_nd = baseline_case_nd,
                                    baseline_case_d = baseline_case_d,
                                    use_probit_link = TRUE,
                                    test_names = test_names,
                                    compute_comparisons = TRUE)
  toc()
}


new_baseline_outs$comparison_results




# plots_DUMMY <- model_summary_and_trace_obj_DUMMY$plot_sROC(test_names = test_names)
##
plots_DUMMY <- R_fn_sROC_plot_NMA(    stan_model_file_name = model_summary_and_trace_obj_DUMMY$internal_obj$outs_stan_model_name$stan_model_file_name,
                                      ##
                                      stan_mod_samples     = model_summary_and_trace_obj_DUMMY$internal_obj$outs_stan_sampling$stan_mod_samples,
                                      tibbles_provided = TRUE,
                                      tibbles = new_baseline_outs,
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

 
 
 




 
# 
# 
# str(outs$Se_baseline)
# 
# outs$Se_summaries

# # Compile and fit model
Stan_model_file_path <- model_summary_and_trace_obj_META_REG_DUMMY$internal_obj$outs_stan_compile$stan_model_file_path
Stan_functions_path <- model_summary_and_trace_obj_META_REG_DUMMY$internal_obj$outs_stan_compile$stan_functions_directory
# ##
# mod <- cmdstan_model(Stan_model_file_path,
#                      include_paths = c(Stan_functions_path))

##
# ## set Stan model file path for your Stan model (replace with your path)
Stan_model_file_path <- system.file("stan_models/basic_logistic.stan", package = "BayesMVPtest")
# make y a matrix first (matrix w/ 1 col)
### input N to use
N <- 1000
latent_results <- rlogis(n = N, location =  0.0) # for logistic regression
y <- ifelse(latent_results > 0, 1, 0)
y <- matrix(data = c(y), ncol = 1)
##
# Stan_init_list <- init_lists_per_chain[[1]]
# init_lists_per_chain <- rep(list(Stan_init_list), n_chains)
##
{
  Stan_data_list <- stan_data_list
  ##
  Stan_data_list$y <- NULL
  Stan_data_list$N <- NULL
  ##
  Stan_data_list$softplus <- 0
  ##r
  ##
  # str(Stan_data_list)
  Stan_data_list$x_with_missings <- NULL
  Stan_data_list$n_covariates_max <- NULL
  ##
  Stan_data_list_old <- Stan_data_list
  Stan_data_list <- c(Stan_data_list_old, priors)
  # str(Stan_data_list)
  Stan_data_list$cts_thr_values <- NULL
  Stan_data_list$cts_thr_values_nd <- NULL
  Stan_data_list$cts_thr_values_d <- NULL
  Stan_data_list$compound_symmetry
}
##
## set bs environment variable (otherwise it'll try downloading it even if already installed...)
bs_path <- bridgestan_path()
Sys.setenv(BRIDGESTAN = bs_path)
##
## For fixed-C:
if (random_thresholds == FALSE) {
  Stan_data_list$prior_kappa_mean <- NULL
  Stan_data_list$prior_kappa_SD <- NULL
  Stan_data_list$kappa_lb <- NULL
  Stan_data_list$n_total_C_if_random <- NULL
  Stan_data_list$prior_dirichlet_phi <- NULL
  Stan_data_list$n_total_summary_cat <- NULL
  Stan_data_list$n_total_cat <- NULL
  Stan_data_list$n_total_cutpoints <- NULL
  Stan_data_list$prior_dirichlet_alpha <- array(1.0, dim = c(n_index_tests, Stan_data_list$n_cat_max))
}

##
str(Stan_data_list)
Stan_data_list$softplus <- NULL  # Remove one instance
Stan_data_list$softplus <- 0     # Set it once
##
stanc_args <- list(paste0("--include-paths=", Stan_functions_path))
##
#  


## get_n_params_main_from_n_nuisance_using_Stan_bs_model
## get_n_params_main_from_n_nuisance_using_Stan_model

outs <- BayesMVPtest:::get_n_params_main_from_n_nuisance_using_Stan_bs_model( n_nuisance = 0,
                                                                           Stan_data_list = stan_data_list,
                                                                           Stan_model_file_path = Stan_model_file_path,
                                                                           stanc_args = stanc_args)


bs_model <- outs$bs_model

bs_model$param_constrain(theta_unc = )

# 
# n_params_main <- outs$n_params_main
# n_params <- outs$n_params
# 
# # ###  -----------  Compile + initialise the model using "MVP_model$new(...)" 
# BayesMVP_model_prep_obj <- BayesMVPtest:::MVP_model$new(   Model_type =  "Stan",
#                                             y = NULL,
#                                             N = NULL,
#                                             ##  model_args_list = model_args_list, # this arg is only needed for BUILT-IN (not Stan) models
#                                             Stan_data_list = Stan_data_list,
#                                             Stan_model_file_path = Stan_model_file_path,
#                                             init_lists_per_chain = init_lists_per_chain,
#                                             sample_nuisance = FALSE,
#                                             n_chains_burnin = 8,
#                                             ##
#                                             n_params_main = n_params_main,
# #                                             n_nuisance = 0,
# #                                             ##
# #                                             stanc_args = stanc_args)
# # 
# model_so_file <- BayesMVP_model_prep_obj$init_object$model_so_file
# json_file_path <- BayesMVP_model_prep_obj$init_object$json_file_path
# 
# # Get trace
# trace_gq <- model_summary_and_trace_obj_META_REG_DUMMY$get_trace_generated_quantities()
# trace_dims <- dim(trace_gq) ; trace_dims
# trace_vec <- as.vector(trace_gq)
# ##
# trace_main <- model_summary_and_trace_obj_META_REG_DUMMY$get_all_traces_list()$traces_as_arrays$trace_main
# trace_dims <- dim(trace_main) ; trace_dims
# trace_vec <- as.vector(trace_main)
# ##
# stan_functions_directory <- MR_Model_dummy$model_prep_obj$internal_obj$outs_stan_compile$stan_functions_directory ; stan_functions_directory
# stan_model_file_path     <- MR_Model_dummy$model_prep_obj$internal_obj$outs_stan_compile$stan_model_file_path ; stan_model_file_path
# stan_model_file_name <- MR_Model_dummy$model_prep_obj$internal_obj$outs_stan_compile$stan_model_file_name ; stan_model_file_name
# # stan_data_list           <- model_summary_and_trace_obj_META_REG_DUMMY$internal_obj$outs_data$stan_data_list 
# stan_data_list <- MR_Model_dummy$model_prep_obj$internal_obj$outs_data$stan_data_list
# str(stan_data_list)
# ##
# stan_mod_samples <- model_summary_and_trace_obj_META_REG_DUMMY$internal_obj$outs_stan_sampling$stan_mod_samples
# str(stan_mod_samples)
# ##
# new_cov_data <- list( baseline_case_nd = baseline_case_nd,
#                       baseline_case_d = baseline_case_d)
# ##
# # new_cov_data <- model_obj$cov_data
# 
# 
# 
# str(trace_main)
# 
# n_par_inc_tp_and_gq <- length(bs_model$param_names(include_tp = TRUE, include_gq = TRUE))
# 
# 
# trace_nuisance <- array(0, dim = c(500, 32, 10))
# trace_main_as_list <- trace_nuisance_as_list <- list()
# 
# for (kk in 1:32) { 
#   trace_main_as_list[[kk]] <- t(trace_main[, kk, ])
#   trace_nuisance_as_list[[kk]] <- t(trace_nuisance[, kk, ])
# }
# #
# 
# ## print(n_params_main)
# pars_indicies_to_track <- 1:n_par_inc_tp_and_gq
# n_params_full <- n_par_inc_tp_and_gq
# all_param_outs_trace <-    (BayesMVPtest:::fn_compute_param_constrain_from_trace_parallel(     unc_params_trace_input_main = trace_main_as_list,
#                                                                                                unc_params_trace_input_nuisance = trace_nuisance_as_list,
#                                                                                                pars_indicies_to_track = pars_indicies_to_track,
#                                                                                                n_params_full = n_params_full,
#                                                                                                n_nuisance = 0,
#                                                                                                n_params_main = n_params_main,
#                                                                                                include_nuisance = FALSE,
#                                                                                                model_so_file = model_so_file,
#                                                                                                json_file_path = json_file_path))
# 
# 
# 
# 
# str(all_param_outs_trace)
# 
# 
# 
# 
# 
# 
# ##
# stan_mod_samples_to_use <- R_fn_GQ_covariates_NMA(  new_cov_data = new_cov_data,
#                                                     stan_mod_samples = stan_mod_samples,
#                                                     stan_data_list = stan_data_list,
#                                                     stan_model_file_path = stan_model_file_path,
#                                                     stan_functions_directory = stan_functions_directory)
# 
# 



  
sum(MR_Model_dummy_X[[1]][[1]][, 3])/n_studies

X_summary_dummy <- R_fn_summarize_covariates_general(X_list = MR_Model_dummy_X[[1]],
                                                     test_names = test_names,
                                                     continuous_vars = c("logit_prev_GAD"),
                                                     binary_vars = c("Ref_test_clean_SCID", "Ref_test_clean_Structured")
                                                     # logit_vars = c("logit_prev_GAD")  # doesnt work!
                                                     )
                                                     
X_summary_dummy
##
# X_summary_real <- R_fn_summarize_covariates_general(X_list = real_data$cov_data$X[[1]])
# X_summary_real
##
beta_mu <- model_summary_and_trace_obj_META_REG_DUMMY$extract_params(params = c("beta_mu")) %>% print(n = 100)
##
beta_tibble
covariate_names
test_names
##


 

outs <- R_fn_map_beta_coefficients(beta_tibble = beta_mu, 
                           covariate_names = dimnames(MR_Model_dummy_X[[1]][[1]])[[2]], 
                           test_names = test_names, 
                           min_magnitude = 0.125,
                           borderline_threshold = 0.10,
                           show_borderline = TRUE)
##
## outs$full_table %>% print(n = 100)
##
outs$full_table %>%
  filter(abs(`50%`) > 0.125) %>%
  filter(is_significant == FALSE, covariate_name != "intercept") %>%
  print(n = 100)
##
outs_latex <- create_metareg_coef_table( beta_tibble = outs$full_table, 
                                         test_names = test_names)
cat(outs_latex)
# \begin{table}[!h]
# \centering
# \caption{Meta-regression coefficients (probit scale)}
# \small
# \begin{tabular}{lcc}
# \toprule
# & \multicolumn{1}{c}{Non-diseased} & \multicolumn{1}{c}{Diseased} \\
# \midrule
# \multicolumn{3}{l}{\textbf{GAD-2}} \\
# Intercept            & $\mathbf{-0.970 ~ \ci{-1.751}{-0.204}}$ & $0.295 ~ \ci{-0.483}{1.057}$ \\
# Prevalence (per SD)  & $0.089 ~ \ci{-0.023}{0.198}$ & $\mathbf{0.329 ~ \ci{0.208}{0.446}}$ \\
# Ref test: SCID       & $0.046 ~ \ci{-0.260}{0.346}$ & $0.233 ~ \ci{-0.066}{0.531}$ \\
# Ref_test_clean_Structured & $\mathbf{0.581 ~ \ci{0.123}{1.023}}$ & $-0.399 ~ \ci{-0.845}{0.040}$ \\
# %%%%
#   \midrule
# %%%%
#   \multicolumn{3}{l}{\textbf{GAD-7}} \\
# Intercept            & $\mathbf{-1.439 ~ \ci{-1.892}{-0.993}}$ & $0.109 ~ \ci{-0.350}{0.557}$ \\
# Prevalence (per SD)  & $\mathbf{0.113 ~ \ci{0.005}{0.221}}$ & $\mathbf{0.283 ~ \ci{0.166}{0.398}}$ \\
# Ref test: SCID       & $0.187 ~ \ci{-0.033}{0.405}$ & $\mathbf{0.302 ~ \ci{0.059}{0.544}}$ \\
# Ref_test_clean_Structured & $0.113 ~ \ci{-0.272}{0.486}$ & $-0.369 ~ \ci{-0.761}{0.019}$ \\
# %%%%
#   \midrule
# %%%%
#   \multicolumn{3}{l}{\textbf{HADS}} \\
# Intercept            & $\mathbf{-1.131 ~ \ci{-1.597}{-0.676}}$ & $0.197 ~ \ci{-0.275}{0.666}$ \\
# Prevalence (per SD)  & $0.012 ~ \ci{-0.158}{0.178}$ & $\mathbf{-0.317 ~ \ci{-0.499}{-0.138}}$ \\
# Ref test: SCID       & $0.030 ~ \ci{-0.249}{0.306}$ & $-0.248 ~ \ci{-0.555}{0.052}$ \\
# Ref_test_clean_Structured & $-0.148 ~ \ci{-0.554}{0.248}$ & $-0.244 ~ \ci{-0.689}{0.201}$ \\
# %%%%
#   \midrule
# %%%%
#   \multicolumn{3}{l}{\textbf{BAI}} \\
# Intercept            & $\mathbf{-1.821 ~ \ci{-2.198}{-1.442}}$ & $-0.232 ~ \ci{-0.651}{0.176}$ \\
# Prevalence (per SD)  & $0.199 ~ \ci{-0.056}{0.451}$ & $-0.116 ~ \ci{-0.399}{0.167}$ \\
# Ref test: SCID       & $0.185 ~ \ci{-0.299}{0.660}$ & $-0.228 ~ \ci{-0.763}{0.293}$ \\
# Ref_test_clean_Structured & $0.001 ~ \ci{-1.005}{0.960}$ & $-0.001 ~ \ci{-1.001}{0.967}$ \\
# %%%%
#   \bottomrule
# %%%%
#   \end{tabular}
# \begin{tablenotes}
# \item Bold indicates significant at 95\% level
# \end{tablenotes}
# \end{table}



##
## ---- Load chosen intercept-only model (i.e. the MAIN model):
##
Model_dummy <- readRDS(
  file.path("application_results_dummy_data", 
            "seed_1_dummy_data_1intercept_only_1N_covs_1param_Xurand_thr_FALSECS_1het_sigma_0_applied_results.RDS"))

model_summary_and_trace_obj_DUMMY <- Model_dummy$model_summary_and_trace_obj

 
 
{
  ##
  tibble_total_SD_inc_C_MR <- model_summary_and_trace_obj_META_REG_DUMMY$extract_params(params = c("total_SD_inc_C")) %>% 
    select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
    print(n = 100)
  ##
  tibble_beta_sigma_MR <- model_summary_and_trace_obj_META_REG_DUMMY$extract_params(params = c("beta_sigma")) %>% 
    select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
    print(n = 100)
  ##
  tibble_beta_tau_MR <- model_summary_and_trace_obj_META_REG_DUMMY$extract_params(params = c("beta_tau")) %>% 
    select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
    print(n = 100)
  ##
  tibble_beta_corr_MR <- model_summary_and_trace_obj_META_REG_DUMMY$extract_params(params = c("beta_corr")) %>% 
    select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
    print(n = 100)
  ##
  tibble_rho12_MR <- model_summary_and_trace_obj_META_REG_DUMMY$extract_params(params = c("rho12")) %>% 
    select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
    print(n = 100)
}
{
  ##
  tibble_total_SD_inc_C <- model_summary_and_trace_obj_DUMMY$extract_params(params = c("total_SD_inc_C")) %>% 
    select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
    print(n = 100)
  ##
  tibble_beta_sigma <- model_summary_and_trace_obj_DUMMY$extract_params(params = c("beta_sigma")) %>% 
    select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
    print(n = 100)
  ##
  tibble_beta_tau <- model_summary_and_trace_obj_DUMMY$extract_params(params = c("beta_tau")) %>% 
    select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
    print(n = 100)
  ##
  tibble_beta_corr <- model_summary_and_trace_obj_DUMMY$extract_params(params = c("beta_corr")) %>% 
    select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
    print(n = 100)
  ##
  tibble_rho12 <- model_summary_and_trace_obj_DUMMY$extract_params(params = c("rho12")) %>% 
    select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
    print(n = 100)
}


outs_latex <- create_heterogeneity_comparison_table( total_SD_intercept = tibble_total_SD_inc_C,
                                       beta_sigma_intercept = tibble_beta_sigma,
                                       beta_tau_intercept = tibble_beta_tau,
                                       beta_corr_intercept = tibble_beta_corr,
                                       rho12_intercept = tibble_rho12,
                                       ##
                                       total_SD_covariate = tibble_total_SD_inc_C_MR,
                                       beta_sigma_covariate = tibble_beta_sigma_MR,
                                       beta_tau_covariate = tibble_beta_tau_MR,
                                       beta_corr_covariate = tibble_beta_corr_MR,
                                       rho12_covariate = tibble_rho12_MR)
cat(outs_latex)
# \begin{table}[!h]
# \centering
# \caption{Impact of covariates on heterogeneity parameters}
# \small
# \begin{tabular}{lcc}
# \toprule
# \textbf{Parameter}  & \textbf{Model A}          & \textbf{Model A + Covariates} \\
# & \textit{(intercept-only)} & \textit{(prevalence, ref test, setting)} \\
# %%%%
#   \midrule
# %%%%
#   \multicolumn{3}{l}{\textit{\textbf{NMA shared between-study SD's:}}} \\
# $\sigma_{\beta}^{(d-)}$ & $0.215 ~ \ci{0.101}{0.304}$ & $0.207 ~ \ci{0.098}{0.292}$ \\
# $\sigma_{\beta}^{(d+)}$ & $0.421 ~ \ci{0.319}{0.531}$ & $0.358 ~ \ci{0.289}{0.439}$ \\
# %%%%
# \midrule
# %%%%
# \multicolumn{3}{l}{\textit{\textbf{NMA test-specific deviation SD's (CS structure):}}} \\
# $\tau^{(d-)}$ (all tests) & $0.287 ~ \ci{0.232}{0.352}$ & $0.271 ~ \ci{0.218}{0.338}$ \\
# $\tau^{(d+)}$ (all tests) & $0.329 ~ \ci{0.264}{0.417}$ & $0.162 ~ \ci{0.118}{0.221}$ \\
# %%%%
#   \midrule
# %%%%
#   \multicolumn{3}{l}{\textit{\textbf{Total within-test, between-study heterogeneity:}}} \\
# $\sigma_{total}^{(d-)}$ (all tests) & $0.360 ~ \ci{0.314}{0.415}$ & $0.344 ~ \ci{0.299}{0.399}$ \\
# $\sigma_{total}^{(d+)}$ (all tests) & $0.537 ~ \ci{0.468}{0.627}$ & $0.395 ~ \ci{0.337}{0.469}$ \\
# %%%%
#   \midrule
# %%%%
#   \multicolumn{3}{l}{\textit{\textbf{Correlation structure:}}} \\
# $\rho_{\beta}$ (NMA shared corr.) & $0.478 ~ \ci{0.096}{0.788}$ & $0.423 ~ \ci{0.045}{0.728}$ \\
# $\rho_{12}$ (within-test corr.) & $0.210 ~ \ci{0.037}{0.383}$ & $0.223 ~ \ci{0.024}{0.406}$ \\
# %%%%
#   \bottomrule
# %%%%
#   \end{tabular}
# %%%%
#   \begin{tablenotes}
# \footnotesize
# \item Both models use fixed cutpoints with compound symmetry correlation (Model A specification).
# \item CS = compound symmetry; all tests share the same variance components.
# \end{tablenotes}
# %%%%
#   \end{table}






 


##
## ---- Extract NMA AUC metrics: ----------------------------------------------------------------------------------------
##
 
 
# ##
# ## ---- For dummy applied data:
# ##
# baseline_case_nd <- model_summary_and_trace_obj_META_REG_DUMMY$cov_data$baseline_case_nd
# baseline_case_d <- model_summary_and_trace_obj_META_REG_DUMMY$cov_data$baseline_case_d
# ##
# baseline_case_d <- baseline_case_nd <- list()
# ##
# baseline_case_d[["GAD-2"]] <- c( 1, 0.0,  ## mean prev_GAD
#                                  0, 0) ## MINI
# baseline_case_nd[["GAD-2"]] <- baseline_case_d[["GAD-2"]]
# ##
# baseline_case_d[["GAD-7"]] <- c( 1, 0.0,  ## mean prev_GAD
#                                  0, 0) ## MINI
# baseline_case_nd[["GAD-7"]] <- baseline_case_d[["GAD-7"]]
# ##
# baseline_case_d[["HADS"]] <- c( 1, 0.0,  ## mean prev_GAD
#                                 0, 0) ## MINI
# baseline_case_nd[["HADS"]] <- baseline_case_d[["HADS"]]
# ##
# baseline_case_d[["BAI"]] <- c( 1, 0.0,  ## mean prev_GAD
#                                0, 0) ## MINI
# baseline_case_nd[["BAI"]] <- baseline_case_d[["BAI"]]
# ##
# ## ---- check:
# ##
# baseline_case_d
# baseline_case_nd

# MR_AUC_outs_DUMMY <- R_fn_using_Rcpp_compute_NMA_AUC_custom_baseline(  model_obj = model_summary_and_trace_obj_META_REG_DUMMY,
#                                                                        baseline_case_nd = baseline_case_nd,
#                                                                        baseline_case_d = baseline_case_d,
#                                                                        n_thr = n_thr,
#                                                                        test_names = test_names)

# 
# MR_AUC_outs_DUMMY <- R_fn_extract_NMA_AUC_summary(  tibble_gq = model_summary_and_trace_obj_META_REG_DUMMY$get_summary_generated_quantities(),
#                                                  n_thr = n_thr,
#                                                  test_names = test_names,
#                                                  network = network)
##
MR_AUC_outs_DUMMY$auc
MR_AUC_outs_DUMMY$auc_diff
MR_AUC_outs_DUMMY$auc_pred
# # A tibble: 4 × 7
# test test_name AUC_mean AUC_median AUC_lower AUC_upper AUC_sd
# <dbl> <chr>        <dbl>      <dbl>     <dbl>     <dbl>  <dbl>
# 1     1 GAD-2        0.785      0.785     0.745     0.821 0.0196
# 2     2 GAD-7        0.809      0.809     0.771     0.843 0.0187
# 3     3 HADS         0.873      0.873     0.843     0.900 0.0147
# 4     4 BAI          0.865      0.866     0.814     0.906 0.0234
# > MR_AUC_outs_DUMMY$auc_diff
# # A tibble: 6 × 10
# test1 test2 test1_name test2_name AUC_diff_mean AUC_diff_median AUC_diff_lower AUC_diff_upper AUC_diff_sd prob_test1_better
# <int> <int> <chr>      <chr>              <dbl>           <dbl>          <dbl>          <dbl>       <dbl>             <dbl>
# 1     1     2 GAD-2      GAD-7           -0.0240         -0.0240         -0.0721        0.0226       0.0245                NA
# 2     1     3 GAD-2      HADS            -0.0880         -0.0876         -0.136        -0.0428       0.0238                NA
# 3     2     3 GAD-7      HADS            -0.0798         -0.0803         -0.137        -0.0211       0.0300                NA
# 4     1     4 GAD-2      BAI             -0.0640         -0.0637         -0.110        -0.0203       0.0231                NA
# 5     2     4 GAD-7      BAI             -0.0559         -0.0568         -0.112         0.00300      0.0298                NA
# 6     3     4 HADS       BAI              0.00815         0.00715        -0.0420        0.0622       0.0269                NA
# MR_AUC_outs_DUMMY$auc_pred
# # A tibble: 4 × 7
# test test_name AUC_mean AUC_median AUC_lower AUC_upper AUC_sd
# <dbl> <chr>        <dbl>      <dbl>     <dbl>     <dbl>  <dbl>
# 1     1 GAD-2        0.763      0.780     0.488     0.941 0.120 
# 2     2 GAD-7        0.782      0.801     0.500     0.952 0.120 
# 3     3 HADS         0.849      0.870     0.602     0.975 0.0981
# 4     4 BAI          0.840      0.862     0.588     0.974 0.102 
##
## ---- Make LaTeX table:
##
AUC_latex_table_DUMMY <- create_AUC_latex_table( auc_outs = MR_AUC_outs_DUMMY, 
                                                 digits = 3)
cat(AUC_latex_table_DUMMY)


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
plots_REAL <- R_fn_sROC_plot_NMA(    stan_model_file_name = model_summary_and_trace_obj_REAL$internal_obj$outs_stan_model_name$stan_model_file_name,
                                     stan_mod_samples     =  model_summary_and_trace_obj_REAL$internal_obj$outs_stan_sampling$stan_mod_samples,
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
                                     conf_region_alpha = 0.30,
                                     pred_region_alpha = 0.15,
                                     ##
                                     base_size = 24,
                                     point_size = 4.0,
                                     linewidth = 1.0,
                                     ##
                                     n_index_tests = n_index_tests,
                                     n_thr = n_thr)
##
plots_REAL$plot_list[[1]]
#### plots_REAL$plot_list[[2]]
#### plots_REAL$plot_list[[3]]
#### plots_REAL$plot_list[[4]]
plots_REAL$plot_list[[5]]
#### plots_REAL$plot_list[[6]]
plots_REAL$plot_list[[7]]
plots_REAL$plot_list[[8]]
plots_REAL$plot_list[[9]]
##
ggsave("REAL_Application_baseline_sROC.png", plots_REAL$plot_list[[1]], width = 16, height = 9, dpi = 400)
ggsave("REAL_Application_baseline_sROC_panel_w_CrI_PrI.png", plots_REAL$plot_list[[6]], width = 16, height = 9, dpi = 400)
##
ggsave("REAL_Application_baseline_Se_Sp_vs_thr.png", plots_REAL$plot_list[[7]], width = 16, height = 16, dpi = 400)
ggsave("REAL_Application_baseline_Se_Sp_vs_thr_w_CrI.png", plots_REAL$plot_list[[8]], width = 16, height = 16, dpi = 400)
ggsave("REAL_Application_baseline_Se_Sp_vs_thr_panel_w_CrI_PrI.png", plots_REAL$plot_list[[9]], width = 16, height = 16, dpi = 400)
##
##
##
# plots_DUMMY <- model_summary_and_trace_obj_DUMMY$plot_sROC(test_names = test_names)
##
plots_DUMMY <- R_fn_sROC_plot_NMA(    stan_model_file_name = model_summary_and_trace_obj_DUMMY$internal_obj$outs_stan_model_name$stan_model_file_name,
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











 
 












 
# 
# 
# 
# 
# {
#   
#   
#   
#   basic_model_options = list(
#     network = NULL,
#     cts = NULL,
#     prior_only = NULL
#   )
#   ####
#   advanced_model_options = list(
#     model_parameterisation = NULL,
#     random_thresholds = NULL,
#     Dirichlet_random_effects_type = NULL,
#     box_cox = NULL,
#     softplus = NULL
#   )
#   ####
#   other_advanced_options = list(
#     advanced_compile = NULL, ## default is false (aka a "basic" compile)
#     ##
#     force_recompile = NULL,
#     quiet = NULL,
#     compile = NULL,
#     ##
#     set_custom_CXX_CPP_flags = NULL,
#     CCACHE_PATH = NULL,
#     custom_cpp_user_header_file_path = NULL,
#     CXX_COMPILER_PATH = NULL,
#     CPP_COMPILER_PATH = NULL,
#     MATH_FLAGS = NULL,
#     FMA_FLAGS = NULL,
#     AVX_FLAGS = NULL,
#     THREAD_FLAGS = NULL
#   )
#   ##
#   ####
#   priors = NULL
#   ####
#   # init_lists_per_chain = NULL,
#   ####
#   MCMC_params = list(
#     seed = NULL,
#     n_superchains = NULL,
#     n_chains = NULL,
#     n_iter = NULL,
#     n_burnin = NULL,
#     adapt_delta = NULL,
#     max_treedepth = NULL,
#     metric_shape = NULL
#   )
#   
#   ##
#   ## ---- Main "internal_obj" list:
#   ##
#   internal_obj = list(
#     outs_data = NULL,
#     outs_stan_model_name = NULL,
#     outs_stan_compile = NULL,
#     outs_stan_init = NULL,
#     outs_stan_sampling = NULL,
#     ##
#     HMC_info = NULL, ## new
#     efficiency_info = NULL, ## new
#     summaries = NULL, ## new
#     traces = NULL ## new
#   )
#   
#   debugging <- TRUE
#   ##
#   n_iter = 500
#   n_burnin = 500
#   ##
#   priors <- NULL
#   ##
#   x = x
#   n_chains = n_chains
#   ##
#   cts = FALSE
#   network = TRUE
#   ##
#   prior_only = FALSE
#   ##
#   softplus = FALSE
#   ##
#   ##
#   model_parameterisation = "R&G"
#   ##
#   random_thresholds = TRUE
#   Dirichlet_random_effects_type = "kappa"
#   ##
#   box_cox <- FALSE ## cannot input ACTUAL NA's into Stan!
#   ##
#   init_lists_per_chain = NULL
#   
#   
#   basic_model_options$network <- network
#   basic_model_options$cts <- cts
#   basic_model_options$prior_only <- prior_only
#   ##
#   advanced_model_options$model_parameterisation <- model_parameterisation
#   advanced_model_options$random_thresholds <- random_thresholds
#   advanced_model_options$Dirichlet_random_effects_type <- Dirichlet_random_effects_type
#   advanced_model_options$box_cox <- box_cox
#   advanced_model_options$softplus <- softplus
#   ##
#   MCMC_params$n_chains <- n_chains
#   
# }
# 
# 
# 
# 
# 
# outs_data = list(
#   stan_data_list = NULL,
#   n_tests = NULL,
#   n_studies = NULL,
#   n_thr = NULL,
#   n_cat = NULL
# )
# ##
# outs_stan_model_name = list(
#   stan_model_file_name = NULL
# )
# ##
# outs_stan_compile = list(
#   stan_model_obj = NULL,
#   stan_model_file_name = NULL,
#   stan_model_file_path = NULL,
#   ##
#   pkg_root_directory = NULL,
#   stan_models_directory = NULL,
#   stan_functions_directory = NULL,
#   ##
#   stan_MA_directory = NULL,
#   stan_MA_prior_directory = NULL,
#   ##
#   stan_NMA_directory = NULL,
#   stan_NMA_directory = NULL
# )
# ##
# outs_stan_init = list(
#   inits_unconstrained_vec_per_chain = NULL,
#   stan_param_names_list = NULL,
#   stan_param_names_main = NULL,
#   stan_init_pseudo_sampling_outs = NULL,
#   stan_model_obj = NULL,
#   json_file_path = NULL,
#   stan_model_file_path = NULL
# )
# ##
# outs_stan_sampling = list(
#   stan_mod_samples = NULL,
#   time_total = NULL
# )
# 
# X <- NULL
# 
# # x <- x_NMA
# 
# MetaOrdDTA:::prep_data_and_model(  
#   debugging = self$debugging,
#   ##
#   x = self$x,
#   ##
#   X = self$X, ## covariates
#   ##
#   indicator_index_test_in_study = self$indicator_index_test_in_study,
#   ##
#   internal_obj = self$internal_obj,
#   ##
#   basic_model_options    = self$basic_model_options,
#   advanced_model_options = self$advanced_model_options,
#   MCMC_params            = self$MCMC_params,
#   other_advanced_options = self$other_advanced_options,
#   ##
#   priors = self$priors,
#   init_lists_per_chain = self$init_lists_per_chain,
#   ##
#   compute_sim_study_metrics = self$compute_sim_study_metrics,
#   vec_index_inner_thr       = self$vec_index_inner_thr)
# 
# 
# 
# 
# 
# 
# 
# x = x
# X = NULL
# model_parameterisation = model_parameterisation
# n_index_tests_per_study = n_index_tests_per_study
# indicator_index_test_in_study = indicator_index_test_in_study
# 
# 
# 
# outs_data <- R_fn_prep_NMA_data(x = x,
#                                 X = NULL,
#                                 model_parameterisation = model_parameterisation,
#                                 n_index_tests_per_study = n_index_tests_per_study,
#                                 indicator_index_test_in_study = indicator_index_test_in_study)
# 
# 
# internal_obj$outs_data <- outs_data
# 



