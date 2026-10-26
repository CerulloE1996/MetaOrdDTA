



##
## ---- NON-meta-regression sROC plots first:  ---------------------------------------------------------------------------------------------
##
{
      rm(create_threshold_plots_flexible)
      rm(create_sROC_plots_flexible)
      # rm(MetaOrd_model)
      ##
      MetaOrd_model <- NULL
      source(file.path("R", "R_R6_class_aaa_init_and_sample.R"))
      ##
      ## Get the object's environment:
      ##
      updated_methods <- c("extract_AUC",
                           "extract_NMA_performance",
                           "extract_NMA_comparisons",
                           "plot_sROC",
                           "plot_threshold",
                           "run_k_fold_CV")
      ##
      ## Get the object's environment:
      ##
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
}
##
## ---- Extract results + plots for FINAL INTERCEPT-ONLY chosen model (i.e. FIXED-cutpoints + CS) ------------------------------------------
##
target_condition <- "AAD"
##
Baseline_Model_best <-  readRDS(
                        file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                  target_condition,
                                  "seed_2_dummy_data_0n_tests_8tar_cond_AADintercept_only_1N_covs_1param_Xurand_thr_FALSECS_TRUEhet_sigma_0_applied_results.RDS"))
##
# model_summary_and_trace_obj_DUMMY <- Model_dummy$model_summary_and_trace_obj
##
# model_summary_and_trace_obj <- Baseline_Model_best$model_summary_and_trace_obj
# model_prep_obj <- Baseline_Model_best$model_prep_obj
##
relevant_thresholds = list( "GAD-2" = c(1:6),
                            "GAD-7" = c(3:18),
                            "HADS"  = c(3:18),
                            ##
                            "BAI"    = c(3:45),
                            "OASIS"  = c(3:18),
                            "STAI_S" = c(3:45),
                            "STAI_T" = c(3:45),
                            "PROMIS" = c(1:32))
##
test_names <- c("GAD-2", "GAD-7", "HADS", "BAI", 
                "OASIS", "STAI_S", "STAI_T", "PROMIS")
names(test_names) <- test_names
##
## # max_scores <- c( 6, 21, 21, 
##                    63, 20, 60, 60, 32)
##
by_scenario_sROC_settings <- list()
by_scenario_sROC_settings$base_size  <- 20
by_scenario_sROC_settings$point_size <- 4
##
# by_scenario_sROC_settings$add_regions <- TRUE ## not recommended (too muh overlapping unless 2-3 tests maybe...)
# by_scenario_sROC_settings$show_conf <- TRUE ## not recommended (too muh overlapping unless 2-3 tests maybe...)
# by_scenario_sROC_settings$conf_region_alpha <- 0.5 ## not recommended (too muh overlapping unless 2-3 tests maybe...)
# by_scenario_sROC_settings$conf_region_colour <- NULL ## not recommended (too muh overlapping unless 2-3 tests maybe...)
##
by_scenario_grid_separate_sROC_settings <- list()
by_scenario_grid_separate_sROC_settings$base_size  <- 20
by_scenario_grid_separate_sROC_settings$point_size <- 4
# by_scenario_sROC_settings$conf_region_colour <- NULL
##
outs_plots_sROC <-  Baseline_Model_best$model_summary_and_trace_obj$plot_sROC(
                    test_names = test_names,
                    relevant_thresholds = relevant_thresholds,
                    by_scenario_sROC_settings = by_scenario_sROC_settings,
                    by_scenario_grid_separate_sROC_settings = by_scenario_grid_separate_sROC_settings)
##
outs_plots_sROC$all_tests
outs_plots_sROC$panels_with_regions
##
## Save plots:
##
ggsave("AAD_REAL_Application_baseline_sROC.png",
       outs_plots_sROC$all_tests,
       width = 16*0.75, height = 9*0.75, dpi = 200)
##
ggsave("AAD_REAL_Application_baseline_sROC_panel_w_CrI_PrI.png",
       outs_plots_sROC$panels_with_regions,
       width = 16*0.75, height = 16*0.75, dpi = 200)








##
## ---- NON-meta-regression Accuracy vs. thr plots first: -----------------------------------------------------------------------------------
##
by_scenario_threshold_settings <- list()
by_scenario_threshold_settings$base_size  <- 20
by_scenario_threshold_settings$point_size <- 4
by_scenario_threshold_settings$show_conf <- F ## not recommended (too muh overlapping unless 2-3 tests maybe...)
by_scenario_threshold_settings$show_pred <- F ## not recommended (too muh overlapping unless 2-3 tests maybe...)
by_scenario_threshold_settings$conf_region_alpha  <- 0.5 ## not recommended (too muh overlapping unless 2-3 tests maybe...)
by_scenario_threshold_settings$conf_region_colour <- NULL ## not recommended (too muh overlapping unless 2-3 tests maybe...)
##
by_scenario_grid_separate_threshold_settings <- list()
by_scenario_grid_separate_threshold_settings$base_size  <- 20
by_scenario_grid_separate_threshold_settings$point_size <- 4
# by_scenario_grid_separate_threshold_settings$conf_region_colour <- NULL
##
outs_plots_thr <- Baseline_Model_best$model_summary_and_trace_obj$plot_threshold( 
                  test_names = test_names,
                  ##
                  #  plot_type = "Se",
                  ## layout = c("by_scenario", "grid_separate"),
                  ##
                  plot_type = "both",
                  show_actual_thresholds = FALSE,
                  ##
                  relevant_thresholds = relevant_thresholds,
                  ##
                  by_scenario_threshold_settings = by_scenario_threshold_settings,
                  by_scenario_grid_separate_threshold_settings = by_scenario_grid_separate_threshold_settings)      
##
# outs_plots_thr$all_tests_combined_Se_and_Sp
# 
##
## ---- If Se and/or Sp only:
##
outs_plots_thr$all_tests_combined_Se
outs_plots_thr$all_tests_combined_Sp
## outs_plots_thr$by_test
## outs_plots_thr$Se_grid
## outs_plots_thr$Sp_grid
## outs_plots_thr$both_grid
## outs_plots_thr$panels_with_regions
outs_plots_thr$panels_with_regions_Se
#
##
## ---- If BOTH:
##
outs_plots_thr$all_tests_combined_Se_and_Sp
##
# outs_plots_thr$panels_with_regions_Se
# outs_plots_thr$panels_with_regions_Sp
outs_plots_thr$panels_with_regions_Se_and_Sp
##
##
# outs_plots_thr$grid_separate$`Se_grid_separate_scenario_Baseline (no covariates`
# outs_plots_thr$grid_separate$`Sp_grid_separate_scenario_Baseline (no covariates`
# ##
# outs_plots_thr$grid_separate$`both_joined_grid_separate_scenario_Baseline (no covariates`
##
## ---- save plots:
##
# ggsave("AAD_REAL_Application_baseline_Se_Sp_vs_thr.png",
#        outs_plots_thr$all_tests_combined_Se_and_Sp,
#        width = 16*0.75, height = 11*0.75, dpi = 200)
##
ggsave("AAD_REAL_Application_baseline_Se_vs_thr.png",
       outs_plots_thr$all_tests_combined_Se,
       width = 16*0.75, height = 11*0.75, dpi = 200)
##
ggsave("AAD_REAL_Application_baseline_Sp_vs_thr.png",
       outs_plots_thr$all_tests_combined_Sp,
       width = 16*0.75, height = 11*0.75, dpi = 200)
##
ggsave("AAD_REAL_Application_baseline_Se_Sp_vs_thr_panel_w_CrI_PrI.png", 
       outs_plots_thr$panels_with_regions_Se_and_Sp,
       width = 16*0.75, height = 12*0.75, dpi = 200)





























































