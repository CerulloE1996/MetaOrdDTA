

require(ggplot2)
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
{
        rm(create_threshold_plots_flexible)
        rm(create_sROC_plots_flexible)
        rm(MetaOrd_model)
        ##
        MetaOrd_model <- NULL
        source(file.path("R", "R_R6_class_aaa_init_and_sample.R"))
        # Get the object's environment
        updated_methods <- c("run_k_fold_CV",
                             "extract_AUC", 
                             "extract_NMA_performance", 
                             "extract_NMA_comparisons",
                             "plot_sROC",
                             "plot_threshold",
                             "plot_pairwise")
        ##
        # Get the object's environment
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
## ---- Extract NMA Se/Sp comparison metrics for tables + plots  ------------------------------------------------------------------------
##
# draws_array <- Baseline_Model_best$model_summary_and_trace_obj$get_all_traces_list()$traces_as_arrays$draws_array
# # trace_gq <-  Baseline_Model_best$model_summary_and_trace_obj$get_all_traces_list()$traces_as_arrays$trace_gq
# # str(trace_gq)
# ##
# outs_NMA_comp <-  R_fn_using_Rcpp_compute_NMA_comparisons_posthoc(  trace_gq = draws_array,
#                                                                           n_thr = n_thr,
#                                                                           test_names = test_names)
##
##
# MetaOrdDTA:::extract_params_from_tibble_batch(debugging = FALSE,
#                                               tibble  = Baseline_Model_best$model_summary_and_trace_obj$get_summary_generated_quantities(), 
#                                               param_strings_vec = c("Se", "Sp"),
#                                               condition = "containing")
                                              
           
                                
# screening_comps <- filter_clinically_relevant_comparisons( tibble_comp = outs_NMA_comp,
#                                                            relevant_thresholds = thr_screening) 
# screening_comps
# ##
# screening_comps_2 <- filter_clinically_relevant_comparisons( tibble_comp = outs_NMA_comp,
#                                                              relevant_thresholds = thr_screening_2)
# screening_comps_2
# ##
# screening_comps <- screening_comps %>%
#   arrange(test1, threshold1,
#           test2, threshold2)
# ##
# screening_comps_2 <- screening_comps_2 %>%
#   arrange(test1, threshold1,
#           test2, threshold2)
##
# screening_comps   




##
## ---- First, extract NMA comparisons using "$extract_NMA_comparisons()" R6 class method:
##
## ---- Note: this calls an internal Rcpp (C++) function to rapidly compute all pairwise comparisons for all test/threshold combos,
## ---- using the full posterior Stan (or BayesMVP) traces - not just naively using posterior medians or means
##
outs_NMA_comp <- Baseline_Model_best$model_summary_and_trace_obj$extract_NMA_comparisons( test_names = test_names)
##
## ---- Make R list w/ thresholds of interest for each test:
##
thr_screening = list(
      "GAD-2" = c(3),                 # 3  is SCREENING standard for GAD
      "GAD-7" = c(8, 9, 10, 11),      # 10 is SCREENING standard
      "HADS"  = c(8, 9, 10, 11)       # 8 for SCREENING ("possible") (11 is "probable")
      ## "BAI"   = c(17, 18, 19)      # 16 is standard
)
# ##
# ## More thresholds here, use this for more detailed (e.g. appendix for most/all of them) pairwise plots (choose(4, 2) = 6 plots):
# ##
# thr_screening_2 = list(
#       "GAD-2" = c(3, 4),      # 3 is SCREENING standard for GAD
#       "GAD-7" = c(8, 9, 10, 11),    # 10 is SCREENING standard for GAD
#       "HADS" = c(8, 9, 10, 11, 12),    # 11 is standard
#       "BAI" = c(16, 17, 18, 19, 20)      # 16 is standard
# )
#  
##
## ---- Create forest plots for Se + Sp (for screening thresholds) using "$plot_pairwise()" R6 class method:
##
# plot_pairwise_differences(tibble_comp = outs_NMA_comp,
#                           comparison_type = "both",
#                           ##
#                           test_names = test_names,
#                           relevant_thresholds = thr_screening,
#                           ##
#                           baseline_plot_settings = baseline_plot_settings,
#                           metareg_plot_settings = NULL)
# ##
# ##
outs_NMA_comp_subset <- outs_NMA_comp %>% filter(test1_name != "BAI", test2_name != "BAI") ## filter out BAI
# outs_NMA_comp_subset %>% print(n = 1000)
##
baseline_plot_settings <- list()
baseline_plot_settings$show_values <- TRUE ## shows posterior medians and 95% CrI's for pairwise diffs
##
outs_plot_pairwise <- Baseline_Model_best$model_summary_and_trace_obj$plot_pairwise( tibble_comp = outs_NMA_comp_subset,
                                                                                     comparison_type = "both",
                                                                                     ##
                                                                                     test_names = c("GAD-2", "GAD-7", "HADS"),
                                                                                     relevant_thresholds = thr_screening,
                                                                                     ##
                                                                                     baseline_plot_settings = baseline_plot_settings,
                                                                                     metareg_plot_settings = NULL)
##
# outs_plot_pairwise$Se
# outs_plot_pairwise$Sp
outs_plot_pairwise$combined
##
## ---- Save plot(s):
##
ggsave( "DUMMY_Application_Se_and_Sp_pairwise_diffs_screening_all.png",
        outs_plot_pairwise$combined,
        width = 12, height = 12, dpi = 200)










 
##
## ---- For meta-reg: ---------------------------------------------------------------------------------------------------------------------
##
##
MR_Model_dummy <- readRDS(
  file.path("application_results_dummy_data", 
            "seed_1_dummy_data_1intercept_only_0N_covs_4param_Xurand_thr_FALSECS_TRUEhet_sigma_0_MR_model_6_applied_results.RDS"))
##

  ##
  rm(create_threshold_plots_flexible)
  rm(create_sROC_plots_flexible)
  rm(MetaOrd_model)
  rm(plot_pairwise_differences)
{
  MetaOrd_model <- NULL
  source(file.path("R", "R_R6_class_aaa_init_and_sample.R"))
  # Get the object's environment
  updated_methods <- c("run_k_fold_CV",
                       "extract_AUC", 
                       "extract_NMA_performance", 
                       "extract_NMA_comparisons",
                       "plot_sROC",
                       "plot_threshold",
                       "plot_pairwise",
                       "recompute_baseline",
                       "compute_pairwise_comparisons")
  
  # Get the object's environment
  obj_env <- environment(MR_Model_dummy$model_summary_and_trace_obj$initialize)
  obj_as_env <- as.environment(MR_Model_dummy$model_summary_and_trace_obj)
  
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
# ##
all_MR_scenarios <- list(
  ##
  ## ---- Ref test = MINI:
  ##
  "Low GAD prev. (~6%), \nMINI"   = list( logit_prev_GAD = -1.0,
                                        Ref_test_clean_SCID = 0,
                                        Ref_test_clean_Structured = 0),
  ##
  "Median GAD prev. (~15%), \nMINI" = list( logit_prev_GAD =  0.0,
                                          Ref_test_clean_SCID = 0,
                                          Ref_test_clean_Structured = 0),
  ##
  "High GAD prev. (~33%), \nMINI"   = list( logit_prev_GAD = +1.0,
                                          Ref_test_clean_SCID = 0,
                                          Ref_test_clean_Structured = 0),
  ##
  ## ---- Ref test = SCID:
  ##
  "Low GAD prev. (~6%), \nSCID"   = list( logit_prev_GAD = -1.0,
                                        Ref_test_clean_SCID = 1,
                                        Ref_test_clean_Structured = 0),
  ##
  "Median GAD prev. (~15%), \nSCID" = list( logit_prev_GAD =  0.0,
                                          Ref_test_clean_SCID = 1,
                                          Ref_test_clean_Structured = 0),
  ##
  "High GAD prev. (~33%), \nSCID"   = list( logit_prev_GAD = +1.0,
                                          Ref_test_clean_SCID = 1,
                                          Ref_test_clean_Structured = 0),
  ##
  ## ---- Ref test = Structured:
  ##
  "Low GAD prev. (~6%), \nStructured"   = list( logit_prev_GAD = -1.0,
                                              Ref_test_clean_SCID = 0,
                                              Ref_test_clean_Structured = 1),
  ##
  "Median GAD prev. (~15%), \nStructured" = list( logit_prev_GAD =  0.0,
                                                Ref_test_clean_SCID = 0,
                                                Ref_test_clean_Structured = 1),
  ##
  "High GAD prev. (~33%), \nStructured"   = list( logit_prev_GAD = +1.0,
                                                Ref_test_clean_SCID = 0,
                                                Ref_test_clean_Structured = 1)
)
##
relevant_thresholds = list( "GAD-2" = c(1:6), 
                            "GAD-7" =  c(3:18), 
                            "HADS" =  c(3:18),
                            "BAI" = c(3:38))
##
## ---- First, extract NMA comparisons using "$extract_NMA_comparisons()" R6 class method:
##
## ---- Note: this calls an internal Rcpp (C++) function to rapidly compute all pairwise comparisons for all test/threshold combos,
## ---- using the full posterior Stan (or BayesMVP) traces - not just naively using posterior medians or means
# ##
# outs_NMA_comp <- MR_Model_dummy$model_summary_and_trace_obj$extract_NMA_comparisons(test_names = test_names)
##
## ---- Make R list w/ thresholds of interest for each test:
##
thr_screening = list(
  "GAD-2" = c(3),              # 3 is SCREENING standard for GAD
  "GAD-7" = c(8, 9, 10, 11),       # 10 is SCREENING standard  
  "HADS"  = c(8, 9, 10, 11)       # 8 for SCREENING ("possible") (11 is "probable")
  ## "BAI"   = c(17, 18, 19)      # 16 is standard
)
# ##
# ## More thresholds here, use this for more detailed (e.g. appendix for most/all of them) pairwise plots (choose(4, 2) = 6 plots):
# ##
# thr_screening_2 = list(
#       "GAD-2" = c(3, 4),      # 3 is SCREENING standard for GAD
#       "GAD-7" = c(8, 9, 10, 11),    # 10 is SCREENING standard for GAD
#       "HADS" = c(8, 9, 10, 11, 12),    # 11 is standard
#       "BAI" = c(16, 17, 18, 19, 20)      # 16 is standard
# )
#  
# ## ---- Subset data (excluding BAI for these plots):
# outs_NMA_comp_subset <- outs_NMA_comp %>% filter(test1_name != "BAI", test2_name != "BAI") ## filter out BAI
##
##
## ---- First, make sure "$recompute_baseline()" is run first (runs internal fast C++ fn):
##
MR_baseline_outs <-  MR_Model_dummy$model_summary_and_trace_obj$recompute_baseline(
  test_names = test_names,
  MR_scenarios = all_MR_scenarios)
##
## ---- Then run "$compute_pairwise_comparisons()" for NMA pairwise Se/Sp estimates (runs internal fast C++ fn):
##
MR_pairwise_baseline_outs <- MR_Model_dummy$model_summary_and_trace_obj$compute_pairwise_comparisons( test_names = test_names,
                                                                                                      recompute_baseline = FALSE,
                                                                                                      MR_scenarios = all_MR_scenarios)
# str(MR_pairwise_baseline_outs$"Low GAD prev. (~6%), MINI")

## ---- Then call "$plot_pairwise()":
##
metareg_plot_settings <- list()
metareg_plot_settings$show_values <- TRUE ## shows posterior medians and 95% CrI's for pairwise diffs
metareg_plot_settings$facet_ncol <- 3 ## default is 2
metareg_plot_settings$facet_scales <- "fixed"
metareg_plot_settings$ratio_plot <- 3.0
metareg_plot_settings$ratio_text <- 2.25
##
outs_plot_pairwise <- MR_Model_dummy$model_summary_and_trace_obj$plot_pairwise(   
                                                                                 comparison_type = "both",
                                                                                 ##
                                                                                 test_names = c("GAD-2", "GAD-7", "HADS"),
                                                                                 relevant_thresholds = thr_screening,
                                                                                 ##
                                                                                 baseline_plot_settings = NULL,
                                                                                 metareg_plot_settings = metareg_plot_settings)
##
outs_plot_pairwise$Se_forest
##
##
## ---- Save plot(s):
##
ggsave("DUMMY_Application_MR_Se_pairwise_diffs_screening_all.png", 
       outs_plot_pairwise$Se_forest, 
       width = 20*0.75, height = 12*0.75, dpi = 200)
##
ggsave("DUMMY_Application_MR_Sp_pairwise_diffs_screening_all.png", 
       outs_plot_pairwise$Sp_forest, 
       width = 20*0.75, height = 12*0.75, dpi = 200)


##
## ---- DEBUGGING:
##
MR_Model_dummy$model_summary_and_trace_obj$pairwise_comparisons_MR

# Just retrieve the stored comparisons - don't recompute!
all_comparisons <- list()
for (scenario_name in names(MR_Model_dummy$model_summary_and_trace_obj$pairwise_comparisons_MR)) {
  comps <- MR_Model_dummy$model_summary_and_trace_obj$pairwise_comparisons_MR[[scenario_name]]
  comps$scenario <- scenario_name
  all_comparisons[[scenario_name]] <- comps
}
tibble_comp <- bind_rows(all_comparisons)




outs_plot_pairwise <- plot_pairwise_differences(  tibble_comp = tibble_comp,
                            comparison_type = "both",
                            ##
                            test_names = c("GAD-2", "GAD-7", "HADS"),
                            relevant_thresholds = thr_screening,
                            ##
                            baseline_plot_settings = NULL,
                            metareg_plot_settings = metareg_plot_settings)
##
outs_plot_pairwise$Se_forest
outs_plot_pairwise$Sp_forest              
outs_plot_pairwise$combined             





str(tibble_comp)





# combined <- outs_plot_pairwise$Se_forest + outs_plot_pairwise$Sp_forest +
#   plot_layout(ncol = 4)
# combined
# 
# # Then modify just the second plot in the combined layout
# combined[[2]] <- combined[[2]] + 
#   theme(
#     axis.text.y = element_blank(),
#     axis.ticks.y = element_blank()
#   )






# ##
# ## --- DEBUGGING (SANITY CHECK FOR PAIRWISE DIFFS):
# ##
# self <- MR_Model_dummy$model_summary_and_trace_obj
# ##
# # Get the first scenario's raw data
# scenario_1 <- self$baseline_outs_MR[[1]]
# 
# # Look at the summaries for GAD-2 (test 1) and GAD-7 (test 2)
# se_summaries <- scenario_1$Se_summaries
# sp_summaries <- scenario_1$Sp_summaries
# 
# str(se_summaries)
# 
# # Parse the variable column to extract test and threshold
# se_summaries_parsed <- se_summaries %>%
#   mutate(
#     test = as.integer(str_extract(variable, "(?<=\\[)\\d+(?=,)")),
#     threshold = as.integer(str_extract(variable, "(?<=,)\\d+(?=\\])"))
#   )
# 
# sp_summaries_parsed <- sp_summaries %>%
#   mutate(
#     test = as.integer(str_extract(variable, "(?<=\\[)\\d+(?=,)")),
#     threshold = as.integer(str_extract(variable, "(?<=,)\\d+(?=\\])"))
#   )
# 
# str(sp_summaries_parsed)
# 
# # Now extract values
# gad2_thr1_se <- se_summaries_parsed %>% 
#   filter(test == 1, threshold == 1) %>% 
#   pull(`50%`)
# 
# gad7_thr1_se <- se_summaries_parsed %>% 
#   filter(test == 2, threshold == 1) %>% 
#   pull(`50%`)
# 
# gad2_thr1_sp <- sp_summaries_parsed %>% 
#   filter(test == 1, threshold == 1) %>% 
#   pull(`50%`)
# 
# gad7_thr1_sp <- sp_summaries_parsed %>% 
#   filter(test == 2, threshold == 1) %>% 
#   pull(`50%`)
# 
# # Calculate differences
# diff_se_manual <- gad2_thr1_se - gad7_thr1_se
# diff_sp_manual <- gad2_thr1_sp - gad7_thr1_sp
# 
# # Compare to tibble_comp
# comp_row <- tibble_comp %>% 
#   filter(test1 == 1, test2 == 2, threshold1 == 1, threshold2 == 1, 
#          scenario == "Low GAD prev. (~6%), MINI")
# 
# cat("GAD-2 thr1 Se:", gad2_thr1_se, "\n")
# cat("GAD-7 thr1 Se:", gad7_thr1_se, "\n")
# cat("Manual diff_Se:", diff_se_manual, "\n")
# cat("Tibble diff_Se:", comp_row$diff_Se_median, "\n\n")
# diff_se_manual - comp_row$diff_Se_median
# 
# cat("GAD-2 thr1 Sp:", gad2_thr1_sp, "\n")
# cat("GAD-7 thr1 Sp:", gad7_thr1_sp, "\n")
# cat("Manual diff_Sp:", diff_sp_manual, "\n")
# cat("Tibble diff_Sp:", comp_row$diff_Sp_median, "\n") 
# diff_sp_manual - comp_row$diff_Sp_median
# 
# 












outs_NMA_comp




plot_screening_Se <- plot_se_differences( tibble = outs_NMA_comp,
                                          relevant_thresholds = thr_screening,
                                          filter = TRUE,
                                          plot_title = "Sensitivity (Se) Differences: Screening Thresholds")
plot_screening_Se
##
plot_screening_Sp <- plot_sp_differences( tibble = outs_NMA_comp_DUMMY,
                                          relevant_thresholds = thr_screening,
                                          filter = TRUE,
                                          plot_title = "Specificity (Sp) Differences: Screening Thresholds")
## 
plot_screening_Se_and_Sp_main <- plot_screening_Se$plot + plot_screening_Sp$plot
plot_screening_Se_and_Sp_main
# nrow(plot_screening_Se$plot_data)
##
ggsave("DUMMY_Application_Se_and_Sp_pairwise_diffs_screening_GAD.png", 
       plot_screening_Se_and_Sp_main, 
       width = 16, height = 16, dpi = 200)
##
data_screening_Se_and_Sp_main <- plot_screening_Se$plot_data
data_screening_Se_and_Sp_main %>% print(n=100)


##
## ---- Se/Sp pairwise diff's subsetted by test-pair combos's (choose(4, 2) = 6 plots in total):
##
test_names
##


##
## ---- GAD-2 vs. GAD-7:
##
test_pair <- c("GAD-2", "GAD-7")
screening_comps_GAD_2_vs_GAD_7 <- dplyr::filter(screening_comps, (test1_name %in% test_pair) & (test2_name %in% test_pair))
##
plot_screening_Se <- plot_se_differences( filtered_tibble = screening_comps_GAD_2_vs_GAD_7,
                                          plot_title = "Sensitivity (Se) Differences: Screening Thresholds")
##
plot_screening_Sp <- plot_sp_differences( filtered_tibble = screening_comps_GAD_2_vs_GAD_7,
                                          plot_title = "Specificity (Sp) Differences: Screening Thresholds")
## 
plot_screening_Se_and_Sp_GAD_2_vs_GAD_7 <- plot_screening_Se$plot + plot_screening_Sp$plot
plot_screening_Se_and_Sp_GAD_2_vs_GAD_7
##
ggsave("DUMMY_Application_Se_and_Sp_pairwise_diffs_screening_GAD_GAD_2_vs_GAD_7.png", 
       plot_screening_Se_and_Sp_GAD_2_vs_GAD_7, 
       width = 16, height = 16, dpi = plot_dpi)
##
data_screening_Se_and_Sp_GAD_2_vs_GAD_7 <- plot_screening_Se$plot_data
data_screening_Se_and_Sp_GAD_2_vs_GAD_7 %>% print(n=100)



##
## ---- GAD-2 vs. HADS
##
test_pair <- c("GAD-2", "HADS")
screening_comps_GAD_2_vs_HADS <- dplyr::filter(screening_comps_2, (test1_name %in% test_pair) & (test2_name %in% test_pair))
##
plot_screening_Se <- plot_se_differences( filtered_tibble = screening_comps_GAD_2_vs_HADS,
                                          plot_title = "Sensitivity (Se) Differences: Screening Thresholds")
##
plot_screening_Sp <- plot_sp_differences( filtered_tibble = screening_comps_GAD_2_vs_HADS,
                                          plot_title = "Specificity (Sp) Differences: Screening Thresholds")
## 
plot_screening_Se_and_Sp_GAD_2_vs_HADS <- plot_screening_Se$plot + plot_screening_Sp$plot
plot_screening_Se_and_Sp_GAD_2_vs_HADS
##
ggsave("DUMMY_Application_Se_and_Sp_pairwise_diffs_screening_GAD_GAD_2_vs_HADS.png", 
       plot_screening_Se_and_Sp_GAD_2_vs_HADS, 
       width = 16, height = 16, dpi = plot_dpi)
##
data_screening_Se_and_Sp_GAD_2_vs_HADS <-  plot_screening_Se$plot_data
data_screening_Se_and_Sp_GAD_2_vs_HADS %>% print(n=100)



##
## ---- GAD-7 vs. HADS:
##
test_pair <- c("GAD-7", "HADS")
screening_comps_GAD_7_vs_HADS <- dplyr::filter(screening_comps_2, (test1_name %in% test_pair) & (test2_name %in% test_pair))
##
plot_screening_Se <- plot_se_differences( filtered_tibble = screening_comps_GAD_7_vs_HADS,
                                          plot_title = "Sensitivity (Se) Differences: Screening Thresholds")
##
plot_screening_Sp <- plot_sp_differences( filtered_tibble = screening_comps_GAD_7_vs_HADS,
                                          plot_title = "Specificity (Sp) Differences: Screening Thresholds")
## 
plot_screening_Se_and_Sp_GAD_7_vs_HADS <- plot_screening_Se$plot + plot_screening_Sp$plot
plot_screening_Se_and_Sp_GAD_7_vs_HADS
##
ggsave("DUMMY_Application_Se_and_Sp_pairwise_diffs_screening_GAD_GAD_7_vs_HADS.png", 
       plot_screening_Se_and_Sp_GAD_7_vs_HADS, 
       width = 16, height = 16, dpi = plot_dpi)
##
data_screening_Se_and_Sp_GAD_7_vs_HADS <-  plot_screening_Se$plot_data
data_screening_Se_and_Sp_GAD_7_vs_HADS %>% print(n=100)


##
## ---- GAD-2 vs. BAI:
##
test_pair <- c("GAD-2", "BAI")
screening_comps_GAD_2_vs_BAI <- dplyr::filter(screening_comps_2, (test1_name %in% test_pair) & (test2_name %in% test_pair))
##
plot_screening_Se <- plot_se_differences( filtered_tibble = screening_comps_GAD_2_vs_BAI,,
                                          plot_title = "Sensitivity (Se) Differences: Screening Thresholds")
##
plot_screening_Sp <- plot_sp_differences( filtered_tibble = screening_comps_GAD_2_vs_BAI,
                                          plot_title = "Specificity (Sp) Differences: Screening Thresholds")
## 
plot_screening_Se_and_Sp_GAD_2_vs_BAI <- plot_screening_Se$plot + plot_screening_Sp$plot
plot_screening_Se_and_Sp_GAD_2_vs_BAI
##
ggsave("DUMMY_Application_Se_and_Sp_pairwise_diffs_screening_GAD_GAD_2_vs_BAI.png", 
       plot_screening_Se_and_Sp_GAD_2_vs_BAI, 
       width = 16, height = 16, dpi = plot_dpi)
##
data_screening_Se_and_Sp_GAD_2_vs_BAI <-  plot_screening_Se$plot_data
data_screening_Se_and_Sp_GAD_2_vs_BAI %>% print(n=100)


##
## ---- GAD-7 vs. BAI:
##
test_pair <- c("GAD-7", "BAI")
screening_comps_GAD_7_vs_BAI <- dplyr::filter(screening_comps_2, (test1_name %in% test_pair) & (test2_name %in% test_pair))
##
plot_screening_Se <- plot_se_differences( filtered_tibble = screening_comps_GAD_7_vs_BAI,
                                          plot_title = "Sensitivity (Se) Differences: Screening Thresholds")
##
plot_screening_Sp <- plot_sp_differences( filtered_tibble = screening_comps_GAD_7_vs_BAI,
                                          plot_title = "Specificity (Sp) Differences: Screening Thresholds")
## 
plot_screening_Se_and_Sp_GAD_7_vs_BAI <- plot_screening_Se$plot + plot_screening_Sp$plot
plot_screening_Se_and_Sp_GAD_7_vs_BAI
##
ggsave("DUMMY_Application_Se_and_Sp_pairwise_diffs_screening_GAD_GAD_7_vs_BAI.png", 
       plot_screening_Se_and_Sp_GAD_7_vs_BAI, 
       width = 16, height = 16, dpi = plot_dpi)
##
data_screening_Se_and_Sp_GAD_7_vs_BAI <-   plot_screening_Se$plot_data
data_screening_Se_and_Sp_GAD_7_vs_BAI %>% print(n=100)


##
## ---- HADS vs. BAI:
##
test_pair <- c("HADS", "BAI")
screening_comps_HADS_vs_BAI <- dplyr::filter(screening_comps_2, (test1_name %in% test_pair) & (test2_name %in% test_pair))
##
plot_screening_Se <- plot_se_differences( filtered_tibble = screening_comps_HADS_vs_BAI,
                                          plot_title = "Sensitivity (Se) Differences: Screening Thresholds")
##
plot_screening_Sp <- plot_sp_differences( filtered_tibble = screening_comps_HADS_vs_BAI,
                                          plot_title = "Specificity (Sp) Differences: Screening Thresholds")
## 
plot_screening_Se_and_Sp_HADS_vs_BAI <- plot_screening_Se$plot + plot_screening_Sp$plot
plot_screening_Se_and_Sp_HADS_vs_BAI
##
ggsave("DUMMY_Application_Se_and_Sp_pairwise_diffs_screening_GAD_HADS_vs_BAI.png", 
       plot_screening_Se_and_Sp_HADS_vs_BAI, 
       width = 16, height = 16, dpi = plot_dpi)
##
data_screening_Se_and_Sp_HADS_vs_BAI <-   plot_screening_Se$plot_data
data_screening_Se_and_Sp_HADS_vs_BAI %>% print(n=100)




## Main forest plot:
main_comps <- data_screening_Se_and_Sp_main %>%
              select(-test1, -test2,
                     -test1_name, -test2_name, 
                     -threshold1, -threshold2,
                     -diff_Se_mean, -diff_Sp_mean,
                     -diff_Se_sd, -diff_Sp_sd)
##
main_comps %>% print(n = 200)

## bind all to explore:
# all_comps <- tibble(rbind(data_screening_Se_and_Sp_GAD_2_vs_GAD_7_data, 
#                           data_screening_Se_and_Sp_GAD_2_vs_HADS_data,
#                           data_screening_Se_and_Sp_GAD_2_vs_BAI_data,
#                           ##
#                           data_screening_Se_and_Sp_GAD_7_vs_HADS_data,
#                           data_screening_Se_and_Sp_GAD_7_vs_BAI_data,
#                           ##
#                           data_screening_Se_and_Sp_HADS_vs_BAI_data) %>%
##
plot_screening_Se_and_Sp_GAD_2_vs_GAD_7
data_screening_Se_and_Sp_GAD_2_vs_GAD_7 %>%  select(comparison, threshold_pair,
                                                   diff_Se_lower, diff_Se_median, diff_Se_upper,
                                                   diff_Sp_lower, diff_Sp_median, diff_Sp_upper) %>% print(n = 100)
##
plot_screening_Se_and_Sp_GAD_2_vs_HADS
data_screening_Se_and_Sp_GAD_2_vs_HADS %>%  select(comparison, threshold_pair,
                                                         diff_Se_lower, diff_Se_median, diff_Se_upper,
                                                         diff_Sp_lower, diff_Sp_median, diff_Sp_upper) %>% print(n = 100)
##
plot_screening_Se_and_Sp_GAD_2_vs_BAI
data_screening_Se_and_Sp_GAD_2_vs_BAI  %>%  select(comparison, threshold_pair,
                                                         diff_Se_lower, diff_Se_median, diff_Se_upper,
                                                         diff_Sp_lower, diff_Sp_median, diff_Sp_upper) %>% print(n = 100)
##
plot_screening_Se_and_Sp_GAD_7_vs_HADS
data_screening_Se_and_Sp_GAD_7_vs_HADS  %>%  select(comparison, threshold_pair,
                                                         diff_Se_lower, diff_Se_median, diff_Se_upper,
                                                         diff_Sp_lower, diff_Sp_median, diff_Sp_upper) %>% print(n = 100)
##
plot_screening_Se_and_Sp_GAD_7_vs_BAI
data_screening_Se_and_Sp_GAD_7_vs_BAI %>%  select(comparison, threshold_pair,
                                                         diff_Se_lower, diff_Se_median, diff_Se_upper,
                                                         diff_Sp_lower, diff_Sp_median, diff_Sp_upper) %>% print(n = 100)
##
plot_screening_Se_and_Sp_HADS_vs_BAI
data_screening_Se_and_Sp_HADS_vs_BAI %>%  select(comparison, threshold_pair,
                                                         diff_Se_lower, diff_Se_median, diff_Se_upper,
                                                         diff_Sp_lower, diff_Sp_median, diff_Sp_upper) %>% print(n = 100)


##
main_comps %>% print(n = 200)







identify_sroc_position(main_comps)
##
identify_sroc_position(data_screening_Se_and_Sp_GAD_2_vs_GAD_7)
identify_sroc_position(data_screening_Se_and_Sp_GAD_2_vs_HADS)  
identify_sroc_position(data_screening_Se_and_Sp_GAD_2_vs_BAI)  
identify_sroc_position(data_screening_Se_and_Sp_GAD_7_vs_HADS) 
identify_sroc_position(data_screening_Se_and_Sp_GAD_7_vs_BAI)  
identify_sroc_position(data_screening_Se_and_Sp_HADS_vs_BAI) 

diff_data = data_screening_Se_and_Sp_GAD_2_vs_HADS

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
























































