


rm(list = ls())
.rs.restartR()  # In RStudio 
##
{
  os <- .Platform$OS.type
  ##
  if (os == "unix") {
    user_root_dir <- Sys.getenv("PWD")
  } else if (os == "windows") {
    user_root_dir <- Sys.getenv("USERPROFILE")
  }
  ##
  local_pkg_dir <- file.path(user_root_dir, "Documents/Work/PhD_work/R_packages/MetaOrdDTA")
  ##
  {
    ## First remove any possible package fragments:
    ## Find user_pkg_ins  tall_dir:
    user_pkg_install_dir <- Sys.getenv("R_LIBS_USER")
    print(paste("user_pkg_install_dir = ", user_pkg_install_dir))
    ##
    ## Find pkg_install_path + pkg_temp_install_path:
    pkg_install_path <- file.path(user_pkg_install_dir, "MetaOrdDTA")
    pkg_temp_install_path <- file.path(user_pkg_install_dir, "00LOCK-MetaOrdDTA")
    ##
    ## Remove any (possible) MetaOrdDTA package fragments:
    remove.packages("MetaOrdDTA")
    unlink(pkg_install_path, recursive = TRUE, force = TRUE)
    unlink(pkg_temp_install_path, recursive = TRUE, force = TRUE)
  }
  ##
  ## ---- ACTUAL (LOCAL) INSTALL:
  ##
  devtools::clean_dll(local_pkg_dir)
  devtools::document(local_pkg_dir)
  ## sroxygen2::roxygenize(local_pkg_dir)
  ##
  ## ---- Install:
  ##
  devtools::clean_dll(local_pkg_dir)
  devtools::install(local_pkg_dir,
                    ## build_vignettes = FALSE,
                    upgrade = "never",
                    quick = TRUE)
  ##
  ## ---- May need to restart R:
  ##
  .rs.restartR() # In RStudio
  ##
  # ?devtools::install
  ##
  require(MetaOrdDTA)
  ##
  # # Remove any lock files
  # unlink(file.path(local_pkg_dir, "00LOCK-MetaOrdDTA"), recursive = TRUE)
  #
  # # Clean compiled objects
  # pkgbuild::clean_dll(local_pkg_dir)
  #
  # # Then just load
  # devtools::load_all(local_pkg_dir)
}
##
require(MetaOrdDTA)
require(RcppParallel)
## require(BayesMVPtest)
##
require(dplyr)
##








##
## ---- NON-meta-regression Pairwise differences plots: ------------------------------------------------------------------------------------------------
##
# cat("Studies per test:", colSums(indicator_index_test_in_study), "\n")
# Studies per test: 20 29 40 4 2 9 7 3 
# >     cat("Total studies:", nrow(cov_data_Klaus_3), "\n")
# Total studies: 138
##
target_condition <- "AAD"
##
Baseline_Model_best <-  readRDS(
                        file.path("application_results_real_data",
                                  target_condition,
                                  "seed_2_dummy_data_0n_tests_8tar_cond_AADintercept_only_1N_covs_1param_Xurand_thr_FALSECS_TRUEhet_sigma_0_applied_results.RDS"))
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
# relevant_thresholds = list( "GAD-2" = c(1:6),
#                             "GAD-7" = c(3:18),
#                             "HADS"  = c(3:18),
#                             ##
#                             "BAI"    = c(3:45),
#                             "OASIS"  = c(3:18),
#                             "STAI_S" = c(3:45),
#                             "STAI_T" = c(3:45),
#                             "PROMIS" = c(1:32))
##
thr_screening_high_study_count_tests = list(
  "GAD-2" = c(3),
  "GAD-7" = c(8, 9, 10, 11),
  "HADS"  = c(8, 9, 10, 11)
)
##
## More thresholds here, use this for more detailed (e.g. appendix for most/all of them) pairwise plots (choose(4, 2) = 6 plots):
##
# thr_screening_2 = list(
#       "GAD-2" = c(3, 4),              # 3  is SCREENING standard for GAD
#       "GAD-7" = c(8, 9, 10, 11),      # 10 is SCREENING standard for GAD
#       "HADS"  = c(8, 9, 10, 11, 12),  # 11 is standard
#       "BAI"   = c(16, 17, 18, 19, 20) # 16 is standard
# )
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
##
# outs_NMA_comp_subset <- outs_NMA_comp %>% filter( test1_name != "BAI",
#                                                   test2_name != "BAI")
##
baseline_plot_settings <- list()
baseline_plot_settings$show_values <- TRUE ## shows posterior medians and 95% CrI's for pairwise diffs
##
outs_plot_pairwise_high_study_count_tests <-  Baseline_Model_best$model_summary_and_trace_obj$plot_pairwise(
                                              comparison_type = "both",
                                              ##
                                              test_names = c("GAD-2", "GAD-7", "HADS"),
                                              relevant_thresholds = thr_screening_high_study_count_tests,
                                              ##
                                              baseline_plot_settings = baseline_plot_settings,
                                              metareg_plot_settings = NULL)
##
# outs_plot_pairwise$Se
# outs_plot_pairwise$Sp
outs_plot_pairwise_high_study_count_tests$combined
##
ggsave("AAD_REAL_Application_Se_Sp_pairwise_diffs_screening_high_study_count_tests.png",
       outs_plot_pairwise_high_study_count_tests$combined,
       width = 12, height = 12, dpi = 200)
##
## -------------------------------------------------------------------------------------------------------------------------------------
##
all_test_names <- c("GAD-2", "GAD-7", "HADS", "BAI", "OASIS", "STAI_S", "STAI_T", "PROMIS")
##
thr_screening_low_study_count_tests = list(
  "BAI"    = c(8, 10, 16),   # 20 ---- 8-10 most common screening cutoffs
  "OASIS"  = c(8, 9, 10, 11, 12),   # 60 ---- 8 is standard screening, 10-12 also used
  "STAI_S" = c(20, 21, 22, 23, 24, 25),  # 60 ---- After -20 transform: screening range
  "STAI_T" = c(20, 21, 22, 23, 24, 25),  # 60 ---- Match STAI-S
  "PROMIS" = c(16, 17, 18, 19, 20, 21, 22)   # 32 ---- T-score ~60-65 equivalent for screening
)
##
# thr_screening_low_study_count_tests = list(
#   "BAI"    = c(8, 9, 10, 11, 12),   # 20 ---- 8-10 most common screening cutoffs
#   "OASIS"  = c(8, 9, 10, 11, 12)   # 60 ---- 8 is standard screening, 10-12 also used
#   # "STAI_S" = c(20, 21, 22, 23, 24, 25),  # 60 ---- After -20 transform: screening range
#   # "STAI_T" = c(20, 21, 22, 23, 24, 25),  # 60 ---- Match STAI-S
#   # "PROMIS" = c(16, 17, 18, 19, 20, 21, 22)   # 32 ---- T-score ~60-65 equivalent for screening
# )
##
outs_plot_pairwise_low_study_count_tests <-   Baseline_Model_best$model_summary_and_trace_obj$plot_pairwise(
                                              comparison_type = "both",
                                              ##
                                              # test_names = c("BAI", "OASIS"),
                                              # test_names = c("BAI", "OASIS", "STAI_S", "STAI_T", "PROMIS"),
                                              test_names = all_test_names,
                                              ##
                                              relevant_thresholds = thr_screening_low_study_count_tests,
                                              ##
                                              baseline_plot_settings = baseline_plot_settings,
                                              metareg_plot_settings = NULL)
##
outs_plot_pairwise_low_study_count_tests$combined
##
## ---- Save plot(s):
##
ggsave("AAD_REAL_Application_Se_Sp_pairwise_diffs_screening_low_study_count_tests.png",
       outs_plot_pairwise_low_study_count_tests$combined,
       width = 12, height = 12, dpi = 200)
##
## -------------------------------------------------------------------------------------------------------------------------------------
##
outs_NMA_comp
##
all_test_names <- c("GAD-2", "GAD-7", "HADS", "BAI", "OASIS", "STAI_S", "STAI_T", "PROMIS")
##
thr_screening_all_tests = list(
  "GAD-2" = c(3), ## 7 ---- 3 is standard screening thr.
  "GAD-7" = c(8, 9, 10, 11), ## 21 ---- 
  "HADS"  = c(8, 9, 10, 11), ## 21 ----
  ##
  "BAI"    = c(12, 14, 15, 16, 17, 18, 20, 24),  # 20 ---- 8-10 most common screening cutoffs
  "OASIS"  = c(6, 7, 8, 9, 10, 12),   # 60 ---- 8 is standard screening, 10-12 also used
  "STAI_S" = c(64, 68, 72, 74, 75, 76, 80) - 20, # 60 ---- After -20 transform: screening range
  "STAI_T" = c(64, 68, 72, 74, 75, 76, 80) - 20,  # 60 ---- Match STAI-S
  "PROMIS" = c(18, 20, 21, 22, 23, 24)   # 32 ---- T-score ~60-65 equivalent for screening
)
# ##
# outs_plot_pairwise_all_tests <-   Baseline_Model_best$model_summary_and_trace_obj$plot_pairwise(
#                                   comparison_type = "both",
#                                   ##
#                                   # test_names = c("BAI", "OASIS"),
#                                   # test_names = c("BAI", "OASIS", "STAI_S", "STAI_T", "PROMIS"),
#                                   test_names = all_test_names,
#                                   ##
#                                   relevant_thresholds = thr_screening_all_tests,
#                                   ##
#                                   baseline_plot_settings = baseline_plot_settings,
#                                   metareg_plot_settings = NULL)
# ##
# # outs_plot_pairwise_high_study_count_tests
# outs_plot_pairwise_all_tests$combined
# ##
# ## ---- Save plot(s):
# ##
# ggsave("AAD_REAL_Application_Se_Sp_pairwise_diffs_screening_all_tests.png",
#        outs_plot_pairwise_all_tests$combined,
#        width = 12, height = 12, dpi = 200)
##
## -------------------------------------------------------------------------------------------------------------------------------------
##
# Get your comparisons
# tibble_comp <- Baseline_Model_best$extract_NMA_comparisons(
#   test_names = c("GAD-2", "GAD-7", "HADS", "BAI", "OASIS", "STAI_S", "STAI_T", "PROMIS")
# )
thr_screening_all_tests = list(
  "GAD-2" = c(2), ## 7 ---- 3 is standard screening thr.
  "GAD-7" = c(4, 5, 6, 7), ## 21 ---- 
  "HADS"  = c(5, 6, 7, 8), ## 21 ----
  ##
  "BAI"    = c(2, 4, 6, 8, 10, 12),  # 20 ---- 8-10 most common screening cutoffs
  "OASIS"  = c(1, 2, 3, 4, 5, 6),   # 60 ---- 8 is standard screening, 10-12 also used
  "STAI_S" = c(52, 54, 56, 58, 60, 62, 64, 66, 68) - 20, # 60 ---- After -20 transform: screening range
  "STAI_T" = c(52, 54, 56, 58, 60, 62, 64, 66, 68) - 20,  # 60 ---- Match STAI-S
  "PROMIS" = c(12, 14, 16, 18, 20, 22)   # 32 ---- T-score ~60-65 equivalent for screening
)
##
tibble_comp <- outs_NMA_comp
##
## Add scenario column
##
tibble_comp$scenario <- "Baseline (no covariates)"
##
## Filter by thresholds
##
tibble_comp_filtered <- filter_clinically_relevant_comparisons(
  tibble_comp = tibble_comp,
  relevant_thresholds = thr_screening_all_tests
)
#
##
## Split by reference test (test1)
##
reference_tests <- c("GAD-2", "GAD-7", "HADS", "BAI", "OASIS", "STAI_S", "STAI_T")
##
plots_by_ref <- subset_comp <- list()
##
for (ref_test in reference_tests) {
  
      # Filter to just comparisons involving this reference test
      subset_comp[[ref_test]] <- tibble_comp_filtered %>% 
                                 filter(test1_name == ref_test)
      
      if (nrow(subset_comp[[ref_test]]) > 0) {
        plots_by_ref[[ref_test]] <- plot_pairwise_differences(  tibble_comp = subset_comp[[ref_test]],
                                                                comparison_type = "both",
                                                                test_names = c("GAD-2", "GAD-7", "HADS",
                                                                               "BAI", "OASIS", "STAI_S", "STAI_T", "PROMIS"),
                                                                relevant_thresholds = NULL,  # Already filtered!
                                                                baseline_plot_settings = baseline_plot_settings)
      }
}
##
## ---- Plots:
##
plots_by_ref$`GAD-2`$combined
plots_by_ref$`GAD-7`$combined
plots_by_ref$HADS$combined
##
plots_by_ref$BAI$combined
plots_by_ref$OASIS$combined
plots_by_ref$STAI_S$combined
plots_by_ref$STAI_T$combined
##
## ---- Tibble's:
##
subset_comp$`GAD-2` %>% print(n = 50)
subset_comp$`GAD-7` %>% print(n = 200)
subset_comp$HADS    %>% print(n = 150)
##
subset_comp$BAI    %>% print(n = 100)
subset_comp$OASIS  %>% print(n = 100)
subset_comp$STAI_S %>% print(n = 100)
subset_comp$STAI_T %>% print(n = 100)
##
## ---- "Roughly equivalent" SCREENING thresholds for all 8 tests (prioritises Se, as screening):
##
# Define your screening thresholds
# thr_screening_reference_tests <- list(
#   "GAD-2"  = c(3),
#   "GAD-7"  = c(8, 9, 10, 11),
#   "HADS"   = c(8, 9, 10, 11),
#   "BAI"    = c(14, 15, 16, 17, 18, 20),
#   "OASIS"  = c(6, 7, 8, 9, 10, 11),
#   "STAI_S" = c(64, 68, 72, 74, 75, 76) - 20,
#   "STAI_T" = c(64, 68, 72, 74, 75, 76) - 20
# )
##
thr_screening_reference_tests = list(
  "GAD-2" = c(2), ## 7 ---- 3 is standard screening thr.
  "GAD-7" = c(4, 5, 6, 7), ## 21 ---- 
  "HADS"  = c(5, 6, 7, 8), ## 21 ----
  ##
  "BAI"    = c(2, 4, 6, 8, 10, 12),  # 20 ---- 8-10 most common screening cutoffs
  "OASIS"  = c(1, 2, 3, 4, 5, 6),   # 60 ---- 8 is standard screening, 10-12 also used
  "STAI_S" = c(52, 54, 56, 58, 60, 62, 64, 66, 68) - 20, # 60 ---- After -20 transform: screening range
  "STAI_T" = c(52, 54, 56, 58, 60, 62, 64, 66, 68) - 20# ,  # 60 ---- Match STAI-S
  # "PROMIS" = c(12, 14, 16, 18, 20, 22)   # 32 ---- T-score ~60-65 equivalent for screening
)
##
# Subset to bracketing thresholds only
subset_comp_bracketed <- subset_to_bracketing_thresholds(
  subset_comp = subset_comp,
  thr_screening_reference_tests = thr_screening_reference_tests
)

# Check the reduction
sapply(subset_comp, nrow)            # Original sizes
sapply(subset_comp_bracketed, nrow)  # Bracketed sizes (should be much smaller!)

# Extract the relevant thresholds for plotting
relevant_thresholds_bracketed <- extract_bracketed_thresholds( subset_comp_bracketed = subset_comp_bracketed)

# Now plot each reference test
plots_bracketed <- list()

for (ref_test in names(subset_comp_bracketed)) {
  
      plots_bracketed[[ref_test]] <-  plot_pairwise_differences(
                                      tibble_comp = subset_comp_bracketed[[ref_test]],
                                      comparison_type = "both",
                                      test_names = c("GAD-2", "GAD-7", "HADS", "BAI", "OASIS", "STAI_S", "STAI_T", "PROMIS"),
                                      relevant_thresholds = NULL,  # Already subsetted!
                                      baseline_plot_settings = baseline_plot_settings)
}

##
## ---- Plots:
##
plots_bracketed$`GAD-2`$combined
plots_bracketed$`GAD-7`$combined
plots_bracketed$HADS$combined
##
plots_bracketed$BAI$combined
plots_bracketed$OASIS$combined
plots_bracketed$STAI_S$combined
plots_bracketed$STAI_T$combined
##
## ---- Tibble's:
##
subset_comp_bracketed$`GAD-2` %>% select(test1_name, threshold1, test2_name, threshold2,
                                         diff_Se_mean, diff_Se_lower, diff_Se_upper,
                                         diff_Sp_mean, diff_Sp_lower, diff_Sp_upper) %>% print(n = 50, width = 9000)
subset_comp_bracketed$`GAD-7` %>% select(test1_name, threshold1, test2_name, threshold2,
                                         diff_Se_mean, diff_Se_lower, diff_Se_upper,
                                         diff_Sp_mean, diff_Sp_lower, diff_Sp_upper) %>% print(n = 50, width = 9000)
subset_comp_bracketed$HADS %>% select(test1_name, threshold1, test2_name, threshold2,
                                         diff_Se_mean, diff_Se_lower, diff_Se_upper,
                                         diff_Sp_mean, diff_Sp_lower, diff_Sp_upper) %>% print(n = 50, width = 9000)
##
subset_comp_bracketed$BAI %>% select(test1_name, threshold1, test2_name, threshold2,
                                        diff_Se_mean, diff_Se_lower, diff_Se_upper,
                                        diff_Sp_mean, diff_Sp_lower, diff_Sp_upper) %>% print(n = 50, width = 9000)
subset_comp_bracketed$OASIS %>% select(test1_name, threshold1, test2_name, threshold2,
                                     diff_Se_mean, diff_Se_lower, diff_Se_upper,
                                     diff_Sp_mean, diff_Sp_lower, diff_Sp_upper) %>% print(n = 50, width = 9000)
subset_comp_bracketed$STAI_S %>% select(test1_name, threshold1, test2_name, threshold2,
                                     diff_Se_mean, diff_Se_lower, diff_Se_upper,
                                     diff_Sp_mean, diff_Sp_lower, diff_Sp_upper) %>% print(n = 50, width = 9000)
subset_comp_bracketed$STAI_T %>% select(test1_name, threshold1, test2_name, threshold2,
                                     diff_Se_mean, diff_Se_lower, diff_Se_upper,
                                     diff_Sp_mean, diff_Sp_lower, diff_Sp_upper) %>% print(n = 50, width = 9000)
##
## ---- Save:
##
scale_factor <- 0.80
##
ggsave( filename = paste0(target_condition, "_pairwise_bracketed_", gsub("-", "", "GAD-2"), ".png"),
        plot = plots_bracketed$`GAD-2`$combined,
        width = 16*scale_factor, height = 10*scale_factor,
        dpi = 200)
##
ggsave( filename = paste0(target_condition, "_pairwise_bracketed_", gsub("-", "", "GAD-7"), ".png"),
        plot = plots_bracketed$`GAD-7`$combined,
        width = 16*scale_factor, height = 16*scale_factor,
        dpi = 200)
##
ggsave( filename = paste0(target_condition, "_pairwise_bracketed_", gsub("-", "", "HADS"), ".png"),
        plot = plots_bracketed$HADS$combined,
        width = 16*scale_factor, height = 16*scale_factor,
        dpi = 200)
##
##
##
ggsave( filename = paste0(target_condition, "_pairwise_bracketed_", gsub("-", "", "BAI"), ".png"),
        plot = plots_bracketed$BAI$combined,
        width = 16*scale_factor, height = 20*scale_factor,
        dpi = 200)
##
##
##
ggsave( filename = paste0(target_condition, "_pairwise_bracketed_", gsub("-", "", "OASIS"), ".png"),
        plot = plots_bracketed$OASIS$combined,
        width = 16*scale_factor, height = 16*scale_factor,
        dpi = 200)
##
##
##
ggsave( filename = paste0(target_condition, "_pairwise_bracketed_", gsub("-", "", "STAI_S"), ".png"),
        plot = plots_bracketed$STAI_S$combined,
        width = 16*scale_factor, height = 12*scale_factor,
        dpi = 200)
##
ggsave( filename = paste0(target_condition, "_pairwise_bracketed_", gsub("-", "", "STAI_T"), ".png"),
        plot = plots_bracketed$STAI_T$combined,
        width = 16*scale_factor, height = 8*scale_factor,
        dpi = 200)
 



































equiv_screening_thr <- list()
##
for (ref_test in reference_tests) {
      equiv_screening_thr[[ref_test]] <- find_equivalent_thresholds_bracketed(  tibble_comp = subset_comp[[ref_test]], 
                                                                                require_CI_cross = TRUE,
                                                                                Sp_warning_threshold = 0.10)
}

equiv_screening_thr$`GAD-2`
##
subset_comp$`GAD-7`
equiv_screening_thr$`GAD-7`

tibble_comp_filtered %>% print(n = 200)
##


generate_equiv_thresholds_latex(  results_tibble = equiv_screening_thr[[ref_test]], 
                                  reference_test = "GAD-2",
                                  reference_threshold = 3,
                                  caption = NULL,
                                  label = "tab:equiv_thresholds")


choose(8, 2)
choose(5, 2)




# plot_data <- outs_plot_pairwise_high_study_count_tests$combined
# 
# max_per_plot = 10
# 
# comparison_labels <- unique(plot_data$comparison_label)
# n_chunks <- ceiling(length(comparison_labels) / max_per_plot)
# 
# for (i in 1:n_chunks) {
#   start_idx <- (i-1) * max_per_plot + 1
#   end_idx <- min(i * max_per_plot, length(comparison_labels))
#   
#   chunk_comparisons <- comparison_labels[start_idx:end_idx]
#   chunk_data <- plot_data %>% filter(comparison_label %in% chunk_comparisons)
#   
#   # plots[[paste0("chunk_", i)]]$Se <- ...
# }
# 
# 
# 
# 
# 















real_data$x[[4]][[1]]



# Check what thresholds actually exist for BAI
test_comp_debug <- Baseline_Model_best$extract_NMA_comparisons(test_names = test_names)

test_comp_debug <- outs_NMA_comp

cat("\n=== Available BAI thresholds ===\n")
test_comp_debug %>% 
  filter(test1_name == "BAI" | test2_name == "BAI") %>%
  pull(threshold1) %>% 
  unique() %>% 
  sort() %>%
  print()



# Debug: Check data BEFORE plotting
cat("\n=== Checking filtered data ===\n")
cat("Total rows before filtering:", nrow(test_comp_debug), "\n")

# Manual filter to see what's happening
filtered_debug <- test_comp_debug %>%
  filter(
    test1_name %in% c("BAI", "OASIS") &
      test2_name %in% c("BAI", "OASIS")
  )

cat("Rows after test filtering:", nrow(filtered_debug), "\n")

# Now check threshold filtering
filtered_debug2 <- filtered_debug %>%
  filter(
    (test1_name == "BAI" & threshold1 %in% c(8, 10, 16)) |
      (test1_name == "OASIS" & threshold1 %in% c(8, 9, 10, 11, 12))
  ) %>%
  filter(
    (test2_name == "BAI" & threshold2 %in% c(8, 10, 16)) |
      (test2_name == "OASIS" & threshold2 %in% c(8, 9, 10, 11, 12))
  )

cat("Rows after threshold filtering:", nrow(filtered_debug2), "\n")

if (nrow(filtered_debug2) == 0) {
  cat("\nPROBLEM: No data left after filtering!\n")
  cat("Checking first few rows of unfiltered:\n")
  filtered_debug %>% select(test1_name, test2_name, threshold1, threshold2) %>% head(20) %>% print()
}







cat("\n=== The 15 rows that survived filtering ===\n")
filtered_debug2 %>% 
  select(test1_name, test2_name, threshold1, threshold2, 
         diff_Se_mean, diff_Sp_mean) %>% 
  print(n = 15)


# Check for any NA or empty columns
cat("\n=== Checking for problematic columns ===\n")
sapply(filtered_debug2, function(x) sum(is.na(x))) %>% 
  .[. > 0] %>% 
  print()

# Check unique values in key columns
cat("\nUnique test1_name:", unique(filtered_debug2$test1_name), "\n")
cat("Unique test2_name:", unique(filtered_debug2$test2_name), "\n")






cat("\n=== All column names ===\n")
colnames(filtered_debug2) %>% print()

if ("scenario" %in% colnames(filtered_debug2)) {
  cat("\nUnique scenario values:", unique(filtered_debug2$scenario), "\n")
}

















