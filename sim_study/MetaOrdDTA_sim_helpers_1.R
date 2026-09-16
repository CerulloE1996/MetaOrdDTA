


# ============================================================================================================
# Process Single Scenario Function - MetaOrdDTA Version
# ============================================================================================================
process_single_scenario_MetaOrd <- function(output_dir,
                                            reduction_factor,
                                            scenario_row,
                                            target_MCSE_RMSE,
                                            max_seed_to_check,
                                            min_runs_needed,
                                            min_ESS_thr,
                                            vec_index_inner,
                                            min_studies_per_threshold
) {
        
        # Extract scenario parameters
        DGM <- scenario_row$DGM
        test <- scenario_row$index_test_MA
        n_studies <- scenario_row$n_studies
        model_parameterisation <- scenario_row$model_parameterisation
        ##
        true_vals <- get_true_values(test = test,
                                     DGM = DGM)
        true_DGM_Se <- true_vals$Se
        true_DGM_Sp <- true_vals$Sp
        
        # Create file name string
        file_name_string <- create_file_name_string( scenario_row = scenario_row,
                                                     reduction_factor = reduction_factor)
        
        # Load results
        complete_results_list <- list()
        for (seed in 1:max_seed_to_check) {
          file_path <- file.path(output_dir, 
                                 paste0("seed_", seed, "_", 
                                        file_name_string, 
                                        "_results_list.RDS"))
          if (file.exists(file_path)) {
            try({
              complete_results_list[[seed]] <- readRDS(file_path)
            }, silent = TRUE)
          }
        }
        
        # Remove NULL entries
        complete_results_list <- complete_results_list[!sapply(complete_results_list, is.null)]
        
        if (length(complete_results_list) < min_runs_needed) {
              
              warning(sprintf("Insufficient results for scenario %d (found %d, need %d)", 
                              scenario_row$scenario_id, length(complete_results_list), min_runs_needed))
              return(NULL)
          
        }
        
        # Compute metrics
        metrics <-  compute_simulation_metrics(
                    test = test,
                    complete_results_list = complete_results_list,
                    output_dir = output_dir,
                    true_DGM_Se = true_DGM_Se,
                    true_DGM_Sp = true_DGM_Sp,
                    min_ESS_thr = min_ESS_thr,
                    vec_index_inner = vec_index_inner,
                    min_studies_per_threshold = min_studies_per_threshold)
        
        mean_min_ESS <- metrics$mean_min_ESS
     
        # Create summary table
        metrics_tables <- create_metrics_summary_table(
                          metrics = metrics,
                          vec_index_inner = vec_index_inner)
        ##
        max_MCSE_Se <- max(as.numeric(metrics_tables$by_threshold$MCSE_RMSE_Se), na.rm = TRUE)
        max_MCSE_Sp <- max(as.numeric(metrics_tables$by_threshold$MCSE_RMSE_Sp), na.rm = TRUE)
        max_MCSE <- max(max_MCSE_Se, max_MCSE_Sp)
        ##
        avg_MCSE_Se <- mean(as.numeric(metrics_tables$by_threshold$MCSE_RMSE_Se), na.rm = TRUE)
        avg_MCSE_Sp <- mean(as.numeric(metrics_tables$by_threshold$MCSE_RMSE_Sp), na.rm = TRUE)
        avg_MCSE <- max(avg_MCSE_Se, avg_MCSE_Sp)
        ##
        MCSE_val <- round((max_MCSE + avg_MCSE)/2.0)
        ##
        mcse_met <- (!is.na(MCSE_val) && MCSE_val <= target_MCSE_RMSE)
        ##
        # Calculate minimum simulations needed
        min_sims_needed <-  calc_min_sim_size_RMSE( metrics = metrics,
                                                    target_MCSE_RMSE = target_MCSE_RMSE,
                                                    vec_index_inner = vec_index_inner)
        
        return(list(
          mean_min_ESS = mean_min_ESS,
          metrics = metrics,
          metrics_tables = metrics_tables,
          n_sims = metrics$n_sims,
          n_excluded = metrics$overall$n_excluded,
          MCSE_val = MCSE_val,
          mcse_met = mcse_met,
          min_sims_needed = min_sims_needed,
          test = test,
          n_studies = n_studies,
          model_parameterisation = model_parameterisation
        ))
  
}




# ============================================================================================================
# Create Summary Row Function - MetaOrdDTA Version
# ============================================================================================================
create_summary_row_MetaOrd <- function(result,
                                       scenario_row
) {
  
        if (is.null(result)) {
          return(tibble())
        }
        
        metrics <- result$metrics
        
        tibble(
          # Scenario info
          scenario_id = scenario_row$scenario_id,
          test = scenario_row$index_test_MA,
          DGM = scenario_row$DGM,
          n_studies = scenario_row$n_studies,
          model_parameterisation = scenario_row$model_parameterisation,
          alpha_pi = scenario_row$alpha_pi,
          softplus = scenario_row$softplus,
          random_thresholds = scenario_row$random_thresholds,
          Dirichlet_random_effects_type = scenario_row$Dirichlet_random_effects_type,
          N_per_study_mean = scenario_row$N_per_study_mean,
          ##
          n_sims = result$n_sims,
          n_excluded = result$n_excluded,
          ##
          min_sims_needed = result$min_sims_needed$n_sims_required,
          ##
          MCSE_val = result$MCSE_val,
          ##
          mcse_met = result$mcse_met,
          ##
          ## ---- RMSE:
          ##
          avg_RMSE_Se = metrics$overall$avg_RMSE_Se,
          avg_RMSE_Sp = metrics$overall$avg_RMSE_Sp,
          avg_RMSE = metrics$overall$avg_RMSE_Se + metrics$overall$avg_RMSE_Sp,
          ##
          max_RMSE_Se = metrics$overall$max_RMSE_Se,
          max_RMSE_Sp = metrics$overall$max_RMSE_Sp,
          ##
          ## ---- MCSE for RMSE:
          ##
          avg_MCSE_RMSE_Se = metrics$overall$avg_MCSE_RMSE_Se,
          avg_MCSE_RMSE_Sp = metrics$overall$avg_MCSE_RMSE_Sp,
          ##
          max_MCSE_RMSE_Se = metrics$overall$max_MCSE_RMSE_Se,
          max_MCSE_RMSE_Sp = metrics$overall$max_MCSE_RMSE_Sp,
          ##
          ## ---- Bias:
          ##
          avg_abs_bias_Se = metrics$overall$avg_abs_bias_Se,
          avg_abs_bias_Sp = metrics$overall$avg_abs_bias_Sp,
          avg_abs_bias = metrics$overall$avg_abs_bias_Se + metrics$overall$avg_abs_bias_Sp,
          ##
          max_abs_bias_Se = metrics$overall$max_abs_bias_Se,
          max_abs_bias_Sp = metrics$overall$max_abs_bias_Sp,
          ##
          ## ---- MCSE for Bias:
          ##
          avg_MCSE_bias_Se = metrics$overall$avg_MCSE_bias_Se,
          avg_MCSE_bias_Sp = metrics$overall$avg_MCSE_bias_Sp,
          ##
          ## ---- Coverage:
          ##
          avg_coverage_Se = metrics$overall$avg_coverage_prob_Se,
          avg_coverage_Sp = metrics$overall$avg_coverage_prob_Sp,
          ##
          ## ---- CI Width:
          ##
          avg_width_Se = metrics$overall$avg_interval_width_Se,
          avg_width_Sp = metrics$overall$avg_interval_width_Sp,
          ##
          mean_min_ESS = metrics$mean_min_ESS
        )
  
}







# ============================================================================================================
# Process All Scenarios - MetaOrdDTA Version
# ============================================================================================================
process_all_scenarios_MetaOrd <- function(  scenario_df,
                                            output_dir,
                                            test_filter = NULL,
                                            reduction_factor = NULL,
                                            n_studies_filter = NULL,
                                            model_filter = NULL,
                                            target_MCSE_RMSE,
                                            max_seed_to_check,
                                            min_runs_needed = 30,
                                            min_ESS_thr = 1000,
                                            vec_index_inner,
                                            batch_size = 100,
                                            max_iterations = 50,
                                            min_studies_per_threshold
                                            
) {
        
        require(dplyr)
  
        if (is.null(reduction_factor)) {
          reduction_factor <- 1.0
        }
        
        # Filter scenarios if needed
        scenarios_to_process <- scenario_df
        
        if (!is.null(test_filter)) {
          scenarios_to_process <- scenarios_to_process %>% 
            filter(index_test_MA %in% test_filter)
        }
        
        if (!is.null(n_studies_filter)) {
          scenarios_to_process <- scenarios_to_process %>% 
            filter(n_studies %in% n_studies_filter)
        }
        
        if (!is.null(model_filter)) {
          scenarios_to_process <- scenarios_to_process %>% 
            filter(model_parameterisation %in% model_filter)
        }
        
        # Store all results
        all_results <- list()
        summary_table <- tibble()
        
        i <- 1
        
        for (i in 1:nrow(scenarios_to_process)) {

                scenario_row <- scenarios_to_process[i, ]

                # Process this scenario
                result <- process_single_scenario_MetaOrd(
                          reduction_factor = reduction_factor,
                          output_dir = output_dir,
                          scenario_row = scenario_row,
                          target_MCSE_RMSE = target_MCSE_RMSE,
                          max_seed_to_check = max_seed_to_check,
                          min_runs_needed = min_runs_needed,
                          min_ESS_thr = min_ESS_thr,
                          vec_index_inner = vec_index_inner,
                          min_studies_per_threshold = min_studies_per_threshold)

                mean_min_ESS <- result$mean_min_ESS

                # Store results
                if (!is.null(result)) {
                  # if (mean_min_ESS > min_ESS_thr) {
                      result$scenario_info <- scenario_row
                      all_results[[paste0("scenario_", scenario_row$scenario_id)]] <- result

                      # Add to summary table
                      summary_row <- create_summary_row_MetaOrd(result,
                                                                scenario_row)
                      summary_table <- bind_rows(summary_table,
                                                 summary_row)
                  # }
                }

        }
        
        summary_table <- summary_table %>% arrange(test, DGM)
        
        return(list(
          all_results = all_results,
          summary_table = summary_table
        ))
  
}








# ---- MetaOrdDTA MODEL SUMMARY ----
create_MetaOrd_summary <- function( data, 
                                    metric = "RMSE"
) {
  
        group_by_alpha <- FALSE
        
        require(dplyr)
        require(knitr)
        
        # Check if all MCSEs are at or below 0.125
        all_mcse_se <- data$avg_MCSE_RMSE_Se
        all_mcse_sp <- data$avg_MCSE_RMSE_Sp
        
        if(any(all_mcse_se > 0.125, na.rm = TRUE) || any(all_mcse_sp > 0.125, na.rm = TRUE)) {
          
              cat("\n================================================================================\n")
              cat("ERROR: Not all MCSE values are at or below 0.125\n")
              cat("Cannot generate summary - some models have insufficient precision\n")
              cat("Models with MCSE > 0.125:\n")
              
              problem_models <- data %>%
                filter(avg_MCSE_RMSE_Se > 0.125 | avg_MCSE_RMSE_Sp > 0.125) %>%
                select(test, DGM, n_studies, model_parameterisation, 
                       avg_MCSE_RMSE_Se, avg_MCSE_RMSE_Sp)
              
              print(kable(problem_models, digits = 3, format = "simple"))
              cat("================================================================================\n")
              return(invisible(NULL))
          
        }
        
        cat("\n================================================================================\n")
        cat("MetaOrdDTA MODEL SUMMARY:\n")
        cat("Note: Using fixed MCSE = 0.125 for all CI calculations\n")
        cat("================================================================================\n\n")
        
        # Helper function to find statistically equivalent models based on CI overlap
        find_equivalent_models <- function(df,
                                           metric_col,
                                           mcse_col = NULL,
                                           find_worst = FALSE,
                                           z_score = 1.96) {  # 95% CI by default
          
                  if(find_worst) {
                    df <- df %>% arrange(desc(!!sym(metric_col)))  # Sort descending for worst
                  } else {
                    df <- df %>% arrange(!!sym(metric_col))  # Sort ascending for best
                  }
                  
                  if(nrow(df) == 0) return(df)
                  
                  if(!is.null(mcse_col)) {
                    # Calculate CIs for all models
                    df <- df %>%
                      mutate(
                        lower_CI = !!sym(metric_col) - z_score * !!sym(mcse_col),
                        upper_CI = !!sym(metric_col) + z_score * !!sym(mcse_col)
                      )
                    
                    best_lower <- df$lower_CI[1]
                    best_upper <- df$upper_CI[1]
                    
                    # Check which models have overlapping CIs with the best
                    equivalent <- df %>%
                      filter(lower_CI < best_upper & upper_CI > best_lower) %>%
                      select(-lower_CI, -upper_CI)
                  } else {
                    # Fallback if no MCSE column provided
                    equivalent <- df %>% filter(row_number() == 1)
                  }
                  
                  equivalent <- equivalent %>%
                    mutate(
                      rank = row_number(),
                      significant = ifelse(rank == 1 & nrow(equivalent) == 1, "*", ""),
                      equiv_group = paste0("Tied (n=", nrow(equivalent), ")")
                    )
                  
                  return(equivalent)
          
        }
        
        # Best performing scenarios per test/n_studies
        cat("\n1. BEST PERFORMING MODELS:\n")
        cat("----------------------------------------\n")
        cat("Note: * indicates statistically significant difference (non-overlapping 95% CIs)\n")
        cat("      Multiple models shown if CIs overlap with best\n\n")
        
        # cat("1a. BEST FOR SENSITIVITY (RMSE_Se in %):\n")
        best_Se <- data %>%
          mutate(RMSE_Se_pct = avg_RMSE_Se * 1,
                 MCSE_Se_pct = avg_MCSE_RMSE_Se * 1,
                 Bias_Se_pct = abs(avg_abs_bias_Se) * 1,
                 MCSE_Bias_Se_pct = avg_MCSE_bias_Se * 1)
        
        # Group by test, n_studies, and optionally alpha_pi
        if (group_by_alpha) {
          best_Se <- best_Se %>% group_by(test, DGM, n_studies, alpha_pi)
        } else {
          best_Se <- best_Se %>% group_by(test, DGM, n_studies)
        }
        
        if (metric == "RMSE") {
          
                best_Se <- best_Se %>% 
                  group_modify(~ find_equivalent_models(.x, "RMSE_Se_pct", "MCSE_Se_pct")) %>%
                  select(test, DGM, n_studies, model_parameterisation, alpha_pi, random_thresholds,
                         RMSE_Se_pct, MCSE_Se_pct, Bias_Se_pct, avg_coverage_Se, 
                         significant, equiv_group) %>%
                  rename(`Test` = test,
                         `N Studies` = n_studies,
                         `Model` = model_parameterisation,
                         `Alpha_pi` = alpha_pi,
                         `Rand_Thr` = random_thresholds,
                         `RMSE_Se(%)` = RMSE_Se_pct,
                         `MCSE(%)` = MCSE_Se_pct,
                         `|Bias_Se|(%)` = Bias_Se_pct,
                         `Coverage_Se(%)` = avg_coverage_Se,
                         Sig = significant,
                         Group = equiv_group)
          
        } else { 
                
                best_Se <- best_Se %>% 
                  group_modify(~ find_equivalent_models(.x, "Bias_Se_pct", "MCSE_Bias_Se_pct")) %>%
                  select(test, DGM, n_studies, model_parameterisation, alpha_pi, random_thresholds,
                         RMSE_Se_pct, MCSE_Se_pct, Bias_Se_pct, avg_coverage_Se, 
                         significant, equiv_group) %>%
                  rename(`Test` = test,
                         `N Studies` = n_studies,
                         `Model` = model_parameterisation,
                         `Alpha_pi` = alpha_pi,
                         `Rand_Thr` = random_thresholds,
                         `RMSE_Se(%)` = RMSE_Se_pct,
                         `MCSE(%)` = MCSE_Se_pct,
                         `|Bias_Se|(%)` = Bias_Se_pct,
                         `Coverage_Se(%)` = avg_coverage_Se,
                         Sig = significant,
                         Group = equiv_group)
          
        }
        
        # print(kable(best_Se, digits = 2, format = "simple"))
        
        # cat("\n1b. BEST FOR SPECIFICITY (RMSE_Sp in %):\n")
        best_Sp <- data %>%
          mutate(RMSE_Sp_pct = avg_RMSE_Sp * 1,
                 MCSE_Sp_pct = avg_MCSE_RMSE_Sp * 1,
                 Bias_Sp_pct = abs(avg_abs_bias_Sp) * 1,
                 MCSE_Bias_Sp_pct = avg_MCSE_bias_Sp * 1)
        
        if (group_by_alpha) {
          best_Sp <- best_Sp %>% group_by(test, DGM, n_studies, alpha_pi)
        } else {
          best_Sp <- best_Sp %>% group_by(test, DGM, n_studies)
        }
        
        if (metric == "RMSE") {
          
          best_Sp <- best_Sp %>% 
            group_modify(~ find_equivalent_models(.x, "RMSE_Sp_pct", "MCSE_Sp_pct")) %>%
            select(test, DGM, n_studies, model_parameterisation, alpha_pi, random_thresholds,
                   RMSE_Sp_pct, MCSE_Sp_pct, Bias_Sp_pct, avg_coverage_Sp, 
                   significant, equiv_group) %>%
            rename(`Test` = test,
                   `N Studies` = n_studies,
                   `Model` = model_parameterisation,
                   `Alpha_pi` = alpha_pi,
                   `Rand_Thr` = random_thresholds,
                   `RMSE_Sp(%)` = RMSE_Sp_pct,
                   `MCSE(%)` = MCSE_Sp_pct,
                   `|Bias_Sp|(%)` = Bias_Sp_pct,
                   `Coverage_Sp(%)` = avg_coverage_Sp,
                   Sig = significant,
                   Group = equiv_group)
          
        } else { 
          
          best_Sp <- best_Sp %>% 
            group_modify(~ find_equivalent_models(.x, "Bias_Sp_pct", "MCSE_Bias_Sp_pct")) %>%
            select(test, DGM, n_studies, model_parameterisation, alpha_pi, random_thresholds,
                   RMSE_Sp_pct, MCSE_Sp_pct, Bias_Sp_pct, avg_coverage_Sp, 
                   significant, equiv_group) %>%
            rename(`Test` = test,
                   `N Studies` = n_studies,
                   `Model` = model_parameterisation,
                   `Alpha_pi` = alpha_pi,
                   `Rand_Thr` = random_thresholds,
                   `RMSE_Sp(%)` = RMSE_Sp_pct,
                   `MCSE(%)` = MCSE_Sp_pct,
                   `|Bias_Sp|(%)` = Bias_Sp_pct,
                   `Coverage_Sp(%)` = avg_coverage_Sp,
                   Sig = significant,
                   Group = equiv_group)
          
        }
        
        # print(kable(best_Sp, digits = 2, format = "simple"))
        
        cat("\n1c. BEST OVERALL (RMSE_Se + RMSE_Sp in %):\n")
        best_overall <- data %>%
          mutate(RMSE_Se_pct = avg_RMSE_Se,
                 RMSE_Sp_pct = avg_RMSE_Sp,
                 RMSE_total_pct = RMSE_Se_pct + RMSE_Sp_pct,
                 MCSE_Se_pct = avg_MCSE_RMSE_Se,
                 MCSE_Sp_pct = avg_MCSE_RMSE_Sp,
                 MCSE_total_pct = sqrt(MCSE_Se_pct^2 + MCSE_Sp_pct^2),
                 Bias_Se_pct = abs(avg_abs_bias_Se),
                 Bias_Sp_pct = abs(avg_abs_bias_Sp),
                 Bias_total_pct = Bias_Se_pct + Bias_Sp_pct)
        
        if (group_by_alpha) {
          best_overall <- best_overall %>% group_by(test, DGM, n_studies, alpha_pi)
        } else {
          best_overall <- best_overall %>% group_by(test, DGM, n_studies)
        }
        
        best_overall <- best_overall %>%
          group_modify(~ find_equivalent_models(.x, "RMSE_total_pct", "MCSE_total_pct")) %>%
          select(test, DGM, n_studies, model_parameterisation, alpha_pi, random_thresholds,
                 RMSE_total_pct, MCSE_total_pct,
                 Bias_total_pct, Bias_Se_pct, Bias_Sp_pct,
                 avg_coverage_Se, avg_coverage_Sp, 
                 significant, equiv_group) %>%
          rename(`Test` = test,
                 `DGM` = DGM,
                 `N Studies` = n_studies,
                 `Model` = model_parameterisation,
                 `Alpha_pi` = alpha_pi,
                 `Rand_Thr` = random_thresholds,
                 `RMSE_Total(%)` = RMSE_total_pct,
                 `MCSE_Total(%)` = MCSE_total_pct,
                 `|Bias_Se|(%)` = Bias_Se_pct,
                 `|Bias_Sp|(%)` = Bias_Sp_pct,
                 `|Bias_Total|(%)` = Bias_total_pct,
                 `Cvg_Se(%)` = avg_coverage_Se,
                 `Cvg_Sp(%)` = avg_coverage_Sp,
                 Sig = significant,
                 Group = equiv_group)
        
        print(kable(best_overall, digits = 2, format = "simple"))
        
        
        # SECTION 2: WORST performing scenarios
        cat("\n\n2. WORST PERFORMING MODELS:\n")
        cat("----------------------------------------\n")
        cat("Note: Shows models NOT in the best group for each scenario\n\n")
        ##
        # # 2a. WORST FOR SENSITIVITY
        # cat("2a. WORST FOR SENSITIVITY (RMSE_Se in %):\n")
        # 
        # # First get the best performers
        # best_Se_names <- data %>%
        #   mutate(RMSE_Se_pct = avg_RMSE_Se,
        #          MCSE_Se_pct = avg_MCSE_RMSE_Se,
        #          Bias_Se_pct = abs(avg_abs_bias_Se),
        #          MCSE_Bias_Se_pct = avg_MCSE_bias_Se)
        # 
        # if (group_by_alpha) {
        #   best_Se_names <- best_Se_names %>% group_by(test, DGM, n_studies, alpha_pi)
        # } else {
        #   best_Se_names <- best_Se_names %>% group_by(test, DGM, n_studies)
        # }
        # 
        # if (metric == "RMSE") {
        #   best_Se_names <- best_Se_names %>%
        #     group_modify(~ find_equivalent_models(.x, "RMSE_Se_pct", "MCSE_Se_pct")) %>%
        #     select(test, DGM, n_studies, model_parameterisation, alpha_pi)
        # } else {
        #   best_Se_names <- best_Se_names %>%
        #     group_modify(~ find_equivalent_models(.x, "Bias_Se_pct", "MCSE_Bias_Se_pct")) %>%
        #     select(test, DGM, n_studies, model_parameterisation, alpha_pi)
        # }
        # 
        # # Now get worst as complement
        # worst_Se <- data %>%
        #   mutate(RMSE_Se_pct = avg_RMSE_Se,
        #          MCSE_Se_pct = avg_MCSE_RMSE_Se,
        #          Bias_Se_pct = abs(avg_abs_bias_Se)) %>%
        #   anti_join(best_Se_names, by = c("test", "DGM", "n_studies", "model_parameterisation", "alpha_pi"))
        # 
        # if (group_by_alpha) {
        #   worst_Se <- worst_Se %>% group_by(test, DGM, n_studies, alpha_pi)
        # } else {
        #   worst_Se <- worst_Se %>% group_by(test, DGM, n_studies)
        # }
        # 
        # worst_Se <- worst_Se %>%
        #   arrange(test, n_studies, desc(RMSE_Se_pct)) %>%
        #   select(test, DGM, n_studies, model_parameterisation, alpha_pi, random_thresholds,
        #          RMSE_Se_pct, MCSE_Se_pct, Bias_Se_pct, avg_coverage_Se) %>%
        #   rename(`Test` = test,
        #          `N Studies` = n_studies,
        #          `Model` = model_parameterisation,
        #          `Alpha_pi` = alpha_pi,
        #          `Rand_Thr` = random_thresholds,
        #          `RMSE_Se(%)` = RMSE_Se_pct,
        #          `MCSE(%)` = MCSE_Se_pct,
        #          `|Bias_Se|(%)` = Bias_Se_pct,
        #          `Coverage_Se(%)` = avg_coverage_Se)
        ##
        # if(nrow(worst_Se) > 0) {
        #   print(kable(worst_Se, digits = 2, format = "simple"))
        # } else {
        #   cat("(All models tied as best for all scenarios)\n")
        # # }
        # 
        # # 2b. WORST FOR SPECIFICITY
        # cat("\n2b. WORST FOR SPECIFICITY (RMSE_Sp in %):\n")
        # 
        # best_Sp_names <- data %>%
        #   mutate(RMSE_Sp_pct = avg_RMSE_Sp * 1,
        #          MCSE_Sp_pct = avg_MCSE_RMSE_Sp * 1,
        #          Bias_Sp_pct = abs(avg_abs_bias_Sp) * 1,
        #          MCSE_Bias_Sp_pct = avg_MCSE_bias_Sp * 1)
        # 
        # if (group_by_alpha) {
        #   best_Sp_names <- best_Sp_names %>% group_by(test, DGM, n_studies, alpha_pi)
        # } else {
        #   best_Sp_names <- best_Sp_names %>% group_by(test, DGM, n_studies)
        # }
        # 
        # if (metric == "RMSE") {
        #   best_Sp_names <- best_Sp_names %>%
        #     group_modify(~ find_equivalent_models(.x, "RMSE_Sp_pct", "MCSE_Sp_pct")) %>%
        #     select(test, DGM, n_studies, model_parameterisation, alpha_pi)
        # } else {
        #   best_Sp_names <- best_Sp_names %>%
        #     group_modify(~ find_equivalent_models(.x, "Bias_Sp_pct", "MCSE_Bias_Sp_pct")) %>%
        #     select(test, DGM, n_studies, model_parameterisation, alpha_pi)
        # }
        # 
        # worst_Sp <- data %>%
        #   mutate(RMSE_Sp_pct = avg_RMSE_Sp * 1,
        #          MCSE_Sp_pct = avg_MCSE_RMSE_Sp * 1,
        #          Bias_Sp_pct = abs(avg_abs_bias_Sp) * 1) %>%
        #   anti_join(best_Sp_names, by = c("test",  "DGM", "n_studies", "model_parameterisation", "alpha_pi"))
        # 
        # if (group_by_alpha) {
        #   worst_Sp <- worst_Sp %>% group_by(test, DGM, n_studies, alpha_pi)
        # } else {
        #   worst_Sp <- worst_Sp %>% group_by(test, DGM, n_studies)
        # }
        # 
        # worst_Sp <- worst_Sp %>%
        #   arrange(test, n_studies, desc(RMSE_Sp_pct)) %>%
        #   select(test, DGM, n_studies, model_parameterisation, alpha_pi, random_thresholds,
        #          RMSE_Sp_pct, MCSE_Sp_pct, Bias_Sp_pct, avg_coverage_Sp) %>%
        #   rename(`Test` = test,
        #          `N Studies` = n_studies,
        #          `Model` = model_parameterisation,
        #          `Alpha_pi` = alpha_pi,
        #          `Rand_Thr` = random_thresholds,
        #          `RMSE_Sp(%)` = RMSE_Sp_pct,
        #          `MCSE(%)` = MCSE_Sp_pct,
        #          `|Bias_Sp|(%)` = Bias_Sp_pct,
        #          `Coverage_Sp(%)` = avg_coverage_Sp)
        
        # if(nrow(worst_Sp) > 0) {
        #   print(kable(worst_Sp, digits = 2, format = "simple"))
        # } else {
        #   cat("(All models tied as best for all scenarios)\n")
        # }
        
        # 2c. WORST OVERALL
        cat("\n2c. WORST OVERALL (RMSE_Se + RMSE_Sp in %):\n")
        
        # Get the list of best model identifiers BEFORE renaming
        best_overall_names <- data %>%
          mutate(RMSE_Se_pct = avg_RMSE_Se,
                 RMSE_Sp_pct = avg_RMSE_Sp,
                 RMSE_total_pct = RMSE_Se_pct + RMSE_Sp_pct,
                 MCSE_Se_pct = avg_MCSE_RMSE_Se,
                 MCSE_Sp_pct = avg_MCSE_RMSE_Sp,
                 MCSE_total_pct = sqrt(MCSE_Se_pct^2 + MCSE_Sp_pct^2)) %>%
          group_by(test, DGM, n_studies) %>%
          group_modify(~ find_equivalent_models(.x, "RMSE_total_pct", "MCSE_total_pct")) %>%
          ungroup() %>%
          select(test, DGM, n_studies, model_parameterisation, alpha_pi, random_thresholds)
        
        ## ---- Now get ALL models NOT in that list
        worst_overall <- data %>%
          anti_join(best_overall_names, by = c("test", "DGM", "n_studies", 
                                               "model_parameterisation", 
                                               "alpha_pi", 
                                               "random_thresholds")) %>%
          mutate(RMSE_Se_pct = avg_RMSE_Se,
                 RMSE_Sp_pct = avg_RMSE_Sp,
                 RMSE_total_pct = RMSE_Se_pct + RMSE_Sp_pct,
                 Bias_Se_pct = abs(avg_abs_bias_Se),
                 Bias_Sp_pct = abs(avg_abs_bias_Sp),
                 Bias_total_pct = Bias_Se_pct + Bias_Sp_pct,
                 MCSE_total_pct = sqrt(avg_MCSE_RMSE_Se^2 + avg_MCSE_RMSE_Sp^2)) %>%
          arrange(desc(RMSE_total_pct)) %>%
          select(test, DGM, n_studies, model_parameterisation, alpha_pi, random_thresholds,
                 RMSE_total_pct, MCSE_total_pct, Bias_total_pct, Bias_Se_pct, Bias_Sp_pct,
                 avg_coverage_Se, avg_coverage_Sp) %>%
          rename(`Test` = test, 
                 `DGM` = DGM,
                 `N Studies` = n_studies, 
                 `Model` = model_parameterisation,
                 `Alpha_pi` = alpha_pi, 
                 `Rand_Thr` = random_thresholds,
                 `RMSE_Total(%)` = RMSE_total_pct, 
                 `MCSE_Total(%)` = MCSE_total_pct,
                 `|Bias_Total|(%)` = Bias_total_pct, 
                 `|Bias_Se|(%)` = Bias_Se_pct,
                 `|Bias_Sp|(%)` = Bias_Sp_pct, 
                 `Cvg_Se(%)` = avg_coverage_Se,
                 `Cvg_Sp(%)` = avg_coverage_Sp)
        
        if(nrow(worst_overall) > 0) {
          print(kable(worst_overall, digits = 2, format = "simple"))
        } else {
          cat("(All models tied as best)\n")
        }
        
}





 






