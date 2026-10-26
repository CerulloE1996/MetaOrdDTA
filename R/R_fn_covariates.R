





#' subset_X_covariates
#' @export
subset_X_covariates <- function( X, 
                                 covariates_to_keep) {
  
      # Get dimensions from X structure
      n_classes <- length(X)
      n_tests <- length(X[[1]])
        
      # Create output with same structure
      X_subset <- list()
        
      # Loop through disease classes
      for (c in 1:n_classes) {
          
          X_subset[[c]] <- list()
            
          # Loop through tests
          for (t in 1:n_tests) {
            
                # Get current matrix
                X_mat <- X[[c]][[t]]
                
                # Subset columns by name in the specified order
                # This handles the case where some matrices might be all zeros
                if (sum(X_mat[,1]) > 0) {  # Check if any non-missing data
                  X_subset[[c]][[t]] <- X_mat[, covariates_to_keep, drop = FALSE]
                } else {
                  # For all-zero matrices, create zero matrix with right columns
                  X_subset[[c]][[t]] <- matrix(0, 
                                               nrow = nrow(X_mat), 
                                               ncol = length(covariates_to_keep))
                  colnames(X_subset[[c]][[t]]) <- covariates_to_keep
                }
              
          }
        
      }

      return(X_subset)
      
}





#' map_beta_coefficients
#' @export 
map_beta_coefficients <- function( beta_tibble, 
                                   covariate_names = NULL,
                                   test_names = NULL,
                                   min_magnitude = 0.125,
                                   borderline_threshold = 0.10,
                                   show_borderline = TRUE
) {

      alpha <- 0.05
      
      library(tidyverse)

      # Extract indices from parameter names
      beta_tibble <- beta_tibble %>%
        mutate(
          # Extract indices using regex
          test_idx = as.numeric(str_extract(variable, "(?<=\\[)\\d+(?=,)")),
          class_idx = as.numeric(str_extract(variable, "(?<=,)\\d+(?=,)")),
          cov_idx = as.numeric(str_extract(variable, "(?<=,)\\d+(?=\\])")),

          # Determine significance at 95% level (using available columns)
          is_significant = !(`2.5%` <= 0 & `97.5%` >= 0),

          # Borderline: CI comes within threshold of excluding 0
          distance_to_zero = pmin(abs(`2.5%`), abs(`97.5%`)),
          is_borderline = !is_significant & distance_to_zero <= borderline_threshold,

          # Check magnitude
          has_min_magnitude = abs(median) >= min_magnitude,

          # Combined criteria
          meets_criteria = is_significant | (is_borderline & has_min_magnitude),

          # Add CI width for precision assessment
          ci_width = `97.5%` - `2.5%`
        )
      
        
        # # Extract indices from parameter names
        # beta_tibble <- beta_tibble %>%
        #   mutate(
        #     # Extract indices using regex
        #     test_idx = as.numeric(str_extract(variable, "(?<=\\[)\\d+(?=,)")),
        #     class_idx = as.numeric(str_extract(variable, "(?<=,)\\d+(?=,)")),
        #     cov_idx = as.numeric(str_extract(variable, "(?<=,)\\d+(?=\\])")),
        #     
        #     # Determine significance at 95% level (using available columns)
        #     is_significant = !(q5 <= 0 & q95 >= 0),
        #     
        #     # Borderline: CI comes within threshold of excluding 0
        #     distance_to_zero = pmin(abs(q5), abs(q95)),
        #     is_borderline = !is_significant & distance_to_zero <= borderline_threshold,
        #     
        #     # Check magnitude
        #     has_min_magnitude = abs(median) >= min_magnitude,
        #     
        #     # Combined criteria
        #     meets_criteria = is_significant | (is_borderline & has_min_magnitude),
        #     
        #     # Add CI width for precision assessment
        #     ci_width = q95 - q5
        #   )
        
        
        # Determine dimensions
        n_tests <- max(beta_tibble$test_idx)
        n_classes <- max(beta_tibble$class_idx)
        n_covariates <- max(beta_tibble$cov_idx)
        
        # Set default names if not provided
        if (is.null(test_names)) {
          test_names <- paste("Test", 1:n_tests)
        }
        if (is.null(covariate_names)) {
          covariate_names <- c("Intercept", paste0("Covariate_", 2:n_covariates))
        }
        
        class_names <- c("Non-diseased", "Diseased")
        
        # Add readable names to tibble
        beta_tibble <- beta_tibble %>%
          mutate(
            test_name = test_names[test_idx],
            class_name = class_names[class_idx],
            covariate_name = covariate_names[cov_idx]
          )
        
        # Create structured output
        results <- list()
        results$full_table <- beta_tibble
        results$dimensions <- list(n_tests = n_tests, n_classes = n_classes, n_covariates = n_covariates)
        results$alpha <- alpha
        results$min_magnitude <- min_magnitude
        
        # Create arrays
        coef_array <- array(NA, dim = c(n_tests, n_classes, n_covariates),
                            dimnames = list(test_names, class_names, covariate_names))
        
        ci_lower_array <- array(NA, dim = c(n_tests, n_classes, n_covariates),
                                dimnames = list(test_names, class_names, covariate_names))
        
        ci_upper_array <- array(NA, dim = c(n_tests, n_classes, n_covariates),
                                dimnames = list(test_names, class_names, covariate_names))
        
        signif_array <- array(FALSE, dim = c(n_tests, n_classes, n_covariates),
                              dimnames = list(test_names, class_names, covariate_names))
        
        meets_criteria_array <- array(FALSE, dim = c(n_tests, n_classes, n_covariates),
                                      dimnames = list(test_names, class_names, covariate_names))
        
        # Fill arrays
        for (i in 1:nrow(beta_tibble)) {
              
              t_idx <- beta_tibble$test_idx[i]
              c_idx <- beta_tibble$class_idx[i]
              j_idx <- beta_tibble$cov_idx[i]
              
              # coef_array[t_idx, c_idx, j_idx] <- beta_tibble$`50%`[i]
              ci_lower_array[t_idx, c_idx, j_idx] <- beta_tibble$`2.5%`[i]
              ci_upper_array[t_idx, c_idx, j_idx] <- beta_tibble$`97.5%`[i]
              ##
              coef_array[t_idx, c_idx, j_idx] <- beta_tibble$median[i]
              # ci_lower_array[t_idx, c_idx, j_idx] <- beta_tibble$q5[i]
              # ci_upper_array[t_idx, c_idx, j_idx] <- beta_tibble$q95[i]
              ##
              signif_array[t_idx, c_idx, j_idx] <- beta_tibble$is_significant[i]
              meets_criteria_array[t_idx, c_idx, j_idx] <- beta_tibble$meets_criteria[i]
          
        }
        
        results$coefficients <- coef_array
        results$ci_lower <- ci_lower_array
        results$ci_upper <- ci_upper_array
        results$significant <- signif_array
        results$meets_criteria <- meets_criteria_array
        
        # Print summary
        cat("\n", rep("=", 70), "\n", sep = "")
        cat("BETA COEFFICIENT SUMMARY\n")
        cat(rep("=", 70), "\n", sep = "")
        cat(sprintf("\nSignificance level: %.1f%% CI\n", (1-alpha)*100))
        cat(sprintf("Minimum magnitude threshold: %.3f\n", min_magnitude))
        
        # Overall significance summary
        n_signif <- sum(beta_tibble$is_significant)
        n_borderline <- sum(beta_tibble$is_borderline & !beta_tibble$is_significant)
        n_large_borderline <- sum(beta_tibble$is_borderline & beta_tibble$has_min_magnitude & !beta_tibble$is_significant)
        n_total <- nrow(beta_tibble)
        
        cat(sprintf("\nTotal significant effects (%.0f%% CI): %d/%d (%.1f%%)\n", 
                    (1-alpha)*100, n_signif, n_total, 100 * n_signif / n_total))
        
        if (show_borderline) {
          cat(sprintf("Total with |β| ≥ %.3f: %d\n", 
                      min_magnitude, sum(beta_tibble$has_min_magnitude)))
          cat(sprintf("Large magnitude but not significant: %d\n", 
                      sum(beta_tibble$has_min_magnitude & !beta_tibble$is_significant)))
        }
        
        # Summary by test
        cat("\n", "Significant Effects by Test:\n", sep = "")
        cat(rep("-", 50), "\n", sep = "")
        
        for (t in 1:n_tests) {
          test_data <- beta_tibble %>% filter(test_idx == t)
          n_sig_test <- sum(test_data$is_significant)
          
          cat(sprintf("\n%s: %d significant effects\n", test_names[t], n_sig_test))
          
          if (n_sig_test > 0) {
            sig_effects <- test_data %>% 
              filter(is_significant) %>%
              arrange(class_idx, cov_idx)
            
            for (i in 1:nrow(sig_effects)) {
              cat(sprintf("  - %s, %s: β = %.3f [%.3f, %.3f]\n",
                          sig_effects$class_name[i],
                          sig_effects$covariate_name[i],
                          sig_effects$median[i],
                          sig_effects$q5[i],
                          sig_effects$q95[i]))
            }
          }
        }
        
        # Summary by covariate
        cat("\n\n", "Significant Effects by Covariate:\n", sep = "")
        cat(rep("-", 50), "\n", sep = "")
        
        for (j in 1:n_covariates) {
          cov_data <- beta_tibble %>% filter(cov_idx == j)
          sig_cov_data <- cov_data %>% filter(is_significant)
          
          if (nrow(sig_cov_data) > 0) {
            cat(sprintf("\n%s:\n", covariate_names[j]))
            
            for (i in 1:nrow(sig_cov_data)) {
              cat(sprintf("  - %s, %s: β = %.3f [%.3f, %.3f]\n",
                          sig_cov_data$test_name[i],
                          sig_cov_data$class_name[i],
                          sig_cov_data$median[i],
                          sig_cov_data$q5[i],
                          sig_cov_data$q95[i]))
            }
          }
        }
        
        # Add summary of large but non-significant effects
        if (show_borderline) {
          cat("\n\n", "Large Effects (|β| ≥ ", min_magnitude, ") Not Reaching Significance:\n", sep = "")
          cat(rep("-", 50), "\n", sep = "")
          
          large_nonsig <- beta_tibble %>%
            filter(has_min_magnitude & !is_significant & is_borderline & covariate_name != "intercept") %>%
            arrange(desc(abs(median)))
          
          if (nrow(large_nonsig) > 0) {
            for (i in 1:min(nrow(large_nonsig), 10)) {  # Show top 10
              cat(sprintf("%s %s, %s: β = %.3f [%.3f, %.3f]\n",
                          large_nonsig$test_name[i],
                          large_nonsig$class_name[i],
                          large_nonsig$covariate_name[i],
                          large_nonsig$median[i],
                          large_nonsig$q5[i],
                          large_nonsig$q95[i]))
            }
          } else {
            cat("None\n")
          }
        }
        
        
        # In the summary, note the borderline threshold
        cat(sprintf("\nBorderline threshold: Effects within %.3f of significance\n", borderline_threshold))
        
        # Significance matrix
        cat("\n\n", "Significance Matrix (* = significant at ", (1-alpha)*100, "% level):\n", sep = "")
        cat(rep("-", 70), "\n", sep = "")
        
        # Create header
        cat(sprintf("%-20s", "Test/Class"))
        for (j in 1:min(n_covariates, 7)) {
          cat(sprintf("%-10s", substr(covariate_names[j], 1, 9)))
        }
        cat("\n")
        cat(rep("-", 70), "\n", sep = "")
        
        # Print significance indicators
        for (t in 1:n_tests) {
          for (c in 1:n_classes) {
            cat(sprintf("%-20s", paste(test_names[t], "/", substr(class_names[c], 1, 3))))
            
            for (j in 1:min(n_covariates, 7)) {
              if (signif_array[t, c, j]) {
                cat(sprintf("%-10s", "*"))
              } else {
                cat(sprintf("%-10s", ""))
              }
            }
            cat("\n")
          }
        }
        
        # Store large non-significant effects
        results$large_nonsignificant <- large_nonsig
        
        invisible(results)
  
}



        
#' Helper function to create a publication-ready table
#' R_fn_create_coefficient_table
#' @export 
R_fn_create_coefficient_table <- function(mapped_results, 
                                     tests_to_include = NULL,
                                     classes_to_include = NULL,
                                     covariates_to_include = NULL,
                                     digits = 3) {
  
  # Filter the full table based on selections
  table_data <- mapped_results$full_table
  
  if (!is.null(tests_to_include)) {
    table_data <- table_data %>% filter(test_name %in% tests_to_include)
  }
  
  if (!is.null(classes_to_include)) {
    table_data <- table_data %>% filter(class_name %in% classes_to_include)
  }
  
  if (!is.null(covariates_to_include)) {
    table_data <- table_data %>% filter(covariate_name %in% covariates_to_include)
  }
  
  # Format for publication
  pub_table <- table_data %>%
    mutate(
      estimate_ci = sprintf("%.3f [%.3f, %.3f]%s", 
                            mean, `2.5%`, `97.5%`,
                            ifelse(is_significant, "*", "")),
      effect_size = case_when(
        abs(mean) < 0.2 ~ "Small",
        abs(mean) < 0.5 ~ "Medium",
        abs(mean) < 0.8 ~ "Large",
        TRUE ~ "Very Large"
      )
    ) %>%
    select(test_name, class_name, covariate_name, estimate_ci, effect_size) %>%
    pivot_wider(names_from = class_name, 
                values_from = c(estimate_ci, effect_size),
                names_sep = " - ")
  
  return(pub_table)

}
        

 


 
#' R_fn_summarize_covariates_general
#' @export 
R_fn_summarize_covariates_general <- function(X_list, 
                                              test_names = NULL,
                                              continuous_vars = NULL,
                                              binary_vars = NULL,
                                              categorical_vars = NULL,
                                              logit_vars = NULL) {
  
        # Auto-detect test names if not provided
        if (is.null(test_names)) {
          test_names <- paste0("Test_", 1:length(X_list))
        }
        
        # Get column names from first non-empty matrix
        col_names <- NULL
        for (X in X_list) {
          if (!is.null(X) && nrow(X) > 0) {
            col_names <- colnames(X)
            break
          }
        }
        
        if (is.null(col_names)) {
          stop("No column names found in X_list")
        }
        
        # Initialize results
        results <- list()
        
        # Process each test
        for (t in 1:length(X_list)) {
          X_test <- X_list[[t]]
          
          if (is.null(X_test) || nrow(X_test) == 0) {
            cat("\nTest", t, ":", test_names[t], "- No data\n")
            next
          }
          
          # IMPORTANT: Keep ALL rows to get correct proportions
          n_total_studies <- nrow(X_test)
          
          # Find which studies actually have this test
          non_zero_rows <- rowSums(abs(X_test)) > 0
          n_studies_with_test <- sum(non_zero_rows)
          
          # Basic summary
          cat("\n", rep("=", 60), "\n", sep = "")
          cat("Test", t, ":", test_names[t], "\n")
          cat(rep("=", 60), "\n", sep = "")
          cat("Studies with this test:", n_studies_with_test, "/", n_total_studies, "\n\n")
          
          # Summarize each covariate
          cat("Covariate Distributions (among ALL", n_total_studies, "studies):\n")
          cat(rep("-", 40), "\n", sep = "")
          
          cov_summary <- list()
          
          # Auto-detect variable types if not specified
          if (is.null(continuous_vars) && is.null(binary_vars)) {
            continuous_vars <- c()
            binary_vars <- c()
            
            for (j in 2:ncol(X_test)) {  # Skip intercept
              col_data <- X_test[, j]
              col_data_nonzero <- col_data[col_data != 0]  # Exclude zeros for type detection
              
              if (length(unique(col_data_nonzero)) > 2) {
                continuous_vars <- c(continuous_vars, col_names[j])
              } else {
                binary_vars <- c(binary_vars, col_names[j])
              }
            }
          }
          
          # Process each covariate (skip intercept)
          for (j in 2:ncol(X_test)) {
            cov_name <- col_names[j]
            values <- X_test[, j]
            
            # Determine variable type
            if (cov_name %in% continuous_vars) {
              # Continuous variable - only summarize non-zero values
              non_zero_vals <- values[values != 0]
              if (length(non_zero_vals) > 0) {
                cat(sprintf("%-30s: Mean = %6.3f, SD = %6.3f (n=%d)\n",
                            cov_name, mean(non_zero_vals), sd(non_zero_vals), length(non_zero_vals)))
                
                # If it's a logit variable, show on probability scale
                if (!is.null(logit_vars) && cov_name %in% logit_vars) {
                  prob_values <- plogis(non_zero_vals)
                  cat(sprintf("%-30s: Mean = %6.1f%%, Range = [%5.1f%%, %5.1f%%]\n",
                              "  (Probability scale)", mean(prob_values)*100, 
                              min(prob_values)*100, max(prob_values)*100))
                }
              }
              
            } else if (cov_name %in% binary_vars) {
              # Binary variable - count 1s among ALL studies
              n_ones <- sum(values == 1)
              prop_ones <- n_ones / n_total_studies
              cat(sprintf("%-30s: %d/%d (%5.1f%%)\n", 
                          cov_name, n_ones, n_total_studies, prop_ones * 100))
            }
          }
          
          # Now show distributions among studies WITH this test
          if (n_studies_with_test > 0 && n_studies_with_test < n_total_studies) {
            cat("\n", "Among studies WITH this test (n=", n_studies_with_test, "):\n", sep = "")
            cat(rep("-", 40), "\n", sep = "")
            
            X_active <- X_test[non_zero_rows, , drop = FALSE]
            
            for (j in 2:ncol(X_active)) {
              cov_name <- col_names[j]
              values_active <- X_active[, j]
              
              if (cov_name %in% continuous_vars) {
                cat(sprintf("%-30s: Mean = %6.3f, SD = %6.3f\n",
                            cov_name, mean(values_active), sd(values_active)))
              } else if (cov_name %in% binary_vars) {
                n_ones_active <- sum(values_active == 1)
                prop_ones_active <- n_ones_active / n_studies_with_test
                cat(sprintf("%-30s: %d/%d (%5.1f%%)\n", 
                            cov_name, n_ones_active, n_studies_with_test, prop_ones_active * 100))
              }
            }
          }
          
          # Store results
          results[[test_names[t]]] <- list(
            n_total_studies = n_total_studies,
            n_studies_with_test = n_studies_with_test,
            cov_summary = cov_summary,
            X_full = X_test,
            X_active = if(n_studies_with_test > 0) X_test[non_zero_rows, , drop = FALSE] else NULL
          )
        }
        
        invisible(results)
  
}




#' Helper function to generate specific scenarios for sROC plots
#' R_fn_summarize_covariates_by_test
#' @keywords internal
#' @export
R_fn_generate_sroc_scenarios <- function( summary_results,
                                          max_scenarios = 10) {
  
            # Extract test names
            test_names <- names(summary_results)[!names(summary_results) %in% "scenario_recommendations"]
            
            # Initialize scenarios list
            scenarios <- list()
            
            # 1. Baseline scenario for each test
            baseline_scenario <- list()
            for (test in test_names) {
              if (!is.null(summary_results[[test]]$baseline_values)) {
                baseline_scenario[[test]] <- summary_results[[test]]$baseline_values
              }
            }
            
            scenarios$baseline <- list(
              name = "Baseline (most common values)",
              values = baseline_scenario,
              description = "Most common or median values for all covariates"
            )
            
            # 2. Scenarios varying key covariates
            if (!is.null(summary_results$scenario_recommendations$key_covariates)) {
              key_covs <- summary_results$scenario_recommendations$key_covariates
              
              for (cov in key_covs) {
                # Determine the type and range of this covariate across tests
                cov_type <- NULL
                cov_values <- list()
                
                for (test in test_names) {
                  if (!is.null(summary_results[[test]]$cov_summary[[cov]])) {
                    info <- summary_results[[test]]$cov_summary[[cov]]
                    cov_type <- info$type
                    
                    if (info$type == "continuous") {
                      cov_values[[test]] <- c(info$q25, info$median, info$q75)
                    } else if (info$type == "binary") {
                      cov_values[[test]] <- info$values
                    }
                  }
                }
                
                # Create scenarios for this covariate
                if (cov_type == "continuous" && length(cov_values) > 0) {
                  # Use quartiles
                  all_values <- sort(unique(unlist(cov_values)))
                  
                  for (val in all_values[c(1, ceiling(length(all_values)/2), length(all_values))]) {
                    scenario_name <- sprintf("%s = %.2f", cov, val)
                    scenario_values <- baseline_scenario
                    
                    # Update the covariate value for each test
                    for (test in test_names) {
                      if (!is.null(scenario_values[[test]][[cov]])) {
                        scenario_values[[test]][[cov]] <- val
                      }
                    }
                    
                    scenarios[[paste0(cov, "_", round(val, 2))]] <- list(
                      name = scenario_name,
                      values = scenario_values,
                      description = sprintf("Varying %s to %.2f, others at baseline", cov, val)
                    )
                  }
                  
                } else if (cov_type == "binary" && length(cov_values) > 0) {
                  # Create scenario for each binary value
                  unique_vals <- sort(unique(unlist(cov_values)))
                  
                  for (val in unique_vals) {
                    scenario_name <- sprintf("%s = %g", cov, val)
                    scenario_values <- baseline_scenario
                    
                    # Update the covariate value for each test
                    for (test in test_names) {
                      if (!is.null(scenario_values[[test]][[cov]])) {
                        scenario_values[[test]][[cov]] <- val
                      }
                    }
                    
                    scenarios[[paste0(cov, "_", val)]] <- list(
                      name = scenario_name,
                      values = scenario_values,
                      description = sprintf("%s set to %g, others at baseline", cov, val)
                    )
                  }
                }
              }
            }
            
            # Limit number of scenarios if requested
            if (length(scenarios) > max_scenarios) {
              scenarios <- scenarios[1:max_scenarios]
              cat(sprintf("\nNote: Limited to first %d scenarios. Total possible: %d\n", 
                          max_scenarios, length(scenarios)))
            }
            
            # Print summary
            cat("\nGenerated Scenarios for sROC Plots:\n")
            cat(rep("=", 50), "\n")
            
            for (i in 1:length(scenarios)) {
              scenario <- scenarios[[i]]
              cat(sprintf("\n%d. %s\n", i, scenario$name))
              cat(sprintf("   %s\n", scenario$description))
            }
            
            # Add usage instructions
            cat("\n", rep("-", 50), "\n", sep = "")
            cat("To use these scenarios in Stan/R:\n")
            cat("1. Extract values: scenario_vals <- scenarios$baseline$values\n")
            cat("2. For test-specific values: scenario_vals$'Test 1'$covariate_name\n")
            cat("3. Create baseline vectors for Stan model accordingly\n")
            
            return(scenarios)
          }
          
          # Function to create a scenario matrix for a specific test
          create_scenario_matrix <- function(scenarios, test_name, covariate_names) {
            # Create a matrix where each row is a scenario
            n_scenarios <- length(scenarios)
            n_covariates <- length(covariate_names)
            
            scenario_matrix <- matrix(NA, nrow = n_scenarios, ncol = n_covariates)
            colnames(scenario_matrix) <- covariate_names
            rownames(scenario_matrix) <- names(scenarios)
            
            for (i in 1:n_scenarios) {
              scenario_values <- scenarios[[i]]$values[[test_name]]
              
              for (j in 1:n_covariates) {
                cov_name <- covariate_names[j]
                if (!is.null(scenario_values[[cov_name]])) {
                  scenario_matrix[i, j] <- scenario_values[[cov_name]]
                }
              }
            }
            
            return(scenario_matrix)
  
}

 





#' R_fn_get_covariate_info_MA
#' @keywords internal
#' @export 
R_fn_get_covariate_info_MA <- function(  intercept_only, 
                                         cov_data, 
                                         model_parameterisation,
                                         n_studies
) {
  
            
            cov_info_list  <- list()
            
            if (intercept_only) { 
              X <- NULL
            } else { 
              X <- cov_data$X
              baseline_case_nd <- cov_data$baseline_case_nd
              baseline_case_d  <- cov_data$baseline_case_d
              baseline_case    <- cov_data$baseline_case
              ##
              cov_info_list$baseline_case_nd <- baseline_case_nd
              cov_info_list$baseline_case_d  <- baseline_case_d
              cov_info_list$baseline_case <- baseline_case
            }
            
            if (is.null(X)) { ## i.e. intercept-only 
                    
                    if (grepl("ord_HSROC|ordinal_HSROC|R&G|HSROC|Gatsonis", model_parameterisation, ignore.case = TRUE)) {
                      
                            X_t <- matrix(1.0, nrow = n_studies, ncol = 1)
                            cov_info_list$X <- X_t
                            ##
                            cov_info_list$baseline_case  <- 1
                            cov_info_list$n_covariates   <- 1
                            cov_info_list$n_covariates_max <- 1
                      
                      
                    } else {   ## For ord_biv or Jones models
                      
                            cov_info_list$n_covariates_nd  <- 1
                            cov_info_list$n_covariates_d   <- 1
                            cov_info_list$n_covariates_max <- 1
                            ##
                            X_nd_t <- matrix(1.0, nrow = n_studies, ncol = 1)
                            cov_info_list$X_nd <- X_nd_t
                            ##
                            X_d_t <- matrix(1.0, nrow = n_studies, ncol = 1)
                            cov_info_list$X_d <- X_d_t
                            ##
                            X <- list(cov_info_list$X_nd, cov_info_list$X_d)
                            ##
                            cov_info_list$baseline_case_nd <- 1
                            cov_info_list$baseline_case_d  <- 1
                      
                      
                    }
                    
            } else {  ## i.e. if covariates
                    
                    if (grepl("ord_HSROC|ordinal_HSROC|R&G|HSROC|Gatsonis", model_parameterisation, ignore.case = TRUE)) {
                      
                            cov_info_list$X <- X
                            cov_info_list$n_covariates  <-  ncol(X)
                      
                      
                    } else {   ## For ord_biv or Jones models
                            
                            cov_info_list$X_nd <- X[[1]] ## n_studies x n_covariates_nd matrix
                            cov_info_list$n_covariates_nd <- ncol(cov_info_list$X_nd)
                            ##
                            cov_info_list$X_d  <- X[[2]] ## n_studies x n_covariates_d matrix
                            cov_info_list$n_covariates_d <- ncol(cov_info_list$X_d)
                            ##
                            cov_info_list$n_covariates_max <- max(cov_info_list$n_covariates_nd, cov_info_list$n_covariates_d)
                      
                      
                    }
              
            }
            
            return(cov_info_list)
  
}




#' R_fn_get_covariate_info_NMA
#' @keywords internal
#' @export
R_fn_get_covariate_info_NMA <- function( intercept_only,
                                         cov_data, 
                                         model_parameterisation,
                                         n_index_tests,
                                         n_studies
) {
  
            cov_info_list  <- list()
            
            if (intercept_only) { 
                    X <- NULL
            } else { 
                    X <- cov_data$X
                    baseline_case_nd <- cov_data$baseline_case_nd
                    baseline_case_d  <- cov_data$baseline_case_d
                    baseline_case    <- cov_data$baseline_case
                    ##
                    cov_info_list$baseline_case_nd <- baseline_case_nd
                    cov_info_list$baseline_case_d  <- baseline_case_d
                    cov_info_list$baseline_case <- baseline_case
            }
            
            if (is.null(X)) { ## i.e. intercept-only 
              
                  if (grepl("ord_HSROC|ordinal_HSROC|R&G|HSROC|Gatsonis", model_parameterisation, ignore.case = TRUE)) {
                    
                          X_t <- matrix(1.0, nrow = n_studies, ncol = 1)
                          cov_info_list$X <- rep(list(X_t), n_index_tests)
                          ##
                          cov_info_list$baseline_case  <- rep(list(1), n_index_tests)
                          cov_info_list$n_covariates   <- rep(1, n_index_tests)
                          cov_info_list$n_covariates_max <- 1
                          
                  } else {   ## For ord_biv or Jones models
                          
                          cov_info_list$n_covariates_nd <- rep(1, n_index_tests)
                          cov_info_list$n_covariates_d  <- rep(1, n_index_tests)
                          cov_info_list$n_covariates_max <- 1
                          ##
                          X_nd_t <- matrix(1.0, nrow = n_studies, ncol = 1)
                          cov_info_list$X_nd <- rep(list(X_nd_t), n_index_tests)
                          ##
                          X_d_t <- matrix(1.0, nrow = n_studies, ncol = 1)
                          cov_info_list$X_d <- rep(list(X_d_t), n_index_tests)
                          ##
                          cov_info_list$X <- list(cov_info_list$X_nd, cov_info_list$X_d)
                          ##
                          cov_info_list$baseline_case_nd <- rep(list(1), n_index_tests)
                          cov_info_list$baseline_case_d  <- rep(list(1), n_index_tests)
                    
                  }
              
            } else { ## i.e. if covariates
              
                  if (grepl("ord_HSROC|ordinal_HSROC|R&G|HSROC|Gatsonis", model_parameterisation, ignore.case = TRUE)) {
                          
                          cov_info_list$X <- X
                          cov_info_list$n_covariates  <-  c() ##  ncol(X)
                          ##
                          for (t in 1:n_index_tests) {
                            cov_info_list$n_covariates[t] <- ncol(X[[t]])
                          }
                          ##
                          cov_info_list$n_covariates_max <- max(cov_info_list$n_covariates)
                    
                  } else {   ## For ord_biv or Jones models
                    
                          cov_info_list$X_nd <- X[[1]] ## n_studies x n_covariates_nd matrix
                          cov_info_list$X_d  <- X[[2]] ## n_studies x n_covariates_nd matrix
                          cov_info_list$X <- X
                          ##
                          cov_info_list$n_covariates_nd <- c()
                          cov_info_list$n_covariates_d  <- c()
                          ##
                          for (t in 1:n_index_tests) {
                            cov_info_list$n_covariates_nd[t] <- ncol(cov_info_list$X_nd[[t]])
                            cov_info_list$n_covariates_d[t]  <- ncol(cov_info_list$X_d[[t]])
                          }
                          ##
                          cov_info_list$n_covariates_max <- max(max(cov_info_list$n_covariates_nd), max(cov_info_list$n_covariates_d))
                    
                  }
              
            }
            
            return(cov_info_list)

}













#' Function to convert multi-level categorical to dummy variables
#' R_fn_create_dummy_variables
#' @export
R_fn_create_dummy_variables <- function( x, 
                                         var_name, 
                                         reference_level = NULL) {
          # Convert to factor if not already
          if (!is.factor(x)) {
            x <- factor(x)
          }
          
          levels_x <- levels(x)
          n_levels <- length(levels_x)
          
          if (n_levels < 2) {
            stop(paste("Variable", var_name, "has fewer than 2 levels"))
          }
          
          # Set reference level (default is first alphabetically)
          if (is.null(reference_level)) {
            reference_level <- levels_x[1]
          }
          
          # Create dummy variables (n_levels - 1)
          dummy_matrix <- matrix(0, nrow = length(x), ncol = n_levels - 1)
          colnames(dummy_matrix) <- paste0(var_name, "_", levels_x[-which(levels_x == reference_level)])
          
          # Fill in dummy values
          for (i in 1:(n_levels - 1)) {
            level_idx <- which(levels_x == levels_x[-which(levels_x == reference_level)][i])
            dummy_matrix[x == levels_x[level_idx], i] <- 1
          }
          
          return(dummy_matrix)
          
}


 




#' R_fn_prepare_NMA_covariates
#' @export
R_fn_prepare_NMA_covariates <- function(cov_data, 
                                        test_data_list,
                                        test_names,
                                        continuous_covs,
                                        binary_covs,
                                        categorical_covs,
                                        standardize_continuous = TRUE,
                                        handle_missing = "median",
                                        baseline_values = NULL) {
  
          require(dplyr)
          
          # Get unique studies and ensure proper ordering
          unique_studies <- cov_data %>% 
            distinct(Study) %>% 
            arrange(Study) %>% 
            pull(Study)
          
          n_studies <- length(unique_studies)
          n_index_tests <- length(test_names)
          
          # Prepare covariate data - one row per study
          cov_data_clean <- cov_data %>%
            distinct(Study, .keep_all = TRUE) %>%
            arrange(Study)
          
          # Initialize the design matrix with intercept
          design_matrix <- data.frame(intercept = rep(1, n_studies))
          
          # Process continuous variables
          for (var in continuous_covs) {
            if (var %in% names(cov_data_clean)) {
              # Convert to numeric
              var_values <- as.numeric(cov_data_clean[[var]])
              
              # Handle missing values
              if (handle_missing == "median") {
                median_val <- median(var_values, na.rm = TRUE)
                var_values[is.na(var_values)] <- median_val
              } else if (handle_missing == "mean") {
                mean_val <- mean(var_values, na.rm = TRUE)
                var_values[is.na(var_values)] <- mean_val
              }
              
              # Standardize if requested
              if (standardize_continuous) {
                mean_val <- mean(var_values, na.rm = TRUE)
                sd_val <- sd(var_values, na.rm = TRUE)
                if (sd_val > 0) {
                  var_values <- (var_values - mean_val) / sd_val
                }
              }
              
              design_matrix[[var]] <- var_values
            }
          }
          
          # Process binary variables
          if (!is.null(binary_covs)) {
            for (var in binary_covs) {
              if (var %in% names(cov_data_clean)) {
                # Convert Yes/No to 1/0
                var_values <- cov_data_clean[[var]]
                var_values <- ifelse(toupper(as.character(var_values)) %in% c("YES", "1") | var_values == 1, 1, 0)
                design_matrix[[var]] <- var_values
              }
            }
          }
          
          # Process categorical variables (create dummy variables)
          if (!is.null(categorical_covs) && length(categorical_covs) > 0) {
            for (var in categorical_covs) {
              if (var %in% names(cov_data_clean)) {
                cat_values <- cov_data_clean[[var]]
                
                # Debug print
                cat("\nProcessing categorical variable:", var, "\n")
                cat("Unique values:", unique(cat_values), "\n")
                
                # Convert to factor
                if (!is.factor(cat_values)) {
                  cat_values <- factor(cat_values)
                }
                
                # Only create dummies if more than 2 levels
                n_levels <- length(levels(cat_values))
                if (n_levels > 1) {
                  # Create dummy variables manually
                  for (i in 2:n_levels) {  # Skip first level (reference)
                    level_name <- levels(cat_values)[i]
                    dummy_name <- paste0(var, "_", level_name)
                    design_matrix[[dummy_name]] <- as.numeric(cat_values == level_name)
                  }
                }
              }
            }
          }
          
          # Convert to numeric matrix
          X_mat <- as.matrix(design_matrix)
          # Ensure all values are numeric
          storage.mode(X_mat) <- "numeric"
          
          # Create X_nd and X_d arrays (4D for Stan)
          n_max_covariates <- ncol(X_mat)
          
          # Initialize the arrays
          X_nd_array <- array(0, dim = c(n_index_tests, n_studies, n_max_covariates))
          X_d_array <- array(0, dim = c(n_index_tests, n_studies, n_max_covariates))
          
          # Fill the arrays
          for (t in 1:n_index_tests) {
            X_nd_array[t, , ] <- X_mat
            X_d_array[t, , ] <- X_mat
          }
          
          # Create lists for compatibility
          X_nd_list <- list()
          X_d_list <- list()
          for (t in 1:n_index_tests) {
            X_nd_list[[t]] <- X_mat
            X_d_list[[t]] <- X_mat
          }
          
          n_covariates_nd <- rep(ncol(X_mat), n_index_tests)
          n_covariates_d <- rep(ncol(X_mat), n_index_tests)
          n_covariates_max <- ncol(X_mat)
          
          # Create baseline case vectors
          if (is.null(baseline_values)) {
            baseline_case_nd <- list()
            baseline_case_d <- list()
            
            baseline_vec <- rep(0, ncol(X_mat))
            baseline_vec[1] <- 1  # intercept
            
            for (t in 1:n_index_tests) {
              baseline_case_nd[[t]] <- baseline_vec
              baseline_case_d[[t]] <- baseline_vec
            }
          } else {
            baseline_case_nd <- baseline_values$nd
            baseline_case_d <- baseline_values$d
          }
          
          # Create indicator matrix for which tests are in which studies
          indicator_index_test_in_study <- matrix(0, nrow = n_studies, ncol = n_index_tests)
          
          for (t in 1:n_index_tests) {
            test_data <- test_data_list[[t]]
            studies_with_test <- unique(test_data$Study)
            
            for (s in 1:n_studies) {
              if (unique_studies[s] %in% studies_with_test) {
                indicator_index_test_in_study[s, t] <- 1
              }
            }
          }
          
          # Create prior matrices
          prior_beta_mu_mean <- list()
          prior_beta_mu_SD <- list()
          
          for (t in 1:n_index_tests) {
            prior_beta_mu_mean[[t]] <- matrix(0, nrow = 2, ncol = n_covariates_max)
            prior_beta_mu_SD[[t]] <- matrix(2, nrow = 2, ncol = n_covariates_max)
            prior_beta_mu_SD[[t]][, 1] <- 1  # Tighter prior for intercept
          }
          
          # Print summary
          cat("\n=== Covariate Preparation Summary ===\n")
          cat("Number of studies:", n_studies, "\n")
          cat("Number of index tests:", n_index_tests, "\n")
          cat("Total number of covariates (including intercept and dummies):", ncol(X_mat), "\n")
          cat("Covariate names:", paste(colnames(X_mat), collapse = ", "), "\n")
          cat("\nFirst few rows of covariate matrix:\n")
          print(head(X_mat, 10))
          cat("\nBaseline case values:", baseline_vec, "\n")
          
          # Return all components
          return(list(
            n_studies = n_studies,
            n_index_tests = n_index_tests,
            n_covariates_nd = n_covariates_nd,
            n_covariates_d = n_covariates_d,
            n_covariates_max = n_covariates_max,
            X_nd = X_nd_list,
            X_d = X_d_list,
            baseline_case_nd = baseline_case_nd,
            baseline_case_d = baseline_case_d,
            indicator_index_test_in_study = indicator_index_test_in_study,
            prior_beta_mu_mean = prior_beta_mu_mean,
            prior_beta_mu_SD = prior_beta_mu_SD,
            study_mapping = unique_studies,
            covariate_names = colnames(X_mat),
            design_matrix = X_mat
          ))
  
}




#' R_fn_pad_covariate_matrices
#' @export
R_fn_pad_covariate_matrices <- function( X_list, 
                                         fill_value = 0 ## 0 won't contribute to Xbeta/lp
) {
  
      # Find the maximum number of rows across all matrices
      max_rows <- max(sapply(X_list, nrow))
      
      # Pad each matrix to have max_rows
      X_list_padded <- lapply(X_list, function(X) {
        current_rows <- nrow(X)
        current_cols <- ncol(X)
        
        if (current_rows < max_rows) {
          # Create padding matrix filled with fill_value
          padding_rows <- max_rows - current_rows
          padding_matrix <- matrix(fill_value, 
                                   nrow = padding_rows, 
                                   ncol = current_cols)
          
          # Set column names if they exist
          if (!is.null(colnames(X))) {
            colnames(padding_matrix) <- colnames(X)
          }
          
          # Combine original matrix with padding
          X_padded <- rbind(X, padding_matrix)
        } else {
          # Matrix already has max rows, no padding needed
          X_padded <- X
        }
        
        return(X_padded)
      })
  
  return(X_list_padded)
  
}




#' expand_covariates_to_all_studies
#' @export
expand_covariates_to_all_studies <- function(X, 
                                             indicator_index_test_in_study, 
                                             study_mappings = NULL,
                                             missing_value = 0 ## 0 won't contribute to Xbeta/lp
) {
  
            # First, check the structure of X
            if (!is.list(X) || length(X) != 2) {
              stop("X must be a list of length 2 (for nd and d)")
            }
            
            n_tests <- length(X[[1]])  # Number of tests
            
            # Verify both disease groups have same number of tests
            if (length(X[[2]]) != n_tests) {
              stop("X[[1]] and X[[2]] must have the same number of tests")
            }
            
            n_studies <- nrow(indicator_index_test_in_study)
            
            # Get number of covariates (check first non-NULL matrix)
            n_covariates <- NULL
            for (disease_group in 1:2) {
              for (t in 1:n_tests) {
                if (!is.null(X[[disease_group]][[t]]) && is.matrix(X[[disease_group]][[t]])) {
                  n_covariates <- ncol(X[[disease_group]][[t]])
                  break
                }
              }
              if (!is.null(n_covariates)) break
            }
            
            if (is.null(n_covariates)) {
              stop("Could not determine number of covariates from X")
            }
            
            # Check if already in expanded format
            is_already_expanded <- TRUE
            for (disease_group in 1:2) {
              for (t in 1:n_tests) {
                # Check if this element exists and is a matrix
                if (!is.null(X[[disease_group]][[t]]) && is.matrix(X[[disease_group]][[t]])) {
                  if (nrow(X[[disease_group]][[t]]) != n_studies) {
                    is_already_expanded <- FALSE
                    break
                  }
                }
              }
              if (!is_already_expanded) break
            }
            
            # If already expanded, just return as is
            if (is_already_expanded) {
              return(X)
            }
            
            # Calculate n_studies_per_test from the indicator matrix
            n_studies_per_test <- colSums(indicator_index_test_in_study)
            
            # If study_mappings not provided, create them
            if (is.null(study_mappings)) {
              study_mappings <- list()
              for (t in 1:n_tests) {
                studies_with_test <- which(indicator_index_test_in_study[, t] == 1)
                study_mappings[[t]] <- studies_with_test
              }
            }
            
            # Create expanded X with all studies
            X_expanded <- list()
            
            for (disease_group in 1:2) {  # nd and d
              X_expanded[[disease_group]] <- list()
              
              for (t in 1:n_tests) {
                # Initialize with missing values for all studies
                X_expanded[[disease_group]][[t]] <- matrix(missing_value, 
                                                           nrow = n_studies, 
                                                           ncol = n_covariates)
                
                # Check if we have data for this test and disease group
                if (!is.null(X[[disease_group]][[t]]) && is.matrix(X[[disease_group]][[t]])) {
                  
                  # Preserve column names if they exist
                  if (!is.null(colnames(X[[disease_group]][[t]]))) {
                    colnames(X_expanded[[disease_group]][[t]]) <- colnames(X[[disease_group]][[t]])
                  }
                  
                  # Get the actual number of studies for this test
                  n_studies_for_test <- n_studies_per_test[t]
                  n_rows_in_data <- nrow(X[[disease_group]][[t]])
                  
                  # Only process up to the minimum of actual studies or rows in data
                  n_rows_to_process <- min(n_studies_for_test, n_rows_in_data)
                  
                  # Fill in data for studies that have this test
                  for (row_idx in 1:n_rows_to_process) {
                    global_study_idx <- study_mappings[[t]][row_idx]
                    
                    # Get the row data
                    row_data <- X[[disease_group]][[t]][row_idx, ]
                    
                    # Check if we have valid data for this row
                    if (!all(is.na(row_data) | row_data == missing_value)) {
                      X_expanded[[disease_group]][[t]][global_study_idx, ] <- row_data
                    }
                  }
                }
              }
            }
            
            return(X_expanded)
  
}


#






 
#' process_covs
#' @export
process_covs <- function(covs, 
                         n_tests) {
  
      if (is.null(covs)) return(rep(list(NULL), n_tests))
      if (!is.list(covs)) {
        # Single vector -> use for all tests
        return(rep(list(covs), n_tests))
      } else if (length(covs) == n_tests) {
        # Already per-test
        return(covs)
      } else {
        stop("Covariate list length must match n_tests")
      }
  
}


#' process_covariates_for_group
#' @export
process_covariates_for_group <- function(cov_data,
                                         test_studies, 
                                         continuous, 
                                         binary, 
                                         categorical,
                                         all_categorical_levels
) {
        # Filter and prepare data
        # cov_data_clean <- cov_data %>%
        #   filter(Study %in% test_studies) %>%
        #   distinct(Study, .keep_all = TRUE) %>%
        #   arrange(Study)
        
        cov_data_clean <- cov_data %>%
          distinct(Study, .keep_all = TRUE) %>%
          arrange(Study)
        
        # Convert missing to 999
        for (col in names(cov_data_clean)) {
          if (col != "Study") {
            if (is.numeric(cov_data_clean[[col]])) {
              cov_data_clean[[col]][is.na(cov_data_clean[[col]])] <- 999
              cov_data_clean[[col]][is.nan(cov_data_clean[[col]])] <- 999
              cov_data_clean[[col]][is.infinite(cov_data_clean[[col]])] <- 999
            } else {
              cov_data_clean[[col]][is.na(cov_data_clean[[col]])] <- 999
              cov_data_clean[[col]][cov_data_clean[[col]] == ""] <- 999
            }
          }
        }
        
        # Create design matrix with intercept
        design_matrix <- data.frame(intercept = rep(1, nrow(cov_data_clean)))
        
        # Process continuous variables
        for (var in continuous) {
          if (var %in% names(cov_data_clean)) {
            var_values <- as.numeric(cov_data_clean[[var]])
            
            # Handle 999 as missing (instead of NA median imputation)
            missing_idx <- which(var_values == 999)
            if (length(missing_idx) > 0) {
              var_values[missing_idx] <- 999 # Set to 999, Stan will estimate
              cat("Note:", length(missing_idx), "missing values (999) in", var, "will be estimated by Stan\n")
            }
            design_matrix[[var]] <- var_values
          }
        }
        
        # Process binary variables
        for (var in binary) {
          if (var %in% names(cov_data_clean)) {
            var_values <- cov_data_clean[[var]]
            
            # First convert any YES/NO strings to 1/0
            var_values <- ifelse(toupper(as.character(var_values)) %in% c("YES", "1"), 1,
                                 ifelse(toupper(as.character(var_values)) %in% c("NO", "0"), 0,
                                        var_values))  # Keep 999 as is
            
            # Now handle 999s
            if (any(var_values == 999, na.rm = TRUE)) {
              # Replace 999 with 0 in main variable
              design_matrix[[var]] <- ifelse(var_values == 999, 0, var_values)
              # Add missing indicator
              design_matrix[[paste0(var, "_missing")]] <- as.numeric(var_values == 999)
              cat("Note: Missing values (999) in", var, "coded as separate indicator\n")
            } else {
              # No missing values - just use as is
              design_matrix[[var]] <- var_values
            }
          }
        }
        
        # Process categorical variables - handle 999 DURING dummy creation
        if (!is.null(categorical) && length(categorical) > 0) {
          
          for (var in categorical) {
            if (var %in% names(cov_data_clean)) {
              cat_values <- cov_data_clean[[var]]
              
              # Check if there are any 999s
              has_missing <- any(cat_values == 999, na.rm = TRUE)
              
              # Get ALL levels from all_categorical_levels
              if (var %in% names(all_categorical_levels)) {
                unique_levels <- factor(all_categorical_levels[[var]])
              } else {
                # This shouldn't happen but fallback just in case
                unique_levels <- unique(cat_values[cat_values != 999])
                if (!is.factor(unique_levels)) {
                  unique_levels <- factor(unique_levels)
                }
              }
              
              # Create dummies for each level (except reference)
              n_levels <- length(levels(unique_levels))
              if (n_levels > 1) {
                for (i in 2:n_levels) {
                  level_name <- levels(unique_levels)[i]
                  dummy_name <- paste0(var, "_", level_name)
                  # When original is 999, dummy gets 0
                  design_matrix[[dummy_name]] <- ifelse(cat_values == 999, 0,
                                                        as.numeric(cat_values == level_name))
                }
              }
              
              # Add missing indicator if needed
              if (has_missing) {
                design_matrix[[paste0(var, "_miss")]] <- as.numeric(cat_values == 999)
              }
            }
          }
        }
        
        # Convert to matrix and return
        X_mat <- as.matrix(design_matrix)
        storage.mode(X_mat) <- "numeric"
        
        return(X_mat)
  
}





#' create_subsetted_covariates_for_tests
#' @export
create_subsetted_covariates_for_tests <- function(     cov_data,
                                                       indicator_index_test_in_study, 
                                                       test_names, 
                                                       continuous_covs,
                                                       binary_covs,
                                                       categorical_covs,
                                                       center_and_scale_cts) {
          require(dplyr)
          
          # Handle single test (non-NMA) case
          if (is.vector(indicator_index_test_in_study)) {
            indicator_index_test_in_study <- matrix(indicator_index_test_in_study, ncol = 1)
          }
          
          n_tests <- ncol(indicator_index_test_in_study)
          n_max_studies <- nrow(indicator_index_test_in_study)
          
          # PARSE INPUT into standardized format
          if (is.data.frame(cov_data)) {
            cov_data_nd <- rep(list(cov_data), n_tests)
            cov_data_d <- rep(list(cov_data), n_tests)
          } else if (is.list(cov_data)) {
            if (length(cov_data) == 2 && is.list(cov_data[[1]])) {
              cov_data_nd <- cov_data[[1]]
              cov_data_d <- cov_data[[2]]
            } else if (length(cov_data) == n_tests && is.data.frame(cov_data[[1]])) {
              cov_data_nd <- cov_data
              cov_data_d <- cov_data
            } else {
              stop("Invalid cov_data format")
            }
          } else {
            stop("cov_data must be tibble or list")
          }
          
          continuous_per_test <- process_covs(continuous_covs, n_tests)
          binary_per_test <- process_covs(binary_covs, n_tests)
          categorical_per_test <- process_covs(categorical_covs, n_tests)
          
          # Collect ALL categorical levels across ALL tests for SHARED covariates
          all_categorical_levels <- list()
          all_categorical_vars <- unique(unlist(categorical_per_test))
          
          for (var in all_categorical_vars[!is.null(all_categorical_vars)]) {
            all_vals <- c()
            
            for (t in 1:n_tests) {
              if (!is.null(cov_data_nd[[t]]) && var %in% names(cov_data_nd[[t]])) {
                vals <- cov_data_nd[[t]][[var]]
                all_vals <- c(all_vals, vals[vals != 999 & !is.na(vals)])
              }
              if (!is.null(cov_data_d[[t]]) && var %in% names(cov_data_d[[t]])) {
                vals <- cov_data_d[[t]][[var]]
                all_vals <- c(all_vals, vals[vals != 999 & !is.na(vals)])
              }
            }
            
            if (length(all_vals) > 0) {
              all_categorical_levels[[var]] <- sort(unique(all_vals))
            }
          }
          
          print("DEBUG - Categorical levels collected:")
          print(all_categorical_levels)
          
          # Process each test
          X_nd_list <- list()
          X_d_list <- list()
          n_studies_per_test <- numeric(n_tests)
          study_mappings_per_test <- list()
 
          
          # for (t in 1:n_tests) {
          #   studies_in_test <- which(indicator_index_test_in_study[, t] == 1)
          #   n_studies_per_test[t] <- length(studies_in_test)
          #   
          #   study_mappings_per_test[[t]] <- data.frame(
          #     original_study = studies_in_test,
          #     new_row_index = 1:length(studies_in_test)
          #   )
          #   
          #   # Process ND
          #   X_nd_list[[t]] <- process_covariates_for_group(
          #     cov_data_nd[[t]],
          #     studies_in_test,
          #     continuous_per_test[[t]],
          #     binary_per_test[[t]],
          #     categorical_per_test[[t]],
          #     all_categorical_levels  # Passed internally
          #   )
          #   
          #   # Process D  
          #   X_d_list[[t]] <- process_covariates_for_group(
          #     cov_data_d[[t]],
          #     studies_in_test,
          #     continuous_per_test[[t]],
          #     binary_per_test[[t]],
          #     categorical_per_test[[t]],
          #     all_categorical_levels  # Passed internally
          #   )
          # }
          
          
          for (t in 1:n_tests) {
            # Get the actual number of studies in this test's data
            n_studies_this_test <- nrow(cov_data_nd[[t]])
            n_studies_per_test[t] <- n_studies_this_test
            
            # Don't need test_studies anymore - just pass NULL
            X_nd_list[[t]] <- process_covariates_for_group(
              cov_data_nd[[t]],
              NULL,  # Don't need this anymore
              continuous_per_test[[t]],
              binary_per_test[[t]],
              categorical_per_test[[t]],
              all_categorical_levels
            )
            
            # Same for D
            X_d_list[[t]] <- process_covariates_for_group(
              cov_data_d[[t]],
              NULL,  # Don't need this anymore
              continuous_per_test[[t]],
              binary_per_test[[t]],
              categorical_per_test[[t]],
              all_categorical_levels
            )
          }

          
          # Get dimensions
          n_covariates_nd <- sapply(X_nd_list, ncol)
          n_covariates_d <- sapply(X_d_list, ncol)
          n_covariates_max <- max(c(n_covariates_nd, n_covariates_d))
          
          # Create baseline vectors
          baseline_vec <- rep(0, n_covariates_max)
          baseline_vec[1] <- 1  # intercept
          
          baseline_case_nd <- rep(list(baseline_vec), n_tests)
          baseline_case_d <- rep(list(baseline_vec), n_tests)
          
          # Center and scale continuous variables
          if (center_and_scale_cts) {
            all_continuous_vars <- unique(unlist(continuous_per_test))
            
            for (var in all_continuous_vars[!is.null(all_continuous_vars)]) {
              # Collect all NON-MISSING values across all tests
              all_values <- c()
              for (t in 1:n_tests) {
                col_idx <- which(colnames(X_nd_list[[t]]) == var)
                if (length(col_idx) > 0) {
                  vals <- X_nd_list[[t]][, col_idx]
                  all_values <- c(all_values, vals[vals != 999])
                }
                col_idx <- which(colnames(X_d_list[[t]]) == var)
                if (length(col_idx) > 0) {
                  vals <- X_d_list[[t]][, col_idx]
                  all_values <- c(all_values, vals[vals != 999])
                }
              }
              
              if (length(all_values) > 0) {
                global_mean <- mean(all_values)
                global_sd <- sd(all_values)
                
                # Apply scaling but KEEP 999s as 999
                for (t in 1:n_tests) {
                  # ND
                  col_idx <- which(colnames(X_nd_list[[t]]) == var)
                  if (length(col_idx) > 0) {
                    vals <- X_nd_list[[t]][, col_idx]
                    non_missing <- vals != 999
                    vals[non_missing] <- (vals[non_missing] - global_mean) / global_sd
                    X_nd_list[[t]][, col_idx] <- vals
                  }
                  
                  # D
                  col_idx <- which(colnames(X_d_list[[t]]) == var)
                  if (length(col_idx) > 0) {
                    vals <- X_d_list[[t]][, col_idx]
                    non_missing <- vals != 999
                    vals[non_missing] <- (vals[non_missing] - global_mean) / global_sd
                    X_d_list[[t]][, col_idx] <- vals
                  }
                }
              }
            }
          }
          
          # Print summary
          cat("\n=== Covariate Preparation Summary ===\n")
          cat("Total studies:", n_max_studies, "\n")
          cat("Studies per test:", n_studies_per_test, "\n")
          cat("Max covariates:", n_covariates_max, "\n")
          
          return(list(
            n_studies = n_max_studies,
            n_index_tests = n_tests,
            n_covariates_nd = n_covariates_nd,
            n_covariates_d = n_covariates_d,
            n_covariates_max = n_covariates_max,
            X_nd = X_nd_list,
            X_d = X_d_list,
            baseline_case_nd = baseline_case_nd,
            baseline_case_d = baseline_case_d,
            indicator_index_test_in_study = indicator_index_test_in_study,
            study_mappings_per_test = study_mappings_per_test,
            n_studies_per_test = n_studies_per_test,
            continuous_covs = continuous_per_test,
            binary_covs = binary_per_test,
            categorical_covs = categorical_per_test
          ))
  
}








#' prepare_NMA_data
#' @export
prepare_NMA_data <- function(test_data_list, 
                             test_names = NULL,
                             max_scores = NULL
) {
  
        n_index_tests <- length(test_data_list)
        
        # Auto-detect test names if not provided
        if (is.null(test_names)) {
          test_names <- names(test_data_list)
          if (is.null(test_names)) test_names <- paste0("Test_", 1:n_index_tests)
        }
        
        min_scores <- rep(0, n_index_tests) ## MUST be zero!
        # Default score ranges if not provided
        if (is.null(max_scores)) {
          # You'll need to specify these or detect them from data
          stop("max_scores must be provided")
        }
        
        # Get unique studies per test
        study_ids_by_test <- lapply(test_data_list, function(x) unique(x$Study))
        all_study_ids <- sort(unique(unlist(study_ids_by_test)))
        n_studies <- length(all_study_ids)
        
        # Create indicator matrix
        indicator_matrix <- matrix(0, nrow = n_studies, ncol = n_index_tests)
        study_id_to_row <- setNames(1:n_studies, all_study_ids)
        
        for (t in 1:n_index_tests) {
          for (s in study_ids_by_test[[t]]) {
            row_idx <- study_id_to_row[as.character(s)]
            indicator_matrix[row_idx, t] <- 1
          }
        }
        
        # Initialize x arrays
        x <- list()
        n_thr <- numeric(n_index_tests)
        
        for (t in 1:n_index_tests) {
          # Transform data
          data_transformed <- transform_questionnaire_data_by_disease(
            data = test_data_list[[t]],
            min_score = min_scores[t],
            max_score = max_scores[t]
          )
          
          n_thr[t] <- ncol(data_transformed$diseased) - 1
          
          # Convert to integers
          x_nd <- convert_to_integers_simple(data_transformed$nondiseased)
          x_d <- convert_to_integers_simple(data_transformed$diseased)
          
          # Initialize matrices
          x[[t]] <- list(
            matrix(-1, nrow = n_studies, ncol = n_thr[t] + 1),
            matrix(-1, nrow = n_studies, ncol = n_thr[t] + 1)
          )
          
          # Fill in data
          unique_studies <- unique(test_data_list[[t]]$Study)
          for (s_local in 1:length(unique_studies)) {
            study_id <- unique_studies[s_local]
            s_row <- study_id_to_row[as.character(study_id)]
            
            x[[t]][[1]][s_row, ] <- x_nd[s_local, ]
            x[[t]][[2]][s_row, ] <- x_d[s_local, ]
          }
        }
        
        # Return clean list
        list(
          n_studies = n_studies,
          n_index_tests = n_index_tests,
          n_thr = n_thr,
          n_cat = n_thr + 1,
          n_thr_max = max(n_thr),
          n_cat_max = max(n_thr + 1),
          n_studies_per_test = sapply(study_ids_by_test, length),
          x = x,
          indicator_index_test_in_study = indicator_matrix,
          test_names = test_names
        )
  
}






#' prepare_covariate_data
#' @export
prepare_covariate_data <- function(X, 
                                   indicator_matrix, 
                                   n_index_tests, 
                                   baseline_case = NULL) {
  
        n_covariates <- ncol(X[[1]][[1]])
        
        # Expand covariates to all studies if needed
        X_expanded <- MetaOrdDTA:::expand_covariates_to_all_studies(
          X, 
          indicator_index_test_in_study = indicator_matrix
        )
        
        # Create baseline if not provided
        if (is.null(baseline_case)) {
          baseline_case <- rep(0, n_covariates)
          baseline_case[1] <- 1  # Intercept
        }
        
        list(
          X = X_expanded,
          baseline_case_nd = rep(list(baseline_case), n_index_tests),
          baseline_case_d = rep(list(baseline_case), n_index_tests)
        )
  
}





#' validate_X_matrices
#' @export
validate_X_matrices <- function(X, 
                                tolerance = 1e-10,
                                test_specific_covs = NULL,
                                indicator_index_test_in_study = NULL,
                                is_X_expanded = NULL
) {
        
        cat("=== VALIDATING X MATRICES ===\n\n")
        
        n_disease_groups <- length(X)
        n_tests <- length(X[[1]])
        
        cat("Structure check:\n")
        cat("- Disease groups:", n_disease_groups, "\n")
        cat("- Tests:", n_tests, "\n")
        
        # Require explicit specification if indicator provided
        if (!is.null(indicator_index_test_in_study) && is.null(is_X_expanded)) {
          stop("When indicator_index_test_in_study is provided, you must specify is_X_expanded=TRUE or FALSE")
        }
        
        if (!is.null(is_X_expanded)) {
          cat("- Format:", ifelse(is_X_expanded, "EXPANDED", "NON-EXPANDED"), "\n")
        }
        cat("\n")
        

        # Check dimensions
        cat("Dimensions per test:\n")
        for (t in 1:n_tests) {
          n_rows_nd <- nrow(X[[1]][[t]])
          n_cols_nd <- ncol(X[[1]][[t]])
          n_rows_d <- nrow(X[[2]][[t]])
          n_cols_d <- ncol(X[[2]][[t]])
          
          cat("  Test", t, ": ND =", n_rows_nd, "x", n_cols_nd, 
              ", D =", n_rows_d, "x", n_cols_d, "\n")
          
          if (n_rows_nd != n_rows_d || n_cols_nd != n_cols_d) {
            cat("    ⚠️ WARNING: ND and D dimensions differ!\n")
          }
        }
        
        covariate_names <- colnames(X[[1]][[1]])
        cat("\nCovariate names:", paste(covariate_names, collapse = ", "), "\n\n")
        
        # Check ND vs D
        cat("=== Checking ND vs D consistency within tests ===\n")
        nd_d_match <- TRUE
        for (t in 1:n_tests) {
          max_diff <- max(abs(X[[1]][[t]] - X[[2]][[t]]))
          if (max_diff > tolerance) {
            cat("  ❌ Test", t, ": ND and D differ (max diff =", max_diff, ")\n")
            nd_d_match <- FALSE
          } else {
            cat("  ✓ Test", t, ": ND and D identical\n")
          }
        }
        
        # Check cross-test consistency if indicator provided
        if (!is.null(indicator_index_test_in_study)) {
          
                cat("\n=== Checking cross-test consistency for shared studies ===\n")
                
                # Derive study mappings from indicator
                study_mappings <- list()
                for (t in 1:n_tests) {
                  study_mappings[[t]] <- which(indicator_index_test_in_study[, t] == 1)
                }
                
                # Determine covariate types
                if (is.null(test_specific_covs)) {
                  test_specific_covs <- character(0)
                }
                study_level_covs <- setdiff(covariate_names, c("intercept", test_specific_covs))
                
                cat("Test-specific covariates (allowed to vary):", 
                    ifelse(length(test_specific_covs) > 0, 
                           paste(test_specific_covs, collapse = ", "), 
                           "none"), "\n")
                cat("Study-level covariates (must be consistent):", 
                    ifelse(length(study_level_covs) > 0,
                           paste(study_level_covs, collapse = ", "),
                           "none"), "\n\n")
                
                # Find studies in multiple tests
                all_studies <- 1:nrow(indicator_index_test_in_study)
                cross_test_valid <- TRUE
                
                for (study_id in all_studies) {
                  tests_with_study <- which(indicator_index_test_in_study[study_id, ] == 1)
                  
                  if (length(tests_with_study) > 1) {
                    # This study appears in multiple tests
                    for (cov in study_level_covs) {
                      ref_test <- tests_with_study[1]
                      
                      if (is_X_expanded) {
                        # Expanded: study_id IS the row number
                        ref_val <- X[[1]][[ref_test]][study_id, cov]
                        
                        for (t in tests_with_study[-1]) {
                          test_val <- X[[1]][[t]][study_id, cov]
                          
                          if (abs(ref_val - test_val) > tolerance) {
                            cat("  ❌ Study", study_id, ", Covariate '", cov, 
                                "': Test", ref_test, "=", ref_val, 
                                ", Test", t, "=", test_val, "\n")
                            cross_test_valid <- FALSE
                          }
                        }
                      } else {
                        # Non-expanded: need to find which row this study is
                        # Build study mappings once outside the loop
                        if (!exists("study_mappings")) {
                          study_mappings <- list()
                          for (t in 1:n_tests) {
                            study_mappings[[t]] <- which(indicator_index_test_in_study[, t] == 1)
                          }
                        }
                        
                        ref_row <- which(study_mappings[[ref_test]] == study_id)
                        ref_val <- X[[1]][[ref_test]][ref_row, cov]
                        
                        for (t in tests_with_study[-1]) {
                          test_row <- which(study_mappings[[t]] == study_id)
                          test_val <- X[[1]][[t]][test_row, cov]
                          
                          if (abs(ref_val - test_val) > tolerance) {
                            cat("  ❌ Study", study_id, ", Covariate '", cov, 
                                "': Test", ref_test, "=", ref_val, 
                                ", Test", t, "=", test_val, "\n")
                            cross_test_valid <- FALSE
                          }
                        }
                      }
                    }
                  }
                }
                
                if (cross_test_valid && length(study_level_covs) > 0) {
                  cat("  ✓ All study-level covariates consistent across tests\n")
                }
                
        }
        
        return(nd_d_match)
  
}


# validate_X_matrices <- function(X, 
#                                 tolerance = 1e-10) {
#   
#           cat("=== VALIDATING X MATRICES ===\n\n")
#           
#           # Basic structure checks
#           if (!is.list(X) || length(X) != 2) {
#             stop("X must be a list of length 2 (for non-diseased and diseased)")
#           }
#           
#           n_disease_groups <- length(X)
#           n_tests <- length(X[[1]])
#           
#           if (length(X[[2]]) != n_tests) {
#             stop("X[[1]] and X[[2]] must have same number of tests")
#           }
#           
#           # Get dimensions
#           n_studies <- nrow(X[[1]][[1]])
#           n_covariates <- ncol(X[[1]][[1]])
#           covariate_names <- colnames(X[[1]][[1]])
#           
#           cat("Structure check:\n")
#           cat("- Disease groups:", n_disease_groups, "\n")
#           cat("- Tests:", n_tests, "\n")
#           cat("- Studies:", n_studies, "\n")
#           cat("- Covariates:", n_covariates, "\n")
#           cat("- Covariate names:", paste(covariate_names, collapse = ", "), "\n\n")
#           
#           # Check dimensions consistency
#           cat("Checking dimensions consistency...\n")
#           for (dg in 1:2) {
#             for (test in 1:n_tests) {
#               if (nrow(X[[dg]][[test]]) != n_studies) {
#                 stop(sprintf("Inconsistent n_studies: DG %d, Test %d has %d rows, expected %d", 
#                              dg, test, nrow(X[[dg]][[test]]), n_studies))
#               }
#               if (ncol(X[[dg]][[test]]) != n_covariates) {
#                 stop(sprintf("Inconsistent n_covariates: DG %d, Test %d has %d cols, expected %d", 
#                              dg, test, ncol(X[[dg]][[test]]), n_covariates))
#               }
#             }
#           }
#           cat("✓ All dimensions consistent\n\n")
#           
#           # CRITICAL CHECK: Study-level covariates must be identical across tests
#           cat("Checking study-level covariate consistency...\n")
#           
#           all_valid <- TRUE
#           inconsistent_studies <- c()
#           inconsistent_covariates <- c()
#           
#           for (study in 1:n_studies) {
#             # Get reference values from first non-zero occurrence
#             reference_vals <- NULL
#             reference_test <- NULL
#             
#             # Find first test where this study has data (intercept != 0)
#             for (test in 1:n_tests) {
#               if (X[[1]][[test]][study, "intercept"] != 0) {
#                 reference_vals <- X[[1]][[test]][study, ]
#                 reference_test <- test
#                 break
#               }
#             }
#             
#             # If no reference found, study doesn't appear in any test
#             if (is.null(reference_vals)) next
#             
#             # Check consistency across all other tests
#             for (test in 1:n_tests) {
#               # Skip if study not in this test
#               if (X[[1]][[test]][study, "intercept"] == 0) next
#               
#               # Compare with reference
#               for (cov in covariate_names) {
#                 # Check both disease groups
#                 for (dg in 1:2) {
#                   val_ref <- reference_vals[cov]
#                   val_test <- X[[dg]][[test]][study, cov]
#                   
#                   if (abs(val_ref - val_test) > tolerance) {
#                     all_valid <- FALSE
#                     inconsistent_studies <- c(inconsistent_studies, study)
#                     inconsistent_covariates <- c(inconsistent_covariates, cov)
#                     
#                     cat(sprintf("  ❌ Study %d, Covariate '%s': Test %d = %.6f, Test %d = %.6f (diff = %.6f)\n",
#                                 study, cov, reference_test, val_ref, test, val_test, 
#                                 abs(val_ref - val_test)))
#                   }
#                 }
#               }
#             }
#           }
#           
#           if (all_valid) {
#             cat("✓ All study-level covariates are consistent across tests!\n\n")
#           } else {
#             cat("\n❌ VALIDATION FAILED!\n")
#             cat("Inconsistent studies:", unique(inconsistent_studies), "\n")
#             cat("Inconsistent covariates:", unique(inconsistent_covariates), "\n\n")
#           }
#           
#           # Additional check: Are diseased and non-diseased X matrices identical?
#           cat("Checking if diseased and non-diseased X matrices are identical...\n")
#           matrices_identical <- TRUE
#           for (test in 1:n_tests) {
#             if (!identical(X[[1]][[test]], X[[2]][[test]])) {
#               matrices_identical <- FALSE
#               cat(sprintf("  Test %d: Diseased and non-diseased differ\n", test))
#             }
#           }
#           
#           if (matrices_identical) {
#             cat("✓ Diseased and non-diseased X matrices are identical (as expected)\n\n")
#           }
#           
#           # Summary statistics
#           cat("Summary statistics:\n")
#           for (cov in covariate_names) {
#             if (cov != "intercept") {
#               # Collect all non-zero values
#               all_vals <- c()
#               for (test in 1:n_tests) {
#                 vals <- X[[1]][[test]][X[[1]][[test]][, "intercept"] != 0, cov]
#                 all_vals <- c(all_vals, vals)
#               }
#               
#               if (length(unique(all_vals)) <= 10) {
#                 cat(sprintf("  %s: %d unique values: %s\n", 
#                             cov, length(unique(all_vals)), 
#                             paste(sort(unique(all_vals)), collapse = ", ")))
#               } else {
#                 cat(sprintf("  %s: mean = %.3f, sd = %.3f, range = [%.3f, %.3f]\n",
#                             cov, mean(all_vals), sd(all_vals), 
#                             min(all_vals), max(all_vals)))
#               }
#             }
#           }
#           
#           return(all_valid)
#   
# }
# 
#  
# 
# 
















check_X_column_order <- function(X, 
                                 expected_order = NULL) {
  
          # Default expected order
          if (is.null(expected_order)) {
            expected_order <- c("intercept", 
                                "logit_prev_GAD", 
                                "low_RoB_QUADAS_clean",
                                "study_setting_2", 
                                "study_setting_3",
                                "Ref_test_clean_SCID", 
                                "Ref_test_clean_Structured")
          }
          
          cat("\n=== CHECKING X COLUMN ORDER ===\n")
          
          all_good <- TRUE
          
          # Check each disease group and test
          for (dg in 1:length(X)) {
            for (test in 1:length(X[[dg]])) {
              if (!is.null(X[[dg]][[test]]) && ncol(X[[dg]][[test]]) > 0) {
                current_cols <- colnames(X[[dg]][[test]])
                
                # Check if columns are in expected order
                expected_subset <- expected_order[expected_order %in% current_cols]
                
                if (!identical(current_cols, expected_subset)) {
                  cat(sprintf("\n❌ WRONG ORDER in X[[%d]][[%d]]:\n", dg, test))
                  cat("  Current: ", paste(current_cols, collapse = ", "), "\n")
                  cat("  Expected:", paste(expected_subset, collapse = ", "), "\n")
                  all_good <- FALSE
                }
              }
            }
          }
          
          if (all_good) {
            # Show what we have
            example_cols <- colnames(X[[1]][[1]])
            cat("\n✅ All columns in correct order!\n")
            cat("Columns present:", paste(example_cols, collapse = ", "), "\n")
            
            # Show positions
            cat("\nColumn positions:\n")
            for (i in 1:length(example_cols)) {
              cat(sprintf("  Position %d: %s\n", i, example_cols[i]))
            }
          }
          
          return(all_good)
  
}







