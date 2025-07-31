


#' resize_init_list
#' @export
resize_init_list <- function(init_lists_per_chain, 
                             n_chains_new) {
  
        n_chains_old <- length(init_lists_per_chain)
        
        if (n_chains_new == n_chains_old) {
          # No change needed
          return(init_lists_per_chain)
        } else if (n_chains_new < n_chains_old) {
          # Fewer chains - just take first n_chains_new
          return(init_lists_per_chain[1:n_chains_new])
        } else {
          # More chains - repeat cyclically
          indices <- rep(1:n_chains_old, length.out = n_chains_new)
          return(init_lists_per_chain[indices])
        }
  
}



#' create_cutpoint_index_matrix
#' @export
create_cutpoint_index_matrix <- function( data_matrix, 
                                          max_score
) {
  
        # Get dimensions of the input matrix
        n_studies <- nrow(data_matrix)
        n_thr <- ncol(data_matrix)
        
        # Create a new matrix of the same dimensions
        cutpoint_index <- matrix(-1, nrow = n_studies, ncol = n_thr)
        
        # Retain the row and column names
        rownames(cutpoint_index) <- rownames(data_matrix)
        colnames(cutpoint_index) <- colnames(data_matrix)
        
        # For each row (study)
        for (i in 1:n_studies) {
          # For each column (threshold)
          for (j in 1:n_thr) {
            # If the data isn't missing (-1), replace with the cutpoint index (j-1)
            if (data_matrix[i, j] != -1) {
              cutpoint_index[i, j] <- j - 1  # 0-based indexing: 0, 1, 2, ..., max_score
            }
          }
        }
        
        return(cutpoint_index)
  
}




#' compare_indicator_matrices
#' @export
compare_indicator_matrices <- function(sim_matrix, 
                                       real_matrix, 
                                       test_names = c("Test1", "Test2", "Test3", "Test4"),
                                       detailed = TRUE) {
  
        # Ensure same number of columns
        if (ncol(sim_matrix) != ncol(real_matrix)) {
          stop("Matrices must have same number of columns (tests)")
        }
        
        # 1. Test frequencies
        sim_freq <- colSums(sim_matrix) / nrow(sim_matrix)
        real_freq <- colSums(real_matrix) / nrow(real_matrix)
        freq_diff <- abs(sim_freq - real_freq)
        
        # 2. Overlap proportions
        sim_overlap <- compute_study_overlap_prop(sim_matrix)
        real_overlap <- compute_study_overlap_prop(real_matrix)
        overlap_diff <- abs(sim_overlap - real_overlap)
        
        # 3. Test combinations
        get_combinations <- function(mat) {
          # Convert each row to a string pattern
          patterns <- apply(mat, 1, function(x) paste(x, collapse = ""))
          table(patterns) / nrow(mat)
        }
        
        sim_combos <- get_combinations(sim_matrix)
        real_combos <- get_combinations(real_matrix)
        
        # 4. Pairwise test associations (correlation)
        sim_cor <- cor(sim_matrix)
        real_cor <- cor(real_matrix)
        cor_diff <- abs(sim_cor - real_cor)
        
        # 5. Number of tests per study distribution
        sim_n_tests <- table(rowSums(sim_matrix)) / nrow(sim_matrix)
        real_n_tests <- table(rowSums(real_matrix)) / nrow(real_matrix)
        
        # 6. Calculate similarity metrics
        # Hellinger distance for combination frequencies
        all_patterns <- unique(c(names(sim_combos), names(real_combos)))
        sim_combo_vec <- rep(0, length(all_patterns))
        real_combo_vec <- rep(0, length(all_patterns))
        names(sim_combo_vec) <- names(real_combo_vec) <- all_patterns
        sim_combo_vec[names(sim_combos)] <- sim_combos
        real_combo_vec[names(real_combos)] <- real_combos
        
        hellinger_dist <- sqrt(0.5 * sum((sqrt(sim_combo_vec) - sqrt(real_combo_vec))^2))
        
        # 7. Overall similarity score (0-100)
        freq_score <- 100 * (1 - mean(freq_diff))
        overlap_score <- 100 * (1 - overlap_diff)
        cor_score <- 100 * (1 - mean(cor_diff[upper.tri(cor_diff)]))
        combo_score <- 100 * (1 - hellinger_dist)
        
        overall_score <- mean(c(freq_score, overlap_score, cor_score, combo_score))
        
        # Print results
        cat("===== INDICATOR MATRIX COMPARISON =====\n\n")
        
        cat("1. TEST FREQUENCIES:\n")
        freq_table <- data.frame(
          Test = test_names,
          Simulated = sprintf("%.1f%%", 100 * sim_freq),
          Real = sprintf("%.1f%%", 100 * real_freq),
          Difference = sprintf("%.1f%%", 100 * freq_diff)
        )
        print(freq_table, row.names = FALSE)
        
        cat("\n2. OVERLAP PROPORTION:\n")
        cat(sprintf("Simulated: %.1f%%\n", 100 * sim_overlap))
        cat(sprintf("Real: %.1f%%\n", 100 * real_overlap))
        cat(sprintf("Difference: %.1f%%\n", 100 * overlap_diff))
        
        cat("\n3. TESTS PER STUDY:\n")
        cat("Simulated distribution:\n")
        print(round(100 * sim_n_tests, 1))
        cat("Real distribution:\n")
        print(round(100 * real_n_tests, 1))
        
        if (detailed) {
          cat("\n4. TOP TEST COMBINATIONS:\n")
          cat("\nSimulated (top 5):\n")
          sim_top <- head(sort(sim_combos, decreasing = TRUE), 5)
          for (i in seq_along(sim_top)) {
            pattern <- names(sim_top)[i]
            tests <- which(strsplit(pattern, "")[[1]] == "1")
            test_combo <- paste(test_names[tests], collapse = "+")
            if (test_combo == "") test_combo <- "None"
            cat(sprintf("  %s: %.1f%%\n", test_combo, 100 * sim_top[i]))
          }
          
          cat("\nReal (top 5):\n")
          real_top <- head(sort(real_combos, decreasing = TRUE), 5)
          for (i in seq_along(real_top)) {
            pattern <- names(real_top)[i]
            tests <- which(strsplit(pattern, "")[[1]] == "1")
            test_combo <- paste(test_names[tests], collapse = "+")
            if (test_combo == "") test_combo <- "None"
            cat(sprintf("  %s: %.1f%%\n", test_combo, 100 * real_top[i]))
          }
          
          cat("\n5. TEST CORRELATIONS:\n")
          cat("Simulated:\n")
          print(round(sim_cor, 2))
          cat("\nReal:\n")
          print(round(real_cor, 2))
        }
        
        cat("\n===== SIMILARITY SCORES =====\n")
        cat(sprintf("Frequency similarity: %.0f%%\n", freq_score))
        cat(sprintf("Overlap similarity: %.0f%%\n", overlap_score))
        cat(sprintf("Correlation similarity: %.0f%%\n", cor_score))
        cat(sprintf("Combination similarity: %.0f%%\n", combo_score))
        cat(sprintf("\nOVERALL SIMILARITY: %.0f%%\n", overall_score))
        
        # Return detailed results invisibly
        invisible(list(
          test_frequencies = list(sim = sim_freq, real = real_freq, diff = freq_diff),
          overlap = list(sim = sim_overlap, real = real_overlap, diff = overlap_diff),
          combinations = list(sim = sim_combos, real = real_combos),
          correlations = list(sim = sim_cor, real = real_cor, diff = cor_diff),
          tests_per_study = list(sim = sim_n_tests, real = real_n_tests),
          scores = list(
            frequency = freq_score,
            overlap = overlap_score,
            correlation = cor_score,
            combination = combo_score,
            overall = overall_score
          ),
          hellinger_distance = hellinger_dist
        ))
        
}

# Helper function if compute_study_overlap_prop isn't already defined
if (!exists("compute_study_overlap_prop")) {
  compute_study_overlap_prop <- function(indicator_matrix) {
    tests_per_study <- rowSums(indicator_matrix)
    studies_with_multiple_tests <- sum(tests_per_study >= 2)
    total_studies <- nrow(indicator_matrix)
    study_overlap_prop <- studies_with_multiple_tests / total_studies
    return(study_overlap_prop)
  }
}





#' compute_study_overlap_prop
#' @export
compute_study_overlap_prop <- function(indicator_matrix) {
  
          # Count number of tests per study (row sums)
          tests_per_study <- rowSums(indicator_matrix)
          
          # Count studies with multiple tests (2 or more)
          studies_with_multiple_tests <- sum(tests_per_study >= 2)
          
          # Total number of studies
          total_studies <- nrow(indicator_matrix)
          
          # Compute proportion
          study_overlap_prop <- studies_with_multiple_tests / total_studies
          
          # Print summary
          cat("Total studies:", total_studies, "\n")
          cat("Studies with 1 test:", sum(tests_per_study == 1), "\n")
          cat("Studies with 2+ tests:", studies_with_multiple_tests, "\n")
          cat("Study overlap proportion:", round(study_overlap_prop, 3), "\n")
          
          return(study_overlap_prop)
          
}

# Usage:
# study_overlap_prop <- compute_study_overlap_prop(indicator_index_test_in_study)


#' create_indicator_index_test_in_study
#' @export
create_indicator_index_test_in_study <- function( x, 
                                                  study_missing_indicator = -1) {
          
          # Helper function to check if a row has valid data
          has_valid_data <- function(row) {
            !all(row == study_missing_indicator) && 
              !(row[1] != study_missing_indicator && all(row[-1] == study_missing_indicator))
          }
          
          # Get number of tests and studies
          n_tests <- length(x)
          n_studies <- nrow(x[[1]][[1]])
          
          # Initialize indicator matrix
          indicator_matrix <- matrix(0, nrow = n_studies, ncol = n_tests)
          
          # For each test
          for (test in 1:n_tests) {
            # Check both diseased and non-diseased groups
            has_data_nd <- apply(x[[test]][[1]], 1, has_valid_data)
            has_data_d <- apply(x[[test]][[2]], 1, has_valid_data)
            
            # Study has test if either group has data
            indicator_matrix[, test] <- as.numeric(has_data_nd | has_data_d)
          }
          
          # Add column names
          colnames(indicator_matrix) <- paste0("Test", 1:n_tests)
          
          return(indicator_matrix)
  
}





#' arrange_missing_values
#' @export
arrange_missing_values <- function(data_list) {
        
        # Create a copy of the input list to avoid modifying the original
        result_list <- list(
          non_diseased = data_list[[1]],
          diseased = data_list[[2]]
        )
        
        # Process each matrix in the list
        for (i in 1:2) {
          matrix_data <- result_list[[i]]
          
          # Process each row
          for (row in 1:nrow(matrix_data)) {
            # Get current row data
            row_data <- matrix_data[row, ]
            
            # Split into regular values and missing values (-1)
            regular_values <- row_data[row_data != -1]
            missing_values <- row_data[row_data == -1]
            
            # Combine regular values first, then missing values
            new_row <- c(regular_values, missing_values)
            
            # Replace the original row
            matrix_data[row, ] <- new_row
          }
          
          result_list[[i]] <- matrix_data
        }
        
        return(result_list)
  
}






#' apply_missingness_pattern
#' @export
apply_missingness_pattern <- function( x_complete, 
                                       x_pattern, 
                                       enforce_consecutive_missingness = FALSE) {
  
        
        # Check if inputs are valid
        if (!is.list(x_complete) || !is.list(x_pattern) || length(x_complete) != 2 || length(x_pattern) != 2) {
          stop("Both inputs must be lists with exactly 2 elements (non-diseased and diseased matrices)")
        }
        
        # Create copies of the complete data to modify
        result <- list(
          non_diseased = x_complete[[1]],
          diseased = x_complete[[2]]
        )
        
        # Get dimensions of the first complete matrix to use as reference
        n_rows_complete <- nrow(x_complete[[1]])
        n_cols_complete <- ncol(x_complete[[1]])
        
        # Calculate the missingness pattern based on the first pattern matrix
        missingness_indices <- list()
        pattern_matrix <- x_pattern[[1]]
        n_rows_pattern <- nrow(pattern_matrix)
        n_cols_pattern <- ncol(pattern_matrix)
        
        # Find missing indices for each column
        for (col in 2:min(n_cols_complete, n_cols_pattern)) {
          # Find rows with -1 in the pattern matrix
          missing_rows_in_pattern <- which(pattern_matrix[, col] == -1)
          
          if (length(missing_rows_in_pattern) > 0) {
            # Calculate proportion of missing values in this column of pattern
            missing_proportion <- length(missing_rows_in_pattern) / n_rows_pattern
            
            # Calculate how many rows should be missing in the complete matrix
            n_missing_in_complete <- round(n_rows_complete * missing_proportion)
            
            if (n_missing_in_complete > 0) {
              # Randomly select rows to make missing
              rows_to_make_missing <- sample(1:n_rows_complete, n_missing_in_complete)
              
              # Store these indices
              missingness_indices[[col]] <- rows_to_make_missing
            } else {
              missingness_indices[[col]] <- integer(0)
            }
          } else {
            missingness_indices[[col]] <- integer(0)
          }
        }
        
        # Apply the same missingness pattern to both populations
        for (pop_idx in 1:2) {
          # For each column in the complete matrix
          for (col in 2:min(n_cols_complete, n_cols_pattern)) {
            if (length(missingness_indices[[col]]) > 0) {
              # Set values to -1 using the same indices
              result[[pop_idx]][missingness_indices[[col]], col] <- -1
            }
          }
          
          # Optionally apply the consecutive missingness rule
          if (enforce_consecutive_missingness) {
            # If there's a -1 in column j, all columns to the right should also be -1
            for (row in 1:nrow(result[[pop_idx]])) {
              missing_found <- FALSE
              for (col in 2:ncol(result[[pop_idx]])) {
                if (missing_found || result[[pop_idx]][row, col] == -1) {
                  missing_found <- TRUE
                  result[[pop_idx]][row, col] <- -1
                }
              }
            }
          }
        }
        
        return(result)
  
}




#' Wrapper for apply_missingness_pattern that handles NMA data
#' 
#' @param x_complete List with 2 matrices - complete simulated NMA data (may have rows of all -999)
#' @param x_pattern List with 2 matrices - pattern data
#' @param enforce_consecutive_missingness Passed to original function
#' @param x_complete_missing_indicator Value indicating "no test" in x_complete (default -999)
#' @param x_pattern_missing_indicator Value indicating "missing threshold" in x_pattern (default -1)
#' 
#' @export
apply_missingness_pattern_NMA <- function(x_complete, 
                                          x_pattern,
                                          enforce_consecutive_missingness = FALSE,
                                          x_complete_missing_indicator = -999,
                                          x_pattern_missing_indicator = -1) {
        
        # Find rows that have actual test data (not all missing)
        has_test_data <- apply(x_complete[[1]], 1, function(row) {
          any(row != x_complete_missing_indicator)
        })
        
        # Get indices of valid rows and missing rows
        valid_row_indices <- which(has_test_data)
        missing_row_indices <- which(!has_test_data)
        
        # If no valid rows, return original
        if (length(valid_row_indices) == 0) {
          warning("No valid rows found in x_complete. Returning original data.")
          return(x_complete)
        }
        
        # Extract only the valid rows
        x_complete_valid <- list(
          x_complete[[1]][valid_row_indices, , drop = FALSE],
          x_complete[[2]][valid_row_indices, , drop = FALSE]
        )
        
        # Apply the original function to valid data only
        x_valid_with_pattern <- apply_missingness_pattern(
          x_complete = x_complete_valid,
          x_pattern = x_pattern,
          enforce_consecutive_missingness = enforce_consecutive_missingness
        )
        
        # Now reconstruct the full matrices with original row order
        result <- list(
          matrix(x_complete_missing_indicator, 
                 nrow = nrow(x_complete[[1]]), 
                 ncol = ncol(x_complete[[1]])),
          matrix(x_complete_missing_indicator, 
                 nrow = nrow(x_complete[[2]]), 
                 ncol = ncol(x_complete[[2]]))
        )
        
        # Put the processed valid rows back in their original positions
        result[[1]][valid_row_indices, ] <- x_valid_with_pattern[[1]]
        result[[2]][valid_row_indices, ] <- x_valid_with_pattern[[2]]
        
        # Ensure rows that originally had no test remain all missing
        if (length(missing_row_indices) > 0) {
          result[[1]][missing_row_indices, ] <- x_complete_missing_indicator
          result[[2]][missing_row_indices, ] <- x_complete_missing_indicator
        }
        
        # Preserve column names if they exist
        if (!is.null(colnames(x_complete[[1]]))) {
          colnames(result[[1]]) <- colnames(x_complete[[1]])
          colnames(result[[2]]) <- colnames(x_complete[[2]])
        }
        
        # Preserve list names if they exist
        if (!is.null(names(x_complete))) {
          names(result) <- names(x_complete)
        }
        
        return(result)
  
}






#' generate_ordinal_dirichlet_prior
#' @export
generate_ordinal_dirichlet_prior <- function(n_cat, 
                                             max_alpha = 5, 
                                             min_alpha = 1,
                                             left_skew_factor = 0.3) {
        
        # n_cat: Number of categories in the ordinal scale
        # max_alpha: Maximum concentration parameter (for early categories)
        # min_alpha: Minimum concentration parameter (for later categories)
        # left_skew_factor: Controls how quickly values decrease (lower = more left-skewed)
        
        # Generate positions from 0 to 1
        positions <- seq(0, 1, length.out = n_cat)
        
        # Apply a transformation to create the left skew
        # Using a power function where lower values of left_skew_factor create more skew
        transformed_positions <- positions^left_skew_factor
        
        # Normalize to [0,1] range
        transformed_positions <- transformed_positions / max(transformed_positions)
        
        # Map to alpha parameter range
        alpha_values <- max_alpha - (max_alpha - min_alpha) * transformed_positions
        
        # Round to create cleaner values (optional)
        alpha_values <- round(alpha_values, 1)
        
        return(alpha_values)
        
}









#' generate_dual_ordinal_dirichlet_priors
#' @export
generate_dual_ordinal_dirichlet_priors <- function(n_cat, 
                                                   nd_peak = 0.15,          # Where non-diseased peaks (as fraction)
                                                   d_peak = 0.30,           # Where diseased peaks (as fraction)
                                                   nd_max_alpha = 5.0,      # Max alpha for non-diseased
                                                   d_max_alpha = 5.0,       # Max alpha for diseased
                                                   nd_min_alpha = 1.0,      # Min alpha for non-diseased
                                                   d_min_alpha = 1.0,
                                                   smoothness = 0.15) {     # Higher = smoother distribution
  
        # Create sequence for category positions
        positions <- seq(0, 1, length.out = n_cat)
        
        # Generate smoother non-diseased alphas
        nd_spread <- smoothness * 1.2  # Slightly wider for non-diseased
        nd_alpha_values <- nd_min_alpha + (nd_max_alpha - nd_min_alpha) * 
          exp(-((positions - nd_peak)^2) / (2 * nd_spread^2))
        
        # Generate smoother diseased alphas
        d_spread <- smoothness * 1.5   # Even wider for diseased
        d_alpha_values <- d_min_alpha + (d_max_alpha - d_min_alpha) * 
          exp(-((positions - d_peak)^2) / (2 * d_spread^2))
        
        # Round for cleaner values
        nd_alpha_values <- round(nd_alpha_values, 1)
        d_alpha_values <- round(d_alpha_values, 1)
        
        return(list(
          non_diseased = nd_alpha_values,
          diseased = d_alpha_values
        ))
  
}



#' divide_studies
#' @export
divide_studies <- function( n_studies, 
                            n_groups) {
  
        # Calculate base size for each group
        base_size <- floor(n_studies / n_groups)
        
        # Calculate remainder
        remainder <- n_studies %% n_groups
        
        # Create list to store group vectors
        groups <- list()
        
        # Start index
        start_idx <- 1
        
        # For each group
        for (i in 1:n_groups) {
            
            # Determine group size (add 1 to base_size if there's still remainder)
            group_size <- base_size + (if (i <= remainder) 1 else 0)
            
            # Calculate end index
            end_idx <- start_idx + group_size - 1
            
            # Store group vector
            groups[[i]] <- start_idx:end_idx
            
            # Update start index for next group
            start_idx <- end_idx + 1
          
        }
        
        return(groups)
        
}
















#' R_fn_apply_thr_missingness
#' @export
R_fn_apply_thr_missingness_to_test <- function(  data_list_of_mat, 
                                                 thr_combo_vec = NULL, 
                                                 apply_missings_to_these_rows
) {

         data_list_of_mat_original <- data_list_of_mat
          
         for (c in 1:2) {
              data_list_of_mat[[c]] <-  R_fn_apply_thr_missingness_to_mat( data_matrix = data_list_of_mat_original[[c]],
                                                                           thr_combo_vec = thr_combo_vec,
                                                                           apply_missings_to_these_rows = apply_missings_to_these_rows)
         }
          
         return(data_list_of_mat)
  
}





 


#' R_fn_apply_missingness_pattern_PHQ_9
#' @export
R_fn_apply_missingness_pattern_PHQ_9 <- function( x_PHQ9, 
                                                  case_list, 
                                                  seed = 123) {
  
          # set.seed(123)
  
          # Create a copy of the original data
          x_PHQ_missing <- x_PHQ9
          
          # Loop through each case
          for (case_name in names(case_list)) {
                  
                  case <- case_list[[case_name]]
                  studies <- case$studies
                
                  if (case_name == "case_3") {
                    
                          # For case 3, we need to generate random thresholds for each study:
                          for (study_idx in studies) {
                            # Select a random threshold between 8 and 15
                            random_thr <- sample(8:15, 1)
                            
                            # Create mask where everything is -1 except the random threshold
                            for (group in 1:2) { # Loop over diseased and non-diseased groups
                              for (col in 1:ncol(x_PHQ9[[group]])) {
                                if (col != random_thr) {
                                  try({  
                                    x_PHQ_missing[[group]][study_idx, col + 1] <- -1
                                  })
                                }
                              }
                            }
                          }
                    
                  } else if (case_name == "case_4") {
                    
                          ## For case 4, we randomly select one of the threshold combinations for each study:
                          # Extract all threshold combination vectors
                          thr_combo_vecs <- list()
                          for (i in 1:12) {
                            vec_name <- paste0("thr_combo_vec_", i)
                            if (vec_name %in% names(case)) {
                              thr_combo_vecs[[i]] <- case[[vec_name]]
                            }
                          }
                          
                          # Apply a randomly selected threshold combination to each study
                          for (study_idx in studies) {
                            # Randomly select one of the threshold combinations
                            selected_combo_idx <- sample(1:length(thr_combo_vecs), 1)
                            selected_thresholds <- thr_combo_vecs[[selected_combo_idx]]
                            
                            # Create mask where everything is -1 except the selected thresholds
                            for (group in 1:2) {
                              for (col in 1:ncol(x_PHQ9[[group]])) {
                                if (!(col %in% selected_thresholds)) {
                                  try({  
                                    x_PHQ_missing[[group]][study_idx, col + 1] <- -1
                                  })
                                }
                              }
                            }
                          }
                    
                  } else if (case_name %in% c("case_1", "case_2", "case_5")) {
                    
                          ## For cases 1, 2, and 5, the threshold pattern is the same for all studies in the case:
                          ## Extract the threshold vector (it's called thr_combo_vec_1 for these cases):
                          thresholds <- case$thr_combo_vec_1
                          
                          ## Apply the missingness pattern to all studies in this case:
                          for (study_idx in studies) {
                            for (group in 1:2) {
                              for (col in 1:ncol(x_PHQ9[[group]])) {
                                if (!(col %in% thresholds)) {
                                  try({  
                                    x_PHQ_missing[[group]][study_idx, col + 1] <- -1
                                  })
                                }
                              }
                            }
                          }
                    
                  }
            
          }
          
          return(x_PHQ_missing)
  
}








#' convert_to_aggregate_counts
#' @export
convert_to_aggregate_counts <- function(y_list, 
                                        n_studies, 
                                        n_tests_inc_ref,
                                        n_thr) {
  
          ## Add BINARY ref test n_thr (=1) to "n_thr":
          n_thr <- c(1, n_thr)
          ##
          n_tests <- n_tests_inc_ref

          n_total_nd <- n_total_d <- numeric(n_studies)
          x_nd_list <- x_d_list <- list()
          Se_per_study_list <- Sp_per_study_list <- list()
          
          ## Initialise lists:
          n_index_tests <- n_tests - 1
          for (t in 2:n_tests) { ## skip the BINARY reference test (test 1)
            
            n_thr_t <- n_thr[t]
            n_cat_t <- n_thr_t + 1
            ##
            x_nd_list[[t - 1]] <- matrix(NA, n_studies, n_cat_t)
            x_d_list[[t - 1]]  <- matrix(NA, n_studies, n_cat_t)
            ##
            Se_per_study_list[[t - 1]] <-  matrix(NA, n_studies, n_thr_t + 1)
            Sp_per_study_list[[t - 1]] <-  matrix(NA, n_studies, n_thr_t + 1)
            
          }
          
          for (s in 1:n_studies) {
            
                    study_data <- y_list[[s]]
                    disease_status <- study_data[, 1] # Col 1 is disease status / reference test
                    
                    ## Get n_total_d and n_total_nd:
                    n_total_d[s]  <- sum(disease_status == 1)
                    n_total_nd[s] <- sum(disease_status == 0)
                
            for (t in 2:n_tests) { ## skip the BINARY reference test (test 1)
              
                    n_thr_t <- n_thr[t]
                    n_cat_t <- n_thr_t + 1
                    index_test_results   <- study_data[, t] ## - 1 ## minus one so {0, .., 27} for e.g. PHQ-9.
                     
                    print(paste("min = ", min(index_test_results)))
                    print(paste("max = ", max(index_test_results)))
                    
                    # x_d_list[[t - 1]][s, 1]  <- n_total_d[s]
                    # x_nd_list[[t - 1]][s, 1] <- n_total_nd[s]
                    # Get counts for each threshold
                    for (k in 1:(n_thr_t + 1)) {
                          ##
                          ## True-positives (TP's):
                          ##
                          x_d_list[[t - 1]][s, k]  <- sum(index_test_results[disease_status == 1] >= k - 1) # Pr(testing POSITIVE at threshold k)
                          Se_per_study_list[[t - 1]][s, k] <-   x_d_list[[t - 1]][s, k]  / n_total_d[s]
                          ##
                          ## False-positives (FP's):
                          ##
                          x_nd_list[[t - 1]][s, k] <- sum(index_test_results[disease_status == 0] >= k - 1) # Pr(testing POSITIVE at threshold k)
                          Sp_per_study_list[[t - 1]][s, k] <-   x_nd_list[[t - 1]][s, k]  / n_total_d[s]
                          
                    }
                    # for (k in 1:(n_thr_t + 1)) {
                    # 
                    #     Se_per_study_list[[t - 1]][s, k] <-        ( sum(index_test_results[disease_status == 1] >= k - 1) / n_total_d[s])
                    #     Sp_per_study_list[[t - 1]][s, k] <- 1.0 -  ( sum(index_test_results[disease_status == 0] >= k - 1)/ n_total_nd[s])
                    #   
                    # }
            }
            
          }

          
          
          return(list(
            n_total_nd = n_total_nd,
            n_total_d  = n_total_d, 
            x_nd_list = x_nd_list,
            x_d_list  = x_d_list,
            Se_per_study_list = Se_per_study_list,
            Sp_per_study_list = Sp_per_study_list
          ))
  
}








#' create_test_in_study_for_NMA - function to create "indicator_test_in_study" indicator matrix
#' @export
create_test_in_study_for_NMA <- function(   seed,
                                            n_studies, 
                                            n_tests_inc_ref, 
                                            prob_present = 0.75, 
                                            min_index_tests_per_study = 1
) {
  
  set.seed(seed, kind = "Mersenne-Twister")
  
  n_tests <- n_tests_inc_ref
  
  ## At least 1 index tests per study -> at least 2 total tests per study:
  min_tests_per_study <- min_index_tests_per_study + 1
  
  ## Initialize the matrix:
  indicator_index_test_in_study <- matrix(0, nrow = n_studies, ncol = n_tests)
  
  ## Ref test is ALWAYS included: 
  reference_test_prob <- 1.0
  
  ## For each study, decide which tests are present:
  for (s in 1:n_studies) {
    
    ## Give reference test higher probability of being included
    present_probs <- rep(prob_present, n_tests)
    ## Assuming test 1 is ref test:
    present_probs[1] <- reference_test_prob  
    
    ## Randomly determine which tests are present (1 = present, 0 = absent)
    test_presence <- rbinom(n_tests, 1, present_probs)
    
    ## Ensure at least min_tests_per_study tests are present
    while (sum(test_presence) < min_tests_per_study) {
      ## If not enough tests, add one more randomly
      missing_tests <- which(test_presence == 0)
      if (length(missing_tests) > 0) {
        add_test <- sample(missing_tests, 1)
        test_presence[add_test] <- 1
      } else {
        break  # All tests already included
      }
    }
    
    indicator_index_test_in_study[s, ] <- test_presence
    
  }
  
  ## Optional: Ensure each test appears in at least one study:
  for (t in 1:n_tests) {
    
    if (sum(indicator_index_test_in_study[, t]) == 0) {
      # If a test doesn't appear in any study, add it to a random study
      random_study <- sample(1:n_studies, 1)
      indicator_index_test_in_study[random_study, t] <- 1
    }
    
  }
  
  indicator_index_test_in_study <- indicator_index_test_in_study[1:n_studies, 2:n_tests]
  
  ## Count number of (observed / non-missing): 
  n_index_tests_per_study <- rowSums(indicator_index_test_in_study)
  
  return(list( indicator_index_test_in_study = indicator_index_test_in_study,
               n_index_tests_per_study = n_index_tests_per_study))
  
}















#' apply_thr_missingness
#' @keywords internal
#' @export
apply_thr_missingness <- function(  agg_data_cumulative, 
                                    studies_subset_vec,
                                    missing_thr_subset_vec,
                                    missing_indicator = -1) { 
  
  agg_data_cumulative$x_diseased[studies_subset_vec, missing_thr_subset_vec] <- missing_indicator
  agg_data_cumulative$x_non_diseased[studies_subset_vec, missing_thr_subset_vec] <- missing_indicator
  
  return(agg_data_cumulative)
  
}
















#' convert_cumulative_to_category
#' @keywords internal
#' @export
convert_cumulative_to_category <- function(cumulative_matrix,
                                           missing_indicator = -1) {
  
  # Get dimensions and create output matrix with same dimensions plus one column
  n_rows <- nrow(cumulative_matrix)
  n_cols <- ncol(cumulative_matrix)
  category_matrix <- matrix(-1, nrow = n_rows, ncol = n_cols + 1)
  
  for (i in 1:n_rows) {
    # Get current row
    row_data <- cumulative_matrix[i, ]
    
    # First category is the first cumulative count
    category_matrix[i, 1] <- row_data[1]
    
    # Second category is second column minus first column
    if (n_cols >= 2) {
      if (!is.na(row_data[1]) && !is.na(row_data[2])) {
        category_matrix[i, 2] <- row_data[2] - row_data[1]
      }
    }
    
    # Process remaining categories by taking differences
    for (j in 2:(n_cols-1)) {
      if (!is.na(row_data[j]) && !is.na(row_data[j+1])) {
        category_matrix[i, j+1] <- row_data[j+1] - row_data[j]
      }
    }
    
    # Last category is the total minus the last cumulative count
    # (Assuming we want all categories to sum to the total)
    if (!is.na(row_data[n_cols])) {
      total <- row_data[n_cols]
      category_matrix[i, n_cols+1] <- 0  # Default to zero
      
      # If we want a true non-zero last category, we would need the actual total N
      # which isn't directly provided in this data format
    }
  }
  
  return(category_matrix[, 1:n_cols])
  
}

















#' categorical_to_individual
#' @keywords internal
#' @export
categorical_to_individual <- function( categorical_matrix,
                                       binary_disease_indicator,
                                       missing_indicator = -1) {
  
  # Initialize list to store results
  results <- list()
  n_studies <- nrow(categorical_matrix)
  n_categories <- ncol(categorical_matrix)
  
  for (study in 1:n_studies) {
    # Get counts for current study
    study_counts <- categorical_matrix[study, ]
    
    # Initialize vectors for this study
    study_id <- vector()
    category_values <- vector()
    
    # For each category
    for (cat in 1:n_categories) {
      count <- study_counts[cat]
      
      if (count == missing_indicator) {
        # For missing data, create one observation with missing_indicator
        study_id <- c(study_id, study)
        category_values <- c(category_values, missing_indicator)
      } else if (count > 0) {
        # For non-missing data, replicate the category value count times
        study_id <- c(study_id, rep(study, count))
        category_values <- c(category_values, rep(cat, count))
      }
    }
    
    # Create data frame for this study
    study_data <- data.frame(
      study_id = study_id,
      group = binary_disease_indicator,  # 1 for diseased, 0 for non-diseased
      value = category_values
    )
    
    results[[study]] <- study_data
  }
  
  # Combine all studies into one data frame
  final_data <- do.call(rbind, results)
  rownames(final_data) <- NULL  # Reset row names
  
  return(tibble(final_data))
  
}










#' box_cox_grid
#' @keywords internal
#' @export
box_cox_grid <- function(x, 
                         lambda) {
  
  if (abs(lambda) < 0.001) {
    return(log(x))
  } else {
    return((x^lambda - 1)/lambda)
  }
  
}
 




# individual_obs_tibble <- x_individual
# group_name <- "non-diseased"
# study_index <- 1
# 
# # 
# # 
# # Check what study_id values exist
# unique_studies <- unique(individual_obs_tibble$study_id)
# print("Available study IDs:")
# print(unique_studies)
# 
# # Check what group values exist
# unique_groups <- unique(individual_obs_tibble$group)
# print("Available group values:")
# print(unique_groups)



## Function for log-logistic density:
#' dloglogistic
#' @keywords internal
#' @export
dloglogistic <- function(x, 
                         alpha, 
                         beta) {
  return((alpha/beta) * (x/beta)^(alpha-1) / (1 + (x/beta)^alpha)^2)
}





## Function for Box-Cox density
#' box_cox_density
#' @keywords internal
#' @export
box_cox_density <- function(x, 
                            lambda,
                            data) {
  # Only process positive values (Box-Cox requires positive values)
  if(any(x <= 0)) return(rep(0, length(x)))
  
  # Transform the data
  transformed_data <- box_cox_transform(data, lambda)
  
  # Fit normal to transformed data
  mean_transformed <- mean(transformed_data)
  sd_transformed <- sd(transformed_data)
  
  # Calculate density with Jacobian adjustment
  density_values <- dnorm(box_cox_transform(x, lambda), 
                          mean = mean_transformed, 
                          sd = sd_transformed) * x^(lambda-1)
  
  return(density_values)
}



 
  ## 
  ## lambda <- -0.5:
  ##
  box_cox_neg05 <- function(x) {
    return((x^(-0.5) - 1.0)/(-0.5))
  }
  
  bc_neg05_jacobian <- function(x) { 
    return(x^(-1.5))
  }
  
  bc_neg05_density <- function(x) {
    
    unajusted_dens <- dnorm(box_cox_neg05(x), 
                            mean = mean(x), 
                            sd   = sd(x))
    
    return(unajusted_dens * bc_05_jacobian(x))
    
  }
  
  ## 
  ## lambda <- +0.5:
  ##
  box_cox_05 <- function(x) {
    return((x^0.5 - 1.0)/0.5)
  }
  bc_05_jacobian <- function(x) { 
    return(x^(-0.5))
  }
  bc_05_density <- function(x) {
    
    unajusted_dens <- dnorm(box_cox_05(x), 
                            mean = mean(x), 
                            sd   = sd(x))
    
    return(unajusted_dens * bc_05_jacobian(x))
    
  }
  ## 
  ## lambda > 0:
  ##
  box_cox_pos <- function(x, 
                          lambda) {
    return((x^lambda - 1.0)/lambda)
  }
  bc_pos_jacobian <- function(x,
                              lambda) { 
    return(x^(1.0 - lambda))
  }
  bc_pos_density <- function(x,
                             lambda) {
    
    trans_x <- box_cox_pos(x = x, 
                           lambda = lambda)
    
    unajusted_dens <- dnorm(trans_x, 
                            mean = mean(x), 
                            sd   = sd(x))
    
    Jacobian <- bc_pos_jacobian(x = x,
                                lambda = lambda)
    
    return(unajusted_dens * Jacobian)
    
  }
  
 


#' box_cox_transform
#' @keywords internal
#' @export
box_cox_transform <- function(x, 
                              lambda) {
  
  if(abs(lambda) < 0.001) return(log(x))
  return((x^lambda - 1)/lambda)
  
}




#' bc_density
#' @keywords internal
#' @export
bc_density <- function(x, 
                       lambda, 
                       observations) {
  
  ## Transform original data:
  transformed_data <- box_cox_transform(observations, 
                                        lambda)
  transformed_mean <- mean(transformed_data)
  transformed_sd   <- sd(transformed_data)
  
  ## Apply Jacobian adjustment:
  jacobian <- x^(lambda-1)
  
  ## Calculate density:
  dnorm(box_cox_transform(x, lambda), 
        mean = transformed_mean, 
        sd = transformed_sd) * jacobian
  
}






# study_indexes <- c(1:n_studies)

##
## Function to fit distributions and create plots:
##
#' plot_distribution_fits
#' @keywords internal
#' @export
plot_distribution_fits <- function(individual_obs_tibble, 
                                   n_studies,
                                   study_indexes = NULL, 
                                   group_name = "non-diseased") {
  
  if (group_name == "non-diseased") { 
    group_value <- 0.0
  } else { 
    group_value <- 1.0
  }
  
  if (is.null(study_indexes)) { 
    study_indexes <- seq(from = 1, to = n_studies)
  }
  
  individual_obs_tibble
  df_filtered <- dplyr::filter(individual_obs_tibble, 
                               study_id %in% study_indexes, 
                               group == as.numeric(group_value)) %>% 
    print(n = 100)
  
  ## Extract test values
  observations <- df_filtered$value
  
  ## Calculate bin width based on data range
  bin_width <- max(1, (max(observations) - min(observations)) / 30)
  
  ## Log-normal (need to handle zeros if present)
  obs_for_log <- observations
  if(any(obs_for_log <= 0)) obs_for_log <- observations + 1 # Shift if needed
  
  log_normal_fit <- try(MASS::fitdistr(obs_for_log, "lognormal"))
  if(!inherits(lognormal_fit, "try-error")) {
    lognormal_meanlog <- lognormal_fit$estimate["meanlog"]
    lognormal_sdlog <- lognormal_fit$estimate["sdlog"]
  } else {
    lognormal_meanlog <- mean(log(obs_for_log))
    lognormal_sdlog <- sd(log(obs_for_log))
  }
  
  ## Log-logistic (need to handle zeros if present)
  # Fit log-logistic parameters
  # For log-logistic, we'll estimate shape and scale parameters
  # using method of moments from log-transformed data
  log_data <- log(obs_for_log)
  mu_log   <- mean(log_data)
  sd_log   <- sd(log_data)
  
  ## Convert to log-logistic parameters (approximation)
  ## Shape parameter (alpha) related to the variability
  alpha <- pi/(sd_log * sqrt(3))
  ## Scale parameter (beta) related to the median
  beta <- exp(mu_log)
  # ##
  # ## Create box-cox densities:
  # ##
  # bc_neg1_density   <- bc_density(lambda = -1.0, x = observations)
  # bc_neg05_density  <- bc_density(lambda = -0.5, x = observations)
  # ## For lambda > 0:
  # bc_05_density  <- bc_density(lambda = 0.5, x = observations)
  # bc_1_density   <- bc_density(lambda = 1.0, x = observations)
  # bc_15_density  <- bc_density(lambda = 1.5, x = observations)
  # bc_2_density   <- bc_density(lambda = 2.0, x = observations)
  # # bc_25_density  <- bc_density(lambda = 2.5, x = observations)
  # bc_3_density   <- bc_density(lambda = 3.0, x = observations)
  ##
  ## Plot title:
  ##
  if (length(study_indexes == n_studies)) { 
    plot_title <- paste("Distribution of observed", group_name, "data for ALL", n_studies, "studies")
  } else { 
    plot_title <- paste("Distribution of observed", group_name, "data for all studies", study_index)
  }
  ##
  ## Create plot:
  ##
  plot <- ggplot(data.frame(value = observations), aes(x = value)) +
    ##
    ## Histogram of OBSERVED data:
    ##
    geom_histogram(aes(y = after_stat(density)), binwidth = bin_width, 
                   fill = "lightblue", color = "darkblue", alpha = 0.8) +
    ##
    ## Add density estimate:
    ##
    geom_density(color = "black", size = 2.0) +
    ##
    ## Add log-normal density (if data is positive):
    ##
    stat_function(fun = function(x) dlnorm(x,
                                           meanlog = lognormal_meanlog, 
                                           sdlog   = lognormal_sdlog),
                  color = "green", size = 1.0, linetype = "dashed") + 
    ##
    ## Add log-logistic density:
    ##
    stat_function(fun = function(x) dloglogistic(x, 
                                                 alpha, 
                                                 beta),
                  color = "blue", size = 1.0, linetype = "dashed") + 
    ##
    ## Add Box-Cox densities:
    ##
    stat_function(fun = function(x) bc_density(x, lambda = -1.0, observations = observations), color = "darkblue", size = 1.5) +
    stat_function(fun = function(x) bc_density(x, lambda = -0.5, observations = observations), color = "darkred",  size = 1.5) +
    ##
    stat_function(fun = function(x) bc_density(x, lambda = +0.5,  observations = observations), color = "orange", size = 1.5) +
    stat_function(fun = function(x) bc_density(x, lambda = +1.0,  observations = observations), color = "green",  size = 1.5) +
    stat_function(fun = function(x) bc_density(x, lambda = +1.5,  observations = observations), color = "purple", size = 1.5) +
    stat_function(fun = function(x) bc_density(x, lambda = +2.0,  observations = observations), color = "red",    size = 1.5) +
    ## 
    ## Formatting:
    ##
    labs(title = plot_title,
         subtitle = paste("n =", length(observations)),
         x = "Test Value", 
         y = "Density") +
    theme_minimal() +
    ##
    ## Add legends:
    ##
    annotate("text", x = max(observations)*0.8, y = max(density(observations)$y)*0.90, label = "Log-normal (dashed line)",   color = "green") +
    annotate("text", x = max(observations)*0.8, y = max(density(observations)$y)*0.85, label = "Log-logistic (dashed line)", color = "blue") + 
    ## Box-cox:
    annotate("text", x = max(observations)*0.8, y = max(density(observations)$y)*0.75, label = "Box-Cox (λ=-1.0)", color = "darkblue") +
    annotate("text", x = max(observations)*0.8, y = max(density(observations)$y)*0.70, label = "Box-Cox (λ=-0.5)", color = "darkred") + 
    ## Box-cox:
    annotate("text", x = max(observations)*0.8, y = max(density(observations)$y)*0.65, label = "Box-Cox (λ=+0.5)", color = "orange") + 
    annotate("text", x = max(observations)*0.8, y = max(density(observations)$y)*0.60, label = "Box-Cox (λ=+1.0)", color = "green") + 
    annotate("text", x = max(observations)*0.8, y = max(density(observations)$y)*0.55, label = "Box-Cox (λ=+1.5)", color = "purple") + 
    annotate("text", x = max(observations)*0.8, y = max(density(observations)$y)*0.50, label = "Box-Cox (λ=+2.0)", color = "red")
  
  print(plot)
  
  return(plot)
  
}










# identify_important_comparisons <- function(tibble_comp, 
#                                            n_per_type = 10) {
#   
#         # Define what makes a comparison "important"
#         important_comps <- tibble_comp %>%
#           mutate(
#             # 1. Statistical significance
#             is_significant_se = (diff_Se_lower > 0 | diff_Se_upper < 0),
#             is_significant_sp = (diff_Sp_lower > 0 | diff_Sp_upper < 0),
#             
#             # 2. Effect size (standardized)
#             effect_size_se = abs(diff_Se_median) / diff_Se_sd,
#             effect_size_sp = abs(diff_Sp_median) / diff_Sp_sd,
#             
#             # 3. Clinical relevance score
#             clinical_score = case_when(
#               # Same test, adjacent thresholds
#               test1 == test2 & abs(threshold1 - threshold2) == 1 ~ 3,
#               # Same test, 2 thresholds apart
#               test1 == test2 & abs(threshold1 - threshold2) == 2 ~ 2,
#               # Different tests, similar purpose
#               test1 != test2 ~ 1,
#               TRUE ~ 0
#             ),
#             
#             # Combined importance score
#             importance_score = clinical_score * pmax(effect_size_se, effect_size_sp)
#           )
#         
#         # Get top comparisons by type
#         within_test <- important_comps %>%
#           filter(test1 == test2) %>%
#           arrange(desc(importance_score)) %>%
#           slice_head(n = n_per_type)
#         
#         between_test <- important_comps %>%
#           filter(test1 != test2) %>%
#           arrange(desc(importance_score)) %>%
#           slice_head(n = n_per_type)
#         
#         bind_rows(
#           within_test %>% mutate(comparison_type = "Within-test"),
#           between_test %>% mutate(comparison_type = "Between-test")
#         )
#   
# }
# 
# 
# 


#' filter_tibble_by_thresholds
#' @export
filter_tibble_by_thresholds <- function(tibble_data, 
                                        relevant_thresholds = list(
                                          "GAD-2" = c(1:7), 
                                          "GAD-7" = c(3:18), 
                                          "HADS" = c(3:18),
                                          "BAI" = c(3:40)
                                        )) {
  
        # Create a dataframe of valid test-threshold combinations
        valid_combos <- bind_rows(
          lapply(names(relevant_thresholds), function(test) {
            data.frame(
              test_name = test,
              threshold = relevant_thresholds[[test]],
              stringsAsFactors = FALSE
            )
          })
        )
        
        # Filter using inner join
        filtered_tibble <- tibble_data %>%
          inner_join(valid_combos, by = c("test_name", "threshold"))
        
        return(filtered_tibble)
        
}

 
 

#' filter_clinically_relevant_comparisons
#' @export
filter_clinically_relevant_comparisons <- function(tibble_comp, 
                                                   relevant_thresholds = list(
                                                     "GAD-2" = c(2, 3),
                                                     "GAD-7" = c(5, 10, 15),
                                                     "HADS" = c(8, 11),
                                                     "BAI" = c(8, 16, 26)
                                                   )
) {
  
        # Create valid combinations
        valid_combos <- bind_rows(
          lapply(names(relevant_thresholds), function(test) {
            data.frame(
              test_name = test,
              threshold = relevant_thresholds[[test]]
            )
          })
        )
        
        # Filter using joins
        tibble_comp %>%
          inner_join(valid_combos, by = c("test1_name" = "test_name", "threshold1" = "threshold")) %>%
          inner_join(valid_combos, by = c("test2_name" = "test_name", "threshold2" = "threshold")) %>%
          mutate(
            same_test = (test1 == test2),
            comparison_type = ifelse(same_test, "Within-test", "Between-test"),
            comparison_label = paste0(
              test1_name, " (≥", threshold1, ") vs ",
              test2_name, " (≥", threshold2, ")"
            )
          )
  
}










#' extract_param_latex
#' @export
extract_param_latex <- function(model_or_tibble, 
                                param_name, 
                                digits = 3) {
  
        # Check what we're dealing with
        if(is.data.frame(model_or_tibble)) {
          # It's already a tibble
          tibble <- model_or_tibble
        } else if("model_summary_and_trace_obj" %in% names(model_or_tibble)) {
          # It's your model structure!
          param_parts <- strsplit(param_name, "\\[")[[1]][1]
          tibble <- model_or_tibble$model_summary_and_trace_obj$extract_params(params = param_parts)
        } else if("extract_params" %in% names(model_or_tibble)) {
          # It's the R6 object directly
          param_parts <- strsplit(param_name, "\\[")[[1]][1]
          tibble <- model_or_tibble$extract_params(params = param_parts)
        } else {
          return("$-0.aaa ~ \\ci{0.aaa}{0.aaa}$")
        }
        
        # Filter using base R to avoid dplyr issues
        row <- tibble[tibble$parameter == param_name, ]
        
        if(nrow(row) == 0 || is.na(row$`50%`)) {
          return("$-0.aaa ~ \\ci{0.aaa}{0.aaa}$")
        }
        
        sprintf("$%.3f ~ \\ci{%.3f}{%.3f}$", 
                round(row$`50%`, digits), 
                round(row$`2.5%`, digits), 
                round(row$`97.5%`, digits))
        
      }
      
      # Function 2: Extract from multiple models and create table row
      extract_param_row <- function(model_list, 
                                    param_name, 
                                    digits = 3) {
        
        values <- sapply(model_list, function(model) {
          if(is.null(model)) {
            return("$-0.aaa ~ \\ci{0.aaa}{0.aaa}$")
          }
          extract_param_latex(model, param_name, digits)
        })
        
        paste(values, collapse = " & ")
  
}







#' Updated function to extract Se or Sp values WITH CIs
#' extract_se_sp_latex
#' @export
extract_se_sp_latex <- function(model_or_tibble, 
                                param_type = "Se_baseline", 
                                n_tests = 4,
                                max_n_thr = 63, 
                                digits = 3) {
  
        
        # Get the tibble
        if(is.data.frame(model_or_tibble)) {
          tibble <- model_or_tibble
        } else if("model_summary_and_trace_obj" %in% names(model_or_tibble)) {
          tibble <- model_or_tibble$model_summary_and_trace_obj$extract_params(params = param_type)
        } else {
          return(NULL)
        }
        
        # Create empty arrays for estimates and CIs
        est_matrix <- matrix(NA, nrow = n_tests, ncol = max_n_thr)
        lower_matrix <- matrix(NA, nrow = n_tests, ncol = max_n_thr)
        upper_matrix <- matrix(NA, nrow = n_tests, ncol = max_n_thr)
        
        # Fill the matrices
        for(i in 1:nrow(tibble)) {
          param <- tibble$parameter[i]
          # Extract indices from "Se_baseline[1,2]" format
          indices <- as.numeric(gsub(".*\\[([0-9]+),([0-9]+)\\]", "\\1,\\2", param) %>% 
                                  strsplit(",") %>% unlist())
          
          if(length(indices) == 2 && !is.na(tibble$`50%`[i]) && tibble$`50%`[i] != -1) {
            est_matrix[indices[1], indices[2]] <- tibble$`50%`[i]
            lower_matrix[indices[1], indices[2]] <- tibble$`2.5%`[i]
            upper_matrix[indices[1], indices[2]] <- tibble$`97.5%`[i]
          }
        }
        
        return(list(est = est_matrix, 
                    lower = lower_matrix,
                    upper = upper_matrix))
  
}






 
#' Updated function to format Se/Sp rows with CIs
#' format_se_sp_rows
#' @export
format_se_sp_rows <- function(model_list, 
                              test_num, 
                              thresholds, 
                              n_tests,
                              max_n_thr,
                              threshold_labels,
                              param_type = "Se_baseline",
                              test_name = "Test",
                              digits = 3) {
  
        lines <- c()
        
        # Header
        param_label <- ifelse(grepl("Se", param_type), "Sensitivity", "Specificity")
        lines <- c(lines, sprintf("\\multicolumn{4}{l}{\\textit{\\textbf{%s (%s):}}} \\\\", 
                                  param_label, test_name))
        
        # Extract matrices for all models
        results_list <- lapply(model_list, function(model) {
          extract_se_sp_latex(model, 
                              n_tests = n_tests,
                              max_n_thr = max_n_thr,
                              param_type = param_type,
                              digits = digits)
        })
        
        # For each threshold
        for(i in seq_along(thresholds)) {
          thr <- thresholds[i]
          label <- threshold_labels[i]
          
          # Get values from each model WITH CIs
          values <- sapply(results_list, function(res) {
            if(is.null(res) || is.na(res$est[test_num, thr])) {
              return("$-0.aaa ~ \\ci{0.aaa}{0.aaa}$")
            } else {
              return(sprintf("$%.3f ~ \\ci{%.3f}{%.3f}$", 
                             res$est[test_num, thr],
                             res$lower[test_num, thr],
                             res$upper[test_num, thr]))
            }
          })
          
          # Format the row
          line <- sprintf("$%s_{%d, ~ \\ge %d}$ %s & %s \\\\", 
                          ifelse(grepl("Se", param_type), "Se", "Sp"),
                          test_num, thr, label,
                          paste(values, collapse = " & "))
          
          lines <- c(lines, line)
        }
        
        return(lines)
  
}

 





#' Full function to do all tests at once:
#' create_all_se_sp_sections
#' @export
create_all_se_sp_sections <- function(model_list,
                                      digits,
                                      n_tests,
                                      max_n_thr
) {
  
        all_lines <- c()
        
        # Test configurations
        test_configs <- list(
          list(name = "GAD-2", num = 1, thresholds = c(2, 3, 4), # 3 is SCREENING standard for GAD
               labels = c(" ",
                          "(standard screening thr.)", 
                          " ")),
          list(name = "GAD-7", num = 2, thresholds = c(8, 9, 10, 11), # 10 is SCREENING standard for GAD (10 = "moderate", 5 = "mild", 15 = "severe")
               labels = c(" ",
                          " ",
                          "(moderate; standard screening thr)",
                          " ")),
          list(name = "HADS", num = 3, thresholds = c(8, 9, 10, 11, 12), ## 8 = "Possible", 11 = "probable"   # 11 is standard screening thr.
               labels = c("(possible; screening thr.)", 
                          " ", 
                          " ",
                          "(probable thr.)", 
                          " ")),
          list(name = "BAI", num = 4, thresholds = c(16, 17, 18, 19, 20), # 16 is standard
               labels = c("(moderate thr.)",
                          " ",
                          " ", 
                          " ", 
                          " "))
        )
        
        for(config in test_configs) {
          # Sensitivity
          se_lines <- format_se_sp_rows(  model_list = model_list, 
                                          test_num = config$num, 
                                          thresholds = config$thresholds,
                                          n_tests = n_tests,
                                          max_n_thr = max_n_thr,
                                          threshold_labels = config$labels, 
                                          param_type = "Se_baseline", 
                                          test_name = config$name,
                                          digits = digits)
          all_lines <- c(all_lines, se_lines)
          ##
          # Specificity
          sp_lines <- format_se_sp_rows(  model_list = model_list, 
                                          test_num = config$num, 
                                          thresholds = config$thresholds,
                                          n_tests = n_tests,
                                          max_n_thr = max_n_thr,
                                          threshold_labels = config$labels, 
                                          param_type = "Sp_baseline", 
                                          test_name = config$name,
                                          digits = digits)
          all_lines <- c(all_lines, sp_lines)
          
          # Add midrule between tests (except last)
          if(config$name != "BAI") {
            all_lines <- c(all_lines, "%%%%", "\\midrule", "%%%%")
          }
        }
        
        return(all_lines)
  
}





  
# R_fn_extract_NMA_AUC_summary
# function (tibble_gq, n_thr, test_names = NULL) 
# {
#   require(dplyr)
#   require(tidyr)
#   auc_summary <- tibble_gq %>% filter(grepl("^AUC\\[", parameter)) %>% 
#     mutate(test = as.numeric(gsub(".*\\[(\\d+)\\]", "\\1", 
#                                   parameter))) %>% mutate(test_name = if (!is.null(test_names)) 
#                                     test_names[test]
#                                     else paste0("Test", test)) %>% select(test, test_name, AUC_mean = mean, 
#                                                                           AUC_median = `50%`, AUC_lower = `2.5%`, AUC_upper = `97.5%`, 
#                                                                           AUC_sd = sd)
#   auc_diff_summary <- tibble_gq %>% filter(grepl("^AUC_diff\\[", 
#                                                  parameter)) %>% mutate(index = as.numeric(gsub(".*\\[(\\d+)\\]", 
#                                                                                                 "\\1", parameter))) %>% left_join(expand.grid(test1 = 1:(length(n_thr) - 
#                                                                                                                                                            1), test2 = 1:length(n_thr)) %>% filter(test2 > test1) %>% 
#                                                                                                                                     mutate(index = row_number(), test1_name = if (!is.null(test_names)) 
#                                                                                                                                       test_names[test1]
#                                                                                                                                       else paste0("Test", test1), test2_name = if (!is.null(test_names)) 
#                                                                                                                                         test_names[test2]
#                                                                                                                                       else paste0("Test", test2)), by = "index") %>% select(test1, 
#                                                                                                                                                                                             test2, test1_name, test2_name, AUC_diff_mean = mean, 
#                                                                                                                                                                                             AUC_diff_median = `50%`, AUC_diff_lower = `2.5%`, AUC_diff_upper = `97.5%`, 
#                                                                                                                                                                                             AUC_diff_sd = sd) %>% mutate(prob_test1_better = NA_real_)
#   return(list(auc = auc_summary, auc_diff = auc_diff_summary))
# }
# <bytecode: 0x560bc37856e8>
#   <environment: namespace:MetaOrdDTA>
#   
#   
  
 
#' Function to handle N/A values based on model type
#' get_param_or_na
#' @export
get_param_or_na <- function(model, 
                            param_name,
                            digits = 3) {
        
        # Determine model type from its position
        model_type <- which(sapply(model_list, identical, model))
        
        # Define which params are N/A for which models
        na_params <- list(
          # Model A (fixed cutpoints, CS)
          "A" = c("beta_sigma", "beta_tau\\[\\d+,", "beta_corr", "rho12\\[1,2\\]", "rho12\\[1,3\\]", "rho12\\[1,4\\]"),
          # Model B (fixed cutpoints, UN)
          "B" = c("tau^\\(d-\\).*all tests", "tau^\\(d\\+\\).*all tests"),
          # Model C (random cutpoints, CS)
          "C" = c("beta_sigma", "beta_tau\\[\\d+,", "beta_corr", "rho12\\[1,2\\]", "rho12\\[1,3\\]", "rho12\\[1,4\\]"),
          # Model D (random cutpoints, UN)
          "D" = c("tau^\\(d-\\).*all tests", "tau^\\(d\\+\\).*all tests")
        )
        
        model_letter <- c("A", "B", "C", "D")[model_type]
        
        # Check if this param should be N/A for this model
        for(pattern in na_params[[model_letter]]) {
          if(grepl(pattern, param_name)) {
            return("N/A")
          }
        }
        
        # Otherwise extract the parameter
        return(extract_param_latex(model, 
                                   param_name,
                                   digits))
  
}








 
#' Function to create table 2
#' create_table2_latex
#' @export
create_table2_latex <- function(model_list, 
                                digits = 3) {
  
        lines <- c()
        
        # Test names
        test_names <- c("GAD-2", "GAD-7", "HADS", "BAI")
        
        # Location parameters (non-diseased)
        lines <- c(lines, "\\multicolumn{4}{l}{\\textit{\\textbf{Location parameters (intercept; non-diseased):}}}  \\\\")
        
        for(i in 1:4) {
          param <- paste0("beta_mu[", i, ",1,1]")
          
          values <- sapply(model_list, function(model) {
            extract_param_latex(model, param, digits)
          })
          
          line <- sprintf("${\\beta_{\\bullet, %d, 0}}^{[d-]}$ \\textit{ (%s) } & %s \\\\", 
                          i, test_names[i], paste(values, collapse = " & "))
          lines <- c(lines, line)
        }
        
        # Location parameters (diseased)
        lines <- c(lines, "\\multicolumn{4}{l}{\\textit{\\textbf{Location parameters (intercept; diseased):}}}  \\\\")
        
        for(i in 1:4) {
          param <- paste0("beta_mu[", i, ",2,1]")
          
          values <- sapply(model_list, function(model) {
            extract_param_latex(model, param, digits)
          })
          
          line <- sprintf("${\\beta_{\\bullet, %d, 0}}^{[d+]}$ \\textit{ (%s) } & %s \\\\", 
                          i, test_names[i], paste(values, collapse = " & "))
          lines <- c(lines, line)
        }
        
        # Midrule
        lines <- c(lines, "%%%%", "\\midrule", "%%%%")
        
        # NMA shared between-study SD's - NO N/A VALUES
        lines <- c(lines, "\\multicolumn{5}{l}{\\textit{\\textbf{NMA shared between-study SD's:}}} \\\\")
        
        # beta_sigma[1] for non-diseased
        values <- sapply(model_list, function(model) {
          extract_param_latex(model, "beta_sigma[1]", digits)
        })
        lines <- c(lines, sprintf("$\\sigma_{\\beta}^{(d-)}$ & %s \\\\", paste(values, collapse = " & ")))
        
        # beta_sigma[2] for diseased  
        values <- sapply(model_list, function(model) {
          extract_param_latex(model, "beta_sigma[2]", digits)
        })
        lines <- c(lines, sprintf("$\\sigma_{\\beta}^{(d+)}$ & %s \\\\", paste(values, collapse = " & ")))
        
        # Midrule
        lines <- c(lines, "%%%%", "\\midrule", "%%%%")
        
        # NMA test-specific deviation SD's
        lines <- c(lines, "\\multicolumn{5}{l}{\\textit{\\textbf{NMA test-specific deviation SD's:}}} \\\\")
        lines <- c(lines, "%%%%")
        
        # Non-diseased
        lines <- c(lines, "\\multicolumn{5}{l}{\\textit{Non-diseased:}} \\\\")
        
        for(i in 1:4) {
          param <- paste0("beta_tau[", i, ",1]")
          
          values <- sapply(model_list, function(model) {
            extract_param_latex(model, param, digits)
          })
          
          line <- sprintf("$\\tau_{%d}^{(d-)}$ (%s) & %s \\\\", i, test_names[i], paste(values, collapse = " & "))
          lines <- c(lines, line)
        }
        
        # Diseased
        lines <- c(lines, "\\multicolumn{5}{l}{\\textit{Diseased:}} \\\\")
        
        for(i in 1:4) {
          param <- paste0("beta_tau[", i, ",2]")
          
          values <- sapply(model_list, function(model) {
            extract_param_latex(model, param, digits)
          })
          
          line <- sprintf("$\\tau_{%d}^{(d+)}$ (%s) & %s \\\\", i, test_names[i], paste(values, collapse = " & "))
          lines <- c(lines, line)
        }
        
        # Midrule
        lines <- c(lines, "%%%%", "%%%%", "\\midrule", "%%%%")
        
        # Total heterogeneity
        lines <- c(lines, "\\multicolumn{5}{l}{\\textit{\\textbf{Total within-test, between-study heterogeneity (includes cutpoint variation for RC models):}}} \\\\")
        lines <- c(lines, "\\multicolumn{5}{l}{\\textit{Non-diseased:}} \\\\")
        
        for(i in 1:4) {
          param <- paste0("total_SD_inc_C[", i, ",1]")
          
          values <- sapply(model_list, function(model) {
            extract_param_latex(model, param, digits)
          })
          
          line <- sprintf("$\\sigma_{total,%d}^{(d-)}$ (%s) & %s \\\\", i, test_names[i], paste(values, collapse = " & "))
          lines <- c(lines, line)
        }
        
        lines <- c(lines, "\\multicolumn{5}{l}{\\textit{Diseased:}} \\\\")
        
        for(i in 1:4) {
          param <- paste0("total_SD_inc_C[", i, ",2]")
          
          values <- sapply(model_list, function(model) {
            extract_param_latex(model, param, digits)
          })
          
          line <- sprintf("$\\sigma_{total,%d}^{(d+)}$ (%s) & %s \\\\", i, test_names[i], paste(values, collapse = " & "))
          lines <- c(lines, line)
        }
        
        
        
        # Midrule
        lines <- c(lines, "%%%%", "\\midrule", "%%%%")
        
        # Correlation structure
        lines <- c(lines, "\\multicolumn{5}{l}{\\textit{\\textbf{Correlation structure:}}} \\\\")
        
        # beta_corr  
        values <- sapply(model_list, function(model) {
          extract_param_latex(model, "beta_corr", digits)
        })
        lines <- c(lines, sprintf("$\\rho_{\\beta}$ (NMA shared corr.) & %s \\\\", paste(values, collapse = " & ")))
        
        # rho12
        lines <- c(lines, "\\multicolumn{5}{l}{\\textit{Between-study, within-test correlations ($\\rho_{12}$):}} \\\\")
        
        for(i in 1:4) {
          param <- paste0("rho12[", i, ",", i, "]")
          
          values <- sapply(model_list, function(model) {
            extract_param_latex(model, param, digits)
          })
          
          line <- sprintf("%s & %s \\\\", test_names[i], paste(values, collapse = " & "))
          lines <- c(lines, line)
        }
        
        return(lines)
  
}



create_AUC_latex_table <- function( auc_outs,
                                    digits = 3) {
  
        
        auc_main <- auc_outs$auc
        auc_diff <- auc_outs$auc_diff
        auc_pred <- auc_outs$auc_pred
        
        lines <- c()
        
        # Header
        lines <- c(lines, 
                   "\\begin{table}[!h]",
                   "\\centering", 
                   "\\caption{AUC values, pairwise differences, and predictive intervals}",
                   "\\small",
                   "\\begin{tabular}{lcc}",
                   "\\toprule")
        
        # Section 1: Main AUC values
        lines <- c(lines, 
                   "\\multicolumn{3}{l}{\\textbf{AUC estimates (95\\% CI):}} \\\\",
                   "\\midrule")
        
        for(i in 1:nrow(auc_main)) {
          line <- sprintf("%s & %.3f & (%.3f, %.3f) \\\\",
                          auc_main$test_name[i],
                          auc_main$AUC_median[i],
                          auc_main$AUC_lower[i],
                          auc_main$AUC_upper[i])
          lines <- c(lines, line)
        }
        
        # Section 2: Pairwise differences
        lines <- c(lines,
                   "\\midrule",
                   "\\multicolumn{3}{l}{\\textbf{Pairwise AUC differences:}} \\\\",
                   "\\midrule")
        
        for(i in 1:nrow(auc_diff)) {
          sig_marker <- ifelse(auc_diff$AUC_diff_lower[i] > 0 | auc_diff$AUC_diff_upper[i] < 0, "*", "")
          line <- sprintf("%s vs %s & %.3f & (%.3f, %.3f)%s \\\\",
                          auc_diff$test1_name[i],
                          auc_diff$test2_name[i],
                          auc_diff$AUC_diff_median[i],
                          auc_diff$AUC_diff_lower[i],
                          auc_diff$AUC_diff_upper[i],
                          sig_marker)
          lines <- c(lines, line)
        }
        
        # Section 3: Predictive intervals
        lines <- c(lines,
                   "\\midrule",
                   "\\multicolumn{3}{l}{\\textbf{AUC predictive intervals:}} \\\\",
                   "\\midrule")
        
        for(i in 1:nrow(auc_pred)) {
          line <- sprintf("%s & %.3f & (%.3f, %.3f) \\\\",
                          auc_pred$test_name[i],
                          auc_pred$AUC_median[i],
                          auc_pred$AUC_lower[i],
                          auc_pred$AUC_upper[i])
          lines <- c(lines, line)
        }
        
        # Footer
        lines <- c(lines,
                   "\\bottomrule",
                   "\\end{tabular}",
                   "\\begin{tablenotes}",
                   "\\footnotesize",
                   "\\item * indicates 95\\% CI excludes zero (significant difference)",
                   "\\item Predictive intervals incorporate between-study heterogeneity",
                   "\\end{tablenotes}",
                   "\\end{table}")
        
        return(paste(lines, collapse = "\n"))
        
}













# Function 1: Create meta-regression coefficients table
create_metareg_coef_table <- function(beta_tibble, 
                                      test_names = c("GAD-2", "GAD-7", "HADS", "BAI")) {
  
          lines <- c()
          
          # Header
          lines <- c(lines,
                     "\\begin{table}[!h]",
                     "\\centering",
                     "\\caption{Meta-regression coefficients (probit scale)}",
                     "\\small",
                     "\\begin{tabular}{lcc}",
                     "\\toprule",
                     "& \\multicolumn{1}{c}{Non-diseased} & \\multicolumn{1}{c}{Diseased} \\\\",
                     "\\midrule"
          )
          
          # Process each test
          for (test in test_names) {
            # Add test header
            lines <- c(lines, sprintf("\\multicolumn{3}{l}{\\textbf{%s}} \\\\", test))
            
            # Get unique covariates for this test
            test_data <- beta_tibble %>% filter(test_name == test)
            covariates <- unique(test_data$covariate_name)
            
            for (cov in covariates) {
              # Get non-diseased and diseased values
              nd_row <- test_data %>% filter(covariate_name == cov, class_name == "Non-diseased")
              d_row <- test_data %>% filter(covariate_name == cov, class_name == "Diseased")
              
              # Format values
              if (nrow(nd_row) > 0) {
                nd_val <- sprintf("%.3f ~ \\ci{%.3f}{%.3f}", nd_row$mean, nd_row$`2.5%`, nd_row$`97.5%`)
                if (nd_row$is_significant) nd_val <- sprintf("\\mathbf{%s}", nd_val)
              } else {
                nd_val <- "--"
              }
              
              if (nrow(d_row) > 0) {
                d_val <- sprintf("%.3f ~ \\ci{%.3f}{%.3f}", d_row$mean, d_row$`2.5%`, d_row$`97.5%`)
                if (d_row$is_significant) d_val <- sprintf("\\mathbf{%s}", d_val)
              } else {
                d_val <- "--"
              }
              
              # Clean covariate name
              cov_display <- case_when(
                cov == "intercept" ~ "Intercept",
                cov == "logit_prev_GAD" ~ "Prevalence (per SD)",
                cov == "Ref_test_clean_SCID" ~ "Ref test: SCID",
                cov == "study_setting_3" ~ "Study setting 3",
                TRUE ~ cov
              )
              
              # Add row
              lines <- c(lines, sprintf("%-20s & $%s$ & $%s$ \\\\", cov_display, nd_val, d_val))
            }
            
            # Add space between tests
            if (test != test_names[length(test_names)]) {
              lines <- c(lines, "%%%%", "\\midrule", "%%%%")
            }
          }
          
          # Footer
          lines <- c(lines,
                     "%%%%",
                     "\\bottomrule",
                     "%%%%",
                     "\\end{tabular}",
                     "\\begin{tablenotes}",
                     "\\item Bold indicates significant at 95\\% level",
                     "\\end{tablenotes}",
                     "\\end{table}"
          )
          
          return(paste(lines, collapse = "\n"))
  
}





# Function 2: Create heterogeneity comparison table
create_heterogeneity_comparison_table <- function(total_SD_intercept, total_SD_covariate,
                                                  beta_sigma_intercept, beta_sigma_covariate,
                                                  beta_tau_intercept, beta_tau_covariate,
                                                  beta_corr_intercept, beta_corr_covariate,
                                                  rho12_intercept, rho12_covariate) {
  
          lines <- c()
          
          # Helper function to format values
          format_val <- function(row) {
            sprintf("%.3f ~ \\ci{%.3f}{%.3f}", row$`50%`, row$`2.5%`, row$`97.5%`)
          }
          
          # Header
          lines <- c(lines,
                     "\\begin{table}[!h]",
                     "\\centering",
                     "\\caption{Impact of covariates on heterogeneity parameters}",
                     "\\small",
                     "\\begin{tabular}{lcc}",
                     "\\toprule",
                     "\\textbf{Parameter}  & \\textbf{Model A}          & \\textbf{Model A + Covariates} \\\\",
                     "                    & \\textit{(intercept-only)} & \\textit{(prevalence, ref test, setting)} \\\\",
                     "%%%%",
                     "\\midrule",
                     "%%%%"
          )
          
          # NMA shared between-study SD's
          lines <- c(lines,
                     "\\multicolumn{3}{l}{\\textit{\\textbf{NMA shared between-study SD's:}}} \\\\",
                     sprintf("$\\sigma_{\\beta}^{(d-)}$ & $%s$ & $%s$ \\\\",
                             format_val(beta_sigma_intercept[1,]),
                             format_val(beta_sigma_covariate[1,])),
                     sprintf("$\\sigma_{\\beta}^{(d+)}$ & $%s$ & $%s$ \\\\",
                             format_val(beta_sigma_intercept[2,]),
                             format_val(beta_sigma_covariate[2,]))
          )
          
          # Test-specific deviation SD's (CS structure - same for all tests)
          lines <- c(lines,
                     "%%%%",
                     "\\midrule",
                     "%%%%",
                     "\\multicolumn{3}{l}{\\textit{\\textbf{NMA test-specific deviation SD's (CS structure):}}} \\\\",
                     sprintf("$\\tau^{(d-)}$ (all tests) & $%s$ & $%s$ \\\\",
                             format_val(beta_tau_intercept[1,]),  # Just use first row since CS
                             format_val(beta_tau_covariate[1,])),
                     sprintf("$\\tau^{(d+)}$ (all tests) & $%s$ & $%s$ \\\\",
                             format_val(beta_tau_intercept[5,]),  # Row 5 is diseased
                             format_val(beta_tau_covariate[5,]))
          )
          
          # Total heterogeneity (CS structure - same for all tests)
          lines <- c(lines,
                     "%%%%",
                     "\\midrule",
                     "%%%%",
                     "\\multicolumn{3}{l}{\\textit{\\textbf{Total within-test, between-study heterogeneity:}}} \\\\",
                     sprintf("$\\sigma_{total}^{(d-)}$ (all tests) & $%s$ & $%s$ \\\\",
                             format_val(total_SD_intercept[1,]),  # Just use first row since CS
                             format_val(total_SD_covariate[1,])),
                     sprintf("$\\sigma_{total}^{(d+)}$ (all tests) & $%s$ & $%s$ \\\\",
                             format_val(total_SD_intercept[5,]),  # Row 5 is diseased
                             format_val(total_SD_covariate[5,]))
          )
          
          # Correlation structure
          lines <- c(lines,
                     "%%%%",
                     "\\midrule",
                     "%%%%",
                     "\\multicolumn{3}{l}{\\textit{\\textbf{Correlation structure:}}} \\\\",
                     sprintf("$\\rho_{\\beta}$ (NMA shared corr.) & $%s$ & $%s$ \\\\",
                             format_val(beta_corr_intercept),
                             format_val(beta_corr_covariate)),
                     sprintf("$\\rho_{12}$ (within-test corr.) & $%s$ & $%s$ \\\\",
                             format_val(rho12_intercept[1,]),  # Just use first non-NA row
                             format_val(rho12_covariate[1,]))
          )
          
          # Footer
          lines <- c(lines,
                     "%%%%",
                     "\\bottomrule",
                     "%%%%",
                     "\\end{tabular}",
                     "%%%%",
                     "\\begin{tablenotes}",
                     "\\footnotesize",
                     "\\item Both models use fixed cutpoints with compound symmetry correlation (Model A specification).",
                     "\\item CS = compound symmetry; all tests share the same variance components.",
                     "\\end{tablenotes}",
                     "%%%%",
                     "\\end{table}"
          )
          
          return(paste(lines, collapse = "\n"))
  
}

 





