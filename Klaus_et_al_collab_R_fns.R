




compute_Se_Sp_by_cutpoint <- function(data,
                                      min_score,
                                      max_score) {
  
          # Ensure data is a data frame
          data <- as.data.frame(data)
          
          # Get unique cutoff values (assuming they're integers from 0 to 6)
          cutpoints <- min_score:max_score
          
          # Initialize result storage
          result <- list(
            cutpoints = cutpoints,
            Se = numeric(length(cutpoints)),
            Sp = numeric(length(cutpoints)),
            n_studies = numeric(length(cutpoints))
          )
          
          # Calculate mean Se and Sp for each cutpoint
          for (i in seq_along(cutpoints)) {
            cutpoint <- cutpoints[i]
            
            # Filter data for this cutpoint
            subset_data <- data[data$Cutoff == cutpoint, ]
            
            # Count number of studies for this cutpoint
            n_studies <- nrow(subset_data)
            
            if (n_studies > 0) {
              # Calculate mean Se and Sp
              mean_Se <- mean(subset_data$Se, na.rm = TRUE)
              mean_Sp <- mean(subset_data$Sp, na.rm = TRUE)
              
              # Store results
              result$Se[i] <- mean_Se
              result$Sp[i] <- mean_Sp
              result$n_studies[i] <- n_studies
            } else {
              # If no studies for this cutpoint, set to NA
              result$Se[i] <- NA
              result$Sp[i] <- NA
              result$n_studies[i] <- 0
            }
          }
          
          # Convert to a data frame for easier viewing
          result_df <- data.frame(
            Cutpoint = result$cutpoints,
            Se = result$Se,
            Sp = result$Sp,
            n_studies = result$n_studies
          )
          
          # Add formatted values for printing
          result_df$Se_formatted <- sprintf("%.3f", result$Se)
          result_df$Sp_formatted <- sprintf("%.3f", result$Sp)
          
          return(result_df)
  
}








transform_questionnaire_data_by_disease <- function(data, 
                                                    min_score,
                                                    max_score) {
  
          # Get unique studies
          unique_studies <- unique(data$Study)
          n_studies <- length(unique_studies)
          
          # Create empty matrices with rows = n_studies and cols = max_score + 1 thresholds (0 to max_score)
          # One for diseased and one for non-diseased populations
          diseased_matrix <- matrix(-1, nrow = n_studies, ncol = max_score + 1)
          nondiseased_matrix <- matrix(-1, nrow = n_studies, ncol = max_score + 1)
          
          colnames(diseased_matrix) <- colnames(nondiseased_matrix) <- paste0(">=", min_score:max_score)
          rownames(diseased_matrix) <- rownames(nondiseased_matrix) <- unique_studies
          
          # Fill the matrices
          for (i in seq_along(unique_studies)) {
            
                  study_id <- unique_studies[i]
                  study_data <- data[data$Study == study_id, ]
                  
                  # For each study, sort by cutoff
                  study_data <- study_data[order(study_data$Cutoff), ]
                  
                  # Get the diseased and non-diseased Ns for this study
                  n_diseased <- unique(study_data$N_diseased)
                  if(length(n_diseased) > 1) {
                    warning(paste("Study", study_id, "has inconsistent N_diseased values. Using the first one."))
                    n_diseased <- n_diseased[1]
                  }
                  
                  total_n <- unique(study_data$N)
                  if(length(total_n) > 1) {
                    warning(paste("Study", study_id, "has inconsistent N values. Using the first one."))
                    total_n <- total_n[1]
                  }
                  
                  n_nondiseased <- total_n - n_diseased
                  
                  # Fill in the >=0 column with total N for each group (everyone scores at least 0)
                  diseased_matrix[i, 1] <- n_diseased
                  nondiseased_matrix[i, 1] <- n_nondiseased
                  
                  # Fill in the other thresholds
                  for (threshold in 1:max_score) {
                    
                          # Find the rows matching this cutoff
                          cutoff_rows <- which(study_data$Cutoff == threshold)
                          
                          if (length(cutoff_rows) > 0) {
                            # If we have data for this threshold
                            row_to_use <- cutoff_rows[1]  # Use the first one if multiple
                            
                            # For diseased population: true positives (TP)
                            # These are the people with disease scoring at or above the threshold
                            tp <- study_data$TP[row_to_use]
                            if (!is.na(tp)) {
                              diseased_matrix[i, threshold + 1] <- tp
                            }
                            
                            # For non-diseased population: false positives (FP)
                            # These are the people without disease scoring at or above the threshold
                            fp <- study_data$FP[row_to_use]
                            if (!is.na(fp)) {
                              nondiseased_matrix[i, threshold + 1] <- fp
                            }
                          }
                          # If no data for this threshold, leave as -1
                          
                  }
            
          }
          
          # Convert to integers and ensure monotonic decrease
          diseased_matrix <- convert_to_integers_simple(diseased_matrix)
          nondiseased_matrix <- convert_to_integers_simple(nondiseased_matrix)
          
          return(list(diseased = diseased_matrix, 
                      nondiseased = nondiseased_matrix))
  
}









# Helper function to convert to integers and ensure monotonic decrease
convert_to_integers_simple <- function(count_matrix) {
  
          integer_matrix <- count_matrix  # Start with a copy
          
          for (i in 1:nrow(count_matrix)) {
            # First pass: round all non-missing values
            for (j in 1:ncol(count_matrix)) {
              if (count_matrix[i, j] != -1) {
                integer_matrix[i, j] <- round(count_matrix[i, j])
              }
            }
            
            # Second pass: ensure monotonic decrease
            for (j in 2:ncol(count_matrix)) {
              if (integer_matrix[i, j] != -1 && integer_matrix[i, j-1] != -1) {
                if (integer_matrix[i, j] > integer_matrix[i, j-1]) {
                  integer_matrix[i, j] <- integer_matrix[i, j-1]
                }
              }
            }
          }
          
          return(integer_matrix)
  
}







create_cutpoint_index_matrix <- function(data_matrix, 
                                         max_score) {
          
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










impute_missing_values <- function(x_list) {
  
        # Create a copy of the input list to avoid modifying the original
        imputed_x_list <- x_list
        
        # Process each disease status group (non-diseased and diseased)
        for (j in 1:length(x_list)) {
          x_matrix <- imputed_x_list[[j]]
          n_rows <- nrow(x_matrix)
          n_cols <- ncol(x_matrix)
          
          # Process each row (study)
          for (i in 1:n_rows) {
            row_data <- x_matrix[i, ]
            
            # First, ensure the total (first column) is not missing
            if (row_data[1] == -1) {
              # If total is missing, find the first non-missing value
              first_non_missing <- min(which(row_data != -1), na.rm = TRUE)
              if (length(first_non_missing) > 0 && first_non_missing > 1) {
                # Set total to a value slightly higher than the next non-missing value
                x_matrix[i, 1] <- ceiling(row_data[first_non_missing] * 1.2)
              } else {
                # If all values are missing (unlikely), set a default
                x_matrix[i, 1] <- 100
              }
            }
            
            # Now impute missing values while ensuring monotonic decrease
            for (k in 2:n_cols) {
              if (row_data[k] == -1) {
                # Find nearest non-missing values before and after
                prev_indices <- which(row_data[1:(k-1)] != -1)
                next_indices <- which(row_data[(k+1):n_cols] != -1)
                
                if (length(prev_indices) > 0 && length(next_indices) > 0) {
                  # Both before and after values exist
                  prev_val <- row_data[max(prev_indices)]
                  next_val <- row_data[k + min(next_indices)]
                  
                  # Linear interpolation
                  prev_pos <- max(prev_indices)
                  next_pos <- k + min(next_indices)
                  proportion <- (k - prev_pos) / (next_pos - prev_pos)
                  imputed_value <- round(prev_val - proportion * (prev_val - next_val))
                  
                } else if (length(prev_indices) > 0) {
                  # Only previous values exist
                  prev_val <- row_data[max(prev_indices)]
                  
                  # Exponential decay
                  steps_from_prev <- k - max(prev_indices)
                  decay_factor <- 0.6  # Adjust as needed
                  imputed_value <- round(prev_val * (decay_factor^steps_from_prev))
                  
                } else if (length(next_indices) > 0) {
                  # Only next values exist
                  next_val <- row_data[k + min(next_indices)]
                  
                  # Create a value greater than the next value
                  imputed_value <- ceiling(next_val * 1.3)
                  
                } else {
                  # No other values to reference
                  # Use a fraction of the total
                  imputed_value <- round(row_data[1] * (0.5^(k-1)))
                }
                
                # Assign the imputed value
                x_matrix[i, k] <- imputed_value
              }
            }
            
            # Final check to ensure values are monotonically decreasing
            for (k in 2:n_cols) {
              if (x_matrix[i, k] > x_matrix[i, k-1]) {
                x_matrix[i, k] <- x_matrix[i, k-1]
              }
            }
          }
          
          imputed_x_list[[j]] <- x_matrix
        }
        
        return(imputed_x_list)
  
}


























