



#' R_fn_avg_matrix_columns
#' @keywords internal
#' @export
R_fn_avg_matrix_columns <- function(matrix, 
                                    n_chains_per_superchain) {
  
        n_params <- nrow(matrix)
        n_chains <- ncol(matrix)
        
        # Calculate the number of full superchains
        n_full_superchains <- floor(n_chains / n_chains_per_superchain)
        
        # Process full superchains
        for (i in 1:n_full_superchains) { 
          start_col <- (i - 1) * n_chains_per_superchain + 1
          end_col <- i * n_chains_per_superchain
          subset <- matrix[, start_col:end_col, drop = FALSE]
          avg <- rowMeans(subset)
          matrix[, start_col:end_col] <- matrix[, rep(start_col, n_chains_per_superchain)]
        }
        
        # Handle remaining columns if any
        if (n_full_superchains * n_chains_per_superchain < n_chains) {
          start_col <- n_full_superchains * n_chains_per_superchain + 1
          subset <- matrix[, start_col:n_chains, drop = FALSE]
          avg <- rowMeans(subset)
          matrix[, start_col:n_chains] <- matrix[, rep(start_col, ncol(subset))]
        }
        
        return(matrix)
  
}



