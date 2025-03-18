 

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




#' R_fn_grp_mat_cols_by_first_col_per_superchain_grp
#' @keywords internal
#' @export
R_fn_grp_mat_cols_by_first_col_per_superchain_grp <- function( matrix,
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
          # Take first column of subset and copy it across
          matrix[, start_col:end_col] <- matrix[, rep(start_col, n_chains_per_superchain)]
        }
        
        # Handle remaining columns if any
        if (n_full_superchains * n_chains_per_superchain < n_chains) {
          start_col <- n_full_superchains * n_chains_per_superchain + 1
          subset <- matrix[, start_col:n_chains, drop = FALSE]
          # Take first column of remaining subset and copy it across
          matrix[, start_col:n_chains] <- matrix[, rep(start_col, ncol(subset))]
        }
        
        return(matrix)
  
}




#' R_fn_post_burnin_prep_for_sampling
#' @keywords internal
#' @export
R_fn_post_burnin_prep_for_sampling <- function(n_chains,
                                               n_superchains, 
                                               n_params_raw,
                                               theta_main_vectors_all_chains_input_from_R
                                               )  {
  
          n_chains_per_superchain <- n_chains / n_superchains ; n_chains_per_superchain
          
          if (n_superchains * n_chains_per_superchain > n_chains) { 
            
                warning("n_superchains * n_chains_per_superchain > n_chains - lowering n_chains as necessary")
                
                
                while  (n_superchains * n_chains_per_superchain > n_chains) {
                  
                  n_chains <- n_chains - 1
                
                } 
            
          }
          
          
          theta_main_vectors_all_chains_input_from_R_sampling <- array(dim = c(n_params_raw, n_chains))
          
          counter = 1
          for (kk in 1:n_chains) { 
              theta_main_vectors_all_chains_input_from_R_sampling[,kk] <-  theta_main_vectors_all_chains_input_from_R[,counter]
        
              if (kk %% n_chains_burnin == 0) { 
                counter = counter + 1
              }
          }
          
          # re-assign
          theta_main_vectors_all_chains_input_from_R <- theta_main_vectors_all_chains_input_from_R_sampling
          
          #n_superchains <- 4  # we will have n_superchains sets of initial values (chains within same superchains share same init's)
          theta_main_vectors_all_chains_input_from_R <- R_fn_grp_mat_cols_by_first_col_per_superchain_grp( theta_main_vectors_all_chains_input_from_R,
                                                                                                           n_chains_per_superchain)

          
          return(list(theta_main_vectors_all_chains_input_from_R = theta_main_vectors_all_chains_input_from_R))
  
  
}







