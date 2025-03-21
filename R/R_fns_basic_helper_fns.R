




#' simplex_to_unconstrained
#' @export
simplex_to_unconstrained <- function(x) {
      
      N <- length(x)
      y <- numeric(N-1)
      
      # Small constant to avoid exact 0 or 1 values
      epsilon <- 1e-10
      
      # Ensure x is a valid simplex and avoid exact 0s and 1s
      x <- pmax(x, epsilon)
      x <- pmin(x, 1-epsilon)
      x <- x / sum(x)  # Re-normalize to ensure sum is 1
      
      # Calculate the stick-breaking proportions
      remaining_mass <- 1.0
      for (i in 1:(N-1)) {
        # Calculate proportion with numerical safeguards
        z_i <- min(max(x[i] / remaining_mass, epsilon), 1-epsilon)
        
        # Transform to unconstrained space
        y[i] <- log(z_i/(1-z_i)) + log(N-i)
        
        # Update remaining mass
        remaining_mass <- remaining_mass - x[i]
        
        # Avoid numerical issues when remaining mass gets close to 0
        if (remaining_mass < epsilon) {
          # Fill remaining elements with reasonable defaults
          if (i < (N-1)) {
            y[(i+1):(N-1)] <- rnorm(N-1-i, 0, 1)
          }
          break
        }
      }
      
      return(y)
  
}


#' generate_inits_for_raw_simplex_vec
#' @export
generate_inits_for_raw_simplex_vec <- function( n_thr, 
                                                seed, 
                                                width = 0.01) {
  
    set.seed(seed, kind = "L'Ecuyer-CMRG")
   
    n_cat_t <- n_thr + 1
    uniform_simplex <- rep(1/n_cat_t, n_cat_t)  # Start with uniform
    # Add some variation while ensuring it remains a valid simplex
    noise <- runif(n_cat_t, -width, width)
    informative_simplex <- uniform_simplex + noise
    # Make sure it's still a valid simplex
    informative_simplex <- informative_simplex / sum(informative_simplex)
    
    dirichlet_cat_means_phi_raw_mat <- simplex_to_unconstrained(informative_simplex)
    
    return(dirichlet_cat_means_phi_raw_mat)
 
  
  
}



#' R_fn_remove_duplicates_from_list
#' @keywords internal
#' @export
R_fn_remove_duplicates_from_list <- function(stan_data_list) {
  
        ##
        ## First, check for duplicates in the data list:
        ##
        all_names <- names(stan_data_list)
        duplicate_names <- all_names[duplicated(all_names)]
        print(duplicate_names)
        ##
        ## Remove duplicates by keeping only the first occurrence:
        ##
        unique_stan_data_list <- stan_data_list[!duplicated(names(stan_data_list))]
        ##
        ## Output unique list:
        ##
        return(unique_stan_data_list)
  
}






#' check_if_BayesMVP_R_pkg_installed
#' @keywords internal
#' @export
check_if_BayesMVP_R_pkg_installed <- function(if_true = TRUE,
                                              if_false = FALSE) {
          
        # Check if package is among installed packages:
        is_installed <- "BayesMVP" %in% installed.packages()[,"Package"]
        
        # Return the appropriate value based on whether the package is installed
        if(is_installed) {
          return(if_true)
        } else {
          return(if_false)
        }
          
}






## Helper R function to remove "log_lik" parameter from trace array 
#' R_fn_remove_log_lik_from_array
#' @keywords internal
#' @export
R_fn_remove_log_lik_from_array <- function(arr) {
  
          format(object.size(arr), units = "MB")
          
          # Find which parameters don't contain "log_lik"
          keep_params <- !grepl("log_lik", dimnames(arr)[[1]])
          
          # Create new array without the log_lik parameters
          arr_filtered <- arr[, , keep_params, drop = FALSE]
          
          return(arr_filtered)
  
}






#' get_BayesMVP_stan_paths
#' @keywords internal
#' @export
get_MetaOrdDTA_stan_paths <- function() {

        Sys.setenv(stan_THREADS = "true")
        
        ## Get package directory paths:
        pkg_root_dir <- system.file(package = "MetaOrdDTA")
        pkg_data_dir <- file.path(pkg_root_dir, "stan_data")  # directory to store data inc. JSON data files
        pkg_stan_dir <- file.path(pkg_root_dir, "stan_models")
        ##
        print(paste("pkg_root_dir = ", pkg_root_dir))
        print(paste("pkg_data_dir = ", pkg_data_dir))
        print(paste("pkg_stan_dir = ", pkg_stan_dir))
        
        return(list(pkg_root_dir = pkg_root_dir,
                    pkg_data_dir = pkg_data_dir,
                    pkg_stan_dir = pkg_stan_dir))
    
}






#' make_array_dims_equal_Bayesplot_order
#' @keywords internal
#' @export
make_array_dims_equal_Bayesplot_order <- function( input_array, 
                                                   iter_dim,
                                                   chain_dim,
                                                   params_dim) {
  
        # Create the permutation vector for aperm
        # The desired order is: Iteration, Chain, Parameter
        desired_array_dims_vec <- c(iter_dim, chain_dim, params_dim)
        
        # Use base::aperm to reorder dimensions
        new_array <- base::aperm(input_array, desired_array_dims_vec)
        ##
        # print(paste("make_array_dims_equal_Bayesplot_order - new_array = "))
        # print(str(new_array))
        ##
        return(new_array)
        
}




#' format_named_array_for_bayesplot
#' @keywords internal
#' @export
format_named_array_for_bayesplot <- function(named_input_array) { 
  
  
        dims_numeric <- dim(named_input_array)
        named_dims <- dimnames(named_input_array)
        names_of_named_dims <- names(named_dims)
        cat("names_of_named_dims = ", names_of_named_dims)
        ##
        first_dim_name  <- names_of_named_dims[1]
        second_dim_name <- names_of_named_dims[2]
        third_dim_name  <- names_of_named_dims[3]
        ##
        old_iter_dim <- which(stringr::str_detect(names_of_named_dims, "iter") == TRUE)
        print(paste("old_iter_dim = ", old_iter_dim))
        ##
        old_chain_dim <- which(stringr::str_detect(names_of_named_dims, "chain") == TRUE)
        print(paste("old_chain_dim = ", old_chain_dim))
        ##
        condition_temp <- (stringr::str_detect(names_of_named_dims, "var") | stringr::str_detect(names_of_named_dims, "par"))
        old_params_dim <- which(condition_temp == TRUE)
        print(paste("old_params_dim = ", old_params_dim))
        ##
        ## For Bayesplot, you need dims in this order: ⁠Iteration, Chain, Parameter
        ##
        new_array <- make_array_dims_equal_Bayesplot_order( input_array = named_input_array, 
                                                            iter_dim = old_iter_dim,
                                                            chain_dim = old_chain_dim,
                                                            params_dim = old_params_dim)
        ##
        print(paste("format_named_array_for_bayesplot - new_array = "))
        print(str(new_array))
        ##
        return(new_array)
  
}








# param_string <- "beta"

#' find_all_stan_param_names
#' @keywords internal
#' @export
find_all_stan_param_names <- function( stan_draws_array, 
                                       exclude_NA = FALSE,
                                       debugging = FALSE) {
          
          ## Get parameter names that contain the string anywhere:
          all_stan_param_names_inc_NAs <- dimnames(stan_draws_array)[[3]]
          
          if (exclude_NA) {
            
                ## Filter out parameters that are all NA:
                param_names_wo_NAs <- character(0)
                for (param in all_stan_param_names_inc_NAs) {
                  if (!all(is.na(stan_draws_array[,,param]))) {
                    param_names_wo_NAs <- c(param_names_wo_NAs, param)
                  }
                }
                
                if (length(param_names_wo_NAs) == 0) {
                  warning(paste("All parameters containing", param_string, "contain only NAs"))
                  return(list())  # Return empty list if no valid parameters
                }
            
          } else { 
            
                param_names_wo_NAs <- all_stan_param_names_inc_NAs
            
          }
          
          if (debugging) {
            print(param_names_wo_NAs)
          }
          
          return(param_names_wo_NAs)
  
}







 
# format_array_for_bayesplot <- function(input_array) { 
#           
#           input_array_2 <- (base::aperm(input_array,   c(2,1,3))) ; str(input_array_2)
#           input_array_3 <- (base::aperm(input_array_2, c(1,3,2))) ; str(input_array_3)
#           return(input_array_3)
#   
# }
# 
# 




#' convert_list_of_arrays_to_one_big_array
#' @keywords internal
#' @export
convert_list_of_arrays_to_one_big_array <- function(draws_array_list) {
  
  
  str(draws_array_list)
  
          if (is.list(draws_array_list)) {
  
                ##
                ## Verify all arrays have same iterations and chains:
                ##
                n_iter   <- dim(draws_array_list[[1]])[1] ; n_iter
                print(paste("convert_list_of_arrays_to_one_big_array - n_iter = ", n_iter))
                ##
                n_chains <- dim(draws_array_list[[1]])[2] ; n_chains
                print(paste("convert_list_of_arrays_to_one_big_array - n_iter = ", n_iter))
                ##
                ## Calculate total number of variables:
                ##
                n_variables_total <- 0
                n_variables_per_list <- c()
                param_names_per_list <- list()
                ##
      
                for (l in 1:length(draws_array_list)) {
                  
                      list_or_mat_l <- draws_array_list[[l]]
                      is_matrix <- is.matrix(list_or_mat_l)
                      ##
                      if (is_matrix) n_variables_list_l <- 1
                      else           n_variables_list_l <- dim(list_or_mat_l)[3]
                      ##
                      n_variables_total <- n_variables_total + n_variables_list_l
                      n_variables_per_list[l] <- n_variables_list_l
                      try({  
                        if (!is_matrix) param_names_per_list[[l]] <-  find_all_stan_param_names(list_or_mat_l)
                        else            param_names_per_list[[l]] <-  names(draws_array_list)[l]
                      })
                      ##
                      # dimnames_per_list[[l]] <- find_all_stan_param_names(named_draws_array)
                     
                      
                }
                ##
                ## Make vector of full dim names:
                ##
                full_dim_names_vec <- c(unlist(param_names_per_list))
                ##
                ## Create empty array to hold all variables:
                ##
                combined_array <- array(
                  NA, 
                  dim = c(n_iter, n_chains, n_variables_total),
                  dimnames = list(
                    iteration = 1:n_iter,
                    chain = 1:n_chains,
                    variable = full_dim_names_vec  # Will set this later
                  )
                )
                ##
                ## Combine all variable names:
                ##
                all_var_names <- c()
                ##
                ## Track where to insert values:
                ##
                var_index <- 1
                ##
                # ## Fill the combined array:
                # str(current_array)
                # ##
                # i = 1
                for (l in 1:length(draws_array_list)) {
                  
                      current_array <- draws_array_list[[l]]
                      ##
                      n_vars <- n_variables_per_list[l]
                      ##
                      current_var_names <- param_names_per_list[l]
                      ##
                      ## Copy values to the combined array:
                      ##
                      combined_array[, , var_index:(var_index + n_vars - 1)] <- current_array
                      ##
                      ## Add variable names
                      ##
                      all_var_names <- c(all_var_names, current_var_names)
                      ##
                      ## Update index:
                      ##
                      var_index <- var_index + n_vars
                  
                }
                # str(combined_array)
                ##
                ## Set the variable names for the combined array:
                # ##
                # str(combined_array)
                # # dimnames(combined_array)[[3]] <- c()
                # dimnames(combined_array)[[3]] <- all_var_names
                # dimnames(combined_array)[[2]]
                ##
                # # Apply the draws_array class if the original arrays had it
                # if ("draws_array" %in% class(draws_array_list[[1]])) {
                #   class(combined_array) <- c("draws_array", class(combined_array))
                # }
                ##
                return(combined_array)
          
          } else { 
            
                return(draws_array_list)
            
          }
  
  
}












#' if_null_then_set_to
#' @keywords internal
#' @export
if_null_then_set_to <- function(x, 
                                set_to_thif_null, 
                                debugging = FALSE)  {
  
        if (is.null(x)) {
            y <- set_to_thif_null
        } else { 
            y <- x
        }
  
        if (debugging) print(paste(y))
        
        return(y)

}






#' unpack_list
#' @keywords internal
#' @export
unpack_list <- function(list_obj, 
                        env = parent.frame()) {
  
      for (name in names(list_obj)) {
        assign(name, list_obj[[name]], envir = env)
      }
      
      invisible(NULL)
  
}






#' pack_list
#' @keywords internal
#' @export
pack_list <- function(list_obj, 
                      env = parent.frame()) {
  
      for (name in names(list_obj)) {
        if (exists(name, envir = env)) {
          list_obj[[name]] <- get(name, envir = env)
        }
      }
      
      return(list_obj)
  
}












#' check_vector_length
#' @keywords internal
#' @export
check_vec_length <- function(   input_list, 
                                param_name, 
                                expected_length) {
 
          result <- list()
          result$valid <- TRUE
          ##
          param_value <- input_list[[param_name]]
          ##
          if (!is.null(param_value)) {
            
            if (!is.vector(param_value) || length(param_value) != expected_length) {
              
              result$valid <- FALSE
              result$messages <- c(result$messages,
                                   sprintf("%s must be a vector of length %d, but is %s of length %d",
                                           param_name, expected_length,
                                           class(param_value)[1], length(param_value)))
              stop(print(result$messages))
              
            }
            
          }
          
          return(result)
          
}






#' check_list_of_vectors_length
#' @keywords internal
#' @export
check_list_of_vectors_length <- function(input_list,
                                         param_name, 
                                         expected_length) {
          
          list_length <- length(input_list)
          ##
          result_per_list_element <- list()
          ##
          for (i in 1:list_length) { 
            result_per_list_element[[i]] <- check_vec_length(input_list = input_list,
                                                             param_name = param_name,
                                                             expected_length = expected_length)
          }
          
          return(result_per_list_element)
  
}






#' check_array_dims
#' @keywords internal
#' @export
check_array_dims <- function(   input_list, 
                                param_name, 
                                expected_dims) {
  
          result <- list()
          result$valid <- TRUE
          ##
          param_value <- input_list[[param_name]]
          ##
          if (!is.null(param_value)) {
            
            inputted_dims <- c(base::dim(param_value))
            print(paste("inputted_dims = ", inputted_dims))
            ##
            n_dims_which_dont_match <- sum(as.numeric(inputted_dims != expected_dims))
            print(paste("n_dims_which_dont_match = ", n_dims_which_dont_match))
            
            if (n_dims_which_dont_match != 0) {
              
              result$valid <- FALSE
              result$messages <- c(result$messages,
                                   sprintf("%s must be an array of dimension %d, but is %s of dimension %d",
                                           param_name, expected_dims,
                                           class(param_value)[1], inputted_dims))
              stop(print(result$messages))
              
            }
            
          }
          
          return(result)
  
}







#' check_scalar
#' Function to check scalar value
#' @keywords internal
#' @export
check_scalar <- function(input_list,
                         param_name) {
  
 
          result <- list()
          result$valid <- TRUE
          ##
          param_value <- input_list[[param_name]]
          ##
          if (!is.null(param_value)) {
            
            if (!is.vector(param_value) || length(param_value) != 1) {
              
              result$valid <- FALSE
              result$messages <- c(result$messages, 
                                   sprintf("%s must be a scalar value, but is %s of length %d", 
                                           param_name, class(param_value)[1], length(param_value)))
              stop(print(result$messages))
            }
            
          }
          
          return(result)
  
}








