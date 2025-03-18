


# param_string <- "beta"

#' find_all_param_names
#' @keywords internal
#' @export
find_all_Stan_param_names <- function( draws_array, 
                                  print_findings = FALSE) {
  
            # Get parameter names that contain the string anywhere
            all_Stan_param_names <- dimnames(draws_array)[[3]]
              
            # Filter out parameters that are all NA
            valid_params <- character(0)
            for (param in param_names) {
              if (!all(is.na(draws_array[,,param]))) {
                valid_params <- c(valid_params, param)
              }
            }

            if (length(valid_params) == 0) {
              warning(paste("All parameters containing", param_string, "contain only NAs"))
              return(list())  # Return empty list if no valid parameters
            }
           
            if (print_findings) {
              print(all_Stan_param_names)
            }
            
            return(all_Stan_param_names)
            
}








#' find_param_names_containing_string
#' @keywords internal
#' @export
find_param_names_containing_string <- function( param_names, 
                                                param_string, 
                                                print_findings = FALSE) {
            
            
            # Get parameter names that contain the string anywhere
            param_names_filtered <- param_names[grep(param_string, param_names, fixed = TRUE)]
            
            if (length(param_names_filtered) == 0) {
              stop(paste("No parameters found containing", param_string))
            }
            
            if (print_findings) { 
              print(param_names_filtered)
            }
            
            return(param_names_filtered)
  
}









#' find_param_names_containing_string_batched
#' @keywords internal
#' @export
find_param_names_containing_string_batched <- function(  all_param_names, 
                                                         param_prefixes) {
            
            
            # Convert single string to list if needed
            if (is.character(param_prefixes) && length(param_prefixes) == 1) {
              param_prefixes <- list(param_prefixes)
            }
            
            # Store plots for each parameter group
            matching_param_names_list <- list()
            
            # Generate plots for each parameter prefix
            for (prefix in param_prefixes) {
                    # 
                    # cat("\nTrying prefix:", prefix, "\n")
                    #### pattern <- paste0("^", prefix, "\\[")  
                    pattern <- paste0("^", prefix, "($|\\[)")
                    # cat("Using pattern:", pattern, "\n")
                    
                    matching_params <- grep(pattern, all_param_names, value = TRUE)
                    # cat("Found matching params:", paste(matching_params, collapse=", "), "\n")
                    
                    if(length(matching_params) == 0) {
                      stop(sprintf("No parameters found starting with %s\nAvailable parameters: %s", 
                                   prefix, paste(all_param_names, collapse=", ")))
                    }
                    
                    params <- find_param_names_containing_string( param_names = all_param_names,
                                                                  param_string = prefix)
                    
                    # plots <- plot_param_group_batched(  draws_array = draws_array, 
                    #                                     param_prefix = prefix, 
                    # 
                    matching_param_names_list[[prefix]] <- params
                    
            }
            
            return(matching_param_names_list)
  
}









#' find_params_containing_string
#' @keywords internal
#' @export
find_params_containing_string <- function(draws_array, 
                                          param_string, 
                                          debugging = FALSE) {
  
              ## Get parameter names that contain the string anywhere:
              param_names <- dimnames(draws_array)[[3]][grep(param_string, dimnames(draws_array)[[3]], fixed = TRUE)]
              ##
              param_names_filtered <- find_param_names_containing_string( param_names = param_names,
                                                                          param_string = param_string,
                                                                          debugging = debugging)
              
              str(draws_array)
              draws_array_filtered <- draws_array[,,param_names_filtered]
              str(draws_array_filtered)
              
              # ## Filter out parameters that are all NA:
              # valid_params <- character(0)
              # for (param in param_names_filtered) {
              #   if (!all(is.na(draws_array[,,param]))) {
              #     valid_params <- c(valid_params, param)
              #   }
              # }
              # if (length(valid_params) == 0) {
              #   warning(paste("All parameters containing", param_string, "contain only NAs"))
              #   return(list())  # Return empty list if no valid parameters
              # }
              
              if (debugging) { 
                print(valid_params)
              }
              
              
              return(draws_array_filtered)
    
}

# 
# all_draws_array <- draws_array
# param_prefixes <- c("beta", "raw_scale")
# 







#' combine_draws_arrays
#' @keywords internal
#' @export
combine_draws_arrays <- function(draws_array_list) {
  
          # Verify all arrays have same iterations and chains
          n_iter   <- dim(draws_array_list[[1]])[1]
          n_chains <- dim(draws_array_list[[1]])[2]
          
          # Calculate total number of variables
          total_vars <- sum(sapply(draws_array_list, function(x) dim(x)[3]))
          
          # Create empty array to hold all variables
          combined_array <- array(
            NA, 
            dim = c(n_iter, n_chains, total_vars),
            dimnames = list(
              iteration = 1:n_iter,
              chain = 1:n_chains,
              variable = NULL  # Will set this later
            )
          )
          
          # Combine all variable names
          all_var_names <- c()
          
          # Track where to insert values
          var_index <- 1
          
          # Fill the combined array
          for (i in 1:length(draws_array_list)) {
            current_array <- draws_array_list[[i]]
            n_vars <- dim(current_array)[3]
            current_var_names <- dimnames(current_array)[[3]]
            
            # Copy values to the combined array
            combined_array[, , var_index:(var_index + n_vars - 1)] <- current_array
            
            # Add variable names
            all_var_names <- c(all_var_names, current_var_names)
            
            # Update index
            var_index <- var_index + n_vars
          }
          
          # Set the variable names for the combined array
          dimnames(combined_array)[[3]] <- all_var_names
          
          # Apply the draws_array class if the original arrays had it
          if ("draws_array" %in% class(draws_array_list[[1]])) {
            class(combined_array) <- c("draws_array", class(combined_array))
          }
          
          return(combined_array)
          
}










#' find_param_names_containing_string_batched
#' @keywords internal
#' @export
find_params_containing_string_batched <- function(  all_draws_array,
                                                    param_prefixes) {

            all_param_names <- dimnames(all_draws_array)$variable
            # str(all_param_names)
            ##
            param_names_filtered_list <- find_param_names_containing_string_batched( all_param_names = all_param_names,
                                                                                     param_prefixes = param_prefixes)
            # param_names_filtered_list
            
            drawws_array_filtered_list <- list()
            for (i in 1:length(param_prefixes)) {
                drawws_array_filtered_list[[i]] <- all_draws_array[,,param_names_filtered_list[[i]]]
            }
            
            drawws_array_filtered <- combine_draws_arrays(drawws_array_filtered_list)
            
            return(list(drawws_array_filtered_list = drawws_array_filtered_list,
                        drawws_array_filtered = drawws_array_filtered))
  
}















#' find_param_names_matching_string
#' @keywords internal
#' @export
find_param_names_matching_string <- function(  draws_array, 
                                               param_prefix) {
              
              # # # Get parameter names that start with the prefix
              # # param_names <- dimnames(draws_array)$parameters[grep(paste0("^", param_prefix, "\\["), dimnames(draws_array)$parameters)]
              # 
              # # Get parameter names that start with the prefix
              # ####  param_names <- dimnames(draws_array)[[3]][grep(paste0("^", param_prefix, "\\."), dimnames(draws_array)[[3]])]
              
              # # Get parameter names that start with the prefix
              # param_names <- dimnames(draws_array)[[3]][grep(paste0("^", param_prefix, "\\["), dimnames(draws_array)[[3]])]
              
              # Modified pattern to match both cases:
              # 1. Exact match of prefix (for parameters like "alpha")
              # 2. Prefix followed by [ (for parameters like "beta[1]")
              param_names <- dimnames(draws_array)[[3]][grep(paste0("^", param_prefix, "($|\\[)"), dimnames(draws_array)[[3]])]
              
              if (length(param_names) == 0) {
                stop(paste("No parameters found starting with", param_prefix))
              }
              
              # Filter out parameters that are all NA
              valid_params <- character(0)
              for (param in param_names) {
                if (!all(is.na(draws_array[,,param]))) {
                  valid_params <- c(valid_params, param)
                }
              }
              
              if (length(valid_params) == 0) {
                warning(paste("All parameters starting with", param_prefix, "contain only NAs"))
                return(list())  # Return empty list if no valid parameters
              }
            
              return(valid_params)
  
}









#' find_param_names_containing_string_batched
#' @keywords internal
#' @export
find_param_names_matching_string_batched <- function(    all_param_names, 
                                                         param_prefixes) {
  
  
            
            # Convert single string to list if needed
            if (is.character(param_prefixes) && length(param_prefixes) == 1) {
              param_prefixes <- list(param_prefixes)
            }
            
            # Store plots for each parameter group
            matching_param_names_list <- list()
            
            # Generate plots for each parameter prefix
            for (prefix in param_prefixes) {
              
              cat("\nTrying prefix:", prefix, "\n")
              #### pattern <- paste0("^", prefix, "\\[")  
              pattern <- paste0("^", prefix, "($|\\[)")
              cat("Using pattern:", pattern, "\n")
              
              matching_params <- grep(pattern, all_param_names, value = TRUE)
              cat("Found matching params:", paste(matching_params, collapse=", "), "\n")
              
              if(length(matching_params) == 0) {
                stop(sprintf("No parameters found starting with %s\nAvailable parameters: %s", 
                             prefix, paste(all_param_names, collapse=", ")))
              }
              
              params <- find_param_names_containing_string( param_names = all_param_names,
                                                            param_string = prefix)
              
              # plots <- plot_param_group_batched(  draws_array = draws_array, 
              #                                     param_prefix = prefix, 
              # 
              matching_param_names_list[[prefix]] <- params
              
            }
            
            return(matching_param_names_list)
  
}









 

#' plot_param_group_batched
#' @keywords internal
#' @export
plot_param_group_batched <- function(draws_array, 
                                     param_prefix, 
                                     plot_type = c("density", "trace"),
                                     batch_size = 9) {

    
              plot_type <- match.arg(plot_type)
              
              ## Get param names:
              valid_params <- find_param_names_matching_string(draws_array, param_prefix)
              
              ## Choose plot function:
              plot_func <- switch(plot_type,
                                  density = bayesplot::mcmc_dens,
                                  trace =   bayesplot::mcmc_trace)
              
              ## Plot in batches:
              n_batches <- ceiling(length(valid_params) / batch_size)
              
              plots <- list()
              for (i in 1:n_batches) {
                
                        start_idx <- (i-1) * batch_size + 1
                        end_idx <- min(i * batch_size, length(valid_params))
                        batch_params <- valid_params[start_idx:end_idx]
                        
                        ### plot 
                        plots[[i]] <- plot_func(draws_array[,,batch_params, drop=FALSE]) +
                                      ggplot2::ggtitle(paste0(param_prefix, " (batch ", i, " of ", n_batches, ")"))
                      
              }
              
              return(plots)
      
}










#' plot_multiple_params_batched
#' @keywords internal
#' @export
plot_multiple_params_batched <- function(draws_array, 
                                         param_prefixes, 
                                         plot_type = c("density", "trace"),
                                         batch_size = 9) {
  
  
              # # Debug prints
              # cat("\nDebugging plot_multiple_params_batched:\n")
              # cat("Array dimensions:", paste(dim(draws_array), collapse=" x "), "\n")
              # cat("Parameter names in array:", paste(dimnames(draws_array)[[3]], collapse=", "), "\n")
              # cat("Looking for prefix:", param_prefixes, "\n")
              
              # Get parameter names from array
              param_names <- dimnames(draws_array)[[3]]
  
              # Convert single string to list if needed
              if (is.character(param_prefixes) && length(param_prefixes) == 1) {
                param_prefixes <- list(param_prefixes)
              }
              
              # Store plots for each parameter group
              all_plots <- list()
              
              # Generate plots for each parameter prefix
              for (prefix in param_prefixes) {
                
                    cat("\nTrying prefix:", prefix, "\n")
                    #### pattern <- paste0("^", prefix, "\\[")  
                    pattern <- paste0("^", prefix, "($|\\[)")
                    cat("Using pattern:", pattern, "\n")
                    
                    matching_params <- grep(pattern, param_names, value = TRUE)
                    cat("Found matching params:", paste(matching_params, collapse=", "), "\n")
                    
                    if(length(matching_params) == 0) {
                      stop(sprintf("No parameters found starting with %s\nAvailable parameters: %s", 
                                   prefix, paste(param_names, collapse=", ")))
                    }
                    
                    
                        plots <- plot_param_group_batched(  draws_array = draws_array, 
                                                            param_prefix = prefix, 
                                                            plot_type = plot_type, 
                                                            batch_size = batch_size)
                        
                        all_plots[[prefix]] <- plots
                    
              }
              
              return(all_plots)
  
}






















