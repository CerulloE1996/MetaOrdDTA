
 


 





#' find_param_names_matching_string
#' @keywords internal
#' @export
find_param_names_matching_string <- function( all_param_names, 
                                              param_string, 
                                              print = FALSE) {
  
            filtered_param_names <- all_param_names[grep(paste0("^", param_string, "($|\\[)"), all_param_names)]
            
            if (length(filtered_param_names) == 0) {
              stop(paste("No parameters found starting with", param_string))
            }
            
            if (print) { 
              print(filtered_param_names)
            }
            
            return(filtered_param_names)
  
}






#' find_param_names_matching_string_batched
#' @keywords internal
#' @export
find_param_names_matching_string_batched <- function(  all_param_names_vec,
                                                       param_strings,
                                                       print = FALSE) {
  
            ## Convert single string to list if needed
            if (is.character(param_strings) && length(param_strings) == 1) {
              param_strings <- list(param_strings)
            }
            
            ## Store in list for each parameter in "all_param_names_vec"
            filtered_param_names_list <- list()
         
            for (string in param_strings) {
                   
                    if (print) { 
                      cat("\nTrying string:", string, "\n")
                    }
                    ##
                    #### pattern <- paste0("^", string, "\\[")  
                    pattern <- paste0("^", string, "($|\\[)")
                    if (print) cat("Using pattern:", pattern, "\n")
                    
                    # matching_params <- grep(pattern, all_param_names_vec, value = TRUE)
                    # ##
                    # if (print) cat("Found matching params:", paste(matching_params, collapse=", "), "\n")
                    # 
                    # if(length(matching_params) == 0) {
                    #   stop(sprintf("No parameters found starting with %s\nAvailable parameters: %s", 
                    #                string, paste(all_param_names_vec, collapse=", ")))
                    # }
                    
                    params <- find_param_names_matching_string(   all_param_names = all_param_names_vec,
                                                                  param_string = string, 
                                                                  print = print)
                
                    filtered_param_names_list[[string]] <- params
              
            }
            
            if (print) { 
              print(filtered_param_names_list)
            }
            
            return(filtered_param_names_list)
  
}





#' find_params_matching_string
#' @keywords internal
#' @export
find_param_draws_matching_string <- function(  all_draws_array, 
                                               param_string, 
                                               exclude_NA = FALSE,
                                               print = FALSE) {
  
            ## Get parameter names that contain the string anywhere:
            param_names <- dimnames(all_draws_array)[[3]][grep(param_string, dimnames(all_draws_array)[[3]], fixed = TRUE)]
            ##
            param_names_filtered <- find_param_names_matching_string(   param_names = param_names,
                                                                        param_string = param_string,
                                                                        print = print)
            
            ## Then filter/subset "all_draws_array":
            draws_array_filtered <- all_draws_array[,,param_names_filtered]
            
            if (!(exclude_NA)) { 
              
                  if (print) { 
                    print(draws_array_filtered)
                  }
                  
                  return(draws_array_filtered)
              
            } else if (exclude_NA == TRUE) { 
              
                  ## Filter out parameters that are all NA:
                  param_names_filtered_wo_NAs <- character(0)
                  for (param in param_names_filtered) {
                    if (!all(is.na(draws_array_filtered[,,param]))) {
                      param_names_filtered_wo_NAs <- c(param_names_filtered_wo_NAs, param)
                    }
                  }
                  if (length(param_names_filtered_wo_NAs) == 0) {
                    warning(paste("All parameters containing", param_string, "contain only NAs"))
                    return(list())  # Return empty list if no valid parameters
                  }
                  draws_array_filtered_wo_NAs <- draws_array_filtered[,,param_names_filtered_wo_NAs]
                  
                  if (print) { 
                    print(draws_array_filtered_wo_NAs)
                  }
                  
                  return(draws_array_filtered_wo_NAs)
              
            }
  
}






#' find_param_draws_matching_string_batched
#' @keywords internal
#' @export
find_param_draws_matching_string_batched <- function(   all_draws_array,
                                                        param_prefixes,
                                                        exclude_NA = FALSE,
                                                        print = FALSE) {
              
              all_param_names <- dimnames(all_draws_array)$variable
              head(all_param_names)
              # str(all_param_names)
              ##
              param_names_filtered_list <- find_param_names_matching_string_batched( all_param_names = all_param_names,
                                                                                       param_prefixes = param_prefixes)
              # param_names_filtered_list
              
              drawws_array_filtered_list <- list()
              for (i in 1:length(param_prefixes)) {
                drawws_array_filtered_list[[i]] <- all_draws_array[,,param_names_filtered_list[[i]]]
              }
              
              drawws_array_filtered <- combine_draws_arrays(drawws_array_filtered_list)
              
              out_list <- list( drawws_array_filtered_list = drawws_array_filtered_list,
                                drawws_array_filtered = drawws_array_filtered)
              
              if (print) { 
                print(out_list)
              }
              
              return(out_list)
  
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






















