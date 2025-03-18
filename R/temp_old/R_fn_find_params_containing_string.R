


# param_string <- "beta"

#' find_all_Stan_param_names
#' @keywords internal
#' @export
find_all_Stan_param_names <- function( stan_draws_array, 
                                       exclude_NA = FALSE,
                                       print = FALSE) {
  
            ## Get parameter names that contain the string anywhere:
            all_Stan_param_names_inc_NAs <- dimnames(stan_draws_array)[[3]]
            
            if (exclude_NA) {
              
                      ## Filter out parameters that are all NA:
                      param_names_wo_NAs <- character(0)
                      for (param in all_Stan_param_names_inc_NAs) {
                        if (!all(is.na(stan_draws_array[,,param]))) {
                          param_names_wo_NAs <- c(param_names_wo_NAs, param)
                        }
                      }
          
                      if (length(param_names_wo_NAs) == 0) {
                        warning(paste("All parameters containing", param_string, "contain only NAs"))
                        return(list())  # Return empty list if no valid parameters
                      }
                
            } else { 
              
                      param_names_wo_NAs <- all_Stan_param_names_inc_NAs
              
            }
           
            if (print) {
              print(param_names_wo_NAs)
            }
            
            return(param_names_wo_NAs)
            
}










#' combine_draws_arrays
#' @keywords internal
#' @export
combine_draws_arrays <- function(draws_array_list) {
  
            ##
            ## Verify all arrays have same iterations and chains:
            ##
            n_iter   <- dim(draws_array_list[[1]])[1] ; n_iter
            n_chains <- dim(draws_array_list[[1]])[2] ; n_chains
            ##
            ## Calculate total number of variables:
            ##
            total_vars <- sum(sapply(draws_array_list, function(x) dim(x)[3]))
            ##
            ## Create empty array to hold all variables:
            ##
            combined_array <- array(
              NA, 
              dim = c(n_iter, n_chains, total_vars),
              dimnames = list(
                iteration = 1:n_iter,
                chain = 1:n_chains,
                variable = NULL  # Will set this later
              )
            )
            #### str(combined_array)
            ##
            ## Combine all variable names:
            ##
            all_var_names <- c()
            ##
            ## Track where to insert values:
            ##
            var_index <- 1
            ##
            ## Fill the combined array:
            ##
            #### i = 1
            for (i in 1:length(draws_array_list)) {
              
                    current_array <- draws_array_list[[i]]
                    n_vars <- dim(current_array)[3] ; n_vars
                    current_var_names <- dimnames(current_array)[[3]]
                    
                    # Copy values to the combined array
                    combined_array[, , var_index:(var_index + n_vars - 1)] <- current_array
                    
                    # Add variable names
                    all_var_names <- c(all_var_names, current_var_names)
                    
                    # Update index
                    var_index <- var_index + n_vars
                  
            }
            ##
            ## Set the variable names for the combined array:
            ##
            dimnames(combined_array)[[3]] <- all_var_names
            ##
            # # Apply the draws_array class if the original arrays had it
            # if ("draws_array" %in% class(draws_array_list[[1]])) {
            #   class(combined_array) <- c("draws_array", class(combined_array))
            # }
            ##
            return(combined_array)
  
}











## ---------------------------------------------------------------------------------------------------

#' filtered_param_names_string
#' @keywords internal
#' @export
filtered_param_names_string <-   function( all_param_names_vec, 
                                         param_string, 
                                         condition,
                                         print = FALSE) {
  
            if (condition == "exact_match") { 
              
                 filtered_param_names_vec <- all_param_names_vec[grep(paste0("^", param_string, "($|\\[)"), all_param_names_vec)]
              
            } else if (condition == "containing") {
              
                 filtered_param_names_vec <- all_param_names_vec[grep(param_string, all_param_names_vec, fixed = TRUE)]
              
            } else { 
              
                 stop("ERROR: 'condition' argument must be set to either 'exact_match' or 'containing'")
              
            }

            if (length(filtered_param_names_vec) == 0) {
              warning(paste("No parameters found for condition:", condition, ", for the string:", param_string))
            }
            
            if (print) { 
              print(filtered_param_names_vec)
            }
            
            return(filtered_param_names_vec)
  
}





 

#' filter_param_names_string_batched
#' @keywords internal
#' @export
filter_param_names_string_batched <- function(  all_param_names_vec, 
                                                param_strings_vec,
                                                condition,
                                                print = FALSE) {
            
            
            # Convert single string to list if needed
            if (is.character(param_strings_vec) && length(param_strings_vec) == 1) {
              param_strings_vec <- list(param_strings_vec)
            }
            
            # Store in list for each parameter in "all_param_names_vec"
            filtered_param_names_list <- list()
          
            for (param_string in param_strings_vec) {
                    
                    if (print) { 
                       cat("\nTrying param_string:", param_string, "\n")
                    }
              
                    # #### pattern <- paste0("^", param_string, "\\[")  
                    # pattern <- paste0("^", param_string, "($|\\[)")
                    # if (print) cat("Using pattern:", pattern, "\n")
                    
                    params <- filtered_param_names_string( all_param_names_vec = all_param_names_vec,
                                                           param_string = param_string,
                                                           condition = condition, 
                                                           print = print)

                   
                    filtered_param_names_list[[param_string]] <- params
                    
            }
            
            return(filtered_param_names_list)
  
}








#' filter_param_draws_string
#' @keywords internal
#' @export
filter_param_draws_string <- function(  all_draws_array, 
                                        param_string, 
                                        condition,
                                        exclude_NA = FALSE,
                                        print = FALSE) {
  
              ## Get parameter names that contain the string ANYWHERE:
              all_param_names_vec <- dimnames(all_draws_array)[[3]] ## [grep(param_string, dimnames(all_draws_array)[[3]], fixed = TRUE)]
              
              
              {
                  draws_array_filtered <- all_draws_array
                  try({  
                 
                      ## Filter parameter names (note this WILL include params which have "ALL NA's":
                      param_names_filtered <- filtered_param_names_string(  all_param_names_vec = all_param_names_vec,
                                                                            param_string = param_string,
                                                                            condition = condition,
                                                                            print = print)
                      
                      ## Then filter/subset "all_draws_array":
                      draws_array_filtered <- all_draws_array[,,param_names_filtered]
                  })
              }
              
              ## Then exclude any draws which are ALL NA's (if chosen):
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


# all_draws_array <- draws_array
# param_strings_vec <- c("beta", "raw_scale")
# 
# 
# 
# 
# 
# all_draws_array = all_draws_array
# param_strings_vec = grouped_params$bs_names_main
# condition = "exact_match"
# exclude_NA = exclude_NA
# print = TRUE


#' filter_param_draws_string_batched
#' @keywords internal
#' @export
filter_param_draws_string_batched <- function(    all_draws_array,
                                                  param_strings_vec,
                                                  condition,
                                                  exclude_NA = FALSE,
                                                  print = FALSE) {

            param_draws_filtered_list <- list()
            
            for (param_string_i in param_strings_vec) {
              
              try({ 
            
                  param_draws_filtered_list[[param_string_i]] <- filter_param_draws_string( all_draws_array = all_draws_array,
                                                                                            param_string = param_string_i, 
                                                                                            condition = condition,
                                                                                            exclude_NA = exclude_NA,
                                                                                            print = print)
              })
              
            }

            
            ## Turn list into one big 3D array:
            param_draws_filtered_array <- NULL
            try({ 
              param_draws_filtered_array <- combine_draws_arrays( draws_array_list = param_draws_filtered_list )
              param_draws_filtered_array
            })
            
            # param_draws_filtered_list <- list()
            # for (i in 1:length(param_strings_vec)) {
            #     param_draws_filtered_list[[i]] <- all_draws_array[,,param_names_filtered_list[[i]]]
            # }
            # 
            # 
            # 
            # ## Turn list into one big 3D array:
            # param_draws_filtered_array <- combine_draws_arrays(param_draws_filtered_list)
            # 
            return(list(param_draws_filtered_list = param_draws_filtered_list,
                        param_draws_filtered_array = param_draws_filtered_array))
  
}











 














