










## ---------------------------------------------------------------------------------------------------

#' filtered_param_names_string
#' @keywords internal
#' @export
filtered_param_names_string <-   function( all_param_names_vec, 
                                           param_string, 
                                           condition,
                                           debugging = FALSE) {
  
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
            
            if (debugging) { 
              print(filtered_param_names_vec)
            }
            
            return(filtered_param_names_vec)
  
}





 

#' filter_param_names_string_batched
#' @keywords internal
#' @export
filter_param_names_string_batched <- function(  all_param_names_vec, 
                                                param_string_vec,
                                                condition,
                                                debugging = FALSE) {
            
            
            # Convert single string to list if needed
            if (is.character(param_string_vec) && length(param_string_vec) == 1) {
              param_string_vec <- list(param_string_vec)
            }
            
            # Store in list for each parameter in "all_param_names_vec"
            filtered_param_names_list <- list()
          
            for (param_string_i in param_string_vec) {
                    
                    # if (debugging) { 
                    #    cat("\nTrying param_string:", param_string, "\n")
                    # }
              
                    # #### pattern <- paste0("^", param_string, "\\[")  
                    # pattern <- paste0("^", param_string, "($|\\[)")
                    # if (debugging) cat("Using pattern:", pattern, "\n")
                    
                    params <- filtered_param_names_string( all_param_names_vec = all_param_names_vec,
                                                           param_string = param_string_i,
                                                           condition = condition, 
                                                           debugging = debugging)

                   
                    filtered_param_names_list[[param_string_i]] <- params
                    
            }
            
            return(filtered_param_names_list)
  
}








#' filter_param_draws_string
#' @keywords internal
#' @export
filter_param_draws_string <- function(  named_draws_array,
                                        param_string, 
                                        condition,
                                        exclude_NA = FALSE,
                                        debugging = FALSE) {
  
  
              named_draws_array_2 <- format_named_array_for_bayesplot(named_draws_array)
  
              ## Get parameter names that contain the string ANYWHERE:
              all_param_names_vec <- dimnames(named_draws_array_2)[[3]] 
              ## [grep(param_string, dimnames(named_draws_array_2)[[3]], fixed = TRUE)]
              
              {
                  draws_array_filtered <- named_draws_array_2
                  try({  
                 
                      ## Filter parameter names (note this WILL include params which have "ALL NA's":
                      param_names_filtered <- filtered_param_names_string(  all_param_names_vec = all_param_names_vec,
                                                                            param_string = param_string,
                                                                            condition = condition,
                                                                            debugging = debugging)
                      
                      ## Then filter/subset "named_draws_array_2":
                      draws_array_filtered <- named_draws_array_2[,,param_names_filtered]
                  })
              }
              
              ## Then exclude any draws which are ALL NA's (if chosen):
              if (!(exclude_NA)) { 
                        
                        if (debugging) { 
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
                        
                        if (debugging) { 
                          print(draws_array_filtered_wo_NAs)
                        }
                        
                        return(draws_array_filtered_wo_NAs)
                
              }
    
}


# named_draws_array <- draws_array
# param_strings_vec <- c("beta", "raw_scale")
# 
# 
# 
# 
# 
# named_draws_array = named_draws_array
# param_strings_vec = grouped_params$bs_names_main
# condition = "exact_match"
# exclude_NA = exclude_NA
# debugging = TRUE


#' filter_param_draws_string_batched
#' @keywords internal
#' @export
filter_param_draws_string_batched <- function(    named_draws_array,
                                                  param_strings_vec,
                                                  condition,
                                                  exclude_NA = FALSE,
                                                  debugging = FALSE) {

            param_draws_filtered_list <- list()
            
            for (param_string_i in param_strings_vec) {
                
                try({ 
              
                    param_draws_filtered_list[[param_string_i]] <- filter_param_draws_string( named_draws_array = named_draws_array,
                                                                                              param_string = param_string_i, 
                                                                                              condition = condition,
                                                                                              exclude_NA = exclude_NA,
                                                                                              debugging = debugging)
                })
              
            }

            dimnames(named_draws_array)[[3]]
            # str(param_draws_filtered_list)
            
            if (is.list(param_draws_filtered_list)) {
              
                      ## Turn list into one big 3D array:
                      param_draws_filtered_array <- NULL
                      try({ 
                        param_draws_filtered_array <- convert_list_of_arrays_to_one_big_array( 
                                                                          draws_array_list = param_draws_filtered_list )
                        param_draws_filtered_array    
                      })
                      
                      # param_draws_filtered_list <- list()
                      # for (i in 1:length(param_strings_vec)) {
                      #     param_draws_filtered_list[[i]] <- named_draws_array[,,param_names_filtered_list[[i]]]
                      # }
                      # 
                      # 
                      # 
                      # ## Turn list into one big 3D array:
                      # param_draws_filtered_array <- convert_list_of_arrays_to_one_big_array(param_draws_filtered_list)
                      # 
                      return(list(param_draws_filtered_list = param_draws_filtered_list,
                                  param_draws_filtered_array = param_draws_filtered_array))
            
            } else { 
                  
                      param_draws_filtered_array <- param_draws_filtered_list
                      return(list(param_draws_filtered_array = param_draws_filtered_array))
              
            }
            
  
}



# 
# 
# 
# data.table::rbindlist(param_draws_filtered_list)
# 
# 
# 
#  














