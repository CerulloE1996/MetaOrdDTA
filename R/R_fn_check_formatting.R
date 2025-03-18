

#' check_format
#' @keywords internal
#' @export
check_format <- function(   input_list, 
                                param_name, 
                                expected_possible_formats) {
  ##
  result <- list()
  result$valid <- TRUE
  result$messages <- character(0) 
  ##
  param_value <- input_list[[param_name]]
  valid_param_values <- expected_possible_formats[[param_name]]
  ##
  if (!is.null(param_value)) {
        
        if  (!(param_value %in% valid_param_values))  {
          
              # Create a nicely formatted message
              if (is.character(valid_param_values)) {
                # For character values, put quotes around each option
                options_str <- paste(sprintf("'%s'", valid_param_values), collapse = ", ")
              } else {
                # For numeric or other types, just convert to string and join
                options_str <- paste(valid_param_values, collapse = ", ")
              }
              
              
              # Format the message
              error_msg <- sprintf("%s must be one of the following: %s", param_name, options_str)
              result$messages <- c(result$messages, error_msg)
              
              # Stop with a nicely formatted error
              stop(error_msg, call. = FALSE)
              
              # result$valid <- FALSE
              # result$messages <- c(result$messages,
              #                      cat(param_name, "must be one of the following:", valid_param_values))
              # stop(print(result$messages))
              
        }
    
  }
  
  return(result)
  
}




#' R_fn_get_possible_values_list
#' @keywords internal
#' @export
R_fn_get_possible_model_args_list <- function() {
  
        
        {
          possible_values_list <- list()
          ##
          ## Basic and/or mandatory model arguments:
          ##
          # possible_values_list$test_type  <- list("continuous", "cts", "ordinal", "ord")
          possible_values_list$network    <- list(TRUE, 1, FALSE, 0)
          possible_values_list$prior_only <- list(TRUE, 1, FALSE, 0)
          ##
          ## "Advanced" options for ordinal models ONLY:
          ##
          possible_values_list$parameterisation  <- list("Xu", "Gatsonis", "Jones")
          possible_values_list$random_thresholds <- list(TRUE, 1, FALSE, 0)
          ##
          ## "Advanced" options for cts (i.e., Jones) models ONLY:
          ##
          possible_values_list$box_cox <- list(TRUE, 1, FALSE, 0)
          ##
          ## "Advanced" options for ALL models (i.e., ANY MODEL):
          ##
          possible_values_list$softplus <- list(TRUE, 1, FALSE, 0)
        }
  
        return(possible_values_list)
  
}




 

 

#' try_silent
#' @keywords internal
#' @export
try_silent <- function(x) { 
  
  try({  x  }, silent = TRUE) 
  
}



#' check_data_format
#' @keywords internal
#' @export
check_data_format <- function(x, 
                              network) { 
                  
          result <- list()
          result$valid <- TRUE
          result$messages <- character(0) 
          ##
          test_sum <- 0
          is_x_list     <- try_silent(is.list(x))
          is_x_length_2 <- try_silent(length(x) == 2)
          ##
          test_sum <- try_silent(sum(is_x_list, is_x_length_2))
          ##
          error_msg_for_x <- "The data (x) MUST be an R list of length 2, where element #1 corresponds \n to the data in the non-diseased population \n 
                         and element #2 the diseased populatio. \n
                         For meta-analysis (MA), both elements should be a matrix with #rows = #{studies} and #cols = #{observed test thresholds}. \n
                         For network-meta-analysis (NMA), both elements should be a a list (hence x is a list of lists), where each inner list corresponds to \n
                         a different test, and then the format within each inner list is the same as for MA (i.e., for NMA x will be a 'list of lists of matrices'"
          ##
          if (test_sum != 2) { 
            result$test_failed_at <- 1
            result$messages <- error_msg_for_x
          }
          ##
          if (network == TRUE) { 
            
                    try({ 
                      x_first_element <- try_silent(x[[1]])
                      is_list <- is.list(x_first_element)
                      if (!(is_list)) { 
                        result$test_failed_at <- 2
                        result$messages <- error_msg_for_x
                      }
                    })
                    ##
                    try({ 
                      x_first_nested_element <- try_silent(x[[1]][[1]])
                      is_matrix <- is.matrix(x_first_nested_element)
                      if (!(is_matrix)) { 
                        result$test_failed_at <- 3
                        result$messages <- error_msg_for_x
                      }
                    })
                    
          } else if (network == FALSE) { 
            
                    try({ 
                      x_first_element <- try_silent(x[[1]])
                      is_matrix <- is.matrix(x_first_element)
                      if (!(is_matrix)) { 
                        result$test_failed_at <- 4
                        result$messages <- error_msg_for_x
                        
                      }
                    })
            
          }
          
          
          try({  
          test_same <- (result$messages == error_msg_for_x)
          if (test_same == TRUE) { 
            stop(error_msg_for_x) # , call. = FALSE)
          }
          }, silent = TRUE)
          
          return(result)
  
}

 


