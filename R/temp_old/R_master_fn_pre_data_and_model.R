




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





#' check_model_args_list_format - Function to check that the "main" and/or "advanced" user-inputted modelling options are in the correct format:
#' @keywords internal
#' @export
check_model_args_list_format <- function(model_args_list) {
        
  
        {
          possible_values_list <- list()
          ##
          ## Basic and/or mandatory model arguments:
          ##
          possible_values_list$test_type  <- list("continuous", "cts", "ordinal", "ord")
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
        ##
        model_args_list_param_names <- list("test_type", "network", "prior_only", 
                                            "parameterisation", "random_thresholds", 
                                            "box_cox",
                                            "softplus")
        n_model_args <- length(model_args_list_param_names)
        ##
        ## Check format of arguments given in "model_args_list":
        ##
        for (i in 1:n_model_args) {
          result <- check_format( input_list = model_args_list, 
                        param_name = model_args_list_param_names[[i]],
                        expected_possible_formats = possible_values_list)
        }
        
        return(list(result = result, 
                    possible_values_list = possible_values_list))
  
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



#' default_network
#' @keywords internal
#' @export
default_network <- function(x) { 
  
  
      network <- NA
      ##
      is_x_list <- try_silent(is.list(x))
      is_x_nested_list <- try_silent(is.list(x[[1]]))
      ##
      if (is_x_list) { 
        network <- FALSE
      } else if (is_x_nested_list) { 
        network <- TRUE
      }
      
      return(network)
  
  
  
}




default_parameterisation <- function() { 
  
}

default_random_thresholds <- function() { 
  
}

default_softplus <- function() { 
  
}

default_box_cox <- function() { 
  
}


prior_only <- NULL
model_args_list <- NULL

#' prep_data_and_model
#' @keywords internal
#' @export
prep_data_and_model <- function(x, 
                                test_type,
                                network,
                                prior_only,
                                ##
                                model_args_list = NULL,
                                ##
                                priors,
                                ##
                                ## "Advanced" options for ordinal models ONLY:
                                ##
                                parameterisation = NULL,
                                random_thresholds = NULL,
                                ##
                                ## "Advanced" options for cts (i.e., Jones) models ONLY:
                                ##
                                box_cox = NULL,
                                ##
                                ## "Advanced" options for ALL models:
                                ##
                                softplus = NULL,
                                ##
                                ## Other stuff:
                                ##
                                init_lists_per_chain,
                                n_chains
) {
  
  
 

          ##
          ## Initialise "model_args_list" object (if not already initialised):
          ##
          model_args_list <- if_null_then_set_to(model_args_list, list())
          ##
          ## Put the (if any) user-inputted model arguments inside "model_args_list" (if/when available):
          ##
          try({  
              model_args_list$test_type  <- if_null_then_set_to(test_type, "ordinal")
              model_args_list$network    <- if_null_then_set_to(network,  default_network(x))
              model_args_list$prior_only <- if_null_then_set_to(prior_only, FALSE)
              ##
              ## "Advanced" options for ordinal models ONLY:
              ##
              model_args_list$model_parameterisation  <- if_null_then_set_to(parameterisation,  default_parameterisation)
              model_args_list$random_thresholds       <- if_null_then_set_to(random_thresholds, default_random_thresholds)
              ##
              ## "Advanced" options for cts (i.e., Jones) models ONLY:
              ##
              model_args_list$use_box_cox             <- if_null_then_set_to(box_cox, default_box_cox)
              ##
              ## "Advanced" options for ALL models:
              ##
              model_args_list$use_softplus <- if_null_then_set_to(softplus, default_softplus)
      
          }, silent = TRUE)
          ##
          ## Check that arguments inside "model_args_list" are in the correct format:
          ##
          check_model_args_list_format_outs <- check_model_args_list_format(model_args_list = model_args_list)
          # check_model_args_list_format_outs$result
          ##
          ## Check DATA (x) format (do this AFTER checking model args):
          ##
          check_data_format(x, model_args_list$network)


          ##
          ## -----------------  Setup data: ---------------------------------------------------------------------------
          ##
          is_network <- model_args_list$network
          model_args_list$is_network <- is_network
          ##
          if (test_type %in% c("cts", "continuous")) { 
            is_cts <- TRUE
          } else { 
            is_cts <- FALSE
          }
          ##
          data_setup_outs <- ifelse(is_network, R_fn_prep_NMA_data(x = x), R_fn_prep_MA_data(x = x))
          ##
          Stan_data_list <- data_setup_outs$Stan_data_list
          n_studies <- data_setup_outs$n_studies
          n_cat <- data_setup_outs$n_cat
          n_thr <- data_setup_outs$n_thr
          ##
          ## update model_args_list:
          ##
          model_args_list$Stan_data_list <- Stan_data_list
          model_args_list$n_studies <- n_studies
          model_args_list$n_cat <- n_cat
          model_args_list$n_thr <- n_thr
          ##
          model_args_list$is_cts <- is_cts
          model_args_list$is_network <- is_network
          ##
          # ##
          # ## ----- Get / set basic model info:
          # ##
          # info_outs <- R_fn_get_and_set_basic_model_info( Model_type = Model_type,
          #                                                 model_args_list = model_args_list)
          # ##
          # # is_cts <- info_outs$is_cts
          # # is_network <- info_outs$is_network
          # ##
          # ## update model_args_list:
          # ##
          # model_args_list <- info_outs$model_args_list
          ##
          ## ----- Get the Stan model file (.stan file) name: ----------------------------------------------------------
          ##
          Stan_model_name_outs <- R_fn_get_Stan_model_file_name(Model_type, 
                                                                model_args_list)
          ##
          Stan_model_file_name <- Stan_model_name_outs$Stan_model_file_name
          model_args_list <- Stan_model_name_outs$model_args_list
          ##
          ## ----------------- Compile Stan model (if necessary): ------------------------------------------------------
          ##
          outs_compile <- R_fn_compile_Stan_model_basic_given_file_name( Stan_model_file_name = Stan_model_file_name,
                                                                         debugging = TRUE,
                                                                         model_args_list = model_args_list)
          stan_model_obj <- outs_compile$stan_model_obj
          stan_functions_dir <- outs_compile$stan_functions_dir
          ##
          model_args_list <- outs_compile$model_args_list
          ##
          ## -----------------  Set priors (using R package defaults): --------------------------------------------------
          ##
          outs_priors <- R_fn_set_priors_MA(Model_type, 
                                            model_args_list, 
                                            priors)
          priors <-  outs_priors$priors
          priors
          ##
          ## ----------------- Set initial values (using R package defaults) - AFTER priors have been set!: -------------
          ##
          if (is.null(init_lists_per_chain)) { 
            init_lists_per_chain <- list()
            for (i in 1:n_chains) {
              init_lists_per_chain[[i]] <- list(1) ## placeholder (NULL doesnt work here)
            }
          } else { 
            if (length(init_lists_per_chain) != n_chains) { 
              stop("init_lists_per_chain must either be set to 'NULL' (default inits)
                           or be a list of lists of length 'n_chains")
            }
          }
          ##
          for (i in 1:n_chains) {
            outs_inits <- R_fn_set_inits_MA(  Model_type = Model_type,
                                              model_args_list = model_args_list,
                                              priors = priors,
                                              inits = init_lists_per_chain[[i]])
            ##
            init_lists_per_chain[[i]] <- outs_inits$inits
            # model_args_list$inits <- inits
          }
          init_lists_per_chain
          ##
          ## ----------------- Output list: -----------------------------------------------------------------------------
          ##
          return(list( model_args_list = model_args_list,
                       # user_input_model_args_list = user_input_model_args_list,
                       ##
                       priors = priors,
                       init_lists_per_chain = init_lists_per_chain,
                       ##
                       Stan_data_list = Stan_data_list,
                       Stan_model_file_name = Stan_model_file_name,
                       stan_model_obj = stan_model_obj,
                       stan_functions_dir = stan_functions_dir,
                       ##
                       is_cts = is_cts,
                       is_network = is_network,
                       ##
                       n_studies = n_studies,
                       n_cat = n_cat,
                       n_thr = n_thr))
          
}




