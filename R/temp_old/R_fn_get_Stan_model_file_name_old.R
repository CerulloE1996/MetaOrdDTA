





#' R_fn_get_basic_model_info
#' @keywords internal
#' @export
R_fn_get_and_set_basic_model_info <- function(  Model_type,
                                                model_args_list
                                       
)  {

          ##
          ## First determine whether model is for MA or NMA:
          ##
          is_network <- stringr::str_detect(Model_type, "network")
          print(paste("is_network = ", is_network))
          ##
          ## Then determine whether model is for ordinal test(s) or continuous test(s):
          ##
          is_cts <- stringr::str_detect(Model_type, "cts")
          print(paste("is_cts = ", is_cts))
          ##
          ## Then update model_args_list (if necessary):
          ##
          model_args_list$prior_only <- if_null_then_set_to(x = model_args_list$prior_only, FALSE)
          model_args_list$use_softplus_for_scales <- if_null_then_set_to(x = model_args_list$use_softplus_for_scales, TRUE)
          ##
          if (is_cts == FALSE) {
            
                model_args_list$ord_model_parameterisation <- if_null_then_set_to(x = model_args_list$ord_model_parameterisation, "Gatsonis")
                model_args_list$random_thresholds <- if_null_then_set_to(x = model_args_list$random_thresholds, FALSE)
                model_args_list$Dirichlet_random_effects_type <- if_null_then_set_to(x = model_args_list$Dirichlet_random_effects_type, "SD")
            
          } else if (is_cts == TRUE) {
            
                try({ 
                  rm(model_args_list$ord_model_parameterisation)
                  rm(model_args_list$random_thresholds)
                  rm(model_args_list$Dirichlet_random_effects_type)
                }, silent = TRUE)
                ##
                model_args_list$use_box_cox <- if_null_then_set_to(model_args_list$use_box_cox, TRUE)
            
          }
          ##
          ## Output list:
          ##
          return(list( model_args_list = model_args_list, 
                       is_cts = is_cts, 
                       is_network = is_network))
  
  
}








#' R_fn_get_Stan_model_file_name
#' @keywords internal
#' @export
R_fn_get_Stan_model_file_name <- function( Model_type,
                                           model_args_list
)  {
  
        ##
        ## First update model_args_list (if necessary) and get basic model info:
        ##
        outs <- R_fn_get_and_set_basic_model_info(Model_type, model_args_list)
        model_args_list <- outs$model_args_list
        print(paste("model_args_list = ", model_args_list))
        ##
        is_cts          <- outs$is_cts
        print(paste("is_cts = ", is_cts))
        ##
        is_network      <- outs$is_network
        print(paste("is_network = ", is_network))
        ##
        is_prior_only <- outs$model_args_list$prior_only
        ## Then get the file name string:
        ##
        if (is_network == TRUE) { 
          file_string_1 <- "DTA_NMA_Nyaga_"
        } else { 
          file_string_1 <- "DTA_MA_"
        }
        
        
        if (is_cts == FALSE) {
          
                  ##
                  ## Now get the Stan model file name string:
                  ##
                  if (model_args_list$ord_model_parameterisation == "Gatsonis") { 
                    file_string_2 <- "Gat_"
                  } else { 
                    file_string_2 <- "Xu_"
                  }
                  ##
                  if (model_args_list$random_thresholds == FALSE) { 
                    file_string_3 <- "FIXEDthr"
                  } else { 
                    file_string_3 <- "RANDthr"
                  }
                  ##
                  if (model_args_list$random_thresholds == TRUE) { 
                      if (model_args_list$Dirichlet_random_effects_type == "SD") { 
                        file_string_4 <- "_SD"
                      } else { 
                        file_string_4 <- "_kappa"
                      }
                  } else { 
                      file_string_4 <- "" ## no string here if fixed thresholds.
                      model_args_list$Dirichlet_random_effects_type <- "none"
                  }
                  ##
                  ## Now construct the string:
                  ##
                  Stan_model_file_name <- paste0(file_string_1, 
                                                 file_string_2, 
                                                 file_string_3,
                                                 file_string_4,
                                                 ".stan")
                  print(paste("Stan_model_file_name = ",Stan_model_file_name))
          
          
        } else if (is_cts == TRUE) {
          
                  ##
                  ## Now get the Stan model file name string:
                  ##
                  file_string_2 <- "Jones"
                  ##
                  ## Now construct the string:
                  ##
                  Stan_model_file_name <- paste0(file_string_1, 
                                                 file_string_2, 
                                                 ".stan")
                  print(paste("Stan_model_file_name = ",Stan_model_file_name))
                  
        }
        ##
        if (is_prior_only) {
          
                 Stan_model_file_name <- paste0("PO_", Stan_model_file_name)
          
        }
          
        
        ##
        ## Output list:
        ##
        return(list( Stan_model_file_name = Stan_model_file_name, 
                     model_args_list = model_args_list))

}










 




















