




#' R_fn_get_stan_model_file_name
#' @keywords internal
#' @export
R_fn_get_stan_model_file_name <- function( cts,
                                           network,
                                           prior_only,
                                           ##
                                           model_parameterisation,
                                           random_thresholds,
                                           Dirichlet_random_effects_type
                                           
)  {
  
        ##
        ## Then get the file name string:
        ##
        if (network == TRUE) { 
          file_string_1 <- "DTA_NMA_Nyaga_"
        } else { 
          file_string_1 <- "DTA_MA_"
        }
        ##
        if (cts == FALSE) {
          
                  ##
                  ## Now get the Stan model file name string:
                  ##
                  if (model_parameterisation == "Gatsonis") { 
                    file_string_2 <- "Gat_"
                  } else { 
                    file_string_2 <- "Xu_"
                  }
                  ##
                  if (random_thresholds == FALSE) { 
                    file_string_3 <- "FIXEDthr"
                  } else { 
                    file_string_3 <- "RANDthr"
                  }
                  ##
                  if (random_thresholds == TRUE) { 
                      if (Dirichlet_random_effects_type == "SD") { 
                        file_string_4 <- "_SD"
                      } else { 
                        file_string_4 <- "_kappa"
                      }
                  } else { 
                      file_string_4 <- "" ## no string here if fixed thresholds.
                      Dirichlet_random_effects_type <- "none"
                  }
                  ##
                  ## Now construct the string:
                  ##
                  stan_model_file_name <- paste0(file_string_1, 
                                                 file_string_2, 
                                                 file_string_3,
                                                 file_string_4,
                                                 ".stan")
           
          
          
        } else if (cts == TRUE) {
          
                  ##
                  ## Now get the Stan model file name string:
                  ##
                  file_string_2 <- "Jones"
                  ##
                  ## Now construct the string:
                  ##
                  stan_model_file_name <- paste0(file_string_1, 
                                                 file_string_2, 
                                                 ".stan")
                  
        }
        ##
        if (prior_only) {
          
                  stan_model_file_name <- paste0("PO_", stan_model_file_name)
          
        }
  
        print(paste("stan_model_file_name = ", stan_model_file_name))
        ##
        ## Output:
        ##
         ### return(stan_model_file_name)
        return(list(stan_model_file_name = stan_model_file_name))

}










 




















