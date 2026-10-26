
#' R_fn_get_stan_model_file_name
#' @keywords internal
#' @export
R_fn_get_stan_model_file_name <- function( custom_file_name,
                                           use_Se_Sp_prior,
                                           ##
                                           cts,
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
                  if (grepl("ord_HSROC|ordinal_HSROC|R&G|HSROC|Gatsonis",
                            model_parameterisation,
                            ignore.case = TRUE)) { 
                    file_string_2 <- "ord_HSROC_"
                  } else { 
                    file_string_2 <- "ord_biv_"
                  }
                  ##
                  if (random_thresholds == FALSE) { 
                    file_string_3 <- "FIXEDthr"
                  } else { 
                    file_string_3 <- "RANDthr"
                  }
                  ##
                  if (random_thresholds == TRUE) { 
                      rand_type <- Dirichlet_random_effects_type
                      if (rand_type == "SD") { 
                        file_string_4 <- "_SD"
                      } else if (rand_type == "kappa") { 
                        file_string_4 <- "_kappa"
                      } else if (rand_type == "alpha") { 
                        file_string_4 <- "_alpha"
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
        # if (prior_only) {
        #   
        #         stan_model_file_name <- paste0("PO_", stan_model_file_name)
        #   
        # }
        ## (prior_only = TRUE runs the model's own Stan file with its likelihood switched off - every study's
        ##  K_fold_CV_indicator = 0, set in initialise_update_run_model() - so that the draws are the model's
        ##  real induced prior; the PO_ files of the prior_only folders are older copies, not the model)
  
        print(paste("stan_model_file_name = ", stan_model_file_name))
        ##
        ## Output:
        ##
         ### return(stan_model_file_name)
        
        if (!(is.null(custom_file_name))) { 
          stan_model_file_name <- custom_file_name
        }
        ##
        ## The model's _Se_Sp_prior file when the user gives direct Se / Sp priors (priors$Se_Sp_prior):
        ##
        if (use_Se_Sp_prior == TRUE && !grepl("_Se_Sp_prior", stan_model_file_name, fixed = TRUE)) {
          stan_model_file_name <- sub("(_NicoStan_fused_AVX)?[.]stan$",
                                      "_Se_Sp_prior\\1.stan",
                                      stan_model_file_name)
        }
        
        return(list(stan_model_file_name = stan_model_file_name))

}










 




















