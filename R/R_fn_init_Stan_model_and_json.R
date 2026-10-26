
  
  

 




#' init_stan_model_external
#' @keywords internal
#' @export
init_stan_model_external <- function( debugging = FALSE,
                                      stan_data_list, 
                                      stan_model_file_path,
                                      stan_model_obj,
                                      n_chains,
                                      init_lists_per_chain
) {
  
        ## Get package directory paths
        outs <- get_MetaOrdDTA_stan_paths()
        pkg_root_dir <- outs$pkg_root_dir
        pkg_data_dir <- outs$pkg_data_dir
        pkg_stan_dir <- outs$pkg_stan_dir
        ##
        ## JSON:
        ##
        json_file_path <- convert_stan_data_list_to_JSON( stan_data_list = stan_data_list,
                                                          pkg_data_dir = pkg_data_dir)
        ##
        validated_json_string <- paste(readLines(json_file_path), collapse="")
        jsonlite::write_json(x = validated_json_string, path = json_file_path)
        ## (stan_model_obj is NULL when the model is not compiled with cmdstanr, e.g. for a NicoStan fit)

        ## Sampling five throw-away iterations used to happen here merely to inspect
        ## the model variables. Apart from doubling work, that run could fail and
        ## leave the real fit in a partially initialised state. CmdStanR exposes the
        ## variable metadata directly from the compiled model, so no pseudo-fit is
        ## required.
        stan_model_file_path <- file.path(stan_model_file_path)
        json_file_path       <- normalizePath(json_file_path)
        stan_model_file_path <- normalizePath(stan_model_file_path)
        ##
        ## Output:
        ##
        return(list(
          stan_init_pseudo_sampling_outs = NULL,
          stan_model_obj = stan_model_obj,
          json_file_path = json_file_path,
          stan_model_file_path = stan_model_file_path))
  
}















#' R_fn_init_stan_inits_external
#' @keywords internal
#' @export
R_fn_init_stan_inits_external <- function( stan_data_list, 
                                           stan_model_file_path,
                                           stan_model_obj,
                                           stan_functions_directory,
                                           n_chains,
                                           init_lists_per_chain
) {
          ##
          ## First, check for duplicates in the data list:
          ##
          stan_data_list <- R_fn_remove_duplicates_from_list(stan_data_list)
          ##
          stan_data_list$priors <- NULL
          ##
          outs_init_stan_model <- init_stan_model_external( stan_data_list = stan_data_list,
                                                            stan_model_file_path = stan_model_file_path,
                                                            stan_model_obj = stan_model_obj,
                                                            n_chains = n_chains,
                                                            init_lists_per_chain = init_lists_per_chain)
          
          stan_init_pseudo_sampling_outs <- outs_init_stan_model$stan_init_pseudo_sampling_outs
          # stan_model_obj <- outs_init_stan_model$stan_model_obj
          json_file_path <- outs_init_stan_model$json_file_path
          stan_model_file_path <- outs_init_stan_model$stan_model_file_path
          ##
          ## Get variable names (from/using ".stan" file):
          ##
          {
            stan_param_names_list <- list()
            ##
            if (!is.null(stan_model_obj)) {
              ##
              stan_param_names_list$parameters <- names(stan_model_obj$variables()$parameters)
              ##
              stan_param_names_list$transformed_parameters <- names(stan_model_obj$variables()$transformed_parameters)
              ##
              stan_param_names_list$generated_quantities <- names(stan_model_obj$variables()$generated_quantities)
              ##
            } else {
              ##
              ## (no cmdstanr model, e.g. a NicoStan fit: the same names from stanc --info, which is what
              ##  cmdstanr's $variables() reads)
              stanc_file_path <- file.path(cmdstanr_path(), "bin", "stanc")
              info_of_the_model <- jsonlite::fromJSON(system2(stanc_file_path,
                                                              c("--info",
                                                                paste0("--include-paths=", stan_functions_directory),
                                                                stan_model_file_path),
                                                              stdout = TRUE,
                                                              stderr = FALSE))
              ##
              stan_param_names_list$parameters <- names(info_of_the_model[["parameters"]])
              ##
              stan_param_names_list$transformed_parameters <- names(info_of_the_model[["transformed parameters"]])
              ##
              stan_param_names_list$generated_quantities <- names(info_of_the_model[["generated quantities"]])
              ##
            }
            ##
            stan_param_names_list$all_params <- c(
              stan_param_names_list$parameters,
              stan_param_names_list$transformed_parameters,
              stan_param_names_list$generated_quantities
            )
          }
          ##
          ## Get "parameters" block param names (will be in the correct order):
          ##
          stan_param_names_main <- stan_param_names_list$parameters ; stan_param_names_main
          ##
          ## This legacy field is retained for compatibility. The vectors were not
          ## used by sampling or summaries, and attempting to infer Stan's
          ## unconstrained representation from JSON was not reliable.
          inits_unconstrained_vec_per_chain <- NULL
          ##
          ## ---- Output:
          ##
          return(list(  inits_unconstrained_vec_per_chain = inits_unconstrained_vec_per_chain,
                        stan_param_names_list = stan_param_names_list,
                        stan_param_names_main = stan_param_names_main,
                        stan_init_pseudo_sampling_outs = stan_init_pseudo_sampling_outs,
                        stan_model_obj = stan_model_obj,
                        json_file_path = json_file_path,
                        stan_model_file_path = stan_model_file_path))

}











#' init_stan_model_internal
#' @keywords internal
#' @export
init_stan_model_internal <- function( stan_data_list, 
                                      stan_model_name,
                                      stan_model_stan_include_paths = NULL,
                                      n_chains,
                                      init_lists_per_chain
) {
  
            ##
            ## Get package directory paths:
            ##
            outs <- get_MetaOrdDTA_stan_paths()
            pkg_root_dir <- outs$pkg_root_dir
            pkg_data_dir <- outs$pkg_data_dir
            pkg_stan_dir <- outs$pkg_stan_dir
            ##
            ## Get Stan model path:
            ##
            stan_model_file_path <- file.path(pkg_stan_dir, stan_model_name)
            ## Compile if necessary:
            if (is.null(stan_model_stan_include_paths)) { 
              stan_model_obj <- cmdstanr::cmdstan_model(stan_model_file_path)
            } else { 
              stan_model_obj <- cmdstanr::cmdstan_model(stan_model_file_path,
                                                        include_paths = stan_model_stan_include_paths) ## , force_recompile = TRUE)
            }
            
            outs_stan_model <- init_stan_model_external(  stan_data_list = stan_data_list, 
                                                          stan_model_file_path = stan_model_file_path,
                                                          stan_model_obj = stan_model_obj,
                                                          n_chains = n_chains,
                                                          init_lists_per_chain = init_lists_per_chain)
            
            return(outs_stan_model)
  
}

















#' R_fn_init_stan_inits_internal
#' @keywords internal
#' @export
R_fn_init_stan_inits_internal <- function( stan_data_list, 
                                           stan_model_name,
                                           stan_model_stan_include_paths = NULL,
                                           n_chains,
                                           init_lists_per_chain
) {
  
          outs_init_stan_model <- init_stan_model_internal( stan_data_list = stan_data_list,
                                                            stan_model_name = stan_model_name,
                                                            stan_model_stan_include_paths = stan_model_stan_include_paths,
                                                            n_chains = n_chains,
                                                            init_lists_per_chain = init_lists_per_chain)
          ##
          # require(bridgestan)
          # ?bridgestan::StanModel
          ##
          # stan_model_file_path$sample()
          # stan_model_file_path
          ##
          stan_init_pseudo_sampling_outs <- outs_init_stan_model$stan_init_pseudo_sampling_outs
          stan_model_obj <- outs_init_stan_model$stan_model_obj
          json_file_path <- outs_init_stan_model$json_file_path
          stan_model_file_path <- outs_init_stan_model$stan_model_file_path
          ##
          ## Get variable names:
          ##
          {
                stan_param_names_list <- list()
                ##
                message(paste("parameters = "))
                stan_param_names_list$parameters             <- print(names(stan_model_obj$variables()$parameters)) ## length = "n_params_main"
                ##
                message(paste("transformed_parameters = "))
                stan_param_names_list$transformed_parameters <- print(names(stan_model_obj$variables()$transformed_parameters))
                ##
                message(paste("generated_quantities = "))
                stan_param_names_list$generated_quantities   <- print(names(stan_model_obj$variables()$generated_quantities))
                ##
                message(paste("All Stan model variables = "))
                stan_param_names_list$all_params <- print(c( stan_param_names_list$parameters, 
                                                             stan_param_names_list$transformed_parameters, 
                                                        stan_param_names_list$generated_quantities))
          }
          ##
          ## Get "parameters" block param names (will be in the correct order):
          ##
          stan_param_names_main <- stan_param_names_list$parameters ; stan_param_names_main
          ##
          ## Set init vectors (in correct order)
          ##
          inits_unconstrained_vec_per_chain <- list()
          ##
          for (kk in 1:n_chains) {
              ##
              ## Get/format json string:
              ##
              json_string_for_inits_chain_kk <- convert_stan_data_list_to_JSON(init_lists_per_chain[[kk]])
              validated_json_string          <- paste(readLines(json_string_for_inits_chain_kk), collapse="")
              ##
              ## Turn the json string (of unconstrained params or "parameters block" params) into a vector of unconstrained 
              ## params or "parameters block" params:
              ##
              inits_unconstrained_vec_per_chain[[kk]] <- convert_JSON_string_to_ordered_R_vector(validated_json_string, stan_param_names_main)
          }
          
          
          return(list(  inits_unconstrained_vec_per_chain = inits_unconstrained_vec_per_chain,
                        stan_param_names_list = stan_param_names_list,
                        stan_param_names_main = stan_param_names_main,
                        stan_init_pseudo_sampling_outs = stan_init_pseudo_sampling_outs,
                        stan_model_obj = stan_model_obj,
                        json_file_path = json_file_path,
                        stan_model_file_path = stan_model_file_path))
  
}












