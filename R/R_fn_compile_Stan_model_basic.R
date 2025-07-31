




#' R_fn_compile_stan_model_basic_given_file_name
#' @keywords internal
#' @export
R_fn_compile_stan_model_basic_given_file_name <- function( stan_model_file_name,
                                                           cts,
                                                           network,
                                                           prior_only,
                                                           debugging,
                                                           force_recompile,
                                                           quiet,
                                                           compile
)  {
  
  
          stan_model_obj <- NULL
          
          if (debugging) {
            
                    os <- .Platform$OS.type
            
                    if (os == "unix") { 
                      user_root_dir <- Sys.getenv("PWD")
                    } else if (os == "windows") { 
                      user_root_dir <- Sys.getenv("USERPROFILE")
                    }
                    local_pkg_dir <- file.path(user_root_dir, "Documents/Work/PhD_work/R_packages/MetaOrdDTA")
                    
                    pkg_root_directory       <- local_pkg_dir
                    stan_models_directory    <- file.path(pkg_root_directory, "inst", "stan_models")
                    stan_functions_directory <- file.path(pkg_root_directory, "inst", "stan_functions")
                    ##
                    stan_MA_directory        <- file.path(stan_models_directory, "stan_models_MA")
                    stan_MA_prior_directory  <- file.path(stan_MA_directory, "prior_only")
                    ##
                    stan_NMA_directory       <- file.path(stan_models_directory, "stan_models_NMA")
                    stan_NMA_prior_directory <- file.path(stan_NMA_directory, "prior_only")
                    ##
                    if (network == TRUE) {
                      if (prior_only) stan_directory_to_use <-  stan_NMA_prior_directory
                      else            stan_directory_to_use <-  stan_NMA_directory
                    } else {
                      if (prior_only) stan_directory_to_use <-  stan_MA_prior_directory
                      else            stan_directory_to_use <-  stan_MA_directory
                    }
                    
                    stan_model_file_path <-  file.path(stan_directory_to_use, stan_model_file_name)
            
          } else {
                    
                    pkg_root_directory       <- system.file(package = "MetaOrdDTA")
                    ##
                    stan_models_directory    <- file.path(pkg_root_directory, "stan_models")
                    stan_functions_directory <- file.path(pkg_root_directory, "stan_functions")
                    ##
                    stan_MA_directory        <- file.path(stan_models_directory, "stan_models_MA")
                    stan_MA_prior_directory  <- file.path(stan_MA_directory, "prior_only")
                    ##
                    stan_NMA_directory       <- file.path(stan_models_directory, "stan_models_NMA")
                    stan_NMA_prior_directory <- file.path(stan_NMA_directory, "prior_only")
                    ##
                    if (network == TRUE) {
                      if (prior_only) stan_directory_to_use <-  stan_NMA_prior_directory
                      else            stan_directory_to_use <-  stan_NMA_directory
                    } else {
                      if (prior_only) stan_directory_to_use <-  stan_MA_prior_directory
                      else            stan_directory_to_use <-  stan_MA_directory
                    }
                    
                    stan_model_file_path <-  file.path(stan_directory_to_use, stan_model_file_name)
                    
          }
          ##
          ## Compile model:
          ##
          try({ 
                if (compile) {
                    
                    try({  
                      
                        stan_model_obj <- cmdstanr::cmdstan_model( stan_model_file_path,
                                                                   force_recompile = force_recompile,
                                                                   quiet = quiet, 
                                                                   include_paths = stan_functions_directory)
                      
                    })
                }
          })
          
          ##
          ## output list:
          ##
          return(list( stan_model_obj = stan_model_obj, 
                       stan_model_file_name = stan_model_file_name,
                       stan_model_file_path = stan_model_file_path,
                       ##
                       pkg_root_directory = pkg_root_directory,
                       stan_models_directory = stan_models_directory,
                       stan_functions_directory = stan_functions_directory,
                       ##
                       stan_MA_directory = stan_MA_directory,
                       stan_MA_prior_directory = stan_MA_prior_directory,
                       ##
                       stan_NMA_directory = stan_NMA_directory,
                       stan_NMA_prior_directory = stan_NMA_prior_directory))
  
}

















