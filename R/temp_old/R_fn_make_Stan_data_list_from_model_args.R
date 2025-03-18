


 




#' R_fn_make_Stan_data_list_from_model_args
#' @keywords internal
#' @export
R_fn_make_Stan_data_list_from_model_args <- function(debugging = FALSE,
                                                     # initial_object,
                                                     Stan_data_list,
                                                     priors,
                                                     prior_only) {
  
                try({  rm(Stan_data_list$priors) }, silent = TRUE)
                Stan_data_list <- c(Stan_data_list, priors)
                
                Stan_data_list$x_with_missings <- Stan_data_list$x
                
                try({  rm(Stan_data_list$test_type) })
                
                return(Stan_data_list)
  
}






 





