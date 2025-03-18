
#source(file.path(pkg_dir, "R/R_fn_EHMC_burnin_v2.R")) ### load R burnin-fn 





# -------- -- R function to prepare either hard-coded or user-supplied Stan mode, makes and prepares lists for the C++ structs  ----------------------------------------------------------------------------- 

#' init_model
#' @keywords internal
#' @export
init_model <- function(Model_type,
                       y,
                       X = NULL,
                       n_params_main = NULL,
                       model_args_list,
                       Stan_model_file_path,
                       Stan_cpp_user_header,
                       Stan_data_list
) {
  
  
  
  models_vec <- c("meta_ord", "meta_cts", "network_meta_ord", "network_meta_cts")
  
  if (!(Model_type %in% models_vec)) { 
    stop("Model_type must be set and be one of the following: 'meta_ord' 'meta_cts', 'network_meta_ord', 'network_meta_cts'")
  }
  
  if (is.null(y)) { 
    stop("ERROR: y (data) not inputted")
  }
  
  #### load hard-coded model args 
  if (Model_type  %in% models_vec) {
    
             outs_init_hard_coded_model <- init_hard_coded_model( Model_type = Model_type,
                                                                  y = y,
                                                                  N = N,
                                                                  model_args_list = model_args_list,
                                                                  ...)
             
             #### load fn args from list
             n_nuisance <- outs_init_hard_coded_model$n_nuisance
             n_params_main <- outs_init_hard_coded_model$n_params_main
             n_params <- outs_init_hard_coded_model$n_params
             ####
             n_covariates_per_outcome_mat <- outs_init_hard_coded_model$n_covariates_per_outcome_mat
             prior_only <- outs_init_hard_coded_model$prior_only
             prior_coeffs_mean_mat <- outs_init_hard_coded_model$prior_coeffs_mean_mat
             prior_coeffs_sd_mat <- outs_init_hard_coded_model$prior_coeffs_sd_mat
             corr_force_positive <- outs_init_hard_coded_model$corr_force_positive
             corr_param <- outs_init_hard_coded_model$corr_param
             lkj_cholesky_eta <- outs_init_hard_coded_model$lkj_cholesky_eta
             corr_prior_norm <- outs_init_hard_coded_model$corr_prior_norm
             corr_prior_beta <- outs_init_hard_coded_model$corr_prior_beta
             ub_corr <- outs_init_hard_coded_model$ub_corr
             lb_corr <- outs_init_hard_coded_model$lb_corr
             prev_prior_a <- outs_init_hard_coded_model$prev_prior_a
             prev_prior_b <- outs_init_hard_coded_model$prev_prior_b
             prior_for_skewed_LKJ_a <- outs_init_hard_coded_model$prior_for_skewed_LKJ_a
             prior_for_skewed_LKJ_b <- outs_init_hard_coded_model$prior_for_skewed_LKJ_b
             prior_for_corr_a <- outs_init_hard_coded_model$prior_for_corr_a
             prior_for_corr_b <- outs_init_hard_coded_model$prior_for_corr_b
             known_values_indicator_list <- outs_init_hard_coded_model$known_values_indicator_list
             known_values_list <- outs_init_hard_coded_model$known_values_list
             LT_b_priors_shape <- outs_init_hard_coded_model$LT_b_priors_shape
             LT_b_priors_scale <- outs_init_hard_coded_model$LT_b_priors_scale
             LT_known_bs_values <- outs_init_hard_coded_model$LT_known_bs_values
             LT_known_bs_indicator <- outs_init_hard_coded_model$LT_known_bs_indicator
             init_coeffs <- outs_init_hard_coded_model$init_coeffs
             init_raw_corrs <- outs_init_hard_coded_model$init_raw_corrs
             init_prev <- outs_init_hard_coded_model$init_prev
             X <- outs_init_hard_coded_model$X
             
             model_args_list <- outs_init_hard_coded_model
             
             ## load other fn args needed from outs_init_hard_coded_model
             # .... code to load other fn args here 
             Model_args_as_Rcpp_List = outs_init_hard_coded_model$Model_args_as_Rcpp_List

             Stan_model <- NULL
    
             json_file_path <- "none"
           
             bs_model <- NULL
             
             
  } else if (Model_type == "Stan") { 

            Stan_model <- NULL
            Stan_cpp_user_header = NULL
            Stan_cpp_flags = NULL
            
            bs_model_outs <- init_bs_model_external(Stan_data_list = Stan_data_list, 
                                                    Stan_model_file_path = Stan_model_file_path)
            
            ##
            bs_model <- bs_model_outs$bs_model
            print(paste("bs_model = "))
            print(bs_model)
            ##
            json_file_path <- bs_model_outs$json_file_path
            print(paste("json_file_path = "))
            print(json_file_path)
            
            ##
            Model_args_as_Rcpp_List <- list()
            Model_args_as_Rcpp_List$N <- N
            Model_args_as_Rcpp_List$n_params_main <- n_params_main
            Model_args_as_Rcpp_List$n_nuisance <- n_nuisance
            Model_args_as_Rcpp_List$json_file_path <- json_file_path
          
    
  }
  
          
          return(list(  bs_model = bs_model,
                        model_args_list = model_args_list,
                        Model_args_as_Rcpp_List = Model_args_as_Rcpp_List, 
                        Model_type = Model_type, 
                        y = y, 
                        X = X, 
                        Stan_model = Stan_model,
                        Stan_data_list = Stan_data_list, 
                        json_file_path = json_file_path))


}





