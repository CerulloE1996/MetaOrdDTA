# source("R_fn_init_model_and_vals.R")

update_model    <- function(init_object, ### this is the object to update 
                            Model_type,
                            y,
                            N,
                            init_lists_per_chain,
                            sample_nuisance,
                            model_args_list,
                            Stan_data_list,
                            n_params_main,
                            n_nuisance,
                            n_chains_burnin,
                            force_recompile,
                            ...) {
  
 
  # ### old arguments 
    Model_type_old <- init_object$Model_type
  # 
  # y_old <- init_object$y
  #  
  # init_lists_per_chain_old <- init_object$init_lists_per_chain
  # 
  # Stan_data_list_old <- init_object$Stan_data_list
  # model_args_list_old <- init_object$model_args_list
 
  
  ### Model_type (can't be updated here)
  if (!(is.null(Model_type))) { 
    if (Model_type != Model_type_old) { 
      stop("To update Model_type, the model needs to be re-initialised (i.e. please run the $new() method again)")
    }
  }
  
  ### update dataset (y)
  if (!(is.null(y))) {
   # message("Updating data (y)")
    init_object$y <- y
    warning("If updating data, if appropriate for your model, be sure to also update n_nuisance")
  } else {  # don't update
    y <-  init_object$y
  }
  
  ### update initial values per chain (init_lists_per_chain)
  if (!(is.null(init_lists_per_chain))) {
   # message("Updating initial values per chain (init_lists_per_chain)")
    if (  length(init_lists_per_chain) != n_chains_burnin) { 
      warning("The length of init_lists_per_chain must be equal to n_chains_burnin")
      message("Please ensure that the length of init_lists_per_chain = n_chains_burnin")
    }
    init_object$init_lists_per_chain <- init_lists_per_chain
  } else {  # don't update
    init_lists_per_chain <-  init_object$init_lists_per_chain
  }
  
  ### update sample_nuisance
  if (!(is.null(sample_nuisance))) {
    init_object$sample_nuisance <- sample_nuisance
  } else {  # don't update
     sample_nuisance <-  init_object$sample_nuisance
  }
  
  ### update model_args_list
  if (!(is.null(model_args_list))) {
    init_object$model_args_list <- model_args_list
  } else {  # don't update
    model_args_list <-  init_object$model_args_list
  }
  
  
  ### update Stan_data_list
  if (!(is.null(Stan_data_list))) {
    init_object$Stan_data_list <- Stan_data_list
  } else {  # don't update
    Stan_data_list <-  init_object$Stan_data_list
  }
  
  ### update n_params_main
  if (!(is.null(n_params_main))) {
    init_object$n_params_main <- n_params_main
  } else {  # don't update
    n_params_main <-  init_object$n_params_main
  }
  
  ### update n_nuisance
  if (!(is.null(n_nuisance))) {
    init_object$n_nuisance <- n_nuisance
  } else {  # don't update
    n_nuisance <-  init_object$n_nuisance
  }
    

  ### now re-initialise but WITHOUT re-compiling the model
    init_object <- BayesMVP:::initialise_model(  Model_type = Model_type,
                                                 compile = FALSE, ## DONT re-compile model
                                                 force_recompile = force_recompile,
                                                 cmdstanr_model_fit_obj = init_object$cmdstanr_model_fit_obj,
                                                 y = y,
                                                 N = N,
                                                 init_lists_per_chain = init_lists_per_chain,
                                                 n_params_main = n_params_main,
                                                 n_nuisance = n_nuisance, 
                                                 sample_nuisance = sample_nuisance,
                                                 n_chains_burnin = n_chains_burnin,
                                                 model_args_list = model_args_list,
                                                 Stan_data_list = Stan_data_list,
                                                 Stan_model_file_path = init_object$Stan_model_file_path,
                                                 Stan_cpp_user_header = init_object$Stan_cpp_user_header,
                                                 Stan_cpp_flags = init_object$Stan_cpp_flags,
                                                 ...)
  
  message(  cat(colourise(     (paste(Model_type, "Model updated!"))         , "green"), "\n") )
  
  
  return(init_object)
  
  
}

















                  
                   
                  
                 
  