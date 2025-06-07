




#' prep_data_and_model
#' @keywords internal
#' @export
prep_data_and_model <- function(  debugging,
                                  ##
                                  x, 
                                  ##
                                  X,
                                  ##
                                  # n_index_tests_per_study = NULL, ## only needed for NMA
                                  indicator_index_test_in_study, ## only needed for NMA
                                  ##
                                  internal_obj,
                                  ##
                                  basic_model_options,
                                  advanced_model_options,
                                  MCMC_params,
                                  other_advanced_options,
                                  ##
                                  priors,
                                  ##
                                  init_lists_per_chain,
                                  ##
                                  compute_sim_study_metrics,
                                  vec_index_inner_thr
) {
          
          if (is.null(other_advanced_options$advanced_compile)) { 
            other_advanced_options$advanced_compile <- FALSE
          }
          ##
          ## ---- Unpack elements from grouped lists:
          ##
          default_fns_list <- list(  network = MetaOrdDTA:::default_network,
                                     cts = MetaOrdDTA:::default_cts,
                                     model_parameterisation = MetaOrdDTA:::default_model_parameterisation,
                                     random_thresholds = MetaOrdDTA:::default_random_thresholds,
                                     Dirichlet_random_effects_type = MetaOrdDTA:::default_Dirichlet_random_effects_type,
                                     ##
                                     box_cox = MetaOrdDTA:::default_box_cox,
                                     ##
                                     softplus = MetaOrdDTA:::default_softplus)
          
          ##
          default_n_chains <- round(parallel::detectCores()/2)
          MCMC_params$n_chains <- MetaOrdDTA:::if_null_then_set_to(MCMC_params$n_chains, default_n_chains)
          print(paste("MCMC_params$n_chains = ", MCMC_params$n_chains))
          ##
          try({  
              ##
              # ## Then we can set the "cts" binary indicator:
              # if (test_type %in% c("cts", "continuous")) { 
              #   cts <- TRUE
              # } else { 
              #   cts <- FALSE
              # }
              basic_model_options$cts <-  MetaOrdDTA:::if_null_then_set_to( basic_model_options$cts,  
                                                                           default_fns_list$cts()[[1]])
              print(paste("basic_model_options$cts = ", basic_model_options$cts))
              ##
              ## Check if NMA or not (and set "network" binary indicator):
              ##
              basic_model_options$network <-  MetaOrdDTA:::if_null_then_set_to( basic_model_options$network,  
                                                                               default_fns_list$network(x)[[1]])
              print(paste("basic_model_options$network = ", basic_model_options$network))
              ##
              ## Check if prior-only model:
              ##
              basic_model_options$prior_only <-  MetaOrdDTA:::if_null_then_set_to( basic_model_options$prior_only, 
                                                                                  FALSE)
              print(paste("basic_model_options$prior_only = ", basic_model_options$prior_only))
              ##
              ## "Advanced" options for ordinal models ONLY:
              ##
              advanced_model_options$model_parameterisation        <-  MetaOrdDTA:::if_null_then_set_to( advanced_model_options$model_parameterisation,  
                                                                                                        default_fns_list$model_parameterisation(
                                                                                                            cts = basic_model_options$cts)[[1]])
              print(paste("advanced_model_options$model_parameterisation = ", advanced_model_options$model_parameterisation))
              ##
              ##
              ##
              advanced_model_options$random_thresholds             <-  MetaOrdDTA:::if_null_then_set_to( advanced_model_options$random_thresholds, 
                                                                                                        default_fns_list$random_thresholds(
                                                                                                            cts = basic_model_options$cts)[[1]])
              print(paste("advanced_model_options$random_thresholds = ", advanced_model_options$random_thresholds))
              ##
              ##
              ##
              advanced_model_options$Dirichlet_random_effects_type <-  MetaOrdDTA:::if_null_then_set_to( advanced_model_options$Dirichlet_random_effects_type, 
                                                                                              default_fns_list$Dirichlet_random_effects_type(
                                                                                              cts = basic_model_options$cts, 
                                                                                              random_thresholds = advanced_model_options$random_thresholds)[[1]])
              print(paste("advanced_model_options$Dirichlet_random_effects_type = ", advanced_model_options$Dirichlet_random_effects_type))
              ##
              ##
              ##
              advanced_model_options$box_cox <-  MetaOrdDTA:::if_null_then_set_to( advanced_model_options$box_cox, 
                                                                                  default_fns_list$box_cox(cts = basic_model_options$cts)[[1]])
              print(paste("advanced_model_options$box_cox = ", advanced_model_options$box_cox))
              ##
              ##
              ##
              advanced_model_options$softplus <-  MetaOrdDTA:::if_null_then_set_to( advanced_model_options$softplus, 
                                                                                   default_fns_list$softplus()[[1]])
              print(paste("advanced_model_options$softplus = ", advanced_model_options$softplus))
          })
          ##
          ## -----------------  Setup data: ---------------------------------------------------------------------------
          ##
          if (basic_model_options$network) {
              try({ 
                n_index_tests_per_study <- rowSums(indicator_index_test_in_study)
              }, silent = TRUE)
              print(paste("n_index_tests_per_study = "))
              cat(n_index_tests_per_study)
          } else { 
              n_index_tests_per_study <- NULL
          }
          
          {
                 if (basic_model_options$network)  { 
                      data_fn <- MetaOrdDTA:::R_fn_prep_NMA_data
                 } else { 
                      data_fn <- MetaOrdDTA:::R_fn_prep_MA_data
                 }
              
                 outs_data <- data_fn( x = x,
                                       ##
                                       X = X,
                                       model_parameterisation = model_parameterisation,
                                       ##
                                       n_index_tests_per_study = n_index_tests_per_study,
                                       indicator_index_test_in_study = indicator_index_test_in_study
                                       ##
                                       # compute_sim_study_metrics = compute_sim_study_metrics,
                                       # vec_index_inner_thr = vec_index_inner_thr
                                       )
          }
           ##
          internal_obj$outs_data <- outs_data
          # ##
          # print(paste("internal_obj$outs_data$stan_data_list = "))
          # print(str(internal_obj$outs_data$stan_data_list))
          ##
          if (basic_model_options$network == FALSE) { 
                 internal_obj$outs_data$n_tests <- 1
          }
          # internal_obj$outs_data$stan_data_list
          ##
          ## ----- Get the Stan model file (.stan file) name: ----------------------------------------------------------
          ##
          internal_obj$outs_stan_model_name  <- MetaOrdDTA:::R_fn_get_stan_model_file_name(
                                                                custom_file_name = advanced_model_options$custom_file_name,
                                                                ##
                                                                network    = basic_model_options$network, 
                                                                cts        = basic_model_options$cts, 
                                                                prior_only = basic_model_options$prior_only, 
                                                                ##
                                                                model_parameterisation        = advanced_model_options$model_parameterisation,
                                                                random_thresholds             = advanced_model_options$random_thresholds,
                                                                Dirichlet_random_effects_type = advanced_model_options$Dirichlet_random_effects_type)
          ##
          message(paste("internal_obj$outs_stan_model_name$stan_model_file_name = "))
          print(internal_obj$outs_stan_model_name)
          ##
          ## ----------------- Compile Stan model (if necessary): ------------------------------------------------------
          ## 
          ## R_fn_compile_stan_advanced_given_file_name
          ## R_fn_compile_stan_model_basic_given_file_name
          if (is.null(other_advanced_options$advanced_compile)) { 
            other_advanced_options$advanced_compile <- FALSE
          }
          
          if (other_advanced_options$advanced_compile) {
                  
                  other_advanced_options$force_recompile <- MetaOrdDTA:::if_null_then_set_to(other_advanced_options$force_recompile, TRUE)
                  other_advanced_options$quiet <- MetaOrdDTA:::if_null_then_set_to(other_advanced_options$quiet, FALSE)
                  ##
                  set_custom_CXX_CPP_flags <- other_advanced_options$set_custom_CXX_CPP_flags
                  CCACHE_PATH <- other_advanced_options$CCACHE_PATH
                  custom_cpp_user_header_file_path <- other_advanced_options$custom_cpp_user_header_file_path
                  CXX_COMPILER_PATH <- other_advanced_options$CXX_COMPILER_PATH
                  CPP_COMPILER_PATH <- other_advanced_options$CPP_COMPILER_PATH
                  MATH_FLAGS <- other_advanced_options$MATH_FLAGS
                  FMA_FLAGS <- other_advanced_options$FMA_FLAGS
                  AVX_FLAGS <- other_advanced_options$AVX_FLAGS
                  THREAD_FLAGS <- other_advanced_options$THREAD_FLAGS
 
                  internal_obj$outs_stan_compile <- MetaOrdDTA:::R_fn_compile_stan_advanced_given_file_name( 
                                                                                 stan_model_file_name = internal_obj$outs_stan_model_name$stan_model_file_name,
                                                                                 ##
                                                                                 cts        = basic_model_options$cts, 
                                                                                 network    = basic_model_options$network, 
                                                                                 prior_only = basic_model_options$prior_only,
                                                                                 ##
                                                                                 debugging = debugging,
                                                                                 force_recompile = other_advanced_options$force_recompile,
                                                                                 quiet = other_advanced_options$quiet,
                                                                                 ##
                                                                                 set_custom_CXX_CPP_flags = set_custom_CXX_CPP_flags, 
                                                                                 CCACHE_PATH = " ", 
                                                                                 custom_cpp_user_header_file_path = custom_cpp_user_header_file_path, ## if wish to include custom C++ files
                                                                                 CXX_COMPILER_PATH = CXX_COMPILER_PATH,
                                                                                 CPP_COMPILER_PATH = CPP_COMPILER_PATH,
                                                                                 MATH_FLAGS = MATH_FLAGS,
                                                                                 FMA_FLAGS = FMA_FLAGS,
                                                                                 AVX_FLAGS = AVX_FLAGS,
                                                                                 THREAD_FLAGS = THREAD_FLAGS)
          } else { 
                  
                  internal_obj$outs_stan_compile <- MetaOrdDTA:::R_fn_compile_stan_model_basic_given_file_name( 
                                                                                stan_model_file_name = internal_obj$outs_stan_model_name$stan_model_file_name,
                                                                                ##
                                                                                cts        = basic_model_options$cts, 
                                                                                network    = basic_model_options$network, 
                                                                                prior_only = basic_model_options$prior_only,
                                                                                ##
                                                                                debugging = debugging,
                                                                                force_recompile = FALSE,
                                                                                quiet   = FALSE,
                                                                                compile = TRUE)
                                                                            
            
          }
          
          # ##
          # message(paste("internal_obj$outs_stan_compile = "))
          # print(internal_obj$outs_stan_compile)
          ##
          ## -----------------  Set priors (using R package defaults): --------------------------------------------------
          ##
          priors <- if_null_then_set_to(x = priors, list())
          ##
          if (basic_model_options$network)  { 
            prior_fn <- R_fn_set_priors_NMA
          } else { 
            prior_fn <- R_fn_set_priors_MA
          }
          ##
          if (!(basic_model_options$network)) { 
            internal_obj$outs_data$n_index_tests <- 1
          }
          ##
          priors <- prior_fn(  
                             priors = priors,
                             ##
                             X = X,
                             ##
                             n_index_tests = internal_obj$outs_data$n_index_tests,
                             n_studies = internal_obj$outs_data$n_studies,
                             ##
                             cts     = basic_model_options$cts,
                             ##
                             model_parameterisation        = advanced_model_options$model_parameterisation,
                             random_thresholds             = advanced_model_options$random_thresholds,
                             Dirichlet_random_effects_type = advanced_model_options$Dirichlet_random_effects_type,
                             softplus                      = advanced_model_options$softplus,
                             ##
                             n_cat = internal_obj$outs_data$n_cat,
                             n_thr = internal_obj$outs_data$n_thr)
          ##
          # message(paste("priors = "))
          # print(priors)
          ##
          ## ----------------- Set initial values (using R package defaults) - AFTER priors have been set!: -------------
          ##
          if (MCMC_params$n_chains < 1) {  MCMC_params$n_chains <- round(parallel::detectCores()/2) }
          # try({ 
          #   MCMC_params$n_chains <- MetaOrdDTA::if_null_then_set_to( MCMC_params$n_chains, length(init_lists_per_chain))
          # })
          # try({ 
          #   MCMC_params$n_chains <- MetaOrdDTA::if_null_then_set_to( MCMC_params$n_chains, round(parallel::detectCores()/2))
          # })
          ##
          # # init_lists_per_chain = NULL
          # print(paste("n_chains = ", MCMC_params$n_chains))
          # ##
          # try({ 
          #     print(paste("init_lists_per_chain[[1]] = "))
          #     print(init_lists_per_chain)
          # })
          # ##
          try({
            if (is.null(init_lists_per_chain)) {
                  init_lists_per_chain <- list()
                  for (i in 1:MCMC_params$n_chains) {
                    init_lists_per_chain[[i]] <- list(one = 1) ## placeholder (NULL doesnt work here)
                  }
            }

            # else {
            #       if (length(init_lists_per_chain) != MCMC_params$n_chains) {
            #         stop("init_lists_per_chain must either be set to 'NULL' (default inits)
            #                      or be a list of lists of length 'MCMC_params$n_chains")
            #       }
            # }
          })
          ##
          # print(paste("advanced_model_options = "))
          # print(advanced_model_options)
          ##
          if (basic_model_options$network)  { 
            inits_fn <- R_fn_set_inits_NMA
          } else { 
            inits_fn <- R_fn_set_inits_MA
          }
          ##
          try({ 
            for (i in 1:MCMC_params$n_chains) {
                 init_lists_per_chain[[i]] <- inits_fn( 
                                                       inits  = init_lists_per_chain[[i]],
                                                       ##
                                                       priors = priors,
                                                       ##
                                                       X = X,
                                                       ##
                                                       cts    = basic_model_options$cts,
                                                       ##
                                                       n_index_tests = internal_obj$outs_data$n_index_tests,
                                                       n_studies     = internal_obj$outs_data$n_studies,
                                                       n_thr         = internal_obj$outs_data$n_thr,
                                                       ##
                                                       model_parameterisation        = advanced_model_options$model_parameterisation,
                                                       random_thresholds             = advanced_model_options$random_thresholds,
                                                       Dirichlet_random_effects_type = advanced_model_options$Dirichlet_random_effects_type,
                                                       softplus                      = advanced_model_options$softplus)
            }
          })
          # ##
          # message(paste("init_lists_per_chain[[1]] = "))
          # print(init_lists_per_chain[[1]])
          # init_lists_per_chain
          # init_lists_per_chain[[1]]$inits$
          # ##
          # ## ---- Pack the modified variables back into the lists:
          # ##
          ##
          ## ----------------- Output list: -----------------------------------------------------------------------------
          ##
          return(list(internal_obj = internal_obj,
                      ##
                      basic_model_options = basic_model_options,
                      advanced_model_options = advanced_model_options,
                      MCMC_params = MCMC_params,
                      other_advanced_options = other_advanced_options,
                      ##
                      init_lists_per_chain = init_lists_per_chain,
                      ##
                      priors = priors,
                      ##
                      X = X))
          
}


 



