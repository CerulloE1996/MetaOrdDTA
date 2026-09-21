
#' @title MetaOrdDTA Model Class
#' @description R6 Class for MetaOrdDTA model initialization and sampling.
#'
#' @details 
#' This class handles:
#' * Model initialization and compilation
#' * MCMC sampling with NUTS-HMC (using Stan).
#' * Parameter updates and diagnostics
#'
#' @section Model Types:
#' * meta_ord: Meta-analysis model for a single ordinal test (e.g., questionnaires). Uses model proposed by Cerullo et al. 
#'   Supports missing thresholds in some or all studies.
#'   
#' * meta_cts: Meta-analysis model for a single continuous test (e.g., biomarkers). Uses the model proposed by Jones et al. 
#'   Supports missing thresholds in some or all studies.
#'   
#' * network_meta_ord: Meta-analysis model for multiple ordinal tests (e.g., questionnaires). Uses the model proposed by Cerullo et al. 
#'   and NMA model proposed by Nyaga et al. Supports missing thresholds in some or all studies.
#'   
#' * network_meta_cts: Meta-analysis model for multiple continuous tests (e.g., biomarkers). Uses the model proposed by Jones et al. and 
#'   Nyaga et al. Supports missing thresholds in some or all studies.
#' 
#' @section EXAMPLE - Typical workflow:
#'  \preformatted{
#'  #### NOTE: Please see the ".R" example files for full details & examples, 
#'  #### as this is NOT a complete working example, but is provided just to show 
#'  #### what order you should call things in. 
#'  ####  -----------  Compile + initialise the model using "MetaOrd_model$new(...)"  ---------------- 
#'    require(MetaOrdDTA)
#'    ## NOTE: The "model_args_list" argument is only needed for BUILT-IN (not Stan) models.
#'    ##
#'    ## --- INITIALISE MODEL:  --------------------------------------------------------------- 
#'    ##
#'    model_obj <- MetaOrdDTA::MetaOrd_model$new(   
#'                    Model_type = "meta_ord",
#'                    x = x, ## data
#'                    model_args_list = model_args_list,
#'                    n_chains = n_chains
#'                    init_lists_per_chain = init_lists_per_chain)
#'    ##
#'    ## --- SAMPLE MODEL: --------------------------------------------------------------------
#'    ##
#'    model_samples <-  model_obj$sample(
#'                    seed = 123,
#'                    n_burnin = 500,
#'                    n_iter = 1000,
#'                    x = x,
#'                    n_chains = 16,
#'                    init_lists_per_chain = init_lists_per_chain,
#'                    model_args_list = model_args_list)
#'    ##
#'    ## --- EXTRACT MODEL RESULTS SUMMARY + DIAGNOSTICS ---------------------------------------  
#'    ##         
#'    model_fit <- model_samples$summary(
#'                    save_log_lik_trace = FALSE, 
#'                    compute_nested_rhat = FALSE,
#'                    compute_transformed_parameters = FALSE) #'    
#'    
#'    ## extract # divergences + % of sampling iterations which have divergences
#'    model_fit$get_divergences()
#'    
#'    #### --- TRACE PLOTS  ----------------------------------------------------------------------
#'    ## trace_plots_all <- model_samples$plot_traces() # if want the trace for all parameters 
#'    trace_plots <- model_fit$plot_traces(
#'                    params = c("beta", "Omega", "p"), 
#'                    batch_size = 12)
#'    
#'    ## you can extract parameters by doing: "trace$param_name()". 
#'    ## For example:
#'    ## display each panel for beta and Omega ("batch_size" controls the # of plots per panel)
#'    trace_plots$beta[[1]] # 1st (and only) panel
#'    
#'    
#'    #### --- POSTERIOR DENSITY PLOTS -----------------------------------------------------------
#'    ## density_plots_all <- model_samples$plot_densities() # if want the densities 
#'    ## for all parameters.
#'    ## Let's plot the densities for: sensitivity, specificity, and prevalence 
#'    density_plots <- model_fit$plot_densities(
#'                      params = c("Se_bin", "Sp_bin", "p"), 
#'                      batch_size = 12)
#'    
#'    ## you can extract parameters by doing: "trace$param_name()". 
#'    ## For example:
#'    ## display each panel for beta and Omega ("batch_size" controls the # of plots per panel)
#'    density_plots$Se[[1]] # Se - 1st (and only) panel
#'    density_plots$Sp[[1]] # Sp - 1st (and only) panel
#'    density_plots$p[[1]] # p (prevelance) - 1st (and only) panel
#'    
#'    
#'    ###### --- EXTRACT PARAMETER SUMMARY TIBBLES: ----------------------------------------------
#'    ## The "model_summary" object (created using the "$summary()" method) contains
#'    ##  many useful objects. 
#'    ## For example:
#'    require(dplyr)
#'    ## nice summary tibble for main parameters, includes ESS/Rhat, etc:
#'    model_fit$get_summary_main() %>% print(n = 50) 
#'    ## nice summary tibble for transformed parameters, includes ESS/Rhat, etc:
#'    model_fit$get_summary_transformed() %>% print(n = 150) 
#'    ## nice summary tibble for generated quantities, includes 
#'    ## ESS/Rhat, etc (for LC-MVP this includes Se/Sp/prevalence):
#'    model_fit$get_summary_generated_quantities () %>% print(n = 150) 
#'    
#' }
#' 
#' 
#'
#' @export
MetaOrd_model <- R6::R6Class("MetaOrd_model",
                    
              public = list(
                    ##
                    ## -------------------------  Store all "important" parameters as class members / define them first:
                    ##
                    debugging = NULL,
                    ##
                    baseline_outs_MR = NULL,
                    pairwise_comparisons_MR = NULL,
                    ##
                    ## ----  Core / internal objects (NOT user-inputted)
                    ##
                    outs_data = list(
                        stan_data_list = NULL,
                        n_tests = NULL,
                        n_studies = NULL,
                        n_studies_per_test = NULL,
                        n_thr = NULL,
                        n_cat = NULL
                    ),
                    ##
                    outs_stan_model_name = list( 
                        stan_model_file_name = NULL
                    ),
                    ##
                    outs_stan_compile = list( 
                        stan_model_obj = NULL, 
                        #### stan_model_file_name = NULL,
                        stan_model_file_path = NULL,
                        ##
                        pkg_root_directory = NULL,
                        stan_models_directory = NULL,
                        stan_functions_directory = NULL,
                        ##
                        stan_MA_directory = NULL,
                        stan_MA_prior_directory = NULL,
                        ##
                        stan_NMA_directory = NULL,
                        stan_NMA_prior_directory = NULL
                    ),
                    ##
                    outs_stan_init = list( 
                        inits_unconstrained_vec_per_chain = NULL,
                        stan_param_names_list = NULL,
                        stan_param_names_main = NULL,
                        stan_init_pseudo_sampling_outs = NULL,
                        stan_model_obj = NULL,
                        json_file_path = NULL,
                        stan_model_file_path = NULL
                    ),
                    ##
                    outs_stan_sampling = list( 
                        stan_mod_samples = NULL,
                        time_total = NULL
                    ),
                    ##
                    ## ---- Main "internal_obj" list:
                    ##
                    internal_obj = list(
                        outs_data = NULL,
                        ##
                        outs_stan_model_name = NULL,
                        outs_stan_compile = NULL,
                        outs_stan_init = NULL,
                        outs_stan_sampling = NULL,
                        ##
                        HMC_info = NULL,
                        efficiency_info = NULL,
                        summaries = NULL,
                        traces = NULL
                    ),
                    ##
                    ## ---- Data (user-inputted - MANDATORY argument)
                    ##
                    x = NULL,
                    ##
                    ## ---- Data (user-inputted - MANDATORY argument)
                    ##
                    intercept_only = NULL,
                    cov_data = NULL,
                    ##
                    indicator_index_test_in_study = NULL,
                    ##
                    ## ---- Basic modelling options (user-inputted element-by-element)
                    ##
                    basic_model_options = list(
                        network = NULL,
                        cts = NULL,
                        prior_only = NULL
                    ),
                    #' @field network Whether to perform meta-analysis (MA) or network-meta-analysis (NMA). The default will depend on the structure of the data (i.e. on \code{x}).
                    #' @field cts Whether to use continuous model (i.e., Jones et al-based) or ordinal. Default is ordinal. Type of test (and hence the type of model to fit) -
                    #' can be either "continuous" (or simply "cts") or "ordinal" (or simply "ord"). The former uses the
                    #' model from Jones et al and assumes, that the test is continuous (e.g., biokarkers) and the latter assumes the data is ordinal (e.g., questionnaires). 
                    #' The latter is generally more flexible so if you are unsure which model type to use, we suggest using "ordinal". 
                    #' Furthermore, if performing network-meta-analysis (NMA) and you have a mix of continuous and ordinal tests, we also suggest using "ordinal". 
                    #' The default is "ordinal".
                    #' @field prior_only Whether to run prior-only (for prior predictive checking) or not. Default is FALSE. 
                    ##
                    ## ---- "advanced" user-inputted modelling options  (OPTIONAL - user-inputted element-by-element)
                    ##
                    advanced_model_options = list(
                        model_parameterisation = NULL,
                        ##
                        random_thresholds = NULL,
                        Dirichlet_random_effects_type = NULL,
                        ##
                        box_cox = NULL,
                        softplus = NULL,
                        custom_file_name = NULL,
                        ##
                        compound_symmetry = NULL,
                        hetero_sigma = NULL
                    ),
                    ##
                    ## ---- other "advanced" user-inputted options (e.g. C++ compiler flags - OPTIONAL - user-inputted element-by-element)
                    ##
                    other_advanced_options = list(
                        advanced_compile = NULL, ## default is false (aka a "basic" compile)
                        ##
                        force_recompile = NULL,
                        quiet = NULL,
                        compile = NULL,
                        ##
                        set_custom_CXX_CPP_flags = NULL,
                        CCACHE_PATH = NULL,
                        custom_cpp_user_header_file_path = NULL,
                        CXX_COMPILER_PATH = NULL,
                        CPP_COMPILER_PATH = NULL,
                        MATH_FLAGS = NULL,
                        FMA_FLAGS = NULL,
                        AVX_FLAGS = NULL,
                        THREAD_FLAGS = NULL
                    ),
                    ##
                    ## ---- priors  (OPTIONAL - user-inputted element-by-element)
                    ##
                    priors = NULL,
                    #' @field priors  An R list of priors. Default will depend on model.
                    ##
                    ## ---- initial values for each chain (OPTIONAL - user-inputted element-by-element)
                    ##
                    init_lists_per_chain = NULL,
                    #'@field init_lists_per_chain A list of dimension n_chains, where each element of the list is another list which contains the #
                    #' initial values for the chain. Note that this is the same format for initial values that Stan (both rstan and cmdstanr) uses. 
                    ##
                    ## ----  (OPTIONAL - user-inputted element-by-element)
                    ##
                    MCMC_params = list(
                      seed = NULL,
                      n_superchains = NULL,
                      n_chains = NULL,
                      n_iter = NULL,
                      n_burnin = NULL,
                      adapt_delta = NULL,
                      max_treedepth = NULL,
                      metric_shape = NULL
                    ),
                    ##
                    # cmdstanr_args = list(),
                    ##
                    ## ---- Internal use only (for sim study):
                    ##
                    compute_sim_study_metrics = NULL,
                    vec_index_inner_thr = NULL,
                    ##
                    ## ---- Other class members:
                    ##
                    ## ---------- constructor - initialize using the prep_data_and_model fn (this wraps the prep_data_and_model function with $new())
                    ## - store all important parameters:
                    #'@description Create a new "MetaOrd_model" object.
                    ##
                    #'@param x The dataset. See class documentation for details.
                    ##
                    #'@param intercept_only Optional covariates and baseline cases. See class documentation for details.
                    #'@param cov_data Optional covariates and baseline cases. See class documentation for details.
                    #'
                    ##
                    #'@param test_type Type of test. See class documentation for details.
                    #'@param network See class documentation for details.
                    #'@param cts Whether to use a continuous ("cts") model or whether to use an ordinal ("ord") model.
                    ##
                    #'@param prior_only See class documentation for details.
                    ##
                    #'@param priors See class documentation for details.
                    ##
                    #'@param model_parameterisation See class documentation for details.
                    #'@param random_thresholds See class documentation for details.
                    #'@param Dirichlet_random_effects_type See class documentation for details.
                    ##
                    #'@param box_cox See class documentation for details.
                    ##
                    #'@param softplus See class documentation for details.
                    ##
                    #'@param n_chains Number of chains used for burnin. See class documentation for details.
                    #'@param init_lists_per_chain List of initial values for each chain. See class documentation for details.
                    #'@return Returns self$model_initial_obj, an object generated from the "MetaOrdDTA::prep_data_and_model" function.
                    ##
                    initialize = function(  debugging = FALSE,
                                            ##
                                            x, ## Mandatory
                                            ##
                                            intercept_only, ## Mandatory
                                            cov_data = self$cov_data,
                                            ##
                                            indicator_index_test_in_study = self$indicator_index_test_in_study,
                                            ##
                                            network = NULL,
                                            compound_symmetry = NULL,
                                            hetero_sigma = NULL,
                                            ##
                                            prior_only = FALSE,
                                            ##
                                            ## "Advanced" options for ordinal models ONLY:
                                            ##
                                            model_parameterisation, ## Mandatory
                                            random_thresholds = NULL,
                                            Dirichlet_random_effects_type = NULL,
                                            ##
                                            ## "Advanced" options for cts (i.e., Jones) models ONLY:
                                            ##
                                            box_cox = NULL, ## only for "Jones" model
                                            ##
                                            ## "Advanced" options for ALL models:
                                            ##
                                            softplus = FALSE, ## exp() is default (for cutpoints and also scale parameters for "Jones" model)
                                            ##
                                            priors = NULL,
                                            ##
                                            n_chains = NULL,
                                            init_lists_per_chain = self$init_lists_per_chain,
                                            ##
                                            ## Other "advanced" options (non-modelling options e.g. compiler / C++ / Stan flags):
                                            ##
                                            advanced_compile = FALSE, ## default is standard/"basic" compile
                                            ##
                                            force_recompile = FALSE,
                                            quiet = FALSE,
                                            compile = TRUE,
                                            ##
                                            set_custom_CXX_CPP_flags = FALSE, 
                                            CCACHE_PATH = " ", 
                                            custom_cpp_user_header_file_path = NULL, ## if wish to include custom C++ files
                                            CXX_COMPILER_PATH = NULL,
                                            CPP_COMPILER_PATH = NULL,
                                            MATH_FLAGS = NULL,
                                            FMA_FLAGS = NULL,
                                            AVX_FLAGS = NULL,
                                            THREAD_FLAGS = NULL,
                                            ##
                                            custom_file_name = NULL
                                            ) {
                      
                              # message(("aaa_class_hello_1"))
                              ##
                              if (model_parameterisation %in% c("Jones", "jones")) {
                                
                                      cts <- TRUE
                                      random_thresholds <- FALSE
                                      Dirichlet_random_effects_type <- "none"
                                      ##
                                      if (is.null(box_cox)) { 
                                        warning("box_cox not specificed - setting to TRUE")
                                        box_cox <- TRUE
                                      }
                                    
                              } else  {
                                
                                    cts <- FALSE
                                    box_cox <- FALSE
                                    ##
                                    if (is.null(random_thresholds)) { 
                                      warning("random_thresholds not specificed - setting to FALSE")
                                      random_thresholds <- FALSE
                                    }
                                    ##
                                    if (random_thresholds == TRUE) {
                                      if (is.null(Dirichlet_random_effects_type)) { 
                                        warning("Dirichlet_random_effects_type not specificed - setting to 'kappa' (recommended)")
                                        Dirichlet_random_effects_type <- "kappa"
                                      }
                                    }
                                    ##
                                    if (grepl("ord_biv|ordinal_biv|biv|Xu|O-", model_parameterisation, ignore.case = TRUE)) {
                                      
                                    } else if (grepl("ord_HSROC|ordinal_HSROC|R&G|HSROC|Gatsonis", model_parameterisation, ignore.case = TRUE)) {
                                      
                                    } else { 
                                      stop("Invalid 'model_parameterisation' selected, must be one of: 'Jones', 'ord_biv', or 'ord_HSROC'")
                                    }
                                    
                              }
                              ##
                              ## ---- Store important parameters as class members:
                              ##
                              self$debugging <- debugging
                              ##
                              ## ---- Data:
                              ##
                              self$x <- x
                              ##
                              self$intercept_only <- intercept_only
                              self$cov_data <- cov_data ## covariates
                              ##
                              self$indicator_index_test_in_study <- indicator_index_test_in_study
                              ##
                              ## ---- make INTERNAL "basic_model_options" list:
                              ##
                              self$basic_model_options$network    <- network
                              self$basic_model_options$cts        <- cts
                              self$basic_model_options$prior_only <- prior_only
                              ##
                              ## ---- make INTERNAL "advanced_model_options" list:
                              ##
                              self$advanced_model_options$model_parameterisation        <- model_parameterisation
                              self$advanced_model_options$random_thresholds             <- random_thresholds
                              self$advanced_model_options$Dirichlet_random_effects_type <- Dirichlet_random_effects_type
                              self$advanced_model_options$box_cox  <- box_cox
                              self$advanced_model_options$softplus <- softplus
                              self$advanced_model_options$custom_file_name <- custom_file_name
                              ##
                              self$advanced_model_options$compound_symmetry <- compound_symmetry
                              self$advanced_model_options$hetero_sigma <- hetero_sigma
                              ##
                              ## ---- make INTERNAL "other_advanced_options" list:
                              ##
                              self$other_advanced_options$advanced_compile <- advanced_compile
                              ##
                              self$other_advanced_options$force_recompile <- force_recompile
                              self$other_advanced_options$quiet <- quiet
                              self$other_advanced_options$compile <- compile
                              ##
                              self$other_advanced_options$set_custom_CXX_CPP_flags <- set_custom_CXX_CPP_flags
                              self$other_advanced_options$CCACHE_PATH              <- CCACHE_PATH
                              self$other_advanced_options$custom_cpp_user_header_file_path <- custom_cpp_user_header_file_path
                              self$other_advanced_options$CXX_COMPILER_PATH  <- CXX_COMPILER_PATH
                              self$other_advanced_options$CXX_COMPILER_PATH  <- CXX_COMPILER_PATH
                              self$other_advanced_options$MATH_FLAGS   <- MATH_FLAGS
                              self$other_advanced_options$FMA_FLAGS    <- FMA_FLAGS
                              self$other_advanced_options$AVX_FLAGS    <- AVX_FLAGS
                              self$other_advanced_options$THREAD_FLAGS <- THREAD_FLAGS
                              ##
                              ## ---- priors:
                              ##
                              self$priors <- priors
                              ##
                              ## ---- initial values (per chain):
                              ##
                              self$init_lists_per_chain <- init_lists_per_chain
                              ##
                              ## ---- MCMC params (only n_chains can be inputted for initialisation):
                              ##
                              self$MCMC_params$n_chains <- n_chains
                              ##
                              # self$cmdstanr_args <- cmdstanr_args
                              ##
                              # message(("aaa_class_hello_2"))
                              ##
                              ## -----------  call initialising fn's: ------------------------------------------------------------------------------------------------------------
                              ##
                              local_prep_data_and_model_outs <-  MetaOrdDTA:::prep_data_and_model(  
                                                                                    debugging = self$debugging,
                                                                                    ##
                                                                                    x = self$x,
                                                                                    ##
                                                                                    intercept_only = self$intercept_only,
                                                                                    cov_data = self$cov_data, ## covariates
                                                                                    ##
                                                                                    indicator_index_test_in_study = self$indicator_index_test_in_study,
                                                                                    ##
                                                                                    internal_obj = self$internal_obj,
                                                                                    ##
                                                                                    basic_model_options    = self$basic_model_options,
                                                                                    advanced_model_options = self$advanced_model_options,
                                                                                    MCMC_params            = self$MCMC_params,
                                                                                    other_advanced_options = self$other_advanced_options,
                                                                                    ##
                                                                                    priors = self$priors,
                                                                                    init_lists_per_chain = self$init_lists_per_chain)
                              
                              ##
                              ## ---- Update R grouped lists:
                              ##
                              self$internal_obj           <- local_prep_data_and_model_outs$internal_obj
                              ##
                              self$basic_model_options    <- local_prep_data_and_model_outs$basic_model_options
                              self$advanced_model_options <- local_prep_data_and_model_outs$advanced_model_options
                              self$MCMC_params            <- local_prep_data_and_model_outs$MCMC_params
                              self$other_advanced_options <- local_prep_data_and_model_outs$other_advanced_options
                              ##
                              self$priors <- local_prep_data_and_model_outs$priors
                              ##
                              self$init_lists_per_chain <- local_prep_data_and_model_outs$init_lists_per_chain
                              
                              return(self)
                              
                              
                     },
                    ##  
                    ## --------  wrap the sample fn -------------------------------------------------------------------------------------------------------------------------------
                    ##
                    #'@description
                    #'Sample from the model
                    #'@param init_lists_per_chain List of initial values for each chain. See class documentation for details.
                    #'@param model_args_list List of model arguments. See class documentation for details.
                    #'@param stan_data_list List of Stan data (optional). See class documentation for details.
                    #'@param x The dataset. See class documentation for details.
                    #'@param n_chains Number of chains used for sampling (i.e., post-burnin). Note that this must match the length of "init_lists_per_chain".
                    #'@param n_superchains  Number of superchains for nested R-hat (nR-hat; see Margossian et al, 2023) computation. 
                    #'@param seed Random MCMC seed.
                    #'@param n_burnin The number of burnin iterations. Default is 500. 
                    #'@param n_iter The number of sampling (i.e., post-burnin) iterations. Default is 1000. 
                    #'@param adapt_delta The Metropolis-Hastings target acceptance rate. Default is 0.80. If there are divergences, sometimes increasing this can help, at the cost 
                    #'of efficiency (since this will decrease the step-size (epsilon) and increase the number of leapfrog steps (L) per iteration).
                    #'@param max_treedepth Proportional to the maximum number of leapfrog steps. The default is 10.
                    #'@param metric_shape The shape of the metric to use for the main parameters, which is adapted during the burnin period. Can be either \code{"diag"} or 
                    #'\code{"dense"}. The default is \code{"diag"} (i.e. a diagonal metric).
                    #'@return Returns self invisibly, allowing for method chaining of this classes (MetaOrd_model) methods. E.g.: model$sample(...)$summary(...). 
                    ##
                    sample = function(  debugging = self$debugging,
                                        ##
                                        x = self$x, ## data is allowed to be updated without re-initialising (as long as "network" and "cts" are the same)
                                        ##
                                        intercept_only = self$intercept_only, ## data is allowed to be updated without re-initialising
                                        cov_data = self$cov_data, ## data is allowed to be updated without re-initialising
                                        ##
                                        indicator_index_test_in_study = self$indicator_index_test_in_study,
                                        ##
                                        # cts = self$basic_model_options$cts,
                                        ##
                                        network = self$basic_model_options$network,
                                        compound_symmetry =  self$advanced_model_options$compound_symmetry,
                                        hetero_sigma =  self$advanced_model_options$hetero_sigma,
                                        ##
                                        prior_only = self$basic_model_options$prior_only,
                                        ##
                                        model_parameterisation        = self$advanced_model_options$model_parameterisation,
                                        random_thresholds             = self$advanced_model_options$random_thresholds,
                                        Dirichlet_random_effects_type = self$advanced_model_options$Dirichlet_random_effects_type,
                                        ##
                                        box_cox  = self$advanced_model_options$box_cox,
                                        softplus = self$advanced_model_options$softplus,
                                        ##
                                        custom_file_name = self$advanced_model_options$custom_file_name,
                                        ##
                                        priors = self$priors,
                                        ##
                                        init_lists_per_chain = self$init_lists_per_chain,
                                        ##
                                        ## MCMC / Stan params:
                                        ##
                                        algorithm = "Stan",
                                        ##
                                        seed          = NULL,
                                        n_superchains = NULL,
                                        n_chains = self$MCMC_params$n_chains,
                                        n_burnin      = 1000,
                                        n_iter        = 1000,
                                        adapt_delta   = 0.80,
                                        max_treedepth = 10,
                                        metric_shape  = "diag_e",
                                        ##
                                        cmdstanr_args = list()
                                        ) {
                                
                                if (model_parameterisation %in% c("Jones", "jones")) {
                                        
                                        cts <- TRUE
                                        random_thresholds <- FALSE
                                        Dirichlet_random_effects_type <- "none"
                                        ##
                                        if (is.null(box_cox)) { 
                                          warning("box_cox not specificed - setting to TRUE")
                                          box_cox <- TRUE
                                        }
                                  
                                } else  {
                                        
                                        cts <- FALSE
                                        box_cox <- FALSE
                                        ##
                                        if (is.null(random_thresholds)) { 
                                          warning("random_thresholds not specificed - setting to FALSE")
                                          random_thresholds <- FALSE
                                        }
                                        ##
                                        if (random_thresholds == TRUE) {
                                          if (is.null(Dirichlet_random_effects_type)) { 
                                            warning("Dirichlet_random_effects_type not specificed - setting to 'kappa' (recommended)")
                                            Dirichlet_random_effects_type <- "kappa"
                                          }
                                        }
                                        ##
                                        if (grepl("ord_biv|ordinal_biv|biv|Xu|O-", model_parameterisation, ignore.case = TRUE)) {
                                          
                                        } else if (grepl("ord_HSROC|ordinal_HSROC|R&G|HSROC|Gatsonis", model_parameterisation, ignore.case = TRUE)) {
                                          
                                        } else { 
                                          stop("Invalid 'model_parameterisation' selected, must be one of: 'Jones', 'ord_biv', or 'ord_HSROC'")
                                        }
                                  
                                }
                                ##
                                ## ---- Validate initialization:
                                ##
                                if (is.null(self$internal_obj)) {
                                  stop("Model was not properly initialize")
                                }
                                ##
                                # seed <- if_null_then_set_to(seed, 123) 
                                # n_superchains <- if_null_then_set_to(n_superchains, n_chains)
                                # ##
                                # ## ---- Extract params which can't be updated (i.e. params which require "$new()" again):
                                # ##
                                # network   <- self$network
                                # cts       <- self$cts
                                ##
                                ## ---- Update MCMC params:
                                ##
                                self$MCMC_params$n_superchains <- n_superchains
                                self$MCMC_params$n_chains <- n_chains
                                self$MCMC_params$seed <- seed
                                self$MCMC_params$n_burnin <- n_burnin
                                self$MCMC_params$n_iter <- n_iter
                                self$MCMC_params$adapt_delta <- adapt_delta
                                self$MCMC_params$max_treedepth <- max_treedepth
                                self$MCMC_params$metric_shape <- metric_shape
                                ##
                                # self$cmdstanr_args <- cmdstanr_args
                                ##
                                ## ---- Update params (if any changed):
                                ##
                                params_same <- 1
                                ##
                                {
                                    ##
                                    ## Update class members ** only if ** new values provided:
                                    ##
                                    # if (!identical(self$cmdstanr_args, cmdstanr_args))  { self$cmdstanr_args <- cmdstanr_args ; params_same <- 0 }
                                    # if (!identical(self$algorithm, algorithm))  { self$algorithm <- algorithm ; params_same <- 0 }
                                    ##
                                    if (!identical(self$x, x))  { self$x <- x ; params_same <- 0 }
                                    ##
                                    if (!identical(self$intercept_only, intercept_only))  { self$intercept_only <- intercept_only ; params_same <- 0 }
                                    if (!identical(self$cov_data, cov_data))  { self$cov_data <- cov_data ; params_same <- 0 }
                                    ##
                                    if (!identical(self$indicator_index_test_in_study, indicator_index_test_in_study))  { self$indicator_index_test_in_study <- indicator_index_test_in_study ; params_same <- 0 }
                                    ##
                                    if (!identical(self$basic_model_options$network, network)) { self$basic_model_options$network <- network ; params_same <- 0 }
                                    if (!identical(self$basic_model_options$cts, cts)) { self$basic_model_options$cts <- cts ; params_same <- 0 }
                                    if (!identical(self$basic_model_options$prior_only, prior_only)) { self$basic_model_options$prior_only <- prior_only ; params_same <- 0 }
                                    ##
                                    if (!identical(self$advanced_model_options$model_parameterisation, model_parameterisation))               { self$advanced_model_options$model_parameterisation        <- model_parameterisation        ; params_same <- 0 }
                                    if (!identical(self$advanced_model_options$random_thresholds, random_thresholds))                         { self$advanced_model_options$random_thresholds             <- random_thresholds             ; params_same <- 0 }
                                    if (!identical(self$advanced_model_options$Dirichlet_random_effects_type, Dirichlet_random_effects_type)) { self$advanced_model_options$Dirichlet_random_effects_type <- Dirichlet_random_effects_type ; params_same <- 0 }
                                    if (!identical(self$advanced_model_options$box_cox, box_cox))   { self$advanced_model_options$box_cox  <- box_cox  ; params_same <- 0 }
                                    if (!identical(self$advanced_model_options$softplus, softplus)) { self$advanced_model_options$softplus <- softplus ; params_same <- 0 }
                                    if (!identical(self$advanced_model_options$custom_file_name, custom_file_name)) { self$advanced_model_options$custom_file_name <- custom_file_name ; params_same <- 0 }
                                    if (!identical(self$advanced_model_options$compound_symmetry, compound_symmetry)) { self$advanced_model_options$compound_symmetry <- compound_symmetry ; params_same <- 0 }
                                    if (!identical(self$advanced_model_options$hetero_sigma, hetero_sigma)) { self$advanced_model_options$hetero_sigma <- hetero_sigma ; params_same <- 0 }
                                    ##
                                    if (!identical(self$priors, priors))         { self$priors     <- priors     ; params_same <- 0 }
                                    ##
                                    if (!identical(self$init_lists_per_chain, init_lists_per_chain))  { self$init_lists_per_chain <- init_lists_per_chain ; params_same <- 0 }
                                    ##
                                    if (!identical(self$MCMC_params$n_chains, n_chains)) { self$MCMC_params$n_chains <- n_chains ; params_same <- 0 }
                                    # ##
                                    # if (!identical(self$compute_sim_study_metrics, compute_sim_study_metrics)) { self$compute_sim_study_metrics <- compute_sim_study_metrics ; params_same <- 0 }
                                    # if (!identical(self$vec_index_inner_thr, vec_index_inner_thr)) { self$vec_index_inner_thr <- vec_index_inner_thr ; params_same <- 0 }
                                    
                                }
                            
                                # then update model if any of needed parameters changed
                                if (params_same == 0) {
                                    ##
                                    ## ---------- call "update_model" fn  -----------------------------------------------------------------------------------------
                                    ##
                                    #     self$model_prep_obj <-   MetaOrdDTA:::update_model(    test_type = test_type,
                                    #                                                       model_prep_obj = model_prep_obj,
                                    #                                                       y = self$y,
                                    #                                                       n_chains = self$n_chains
                                    #                                                       init_lists_per_chain = self$init_lists_per_chain,
                                    #                                                       model_args_list = self$model_args_list,
                                    #                                                       stan_data_list =  self$stan_data_list)
    
                                }
                                #
                                if (!is.null(adapt_delta) && (adapt_delta <= 0 || adapt_delta >= 1)) {
                                  stop("adapt_delta must be between 0 and 1")
                                }
                                
                                ##
                                ## -----------  call "initialise_update_run_model" fn ----------------------------------------------------------------------------
                                ##
                                local_model_samples_obj <-       MetaOrdDTA:::initialise_update_run_model(     
                                                                                               debugging = self$debugging,
                                                                                               ##
                                                                                               x = self$x,
                                                                                               ##
                                                                                               intercept_only = self$intercept_only,
                                                                                               cov_data = self$cov_data, ## covariates
                                                                                               ##
                                                                                               indicator_index_test_in_study = self$indicator_index_test_in_study,
                                                                                               ##
                                                                                               internal_obj = self$internal_obj,
                                                                                               ##
                                                                                               basic_model_options = self$basic_model_options,
                                                                                               advanced_model_options = self$advanced_model_options,
                                                                                               MCMC_params = self$MCMC_params,
                                                                                               other_advanced_options = self$other_advanced_options,
                                                                                               ##
                                                                                               priors = self$priors,
                                                                                               ##
                                                                                               init_lists_per_chain = self$init_lists_per_chain,
                                                                                               ##
                                                                                               cmdstanr_args = cmdstanr_args)

                                self$internal_obj <- local_model_samples_obj$internal_obj
                                ##
                                self$basic_model_options    <- local_model_samples_obj$basic_model_options
                                self$advanced_model_options <- local_model_samples_obj$advanced_model_options
                                self$MCMC_params            <- local_model_samples_obj$MCMC_params
                                ##
                                self$priors <- local_model_samples_obj$priors
                                ##
                                self$init_lists_per_chain <- local_model_samples_obj$init_lists_per_chain
                                                                                     
                                return(invisible(self))
                          
                                        
                      
                    },
                    ##        
                    ## --------  wrap the create_summary_and_traces fn + call the "BayesMVP::MVP_class_plot" R6 class  ----------------------------------------------------------
                    ##
                    #'@description
                    #'Create and compute summary statistics, traces and model diagnostics. 
                    #'@param compute_main_params Whether to compute the main parameter summaries. Default is TRUE.
                    #'@param compute_transformed_parameters Whether to compute transformed parameter summaries. Default is TRUE. For Stan models, this will be for 
                    #'all of the parameters defined in the "transformed parameters" block, EXCEPT for the (transformed) nuisance parameters and log_lik (see "save_log_lik_trace"
                    #'for more information on log_lik). 
                    #'@param compute_generated_quantities Whether to compute the summaries for generated quantities. Default is TRUE.
                    #'@param save_log_lik_trace Whether to save the log-likelihood (log_lik) trace. Default is FALSE. For Stan models, this will only work 
                    #'if there is a "log_lik" parameter defined in the "transformed parameters" model block.
                    #'@param compute_nested_rhat Whether to compute the nested rhat diagnostic (nR-hat) (Margossian et al, 2023). This is useful when
                    #'running many (usually short) chains. Also see "n_superchains" argument. 
                    #'@param n_superchains The number of superchains to use for the computation of nR-hat. 
                    #'Only relevant if \code{compute_nested_rhat = TRUE}.
                    #'@param save_trace_tibbles Whether to save the trace as tibble dataframes as well as 3D arrays. Default is FALSE. 
                    #'@return Returns a new MetaOrd_plot_and_diagnose object (from the "MetaOrd_plot_and_diagnose" R6 class) for creating MCMC diagnostics and plots.
                    ##
                    summary = function(       debugging = self$debugging,
                                              ##
                                              compute_main_params = TRUE,
                                              compute_transformed_parameters = TRUE,
                                              compute_generated_quantities = TRUE,
                                              save_log_lik_trace = TRUE,
                                              compute_nested_rhat = NULL,
                                              n_superchains = NULL,
                                              ##
                                              save_trace_tibbles = FALSE,
                                              ##
                                              use_BayesMVP_for_faster_summaries = NULL
                                              ##
                                              # other_args = list()
                                              ) {
                        
                            ## validate initialization:
                            if (is.null(self$internal_obj)) {
                              stop("Model was not properly initialize")
                            }
                            ##
                            if (debugging) {
                              message("Summarising fitted MetaOrdDTA model")
                              str(self$internal_obj, max.level = 1L)
                            }
                            ##
                            use_BayesMVP_for_faster_summaries <- if_null_then_set_to( use_BayesMVP_for_faster_summaries, 
                                                                                      MetaOrdDTA:::check_if_BayesMVP_R_pkg_installed(TRUE, FALSE))
                            ##
                            # initial_object <- self$initial_object
                            ##
                            # stan_model_sample_output <- self$model_samples_obj$full_outputs$stan_model_sample_output ## add "stan_model_sample_output" to "internal_obj" ???
                            # stan_param_names_list    <- self$internal_obj$stan_param_names_list
                            ##
                            ## create model fit object (includes model summary tables + traces + divergence info) by calling "BayesMVP::create_summary_and_traces" -----------------
                            ##
                            local_model_summary_and_trace_obj <- MetaOrdDTA:::create_summary_and_traces(  
                                                                                  #### debugging = self$debugging,
                                                                                  ##
                                                                                  package = "MetaOrdDTA",
                                                                                  use_bridgestan = FALSE,
                                                                                  use_BayesMVP_for_faster_summaries = use_BayesMVP_for_faster_summaries,
                                                                                  ##
                                                                                  # stan_model_sample_output = stan_model_sample_output,
                                                                                  # stan_param_names_list = stan_param_names_list,
                                                                                  internal_obj = self$internal_obj,
                                                                                  MCMC_params  = self$MCMC_params,
                                                                                  ##
                                                                                  model_results  = NULL,
                                                                                  init_object = NULL,
                                                                                  ##
                                                                                  n_nuisance = 0,
                                                                                  ##
                                                                                  compute_main_params = compute_main_params,
                                                                                  compute_transformed_parameters = compute_transformed_parameters,
                                                                                  compute_generated_quantities = compute_generated_quantities,
                                                                                  ##
                                                                                  save_log_lik_trace = save_log_lik_trace,
                                                                                  ##
                                                                                  compute_nested_rhat = compute_nested_rhat,
                                                                                  n_superchains = n_superchains,
                                                                                  save_trace_tibbles = save_trace_tibbles
                                                                                  ##
                                                                                  # other_args = other_args
                                                                                  )
                            
                            self$internal_obj$HMC_info <- local_model_summary_and_trace_obj$HMC_info
                            self$internal_obj$efficiency_info <- local_model_summary_and_trace_obj$efficiency_info
                            self$internal_obj$summaries <- local_model_summary_and_trace_obj$summaries
                            self$internal_obj$traces <- local_model_summary_and_trace_obj$traces
                            
                            return(self)
                      
                   },
                   ##
                   ## ------------- plots + MCMC & efficiency diagnostics:   ------------------------------------------------------------------------------------------------------
                   ##
                   ##
                   ## ---- Get divergence info:
                   ##
                   #'@return Returns self$internal_obj$summaries$divergences obtained from the main "MetaOrd_model" R6 class. 
                   #'self$internal_obj$summaries$divergences contains divergence information (e.g. % of transitions which diverged). 
                   get_divergences = function() {
                     return(self$internal_obj$HMC_info$HMC_diagnostic_info$divergences)
                   },
                   ##
                   ## ---- HMC diagnostic info:
                   ##
                   #'@return Returns self$internal_obj$summaries$divergences obtained from the main "MetaOrd_model" R6 class. 
                   #'self$internal_obj$summaries$divergences contains divergence information (e.g. % of transitions which diverged). 
                   get_HMC_diagnostic_info = function() { 
                     return(self$internal_obj$HMC_info$HMC_diagnostic_info) ## Extract the "HMC_info" list
                   },
                   ##
                   ## --------- Convenience method for HMC algorithm metrics:
                   ##
                   #' @return Returns a list containing HMC algorithm info (including timings and some efficiency metrics)
                   get_HMC_info = function() { 
                     return(self$internal_obj$HMC_info) ## Extract the "HMC_info" list
                   },
                   ##
                   ## ---- Summary tibble methods:
                   ##
                   #'@return Returns main parameter summaries as a tibble.
                   ##
                   get_summary_main = function() {  #### --------------------------------------------------- WORKING
                     return(self$internal_obj$summaries$summary_tibbles$summary_tibble_main_params)
                   },
                   #'@return Returns transformed parameter summaries as a tibble.
                   ##
                   get_summary_transformed = function() {
                     return(self$internal_obj$summaries$summary_tibbles$summary_tibble_transformed_parameters)
                   },
                   #'@return Returns generated quantities as a tibble.
                   ##
                   get_summary_generated_quantities = function() {  #### --------------------------------------------------- WORKING
                     return(self$internal_obj$summaries$summary_tibbles$summary_tibble_generated_quantities)
                   },
                   ##
                   ## ---- Posterior draws / 3D traces methods:
                   ##
                   #'Convenience method for ALL traces
                   #'@return Returns all traces including as tibbles (if enabled) and as 3D array's (iterations, chains, parameters). 
                   ##
                   get_all_traces_list = function() {  #### --------------------------------------------------- WORKING
                     return(self$internal_obj$traces)
                   },
                   #'Convenience method for log_lik trace
                   #'@return Returns  log_lik trace as 3D array (iterations, chains, parameters).
                   get_trace_main = function() { #### --------------------------------------------------- WORKING
                     return(self$internal_obj$traces$traces_as_arrays$trace_params_main)
                   },
                   #'Convenience method for log_lik trace
                   #'@return Returns  log_lik trace as 3D array (iterations, chains, parameters).
                   get_trace_transformed = function() {  #### --------------------------------------------------- WORKING
                     return(self$internal_obj$traces$traces_as_arrays$trace_tp_wo_log_lik)
                   },
                   #'Convenience method for log_lik trace
                   #'@return Returns  log_lik trace as 3D array (iterations, chains, parameters).
                   get_trace_generated_quantities = function() { #### --------------------------------------------------- WORKING
                     return(self$internal_obj$traces$traces_as_arrays$trace_gq)
                   },
                   #'Convenience method for log_lik trace
                   #'@return Returns  log_lik trace as 3D array (iterations, chains, parameters).
                   get_trace_log_lik = function() {  #### --------------------------------------------------- WORKING
                     return(self$internal_obj$traces$traces_as_arrays$trace_log_lik)
                   },
                   #'Convenience method for posterior draws (as tibbles)
                   #'@return Returns trace as tibbles. The trace gets returned as 3 seperate tibbles (one tibble for the 
                   #'main parameters, one for the transformed parameters, and one for generated quantities). 
                   ##
                   get_traces_as_tibbles = function() {  #### --------------------------------------------------- WORKING
                     list(
                       trace_as_tibble_main_params          =  self$internal_obj$traces$traces_as_tibbles$trace_params_main_tibble,
                       trace_as_tibble_transformed_params   =  self$internal_obj$traces$traces_as_tibbles$trace_transformed_params_tibble,
                       trace_as_tibble_generated_quantities =  self$internal_obj$traces$traces_as_tibbles$trace_generated_quantities_tibble
                     )
                   },
                   ##
                   ## --------- Convenience method for efficiency metrics:
                   ##
                   #' @return Returns a list containing efficiency metrics including:
                   #' \itemize{
                   #'   \item time_burnin: Time spent in burnin phase
                   #'   \item time_sampling: Time spent sampling
                   #'   \item time_total_MCMC: Total time spent doing MCMC/HMC sampling, excluding time to compute summaries and diagnostics.
                   #'   \item time_total_inc_summaries: Total time spent doing MCMC/HMC sampling, including time to compute summaries and diagnostics.
                   #'   \item Min_ESS_main: The minimum ESS of the main model parameters. 
                   #'   \item Min_ESS_per_sec_sampling: The minimum ESS per second for the sampling phase. 
                   #'   \item Min_ESS_per_sec_total: The minimum ESS per second for the total model run time, including any time spent computing
                   #'   summaries and diagnostics. 
                   #'   \item Min_ESS_per_grad_sampling:  The minimum ESS per gradient evaluation for the sampling phase. 
                   #'   \item grad_evals_per_sec: The number of gradient evaluations performed per second. 
                   #'   \item est_time_to_100_ESS_sampling: The estimated sampling time to reach a minimum ESS of 100.
                   #'   \item est_time_to_1000_ESS_sampling: The estimated sampling time to reach a minimum ESS of 1000.
                   #'   \item est_time_to_10000_ESS_sampling: The estimated sampling time to reach a minimum ESS of 10,000.
                   #'   \item est_time_to_100_ESS_wo_summaries: The estimated total time (expluding time to compute model summaries and diagnostics)
                   #'   to reach a minimum ESS of 100.
                   #'   \item est_time_to_1000_ESS_wo_summaries: The estimated total time (expluding time to compute model summaries and diagnostics)
                   #'   to reach a minimum ESS of 1000.
                   #'   \item est_time_to_10000_ESS_wo_summaries: The estimated total time (expluding time to compute model summaries and diagnostics)
                   #'   to reach a minimum ESS of 10,000.
                   #'   \item est_time_to_100_ESS_inc_summaries: The estimated total time (including time spent computing model summaries and 
                   #'   diagnostics) to reach a minimum ESS of 100.
                   #'   \item est_time_to_1000_ESS_inc_summaries: The estimated total time (including time spent computing model summaries and 
                   #'   diagnostics) to reach a minimum ESS of 1000.
                   #'   \item est_time_to_10000_ESS_inc_summaries: The estimated total time (including time spent computing model summaries and 
                   #'   diagnostics) to reach a minimum ESS of 10,000.
                   #' }
                   ##
                   get_efficiency_metrics = function() {
                     return(self$internal_obj$efficiency_info) ## Extract the "efficiency_info" list
                   },
                   ##
                   ## --------- Convenience fn to compute "time to X ESS" using the "$time_to_target_ESS()" method
                   ##
                   #'@param target_ESS The target ESS. 
                   #'@return Returns a list called "target_ESS_times" which contains the estimated sampling time to reach the target ESS 
                   #'("sampling_time_to_target_ESS"), the estimated total time excluding the time it takes to compute model summaries 
                   #'("total_time_to_target_ESS_wo_summaries"), and the total estimated time to reach the 
                   #'target ESS ("total_time_to_target_ESS_with_summaries").
                   time_to_target_ESS = function(target_ESS) {
                     
                     if (is.null(self$internal_obj)) {
                       stop("No summary object available")   
                     }
                     
                     ## get required efficiency info first
                     Min_ESS_main   <- self$internal_obj$efficiency_info$Min_ESS_main
                     time_burnin    <- self$internal_obj$efficiency_info$time_burnin
                     time_sampling  <- self$internal_obj$efficiency_info$time_sampling
                     time_summaries <- self$internal_obj$efficiency_info$time_summaries
                     
                     sampling_time_to_Min_ESS <- time_sampling
                     
                     ## est. sampling time to target_ESS
                     sampling_time_to_target_ESS <- (target_ESS / Min_ESS_main) * sampling_time_to_Min_ESS
                     
                     ## total w/o summary time
                     total_time_to_target_ESS_wo_summaries <-   time_burnin + sampling_time_to_target_ESS
                     
                     ## total w/ summary time
                     summary_time_to_Min_ESS <- time_summaries
                     summary_time_to_target_ESS <- (target_ESS / Min_ESS_main) * summary_time_to_Min_ESS
                     
                     total_time_to_target_ESS_with_summaries <-   time_burnin + sampling_time_to_target_ESS   + summary_time_to_target_ESS
                     
                     target_ESS_times <- list(sampling_time_to_target_ESS = sampling_time_to_target_ESS,
                                              total_time_to_target_ESS_wo_summaries = total_time_to_target_ESS_wo_summaries,
                                              total_time_to_target_ESS_with_summaries = total_time_to_target_ESS_with_summaries)
                     
                     return(target_ESS_times)
                     
                     
                   },
                   ##
                   ## --------- Convenience fn to compute "n_iter to X ESS" using the "$iter_to_target_ESS()" method
                   ##
                   #'@param target_ESS The target ESS. 
                   #'@return Returns a list called "target_ESS_iter" which contains the estimated number of sampling iterations (n_iter) to 
                   #'reach the target ESS ("sampling_iter_to_target_ESS").
                   iter_to_target_ESS = function(target_ESS) {
                     
                               if (is.null(self$internal_obj)) {
                                 stop("No internal_obj object available")   
                               }
                               
                               ## get required efficiency info first:
                               Min_ESS_main   <- self$internal_obj$efficiency_info$Min_ESS_main
                               n_iter         <- self$internal_obj$efficiency_info$n_iter
                               n_iter_to_min_ESS <- n_iter
                               
                               ## est. sampling time to target_ESS
                               sampling_iter_to_target_ESS <- (target_ESS / Min_ESS_main) * n_iter_to_min_ESS
                               
                               
                               target_ESS_iter <- list(sampling_iter_to_target_ESS = sampling_iter_to_target_ESS)
                               
                               return(target_ESS_iter)
                     
                   },
                   ##
                   ## --------- MCMC trace plots method
                   ##
                   #'@param debugging Whether debugging mode is turned on (\code{TRUE}) or off (\code{FALSE}). The default is FALSE.
                   #'@param param_string The parameters to generate trace plots for. This is a character vector - e.g. to plot trace plots for 
                   #'beta: param_string = c("beta"). 
                   #'The default is NULL and the trace plots for all model parameters (which have a trace array) will be plotted. 
                   #'@param batch_size The number of trace plots to display per panel. Default is 9. 
                   #'@return If no parameters specified (i.e. param_string = NULL), then this will return an object containing the trace plots for all
                   #'model parameters  which have a trace array.
                   plot_traces = function( debugging = FALSE,
                                           param_string = NULL, 
                                           condition = "exact_match",
                                           batch_size = 9
                   ) {
                     
                               if (is.null(self$internal_obj)) {
                                 stop("No internal_obj object available")  
                               }
                               if (!is.null(param_string) && !is.character(param_string)) {
                                 stop("param_string must be NULL or a character vector")
                               }
                               if (!is.numeric(batch_size) || batch_size < 1) {
                                 stop("batch_size must be a positive integer")
                               }
                               if (debugging) {
                                 cat("Summary object structure:\n")
                                 str(self$internal_obj)
                               }
                               # ##
                               # ## ---- Format the array for Bayesplot R package:
                               # ##
                               # draws_array <- MetaOrdDTA:::format_named_array_for_bayesplot(self$internal_obj$traces$traces_as_arrays$draws_array)
                               ##
                               if (is.null(param_string)) { # plot all param_string
                                 
                                     bayesplot::mcmc_trace(draws_array)
                                 
                               } else {   # plot specific param_string using custom "plot_multiple_params_batched" fn - mimics Stan's method but uses bayesplot 
                                 
                                     MetaOrdDTA:::plot_param_group_batched(       debugging = debugging,
                                                                                  draws_array = self$internal_obj$traces$traces_as_arrays$draws_array,
                                                                                  plot_type = "trace",
                                                                                  param_string = param_string,
                                                                                  condition = condition,
                                                                                  batch_size = batch_size)
                                                                        
                               }
                               
                   },
                   ##
                   ## --------- MCMC density plots method
                   ##
                   #'@param debugging Whether debugging mode is turned on (\code{TRUE}) or off (\code{FALSE}). The default is FALSE.
                   #'@param param_string The parameters to generate posterior density plots for. This is a character vector - e.g. to plot the 
                   #'posterior density plots for beta: param_string = c("beta"). 
                   #'The default is NULL and the posterior density plots for all model parameters (which have a trace array) will be plotted. 
                   #'@param batch_size The number of posterior density plots to display per panel. Default is 9. 
                   #'@return If no parameters specified (i.e. param_string = NULL), then this will return an object containing the posterior density
                   #'plots for all model parameters which have a trace array.
                   plot_densities = function( debugging = FALSE,
                                              param_string = NULL, 
                                              condition = "exact_match",
                                              batch_size =  9
                   ) {
                     
                               if (is.null(self$internal_obj)) {
                                 stop("No internal_obj object available")  
                               }
                               if (!is.null(param_string) && !is.character(param_string)) {
                                 stop("param_string must be NULL or a character vector")
                               }
                               if (!is.numeric(batch_size) || batch_size < 1) {
                                 stop("batch_size must be a positive integer")
                               }
                               if (debugging) {
                                 cat("Summary object structure:\n")
                                 str(self$internal_obj)
                               }
                               # ##
                               # ## ---- Format the array for Bayesplot R package:
                               # ##
                               # draws_array <- MetaOrdDTA:::format_named_array_for_bayesplot(self$internal_obj$traces$traces_as_arrays$draws_array)
                               ##
                               if (is.null(param_string)) { # plot all param_string
                                 
                                     bayesplot::mcmc_dens(draws_array)
                                 
                               } else {   # plot specific param_string using custom "plot_multiple_params_batched" fn - mimics Stan's method but uses bayesplot
                                 
                                     MetaOrdDTA:::plot_param_group_batched(    plot_type = "density",
                                                                               debugging = debugging,
                                                                               draws_array = self$internal_obj$traces$traces_as_arrays$draws_array,
                                                                               param_string = param_string,
                                                                               condition = condition,
                                                                               batch_size = batch_size)
                               }
                     
                   },
                   ##
                   ## --------- Extract Params from summary tibble's
                   ##
                   #'@param params ...
                   #'@param debugging ...
                   #'@param tibble ...
                   #'@param condition ...
                   #'@param combine ...
                   #'@return Extract params
                   extract_params = function( params,
                                              debugging = FALSE,
                                              tibble = NULL,
                                              condition = "exact_match",
                                              combine = TRUE
                   ) {
                     
                               if (is.null(self$internal_obj)) {
                                 stop("No model object ('internal_obj') available")
                               }
                               ##
                               if (debugging) {
                                 cat("Summary object structure:\n")
                                 str(self$internal_obj)
                               }
                               ##
                               if (is.null(tibble)) { 
                                   tibble_main <- self$internal_obj$summaries$summary_tibbles$summary_tibble_main_params
                                   tibble_tp   <- self$internal_obj$summaries$summary_tibbles$summary_tibble_transformed_parameters
                                   tibble_gq   <- self$internal_obj$summaries$summary_tibbles$summary_tibble_generated_quantities
                                   ##
                                   tibble <- rbind(tibble_main, tibble_tp, tibble_gq)
                               }
                               ##
                               outs <- MetaOrdDTA:::extract_params_from_tibble_batch(   debugging = debugging,
                                                                                        tibble = tibble,
                                                                                        param_strings_vec = params,
                                                                                        condition = condition,
                                                                                        combine = combine)
                               ##
                               return(outs)
                     
                   },
                   ##
                   ## --------- Run full K-fold cross-validation
                   ##
                   #'@param debugging ...
                   #'@param seed ...
                   #'@param K ...
                   #'
                   make_folds_k_fold_CV = function(  debugging = FALSE,
                                                     seed,
                                                     K
                                                  
                   ) {
                     
                               set.seed(seed, kind = "Mersenne-Twister")
                               ##
                               network <- self$basic_model_options$network
                               n_studies <- self$internal_obj$outs_data$stan_data_list$n_studies
                               indicator_index_test_in_study <- self$indicator_index_test_in_study
                               ##
                               fold_assignments <- create_folds(
                                 indicator_index_test_in_study = indicator_index_test_in_study,
                                 K = K,
                                 seed = seed,
                                 network = network,
                                 n_studies = n_studies
                               )
                               ## Inspect folds:
                               print(table(fold_assignments))
                               for (k in 1:K) {
                                 cat("Fold", k, ":", which(fold_assignments == k), "\n")
                               }
                               ##
                               return(fold_assignments)
                           
                   },
                   ##
                   ## --------- Run full K-fold cross-validation
                   ##
                   #'@param debugging ...
                   #'@param K ...
                   #'@param n_burnin ...
                   #'@param n_iter ...
                   #'@param n_chains ...
                   #'@param adapt_delta ...
                   #'@param max_treedepth ...
                   #'@param parallel ...
                   #'@param n_workers ...
                   #'@param output_dir ...
                   #'@param use_BayesMVP_for_faster_summaries ...
                   #'@return Extract params
                   #'
                   run_k_fold_CV = function(  debugging = FALSE,
                                              ##
                                              seed,
                                              fold_assignments = NULL,
                                              ##
                                              cmdstanr_args = NULL,
                                              ##
                                              init_lists_per_chain = NULL,
                                              priors = self$priors,
                                              ##
                                              n_burnin = 1000,
                                              n_iter = 1000,
                                              n_chains = NULL,
                                              adapt_delta = 0.80,
                                              max_treedepth = 10,
                                              ##
                                              parallel = NULL,
                                              n_workers = NULL,
                                              output_dir = metaorddta_output_path("cv_results"),
                                              ##
                                              use_BayesMVP_for_faster_summaries = NULL
                   ) {
                     
                               if (is.null(self$internal_obj)) {
                                 stop("No model samples object ('internal_obj') available")
                               }
                               use_BayesMVP_for_faster_summaries <- if_null_then_set_to( 
                                                                    use_BayesMVP_for_faster_summaries, 
                                                                    MetaOrdDTA:::check_if_BayesMVP_R_pkg_installed(TRUE, FALSE))
                               ##
                               model_prep_obj <- self
                               stan_data_list <- self$internal_obj$outs_data$stan_data_list
                               ##
                               if (is.null(fold_assignments)) { 
                                 stop("Please run the '$make_folds_k_fold_CV()' R6 class method to make folds before running K-fold CV")
                               }
                               ##
                               K <- length(unique(fold_assignments))
                               ##
                               if (is.null(n_chains)) { 
                                 n_chains <- 2
                               }
                               if (is.null(init_lists_per_chain)) {
                                 init_lists_per_chain <- self$init_lists_per_chain
                               }
                               init_lists_per_chain <- resize_init_list(
                                 init_lists_per_chain = init_lists_per_chain,
                                 n_chains_new = n_chains
                               )
                               ##
                               ## ---- Decide whether to use parallel version:
                               ##
                               if (is.null(parallel)) {
                                 # Use parallel if we have enough cores for at least 2 workers
                                 available_cores <- parallel::detectCores()
                                 potential_workers <- floor(available_cores / n_chains)
                                 use_parallel <- potential_workers >= 2 && K >= 4
                               } else {
                                 use_parallel <- parallel
                               }
                               ##
                               if (use_parallel) {
                                 cat("Using PARALLEL k-fold CV\n")
                                 outs <- R_fn_run_kfold_cv_parallel(   debugging = debugging,
                                                                       ##
                                                                       seed = seed,
                                                                       K = K,
                                                                       fold_assignments = fold_assignments,
                                                                       ##
                                                                       model_prep_obj = model_prep_obj,
                                                                       stan_data_list = stan_data_list,
                                                                       ##
                                                                       cmdstanr_args = cmdstanr_args,
                                                                       ##
                                                                       priors = priors,
                                                                       init_lists_per_chain = init_lists_per_chain,
                                                                       ##
                                                                       n_burnin = n_burnin,
                                                                       n_iter = n_iter,
                                                                       n_chains = n_chains,
                                                                       adapt_delta = adapt_delta,
                                                                       max_treedepth = max_treedepth,
                                                                       ##
                                                                       n_workers = n_workers,
                                                                       output_dir = output_dir,
                                                                       ##
                                                                       use_BayesMVP_for_faster_summaries = use_BayesMVP_for_faster_summaries)
                               } else { 
                                 cat("Using SEQUENTIAL k-fold CV\n")
                                 outs <- R_fn_run_kfold_cv(    debugging = debugging,
                                                               ##
                                                               seed = seed,
                                                               K = K,
                                                               fold_assignments = fold_assignments,
                                                               ##
                                                               model_prep_obj = model_prep_obj,
                                                               stan_data_list = stan_data_list,
                                                               ##
                                                               cmdstanr_args = cmdstanr_args,
                                                               ##
                                                               priors = priors,
                                                               init_lists_per_chain = init_lists_per_chain,
                                                               ##
                                                               n_burnin = n_burnin,
                                                               n_iter = n_iter,
                                                               n_chains = n_chains,
                                                               adapt_delta = adapt_delta,
                                                               max_treedepth = max_treedepth,
                                                               ##
                                                               output_dir = output_dir,
                                                               ##
                                                               use_BayesMVP_for_faster_summaries = use_BayesMVP_for_faster_summaries)
                               }
                               
                               ##
                               return(outs)
                     
                   },
                   ##
                   ## --------- Test accuracy-specific plots method: sROC plot (w/ confidence + prediction regions/curves):
                   ##
                   #'@param debugging ...
                   #'@param test_names ...
                   #'@param relevant_thresholds ...
                   #'@param df_true ... 
                   #'@param layout ...
                   #'@param separate_sROC_settings ...
                   #'@param by_test_sROC_settings ...
                   #'@param by_scenario_sROC_settings ...
                   #'@param by_test_grid_sROC_settings ...
                   #'@param by_scenario_grid_sROC_settings ...
                   #'@param by_test_grid_separate_sROC_settings ...
                   #'@param by_scenario_grid_separate_sROC_settings ...
                   #'@param aesthetic_mapping ...
                   #'@return sROC plot(s).
                   #'
                   plot_sROC = function( debugging = FALSE,
                                         ##
                                         test_names = NULL,
                                         relevant_thresholds = NULL,
                                         ##
                                         df_true = NULL,
                                         ##
                                         ## ---- SOME of these args are for MR-only:
                                         ##
                                         layout = c("separate", "by_test", "by_scenario", "grid", "grid_separate"),
                                         ##
                                         separate_sROC_settings = NULL,
                                         by_test_sROC_settings = NULL,
                                         by_scenario_sROC_settings = NULL, ## this is also used for non-MR
                                         by_test_grid_sROC_settings = NULL,
                                         by_scenario_grid_sROC_settings = NULL,
                                         by_test_grid_separate_sROC_settings = NULL,
                                         by_scenario_grid_separate_sROC_settings = NULL, ## this is also used for non-MR
                                         ##
                                         aesthetic_mapping = NULL
                   ) {
                     
                               if (is.null(self$internal_obj)) {
                                 stop("No model samples object ('internal_obj') available")
                               }
                               # ##
                               if (debugging) {
                                 cat("Summary object structure:\n")
                                 str(self$internal_obj)
                               }
                               ##
                               network <- self$basic_model_options$network
                               ##
                               if (network) { 
                                 n_index_tests <- self$internal_obj$outs_data$stan_data_list$n_index_tests
                                 print(paste("n_index_tests = ", n_index_tests))
                                 ##
                                 n_thr         <- self$internal_obj$outs_data$stan_data_list$n_thr
                                 print(paste("n_thr = ", n_thr))
                                 ##
                               } else { 
                                 n_index_tests <- 1     ## dummy arg - needed for compatability / R6 class
                                 n_thr         <- NULL  ## dummy arg - needed for compatability / R6 class
                               }
                               ##
                               n_covariates <- self$internal_obj$outs_data$stan_data_list$n_covariates_max
                               print(paste("n_covariates = ", n_covariates))
                               ##
                               if (n_covariates == 1) {  # No covariates - create simple baseline
                                       
                                       covariates <- FALSE
                                       use_MR_scenarios <- FALSE
                                       ##
                                       baseline_case_nd <- baseline_case_d <- list()
                                       for (t in 1:n_index_tests) {
                                         tn <- ifelse(!is.null(test_names), test_names[t], paste0("Test", t))
                                         baseline_case_nd[[tn]] <- 1
                                         baseline_case_d[[tn]] <- 1
                                       }
                                       ##
                                       scenario_results <- R_fn_compute_MR_multiple_scenarios(
                                                           model_prep_obj = self,
                                                           model_summary_and_trace_obj = self,
                                                           baseline_case_nd_list = list("Baseline (no covariates)" = baseline_case_nd),
                                                           baseline_case_d_list = list("Baseline (no covariates)" = baseline_case_d),
                                                           scenario_names = "Baseline (no covariates)",
                                                           test_names = test_names,
                                                           compute_comparisons = FALSE,
                                                           save_results = FALSE,
                                                           parallel = FALSE)
                                     
                                 
                               } else {   ## Meta-regression - MUST use stored scenarios
                                     
                                       covariates <- TRUE
                                       use_MR_scenarios <- TRUE
                                       ##
                                       if (is.null(self$baseline_outs_MR)) {
                                         stop("No MR scenarios computed. Run '$recompute_baseline()' first, 
                                              then try running '$extract_AUC()' again.")
                                       }
                                       ##
                                       scenario_results <- self$baseline_outs_MR
                                   
                               }
                               ##
                               combined_data <- combine_sroc_data_scenarios(scenario_results)
                               ##
                               plots <- create_sROC_plots_flexible(
                                 combined_data = combined_data,
                                 ##
                                 test_names = test_names,
                                 relevant_thresholds = relevant_thresholds,
                                 ##
                                 layout = layout,
                                 ##
                                 by_test_sROC_settings = by_test_sROC_settings,
                                 by_scenario_sROC_settings = by_scenario_sROC_settings,
                                 separate_sROC_settings = separate_sROC_settings,
                                 by_test_grid_sROC_settings = by_test_grid_sROC_settings,
                                 by_scenario_grid_sROC_settings = by_scenario_grid_sROC_settings,
                                 by_test_grid_separate_sROC_settings = by_test_grid_separate_sROC_settings,
                                 by_scenario_grid_separate_sROC_settings = by_scenario_grid_separate_sROC_settings)
                               ##
                               return(plots)
                     
                   },
                   ##
                   ## --------- Test accuracy-specific plots method: Accuracy (Se/Sp) vs. threshold plots
                   ##
                   #'@param debugging ...
                   #'@param test_names ...
                   #'@param relevant_thresholds ...
                   #'@param plot_type ...
                   #'@param df_true ... 
                   #'@param new_cov_data ...
                   #'@param layout ...
                   #'@param separate_threshold_settings ...
                   #'@param by_test_threshold_settings ...
                   #'@param by_scenario_threshold_settings ...
                   #'@param by_test_grid_threshold_settings ...
                   #'@param by_scenario_grid_threshold_settings ...
                   #'@param by_test_grid_separate_threshold_settings ...
                   #'@param by_scenario_grid_separate_threshold_settings ...
                   #'@param aesthetic_mapping ...
                   #'@param show_actual_thresholds ...
                   #'@return Accuracy vs. threshold plot(s).
                   plot_threshold = function( debugging = FALSE,
                                         ##
                                         test_names = NULL,
                                         relevant_thresholds = NULL,
                                         ##
                                         plot_type = c("Se", "Sp", "both"),
                                         ##
                                         df_true = NULL,
                                         ##
                                         ## ---- SOME of these args are for MR-only:
                                         ##
                                         new_cov_data = NULL,
                                         layout = c("separate", "by_test", "by_scenario", "grid", "grid_separate"),
                                         ##
                                         separate_threshold_settings = NULL,
                                         by_test_threshold_settings = NULL,
                                         by_scenario_threshold_settings = NULL, ## this is also used for non-MR
                                         by_test_grid_threshold_settings = NULL,
                                         by_scenario_grid_threshold_settings = NULL,
                                         by_test_grid_separate_threshold_settings = NULL,
                                         by_scenario_grid_separate_threshold_settings = NULL, ## this is also used for non-MR
                                         ##
                                         aesthetic_mapping = NULL,
                                         show_actual_thresholds = FALSE
                   ) { 
                     
                             if (is.null(self$internal_obj)) {
                               stop("No model samples object ('internal_obj') available")
                             }
                             # ##
                             if (debugging) {
                               cat("Summary object structure:\n")
                               str(self$internal_obj)
                             }
                             ##
                             network <- self$basic_model_options$network
                             ##
                             if (network) { 
                               n_index_tests <- self$internal_obj$outs_data$stan_data_list$n_index_tests
                               print(paste("n_index_tests = ", n_index_tests))
                               ##
                               n_thr         <- self$internal_obj$outs_data$stan_data_list$n_thr
                               print(paste("n_thr = ", n_thr))
                               ##
                             } else { 
                               n_index_tests <- 1     ## dummy arg - needed for compatability / R6 class
                               n_thr         <- NULL  ## dummy arg - needed for compatability / R6 class
                             }
                             ##
                             n_covariates <- self$internal_obj$outs_data$stan_data_list$n_covariates_max
                             print(paste("n_covariates = ", n_covariates))
                             ##
                             if (n_covariates == 1) {  # No covariates - create simple baseline
                                     
                                     covariates <- FALSE
                                     use_MR_scenarios <- FALSE
                                     ##
                                     baseline_case_nd <- baseline_case_d <- list()
                                     for (t in 1:n_index_tests) {
                                       tn <- ifelse(!is.null(test_names), test_names[t], paste0("Test", t))
                                       baseline_case_nd[[tn]] <- 1
                                       baseline_case_d[[tn]] <- 1
                                     }
                                     ##
                                     scenario_results <- R_fn_compute_MR_multiple_scenarios(
                                       model_prep_obj = self,
                                       model_summary_and_trace_obj = self,
                                       baseline_case_nd_list = list("Baseline (no covariates)" = baseline_case_nd),
                                       baseline_case_d_list = list("Baseline (no covariates)" = baseline_case_d),
                                       scenario_names = "Baseline (no covariates)",
                                       test_names = test_names,
                                       compute_comparisons = FALSE,
                                       save_results = FALSE,
                                       parallel = FALSE)
                               
                               
                             } else {   ## Meta-regression - MUST use stored scenarios
                               
                                     covariates <- TRUE
                                     use_MR_scenarios <- TRUE
                                     ##
                                     if (is.null(self$baseline_outs_MR)) {
                                       stop("No MR scenarios computed. Run '$recompute_baseline()' first, 
                                                    then try running '$extract_AUC()' again.")
                                     }
                                     ##
                                     scenario_results <- self$baseline_outs_MR
                               
                             }
                             ##
                             # Convert to combined data format
                             combined_data <- extract_threshold_data_scenarios(scenario_results)
                             ##
                             plots <- create_threshold_plots_flexible(
                               combined_data = combined_data,
                               ##
                               test_names = test_names,
                               relevant_thresholds = relevant_thresholds,
                               ##
                               plot_type = plot_type,
                               layout = layout,
                               ##
                               separate_threshold_settings = separate_threshold_settings,
                               by_test_threshold_settings = by_test_threshold_settings,
                               by_scenario_threshold_settings = by_scenario_threshold_settings,
                               by_test_grid_threshold_settings = by_test_grid_threshold_settings,
                               by_scenario_grid_threshold_settings = by_scenario_grid_threshold_settings,
                               by_test_grid_separate_threshold_settings = by_test_grid_separate_threshold_settings,
                               by_scenario_grid_separate_threshold_settings = by_scenario_grid_separate_threshold_settings,
                               ##
                               aesthetic_mapping = aesthetic_mapping,
                               ##
                               show_actual_thresholds = show_actual_thresholds)
                             ##
                             return(plots)
                             
                   },
                   ##
                   ## --------- Test accuracy-specific plots method: Accuracy (Se/Sp) vs. threshold plots
                   ##
                   #'@param debugging ...
                   #'@param comparison_type ...
                   #'@param test_names ...
                   #'@param relevant_thresholds ...
                   #'@param baseline_plot_settings ...
                   #'@param metareg_plot_settings ...
                   #'@param aesthetic_mapping ...
                   #'@return Accuracy vs. threshold plot(s).
                   plot_pairwise = function(  debugging = FALSE,
                                              ##
                                              comparison_type = "both",
                                              ##
                                              test_names = NULL,
                                              relevant_thresholds = NULL,
                                              ##
                                              baseline_plot_settings = NULL, ## NON-MR-only
                                              metareg_plot_settings = NULL, ## MR-only
                                              ##
                                              aesthetic_mapping = NULL
                   ) { 
                             
                             if (is.null(self$internal_obj)) {
                               stop("No model samples object ('internal_obj') available")
                             }
                             ##
                             if (debugging) {
                               cat("Summary object structure:\n")
                               str(self$internal_obj)
                             }
                             ##
                             network <- self$basic_model_options$network
                             ##
                             if (network) { 
                               n_index_tests <- self$internal_obj$outs_data$stan_data_list$n_index_tests
                               print(paste("n_index_tests = ", n_index_tests))
                               n_thr <- self$internal_obj$outs_data$stan_data_list$n_thr
                               print(paste("n_thr = ", n_thr))
                             } else { 
                               n_index_tests <- 1
                               n_thr <- NULL
                             }
                             ##
                             n_covariates <- self$internal_obj$outs_data$stan_data_list$n_covariates_max
                             print(paste("n_covariates = ", n_covariates))
                             ##
                             ## ---- Extract comparisons
                             ##
                             if (n_covariates == 1) {   # No covariates - extract directly
                               
                                     covariates <- FALSE
                                     tibble_comp <- self$extract_NMA_comparisons(test_names = test_names)
                               
                             } else { # Meta-regression - use STORED comparisons
                               
                                     covariates <- TRUE
                                     
                                     if (is.null(self$pairwise_comparisons_MR)) {
                                       stop("No pairwise comparisons computed. Run '$compute_pairwise_comparisons()' first.")
                                     }
                                     
                                     # Just retrieve the stored comparisons - don't recompute!
                                     all_comparisons <- list()
                                     for (scenario_name in names(self$pairwise_comparisons_MR)) {
                                       comps <- self$pairwise_comparisons_MR[[scenario_name]]
                                       comps$scenario <- scenario_name
                                       all_comparisons[[scenario_name]] <- comps
                                     }
                                     tibble_comp <- bind_rows(all_comparisons)
                                     
                             }
                             ##
                             ## ---- Call the unified plotting function:
                             ##
                             plots <-  plot_pairwise_differences(
                                       tibble_comp = tibble_comp,
                                       comparison_type = comparison_type,
                                       test_names = test_names,
                                       relevant_thresholds = relevant_thresholds,
                                       baseline_plot_settings = baseline_plot_settings,
                                       metareg_plot_settings = metareg_plot_settings)
                             ##
                             return(plots)
                             
                   },
                   ##
                   ## --------- Extract NMA-specific COMPARISONS summaries: Se/Sp diffs + ratio's
                   ##
                   #'@param debugging ...
                   #'@param test_names ...
                   #'@return NMA comparisons summary
                   #'
                   extract_NMA_comparisons = function( debugging = FALSE,
                                                       test_names = NULL,
                                                       MR_scenario = NULL
                   ) {
                     
                             if (is.null(self$internal_obj)) {
                               stop("No model samples object ('internal_obj') available")
                             }
                             ##
                             if (debugging) {
                               cat("Summary object structure:\n")
                               str(self$internal_obj)
                             }
                             ##
                             network <- self$basic_model_options$network
                             ##
                             if (network) { 
                               n_index_tests <- self$internal_obj$outs_data$stan_data_list$n_index_tests
                               print(paste("n_index_tests = ", n_index_tests))
                               n_thr <- self$internal_obj$outs_data$stan_data_list$n_thr
                               print(paste("n_thr = ", n_thr))
                             } else { 
                               n_index_tests <- 1
                               n_thr <- NULL
                             }
                             ##
                             n_covariates <- self$internal_obj$outs_data$stan_data_list$n_covariates_max
                             print(paste("n_covariates = ", n_covariates))
                             ##
                             if (n_covariates == 1) {
                               covariates <- FALSE
                             } else {
                               covariates <- TRUE
                             }
                             ##
                             if (covariates) { 
                               
                                     if (is.null(MR_scenario)) { 
                                       stop("You must supply 'MR_scenario' for meta-regression models")
                                     }
                                     
                                     # Convert scenario to baseline_case format (same as recompute_baseline does)
                                     baseline_vec <- c(1, unlist(MR_scenario))
                                     baseline_case_nd <- rep(list(baseline_vec), length(test_names))
                                     baseline_case_d <- rep(list(baseline_vec), length(test_names))
                                     
                                     stan_mod_samples_to_use <- R_fn_GQ_covariates_NMA(
                                       baseline_case_nd = baseline_case_nd,
                                       baseline_case_d = baseline_case_d,
                                       model_summary_and_trace_obj = self,
                                       model_prep_obj = self
                                     )
                                     
                                     trace_gq <- stan_mod_samples_to_use$draws()
                               
                             } else { 
                               
                                     stan_mod_samples_to_use <- self$internal_obj$outs_stan_sampling$stan_mod_samples
                                     trace_gq <- self$internal_obj$traces$traces_as_arrays$trace_gq
                                     
                             }
                             ##
                             ## ---- Use fast Rcpp fn:
                             ##
                             outs <- MetaOrdDTA:::R_fn_using_Rcpp_compute_NMA_comparisons_posthoc(
                               trace_gq = trace_gq,
                               n_thr = n_thr,
                               test_names = test_names
                             )
                             ##
                             return(outs)
                     
                   },
                   ##
                   ## --------- Extract NMA-specific PERFORMANCE summaries (e.g., DOR)
                   ##
                   #'@param debugging ...
                   #'@param test_names ...
                   #'@param new_cov_data ...
                   #'@return NMA performance metrics (e.g., DOR) summary
                   #'
                   extract_NMA_performance = function( debugging = FALSE,
                                                       test_names = NULL,
                                                       new_cov_data = NULL
                   ) {
                     
                               if (is.null(self$internal_obj)) {
                                 stop("No model samples object ('internal_obj') available")
                               }
                               if (is.null(self$internal_obj)) {
                                 stop("No model summary + trace object ('internal_obj') available")
                               }
                               ##
                               if (debugging) {
                                 cat("Summary object structure:\n")
                                 str(self$internal_obj)
                               }
                               ##
                               network <- self$basic_model_options$network
                               ##
                               if (network) { 
                                 n_index_tests <- self$internal_obj$outs_data$stan_data_list$n_index_tests
                                 print(paste("n_index_tests = ", n_index_tests))
                                 ##
                                 n_thr         <- self$internal_obj$outs_data$stan_data_list$n_thr
                                 print(paste("n_thr = ", n_thr))
                                 ##
                               } else { 
                                 n_index_tests <- 1     ## dummy arg - needed for compatability / R6 class
                                 n_thr         <- NULL  ## dummy arg - needed for compatability / R6 class
                               }
                               ##
                               n_covariates <- self$internal_obj$outs_data$stan_data_list$n_covariates_max
                               print(paste("n_covariates = ", n_covariates))
                               ##
                               if (n_covariates == 1) {
                                 covariates <- FALSE
                               } else {
                                 covariates <- TRUE
                               }
                               ##
                               if (covariates) { 
                                 if (is.null(new_cov_data)) { 
                                   stop("You must supply 'new_cov_data' (an R list containing 'baseline_case_nd' and 
                                                        'baseline_case_d' for Xu/bivariate model or Jones model, or an R list 
                                                        containing 'baseline_case' for R&G model")
                                 }
                                 stan_mod_samples_to_use <- R_fn_GQ_covariates_NMA(
                                   baseline_case_nd = new_cov_data$baseline_case_nd,
                                   baseline_case_d = new_cov_data$baseline_case_d,
                                   model_summary_and_trace_obj = self,
                                   model_prep_obj = self
                                 )
                                 ##
                                 trace_gq <- stan_mod_samples_to_use$draws()
                               } else { 
                                 stan_mod_samples_to_use <- self$internal_obj$outs_stan_sampling$stan_mod_samples
                                 ##
                                 trace_gq <- self$internal_obj$traces$traces_as_arrays$trace_gq
                               }
                               # Use fast Rcpp version
                               outs <-  R_fn_using_Rcpp_compute_NMA_performance_posthoc(  trace_gq = trace_gq,
                                                                                                       n_thr = n_thr,
                                                                                                       test_names = test_names)
                               ##
                               return(outs)
                     
                   },
                   ##
                   ## --------- Meta-regression only: recompute statistics (Se/Sp/AUC) for new baseline:
                   ##
                   #'@param debugging ...
                   #'@param test_names ...
                   #'@param MR_scenarios ...
                   #'@param use_probit_link ..
                   #'@return new_baseline_outs
                   #'
                   recompute_baseline = function(  debugging = FALSE,
                                                   test_names,
                                                   MR_scenarios = NULL,
                                                   use_probit_link = TRUE
                   ) {
                             
                               if (debugging) {
                                 cat("Summary object structure:\n")
                                 str(self$internal_obj)
                               }
                               if (is.null(MR_scenarios)) { 
                                 stop("Must provide the R list 'MR_scenarios'")
                               }
                               ##
                               network <- self$basic_model_options$network
                               print(paste("network = ", network))
                               ##
                               if (network) { 
                                 n_index_tests <- self$internal_obj$outs_data$stan_data_list$n_index_tests
                                 print(paste("n_index_tests = ", n_index_tests))
                                 ##
                                 n_thr         <- self$internal_obj$outs_data$stan_data_list$n_thr
                                 print(paste("n_thr = ", n_thr))
                                 ##
                               } else { 
                                 n_index_tests <- 1     ## dummy arg - needed for compatability / R6 class
                                 n_thr         <- NULL  ## dummy arg - needed for compatability / R6 class
                               }
                               ##
                               n_covariates <- self$internal_obj$outs_data$stan_data_list$n_covariates_max
                               print(paste("n_covariates = ", n_covariates))
                               ##
                               if (n_covariates == 1) {
                                 stop("The 'recompute_baseline' R6 method only works for meta-regression models!")
                               }
                               ## ---- If MR_scenarios provided, loop through them:
                               if ((!is.null(MR_scenarios)) && (n_covariates > 1)) {
                                   
                                     baseline_outs_list <- list()
                                     
                                     for (k in 1:length(MR_scenarios)) {
                                       
                                           baseline_vec_k <- c(1, unlist(MR_scenarios[[k]]))
                                           baseline_case_nd_k <- rep(list(baseline_vec_k), length(test_names))
                                           baseline_case_d_k <- rep(list(baseline_vec_k), length(test_names))
                                             
                                           baseline_outs_list[[k]] <- R_fn_using_Rcpp_compute_MR_Se_Sp_AUC_baseline(
                                                                                 model_prep_obj = self,
                                                                                 model_summary_and_trace_obj = self,
                                                                                 baseline_case_nd = baseline_case_nd_k,
                                                                                 baseline_case_d = baseline_case_d_k,
                                                                                 use_probit_link = use_probit_link,
                                                                                 test_names = test_names)
                                           ##
                                           # ADD THIS: Store the original scenario definition
                                           baseline_outs_list[[k]]$original_scenario <- MR_scenarios[[k]]
                                         
                                     }
                                     
                                     names(baseline_outs_list) <- names(MR_scenarios)
                                     self$baseline_outs_MR <- baseline_outs_list  # Store for extract_AUC
                                     ##
                                     return(baseline_outs_list)
                            
                               } 
                               
                   },
                   ##
                   ## --------- Meta-regression only: recompute statistics (Se/Sp/AUC) for new baseline:
                   ##
                   #'@param debugging ...
                   #'@param test_names ...
                   #'@param MR_scenarios ...
                   #'@param use_probit_link ..
                   #'@return new_baseline_outs
                   #'
                   compute_pairwise_comparisons = function(debugging = FALSE, 
                                                           test_names,
                                                           MR_scenarios = NULL,
                                                           recompute_baseline = TRUE,
                                                           use_probit_link = TRUE
                   ) {
                     
                             n_covariates <- self$internal_obj$outs_data$stan_data_list$n_covariates_max
                             ##
                             if (n_covariates == 1) {
                               stop("Only works for meta-regression models!")
                             }
                             ##
                             n_thr <- self$internal_obj$outs_data$stan_data_list$n_thr
                             ##
                             if (recompute_baseline) {  ## OPTION 1: User provides MR_scenarios - compute fresh
                               
                                     if (is.null(MR_scenarios)) { 
                                       stop("If recompute_baseline == TRUE, you must provide 'MR_scenarios' R list")
                                     }
                                     ##
                                     if (debugging) cat("Computing Se/Sp arrays fresh from MR_scenarios...\n")
                                     ##
                                     Se_Sp_list <- list()
                                     ##
                                     for (k in 1:length(MR_scenarios)) {
                                       baseline_vec_k <- c(1, unlist(MR_scenarios[[k]]))
                                       baseline_case_nd_k <- rep(list(baseline_vec_k), length(test_names))
                                       baseline_case_d_k <- rep(list(baseline_vec_k), length(test_names))
                                       
                                       Se_Sp_list[[k]] <- R_fn_using_Rcpp_compute_MR_Se_Sp_AUC_baseline(
                                         model_prep_obj = self,
                                         model_summary_and_trace_obj = self,
                                         baseline_case_nd = baseline_case_nd_k,
                                         baseline_case_d = baseline_case_d_k,
                                         use_probit_link = use_probit_link,
                                         test_names = test_names
                                       )
                                     }
                                     names(Se_Sp_list) <- names(MR_scenarios)
                        
                             } else if (!is.null(self$baseline_outs_MR)) { ## OPTION 2: Use stored results from recompute_baseline
                               
                                     if (debugging) cat("Using stored Se/Sp arrays from recompute_baseline...\n")
                                     Se_Sp_list <- self$baseline_outs_MR
                             
                             } else { ## ERROR: Neither option available
                               
                                     stop("Must either provide 'MR_scenarios' OR run 'recompute_baseline()' first")
                               
                             }
                             ## ---- Now compute pairwise comparisons from the Se/Sp arrays (SHARED CODE)
                             pairwise_comparisons_list <- list()
                             ##
                             for (k in 1:length(Se_Sp_list)) {
                                     
                                     Se_baseline <- Se_Sp_list[[k]]$Se_baseline
                                     Sp_baseline <- Se_Sp_list[[k]]$Sp_baseline
                                     array_dims <- as.integer(dim(Se_baseline))
                                     
                                     # Call Rcpp function
                                     rcpp_result <- MetaOrdDTA:::Rcpp_compute_NMA_comparisons_from_arrays(
                                       Se_array = Se_baseline,
                                       Sp_array = Sp_baseline,
                                       array_dims = array_dims,
                                       n_thr = n_thr
                                     )
                                     
                                     # Add metadata - CRITICAL for plotting!
                                     comp_metadata <- generate_comparison_metadata(n_thr, test_names)
                                     
                                     # Combine into tibble
                                     pairwise_comparisons_list[[k]] <- tibble(cbind(comp_metadata,
                                                                                    rcpp_result$summaries))
                               
                             }
                             
                             
                             # # Now compute pairwise comparisons from the Se/Sp arrays
                             # pairwise_comparisons_list <- list()
                             # 
                             # for (k in 1:length(Se_Sp_list)) {
                             #         
                             #         Se_baseline <- Se_Sp_list[[k]]$Se_baseline
                             #         Sp_baseline <- Se_Sp_list[[k]]$Sp_baseline
                             #         array_dims <- as.integer(dim(Se_baseline))
                             #         
                             #         pairwise_comparisons_list[[k]] <- MetaOrdDTA:::Rcpp_compute_NMA_comparisons_from_arrays(
                             #           Se_array = Se_baseline,
                             #           Sp_array = Sp_baseline,
                             #           array_dims = array_dims,
                             #           n_thr = n_thr)
                             #  
                             # }
                             
                             names(pairwise_comparisons_list) <- names(Se_Sp_list)
                             self$pairwise_comparisons_MR <- pairwise_comparisons_list
                             
                             return(pairwise_comparisons_list)
                     
                   },
                   ##
                   ## --------- Extract AUC summaries:
                   ##
                   #'@param debugging ...
                   #'@param test_names ...
                   #'@param use_BayesMVP_for_faster_summaries ...
                   #'@return AUC summary
                   #'
                   extract_AUC = function( debugging = FALSE,
                                           test_names,
                                           use_BayesMVP_for_faster_summaries = NULL
                   ) {
                     
                               if (is.null(self$internal_obj)) {
                                 stop("No model samples object ('internal_obj') available")
                               }
                               if (is.null(self$internal_obj)) {
                                 stop("No model summary + trace object ('internal_obj') available")
                               }
                               ##
                               if (debugging) {
                                 cat("Summary object structure:\n")
                                 str(self$internal_obj)
                               }
                               ##
                               network <- self$basic_model_options$network
                               print(paste("network = ", network))
                               ##
                               if (network) { 
                                 n_index_tests <- self$internal_obj$outs_data$stan_data_list$n_index_tests
                                 print(paste("n_index_tests = ", n_index_tests))
                                 ##
                                 n_thr         <- self$internal_obj$outs_data$stan_data_list$n_thr
                                 print(paste("n_thr = ", n_thr))
                                 ##
                               } else { 
                                 n_index_tests <- 1     ## dummy arg - needed for compatability / R6 class
                                 n_thr         <- NULL  ## dummy arg - needed for compatability / R6 class
                               }
                               ##
                               n_covariates <- self$internal_obj$outs_data$stan_data_list$n_covariates_max
                               print(paste("n_covariates = ", n_covariates))
                               ##
                               if (n_covariates == 1) {
                                 covariates <- FALSE
                                 use_MR_scenarios <- FALSE
                               } else {
                                 covariates <- TRUE
                                 use_MR_scenarios <- TRUE
                                 if (is.null(self$baseline_outs_MR)) {
                                   stop("No MR scenarios computed. Run compute_MR_scenarios() first, then try running '$extract_AUC()' again.")
                                 }
                               }
                               # ##
                               if (covariates) { 
                                         
                                               # Extract AUC from all stored scenarios
                                               outs_AUC <- list()
                                               outs_AUC$scenarios <- names(self$baseline_outs_MR)
                                               outs_AUC$auc_by_scenario <- lapply(self$baseline_outs_MR, function(x) x$AUC_summary)
                                               outs_AUC$auc_pred_by_scenario <- lapply(self$baseline_outs_MR, function(x) x$AUC_pred_summary)
                                               outs_AUC$auc_diff_by_scenario <- lapply(self$baseline_outs_MR, function(x) x$AUC_diff_summary)
                                 
                               } else {  ## No covariates case
                                 
                                         tibble_gq <- self$internal_obj$summaries$summary_tibbles$summary_tibble_generated_quantities
                                         ##
                                         outs_AUC <- R_fn_extract_NMA_AUC_summary( tibble_gq = tibble_gq,
                                                                                   network = network,
                                                                                   n_thr = n_thr,
                                                                                   test_names = test_names)
                                 
                                 
                               }
                               ##
                               outs <- list(outs_AUC = outs_AUC)
                               ##
                               return(outs)
                     
                   },
                   ##
                   ## --------- Extract Se, Sp summaries:
                   ##
                   #'@param debugging ...
                   #'@param test_names ...
                   #'@param use_BayesMVP_for_faster_summaries ...
                   #'@return Se/Sp summary
                   #'
                   extract_MR_Se_Sp = function( debugging = FALSE,
                                                test_names,
                                                use_BayesMVP_for_faster_summaries = NULL
                   ) {
                     
                     if (is.null(self$internal_obj)) {
                       stop("No model samples object ('internal_obj') available")
                     }
                     if (is.null(self$internal_obj)) {
                       stop("No model summary + trace object ('internal_obj') available")
                     }
                     ##
                     if (debugging) {
                       cat("Summary object structure:\n")
                       str(self$internal_obj)
                     }
                     ##
                     network <- self$basic_model_options$network
                     print(paste("network = ", network))
                     ##
                     if (network) { 
                       n_index_tests <- self$internal_obj$outs_data$stan_data_list$n_index_tests
                       print(paste("n_index_tests = ", n_index_tests))
                       ##
                       n_thr         <- self$internal_obj$outs_data$stan_data_list$n_thr
                       print(paste("n_thr = ", n_thr))
                       ##
                     } else { 
                       n_index_tests <- 1     ## dummy arg - needed for compatability / R6 class
                       n_thr         <- NULL  ## dummy arg - needed for compatability / R6 class
                     }
                     ##
                     n_covariates <- self$internal_obj$outs_data$stan_data_list$n_covariates_max
                     print(paste("n_covariates = ", n_covariates))
                     ##
                     if (n_covariates == 1) {
                       covariates <- FALSE
                       use_MR_scenarios <- FALSE
                     } else {
                       covariates <- TRUE
                       use_MR_scenarios <- TRUE
                       if (is.null(self$baseline_outs_MR)) {
                         stop("No MR scenarios computed. Run compute_MR_scenarios() first, then try running '$extract_AUC()' again.")
                       }
                     }
                     # ##
                     if (covariates) { 
                           
                           # Extract AUC from all stored scenarios
                           outs_Se_Sp <- list()
                           outs_Se_Sp$scenarios <- names(self$baseline_outs_MR)
                           outs_Se_Sp$Se_by_scenario <- lapply(self$baseline_outs_MR, function(x) x$Se_summaries)
                           outs_Se_Sp$Sp_by_scenario <- lapply(self$baseline_outs_MR, function(x) x$Sp_summaries)
                           outs_Se_Sp$Se_pred_by_scenario <- lapply(self$baseline_outs_MR, function(x) x$Se_pred_summaries)
                           outs_Se_Sp$Sp_pred_by_scenario <- lapply(self$baseline_outs_MR, function(x) x$Sp_pred_summaries)
                       
                     } else {  ## No covariates case
                       
                          stop("This function is only relevant for meta-regression")
                       
                     }
                     ##
                     outs <- list(outs_Se_Sp = outs_Se_Sp)
                     ##
                     return(outs)
                     
                   }
                   
  

  
  


)
)







