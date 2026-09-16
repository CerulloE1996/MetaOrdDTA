# ============================================================================================================
# Run Single MetaOrdDTA Simulation
# ============================================================================================================

run_single_MetaOrd_simulation <- function(seed,
                                          scenario,
                                          config,
                                          local_pkg_dir,
                                          output_dir,
                                          force_autodiff_for_metric = TRUE,
                                          check_mcse_every = 5) {
  
        set.seed(seed)
        
        # Extract scenario parameters
        list2env(scenario, environment())
        list2env(config, environment())
        
        # Get threshold configuration
        threshold_config <- config$threshold_configs[[index_test_MA]]
        vec_index_inner <- threshold_config$vec_index_inner
        n_thr <- threshold_config$n_thr
        
        # Load data
        source(file.path(local_pkg_dir, "R", "R_fn_sim_data_ord_MA_v2.R"), local = TRUE)
        source(file.path(local_pkg_dir, "inst", "examples", "MA_missing_thr_prep_data.R"), local = TRUE)
        
        # Select appropriate dataset
        if (missing_thr == TRUE) {
          x <- get(paste0("x_", gsub("_", "", index_test_MA), "_w_missing_thr"))
        } else {
          x <- get(paste0("x_", gsub("_", "", index_test_MA)))
        }
        
        n_studies <- nrow(x[[1]])
        
        # Initialize model
        model_prep_obj <- MetaOrdDTA::MetaOrd_model$new(
          debugging = FALSE,
          ##
          x = x,
          intercept_only = TRUE,
          ##
          n_chains = config$n_chains,
          ##
          network = FALSE,
          ##
          model_parameterisation = model_parameterisation,
          random_thresholds = random_thresholds,
          Dirichlet_random_effects_type = Dirichlet_random_effects_type
        )
        
        # Set up priors
        priors <- setup_priors_MetaOrd(
          model_prep_obj = model_prep_obj,
          alpha_pi = alpha_pi,
          n_thr = n_thr,
          model_parameterisation = model_parameterisation
        )
        
        # Set up initial values
        init_lists_per_chain <- setup_inits_MetaOrd(
          model_parameterisation = model_parameterisation,
          n_chains = config$n_chains,
          n_thr = n_thr
        )
        
        # Sample model
        n_iter <- round(config$n_iter_base / (config$n_chains / 4))
        
        model_samples_obj <- model_prep_obj$sample(
          n_burnin = config$n_burnin,
          n_iter = n_iter,
          adapt_delta = config$adapt_delta,
          max_treedepth = config$max_treedepth,
          metric_shape = "diag_e",
          priors = priors,
          init_lists_per_chain = init_lists_per_chain,
          n_chains = config$n_chains
        )
        
        # Get summaries
        model_summary_obj <- model_samples_obj$summary(
          compute_main_params = TRUE,
          compute_generated_quantities = TRUE,
          use_BayesMVP_for_faster_summaries = TRUE
        )
        
        # Extract results
        results <- extract_results_MetaOrd(
          model_summary_obj = model_summary_obj,
          seed = seed,
          scenario = scenario,
          vec_index_inner = vec_index_inner
        )
        
        # Check ESS
        min_ESS <- min(results$tibble_gq$n_eff, na.rm = TRUE)
        
        if (min_ESS > config$min_ESS_threshold) {    # Save results
             
                save_path <- file.path(output_dir,
                                       paste0("seed_", seed, "_",
                                              create_file_name_string(scenario),
                                              "_results_list.RDS"))
                saveRDS(results, save_path)
                
                return(list(success = TRUE, min_ESS = min_ESS))
          
        } else {
          
                return(list(success = FALSE, min_ESS = min_ESS))
          
        }
  
}























