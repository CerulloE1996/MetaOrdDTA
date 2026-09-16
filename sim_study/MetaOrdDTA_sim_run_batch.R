# ============================================================================================================
# Batch Runner for MetaOrdDTA Simulations
# ============================================================================================================

run_batch_MetaOrd_simulations <- function(scenario,
                                          config,
                                          n_sims_total,
                                          n_workers,
                                          sims_per_batch,
                                          local_pkg_dir,
                                          output_dir,
                                          starting_seed = 1,
                                          stop_early_if_target_met = TRUE) {
  
        # Pre-compile Stan models in main process
        precompile_stan_models(scenario, config, local_pkg_dir)
        
        # Calculate batches
        n_batches <- ceiling(n_sims_total / sims_per_batch)
        all_batch_results <- list()
        
        for (batch_id in 1:n_batches) {
          
          cat(sprintf("\n========== BATCH %d/%d ==========\n", batch_id, n_batches))
          
          # Check MCSE status before batch
          if (batch_id > 1) {
            mcse_status <- check_MCSE_status_MetaOrd(
              output_dir = output_dir,
              file_name_string = create_file_name_string(scenario),
              target_MCSE_RMSE = config$target_MCSE_RMSE,
              vec_index_inner = config$threshold_configs[[scenario$index_test_MA]]$vec_index_inner,
              true_DGM_Se = true_DGM_Se,
              true_DGM_Sp = true_DGM_Sp
            )
            
            if (mcse_status$mcse_met && stop_early_if_target_met) {
              cat("✓ TARGET MCSE REACHED! Stopping early.\n")
              break
            }
          }
          
          # Calculate seed range for this batch
          start_seed <- starting_seed + (batch_id - 1) * sims_per_batch
          end_seed <- min(starting_seed + batch_id * sims_per_batch - 1, 
                          starting_seed + n_sims_total - 1)
          batch_seeds <- start_seed:end_seed
          
          # Distribute seeds to workers
          seeds_per_worker <- split(batch_seeds, 
                                    cut(seq_along(batch_seeds), n_workers, labels = FALSE))
          
          # Create cluster
          cl <- makeCluster(n_workers, type = "PSOCK", outfile = "")
          registerDoParallel(cl)
          
          # Export necessary objects
          clusterExport(cl, c("scenario", "config", "local_pkg_dir", "output_dir"),
                        envir = environment())
          
          # Run parallel simulations
          batch_results <- foreach(worker_id = 1:n_workers,
                                   .combine = 'c',
                                   .errorhandling = "pass",
                                   .packages = c("MetaOrdDTA", "cmdstanr", "dplyr")) %dopar% {
                                     
                                     # Set up worker environment
                                     cmdstanr::set_cmdstan_path()
                                     worker_seeds <- seeds_per_worker[[worker_id]]
                                     worker_results <- list()
                                     
                                     for (seed in worker_seeds) {
                                       result <- run_single_MetaOrd_simulation(
                                         seed = seed,
                                         scenario = scenario,
                                         config = config,
                                         local_pkg_dir = local_pkg_dir,
                                         output_dir = output_dir
                                       )
                                       worker_results[[as.character(seed)]] <- result
                                     }
                                     
                                     return(worker_results)
                                   }
          
          # Stop cluster
          stopCluster(cl)
          
          # Store batch results
          all_batch_results[[batch_id]] <- batch_results
          
          # Clean up
          gc(reset = TRUE, full = TRUE)
          
          cat(sprintf("Batch %d completed\n", batch_id))
        }
        
        return(all_batch_results)
  
}




# Pre-compile Stan models
precompile_stan_models <- function(scenario, 
                                   config, 
                                   local_pkg_dir) {
        
        cat("Pre-compiling Stan models...\n")
        
        # Create dummy data for compilation
        x_dummy <- list(
          matrix(1, nrow = 5, ncol = 5),
          matrix(1, nrow = 5, ncol = 5)
        )
        
        # Compile model
        dummy_model <- MetaOrdDTA::MetaOrd_model$new(
          x = x_dummy,
          model_parameterisation = scenario$model_parameterisation,
          random_thresholds = scenario$random_thresholds,
          n_chains = 1
        )
        
        cat("Models pre-compiled successfully\n")
  
}




