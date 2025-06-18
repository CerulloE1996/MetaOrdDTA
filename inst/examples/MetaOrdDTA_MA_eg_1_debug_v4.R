






{
    
    os <- .Platform$OS.type
    
    
    if (os == "unix") {
      user_root_dir <- Sys.getenv("PWD")
    } else if (os == "windows") {
      user_root_dir <- Sys.getenv("USERPROFILE")
    }
    local_pkg_dir <- file.path(user_root_dir, "Documents/Work/PhD_work/R_packages/MetaOrdDTA")
    
    
    setwd(local_pkg_dir)
    
    {
      require(MetaOrdDTA)
      require(RcppParallel)
      require(BayesMVP)
      require(dplyr)
    }
    
    source(file.path(getwd(), "inst", "examples", "thr_sim_study_helper_fns.R"))
    
    
    
    
    sim_index <- 1
    # results_list <- c()
    
    seed <- 1
    
    
    # total_target_sims <- n_sims  # Your final target
    # 
    # # Calculate how many more sims to run
    # n_sims <- total_target_sims - start_sim_index + 1  # This would be 470 in your example
    
    iii <- 1
    
 
    file_name_string <- paste0( "test_", index_test_MA, 
                                "N_", N_per_study_mean,
                                "n_studies_", n_studies, 
                                ##
                                "param_", model_parameterisation, 
                                "soft_", softplus,
                                "rand_thr_", random_thresholds,
                                "alpha_pi_", alpha_pi
    )
    
    results_list_file_name <- paste0(file_name_string, "_results_list.RDS")
    # summary_df_file_name   <- paste0(file_name_string, "_summary_df.RDS")
    
    results_list_path <- file.path(output_dir, results_list_file_name)
    # summary_df_path   <- file.path(output_dir, summary_df_file_name)
    
    # Store the total target before we modify n_sims
    total_target_sims <- n_sims  
    
    # Calculate how many more sims to run
    n_sims_to_run <- total_target_sims - start_sim_index + 1
    
    # Before the simulation loop, initialize results_list
    if (start_sim_index == 1) {
      # Fresh start
      results_list <- vector("list", total_target_sims)
      # summary_df <- create_summary_df(n_sims = total_target_sims,
      #                                 min_k = vec_index_inner[1], 
      #                                 max_k = tail(vec_index_inner, 1),
      #                                 start_index = 1)
    } else {
      
      # Initialize with proper size
      results_list <- vector("list", total_target_sims)
      # summary_df <- create_summary_df(n_sims = total_target_sims,
      #                                 min_k = vec_index_inner[1], 
      #                                 max_k = tail(vec_index_inner, 1),
      #                                 start_index = 1)
      
      # Load individual seed results that were already completed
      cat("Loading previous results...\n")
      n_loaded <- 0
      
      for (s in 1:(start_sim_index - 1)) {
        
                seed_file <- file.path(output_dir, paste0("seed_", s, "_", results_list_file_name))
                ##
                if (file.exists(seed_file)) {
                  results_list[[s]] <- readRDS(seed_file)
                  
                  if (!is.null(results_list[[s]])) {
                    # # Create row data
                    # row_data <- list(
                    #   seed = s,  # or results_list[[s]]$seed
                    #   # sim_index = s,  
                    #   n_studies = results_list[[s]]$n_studies,
                    #   model_parameterisation = results_list[[s]]$model_parameterisation,
                    #   box_cox = results_list[[s]]$box_cox,
                    #   softplus = results_list[[s]]$softplus,
                    #   random_thresholds = results_list[[s]]$random_thresholds
                    # )
                    # ##
                    # # Add Se and Sp values
                    # for (k in vec_index_inner) {
                    #   row_data[[paste0("Se_mean_", k)]] <- results_list[[s]]$Se_median_vec[k]
                    # }
                    # for (k in vec_index_inner) {
                    #   row_data[[paste0("Sp_mean_", k)]] <- results_list[[s]]$Sp_median_vec[k]
                    # }
                    ##
                    # Update the row using s (which should equal sim_index for loaded data)
                    # summary_df[s, ] <- row_data
                    # n_loaded <- n_loaded + 1
                  }
                }
        
      }
      
      cat(sprintf("Loaded %d previous results. Starting from simulation %d.\n", n_loaded, start_sim_index))
    }



}










{
    
  
# Run simulations with different seeds
for (iii in 1:n_sims_to_run) {
  
  sim_index <- start_sim_index + iii - 1

  
  require(MetaOrdDTA)
  
  seed <- sim_index
  # results_list_file_w_seed <- file.path(output_dir, paste0("seed_", seed, "_", results_list_file))
  # results_list_file_w_seed <- file.path(output_dir, paste0("seed_", seed, "_", results_list_file))
  ##
  cat(sprintf("\n\n--- Running simulation %d/%d with seed = %d ---\n\n", sim_index, n_sims, seed))
  
  # Set seed for reproducibility
  set.seed(seed)
  
  
  {
    # ##
    # ## ---- Load NMA data:
    # ##
    setwd(local_pkg_dir)
    ##
    source(file.path(getwd(), "R", "R_fn_sim_data_ord_MA.R"))
    source(file.path(getwd(), "inst", "examples", "NMA_missing_thr_prep_data.R"))
    ##
    # data <- readRDS(file.path(getwd(), "inst", "examples", "data_example_1_NMA_list.RDS"))
    ##
    # x_NMA <- data$x_NMA
    # ##
    # x_MA <- list()
    # for (c in 1:2) {
    #   x_MA[[c]] <- x_NMA[[c]][[index_test_chosen_index - 1]]
    # }
    # ##
    # ## Select which data ("x") to use:
    # ##
    # str(x_NMA)
    ##
    ## ---- Select x to use:
    ##
    if (index_test_MA == "GAD_7") {
      if (missing_thr == TRUE) x <- x_GAD7_w_missing_thr
      else                     x <- x_GAD7
    } else if (index_test_MA == "GAD_2") { 
      if (missing_thr == TRUE) x <- x_GAD2_w_missing_thr
      else                     x <- x_GAD2
    } else if (index_test_MA == "HADS") {
      if (missing_thr == TRUE) x <- x_HADS_w_missing_thr
      else                     x <- x_HADS
    }
    ##
    n_thr <- ncol(x[[1]]) - 1
    ##
    n_studies <- nrow(x[[1]])
    n_studies
    ##
    # x <- arrange_missing_values(x)
  }
  
  
  try({  
  

  ## true_DGM_Fp

##
## ----  Initialise / select model: --------------------------------------------------------------------------- NMA:
##

# internal_obj$outs_data$stan_data_list$x
# internal_obj$outs_data$stan_data_list$n
# 
# any(is.na(unlist(internal_obj$outs_data$stan_data_list$n)))
init_lists_per_chain <- NULL
##
##
##
##


if (model_parameterisation == "Bivariate") {

    custom_file_name <- "DTA_MA_strat_bivariate.stan"
    
} else { 
  {
    if (use_custom_file == TRUE) custom_file_name <- file_name
    else  custom_file_name <- NULL
  }
  
}



# BASE_FLAGS <- "-O3  -march=native  -mtune=native -mfma -mavx -mavx2   -mavx512vl -mavx512dq  -mavx512f -fno-math-errno  -fno-signed-zeros -fno-trapping-math -fPIC -DNDEBUG -fpermissive  -DBOOST_DISABLE_ASSERTS -Wno-deprecated-declarations -Wno-sign-compare -Wno-ignored-attributes -Wno-class-memaccess -Wno-class-varargs"
# 
# Sys.setenv(PKG_CPPFLAGS = BASE_FLAGS)
# Sys.setenv(PKG_CXXFLAGS = BASE_FLAGS)
# Sys.setenv(CPPFLAGS = BASE_FLAGS)
# # Sys.setenv(CXXFLAGS = BASE_FLAGS)
# Sys.unsetenv("PKG_CPPFLAGS")
# Sys.unsetenv("PKG_CXXFLAGS")
# Sys.unsetenv("CPPFLAGS")
# Sys.unsetenv("CXXFLAGS")
##
model_prep_obj <- MetaOrdDTA::MetaOrd_model$new(  
                  debugging = TRUE,
                  ##
                  x = x, 
                  indicator_index_test_in_study = NULL,
                  # n_index_tests_per_study = NULL,
                  ##
                  # priors = priors,
                  ##
                  n_chains = n_chains,
                  ##
                  cts = cts,
                  network = network,
                  ##
                  prior_only = FALSE,
                  ##
                  softplus = softplus,
                  ##
                  box_cox = box_cox,
                  ##
                  model_parameterisation = model_parameterisation,
                  random_thresholds = random_thresholds,
                  Dirichlet_random_effects_type = Dirichlet_random_effects_type, ## only used if random cutpoints
                  ##
                  init_lists_per_chain = init_lists_per_chain,
                  ##
                  advanced_compile = TRUE, ## default is standard/"basic" compile
                  ##
                  force_recompile = F,
                  ##
                  set_custom_CXX_CPP_flags = TRUE, 
                  CCACHE_PATH = "/usr/bin/ccache", 
                  custom_cpp_user_header_file_path = NULL, ## if wish to include custom C++ files
                  CXX_COMPILER_PATH = "/opt/AMD/aocc-compiler-5.0.0/bin/clang++", ## g++
                  CPP_COMPILER_PATH = "/opt/AMD/aocc-compiler-5.0.0/bin/clang", ## gcc
                  # CXX_COMPILER_PATH = "g++",
                  # CPP_COMPILER_PATH = "gcc",
                  MATH_FLAGS  = "-fno-math-errno  -fno-signed-zeros -fno-trapping-math", ## fine
                  FMA_FLAGS = "-mfma", ## NOT fine !!!#
                  # AVX_FLAGS =  "-mavx", ## NOT working
                  AVX_FLAGS = "-mavx -mavx2 -mavx512vl -mavx512dq  -mavx512f",
                  THREAD_FLAGS = "-pthread -D_REENTRANT", ## fine
                  ##
                  # compute_sim_study_metrics = compute_sim_study_metrics,
                  # vec_index_inner_thr = vec_index_inner,
                  ##
                  custom_file_name = custom_file_name
                  )


{
      # model_prep_obj$internal_obj$outs_data$stan_data_list
      init_lists_per_chain <- model_prep_obj$init_lists_per_chain
      priors <- model_prep_obj$priors
      ##
      n_studies <- nrow(x[[1]])
      ##  n_index_tests <- 1# length(x[[1]])
      ##
      stan_data_list <- model_prep_obj$internal_obj$outs_data$stan_data_list
}




      
      if (alpha_pi == "flat") {
        
              prior_dirichlet_alpha <- rep(2.0, n_thr + 1)
              ##
              if (model_parameterisation == "Xu") {
                 priors$prior_dirichlet_alpha <- list( prior_dirichlet_alpha, prior_dirichlet_alpha)
              } else {
                 priors$prior_dirichlet_alpha <- prior_dirichlet_alpha
              }
              
      } else {
        
              if (model_parameterisation == "Xu") {
              
                        alpha_priors_outs <- MetaOrdDTA:::generate_dual_ordinal_dirichlet_priors(n_cat = n_thr + 1,
                                                                                    nd_peak = 0.15,
                                                                                    d_peak = 0.30,
                                                                                    nd_max_alpha = 5,
                                                                                    d_max_alpha = 5,
                                                                                    nd_min_alpha = 1,
                                                                                    d_min_alpha = 1,
                                                                                    smoothness = 0.15)
            
                        priors$prior_dirichlet_alpha <- list( alpha_priors_outs$non_diseased, 
                                                              alpha_priors_outs$diseased)
              } else { 
                
                        priors$prior_dirichlet_alpha <- generate_ordinal_dirichlet_prior( n_cat = n_thr + 1,
                                                                                          max_alpha = 5,
                                                                                          min_alpha = 1,
                                                                                          left_skew_factor = 0.5)
              }
        
      }




 


priors$use_probit_link <- 1


priors

 stan_data_list$cutpoint_index[1,,]
# 
stan_data_list$x[1,,]


if (model_parameterisation == "Bivariate") {
  
      init_lists_per_chain <- list()
      mu <- list()
      sigma <- list()
      for (k in 1:n_thr) { 
        mu[[k]] <- c(-1, +1)
        sigma[[k]] <- c(0.01, 0.01)
      }
      
      for (kk in 1:n_chains) { 
        init_lists_per_chain[[kk]] <- list()
      }
      
      for (kk in 1:n_chains) {
        init_lists_per_chain[[kk]]$mu <- mu
        init_lists_per_chain[[kk]]$sigma <- sigma
      }
  
} else {
  
      init_lists_per_chain <- NULL
  
}


##
## ----  Sample model: ------------------------------------------------------------------------------
##
model_samples_obj <-  model_prep_obj$sample(   
                             n_burnin = n_burnin,
                             n_iter   = n_iter,
                             adapt_delta = adapt_delta, 
                             max_treedepth = max_treedepth,
                             metric_shape = "diag_e",
                             ##
                             priors = priors,
                             init_lists_per_chain = init_lists_per_chain,
                             ##
                             n_chains = n_chains
                             )

 
##
## ----  Summarise + output results: ----------------------------------------------------------------
##
model_summary_and_trace_obj <- model_samples_obj$summary(
                                          compute_main_params = TRUE,
                                          compute_transformed_parameters = FALSE, 
                                          compute_generated_quantities = TRUE,
                                          ##
                                          save_log_lik_trace = TRUE,
                                          ##
                                          use_BayesMVP_for_faster_summaries = TRUE)

tibble_main <- model_summary_and_trace_obj$get_summary_main() %>% print(n = 200)
# tibble_tp <- model_summary_and_trace_obj$get_summary_transformed() %>% print(n = 100)
tibble_gq   <- model_summary_and_trace_obj$get_summary_generated_quantities() %>% print(n = 1000)
##
  tibble_all <- rbind(tibble_main, tibble_gq)
   # tibble_all <- rbind(tibble_main, tibble_tp, tibble_gq)
##
dplyr::filter(tibble_all, (stringr::str_detect(parameter, "Se"))) %>% print(n = 1000)
true_DGM_Se
mean(dplyr::filter(tibble_all, (stringr::str_detect(parameter, "Se_baseline\\[")))$mean*100 - true_DGM_Se)
##
dplyr::filter(tibble_all, (stringr::str_detect(parameter, "Sp"))) %>% print(n = 1000)
true_DGM_Sp
mean(dplyr::filter(tibble_all, (stringr::str_detect(parameter, "Sp_baseline\\[")))$mean*100 - true_DGM_Sp)
##
try({ 
    dplyr::filter(tibble_all, (stringr::str_detect(parameter, "kappa"))) %>% print(n = 1000)
    dplyr::filter(tibble_all, (stringr::str_detect(parameter, "dirichlet_phi"))) %>% print(n = 1000)
})
try({ 
  dplyr::filter(tibble_all, (stringr::str_detect(parameter, "alpha"))) %>% print(n = 1000)
})



try({ 
  dplyr::filter(tibble_all, (stringr::str_detect(parameter, "dirichlet_cat_SDs_sigma"))) %>% print(n = 1000)
  dplyr::filter(tibble_all, (stringr::str_detect(parameter, "target_sds"))) %>% print(n = 1000)
})

dplyr::filter(tibble_all, (stringr::str_detect(parameter, "C_raw"))) %>% print(n = 1000)


dplyr::filter(tibble_all, (stringr::str_detect(parameter, "beta_mu"))) %>% print(n = 1000)
dplyr::filter(tibble_all, (stringr::str_detect(parameter, "raw_scale_mu"))) %>% print(n = 1000)
dplyr::filter(tibble_all, (stringr::str_detect(parameter, "beta_SD"))) %>% print(n = 1000)
dplyr::filter(tibble_all, (stringr::str_detect(parameter, "raw_scale_SD"))) %>% print(n = 1000)

dplyr::filter(tibble_all, (stringr::str_detect(parameter, "kappa"))) %>% print(n = 1000)
dplyr::filter(tibble_all, (stringr::str_detect(parameter, "dirichlet_cat_SDs_sigma"))) %>% print(n = 1000)

# kappa <- 5
# phi <- rep(2.0, 7)
# ##
# 
# 
# require(truncnorm)
# quantile(truncnorm::rtruncnorm(n = 10000, a = 0.0, mean = 0.0, sd =  priors$prior_kappa_SD), probs = c(0.025, 0.50, 0.975))
# 
# alpha <- phi*kappa
# alpha_0 <- sum(alpha);
# dirichlet_cat_SDs_sigma <- sqrt((alpha * (alpha_0 - alpha) ) / ((alpha_0)^2 * (alpha_0 + 1.0))); 
# dirichlet_cat_SDs_sigma

 

tibble_all ## save this! 
as.numeric(object.size(tibble_all)) / 1024^2 ## < 1 MB

traces_all <- model_summary_and_trace_obj$get_all_traces_list()
trace <- traces_all$traces_as_arrays
trace_gq <- trace$trace_gq ## could save this trace 
as.numeric(object.size(trace_gq)) / 1024^2 ## ~ 10 MB



   model_summary_and_trace_obj$get_efficiency_metrics()$time_sampling
   model_summary_and_trace_obj$get_HMC_info()$stan_HMC_info_list$out_list_sampling$mean_L_sampling
   
   100 * model_summary_and_trace_obj$get_HMC_info()$HMC_diagnostic_info$divergences$n_divs /  (n_chains*(n_iter))
 



##
## 
## ----  Plots: ----------------------------------------------------------------------
# ##
# require(Cairo)
# options(bitmapType = "cairo")

# ##
## Try setting X11 as the device with specific font options
# ## (solves bug on linux for font type when usijng ggplot)
# Sys.setenv("R_GSCMD" = "gs")
# options(device = "X11")
# X11.options(family = "sans")
# ##
# options(bitmapType = "Xlib")
# 
# stan_model_file_name <- model_summary_and_trace_obj$internal_obj$outs_stan_model_name$stan_model_file_name
# stan_mod_samples <- model_summary_and_trace_obj$internal_obj$outs_stan_sampling$stan_mod_samples
# # ##
# R_fn_sROC_plot_MA(stan_model_file_name = stan_model_file_name,
#                                stan_mod_samples = stan_mod_samples)
# ##
# plots <- model_summary_and_trace_obj$plot_sROC()
# plots$plot_list[[1]]
# plots$plot_list[[2]]

# model_summary_and_trace_obj$plot_traces(param_string = "alpha")


{
        n_thr_max <- max(model_prep_obj$internal_obj$outs_data$n_thr)
        ##
        Se_median_vec <- Sp_median_vec <- Fp_median_vec <- c()
        Se_mean_vec   <- Sp_mean_vec   <- Fp_mean_vec   <- c()
        ##
        Se_lower_vec <- Sp_lower_vec <- Fp_lower_vec <- c()
        Se_upper_vec <- Sp_upper_vec <- Fp_upper_vec <- c()
        ##
        Se_pred_lower_vec <- Sp_pred_lower_vec <- Fp_pred_lower_vec <- c()
        Se_pred_upper_vec <- Sp_pred_upper_vec <- Fp_pred_upper_vec <- c()

        # str(Se_array)

        Se <- dplyr::filter(tibble_gq, (stringr::str_detect(parameter, "Se_baseline")) &
                                     (!(stringr::str_detect(parameter, "Se_baseline_pred"))))
        Sp <- dplyr::filter(tibble_gq, (stringr::str_detect(parameter, "Sp_baseline")) &
                              (!(stringr::str_detect(parameter, "Sp_baseline_pred"))))
        Fp  <- dplyr::filter(tibble_gq, (stringr::str_detect(parameter, "Fp_baseline")) &
                              (!(stringr::str_detect(parameter, "Fp_baseline_pred"))))
        ##
        Se_pred <- dplyr::filter(tibble_gq, (stringr::str_detect(parameter, "Se_baseline_pred")))
        Sp_pred <- dplyr::filter(tibble_gq, (stringr::str_detect(parameter, "Sp_baseline_pred")))
        Fp_pred <- dplyr::filter(tibble_gq, (stringr::str_detect(parameter, "Fp_baseline_pred")))

        counter <- 1
        for (k in 1:n_thr_max) {
          for (t in 1:1) {
              # test_vec[counter] <- t
              ## Posterior medians (of pooled estimates):
              Se_median_vec[k] <- 100 * Se$`50%`[counter]
              Sp_median_vec[k] <- 100 * Sp$`50%`[counter]
              Fp_median_vec[k] <- 100 * Fp$`50%`[counter]
              ## Posterior means (of pooled estimates):
              Se_mean_vec[k] <- 100 * Se$mean[counter]
              Sp_mean_vec[k] <- 100 * Sp$mean[counter]
              Fp_mean_vec[k] <- 100 * Fp$mean[counter]
              ## Posterior lower 95% (of pooled estimates):
              Se_lower_vec[k] <- 100 * Se$`2.5%`[counter]
              Sp_lower_vec[k] <- 100 * Sp$`2.5%`[counter]
              Fp_lower_vec[k] <- 100 * Fp$`2.5%`[counter]
              ## Posterior upper 95% (of pooled estimates):
              Se_upper_vec[k] <- 100 * Se$`97.5%`[counter]
              Sp_upper_vec[k] <- 100 * Sp$`97.5%`[counter]
              Fp_upper_vec[k] <- 100 * Fp$`97.5%`[counter]
              ## Posterior lower prediction 95% (of pooled estimates):
              Se_pred_lower_vec[k] <- 100 * Se_pred$`2.5%`[counter]
              Sp_pred_lower_vec[k] <- 100 * Sp_pred$`2.5%`[counter]
              Fp_pred_lower_vec[k] <- 100 * Fp_pred$`2.5%`[counter]
              ## Posterior upper prediction 95% (of pooled estimates):
              Se_pred_upper_vec[k] <- 100 * Se_pred$`97.5%`[counter]
              Sp_pred_upper_vec[k] <- 100 * Sp_pred$`97.5%`[counter]
              Fp_pred_upper_vec[k] <- 100 * Fp_pred$`97.5%`[counter]
              ##
              counter <- counter + 1
          }
        }
}
   
 
 
# Se %>% print(n = 100) # Se_abs_diffs_sim <- Sp_abs_diffs_sim <- list()
# Se_mean_of_diffs_sim <- Sp_mean_of_diffs_sim <- c()
# Se_sum_of_diffs_sim <- Sp_sum_of_diffs_sim <- c()
# Se_max_of_diffs_sim <- Sp_max_of_diffs_sim <- c()

 ESS_vec <- c(Se$n_eff, Sp$n_eff)
 min_ESS <- min(Se$n_eff)
 min_ESS
  

  ##------------------------
  {
     vec_index_all         <-  1:n_thr
        {
          # x[[1]][, 10]
          print(paste("seed = ", seed))
          print(paste("model_parameterisation = ", model_parameterisation))
          print(paste("box_cox = ", box_cox))
          print(paste("softplus = ", softplus))
          print(paste("n_studies = ", n_studies))
          # print(paste("bs_het_C = ", bs_het_C))
        }
  }


 
 if (min_ESS > 100) {
   
   # Store results for this simulation
   results_list[[sim_index]] <- list(
     seed = seed,
     n_studies = n_studies,
     model_parameterisation = model_parameterisation,
     box_cox = box_cox,
     softplus = softplus,
     random_thresholds = random_thresholds,
     ##
     tibble_all, # added this - has posterior means, medians, posterior SD, 2.5 and 97.5% quintiles, n_eff and Rhat
     ##
     true_DGM_Se = true_DGM_Se,
     true_DGM_Sp = true_DGM_Sp,
     ##
     Se_median_vec = Se_median_vec, # added this
     Sp_median_vec = Sp_median_vec, # added this
     ##
     Se_mean_vec = Se_mean_vec, # added this
     Sp_mean_vec = Sp_mean_vec, # added this
     ##
     Se_lower_vec = Se_lower_vec, # added this
     Sp_lower_vec = Sp_lower_vec, # added this
     ##
     Se_upper_vec = Se_upper_vec, # added this
     Sp_upper_vec = Sp_upper_vec, # added this
     ##
     Se_pred_lower_vec = Se_pred_lower_vec, # added this
     Sp_pred_lower_vec = Sp_pred_lower_vec, # added this
     ##
     Se_pred_upper_vec = Se_pred_upper_vec, # added this
     Sp_pred_upper_vec = Sp_pred_upper_vec # added this
   )
   
   # # Update summary dataframe
   # row_data <- list(
   #   seed = seed,
   #   n_studies = n_studies,
   #   model_parameterisation = model_parameterisation,
   #   box_cox = box_cox,
   #   softplus = softplus,
   #   random_thresholds = random_thresholds
   # )
   # 
   # for (k in 1:n_thr) {
   #   row_data[[paste0("Se_mean_", k)]] <- Se_median_vec[k]
   # }
   # for (k in 1:n_thr) {
   #   row_data[[paste0("Sp_mean_", k)]] <- Sp_median_vec[k]
   # }
   # 
   # summary_df[sim_index, ] <- row_data
   
   # Save with overwrite check
   results_list_for_current_seed <- paste0("seed_", seed, "_", results_list_file_name)
   results_list_for_current_seed_path <- file.path(output_dir, results_list_for_current_seed)
   
   ## Remove if exists
   if (file.exists(results_list_for_current_seed_path)) {
     cat(sprintf("Removing existing file: %s\n", results_list_for_current_seed_path))
     file.remove(results_list_for_current_seed_path)
   }
   
   ## Save
   cat(sprintf("Saving to: %s\n", results_list_for_current_seed_path))
   saveRDS(results_list[[sim_index]], results_list_for_current_seed_path)
   
   ## Verify
   if (!file.exists(results_list_for_current_seed_path)) {
     stop(sprintf("Failed to save file: %s", results_list_for_current_seed_path))
   } else {
     cat(sprintf("File saved successfully (size: %d bytes)\n", 
                 file.info(results_list_for_current_seed_path)$size))
   }
   
   cat(sprintf("\n--- Simulation %d/%d completed ---\n", seed, n_sims))
   
 }
 
  
    
  })
  
  } ##  -----------------------  end of N_{sims} loop  --------------------------------------------------------------------------------------


  
  # # Calculate summary statistics across ALL simulations
  # summary_stats <- list(
  #   n_sims = total_target_sims,  # This should be 5 in your case
  #   n_studies = n_studies,
  #   model_parameterisation = model_parameterisation,
  #   box_cox = box_cox,
  #   softplus = softplus,
  #   random_thresholds = random_thresholds,
  #   summary_df = summary_df
  # )
  
  # Also rebuild the complete results_list from individual seed files
  complete_results_list <- vector("list", total_target_sims)
  ##
  for (seed_idx in 1:total_target_sims) {
    # Use the same filename pattern as when you save individual seeds
    seed_file <- file.path(output_dir, paste0("seed_", seed_idx, "_", file_name_string, "_results_list.RDS"))
    if (file.exists(seed_file)) {
      complete_results_list[[seed_idx]] <- readRDS(seed_file)
    } else {
      cat(sprintf("Warning: seed file not found: %s\n", seed_file))
    }
  }

  # # Save overall summary
  # try({ 
  #   summary_file <- file.path(output_dir, paste0(file_name_string, ".RDS"))
  #   saveRDS(summary_stats, summary_file)
  # })
  
  # Save the complete results list
  try({ 
    results_list_file <- file.path(output_dir, paste0(file_name_string, "_complete_results_list.RDS"))
    saveRDS(complete_results_list, results_list_file)
    cat(sprintf("Saved complete results list with %d non-NULL entries to: %s\n", 
                sum(!sapply(complete_results_list, is.null)), results_list_file))
  })
  
  # # Save as CSV for easy viewing
  # try({ 
  #   summary_csv <- file.path(output_dir, paste0(file_name_string, ".csv"))
  #   write.csv(summary_stats$summary_df, summary_csv, row.names = FALSE)
  # })
# 
#   
#   print(summary_stats)
  
  
  require(beepr)
  beepr::beep(sound = 2)
  ##
  gc(reset = TRUE)
  ##
  .rs.restartR()  # In RStudio
  
  
  # 
  # complete_results_list = complete_results_list
  # true_DGM_Se = true_DGM_Se
  # true_DGM_Sp = true_DGM_Sp
  # min_k = vec_index_inner[1]
  # max_k = tail(vec_index_inner, 1)
  
  
  # Usage example:
  metrics <- compute_simulation_metrics(complete_results_list = complete_results_list,
                                        true_DGM_Se = true_DGM_Se, 
                                        true_DGM_Sp = true_DGM_Sp,
                                        min_k = vec_index_inner[1],
                                        max_k = tail(vec_index_inner, 1))
  
  # Create nice tables
  tables <- create_metrics_summary_table(metrics = metrics,
                                         min_k = vec_index_inner[1], 
                                         max_k = tail(vec_index_inner, 1))
  
  # Print results
  cat("\n========== SIMULATION METRICS SUMMARY ==========\n")
  cat(sprintf("Number of simulations: %d\n", metrics$n_sims))
  cat(sprintf("Thresholds analyzed: %d to %d\n\n", vec_index_inner[1], tail(vec_index_inner, 1)))
  
  cat("--- Overall Performance (averaged across thresholds) ---\n")
  print(tables$overall, row.names = FALSE)
  
  cat("\n--- Performance by Threshold ---\n")
  print(tables$by_threshold, row.names = FALSE)
  
  # Save results
  saveRDS(metrics, file.path(output_dir, paste0(file_name_string, "_metrics.RDS")))
  write.csv(tables$by_threshold, 
            file.path(output_dir, paste0(file_name_string, "_metrics_by_threshold.csv")), 
            row.names = FALSE)
  write.csv(tables$overall, 
            file.path(output_dir, paste0(file_name_string, "_metrics_overall.csv")), 
            row.names = FALSE)
  
  
  
  
  
  } ##  -----------------------  end of sim study  --------------------------------------------------------------------------------------


  # for (k in vec_index_inner) {
  #   summary_stats$summary_df[, paste0("Bias_Se_", k)] <- summary_stats$summary_df[, paste0("Se_mean_", k)] - true_DGM_Se[k]
  # }
  # 

source(file.path(getwd(), "inst", "examples", "thr_sim_study_helper_fns.R"))

# summary_stats$summary_df


message(paste("--------------------------------------------------"))

# # For crayon
# library(crayon)
# options(crayon.enabled = TRUE)
# 
# cat(red("Does this show in red?"))
#  


  ##
  ## ---- 
  ##
  {
        # ##
        # ## ---- Compute average bias for thresholds in "vec_index_inner":
        # ##
        # {
        #     result <- compute_average_bias(summary_df =   summary_stats$summary_df,
        #                                    true_DGM_Se = true_DGM_Se,
        #                                    true_DGM_Sp = true_DGM_Sp,
        #                                    min_k = vec_index_inner[1],
        #                                    max_k = tail(vec_index_inner, 1))
        #     ##
        #     cat(paste("\n--------------------------------------------------\n"))
        #     ##
        #     cat(sprintf("Mean Se bias: %s\n", formatC(result$mean_bias_Se, digits=3, format="g")))
        #     cat(sprintf("Max Se bias: %s\n", formatC(result$max_bias_Se, digits=3, format="g")))
        #     ##
        #     # cat(sprintf("Mean Se SD: %s\n", formatC(result$mean_SD_Se, digits=3, format="g")))
        #     # cat(sprintf("Max Se SD: %s\n", formatC(result$max_SD_Se, digits=3, format="g")))
        #     ##
        #     cat(sprintf("Mean Se MCSE(bias): %s\n", formatC(result$mean_MCSE_bias_Se, digits=3, format="g")))
        #     cat(sprintf("Max Se MCSE(bias): %s\n", formatC(result$max_MCSE_bias_Se, digits=3, format="g")))
        #     ##
        #     cat(paste("\n--------------------------------------------------\n"))
        #     ##
        #     cat(sprintf("Mean Sp bias: %s\n", formatC(result$mean_bias_Sp, digits=3, format="g")))
        #     cat(sprintf("Max Sp bias: %s\n", formatC(result$max_bias_Sp, digits=3, format="g")))
        #     ##
        #     # cat(sprintf("Mean Sp SD: %s\n", formatC(result$mean_SD_Sp, digits=3, format="g")))
        #     # cat(sprintf("Max Sp SD: %s\n", formatC(result$max_SD_Sp, digits=3, format="g")))
        #     ##
        #     cat(sprintf("Mean Sp MCSE(bias): %s\n", formatC(result$mean_MCSE_bias_Sp, digits=3, format="g")))
        #     cat(sprintf("Max Sp MCSE(bias): %s\n", formatC(result$max_MCSE_bias_Sp, digits=3, format="g")))
        #     ##
        #     ##
        #     ##
        #     invisible(cat(black(paste("\n model_parameterisation = ", model_parameterisation))))
        #     invisible(cat(black(paste("\n Dirichlet_random_effects_type = ", Dirichlet_random_effects_type))))
        #     cat(paste("\n--------------------------------------------------\n"))
        #     ##
        #     ##
        #     # The function also returns the updated summary_df with all the Bias values
        #     summary_df <- result$summary_df
        #     ##
        #     cat(paste("\n--------------------------------------------------\n"))
        #     ##
        #     invisible(cat(yellow(paste("\n n_sims = ", n_sims))))
        #     cat(paste("\n--------------------------------------------------\n"))
        #     ##
        #     invisible(cat(blue(paste("\n index_test_MA = ", index_test_MA))))
        #     cat(paste("\n--------------------------------------------------\n"))
        #     ##
        #     invisible(cat(red(paste("\n missing_thr = ", missing_thr))))
        #     cat(paste("\n--------------------------------------------------\n"))
        #     ##
        #     ##
        #     invisible(cat(red(paste("\n n_studies = ", n_studies))))
        #     invisible(cat(red(paste("\n N_per_study_mean = ", N_per_study_mean))))
        #     cat(paste("\n--------------------------------------------------\n"))
        #     invisible(cat(blue(paste("\n N_per_study_SD = ", N_per_study_SD))))
        #     cat(paste("\n--------------------------------------------------\n"))
        #     ##
        #     # invisible(cat(green(paste("\n bivariate_locations_nd_bs_het_SD = ", scale_ALL_bs_SDs_by*bivariate_locations_nd_bs_het_SD))))
        #     # invisible(cat(green(paste("\n bivariate_locations_d_bs_het_SD = ", scale_ALL_bs_SDs_by*bivariate_locations_d_bs_het_SD))))
        #     # cat(paste("\n--------------------------------------------------\n"))
        #     ##
        #     # invisible(cat(red(paste("\n bs_het_C_nd_prob_scale = ", scale_ALL_bs_SDs_by*bs_het_C_nd_prob_scale))))
        #     # invisible(cat(red(paste("\n bs_het_C_d_prob_scale = ", scale_ALL_bs_SDs_by*bs_het_C_d_prob_scale))))
        #     # cat(paste("\n--------------------------------------------------\n"))
        #     ##
        #     invisible(cat(green(paste("\n alpha_pi = ", alpha_pi))))
        #     cat(paste("\n--------------------------------------------------\n"))
        #     ##
        #     invisible(cat(green(paste("\n true_Mean_prev = ", true_Mean_prev))))
        #     invisible(cat(green(paste("\n true_SD_probit_prev = ", true_SD_probit_prev))))
        #     cat(paste("\n--------------------------------------------------\n"))
        #     # ##
        #     # invisible(cat(red(paste("bs_het_C_nd_prob_scale = ", bs_het_C_nd_prob_scale))))
        #     # invisible(cat(red(paste("bs_het_C_d_prob_scale = ", bs_het_C_d_prob_scale))))
        #     # cat(paste("\n--------------------------------------------------\n"))
        #     ##
        #     ##
        #     ##
        #     ##
        #     # cat(paste("\n--------------------------------------------------\n"))
        #     # cat(sprintf("Mean Se bias - ALL: %s\n", formatC(result$bias_Se_vec, digits=3, format="g")))
        #     # cat(paste("\n--------------------------------------------------\n"))
        #     # ##
        #     # cat(paste("\n--------------------------------------------------\n"))
        #     # cat(sprintf("Mean Sp bias - ALL: %s\n", formatC(result$bias_Sp_vec, digits=3, format="g")))
        #     # cat(paste("\n--------------------------------------------------\n"))
        # }
        # 
        # 
        # ##
        # ## ----
        # ##
        # ##
        {

          # cat(paste("\n--------------------------------------------------\n"))
          # message("if min_MCSE_pct = 0.5%:")
          # print(calc_min_sim_size(SD_Se_or_Sp_vec = result$SD_Se_vec,
          #                         min_MCSE_pct = 0.50,
          #                         min_k = vec_index_inner[1],
          #                         max_k = tail(vec_index_inner, 1)))
          # cat(paste("\n--------------------------------------------------\n"))
          cat(paste("\n--------------------------------------------------\n"))
          message("if min_MCSE_pct = 0.25%:")
          print(calc_min_sim_size(
            SD_Se_or_Sp_vec = metrics$bias_vectors$SD_Se_vec,
            min_MCSE_pct = 0.25, 
            min_k = vec_index_inner[1],
            max_k = tail(vec_index_inner, 1)
          ))
          
          cat(paste("\n--------------------------------------------------\n"))
          ##
          ##
          ##
          cat(paste("\n--------------------------------------------------\n"))
          message("if min_MCSE_pct = 0.125%:")
          print(calc_min_sim_size(
            SD_Se_or_Sp_vec = metrics$bias_vectors$SD_Se_vec,
            min_MCSE_pct = 0.125, 
            min_k = vec_index_inner[1],
            max_k = tail(vec_index_inner, 1)
          ))
          cat(paste("\n--------------------------------------------------\n"))
          ##
          ##
          ##
          cat(paste("\n--------------------------------------------------\n"))
          message("if min_MCSE_pct = 0.1%:")
          print(calc_min_sim_size(
            SD_Se_or_Sp_vec = metrics$bias_vectors$SD_Se_vec,
            min_MCSE_pct = 0.10, 
            min_k = vec_index_inner[1],
            max_k = tail(vec_index_inner, 1)
          ))
          cat(paste("\n--------------------------------------------------\n"))
          ##
          ##
          ##
          cat(paste("\n--------------------------------------------------\n"))
          message("if min_MCSE_pct = 0.05%:")
          print(calc_min_sim_size(
            SD_Se_or_Sp_vec = metrics$bias_vectors$SD_Se_vec,
            min_MCSE_pct = 0.05, 
            min_k = vec_index_inner[1],
            max_k = tail(vec_index_inner, 1)
          ))
          
          cat(paste("\n--------------------------------------------------\n"))
        }
    
    
    # ##
    {
      cat(paste("\n--------------------------------------------------\n"))
      message("if target_MCSE_RMSE = 0.25%:")
      print(calc_min_sim_size_RMSE(
        metrics = metrics,
        target_MCSE_RMSE = 0.25, 
        min_k = vec_index_inner[1],
        max_k = tail(vec_index_inner, 1)
      ))
      
      cat(paste("\n--------------------------------------------------\n"))
      ##
      ##
      ##
      cat(paste("\n--------------------------------------------------\n"))
      message("if target_MCSE_RMSE = 0.125%:")
      print(calc_min_sim_size_RMSE(
        metrics = metrics,
        target_MCSE_RMSE = 0.125, 
        min_k = vec_index_inner[1],
        max_k = tail(vec_index_inner, 1)
      ))
      cat(paste("\n--------------------------------------------------\n"))
      ##
      ##
      ##
      cat(paste("\n--------------------------------------------------\n"))
      message("if target_MCSE_RMSE = 0.1%:")
      print(calc_min_sim_size_RMSE(
        metrics = metrics,
        target_MCSE_RMSE = 0.10, 
        min_k = vec_index_inner[1],
        max_k = tail(vec_index_inner, 1)
      ))
      cat(paste("\n--------------------------------------------------\n"))
      ##
      ##
      ##
      cat(paste("\n--------------------------------------------------\n"))
      message("if target_MCSE_RMSE = 0.05%:")
      print(calc_min_sim_size_RMSE(
        metrics = metrics,
        target_MCSE_RMSE = 0.05, 
        min_k = vec_index_inner[1],
        max_k = tail(vec_index_inner, 1)
      ))
      
      cat(paste("\n--------------------------------------------------\n"))
    }
    
    
 
    
    
    
    
    
  }





# 
#    n_sims
#    priors$prior_dirichlet_alpha
#    priors
#    
#    
#    
#    model_summary_and_trace_obj$get_efficiency_metrics()$time_sampling
#    model_summary_and_trace_obj$get_HMC_info()$stan_HMC_info_list$out_list_sampling$mean_L_sampling
#    
#    
#  
#    
#    
#    model_summary_and_trace_obj$get_HMC_info()$HMC_diagnostic_info$divergences$pct_divs
#    sum(model_summary_and_trace_obj$get_HMC_info()$HMC_diagnostic_info$n_max_trees_per_chain) / (n_iter * n_chains)
#    
#    dplyr::filter(tibble_all, (stringr::str_detect(parameter, "beta_mu"))) %>% print(n = 1000)
#    
#    dplyr::filter(tibble_all, (stringr::str_detect(parameter, "kappa"))) %>% print(n = 1000)
#    
#    
#    dplyr::filter(tibble_all, (stringr::str_detect(parameter, "pcp_weights"))) %>% print(n = 1000)
#    
#    
#    Se_and_Sp_summaries <- dplyr::filter(tibble_all, 
#                  # stringr::str_detect(parameter, "Fp\\[") |
#                  stringr::str_detect(parameter, "Fp\\[") |
#                  stringr::str_detect(parameter, "Se\\[") ) %>% print(n = 1000)
#    
#    dplyr::filter(tibble_all, 
#                  # stringr::str_detect(parameter, "Fp\\[") |
#                  stringr::str_detect(parameter, "C_pred\\[") |
#                    stringr::str_detect(parameter, "C_mu\\[") ) %>% print(n = 1000)
#    
#    
#    
#    true_DGM_Fp
#    true_DGM_Se
#    
#    # true_DGM_Se
#    # true_Se_OVERALL[[1]]
#    # true_Se_OVERALL_weighted[[1]]
#    # 
#    # true_DGM_Sp
#    # true_Sp_OVERALL[[1]]
#    # true_Sp_OVERALL_weighted[[1]]
#    # 
#    # true_DGM_Fp
#    # true_Sp_OVERALL[[1]]
#    # true_Sp_OVERALL_weighted[[1]]
#    
#    min(Se_and_Sp_summaries$n_eff)
#    min(Se_and_Sp_summaries$n_eff) / model_summary_and_trace_obj$get_efficiency_metrics()$time_sampling
#    
# 
#    
#    
# 
#    stan_data_list$cutpoint_index
#    stan_data_list$x_with_missings
#    stan_data_list$x
#    stan_data_list$n_obs_cutpoints
#    
#    
#    sim_results$N_per_study_vec
#    
#    
#    dplyr::filter(tibble_all, (stringr::str_detect(parameter, "beta_mu"))) %>% print(n = 1000)
#    
#    dplyr::filter(tibble_all, (stringr::str_detect(parameter, "C"))) %>% print(n = 1000)
#    
#    
#    
#    summary_stats <- readRDS(file.path(getwd(), "simulation_results", "N_500n_studies_10param_Jonessoft_FALSErand_thr_FALSEsigma_C_nd0.1225sigma_C_d0.1875.RDS"))
#    
#    
#    
#  #   4.627 / 194.4  ## 0.02380144 ---- time per L
#  #   
#  #   0.24 / 30.852 ## 0.007779074
#  #   
#  #   0.02380144 / 0.007779074
#  #   
#  #   
#  #   0.435  / 30.852 ## = 0.01409957
#  #   
#  #   0.01409957 / 0.007779074
#  #   
#  #   0.02380144 / 0.01409957
#  #    
#  #   
#  #   3.775 / 254.936 ## = 0.01480764
#  #   
#  #   
#  #   6.124 / 261.336 ## 0.02343343
#  #   
#  #   
#  #   3.556 / 190.68
#  #  
#  #  
#  #   0.02380144 /  0.01864905
#  #     
#  #   3.623 / 190.68
#  #   
#  #   0.02380144 / 0.01900042
#  #   
#  #   3.692 / 228.568
#  #   
#  #   0.02380144 / 0.01615274#
#  #   
#  #   
#  #   
#  #   3.598 / 215.256
#  # 
#  #   
#  #   0.02380144 / 0.01671498
#  #   
#  #   3.384 / 255.192
#  #   
#  #   0.02380144 / 0.0132606
#  #   
#  # # true_Sp_OVERALL_weighted[[t]]
#  # # true_DGM_Sp
#  # # Sp_mean_array
#  # # 
#  # 
#  # t
#  # true_Se_OVERALL_weighted[[t]]
#  # true_DGM_Se
#  # Se_mean_array
#  # 
#  # 
#  # 
#  # dplyr::filter(tibble_all, (stringr::str_detect(parameter, "surv_prob"))) %>% print(n = 1000)
#  # 
#  # 
#  # 
#  # 
#  # stan_data_list$x_with_missings
#  # stan_data_list$x
#  # stan_data_list$cutpoint_index
#  # stan_data_list$n_obs_cutpoints
#  # ##
#  # stan_data_list$cts_thr_values
#  # 
#  # 
#  
#  
#  
#  
#  
#  
#   
 
 
 
 