




run_single_simulation <- function(start_sim_index,
                                  n_sims,
                                  model_parameterisation,
                                  softplus,
                                  box_cox,
                                  cts,
                                  random_thresholds,
                                  Dirichlet_random_effects_type,
                                  local_pkg_dir,
                                  ##
                                  sim_config,
                                  model_config
) {
  
  
  {
 
        # Extract all needed variables
        list2env(sim_config, environment())
        list2env(model_config, environment())
        
        # # Get the base directory (where your main script is)
        # if (exists("local_pkg_dir")) {
        #   base_dir <- local_pkg_dir
        # } else {
        #   base_dir <- getwd()
        # }
        
        # Create output directory with absolute path
        output_dir <- file.path(local_pkg_dir, "simulation_results")
        if (!dir.exists(output_dir)) {
          dir.create(output_dir, recursive = TRUE)
        }
        output_dir <- normalizePath(output_dir, mustWork = TRUE)
        
        cat(sprintf("Worker: Using output directory: %s\n", output_dir))
        
        source(file.path(local_pkg_dir, "inst", "examples", "thr_sim_study_helper_fns.R"))
        
      
      {
      
       
          # 
          {
            require(MetaOrdDTA)
            require(RcppParallel)
            require(BayesMVP)
            require(dplyr)
          }
      
      
      
      
          file_name_string <- paste0( "test_", index_test_MA,
                                      "N_", N_per_study_mean,
                                      "n_studies_", n_studies,
                                      ##
                                      "param_", model_parameterisation,
                                      "soft_", softplus,
                                      "rand_thr_", random_thresholds,
                                      "alpha_pi_", alpha_pi)
      
          results_list_file_name <- paste0(file_name_string, "_results_list.RDS")
          results_list_path <- file.path(output_dir, results_list_file_name)
      
          # Store the total target before we modify n_sims
          total_target_sims <- n_sims
      
          # Calculate how many more sims to run
          n_sims_to_run <- total_target_sims - start_sim_index + 1
          
          results_list <- vector("list", total_target_sims)
      
      
      
      
      }
      
       
      
        iii <- 1
  
  
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
  
  covariates <- FALSE

  
  {
    
    
      # ##
      # ## ---- Load NMA data:
      # ##
       # setwd(local_pkg_dir)
      # setwd("/")
      ## print(paste(getwd()))
      ##  
      try({ 
        rm(true_DGM_Se)
        rm(true_DGM_Sp)
        })
      source(file.path(local_pkg_dir, "R", "R_fn_sim_data_ord_MA.R"), local = TRUE)
      source(file.path(local_pkg_dir, "inst", "examples", "NMA_missing_thr_prep_data.R"), local = TRUE)
      source(file.path(local_pkg_dir, "inst", "examples", "thr_sim_study_helper_fns.R"), , local = TRUE)
      ##
      # data <- readRDS(file.path(path, "inst", "examples", "data_example_1_NMA_list.RDS"))
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
      } else if (index_test_MA == "BAI") {
        if (missing_thr == TRUE) x <- x_BAI_w_missing_thr
        else                     x <- x_BAI
        
      
        
        
      }
      ##
      n_thr <- ncol(x[[1]]) - 1
      ##
      n_studies <- nrow(x[[1]])
      n_studies
      ##
      # x <- arrange_missing_values(x)
  }
  # 
  
  
  100 * sum(x_BAI_w_missing_thr[[1]] == -1) / (n_studies * (ncol(x_BAI_w_missing_thr[[1]]) - 1))
  
  try({  
  

  ## true_DGM_Fp

##
## ----  Initialise / select model: ---------------------------------------------------------------------------
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

use_custom_file  <- TRUE

# 
if (model_parameterisation == "Bivariate") {

    custom_file_name <- "DTA_MA_strat_bivariate.stan"

} else {
  # {
  #   if (use_custom_file == TRUE) custom_file_name <- file_name
  #   else  custom_file_name <- NULL
  # }
  
  
    if ((model_parameterisation == "R&G") && (random_thresholds == FALSE)) {
      
      custom_file_name <- "DTA_MA_Gat_FIXEDthr.stan"
      
    } else if ((model_parameterisation == "R&G") && (random_thresholds == TRUE)) {
      
      custom_file_name <- "DTA_MA_Gat_RANDthr_kappa.stan"
      
    } else if ((model_parameterisation == "Xu") && (random_thresholds == FALSE)) {
      
      custom_file_name <- "DTA_MA_Xu_FIXEDthr.stan"
      
    } else if ((model_parameterisation == "Xu") && (random_thresholds == TRUE)) {
      
      custom_file_name <- "DTA_MA_Xu_RANDthr_kappa.stan"
      
    } else if ((model_parameterisation == "Jones"))  {
      
      custom_file_name <- "DTA_MA_Jones.stan"
      
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



if (covariates) {
      
      X <- list(cov_data_list$X_nd[[1]], 
                cov_data_list$X_d[[1]])
      ##
      cov_data <- list()
      ##
      cov_data$X <- X
      ##
      
      if (model_parameterisation == "R&G") { 
        cov_data$baseline_case    <- cov_data_list$baseline_case
      } else { 
        cov_data$baseline_case_nd <- cov_data_list$baseline_case_nd
        cov_data$baseline_case_d  <- cov_data_list$baseline_case_d
      }

}


# mean(X[[1]][, 2])


#
model_prep_obj <- MetaOrdDTA::MetaOrd_model$new(
                  debugging = TRUE,
                  ##
                  x = x,
                  indicator_index_test_in_study = NULL,
                  # n_index_tests_per_study = NULL,
                  ##
                  intercept_only = TRUE,
                  cov_data = cov_data,
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



      if (alpha_pi == "one")  {
        
              prior_dirichlet_alpha <- rep(1.0, n_thr + 1)
              ##
              if (model_parameterisation == "Xu") {
                priors$prior_dirichlet_alpha <- list( prior_dirichlet_alpha, prior_dirichlet_alpha)
              } else {
                priors$prior_dirichlet_alpha <- prior_dirichlet_alpha
              }
      
      } else if (alpha_pi == "flat") {

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

                        priors$prior_dirichlet_alpha <- MetaOrdDTA:::generate_ordinal_dirichlet_prior( n_cat = n_thr + 1,
                                                                                          max_alpha = 5,
                                                                                          min_alpha = 1,
                                                                                          left_skew_factor = 0.5)
              }

      }
      # 



 


priors$use_probit_link <- 1


# priors
# 
#  stan_data_list$cutpoint_index[1,,]
# #
# stan_data_list$x[1,,]


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

if (model_parameterisation == "R&G") {

    priors$prior_raw_scale_mu_SD <- 1
    priors$prior_raw_scale_SD_SD <- 1

}
# 
# exp(0 - 2*priors$prior_raw_scale_mu_SD)
# exp(0 + 2*priors$prior_raw_scale_mu_SD)

priors

priors$K_fold_CV_indicator <- rep(1, n_studies)

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
mean(abs(dplyr::filter(tibble_all, (stringr::str_detect(parameter, "Se_baseline\\[")))$mean*100 - true_DGM_Se), na.rm = TRUE)
##
dplyr::filter(tibble_all, (stringr::str_detect(parameter, "Sp"))) %>% print(n = 1000)
true_DGM_Sp
mean(abs(dplyr::filter(tibble_all, (stringr::str_detect(parameter, "Sp_baseline\\[")))$mean*100 - true_DGM_Sp), na.rm = TRUE)
##
try({
    dplyr::filter(tibble_all, (stringr::str_detect(parameter, "kappa"))) %>% print(n = 1000)
    dplyr::filter(tibble_all, (stringr::str_detect(parameter, "dirichlet_phi"))) %>% print(n = 1000)
})
##
try({
    dplyr::filter(tibble_all, (stringr::str_detect(parameter, "alpha"))) %>% print(n = 1000)
})
##
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



dplyr::filter(tibble_all, (stringr::str_detect(parameter, "Omega"))) %>% print(n = 1000)
dplyr::filter(tibble_all, (stringr::str_detect(parameter, "Sigma"))) %>% print(n = 1000)




tibble_all ## save this!
as.numeric(object.size(tibble_all)) / 1024^2 ## < 1 MB

traces_all <- model_summary_and_trace_obj$get_all_traces_list()
trace <- traces_all$traces_as_arrays
trace_gq <- trace$trace_gq ## could save this trace
as.numeric(object.size(trace_gq)) / 1024^2 ## ~ 10 MB



   model_summary_and_trace_obj$get_efficiency_metrics()$time_sampling
   model_summary_and_trace_obj$get_HMC_info()$stan_HMC_info_list$out_list_sampling$mean_L_sampling

   100 * model_summary_and_trace_obj$get_HMC_info()$HMC_diagnostic_info$divergences$n_divs /  (n_chains*(n_iter))


   
   
   mean(X[[1]][, 2])
   
   mean(X[, 2])

# 
# ##
# ## 
# ## ----  Plots: ----------------------------------------------------------------------
##
dplyr::filter(tibble_all, (stringr::str_detect(parameter, "beta_mu"))) %>% print(n = 1000)
##
new_cov_data <- list()
##
baseline_case <- c(1,              ## intercept
                  ##
                  0.0,   ## prev_GAD
                  ##
                  0, ## low_RoB_QUADAS
                  ##
                  0, ## Ref_test_SCID
                  0, ## Ref_test_Structured
                  ##
                  0, ## study_setting_2
                  0) ## study_setting_3
##
new_cov_data$baseline_case_nd <- baseline_case
new_cov_data$baseline_case_d  <- baseline_case
##
##
if (model_parameterisation == "R&G") { 
 new_cov_data$baseline_case    <- baseline_case
} else { 
 new_cov_data$baseline_case_nd <- baseline_case
 new_cov_data$baseline_case_d  <- baseline_case
}
   
plots <- model_summary_and_trace_obj$plot_sROC(new_cov_data = new_cov_data)

plots$plot_list[[1]]
plots$plot_list[[2]]

str(plots)
   
   
# # ##
# # require(Cairo)
# # options(bitmapType = "cairo")
# 
# # ##
# ## Try setting X11 as the device with specific font options
# # ## (solves bug on linux for font type when usijng ggplot)
# # Sys.setenv("R_GSCMD" = "gs")
# # options(device = "X11")
# # X11.options(family = "sans")
# # ##
# # options(bitmapType = "Xlib")
# # 
# # stan_model_file_name <- model_summary_and_trace_obj$internal_obj$outs_stan_model_name$stan_model_file_name
# # stan_mod_samples <- model_summary_and_trace_obj$internal_obj$outs_stan_sampling$stan_mod_samples
# # # ##
# # R_fn_sROC_plot_MA(stan_model_file_name = stan_model_file_name,
# #                                stan_mod_samples = stan_mod_samples)
# # ##
# # plots <- model_summary_and_trace_obj$plot_sROC()
# # plots$plot_list[[1]]
# # plots$plot_list[[2]]
# 
# # model_summary_and_trace_obj$plot_traces(param_string = "alpha")


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





 ESS_vec <- c(Se$n_eff, Sp$n_eff)
 min_ESS <- min(ESS_vec, na.rm = TRUE)
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
  # 

 # 
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
           tibble_all = tibble_all, # added this - has posterior means, medians, posterior SD, 2.5 and 97.5% quintiles, n_eff and Rhat
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
      
      
      
         # Save with overwrite check
         results_list_file <- file.path(output_dir, paste0("seed_", seed, "_", file_name_string, "_results_list.RDS"))
      
         ## Remove if exists
         try({
           if (file.exists(results_list_file)) {
             cat(sprintf("Removing existing file: %s\n", results_list_file))
             file.remove(results_list_file)
           }
         })
         # Save:
         try({
           results_list_file <- file.path(output_dir, paste0("seed_", seed, "_", file_name_string, "_results_list.RDS"))
           saveRDS(results_list[[sim_index]], results_list_file)
         })
         
         ## Verify
         try({
           if (!file.exists(results_list_file)) {
             stop(sprintf("Failed to save file: %s", results_list_file))
           } else {
             cat(sprintf("File saved successfully (size: %d bytes)\n",
                         file.info(results_list_file)$size))
           }
         })
      
         cat(sprintf("\n--- Simulation %d/%d completed ---\n", seed, n_sims))

  }
 
  
    
  })
  
  } ##  -----------------------  end of N_{sims} loop  --------------------------------------------------------------------------------------


  # 
  # Also rebuild the complete results_list from individual seed files
  complete_results_list <- vector("list", total_target_sims)
  ##
  for (seed_idx in 1:total_target_sims) {
    # Use the same filename pattern as when you save individual seeds
    seed_file <- file.path(output_dir,
                           paste0("seed_", seed_idx, "_", file_name_string, "_results_list.RDS"))
    if (file.exists(seed_file)) {
      complete_results_list[[seed_idx]] <- readRDS(seed_file)
    } else {
      cat(sprintf("Warning: seed file not found: %s\n", seed_file))
    }
  }



  # Save the complete results list
  try({
    results_list_file <- file.path(output_dir, 
                                   paste0(file_name_string, "_complete_results_list.RDS"))
    saveRDS(complete_results_list, results_list_file)
    cat(sprintf("Saved complete results list with %d non-NULL entries to: %s\n",
                sum(!sapply(complete_results_list, is.null)), results_list_file))
  })

  

   
  
  # return(complete_results_list)

  ##
  # .rs.restartR()  # In RStudio
  source(file.path(local_pkg_dir, "inst", "examples", "thr_sim_study_helper_fns.R"), local = TRUE)
  
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



  require(beepr)
  beepr::beep(sound = 2)

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
  ##
  write.csv(tables$by_threshold,
            file.path(output_dir, paste0(file_name_string, "_metrics_by_threshold.csv")),
            row.names = FALSE)
  ##
  write.csv(tables$overall,
            file.path(output_dir, paste0(file_name_string, "_metrics_overall.csv")),
            row.names = FALSE)
  # 
  # 
  ##
  gc(reset = TRUE)

  
  
  
  
  } ##  -----------------------  end of sim study  --------------------------------------------------------------------------------------



  
  
  


message(paste("--------------------------------------------------"))
 

  ##
  ## ---- 
  ##
  {

    # 
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

}








 
 