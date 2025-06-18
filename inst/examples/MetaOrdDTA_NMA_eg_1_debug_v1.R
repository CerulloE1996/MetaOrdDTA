







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
       
 
}







sim_index <- 1
results_list <- c()



{
  
  
# Run simulations with different seeds
for (sim_index in 1:n_sims) {
  
  try({  
  
  require(MetaOrdDTA)
  
  seed <- sim_index
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
          if (missing_thr == TRUE) { 
            x <- x_NMA_w_missing_thr
          } else { 
            x <- x_NMA
          }
          ##
          # n_thr <- ncol(x[[1]]) - 1
          ##
          n_studies <- nrow(x[[1]][[1]])
          n_studies
          ##
          # x <- arrange_missing_values(x)
  }
  
  # true_DGM_Fp

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
{
  if (use_custom_file == TRUE) custom_file_name <- file_name
  else  custom_file_name <- NULL
}
##
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
                  indicator_index_test_in_study = indicator_index_test_in_study,
                  # n_index_tests_per_study = NULL,
                  ##
                  n_chains = n_chains,
                  ##
                  cts = cts,
                  network = TRUE,
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



init_lists_per_chain <- model_prep_obj$init_lists_per_chain
init_lists_per_chain[[1]]



##
## ----  Sample model: ----------------------------------------------------------------
##
model_samples_obj <-  model_prep_obj$sample(   
                             n_burnin = n_burnin,
                             n_iter   = n_iter,
                             adapt_delta = adapt_delta, 
                             max_treedepth = max_treedepth,
                             metric_shape = "diag_e",
                             ##
                          #   priors = priors,
                             ##
                             n_chains = n_chains
                             ##
                             # init_lists_per_chain = init_lists_per_chain
                             )

##
## ----  Summarise + output results: -------------------------------------------------
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
dplyr::filter(tibble_all, (stringr::str_detect(parameter, "Se_baseline"))) %>% print(n = 1000)
true_DGM_Se
mean(dplyr::filter(tibble_all, (stringr::str_detect(parameter, "Se_baseline\\[")))$mean*100 - true_DGM_Se)
##
dplyr::filter(tibble_all, (stringr::str_detect(parameter, "Sp_baseline"))) %>% print(n = 1000)
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
# dplyr::filter(tibble_all, (stringr::str_detect(parameter, "beta_SD"))) %>% print(n = 100)


# dplyr::filter(tibble_all, (stringr::str_detect(parameter, "beta_SD"))) %>% print(n = 1000)

try({ 
  dplyr::filter(tibble_all,
                (stringr::str_detect(parameter, "kappa"))) %>% print(n = 1000)
  dplyr::filter(tibble_all,
                (stringr::str_detect(parameter, "raw_dirichlet_phi_vec"))) %>% print(n = 1000)
})


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




   model_summary_and_trace_obj$get_efficiency_metrics()$time_sampling
   model_summary_and_trace_obj$get_HMC_info()$stan_HMC_info_list$out_list_sampling$mean_L_sampling
   model_summary_and_trace_obj$get_HMC_info()$HMC_diagnostic_info$divergences
 



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
        
        Se <- dplyr::filter(tibble_gq, (stringr::str_detect(parameter, "Se")) &
                                     (!(stringr::str_detect(parameter, "Se_pred"))))
        Sp <- dplyr::filter(tibble_gq, (stringr::str_detect(parameter, "Sp")) &
                              (!(stringr::str_detect(parameter, "Sp_pred"))))
        Fp  <- dplyr::filter(tibble_gq, (stringr::str_detect(parameter, "Fp")) &
                              (!(stringr::str_detect(parameter, "Fp_pred"))))
        ##
        Se_pred <- dplyr::filter(tibble_gq, (stringr::str_detect(parameter, "Se_pred")))
        Sp_pred <- dplyr::filter(tibble_gq, (stringr::str_detect(parameter, "Sp_pred")))
        Fp_pred <- dplyr::filter(tibble_gq, (stringr::str_detect(parameter, "Fp_pred")))
        
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
   
          ##
          ## ---- Store results for this simulation:
          ##
          results_list[[sim_index]] <- list(
            seed = seed,
            n_studies = n_studies,
            model_parameterisation = model_parameterisation,
            box_cox = box_cox,
            softplus = softplus,
            random_thresholds = random_thresholds,
            ##
            true_DGM_Se = true_DGM_Se,
            true_DGM_Sp = true_DGM_Sp,
            ##
            Se_mean_vec = Se_mean_vec,
            Sp_mean_vec = Sp_mean_vec
            ##
           ##  model_summary = model_summary_and_trace_obj  # Optional: save full summary object
          )
          
          ##
          ## ---- Update summary dataframe:
          ##
          ##
          # Create base values for the row
          row_data <- list(
            seed = seed,
            n_studies = n_studies,
            model_parameterisation = model_parameterisation,
            box_cox = box_cox,
            softplus = softplus,
            random_thresholds = random_thresholds
          )
          ##
          ## Add threshold-specific columns in a loop:
          ##
          for (k in vec_index_inner) {
            row_data[[paste0("Se_mean_", k)]] <- Se_mean_vec[k]
          }
          ##
          for (k in vec_index_inner) {
            row_data[[paste0("Sp_mean_", k)]] <- Sp_mean_vec[k]
          }
          ##
          ## Assign all data to the row:
          ##
          summary_df[sim_index, ] <- row_data
          ##
          ## ---- Save results for this simulation:
          ##
          file_name_string <- paste0( "test_", index_test_MA, 
                                      "N_", N_per_study_mean,
                                      "n_studies_", n_studies, 
                                      ##
                                      "param_", model_parameterisation, 
                                      "soft_", softplus,
                                      "rand_thr_", random_thresholds,
                                      "alpha_pi_", alpha_pi,
                                      ##
                                      "sigma_C_nd", bs_het_C_nd_prob_scale,
                                      "sigma_C_d", bs_het_C_d_prob_scale
                                    )
          
          sim_result_file <- file.path(output_dir, paste0("seed_", seed, "_", file_name_string))
          saveRDS(results_list[[sim_index]], sim_result_file)
          ##
          ## ---- Report progress:
          ##
          cat(sprintf("\n--- Simulation %d/%d completed ---\n", sim_index, n_sims))
          # cat(sprintf("Se difference from true value: %.3f\n", Se_diff))
          # cat(sprintf("Sp difference from true value: %.3f\n", Sp_diff))
      
  }
  
    
  })
  
  } ##  -----------------------  end of N_{sims} loop  --------------------------------------------------------------------------------------

  # Calculate summary statistics across ALL simulations
  summary_stats <- list(
    n_sims = n_sims,
    n_studies = n_studies,
    model_parameterisation = model_parameterisation,
    box_cox = box_cox,
    softplus = softplus,
    random_thresholds = random_thresholds,
    ##
    summary_df = summary_df
  )
  
  
  try({ 
    # Save overall summary
    summary_file <- file.path(getwd(), output_dir, paste0(file_name_string, ".RDS"))
    saveRDS(summary_stats, summary_file)
  })
  
  try({ 
    # Also save as CSV for easy viewing
    summary_csv <- file.path(getwd(), output_dir, paste0(file_name_string, ".csv"))
    write.csv(summary_df, summary_csv, row.names = FALSE)
  })
  

  
  print(summary_stats)
  
  
  
  require(beepr)
  beepr::beep(sound = 2)
  gc(reset = TRUE)
  .rs.restartR()  # In RStudio
  
  
  } ##  -----------------------  end of sim study  --------------------------------------------------------------------------------------


  for (k in vec_index_inner) {
    summary_stats$summary_df[, paste0("Bias_Se_", k)] <- summary_stats$summary_df[, paste0("Se_mean_", k)] - true_DGM_Se[k]
  }
  

source(file.path(getwd(), "inst", "examples", "thr_sim_study_helper_fns.R"))

summary_stats$summary_df


message(paste("--------------------------------------------------"))

# For crayon
library(crayon)
options(crayon.enabled = TRUE)

cat(red("Does this show in red?"))
 


  ##
  ## ---- 
  ##
  {
        ##
        ## ---- Compute average bias for thresholds in "vec_index_inner":
        ##
        {
            result <- compute_average_bias(summary_df =   summary_stats$summary_df,
                                           true_DGM_Se = true_DGM_Se,
                                           true_DGM_Sp = true_DGM_Sp,
                                           min_k = vec_index_inner[1],
                                           max_k = tail(vec_index_inner, 1))
            ##
            cat(paste("\n--------------------------------------------------\n"))
            ##
            cat(sprintf("Mean Se bias: %s\n", formatC(result$mean_bias_Se, digits=3, format="g")))
            cat(sprintf("Max Se bias: %s\n", formatC(result$max_bias_Se, digits=3, format="g")))
            ##
            # cat(sprintf("Mean Se SD: %s\n", formatC(result$mean_SD_Se, digits=3, format="g")))
            # cat(sprintf("Max Se SD: %s\n", formatC(result$max_SD_Se, digits=3, format="g")))
            ##
            cat(sprintf("Mean Se MCSE(bias): %s\n", formatC(result$mean_MCSE_bias_Se, digits=3, format="g")))
            cat(sprintf("Max Se MCSE(bias): %s\n", formatC(result$max_MCSE_bias_Se, digits=3, format="g")))
            ##
            cat(paste("\n--------------------------------------------------\n"))
            ##
            cat(sprintf("Mean Sp bias: %s\n", formatC(result$mean_bias_Sp, digits=3, format="g")))
            cat(sprintf("Max Sp bias: %s\n", formatC(result$max_bias_Sp, digits=3, format="g")))
            ##
            # cat(sprintf("Mean Sp SD: %s\n", formatC(result$mean_SD_Sp, digits=3, format="g")))
            # cat(sprintf("Max Sp SD: %s\n", formatC(result$max_SD_Sp, digits=3, format="g")))
            ##
            cat(sprintf("Mean Sp MCSE(bias): %s\n", formatC(result$mean_MCSE_bias_Sp, digits=3, format="g")))
            cat(sprintf("Max Sp MCSE(bias): %s\n", formatC(result$max_MCSE_bias_Sp, digits=3, format="g")))
            ##
            ##
            ##
            invisible(cat(black(paste("\n model_parameterisation = ", model_parameterisation))))
            invisible(cat(black(paste("\n Dirichlet_random_effects_type = ", Dirichlet_random_effects_type))))
            cat(paste("\n--------------------------------------------------\n"))
            ##
            ##
            # The function also returns the updated summary_df with all the Bias values
            summary_df <- result$summary_df
            ##
            cat(paste("\n--------------------------------------------------\n"))
            ##
            invisible(cat(yellow(paste("\n n_sims = ", n_sims))))
            cat(paste("\n--------------------------------------------------\n"))
            ##
            invisible(cat(blue(paste("\n index_test_MA = ", index_test_MA))))
            cat(paste("\n--------------------------------------------------\n"))
            ##
            invisible(cat(red(paste("\n missing_thr = ", missing_thr))))
            cat(paste("\n--------------------------------------------------\n"))
            ##
            ##
            invisible(cat(red(paste("\n n_studies = ", n_studies))))
            invisible(cat(red(paste("\n N_per_study_mean = ", N_per_study_mean))))
            cat(paste("\n--------------------------------------------------\n"))
            invisible(cat(blue(paste("\n N_per_study_SD = ", N_per_study_SD))))
            cat(paste("\n--------------------------------------------------\n"))
            ##
            invisible(cat(green(paste("\n bivariate_locations_nd_bs_het_SD = ", scale_ALL_bs_SDs_by*bivariate_locations_nd_bs_het_SD))))
            invisible(cat(green(paste("\n bivariate_locations_d_bs_het_SD = ", scale_ALL_bs_SDs_by*bivariate_locations_d_bs_het_SD))))
            cat(paste("\n--------------------------------------------------\n"))
            ##
            invisible(cat(red(paste("\n bs_het_C_nd_prob_scale = ", scale_ALL_bs_SDs_by*bs_het_C_nd_prob_scale))))
            invisible(cat(red(paste("\n bs_het_C_d_prob_scale = ", scale_ALL_bs_SDs_by*bs_het_C_d_prob_scale))))
            cat(paste("\n--------------------------------------------------\n"))
            ##
            invisible(cat(green(paste("\n true_Mean_prev = ", true_Mean_prev))))
            invisible(cat(green(paste("\n true_SD_probit_prev = ", true_SD_probit_prev))))
            cat(paste("\n--------------------------------------------------\n"))
            # ##
            # invisible(cat(red(paste("bs_het_C_nd_prob_scale = ", bs_het_C_nd_prob_scale))))
            # invisible(cat(red(paste("bs_het_C_d_prob_scale = ", bs_het_C_d_prob_scale))))
            # cat(paste("\n--------------------------------------------------\n"))
            ##
            ##
            ##
            ##
            cat(paste("\n--------------------------------------------------\n"))
            cat(sprintf("Mean Se bias - ALL: %s\n", formatC(result$bias_Se_vec, digits=3, format="g")))
            cat(paste("\n--------------------------------------------------\n"))
            ##
            cat(paste("\n--------------------------------------------------\n"))
            cat(sprintf("Mean Sp bias - ALL: %s\n", formatC(result$bias_Sp_vec, digits=3, format="g")))
            cat(paste("\n--------------------------------------------------\n"))
        }
        
        
        ##
        ## ----
        ##
        ##
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
          print(calc_min_sim_size(SD_Se_or_Sp_vec = result$SD_Se_vec, 
                            min_MCSE_pct = 0.25, 
                            min_k = vec_index_inner[1],
                            max_k = tail(vec_index_inner, 1))) 
          cat(paste("\n--------------------------------------------------\n"))
          ##
          ##
          ##
          cat(paste("\n--------------------------------------------------\n"))
          message("if min_MCSE_pct = 0.125%:")
          print(calc_min_sim_size(SD_Se_or_Sp_vec = result$SD_Se_vec, 
                            min_MCSE_pct = 0.125, 
                            min_k = vec_index_inner[1],
                            max_k = tail(vec_index_inner, 1)))
          cat(paste("\n--------------------------------------------------\n"))
          ##
          ##
          ##
          cat(paste("\n--------------------------------------------------\n"))
          message("if min_MCSE_pct = 0.1%:")
          print(calc_min_sim_size(SD_Se_or_Sp_vec = result$SD_Se_vec, 
                                  min_MCSE_pct = 0.1, 
                                  min_k = vec_index_inner[1],
                                  max_k = tail(vec_index_inner, 1)))
          cat(paste("\n--------------------------------------------------\n"))
          ##
          ##
          ##
          cat(paste("\n--------------------------------------------------\n"))
          message("if min_MCSE_pct = 0.05%:")
          print(calc_min_sim_size(SD_Se_or_Sp_vec = result$SD_Se_vec, 
                            min_MCSE_pct = 0.05, 
                            min_k = vec_index_inner[1],
                            max_k = tail(vec_index_inner, 1)))
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
 
 
 
 