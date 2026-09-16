





# # # # # # # # ##
rm(list = ls())

# 

{
      os <- .Platform$OS.type


      if (os == "unix") {
        user_root_dir <- Sys.getenv("PWD")
      } else if (os == "windows") {
        user_root_dir <- Sys.getenv("USERPROFILE")
      }
      local_pkg_dir <- file.path(user_root_dir, "Documents/Work/PhD_work/R_packages/MetaOrdDTA")
      #


      {
        ## First remove any possible package fragments:
        ## Find user_pkg_install_dir:
        user_pkg_install_dir <- Sys.getenv("R_LIBS_USER")
        print(paste("user_pkg_install_dir = ", user_pkg_install_dir))
        ##
        ## Find pkg_install_path + pkg_temp_install_path:
        pkg_install_path <- file.path(user_pkg_install_dir, "MetaOrdDTA")
        pkg_temp_install_path <- file.path(user_pkg_install_dir, "00LOCK-MetaOrdDTA")
        ##
        ## Remove any (possible) MetaOrdDTA package fragments:
        remove.packages("MetaOrdDTA")
        unlink(pkg_install_path, recursive = TRUE, force = TRUE)
        unlink(pkg_temp_install_path, recursive = TRUE, force = TRUE)
      }

      #


      #
      #### -------- ACTUAL (LOCAL) INSTALL:
      ## Document:Q
      devtools::clean_dll(local_pkg_dir)
      roxygen2::roxygenize(local_pkg_dir)
      ##
      ## Install (outer pkg):
      ##
      devtools::clean_dll(local_pkg_dir)
      devtools::install(local_pkg_dir,
                        upgrade = "never",
                        quick = TRUE)
      ##
      ## May need to restart R:
      ##
      .rs.restartR()  # In RStudio

      # ?devtools::install

      require(MetaOrdDTA)
}
#


# use_BayesMVP_for_faster_summaries <- TRUE
# output_dir = "cv_results"
# K <- 10
# n_workers <- 10
# debugging <- FALSE

# #### ----
# {
#   model_parameterisation = "Jones"
#   ##
#   softplus <- FALSE
#   ##
#   box_cox <- TRUE
#   cts <- TRUE
#   ##
#   random_thresholds <-  FALSE
#   Dirichlet_random_effects_type <- "none"
# }
# 

#
{
  model_parameterisation = "Xu"
  box_cox <- FALSE
  cts <- FALSE
  ##
  random_thresholds <-  FALSE
  Dirichlet_random_effects_type <- "fixed"
}


{
  model_parameterisation = "Xu"
  box_cox <- FALSE
  cts <- FALSE
  ##
  random_thresholds <-  TRUE
  ##"
  Dirichlet_random_effects_type <- "kappa"
}
# # #
#
# 



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
  }
  
  source(file.path(getwd(), "inst", "examples", "thr_sim_study_helper_fns.R"))
  
  
  
  {
    
    n_chains <-  8 ##  parallel::detectCores() / 2
    ##
    if (random_thresholds)    n_iter   <- ceiling(1000/(n_chains/4))
    if (!(random_thresholds)) n_iter   <- ceiling(1000/(n_chains/4))
    n_iter
    ##
    n_burnin <- 500
    ##
    #   max_treedepth <- 9
    max_treedepth <- 10
    ##
    adapt_delta <- 0.65
    ##  adapt_delta <- 0.80
    
  }
  
  
  ##
  ## n_studies <- 10
  ## n_studies <- 50
  ########  n_studies <- 5
  ## n_studies <- 25
  ######## n_studies <- 100
  ##
  ## N_per_study_mean <- 500
  ######## N_per_study_mean <- 5000
  ######## N_per_study_mean <- 10000
  ##
   ##  N_per_study_mean <- 500  ; n_studies <- 10
  ##  N_per_study_mean <- 500  ; n_studies <- 50
  
    N_per_study_mean <- 500  ; n_studies <- 25
  
  ##
  N_per_study_SD <- N_per_study_mean
  ## N_per_study_SD <- 1 ## N_per_study_mean
  
  
  
  missing_thr <- TRUE
  # missing_thr <- FALSE
  
  
  
  
  
  
  alpha_pi <- "flat"
  # alpha_pi <- "weak_info"
  
  
  {
    
    true_Mean_prev <- 0.25 ;  true_SD_probit_prev <- 0.25
    ##
    scale_ALL_bs_SDs_by <- 0.50
    ##
    bs_het_C_nd_prob_scale <- 0.001      ; bs_het_C_d_prob_scale <- 0.001
    # bs_het_C_nd_prob_scale <- 0.03932276 ; bs_het_C_d_prob_scale <- 0.07712564 ## ---- based on "alpha" Xu model fitted to real GAD-2 data
    ##
    bivariate_locations_nd_bs_het_SD <- 0.45 ; bivariate_locations_d_bs_het_SD <- 0.60 ## ---- based on "alpha" Xu model fitted to real GAD-2 data
    ######## bivariate_locations_nd_bs_het_SD <- 0.25 ; bivariate_locations_d_bs_het_SD <- 0.40  ## ---- made up / arbitrary (ish)
    
  }
  
  
  
  index_test_chosen_index <- 1 + 1
  

  
  
  ## 
  if (n_studies < 50) { 
    n_sims <- 100
  } else { 
    n_sims <- 25
  }
  
  ##
  
  
  
  {
    print(paste("--------------------------------------------------------------"))
    print(paste("n_sims = ", n_sims))
    print(paste("n_studies = ", n_studies))
    print(paste("N_per_study_mean = ", N_per_study_mean))
    ##
    print(paste("--------------------------------------------------------------"))
    print(paste("bivariate_locations_nd_bs_het_SD = ", bivariate_locations_nd_bs_het_SD))
    print(paste("bivariate_locations_d_bs_het_SD = ", bivariate_locations_d_bs_het_SD))
    ##
    print(paste("--------------------------------------------------------------"))
    print(paste("bs_het_C_nd_prob_scale = ", bs_het_C_nd_prob_scale))
    print(paste("bs_het_C_d_prob_scale = ", bs_het_C_d_prob_scale))
    ##
    print(paste("--------------------------------------------------------------"))
    print(paste("model_parameterisation = ", model_parameterisation))
    print(paste("random_thresholds = ", random_thresholds))
    print(paste("Dirichlet_random_effects_type = ", Dirichlet_random_effects_type))
    print(paste("alpha_pi = ", alpha_pi))
    ##
    print(paste("--------------------------------------------------------------"))
  }
  
  
  
  
  
  
  {
    output_dir = "simulation_results"
    
    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # # Initialize storage for results
    # # Set up summary storage
    # summary_df <- create_summary_df(n_sims = n_sims,
    #                                 min_k = vec_index_inner[1],
    #                                 max_k = tail(vec_index_inner, 1))
    
    network <- TRUE
    ##
    softplus <- FALSE
    # softplus <- TRUE
  }
  
  
  
  # file_name <- "DTA_MA_Xu_FIXEDthr_daig_J.stan"
  ##
  # use_custom_file <- TRUE
  use_custom_file <- FALSE
  
  
  
  
}









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






# # Check if your real data X was centered:
# mean(real_data$cov_data$X[[1]][[1]][, "logit_prev_GAD"])  # Should be ~0 if centered
# 
# # Check your simulated:
# mean(X[[1]][[1]][, "logit_prev_GAD"])  # Probably not 0
# 




sim_index <- 1
results_list <- c()




  intercept_only <- TRUE ; MR_model <- "intercept_only"
##
# intercept_only <- FALSE ;  MR_model <- "MR_model_1" ## intercept + logit_prev_GAD
# intercept_only <- FALSE ;  MR_model <- "MR_model_2" ## intercept + logit_prev_GAD + ref_test
# intercept_only <- FALSE ;  MR_model <- "MR_model_3" ## intercept + logit_prev_GAD + study_setting
##
# {
#   
#   
# # Run simulations with different seeds
# for (sim_index in 1:n_sims) {
#   
#   try({  
#   
  require(MetaOrdDTA)
  
  seed <- sim_index
  cat(sprintf("\n\n--- Running simulation %d/%d with seed = %d ---\n\n", sim_index, n_sims, seed))
  
  # Set seed for reproducibility
  seed <- 1
  ##
  set.seed(seed)
  

    {
          # ##
          # ## ---- Load NMA data:
          # ##
          setwd(local_pkg_dir)
          ##
          intercept_only_for_sim_DGM <- FALSE
          MR_model_for_sim_DGM <- "MR_model_5"
          ##
          source(file.path(getwd(), "R", "R_fn_sim_data_ord_NMA.R"))
          source(file.path(getwd(), "inst", "examples", "NMA_missing_thr_prep_data_v2.R"))
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
          x_old <- x
          x <- list()
          ##
          for (t in 1:n_index_tests) { 
            x[[t]] <- list()
          }
          for (t in 1:n_index_tests) { 
            for (c in 1:2) {
              x[[t]][[c]] <- x_old[[c]][[t]]
            }
          }
          indicator_index_test_in_study <- create_indicator_index_test_in_study(x)
  }
  
 x
 X
 
 # X <- MetaOrdDTA:::R_fn_expand_covariates_to_all_studies( 
 #   X = X,
 #   indicator_index_test_in_study = indicator_index_test_in_study)
 # ##
 # X <- apply_test_indicators_to_X(X = X,
 #                                     indicator_matrix = indicator_index_test_in_study)
 # ##
 # X
 
 
 # Debug: What's the actual structure?
 which(rowSums(X[[1]][[4]] != 0) > 0)  # Which rows have data?
 which(rowSums(original_cov_data$X[[1]][[4]] != 0) > 0)  # Which rows have data?
 which(indicator_index_test_in_study[, 4] == TRUE)  # Which studies should have test 4?
 
 ifelse(sim_results$indicator_index_test_in_study, 1, 0)


##
## ----  Initialise / select model: --------------------------------------------------------------------------- NMA:
##
init_lists_per_chain <- NULL
##
##
##
##
{
  if (use_custom_file == TRUE) custom_file_name <- file_name
  else  custom_file_name <- NULL
}

if (intercept_only) { 
  
        X_nd <- rep(list(array(1, dim = c(n_studies, 1))), n_index_tests)
        X_d  <- X_nd 
        
        X <- list( X_nd, X_d)
        ##
        X <- apply_test_indicators_to_X(X = X,
                                        indicator_matrix = indicator_index_test_in_study)
        ##
        cov_data <- list()
        ##
        cov_data$X <- X
        ##
        
        if (model_parameterisation == "R&G") { 
          baseline_case <- 1
          cov_data$baseline_case    <- rep(list(baseline_case), n_index_tests)
        } else { 
          baseline_case_nd <- 1
          baseline_case_d  <- 1
          cov_data$baseline_case_nd <- rep(list(baseline_case_nd), n_index_tests)
          cov_data$baseline_case_d  <- rep(list(baseline_case_d), n_index_tests)
        }
  
  
} else { 
  
        # baseline_case <- c(1,              ## intercept
        #                    ##
        #                    0.0,   ## prev_GAD
        #                    # mean(X[, 3]),   ## prev_AAD
        #                    ##
        #                    0, ## low_RoB_QUADAS
        #                    # 0, ## low_RoB_liberal
        #                    ##
        #                    0, ## Ref_test_SCID
        #                    0, ## Ref_test_Structured
        #                    ##
        #                    0, ## study_setting_2
        #                    1) ## study_setting_3 (most common)
        ##
        if (MR_model == "MR_model_1") { 
            X <- subset_X_covariates( X = sim_results$X_mat, 
                                      covariates_to_keep = c("intercept", 
                                                             "logit_prev_GAD"))
            baseline_case <- c(1, 0.0)
        } else if (MR_model == "MR_model_2") {
            X <- subset_X_covariates( X = sim_results$X_mat, 
                                      covariates_to_keep = c("intercept", 
                                                             "Ref_test_clean_SCID", "Ref_test_clean_Structured"))
            baseline_case <- c(1, 0.0, 0, 0)
        } else if (MR_model == "MR_model_3") {
            X <- subset_X_covariates( X = sim_results$X_mat, 
                                      covariates_to_keep = c("intercept", 
                                                             "study_setting_2", "study_setting_3"))
            baseline_case <- c(1, 0.0, 0, 1)
        } else if (MR_model == "MR_model_4") {
            # X <- subset_X_covariates( X = sim_results$X_mat, 
            #                           covariates_to_keep = c("intercept",
            #                                                  "logit_prev_GAD",
            #                                                  "Ref_test_clean_SCID", "Ref_test_clean_Structured",
            #                                                  "study_setting_2", "study_setting_3"))
            # baseline_case <- c(1, 0.0, 0, 0, 0, 1)
        } else if (MR_model == "MR_model_5") {
            X <- subset_X_covariates( X = sim_results$X_mat, 
                                      covariates_to_keep = c("intercept", 
                                                             "logit_prev_GAD",
                                                             "low_RoB_QUADAS_clean",
                                                             "Ref_test_clean_SCID", "Ref_test_clean_Structured",
                                                             "study_setting_2", "study_setting_3"))
            baseline_case <- c(1, 0.0, 0, 0, 0, 0, 1)
        }
        ##
        X <- apply_test_indicators_to_X(X = X,
                                        indicator_matrix = indicator_index_test_in_study)
        ##
        n_covariates <- ncol(X[[1]][[1]])
        print(paste("n_covariates = ", n_covariates))
        ##
        X <- MetaOrdDTA:::R_fn_expand_covariates_to_all_studies( 
                            X,
                            indicator_index_test_in_study = indicator_index_test_in_study)
        ##
        cov_data <- list()
        ##
        cov_data$X <- X
        ##
        cov_data$baseline_case_nd <- rep(list(baseline_case), n_index_tests)
        cov_data$baseline_case_d  <- rep(list(baseline_case), n_index_tests)
        
}

X


# compute_study_overlap_prop(indicator_index_test_in_study)
# ##
# indicator_index_test_in_study
# ##
# ##
# X[[1]][[4]]
# X[[2]][[4]]
# sum( (X[[1]][[4]] == X[[2]][[4]]) == FALSE)
# ##
# x[[4]][[1]][1:n_studies, 1:10] ## 16, 19, 24 , 35
# x[[4]][[2]][1:n_studies, 1:10] ## 16, 19, 24 , 35
# ##
# real_data$x[[4]][[1]][, 1:10]
# real_data$x[[4]][[2]][, 1:10]
##
# compare_indicator_matrices( sim_matrix = indicator_index_test_in_study,
#                             real_matrix = real_data$indicator_index_test_in_study,
#                             test_names = test_names, 
#                             detailed = TRUE)
# compute_study_overlap_prop(indicator_index_test_in_study)
# ##
# Total studies: 40 
# Studies with 1 test: 23 
# Studies with 2+ tests: 17 
# Study overlap proportion: 0.425 
# Total studies: 76 
# Studies with 1 test: 46 
# Studies with 2+ tests: 30 
# Study overlap proportion: 0.395 
# ===== INDICATOR MATRIX COMPARISON =====
#   
#   1. TEST FREQUENCIES:
#   Test Simulated  Real Difference
# GAD-2     50.0% 38.2%      11.8%
# GAD-7     50.0% 51.3%       1.3%
# HADS     45.0% 47.4%       2.4%
# BAI     10.0%  9.2%       0.8%
# 
# 2. OVERLAP PROPORTION:
#   Simulated: 42.5%
# Real: 39.5%
# Difference: 3.0%
# 
# 3. TESTS PER STUDY:
#   Simulated distribution:
#   
#   1    2    3 
# 57.5 30.0 12.5 
# Real distribution:
#   
#   1    2    3 
# 60.5 32.9  6.6 
# 
# 4. TOP TEST COMBINATIONS:
#   
#   Simulated (top 5):
#   HADS: 32.5%
# GAD-2+GAD-7: 30.0%
# GAD-2+GAD-7+HADS: 12.5%
# BAI: 10.0%
# GAD-7: 7.5%
# 
# Real (top 5):
#   HADS: 38.2%
# GAD-2+GAD-7: 27.6%
# GAD-7: 13.2%
# BAI: 5.3%
# GAD-2+GAD-7+HADS: 5.3%
# 
# 5. TEST CORRELATIONS:
#   Simulated:
#   Test1 Test2 Test3 Test4
# Test1  1.00  0.70  -0.4 -0.33
# Test2  0.70  1.00  -0.4 -0.33
# Test3 -0.40 -0.40   1.0 -0.30
# Test4 -0.33 -0.33  -0.3  1.00
# 
# Real:
#   [,1]  [,2]  [,3]  [,4]
# [1,]  1.00  0.60 -0.53 -0.16
# [2,]  0.60  1.00 -0.66 -0.14
# [3,] -0.53 -0.66  1.00 -0.21
# [4,] -0.16 -0.14 -0.21  1.00
# 
# ===== SIMILARITY SCORES =====
#   Frequency similarity: 96%
# Overlap similarity: 97%
# Correlation similarity: 84%
# Combination similarity: 77%
# 
# OVERALL SIMILARITY: 89%
# > compute_study_overlap_prop(indicator_index_test_in_study)
# Total studies: 40 
# Studies with 1 test: 23 
# Studies with 2+ tests: 17 
# Study overlap proportion: 0.425 
# [1] 0.425

# # Check all your objects
# object.size(stan_data_list)
# object.size(priors)
# object.size(init_lists_per_chain)
# ##
# format(object.size(stan_data_list), units = "MB")
# format(object.size(priors), units = "MB")
# format(object.size(init_lists_per_chain), units = "MB")
# ##
# format(object.size(sim_results), units = "MB")
# format(object.size(tibble_all), units = "MB")
# format(object.size(tibble_gq), units = "MB")
##
cmdstanr_args <- list()
##
n_chains <- 16
##
model_prep_obj <- MetaOrdDTA::MetaOrd_model$new(  
                  debugging = TRUE,
                  ##
                  x = x,
                  indicator_index_test_in_study = indicator_index_test_in_study,
                  ##
                  intercept_only = intercept_only,
                  cov_data = cov_data,
                  ##
                  n_chains = n_chains,
                  ##
                  cts = cts,
                  ##
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
                  # init_lists_per_chain = init_lists_per_chain,
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
  init_lists_per_chain <- model_prep_obj$init_lists_per_chain
  priors <- model_prep_obj$priors
  stan_data_list <- model_prep_obj$internal_obj$outs_data$stan_data_list
}
##
n_total_obs = 0;
for (t in 1:n_index_tests) {
  for (s in 1:n_studies) {
    if (indicator_index_test_in_study[s, t] == 1) {
      n_total_obs  = n_total_obs +  1;
    }
  }
}
##
n_nuisance_z <- 2*n_total_obs + 2*n_studies; n_nuisance_z
# ##
n_thr <- stan_data_list$n_thr
##
model_parameterisation
random_thresholds
##
priors$compound_symmetry <- 0
priors$hetero_sigma      <- 0


##
## ----  Sample model: ----------------------------------------------------------------
##
if (random_thresholds == FALSE) {
    n_iter <- 500
} else { 
    n_iter <- 1000
}
##
n_iter
# ##
  # model_samples_obj <-  model_prep_obj$sample(
  #                            n_burnin = n_burnin,
  #                            n_iter   = n_iter,
  #                            adapt_delta = adapt_delta,
  #                            max_treedepth = max_treedepth,
  #                            metric_shape = "diag_e",
  #                            ##
  #                            priors = priors,
  #                            ##
  #                            n_chains = n_chains,
  #                            ##
  #                            init_lists_per_chain = init_lists_per_chain
  #                            )

##
## ---- Using BayesMVP:
##
# Compile and fit model
Stan_model_file_path <- model_prep_obj$internal_obj$outs_stan_compile$stan_model_file_path
Stan_functions_path <- model_prep_obj$internal_obj$outs_stan_compile$stan_functions_directory
##
mod <- cmdstan_model(Stan_model_file_path,
                     include_paths = c(Stan_functions_path))
                     
##
# ## set Stan model file path for your Stan model (replace with your path)
# Stan_model_file_path <- system.file("stan_models/basic_logistic.stan", package = "BayesMVPtest")
## make y a matrix first (matrix w/ 1 col)
### input N to use 
N <- 1000
latent_results <- rlogis(n = N, location =  0.0) # for logistic regression
y <- ifelse(latent_results > 0, 1, 0)
y <- matrix(data = c(y), ncol = 1)
##
# Stan_init_list <- init_lists_per_chain[[1]]
# init_lists_per_chain <- rep(list(Stan_init_list), n_chains) 
##
{
  stan_data <- stan_data_list
  stan_data$y <- y ## dummy data
  ##
  model_prep_obj$internal_obj$outs_stan_compile$stan_model_obj$variables()$parameters
  model_prep_obj$internal_obj$outs_stan_compile$stan_model_obj
  ##
  
  stan_data$N <- N
  stan_data$y <- y
  ##
  Stan_data_list = stan_data
  ##
  Stan_data_list$y <- NULL
  Stan_data_list$N <- NULL
  ##
  Stan_data_list$softplus <- 0
  ##
  # outs <- get_BayesMVP_Stan_paths()
  # pkg_dir <- outs$pkg_dir
  # data_dir <- outs$data_dir
  # stan_dir <- outs$stan_dir
  ##
  # str(Stan_data_list)
  Stan_data_list$x_with_missings <- NULL
  Stan_data_list$n_covariates_max <- NULL
  ##
  Stan_data_list_old <- Stan_data_list
  Stan_data_list <- c(Stan_data_list_old, priors)
  # str(Stan_data_list)
  Stan_data_list$cts_thr_values <- NULL
  Stan_data_list$cts_thr_values_nd <- NULL
  Stan_data_list$cts_thr_values_d <- NULL
  Stan_data_list$compound_symmetry
}
##
## set bs environment variable (otherwise it'll try downloading it even if already installed...)
bs_path <- bridgestan_path()
Sys.setenv(BRIDGESTAN = bs_path)
##
## For fixed-C:
if (random_thresholds == FALSE) {
    Stan_data_list$prior_kappa_mean <- NULL
    Stan_data_list$prior_kappa_SD <- NULL
    Stan_data_list$kappa_lb <- NULL
    Stan_data_list$n_total_C_if_random <- NULL
    Stan_data_list$prior_dirichlet_phi <- NULL
    Stan_data_list$n_total_summary_cat <- NULL
    Stan_data_list$n_total_cat <- NULL
    Stan_data_list$n_total_cutpoints <- NULL
    Stan_data_list$prior_dirichlet_alpha <- array(1.0, dim = c(n_index_tests, Stan_data_list$n_cat_max))
}  

##
# str(Stan_data_list)
##
stanc_args <- list(paste0("--include-paths=", Stan_functions_path))
##
n_nuisance <-  n_nuisance_z 
outs <- get_n_params_main_from_n_nuisance_using_Stan_model( n_nuisance = n_nuisance,
                                                    Stan_data_list = Stan_data_list,
                                                    Stan_model_file_path = Stan_model_file_path,
                                                    stanc_args = stanc_args)
n_params <- outs$n_params ; n_params
n_params_main <- outs$n_params_main ; n_params_main
##
n_chains_burnin <- 8
init_lists_per_chain <- resize_init_list( init_lists_per_chain = init_lists_per_chain, 
                                          n_chains_new = n_chains_burnin)
##
sample_nuisance <- TRUE
##
Stan_data_list$n_studies

for (kk in 1:length(init_lists_per_chain)) {
    init_lists_per_chain[[kk]]$beta_eta_z <- array(0.001, dim = c(2, n_studies))
}

###  -----------  Compile + initialise the model using "MVP_model$new(...)" 
model_prep_obj <- BayesMVPtest:::MVP_model$new(   Model_type =  "Stan",
                                            y = y,
                                            N = N,
                                            ##  model_args_list = model_args_list, # this arg is only needed for BUILT-IN (not Stan) models
                                            Stan_data_list = Stan_data_list,
                                            Stan_model_file_path = Stan_model_file_path,
                                            init_lists_per_chain = init_lists_per_chain,
                                            sample_nuisance = sample_nuisance,
                                            n_chains_burnin = n_chains_burnin,
                                            ##
                                            n_params_main = n_params_main,
                                            n_nuisance = n_nuisance,
                                            ##
                                            stanc_args = stanc_args)



## ----------- Set basic sampler settings
{
  seed <- 1
  # n_chains_sampling <- max(64, parallel::detectCores() / 2)
  n_chains_sampling <- 64
  n_superchains <- min(8, parallel::detectCores() / 2)  ## round(n_chains_sampling / n_chains_burnin) # Each superchain is a "group" or "nest" of chains. If using ~8 chains or less, set this to 1. 
  n_iter <- 1000                                 
  n_burnin <- 500
  ## n_nuisance_to_track <- n_nuisance # set to some small number (< 10) if don't care about making inference on nuisance params (which is most of the time!)
}


## Since this model (i.e., basic logistic reg.) does NOT have any nuisance parameters, we will set:
## partitioned_HMC = FALSE - this is because partitioned HMC is only available if sample_nuisance == TRUE
## (since we sample the nuisance params and main params seperately)
##
## ---- DEFAULT - SNAPER-HMC for burnin + standard-HMC for sampling - NOT WORKING PROPERLY (L gets too low/ adapts to L = 1 !!!!!!!!)
##
partitioned_HMC <- FALSE ;    diffusion_HMC <- FALSE ## this uses SNAPER-HMC for burnin + standard HMC forW sampling
# metric_shape_main <- "diag" ; metric_type_main <- "Hessian"
metric_shape_main <- "diag" ; metric_type_main <- "Empirical"
##
## ---- WITH diffusion-pathspace-HMC:
# ##
partitioned_HMC <- TRUE ;    diffusion_HMC <- TRUE ## this uses SNAPER-HMC for burnin + DIFFUSION-PATHSPACE HMC for sampling (on models w/ nuisance params)
# metric_shape_main <- "diag" ; metric_type_main <- "Hessian"
metric_shape_main <- "diag" ; metric_type_main <- "Empirical"
##
# ## ---- PARTITIONED HMC but WITHOUT diffusion-pathspace-HMC:
# ##
partitioned_HMC <- TRUE ;    diffusion_HMC <- FALSE ## this uses SNAPER-HMC for burnin + standard HMC forW sampling
 # metric_shape_main <- "diag" ; metric_type_main <- "Hessian"
metric_shape_main <- "diag" ; metric_type_main <- "Empirical"
# # # # # 

##
model_samples <-  model_prep_obj$sample(  partitioned_HMC = partitioned_HMC,
                                     diffusion_HMC = diffusion_HMC,
                                     seed = seed,
                                     n_burnin = n_burnin,
                                     n_iter = n_iter,
                                     n_chains_sampling = n_chains_sampling,
                                     n_superchains = n_superchains,
                                     # ## Some other arguments:
                                     # y = y,
                                     # N = N,
                                     # n_params_main = n_params_main,
                                     # n_nuisance = n_nuisance,
                                     init_lists_per_chain = init_lists_per_chain,
                                     n_chains_burnin = n_chains_burnin,
                                     # model_args_list = model_args_list,
                                     ## Some other SAMPLER / MCMC arguments:
                                     # sample_nuisance = sample_nuisance,
                                     ##
                                     adapt_delta = 0.80,
                                     learning_rate = 0.075,
                                     ##
                                     metric_shape_main = metric_shape_main,
                                     metric_type_main = metric_type_main,
                                     ##
                                     tau_mult = 1.6,
                                     ##
                                     clip_iter = 50,
                                     clip_iter_tau = 150,
                                     ##
                                     interval_width_main = 10,
                                     interval_width_nuisance = 10,
                                     ratio_M_us = 0.50,
                                     ratio_M_main = 0.50,
                                     ##
                                     beta1_adam = 0.0,
                                     beta2_adam = 0.95,
                                     eps_adam = 1e-8,
                                     ##
                                     parallel_method = "RcppParallel")
 

require(bridgestan)
model_fit <- model_samples$summary(save_log_lik_trace = FALSE, 
                                   ##
                                   compute_nested_rhat = FALSE,
                                   ##
                                   compute_transformed_parameters = TRUE,
                                   compute_generated_quantities = TRUE)


model_fit$get_HMC_info()
model_fit$get_summary_main() %>% print(n=10)
model_fit$get_summary_transformed()
model_fit$get_summary_generated_quantities()






MetaOrdDTA:::extract_params_from_tibble_batch(debugging = FALSE,
                                              tibble = tibble_gq, 
                                              param_strings_vec = c("Se"),
                                              condition = "containing") %>% print(n=100)





# 
# A tibble: 8 × 6
# parameter      `2.5%`  `50%` `97.5%` n_eff  Rhat
# <chr>           <dbl>  <dbl>   <dbl> <dbl> <dbl>
# 1 beta_mu[1,1,1] -1.92  -1.14   -0.428   327  1.05
# 2 beta_mu[2,1,1] -2.04  -1.57   -1.14    203  1.06
# 3 beta_mu[3,1,1] -1.61  -1.15   -0.735    97  1.12
# 4 beta_mu[4,1,1] -2.46  -1.92   -1.10    834  1.01
# 5 beta_mu[1,2,1] -0.353  0.272   0.890   793  1.01
# 6 beta_mu[2,2,1] -0.891 -0.370   0.121   594  1.03
# 7 beta_mu[3,2,1] -0.795 -0.285   0.195   205  1.07
# 8 beta_mu[4,2,1] -1.45  -0.712   0.161   237  1.05
#  
# 
# 
# model_samples$n_nuisance <- n_nuisance
#  model_samples$n_params_main <- n_params - n_nuisance


 
A tibble: 2,589 × 9
parameter           mean     sd  `2.5%`   `50%`  `97.5%` n_eff  Rhat n_Rhat
<chr>              <dbl>  <dbl>   <dbl>   <dbl>    <dbl> <dbl> <dbl> <lgl> 
1 beta_mu.1.1.1    -1.13   0.377  -1.89   -1.13   -0.406    1057  1.05 NA    
2 beta_mu.2.1.1    -1.57   0.241  -2.06   -1.56   -1.12      531  1.08 NA    
3 beta_mu.3.1.1    -1.16   0.228  -1.61   -1.16   -0.728     372  1.13 NA    
4 beta_mu.4.1.1    -1.91   0.326  -2.47   -1.94   -1.19      456  1.09 NA    
5 beta_mu.1.2.1     0.283  0.320  -0.352   0.282   0.901    1734  1.03 NA    
6 beta_mu.2.2.1    -0.381  0.238  -0.852  -0.371   0.0672    324  1.14 NA    
7 beta_mu.3.2.1    -0.281  0.256  -0.785  -0.275   0.194     352  1.13 NA    
8 beta_mu.4.2.1    -0.728  0.404  -1.47   -0.742   0.100     671  1.07 NA    
## min_ESS (snaper) = 324 (for coeffs/beta). 
## Time w/ 64 chains + 1000 iter: 274.709 seconds / 4.578483 mins
## ESS/sec (sampling) =  1.17943
##
## Overall Ess/sec =  0.5637636

A tibble: 8 × 6
parameter      `2.5%`  `50%` `97.5%` n_eff  Rhat
<chr>           <dbl>  <dbl>   <dbl> <dbl> <dbl>
1 beta_mu[1,1,1] -1.90  -1.13   -0.418  2348  1.03
2 beta_mu[2,1,1] -2.09  -1.59   -1.13   1438  1.03
3 beta_mu[3,1,1] -1.68  -1.17   -0.728   789  1.06
4 beta_mu[4,1,1] -2.51  -1.93   -1.11   4237  1.01
5 beta_mu[1,2,1] -0.343  0.282   0.910  5289  1.01
6 beta_mu[2,2,1] -0.888 -0.362   0.139  3887  1.02
7 beta_mu[3,2,1] -0.823 -0.267   0.247   776  1.06
8 beta_mu[4,2,1] -1.48  -0.756   0.104  1910  1.03
## min_ESS (stan) = 776 (for coeffs/beta)
## Time w/ 64 chains + 1000 iter: 1209.98 secs /  20.16633 mins
## ESS/sec (sampling) = 0.6413329
##
## time_sampling = 1209.98
## time_total_wo_summaries = 2693.472 / 44 mins
## Overall Ess/sec = 0.288104

# A tibble: 8 × 6
# parameter      `2.5%`  `50%` `97.5%` n_eff  Rhat
# <chr>           <dbl>  <dbl>   <dbl> <dbl> <dbl>
# 1 beta_mu[1,1,1] -1.92  -1.14   -0.428   327  1.05
# 2 beta_mu[2,1,1] -2.04  -1.57   -1.14    203  1.06
# 3 beta_mu[3,1,1] -1.61  -1.15   -0.735    97  1.12
# 4 beta_mu[4,1,1] -2.46  -1.92   -1.10    834  1.01
# 5 beta_mu[1,2,1] -0.353  0.272   0.890   793  1.01
# 6 beta_mu[2,2,1] -0.891 -0.370   0.121   594  1.03
# 7 beta_mu[3,2,1] -0.795 -0.285   0.195   205  1.07
# 8 beta_mu[4,2,1] -1.45  -0.712   0.161   237  1.05
#  


model_summary_and_trace_obj$get_HMC_info()
model_summary_and_trace_obj$get_efficiency_metrics()


##
## ----  Summarise + output results: -------------------------------------------------
##
RcppParallel::setThreadOptions(numThreads = n_chains);
##
model_summary_and_trace_obj <- model_samples_obj$summary(
                                          compute_main_params = TRUE,
                                          compute_transformed_parameters = TRUE, 
                                          compute_generated_quantities = TRUE,
                                          ##
                                          save_log_lik_trace = TRUE,
                                          ##
                                          use_BayesMVP_for_faster_summaries = TRUE,
                                          ##
                                          compute_nested_rhat = FALSE)
##
test_names = c("GAD-2", "GAD-7", "HADS", "BAI")
##
tibble_main <- model_summary_and_trace_obj$get_summary_main() %>% print(n = 200)
tibble_tp <- model_summary_and_trace_obj$get_summary_transformed() %>% print(n = 100)
tibble_gq   <- model_summary_and_trace_obj$get_summary_generated_quantities() %>% print(n = 1000)
##
## tibble_all <- rbind(tibble_main, tibble_gq)
tibble_all <- rbind(tibble_main, tibble_tp, tibble_gq)
##
## ---- Save full model output:
##
{
    if (length(test_names) == 4) { 
      dummy_data <- 1
      output_dir <- "application_results_dummy_data"
    } else { 
      dummy_data <- 0
      output_dir <- "application_results_real_data"
    }
    ##
    if (stan_data_list$n_covariates_max == 1) { 
      intercept_only <- 1
    } else { 
      intercept_only <- 0
    }
    ##
    file_name_string <- paste0( "dummy_data_", dummy_data,
                                "intercept_only_", intercept_only,
                                "N_covs_", stan_data_list$n_covariates_max,
                                "param_", model_parameterisation,
                                "rand_thr_", random_thresholds,
                                "CS_", priors$compound_symmetry,
                                "het_sigma_", priors$hetero_sigma)
    ##
    # results_list_file_name <- paste0(file_name_string, "_applied_results.RDS")
    # results_list_path <- file.path(output_dir, results_list_file_name)
    ##
    # Save with overwrite check
    results_list_file <- file.path(output_dir, paste0("seed_", seed, "_", file_name_string, "_applied_results.RDS"))
    ##
    ## Remove if exists
    try({
      if (file.exists(results_list_file)) {
        cat(sprintf("Removing existing file: %s\n", results_list_file))
        file.remove(results_list_file)
      }
    })
    # Save:
    try({
      results_list_file <- file.path(output_dir, paste0("seed_", seed, "_", file_name_string, "_applied_results.RDS"))
      saveRDS(list("model_prep_obj" = model_prep_obj,
                   "model_summary_and_trace_obj" = model_summary_and_trace_obj
                   # "model_samples_obj" = model_samples_obj
                   ), 
              results_list_file)
    })
    # ##
    # ## Verify
    # try({
    #   if (!file.exists(results_list_file)) {
    #     stop(sprintf("Failed to save file: %s", results_list_file))
    #   } else {
    #     cat(sprintf("File saved successfully (size: %d bytes)\n",
    #                 file.info(results_list_file)$size))
    #   }
    # })
}

##
## ----------- Plot sROC curve(s) (NMA - if WITHOUT covariates):
##
plots <- model_summary_and_trace_obj$plot_sROC(test_names = test_names)
##
plots$plot_list[[1]]
# plots$plot_list[[2]]
plots$plot_list[[3]]
# plots$plot_list[[4]]
plots$plot_list[[5]]


 



# 1. Run k-fold for each model
##
# 2. Create table:
#    Model               | K-fold ELPD | SE    | Δ from intercept-only
#    Intercept-only      | -15,554     | 3,668 | --
#    Model 1 (prev)      | ?           | ?     | ?
#    Model 2 (prev+ref)  | ?           | ?     | ?
#    Model 3 (prev+sett) | ?           | ?     | ?
##
# 3. RoB subgroup analysis separately (no covariates)



## bookmark:
# We used a hierarchical approach to meta-regression, starting with disease 
# prevalence and systematically adding reference test type and study setting. 
# Model selection was based on k-fold cross-validation. Risk of bias was examined 
# through pre-specified subgroup analysis.

##
## 
## -------------------------------------------------------------------------
# ##
 

Dirichlet_random_effects_type
priors$compound_symmetry



# ##
# dplyr::filter(tibble_all, (stringr::str_detect(parameter, "beta_Sigma"))) %>% print(n = 100)
# dplyr::filter(tibble_all, (stringr::str_detect(parameter, "beta_Omega"))) %>% print(n = 100)
# ##
# dplyr::filter(tibble_all, (stringr::str_detect(parameter, "rho12"))) %>% print(n = 100)
# dplyr::filter(tibble_all, (stringr::str_detect(parameter, "rho"))) %>% print(n = 100)
# ##
##
## ---- For 1st table:
##
{
    # dplyr::filter(tibble_all, (stringr::str_detect(parameter, "AUC"))) %>% print(n = 100)
    ##
      tibble_AUC <- model_summary_and_trace_obj$extract_params(params = c("AUC"))  %>% 
      select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
      print(n = 100)
    ##
      tibble_AUC_pred <- model_summary_and_trace_obj$extract_params(params = c("AUC_pred"))  %>% 
      select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
      print(n = 100)
    ##
      tibble_AUC_diff <- model_summary_and_trace_obj$extract_params(params = c("AUC_diff"))  %>% 
      select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
      print(n = 100)
      ##
      ##
    tibble_Se_baseline  <- model_summary_and_trace_obj$extract_params(params = c("Se_baseline"))  %>% 
      select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
      print(n = 100)
    ##
    tibble_Sp_baseline  <- model_summary_and_trace_obj$extract_params(params = c("Sp_baseline"))  %>% 
      select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
      print(n = 100)
}
##
## ---- For 2nd table:
##
{
      tibble_beta_mu <- model_summary_and_trace_obj$extract_params(params = c("beta_mu"))  %>% 
      select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
      print(n = 100)
    ##
      tibble_total_SD_inc_C <- model_summary_and_trace_obj$extract_params(params = c("total_SD_inc_C")) %>% 
      select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
      print(n = 100)
    ##
      tibble_beta_sigma <- model_summary_and_trace_obj$extract_params(params = c("beta_sigma")) %>% 
      select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
      print(n = 100)
    ##
      tibble_beta_tau <- model_summary_and_trace_obj$extract_params(params = c("beta_tau")) %>% 
      select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
      print(n = 100)
    ##
      tibble_beta_corr <- model_summary_and_trace_obj$extract_params(params = c("beta_corr")) %>% 
      select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
      print(n = 100)
    ##
      tibble_rho12 <- model_summary_and_trace_obj$extract_params(params = c("rho12")) %>% 
      select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
      print(n = 100)
}



# X_summary <- R_fn_summarize_covariates_by_test(X_list = X[[1]])
# X_summary
##
beta_mu <- model_summary_and_trace_obj$extract_params(params = c("beta_mu")) %>% print(n = 100)
##
R_fn_map_beta_coefficients(beta_tibble = beta_mu, 
                           covariate_names = dimnames(X[[1]][[1]])[[2]], 
                           test_names = test_names, 
                           alpha = 0.05)




 


Model_A <- readRDS(file.path(
           "application_results_dummy_data",
           "seed_1_dummy_data_1intercept_only_1N_covs_1param_Xurand_thr_FALSECS_1het_sigma_0_applied_results.RDS"))
##
Model_B <- readRDS(file.path(
           "application_results_dummy_data",
           "seed_1_dummy_data_1intercept_only_1N_covs_1param_Xurand_thr_FALSECS_0het_sigma_0_applied_results.RDS"))
##
Model_C <- readRDS(file.path(
           "application_results_dummy_data",
           "seed_1_dummy_data_1intercept_only_1N_covs_1param_Xurand_thr_TRUECS_1het_sigma_0_applied_results.RDS"))
##
Model_D <- readRDS(file.path(
           "application_results_dummy_data",
           "seed_1_dummy_data_1intercept_only_1N_covs_1param_Xurand_thr_TRUECS_0het_sigma_0_applied_results.RDS"))
# ##
# ## ---- Model w/ covariates:
# ##
# Model_E <- readRDS(file.path(
#            "application_results_dummy_data",
#            "seed_1_dummy_data_1intercept_only_0N_covs_7param_Xurand_thr_TRUECS_0het_sigma_0_applied_results.RDS"))
##
# model_summary_and_trace_obj <- load_list$model_summary_and_trace_obj
# str(trace_gq)


summary_main <- Model_D$model_summary_and_trace_obj$get_summary_main()
summary_main

Model_D$model_summary_and_trace_obj$get_HMC_info()$L_main_during_sampling ## 166.5072
Model_D$model_summary_and_trace_obj$get_HMC_info()$eps ## 0.02857269

# Function 1: Extract single parameter in LaTeX format
extract_param_latex <- function(model_or_tibble, 
                                param_name, 
                                digits = 3) {
  
        # Check what we're dealing with
        if(is.data.frame(model_or_tibble)) {
          # It's already a tibble
          tibble <- model_or_tibble
        } else if("model_summary_and_trace_obj" %in% names(model_or_tibble)) {
          # It's your model structure!
          param_parts <- strsplit(param_name, "\\[")[[1]][1]
          tibble <- model_or_tibble$model_summary_and_trace_obj$extract_params(params = param_parts)
        } else if("extract_params" %in% names(model_or_tibble)) {
          # It's the R6 object directly
          param_parts <- strsplit(param_name, "\\[")[[1]][1]
          tibble <- model_or_tibble$extract_params(params = param_parts)
        } else {
          return("$-0.aaa ~ \\ci{0.aaa}{0.aaa}$")
        }
        
        # Filter using base R to avoid dplyr issues
        row <- tibble[tibble$parameter == param_name, ]
        
        if(nrow(row) == 0 || is.na(row$`50%`)) {
          return("$-0.aaa ~ \\ci{0.aaa}{0.aaa}$")
        }
        
        sprintf("$%.3f ~ \\ci{%.3f}{%.3f}$", 
                round(row$`50%`, digits), 
                round(row$`2.5%`, digits), 
                round(row$`97.5%`, digits))
  
}

# Function 2: Extract from multiple models and create table row
extract_param_row <- function(model_list, 
                              param_name, 
                              digits = 3) {
  
        values <- sapply(model_list, function(model) {
          if(is.null(model)) {
            return("$-0.aaa ~ \\ci{0.aaa}{0.aaa}$")
          }
          extract_param_latex(model, param_name, digits)
        })
        
        paste(values, collapse = " & ")
  
}

# Now your code should work:
model_list <- list(Model_A, 
                   Model_B, 
                   Model_C, 
                   Model_D)
##
cat(extract_param_row(model_list, "AUC[1]", digits = 3))
##
## For all AUCs:
##
for (i in 1:length(test_names)) {
  cat(test_names[i], "&", extract_param_row(model_list, paste0("AUC[", i, "]")), "\\\\\n")
}
##
## For all AUCs:
##
for (i in 1:length(test_names)) {
  cat(test_names[i], "&", extract_param_row(model_list, paste0("Se_baseline[", i, "]")), "\\\\\n")
}




##
##
##
# For beta_mu parameters:
for (i in 1:length(test_names)) {
  cat("Beta_mu", i, "&", extract_param_row(model_list, paste0("beta_mu[", i, ",1,1]")), "\\\\\n")
}
 
 
 


# Updated function to extract Se or Sp values WITH CIs
extract_se_sp_latex <- function(model_or_tibble, 
                                param_type = "Se_baseline", 
                                n_tests = 4,
                                max_n_thr = 63, 
                                digits = 3) {
        
        
        # Get the tibble
        if(is.data.frame(model_or_tibble)) {
          tibble <- model_or_tibble
        } else if("model_summary_and_trace_obj" %in% names(model_or_tibble)) {
          tibble <- model_or_tibble$model_summary_and_trace_obj$extract_params(params = param_type)
        } else {
          return(NULL)
        }
        
        # Create empty arrays for estimates and CIs
        est_matrix <- matrix(NA, nrow = n_tests, ncol = max_n_thr)
        lower_matrix <- matrix(NA, nrow = n_tests, ncol = max_n_thr)
        upper_matrix <- matrix(NA, nrow = n_tests, ncol = max_n_thr)
        
        # Fill the matrices
        for(i in 1:nrow(tibble)) {
          param <- tibble$parameter[i]
          # Extract indices from "Se_baseline[1,2]" format
          indices <- as.numeric(gsub(".*\\[([0-9]+),([0-9]+)\\]", "\\1,\\2", param) %>% 
                                  strsplit(",") %>% unlist())
          
          if(length(indices) == 2 && !is.na(tibble$`50%`[i]) && tibble$`50%`[i] != -1) {
            est_matrix[indices[1], indices[2]] <- tibble$`50%`[i]
            lower_matrix[indices[1], indices[2]] <- tibble$`2.5%`[i]
            upper_matrix[indices[1], indices[2]] <- tibble$`97.5%`[i]
          }
        }
        
        return(list(est = est_matrix, lower = lower_matrix, upper = upper_matrix))
  
}

# Updated function to format Se/Sp rows with CIs
format_se_sp_rows <- function(model_list, 
                              test_num, 
                              thresholds, 
                              n_tests,
                              max_n_thr,
                              threshold_labels,
                              param_type = "Se_baseline",
                              test_name = "Test",
                              digits = 3) {
  
        lines <- c()
        
        # Header
        param_label <- ifelse(grepl("Se", param_type), "Sensitivity", "Specificity")
        lines <- c(lines, sprintf("\\multicolumn{4}{l}{\\textit{\\textbf{%s (%s):}}} \\\\", 
                                  param_label, test_name))
        
        # Extract matrices for all models
        results_list <- lapply(model_list, function(model) {
          extract_se_sp_latex(model, 
                              n_tests = n_tests,
                              max_n_thr = max_n_thr,
                              param_type = param_type,
                              digits = digits)
        })
        
        # For each threshold
        for(i in seq_along(thresholds)) {
          thr <- thresholds[i]
          label <- threshold_labels[i]
          
          # Get values from each model WITH CIs
          values <- sapply(results_list, function(res) {
            if(is.null(res) || is.na(res$est[test_num, thr])) {
              return("$-0.aaa ~ \\ci{0.aaa}{0.aaa}$")
            } else {
              return(sprintf("$%.3f ~ \\ci{%.3f}{%.3f}$", 
                             res$est[test_num, thr],
                             res$lower[test_num, thr],
                             res$upper[test_num, thr]))
            }
          })
          
          # Format the row
          line <- sprintf("$%s_{%d, ~ \\ge %d}$ %s & %s \\\\", 
                          ifelse(grepl("Se", param_type), "Se", "Sp"),
                          test_num, thr, label,
                          paste(values, collapse = " & "))
          
          lines <- c(lines, line)
        }
        
        return(lines)
  
}

# Full function to do all tests at once:
create_all_se_sp_sections <- function(model_list,
                                      digits) {
        
        all_lines <- c()
        ##
        n_tests <- 4
        max_n_thr <- 63
        
        # Test configurations
        test_configs <- list(
          list(name = "GAD-2", num = 1, thresholds = c(2, 3, 4),
               labels = c("(screening thr.)", "(standard clinical thr.)", "(higher Sp thr.)")),
          list(name = "GAD-7", num = 2, thresholds = c(5, 10, 15),
               labels = c("(mild)", "(moderate)", "(severe)")),
          list(name = "HADS", num = 3, thresholds = c(8, 11, 15),
               labels = c("(possible; screening thr.)", "(probable; standard thr.)", "(severe)")),
          list(name = "BAI", num = 4, thresholds = c(8, 16, 22, 26, 36),
               labels = c("(mild thr.)", "(moderate thr.)", "(moderate-severe thr.)", 
                          "(severe thr.)", "(very severe thr.)"))
        )
        
        for(config in test_configs) {
          # Sensitivity
          se_lines <- format_se_sp_rows(  model_list = model_list, 
                                          test_num = config$num, 
                                          thresholds = config$thresholds,
                                          n_tests = n_tests,
                                          max_n_thr = max_n_thr,
                                          threshold_labels = config$labels, 
                                          param_type = "Se_baseline", 
                                          test_name = config$name,
                                          digits = digits)
          all_lines <- c(all_lines, se_lines)
          ##
          # Specificity
          sp_lines <- format_se_sp_rows(  model_list = model_list, 
                                          test_num = config$num, 
                                          thresholds = config$thresholds,
                                          n_tests = n_tests,
                                          max_n_thr = max_n_thr,
                                          threshold_labels = config$labels, 
                                          param_type = "Sp_baseline", 
                                          test_name = config$name,
                                          digits = digits)
          all_lines <- c(all_lines, sp_lines)
          
          # Add midrule between tests (except last)
          if(config$name != "BAI") {
            all_lines <- c(all_lines, "%%%%", "\\midrule", "%%%%")
          }
        }
        
        return(all_lines)
  
}



# The rest of the usage stays the same, but now includes CIs:
model_list <- list(Model_A, Model_B, Model_C, Model_D)


# GAD-2 Sensitivity
gad2_se <- format_se_sp_rows(
  model_list = model_list,
  n_tests = 4,
  max_n_thr = n_thr_max,
  test_num = 1,
  thresholds = c(2, 3, 4),
  threshold_labels = c("(screening thr.)", "(standard clinical thr.)", "(higher Sp thr.)"),
  param_type = "Se_baseline",
  test_name = "GAD-2"
)

# GAD-2 Specificity  
gad2_sp <- format_se_sp_rows(
  model_list = model_list,
  n_tests = 4,
  max_n_thr = n_thr_max,
  test_num = 1,
  thresholds = c(2, 3, 4),
  threshold_labels = c("(screening thr.)", "(standard clinical thr.)", "(higher Sp thr.)"),
  param_type = "Sp_baseline",
  test_name = "GAD-2"
)

# Print them
cat(paste(gad2_se, collapse = "\n"), "\n")
cat(paste(gad2_sp, collapse = "\n"), "\n")


# Get all Se/Sp sections at once:
all_se_sp <- create_all_se_sp_sections(model_list, digits = 3)
cat(paste(all_se_sp, collapse = "\n"))
 











 

##
Dirichlet_random_effects_type
priors$compound_symmetry
##
## ---- Run k-fold CV ----------------------------------------------------------------------------------------
##
# outs_kfold <- model_summary_and_trace_obj$run_k_fold_CV(  debugging = FALSE,
#                                                           ##
#                                                           K = 10,
#                                                           ##
#                                                           n_chains = 4,
#                                                           n_iter = 500,
#                                                           ##
#                                                           n_workers = 5,
#                                                           parallel = TRUE,
#                                                           ##
#                                                           use_BayesMVP_for_faster_summaries = TRUE)
# outs_kfold
# ##
# ##
# ##
# ##
source(file.path(getwd(), "R", "R_fns_misc.R"))
source(file.path(getwd(), "R", "R_fn_K_fold_CV.R"))
##
priors$compound_symmetry
##
gc() ; gc()
##
outs_kfold <- R_fn_run_kfold_cv_parallel(  debugging = FALSE,
                                           ##
                                           K = 10,
                                           ##
                                           model_prep_obj = model_prep_obj,
                                           stan_data_list = stan_data_list,
                                           ##
                                           priors = priors,
                                           init_lists_per_chain = init_lists_per_chain,
                                           ##
                                           n_burnin = n_burnin,
                                           # n_iter = 2000,
                                           n_iter = 500,
                                           n_chains = 8,
                                           adapt_delta = adapt_delta,
                                           max_treedepth = max_treedepth,
                                           ##
                                           n_workers = 10,
                                           output_dir = "cv_results",
                                           ##
                                           use_BayesMVP_for_faster_summaries = TRUE)
##
outs_kfold
str(outs_kfold)
                    




# 
# debugging = FALSE
# ##
# K = 10
# ##
# # model_prep_obj
# # stan_data_list
# ##
# priors = priors
# init_lists_per_chain = init_lists_per_chain
# ##
# n_burnin = n_burnin
# n_iter = 500
# n_chains = 4
# adapt_delta = adapt_delta
# max_treedepth = max_treedepth
# ##
# n_workers = 10
# output_dir = "cv_results"
# ##
# use_BayesMVP_for_faster_summaries = TRUE



##
## ---- Save K-fold results:
##
{
  if (length(test_names) == 4) { 
    dummy_data <- 1
    output_dir <- "application_results_dummy_data"
  } else { 
    dummy_data <- 0
    output_dir <- "application_results_real_data"
  }
  ##
  if (stan_data_list$n_covariates_max == 1) { 
    intercept_only <- 1
  } else { 
    intercept_only <- 0
  }
  ##
  file_name_string <- paste0( "K_fold_",  
                              "dummy_data_", dummy_data,
                              "intercept_only_", intercept_only,
                              "N_covs_", stan_data_list$n_covariates_max,
                              "param_", model_parameterisation,
                              "rand_thr_", random_thresholds,
                              "CS_", priors$compound_symmetry,
                              "het_sigma_", priors$hetero_sigma)
  ##
  # results_list_file_name <- paste0(file_name_string, "_applied_results.RDS")
  # results_list_path <- file.path(output_dir, results_list_file_name)
  ##
  # Save with overwrite check
  results_list_file <- file.path(output_dir, paste0("seed_", seed, "_", file_name_string, "_applied_results.RDS"))
  ##
  ## Remove if exists
  try({
    if (file.exists(results_list_file)) {
      cat(sprintf("Removing existing file: %s\n", results_list_file))
      file.remove(results_list_file)
    }
  })
  # Save:
  try({
    results_list_file <- file.path(output_dir, paste0("seed_", seed, "_", file_name_string, "_applied_results.RDS"))
    saveRDS(list("outs_kfold" = outs_kfold),
            results_list_file)
  })
  ##
  ## Verify
  try({
    if (!file.exists(results_list_file)) {
      stop(sprintf("Failed to save file: %s", results_list_file))
    } else {
      cat(sprintf("File saved successfully (size: %d bytes)\n",
                  file.info(results_list_file)$size))
    }
  })
}


##
## ---- Load K-fold CV results:
##
##
## ---- Model A: CS + FIXED C (DONE - 3000 iter):
##
file <- readRDS(file.path(output_dir, "seed_1_K_fold_dummy_data_1intercept_only_1N_covs_1param_Xurand_thr_FALSECS_1het_sigma_0_applied_results.RDS"))
file$outs_kfold$elpd_kfold
file$outs_kfold$se_elpd
##
# $outs_kfold
# $outs_kfold$elpd_kfold
# [1] -12263.17
# 
# $outs_kfold$se_elpd
# [1] 1712.946
# 
# $outs_kfold$K
# [1] 10
##
## ---- Model B: UN + FIXED C (currently running w/ 3000 iter):
##
file <- readRDS(file.path(output_dir, "seed_1_K_fold_dummy_data_1intercept_only_1N_covs_1param_Xurand_thr_FALSECS_0het_sigma_0_applied_results.RDS"))
file$outs_kfold$elpd_kfold
file$outs_kfold$se_elpd
# $outs_kfold
# $outs_kfold$elpd_kfold
# [1] -13061.45
# 
# $outs_kfold$se_elpd
# [1] 1527.98
# 
# $outs_kfold$K
# [1] 10
##
## ---- Model C: CS + RANDOM C (DONE - 3000 iter):
##
file <- readRDS(file.path(output_dir, "seed_1_K_fold_dummy_data_1intercept_only_1N_covs_1param_Xurand_thr_TRUECS_1het_sigma_0_applied_results.RDS"))
file$outs_kfold$elpd_kfold
file$outs_kfold$se_elpd
##
# $outs_kfold
# $outs_kfold$elpd_kfold
# [1] -13623.94
# 
# $outs_kfold$se_elpd
# [1] 1482.197
# 
# $outs_kfold$K
# [1] 10
##
## ---- Model D: UN + RANDOM C (DONE - 3000 iter):
##
file <- readRDS(file.path(output_dir, "seed_1_K_fold_dummy_data_1intercept_only_1N_covs_1param_Xurand_thr_TRUECS_0het_sigma_0_applied_results.RDS"))
file$outs_kfold$elpd_kfold
file$outs_kfold$se_elpd
# $outs_kfold
# $outs_kfold$elpd_kfold
# [1] -11030.59
# 
# $outs_kfold$se_elpd
# [1] 1936.577
# 
# $outs_kfold$K
# [1] 10
##
## ----  Best model = Model D (UN + RANDOM C)! 
##
##
## ----  Model E: Model D, but w/ all 4 (i.e. really 7!) covariates (currently running w/ 3000 iter):
##
file <- readRDS(file.path(output_dir, "seed_1_K_fold_dummy_data_1intercept_only_0N_covs_7param_Xurand_thr_TRUECS_0het_sigma_0_applied_results.RDS"))
file$outs_kfold$elpd_kfold
file$outs_kfold$se_elpd
##
## elpd_kfold = -18076.14
## se_elpd = 2283.006






# dplyr::filter(tibble_all, parameter %in% c("beta_sigma", "beta_corr")) %>% print()
# dplyr::filter(tibble_all, stringr::str_detect(parameter, "beta_tau")) %>% print(n=20)
# 
# 
# 
# 
# 
# dplyr::filter(tibble_all, stringr::str_detect(parameter, "beta_sigma")) %>% print(n=20)
# dplyr::filter(tibble_all, stringr::str_detect(parameter, "log_beta_sigma_MU")) %>% print(n=20)
# dplyr::filter(tibble_all, stringr::str_detect(parameter, "log_beta_sigma_SD")) %>% print(n=20)

# dplyr::filter(tibble_all, (stringr::str_detect(parameter, "diff_Se_flat"))) %>% print(n = 100)
# 
# ##
# ## Debugging:
# ##
# outs <- R_fn_sROC_plot_NMA(stan_model_file_name = stan_model_file_name,
#                    stan_mod_samples = stan_mod_samples,
#                    ##
#                    n_index_tests = n_index_tests,
#                    n_thr = n_thr_vec,
#                    ##
#                    test_names = test_names)
# ##
# outs$tibble_NMA %>% print(n = 100)
# outs$plot_list[[1]]
# outs$plot_list[[2]]
# outs$plot_list[[3]]
# outs$plot_list[[4]]
# outs$plot_list[[5]]
# 
# 
# 
stan_mod_samples <- model_summary_and_trace_obj$internal_obj$outs_stan_sampling$stan_mod_samples
fitted_model <- stan_mod_samples
stan_data_list <- model_prep_obj$internal_obj$outs_data$stan_data_list


# 
# # 
stan_functions_directory <- c("/home/enzocerullo/Documents/Work/PhD_work/R_packages/MetaOrdDTA/inst/stan_functions")
stan_model_file_path <- "/home/enzocerullo/Documents/Work/PhD_work/R_packages/MetaOrdDTA/inst/stan_models/stan_models_NMA/DTA_NMA_Nyaga_Xu_RANDthr_kappa.stan"
# # ##
new_cov_data <- list()
##
##
##
GQ_cov_results <- R_fn_GQ_covariates_NMA(  new_cov_data = new_cov_data,
                                           stan_mod_samples = stan_mod_samples,
                                           stan_data_list = stan_data_list,
                                           stan_model_file_path = stan_model_file_path,
                                           stan_functions_directory = stan_functions_directory)

str(GQ_cov_results)

draws <- GQ_cov_results$draws()

tibble_summary <- GQ_cov_results$summary()
tibble_summary
tibble_gq


str(draws)
str(trace_gq)
##
# Use fast Rcpp version
outs <-  R_fn_using_Rcpp_compute_NMA_comparisons_posthoc(  trace_gq = draws,
                                                                        n_thr = n_thr,
                                                                        test_names = test_names)
##
str(outs)
outs
##
## Then plot the results for these user-supplied values
##
outs <- R_fn_sROC_plot_NMA(stan_model_file_name = stan_model_file_name,
                           stan_mod_samples = stan_mod_samples,
                           ##
                           n_index_tests = n_index_tests,
                           n_thr = n_thr_vec,
                           ##
                           test_names = test_names)
##
outs$tibble_NMA %>% print(n = 100)
outs$plot_list[[1]]
# outs$plot_list[[2]]
# outs$plot_list[[3]]
# outs$plot_list[[4]]
# outs$plot_list[[5]]
# 
# model_prep_obj$internal_obj$outs_stan_compile
# model_summary_and_trace_obj$internal_obj$outs_stan_compile$stan_model_file_path# 
# model_summary_and_trace_obj$advanced_model_options$model_parameterisation

##
## ---- Extract NMA performance metrics:
##
tibble_NMA_performance_metrics <- model_summary_and_trace_obj$extract_NMA_performance()
tibble_NMA_performance_metrics
##
## ---- Extract NMA Se/Sp comparison metrics:
##
tibble_NMA_comparisons_metrics <- model_summary_and_trace_obj$extract_NMA_comparisons()
tibble_NMA_comparisons_metrics
##
## ---- Extract NMA AUC metrics:
##
tibble_NMA_AUC_metrics <- model_summary_and_trace_obj$extract_AUC()
tibble_NMA_AUC_metrics
##
## ----------- Plot sROC curve(s) (NMA - if WITHOUT covariates):
##
plots <- model_summary_and_trace_obj$plot_sROC(test_names = test_names)
##
plots$plot_list[[1]]
# plots$plot_list[[2]]
plots$plot_list[[3]]
# plots$plot_list[[4]]
plots$plot_list[[5]]
##
ggsave("Application_sim_baseline_sROC.png", plots$plot_list[[1]], width = 16, height = 9, dpi = 400)
ggsave("Application_sim_baseline_sROC_panel_w_CrI_PrI.png", plots$plot_list[[5]], width = 16, height = 9, dpi = 400)
# ggsave("Application_sim_baseline_sROC_w_CrI_PrI.png", coverage_plots[[6]], width = 16, height = 9, dpi = 400)
##
## ---- Plot sROC curve(s) (NMA - if WITH covariates):
##
X_summary <- R_fn_summarize_covariates_by_test(X_list = real_data$cov_data$X[[1]])
X_summary
##
beta_mu <- model_summary_and_trace_obj$extract_params(params = c("beta_mu")) %>% print(n = 100)
##
R_fn_map_beta_coefficients(beta_tibble = beta_mu, 
                           covariate_names = dimnames(X[[1]][[1]])[[2]], 
                           test_names = test_names, 
                           alpha = 0.05)
##
new_cov_data <- list()
##
baseline_case <- c(1,   ## intercept
                   ##
                   0.0, ## prev_GAD
                   ##
                   0,   ## low_RoB_QUADAS
                   ##
                   0,   ## Ref_test_SCID
                   0,   ## Ref_test_Structured
                   ##
                   0,   ## study_setting_2
                   1)   ## study_setting_3
##
new_cov_data$baseline_case_nd <- baseline_case
new_cov_data$baseline_case_d  <- baseline_case
##
##
if (model_parameterisation == "R&G") { 
  new_cov_data$baseline_case    <- rep(list(baseline_case), n_index_tests)
} else { 
  new_cov_data$baseline_case_nd <- rep(list(new_cov_data$baseline_case_nd), n_index_tests)
  new_cov_data$baseline_case_d  <- rep(list(new_cov_data$baseline_case_d), n_index_tests)
}
##
plots <- model_summary_and_trace_obj$plot_sROC(new_cov_data = new_cov_data, 
                                               test_names = test_names)
##
plots$plot_list[[1]]
plots$plot_list[[2]]
plots$plot_list[[3]]
plots$plot_list[[4]]
plots$plot_list[[5]]
##
## ---- Deviance:
##
# tibble_deviance <- model_summary_and_trace_obj$extract_params(tibble = tibble_all, 
#                                            params = c("deviance")) %>% print(n = 100)
# sum(tibble_deviance$`50%`, na.rm = TRUE)
# dplyr::filter(tibble_all, stringr::str_detect(parameter, "^deviance\\[\\d+\\]$")) %>% print(n = 100)
# dplyr::filter(tibble_all, stringr::str_detect(parameter, "deviance")) %>% 
#   pull(parameter) %>% 
#   unique()
tibble_deviance <- R_fn_deviance_summary( summary_tibble = tibble_gq, 
                               test_names = test_names)
tibble_deviance %>% print(n = 100)
##
sum(tibble_deviance$deviance_mean, na.rm = TRUE)

##
## ---- Study-specific accuracy estimates:
##
outs <- R_fn_extract_study_accuracy_summary(summary_tibble = tibble_gq, 
                                    n_studies = n_studies,
                                    n_thr = n_thr,
                                    test_names = test_names)
outs %>% print(n = 100)
##
## ---- Log-lik:
##
outs <- R_fn_log_lik_summary( summary_tibble = tibble_gq, 
                              n_studies = n_studies,
                              n_thr = n_thr,
                              test_names = test_names)
outs %>% print(n = 100)




extract_params_from_tibble_batch(tibble_all = tibble_all,
                                 param_strings_vec = c("deviance"))








trace_log_lik <- model_summary_and_trace_obj$get_trace_log_lik()

str(trace_log_lik)



posterior_draws <- model_summary_and_trace_obj$get_all_traces_list()
trace_gq <- posterior_draws$traces_as_arrays$trace_gq


log_lik_array <- filter_param_draws_string(named_draws_array = trace_gq,
                                           param_string = "log_lik_study",
                                  condition = "exact_match",
                                  exclude_NA = TRUE)

log_lik_array

str(log_lik_array)

require(loo)

loo::loo(x = log_lik_array)


waic(x = log_lik_array)





outs <- R_fn_using_Rcpp_compute_NMA_comparisons_posthoc(trace_gq = trace_gq,
                                     n_thr = n_thr_vec,
                                     test_names = test_names)

str(outs)

outs





outs_kfold <- model_summary_and_trace_obj$run_k_fold_CV(K = 10,
                                                        n_chains = 4,
                                                        n_workers = 10,
                                                        parallel = TRUE
                                                        )
outs_kfold






# 
# 
# # stan_fit <- stan_mod_samples
# # samples <- stan_mod_samples$draws()
# # 
# # 
# # str(stan_mod_samples)
# # str(samples)
# # 
# #  
# # 
# #   
# # test_names <- c("GAD-2", "BAI", "HADS", "GAD-7")
# # 
# # 
# # dplyr::filter(tibble_all, (stringr::str_detect(parameter, "diff_Se_flat"))) %>% print(n = 100)
# # nma_summary$comparisons_summary %>% print(n = 100)
# # 
# # 
# # 
# # comparisons_tibble %>% print(n = 100)
# # performance_tibble %>% print(n = 50)
# # 
# # 
# # nma_summary$comparisons_summary %>% print(n = 500)
# # nma_summary$performance_summary %>% print(n = 50)
# # 
# # 
# # 
# # NMA_summary <- extract_NMA_summary_direct( tibble_gq = tibble_gq, 
# #                                            n_thr = n_thr, 
# #                                            test_names = test_names)
# # 
# # NMA_summary$comparisons %>% print(width = 10000)
# # 
# # 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# sim_results$true_DGM_Fp
# sim_results$true_DGM_Sp
# sim_results$true_DGM_Se
# 
# 
# filter(outs$tibble_NMA, test == 1) %>% print(n = 7)
# sim_results$true_DGM_Se[[1]]
# sim_results$true_DGM_Sp[[1]]
# ##
# filter(outs$tibble_NMA, test == 2) %>% print(n = 22)
# sim_results$true_DGM_Se[[3]]
# sim_results$true_DGM_Sp[[3]]
# ##
# filter(outs$tibble_NMA, test == 3) %>% print(n = 22)
# sim_results$true_DGM_Se[[4]]
# sim_results$true_DGM_Sp[[4]]
# 
# ##
# filter(outs$tibble_NMA, test == 4) %>% print(n = 63)
# sim_results$true_DGM_Se[[2]]
# sim_results$true_DGM_Sp[[2]]
# 
# 
# 
# 
# ##
# 
# # # ##
# # R_fn_sROC_plot_MA(stan_model_file_name = stan_model_file_name,
# #                                stan_mod_samples = stan_mod_samples)
# # ##
# # plots <- model_summary_and_trace_obj$plot_sROC()
# # plots$plot_list[[1]]
# # plots$plot_list[[2]]
# 
# # model_summary_and_trace_obj$plot_traces(param_string = "alpha")
# 
# 
# {
#         n_thr_max <- max(model_prep_obj$internal_obj$outs_data$n_thr)
#         ##
#         Se_median_vec <- Sp_median_vec <- Fp_median_vec <- c()
#         Se_mean_vec   <- Sp_mean_vec   <- Fp_mean_vec   <- c()
#         ##
#         Se_lower_vec <- Sp_lower_vec <- Fp_lower_vec <- c()
#         Se_upper_vec <- Sp_upper_vec <- Fp_upper_vec <- c()
#         ##
#         Se_pred_lower_vec <- Sp_pred_lower_vec <- Fp_pred_lower_vec <- c()
#         Se_pred_upper_vec <- Sp_pred_upper_vec <- Fp_pred_upper_vec <- c()
#         
#         # str(Se_array)
#         
#         Se <- dplyr::filter(tibble_gq, (stringr::str_detect(parameter, "Se")) &
#                                      (!(stringr::str_detect(parameter, "Se_pred"))))
#         Sp <- dplyr::filter(tibble_gq, (stringr::str_detect(parameter, "Sp")) &
#                               (!(stringr::str_detect(parameter, "Sp_pred"))))
#         Fp  <- dplyr::filter(tibble_gq, (stringr::str_detect(parameter, "Fp")) &
#                               (!(stringr::str_detect(parameter, "Fp_pred"))))
#         ##
#         Se_pred <- dplyr::filter(tibble_gq, (stringr::str_detect(parameter, "Se_pred")))
#         Sp_pred <- dplyr::filter(tibble_gq, (stringr::str_detect(parameter, "Sp_pred")))
#         Fp_pred <- dplyr::filter(tibble_gq, (stringr::str_detect(parameter, "Fp_pred")))
#         
#         counter <- 1
#         for (k in 1:n_thr_max) {
#           for (t in 1:1) {
#               # test_vec[counter] <- t
#               ## Posterior medians (of pooled estimates):
#               Se_median_vec[k] <- 100 * Se$`50%`[counter]
#               Sp_median_vec[k] <- 100 * Sp$`50%`[counter]
#               Fp_median_vec[k] <- 100 * Fp$`50%`[counter]
#               ## Posterior means (of pooled estimates):
#               Se_mean_vec[k] <- 100 * Se$mean[counter]
#               Sp_mean_vec[k] <- 100 * Sp$mean[counter]
#               Fp_mean_vec[k] <- 100 * Fp$mean[counter]
#               ## Posterior lower 95% (of pooled estimates):
#               Se_lower_vec[k] <- 100 * Se$`2.5%`[counter]
#               Sp_lower_vec[k] <- 100 * Sp$`2.5%`[counter]
#               Fp_lower_vec[k] <- 100 * Fp$`2.5%`[counter]
#               ## Posterior upper 95% (of pooled estimates):
#               Se_upper_vec[k] <- 100 * Se$`97.5%`[counter]
#               Sp_upper_vec[k] <- 100 * Sp$`97.5%`[counter]
#               Fp_upper_vec[k] <- 100 * Fp$`97.5%`[counter]
#               ## Posterior lower prediction 95% (of pooled estimates):
#               Se_pred_lower_vec[k] <- 100 * Se_pred$`2.5%`[counter]
#               Sp_pred_lower_vec[k] <- 100 * Sp_pred$`2.5%`[counter]
#               Fp_pred_lower_vec[k] <- 100 * Fp_pred$`2.5%`[counter]
#               ## Posterior upper prediction 95% (of pooled estimates):
#               Se_pred_upper_vec[k] <- 100 * Se_pred$`97.5%`[counter]
#               Sp_pred_upper_vec[k] <- 100 * Sp_pred$`97.5%`[counter]
#               Fp_pred_upper_vec[k] <- 100 * Fp_pred$`97.5%`[counter]
#               ##
#               counter <- counter + 1
#           }
#         }
# }
#  
# # Se %>% print(n = 100) # Se_abs_diffs_sim <- Sp_abs_diffs_sim <- list()
# # Se_mean_of_diffs_sim <- Sp_mean_of_diffs_sim <- c()
# # Se_sum_of_diffs_sim <- Sp_sum_of_diffs_sim <- c()
# # Se_max_of_diffs_sim <- Sp_max_of_diffs_sim <- c()
# 
#  ESS_vec <- c(Se$n_eff, Sp$n_eff)
#  min_ESS <- min(Se$n_eff)
#   
# 
#   ##------------------------
#   {
#      vec_index_all         <-  1:n_thr
#         {
#           # x[[1]][, 10]
#           print(paste("seed = ", seed))
#           print(paste("model_parameterisation = ", model_parameterisation))
#           print(paste("box_cox = ", box_cox))
#           print(paste("softplus = ", softplus))
#           print(paste("n_studies = ", n_studies))
#           # print(paste("bs_het_C = ", bs_het_C))
#         }
#   }
# 
#  
#   if (min_ESS > 100) {
#    
#           ##
#           ## ---- Store results for this simulation:
#           ##
#           results_list[[sim_index]] <- list(
#             seed = seed,
#             n_studies = n_studies,
#             model_parameterisation = model_parameterisation,
#             box_cox = box_cox,
#             softplus = softplus,
#             random_thresholds = random_thresholds,
#             ##
#             true_DGM_Se = true_DGM_Se,
#             true_DGM_Sp = true_DGM_Sp,
#             ##
#             Se_mean_vec = Se_mean_vec,
#             Sp_mean_vec = Sp_mean_vec
#             ##
#            ##  model_summary = model_summary_and_trace_obj  # Optional: save full summary object
#           )
#           
#           ##
#           ## ---- Update summary dataframe:
#           ##
#           ##
#           # Create base values for the row
#           row_data <- list(
#             seed = seed,
#             n_studies = n_studies,
#             model_parameterisation = model_parameterisation,
#             box_cox = box_cox,
#             softplus = softplus,
#             random_thresholds = random_thresholds
#           )
#           ##
#           ## Add threshold-specific columns in a loop:
#           ##
#           for (k in vec_index_inner) {
#             row_data[[paste0("Se_mean_", k)]] <- Se_mean_vec[k]
#           }
#           ##
#           for (k in vec_index_inner) {
#             row_data[[paste0("Sp_mean_", k)]] <- Sp_mean_vec[k]
#           }
#           ##
#           ## Assign all data to the row:
#           ##
#           summary_df[sim_index, ] <- row_data
#           ##
#           ## ---- Save results for this simulation:
#           ##
#           file_name_string <- paste0( "test_", index_test_MA, 
#                                       "N_", N_per_study_mean,
#                                       "n_studies_", n_studies, 
#                                       ##
#                                       "param_", model_parameterisation, 
#                                       "soft_", softplus,
#                                       "rand_thr_", random_thresholds,
#                                       "alpha_pi_", alpha_pi,
#                                       ##
#                                       "sigma_C_nd", bs_het_C_nd_prob_scale,
#                                       "sigma_C_d", bs_het_C_d_prob_scale
#                                     )
#           
#           sim_result_file <- file.path(output_dir, paste0("seed_", seed, "_", file_name_string))
#           saveRDS(results_list[[sim_index]], sim_result_file)
#           ##
#           ## ---- Report progress:
#           ##
#           cat(sprintf("\n--- Simulation %d/%d completed ---\n", sim_index, n_sims))
#           # cat(sprintf("Se difference from true value: %.3f\n", Se_diff))
#           # cat(sprintf("Sp difference from true value: %.3f\n", Sp_diff))
#       
#   }
#   
#     
#   })
#   
#   } ##  -----------------------  end of N_{sims} loop  --------------------------------------------------------------------------------------
# 
#   # Calculate summary statistics across ALL simulations
#   summary_stats <- list(
#     n_sims = n_sims,
#     n_studies = n_studies,
#     model_parameterisation = model_parameterisation,
#     box_cox = box_cox,
#     softplus = softplus,
#     random_thresholds = random_thresholds,
#     ##
#     summary_df = summary_df
#   )
#   
#   
#   try({ 
#     # Save overall summary
#     summary_file <- file.path(getwd(), output_dir, paste0(file_name_string, ".RDS"))
#     saveRDS(summary_stats, summary_file)
#   })
#   
#   try({ 
#     # Also save as CSV for easy viewing
#     summary_csv <- file.path(getwd(), output_dir, paste0(file_name_string, ".csv"))
#     write.csv(summary_df, summary_csv, row.names = FALSE)
#   })
#   
# 
#   
#   print(summary_stats)
#   
#   
#   
#   require(beepr)
#   beepr::beep(sound = 2)
#   gc(reset = TRUE)
#   .rs.restartR()  # In RStudio
#   
#   
#   } ##  -----------------------  end of sim study  --------------------------------------------------------------------------------------
# 
# 
#   for (k in vec_index_inner) {
#     summary_stats$summary_df[, paste0("Bias_Se_", k)] <- summary_stats$summary_df[, paste0("Se_mean_", k)] - true_DGM_Se[k]
#   }
#   
# 
# source(file.path(getwd(), "inst", "examples", "thr_sim_study_helper_fns.R"))
# 
# summary_stats$summary_df
# 
# 
# message(paste("--------------------------------------------------"))
# 
# # For crayon
# library(crayon)
# options(crayon.enabled = TRUE)
# 
# cat(red("Does this show in red?"))
#  
# 
# 
#   ##
#   ## ---- 
#   ##
#   {
#         ##
#         ## ---- Compute average bias for thresholds in "vec_index_inner":
#         ##
#         {
#             result <- compute_average_bias(summary_df =   summary_stats$summary_df,
#                                            true_DGM_Se = true_DGM_Se,
#                                            true_DGM_Sp = true_DGM_Sp,
#                                            min_k = vec_index_inner[1],
#                                            max_k = tail(vec_index_inner, 1))
#             ##
#             cat(paste("\n--------------------------------------------------\n"))
#             ##
#             cat(sprintf("Mean Se bias: %s\n", formatC(result$mean_bias_Se, digits=3, format="g")))
#             cat(sprintf("Max Se bias: %s\n", formatC(result$max_bias_Se, digits=3, format="g")))
#             ##
#             # cat(sprintf("Mean Se SD: %s\n", formatC(result$mean_SD_Se, digits=3, format="g")))
#             # cat(sprintf("Max Se SD: %s\n", formatC(result$max_SD_Se, digits=3, format="g")))
#             ##
#             cat(sprintf("Mean Se MCSE(bias): %s\n", formatC(result$mean_MCSE_bias_Se, digits=3, format="g")))
#             cat(sprintf("Max Se MCSE(bias): %s\n", formatC(result$max_MCSE_bias_Se, digits=3, format="g")))
#             ##
#             cat(paste("\n--------------------------------------------------\n"))
#             ##
#             cat(sprintf("Mean Sp bias: %s\n", formatC(result$mean_bias_Sp, digits=3, format="g")))
#             cat(sprintf("Max Sp bias: %s\n", formatC(result$max_bias_Sp, digits=3, format="g")))
#             ##
#             # cat(sprintf("Mean Sp SD: %s\n", formatC(result$mean_SD_Sp, digits=3, format="g")))
#             # cat(sprintf("Max Sp SD: %s\n", formatC(result$max_SD_Sp, digits=3, format="g")))
#             ##
#             cat(sprintf("Mean Sp MCSE(bias): %s\n", formatC(result$mean_MCSE_bias_Sp, digits=3, format="g")))
#             cat(sprintf("Max Sp MCSE(bias): %s\n", formatC(result$max_MCSE_bias_Sp, digits=3, format="g")))
#             ##
#             ##
#             ##
#             invisible(cat(black(paste("\n model_parameterisation = ", model_parameterisation))))
#             invisible(cat(black(paste("\n Dirichlet_random_effects_type = ", Dirichlet_random_effects_type))))
#             cat(paste("\n--------------------------------------------------\n"))
#             ##
#             ##
#             # The function also returns the updated summary_df with all the Bias values
#             summary_df <- result$summary_df
#             ##
#             cat(paste("\n--------------------------------------------------\n"))
#             ##
#             invisible(cat(yellow(paste("\n n_sims = ", n_sims))))
#             cat(paste("\n--------------------------------------------------\n"))
#             ##
#             invisible(cat(blue(paste("\n index_test_MA = ", index_test_MA))))
#             cat(paste("\n--------------------------------------------------\n"))
#             ##
#             invisible(cat(red(paste("\n missing_thr = ", missing_thr))))
#             cat(paste("\n--------------------------------------------------\n"))
#             ##
#             ##
#             invisible(cat(red(paste("\n n_studies = ", n_studies))))
#             invisible(cat(red(paste("\n N_per_study_mean = ", N_per_study_mean))))
#             cat(paste("\n--------------------------------------------------\n"))
#             invisible(cat(blue(paste("\n N_per_study_SD = ", N_per_study_SD))))
#             cat(paste("\n--------------------------------------------------\n"))
#             ##
#             invisible(cat(green(paste("\n bivariate_locations_nd_bs_het_SD = ", scale_ALL_bs_SDs_by*bivariate_locations_nd_bs_het_SD))))
#             invisible(cat(green(paste("\n bivariate_locations_d_bs_het_SD = ", scale_ALL_bs_SDs_by*bivariate_locations_d_bs_het_SD))))
#             cat(paste("\n--------------------------------------------------\n"))
#             ##
#             invisible(cat(red(paste("\n bs_het_C_nd_prob_scale = ", scale_ALL_bs_SDs_by*bs_het_C_nd_prob_scale))))
#             invisible(cat(red(paste("\n bs_het_C_d_prob_scale = ", scale_ALL_bs_SDs_by*bs_het_C_d_prob_scale))))
#             cat(paste("\n--------------------------------------------------\n"))
#             ##
#             invisible(cat(green(paste("\n true_Mean_prev = ", true_Mean_prev))))
#             invisible(cat(green(paste("\n true_SD_probit_prev = ", true_SD_probit_prev))))
#             cat(paste("\n--------------------------------------------------\n"))
#             # ##
#             # invisible(cat(red(paste("bs_het_C_nd_prob_scale = ", bs_het_C_nd_prob_scale))))
#             # invisible(cat(red(paste("bs_het_C_d_prob_scale = ", bs_het_C_d_prob_scale))))
#             # cat(paste("\n--------------------------------------------------\n"))
#             ##
#             ##
#             ##
#             ##
#             cat(paste("\n--------------------------------------------------\n"))
#             cat(sprintf("Mean Se bias - ALL: %s\n", formatC(result$bias_Se_vec, digits=3, format="g")))
#             cat(paste("\n--------------------------------------------------\n"))
#             ##
#             cat(paste("\n--------------------------------------------------\n"))
#             cat(sprintf("Mean Sp bias - ALL: %s\n", formatC(result$bias_Sp_vec, digits=3, format="g")))
#             cat(paste("\n--------------------------------------------------\n"))
#         }
#         
#         
#         ##
#         ## ----
#         ##
#         ##
#         {
#           # cat(paste("\n--------------------------------------------------\n"))
#           # message("if min_MCSE_pct = 0.5%:")
#           # print(calc_min_sim_size(SD_Se_or_Sp_vec = result$SD_Se_vec, 
#           #                         min_MCSE_pct = 0.50, 
#           #                         min_k = vec_index_inner[1],
#           #                         max_k = tail(vec_index_inner, 1))) 
#           # cat(paste("\n--------------------------------------------------\n"))
#           cat(paste("\n--------------------------------------------------\n"))
#           message("if min_MCSE_pct = 0.25%:")
#           print(calc_min_sim_size(SD_Se_or_Sp_vec = result$SD_Se_vec, 
#                             min_MCSE_pct = 0.25, 
#                             min_k = vec_index_inner[1],
#                             max_k = tail(vec_index_inner, 1))) 
#           cat(paste("\n--------------------------------------------------\n"))
#           ##
#           ##
#           ##
#           cat(paste("\n--------------------------------------------------\n"))
#           message("if min_MCSE_pct = 0.125%:")
#           print(calc_min_sim_size(SD_Se_or_Sp_vec = result$SD_Se_vec, 
#                             min_MCSE_pct = 0.125, 
#                             min_k = vec_index_inner[1],
#                             max_k = tail(vec_index_inner, 1)))
#           cat(paste("\n--------------------------------------------------\n"))
#           ##
#           ##
#           ##
#           cat(paste("\n--------------------------------------------------\n"))
#           message("if min_MCSE_pct = 0.1%:")
#           print(calc_min_sim_size(SD_Se_or_Sp_vec = result$SD_Se_vec, 
#                                   min_MCSE_pct = 0.1, 
#                                   min_k = vec_index_inner[1],
#                                   max_k = tail(vec_index_inner, 1)))
#           cat(paste("\n--------------------------------------------------\n"))
#           ##
#           ##
#           ##
#           cat(paste("\n--------------------------------------------------\n"))
#           message("if min_MCSE_pct = 0.05%:")
#           print(calc_min_sim_size(SD_Se_or_Sp_vec = result$SD_Se_vec, 
#                             min_MCSE_pct = 0.05, 
#                             min_k = vec_index_inner[1],
#                             max_k = tail(vec_index_inner, 1)))
#           cat(paste("\n--------------------------------------------------\n"))
#         }
#   }

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
 
 
 
 




{
  
  
  
  basic_model_options = list(
    network = NULL,
    cts = NULL,
    prior_only = NULL
  )
  ####
  advanced_model_options = list(
    model_parameterisation = NULL,
    random_thresholds = NULL,
    Dirichlet_random_effects_type = NULL,
    box_cox = NULL,
    softplus = NULL
  )
  ####
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
  )
  ##
  ####
  priors = NULL
  ####
  # init_lists_per_chain = NULL,
  ####
  MCMC_params = list(
    seed = NULL,
    n_superchains = NULL,
    n_chains = NULL,
    n_iter = NULL,
    n_burnin = NULL,
    adapt_delta = NULL,
    max_treedepth = NULL,
    metric_shape = NULL
  )
  
  ##
  ## ---- Main "internal_obj" list:
  ##
  internal_obj = list(
    outs_data = NULL,
    outs_stan_model_name = NULL,
    outs_stan_compile = NULL,
    outs_stan_init = NULL,
    outs_stan_sampling = NULL,
    ##
    HMC_info = NULL, ## new
    efficiency_info = NULL, ## new
    summaries = NULL, ## new
    traces = NULL ## new
  )
  
  debugging <- TRUE
  ##
  n_iter = 500
  n_burnin = 500
  ##
  priors <- NULL
  ##
  x = x
  n_chains = n_chains
  ##
  cts = FALSE
  network = TRUE
  ##
  prior_only = FALSE
  ##
  softplus = FALSE
  ##
  ##
  model_parameterisation = "R&G"
  ##
  random_thresholds = TRUE
  Dirichlet_random_effects_type = "kappa"
  ##
  box_cox <- FALSE ## cannot input ACTUAL NA's into Stan!
  ##
  init_lists_per_chain = NULL
  
  
  basic_model_options$network <- network
  basic_model_options$cts <- cts
  basic_model_options$prior_only <- prior_only
  ##
  advanced_model_options$model_parameterisation <- model_parameterisation
  advanced_model_options$random_thresholds <- random_thresholds
  advanced_model_options$Dirichlet_random_effects_type <- Dirichlet_random_effects_type
  advanced_model_options$box_cox <- box_cox
  advanced_model_options$softplus <- softplus
  ##
  MCMC_params$n_chains <- n_chains
  
}





outs_data = list(
  stan_data_list = NULL,
  n_tests = NULL,
  n_studies = NULL,
  n_thr = NULL,
  n_cat = NULL
)
##
outs_stan_model_name = list(
  stan_model_file_name = NULL
)
##
outs_stan_compile = list(
  stan_model_obj = NULL,
  stan_model_file_name = NULL,
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
  stan_NMA_directory = NULL
)
##
outs_stan_init = list(
  inits_unconstrained_vec_per_chain = NULL,
  stan_param_names_list = NULL,
  stan_param_names_main = NULL,
  stan_init_pseudo_sampling_outs = NULL,
  stan_model_obj = NULL,
  json_file_path = NULL,
  stan_model_file_path = NULL
)
##
outs_stan_sampling = list(
  stan_mod_samples = NULL,
  time_total = NULL
)

X <- NULL

# x <- x_NMA

MetaOrdDTA:::prep_data_and_model(  
  debugging = self$debugging,
  ##
  x = self$x,
  ##
  X = self$X, ## covariates
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
  init_lists_per_chain = self$init_lists_per_chain,
  ##
  compute_sim_study_metrics = self$compute_sim_study_metrics,
  vec_index_inner_thr       = self$vec_index_inner_thr)







x = x
X = NULL
model_parameterisation = model_parameterisation
n_index_tests_per_study = n_index_tests_per_study
indicator_index_test_in_study = indicator_index_test_in_study



outs_data <- R_fn_prep_NMA_data(x = x,
                                X = NULL,
                                model_parameterisation = model_parameterisation,
                                n_index_tests_per_study = n_index_tests_per_study,
                                indicator_index_test_in_study = indicator_index_test_in_study)


internal_obj$outs_data <- outs_data




