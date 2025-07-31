



rm(list = ls())
.rs.restartR()  # In RStudio


# # # # # # # # # ##
# rm(list = ls())
# 
# # 
# 
# {
#       os <- .Platform$OS.type
# 
# 
#       if (os == "unix") {
#         user_root_dir <- Sys.getenv("PWD")
#       } else if (os == "windows") {
#         user_root_dir <- Sys.getenv("USERPROFILE")
#       }
#       local_pkg_dir <- file.path(user_root_dir, "Documents/Work/PhD_work/R_packages/MetaOrdDTA")
#       #
# 
# 
#       {
#         ## First remove any possible package fragments:
#         ## Find user_pkg_install_dir:
#         user_pkg_install_dir <- Sys.getenv("R_LIBS_USER")
#         print(paste("user_pkg_install_dir = ", user_pkg_install_dir))
#         ##
#         ## Find pkg_install_path + pkg_temp_install_path:
#         pkg_install_path <- file.path(user_pkg_install_dir, "MetaOrdDTA")
#         pkg_temp_install_path <- file.path(user_pkg_install_dir, "00LOCK-MetaOrdDTA")
#         ##
#         ## Remove any (possible) MetaOrdDTA package fragments:
#         remove.packages("MetaOrdDTA")
#         unlink(pkg_install_path, recursive = TRUE, force = TRUE)
#         unlink(pkg_temp_install_path, recursive = TRUE, force = TRUE)
#       }
# 
#       #
# 
# 
#       #
#       #### -------- ACTUAL (LOCAL) INSTALL:
#       ## document:Q
#       devtools::clean_dll(local_pkg_dir)
#       devtools::document(local_pkg_dir)
#       roxygen2::roxygenize(local_pkg_dir)
#       ##
#       ## Install (outer pkg):
#       ##
#       devtools::clean_dll(local_pkg_dir)
#       devtools::install(local_pkg_dir,
#                         upgrade = "never",
#                         quick = TRUE)
#       ##
#       ## May need to restart R:
#       ##
#       .rs.restartR()  # In RStudio
# 
#       # ?devtools::install
# 
#       require(MetaOrdDTA)
# }
# #


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


# {
#   model_parameterisation = "Xu"
#   box_cox <- FALSE
#   cts <- FALSE
#   ##
#   random_thresholds <-  TRUE
#   ##"
#   Dirichlet_random_effects_type <- "kappa"
# }
# # # # #
# # #
# # 



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
## 1-way (univariate) MR models:
##
# intercept_only <- FALSE ;  MR_model <- "MR_model_1" ## intercept + logit_prev_GAD
# intercept_only <- FALSE ;  MR_model <- "MR_model_2" ## intercept + ref_test
# intercept_only <- FALSE ;  MR_model <- "MR_model_3" ## intercept + study_setting
# intercept_only <- FALSE ;  MR_model <- "MR_model_4" ## intercept + RoB
# ##
# ## FULL MR model:
# ##
# intercept_only <- FALSE ;  MR_model <- "MR_model_5_FULL"  
####################### intercept_only <- FALSE ;  MR_model <- "MR_model_5_excl_RoB"   ### Same as MR model 12 !!!!!!!
# ##
# ## 2-way MR models:
# ##
intercept_only <- FALSE ;  MR_model <- "MR_model_6"   ## (intercept + prev_GAD + ref_test)   #### best model according to K-fold for GAD data !!!!!!!!!!!
# intercept_only <- FALSE ;  MR_model <- "MR_model_7"   ## (intercept + prev_GAD + study_setting)
# intercept_only <- FALSE ;  MR_model <- "MR_model_8"   ## (intercept + prev_GAD + RoB)
# intercept_only <- FALSE ;  MR_model <- "MR_model_9"   ## (intercept + ref_test + study_setting) ##  same as  "MR_model_5_excl_RoB"  !!!!
# intercept_only <- FALSE ;  MR_model <- "MR_model_10"  ## (intercept + ref_test + RoB)
# intercept_only <- FALSE ;  MR_model <- "MR_model_11"  ## (intercept + study_setting + RoB)
# ##
# ## 3-way MR models:
# ##
# intercept_only <- FALSE ;  MR_model <- "MR_model_12"  ## (intercept + prev_GAD + ref_test + study_setting)
# intercept_only <- FALSE ;  MR_model <- "MR_model_13"  ## (intercept + prev_GAD + ref_test + RoB)
# intercept_only <- FALSE ;  MR_model <- "MR_model_14"  ## (intercept + prev_GAD + study_setting + RoB)
# intercept_only <- FALSE ;  MR_model <- "MR_model_15"  ## (intercept + ref_test + study_setting + RoB)


# 
# intercept_only <- FALSE
# ##
# MR_model_set_1 <- c("MR_model_1", 
#                     "MR_model_2")
#  
# ##
# MR_model_set_2 <- c("MR_model_3", 
#                     "MR_model_4")
#                  
# ##
# MR_model_set_3 <- c("MR_model_5_FULL", 
#                     "MR_model_6")
#                    
# ##
# MR_model_set_4 <- c("MR_model_7", 
#                     "MR_model_8")
#                    
# ##
# MR_model_set_5 <- c("MR_model_9", 
#                     "MR_model_10")
#                    
# ##
# MR_model_set_6 <- c("MR_model_11", 
#                     "MR_model_12")
# ##
# MR_model_set_7 <- c("MR_model_13", 
#                     "MR_model_14")
# ##
# MR_model_set_8 <- c("MR_model_15")
# 




# for (MR_model in MR_model_set_8) {
# {
  
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
          set.seed(seed)
          ##
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
  
 source(file.path(getwd(), "R", "R_fn_covariates.R"))
 #x
 #X
  str(real_data$cov_data$X)
 validate_X_matrices(real_data$cov_data$X)
 validate_X_matrices(X)
 ##
 check_X_column_order(X)
 check_X_column_order(real_data$cov_data$X)
 
 # X[[1]][[1]][1:76, ] == real_data$cov_data$X[[1]][[1]]

 
 # # Debug: What's the actual structure?
 # which(rowSums(X[[1]][[4]] != 0) > 0)  # Which rows have data?
 # which(rowSums(original_cov_data$X[[1]][[4]] != 0) > 0)  # Which rows have data?
 # which(indicator_index_test_in_study[, 4] == TRUE)  # Which studies should have test 4?
 # 
 # ifelse(sim_results$indicator_index_test_in_study, 1, 0)
 # Print column names from your real data
 
 
 # 
 # print(colnames(real_data$cov_data$X[[1]][[1]]))
 # 
 # # Print column names from your simulated data  
 # print(colnames(X[[1]][[1]]))
 # 


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
  
        intercept_only <- TRUE
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
  
        X <- sim_results$X_mat
      
        if (MR_model == "MR_model_1") { 
          
                intercept_only <- FALSE
                X <- subset_X_covariates( X = X, 
                                          covariates_to_keep = c("intercept", 
                                                                 "logit_prev_GAD"))
                baseline_case <- c(1, 
                                   0.0)
            
        } else if (MR_model == "MR_model_2") {
          
                intercept_only <- FALSE
                X <- subset_X_covariates( X = X, 
                                          covariates_to_keep = c("intercept", 
                                                                 "Ref_test_clean_SCID", "Ref_test_clean_Structured"))
                baseline_case <- c(1, 
                                   0, 0)
            
        } else if (MR_model == "MR_model_3") {
          
                intercept_only <- FALSE
                X <- subset_X_covariates( X = X, 
                                          covariates_to_keep = c("intercept", 
                                                                 "study_setting_2", "study_setting_3"))
                baseline_case <- c(1, 
                                   0, 1)
            
        } else if (MR_model == "MR_model_4") {
          
                intercept_only <- FALSE
                X <- subset_X_covariates( X = X,
                                          covariates_to_keep = c("intercept",
                                                                 "low_RoB_QUADAS_clean"))
                baseline_case <- c(1,
                                   0)
            
        } else if (MR_model == "MR_model_5") {
          
                intercept_only <- FALSE
                X <- subset_X_covariates( X = X,
                                          covariates_to_keep = c("intercept",
                                                                 "logit_prev_GAD",
                                                                 "low_RoB_QUADAS_clean",
                                                                 "study_setting_2", "study_setting_3",
                                                                 "Ref_test_clean_SCID", "Ref_test_clean_Structured"))
                      baseline_case <- c(1, 
                                         0.0, 
                                         0,
                                         0, 1,
                                         0, 0)
            
        } else if (MR_model == "MR_model_5_excl_RoB") {
          
                intercept_only <- FALSE
                X <- subset_X_covariates( X = X, 
                                          covariates_to_keep = c("intercept", 
                                                                 "logit_prev_GAD",
                                                                 "study_setting_2", "study_setting_3",
                                                                 "Ref_test_clean_SCID", "Ref_test_clean_Structured"))
                baseline_case <- c(1,
                                   0.0, 
                                   0, 1,
                                   0, 0)
                
        } else if (MR_model == "MR_model_6") { ## (intercept + logit_prev_GAD + ref_test)
          
                intercept_only <- FALSE
                X <- subset_X_covariates( X = X,
                                          covariates_to_keep = c("intercept",
                                                                 "logit_prev_GAD",
                                                                 "Ref_test_clean_SCID", "Ref_test_clean_Structured"))
                
                baseline_case <- c(1, 
                                   0.0,
                                   0, 0)
          
        } else if (MR_model == "MR_model_7") { ## (intercept + logit_prev_GAD + study_setting)
          
                intercept_only <- FALSE
                X <- subset_X_covariates( X = X,
                                          covariates_to_keep = c("intercept",
                                                                 "logit_prev_GAD",
                                                                 "study_setting_2", "study_setting_3"))
                baseline_case <- c(1,
                                   0.0,
                                   0, 1)
          
        } else if (MR_model == "MR_model_8") {  ## (intercept + logit_prev_GAD + RoB)
          
                intercept_only <- FALSE
                X <- subset_X_covariates( X = X,
                                          covariates_to_keep = c("intercept",
                                                                 "logit_prev_GAD",
                                                                 "low_RoB_QUADAS_clean"))
                baseline_case <- c(1, 
                                   0.0, 
                                   0)
            
        } else if (MR_model == "MR_model_9") { ## (intercept + ref_test + study_setting)
          
                intercept_only <- FALSE
                X <- subset_X_covariates( X = X,
                                          covariates_to_keep = c("intercept",
                                                                 "study_setting_2", "study_setting_3",
                                                                 "Ref_test_clean_SCID", "Ref_test_clean_Structured"))
                baseline_case <- c(1, 
                                   0, 1,
                                   0, 0)
                
        } else if (MR_model == "MR_model_10") { ## (intercept + ref_test + RoB)
          
                intercept_only <- FALSE
                X <- subset_X_covariates( X = X,
                                          covariates_to_keep = c("intercept",
                                                                 "low_RoB_QUADAS_clean",
                                                                 "Ref_test_clean_SCID", "Ref_test_clean_Structured"))
                baseline_case <- c(1,
                                   0, 
                                   0, 0)
          
        } else if (MR_model == "MR_model_11") { ## (intercept + study_setting + RoB)
          
          intercept_only <- FALSE
                X <- subset_X_covariates( X = X,
                                          covariates_to_keep = c("intercept",
                                                                 "low_RoB_QUADAS_clean",
                                                                 "study_setting_2", "study_setting_3"))
                baseline_case <- c(1,
                                   0,
                                   0, 1)
          
        } else if (MR_model == "MR_model_12") { ## (intercept + logit_prev_GAD + ref_test + study_setting)
          
                intercept_only <- FALSE
                X <- subset_X_covariates( X = X,
                                          covariates_to_keep = c("intercept",
                                                                 "logit_prev_GAD",
                                                                 "study_setting_2", "study_setting_3",
                                                                 "Ref_test_clean_SCID", "Ref_test_clean_Structured"))
                baseline_case <- c(1, 
                                   0.0, 
                                   0, 1, 
                                   0, 0)
          
        } else if (MR_model == "MR_model_13") { ## (intercept + logit_prev_GAD + ref_test + RoB)
          
                intercept_only <- FALSE
                X <- subset_X_covariates( X = X,
                                          covariates_to_keep = c("intercept",
                                                                 "logit_prev_GAD",
                                                                 "low_RoB_QUADAS_clean",
                                                                 "Ref_test_clean_SCID", "Ref_test_clean_Structured"))
                baseline_case <- c(1, 
                                   0.0, 
                                   0, 
                                   0, 0)
          
        } else if (MR_model == "MR_model_14") { ## (intercept + logit_prev_GAD + study_setting + RoB)
          
                intercept_only <- FALSE
                X <- subset_X_covariates( X = X,
                                          covariates_to_keep = c("intercept",
                                                                 "logit_prev_GAD",
                                                                 "low_RoB_QUADAS_clean",
                                                                 "study_setting_2", "study_setting_3"))
                baseline_case <- c(1, 
                                   0.0,
                                   0,
                                   0, 1)
          
        } else if (MR_model == "MR_model_15") { ## (intercept + ref_test + study_setting + RoB)
          
                intercept_only <- FALSE
                X <- subset_X_covariates( X = X,
                                          covariates_to_keep = c("intercept",
                                                                 "low_RoB_QUADAS_clean",
                                                                 "study_setting_2", "study_setting_3",
                                                                 "Ref_test_clean_SCID", "Ref_test_clean_Structured"))
                baseline_case <- c(1,
                                   0,
                                   0, 1, 
                                   0, 0)
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
## validate_X_matrices(X)
validate_X_matrices(sim_results$X_mat)
X
##)
check_X_column_order(real_data$cov_data$X)



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
# # ##
# Total studies: 80 
# Studies with 1 test: 49 
# Studies with 2+ tests: 31 
# Study overlap proportion: 0.388 
# Total studies: 76 
# Studies with 1 test: 46 
# Studies with 2+ tests: 30 
# Study overlap proportion: 0.395 
# ===== INDICATOR MATRIX COMPARISON =====
#   
#   1. TEST FREQUENCIES:
#   Test Simulated  Real Difference
# GAD-2     31.2% 38.2%       6.9%
# GAD-7     50.0% 51.3%       1.3%
# HADS     50.0% 47.4%       2.6%
# BAI     21.2%  9.2%      12.0%
# 
# 2. OVERLAP PROPORTION:
#   Simulated: 38.8%
# Real: 39.5%
# Difference: 0.7%
# 
# 3. TESTS PER STUDY:
#   Simulated distribution:
#   
#   1    2    3 
# 61.3 25.0 13.8 
# Real distribution:
#   
#   1    2    3 
# 60.5 32.9  6.6 
# 
# 4. TOP TEST COMBINATIONS:
#   
#   Simulated (top 5):
#   HADS: 23.8%
# BAI: 21.2%
# GAD-2+GAD-7+HADS: 13.8%
# GAD-7+HADS: 12.5%
# GAD-2+GAD-7: 12.5%
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
# Test1  1.00  0.46 -0.08 -0.35
# Test2  0.46  1.00  0.05 -0.52
# Test3 -0.08  0.05  1.00 -0.52
# Test4 -0.35 -0.52 -0.52  1.00
# 
# Real:
#   [,1]  [,2]  [,3]  [,4]
# [1,]  1.00  0.60 -0.53 -0.16
# [2,]  0.60  1.00 -0.66 -0.14
# [3,] -0.53 -0.66  1.00 -0.21
# [4,] -0.16 -0.14 -0.21  1.00
# 
# ===== SIMILARITY SCORES =====
#   Frequency similarity: 94%
# Overlap similarity: 99%
# Correlation similarity: 64%
# Combination similarity: 69%
# 
# OVERALL SIMILARITY: 81%







##
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##
# 
# 
# ##
# cts <- box_cox <- softplus <- FALSE
# custom_file_name <- NULL ##  "DTA_NMA_Nyaga_Xu_RANDthr_kappa.stan"
# # ##
# model_parameterisation <- "Xu"
# Dirichlet_random_effects_type <- "kappa"
# random_thresholds <- TRUE
# ##
# # model_parameterisation <- "Xu"
# # Dirichlet_random_effects_type <- "fixed"
# # random_thresholds <- FALSE
# ##
# network <- TRUE
# ##
# 
# 
#  

 
 
 


##
## --- K-fold -----------------------------------------------------------------------------------------------

source(file.path(getwd(), "R", "R_fns_misc.R"))
source(file.path(getwd(), "R", "R_fn_K_fold_CV.R"))
##
##
##
## random_thresholds <- TRUE; Dirichlet_random_effects_type <- "kappa"
##
random_thresholds <- FALSE; Dirichlet_random_effects_type <- "fixed"
## 
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
  init_lists_per_chain = NULL
)



{
  init_lists_per_chain <- model_prep_obj$init_lists_per_chain
  priors <- model_prep_obj$priors
  stan_data_list <- model_prep_obj$internal_obj$outs_data$stan_data_list
}
##
init_lists_per_chain <- resize_init_list( init_lists_per_chain = init_lists_per_chain, 
                                          n_chains_new = n_chains)
length(init_lists_per_chain)
##
priors$compound_symmetry <- 0
model_prep_obj$priors$compound_symmetry <- priors$compound_symmetry
##
random_thresholds
priors$compound_symmetry 
intercept_only
##
## gc() ; gc()
RcppParallel::setThreadOptions(numThreads = 96);
##
priors
priors$compound_symmetry 
random_thresholds
##
## DO THIS ONCE AND SAVE IT!
set.seed(123, kind = "Mersenne-Twister")  # YOUR SEED
fold_assignments <- create_folds(
  K = 5,
  study_test_matrix = indicator_index_test_in_study,
  seed = 123
)

# Inspect it!
print(table(fold_assignments))
for (k in 1:5) {
  cat("Fold", k, ":", which(fold_assignments == k), "\n")
}
# Fold 1 : 3 8 13 15 19 21 26 29 30 34 41 42 48 53 56 58 64 65 
# Fold 2 : 7 12 14 23 31 44 52 54 59 60 63 66 68 69 72 73 79 
# Fold 3 : 1 9 10 11 17 33 36 39 40 45 47 51 71 74 75 77 
# Fold 4 : 4 5 16 18 22 24 27 28 32 37 38 49 55 61 67 78 
# Fold 5 : 2 6 20 25 35 43 46 50 57 62 70 76 80 
 

##
## ---- THIS MUST BE EQUAL TO 1:
##
sum(fold_assignments == c(3L, 5L, 1L, 4L, 4L, 5L, 2L, 1L, 3L, 3L, 3L, 2L, 1L, 2L, 1L, 
                          4L, 3L, 4L, 1L, 5L, 1L, 4L, 2L, 4L, 5L, 1L, 4L, 4L, 1L, 1L, 2L, 
                          4L, 3L, 1L, 5L, 3L, 4L, 4L, 3L, 3L, 1L, 1L, 5L, 2L, 3L, 5L, 3L, 
                          1L, 4L, 5L, 3L, 2L, 1L, 2L, 4L, 1L, 5L, 1L, 2L, 2L, 4L, 5L, 2L, 
                          1L, 1L, 2L, 4L, 2L, 2L, 5L, 3L, 2L, 2L, 3L, 3L, 5L, 3L, 4L, 2L, 
                          5L))/n_studies

# ##
# ## ---- Sample:
# ##
# n_chains <- 16
# ##
# if (random_thresholds == FALSE) { 
#   n_burnin <- 500
#   n_iter <- 500
# } else { 
#   n_burnin <- 1000
#   n_iter <- 1000
# }
# adapt_delta <- 0.65
# max_treedepth <- 10
# ##
# n_chains*n_iter
# ##
# init_lists_per_chain <- resize_init_list( init_lists_per_chain = init_lists_per_chain,
#                                           n_chains_new = n_chains)
# ##
# model_samples_obj <-  model_prep_obj$sample(
#   n_burnin = n_burnin,
#   n_iter   = n_iter,
#   adapt_delta = adapt_delta,
#   max_treedepth = max_treedepth,
#   metric_shape = "diag_e",
#   ##
#   priors = priors,
#   ##
#   n_chains = n_chains,
#   ##
#   init_lists_per_chain = init_lists_per_chain
#   ##
#   # cmdstanr_args = final_stan_args
# )
# 
# # #  
# # 


##    Sys.sleep(15*60)  # Blocks R for 10 minutes

# ##
# ## ---- And/or: run k-fold:
# ##
# source(file.path(getwd(), "R", "R_fns_misc.R"))
# source(file.path(getwd(), "R", "R_fn_K_fold_CV.R"))
# ##
# n_chains <- 2
# ##
# if (random_thresholds == FALSE) { 
#   n_burnin <- 500
#   n_iter <- 500
# } else { 
#   n_burnin <- 1000
#   n_iter <- 1000
# }
# ## Make sure have >= 2000 total iter for random-effects, and >=1000 total iter for fixed:
# ##
# random_thresholds
# n_iter*n_chains
# ##
# outs_kfold <- R_fn_run_kfold_cv_parallel(  debugging = FALSE,
#                                            ##
#                                            K = 5,
#                                            seed = 123,
#                                            ##
#                                            model_prep_obj = model_prep_obj,
#                                            stan_data_list = stan_data_list,
#                                            ##
#                                            fold_assignments = fold_assignments,
#                                            ##
#                                            cmdstanr_args = NULL,
#                                            ##
#                                            priors = priors,
#                                            init_lists_per_chain = init_lists_per_chain,
#                                            ##
#                                            n_burnin = n_burnin,
#                                            n_iter = n_iter,
#                                            ##
#                                            n_chains = n_chains,
#                                            ##
#                                            adapt_delta = 0.65,
#                                            max_treedepth = 10,
#                                            ##
#                                            n_workers = 5,
#                                            output_dir = "cv_results",
#                                            ##
#                                            use_BayesMVP_for_faster_summaries = TRUE)
#  
#  
# 
# ##
# ## ---- Save K-fold results:
# ##
# dummy_data <- 1
# ##
# {
#   if (dummy_data == 1) { 
#     output_dir <- "application_results_dummy_data"
#   } else { 
#     output_dir <- "application_results_real_data"
#   }
#   ##
#   if (stan_data_list$n_covariates_max == 1) { 
#     intercept_only <- 1
#   } else { 
#     intercept_only <- 0
#   }
#   ##
#   file_name_string <- paste0( "K_fold_",  
#                               "dummy_data_", dummy_data,
#                               "intercept_only_", intercept_only,
#                               "N_covs_", stan_data_list$n_covariates_max,
#                               "param_", model_parameterisation,
#                               "rand_thr_", random_thresholds,
#                               "CS_", priors$compound_symmetry,
#                               "het_sigma_", 0)
#   if (intercept_only == FALSE) {
#     file_name_string <- paste0(file_name_string, "_", MR_model)
#   }
#   ##
#   # results_list_file_name <- paste0(file_name_string, "_applied_results.RDS")
#   # results_list_path <- file.path(output_dir, results_list_file_name)
#   ##
#   # Save with overwrite check
#   results_list_file <- file.path(output_dir, paste0("seed_", seed, "_", file_name_string, "_applied_results.RDS"))
#   ##
#   ## Remove if exists
#   try({
#     if (file.exists(results_list_file)) {
#       cat(sprintf("Removing existing file: %s\n", results_list_file))
#       file.remove(results_list_file)
#     }
#   })
#   # Save:
#   try({
#     results_list_file <- file.path(output_dir, paste0("seed_", seed, "_", file_name_string, "_applied_results.RDS"))
#     saveRDS(list("outs_kfold" = outs_kfold),
#             results_list_file)
#   })
#   ##
#   ## Verify
#   try({
#     if (!file.exists(results_list_file)) {
#       stop(sprintf("Failed to save file: %s", results_list_file))
#     } else {
#       cat(sprintf("File saved successfully (size: %d bytes)\n",
#                   file.info(results_list_file)$size))
#     }
#   })
# }
# 
# ##
# gc() ; gc()
# #

# 
# # }
# 
# outs_kfold$min_ESS_vec
# 
# stan_data_list$n_covariates_max
# 
# ##
# priors$compound_symmetry
# random_thresholds
# 
# 



# model_names <- c("Intercept_only", "Prevalence", "Ref_test", "Setting")

## Note: accidently saved these w/ hetero_sigma = 1, eventhough it's 0 (+ redundant now anyway as removed from Stan models).
kfold_Model_A <- readRDS(file.path("application_results_dummy_data",
                                   "seed_1_K_fold_dummy_data_1intercept_only_1N_covs_1param_Xurand_thr_FALSECS_1het_sigma_0_applied_results.RDS"))
kfold_Model_B <- readRDS(file.path("application_results_dummy_data",
                                   "seed_1_K_fold_dummy_data_1intercept_only_1N_covs_1param_Xurand_thr_FALSECS_0het_sigma_0_applied_results.RDS"))
kfold_Model_C <- readRDS(file.path("application_results_dummy_data",
                                   "seed_1_K_fold_dummy_data_1intercept_only_1N_covs_1param_Xurand_thr_TRUECS_1het_sigma_0_applied_results.RDS"))
kfold_Model_D <- readRDS(file.path("application_results_dummy_data",
                                   "seed_1_K_fold_dummy_data_1intercept_only_1N_covs_1param_Xurand_thr_TRUECS_0het_sigma_0_applied_results.RDS"))


kfold_results_list <- list(
  kfold_Model_A,
  kfold_Model_B,
  kfold_Model_C,
  kfold_Model_D
)

model_names <- c(
                 "Model_A (FIXED-C + CS)", 
                 "Model_B (FIXED-C + UN)",
                 "Model_C (RAND-C + CS)", 
                 "Model_D (RAND-C + UN)")


kfold_Model_C$outs_kfold$min_ESS_vec
kfold_Model_D$outs_kfold$min_ESS_vec

# # For very conservative threshold (exclude only the terrible fold 9):
# filtered_100 <- filter_kfold_by_ess(outs_kfold, min_ess_threshold = 100) 
# 
# # For more aggressive filtering:
# filtered_400 <- filter_kfold_by_ess(outs_kfold, min_ess_threshold = 400)
# 
# 
# # Run comparison:
# compare_kfold_thresholds(outs_kfold)
# 



# 
# kfold_results_list[[1]]$outs_kfold$fold_assignments

comparison <- compare_kfold_models(kfold_results_list = kfold_results_list,
                                   model_names = model_names,
                                   min_ess_threshold = 100)
##
summarize_model_comparison(comparison)


##
## ---- Results (DUMMY data w/ 4 tests):
##
# +                                    min_ess_threshold = 100)
# Checking fold consistency across models...
# [1] 5
# [1] 5
# [1] 5
# ✓ All models used identical fold assignments
# 
# ESS summary across models:
#   Model_A (FIXED-C + CS) Model_B (FIXED-C + UN) Model_C (RAND-C + CS) Model_D (RAND-C + UN)
# [1,]                  364.1                  391.7                 149.8                 211.1
# [2,]                  470.5                  363.2                 208.1                 299.6
# [3,]                  474.5                  430.0                 355.8                  32.2
# [4,]                  450.0                  443.4                 350.3                 196.8
# [5,]                  513.7                  380.8                 127.1                 202.2
# 
# Minimum ESS per fold: 149.8 208.1 32.2 196.8 127.1 
# 
# Excluding 1 fold(s) with ESS < 100 in ANY model:
#   Fold 3: min ESS = 32.2
# 
# 
# Final comparison using 4/5 folds:
#   Model   ELPD   SE Delta_ELPD Delta_SE
# Model_A (FIXED-C + CS) Model_A (FIXED-C + CS) -19243 2222         NA       NA
# Model_B (FIXED-C + UN) Model_B (FIXED-C + UN) -20053 2225       -810     3145
# Model_C (RAND-C + CS)   Model_C (RAND-C + CS) -21555 2687      -2312     3487
# Model_D (RAND-C + UN)   Model_D (RAND-C + UN) -22029 2571      -2787     3398
# > ##
#   > summarize_model_comparison(comparison)
# 
# Model ranking by ELPD:
# 1. Model_A (FIXED-C + CS): ELPD = -19243 (2222) [BEST]
# 2. Model_B (FIXED-C + UN): ELPD = -20053 (2225), Δ = -810 (3145) 
# 3. Model_C (RAND-C + CS): ELPD = -21555 (2687), Δ = -2312 (3487) 
# 4. Model_D (RAND-C + UN): ELPD = -22029 (2571), Δ = -2787 (3398) 
# 
# * = Significantly different from best model (|Δ| > 2*SE)


##
## ---- Meta-regression k-fold:
##
kfold_MR_Model_1 <- readRDS(file.path("application_results_dummy_data",
                                      "seed_1_K_fold_dummy_data_1intercept_only_0N_covs_2param_Xurand_thr_FALSECS_1het_sigma_0_MR_model_1_applied_results.RDS"))
kfold_MR_Model_2 <- readRDS(file.path("application_results_dummy_data",
                                      "seed_1_K_fold_dummy_data_1intercept_only_0N_covs_3param_Xurand_thr_FALSECS_1het_sigma_0_MR_model_2_applied_results.RDS"))
kfold_MR_Model_3 <- readRDS(file.path("application_results_dummy_data",
                                      "seed_1_K_fold_dummy_data_1intercept_only_0N_covs_3param_Xurand_thr_FALSECS_1het_sigma_0_MR_model_3_applied_results.RDS"))
kfold_MR_Model_4 <- readRDS(file.path("application_results_dummy_data",
                                      "seed_1_K_fold_dummy_data_1intercept_only_0N_covs_2param_Xurand_thr_FALSECS_1het_sigma_0_MR_model_4_applied_results.RDS"))
##
kfold_MR_Model_5_FULL <- readRDS(file.path("application_results_dummy_data",
                                           "seed_1_K_fold_dummy_data_1intercept_only_0N_covs_7param_Xurand_thr_FALSECS_1het_sigma_0_MR_model_5_applied_results.RDS"))
# kfold_MR_Model_5_FULL_wo_RoB <- readRDS(file.path("application_results_dummy_data",
#                                            "seed_1_K_fold_dummy_data_1intercept_only_0N_covs_6param_Xurand_thr_FALSECS_1het_sigma_0_MR_model_5_applied_results.RDS"))
##
kfold_MR_Model_6 <- readRDS(file.path("application_results_dummy_data",
                                      "seed_1_K_fold_dummy_data_1intercept_only_0N_covs_4param_Xurand_thr_FALSECS_1het_sigma_0_MR_model_6_applied_results.RDS"))
kfold_MR_Model_7 <- readRDS(file.path("application_results_dummy_data",
                                      "seed_1_K_fold_dummy_data_1intercept_only_0N_covs_4param_Xurand_thr_FALSECS_1het_sigma_0_MR_model_7_applied_results.RDS"))
kfold_MR_Model_8 <- readRDS(file.path("application_results_dummy_data",
                                      "seed_1_K_fold_dummy_data_1intercept_only_0N_covs_3param_Xurand_thr_FALSECS_1het_sigma_0_MR_model_8_applied_results.RDS"))
kfold_MR_Model_9 <- readRDS(file.path("application_results_dummy_data",
                                      "seed_1_K_fold_dummy_data_1intercept_only_0N_covs_5param_Xurand_thr_FALSECS_1het_sigma_0_MR_model_9_applied_results.RDS"))
kfold_MR_Model_10 <- readRDS(file.path("application_results_dummy_data",
                                       "seed_1_K_fold_dummy_data_1intercept_only_0N_covs_4param_Xurand_thr_FALSECS_1het_sigma_0_MR_model_10_applied_results.RDS"))
kfold_MR_Model_11 <- readRDS(file.path("application_results_dummy_data",
                                       "seed_1_K_fold_dummy_data_1intercept_only_0N_covs_4param_Xurand_thr_FALSECS_1het_sigma_0_MR_model_11_applied_results.RDS"))
##
kfold_MR_Model_12 <- readRDS(file.path("application_results_dummy_data",
                                       "seed_1_K_fold_dummy_data_1intercept_only_0N_covs_6param_Xurand_thr_FALSECS_1het_sigma_0_MR_model_12_applied_results.RDS"))
kfold_MR_Model_13 <- readRDS(file.path("application_results_dummy_data",
                                       "seed_1_K_fold_dummy_data_1intercept_only_0N_covs_5param_Xurand_thr_FALSECS_1het_sigma_0_MR_model_13_applied_results.RDS"))
kfold_MR_Model_14 <- readRDS(file.path("application_results_dummy_data",
                                       "seed_1_K_fold_dummy_data_1intercept_only_0N_covs_5param_Xurand_thr_FALSECS_1het_sigma_0_MR_model_14_applied_results.RDS"))
kfold_MR_Model_15 <- readRDS(file.path("application_results_dummy_data",
                                       "seed_1_K_fold_dummy_data_1intercept_only_0N_covs_6param_Xurand_thr_FALSECS_1het_sigma_0_MR_model_15_applied_results.RDS"))



##
# # Model 1 (baseline meta-reg): 2 covariates
# intercept + prevalence
# 
# # Model 2: 3 covariates  
# intercept +  SCID + Structured
# 
# # Model 3: 3 covariates
# intercept +  setting_2 + setting_3
##
##

kfold_results_list <- list(
  kfold_Model_A,
  ##
  kfold_MR_Model_1,
  kfold_MR_Model_2,
  kfold_MR_Model_3,
  kfold_MR_Model_4,
  ##
  kfold_MR_Model_5_FULL,
  ##
  ## 2-way:
  ##
  kfold_MR_Model_6,
  kfold_MR_Model_7,
  kfold_MR_Model_8,
  kfold_MR_Model_9,
  kfold_MR_Model_10,
  kfold_MR_Model_11,
  ##
  ## 3-way:
  ##
  kfold_MR_Model_12,
  kfold_MR_Model_13,
  kfold_MR_Model_14,
  kfold_MR_Model_15
)


model_names <- c("MR_Model_0 (intercept-only - FIXED-C + CS)",
                 ##
                 "MR_Model_1 (intercept + prev_GAD)", 
                 "MR_Model_2 (intercept + ref_test)",
                 "MR_Model_3 (intercept + study_setting)",
                 "MR_Model_4 (intercept + RoB)",
                 ##
                 "MR_Model_FULL (int. + prev_GAD + ref_test + study_setting + RoB)",
                 # "MR_Model_FULL_wo_RoB (intercept + prev_GAD + ref_test + study_setting)", ## same as model 12
                 ##
                 ## first 2-way MR model as prev_GAD and ref_test gave best K-fold's out of univariate MR's:
                 ##
                 "MR_Model_6  (intercept + prev_GAD + ref_test)",
                 "MR_Model_7  (intercept + prev_GAD + study_setting)",
                 "MR_Model_8  (intercept + prev_GAD + RoB)",
                 "MR_Model_9  (intercept + ref_test + study_setting)",
                 "MR_Model_10 (intercept + ref_test + RoB)",
                 "MR_Model_11 (intercept + study_setting + RoB)",
                 ##
                 ## 3-way models:
                 ##
                 "MR_Model_12 (intercept + prev_GAD + ref_test + study_setting)",
                 "MR_Model_13 (intercept + prev_GAD + ref_test + RoB)",
                 "MR_Model_14 (intercept + prev_GAD + study_setting + RoB)",
                 "MR_Model_15 (intercept + ref_test + study_setting + RoB)"
)


# # For very conservative threshold (exclude only the terrible fold 9):
# filtered_100 <- filter_kfold_by_ess(kfold_results_list, min_ess_threshold = 100)
# 
# # For more aggressive filtering:
# filtered_400 <- filter_kfold_by_ess(kfold_results_list, min_ess_threshold = 400)
# 
# 
# # Run comparison:
# compare_kfold_thresholds(kfold_results_list$)

 

comparison <- compare_kfold_models(kfold_results_list = kfold_results_list,
                                   model_names = model_names,
                                   min_ess_threshold = 100)
##
summarize_model_comparison(comparison)
##
##
# Minimum ESS per fold: 326.4 379.9 346.8 313.5 319.1 
# 
# # 
# Final comparison using 5/5 folds:
# Model   ELPD   SE Delta_ELPD Delta_SE
# MR_Model_0 (intercept-only - FIXED-C + CS)                                             MR_Model_0 (intercept-only - FIXED-C + CS) -23640 2191         NA       NA
# MR_Model_1 (intercept + prev_GAD)                                                               MR_Model_1 (intercept + prev_GAD) -22111 2137       1529     3060
# MR_Model_2 (intercept + ref_test)                                                               MR_Model_2 (intercept + ref_test) -22655 1897        985     2898
# MR_Model_3 (intercept + study_setting)                                                     MR_Model_3 (intercept + study_setting) -24328 2232       -688     3128
# MR_Model_4 (intercept + RoB)                                                                         MR_Model_4 (intercept + RoB) -24091 2524       -450     3342
# MR_Model_FULL (int. + prev_GAD + ref_test + study_setting + RoB) MR_Model_FULL (int. + prev_GAD + ref_test + study_setting + RoB) -22244 1988       1396     2958
# MR_Model_6  (intercept + prev_GAD + ref_test)                                       MR_Model_6  (intercept + prev_GAD + ref_test) -21430 1969       2211     2945
# MR_Model_7  (intercept + prev_GAD + study_setting)                             MR_Model_7  (intercept + prev_GAD + study_setting) -22705 2247        936     3138
# MR_Model_8  (intercept + prev_GAD + RoB)                                                 MR_Model_8  (intercept + prev_GAD + RoB) -22600 2469       1041     3301
# MR_Model_9  (intercept + ref_test + study_setting)                             MR_Model_9  (intercept + ref_test + study_setting) -23236 1778        404     2822
# MR_Model_10 (intercept + ref_test + RoB)                                                 MR_Model_10 (intercept + ref_test + RoB) -23120 2276        520     3159
# MR_Model_11 (intercept + study_setting + RoB)                                       MR_Model_11 (intercept + study_setting + RoB) -24857 2540      -1217     3355
# MR_Model_12 (intercept + prev_GAD + ref_test + study_setting)       MR_Model_12 (intercept + prev_GAD + ref_test + study_setting) -21824 1880       1816     2887
# MR_Model_13 (intercept + prev_GAD + ref_test + RoB)                           MR_Model_13 (intercept + prev_GAD + ref_test + RoB) -21680 2178       1960     3089
# MR_Model_14 (intercept + prev_GAD + study_setting + RoB)                 MR_Model_14 (intercept + prev_GAD + study_setting + RoB) -23228 2584        412     3388
# MR_Model_15 (intercept + ref_test + study_setting + RoB)                 MR_Model_15 (intercept + ref_test + study_setting + RoB) -23746 2285       -105     3166
# > ##
#   > summarize_model_comparison(comparison)
# 
# Model ranking by ELPD:
# 1. MR_Model_6  (intercept + prev_GAD + ref_test): ELPD = -21430 (1969) [BEST]
# 2. MR_Model_13 (intercept + prev_GAD + ref_test + RoB): ELPD = -21680 (2178), Δ = 1960 (3089) 
# 3. MR_Model_12 (intercept + prev_GAD + ref_test + study_setting): ELPD = -21824 (1880), Δ = 1816 (2887) 
# 4. MR_Model_1 (intercept + prev_GAD): ELPD = -22111 (2137), Δ = 1529 (3060) 
# 5. MR_Model_FULL (int. + prev_GAD + ref_test + study_setting + RoB): ELPD = -22244 (1988), Δ = 1396 (2958) 
# 6. MR_Model_8  (intercept + prev_GAD + RoB): ELPD = -22600 (2469), Δ = 1041 (3301) 
# 7. MR_Model_2 (intercept + ref_test): ELPD = -22655 (1897), Δ = 985 (2898) 
# 8. MR_Model_7  (intercept + prev_GAD + study_setting): ELPD = -22705 (2247), Δ = 936 (3138) 
# 9. MR_Model_10 (intercept + ref_test + RoB): ELPD = -23120 (2276), Δ = 520 (3159) 
# 10. MR_Model_14 (intercept + prev_GAD + study_setting + RoB): ELPD = -23228 (2584), Δ = 412 (3388) 
# 11. MR_Model_9  (intercept + ref_test + study_setting): ELPD = -23236 (1778), Δ = 404 (2822) 
# 12. MR_Model_0 (intercept-only - FIXED-C + CS): ELPD = -23640 (2191), Δ = NA (NA) NA
# 13. MR_Model_15 (intercept + ref_test + study_setting + RoB): ELPD = -23746 (2285), Δ = -105 (3166) 
# 14. MR_Model_4 (intercept + RoB): ELPD = -24091 (2524), Δ = -450 (3342) 
# 15. MR_Model_3 (intercept + study_setting): ELPD = -24328 (2232), Δ = -688 (3128) 
# 16. MR_Model_11 (intercept + study_setting + RoB): ELPD = -24857 (2540), Δ = -1217 (3355) 
# 
# * = Significantly different from best model (|Δ| > 2*SE)





##
## ---- For REAL data (NOTE: THIS MAY CHANGE WHEN INCLUDE ALL 8 TESTS!!! AS MORE STUDIES TOO!!!)  -------------------------------------------------
##
# "While study setting is often considered in diagnostic meta-analyses, our k-fold cross-validation 
# revealed it actually degraded predictive performance (ΔELPD = -973), suggesting these broad categorizations
# may not capture meaningful heterogeneity."
##
# Of four pre-specified covariates, only reference test type improved predictive performance (ΔELPD = +416). 
# Disease prevalence had negligible effect (+23), while study setting (-973) and risk of bias (-288) substantially
# degraded predictions, suggesting these categorical variables introduce noise rather than explain
# meaningful heterogeneity
##
## ---- How to explain sim study reasoning in paper:
##
# "In the simulated example (generated from a model with all covariates having non-zero effects), k-fold CV correctly
# identified that disease prevalence improved predictive performance (ΔELPD = 3557), unlike in real data where it provided 
# negligible benefit. This demonstrates our methodology's ability to detect dataset-specific covariate importance."
##
##
# You're exactly right about the real data - stop there and use intercept + reference test as your main model!
# For the dummy data, I'd recommend this strategy:
#   Skip RoB (only 128 ELPD) and work with the 3 useful covariates:
#   Option 1: Incremental approach (I recommend this)
# 
# You already have the univariate results
# Fit one model: intercept + prev_GAD + reference_test + study_setting
# See if it beats prev_GAD alone (which was your best univariate)
# 
# My recommendation: Go with Option 1 because:
#   
#   K-fold already told you which covariates help individually
# You're demonstrating the methodology, not doing exhaustive model selection
# It mirrors what you'd do in practice: "these 3 helped, let's try them together"
# 
# Then for your final models:
#   Real data: intercept + reference_test
# 
# SROC curves by reference test type
# Se/Sp tables/plots
# AUC differences between tests
# 
# Dummy data: intercept + prev_GAD + reference_test + study_setting (if it wins)
# 
# Or just prev_GAD if the combined model doesn't improve
# Same outputs as above
# 
# This gives you a clean story: "We used k-fold to identify useful covariates, combined those that helped,
## and present results for the best model." Perfect for a thesis demonstration! 🎯
# 


##
## ----  Summarise + output results: -----------------------------------------------------------------------------------------------------------------------
##
model_summary_and_trace_obj <- model_samples_obj$summary(
  compute_main_params = TRUE,
  compute_transformed_parameters = TRUE, 
  compute_generated_quantities = TRUE,
  ##
  save_log_lik_trace = FALSE,
  ##
  use_BayesMVP_for_faster_summaries = TRUE,
  ##
  compute_nested_rhat = FALSE)
##
test_names = c("GAD-2", "GAD-7", "HADS", "BAI")
##
tibble_main <- model_summary_and_trace_obj$get_summary_main() %>% print(n = 200)
tibble_tp   <- model_summary_and_trace_obj$get_summary_transformed() %>% print(n = 100)
tibble_gq   <- model_summary_and_trace_obj$get_summary_generated_quantities() %>% print(n = 1000)
##
## tibble_all <- rbind(tibble_main, tibble_gq)
tibble_all <- rbind(tibble_main, tibble_tp, tibble_gq)
##
Se <- model_summary_and_trace_obj$extract_params(params = c("Se_baseline")) %>% print(n = 20)
Sp <- model_summary_and_trace_obj$extract_params(params = c("Sp_baseline")) %>% print(n = 20)
Fp <- model_summary_and_trace_obj$extract_params(params = c("Fp_baseline")) %>% print(n = 20)
min(c(Se$n_eff, Sp$n_eff), na.rm = TRUE)
##
## ---- Extract NMA AUC metrics:
##
# tibble_NMA_AUC_metrics <- model_summary_and_trace_obj$extract_AUC(test_names = test_names)
# tibble_NMA_AUC_metrics
# ##
# dplyr::filter(tibble_all, (stringr::str_detect(parameter, "AUC"))) %>% print(n = 100)
# ##
# model_summary_and_trace_obj$extract_params(params = c("kappa")) %>% print(n = 20)
# ##
# intercept_only
# 
# 
# 
# model_summary_and_trace_obj$get_efficiency_metrics()
# 


## w/ log(200/50), df = 15, log_sd = 0.75 priors:
# A tibble: 8 × 10
# param_group parameter   mean     sd `2.5%` `50%` `97.5%` n_eff  Rhat n_Rhat
# <chr>       <chr>      <dbl>  <dbl>  <dbl> <dbl>   <dbl> <dbl> <dbl> <lgl> 
#   1 kappa       kappa[1,1] 105.   24.8    65.9 102.    160.   1179  1.01 NA    
# 2 kappa       kappa[2,1]  25.0   5.51   16.2  24.2    37.3  1379  1.01 NA    
# 3 kappa       kappa[1,2] 202.   31.7   149.  199.    273.    169  1.07 NA    
# 4 kappa       kappa[2,2]  97.5  20.3    66.1  94.5   144.    291  1.04 NA    
# 5 kappa       kappa[1,3] 374.   71.1   262.  364.    540.    147  1.08 NA    
# 6 kappa       kappa[2,3] 260.   97.6   131.  237.    502.    101  1.11 NA    
# 7 kappa       kappa[1,4] 688.  240.    376.  636.   1293.    180  1.07 NA    
# 8 kappa       kappa[2,4] 974.  560.    347.  814.   2412.     35  1.38 NA 









##
## ---- Save full model output:
##
dummy_data <- 1
##
{
  if (dummy_data == 1) {
    output_dir <- "application_results_dummy_data"
  } else {
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
                              "het_sigma_", 0)
  if (intercept_only == FALSE) {
    file_name_string <- paste0(file_name_string, "_", MR_model)
  }
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
    results_list_file <- file.path(getwd(), output_dir, paste0("seed_", seed, "_", file_name_string, "_applied_results.RDS"))
    saveRDS(object = list("model_prep_obj" = model_prep_obj,
                          "model_summary_and_trace_obj" = model_summary_and_trace_obj),
            file = results_list_file)
  })
}




# 
# MetaOrdDTA:::extract_params_from_tibble_batch(debugging = FALSE,
#                                               tibble = tibble_gq,
#                                               param_strings_vec = c("Se"),
#                                               condition = "containing") %>% print(n=100)
# 
# 
# 
# 
# MetaOrdDTA:::extract_params_from_tibble_batch(debugging = FALSE,
#                                               tibble = tibble_all,
#                                               param_strings_vec = c("beta"),
#                                               condition = "containing") %>% print(n=100)
# 
# 



# ##
# ## ----------- Plot sROC curve(s) (NMA - if WITHOUT covariates):
# ##
# plots <- model_summary_and_trace_obj$plot_sROC(test_names = test_names)
# ##
# plots$plot_list[[1]]
# ##
# ##
# ## ---- Plot sROC curve(s) (NMA - if WITH covariates):
# ##
# new_cov_data <- list()
# ##
# baseline_case <- c(1,   ## intercept
#                    ##
#                    0.0, ## prev_GAD
#                    ##
#                    0,   ## low_RoB_QUADAS
#                    ##
#                    0,   ## Ref_test_SCID
#                    0,   ## Ref_test_Structured
#                    ##
#                    0,   ## study_setting_2
#                    1)   ## study_setting_3
# ##
# new_cov_data$baseline_case_nd <- baseline_case
# new_cov_data$baseline_case_d  <- baseline_case
##
##
# if (model_parameterisation == "R&G") { 
# #   new_cov_data$baseline_case    <- rep(list(new_cov_data$baseline_case), n_index_tests)
# # } else { 
#   new_cov_data$baseline_case_nd <- rep(list(new_cov_data$baseline_case_nd), n_index_tests)
#   new_cov_data$baseline_case_d  <- rep(list(new_cov_data$baseline_case_d), n_index_tests)
# # }
# ##
# plots <- model_summary_and_trace_obj$plot_sROC(new_cov_data = new_cov_data, 
#                                                test_names = test_names)
# ##
# plots$plot_list[[1]]
# plots$plot_list[[2]]
# plots$plot_list[[3]]
# plots$plot_list[[4]]
# plots$plot_list[[5]]
# 
# 
# X_summary <- R_fn_summarize_covariates_by_test(X_list = real_data$cov_data$X[[1]])
# X_summary
# ##
# R_fn_generate_sroc_scenarios(summary_results = X_summary, 
#                              max_scenarios = 10)



model_summary_and_trace_obj$get_efficiency_metrics()$time_total




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
                  Dirichlet_random_effects_type = Dirichlet_random_effects_type)
 


{
  init_lists_per_chain <- model_prep_obj$init_lists_per_chain
  priors <- model_prep_obj$priors
  stan_data_list <- model_prep_obj$internal_obj$outs_data$stan_data_list
}
##
# n_total_obs = 0;
# for (t in 1:n_index_tests) {
#   for (s in 1:n_studies) {
#     if (indicator_index_test_in_study[s, t] == 1) {
#       n_total_obs  = n_total_obs +  1;
#     }
#   }
# }
# ##
# n_nuisance_z <- 2*n_total_obs + 2*n_studies; n_nuisance_z
# # ##
# n_thr <- stan_data_list$n_thr
# ##
# model_parameterisation
# random_thresholds
# ##
# priors$compound_symmetry <- 0
# priors$hetero_sigma      <- 0
# 
# 
# ##
# ## ----  Sample model: ----------------------------------------------------------------
# ##
# if (random_thresholds == FALSE) {
#     n_iter <- 500
# } else { 
#     n_iter <- 1000
# }
# ##
# n_iter
# # ##
#   # model_samples_obj <-  model_prep_obj$sample(
#   #                            n_burnin = n_burnin,
#   #                            n_iter   = n_iter,
#   #                            adapt_delta = adapt_delta,
#   #                            max_treedepth = max_treedepth,
#   #                            metric_shape = "diag_e",
#   #                            ##
  #                            priors = priors,
  #                            ##
  #                            n_chains = n_chains,
  #                            ##
  #                            init_lists_per_chain = init_lists_per_chain
  #                            )

 
 

 
 

#  
# ##
# ## ----------- Plot sROC curve(s) (NMA - if WITHOUT covariates):
# ##
# plots <- model_summary_and_trace_obj$plot_sROC(test_names = test_names)
# ##
# plots$plot_list[[1]]
# # plots$plot_list[[2]]
# plots$plot_list[[3]]
# # plots$plot_list[[4]]
# plots$plot_list[[5]]
# 
# 
#  
# 
# 
#  
 
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




# The rest of the usage stays the same, but now includes CIs:
model_list <- list(Model_A,
                   Model_B, 
                   Model_C, 
                   Model_D)

 

 
##
## For all AUCs:
##
for (i in 1:length(test_names)) {
  cat(test_names[i], "&", extract_param_row(model_list, paste0("AUC[", i, "]")), "\\\\\n")
}
 
##
## Get all Se/Sp sections at once:
##
all_se_sp <- create_all_se_sp_sections(model_list, digits = 3)
cat(paste(all_se_sp, collapse = "\n"))

##
## Get all estimates for table 2:
##
table2_lines <- create_table2_latex(model_list, digits = 3)
cat(paste(table2_lines, collapse = "\n"))


##
# ##
# ##
# # For beta_mu parameters:
# for (i in 1:length(test_names)) {
#   cat("Beta_mu", i, "&", extract_param_row(model_list, paste0("beta_mu[", i, ",1,1]")), "\\\\\n")
# }
# 
# 
# 






 

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
# source(file.path(getwd(), "R", "R_fns_misc.R"))
# source(file.path(getwd(), "R", "R_fn_K_fold_CV.R"))
# ##
# priors$compound_symmetry
# ##
# gc() ; gc()
# ##
# outs_kfold <- R_fn_run_kfold_cv_parallel(  debugging = FALSE,
#                                            ##
#                                            K = 10,
#                                            ##
#                                            model_prep_obj = model_prep_obj,
#                                            stan_data_list = stan_data_list,
#                                            ##
#                                            priors = priors,
#                                            init_lists_per_chain = init_lists_per_chain,
#                                            ##
#                                            n_burnin = n_burnin,
#                                            # n_iter = 2000,
#                                            n_iter = 500,
#                                            n_chains = 8,
#                                            adapt_delta = adapt_delta,
#                                            max_treedepth = max_treedepth,
#                                            ##
#                                            n_workers = 10,
#                                            output_dir = "cv_results",
#                                            ##
#                                            use_BayesMVP_for_faster_summaries = TRUE)
# ##
# outs_kfold
# str(outs_kfold)
#                     
# 
#  
# 
# 
# ##
# ## ---- Save K-fold results:
# ##
# {
#   if (length(test_names) == 4) { 
#     dummy_data <- 1
#     output_dir <- "application_results_dummy_data"
#   } else { 
#     dummy_data <- 0
#     output_dir <- "application_results_real_data"
#   }
#   ##
#   if (stan_data_list$n_covariates_max == 1) { 
#     intercept_only <- 1
#   } else { 
#     intercept_only <- 0
#   }
#   ##
#   file_name_string <- paste0( "K_fold_",  
#                               "dummy_data_", dummy_data,
#                               "intercept_only_", intercept_only,
#                               "N_covs_", stan_data_list$n_covariates_max,
#                               "param_", model_parameterisation,
#                               "rand_thr_", random_thresholds,
#                               "CS_", priors$compound_symmetry,
#                               "het_sigma_", priors$hetero_sigma)
#   ##
#   # results_list_file_name <- paste0(file_name_string, "_applied_results.RDS")
#   # results_list_path <- file.path(output_dir, results_list_file_name)
#   ##
#   # Save with overwrite check
#   results_list_file <- file.path(output_dir, paste0("seed_", seed, "_", file_name_string, "_applied_results.RDS"))
#   ##
#   ## Remove if exists
#   try({
#     if (file.exists(results_list_file)) {
#       cat(sprintf("Removing existing file: %s\n", results_list_file))
#       file.remove(results_list_file)
#     }
#   })
#   # Save:
#   try({
#     results_list_file <- file.path(output_dir, paste0("seed_", seed, "_", file_name_string, "_applied_results.RDS"))
#     saveRDS(list("outs_kfold" = outs_kfold),
#             results_list_file)
#   })
#   ##
#   ## Verify
#   try({
#     if (!file.exists(results_list_file)) {
#       stop(sprintf("Failed to save file: %s", results_list_file))
#     } else {
#       cat(sprintf("File saved successfully (size: %d bytes)\n",
#                   file.info(results_list_file)$size))
#     }
#   })
# }
# 
# 
#  








Model_RAND <- readRDS(
  file.path("application_results_dummy_data", 
            "seed_1_dummy_data_1intercept_only_1N_covs_1param_Xurand_thr_TRUECS_0het_sigma_0_applied_results.RDS"))




model_summary_and_trace_obj <- Model_RAND$model_summary_and_trace_obj
model_prep_obj <- Model_RAND$model_prep_obj

##
## ---- Extract results + plots for FINAL INTERCEPT-ONLY chosen model (i.e. FIXED-cutpoints + CS) -------------------------------------------------------------------------
##
Model_dummy <- readRDS(
  file.path("application_results_dummy_data", 
            "seed_1_dummy_data_1intercept_only_1N_covs_1param_Xurand_thr_FALSECS_1het_sigma_0_applied_results.RDS"))

model_summary_and_trace_obj_DUMMY <- Model_dummy$model_summary_and_trace_obj


test_names <- c("GAD-2", "GAD-7", "HADS", "BAI")



# 
# # # 
# stan_mod_samples <- model_summary_and_trace_obj_DUMMY$internal_obj$outs_stan_sampling$stan_mod_samples
# fitted_model <- stan_mod_samples
# stan_data_list <- model_summary_and_trace_obj_DUMMY$internal_obj$outs_data$stan_data_list
# 



Model_real <- readRDS(
  file.path("application_results_real_data", 
            "seed_1_dummy_data_0intercept_only_1N_covs_1param_Xurand_thr_FALSECS_1het_sigma_0_applied_results.RDS"))

model_summary_and_trace_obj_REAL <- Model_real$model_summary_and_trace_obj






 


##
## ---- Extract NMA AUC metrics: ----------------------------------------------------------------------------------------
##
# tibble_NMA_AUC_metrics <- model_summary_and_trace_obj_REAL$extract_AUC(new_cov_data = NULL,
#                                                                        test_names = test_names,
#                                                                        use_BayesMVP_for_faster_summaries = TRUE)
# tibble_NMA_AUC_metrics
# ##
# tibble_NMA_AUC_metrics <- model_summary_and_trace_obj_DUMMY$extract_AUC(new_cov_data = NULL,
#                                                                         test_names = test_names,
#                                                                         use_BayesMVP_for_faster_summaries = TRUE)
# tibble_NMA_AUC_metrics
# # # tibble_NMA_AUC_metrics
# ##
# AUC_outs_REAL <- R_fn_extract_NMA_AUC_summary( tibble_gq = model_summary_and_trace_obj_REAL$get_summary_generated_quantities(),
#                                                n_thr = n_thr,
#                                                test_names = test_names
#                                                # network = network
#                                                )
# ##
# AUC_outs_REAL$auc
# AUC_outs_REAL$auc_diff
# AUC_outs_REAL$auc_pred
# # # A tibble: 4 × 7
# test test_name AUC_mean AUC_median AUC_lower AUC_upper AUC_sd
# <dbl> <chr>        <dbl>      <dbl>     <dbl>     <dbl>  <dbl>
# 1     1 GAD-2        0.834      0.834     0.800     0.865 0.0165
# 2     2 GAD-7        0.858      0.858     0.824     0.887 0.0159
# 3     3 HADS         0.846      0.846     0.813     0.876 0.0162
# 4     4 BAI          0.836      0.839     0.757     0.895 0.0352
# > AUC_outs_REAL$auc_diff
# # A tibble: 6 × 10
# test1 test2 test1_name test2_name AUC_diff_mean AUC_diff_median AUC_diff_lower AUC_diff_upper AUC_diff_sd prob_test1_better
# <int> <int> <chr>      <chr>              <dbl>           <dbl>          <dbl>          <dbl>       <dbl>             <dbl>
# 1     1     2 GAD-2      GAD-7           -0.0241         -0.0239         -0.0613         0.0126      0.0187                NA
# 2     1     3 GAD-2      HADS            -0.0119         -0.0121         -0.0548         0.0314      0.0218                NA
# 3     2     3 GAD-7      HADS            -0.00196        -0.00451        -0.0693         0.0801      0.0382                NA
# 4     1     4 GAD-2      BAI              0.0122          0.0122         -0.0305         0.0538      0.0214                NA
# 5     2     4 GAD-7      BAI              0.0221          0.0199         -0.0450         0.104       0.0380                NA
# 6     3     4 HADS       BAI              0.00995         0.00697        -0.0553         0.0912      0.0376                NA
# AUC_outs_REAL$auc_pred
# # A tibble: 4 × 7
# test test_name AUC_mean AUC_median AUC_lower AUC_upper AUC_sd
# <dbl> <chr>        <dbl>      <dbl>     <dbl>     <dbl>  <dbl>
# 1     1 GAD-2        0.808      0.823     0.588     0.946 0.0929
# 2     2 GAD-7        0.830      0.848     0.586     0.966 0.0989
# 3     3 HADS         0.826      0.840     0.607     0.954 0.0892
# 4     4 BAI          0.817      0.835     0.574     0.955 0.0986
##
# ## ---- Make LaTeX table:
# ##
# AUC_latex_table_REAL <- create_AUC_latex_table( auc_outs = AUC_outs_REAL, 
#                                                 digits = 3)
# cat(AUC_latex_table_REAL)


str(model_summary_and_trace_obj_DUMMY$get_trace_generated_quantities())
##
## ---- For dummy applied data:
##
tibble_gq = model_summary_and_trace_obj_DUMMY$get_summary_generated_quantities()

AUC_outs_DUMMY <- R_fn_extract_NMA_AUC_summary(  tibble_gq = tibble_gq,
                                                   n_thr = n_thr,
                                                 test_names = test_names,
                                                 network = network
                                                 )
##
AUC_outs_DUMMY$auc
AUC_outs_DUMMY$auc_diff
AUC_outs_DUMMY$auc_pred
 
# # A tibble: 4 × 7
# test test_name AUC_mean AUC_median AUC_lower AUC_upper AUC_sd
# <dbl> <chr>        <dbl>      <dbl>     <dbl>     <dbl>  <dbl>
#   1     1 GAD-2        0.844      0.845     0.808     0.878 0.0181
# 2     2 GAD-7        0.893      0.893     0.868     0.914 0.0118
# 3     3 HADS         0.855      0.855     0.824     0.882 0.0148
# 4     4 BAI          0.867      0.868     0.821     0.905 0.0217
# > AUC_outs_DUMMY$auc_diff
# # A tibble: 6 × 10
# test1 test2 test1_name test2_name AUC_diff_mean AUC_diff_median AUC_diff_lower AUC_diff_upper AUC_diff_sd prob_test1_better
# <int> <int> <chr>      <chr>              <dbl>           <dbl>          <dbl>          <dbl>       <dbl>             <dbl>
#   1     1     2 GAD-2      GAD-7            -0.0482         -0.0479       -0.0856         -0.0138      0.0186                NA
# 2     1     3 GAD-2      HADS             -0.0105         -0.0102       -0.0534          0.0296      0.0212                NA
# 3     2     3 GAD-7      HADS             -0.0225         -0.0231       -0.0768          0.0334      0.0282                NA
# 4     1     4 GAD-2      BAI               0.0377          0.0372        0.00627         0.0701      0.0165                NA
# 5     2     4 GAD-7      BAI               0.0257          0.0252       -0.0202          0.0746      0.0244                NA
# 6     3     4 HADS       BAI              -0.0120         -0.0131       -0.0618          0.0404      0.0262                NA
# > AUC_outs_DUMMY$auc_pred
# # A tibble: 4 × 7
# test test_name AUC_mean AUC_median AUC_lower AUC_upper AUC_sd
# <dbl> <chr>        <dbl>      <dbl>     <dbl>     <dbl>  <dbl>
#   1     1 GAD-2        0.823      0.840     0.592     0.961 0.0974
# 2     2 GAD-7        0.870      0.890     0.654     0.979 0.0853
# 3     3 HADS         0.832      0.852     0.583     0.970 0.102 
# 4     4 BAI          0.846      0.867     0.604     0.972 0.0961
##
## ---- Make LaTeX table:
##
AUC_latex_table_DUMMY <- create_AUC_latex_table( auc_outs = AUC_outs_DUMMY, 
                                                 digits = 3)
cat(AUC_latex_table_DUMMY)
# \begin{table}[!h]
# \centering
# \caption{AUC values, pairwise differences, and predictive intervals}
# \small
# \begin{tabular}{lcc}
# \toprule
# \multicolumn{3}{l}{\textbf{AUC estimates (95\% CI):}} \\
# \midrule
# GAD-2 & 0.845 & (0.808, 0.878) \\
# GAD-7 & 0.893 & (0.868, 0.914) \\
# HADS & 0.855 & (0.824, 0.882) \\
# BAI & 0.868 & (0.821, 0.905) \\
# \midrule
# \multicolumn{3}{l}{\textbf{Pairwise AUC differences:}} \\
# \midrule
# GAD-2 vs GAD-7 & -0.048 & (-0.086, -0.014)* \\
# GAD-2 vs HADS & -0.010 & (-0.053, 0.030) \\
# GAD-7 vs HADS & -0.023 & (-0.077, 0.033) \\
# GAD-2 vs BAI & 0.037 & (0.006, 0.070)* \\
# GAD-7 vs BAI & 0.025 & (-0.020, 0.075) \\
# HADS vs BAI & -0.013 & (-0.062, 0.040) \\
# \midrule
# \multicolumn{3}{l}{\textbf{AUC predictive intervals:}} \\
# \midrule
# GAD-2 & 0.840 & (0.592, 0.961) \\
# GAD-7 & 0.890 & (0.654, 0.979) \\
# HADS & 0.852 & (0.583, 0.970) \\
# BAI & 0.867 & (0.604, 0.972) \\
# \bottomrule
# \end{tabular}
# \begin{tablenotes}
# \footnotesize
# \item * indicates 95\% CI excludes zero (significant difference)
# \item Predictive intervals incorporate between-study heterogeneity
# \end{tablenotes}
# \end{table}

##
## ----------- Plot sROC curve(s) (NMA - for best model WITHOUT covariates): ------------------------------------------
##
relevant_thresholds = list( "GAD-2" = c(1:6), 
                            "GAD-7" =  c(3:18), 
                            "HADS" =  c(3:18),
                            "BAI" = c(3:38))
# stan_model_file_name = model_summary_and_trace_obj_DUMMY$internal_obj$outs_stan_model_name$stan_model_file_name
# stan_mod_samples     = model_summary_and_trace_obj_DUMMY$internal_obj$outs_stan_sampling$stan_mod_samples
# ##
# test_names = test_names
# # covariates = covariates,
# # new_cov_data = new_cov_data,
# ##
# # df_true = df_true
# ##
# conf_region_colour = "blue"
# pred_region_colour = "green"
# ##
# base_size = 24
# point_size = 4.0
# linewidth = 1.0
# ##
# n_index_tests = n_index_tests
# n_thr = n_thr
#
#
##
# plots_REAL <- model_summary_and_trace_obj_REAL$plot_sROC(test_names = test_names)
##
# plots_REAL <- R_fn_sROC_plot_NMA(     stan_model_file_name = model_summary_and_trace_obj_REAL$internal_obj$outs_stan_model_name$stan_model_file_name,
#                                       stan_mod_samples     =  model_summary_and_trace_obj_REAL$internal_obj$outs_stan_sampling$stan_mod_samples,
#                                       ##
#                                       test_names = test_names,
#                                       ##
#                                       relevant_thresholds = relevant_thresholds,
#                                       ##
#                                       # covariates = covariates,
#                                       # new_cov_data = new_cov_data,
#                                       ##
#                                       df_true = NULL,
#                                       ##
#                                       conf_region_colour = "blue", 
#                                       pred_region_colour = "green",
#                                       ##
#                                       conf_region_alpha = 0.30,
#                                       pred_region_alpha = 0.15,
#                                       ##
#                                       base_size = 24,
#                                       point_size = 4.0,
#                                       linewidth = 1.0,
#                                       ##
#                                       n_index_tests = n_index_tests,
#                                       n_thr = n_thr)
# ##
# plots_REAL$plot_list[[1]]
# #### plots_REAL$plot_list[[2]]
# #### plots_REAL$plot_list[[3]]
# #### plots_REAL$plot_list[[4]]
# plots_REAL$plot_list[[5]]
# #### plots_REAL$plot_list[[6]]
# plots_REAL$plot_list[[7]]
# plots_REAL$plot_list[[8]]
# plots_REAL$plot_list[[9]]
# ##
# ggsave("REAL_Application_baseline_sROC.png", plots_REAL$plot_list[[1]], width = 16, height = 9, dpi = 400)
# ggsave("REAL_Application_baseline_sROC_panel_w_CrI_PrI.png", plots_REAL$plot_list[[6]], width = 16, height = 9, dpi = 400)
# ##
# ggsave("REAL_Application_baseline_Se_Sp_vs_thr.png", plots_REAL$plot_list[[7]], width = 16, height = 16, dpi = 400)
# ggsave("REAL_Application_baseline_Se_Sp_vs_thr_w_CrI.png", plots_REAL$plot_list[[8]], width = 16, height = 16, dpi = 400)
# ggsave("REAL_Application_baseline_Se_Sp_vs_thr_panel_w_CrI_PrI.png", plots_REAL$plot_list[[9]], width = 16, height = 16, dpi = 400)
# ##
##
##
# plots_DUMMY <- model_summary_and_trace_obj_DUMMY$plot_sROC(test_names = test_names)
##
n_thr <- c(6, 21, 21, 63)
##
plots_DUMMY <- R_fn_sROC_plot_NMA(     stan_model_file_name = model_summary_and_trace_obj_DUMMY$internal_obj$outs_stan_model_name$stan_model_file_name,
                                       stan_mod_samples     = model_summary_and_trace_obj_DUMMY$internal_obj$outs_stan_sampling$stan_mod_samples,
                                       ##
                                       test_names = test_names,
                                       ##
                                       relevant_thresholds = relevant_thresholds,
                                       ##
                                       # covariates = covariates,
                                       # new_cov_data = new_cov_data,
                                       ##
                                       df_true = NULL,
                                       ##
                                       conf_region_colour = "blue", 
                                       pred_region_colour = "green",
                                       ##
                                       conf_region_alpha = 0.20,
                                       pred_region_alpha = 0.10,
                                       ##
                                       base_size = 24,
                                       point_size = 4.0,
                                       linewidth = 1.0,
                                       ##
                                       n_index_tests = n_index_tests,
                                       n_thr = n_thr)
##
plots_DUMMY$plot_list[[1]]
#### plots_DUMMY$plot_list[[2]]
#### plots_DUMMY$plot_list[[3]]
#### plots_DUMMY$plot_list[[4]]
plots_DUMMY$plot_list[[5]]
#### plots_DUMMY$plot_list[[6]]
plots_DUMMY$plot_list[[7]]
plots_DUMMY$plot_list[[8]]
plots_DUMMY$plot_list[[9]]
##
ggsave("DUMMY_Application_baseline_sROC.png", plots_DUMMY$plot_list[[1]], width = 16, height = 9, dpi = 400)
ggsave("DUMMY_Application_baseline_sROC_panel_w_CrI_PrI.png", plots_DUMMY$plot_list[[5]], width = 16, height = 16, dpi = 400)
##
ggsave("DUMMY_Application_baseline_Se_Sp_vs_thr.png", plots_DUMMY$plot_list[[7]], width = 16, height = 16, dpi = 400)
ggsave("DUMMY_Application_baseline_Se_Sp_vs_thr_w_CrI.png", plots_DUMMY$plot_list[[8]], width = 16, height = 16, dpi = 400)
ggsave("DUMMY_Application_baseline_Se_Sp_vs_thr_panel_w_CrI_PrI.png", plots_DUMMY$plot_list[[9]], width = 16, height = 16, dpi = 400)
# ##
# ## ---- Compare real to simulated/dummy data (for personal interest NOT FOR PAPER/THESIS):
# ##
# require(patchwork)
# plots_DUMMY$plot_list[[1]] + plots_REAL$plot_list[[1]]




##
## ---- Extract NMA Se/Sp comparison metrics for tables + plots  ------------------------------------------------------------------------
##




draws_array <- model_summary_and_trace_obj_DUMMY$get_all_traces_list()$traces_as_arrays$draws_array
# trace_gq <-  model_summary_and_trace_obj_DUMMY$get_all_traces_list()$traces_as_arrays$trace_gq
# str(trace_gq)
##
outs_NMA_comp_DUMMY <-  R_fn_using_Rcpp_compute_NMA_comparisons_posthoc(  trace_gq = draws_array,
                                                                          n_thr = n_thr,
                                                                          test_names = test_names)
##






# Create separate plots for each clinical context

# Define clinical threshold groups
get_threshold_groups <- function() {
  
  list(
  
    screening = list(
      "GAD-2" = c(2, 3),         # 3 is SCREENING standard for GAD
      "GAD-7" = c(9, 10),    # 10 is SCREENING standard  
      # "HADS" = c(7, 8, 9),       # 8 for SCREENING ("possible") (11 is "probable")
      "HADS" = c(9, 10, 11),       # 8 for SCREENING ("possible") (11 is "probable")
      # "BAI" = c(7, 8, 9)         # 8 for mild (screening??)
        "BAI" = c(17, 18, 19)      # 16 is standard
    ),
    
    ## More thresholds here, use this for more detailed (e.g. appendix for most/all
    ## of them) pairwise plots (choose(4, 2) = 6 plots):
    screening_2 = list(
      "GAD-2" = c(2, 3, 4),      # 3 is SCREENING standard for GAD
      "GAD-7" = c(8, 9, 10, 11),    # 10 is SCREENING standard for GAD
      "HADS" = c(8, 9, 10, 11, 12),    # 11 is standard
      "BAI" = c(16, 17, 18, 19, 20)      # 16 is standard
    )
 
  )
  
}

##
threshold_groups <- get_threshold_groups()
threshold_groups
##
threshold_groups$screening
##
screening_comps <- filter_clinically_relevant_comparisons( tibble_comp = outs_NMA_comp_DUMMY,
                                                           relevant_thresholds = threshold_groups$screening) ## %>%
                   # identify_important_comparisons(n_per_type = 100)
screening_comps %>% print(n = 1000)



##
## ---- Create forest plots for Se + Sp - screening thresholds for GAD:
##
plot_dpi <- 100
##
plot_screening_Se <- plot_se_differences( filtered_tibble = screening_comps,
                                          show_top_n = 1000,
                                          order_by = "default",
                                          plot_title = "Sensitivity (Se) Differences: Screening Thresholds")
##
plot_screening_Sp <- plot_sp_differences( filtered_tibble = screening_comps,
                                          show_top_n = 1000,
                                          order_by = "default",
                                          plot_title = "Specificity (Sp) Differences: Screening Thresholds")
## 
plot_screening_Se_and_Sp <- plot_screening_Se + plot_screening_Sp
plot_screening_Se_and_Sp
##
ggsave("DUMMY_Application_Se_and_Sp_pairwise_diffs_screening_GAD.png", 
       plot_screening_Se_and_Sp, 
       width = 16, height = 16, dpi = plot_dpi)


##
## ---- Se/Sp pairwise diff's subsetted by test-pair combos's (choose(4, 2) = 6 plots in total):
##
test_names
##
screening_comps_2 <- filter_clinically_relevant_comparisons( tibble_comp = outs_NMA_comp_DUMMY,
                                                           relevant_thresholds = threshold_groups$screening_2) ## %>%




##
## ---- GAD-2 vs. GAD-7:
##
test_pair <- c("GAD-2", "GAD-7")
screening_comps_GAD_2_vs_GAD_7 <- dplyr::filter(screening_comps_2, (test1_name %in% test_pair) & (test2_name %in% test_pair))
##
plot_screening_Se <- plot_se_differences( filtered_tibble = screening_comps_GAD_2_vs_GAD_7,
                                          show_top_n = 1000,
                                          order_by = "default",
                                          plot_title = "Sensitivity (Se) Differences: Screening Thresholds")
##
plot_screening_Sp <- plot_sp_differences( filtered_tibble = screening_comps_GAD_2_vs_GAD_7,
                                          show_top_n = 1000,
                                          order_by = "default",
                                          plot_title = "Specificity (Sp) Differences: Screening Thresholds")
## 
plot_screening_Se_and_Sp_GAD_2_vs_GAD_7 <- plot_screening_Se + plot_screening_Sp
plot_screening_Se_and_Sp_GAD_2_vs_GAD_7
##
ggsave("DUMMY_Application_Se_and_Sp_pairwise_diffs_screening_GAD_GAD_2_vs_GAD_7.png", 
       plot_screening_Se_and_Sp_GAD_2_vs_GAD_7, 
       width = 16, height = 16, dpi = plot_dpi)



##
## ---- GAD-2 vs. HADS
##
test_pair <- c("GAD-2", "HADS")
screening_comps_GAD_2_vs_HADS <- dplyr::filter(screening_comps_2, (test1_name %in% test_pair) & (test2_name %in% test_pair))
##
plot_screening_Se <- plot_se_differences( filtered_tibble = screening_comps_GAD_2_vs_HADS,
                                          show_top_n = 1000,
                                          order_by = "default",
                                          plot_title = "Sensitivity (Se) Differences: Screening Thresholds")
##
plot_screening_Sp <- plot_sp_differences( filtered_tibble = screening_comps_GAD_2_vs_HADS,
                                          show_top_n = 1000,
                                          order_by = "default",
                                          plot_title = "Specificity (Sp) Differences: Screening Thresholds")
## 
plot_screening_Se_and_Sp_GAD_2_vs_HADS <- plot_screening_Se + plot_screening_Sp
plot_screening_Se_and_Sp_GAD_2_vs_HADS
##
ggsave("DUMMY_Application_Se_and_Sp_pairwise_diffs_screening_GAD_GAD_2_vs_HADS.png", 
       plot_screening_Se_and_Sp_GAD_2_vs_HADS, 
       width = 16, height = 16, dpi = plot_dpi)



##
## ---- GAD-7 vs. HADS:
##
test_pair <- c("GAD-7", "HADS")
screening_comps_GAD_7_vs_HADS <- dplyr::filter(screening_comps_2, (test1_name %in% test_pair) & (test2_name %in% test_pair))
##
plot_screening_Se <- plot_se_differences( filtered_tibble = screening_comps_GAD_7_vs_HADS,
                                          show_top_n = 1000,
                                          order_by = "default",
                                          plot_title = "Sensitivity (Se) Differences: Screening Thresholds")
##
plot_screening_Sp <- plot_sp_differences( filtered_tibble = screening_comps_GAD_7_vs_HADS,
                                          show_top_n = 1000,
                                          order_by = "default",
                                          plot_title = "Specificity (Sp) Differences: Screening Thresholds")
## 
plot_screening_Se_and_Sp_GAD_7_vs_HADS <- plot_screening_Se + plot_screening_Sp
plot_screening_Se_and_Sp_GAD_7_vs_HADS
##
ggsave("DUMMY_Application_Se_and_Sp_pairwise_diffs_screening_GAD_GAD_7_vs_HADS.png", 
       plot_screening_Se_and_Sp_GAD_7_vs_HADS, 
       width = 16, height = 16, dpi = plot_dpi)



##
## ---- GAD-2 vs. BAI:
##
test_pair <- c("GAD-2", "BAI")
screening_comps_GAD_2_vs_BAI <- dplyr::filter(screening_comps_2, (test1_name %in% test_pair) & (test2_name %in% test_pair))
##
plot_screening_Se <- plot_se_differences( filtered_tibble = screening_comps_GAD_2_vs_BAI,
                                          show_top_n = 1000,
                                          order_by = "default",
                                          plot_title = "Sensitivity (Se) Differences: Screening Thresholds")
##
plot_screening_Sp <- plot_sp_differences( filtered_tibble = screening_comps_GAD_2_vs_BAI,
                                          show_top_n = 1000,
                                          order_by = "default",
                                          plot_title = "Specificity (Sp) Differences: Screening Thresholds")
## 
plot_screening_Se_and_Sp_GAD_2_vs_BAI <- plot_screening_Se + plot_screening_Sp
plot_screening_Se_and_Sp_GAD_2_vs_BAI
##
ggsave("DUMMY_Application_Se_and_Sp_pairwise_diffs_screening_GAD_GAD_2_vs_BAI.png", 
       plot_screening_Se_and_Sp_GAD_2_vs_BAI, 
       width = 16, height = 16, dpi = plot_dpi)



##
## ---- GAD-7 vs. BAI:
##
test_pair <- c("GAD-7", "BAI")
screening_comps_GAD_7_vs_BAI <- dplyr::filter(screening_comps_2, (test1_name %in% test_pair) & (test2_name %in% test_pair))
##
plot_screening_Se <- plot_se_differences( filtered_tibble = screening_comps_GAD_7_vs_BAI,
                                          show_top_n = 1000,
                                          order_by = "default",
                                          plot_title = "Sensitivity (Se) Differences: Screening Thresholds")
##
plot_screening_Sp <- plot_sp_differences( filtered_tibble = screening_comps_GAD_7_vs_BAI,
                                          show_top_n = 1000,
                                          order_by = "default",
                                          plot_title = "Specificity (Sp) Differences: Screening Thresholds")
## 
plot_screening_Se_and_Sp_GAD_7_vs_BAI <- plot_screening_Se + plot_screening_Sp
plot_screening_Se_and_Sp_GAD_7_vs_BAI
##
ggsave("DUMMY_Application_Se_and_Sp_pairwise_diffs_screening_GAD_GAD_7_vs_BAI.png", 
       plot_screening_Se_and_Sp_GAD_7_vs_BAI, 
       width = 16, height = 16, dpi = plot_dpi)



##
## ---- HADS vs. BAI:
##
test_pair <- c("HADS", "BAI")
screening_comps_HADS_vs_BAI <- dplyr::filter(screening_comps_2, (test1_name %in% test_pair) & (test2_name %in% test_pair))
##
plot_screening_Se <- plot_se_differences( filtered_tibble = screening_comps_HADS_vs_BAI,
                                          show_top_n = 1000,
                                          order_by = "default",
                                          plot_title = "Sensitivity (Se) Differences: Screening Thresholds")
##
plot_screening_Sp <- plot_sp_differences( filtered_tibble = screening_comps_HADS_vs_BAI,
                                          show_top_n = 1000,
                                          order_by = "default",
                                          plot_title = "Specificity (Sp) Differences: Screening Thresholds")
## 
plot_screening_Se_and_Sp_HADS_vs_BAI <- plot_screening_Se + plot_screening_Sp
plot_screening_Se_and_Sp_HADS_vs_BAI
##
ggsave("DUMMY_Application_Se_and_Sp_pairwise_diffs_screening_GAD_HADS_vs_BAI.png", 
       plot_screening_Se_and_Sp_HADS_vs_BAI, 
       width = 16, height = 16, dpi = plot_dpi)






# # Check what comparisons exist in the original data
# outs_NMA_comp_DUMMY %>%
#   filter(
#     test1_name %in% c("GAD-2", "GAD-7", "HADS", "BAI") &
#       test2_name %in% c("GAD-2", "GAD-7", "HADS", "BAI")
#   ) %>%
#   count(test1_name, test2_name) %>%
#   print(n = 20)
# 
# # See which specific threshold pairs exist
# screening_comps %>%
#   count(test1_name, threshold1, test2_name, threshold2) %>%
#   print(n = 100)
# 


# ##
# ##
# ## ---- Extract NMA performance metrics  -----------------------------------------------------------------------------
# ## NOTE: (not doing these for dummy analysis - only do for real if Klaus et al insist):
# ##
# # tibble_NMA_performance_metrics <- model_summary_and_trace_obj_REAL$extract_NMA_performance()
# # tibble_NMA_performance_metrics
# # ##
# # tibble_NMA_performance_metrics <- model_summary_and_trace_obj_DUMMY$extract_NMA_performance()
# # tibble_NMA_performance_metrics
# ##
# outs <- R_fn_using_Rcpp_compute_NMA_performance_posthoc(  trace_gq = model_summary_and_trace_obj_REAL$get_trace_generated_quantities(),
#                                                           test_names = test_names,
#                                                           n_thr = n_thr)
# outs
# ##
# ##
# outs <- R_fn_using_Rcpp_compute_NMA_performance_posthoc(  trace_gq = model_summary_and_trace_obj_DUMMY$get_trace_generated_quantities(),
#                                                           test_names = test_names,
#                                                           n_thr = n_thr)
# outs
















##
## ----------- Meta-regression results tables + plots: --------------------------------------------------------------------------------------------------------------------------
##









# 
#  
##
all_MR_scenarios <- list(
    ##
    ## ---- Ref test = MINI:
    ##
    "Low GAD prevalence (~5%),   Ref. test = MINI"   = list( logit_prev_GAD = -1.0, 
                                                             Ref_test_clean_SCID = 0, 
                                                             Ref_test_clean_Structured = 0),
    ##
    "Median GAD prevalence (~10%), Ref. test = MINI" = list( logit_prev_GAD =  0.0, 
                                                             Ref_test_clean_SCID = 0, 
                                                             Ref_test_clean_Structured = 0),
    ##
    "High GAD prevalence (~22%), Ref. test = MINI"   = list( logit_prev_GAD = +1.0,
                                                             Ref_test_clean_SCID = 0, 
                                                             Ref_test_clean_Structured = 0),
    ##
    ## ---- Ref test = SCID:
    ##
    "Low GAD prevalence (~5%),   Ref. test = SCID"   = list( logit_prev_GAD = -1.0, 
                                                             Ref_test_clean_SCID = 1, 
                                                             Ref_test_clean_Structured = 0),
    ##
    "Median GAD prevalence (~10%), Ref. test = SCID" = list( logit_prev_GAD =  0.0, 
                                                             Ref_test_clean_SCID = 1, 
                                                             Ref_test_clean_Structured = 0),
    ##
    "High GAD prevalence (~22%), Ref. test = SCID"   = list( logit_prev_GAD = +1.0,
                                                             Ref_test_clean_SCID = 1, 
                                                             Ref_test_clean_Structured = 0),
    ##
    ## ---- Ref test = Structured:
    ##
    "Low GAD prevalence (~5%),   Ref. test = Structured"   = list( logit_prev_GAD = -1.0, 
                                                                   Ref_test_clean_SCID = 0, 
                                                                   Ref_test_clean_Structured = 1),
    ##
    "Median GAD prevalence (~10%), Ref. test = Structured" = list( logit_prev_GAD =  0.0, 
                                                                   Ref_test_clean_SCID = 0, 
                                                                   Ref_test_clean_Structured = 1),
    ##
    "High GAD prevalence (~22%), Ref. test = Structured"   = list( logit_prev_GAD = +1.0,
                                                                   Ref_test_clean_SCID = 0, 
                                                                   Ref_test_clean_Structured = 1)
)
##
baseline_case_nd_list <- baseline_case_d_list <- list()
##
for (k in 1:length(all_MR_scenarios)) {
  
  baseline_case_d <- baseline_case_nd <- list()
  ##
  baseline_case_d[["GAD-2"]] <- c( 1, 
                                   all_MR_scenarios[[k]]$logit_prev_GAD,  ## mean prev_GAD
                                   all_MR_scenarios[[k]]$Ref_test_clean_SCID,
                                   all_MR_scenarios[[k]]$Ref_test_clean_Structured) 
  baseline_case_nd[["GAD-2"]] <- baseline_case_d[["GAD-2"]]
  ##
  baseline_case_d[["GAD-7"]] <- c( 1, 
                                   all_MR_scenarios[[k]]$logit_prev_GAD,  ## mean prev_GAD
                                   all_MR_scenarios[[k]]$Ref_test_clean_SCID,
                                   all_MR_scenarios[[k]]$Ref_test_clean_Structured) 
  baseline_case_nd[["GAD-7"]] <- baseline_case_d[["GAD-7"]]
  ##
  baseline_case_d[["HADS"]] <- c( 1,
                                  all_MR_scenarios[[k]]$logit_prev_GAD,  ## mean prev_GAD
                                  all_MR_scenarios[[k]]$Ref_test_clean_SCID,
                                  all_MR_scenarios[[k]]$Ref_test_clean_Structured) 
  baseline_case_nd[["HADS"]] <- baseline_case_d[["HADS"]]
  ##
  baseline_case_d[["BAI"]] <- c( 1,
                                 all_MR_scenarios[[k]]$logit_prev_GAD,  ## mean prev_GAD
                                 all_MR_scenarios[[k]]$Ref_test_clean_SCID,
                                 all_MR_scenarios[[k]]$Ref_test_clean_Structured) 
  baseline_case_nd[["BAI"]] <- baseline_case_d[["BAI"]]
  ##
  
  baseline_case_nd_list[[k]] <-  baseline_case_nd
  baseline_case_d_list[[k]]  <- baseline_case_d
}

str(baseline_case_d_list)





# Example usage:
# First compute results for all scenarios
all_scenario_results <- R_fn_compute_MR_multiple_scenarios(
  model_prep_obj = MR_Model_dummy$model_prep_obj,
  model_summary_and_trace_obj = MR_Model_dummy$model_summary_and_trace_obj,
  baseline_case_nd_list = baseline_case_nd_list,
  baseline_case_d_list = baseline_case_d_list,
  use_probit_link = TRUE,
  scenario_names = names(all_MR_scenarios),
  test_names = test_names,
  compute_comparisons = TRUE,
  save_results = TRUE,
  output_dir = "MR_results"
)


# str(all_scenario_results)
# 
# all_scenario_results$`Low GAD prevalence (~5%),   Ref. test = MINI`$comparison_results
# 
# 
# # Create sROC plots
# sroc_plots <- R_fn_plot_MR_scenarios_sROC(
#   scenario_results = all_scenario_results,
#   test_names = test_names,
#   plot_type = c("tests_by_scenario", "scenarios_by_test", "individual"),
#   save_plots = TRUE,
#   output_dir = "MR_plots/sROC"
# )
# 
# sroc_plots
# 
# 
# # Create difference plots
# se_diff_plots <- R_fn_plot_MR_scenarios_differences(
#   scenario_results = all_scenario_results,
#   comparison_type = "Se",
#   show_top_n = 30,
#   save_plots = TRUE,
#   output_dir = "MR_plots/differences"
# )
# 
# 
# se_diff_plots
# 
# 
# 
# # Example usage:
# # Compute all comparisons
# all_comparisons <- R_fn_compute_MR_comparisons_all_scenarios(
#   scenario_results = all_scenario_results,
#   test_names = test_names
# )
# 
# 
# # Create forest plots for specific comparisons
# # You can specify which comparisons to show, or let the function choose
# forest_plot_se <- R_fn_plot_MR_forest_comparisons(
#   all_comparisons_data = all_comparisons,
#   selected_comparisons = c(
#     "GAD-2 (3) vs GAD-7 (8)",  # Commonly used thresholds
#     "GAD-2 (3) vs GAD-7 (9)",   # Short vs long form
#     "GAD-2 (3) vs GAD-7 (10)"       # Different scales
#   ),
#   comparison_type = "Se",
#   save_plots = TRUE
# )
# 
# forest_plot_se
# 
# forest_plot_se %>% facet_wrap(ncol = 1)




# Example usage that gives you full control:

# 1. Get the data
combined_data <- combine_sroc_data_scenarios(all_scenario_results)

# 2. Create basic plots
plots <- create_sroc_plots_flexible(
  scenario_results = all_scenario_results,
  test_names = test_names,
  layout = c("by_test", "by_scenario", "grid"),
  base_size = 14
)

plots

# # 3. Now you can modify any plot yourself:
# # Example: Take the GAD-2 plot and customize it
# gad2_plot <- plots[["by_test_GAD-2"]]
# 
# # Add your own customizations
# gad2_plot_custom <- gad2_plot +
#   scale_color_brewer(palette = "Set1") +
#   theme(legend.position = "bottom",
#         legend.text = element_text(size = 10),
#         plot.title = element_text(face = "bold")) +
#   guides(color = guide_legend(nrow = 3))
# 
# gad2_plot_custom
# 
# str(plots)
# 
# plots$`by_scenario_Low GAD prevalence (~5%),   Ref. test = MINI`



# Add confidence intervals to a specific scenario plot
scenario1_plot <- plots$`by_scenario_Low GAD prevalence (~5%),   Ref. test = MINI`

scenario1_plot




scenario1_with_ci <- add_regions_to_sroc(
  plot = scenario1_plot,
  data = combined_data %>% 
    filter(scenario_index == 1) ,
  show_ci = TRUE,
  show_pi = TRUE
)

scenario1_with_ci ## this doesnt work + not using your own ci/pi region fns.... 




##
## ---- sROC plots by test (2x2 panel if 4 tests)
##
# # Create your own custom grid
# library(patchwork)
# 
# # Customize the shared legend
# custom_grid <- (plots[["by_test_GAD-2"]] + plots[["by_test_GAD-7"]]) /
#   (plots[["by_test_HADS"]] + plots[["by_test_BAI"]]) +
#   plot_layout(guides = "collect") +
#   plot_annotation(
#     title = "Meta-regression sROC curves across all scenarios",
#     subtitle = "Each panel shows one test with 9 scenario curves"
#   ) & 
#   theme(legend.position = "bottom",
#         legend.box = "horizontal",
#         legend.title = element_text(size = 12, face = "bold"),
#         legend.text = element_text(size = 10)) &
#   guides(color = guide_legend(nrow = 3, byrow = TRUE))  # 3 rows of legend items
# 
# custom_grid
# custom_grid 
# 
# # Save any plot you want
# ggsave("my_custom_plot.pdf", custom_grid, width = 12, height = 10)
# 




##
## ---- sROC plots by test (2x2 panel if 4 tests)
##
# First, add the prevalence and ref test variables to your combined data
combined_data_enhanced <- combined_data %>%
  mutate(
    prev_level = case_when(
      grepl("Low.*prevalence", scenario) ~ "Low (5%)",
      grepl("Median.*prevalence", scenario) ~ "Medium (10%)",
      grepl("High.*prevalence", scenario) ~ "High (22%)",
      TRUE ~ "Unknown"
    ),
    ref_test = case_when(
      grepl("MINI", scenario) ~ "MINI",
      grepl("SCID", scenario) ~ "SCID",
      grepl("Structured", scenario) ~ "Structured",
      TRUE ~ "Unknown"
    )
  )


# Now recreate the plots with color AND linetype aesthetics
plots_enhanced <- list()

for (t in unique(combined_data_enhanced$test)) {
  
      plot_data <- combined_data_enhanced %>%
        filter(test == t)
      
      test_name <- ifelse(!is.null(test_names), test_names[t], paste0("Test ", t))
      
      p <- ggplot(plot_data, 
                  aes(x = Fp_median, y = Se_median, 
                      color = prev_level,      # Color by prevalence
                      linetype = ref_test,     # Linetype by ref test
                      group = scenario)) +     # Still group by scenario
        geom_line(size = 0.8) +
        geom_point(aes(shape = ref_test), size = 2) +  # Also use shape for ref test
        scale_color_manual(values = c("Low (5%)" = "#1f78b4",
                                      "Medium (10%)" = "#33a02c", 
                                      "High (22%)" = "#e31a1c"),
                           name = "GAD Prevalence") +
        scale_linetype_manual(values = c("MINI" = "solid", 
                                         "SCID" = "dashed", 
                                         "Structured" = "dotted"),
                              name = "Reference Test") +
        scale_shape_manual(values = c("MINI" = 16, 
                                      "SCID" = 17, 
                                      "Structured" = 15),
                           name = "Reference Test") +
        scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                           limits = c(0, 1)) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                           limits = c(0, 1)) +
        labs(x = "False positive rate",
             y = "Sensitivity",
             title = test_name) +
        theme_bw(base_size = 14)
      
      plots_enhanced[[paste0("by_test_", test_name)]] <- p
  
}


# Create grid with the enhanced plots
custom_grid <- (plots_enhanced[["by_test_GAD-2"]] + plots_enhanced[["by_test_GAD-7"]]) /
  (plots_enhanced[["by_test_HADS"]] + plots_enhanced[["by_test_BAI"]]) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Meta-regression sROC curves across all scenarios",
    subtitle = "Colored by prevalence, line type by reference test"
  ) & 
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10))

custom_grid






##
## ---- sROC plots by scenario (3x3 panel if 9 scenarios)
##
# Create your own custom grid
library(patchwork)

# Customize the shared legend
custom_grid <- ( (plots$`by_scenario_Low GAD prevalence (~5%),   Ref. test = MINI` +
                  plots$`by_scenario_Median GAD prevalence (~10%), Ref. test = MINI` +  
                  plots$`by_scenario_High GAD prevalence (~22%), Ref. test = MINI`) /
                  (plots$`by_scenario_Low GAD prevalence (~5%),   Ref. test = SCID`+ 
                  plots$`by_scenario_Median GAD prevalence (~10%), Ref. test = SCID` +
                  plots$`by_scenario_High GAD prevalence (~22%), Ref. test = SCID`) /
                  (plots$`by_scenario_Low GAD prevalence (~5%),   Ref. test = Structured` + 
                  plots$`by_scenario_Median GAD prevalence (~10%), Ref. test = Structured` + 
                  plots$`by_scenario_High GAD prevalence (~22%), Ref. test = Structured`)) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Meta-regression sROC curves across all scenarios",
    subtitle = "Each panel shows one test with 9 scenario curves"
  ) & 
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10)) &
  guides(color = guide_legend(nrow = 3, byrow = TRUE))  # 3 rows of legend items

custom_grid
custom_grid ## need to make it so all 4 panels share just 1 legend as legend is identical always!!!!

# Save any plot you want
ggsave("my_custom_plot.pdf", custom_grid, width = 12, height = 10)








# Example usage for your specific case:
# First, add the group variables to your data if they don't exist
all_comparisons_enhanced <- all_comparisons %>%
  mutate(
    # Parse scenario names - this part is specific to your naming convention
    prev_level = case_when(
      grepl("Low.*prevalence", scenario) ~ "Low (5%)",
      grepl("Median.*prevalence", scenario) ~ "Medium (10%)",
      grepl("High.*prevalence", scenario) ~ "High (22%)",
      TRUE ~ "Unknown"
    ),
    ref_test = case_when(
      grepl("MINI", scenario) ~ "MINI",
      grepl("SCID", scenario) ~ "SCID",
      grepl("Structured", scenario) ~ "Structured",
      TRUE ~ "Unknown"
    )
  )

# Now use the general function
forest_plot <- R_fn_plot_forest_comparisons_general(
  all_comparisons_data = all_comparisons_enhanced,
  comparison_type = "Se",
  scenario_var = "scenario",
  group_var = "prev_level",    # Color by prevalence
  shape_var = "ref_test",       # Shape by reference test
  n_auto_select = 6
)

# Customize colors and shapes after the fact
forest_plot_custom <- forest_plot +
  scale_color_manual(values = c("Low (5%)" = "#1f78b4",
                                "Medium (10%)" = "#33a02c", 
                                "High (22%)" = "#e31a1c"),
                     name = "GAD Prevalence") +
  scale_shape_manual(values = c("MINI" = 16, 
                                "SCID" = 17, 
                                "Structured" = 15),
                     name = "Reference Test")


forest_plot_custom





# Get the base plots
threshold_plots <- create_threshold_plots_flexible(
  scenario_results = all_scenario_results,
  test_names = c("GAD-2", "GAD-7", "HADS", "BAI"),
  plot_type = c("Se", "Sp"),
  layout = c("by_test", "by_scenario")
)


# Get the combined data for customization
combined_data <- extract_threshold_data_scenarios(all_scenario_results)


# Add your specific variables AFTER
combined_data_enhanced <- combined_data %>%
  mutate(
    prev_level = case_when(
      grepl("Low.*prevalence", scenario) ~ "Low (5%)",
      grepl("Median.*prevalence", scenario) ~ "Medium (10%)",
      grepl("High.*prevalence", scenario) ~ "High (22%)"
    ),
    ref_test = case_when(
      grepl("MINI", scenario) ~ "MINI",
      grepl("SCID", scenario) ~ "SCID",  
      grepl("Structured", scenario) ~ "Structured"
    ),
    test_name = c("GAD-2", "GAD-7", "HADS", "BAI")[test]
  )



# Now customize a specific plot
custom_se_gad2 <- ggplot(combined_data_enhanced %>% filter(test == 1),
                         aes(x = threshold_norm, y = Se_median,
                             color = prev_level, 
                             linetype = ref_test,
                             group = scenario)) +
  geom_line(size = 0.8) +
  geom_point(aes(shape = ref_test), size = 2) +
  scale_color_manual(values = c("Low (5%)" = "#1f78b4",
                                "Medium (10%)" = "#33a02c", 
                                "High (22%)" = "#e31a1c"),
                     name = "GAD Prevalence") +
  scale_linetype_manual(values = c("MINI" = "solid", 
                                   "SCID" = "dashed", 
                                   "Structured" = "dotted"),
                        name = "Reference Test") +
  scale_shape_manual(values = c("MINI" = 16, "SCID" = 17, "Structured" = 15),
                     name = "Reference Test") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  labs(x = "Normalized threshold", 
       y = "Sensitivity",
       title = "GAD-2 Sensitivity across Meta-regression Scenarios") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom")


custom_se_gad2





# # In your paper:
# "All plots were created using the MetaOrdDTA R package. Base plots were 
# generated using the following R6 class methods:
# 
# - sROC curves: `model$plot$sroc()`
# - Se/Sp vs threshold: `model$plot$threshold()`  
# - Forest plots: `model$plot$forest()`
# - Pairwise comparisons: `model$plot$comparisons()`
# 
# Custom layouts were created using patchwork. Full code for all figures 
# is available at: https://github.com/yourname/MetaOrdDTA-paper"






# # Create heatmap showing all effects
# heatmap_se <- R_fn_plot_MR_comparison_heatmap(
#   all_comparisons_data = all_comparisons,
#   comparison_type = "Se",
#   n_comparisons = 30,
#   save_plot = TRUE
# )
# 
# heatmap_se





# 
 
# 
# ## If WTE reduction kicks in 1st September 2025:
# 4*(0.3/0.2) ## = 6, so get extra 2 months (so ends 28th feb rather than 31st december)
# ##
# ## If WTE reduction kicks in 1st October 2025:
# 3*(0.3/0.2) ## = 4.5 ## so get extra 14th feb
# 
# 

# 
# 
# MR_Model_dummy$model_prep_obj$internal_obj$outs_stan_compile$stan_model_file_path
# 
# 
# str(stan_mod_samples)
# 
# 
# # Method 3: Check summary
# summary_df <- stan_mod_samples$summary()
# print(head(summary_df))
# print(unique(summary_df$variable)[1:20])  # First 20 unique variable names
# 
# # Method 4: Try to see the actual error from failed chain
# # This is important - it will tell us exactly what went wrong
# gq_model$generate_quantities(
#   fitted_params = fitted_draws[1:10, ],  # Just 10 iterations
#   data = gq_data,
#   seed = 123
#   # chains = 1,
#   # show_messages = TRUE  # This might show more info
# )
# 
# # Then immediately check the output
# gq_model$output(1)
# 
# 
# 
main_trace =trace_main_as_list
nuisance_trace = trace_nuisance_as_list
str(nuisance_trace)
n_nuisance <- 1
include_nuisance

summarise <- model_summary_and_trace_obj_META_REG_DUMMY$summary(  compute_main_params = TRUE,
                                                     compute_transformed_parameters = FALSE, 
                                                     compute_generated_quantities = TRUE,
                                                     ##
                                                     save_log_lik_trace = FALSE,
                                                     ##
                                                     use_BayesMVP_for_faster_summaries = TRUE,
                                                     ##
                                                     compute_nested_rhat = FALSE)
##
summarise

##
## ----------- Plot sROC curve(s) (NMA - for best model WITH covariates): ---------------------------------------------
##
MR_Model_dummy <- readRDS(
  file.path("application_results_dummy_data", 
            "seed_1_dummy_data_1intercept_only_0N_covs_4param_Xurand_thr_FALSECS_1het_sigma_0_MR_model_6_applied_results.RDS"))

model_summary_and_trace_obj_META_REG_DUMMY <- MR_Model_dummy$model_summary_and_trace_obj
 
MR_Model_dummy_X <- subset_X_covariates( X = sim_results$X_mat,
                          covariates_to_keep = c("intercept",
                                                 "logit_prev_GAD",
                                                 "Ref_test_clean_SCID", "Ref_test_clean_Structured"))

baseline_case <- c(1, 
                   0.0,
                   0, 0)


 

stan_mod_samples     = model_summary_and_trace_obj_META_REG_DUMMY$internal_obj$outs_stan_sampling$stan_mod_samples


draws_array <- model_summary_and_trace_obj_META_REG_DUMMY$get_all_traces_list()$traces_as_arrays$draws_array
# trace_gq <-  model_summary_and_trace_obj_META_REG_DUMMY$get_all_traces_list()$traces_as_arrays$trace_gq
# str(trace_gq)
##
outs_NMA_comp_DUMMY <-  R_fn_using_Rcpp_compute_NMA_comparisons_posthoc(  trace_gq = draws_array,
                                                                          n_thr = n_thr,
                                                                          test_names = test_names)
##
str(outs_NMA_comp_DUMMY)
outs_NMA_comp_DUMMY

 

n_params_full <- n_par_inc_tp_and_gq
all_param_outs_trace <-    (BayesMVPtest:::fn_compute_param_constrain_from_trace_parallel(     unc_params_trace_input_main = trace_main_as_list,
                                                                                               unc_params_trace_input_nuisance = trace_nuisance_as_list,
                                                                                               pars_indicies_to_track = pars_indicies_to_track,
                                                                                               n_params_full = n_params_full,
                                                                                               n_nuisance = 0,
                                                                                               n_params_main = n_params_main,
                                                                                               include_nuisance = FALSE,
                                                                                               model_so_file = model_so_file,
                                                                                               json_file_path = json_file_path))




str(all_param_outs_trace)


## R_fn_using_Rcpp_compute_NMA_comparisons_posthoc
 

 
# 
# # 

# 
# # bs_model$param_constrain(theta_unc = )
# 
# prev_val <- -1.0
# 
# baseline_case_d <- baseline_case_nd <- list()
# ##
# baseline_case_d[["GAD-2"]] <- c( 1, prev_val,  ## mean prev_GAD
#                                  0, 0) ## MINI
# baseline_case_nd[["GAD-2"]] <- baseline_case_d[["GAD-2"]]
# ##
# baseline_case_d[["GAD-7"]] <- c( 1, prev_val,  ## mean prev_GAD
#                                  0, 0) ## MINI
# baseline_case_nd[["GAD-7"]] <- baseline_case_d[["GAD-7"]]
# ##
# baseline_case_d[["HADS"]] <- c( 1, prev_val,  ## mean prev_GAD
#                                 0, 0) ## MINI
# baseline_case_nd[["HADS"]] <- baseline_case_d[["HADS"]]
# ##
# baseline_case_d[["BAI"]] <- c( 1, prev_val,  ## mean prev_GAD
#                                0, 0) ## MINI
# baseline_case_nd[["BAI"]] <- baseline_case_d[["BAI"]]
# ##
# ## ---- check:
# ##
# baseline_case_d
# baseline_case_nd


model_prep_obj <- MR_Model_dummy$model_prep_obj
model_summary_and_trace_obj <- MR_Model_dummy$model_summary_and_trace_obj




##
stan_mod_samples_to_use <- R_fn_GQ_covariates_NMA(  baseline_case_nd = baseline_case_nd,
                                                    baseline_case_d = baseline_case_d,
                                                    model_summary_and_trace_obj = model_summary_and_trace_obj,
                                                    model_prep_obj = model_prep_obj)



##
## ------------- Compute new Se/Sp baselines + PI's + AUC's: -------------------------------------------------------
##
new_baseline_outs <- R_fn_using_Rcpp_compute_MR_Se_Sp_AUC_baseline(  
                                      model_prep_obj = MR_Model_dummy$model_prep_obj,
                                      model_summary_and_trace_obj = MR_Model_dummy$model_summary_and_trace_obj,
                                      baseline_case_nd = baseline_case_nd,
                                      baseline_case_d = baseline_case_d,
                                      use_probit_link = TRUE,
                                      test_names = test_names)

##
## ---- New AUC outs for new MR covariate baseline:
##
new_baseline_outs$AUC_summary
new_baseline_outs$AUC_pred_summary
new_baseline_outs$AUC_diff_summary
##
new_baseline_outs$Se_summaries %>% print(n=25)
new_baseline_outs$Sp_summaries %>% print(n=25)
##
new_baseline_outs$Se_pred_summaries %>% print(n=25)
new_baseline_outs$Sp_pred_summaries %>% print(n=25)
# ##
# ## sanity check to make sure matches stan -- ONLY WILL MATCH IF USING SAME MR-COVARIATE BASELINE!!!!
# ##
# options(scipen = 999)
# ##
# outs$Se_summaries %>% print(n=25)
# outs_Stan <- stan_mod_samples$summary(c("Se_baseline"),
#                                       mean, 
#                                       quantiles = ~ quantile(., na.rm = TRUE, probs = c(0.025, 0.50, 0.975))) %>%  print(n=25)
# ##
# outs$Se_summaries$`50%` - outs_Stan$`50%`
# ##
# outs$Se_pred_summaries %>%  print(n=25)
# outs_Stan_pred <- stan_mod_samples$summary(c("Se_baseline_pred"), 
#                                            mean, 
#                                            quantiles = ~ quantile(., na.rm = TRUE, probs = c(0.025, 0.50, 0.975))) %>%  print(n=25)
# # Compare predictive values
# outs$Se_pred_summaries$`50%` - outs_Stan_pred$`50%`

##
require(tictoc)
{
  tic()
  new_baseline_outs_orig <- R_fn_compute_MR_complete_new_baseline_analysis_v1(
                                    debugging = FALSE,
                                    model_prep_obj = MR_Model_dummy$model_prep_obj,
                                    model_summary_and_trace_obj = MR_Model_dummy$model_summary_and_trace_obj,
                                    baseline_case_nd = baseline_case_nd,
                                    baseline_case_d = baseline_case_d,
                                    use_probit_link = TRUE,
                                    test_names = test_names,
                                    compute_comparisons = TRUE)
  toc()
}

new_baseline_outs_orig$comparison_results
##
{
  tic()
  new_baseline_outs <- R_fn_compute_MR_complete_new_baseline_analysis_v2(
                                    debugging = FALSE,
                                    model_prep_obj = MR_Model_dummy$model_prep_obj,
                                    model_summary_and_trace_obj = MR_Model_dummy$model_summary_and_trace_obj,
                                    baseline_case_nd = baseline_case_nd,
                                    baseline_case_d = baseline_case_d,
                                    use_probit_link = TRUE,
                                    test_names = test_names,
                                    compute_comparisons = TRUE)
  toc()
}


new_baseline_outs$comparison_results




# plots_DUMMY <- model_summary_and_trace_obj_DUMMY$plot_sROC(test_names = test_names)
##
plots_DUMMY <- R_fn_sROC_plot_NMA(    stan_model_file_name = model_summary_and_trace_obj_DUMMY$internal_obj$outs_stan_model_name$stan_model_file_name,
                                      ##
                                      stan_mod_samples     = model_summary_and_trace_obj_DUMMY$internal_obj$outs_stan_sampling$stan_mod_samples,
                                      tibbles_provided = TRUE,
                                      tibbles = new_baseline_outs,
                                      ##
                                      test_names = test_names,
                                      ##
                                      relevant_thresholds = relevant_thresholds,
                                      ##
                                      # covariates = covariates,
                                      # new_cov_data = new_cov_data,
                                      ##
                                      df_true = NULL,
                                      ##
                                      conf_region_colour = "blue", 
                                      pred_region_colour = "green",
                                      ##
                                      conf_region_alpha = 0.20,
                                      pred_region_alpha = 0.10,
                                      ##
                                      base_size = 24,
                                      point_size = 4.0,
                                      linewidth = 1.0,
                                      ##
                                      n_index_tests = n_index_tests,
                                      n_thr = n_thr)
##
plots_DUMMY$plot_list[[1]]
#### plots_DUMMY$plot_list[[2]]
#### plots_DUMMY$plot_list[[3]]
#### plots_DUMMY$plot_list[[4]]
plots_DUMMY$plot_list[[5]]
#### plots_DUMMY$plot_list[[6]]
plots_DUMMY$plot_list[[7]]
plots_DUMMY$plot_list[[8]]
plots_DUMMY$plot_list[[9]]
##
ggsave("DUMMY_Application_baseline_sROC.png", plots_DUMMY$plot_list[[1]], width = 16, height = 9, dpi = 400)
ggsave("DUMMY_Application_baseline_sROC_panel_w_CrI_PrI.png", plots_DUMMY$plot_list[[5]], width = 16, height = 16, dpi = 400)
##
ggsave("DUMMY_Application_baseline_Se_Sp_vs_thr.png", plots_DUMMY$plot_list[[7]], width = 16, height = 16, dpi = 400)
ggsave("DUMMY_Application_baseline_Se_Sp_vs_thr_w_CrI.png", plots_DUMMY$plot_list[[8]], width = 16, height = 16, dpi = 400)
ggsave("DUMMY_Application_baseline_Se_Sp_vs_thr_panel_w_CrI_PrI.png", plots_DUMMY$plot_list[[9]], width = 16, height = 16, dpi = 400)
# ##
# ## ---- Compare real to simulated/dummy data (for personal interest NOT FOR PAPER/THESIS):
# ##
# require(patchwork)
# plots_DUMMY$plot_list[[1]] + plots_REAL$plot_list[[1]]

 
 
 




 
# 
# 
# str(outs$Se_baseline)
# 
# outs$Se_summaries

# # Compile and fit model
Stan_model_file_path <- model_summary_and_trace_obj_META_REG_DUMMY$internal_obj$outs_stan_compile$stan_model_file_path
Stan_functions_path <- model_summary_and_trace_obj_META_REG_DUMMY$internal_obj$outs_stan_compile$stan_functions_directory
# ##
# mod <- cmdstan_model(Stan_model_file_path,
#                      include_paths = c(Stan_functions_path))

##
# ## set Stan model file path for your Stan model (replace with your path)
Stan_model_file_path <- system.file("stan_models/basic_logistic.stan", package = "BayesMVPtest")
# make y a matrix first (matrix w/ 1 col)
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
  Stan_data_list <- stan_data_list
  ##
  Stan_data_list$y <- NULL
  Stan_data_list$N <- NULL
  ##
  Stan_data_list$softplus <- 0
  ##r
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
str(Stan_data_list)
Stan_data_list$softplus <- NULL  # Remove one instance
Stan_data_list$softplus <- 0     # Set it once
##
stanc_args <- list(paste0("--include-paths=", Stan_functions_path))
##
#  



outs <- BayesMVPtest:::get_n_params_main_from_n_nuisance_using_Stan_model( n_nuisance = 0,
                                                                           Stan_data_list = stan_data_list,
                                                                           Stan_model_file_path = Stan_model_file_path,
                                                                           stanc_args = stanc_args)


bs_model <- outs$bs_model

bs_model$param_constrain(theta_unc = )

# 
# n_params_main <- outs$n_params_main
# n_params <- outs$n_params
# 
# # ###  -----------  Compile + initialise the model using "MVP_model$new(...)" 
# BayesMVP_model_prep_obj <- BayesMVPtest:::MVP_model$new(   Model_type =  "Stan",
#                                             y = NULL,
#                                             N = NULL,
#                                             ##  model_args_list = model_args_list, # this arg is only needed for BUILT-IN (not Stan) models
#                                             Stan_data_list = Stan_data_list,
#                                             Stan_model_file_path = Stan_model_file_path,
#                                             init_lists_per_chain = init_lists_per_chain,
#                                             sample_nuisance = FALSE,
#                                             n_chains_burnin = 8,
#                                             ##
#                                             n_params_main = n_params_main,
# #                                             n_nuisance = 0,
# #                                             ##
# #                                             stanc_args = stanc_args)
# # 
# model_so_file <- BayesMVP_model_prep_obj$init_object$model_so_file
# json_file_path <- BayesMVP_model_prep_obj$init_object$json_file_path
# 
# # Get trace
# trace_gq <- model_summary_and_trace_obj_META_REG_DUMMY$get_trace_generated_quantities()
# trace_dims <- dim(trace_gq) ; trace_dims
# trace_vec <- as.vector(trace_gq)
# ##
# trace_main <- model_summary_and_trace_obj_META_REG_DUMMY$get_all_traces_list()$traces_as_arrays$trace_main
# trace_dims <- dim(trace_main) ; trace_dims
# trace_vec <- as.vector(trace_main)
# ##
# stan_functions_directory <- MR_Model_dummy$model_prep_obj$internal_obj$outs_stan_compile$stan_functions_directory ; stan_functions_directory
# stan_model_file_path     <- MR_Model_dummy$model_prep_obj$internal_obj$outs_stan_compile$stan_model_file_path ; stan_model_file_path
# stan_model_file_name <- MR_Model_dummy$model_prep_obj$internal_obj$outs_stan_compile$stan_model_file_name ; stan_model_file_name
# # stan_data_list           <- model_summary_and_trace_obj_META_REG_DUMMY$internal_obj$outs_data$stan_data_list 
# stan_data_list <- MR_Model_dummy$model_prep_obj$internal_obj$outs_data$stan_data_list
# str(stan_data_list)
# ##
# stan_mod_samples <- model_summary_and_trace_obj_META_REG_DUMMY$internal_obj$outs_stan_sampling$stan_mod_samples
# str(stan_mod_samples)
# ##
# new_cov_data <- list( baseline_case_nd = baseline_case_nd,
#                       baseline_case_d = baseline_case_d)
# ##
# # new_cov_data <- model_obj$cov_data
# 
# 
# 
# str(trace_main)
# 
# n_par_inc_tp_and_gq <- length(bs_model$param_names(include_tp = TRUE, include_gq = TRUE))
# 
# 
# trace_nuisance <- array(0, dim = c(500, 32, 10))
# trace_main_as_list <- trace_nuisance_as_list <- list()
# 
# for (kk in 1:32) { 
#   trace_main_as_list[[kk]] <- t(trace_main[, kk, ])
#   trace_nuisance_as_list[[kk]] <- t(trace_nuisance[, kk, ])
# }
# #
# 
# ## print(n_params_main)
# pars_indicies_to_track <- 1:n_par_inc_tp_and_gq
# n_params_full <- n_par_inc_tp_and_gq
# all_param_outs_trace <-    (BayesMVPtest:::fn_compute_param_constrain_from_trace_parallel(     unc_params_trace_input_main = trace_main_as_list,
#                                                                                                unc_params_trace_input_nuisance = trace_nuisance_as_list,
#                                                                                                pars_indicies_to_track = pars_indicies_to_track,
#                                                                                                n_params_full = n_params_full,
#                                                                                                n_nuisance = 0,
#                                                                                                n_params_main = n_params_main,
#                                                                                                include_nuisance = FALSE,
#                                                                                                model_so_file = model_so_file,
#                                                                                                json_file_path = json_file_path))
# 
# 
# 
# 
# str(all_param_outs_trace)
# 
# 
# 
# 
# 
# 
# ##
# stan_mod_samples_to_use <- R_fn_GQ_covariates_NMA(  new_cov_data = new_cov_data,
#                                                     stan_mod_samples = stan_mod_samples,
#                                                     stan_data_list = stan_data_list,
#                                                     stan_model_file_path = stan_model_file_path,
#                                                     stan_functions_directory = stan_functions_directory)
# 
# 



  
sum(MR_Model_dummy_X[[1]][[1]][, 3])/n_studies

X_summary_dummy <- R_fn_summarize_covariates_general(X_list = MR_Model_dummy_X[[1]],
                                                     test_names = test_names,
                                                     continuous_vars = c("logit_prev_GAD"),
                                                     binary_vars = c("Ref_test_clean_SCID", "Ref_test_clean_Structured")
                                                     # logit_vars = c("logit_prev_GAD")  # doesnt work!
                                                     )
                                                     
X_summary_dummy
##
# X_summary_real <- R_fn_summarize_covariates_general(X_list = real_data$cov_data$X[[1]])
# X_summary_real
##
beta_mu <- model_summary_and_trace_obj_META_REG_DUMMY$extract_params(params = c("beta_mu")) %>% print(n = 100)
##
beta_tibble
covariate_names
test_names
##


 

outs <- R_fn_map_beta_coefficients(beta_tibble = beta_mu, 
                           covariate_names = dimnames(MR_Model_dummy_X[[1]][[1]])[[2]], 
                           test_names = test_names, 
                           min_magnitude = 0.125,
                           borderline_threshold = 0.10,
                           show_borderline = TRUE)
##
## outs$full_table %>% print(n = 100)
##
outs$full_table %>%
  filter(abs(`50%`) > 0.125) %>%
  filter(is_significant == FALSE, covariate_name != "intercept") %>%
  print(n = 100)
##
outs_latex <- create_metareg_coef_table( beta_tibble = outs$full_table, 
                                         test_names = test_names)
cat(outs_latex)
# \begin{table}[!h]
# \centering
# \caption{Meta-regression coefficients (probit scale)}
# \small
# \begin{tabular}{lcc}
# \toprule
# & \multicolumn{1}{c}{Non-diseased} & \multicolumn{1}{c}{Diseased} \\
# \midrule
# \multicolumn{3}{l}{\textbf{GAD-2}} \\
# Intercept            & $\mathbf{-0.970 ~ \ci{-1.751}{-0.204}}$ & $0.295 ~ \ci{-0.483}{1.057}$ \\
# Prevalence (per SD)  & $0.089 ~ \ci{-0.023}{0.198}$ & $\mathbf{0.329 ~ \ci{0.208}{0.446}}$ \\
# Ref test: SCID       & $0.046 ~ \ci{-0.260}{0.346}$ & $0.233 ~ \ci{-0.066}{0.531}$ \\
# Ref_test_clean_Structured & $\mathbf{0.581 ~ \ci{0.123}{1.023}}$ & $-0.399 ~ \ci{-0.845}{0.040}$ \\
# %%%%
#   \midrule
# %%%%
#   \multicolumn{3}{l}{\textbf{GAD-7}} \\
# Intercept            & $\mathbf{-1.439 ~ \ci{-1.892}{-0.993}}$ & $0.109 ~ \ci{-0.350}{0.557}$ \\
# Prevalence (per SD)  & $\mathbf{0.113 ~ \ci{0.005}{0.221}}$ & $\mathbf{0.283 ~ \ci{0.166}{0.398}}$ \\
# Ref test: SCID       & $0.187 ~ \ci{-0.033}{0.405}$ & $\mathbf{0.302 ~ \ci{0.059}{0.544}}$ \\
# Ref_test_clean_Structured & $0.113 ~ \ci{-0.272}{0.486}$ & $-0.369 ~ \ci{-0.761}{0.019}$ \\
# %%%%
#   \midrule
# %%%%
#   \multicolumn{3}{l}{\textbf{HADS}} \\
# Intercept            & $\mathbf{-1.131 ~ \ci{-1.597}{-0.676}}$ & $0.197 ~ \ci{-0.275}{0.666}$ \\
# Prevalence (per SD)  & $0.012 ~ \ci{-0.158}{0.178}$ & $\mathbf{-0.317 ~ \ci{-0.499}{-0.138}}$ \\
# Ref test: SCID       & $0.030 ~ \ci{-0.249}{0.306}$ & $-0.248 ~ \ci{-0.555}{0.052}$ \\
# Ref_test_clean_Structured & $-0.148 ~ \ci{-0.554}{0.248}$ & $-0.244 ~ \ci{-0.689}{0.201}$ \\
# %%%%
#   \midrule
# %%%%
#   \multicolumn{3}{l}{\textbf{BAI}} \\
# Intercept            & $\mathbf{-1.821 ~ \ci{-2.198}{-1.442}}$ & $-0.232 ~ \ci{-0.651}{0.176}$ \\
# Prevalence (per SD)  & $0.199 ~ \ci{-0.056}{0.451}$ & $-0.116 ~ \ci{-0.399}{0.167}$ \\
# Ref test: SCID       & $0.185 ~ \ci{-0.299}{0.660}$ & $-0.228 ~ \ci{-0.763}{0.293}$ \\
# Ref_test_clean_Structured & $0.001 ~ \ci{-1.005}{0.960}$ & $-0.001 ~ \ci{-1.001}{0.967}$ \\
# %%%%
#   \bottomrule
# %%%%
#   \end{tabular}
# \begin{tablenotes}
# \item Bold indicates significant at 95\% level
# \end{tablenotes}
# \end{table}



##
## ---- Load chosen intercept-only model (i.e. the MAIN model):
##
Model_dummy <- readRDS(
  file.path("application_results_dummy_data", 
            "seed_1_dummy_data_1intercept_only_1N_covs_1param_Xurand_thr_FALSECS_1het_sigma_0_applied_results.RDS"))

model_summary_and_trace_obj_DUMMY <- Model_dummy$model_summary_and_trace_obj

 
 
{
  ##
  tibble_total_SD_inc_C_MR <- model_summary_and_trace_obj_META_REG_DUMMY$extract_params(params = c("total_SD_inc_C")) %>% 
    select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
    print(n = 100)
  ##
  tibble_beta_sigma_MR <- model_summary_and_trace_obj_META_REG_DUMMY$extract_params(params = c("beta_sigma")) %>% 
    select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
    print(n = 100)
  ##
  tibble_beta_tau_MR <- model_summary_and_trace_obj_META_REG_DUMMY$extract_params(params = c("beta_tau")) %>% 
    select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
    print(n = 100)
  ##
  tibble_beta_corr_MR <- model_summary_and_trace_obj_META_REG_DUMMY$extract_params(params = c("beta_corr")) %>% 
    select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
    print(n = 100)
  ##
  tibble_rho12_MR <- model_summary_and_trace_obj_META_REG_DUMMY$extract_params(params = c("rho12")) %>% 
    select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
    print(n = 100)
}
{
  ##
  tibble_total_SD_inc_C <- model_summary_and_trace_obj_DUMMY$extract_params(params = c("total_SD_inc_C")) %>% 
    select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
    print(n = 100)
  ##
  tibble_beta_sigma <- model_summary_and_trace_obj_DUMMY$extract_params(params = c("beta_sigma")) %>% 
    select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
    print(n = 100)
  ##
  tibble_beta_tau <- model_summary_and_trace_obj_DUMMY$extract_params(params = c("beta_tau")) %>% 
    select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
    print(n = 100)
  ##
  tibble_beta_corr <- model_summary_and_trace_obj_DUMMY$extract_params(params = c("beta_corr")) %>% 
    select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
    print(n = 100)
  ##
  tibble_rho12 <- model_summary_and_trace_obj_DUMMY$extract_params(params = c("rho12")) %>% 
    select(parameter, `2.5%`, `50%`, `97.5%`, n_eff, Rhat) %>% 
    print(n = 100)
}


outs_latex <- create_heterogeneity_comparison_table( total_SD_intercept = tibble_total_SD_inc_C,
                                       beta_sigma_intercept = tibble_beta_sigma,
                                       beta_tau_intercept = tibble_beta_tau,
                                       beta_corr_intercept = tibble_beta_corr,
                                       rho12_intercept = tibble_rho12,
                                       ##
                                       total_SD_covariate = tibble_total_SD_inc_C_MR,
                                       beta_sigma_covariate = tibble_beta_sigma_MR,
                                       beta_tau_covariate = tibble_beta_tau_MR,
                                       beta_corr_covariate = tibble_beta_corr_MR,
                                       rho12_covariate = tibble_rho12_MR)
cat(outs_latex)
# \begin{table}[!h]
# \centering
# \caption{Impact of covariates on heterogeneity parameters}
# \small
# \begin{tabular}{lcc}
# \toprule
# \textbf{Parameter}  & \textbf{Model A}          & \textbf{Model A + Covariates} \\
# & \textit{(intercept-only)} & \textit{(prevalence, ref test, setting)} \\
# %%%%
#   \midrule
# %%%%
#   \multicolumn{3}{l}{\textit{\textbf{NMA shared between-study SD's:}}} \\
# $\sigma_{\beta}^{(d-)}$ & $0.215 ~ \ci{0.101}{0.304}$ & $0.207 ~ \ci{0.098}{0.292}$ \\
# $\sigma_{\beta}^{(d+)}$ & $0.421 ~ \ci{0.319}{0.531}$ & $0.358 ~ \ci{0.289}{0.439}$ \\
# %%%%
# \midrule
# %%%%
# \multicolumn{3}{l}{\textit{\textbf{NMA test-specific deviation SD's (CS structure):}}} \\
# $\tau^{(d-)}$ (all tests) & $0.287 ~ \ci{0.232}{0.352}$ & $0.271 ~ \ci{0.218}{0.338}$ \\
# $\tau^{(d+)}$ (all tests) & $0.329 ~ \ci{0.264}{0.417}$ & $0.162 ~ \ci{0.118}{0.221}$ \\
# %%%%
#   \midrule
# %%%%
#   \multicolumn{3}{l}{\textit{\textbf{Total within-test, between-study heterogeneity:}}} \\
# $\sigma_{total}^{(d-)}$ (all tests) & $0.360 ~ \ci{0.314}{0.415}$ & $0.344 ~ \ci{0.299}{0.399}$ \\
# $\sigma_{total}^{(d+)}$ (all tests) & $0.537 ~ \ci{0.468}{0.627}$ & $0.395 ~ \ci{0.337}{0.469}$ \\
# %%%%
#   \midrule
# %%%%
#   \multicolumn{3}{l}{\textit{\textbf{Correlation structure:}}} \\
# $\rho_{\beta}$ (NMA shared corr.) & $0.478 ~ \ci{0.096}{0.788}$ & $0.423 ~ \ci{0.045}{0.728}$ \\
# $\rho_{12}$ (within-test corr.) & $0.210 ~ \ci{0.037}{0.383}$ & $0.223 ~ \ci{0.024}{0.406}$ \\
# %%%%
#   \bottomrule
# %%%%
#   \end{tabular}
# %%%%
#   \begin{tablenotes}
# \footnotesize
# \item Both models use fixed cutpoints with compound symmetry correlation (Model A specification).
# \item CS = compound symmetry; all tests share the same variance components.
# \end{tablenotes}
# %%%%
#   \end{table}






 


##
## ---- Extract NMA AUC metrics: ----------------------------------------------------------------------------------------
##
 
 
# ##
# ## ---- For dummy applied data:
# ##
# baseline_case_nd <- model_summary_and_trace_obj_META_REG_DUMMY$cov_data$baseline_case_nd
# baseline_case_d <- model_summary_and_trace_obj_META_REG_DUMMY$cov_data$baseline_case_d
# ##
# baseline_case_d <- baseline_case_nd <- list()
# ##
# baseline_case_d[["GAD-2"]] <- c( 1, 0.0,  ## mean prev_GAD
#                                  0, 0) ## MINI
# baseline_case_nd[["GAD-2"]] <- baseline_case_d[["GAD-2"]]
# ##
# baseline_case_d[["GAD-7"]] <- c( 1, 0.0,  ## mean prev_GAD
#                                  0, 0) ## MINI
# baseline_case_nd[["GAD-7"]] <- baseline_case_d[["GAD-7"]]
# ##
# baseline_case_d[["HADS"]] <- c( 1, 0.0,  ## mean prev_GAD
#                                 0, 0) ## MINI
# baseline_case_nd[["HADS"]] <- baseline_case_d[["HADS"]]
# ##
# baseline_case_d[["BAI"]] <- c( 1, 0.0,  ## mean prev_GAD
#                                0, 0) ## MINI
# baseline_case_nd[["BAI"]] <- baseline_case_d[["BAI"]]
# ##
# ## ---- check:
# ##
# baseline_case_d
# baseline_case_nd

# MR_AUC_outs_DUMMY <- R_fn_using_Rcpp_compute_NMA_AUC_custom_baseline(  model_obj = model_summary_and_trace_obj_META_REG_DUMMY,
#                                                                        baseline_case_nd = baseline_case_nd,
#                                                                        baseline_case_d = baseline_case_d,
#                                                                        n_thr = n_thr,
#                                                                        test_names = test_names)

# 
# MR_AUC_outs_DUMMY <- R_fn_extract_NMA_AUC_summary(  tibble_gq = model_summary_and_trace_obj_META_REG_DUMMY$get_summary_generated_quantities(),
#                                                  n_thr = n_thr,
#                                                  test_names = test_names,
#                                                  network = network)
##
MR_AUC_outs_DUMMY$auc
MR_AUC_outs_DUMMY$auc_diff
MR_AUC_outs_DUMMY$auc_pred
# # A tibble: 4 × 7
# test test_name AUC_mean AUC_median AUC_lower AUC_upper AUC_sd
# <dbl> <chr>        <dbl>      <dbl>     <dbl>     <dbl>  <dbl>
# 1     1 GAD-2        0.785      0.785     0.745     0.821 0.0196
# 2     2 GAD-7        0.809      0.809     0.771     0.843 0.0187
# 3     3 HADS         0.873      0.873     0.843     0.900 0.0147
# 4     4 BAI          0.865      0.866     0.814     0.906 0.0234
# > MR_AUC_outs_DUMMY$auc_diff
# # A tibble: 6 × 10
# test1 test2 test1_name test2_name AUC_diff_mean AUC_diff_median AUC_diff_lower AUC_diff_upper AUC_diff_sd prob_test1_better
# <int> <int> <chr>      <chr>              <dbl>           <dbl>          <dbl>          <dbl>       <dbl>             <dbl>
# 1     1     2 GAD-2      GAD-7           -0.0240         -0.0240         -0.0721        0.0226       0.0245                NA
# 2     1     3 GAD-2      HADS            -0.0880         -0.0876         -0.136        -0.0428       0.0238                NA
# 3     2     3 GAD-7      HADS            -0.0798         -0.0803         -0.137        -0.0211       0.0300                NA
# 4     1     4 GAD-2      BAI             -0.0640         -0.0637         -0.110        -0.0203       0.0231                NA
# 5     2     4 GAD-7      BAI             -0.0559         -0.0568         -0.112         0.00300      0.0298                NA
# 6     3     4 HADS       BAI              0.00815         0.00715        -0.0420        0.0622       0.0269                NA
# MR_AUC_outs_DUMMY$auc_pred
# # A tibble: 4 × 7
# test test_name AUC_mean AUC_median AUC_lower AUC_upper AUC_sd
# <dbl> <chr>        <dbl>      <dbl>     <dbl>     <dbl>  <dbl>
# 1     1 GAD-2        0.763      0.780     0.488     0.941 0.120 
# 2     2 GAD-7        0.782      0.801     0.500     0.952 0.120 
# 3     3 HADS         0.849      0.870     0.602     0.975 0.0981
# 4     4 BAI          0.840      0.862     0.588     0.974 0.102 
##
## ---- Make LaTeX table:
##
AUC_latex_table_DUMMY <- create_AUC_latex_table( auc_outs = MR_AUC_outs_DUMMY, 
                                                 digits = 3)
cat(AUC_latex_table_DUMMY)


##
## ----------- Plot sROC curve(s) (NMA - for best model WITHOUT covariates): ------------------------------------------
##
relevant_thresholds = list( "GAD-2" = c(1:6), 
                            "GAD-7" =  c(3:18), 
                            "HADS" =  c(3:18),
                            "BAI" = c(3:38))
# stan_model_file_name = model_summary_and_trace_obj_DUMMY$internal_obj$outs_stan_model_name$stan_model_file_name
# stan_mod_samples     = model_summary_and_trace_obj_DUMMY$internal_obj$outs_stan_sampling$stan_mod_samples
# ##
# test_names = test_names
# # covariates = covariates,
# # new_cov_data = new_cov_data,
# ##
# # df_true = df_true
# ##
# conf_region_colour = "blue"
# pred_region_colour = "green"
# ##
# base_size = 24
# point_size = 4.0
# linewidth = 1.0
# ##
# n_index_tests = n_index_tests
# n_thr = n_thr
#
#
##
# plots_REAL <- model_summary_and_trace_obj_REAL$plot_sROC(test_names = test_names)
##
plots_REAL <- R_fn_sROC_plot_NMA(    stan_model_file_name = model_summary_and_trace_obj_REAL$internal_obj$outs_stan_model_name$stan_model_file_name,
                                     stan_mod_samples     =  model_summary_and_trace_obj_REAL$internal_obj$outs_stan_sampling$stan_mod_samples,
                                     ##
                                     test_names = test_names,
                                     ##
                                     relevant_thresholds = relevant_thresholds,
                                     ##
                                     # covariates = covariates,
                                     # new_cov_data = new_cov_data,
                                     ##
                                     df_true = NULL,
                                     ##
                                     conf_region_colour = "blue", 
                                     pred_region_colour = "green",
                                     ##
                                     conf_region_alpha = 0.30,
                                     pred_region_alpha = 0.15,
                                     ##
                                     base_size = 24,
                                     point_size = 4.0,
                                     linewidth = 1.0,
                                     ##
                                     n_index_tests = n_index_tests,
                                     n_thr = n_thr)
##
plots_REAL$plot_list[[1]]
#### plots_REAL$plot_list[[2]]
#### plots_REAL$plot_list[[3]]
#### plots_REAL$plot_list[[4]]
plots_REAL$plot_list[[5]]
#### plots_REAL$plot_list[[6]]
plots_REAL$plot_list[[7]]
plots_REAL$plot_list[[8]]
plots_REAL$plot_list[[9]]
##
ggsave("REAL_Application_baseline_sROC.png", plots_REAL$plot_list[[1]], width = 16, height = 9, dpi = 400)
ggsave("REAL_Application_baseline_sROC_panel_w_CrI_PrI.png", plots_REAL$plot_list[[6]], width = 16, height = 9, dpi = 400)
##
ggsave("REAL_Application_baseline_Se_Sp_vs_thr.png", plots_REAL$plot_list[[7]], width = 16, height = 16, dpi = 400)
ggsave("REAL_Application_baseline_Se_Sp_vs_thr_w_CrI.png", plots_REAL$plot_list[[8]], width = 16, height = 16, dpi = 400)
ggsave("REAL_Application_baseline_Se_Sp_vs_thr_panel_w_CrI_PrI.png", plots_REAL$plot_list[[9]], width = 16, height = 16, dpi = 400)
##
##
##
# plots_DUMMY <- model_summary_and_trace_obj_DUMMY$plot_sROC(test_names = test_names)
##
plots_DUMMY <- R_fn_sROC_plot_NMA(    stan_model_file_name = model_summary_and_trace_obj_DUMMY$internal_obj$outs_stan_model_name$stan_model_file_name,
                                      stan_mod_samples     = model_summary_and_trace_obj_DUMMY$internal_obj$outs_stan_sampling$stan_mod_samples,
                                      ##
                                      test_names = test_names,
                                      ##
                                      relevant_thresholds = relevant_thresholds,
                                      ##
                                      # covariates = covariates,
                                      # new_cov_data = new_cov_data,
                                      ##
                                      df_true = NULL,
                                      ##
                                      conf_region_colour = "blue", 
                                      pred_region_colour = "green",
                                      ##
                                      conf_region_alpha = 0.20,
                                      pred_region_alpha = 0.10,
                                      ##
                                      base_size = 24,
                                      point_size = 4.0,
                                      linewidth = 1.0,
                                      ##
                                      n_index_tests = n_index_tests,
                                      n_thr = n_thr)
##
plots_DUMMY$plot_list[[1]]
#### plots_DUMMY$plot_list[[2]]
#### plots_DUMMY$plot_list[[3]]
#### plots_DUMMY$plot_list[[4]]
plots_DUMMY$plot_list[[5]]
#### plots_DUMMY$plot_list[[6]]
plots_DUMMY$plot_list[[7]]
plots_DUMMY$plot_list[[8]]
plots_DUMMY$plot_list[[9]]
##
ggsave("DUMMY_Application_baseline_sROC.png", plots_DUMMY$plot_list[[1]], width = 16, height = 9, dpi = 400)
ggsave("DUMMY_Application_baseline_sROC_panel_w_CrI_PrI.png", plots_DUMMY$plot_list[[5]], width = 16, height = 16, dpi = 400)
##
ggsave("DUMMY_Application_baseline_Se_Sp_vs_thr.png", plots_DUMMY$plot_list[[7]], width = 16, height = 16, dpi = 400)
ggsave("DUMMY_Application_baseline_Se_Sp_vs_thr_w_CrI.png", plots_DUMMY$plot_list[[8]], width = 16, height = 16, dpi = 400)
ggsave("DUMMY_Application_baseline_Se_Sp_vs_thr_panel_w_CrI_PrI.png", plots_DUMMY$plot_list[[9]], width = 16, height = 16, dpi = 400)
# ##
# ## ---- Compare real to simulated/dummy data (for personal interest NOT FOR PAPER/THESIS):
# ##
# require(patchwork)
# plots_DUMMY$plot_list[[1]] + plots_REAL$plot_list[[1]]











 
 












 
# 
# 
# 
# 
# {
#   
#   
#   
#   basic_model_options = list(
#     network = NULL,
#     cts = NULL,
#     prior_only = NULL
#   )
#   ####
#   advanced_model_options = list(
#     model_parameterisation = NULL,
#     random_thresholds = NULL,
#     Dirichlet_random_effects_type = NULL,
#     box_cox = NULL,
#     softplus = NULL
#   )
#   ####
#   other_advanced_options = list(
#     advanced_compile = NULL, ## default is false (aka a "basic" compile)
#     ##
#     force_recompile = NULL,
#     quiet = NULL,
#     compile = NULL,
#     ##
#     set_custom_CXX_CPP_flags = NULL,
#     CCACHE_PATH = NULL,
#     custom_cpp_user_header_file_path = NULL,
#     CXX_COMPILER_PATH = NULL,
#     CPP_COMPILER_PATH = NULL,
#     MATH_FLAGS = NULL,
#     FMA_FLAGS = NULL,
#     AVX_FLAGS = NULL,
#     THREAD_FLAGS = NULL
#   )
#   ##
#   ####
#   priors = NULL
#   ####
#   # init_lists_per_chain = NULL,
#   ####
#   MCMC_params = list(
#     seed = NULL,
#     n_superchains = NULL,
#     n_chains = NULL,
#     n_iter = NULL,
#     n_burnin = NULL,
#     adapt_delta = NULL,
#     max_treedepth = NULL,
#     metric_shape = NULL
#   )
#   
#   ##
#   ## ---- Main "internal_obj" list:
#   ##
#   internal_obj = list(
#     outs_data = NULL,
#     outs_stan_model_name = NULL,
#     outs_stan_compile = NULL,
#     outs_stan_init = NULL,
#     outs_stan_sampling = NULL,
#     ##
#     HMC_info = NULL, ## new
#     efficiency_info = NULL, ## new
#     summaries = NULL, ## new
#     traces = NULL ## new
#   )
#   
#   debugging <- TRUE
#   ##
#   n_iter = 500
#   n_burnin = 500
#   ##
#   priors <- NULL
#   ##
#   x = x
#   n_chains = n_chains
#   ##
#   cts = FALSE
#   network = TRUE
#   ##
#   prior_only = FALSE
#   ##
#   softplus = FALSE
#   ##
#   ##
#   model_parameterisation = "R&G"
#   ##
#   random_thresholds = TRUE
#   Dirichlet_random_effects_type = "kappa"
#   ##
#   box_cox <- FALSE ## cannot input ACTUAL NA's into Stan!
#   ##
#   init_lists_per_chain = NULL
#   
#   
#   basic_model_options$network <- network
#   basic_model_options$cts <- cts
#   basic_model_options$prior_only <- prior_only
#   ##
#   advanced_model_options$model_parameterisation <- model_parameterisation
#   advanced_model_options$random_thresholds <- random_thresholds
#   advanced_model_options$Dirichlet_random_effects_type <- Dirichlet_random_effects_type
#   advanced_model_options$box_cox <- box_cox
#   advanced_model_options$softplus <- softplus
#   ##
#   MCMC_params$n_chains <- n_chains
#   
# }
# 
# 
# 
# 
# 
# outs_data = list(
#   stan_data_list = NULL,
#   n_tests = NULL,
#   n_studies = NULL,
#   n_thr = NULL,
#   n_cat = NULL
# )
# ##
# outs_stan_model_name = list(
#   stan_model_file_name = NULL
# )
# ##
# outs_stan_compile = list(
#   stan_model_obj = NULL,
#   stan_model_file_name = NULL,
#   stan_model_file_path = NULL,
#   ##
#   pkg_root_directory = NULL,
#   stan_models_directory = NULL,
#   stan_functions_directory = NULL,
#   ##
#   stan_MA_directory = NULL,
#   stan_MA_prior_directory = NULL,
#   ##
#   stan_NMA_directory = NULL,
#   stan_NMA_directory = NULL
# )
# ##
# outs_stan_init = list(
#   inits_unconstrained_vec_per_chain = NULL,
#   stan_param_names_list = NULL,
#   stan_param_names_main = NULL,
#   stan_init_pseudo_sampling_outs = NULL,
#   stan_model_obj = NULL,
#   json_file_path = NULL,
#   stan_model_file_path = NULL
# )
# ##
# outs_stan_sampling = list(
#   stan_mod_samples = NULL,
#   time_total = NULL
# )
# 
# X <- NULL
# 
# # x <- x_NMA
# 
# MetaOrdDTA:::prep_data_and_model(  
#   debugging = self$debugging,
#   ##
#   x = self$x,
#   ##
#   X = self$X, ## covariates
#   ##
#   indicator_index_test_in_study = self$indicator_index_test_in_study,
#   ##
#   internal_obj = self$internal_obj,
#   ##
#   basic_model_options    = self$basic_model_options,
#   advanced_model_options = self$advanced_model_options,
#   MCMC_params            = self$MCMC_params,
#   other_advanced_options = self$other_advanced_options,
#   ##
#   priors = self$priors,
#   init_lists_per_chain = self$init_lists_per_chain,
#   ##
#   compute_sim_study_metrics = self$compute_sim_study_metrics,
#   vec_index_inner_thr       = self$vec_index_inner_thr)
# 
# 
# 
# 
# 
# 
# 
# x = x
# X = NULL
# model_parameterisation = model_parameterisation
# n_index_tests_per_study = n_index_tests_per_study
# indicator_index_test_in_study = indicator_index_test_in_study
# 
# 
# 
# outs_data <- R_fn_prep_NMA_data(x = x,
#                                 X = NULL,
#                                 model_parameterisation = model_parameterisation,
#                                 n_index_tests_per_study = n_index_tests_per_study,
#                                 indicator_index_test_in_study = indicator_index_test_in_study)
# 
# 
# internal_obj$outs_data <- outs_data
# 



