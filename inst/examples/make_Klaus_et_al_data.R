


require(dplyr)
##
Klaus_dir <-  "/home/enzocerullo/Documents/Work/PhD_work/Klaus_collab_NMA_work"
##
source(file.path(local_pkg_dir, "Klaus_et_al_collab_R_fns.R"))
##
## read in data
# data_Klaus <- tibble(read.csv(file.path(Klaus_dir, "DIAQANDI NWMA database final 2025_04_23.csv"))) ## older data (used this for sim study!)
##
# data_Klaus <- tibble(read.csv(file.path(Klaus_dir, "DIAQANDI NWMA database final 2025_05_28.csv")))
data_Klaus <- tibble(read.csv(file.path(Klaus_dir, "DIAQANDI NWMA database updated 2025_11_25_2x2_tables.csv")))
data_Klaus
##
data <- data_Klaus %>%
        dplyr::mutate(id = 1:nrow(data_Klaus)) %>%
        dplyr::select(Study..,
                      Author,
                      Year,
                      Dx,
                      Index.test,
                      Cutoff,
                      n.with.dis.,
                      n.analyzed,
                      disease.prevalence..proportion.,
                      Sens,
                      Spez,
                      TP..a.,
                      FN...c.,
                      FP..b.,
                      TN..d.,
                      published.cut.off.in.analysis.,
                      Cut.off.published.,
                      Remark) %>%
        dplyr::rename(Study = Study..,
               Diagnosis = Dx,
               Index_test = Index.test,
               N_diseased = n.with.dis.,
               N = n.analyzed,
               prev = disease.prevalence..proportion.,
               Se = Sens,
               Sp = Spez,
               TP = TP..a.,
               FN = FN...c.,
               FP = FP..b.,
               TN = TN..d.)
##
test <- "GAD-2"
#### test <- "GAD-7"
test <- "HADS"
test <- "BAI"
##
# unique(data$Diagnosis)
# 
# 
# data$Cut.off.published.
# 
# data$id
# 
# sum(data$Index_test == "GAD-2 ") + 
# sum(data$Index_test == "GAD-2")
# ##
# sum(data$Diagnosis == "GAD") + 
# sum(data$Diagnosis == "GAD ")
##
data_GAD <- filter(data, stringr::str_detect(Diagnosis, "GAD"))
nrow(data_GAD)
##
data_gad2 <- filter(data_GAD, (stringr::str_detect(Index_test, "GAD-2"))) %>% print(n = 200)
saveRDS(object = data_gad2, file = file.path(local_pkg_dir, "temp_old", "Klaus_et_al_data", "real_data_GAD2.RDS"))
##
data_gad7 <- filter(data_GAD, (stringr::str_detect(Index_test, "GAD-7")) ) %>% print(n = 200)
saveRDS(object = data_gad7, file = file.path(local_pkg_dir, "temp_old", "Klaus_et_al_data", "real_data_GAD7.RDS"))
##
data_HADS <- filter(data_GAD, (stringr::str_detect(Index_test, "HADS")) ) %>% print(n = 200)
saveRDS(object = data_HADS, file = file.path(local_pkg_dir, "temp_old", "Klaus_et_al_data", "real_data_HADS.RDS"))
##
data_BAI <- filter(data_GAD, (stringr::str_detect(Index_test, "BAI")) ) %>% print(n = 200)
saveRDS(object = data_BAI, file = file.path(local_pkg_dir, "temp_old", "Klaus_et_al_data", "real_data_BAI.RDS"))
##
data_BAI$Cutoff
##
if (test == "GAD-7") { 
      data_to_use <- data_gad7
      min_score <- 0
      max_score <- 21
      file_name <- "real_x_GAD7.RDS"
} else if (test == "GAD-2") {
      data_to_use <- data_gad2
      min_score <- 0
      max_score <- 6
      file_name <- "real_x_GAD2.RDS"
} else if (test == "HADS") { 
      data_to_use <- data_HADS
      min_score <- 0
      max_score <- 21
      file_name <- "real_x_HADS.RDS"
} else if (test == "BAI") { 
      data_to_use <- data_BAI
      min_score <- 0
      max_score <- 63
      file_name <- "real_x_BAI.RDS"
}
##
# data_gad2
# data_gad7
# data_HADS
# data_BAI
##
##
## ---- Also extract covariate info (NOT for sim study - just for application section and/or Klaus' analysis):
##
cov_data_Klaus <- tibble(read.csv(file.path(Klaus_dir, "DIAQANDI NWMA database final 2025_05_28_descriptive_short.csv")))
cov_data_Klaus %>% print(width = Inf)
##
cov_data_Klaus_2 <- cov_data_Klaus %>% 
  ##
  dplyr::mutate(id = 1:nrow(cov_data_Klaus)) %>%
  ##
  dplyr::rename(Study = Study..,
                study_setting = setting_new,
                Ref_test = ReferenceBrookeLevis,
                DSM_type = DSM.IV.or.DSM.5,
                prev_AAD = Prevalence_AAD,
                prev_GAD = Prevalence_GAD,
                Female_pct = X..female,
                index_test_admin_type = Index.test.admin..mode,
                low_RoB_QUADAS = Low.Overall.RoB.strict...orig.QUADAS.,
                low_RoB_liberal = Low.Overall.RoB...liberal.own.) %>%
  ##
  dplyr::select(Study, 
                Author, 
                Year, 
                Country,
                study_setting, ## CATEGORICAL (3 categories)
                Ref_test,      ## CATEGORICAL (3 categories)
                DSM_type,
                prev_AAD, ## CONTINUOUS
                prev_GAD, ## CONTINUOUS
                Female_pct,
                index_test_admin_type,
                low_RoB_QUADAS,
                low_RoB_liberal)
##
cov_data_Klaus_2 %>% print(n = 250)

# unique(cov_data_Klaus_2$study_setting) ## CATEGORICAL (3 categories + NO MISSING)
# unique(cov_data_Klaus_2$Ref_test)      ## CATEGORICAL (3 categories + NO MISSING)
# ##
# unique(cov_data_Klaus_2$prev_AAD) ## CONTINUOUS (small % missing)
# unique(cov_data_Klaus_2$prev_GAD) ## CONTINUOUS (small % missing)
# ##
# unique(cov_data_Klaus_2$low_RoB_QUADAS)  ## BINARY (+ NO MISSING)
# unique(cov_data_Klaus_2$low_RoB_liberal) ## BINARY (+ NO MISSING)


cov_data_Klaus_3 <- cov_data_Klaus_2 %>% 
  ##
  dplyr::mutate(Ref_test_clean = case_when(Ref_test %in% c("MINI", "MINI ") ~ "MINI",
                                           Ref_test %in% c("SCID", "SCID ") ~ "SCID",
                                           Ref_test %in% c("Structured", "Structured ") ~ "Structured"),
                ##
                low_RoB_QUADAS_clean  = case_when(low_RoB_QUADAS %in% c("No", "No ") ~ "No",
                                                  low_RoB_QUADAS %in% c("Yes", "Yes ") ~ "Yes"),
                ##
                low_RoB_liberal_clean  = case_when(low_RoB_liberal %in% c("No", "No ") ~ "No",
                                                   low_RoB_liberal %in% c("Yes", "Yes ") ~ "Yes"),
                ##
                logit_prev_AAD = qlogis(prev_AAD),
                logit_prev_GAD = qlogis(prev_GAD))
                
# ## Check:
# unique(cov_data_Klaus_3$Ref_test_clean)
# ##
# unique(cov_data_Klaus_3$low_RoB_QUADAS_clean)
# unique(cov_data_Klaus_3$low_RoB_liberal_clean)
# ##
# cov_data_Klaus_3$logit_prev_AAD
# cov_data_Klaus_3$logit_prev_GAD

test_data_list <- list(data_gad2, 
                       data_gad7, 
                       data_HADS, 
                       data_BAI)
##
# continuous_covs <- c("logit_prev_GAD", 
#                      "logit_prev_AAD")
continuous_covs <- c("logit_prev_GAD")
##
# binary_covs <- c(  "low_RoB_QUADAS_clean",
#                    "low_RoB_liberal_clean")
binary_covs <- c("low_RoB_QUADAS_clean")
##
categorical_covs <- c( "study_setting",
                       "Ref_test_clean")


# cov_data = cov_data_Klaus_3
# test_data_list = test_data_list
# test_names = c("GAD-2", "GAD-7", "HADS", "BAI")
# continuous_vars = continuous_covs
# binary_vars = binary_covs
# categorical_vars = categorical_covs
# standardize_continuous = TRUE
# baseline_values = NULL

test_names = c("GAD-2", "GAD-7", "HADS", "BAI")
                     

outs_covs <- MetaOrdDTA::R_fn_prepare_NMA_covariates( cov_data = cov_data_Klaus_3,
                                                      test_data_list = test_data_list,
                                                      test_names = test_names,
                                                      binary_vars = binary_covs,
                                                      continuous_vars = continuous_covs,
                                                      categorical_vars = categorical_covs,
                                                      standardize_continuous = TRUE, ## makes priors for coeffs much easier! (e.g. N(0, 1) or N(0, 2))
                                                      baseline_values = NULL)
# 
# outs_covs$X_nd
# outs_covs$n_covariates_nd
# outs_covs$n_covariates_d
# outs_covs$n_covariates_max


# outs_covs$X_nd[[1]][1:50, ]



outs_covs_subset <- MetaOrdDTA:::create_subsetted_covariates_for_tests(  cov_data = cov_data_Klaus_3,
                                                                              test_data_list = test_data_list,
                                                                              test_names = test_names,
                                                                              continuous_covs = continuous_covs,
                                                                              binary_covs = binary_covs,
                                                                              categorical_covs = categorical_covs,
                                                                              center_cts = TRUE,
                                                                              scale_cts = TRUE)


# # Center your original data before using it
# for (t in 1:4) {
#   # Center continuous variables
#   for (var in covariate_settings$continuous_vars) {
#       col_idx <- which(colnames(outs_covs_subset$X_nd[[t]]) == var)
#       if (length(col_idx) > 0) {
#         outs_covs_subset$X_nd[[t]][, col_idx] <- scale(outs_covs_subset$X_nd[[t]][, col_idx], 
#                                                           center = TRUE, scale = TRUE)[, 1]
#         ##
#         outs_covs_subset$X_d[[t]][, col_idx] <- scale(outs_covs_subset$X_nd[[t]][, col_idx], 
#                                                           center = TRUE, scale = TRUE)[, 1]
#       }
#   }
# }
# 


# outs_covs_subset$X_nd
# 
# 
# 
# 
# 
# unique(data$Diagnosis)
# 
# 
# data$Cut.off.published.
# 
# data$id
# 
# sum(data$Index_test == "GAD-2 ") + 
#   sum(data$Index_test == "GAD-2")
# ##
# sum(data$Diagnosis == "GAD") + 
#   sum(data$Diagnosis == "GAD ")
# 
# 
# data_GAD <- filter(data, stringr::str_detect(Diagnosis, "GAD"))
# nrow(data_GAD)
# 
# data_gad2 <- filter(data_GAD, (stringr::str_detect(Index_test, "GAD-2"))) %>% print(n = 200)
# saveRDS(object = data_gad2, file = file.path(local_pkg_dir, "temp_old", "Klaus_et_al_data", "real_data_GAD2.RDS"))
# ##
# data_gad7 <- filter(data_GAD, (stringr::str_detect(Index_test, "GAD-7")) ) %>% print(n = 200)
# saveRDS(object = data_gad7, file = file.path(local_pkg_dir, "temp_old", "Klaus_et_al_data", "real_data_GAD7.RDS"))
# ##
# data_HADS <- filter(data_GAD, (stringr::str_detect(Index_test, "HADS")) ) %>% print(n = 200)
# saveRDS(object = data_HADS, file = file.path(local_pkg_dir, "temp_old", "Klaus_et_al_data", "real_data_HADS.RDS"))
# ##
# data_BAI <- filter(data_GAD, (stringr::str_detect(Index_test, "BAI")) ) %>% print(n = 200)
# saveRDS(object = data_BAI, file = file.path(local_pkg_dir, "temp_old", "Klaus_et_al_data", "real_data_BAI.RDS"))
# 
# # test <- "GAD-2"
# # test <- "GAD-7"
# # test <- "HADS"




# 
#  
#  
# mean(unique(data_gad2$N)) ; sd(unique(data_gad2$N))
# mean(unique(data_gad7$N)) ; sd(unique(data_gad7$N))
# mean(unique(data_HADS$N)) ; sd(unique(data_HADS$N))
mean(unique(data_BAI$N)) ; sd(unique(data_BAI$N))

# Assuming your data frame is named 'df'
Se_Sp_summary <- compute_Se_Sp_by_cutpoint(data = data_to_use,
                                           min_score = min_score,
                                           max_score = max_score)
print(Se_Sp_summary)




# data_to_use %>% print(n = 32)
data_to_use_2 <- transform_questionnaire_data_by_disease(data = data_to_use,
                                                  min_score = min_score,
                                                  max_score = max_score)
data_to_use_2
##
x_nd <- convert_to_integers_simple(data_to_use_2$nondiseased) 
x_d  <- convert_to_integers_simple(data_to_use_2$diseased)
##
x_nd
x_d
##
saveRDS(object = data_to_use_2, 
        file = file.path(local_pkg_dir, file_name))



 
data_to_use  %>% arrange(Study, Cutoff) %>% print(n = 200)
 
 
# ## SAVE:
# ##
# saveRDS(object = y_tibble, file = "tibble_Baron_et_al.R")
# saveRDS(object = y_tibble_nd, file = "tibble_Baron_et_al_nd.R")
# saveRDS(object = y_tibble_d, file = "tibble_Baron_et_al_d.R")
# 
#  

Klaus_data_list <- list()
##
Klaus_data_list$data_to_use_2 <- list(x_nd,
                                      x_d)
##
# Klaus_data_list$Se_true_gad2 <- Se_true * 100
# Klaus_data_list$Sp_true_gad2 <- Sp_true * 100
##
print(Klaus_data_list)


 


    ##
    ## ---- Add covariate data (OPTIONAL):
    ##
    str(outs_covs)
    Klaus_data_list$n_studies
    
    data_BAI %>% print(n=500)
    
    
    
    cov_data_Klaus_3$Ref_test
    
    baseline_case <- c(1, 
                       0.0, ## centered logit_prev_GAD
                       # 0.0, ## centered logit_prev_AAD
                       0, ## low_RoB_QUADAS_clean
                       # 0, ## low_RoB_liberal_clean
                       0, ## study_setting_2
                       0, ## study_setting_3
                       0, ## Ref_test_clean_SCID
                       0) ## Ref_test_clean_Structured
    
    
    n_covariates_max <- outs_covs_subset$n_covariates_max
    ##
    if (test == "GAD-2") { 
       index <- 1
    } else if (test == "GAD-7") { 
      index <- 2
    } else if (test == "HADS") { 
      index <- 3
    } else if (test == "BAI") {
      index <- 4
    }
    
    covs_given_test <- list( X_nd = outs_covs_subset$X_nd[[index]],
                      X_d  = outs_covs_subset$X_d[[index]],
                      n_covariates_nd = outs_covs_subset$n_covariates_nd[index],
                      n_covariates_d = outs_covs_subset$n_covariates_d[index],
                      n_covariates_max = n_covariates_max,
                      baseline_case_nd = baseline_case,
                      baseline_case_d = baseline_case)
    
    try({ 
      Klaus_data_list$X_nd <- NULL
      Klaus_data_list$X_d <- NULL
      Klaus_data_list$n_covariates_nd <- NULL
      Klaus_data_list$n_covariates_d <- NULL
      Klaus_data_list$n_covariates_max <- NULL
      Klaus_data_list$baseline_case_nd <- NULL
      Klaus_data_list$baseline_case_d <- NULL
    })
    
    Klaus_data_list <- c(Klaus_data_list, covs_given_test)
    
    
    str(Klaus_data_list)
    
    


{
        ##
        Klaus_data_list$n_studies <- nrow(Klaus_data_list$data_to_use_2[[1]])
        Klaus_data_list$n_thr <- ncol(Klaus_data_list$data_to_use_2[[1]]) - 1
        n_thr <- Klaus_data_list$n_thr
        ##
        Klaus_data_list$prior_raw_scale_mu_mean <- 0.0
        Klaus_data_list$prior_raw_scale_mu_SD   <- 1.0
        Klaus_data_list$prior_raw_scale_SD_mean <- 0.0
        Klaus_data_list$prior_raw_scale_SD_SD   <- 1.0
        ##
        Klaus_data_list$prior_dirichlet_alpha <- rep(1.0, Klaus_data_list$n_thr + 1)
        # ##
        # Klaus_data_list$prior_dirichlet_cat_means_alpha <- rep(1.0, Klaus_data_list$n_thr + 1)
        # Klaus_data_list$prior_dirichlet_cat_SDs_mean    <- rep(0.0, Klaus_data_list$n_thr + 1)
        # Klaus_data_list$prior_dirichlet_cat_SDs_SD      <- rep(0.10, Klaus_data_list$n_thr + 1)
        ##
        Klaus_data_list$prior_only <- 0
        Klaus_data_list$kappa_lb <- 0.0
        ##
        Klaus_data_list$softplus <- 0
        ##
        Klaus_data_list$beta_corr_lb <- -1.0
        Klaus_data_list$beta_corr_ub <- +1.0
        Klaus_data_list$prior_beta_corr_LKJ <- 2.0
        ##
        ##
        Klaus_data_list$n_covariates     <- 1
        Klaus_data_list$n_covariates_nd  <- 1
        Klaus_data_list$n_covariates_d   <- 1
        Klaus_data_list$n_covariates_max <- 1
        # ##
        Klaus_data_list$X    <- array(1, dim = c(Klaus_data_list$n_studies, 1))
        Klaus_data_list$X_nd <- array(1, dim = c(Klaus_data_list$n_studies, 1))
        Klaus_data_list$X_d  <- array(1, dim = c(Klaus_data_list$n_studies, 1))
        # ##
        Klaus_data_list$baseline_case <- 1
        Klaus_data_list$baseline_case_nd <- 1
        Klaus_data_list$baseline_case_d  <- 1
        ##
        Klaus_data_list$prior_dirichlet_prob_SDs  <- c(0.05, 0.05)
}
 

 

 

##
## --------------------  Run Stan model to get "true" parameter values for the PHQ-9 test based on real-world data: ---------------------------
##
##
## COMPILE Stan model:
##
{
  
        model_type <- "Jones"
        #  model_type <- "ord_HSROC"
      # model_type <- "ord_bivariate"
      
      if (model_type == "ord_bivariate") { 
        
          stan_model_file_name <-  "DTA_MA_Xu_RANDthr_kappa.stan"
            # stan_model_file_name <-  "DTA_MA_Xu_FIXEDthr.stan"
        
            Klaus_data_list$prior_beta_mu_mean <- array(0.0, dim = c(2, Klaus_data_list$n_covariates_max))
            Klaus_data_list$prior_beta_mu_SD   <- array(1.0, dim = c(2, Klaus_data_list$n_covariates_max))
            ##
            Klaus_data_list$prior_beta_SD_mean <- rep(0.0, 2)
            Klaus_data_list$prior_beta_SD_SD   <- rep(1.0, 2)
            ####
            Klaus_data_list$prior_dirichlet_alpha <- array(1.0, dim = c(2, Klaus_data_list$n_thr + 1))
        
      } else if (model_type == "Jones") { 
        
           stan_model_file_name <-  "DTA_MA_Jones.stan" 
           
           Klaus_data_list$prior_beta_mu_mean <- array(0.0, dim = c(2, Klaus_data_list$n_covariates_max))
           Klaus_data_list$prior_beta_mu_SD   <- array(5.0, dim = c(2, Klaus_data_list$n_covariates_max))
           ##
           Klaus_data_list$prior_beta_SD_mean <- rep(0.0, 2)
           Klaus_data_list$prior_beta_SD_SD   <- rep(1.0, 2)
        
      } else if (model_type == "ord_HSROC") { 
        
           stan_model_file_name <- "DTA_MA_Gat_RANDthr_kappa.stan"
           
           Klaus_data_list$prior_beta_mu_mean <- 0.0
           Klaus_data_list$prior_beta_mu_SD   <- 1.0
           ##
           Klaus_data_list$prior_beta_SD_mean <- 0.0
           Klaus_data_list$prior_beta_SD_SD   <- 1.0
           ##
           Klaus_data_list$prior_dirichlet_alpha <- rep(1.0, Klaus_data_list$n_thr + 1)
           
      }
      
 
       #  
       # # stan_model_file_name <- "DTA_MA_Gat_FIXEDthr.stan"
   
       outs <- MetaOrdDTA:::R_fn_compile_stan_model_basic_given_file_name( 
                                                                           stan_model_file_name =  stan_model_file_name,
                                                                           ##
                                                                           cts        = FALSE, 
                                                                           network    = FALSE,
                                                                           prior_only = FALSE,
                                                                           ##
                                                                           compile = TRUE,
                                                                           debugging = FALSE,
                                                                           force_recompile = FALSE,
                                                                           quiet = FALSE)
      # mod <- outs$mod
      mod <- outs$stan_model_obj

}


##
## Initial values:
##
{
      
      Klaus_init_list <- list()
      ##
      if (model_type == "ord_bivariate") { 
        
            ##
            Klaus_init_list$beta_mu <- array(0.001, dim = c(2, Klaus_data_list$n_covariates_max))
            ##
            Klaus_init_list$beta_SD <- rep(0.001, 2)
            Klaus_init_list$beta_z  <- array(0.001, dim = c(Klaus_data_list$n_studies, 2))
            ##
            Klaus_init_list$C_raw_vec <-   array(-3.0, dim = c(2, Klaus_data_list$n_thr)) ##  Klaus_data_list$n_thr)
            
            ##
            prior_dirichlet_phi <- rep(1.0, n_thr + 1)
            prior_dirichlet_phi[] <- 1.0
            Klaus_data_list$prior_dirichlet_phi <- list(prior_dirichlet_phi, prior_dirichlet_phi)
            ##
            #
            Klaus_data_list$kappa_lb <- 0.001
            Klaus_data_list$prior_kappa_mean <- c(log(100), log(100))
            Klaus_data_list$prior_kappa_SD <- rep(1.0, 2)
        
        
      } else if (model_type == "ord_HSROC")  { 
        
            ##
            Klaus_init_list$beta_mu <- 0.0
            Klaus_init_list$beta_SD <- 0.001
            Klaus_init_list$beta_z <-  rep(0.001, Klaus_data_list$n_studies)
            ##
            Klaus_init_list$raw_scale_mu <-  0.0
            Klaus_init_list$raw_scale_SD <-  0.001
            Klaus_init_list$raw_scale_z  <-  rep(0.001, Klaus_data_list$n_studies)
            ##
            Klaus_init_list$C_raw_vec <-   rep(-3.0,  Klaus_data_list$n_thr)
            ##
            Klaus_data_list$kappa_lb <- 0.001
            Klaus_data_list$prior_kappa_mean <- log(100) ## c(log(200), log(50))
            Klaus_data_list$prior_kappa_SD <- 1.0
            
          
      } else if (model_type == "Jones")  { 
        
            Klaus_init_list$beta_mu <- rep(0.0, 2)
            Klaus_init_list$beta_SD <- rep(0.001, 2)
            Klaus_init_list$beta_z  <- array(0.001, dim = c(Klaus_data_list$n_studies, 2))
            ##
            Klaus_init_list$raw_scale_mu <- rep(0.0, 2)
            Klaus_init_list$raw_scale_SD <- rep(0.001, 2)
            Klaus_init_list$raw_scale_z  <- array(0.001, dim = c(Klaus_data_list$n_studies, 2))
            
            
            Klaus_data_list$cts_thr_values <- seq(from = 1, to = Klaus_data_list$n_thr, by = 1)
            ##
            Klaus_data_list$box_cox <- 1
            Klaus_data_list$prior_boxcox_lambda_mean <- rep(0.0, 1)
            Klaus_data_list$prior_boxcox_lambda_SD   <- rep(1.0, 1)
            ##
            Klaus_data_list$raw_scale_corr_lb <- -1.0
            Klaus_data_list$raw_scale_corr_ub <- +1.0
            Klaus_data_list$prior_raw_scale_corr_LKJ <- 2.0
            ##
            Klaus_data_list$prior_raw_scale_mu_mean <- matrix(rep(0.0, 2))
            Klaus_data_list$prior_raw_scale_mu_SD   <- matrix(rep(1.0, 2))
            Klaus_data_list$prior_raw_scale_SD_mean <- (rep(0.0, 2))
            Klaus_data_list$prior_raw_scale_SD_SD   <- (rep(1.0, 2))
        
      }
  
  
  # Klaus_data_list$compute_sim_study_metrics <- 0
  # Klaus_data_list$n_inner_thr <- 1
  # Klaus_data_list$n_outer_thr <- 1
  # Klaus_data_list$vec_index_inner_thr <- 1
  # Klaus_data_list$vec_index_outer_thr <- 1
  
}


{
    ##
    ## ---- For GAD-2 single-test analysis:
    ##
    Klaus_data_list$x <- Klaus_data_list$data_to_use_2
    ##
    Klaus_data_list$n_obs_cutpoints <- rowSums(Klaus_data_list$x[[1]][, c(-1)] != -1)
    ##
    cutpoint_index <- create_cutpoint_index_matrix(data_matrix = Klaus_data_list$x[[1]], 
                                                                   max_score = Klaus_data_list$n_thr)
    ##
    Klaus_data_list$cutpoint_index <- list(cutpoint_index, 
                                           cutpoint_index)
    
    
    Klaus_data_list
    Klaus_data_list$n_thr
    
    mean(Klaus_data_list$x[[1]][,1] + Klaus_data_list$x[[2]][,1])
    
    sd(Klaus_data_list$x[[1]][,1] + Klaus_data_list$x[[2]][,1])
    
    
    Klaus_data_list$x
    ##
    # Klaus_data_list$x <- MetaOrdDTA:::arrange_missing_values(Klaus_data_list$x)
    # Klaus_data_list$cutpoint_index <- MetaOrdDTA:::arrange_missing_values(Klaus_data_list$cutpoint_index)
    ##
    rearranged <-  MetaOrdDTA::arrange_missing_values_paired(Klaus_data_list$x, Klaus_data_list$cutpoint_index)
    Klaus_data_list$x <- rearranged$x
    Klaus_data_list$cutpoint_index <- rearranged$cutpoint_index
    ##
    Klaus_data_list$use_probit_link <- 1
    
    
    # After creating your list of matrices, convert to 3D array:
    x_array <- array(NA, dim = c(2, Klaus_data_list$n_studies, Klaus_data_list$n_thr + 1))
    for (c in 1:2) {
      x_array[c, , ] <- as.integer(Klaus_data_list$x[[c]])
    }
    Klaus_data_list$x <- x_array
    ##
    # Same for cutpoint_index:
    cutpoint_index_array <- array(NA, dim = c(2, Klaus_data_list$n_studies, Klaus_data_list$n_thr + 1))
    for (c in 1:2) {
      cutpoint_index_array[c, , ] <- as.integer(Klaus_data_list$cutpoint_index[[c]])
    }
    Klaus_data_list$cutpoint_index <- cutpoint_index_array
    
    Klaus_data_list$use_probit_link <- 1
    
    Klaus_data_list$K_fold_CV_indicator <- rep(1, Klaus_data_list$n_studies)
    
 



}



str(Klaus_data_list$x)

## Correction for BAI data:
if (test == "BAI") {
  Klaus_data_list$x[2, 1, 5] <- 23
  # Klaus_data_list$x[2, 1, ]
}


Klaus_data_list$X_nd


Klaus_data_list

Klaus_data_list$prior_beta_mu_SD
Klaus_data_list$prior_beta_SD_SD


Klaus_data_list$box_cox


##
## Sample the model:
##
{
        
        seed <- 123
        
        n_chains <- 8
        init_lists_per_chain <- rep(list(Klaus_init_list), n_chains) 
        
        n_burnin <- 500
        n_iter   <- 1000
        
        Stan_mod_sample <- mod$sample( seed = seed,
                                       data = Klaus_data_list,
                                       init =   0.01,## init_lists_per_chain, 
                                       chains = n_chains,
                                       parallel_chains = n_chains, 
                                       refresh = 50,
                                       iter_sampling = n_iter,
                                       iter_warmup = n_burnin,
                                       max_treedepth = 10, 
                                       adapt_delta = 0.80,
                                       metric = "diag_e")
  
}




# Klaus_data_list$data_to_use_2
# Klaus_data_list$X_nd



##
## Outputs to obtain "true" raw parameter values for simulation study:
##
{
        try({
          Stan_mod_sample$summary(c("beta_mu")) %>% print(n = 100)
          Stan_mod_sample$summary(c("beta_SD")) %>% print(n = 100)
          try({
             Stan_mod_sample$summary(c("beta_Omega", "beta_Sigma")) %>% print(n = 100)
          })
          ##
          try({
          Stan_mod_sample$summary(c("raw_scale_mu")) %>% print(n = 100)
          Stan_mod_sample$summary(c("raw_scale_SD")) %>% print(n = 100)
          })
          ##
          try({
          Stan_mod_sample$summary(c("lambda")) %>% print(n = 100)
          })
          
        })
        try({
           C <- Stan_mod_sample$summary(c("C")) %>% print(n = 100)
          Stan_mod_sample$summary(c("C_MU")) %>% print(n = 100)
        })
        # try({
        #   Stan_mod_sample$summary(c("Se")) %>% print(n = 100)
        #   Stan_mod_sample$summary(c("Sp")) %>% print(n = 100)
        # })
        try({
          Se <- Stan_mod_sample$summary(c("Se_baseline"), quantiles = ~ quantile(., probs = c(0.025, 0.50, 0.975))) %>% print(n = 100)
          # Se_pred <- Stan_mod_sample$summary(c("Se_pred"), quantiles = ~ quantile(., probs = c(0.025, 0.50, 0.975))) %>% print(n = 100)
          ##
          Sp <- Stan_mod_sample$summary(c("Sp_baseline"), quantiles = ~ quantile(., probs = c(0.025, 0.50, 0.975))) %>% print(n = 100)
          # Sp_pred <- Stan_mod_sample$summary(c("Sp_pred"), quantiles = ~ quantile(., probs = c(0.025, 0.50, 0.975))) %>% print(n = 100)
          ##
          Fp <- Stan_mod_sample$summary(c("Fp_baseline"), quantiles = ~ quantile(., probs = c(0.025, 0.50, 0.975))) %>% print(n = 100)
          # Fp_pred <- Stan_mod_sample$summary(c("Fp_pred"), quantiles = ~ quantile(., probs = c(0.025, 0.50, 0.975))) %>% print(n = 100)
        })
  
}




# Se_true_PHQ_9 <- Baron_data_list$Se_true_PHQ_9
# Sp_true_PHQ_9 <- Baron_data_list$Sp_true_PHQ_9

Se_medians <- Se$`50%` * 100
Sp_medians <- Sp$`50%` * 100

# 
# Se_medians
# Sp_medians


## NOTE: ord-bivariate fits data better than the more restricted ord-HSROC model (neglidable
## diff in estimates all very close to 0, but Se estimates of ord-HSROC model are average 
## of 1% off - which is quite notable! however Sp is fine)
## -> Therefore, decided to use ord-bivariate for DGM for PHQ-9 test as it can generate the 
## real data (virtually) perfectly


# dput(C$median)
# 


dirichlet_cat_SDs_sigma <- Stan_mod_sample$summary(c("dirichlet_cat_SDs_sigma")) %>% print(n = 100)
##
C <- Stan_mod_sample$summary(c("C_mu")) %>% print(n = 100)
 C <- Stan_mod_sample$summary(c("C")) %>% print(n = 100)
##
alpha <- Stan_mod_sample$summary(c("alpha")) %>% print(n = 100)
beta_mu <- Stan_mod_sample$summary(c("beta_mu")) %>% print(n = 100)
##
Se <- Stan_mod_sample$summary(c("Se_baseline")) %>% print(n = 100)
Sp <- Stan_mod_sample$summary(c("Sp_baseline")) %>% print(n = 100)
Fp <- Stan_mod_sample$summary(c("Fp_baseline")) %>% print(n = 100)
##
Se_vec <- Sp_vec <- Fp_vec <- c()
##
for (i in 1:Klaus_data_list$n_thr) {
  Se_vec[i] <- Se$median[i]
  Sp_vec[i] <- Sp$median[i]
  Fp_vec[i] <- Fp$median[i]
}

if (model_type == "ord_HSROC")  { 
  
      C_both <- c()
      dirichlet_cat_SDs_sigma_both <- c()
      beta_mu_both <- c()
      ##
      for (i in 1:Klaus_data_list$n_covariates_max) {
        beta_mu_both[i] <- beta_mu$median[i]
      }
      ##
      for (i in 1:Klaus_data_list$n_thr) {
        
        C_both[i]  <- C$median[i]
        ##
        try({ 
          dirichlet_cat_SDs_sigma_both[i] <- dirichlet_cat_SDs_sigma$median[i]
        })
        
      }
      ## 
      alpha_both <- c()
      ##
      for (i in 1:(Klaus_data_list$n_thr + 1)) {
        alpha_both[i]  <- alpha$median[i]
      }
  
} else { 
   
      C_nd <- C_d <- c()
      dirichlet_cat_SDs_sigma_nd <- dirichlet_cat_SDs_sigma_d <- c()
      beta_mu_nd <- beta_mu_d <- c()
      ##
      for (i in 1:Klaus_data_list$n_covariates_max) {
          beta_mu_nd[i] <- beta_mu$median[2*(i - 1) + 1]
          beta_mu_d[i]  <- beta_mu$median[2*i]
      }
      ##
      for (i in 1:Klaus_data_list$n_thr) {
        
            C_nd[i] <- C$median[2*(i - 1) + 1]
            C_d[i]  <- C$median[2*i]
            ##
            try({ 
              dirichlet_cat_SDs_sigma_nd[i] <- dirichlet_cat_SDs_sigma$median[2*(i - 1) + 1]
              dirichlet_cat_SDs_sigma_d[i]  <- dirichlet_cat_SDs_sigma$median[2*i]
            })
      
      }
      ## 
      alpha_nd <- alpha_d <- c()
      ##
      for (i in 1:(Klaus_data_list$n_thr + 1)) {
        alpha_nd[i] <- alpha$median[2*(i - 1) + 1]
        alpha_d[i]  <- alpha$median[2*i]
      }
    
}
##
##
# ##
(dput(round(C_nd, 7)))
(dput(round(C_d, 7)))
# ##
(dput(round(dirichlet_cat_SDs_sigma_nd, 7)))
(dput(round(dirichlet_cat_SDs_sigma_d, 7)))
##
(dput(round(alpha_nd, 7)))
(dput(round(alpha_d, 7)))
##
##
(dput(round(C_both, 7)))
# ##
(dput(round(dirichlet_cat_SDs_sigma_both, 7)))
##
(dput(round(alpha_both, 7)))
##
(dput(round(100*Se_vec, 3)))
(dput(round(100*Fp_vec, 3)))
(dput(round(100*Sp_vec, 3)))
##
##
##
(dput(signif(beta_mu_nd, 3)))
(dput(signif(beta_mu_d, 3)))
##

test ## -------------
model_type
# Klaus_data_list$X_nd
##
##
Stan_mod_sample$summary(c("beta_mu")) %>% print(n = 100)
Stan_mod_sample$summary(c("beta_SD")) %>% print(n = 100)
Stan_mod_sample$summary(c("beta_Omega")) %>% print(n = 100)
##
Stan_mod_sample$summary(c("raw_scale_mu")) %>% print(n = 100)
Stan_mod_sample$summary(c("raw_scale_SD")) %>% print(n = 100)
Stan_mod_sample$summary(c("raw_scale_Omega")) %>% print(n = 100)
##
test
##
# Stan_mod_sample$summary(c("Se_baseline")) %>% print(n = 100)
# Stan_mod_sample$summary(c("Sp_baseline")) %>% print(n = 100)
# Stan_mod_sample$summary(c("Fp_baseline")) %>% print(n = 100)
# 
# 
# Stan_mod_sample$summary(c("C_mu")) %>% print(n = 100)
# Stan_mod_sample$summary(c("C_raw")) %>% print(n = 100)
# 
# 



signif(mean(dirichlet_cat_SDs_sigma_nd), 3)
signif(max(dirichlet_cat_SDs_sigma_nd), 3)
signif(sd(dirichlet_cat_SDs_sigma_nd), 3)
##
signif(mean(dirichlet_cat_SDs_sigma_d), 3)
signif(max(dirichlet_cat_SDs_sigma_d), 3)
signif(sd(dirichlet_cat_SDs_sigma_d), 3)




median(dirichlet_cat_SDs_sigma_nd)
median(dirichlet_cat_SDs_sigma_d)
##
Stan_mod_sample$summary(c("kappa")) %>% print(n = 100)


mean(exp(rnorm(n = 10000, mean = log(100), sd = 1.5)))
median(exp(rnorm(n = 10000, mean = log(100), sd = 1.5)))

quantile(exp(rnorm(n = 10000, mean = log(100), sd = 1.5)), probs = c(0.025, 0.50, 0.975))

dirichlet_cat_SDs_sigma



Stan_mod_sample$summary(c("kappa")) %>% print(n = 100)


Stan_mod_sample$summary(c("Se")) %>% print(n = 100)
Stan_mod_sample$summary(c("Sp")) %>% print(n = 100)

# ##
# ## Store + save "true" parameter values to use for simulation:
# ##
# {
#     true_mean_PHQ_9_Cerullo_Gat_params_list <- list()
#     true_mean_PHQ_9_Cerullo_Gat_params_list$Ord_HSROC <- list()
#     true_mean_PHQ_9_Cerullo_Gat_params_list$Ord_bivariate <- list()
#     ##
#     ## Params for ord-HSROC DGP:
#     ##
#     true_mean_PHQ_9_Cerullo_Gat_params_list$Ord_HSROC$beta_mu <-  0.837
#     true_mean_PHQ_9_Cerullo_Gat_params_list$Ord_HSROC$raw_scale_mu  <- -0.0397 
#     true_mean_PHQ_9_Cerullo_Gat_params_list$Ord_HSROC$C <-  c(-5.057555, -2.012125, -1.54322, -1.161175, -0.8611485, 
#                                                     -0.578003, -0.327348, -0.09665075, 0.09790845, 0.272485, 0.4688355, 
#                                                     0.6720595, 0.8375195, 0.958285, 1.031275, 1.265085, 1.347275, 
#                                                     1.54239, 1.6897, 2.010015, 2.132755, 2.23139, 2.48997, 2.57877, 
#                                                     2.74645, 3.00994, 3.01703)
#     ##
#     ## Params for ord-bivariate DGP:
#     ##
#     true_mean_PHQ_9_Cerullo_Gat_params_list$Ord_bivariate$beta_mu <-  c(-1.72,  -0.159)
#     ##
#     true_mean_PHQ_9_Cerullo_Gat_params_list$Ord_bivariate$C_nd <- c(-2.84831821, -2.3991021125, -2.03067032625, -1.7397899325, 
#                                                                     -1.466111744125, -1.229449821875, -1.011687069, -0.831343162875, 
#                                                                     -0.669146864435, -0.473289126897625, -0.299543294298375, -0.1337316803615, 
#                                                                     -0.01839557547548, 0.078408769757035, 0.344478924181375, 0.399385793905875, 
#                                                                     0.609455893075, 0.843181055, 1.052238411, 1.237164372875, 1.25896475775, 
#                                                                     1.28145171775, 1.30471162275, 1.83039288375, 1.96358992625, 2.139608805, 
#                                                                     2.43692291875)
#     ##
#     true_mean_PHQ_9_Cerullo_Gat_params_list$Ord_bivariate$C_d  <- c(-3.32756737125, -2.42484669625, -2.16246129, -1.99281247125, 
#                                                                    -1.76690405, -1.419330548, -1.139900163625, -0.889595967125, 
#                                                                    -0.6783053727825, -0.514200418774488, -0.235775588161387, -0.0867989636942887, 
#                                                                    0.0385369718472375, 0.065950386820325, 0.248325709728404, 0.35965754351225, 
#                                                                    0.538172728131, 0.63804592431375, 0.994727939625, 1.100288967375, 
#                                                                    1.22190518575, 1.54635508675, 1.673320262875, 1.839718557375, 
#                                                                    2.09151451375, 2.12236033875, 2.15612568625)
#     ##
#     ## SAVE file:
#     ##
#     saveRDS(object = true_mean_PHQ_9_Cerullo_Gat_params_list, file = "true_mean_PHQ_9_Cerullo_Gat_params_list.RDS")
#     
# }
# 
# # true_mean_PHQ_9_Cerullo_Gat_params_list <- readRDS(file = "true_mean_PHQ_9_Cerullo_Gat_params_list.RDS")
# 
# 
#  
# 



















