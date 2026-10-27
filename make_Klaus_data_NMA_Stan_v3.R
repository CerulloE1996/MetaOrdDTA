



rm(list = ls())
.rs.restartR()  # In RStudio 
##
{
        os <- .Platform$OS.type
        ##
        if (os == "unix") {
          user_root_dir <- Sys.getenv("PWD")
        } else if (os == "windows") {
          user_root_dir <- Sys.getenv("USERPROFILE")
        }
        ##
        local_pkg_dir <- file.path(user_root_dir, "Documents/Work/PhD_work/R_packages/MetaOrdDTA")
        ##
        {
            ## First remove any possible package fragments:
            ## Find user_pkg_ins  tall_dir:
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
        ##
        ## ---- ACTUAL (LOCAL) INSTALL:
        ##
        devtools::clean_dll(local_pkg_dir)
        devtools::document(local_pkg_dir)
        ## sroxygen2::roxygenize(local_pkg_dir)
        ##
        ## ---- Install:
        ##
        devtools::clean_dll(local_pkg_dir)
        devtools::install(local_pkg_dir,
                          ## build_vignettes = FALSE,
                          upgrade = "never",
                          quick = TRUE)
        ##
        ## ---- May need to restart R:
        ##
        .rs.restartR() # In RStudio
        ##
        # ?devtools::install
        ##
        require(MetaOrdDTA)
}
##
require(MetaOrdDTA)
require(RcppParallel)
##
## ---- NicoStan and NicoSIMD from the frozen library (the K-fold workers inherit R_LIBS) -----------------------------------------
##
frozen_library <- "/home/enzocerullo/R/frozen_library_of_the_MetaOrdDTA_sim"
Sys.setenv(R_LIBS = frozen_library)
.libPaths(c(frozen_library, .libPaths()))
##
## ---- The sampler: the K-fold is always NicoStan; the main model NicoStan or cmdstanr ("Stan"), for the z-score check -------
##
algorithm_of_the_main_model <- "NicoStan"
# algorithm_of_the_main_model <- "Stan" ## (cmdstanr: the main model only)
##
## (the sim study's NicoStan settings, with auto_chunking = TRUE)
source(file.path("/home/enzocerullo/Documents/Work/PhD_work/R_packages/MetaOrdDTA", "Klaus_sampler_fns.R"))
## require(BayesMVPtest)
##
require(dplyr)
##
#   # n_index_tests <- 8
##

    target_condition <- "AAD"  ; n_index_tests <- 8 ## REAL analysis
    # target_condition <- "GAD"  ; n_index_tests <- 8 ## REAL analysis
    ##
    inform_dir <- FALSE
    ##
    run_single_model <- FALSE
    ##
    # target_condition <- "GAD"  ; n_index_tests <- 4 ## FAKE analysis (only 4/8 tests)
    ##
    # # MR_model_for_real_data <- "MR_model_6"
    # ##
    # # # Model 1 (baseline meta-reg): 2 covariates
    # # intercept + prevalence
    # ##
    # # # Model 2: 3 covariates  
    # # intercept +  SCID + Structured
    # ##
    # # # Model 3: 3 covariates
    # # intercept +  setting_2 + setting_3
    # ##
    # # MR_model_for_real_data <- "MR_model_1"
    # # MR_model_for_real_data <- "MR_model_2"
    # # MR_model_for_real_data <- "MR_model_3"
    # # MR_model_for_real_data <- "MR_model_4"
    # ##
    # ## FULL MR model:
    # ##
    # # MR_model_for_real_data <- "MR_model_5"
    # # MR_model_for_real_data <- "MR_model_5_excl_RoB"  
    # # ##
    # # ## 2-way MR models:
    # # ##
    # MR_model_for_real_data <- "MR_model_6"  ## (intercept + prev_GAD + ref_test)    ### Best model using K-fold, on GAD data w/ 4 tests
    # # MR_model_for_real_data <- "MR_model_7"  ## (intercept + prev_GAD + study_setting)
    # # MR_model_for_real_data <- "MR_model_8"  ## (intercept + prev_GAD + RoB)
    # # MR_model_for_real_data <- "MR_model_9"  ## (intercept + ref_test + study_setting)
    # # MR_model_for_real_data <- "MR_model_10" ## (intercept + ref_test + RoB)
    # # MR_model_for_real_data <- "MR_model_11" ## (intercept + study_setting + RoB)
    # # # ##
    # # # ## 3-way MR models:
    # # # ##
    # # MR_model_for_real_data <- "MR_model_12" ## (intercept + prev_GAD + ref_test + study_setting)
    # MR_model_for_real_data <- "MR_model_13" ## (intercept + prev_GAD + ref_test + RoB)
    # # MR_model_for_real_data <- "MR_model_14" ## (intercept + prev_GAD + study_setting + RoB)
    # # MR_model_for_real_data <- "MR_model_15" ## (intercept + ref_test + study_setting + RoB)
    # # ##
    # # ##
    # ##
    # MR_model_for_real_data_vec_batch <- c("MR_model_1", "MR_model_2",   "MR_model_3")
    # MR_model_for_real_data_vec_batch <- c("MR_model_4", "MR_model_5",   "MR_model_6")
    # MR_model_for_real_data_vec_batch <- c("MR_model_7", "MR_model_8",   "MR_model_9")
    # MR_model_for_real_data_vec_batch <- c("MR_model_10", "MR_model_11", "MR_model_12")
    # MR_model_for_real_data_vec_batch <- c("MR_model_13", "MR_model_14", "MR_model_15") ## 5
    ##
    # MR_model_for_real_data <- c("intercept_only")
    # MR_model_for_real_data_vec_batch <- c("MR_model_1",  "MR_model_2",  "MR_model_3",  "MR_model_4",  "MR_model_5")
    # MR_model_for_real_data_vec_batch <- c("MR_model_6",  "MR_model_7",  "MR_model_8",  "MR_model_9",  "MR_model_10")
    # MR_model_for_real_data_vec_batch <- c("MR_model_11", "MR_model_12", "MR_model_13", "MR_model_14", "MR_model_15")
    ##
    # MR_model_for_real_data_vec_batch <- c("MR_model_1",   "MR_model_2",   "MR_model_3",   "MR_model_4")
    # MR_model_for_real_data_vec_batch <- c("MR_model_5",   "MR_model_6",   "MR_model_7",   "MR_model_8")
    # MR_model_for_real_data_vec_batch <- c("MR_model_9",   "MR_model_10",  "MR_model_11",  "MR_model_12")
    # MR_model_for_real_data_vec_batch <- c("intercept_only", "MR_model_13",  "MR_model_14",  "MR_model_15")
    ##
    # MR_model_for_real_data_vec_batch <- "intercept_only"
    ##
    # # MR_model_for_real_data_vec_batch <- c("MR_model_3",   "MR_model_4")
    # # MR_model_for_real_data_vec_batch <- c( "MR_model_7",   "MR_model_8")
    # # MR_model_for_real_data_vec_batch <- c( "MR_model_11",  "MR_model_12")
    # MR_model_for_real_data_vec_batch <- c( "MR_model_14",  "MR_model_15")
    ##
    fixed_C_indicator <- TRUE ; CS_indicator <- TRUE  ## Model A
    # fixed_C_indicator <- TRUE ; CS_indicator <- FALSE ## Model B
    # fixed_C_indicator <- FALSE ; CS_indicator <- TRUE  ## Model C
    # fixed_C_indicator <- FALSE ; CS_indicator <- FALSE ## Model D
    ##
    # MR_model_for_real_data_vec_batch <- c("MR_model_1",     "MR_model_2",   "MR_model_3",   "MR_model_4",
    #                                       "MR_model_5",     "MR_model_6",   "MR_model_7",   "MR_model_8")
    # MR_model_for_real_data_vec_batch <- c("MR_model_9",     "MR_model_10",  "MR_model_11",  "MR_model_12",
    #                                       "intercept_only", "MR_model_13",  "MR_model_14",  "MR_model_15")
    ##
    MR_model_for_real_data_vec_batch <- "intercept_only"
    #### MR_model_for_real_data_vec_batch <- c("MR_model_4", "MR_model_14")
    # MR_model_for_real_data_vec_batch <- c("MR_model_4", "MR_model_14", "MR_model_7", "intercept_only")
    ##
    # MR_model_for_real_data_vec_batch <- c("MR_model_1",     "MR_model_2",   "MR_model_3",   "MR_model_4",
    #                                       "MR_model_5",     "MR_model_6",   "MR_model_7",   "MR_model_8",
    #                                       "MR_model_9",     "MR_model_10",  "MR_model_11",  "MR_model_12",
    #                                       "intercept_only", "MR_model_13",  "MR_model_14",  "MR_model_15")
    ##
    # MR_model_for_real_data_vec_batch <- "MR_model_11_b"
    # MR_model_for_real_data_vec_batch <- "MR_model_2_b"
    # MR_model_for_real_data_vec_batch <- "MR_model_2"
    # MR_model_for_real_data_vec_batch <- "MR_model_13"
    ##
    # MR_model_for_real_data_vec_batch <- c( "MR_model_4_b", "MR_model_11_b")
    # MR_model_for_real_data_vec_batch <- c( "MR_model_9_b", "MR_model_2_b")
    # MR_model_for_real_data_vec_batch <- c( "MR_model_7_b", "MR_model_7_c")
    # ##
    # MR_model_for_real_data_vec_batch <- c("MR_model_6",  "MR_model_7",  "MR_model_8",  "MR_model_9",  "MR_model_10")
    # MR_model_for_real_data_vec_batch <- c("MR_model_11", "MR_model_12", "MR_model_13", "MR_model_14", "MR_model_15")
    # ##
    # MR_model_for_real_data_vec_batch <- c("MR_model_3", "MR_model_11", 
    #                                       "MR_model_9", "MR_model_1",
    #                                       "MR_model_4", "MR_model_7")
    ##
    # MR_model_for_real_data_vec_batch <- c("MR_model_15")
    ##
    # MR_model_for_real_data <- "MR_model_5"
    ##
    for (MR_model_for_real_data in MR_model_for_real_data_vec_batch) {
    
{
  
    seed <- 1

    {
        os <- .Platform$OS.type
        
        if (os == "unix") {
          user_root_dir <- Sys.getenv("PWD")
        } else if (os == "windows") {
          user_root_dir <- Sys.getenv("USERPROFILE")
        }
        local_pkg_dir <- file.path(user_root_dir, "Documents/Work/PhD_work/R_packages/MetaOrdDTA")
        ##
        setwd(local_pkg_dir)
    }
    
    Klaus_dir <-  "/home/enzocerullo/Documents/Work/PhD_work/Klaus_collab_NMA_work"
    
    source(file.path(local_pkg_dir, "Klaus_et_al_collab_R_fns.R"))
    
    ## read in data
    # data_Klaus <- tibble(read.csv(file.path(Klaus_dir, "DIAQANDI NWMA database final 2025_04_23.csv"))) ## older data (used this for sim study!)
    
    ######## data_Klaus <- tibble(read.csv(file.path(Klaus_dir, "DIAQANDI NWMA database final 2025_05_28.csv")))
    data_Klaus <- tibble(read.csv(file.path(Klaus_dir, "DIAQANDI NWMA database updated 2026_05_11_2x2_tables.csv")))
    data_Klaus
    
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
     
    unique(data$Diagnosis)
    
    {
        data_GAD <- filter(data, stringr::str_detect(Diagnosis, "GAD"))
        nrow(data_GAD)
        saveRDS(object = data_GAD, file = file.path(local_pkg_dir, "temp_old", "Klaus_et_al_data", "real_data_GAD.RDS"))
        ##
        data_AAD <- filter(data, stringr::str_detect(Diagnosis, "AAD"))
        nrow(data_AAD)
        saveRDS(object = data_AAD, file = file.path(local_pkg_dir, "temp_old", "Klaus_et_al_data", "real_data_AAD.RDS"))
        ##
        if (target_condition == "AAD") { 
          data_to_use <- data_AAD
        } else if (target_condition == "GAD") { 
          data_to_use <- data_GAD
        }
        ##
        ## ----data:
        ##
        data_gad2 <- filter(data_to_use, (stringr::str_detect(Index_test, "GAD-2"))) %>% print(n = 100)
        # saveRDS(object = data_gad2, file = file.path(local_pkg_dir, "temp_old", "Klaus_et_al_data", "real_data_GAD2.RDS"))
        ##
        data_gad7 <- filter(data_to_use, (stringr::str_detect(Index_test, "GAD-7")) ) %>% print(n = 100)
        # saveRDS(object = data_GAD_gad7, file = file.path(local_pkg_dir, "temp_old", "Klaus_et_al_data", "real_data_GAD7.RDS"))
        ##
        data_HADS <- filter(data_to_use, (stringr::str_detect(Index_test, "HADS")) ) %>% print(n = 100)
        # saveRDS(object = data_GAD_HADS, file = file.path(local_pkg_dir, "temp_old", "Klaus_et_al_data", "real_data_HADS.RDS"))
        ##
        data_BAI <- filter(data_to_use, (stringr::str_detect(Index_test, "BAI")) ) %>% print(n = 100)
        # saveRDS(object = data_BAI, file = file.path(local_pkg_dir, "temp_old", "Klaus_et_al_data", "real_data_BAI.RDS"))
        
        # # Remove Eack 2008 - no covariate data available
        # data_BAI <- data_BAI %>%
        #   filter(!(Author == "Eack" & Year == 2008))
        
        ##
        data_OASIS <- filter(data_to_use, (stringr::str_detect(Index_test, "OASIS")) ) %>% print(n = 100)
        # saveRDS(object = data_OASIS, file = file.path(local_pkg_dir, "temp_old", "Klaus_et_al_data", "real_data_OASIS.RDS"))
        ##
        data_STAI_S <- filter(data_to_use, (stringr::str_detect(Index_test, "STAI-S")) ) %>% print(n = 100)
        data_STAI_S_adj <- data_STAI_S %>% 
          dplyr::mutate(Cutoff = Cutoff - 20)
        # saveRDS(object = data_STAI_S_adj, file = file.path(local_pkg_dir, "temp_old", "Klaus_et_al_data", "real_data_STAI_S.RDS"))
        ##
        data_STAI_T <- filter(data_to_use, (stringr::str_detect(Index_test, "STAI-T")) ) %>% print(n = 100)
        data_STAI_T_adj <- data_STAI_T %>% 
          dplyr::mutate(Cutoff = Cutoff - 20)
        # saveRDS(object = data_STAI_T_adj, file = file.path(local_pkg_dir, "temp_old", "Klaus_et_al_data", "real_data_STAI_T.RDS"))
        ##
        data_PROMIS <- filter(data_to_use, (stringr::str_detect(Index_test, "PROMIS")) ) %>% print(n = 100)
        data_PROMIS_adj <- data_PROMIS %>% 
          dplyr::mutate(Cutoff = Cutoff - 8)
        # saveRDS(object = data_GAD_PROMIS_adj, file = file.path(local_pkg_dir, "temp_old", "Klaus_et_al_data", "real_data_PROMIS.RDS"))
        ##
        if (n_index_tests == 4) {
            test_names = c("GAD-2", "GAD-7", "HADS", "BAI")
            ##
            test_data_list <- list(data_gad2, 
                                   data_gad7, 
                                   data_HADS, 
                                   data_BAI)
            ##
            n_thr = c(6, 21, 21, 63)
        } else { 
            test_names = c("GAD-2", "GAD-7", "HADS", "BAI",
                           "OASIS", "STAI_S", "STAI_T", "PROMIS")
            ##
            test_data_list <- list(data_gad2, 
                                   data_gad7, 
                                   data_HADS, 
                                   data_BAI,
                                   data_OASIS,
                                   data_STAI_S_adj,
                                   data_STAI_T_adj,
                                   data_PROMIS_adj)
            ##
            n_thr = c(6, 21, 21, 63, 20, 60, 60, 32)
        }
   
    }
     


}
    
    
    ## fix prevs:
    # test_data_list <- fix_prevalences(test_data_list)
      
      
      
      # data_gad2 %>% print(n = 50)
    
      # y_list = test_data_list
      # n_studies = 8
      # n_tests_inc_ref = n_index_tests + 1
      # 
      # convert_to_aggregate_counts( y_list = test_data_list,
      #                              n_studies = 8,
      #                              n_tests_inc_ref = n_index_tests + 1,
      #                              n_thr = c(6, 21, 21, 63, 20, 60, 60, 32))
      
      
      # 
      # # Use it:
      # summarise_crude_accuracy(data_gad2, "GAD-2")
      # target_condition
      # 
      # # Or for all tests at once:
      # crude_summaries <- map2_dfr(
      #   test_data_list, 
      #   test_names,
      #   ~summarise_crude_accuracy(.x, .y)
      # )
      # 
      # # View specific cutoffs (e.g., common published cutoffs)
      # crude_summaries %>% 
      #   filter(Cutoff %in% c(2, 3, 10)) %>%
      #   print(n = 50)
      # 
      # 
      # 
      # require(tidyr)
      # require(ggplot2)
      # require(purrr)
      # 
      # # Usage:
      # plot_score_distribution(data_gad2, "GAD-2")
      # plot_all_distributions(test_data_list, test_names, ncol = 4)
      # 
      
      
      
      

  
# quantile(data_GAD$prev, probs = c(0.05, 0.25, 0.50, 0.75, 0.95))
# 


{
    ##
    ## ---- Also extract covariate info (NOT for sim study - just for application section and/or Klaus' analysis):
    ##
    # cov_data_Klaus <- tibble(read.csv(file.path(Klaus_dir, "DIAQANDI NWMA database final 2025_05_28_descriptive_short.csv")))
    ## DIAQANDI NWMA database updated 2025_11_25_descriptive_short.csv
    cov_data_Klaus <- tibble(read.csv(file.path(Klaus_dir, "DIAQANDI NWMA database updated 2026_05_11_descriptive_short.csv")))
    cov_data_Klaus %>% print(width = Inf)
    ##
    unique(cov_data_Klaus$ReferenceBrookeLevis)
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
    
    cov_data_Klaus_2 %>% print(n = 250)
    
    # unique(cov_data_Klaus_2$study_setting) ## CATEGORICAL (3 categories + NO MISSING)
    # unique(cov_data_Klaus_2$Ref_test)      ## CATEGORICAL (3 categories + NO MISSING)
    # ##
    # unique(cov_data_Klaus_2$prev_AAD) ## CONTINUOUS (small % missing)
    # unique(cov_data_Klaus_2$prev_GAD) ## CONTINUOUS (small % missing)
    # ##
    # unique(cov_data_Klaus_2$low_RoB_QUADAS)  ## BINARY (+ NO MISSING)
    # unique(cov_data_Klaus_2$low_RoB_liberal) ## BINARY (+ NO MISSING)
    {
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
                        
          ## Check:
          unique(cov_data_Klaus_3$Ref_test)
          unique(cov_data_Klaus_3$Ref_test_clean)
          ##
          unique(cov_data_Klaus_3$low_RoB_QUADAS)
          unique(cov_data_Klaus_3$low_RoB_QUADAS_clean)
          ##
          unique(cov_data_Klaus_3$low_RoB_liberal)
          unique(cov_data_Klaus_3$low_RoB_liberal_clean)
          ##
          cov_data_Klaus_3$logit_prev_AAD
          cov_data_Klaus_3$logit_prev_GAD
          ##
          require(MetaOrdDTA)
    }
    
    cov_data_Klaus_3 <- cov_data_Klaus_3 %>%
      mutate(Study = row_number())
    
    ## Remove duplicated studies (5 studies appear twice in descriptive file)
    cov_data_Klaus_3 <- cov_data_Klaus_3 %>%
      filter(!duplicated(paste(trimws(Author), Year)))
    
    # Re-assign sequential Study IDs
    cov_data_Klaus_3 <- cov_data_Klaus_3 %>%
      mutate(Study = row_number())

    
    

}

 

     source(file.path(getwd(), "Klaus_fns.R"))
    
    
    {
 
        
        
        # CRITICAL: Remove old (wrong) Study IDs before fuzzy matching
        for (i in 1:length(test_data_list)) {
          if ("Study" %in% names(test_data_list[[i]])) {
            test_data_list[[i]] <- test_data_list[[i]] %>% select(-Study)
          }
        }
        
        # 2. Match Study IDs (now using corrected cov_data_Klaus_3)
        test_data_list <- Klaus_fn_create_study_ids_fuzzy(test_data_list,
                                                          cov_data_Klaus_3)
        test_data_list <- Klaus_fn_fix_prevalences(test_data_list)
        
        
        
        
        
        # Should be 134
        nrow(cov_data_Klaus_3)
        
        # After fuzzy matching, should be 0 missing
        all_test_studies <- unique(unlist(lapply(test_data_list, function(x) unique(x$Study))))
        missing <- setdiff(all_test_studies, cov_data_Klaus_3$Study)
        cat("Missing:", length(missing), "\n")
        
        
        
        # Remove duplicate rows in test_data_list
        for (i in 1:length(test_data_list)) {
          test_data_list[[i]] <- test_data_list[[i]] %>% distinct()
        }
        
        
        
        
        
        # Check for duplicates in matched data
        test_data_list[[1]] %>%
          select(Author, Year) %>%
          distinct() %>%
          group_by(Author, Year) %>%
          filter(n() > 1)
        
        # Check cov_data_Klaus_3 for duplicates
        cov_data_Klaus_3 %>%
          group_by(Author, Year) %>%
          filter(n() > 1) %>%
          select(Study, Author, Year)
        
        
        
        # Simple check - do Study IDs match to correct authors?
        for (t in 1:n_index_tests) {
          check <- test_data_list[[t]] %>%
            select(Study, Author, Year) %>%
            distinct() %>%
            left_join(cov_data_Klaus_3 %>% select(Study, Author_cov = Author), by = "Study")
          
          # Just print first few to eyeball
          cat("\n=== Test", t, "===\n")
          print(check %>% select(Study, Author, Author_cov) %>% head(5))
        }
        
        
        
        
        cov_data_Klaus_3 %>%
          filter(Author %in% c("Chan ", "Gitaari", "Hébert", "Liu C")) %>%
          select(Study, Author, Year, study_setting, Ref_test_clean, prev_GAD, prev_AAD) %>%
          arrange(Author, Year)
        
        
        # Is Liu C in test_data_list?
        for (t in 1:n_index_tests) {
          liu <- test_data_list[[t]] %>% filter(grepl("Liu C", Author))
          if (nrow(liu) > 0) {
            cat("Test", t, "- Liu C found:\n")
            print(liu %>% select(Study, Author, Year) %>% distinct())
          }
        }
        
        
        # What's the actual test data for both?
        test_data_list[[2]] %>%
          filter(grepl("Liu C", Author)) %>%
          select(Study, Author, Year, Cutoff, N, Se, Sp) %>%
          arrange(Study, Cutoff)
        
        
        # # In the data cleaning script, BEFORE Klaus_fn_create_study_ids_fuzzy,
        # # remove the definite duplicates from cov_data_Klaus_3:
        # cov_data_Klaus_3 <- cov_data_Klaus_3 %>%
        #   filter(!Study %in% c(134, 135, 136))  # Keep 137 for now until confirmed
        # 
        # # Re-assign sequential Study IDs
        # cov_data_Klaus_3 <- cov_data_Klaus_3 %>%
        #   mutate(Study = row_number())
        
        cat("After removing duplicates:", nrow(cov_data_Klaus_3), "studies\n")
        
        
        # Who's missing?
        138 - nrow(cov_data_Klaus_3)  # Should tell you how many removed
        
        
        
        # # After you find out which is correct (e.g., if MINI is correct):
        # cov_data_Klaus_3 <- cov_data_Klaus_3 %>%
        #   mutate(Ref_test_clean = ifelse(Author == "Liu C" & Year == 2025, "MINI", Ref_test_clean))
        
        
        
        # How many now?
        nrow(cov_data_Klaus_3)
        
        # What Study IDs exist?
        range(cov_data_Klaus_3$Study)
        
        # Any gaps in Study IDs?
        setdiff(1:134, cov_data_Klaus_3$Study)
        
        
        # Get all Study IDs used in test_data_list
        all_test_studies <- unique(unlist(lapply(test_data_list, function(x) unique(x$Study))))
        
        # Are any missing from cov_data_Klaus_3?
        missing <- setdiff(all_test_studies, cov_data_Klaus_3$Study)
        
        cat("Studies in tests:", length(all_test_studies), "\n")
        cat("Missing from covariates:", length(missing), "\n")
        if (length(missing) > 0) print(missing)
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

        
        # # Is Eack in any test?
        # for (t in 1:n_index_tests) {
        #   eack <- test_data_list[[t]] %>% filter(grepl("Eack", Author))
        #   if (nrow(eack) > 0) {
        #     cat("Test", t, "- Eack found\n")
        #   }
        # }
        # 
        # 
        # cov_data_Klaus_3 %>% 
        #   filter(grepl("Eack", Author)) %>%
        #   select(Study, Author, Year, prev_GAD, prev_AAD, study_setting, Ref_test_clean)
        
        
        
        
        
        
        
        
        
        
        
        # 3. Create derived objects
        test_prev <- Klaus_fn_get_test_prevalences(test_data_list)
        cov_data_per_test <- Klaus_fn_update_covariates(cov_data_Klaus_3, 
                                                        test_data_list)
        
        # # Rebuild cov_data_per_test using Study IDs (not Author/Year)
        # cov_data_per_test <- list()
        # for (i in 1:length(test_data_list)) {
        #   study_ids <- test_data_list[[i]] %>%
        #     filter(!is.na(Study)) %>%
        #     pull(Study) %>%
        #     unique()
        #   
        #   cov_data_per_test[[i]] <- cov_data_Klaus_3 %>%
        #     filter(Study %in% study_ids)
        # }
        
        # Rebuild cov_data_per_test using Study IDs AND compute prevalence from 2x2 tables
        cov_data_per_test <- list()
        for (i in 1:length(test_data_list)) {
          study_ids <- test_data_list[[i]] %>%
            filter(!is.na(Study)) %>%
            pull(Study) %>%
            unique()
          
          # Get prevalence from 2x2 tables (ALWAYS available)
          prev_from_test <- test_data_list[[i]] %>%
            group_by(Study) %>%
            summarise(prev_from_2x2 = first(prev), .groups = 'drop') %>%
            mutate(logit_prev_GAD = log(prev_from_2x2 / (1 - prev_from_2x2)))
          
          cov_data_per_test[[i]] <- cov_data_Klaus_3 %>%
            filter(Study %in% study_ids) %>%
            select(-prev_GAD, -logit_prev_GAD, -prev_AAD, -logit_prev_AAD) %>%  # Remove file-based prev
            left_join(prev_from_test %>% select(Study, logit_prev_GAD), by = "Study")
        }
        
        
        # 4. Create indicator FROM STUDY IDs (not Author/Year matching)
        indicator_index_test_in_study <- Klaus_fn_create_indicator_from_study_ids(cov_data_Klaus_3, test_data_list)
    }
  
 
##
source(file.path(getwd(), "R", "R_fn_covariates.R"))
##
if (n_index_tests == 4) {
  test_names <- c("GAD-2", "GAD-7", "HADS", "BAI")
  max_scores <- c(6, 21, 21, 63)
} else {
  test_names <- c("GAD-2", "GAD-7", "HADS", "BAI", 
                  "OASIS", "STAI_S", "STAI_T", "PROMIS")
  max_scores <- c(6, 21, 21, 63, 
                  20, 60, 60, 32)
}
 
# # ============================================================
# # CRITICAL: Filter test_data_list to only studies with covariate data
# # ============================================================
# cat("\n=== Filtering to studies with covariate data ===\n")
# for (i in 1:length(test_data_list)) {
#   before_n <- n_distinct(test_data_list[[i]]$Study)
#   valid_studies <- cov_data_per_test[[i]]$Study
#   test_data_list[[i]] <- test_data_list[[i]] %>%
#     filter(Study %in% valid_studies)
#   after_n <- n_distinct(test_data_list[[i]]$Study)
#   
#   if (before_n != after_n) {
#     cat("Test", i, ": FILTERED", before_n, "→", after_n, "studies\n")
#   } else {
#     cat("Test", i, ": OK (", after_n, "studies)\n")
#   }
# }
##
# NOW prepare the NMA data
NMA_data_list <- prepare_NMA_data(
  test_data_list = test_data_list,
  test_names = test_names,
  max_scores = max_scores)
##
indicator_index_test_in_study <- NMA_data_list$indicator_index_test_in_study

# ##
# ## correct data entry errors (PERSONAL CODE):
# if (target_condition == "GAD") {
#       {
#           ##
#           ## ---- Correction for data entry error:
#           ##
#           for (t in 1:4) { 
#             for (c in 1:2) {
#               for (s in 1:76) {
#                 if (any(NMA_data_list$x[[t]][[c]][s, ] == 33)) {
#                     print(paste("t = ", t))
#                     print(paste("c = ", c))
#                     print(paste("s = ", s))
#                     ##
#                     print(NMA_data_list$x[[t]][[c]][s, ])
#                 }
#               }
#             }
#           }
#           ## at: t = 4, c = 2, s = 1 + position/threshold 20:
#           NMA_data_list$x[[4]][[2]][1, 20]
#           ## Replace with 23:
#           NMA_data_list$x[[4]][[2]][1, 20] <- 23
#           ##
#           NMA_data_list$x[[4]][[2]][1, 20]
#       }
# } else { 
#   
#   ##
#   ## ---- Correction for data entry error:
#   ##
#   for (t in 1:8) { 
#     for (c in 1:2) {
#       for (s in 1:80) {
#         if (any(NMA_data_list$x[[t]][[c]][s, ] == 87)) {
#               print(paste("t = ", t))
#               print(paste("c = ", c))
#               print(paste("s = ", s))
#               ##
#               print(NMA_data_list$x[[t]][[c]][s, ])
#         }
#       }
#     }
#   }
# }

{
  ##
  ## ---- MR-Model-1 (intercept + prev_GAD) - 2 covariates total
  ##
  continuous_covs <- c("logit_prev_GAD")
  binary_covs <- c()
  categorical_covs <- c()
  ##
  outs_covs_subset_MR_model_1 <-  MetaOrdDTA:::create_subsetted_covariates_for_tests(    cov_data = cov_data_per_test,
                                                                                         indicator_index_test_in_study = indicator_index_test_in_study,
                                                                                         test_names = test_names,
                                                                                         continuous_covs = continuous_covs,
                                                                                         binary_covs = binary_covs,
                                                                                         categorical_covs = categorical_covs,
                                                                                         center_and_scale_cts = TRUE)
  ##
  outs_covs_subset_MR_model_1$X_nd
  outs_covs_subset_MR_model_1$X_d
  outs_covs_subset_MR_model_1$continuous_covs
  ##
  ## ---- MR-Model-2 (intercept + reference_test) - 4 covariates total (inc. bin ind.)
  ##
  continuous_covs <- c()
  binary_covs <- c()
  categorical_covs <- c("Ref_test_clean")
  ##
  outs_covs_subset_MR_model_2 <-  create_subsetted_covariates_for_tests(    cov_data = cov_data_per_test,
                                                                            indicator_index_test_in_study = indicator_index_test_in_study,
                                                                            test_names = test_names,
                                                                            continuous_covs = continuous_covs,
                                                                            binary_covs = binary_covs,
                                                                            categorical_covs = categorical_covs,
                                                                            center_and_scale_cts = TRUE)
  ##
  outs_covs_subset_MR_model_2$X_nd
  outs_covs_subset_MR_model_2$X_d
  ##
  ## ---- MR-Model-3 (intercept + study_setting) - 4 covariates total (inc. bin ind.)
  ##
  continuous_covs <- c()
  binary_covs <- c()
  categorical_covs <- c("study_setting")
  ##
  outs_covs_subset_MR_model_3 <-  create_subsetted_covariates_for_tests(    cov_data = cov_data_per_test,
                                                                            indicator_index_test_in_study = indicator_index_test_in_study,
                                                                            test_names = test_names,
                                                                            continuous_covs = continuous_covs,
                                                                            binary_covs = binary_covs,
                                                                            categorical_covs = categorical_covs,
                                                                            center_and_scale_cts = TRUE)
  ##
  outs_covs_subset_MR_model_3$X_nd
  outs_covs_subset_MR_model_3$X_d
  ##
  ## ---- MR-Model-4 (intercept + prev_GAD + reference_test + study_setting) - 6 covariates total (inc. bin ind.)
  ##
  continuous_covs <- c()
  binary_covs <- c("low_RoB_QUADAS_clean")
  categorical_covs <- c()
  ##
  outs_covs_subset_MR_model_4 <-  create_subsetted_covariates_for_tests(    cov_data = cov_data_per_test,
                                                                            indicator_index_test_in_study = indicator_index_test_in_study,
                                                                            test_names = test_names,
                                                                            continuous_covs = continuous_covs,
                                                                            binary_covs = binary_covs,
                                                                            categorical_covs = categorical_covs,
                                                                            center_and_scale_cts = TRUE)
  ##
  outs_covs_subset_MR_model_4$X_nd
  outs_covs_subset_MR_model_4$X_d
  ##
  ## ---- MR-Model-5 (intercept + prev_GAD + reference_test + study_setting + RoB_QUADAS) - 7 covariates total (inc. bin ind.)
  ##
  continuous_covs <- c("logit_prev_GAD")
  binary_covs <- c("low_RoB_QUADAS_clean")
  categorical_covs <- c("study_setting", 
                        "Ref_test_clean")
  ##
  outs_covs_subset_MR_model_5 <-  create_subsetted_covariates_for_tests(    cov_data = cov_data_per_test,
                                                                            indicator_index_test_in_study = indicator_index_test_in_study,
                                                                            test_names = test_names,
                                                                            continuous_covs = continuous_covs,
                                                                            binary_covs = binary_covs,
                                                                            categorical_covs = categorical_covs,
                                                                            center_and_scale_cts = TRUE)
  ##
  outs_covs_subset_MR_model_5$X_nd
  outs_covs_subset_MR_model_5$X_d
  ##
  ## ---- MR-Model-6 (intercept + prev_GAD + reference_test) - 3 covariates total
  ##
  continuous_covs <- c("logit_prev_GAD")
  binary_covs <- c()
  categorical_covs <- c("Ref_test_clean")
  ##
  outs_covs_subset_MR_model_6 <-  create_subsetted_covariates_for_tests(    cov_data = cov_data_per_test,
                                                                            indicator_index_test_in_study = indicator_index_test_in_study,
                                                                            test_names = test_names,
                                                                            continuous_covs = continuous_covs,
                                                                            binary_covs = binary_covs,
                                                                            categorical_covs = categorical_covs,
                                                                            center_and_scale_cts = TRUE)
  ##
  outs_covs_subset_MR_model_6$X_nd
  outs_covs_subset_MR_model_6$X_d
  ##
  ## ---- Save to list:
  ##
  outs_covs <- list("MR_model_1" = outs_covs_subset_MR_model_1,
                    "MR_model_2" = outs_covs_subset_MR_model_2,
                    "MR_model_3" = outs_covs_subset_MR_model_3,
                    "MR_model_4" = outs_covs_subset_MR_model_4,
                    "MR_model_5" = outs_covs_subset_MR_model_5,
                    "MR_model_6" = outs_covs_subset_MR_model_6)
  ##
  if (n_index_tests == 4) {
    saveRDS(object = outs_covs, 
            file = file.path(local_pkg_dir, "temp_old", "Klaus_et_al_data", 
                             paste0("real_data_cov_data_4_tests", "_", target_condition)))
  } else { 
    saveRDS(object = outs_covs, 
            file = file.path(local_pkg_dir, "temp_old", "Klaus_et_al_data", paste0("real_data_cov_data_8_tests", "_", target_condition)))
  }
  
}





##
## ---------------------------------------------------------------------------------------------------------------------------------------------------
##
outs_covs_subset <- outs_covs$MR_model_5
X <- list(outs_covs_subset$X_nd, outs_covs_subset$X_d)
##
cov_data <- NULL
##
# MR_model_for_real_data_vec_batch <- c( "MR_model_4_b")
# MR_model_for_real_data_vec_batch <- c( "MR_model_11_b")
# MR_model_for_real_data_vec_batch <- c( "MR_model_9_b")
# MR_model_for_real_data_vec_batch <- c( "MR_model_2_b")
# MR_model_for_real_data_vec_batch <- c( "MR_model_7_b")
# MR_model_for_real_data_vec_batch <- c( "MR_model_7_c")
##
{
        ##
        if (MR_model_for_real_data == "intercept_only") {
          
          print(paste("0"))
          intercept_only <- TRUE
          X <- subset_X_covariates( X = X,
                                    covariates_to_keep = c("intercept"))
          baseline_case <- c(1)
          
        } else if (MR_model_for_real_data %in% c("MR_model_1", "MR_model_1_b")) {
          
          print(paste("1"))
          intercept_only <- FALSE
          X <- subset_X_covariates( X = X, 
                                    covariates_to_keep = c("intercept", 
                                                           "logit_prev_GAD"))
          baseline_case <- c(1, 
                             0.0)
          
        } else if (MR_model_for_real_data %in% c("MR_model_2", "MR_model_2_b")) {
          
          print(paste("2"))
          intercept_only <- FALSE
          X <- subset_X_covariates( X = X, 
                                    covariates_to_keep = c("intercept", 
                                                           "Ref_test_clean_SCID", "Ref_test_clean_Structured"))
          baseline_case <- c(1, 
                             0, 0)
          
        } else if (MR_model_for_real_data %in% c("MR_model_3", "MR_model_3_b")) {
          
          print(paste("3"))
          intercept_only <- FALSE
          X <- subset_X_covariates( X = X, 
                                    covariates_to_keep = c("intercept", 
                                                           "study_setting_2", "study_setting_3"))
          baseline_case <- c(1, 
                             0, 1)
          
        } else if (MR_model_for_real_data %in% c("MR_model_4", "MR_model_4_b")) {
          
          print(paste("4"))
          intercept_only <- FALSE
          X <- subset_X_covariates( X = X,
                                    covariates_to_keep = c("intercept",
                                                           "low_RoB_QUADAS_clean"))
          baseline_case <- c(1, 
                             0)
          
        } else if (MR_model_for_real_data %in% c("MR_model_5", "MR_model_5_b")) {
          
          print(paste("5"))
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
                             0, 0,
                             0, 1)
          
        # } else if (MR_model_for_real_data == "MR_model_5_excl_RoB") { 
        #   
        #   intercept_only <- FALSE
        #   X <- subset_X_covariates( X = X, 
        #                             covariates_to_keep = c("intercept", 
        #                                                    "logit_prev_GAD",
        #                                                    "study_setting_2", "study_setting_3",
        #                                                    "Ref_test_clean_SCID", "Ref_test_clean_Structured"))
        #   baseline_case <- c(1,
        #                      0.0, 
        #                      0, 0,
        #                      0, 1)
          
        } else if (MR_model_for_real_data %in% c("MR_model_6", "MR_model_6_b")) { ## (intercept + logit_prev_GAD + ref_test)
          
          print(paste("6"))
          intercept_only <- FALSE
          X <- subset_X_covariates( X = X,
                                    covariates_to_keep = c("intercept",
                                                           "logit_prev_GAD",
                                                           "Ref_test_clean_SCID", "Ref_test_clean_Structured"))
          
          baseline_case <- c(1, 
                             0.0,
                             0, 0)
          
        } else if (MR_model_for_real_data %in% c("MR_model_7", "MR_model_7_b", "MR_model_7_c")) { ## (intercept + logit_prev_GAD + study_setting)
          
          print(paste("7"))
          intercept_only <- FALSE
          X <- subset_X_covariates( X = X,
                                    covariates_to_keep = c("intercept",
                                                           "logit_prev_GAD",
                                                           "study_setting_2", "study_setting_3"))
          baseline_case <- c(1,
                             0.0, 
                             0, 1)
          
        } else if (MR_model_for_real_data %in% c("MR_model_8", "MR_model_8_b")) {  ## (intercept + logit_prev_GAD + RoB)
          
          print(paste("8"))
          intercept_only <- FALSE
          X <- subset_X_covariates( X = X,
                                    covariates_to_keep = c("intercept",
                                                           "logit_prev_GAD",
                                                           "low_RoB_QUADAS_clean"))
          baseline_case <- c(1, 
                             0.0, 
                             0)
          
        } else if (MR_model_for_real_data %in% c("MR_model_9", "MR_model_9_b")) { ## (intercept + ref_test + study_setting)
          
          print(paste("9"))
          intercept_only <- FALSE
          X <- subset_X_covariates( X = X,
                                    covariates_to_keep = c("intercept",
                                                           "study_setting_2", "study_setting_3",
                                                           "Ref_test_clean_SCID", "Ref_test_clean_Structured"))
          baseline_case <- c(1, 
                             0, 0,
                             0, 1)
          
        } else if (MR_model_for_real_data %in% c("MR_model_10", "MR_model_10_b")) { ## (intercept + ref_test + RoB)
          
          print(paste("10"))
          intercept_only <- FALSE
          X <- subset_X_covariates( X = X,
                                    covariates_to_keep = c("intercept",
                                                           "low_RoB_QUADAS_clean",
                                                           "Ref_test_clean_SCID", "Ref_test_clean_Structured"))
          baseline_case <- c(1,
                             0, 
                             0, 0)
          
        } else if (MR_model_for_real_data %in% c("MR_model_11", "MR_model_11_b")) {  ## (intercept + study_setting + RoB)
          
          print(paste("11"))
          intercept_only <- FALSE
          X <- subset_X_covariates( X = X,
                                    covariates_to_keep = c("intercept",
                                                           "low_RoB_QUADAS_clean",
                                                           "study_setting_2", "study_setting_3"))
          baseline_case <- c(1,
                             0,
                             0, 1)
          
        } else if (MR_model_for_real_data %in% c("MR_model_12", "MR_model_12_b")) {  ## (intercept + logit_prev_GAD + ref_test + study_setting)
          
          print(paste("12"))
          intercept_only <- FALSE
          X <- subset_X_covariates( X = X,
                                    covariates_to_keep = c("intercept",
                                                           "logit_prev_GAD",
                                                           "study_setting_2", "study_setting_3",
                                                           "Ref_test_clean_SCID", "Ref_test_clean_Structured"))
          baseline_case <- c(1, 
                             0.0, 
                             0, 0, 
                             0, 1)
          
        } else if (MR_model_for_real_data %in% c("MR_model_13", "MR_model_13_b")) {  ## (intercept + logit_prev_GAD + ref_test + RoB)
          
          print(paste("13"))
          intercept_only <- FALSE
          X <- subset_X_covariates( X = X,
                                    covariates_to_keep = c("intercept",
                                                           "logit_prev_GAD",
                                                           "low_RoB_QUADAS_clean",
                                                           "Ref_test_clean_SCID", "Ref_test_clean_Structured"))
          baseline_case <- c(1, 
                             0.0, 
                             0, 0, 
                             0)
          
        } else if (MR_model_for_real_data %in% c("MR_model_14", "MR_model_14_b")) { ## (intercept + logit_prev_GAD + study_setting + RoB)
          
          print(paste("14"))
          intercept_only <- FALSE
          X <- subset_X_covariates( X = X,
                                    covariates_to_keep = c("intercept",
                                                           "logit_prev_GAD",
                                                           "low_RoB_QUADAS_clean",
                                                           "study_setting_2", "study_setting_3"))
          baseline_case <- c(1, 
                             0.0,
                             0, 1,
                             0)
          
        } else if (MR_model_for_real_data %in% c("MR_model_15", "MR_model_15_b")) {  ## (intercept + ref_test + study_setting + RoB)
          
          print(paste("15"))
          intercept_only <- FALSE
          X <- subset_X_covariates( X = X,
                                    covariates_to_keep = c("intercept",
                                                           "low_RoB_QUADAS_clean",
                                                           "study_setting_2", "study_setting_3",
                                                           "Ref_test_clean_SCID", "Ref_test_clean_Structured"))
          baseline_case <- c(1,
                             0, 0, 
                             0, 1,
                             0)
        }
        ##
        n_covariates <- ncol(X[[1]][[1]])
        print(paste("n_covariates = ", n_covariates))
        ##
        X_not_expanded <- X
        ##
        # X <- MetaOrdDTA:::R_fn_expand_covariates_to_all_studies(
        #   X,
        #   indicator_index_test_in_study = indicator_index_test_in_study)
        ##
}
##
## Build real_data and cov_data TOGETHER at the end:
##
{
      real_data <- list()
      real_data$x <- NMA_data_list$x
      real_data$indicator_index_test_in_study <- NMA_data_list$indicator_index_test_in_study
      ##
      cov_data <- list()
      n_covariates <- ncol(X[[1]][[1]])
      ##
      X_expanded <- MetaOrdDTA:::expand_covariates_to_all_studies( 
        X_not_expanded, 
        indicator_index_test_in_study = real_data$indicator_index_test_in_study)
      ##
      cov_data$X <- X_expanded
      cov_data$baseline_case_nd <- rep(list(baseline_case), n_index_tests)
      cov_data$baseline_case_d  <- rep(list(baseline_case), n_index_tests)
      real_data$cov_data <- cov_data
      ##
      if (n_index_tests == 4) {
        saveRDS(object = real_data,
                file = file.path(local_pkg_dir, "temp_old", "Klaus_et_al_data",
                                 paste0("real_data_4_tests", "_", target_condition)))
      } else { 
        saveRDS(object = real_data,
                file = file.path(local_pkg_dir, "temp_old", "Klaus_et_al_data",
                                 paste0("real_data_8_tests", "_", target_condition)))
      }
}
##
check_X_column_order(X)
check_X_column_order(cov_data$X)
##
X
validate_X_matrices(cov_data$X)
##
## Check non-expanded X:
##
validate_X_matrices( X = X_not_expanded, 
                     is_X_expanded = FALSE,
                     test_specific_covs = c("logit_prev_GAD"),  # This can vary between tests (in our example)
                     indicator_index_test_in_study = real_data$indicator_index_test_in_study)
##
## Check expanded X:
##
validate_X_matrices( X = X_expanded, 
                     is_X_expanded = TRUE,
                     test_specific_covs = c("logit_prev_GAD"),  # This can vary between tests (in our example)
                     indicator_index_test_in_study = real_data$indicator_index_test_in_study)
# ##
# ## Study 7 appears in tests 1, 2, 3:
# ##
# cov_data$X[[1]][[1]][7, ]  # Study 7 in test 1
# cov_data$X[[1]][[2]][7, ]  # Study 7 in test 2  
# cov_data$X[[1]][[3]][7, ]  # Study 7 in test 3
# # ##
# # ## Original data for Study 7:
# # ##
# # cov_data_Klaus_3 %>% filter(Study == 7)
# # ##
# # sum(X[[2]][[1]][, 3]) + 
# # sum(X[[2]][[2]][, 3]) + 
# # sum(X[[2]][[3]][, 3]) + 
# # sum(X[[2]][[4]][, 3]) 
# # ##
# # sum(X[[2]][[1]][, 4]) + 
# # sum(X[[2]][[2]][, 4]) + 
# # sum(X[[2]][[3]][, 4]) + 
# # sum(X[[2]][[4]][, 4]) 
# # ##
# # 30/76
# # ##
# #  
# # 
# #  
##
# # Now check again
# model_indicator <- create_indicator_index_test_in_study(real_data$x)
# cat("Model indicator studies per test:", colSums(model_indicator), "\n")
# 
# 
# 
# # Count ACTUAL studies (non -1) for Test 4
# sum(NMA_data_list$x[[4]][[1]][,1] != -1)
# 
# # Or use the indicator directly
# sum(NMA_data_list$indicator_index_test_in_study[,4])
# 
# 
# 
# real_data <- list()
# real_data$x <- NMA_data_list$x
# real_data$indicator_index_test_in_study <- NMA_data_list$indicator_index_test_in_study
# 
# model_indicator <- create_indicator_index_test_in_study(real_data$x)
# cat("Test 4:", sum(model_indicator[,4]), "\n")
# 
# model_indicator <- NMA_data_list$indicator_index_test_in_study
# cat("Test 4:", sum(model_indicator[,4]), "\n")
# 
# 
# 
# 
# # Which study is being miscounted?
# nma_indicator <- NMA_data_list$indicator_index_test_in_study[,4]
# fn_indicator <- create_indicator_index_test_in_study(real_data$x)[,4]
# 
# # The problem row(s):
# which(nma_indicator != fn_indicator)
# 
# 
# problem_row <- which(nma_indicator != fn_indicator)[1]
# cat("Row", problem_row, ":\n")
# cat("ND:", real_data$x[[4]][[1]][problem_row, 1:5], "\n")
# cat("D:", real_data$x[[4]][[2]][problem_row, 1:5], "\n")
# 
# 
# 
# 
# # Check ALL columns, not just first 5
# all(real_data$x[[4]][[1]][1, ] == -1)
# all(real_data$x[[4]][[2]][1, ] == -1)
# 
# # Any NA or other weird values?
# sum(is.na(real_data$x[[4]][[1]][1, ]))
# unique(real_data$x[[4]][[1]][1, ])
# 
# 
# 
# # What's in the D row?
# unique(real_data$x[[4]][[2]][1, ])
# 
# # Where is it?
# which(real_data$x[[4]][[2]][1, ] != -1)
# 



# # Which study is in model but NOT in cov_data_per_test for Test 4?
# model_studies_test4 <- test_data_list[[4]] %>%
#   select(Author, Year, Study) %>%
#   distinct()
# 
# cov_studies_test4 <- cov_data_per_test[[4]] %>%
#   select(Author, Year, Study) %>%
#   distinct()
# 
# # The problem study:
# model_studies_test4 %>%
#   anti_join(cov_studies_test4, by = "Study")
# 
# 
# 
# 
# # Is Eack still there?
# test_data_list[[4]] %>%
#   filter(grepl("Eack", Author, ignore.case = TRUE)) %>%
#   select(Author, Year) %>%
#   distinct()
# 
# 
# 






 
# ##
# ## Add covariates:
# ##
# {
#   cov_data <- list()
#   n_covariates <- ncol(X[[1]][[1]])
#   ##
#   X_expanded <- MetaOrdDTA:::expand_covariates_to_all_studies( 
#     X_not_expanded, 
#     indicator_index_test_in_study = real_data$indicator_index_test_in_study)
#   ##
#   cov_data$X <- X_expanded
#   cov_data$baseline_case_nd <- rep(list(baseline_case), n_index_tests)
#   cov_data$baseline_case_d  <- rep(list(baseline_case), n_index_tests)
#   real_data$cov_data <- cov_data
# }
# 

 







##
## --------------------------------------------------------------------------------------------------------------------------------------------
##
# ============================================================
# STEP 1: Basic counts
# ============================================================
cat("=== STEP 1: Basic counts ===\n")
cat("target_condition:", target_condition, "\n")
cat("n_studies in indicator:", nrow(indicator_index_test_in_study), "\n")
cat("n_studies in cov_data_Klaus_3:", nrow(cov_data_Klaus_3), "\n")

cat("\nStudies per test:\n")
cat("From indicator:", colSums(indicator_index_test_in_study), "\n")
cat("From cov_data_per_test:", sapply(cov_data_per_test, nrow), "\n")
cat("From test_data_list:", sapply(test_data_list, function(x) n_distinct(x$Study)), "\n")

# ============================================================
# STEP 2: Check a study that EXISTS - use first study from test 1
# ============================================================
cat("\n=== STEP 2: Spot check first study in test 1 ===\n")

first_study_id <- test_data_list[[1]] %>% pull(Study) %>% unique() %>% first()
cat("Checking Study ID:", first_study_id, "\n")

# What's this study in the covariate file?
cov_data_Klaus_3 %>% 
  filter(Study == first_study_id) %>% 
  select(Study, Author, Year, study_setting, Ref_test_clean, logit_prev_GAD, logit_prev_AAD)

# What's this study in test_data_list?
test_data_list[[1]] %>% 
  filter(Study == first_study_id) %>%
  select(Author, Year, Diagnosis, prev) %>%
  distinct()

# ============================================================
# STEP 3: Check prevalence - is it using correct column (GAD vs AAD)?
# ============================================================
cat("\n=== STEP 3: Prevalence source check ===\n")

for (t in 1:n_index_tests) {
  
    # For test 1, compare prevalence from 2x2 vs what's in cov_data_per_test
    test1_prev <- test_data_list[[t]] %>%
      group_by(Study, Author) %>%
      slice(1) %>%
      ungroup() %>%
      select(Study, Author, prev_from_2x2 = prev) %>%
      mutate(logit_prev_from_2x2 = qlogis(prev_from_2x2))
    
    cov_test1_prev <- cov_data_per_test[[t]] %>%
      select(Study, Author, logit_prev_GAD)
    
    comparison <- test1_prev %>%
      left_join(cov_test1_prev, by = c("Study", "Author")) %>%
      mutate(
        diff = abs(logit_prev_from_2x2 - logit_prev_GAD),
        match = diff < 0.01
      )
    
    cat("Prevalence matches:\n")
    print(table(comparison$match))
    
    cat("\nMismatches (if any):\n")
    comparison %>% filter(!match) %>% print(n = 10)

}


# ============================================================
# STEP 5: Check X matrix columns (for GAD):
# ============================================================
# Get prevalence from 2x2 tables
prev_from_2x2 <- test_data_list[[1]] %>%
  group_by(Study) %>%
  slice(1) %>%
  ungroup() %>%
  arrange(Study) %>%
  mutate(logit_prev = qlogis(prev)) %>%
  pull(logit_prev)

# Get prevalence from X matrix (centered)
prev_from_X <- outs_covs$MR_model_1$X_nd[[1]][, "logit_prev_GAD"]

# Correlation should be 1.0 if mapping is correct
cat("Correlation:", cor(prev_from_2x2, prev_from_X, use = "complete.obs"), "\n")

# Also check they're in same order
plot(prev_from_2x2, prev_from_X, 
     xlab = "logit_prev from 2x2", 
     ylab = "logit_prev from X (centered)",
     main = paste("Correlation:", round(cor(prev_from_2x2, prev_from_X, use = "complete.obs"), 4)))
abline(lm(prev_from_X ~ prev_from_2x2), col = "red")



# Are you accidentally using GAD prevalence for AAD analysis?
cat("target_condition:", target_condition, "\n")
cat("Column name in X:", colnames(outs_covs$MR_model_1$X_nd[[1]])[2], "\n")

# Check if prev in test_data_list matches the CORRECT target
# (i.e., for AAD target, prev should be AAD prevalence)
test_data_list[[1]] %>%
  group_by(Study) %>%
  slice(1) %>%
  select(Study, Author, Diagnosis, prev) %>%
  head(5)


# ============================================================
# STEP 5: Check X matrix columns (for AAD):
# ============================================================
for (t in 1:n_index_tests) {
  
    cat("=== AAD CHECK ===\n")
    cat("target_condition:", target_condition, "\n")
    
    # Correlation check for AAD
    prev_from_2x2 <- test_data_list[[t]] %>%
      group_by(Study) %>%
      slice(1) %>%
      ungroup() %>%
      arrange(Study) %>%
      mutate(logit_prev = qlogis(prev)) %>%
      pull(logit_prev)
    
    prev_from_X <- outs_covs$MR_model_1$X_nd[[t]][, "logit_prev_GAD"]
    
    cat("Correlation:", cor(prev_from_2x2, prev_from_X, use = "complete.obs"), "\n")
    
    # Check for 999s (missing) in X
    cat("\nNumber of 999s in X:\n")
    for (i in 1:length(test_names)) {
      n_999 <- sum(outs_covs$MR_model_1$X_nd[[i]][, "logit_prev_GAD"] == 999, na.rm = TRUE)
      cat("Test", i, ":", n_999, "missing\n")
    }
    
    # Check diagnosis
    test_data_list[[t]] %>%
      select(Diagnosis) %>%
      distinct() %>%
      print(n = 100)
    
}

 




cat("GAD studies per test:\n")
# Run with target_condition = "GAD"
sapply(test_data_list, function(x) n_distinct(x$Study))

cat("AAD studies per test:\n")
# Run with target_condition = "AAD"
sapply(test_data_list, function(x) n_distinct(x$Study))





# Check categorical covariates match between cov_data_per_test and X matrices
for (t in 1:n_index_tests) {
  
  cat("\n=== Test", t, "===\n")
  
  # From cov_data_per_test
  cov_cats <- cov_data_per_test[[t]] %>%
    arrange(Study) %>%
    select(Study, study_setting, Ref_test_clean, low_RoB_QUADAS_clean)
  
  # From X matrix (MR_model_5 has all covariates)
  X_mat <- outs_covs$MR_model_5$X_nd[[t]]
  
  # Check dimensions match
  cat("cov_data rows:", nrow(cov_cats), ", X rows:", nrow(X_mat), "\n")
  
  # Check study_setting encoding
  # study_setting_2 = 1 if setting == 2, study_setting_3 = 1 if setting == 3
  if ("study_setting_2" %in% colnames(X_mat)) {
    setting_from_X <- case_when(
      X_mat[, "study_setting_2"] == 1 ~ 2,
      X_mat[, "study_setting_3"] == 1 ~ 3,
      TRUE ~ 1
    )
    setting_match <- all(setting_from_X == cov_cats$study_setting)
    cat("study_setting match:", setting_match, "\n")
  }
  
  # Check Ref_test encoding
  if ("Ref_test_clean_SCID" %in% colnames(X_mat)) {
    ref_from_X <- case_when(
      X_mat[, "Ref_test_clean_SCID"] == 1 ~ "SCID",
      X_mat[, "Ref_test_clean_Structured"] == 1 ~ "Structured",
      TRUE ~ "MINI"
    )
    ref_match <- all(ref_from_X == cov_cats$Ref_test_clean)
    cat("Ref_test match:", ref_match, "\n")
  }
  
  # Check RoB encoding
  if ("low_RoB_QUADAS_clean" %in% colnames(X_mat)) {
    rob_from_X <- ifelse(X_mat[, "low_RoB_QUADAS_clean"] == 1, "Yes", "No")
    rob_match <- all(rob_from_X == cov_cats$low_RoB_QUADAS_clean)
    cat("RoB match:", rob_match, "\n")
  }
}


##
## --------------------------------------------------------------------------------------------------------------------------------------------
##
source(file.path(getwd(), "R", "R_fn_covariates.R"))
validate_X_matrices(cov_data$X)
##
X
##
check_X_column_order(X)
check_X_column_order(cov_data$X)


# ##
# ## --------------------------------------------------------------------------------------------------------------------------------------------
# ##
if (run_single_model) {
  
    intercept_only <- TRUE
    # intercept_only <- FALSE
    ##
    MR_model_for_real_data
    ##
    # cts <- box_cox <- softplus <- FALSE
    # custom_file_name <- NULL ##  "DTA_NMA_Nyaga_ord_biv_RANDthr_kappa.stan"
    # # ##
    model_parameterisation <- "ord_bivariate"
    ##
    source(file.path(getwd(), "R", "R_fns_misc.R"))
    source(file.path(getwd(), "R", "R_fn_K_fold_CV.R"))
    ##
    ##
    # random_thresholds <- FALSE ; Dirichlet_random_effects_type <- "fixed"  ; compound_symmetry <- TRUE  ## Model A
    # random_thresholds <- FALSE ; Dirichlet_random_effects_type <- "fixed"  ; compound_symmetry <- FALSE  ## Model B
    # random_thresholds <- TRUE  ; Dirichlet_random_effects_type <- "kappa"  ; compound_symmetry <- TRUE  ## Model C
    random_thresholds <- TRUE  ; Dirichlet_random_effects_type <- "kappa"  ; compound_symmetry <- FALSE  ## Model D
    ##
    model_prep_obj <- MetaOrdDTA::MetaOrd_model$new(
                      debugging = TRUE,
                      ##
                      x = real_data$x,
                      indicator_index_test_in_study = create_indicator_index_test_in_study(real_data$x),
                      ##
                      intercept_only = intercept_only,
                      cov_data = real_data$cov_data,
                      ##
                      n_chains = 8,
                      ##
                      compound_symmetry = compound_symmetry,
                      ##
                      network = TRUE,
                      ##
                      prior_only = FALSE,
                      softplus = FALSE,
                      box_cox = FALSE,
                      ##
                      model_parameterisation = model_parameterisation,
                      random_thresholds = random_thresholds,
                      Dirichlet_random_effects_type = Dirichlet_random_effects_type, ## only used if random cutpoints
                      ##
                      init_lists_per_chain = NULL,
                      ##
                      ## (cmdstanr compiles only when the main model runs on cmdstanr)
                      compile = identical(algorithm_of_the_main_model, "Stan"))
    ##
    {
      init_lists_per_chain <- model_prep_obj$init_lists_per_chain
      priors <- model_prep_obj$priors
      stan_data_list <- model_prep_obj$internal_obj$outs_data$stan_data_list
    }
    ##
    stan_data_list$n_covariates_max
    ##
    length(init_lists_per_chain)
    ##
    # priors$compound_symmetry <- compound_symmetry
    # model_prep_obj$priors$compound_symmetry <- compound_symmetry
    ##
    random_thresholds
    priors$compound_symmetry
    stan_data_list$compound_symmetry
    intercept_only
    ##
    ## gc() ; gc()
    RcppParallel::setThreadOptions(numThreads = 96);
    ##
    priors
    priors$compound_symmetry
    random_thresholds
    ##
    ##
    ## ---- Set priors for the locations ("beta"):
    ##
    colSums(indicator_index_test_in_study)
    ##
    {
            
            outs_priors <- Klaus_fn_set_priors(priors = priors,
                                               stan_data_list = stan_data_list,
                                               n_index_tests = n_index_tests,
                                               n_thr = n_thr,
                                               inform_dir = inform_dir)
            priors <- outs_priors$priors
            prior_beta_mu_SD_first_3_tests <- outs_priors$prior_beta_mu_SD_first_3_tests
            prior_beta_mu_SD_last_5_tests <- outs_priors$prior_beta_mu_SD_last_5_tests
            priors
            ##
            quick_kappa_prior(log_mean = log(1000),
                              log_sd = 1,
                              df = 10)
            ##
            quick_kappa_prior(log_mean = log(100),
                              log_sd = 1,
                              df = 10)
            ##
            # ## New prior; D- group:
            # print(paste("New prior; d- group:"))
            # quick_kappa_prior(log_mean = log(100),
            #                   log_sd = 1,
            #                   df = 5)
            # quick_kappa_prior(log_mean = log(500),
            #                   log_sd = 1,
            #                   df = 5)
            # ##
            # ## Old prior; D+ group:
            # print(paste("Old prior; d+ group:"))
            # quick_kappa_prior(log_mean = log(100),
            #                   log_sd = 1,
            #                   df = 3)
            # quick_kappa_prior(log_mean = log(500),
            #                   log_sd = 1,
            #                   df = 3)
            # ## New prior; D+ group:
            # print(paste("New prior; d+ group:"))
            # quick_kappa_prior(log_mean = log(50),
            #                   log_sd = 1,
            #                   df = 15)
    }
    {
      for (kk in 1:8) {
        for (c in 1:2) {
          init_lists_per_chain[[kk]]$kappa[[c]][1:3] <- 100
          init_lists_per_chain[[kk]]$kappa[[c]][4:n_index_tests] <- 1000
        }
      }
    }
    ##
    priors$prior_dirichlet_phi <-   priors$prior_dirichlet_alpha
    ##
    priors
    ##
    ## ---- Sample:  ---------------------------------------------------------------------------
    ##
    {
        tictoc::tic("timer")
        ##
        seed <- 1
        ##
        ## ---- The sim study's settings (run_MULTIPLE_MetaOrd_simulation.R lines 344-357, 409, 412, 636-640;
        ##      more than 10 studies, ord_biv):
        ##
        n_chains <- 4
        if (random_thresholds == TRUE) {
             n_iter_for_4_chains <-  1000 ## 5000
        } else {
             n_iter_for_4_chains <-  2000 ## 1000 ## 5000
        }
        n_iter <- ceiling(n_iter_for_4_chains / (n_chains/4))
        n_burnin <- 500
        ##
        adapt_delta <- 0.80
        max_treedepth <- 10
        ##
        ## (the nuisance blocks and parameter families of interest come from MetaOrdDTA itself)
        name_of_the_NicoStan_model <- sub("(_NicoStan_fused_AVX)?\\.stan$", "_NicoStan_fused_AVX",
                                          basename(model_prep_obj$internal_obj$outs_stan_model_name$stan_model_file_name))
        nuisance_params <- MetaOrdDTA::fn_nuisance_params_of_a_MetaOrdDTA_model(
                                name_of_the_model      = name_of_the_NicoStan_model,
                                model_parameterisation = "NMA",
                                random_thresholds      = random_thresholds)
        interest_only_params <- MetaOrdDTA::fn_interest_only_params_of_a_MetaOrdDTA_model(
                                name_of_the_model      = name_of_the_NicoStan_model,
                                model_parameterisation = "NMA",
                                random_thresholds      = random_thresholds)
        NicoStan_args <- fn_NicoStan_args_of_the_sim_study(n_burnin = n_burnin,
                                                           nuisance_params = nuisance_params,
                                                           interest_only_params = interest_only_params)
        ## (the draws also hold the transformed parameters, as a cmdstanr fit does: e.g. beta_tau)
        NicoStan_args$compute_transformed_parameters <- TRUE
        ##
        n_chains*n_iter
        ##
        random_thresholds
        priors$compound_symmetry
        intercept_only
        ##
        init_lists_per_chain <- resize_init_list( init_lists_per_chain = init_lists_per_chain,
                                                  n_chains_new = n_chains)
        ##
        model_samples_obj <-  model_prep_obj$sample(
          seed = seed,
          ##
          n_burnin = n_burnin,
          n_iter   = n_iter,
          adapt_delta = adapt_delta,
          max_treedepth = max_treedepth,
          metric_shape = "diag_e",
          ##
          priors = priors,
          ##
          n_chains = n_chains,
          ##
          init_lists_per_chain = init_lists_per_chain,
          ##
          algorithm = algorithm_of_the_main_model,
          NicoStan_args = NicoStan_args
        )
    }
    ##
    {
        ##
        ## ----  Summarise + output results: ---------------------------------------------------
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
    
        try({
          print(tictoc::toc(log = TRUE))
          log.txt <- tictoc::tic.log(format = TRUE)
          tictoc::tic.clearlog()
          time_total <- unlist(log.txt)
          ##
          extract_numeric_string <-  stringr::str_extract(time_total, "\\d+\\.\\d+")
          time_total <- as.numeric(extract_numeric_string)
          ##
          print(paste("time_total = ",  time_total))
        })
    }
    ##
    if (n_index_tests > 4) {
       test_names = c("GAD-2", "GAD-7", "HADS", "BAI",
                      "OASIS", "STAI_S", "STAI_T", "PROMIS")
    } else {
       test_names <- c("GAD-2", "GAD-7", "HADS", "BAI")
    }
    ##
    tibble_main <- model_summary_and_trace_obj$get_summary_main()                 %>% print(n = 200)
    tibble_tp   <- model_summary_and_trace_obj$get_summary_transformed()          %>% print(n = 100)
    tibble_gq   <- model_summary_and_trace_obj$get_summary_generated_quantities() %>% print(n = 1000)
    ##
    ## tibble_all <- rbind(tibble_main, tibble_gq)
    tibble_all <- rbind(tibble_main, tibble_tp, tibble_gq)
    ##
    Se <- model_summary_and_trace_obj$extract_params(params = c("Se_baseline")) %>% print(n = 50)
    Sp <- model_summary_and_trace_obj$extract_params(params = c("Sp_baseline")) %>% print(n = 20)
    Fp <- model_summary_and_trace_obj$extract_params(params = c("Fp_baseline")) %>% print(n = 20)
    min(c(Se$n_eff, Sp$n_eff), na.rm = TRUE) ; round(max(c(Se$Rhat, Sp$Rhat), na.rm = TRUE), 2)
    ##
    kappa <- model_summary_and_trace_obj$extract_params(params = c("kappa")) %>% print(n = 20)
    min(kappa$n_eff) ;  round(max(kappa$Rhat), 2)
    ##
    model_summary_and_trace_obj$get_divergences()
    model_summary_and_trace_obj$get_HMC_diagnostic_info()
    model_summary_and_trace_obj$get_efficiency_metrics()$time_total
    ##
    model_summary_and_trace_obj$extract_params(params = c("beta_tau")) %>% print(n = 20)
    model_summary_and_trace_obj$extract_params(params = c("beta_sigma")) %>% print(n = 20)
    ##
    log_lik_study <- model_summary_and_trace_obj$extract_params(params = c("log_lik_study")) %>% print(n = 20)
    min(c(log_lik_study$n_eff), na.rm = TRUE) ; round(max(c(log_lik_study$Rhat), na.rm = TRUE), 2)
    ## model_summary_and_trace_obj$get_HMC_info()
    ##
    ## log_kappa[c][t] ~ normal(0.0, 2.5) is WORSE than:  log_kappa[c][t] ~ normal(0.0, 5.0);
    ##
    ## ---- kappa[c][t] ~ normal(0.0, 5000):
##
    ## Min_ESS(Se, Sp) = 181, Max_rhat(Se, Sp) = 1.04, Min_ESS(kappa) = 13, Max_rhat(kappa) = 1.73
    ##
    ## ---- kappa[c][t] ~ normal(0.0, 1000):
##
    ## Min_ESS(Se, Sp) = 210, Max_rhat(Se, Sp) = 1.03, Min_ESS(kappa) = 40, Max_rhat(kappa) = 1.14
    ##
    ## ---- w/ student_t priors: log_100_100_df_5_5_sd_1_1:
##
    ## Min_ESS(Se, Sp) = 328, Max_rhat(Se, Sp) = 1.03, Min_ESS(kappa) = 9, Max_rhat(kappa) = 3.14
}
# # ##
# # ## ---- Extract NMA AUC metrics:
# # ##
# # # tibble_NMA_AUC_metrics <- model_summary_and_trace_obj$extract_AUC(test_names = test_names)
# # # tibble_NMA_AUC_metrics
# # # ##
# # # dplyr::filter(tibble_all, (stringr::str_detect(parameter, "AUC"))) %>% print(n = 100)
# # ##
# # # model_summary_and_trace_obj$extract_params(params = c("kappa")) %>% print(n = 20)
# # # ##
# # # intercept_only
# # # 
# # # min(model_summary_and_trace_obj$extract_params(params = c("log_lik_study"))$n_eff)
# # # 
# model_summary_and_trace_obj$get_HMC_info()
# # 
# # model_summary_and_trace_obj$extract_params(params = c("beta_mu")) %>% print(n = 20)
# # 
# # ## w/ log(200/50), df = 15, log_sd = 0.75 priors:
# # # A tibble: 8 × 10
# # # param_group parameter   mean     sd `2.5%` `50%` `97.5%` n_eff  Rhat n_Rhat
# # # <chr>       <chr>      <dbl>  <dbl>  <dbl> <dbl>   <dbl> <dbl> <dbl> <lgl> 
# # #   1 kappa       kappa[1,1] 105.   24.8    65.9 102.    160.   1179  1.01 NA    
# # # 2 kappa       kappa[2,1]  25.0   5.51   16.2  24.2    37.3  1379  1.01 NA    
# # # 3 kappa       kappa[1,2] 202.   31.7   149.  199.    273.    169  1.07 NA    
# # # 4 kappa       kappa[2,2]  97.5  20.3    66.1  94.5   144.    291  1.04 NA    
# # # 5 kappa       kappa[1,3] 374.   71.1   262.  364.    540.    147  1.08 NA    
# # # 6 kappa       kappa[2,3] 260.   97.6   131.  237.    502.    101  1.11 NA    
# # # 7 kappa       kappa[1,4] 688.  240.    376.  636.   1293.    180  1.07 NA    
# # # 8 kappa       kappa[2,4] 974.  560.    347.  814.   2412.     35  1.38 NA 
# #
# #
# ##
# ## ---- Save full model output:
# ##
# dummy_data <- 0
# ##
# {
#     if (dummy_data == 1) {
#       output_dir <- MetaOrdDTA::metaorddta_output_path("application_results_dummy_data")
#     } else {
#       output_dir <- MetaOrdDTA::metaorddta_output_path("application_results_real_data")
#     }
#     ##
#     if (stan_data_list$n_covariates_max == 1) {
#       intercept_only <- 1
#     } else {
#       intercept_only <- 0
#     }
#     ##
#     file_name_string <- paste0( "dummy_data_", dummy_data,
#                                 ##
#                                 "n_tests_", n_index_tests, ## new
#                                 "tar_cond_",  target_condition, ## new
#                                 ##
#                                 "intercept_only_", intercept_only,
#                                 "N_covs_", stan_data_list$n_covariates_max,
#                                 "param_", model_parameterisation,
#                                 "rand_thr_", random_thresholds,
#                                 "CS_", compound_symmetry,
#                                 "het_sigma_", 0)
#     if (intercept_only == FALSE) {
#       file_name_string <- paste0(file_name_string, 
#                                  "_", MR_model_for_real_data)
#     }
#     ##
#     # Save with overwrite check
#     results_list_file <- file.path(output_dir, target_condition, paste0("seed_", seed, "_", file_name_string, "_applied_results.RDS"))
#     ##
#     ## Remove if exists
#     try({
#       if (file.exists(results_list_file)) {
#         cat(sprintf("Removing existing file: %s\n", results_list_file))
#         file.remove(results_list_file)
#       }
#     })
#     # Save:
#     try({
#       results_list_file <- file.path(getwd(), output_dir, target_condition, paste0("seed_", seed, "_", file_name_string, "_applied_results.RDS"))
#       saveRDS(object = list("model_prep_obj" = model_prep_obj,
#                             "model_summary_and_trace_obj" = model_summary_and_trace_obj),
#       file = results_list_file)
#     })
# }
# 
# 
# 
# 
# # 
# # MetaOrdDTA:::extract_params_from_tibble_batch(debugging = FALSE,
# #                                               tibble = tibble_gq,
# #                                               param_strings_vec = c("Se"),
# #                                               condition = "containing") %>% print(n=100)
# # 
# # 
# # 
# # 
# # MetaOrdDTA:::extract_params_from_tibble_batch(debugging = FALSE,
# #                                               tibble = tibble_all,
# #                                               param_strings_vec = c("beta"),
# #                                               condition = "containing") %>% print(n=100)
# # 
# # 
# #  
# 
# 
# # ##
# # ## ----------- Plot sROC curve(s) (NMA - if WITHOUT covariates):
# # ##
# # plots <- model_summary_and_trace_obj$plot_sROC(test_names = test_names)
# # ##
# # plots$plot_list[[1]]
# # ##
# # ##
# # ## ---- Plot sROC curve(s) (NMA - if WITH covariates):
# # ##
# # new_cov_data <- list()
# # ##
# # baseline_case <- c(1,   ## intercept
# #                    ##
# #                    0.0, ## prev_GAD
# #                    ##
# #                    0,   ## low_RoB_QUADAS
# #                    ##
# #                    0,   ## Ref_test_SCID
# #                    0,   ## Ref_test_Structured
# #                    ##
# #                    0,   ## study_setting_2
# #                    1)   ## study_setting_3
# # ##
# # new_cov_data$baseline_case_nd <- baseline_case
# # new_cov_data$baseline_case_d  <- baseline_case
# ##
# ##
# # if (model_parameterisation == "R&G") { 
# # #   new_cov_data$baseline_case    <- rep(list(new_cov_data$baseline_case), n_index_tests)
# # # } else { 
# #   new_cov_data$baseline_case_nd <- rep(list(new_cov_data$baseline_case_nd), n_index_tests)
# #   new_cov_data$baseline_case_d  <- rep(list(new_cov_data$baseline_case_d), n_index_tests)
# # # }
# # ##
# # plots <- model_summary_and_trace_obj$plot_sROC(new_cov_data = new_cov_data, 
# #                                                test_names = test_names)
# # ##
# # plots$plot_list[[1]]
# # plots$plot_list[[2]]
# # plots$plot_list[[3]]
# # plots$plot_list[[4]]
# # plots$plot_list[[5]]
# # 
# # 
# # X_summary <- R_fn_summarize_covariates_by_test(X_list = real_data$cov_data$X[[1]])
# # X_summary
# # ##
# # R_fn_generate_sroc_scenarios(summary_results = X_summary, 
# #                              max_scenarios = 10)
# 
#  
# model_summary_and_trace_obj$get_HMC_diagnostic_info()
# 
# 
# model_summary_and_trace_obj$get_efficiency_metrics()$time_total
# 
# 
# compute_study_overlap_prop(indicator_index_test_in_study)
# # Total studies: 76 
# # Studies with 1 test: 38 
# # Studies with 2+ tests: 26 
# # Study overlap proportion: 0.342 
# # [1] 0.3421053
# 
# priors$prior_beta_mu_SD
# 
# 
# 
# 
# 
# 
# model_summary_and_trace_obj$extract_params(params = c("kappa")) %>% print(n = 20)
# beta_mu <- model_summary_and_trace_obj$extract_params(params = c("beta_mu")) %>% print(n = 100)
# ## Min ESS w/ new priors after 500 iter + 16 chains for coeffs = 156 (~ 0)
# ## Min ESS w/ new priors after 500 iter + 16 chains for coeffs = 91 (but >1% divs)
# 
# 
# R_fn_map_beta_coefficients(beta_tibble = beta_mu, 
#                            covariate_names = dimnames(real_data$cov_data$X[[1]][[1]])[[2]], 
#                            test_names = test_names, 
#                            alpha = 0.05)
# 
# 
# 
# min(model_summary_and_trace_obj$extract_params(params = c("log_lik_study"))$n_eff)
# ##
# ## seed_1_dummy_data_0n_tests_8tar_cond_AADintercept_only_1N_covs_1param_ord_bivariaterand_thr_FALSECS_1het_sigma_0_applied_results.RDS
# ## seed_2_dummy_data_0n_tests_8tar_cond_GADintercept_only_1N_covs_1param_ord_bivariaterand_thr_FALSECS_TRUEhet_sigma_0_applied_results.RDS
# Model_A <- readRDS(file.path(
#   MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
#   "old",
#   # target_condition,
#   "seed_1_dummy_data_0n_tests_8tar_cond_AADintercept_only_1N_covs_1param_ord_bivariaterand_thr_FALSECS_1het_sigma_0_applied_results.RDS"))
# ##
# Model_B <- readRDS(file.path(
#   MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
#   target_condition,
#   "seed_1_dummy_data_0n_tests_8tar_cond_AADintercept_only_1N_covs_1param_ord_bivariaterand_thr_FALSECS_0het_sigma_0_applied_results.RDS"))
# ##
# Model_C <- readRDS(file.path(
#   MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
#   target_condition,
#   "seed_1_dummy_data_0n_tests_8tar_cond_AADintercept_only_1N_covs_1param_ord_bivariaterand_thr_TRUECS_1het_sigma_0_applied_results.RDS"))
# ##
# Model_D <- readRDS(file.path(
#   MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
#   target_condition,
#   "seed_1_dummy_data_0n_tests_8tar_cond_AADintercept_only_1N_covs_1param_ord_bivariaterand_thr_TRUECS_0het_sigma_0_applied_results.RDS"))
# ##
# model_list <- list(Model_A)
# 
# 
# 
# 
# ##
# ## For all AUCs:
# ##
# for (i in 1:length(test_names)) {
#   cat(test_names[i], "&", extract_param_row(model_list, paste0("AUC[", i, "]")), "\\\\\n")
# }
# 
# 
# ##
# ## ---- Se/Sp table 1: -----------------------------------------------------------------------------------------------------------
# ##
# ##
# ## Get all Se/Sp sections at once:
# ##
# all_se_sp <- create_all_se_sp_sections( model_list, 
#                                         n_tests = 8, 
#                                         max_n_thr = 63, 
#                                         digits = 3)
# cat(paste(all_se_sp, collapse = "\n"))
# ##
# ## Define thresholds
# ##
# thr_screening_reference_tests <- list(
#   "GAD-2"  = c(1, 2, 3),
#   "GAD-7"  = c(8, 9, 10, 11),
#   "HADS"   = c(8, 9, 10, 11),
#   "BAI"    = c(14, 15, 16, 17, 18, 20),
#   "OASIS"  = c(6, 7, 8, 9, 10, 11),
#   "STAI_S" = c(64, 68, 72, 74, 75, 76) - 20,
#   "STAI_T" = c(64, 68, 72, 74, 75, 76) - 20,
#   "PROMIS" = c(18, 20, 21, 22, 23, 24)
# )
# ##
# ## Optional: Define labels for key thresholds
# ##
# thr_labels <- list(
#   "GAD-2"  = c("(standard screening thr.)"),
#   "GAD-7"  = c(" ", " ", "(moderate; standard)", " "),
#   "HADS"   = c("(possible)", " ", " ", "(probable)"),
#   "BAI"    = c(" ", " ", "(moderate)", " ", " ", " "),
#   "OASIS"  = c(" ", " ", "(standard)", " ", " ", " "),
#   "STAI_S" = c(" ", " ", " ", " ", " ", " "),
#   "STAI_T" = c(" ", " ", " ", " ", " ", " "),
#   "PROMIS" = c(" ", " ", " ", " ", " ", " ")
# )
# ##
# ## Define test number mapping (if different from list order)
# ##
# test_numbers <- c(
#   "GAD-2"  = 1,
#   "GAD-7"  = 2,
#   "HADS"   = 3,
#   "BAI"    = 4,
#   "OASIS"  = 5,
#   "STAI_S" = 6,
#   "STAI_T" = 7,
#   "PROMIS" = 8
# )
# ##
# ## Generate LaTeX lines
# ##
# latex_lines <- create_all_se_sp_sections(
#   model_list = model_list,
#   thresholds_list = thr_screening_reference_tests,
#   threshold_labels = thr_labels,  # Optional - can be NULL
#   digits = 2,
#   n_tests = 8,
#   max_n_thr = 63,
#   test_numbers = test_numbers,  # Optional - can be NULL if order matches
#   add_midrule_between_tests = TRUE
# )
# ##
# cat(latex_lines, sep = "\n")
# 
# 
# ##
# ## ---- Se/Sp table 2 (single/best model): ---------------------------------------------------------------------------------------------------
# ##
# # Define thresholds and labels
# # thr_screening_reference_tests <- list(
# #   "GAD-2"  = c(1, 2, 3),
# #   "GAD-7"  = c(8, 9, 10, 11),
# #   "HADS"   = c(8, 9, 10, 11),
# #   "BAI"    = c(14, 15, 16, 17, 18, 20),
# #   "OASIS"  = c(6, 7, 8, 9, 10, 11),
# #   "STAI_S" = c(44, 48, 52, 54, 55, 56),
# #   "STAI_T" = c(44, 48, 52, 54, 55, 56),
# #   "PROMIS" = c(18, 20, 21, 22, 23, 24)
# # )
# ##
# thr_screening_reference_tests <- list(
#   "GAD-2"  = c(1, 2, 3),
#   "GAD-7"  = c(4, 5, 6, 7),
#   "HADS"   = c(5, 6, 7, 8),
#   "BAI"    = c(2, 4, 6, 8, 10, 12),
#   "OASIS"  = c(1, 2, 3, 4, 5, 6),
#   "STAI_S" = c(52, 54, 56, 58, 60, 62, 64, 66, 68) - 20,
#   "STAI_T" = c(52, 54, 56, 58, 60, 62, 64, 66, 68) - 20,
#   "PROMIS" = c(12, 14, 16, 18, 20, 22)
# )
# ##
# # thr_labels <- list(
# #   "GAD-2"  = c("", "", "(screening)"),
# #   "GAD-7"  = c("", "", "(moderate)", ""),
# #   "HADS"   = c("(possible)", "", "", "(probable)"),
# #   "BAI"    = c("", "", "(moderate)", "", "", ""),
# #   "OASIS"  = c("", "", "(standard)", "", "", ""),
# #   "STAI_S" = c("", "", "", "", "", ""),
# #   "STAI_T" = c("", "", "", "", "", ""),
# #   "PROMIS" = c("", "", "", "", "", "")
# # )
# thr_labels <- list(
#   "GAD-2"  = c("", "", ""),
#   "GAD-7"  = c("", "", "", ""),
#   "HADS"   = c("", "", "", ""),
#   "BAI"    = c("", "", "", "", "", ""),
#   "OASIS"  = c("", "", "", "", "", ""),
#   "STAI_S" = c("", "", "", "", "", ""),
#   "STAI_T" = c("", "", "", "", "", ""),
#   "PROMIS" = c("", "", "", "", "", "")
# )
# ##
# # Generate table
# latex_table <-  create_se_sp_table_compact(
#                 model_list = model_list,
#                 thresholds_list = thr_screening_reference_tests,
#                 threshold_labels = thr_labels,
#                 digits = 2,
#                 n_tests = 8,
#                 max_n_thr = 63,
#                 caption = "Sensitivity and specificity estimates from baseline NMA-DTA model",
#                 label = "tab:se_sp_baseline")
# ##
# writeLines(latex_table)
# 
#  
#  
# 
# 
# ##
# ## Get all estimates for table 2:
# ##
# table2_lines <- create_table2_latex(model_list, digits = 3)
# cat(paste(table2_lines, collapse = "\n"))
# 
# 
# 
# 
# Model_A$model_summary_and_trace_obj$get_HMC_info()
# 
#  
# Model_real <-  Model_A
# 
# # model_summary_and_trace_obj_REAL <- Model_real$model_summary_and_trace_obj
# 
# source(file.path("R", "R_fn_NMA_computations.R"))
# ##
# AUC_outs_REAL <- R_fn_extract_NMA_AUC_summary( tibble_gq = Model_real$model_summary_and_trace_obj$get_summary_generated_quantities(),
#                                                n_thr = n_thr,
#                                                test_names = test_names,
#                                                 network = TRUE
#                                                )
# 
# 
# AUC_outs_REAL$auc
# AUC_outs_REAL$auc_diff
# AUC_outs_REAL$auc_pred
# 
# # AUC_outs_REAL$auc
# # # A tibble: 8 × 7
# # test test_name AUC_mean AUC_median AUC_lower AUC_upper AUC_sd
# # <dbl> <chr>        <dbl>      <dbl>     <dbl>     <dbl>  <dbl>
# #   1     1 GAD-2        0.827      0.828     0.785     0.864 0.0202
# # 2     2 GAD-7        0.857      0.858     0.823     0.887 0.0161
# # 3     3 HADS         0.859      0.860     0.824     0.891 0.0170
# # 4     4 BAI          0.847      0.850     0.768     0.909 0.0362
# # 5     5 OASIS        0.790      0.796     0.659     0.891 0.0595
# # 6     6 STAI_S       0.781      0.787     0.636     0.888 0.0639
# # 7     7 STAI_T       0.773      0.780     0.612     0.889 0.0713
# # 8     8 PROMIS       0.822      0.827     0.701     0.910 0.0538
# # > AUC_outs_REAL$auc_diff
# # # A tibble: 28 × 10
# # test1 test2 test1_name test2_name AUC_diff_mean AUC_diff_median AUC_diff_lower AUC_diff_upper AUC_diff_sd prob_test1_better
# # <int> <int> <chr>      <chr>              <dbl>           <dbl>          <dbl>          <dbl>       <dbl>             <dbl>
# #   1     1     2 GAD-2      GAD-7           -0.0302        -0.0300          -0.0709        0.00823      0.0201                NA
# # 2     1     3 GAD-2      HADS            -0.0326        -0.0325          -0.0828        0.0165       0.0252                NA
# # 3     2     3 GAD-7      HADS            -0.0200        -0.0221          -0.0925        0.0633       0.0400                NA
# # 4     1     4 GAD-2      BAI              0.0367         0.0314          -0.0682        0.169        0.0606                NA
# # 5     2     4 GAD-7      BAI              0.0459         0.0404          -0.0693        0.193        0.0670                NA
# # 6     3     4 HADS       BAI              0.0542         0.0470          -0.0702        0.218        0.0740                NA
# # 7     1     5 GAD-2      OASIS            0.00494       -0.000336        -0.0909        0.126        0.0556                NA
# # 8     2     5 GAD-7      OASIS           -0.00245       -0.00264         -0.0456        0.0413       0.0220                NA
# # 9     3     5 HADS       OASIS            0.0102         0.00784         -0.0586        0.0916       0.0384                NA
# # 10     4     5 BAI        OASIS            0.0669         0.0618          -0.0356        0.198        0.0597                NA
# # # ℹ 18 more rows
# # # ℹ Use `print(n = ...)` to see more rows
# # > AUC_outs_REAL$auc_pred
# # # A tibble: 8 × 7
# # test test_name AUC_mean AUC_median AUC_lower AUC_upper AUC_sd
# # <dbl> <chr>        <dbl>      <dbl>     <dbl>     <dbl>  <dbl>
# #   1     1 GAD-2        0.804      0.823     0.544     0.959  0.109
# # 2     2 GAD-7        0.831      0.851     0.573     0.970  0.104
# # 3     3 HADS         0.837      0.858     0.585     0.972  0.102
# # 4     4 BAI          0.824      0.845     0.550     0.972  0.110
# # 5     5 OASIS        0.768      0.792     0.446     0.960  0.135
# # 6     6 STAI_S       0.754      0.780     0.411     0.955  0.143
# # 7     7 STAI_T       0.746      0.771     0.405     0.950  0.144
# # 8     8 PROMIS       0.801      0.826     0.496     0.970  0.125
# # > 
# 
# 
# ##
# ## ---- Make LaTeX  AUC table:
# ##
# AUC_latex_table_REAL <- create_AUC_latex_table( auc_outs = AUC_outs_REAL, 
#                                                 digits = 3)
# cat(AUC_latex_table_REAL)


# ## sROC plot:
# relevant_thresholds = list( "GAD-2" = c(1:6), 
#                             "GAD-7" =  c(3:18), 
#                             "HADS" =  c(3:18),
#                             "BAI" = c(3:38),
#                             "OASIS" = c(1:20),
#                             "STAI_S" = c(3:40),
#                             "STAI_I" = c(3:40),
#                             "PROMIS" = c(1:32))
# 
# 
# # relevant_thresholds <- NULL
# tibbles_provided = FALSE
# tibbles = NULL
# ##
# base_size = 16
# point_size = 3
# linewidth = 0.5
# ##
# ##
# conf_region_colour = "blue"
# pred_region_colour = "green"
# ##
# conf_region_alpha = 0.15
# pred_region_alpha = 0.3
# 
# stan_model_file_name = model_summary_and_trace_obj_REAL$internal_obj$outs_stan_model_name$stan_model_file_name
# stan_mod_samples     =  model_summary_and_trace_obj_REAL$internal_obj$outs_stan_sampling$stan_mod_samples
# 
# ##
# source(file.path("R", "R_fn_test_accuracy_plots_NMA.R"))
# ##
# plots_REAL <- R_fn_sROC_plot_NMA(     stan_model_file_name = model_summary_and_trace_obj_REAL$internal_obj$outs_stan_model_name$stan_model_file_name,
#                                       stan_mod_samples     =  model_summary_and_trace_obj_REAL$internal_obj$outs_stan_sampling$stan_mod_samples,
#                                       ##
#                                       test_names = test_names,
#                                       ##
#                                       relevant_thresholds = relevant_thresholds,,
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
# 
# 
# 
# 
# 
# 
# # Load JUST the compiled library without attaching the package
# library.dynam("MetaOrdDTA", package = "MetaOrdDTA")
# 
# # Verify it's loaded
# "MetaOrdDTA" %in% names(getLoadedDLLs())
# 
# # Now check what's available
# getDLLRegisteredRoutines("MetaOrdDTA")$.Call

##
## log(100) = 4.60517
##
#
##
## ---- K-fold ---------------------------------------------------------------------------------------------------------------------------
##
source(file.path(getwd(), "R", "R_fns_misc.R"))
source(file.path(getwd(), "R", "R_fn_K_fold_CV.R"))
##
inform_dir <- FALSE
# inform_dir <- TRUE
##
MR_model_for_real_data
##
colSums(indicator_index_test_in_study)
##
if (fixed_C_indicator == TRUE) {
    if (CS_indicator == TRUE) {
        random_thresholds <- FALSE ; Dirichlet_random_effects_type <- "fixed"  ; compound_symmetry <- TRUE  ## Model A
    } else {
        random_thresholds <- FALSE ; Dirichlet_random_effects_type <- "fixed"  ; compound_symmetry <- FALSE  ## Model B  ## Best according to K-fold
    }
} else { 
    if (CS_indicator == TRUE) {
        random_thresholds <- TRUE ; Dirichlet_random_effects_type <- "kappa" ; compound_symmetry <- TRUE  ## Model C
    } else {
        random_thresholds <- TRUE ; Dirichlet_random_effects_type <- "kappa"  ; compound_symmetry <- FALSE  ## Model D
    }
}
##
random_thresholds
compound_symmetry
##
model_prep_obj <- MetaOrdDTA::MetaOrd_model$new(
                  debugging = TRUE,
                  ##
                  x = real_data$x,
                  indicator_index_test_in_study = create_indicator_index_test_in_study(real_data$x),
                  ##
                  intercept_only = intercept_only,
                  cov_data = real_data$cov_data,
                  ##
                  n_chains = 8,
                  ##
                  network = TRUE,
                  ##
                  prior_only = FALSE,
                  softplus = FALSE,
                  box_cox = FALSE,
                  ##
                  compound_symmetry = compound_symmetry,
                  ##
                  model_parameterisation = model_parameterisation,
                  random_thresholds = random_thresholds,
                  Dirichlet_random_effects_type = Dirichlet_random_effects_type, ## only used if random cutpoints
                  ##
                  init_lists_per_chain = NULL,
                  ##
                  ## (the K-fold is NicoStan only: no cmdstanr compile)
                  compile = FALSE)
##
{
    init_lists_per_chain <- model_prep_obj$init_lists_per_chain
    priors <- model_prep_obj$priors
    stan_data_list <- model_prep_obj$internal_obj$outs_data$stan_data_list
}
##
stan_data_list$n_covariates_max
stan_data_list$n_index_tests
##
stan_data_list$baseline_case_d
##
length(init_lists_per_chain)
##
priors$compound_symmetry
model_prep_obj$priors$compound_symmetry 
stan_data_list$compound_symmetry
##
random_thresholds
priors$compound_symmetry 
intercept_only
##
## gc() ; gc()
RcppParallel::setThreadOptions(numThreads = 192);
##
priors
priors$compound_symmetry 
random_thresholds
##
## ---- Set priors for the locations ("beta"):
##
colSums(indicator_index_test_in_study)
##
# MR_model_for_real_data_vec_batch <- c( "MR_model_4_b")
# MR_model_for_real_data_vec_batch <- c( "MR_model_11_b")
# MR_model_for_real_data_vec_batch <- c( "MR_model_9_b")
# MR_model_for_real_data_vec_batch <- c( "MR_model_2_b")
# MR_model_for_real_data_vec_batch <- c( "MR_model_7_b")
# MR_model_for_real_data_vec_batch <- c( "MR_model_7_c")
##
{
      outs_priors <- Klaus_fn_set_priors(priors = priors,
                                         stan_data_list = stan_data_list,
                                         n_index_tests = n_index_tests,
                                         n_thr = n_thr,
                                         inform_dir = inform_dir)
      priors <- outs_priors$priors
      prior_beta_mu_SD_first_3_tests <- outs_priors$prior_beta_mu_SD_first_3_tests
      prior_beta_mu_SD_last_5_tests <- outs_priors$prior_beta_mu_SD_last_5_tests
}
##
priors
##
## ---- Set default priors for the between-study corr matricex for beta ("beta_L_Omega") and "raw_scale_L_Omega"):
##
priors$beta_corr_lb <- if_null_then_set_to( priors$beta_corr_lb, -1.0)
priors$beta_corr_ub <- if_null_then_set_to( priors$beta_corr_ub, +1.0)
priors$prior_beta_corr_LKJ <- if_null_then_set_to( priors$prior_beta_corr_LKJ, 2.0)
##
priors
##
## ---- Always use same seed for every model:
##
set.seed(123, kind = "Mersenne-Twister") # set seed kind to ensure consistency
##
K <- 24
##
real_data$indicator_index_test_in_study
##
fold_assignments <- create_folds_new(
                    K = K,
                    indicator_index_test_in_study = real_data$indicator_index_test_in_study,
                    seed = 123,
                    network = TRUE)
##
## ---- Inspect pattern to make sure its the same:
##
print(table(fold_assignments))
for (k in 1:K) {
  cat("Fold", k, ":", which(fold_assignments == k), "\n")
}
##
# Fold 1 : 2 5 11 14 18 19 20 25 28 34 37 47 49 50 53 55 58 60 65 72 74 75 78 79 81 84 87 91 
# Fold 2 : 3 12 15 22 23 24 27 32 45 48 52 59 61 63 67 73 82 86 89 90 93 96 
# Fold 3 : 1 6 7 9 26 30 35 39 40 41 42 64 68 69 88 92 94 
# Fold 4 : 8 13 17 21 29 31 33 36 38 44 51 54 57 62 70 77 
# Fold 5 : 4 10 16 43 46 56 66 71 76 80 83 85 95 97 98 5 
##
n_studies <- nrow(real_data$indicator_index_test_in_study)
##
## ---- THIS MUST BE EQUAL TO 1:
##
#### dput(fold_assignments)
if (n_index_tests == 8) {
    # sum(fold_assignments ==  c(3L, 5L, 3L, 4L, 7L, 6L, 9L, 9L, 3L, 1L, 8L, 1L, 2L, 5L, 1L, 
    #                            2L, 2L, 2L, 10L, 1L, 2L, 5L, 7L, 8L, 2L, 6L, 1L, 1L, 3L, 8L, 
    #                            10L, 1L, 1L, 2L, 10L, 1L, 8L, 8L, 4L, 5L, 1L, 2L, 9L, 7L, 1L, 
    #                            3L, 10L, 5L, 1L, 6L, 1L, 2L, 3L, 5L, 1L, 4L, 1L, 10L, 3L, 1L, 
    #                            2L, 2L, 6L, 5L, 7L, 2L, 1L, 4L, 7L, 9L, 6L, 4L, 4L, 4L, 9L, 3L, 
    #                            3L, 6L, 1L, 8L, 1L, 1L, 7L))/n_studies
  sum(fold_assignments ==  c(3L, 6L, 2L, 2L, 3L, 1L, 1L, 7L, 6L, 1L, 4L, 3L, 1L, 1L, 1L, 
                             7L, 4L, 3L, 4L, 6L, 4L, 2L, 1L, 1L, 3L, 2L, 1L, 1L, 4L, 2L, 1L, 
                             2L, 7L, 8L, 1L, 6L, 3L, 2L, 2L, 1L, 5L, 4L, 5L, 1L, 1L, 8L, 5L, 
                             5L, 1L, 4L, 2L, 2L, 7L, 5L, 4L, 2L, 6L, 3L, 1L, 1L, 4L, 6L, 2L, 
                             4L, 5L, 8L, 7L, 1L, 3L, 7L, 2L, 1L, 4L, 3L, 1L, 2L, 8L, 1L, 5L, 
                             3L, 5L, 6L, 4L, 3L, 2L, 5L, 1L, 8L, 8L, 2L, 4L, 6L, 2L, 7L, 3L, 
                             2L, 1L, 3L, 5L, 8L, 1L, 3L))/n_studies
} else { 
    sum(fold_assignments ==  c(3L, 5L, 3L, 4L, 7L, 6L, 9L, 9L, 3L, 1L, 8L, 2L, 5L, 2L, 1L, 
                               5L, 7L, 2L, 1L, 8L, 10L, 2L, 3L, 10L, 2L, 8L, 1L, 10L, 4L, 1L, 
                               4L, 1L, 8L, 5L, 4L, 1L, 9L, 1L, 2L, 1L, 10L, 7L, 3L, 3L, 2L, 
                               2L, 6L, 1L, 4L, 3L, 5L, 6L, 9L, 1L, 4L, 4L, 2L, 2L, 7L, 5L, 5L, 
                               3L, 1L, 4L, 7L, 10L, 6L, 6L, 1L, 9L, 9L, 3L, 8L, 6L, 2L, 8L, 
                               1L, 2L, 7L))/n_studies
}
length(fold_assignments)
##
#
dput(fold_assignments)
##
## ---- And/or: run k-fold:
##
source(file.path(getwd(), "R", "R_fns_misc.R"))
source(file.path(getwd(), "R", "R_fn_K_fold_CV.R"))
##
## ---- The sim study's settings (run_MULTIPLE_MetaOrd_simulation.R lines 344-357, 409, 412, 636-640;
##      more than 10 studies, ord_biv), every fold on NicoStan:
##
n_chains <- 4
if (random_thresholds == TRUE) {
     n_iter_for_4_chains <-  1000 ## 5000
} else {
     n_iter_for_4_chains <-  2000 ## 1000 ## 5000
}
n_iter <- ceiling(n_iter_for_4_chains / (n_chains/4))
n_burnin <- 500
##
name_of_the_NicoStan_model <- sub("(_NicoStan_fused_AVX)?\\.stan$", "_NicoStan_fused_AVX",
                                  basename(model_prep_obj$internal_obj$outs_stan_model_name$stan_model_file_name))
nuisance_params <- MetaOrdDTA::fn_nuisance_params_of_a_MetaOrdDTA_model(
                        name_of_the_model      = name_of_the_NicoStan_model,
                        model_parameterisation = "NMA",
                        random_thresholds      = random_thresholds)
interest_only_params <- MetaOrdDTA::fn_interest_only_params_of_a_MetaOrdDTA_model(
                        name_of_the_model      = name_of_the_NicoStan_model,
                        model_parameterisation = "NMA",
                        random_thresholds      = random_thresholds)
NicoStan_args <- fn_NicoStan_args_of_the_sim_study(n_burnin = n_burnin,
                                                   nuisance_params = nuisance_params,
                                                   interest_only_params = interest_only_params)
##
random_thresholds
##
n_burnin
n_iter
n_iter*n_chains
n_chains*K
##
# n_burnin <- 10
# n_iter <- 10
##
# # Close all R connections
# closeAllConnections()
# 
# # Stop any clusters
# try(stopCluster(cl), silent = TRUE)
# gc()
##
try({
  for (kk in 1:n_chains) {
    for (c in 1:2) {
      init_lists_per_chain[[kk]]$kappa[[c]][1:3] <- 100
      init_lists_per_chain[[kk]]$kappa[[c]][4:n_index_tests] <- 1000
    }
  }
})
##
{
    outs_kfold <- R_fn_run_kfold_cv_parallel(  debugging = FALSE,
                                               ##
                                               K = K,
                                               seed = 123,
                                               ##
                                               model_prep_obj = model_prep_obj,
                                               stan_data_list = stan_data_list,
                                               ##
                                               fold_assignments = fold_assignments,
                                               ##
                                               priors = priors,
                                               init_lists_per_chain = init_lists_per_chain,
                                               ##
                                               n_burnin = n_burnin,
                                               n_iter = n_iter,
                                               ##
                                               n_chains = n_chains,
                                               ##
                                               adapt_delta = 0.80,
                                               max_treedepth = 10,
                                               ##
                                               n_workers = K,
                                               output_dir = MetaOrdDTA::metaorddta_output_path("cv_results"),
                                               ##
                                               use_BayesMVP_for_faster_summaries = TRUE,
                                               ##
                                               algorithm = "NicoStan",
                                               NicoStan_args = NicoStan_args)
    
    
    # beepr::beep(sound = 7)
    beepr::beep(sound = 9)
}
# ##
# outs_kfold$fold_results[[1]]$tibble_all %>% print(n=1000)
# ##
# tibble_fixed <- outs_kfold$fold_results[[1]]$tibble_all %>%
#   rename(parameter = variable)
# ##
# MetaOrdDTA:::extract_params_from_tibble_batch( tibble = tibble_fixed, 
#                                                condition = "exact_match",
#                                                param_strings_vec = c("beta_tau"))
# outs_kfold$fold_results[[1]]$beta_mu %>% print(n = 100)
# outs_kfold$fold_results[[1]]$beta_sigma
# outs_kfold$fold_results[[1]]$beta_tau
# outs_kfold
# #str(outs_kfold)
# outs_kfold$elpd_kfold
# outs_kfold$se_elpd
# quantile(outs_kfold$min_ESS_vec)
# quantile(outs_kfold$max_rhat_vec)
# ##
# intercept_only
# ##
##
## ---- Save K-fold results:
##
dummy_data <- 0
##
{
      if (dummy_data == 1) {
        output_dir <- MetaOrdDTA::metaorddta_output_path("application_results_dummy_data")
      } else {
        output_dir <- MetaOrdDTA::metaorddta_output_path("application_results_real_data")
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
                                  ##
                                  "n_tests_", n_index_tests, ## new
                                  "tar_cond_",  target_condition, ## new
                                  ##
                                  "intercept_only_", intercept_only,
                                  "N_covs_", stan_data_list$n_covariates_max,
                                  "param_", model_parameterisation,
                                  "rand_thr_", random_thresholds,
                                  "CS_", compound_symmetry,
                                  "prior_beta_mu_SD_", prior_beta_mu_SD_first_3_tests, "_", prior_beta_mu_SD_last_5_tests,
                                  "inform_dir_", inform_dir)
      ##
      if (intercept_only == FALSE) {
        file_name_string <- paste0(file_name_string, "_", MR_model_for_real_data)
      }
      ##
      # file_name_string <- paste0(file_name_string, "_", "normal_0_1000")
      ##
      # if (stringr::str_detect(MR_model_for_real_data, "_b")) { 
      #   file_name_string <- paste0(file_name_string, "_b")
      # }
      ##
      # results_list_file_name <- paste0(file_name_string, "_applied_results.RDS")
      # results_list_path <- file.path(output_dir, results_list_file_name)
      ##
      # Save with overwrite check
      results_list_file <- file.path( output_dir, target_condition, paste0("seed_", seed, "_", file_name_string, "_applied_results.RDS"))
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
        results_list_file <- file.path( output_dir, target_condition, paste0("seed_", seed, "_", file_name_string, "_applied_results.RDS"))
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
gc() ; gc()


    

    }
    
    
    
outs_kfold
    
    
    
##
outs_kfold$min_ESS_vec
##
stan_data_list$n_covariates_max
##
compound_symmetry
random_thresholds
##
target_condition
##
##
# model_names <- c("Intercept_only", "Prevalence", "Ref_test", "Setting")
##
## Note: accidently saved these w/ hetero_sigma = 1, eventhough it's 0 (+ redundant now anyway as removed from Stan models).
## seed_2_K_fold_dummy_data_0n_tests_8tar_cond_AADintercept_only_1N_covs_1param_ord_bivariaterand_thr_TRUECS_1het_sigma_0_applied_results.RDS
##
kfold_Model_A <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                   target_condition,
                                   paste0( "seed_1_K_fold_dummy_data_0n_tests_8tar_cond_",
                                           target_condition, 
                                           "intercept_only_1N_covs_1param_ord_bivariaterand_thr_FALSECS_TRUEprior_beta_mu_SD_0.5_0.25inform_dir_TRUE_applied_results.RDS")))
kfold_Model_B <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                   target_condition,
                                   paste0( "seed_1_K_fold_dummy_data_0n_tests_8tar_cond_", 
                                           target_condition, 
                                           "intercept_only_1N_covs_1param_ord_bivariaterand_thr_FALSECS_FALSEprior_beta_mu_SD_0.5_0.25inform_dir_TRUE_applied_results.RDS")))
kfold_Model_C <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                   target_condition,
                                   paste0( "seed_1_K_fold_dummy_data_0n_tests_8tar_cond_", 
                                           target_condition, 
                                           "intercept_only_1N_covs_1param_ord_bivariaterand_thr_TRUECS_TRUEprior_beta_mu_SD_0.5_0.25inform_dir_TRUE_applied_results.RDS")))
kfold_Model_D <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                   target_condition,
                                   paste0( "seed_1_K_fold_dummy_data_0n_tests_8tar_cond_", 
                                           target_condition,
                                           "intercept_only_1N_covs_1param_ord_bivariaterand_thr_TRUECS_FALSEprior_beta_mu_SD_0.5_0.25inform_dir_TRUE_applied_results.RDS")))
##
kfold_Model_C_normal <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                   target_condition,
                                   "seed_1_K_fold_dummy_data_0n_tests_8tar_cond_AADintercept_only_1N_covs_1param_ord_bivariaterand_thr_TRUECS_TRUEprior_beta_mu_SD_0.5_0.25inform_dir_TRUE_normal_0_1000_applied_results.RDS"))
kfold_Model_D_normal <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                   target_condition,
                                   "seed_1_K_fold_dummy_data_0n_tests_8tar_cond_AADintercept_only_1N_covs_1param_ord_bivariaterand_thr_TRUECS_FALSEprior_beta_mu_SD_0.5_0.25inform_dir_TRUE_normal_0_1000_applied_results.RDS"))
##
# kfold_Model_A2 <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
#                                    target_condition,
#                                    "seed_1_K_fold_dummy_data_0n_tests_8tar_cond_AADintercept_only_1N_covs_1param_ord_bivariaterand_thr_FALSECS_TRUEprior_beta_mu_SD_1_0.5inform_dir_FALSE_applied_results.RDS"))

##
kfold_results_list <- list(
  kfold_Model_A,
  kfold_Model_B,
  kfold_Model_C,
  kfold_Model_D
  # ##
  # kfold_Model_C_normal,
  # kfold_Model_D_normal
)
##
# kfold_Model_A$outs_kfold$K
# kfold_Model_B$outs_kfold$K
# kfold_Model_C$outs_kfold$K
# kfold_Model_D$outs_kfold$K
##
model_names <- c(
  "Model_A  (FIXED-C + CS)",
  "Model_B (FIXED-C + UN)",
  "Model_C (RAND-C + CS)",
  "Model_D (RAND-C + UN)"
  # ##
  # "Model_C (RAND-C + CS) - normal",
  # "Model_D (RAND-C + UN) - normal"
)
##
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
##
kfold_results_list[[1]]$outs_kfold$fold_assignments
##
comparison <- compare_models_k_fold_CV( kfold_results_list = kfold_results_list,
                                        model_names = model_names,
                                        min_ess_threshold = 100)
##
summarise_model_comparison_k_fold_CV(comparison)
##
## ---- Results (REAL data w/ 8 tests):
##
# comparison <- compare_kfold_models(kfold_results_list = kfold_results_list,
#                                    +                                    model_names = model_names,
#                                    +                                    min_ess_threshold = 40)
# Checking fold consistency across models...
# [1] 5
# [1] 5
# [1] 5
# ✓ All models used identical fold assignments
# 
# ESS summary across models:
#   Model_A (FIXED-C + CS) Model_B (FIXED-C + UN) Model_C (RAND-C + CS) Model_D (RAND-C + UN)
# [1,]                  405.3                  417.4                 246.6                 126.6
# [2,]                  362.2                  392.7                 188.7                 356.3
# [3,]                  369.4                  369.6                 240.8                 284.5
# [4,]                  413.0                  374.9                 219.0                 200.6
# [5,]                  407.9                  426.6                 275.2                 177.7
# 
# Minimum ESS per fold: 126.6 188.7 240.8 200.6 177.7 
# 
# 
# Final comparison using 5/5 folds:
#   Model   ELPD   SE Delta_ELPD Delta_SE
# Model_A (FIXED-C + CS) Model_A (FIXED-C + CS) -14628 3745         NA       NA
# Model_B (FIXED-C + UN) Model_B (FIXED-C + UN) -14868 3904       -240     5410
# Model_C (RAND-C + CS)   Model_C (RAND-C + CS) -15882 3947      -1254     5441
# Model_D (RAND-C + UN)   Model_D (RAND-C + UN) -16113 4128      -1485     5573
# > ##
#   > summarize_model_comparison(comparison)
# 
# Model ranking by ELPD:
#   1. Model_A (FIXED-C + CS): ELPD = -14628 (3745) [BEST]
# 2. Model_B (FIXED-C + UN): ELPD = -14868 (3904), Δ = -240 (5410) 
# 3. Model_C (RAND-C + CS): ELPD = -15882 (3947), Δ = -1254 (5441) 
# 4. Model_D (RAND-C + UN): ELPD = -16113 (4128), Δ = -1485 (5573) 
# 
# * = Significantly different from best model (|Δ| > 2*SE)










##
## seed_1_K_fold_dummy_data_0n_tests_4tar_cond_GADintercept_only_0N_covs_2param_ord_bivariaterand_thr_FALSECS_TRUEhet_sigma_0_MR_model_1_applied_results.RDS



##
## ---- Meta-regression k-fold:
##
# prior_1 <- 1.0
prior_1 <- 0.5
##
# prior_2 <- 1.0
# prior_2 <- 0.5
prior_2 <- 0.25
# prior_2 <- 0.125
## prior_2 <- 0.01 ## NOT good (tested)
##
# CS <- FALSE
CS <- TRUE
##
kfold_Model_A <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                   target_condition,
                                   paste0("seed_1_K_fold_dummy_data_0n_tests_8tar_cond_",
                                          target_condition,
                                          "intercept_only_1N_covs_1param_ord_bivariaterand_thr_FALSECS_", CS, "prior_beta_mu_SD_", 
                                   prior_1, "_",
                                   prior_2,
                                   "inform_dir_TRUE_applied_results.RDS")))
##
kfold_MR_Model_1 <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                      target_condition,
                                      paste0("seed_1_K_fold_dummy_data_0n_tests_8tar_cond_", 
                                             target_condition,
                                             "intercept_only_0N_covs_2param_ord_bivariaterand_thr_FALSECS_", CS, "prior_beta_mu_SD_",
                                             prior_1, "_",
                                             prior_2,
                                             "inform_dir_TRUE_MR_model_1_applied_results.RDS")))
kfold_MR_Model_2 <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                      target_condition,
                                      paste0("seed_1_K_fold_dummy_data_0n_tests_8tar_cond_", 
                                             target_condition,
                                             "intercept_only_0N_covs_3param_ord_bivariaterand_thr_FALSECS_", CS, "prior_beta_mu_SD_", 
                                             prior_1, "_",
                                             prior_2, 
                                             "inform_dir_TRUE_MR_model_2_applied_results.RDS")))
kfold_MR_Model_3 <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                      target_condition,
                                      paste0("seed_1_K_fold_dummy_data_0n_tests_8tar_cond_", 
                                             target_condition,
                                             "intercept_only_0N_covs_3param_ord_bivariaterand_thr_FALSECS_", CS, "prior_beta_mu_SD_", 
                                             prior_1, "_",
                                             prior_2, 
                                             "inform_dir_TRUE_MR_model_3_applied_results.RDS")))
kfold_MR_Model_4 <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                      target_condition,
                                      paste0("seed_1_K_fold_dummy_data_0n_tests_8tar_cond_", 
                                             target_condition,
                                             "intercept_only_0N_covs_2param_ord_bivariaterand_thr_FALSECS_", CS, "prior_beta_mu_SD_", 
                                             prior_1, "_",
                                             prior_2, 
                                             "inform_dir_TRUE_MR_model_4_applied_results.RDS")))
##
kfold_MR_Model_5 <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                      target_condition,
                                      paste0("seed_1_K_fold_dummy_data_0n_tests_8tar_cond_", 
                                             target_condition,
                                             "intercept_only_0N_covs_7param_ord_bivariaterand_thr_FALSECS_", CS, "prior_beta_mu_SD_", 
                                             prior_1, "_",
                                             prior_2, 
                                             "inform_dir_TRUE_MR_model_5_applied_results.RDS")))
##
kfold_MR_Model_6 <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                      target_condition,
                                      paste0("seed_1_K_fold_dummy_data_0n_tests_8tar_cond_", 
                                             target_condition,
                                             "intercept_only_0N_covs_4param_ord_bivariaterand_thr_FALSECS_", CS, "prior_beta_mu_SD_",
                                             prior_1, "_",
                                             prior_2, 
                                             "inform_dir_TRUE_MR_model_6_applied_results.RDS")))
kfold_MR_Model_7 <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                      target_condition,
                                      paste0("seed_1_K_fold_dummy_data_0n_tests_8tar_cond_", 
                                             target_condition,
                                             "intercept_only_0N_covs_4param_ord_bivariaterand_thr_FALSECS_", CS, "prior_beta_mu_SD_", 
                                             prior_1, "_",
                                             prior_2,  
                                             "inform_dir_TRUE_MR_model_7_applied_results.RDS")))
kfold_MR_Model_8 <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                      target_condition,
                                      paste0("seed_1_K_fold_dummy_data_0n_tests_8tar_cond_", 
                                             target_condition,
                                             "intercept_only_0N_covs_3param_ord_bivariaterand_thr_FALSECS_", CS, "prior_beta_mu_SD_", 
                                             prior_1, "_",
                                             prior_2, 
                                             "inform_dir_TRUE_MR_model_8_applied_results.RDS")))
kfold_MR_Model_9 <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                      target_condition,
                                      paste0("seed_1_K_fold_dummy_data_0n_tests_8tar_cond_", 
                                             target_condition,
                                             "intercept_only_0N_covs_5param_ord_bivariaterand_thr_FALSECS_", CS, "prior_beta_mu_SD_",  
                                             prior_1, "_",
                                             prior_2, 
                                             "inform_dir_TRUE_MR_model_9_applied_results.RDS")))
kfold_MR_Model_10 <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                       target_condition,
                                       paste0("seed_1_K_fold_dummy_data_0n_tests_8tar_cond_", 
                                              target_condition,
                                              "intercept_only_0N_covs_4param_ord_bivariaterand_thr_FALSECS_", CS, "prior_beta_mu_SD_",
                                              prior_1, "_",
                                              prior_2, 
                                              "inform_dir_TRUE_MR_model_10_applied_results.RDS")))
kfold_MR_Model_11 <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                       target_condition,
                                       paste0("seed_1_K_fold_dummy_data_0n_tests_8tar_cond_", 
                                              target_condition,
                                              "intercept_only_0N_covs_4param_ord_bivariaterand_thr_FALSECS_", CS, "prior_beta_mu_SD_",
                                              prior_1, "_",
                                              prior_2, 
                                              "inform_dir_TRUE_MR_model_11_applied_results.RDS")))
##
kfold_MR_Model_12 <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                       target_condition,
                                       paste0("seed_1_K_fold_dummy_data_0n_tests_8tar_cond_", 
                                              target_condition,
                                              "intercept_only_0N_covs_6param_ord_bivariaterand_thr_FALSECS_", CS, "prior_beta_mu_SD_", 
                                              prior_1, "_",
                                              prior_2,
                                              "inform_dir_TRUE_MR_model_12_applied_results.RDS")))
kfold_MR_Model_13 <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                       target_condition,
                                       paste0("seed_1_K_fold_dummy_data_0n_tests_8tar_cond_",
                                              target_condition,
                                              "intercept_only_0N_covs_5param_ord_bivariaterand_thr_FALSECS_", CS, "prior_beta_mu_SD_",
                                              prior_1, "_",
                                              prior_2, 
                                              "inform_dir_TRUE_MR_model_13_applied_results.RDS")))
kfold_MR_Model_14 <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                       target_condition,
                                       paste0("seed_1_K_fold_dummy_data_0n_tests_8tar_cond_",
                                              target_condition,
                                              "intercept_only_0N_covs_5param_ord_bivariaterand_thr_FALSECS_", CS, "prior_beta_mu_SD_",
                                              prior_1, "_",
                                              prior_2, 
                                              "inform_dir_TRUE_MR_model_14_applied_results.RDS")))
kfold_MR_Model_15 <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                       target_condition,
                                       paste0("seed_1_K_fold_dummy_data_0n_tests_8tar_cond_", 
                                              target_condition, 
                                              "intercept_only_0N_covs_6param_ord_bivariaterand_thr_FALSECS_", CS, "prior_beta_mu_SD_",
                                              prior_1, "_",
                                              prior_2, 
                                              "inform_dir_TRUE_MR_model_15_applied_results.RDS")))
##
# kfold_MR_Model_7c <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
#                                        target_condition,
#                                        "seed_1_K_fold_dummy_data_0n_tests_8tar_cond_AADintercept_only_0N_covs_4param_ord_bivariaterand_thr_FALSECS_TRUEprior_beta_mu_SD_0.5_0.25inform_dir_TRUE_MR_model_7_c_applied_results.RDS"))
# 
# kfold_MR_Model_2b <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
#                                        target_condition,
#                                        "seed_1_K_fold_dummy_data_0n_tests_8tar_cond_AADintercept_only_0N_covs_3param_ord_bivariaterand_thr_FALSECS_TRUEprior_beta_mu_SD_0.5_0.25inform_dir_TRUE_MR_model_2_b_applied_results.RDS"))
# 
# kfold_MR_Model_11b <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
#                                        target_condition,
#                                        "seed_1_K_fold_dummy_data_0n_tests_8tar_cond_AADintercept_only_0N_covs_4param_ord_bivariaterand_thr_FALSECS_TRUEprior_beta_mu_SD_0.5_0.25inform_dir_TRUE_MR_model_11_b_applied_results.RDS"))
# 
# kfold_MR_Model_7b <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
#                                        target_condition,
#                                        "seed_1_K_fold_dummy_data_0n_tests_8tar_cond_AADintercept_only_0N_covs_4param_ord_bivariaterand_thr_FALSECS_TRUEprior_beta_mu_SD_0.5_0.25inform_dir_TRUE_MR_model_7_b_applied_results.RDS"))
# 
# kfold_MR_Model_9b <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
#                                        target_condition,
#                                        "seed_1_K_fold_dummy_data_0n_tests_8tar_cond_AADintercept_only_0N_covs_5param_ord_bivariaterand_thr_FALSECS_TRUEprior_beta_mu_SD_0.5_0.25inform_dir_TRUE_MR_model_9_b_applied_results.RDS"))
# 
# kfold_MR_Model_4b <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
#                                        target_condition,
#                                        "seed_1_K_fold_dummy_data_0n_tests_8tar_cond_AADintercept_only_0N_covs_2param_ord_bivariaterand_thr_FALSECS_TRUEprior_beta_mu_SD_0.5_0.25inform_dir_TRUE_MR_model_4_b_applied_results.RDS"))



{
  kfold_Model_A$outs_kfold$fold_assignments
  ##
  kfold_MR_Model_1$outs_kfold$fold_assignments
  kfold_MR_Model_2$outs_kfold$fold_assignments
  kfold_MR_Model_3$outs_kfold$fold_assignments
  kfold_MR_Model_4$outs_kfold$fold_assignments
  ##
  kfold_MR_Model_5$outs_kfold$fold_assignments
  ##
  kfold_MR_Model_6$outs_kfold$fold_assignments
  kfold_MR_Model_7$outs_kfold$fold_assignments
  kfold_MR_Model_8$outs_kfold$fold_assignments
  kfold_MR_Model_9$outs_kfold$fold_assignments
  kfold_MR_Model_10$outs_kfold$fold_assignment
  kfold_MR_Model_11$outs_kfold$fold_assignments
  ##
  kfold_MR_Model_12$outs_kfold$fold_assignments
  kfold_MR_Model_13$outs_kfold$fold_assignments
  kfold_MR_Model_14$outs_kfold$fold_assignments
  kfold_MR_Model_15$outs_kfold$fold_assignments
}
##
kfold_Model_A$outs_kfold$K
kfold_MR_Model_1$outs_kfold$K
kfold_MR_Model_4b$outs_kfold$K
#
# kfold_MR_Model_1$outs_kfold$fold_results <- kfold_MR_Model_1$outs_kfold$fold_results[
#   order(sapply(kfold_MR_Model_1$outs_kfold$fold_results, function(x) x$fold))
# ]
# 
# kfold_MR_Model_2$outs_kfold$fold_results <- kfold_MR_Model_2$outs_kfold$fold_results[
#   order(sapply(kfold_MR_Model_2$outs_kfold$fold_results, function(x) x$fold))
# ]
# 
# # How many folds completed?
# length(kfold_MR_Model_1$outs_kfold$fold_results)  # 24
# length(kfold_MR_Model_2$outs_kfold$fold_results)  # 3
# length(kfold_MR_Model_13$outs_kfold$fold_results)
# 
# # Check for errors
# kfold_MR_Model_3$outs_kfold$errors
# 
# 
# # What's stored in the fold results?
# kfold_MR_Model_1$fold_results[[1]]$holdout_studies
# kfold_MR_Model_3$fold_results[[1]]$holdout_studies
# 
# # Or check via outs_kfold
# kfold_MR_Model_1$outs_kfold$fold_results[[1]]$holdout_studies
# kfold_MR_Model_3$outs_kfold$fold_results[[1]]$holdout_studies
# 
# 
# # What fold numbers are in each position?
# sapply(kfold_MR_Model_1$outs_kfold$fold_results, function(x) x$fold)
# sapply(kfold_MR_Model_3$outs_kfold$fold_results, function(x) x$fold)
##
colSums(indicator_index_test_in_study)
##
kfold_results_list <- list( 
  kfold_Model_A,
  kfold_MR_Model_1,
  kfold_MR_Model_2,
  kfold_MR_Model_3,
  kfold_MR_Model_4,
  ##
  kfold_MR_Model_5,
  # ## 2-way:
  kfold_MR_Model_6,
  kfold_MR_Model_7,
  kfold_MR_Model_8,
  kfold_MR_Model_9,
  kfold_MR_Model_10,
  kfold_MR_Model_11,
  ## 3-way:
  kfold_MR_Model_12,
  kfold_MR_Model_13,
  kfold_MR_Model_14,
  kfold_MR_Model_15
  # ##
  # kfold_MR_Model_7c,
  # kfold_MR_Model_2b,
  # kfold_MR_Model_11b,
  # kfold_MR_Model_7b,
  # kfold_MR_Model_9b,
  # kfold_MR_Model_4b
)
##
 model_names <- c(
                  "MR_Model_0 (intercept-only - all MR models using FIXED-C + CS)",
                  ##
                  "MR_Model_1 (intercept + prev_GAD)",
                  "MR_Model_2 (intercept + ref_test)",
                  "MR_Model_3 (intercept + study_setting)",
                  "MR_Model_4 (intercept + RoB)",
                  ##
                  "MR_Model_FULL (intercept + prev_GAD + ref_test + study_setting + RoB)",
                  ##
                  ## First 2-way MR model as prev_GAD and ref_test gave best K-fold's out of univariate MR's:
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
                  # ##
                  # "MR_Model_7c (intercept + prev_GAD [D- only] + study_setting",
                  # "MR_Model_2b (intercept + ref_test [D- only]",
                  # "MR_Model_11b (intercept + ...",
                  # "MR_Model_7b (intercept + ...",
                  # "MR_Model_9b (intercept + ...",
                  # "MR_Model_4b (intercept + ..."
)
##
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
# 
# 
# 
# kfold_results_list[[1]]$outs_kfold$fold_assignments
##
comparison <- MetaOrdDTA:::compare_models_k_fold_CV(  kfold_results_list = kfold_results_list,
                                                      model_names = model_names,
                                                      min_ess_threshold = 100)
##
MetaOrdDTA:::summarise_model_comparison_k_fold_CV(comparison)
##
## ---- Latex table for MR:
##
MetaOrdDTA:::make_latex_table_k_fold_CV( comparison_results = comparison)
##
prior_1
prior_2
##
# Minimum ESS per fold: 339 153.3 220.9 177 248.7 
# 
# 
# Final comparison using 5/5 folds:
#   Model   ELPD   SE Delta_ELPD Delta_SE
# MR_Model_0 (intercept-only - all MR models using FIXED-C + CS)               MR_Model_0 (intercept-only - all MR models using FIXED-C + CS) -15336 2149         NA       NA
# MR_Model_1 (intercept + prev_GAD)                                                                         MR_Model_1 (intercept + prev_GAD) -15241 2186         96     3066
# MR_Model_2 (intercept + ref_test)                                                                         MR_Model_2 (intercept + ref_test) -14905 2160        431     3047
# MR_Model_3 (intercept + study_setting)                                                               MR_Model_3 (intercept + study_setting) -16295 2269       -958     3125
# MR_Model_4 (intercept + RoB)                                                                                   MR_Model_4 (intercept + RoB) -15739 2334       -403     3173
# MR_Model_FULL (intercept + prev_GAD + ref_test + study_setting + RoB) MR_Model_FULL (intercept + prev_GAD + ref_test + study_setting + RoB) -15790 2216       -454     3087
# MR_Model_6  (intercept + prev_GAD + ref_test)                                                 MR_Model_6  (intercept + prev_GAD + ref_test) -14353 2153        984     3042
# MR_Model_7  (intercept + prev_GAD + study_setting)                                       MR_Model_7  (intercept + prev_GAD + study_setting) -16360 2123      -1024     3021
# MR_Model_8  (intercept + prev_GAD + RoB)                                                           MR_Model_8  (intercept + prev_GAD + RoB) -15635 2153       -299     3042
# MR_Model_9  (intercept + ref_test + study_setting)                                       MR_Model_9  (intercept + ref_test + study_setting) -15909 2238       -572     3103
# MR_Model_10 (intercept + ref_test + RoB)                                                           MR_Model_10 (intercept + ref_test + RoB) -15198 2217        138     3088
# MR_Model_11 (intercept + study_setting + RoB)                                                 MR_Model_11 (intercept + study_setting + RoB) -16599 2268      -1263     3125
# MR_Model_12 (intercept + prev_GAD + ref_test + study_setting)                 MR_Model_12 (intercept + prev_GAD + ref_test + study_setting) -15534 2306       -198     3153
# MR_Model_13 (intercept + prev_GAD + ref_test + RoB)                                     MR_Model_13 (intercept + prev_GAD + ref_test + RoB) -14803 2273        534     3128
# MR_Model_14 (intercept + prev_GAD + study_setting + RoB)                           MR_Model_14 (intercept + prev_GAD + study_setting + RoB) -16805 2191      -1468     3069
# MR_Model_15 (intercept + ref_test + study_setting + RoB)                           MR_Model_15 (intercept + ref_test + study_setting + RoB) -16068 2213       -732     3085
# 
# 
# > summarize_model_comparison(comparison)
# 
# Model ranking by ELPD:
#   1. MR_Model_6  (intercept + prev_GAD + ref_test): ELPD = -14353 (2153) [BEST]
# 2. MR_Model_13 (intercept + prev_GAD + ref_test + RoB): ELPD = -14803 (2273), Δ = 534 (3128) 
# 3. MR_Model_2 (intercept + ref_test): ELPD = -14905 (2160), Δ = 431 (3047) 
# 4. MR_Model_10 (intercept + ref_test + RoB): ELPD = -15198 (2217), Δ = 138 (3088) 
# 5. MR_Model_1 (intercept + prev_GAD): ELPD = -15241 (2186), Δ = 96 (3066) 
# 6. MR_Model_0 (intercept-only - all MR models using FIXED-C + CS): ELPD = -15336 (2149), Δ = NA (NA) NA
# 7. MR_Model_12 (intercept + prev_GAD + ref_test + study_setting): ELPD = -15534 (2306), Δ = -198 (3153) 
# 8. MR_Model_8  (intercept + prev_GAD + RoB): ELPD = -15635 (2153), Δ = -299 (3042) 
# 9. MR_Model_4 (intercept + RoB): ELPD = -15739 (2334), Δ = -403 (3173) 
# 10. MR_Model_FULL (intercept + prev_GAD + ref_test + study_setting + RoB): ELPD = -15790 (2216), Δ = -454 (3087) 
# 11. MR_Model_9  (intercept + ref_test + study_setting): ELPD = -15909 (2238), Δ = -572 (3103) 
# 12. MR_Model_15 (intercept + ref_test + study_setting + RoB): ELPD = -16068 (2213), Δ = -732 (3085) 
# 13. MR_Model_3 (intercept + study_setting): ELPD = -16295 (2269), Δ = -958 (3125) 
# 14. MR_Model_7  (intercept + prev_GAD + study_setting): ELPD = -16360 (2123), Δ = -1024 (3021) 
# 15. MR_Model_11 (intercept + study_setting + RoB): ELPD = -16599 (2268), Δ = -1263 (3125) 
# 16. MR_Model_14 (intercept + prev_GAD + study_setting + RoB): ELPD = -16805 (2191), Δ = -1468 (3069) 
# 
# * = Significantly different from best model (|Δ| > 2*SE)




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
# For the simulation study, we generated data using parameters estimated from the best-performing model 
# identified via k-fold cross-validation on real data (fixed cutpoints, compound symmetry, reference test covariate only). 
# While this approach necessarily leads to the same model ranking in simulated data, it provides a realistic demonstration
# of our methodology under data-generating conditions supported by empirical evidence.



##
## ---- intercept + ref_test
##
beta_mu <- kfold_MR_Model_2$outs_kfold$fold_results[[1]]$beta_mu 
beta_mu %>% print(n = 100)
##
## ---- intercept + prev_GAD + ref_test
##
beta_mu <- kfold_MR_Model_6$outs_kfold$fold_results[[1]]$beta_mu 
beta_mu %>% print(n = 100)
##
## ---- intercept + ref_test + RoB
##
beta_mu <- kfold_MR_Model_10$outs_kfold$fold_results[[1]]$beta_mu 
beta_mu %>% print(n = 100)
##
## ---- intercept + prev_GAD + ref_test + RoB
##)
beta_mu <- kfold_MR_Model_13$outs_kfold$fold_results[[1]]$beta_mu 
beta_mu %>% print(n = 100)
##
##
source(file.path("R", "R_fns_other.R"))
##
summarise_beta_mu(  beta_mu = kfold_MR_Model_2$outs_kfold$fold_results[[1]]$beta_mu,
                    ##
                    test_names = c("GAD-2", "GAD-7", "HADS", "BAI",
                                   "OASIS", "STAI_S", "STAI_T", "PROMIS"),
                    ##
                    covariate_names = c("intercept", 
                                        "Ref_test_clean_SCID", "Ref_test_clean_Structured"),
                    ##
                    magnitude_threshold = 0.25,
                    skip_intercept = TRUE)
##
summarise_beta_mu(  beta_mu = kfold_MR_Model_6$outs_kfold$fold_results[[1]]$beta_mu,
                    ##
                    test_names = c("GAD-2", "GAD-7", "HADS", "BAI",
                                   "OASIS", "STAI_S", "STAI_T", "PROMIS"),
                    ##
                    covariate_names = c("intercept", 
                                        "logit_prev_GAD", 
                                        "Ref_test_clean_SCID", "Ref_test_clean_Structured"),
                    ##
                    magnitude_threshold = 0.25,
                    skip_intercept = TRUE)
##
summarise_beta_mu(  beta_mu = kfold_MR_Model_10$outs_kfold$fold_results[[1]]$beta_mu,
                    ##
                    test_names = c("GAD-2", "GAD-7", "HADS", "BAI",
                                   "OASIS", "STAI_S", "STAI_T", "PROMIS"),
                    ##
                    covariate_names = c("intercept", 
                                        "low_RoB_QUADAS_clean",
                                        "Ref_test_clean_SCID", "Ref_test_clean_Structured"),
                    ##
                    magnitude_threshold = 0.25,
                    skip_intercept = TRUE)
##
summarise_beta_mu(  beta_mu = kfold_MR_Model_13$outs_kfold$fold_results[[1]]$beta_mu,
                    ##
                    test_names = c("GAD-2", "GAD-7", "HADS", "BAI",
                                   "OASIS", "STAI_S", "STAI_T", "PROMIS"),
                    ##
                    covariate_names = c("intercept", 
                                        "logit_prev_GAD",
                                        "low_RoB_QUADAS_clean",
                                        "Ref_test_clean_SCID", "Ref_test_clean_Structured"),
                    ##
                    magnitude_threshold = 0.25,
                    skip_intercept = TRUE)
                  
# 
# model_summary_and_trace_obj$extract_params(params = c("kappa")) %>% print(n = 20)
# beta_mu <- model_summary_and_trace_obj$extract_params(params = c("beta_mu")) %>% print(n = 100)
## Min ESS w/ new priors after 500 iter + 16 chains for coeffs = 156 (~ 0)
## Min ESS w/ new priors after 500 iter + 16 chains for coeffs = 91 (but >1% divs)





outs_kfold$fold_results[[1]]$tibble_all %>% print(n=1000)
##
tibble_fixed <- kfold_MR_Model_2$outs_kfold$fold_results[[1]]$tibble_all %>%
  rename(parameter = variable)
# ##
MetaOrdDTA:::extract_params_from_tibble_batch(tibble = tibble_fixed,
                                              param_strings_vec = "log_lik_study")







# map_beta_coefficients( beta_tibble = beta_mu,
#                        covariate_names = dimnames(X[[1]][[1]])[[2]],
#                        test_names = test_names)
# 
# 
# beta_tibble = beta_mu
# covariate_names = dimnames(X[[1]][[1]])[[2]] ; covariate_names
# test_names = test_names
# 
# 
# beta_tibble
# 
# 
# beta_tibble <- fit$summary(
#   variables = "beta_mu",
#   mean, median, sd, mad,
#   ~quantile(.x, probs = c(0.025, 0.975)),
#   rhat, ess_bulk, ess_tail
# )
# 
# 
# colnames(X[[1]][[1]])
# 
# 
# MR_model_for_real_data
# 





require(ggplot2)

# Combine all test data with covariates
check_data <- bind_rows(
  lapply(seq_along(test_data_list), function(i) {
    test_data_list[[i]] %>%
      mutate(test = test_names[i]) %>%
      left_join(cov_data_per_test[[i]] %>% select(Study, study_setting, Ref_test_clean, low_RoB_QUADAS_clean), 
                by = "Study")
  })
)

# 1. Correlations: Se/Sp vs continuous covariates
check_data %>%
  group_by(Study, test) %>%
  slice(1) %>%  # One row per study-test
  ungroup() %>%
  summarise(
    Se_prev_cor = cor(Se, qlogis(prev), use = "complete.obs"),
    Sp_prev_cor = cor(Sp, qlogis(prev), use = "complete.obs")
  )

# # 2. Boxplots: Se/Sp by categorical covariates
# check_data %>%
#   group_by(Study, test) %>%
#   slice(1) %>%
#   ggplot(aes(x = Ref_test_clean, y = Se)) +
#   geom_boxplot() +
#   labs(title = "Se by Reference Test")
# 
# check_data %>%
#   group_by(Study, test) %>%
#   slice(1) %>%
#   ggplot(aes(x = study_setting, y = Se)) +
#   geom_boxplot() +
#   labs(title = "Se by Study Setting")

# 3. Quick t-tests / ANOVA
check_data %>%
  group_by(Study, test) %>%
  slice(1) %>%
  ungroup() %>%
  {t.test(Se ~ low_RoB_QUADAS_clean, data = .)}

# 4. Simple summary table
check_data %>%
  group_by(Study, test) %>%
  slice(1) %>%
  group_by(Ref_test_clean) %>%
  summarise(
    n = n(),
    Se_mean = mean(Se, na.rm = TRUE),
    Sp_mean = mean(Sp, na.rm = TRUE),
    .groups = 'drop'
  )


check_data %>%
  group_by(Study, test) %>%
  slice(1) %>%
  group_by(low_RoB_QUADAS_clean) %>%
  summarise(
    n = n(),
    Se_mean = mean(Se, na.rm = TRUE),
    Sp_mean = mean(Sp, na.rm = TRUE),
    .groups = 'drop'
  )


check_data %>%
  group_by(Study, test) %>%
  slice(1) %>%
  group_by(study_setting) %>%
  summarise(
    n = n(),
    Se_mean = mean(Se, na.rm = TRUE),
    Sp_mean = mean(Sp, na.rm = TRUE),
    .groups = 'drop'
  )




check_data %>%
  group_by(Study, test) %>%
  slice(1) %>%
  group_by(study_setting) %>%
  summarise(
    n = n(),
    prev_mean = mean(prev, na.rm = TRUE),
    prev_median = median(prev, na.rm = TRUE),
    .groups = 'drop'
  )




check_data %>%
  group_by(Study, test) %>%
  slice(1) %>%
  group_by(study_setting) %>%
  summarise(
    n = n(),
    Se_mean = mean(Se, na.rm = TRUE),
    Se_sd = sd(Se, na.rm = TRUE),
    Sp_mean = mean(Sp, na.rm = TRUE),
    Sp_sd = sd(Sp, na.rm = TRUE),
    .groups = 'drop'
  )




 

mean_prev_2x2 <- test_data_list[[1]] %>%
  group_by(Study) %>%
  slice(1) %>%
  pull(prev) %>%
  mean()

mean_prev_cov <- mean(plogis(cov_data_per_test[[1]]$logit_prev_GAD), na.rm = TRUE)

cat("target_condition:", target_condition, "\n")
cat("Mean prev from 2x2:", mean_prev_2x2, "\n")
cat("Mean prev from cov:", mean_prev_cov, "\n")










