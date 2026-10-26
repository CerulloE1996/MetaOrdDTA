



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
        ##
        # # Remove any lock files
        # unlink(file.path(local_pkg_dir, "00LOCK-MetaOrdDTA"), recursive = TRUE)
        #
        # # Clean compiled objects
        # pkgbuild::clean_dll(local_pkg_dir)
        #
        # # Then just load
        # devtools::load_all(local_pkg_dir)
}
##
require(MetaOrdDTA)
require(RcppParallel)
## require(BayesMVPtest)
##
require(dplyr)
##
# ## target_condition <- "AAD"
# 
#  ##
#  
#   
#   # n_index_tests <- 8
##
    # target_condition <- "AAD"  ; n_index_tests <- 8 ## REAL analysis
    target_condition <- "GAD"  ; n_index_tests <- 8 ## REAL analysis
    # target_condition <- "GAD" ;  n_index_tests <- 4
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
    data_Klaus <- tibble(read.csv(file.path(Klaus_dir, "DIAQANDI NWMA database updated 2025_11_25_2x2_tables.csv")))
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
        
        # Remove Eack 2008 - no covariate data available
        data_BAI <- data_BAI %>%
          filter(!(Author == "Eack" & Year == 2008))
        
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
        } else { 
            test_names = c("GAD-2", "GAD-7", "HADS", "BAI",
                           "OASIS", "STAI_S", "STAI_T", "PROMIS")
            ##
            test_data_list <- list(data_gad2, 
                                   data_gad7, 
                                   data_HADS, 
                                   data_BAI,
                                   data_OASIS,
                                   data_STAI_S,
                                   data_STAI_T,
                                   data_PROMIS_adj)
        }
   
    }
     


}
    
    
    ## fix prevs:
    # test_data_list <- fix_prevalences(test_data_list)
    
    
    
    

  
# quantile(data_GAD$prev, probs = c(0.05, 0.25, 0.50, 0.75, 0.95))
# 


{
    ##
    ## ---- Also extract covariate info (NOT for sim study - just for application section and/or Klaus' analysis):
    ##
    # cov_data_Klaus <- tibble(read.csv(file.path(Klaus_dir, "DIAQANDI NWMA database final 2025_05_28_descriptive_short.csv")))
    ## DIAQANDI NWMA database updated 2025_11_25_descriptive_short.csv
    cov_data_Klaus <- tibble(read.csv(file.path(Klaus_dir, "DIAQANDI NWMA database updated 2025_11_25_descriptive_short.csv")))
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
          unique(cov_data_Klaus_3$Ref_test_clean)
          ##
          unique(cov_data_Klaus_3$low_RoB_QUADAS_clean)
          unique(cov_data_Klaus_3$low_RoB_liberal_clean)
          ##
          cov_data_Klaus_3$logit_prev_AAD
          cov_data_Klaus_3$logit_prev_GAD
          ##
          require(MetaOrdDTA)
    }
    
    cov_data_Klaus_3 <- cov_data_Klaus_3 %>%
      mutate(Study = row_number())
    
    

}

 
    
    # GENERAL SOLUTION - NO HARDCODED NAMES
    {
      
         Klaus_fn_fix_prevalences <- function(test_data_list) {
            
            # Fix prevalence calculations for all tests in the list
            test_data_list_fixed <- lapply(test_data_list, 
                                           function(test_data) {
                                             test_data %>%
                                               select(-prev) %>%  # Remove incorrect prev column
                                               group_by(Study) %>%
                                               mutate(
                                                 # Calculate correct prevalence from 2x2 table
                                                 # Use first row's values since total N and diseased are same for all cutoffs
                                                 total_diseased = first(TP + FN),
                                                 total_n = first(TP + FP + TN + FN),
                                                 prev = total_diseased / total_n
                                               ) %>%
                                               select(-total_diseased, -total_n) %>%  # Remove helper columns
                                               ungroup()
                                           })
            
            return(test_data_list_fixed)
            
          }
          
          # Step 1: Get test-specific prevalences from test data
          Klaus_fn_get_test_prevalences <- function(test_data_list) {
            
                test_prev <- list()
                for (i in 1:length(test_data_list)) {
                  test_prev[[i]] <- test_data_list[[i]] %>%
                    select(Author, Year, prev) %>%
                    distinct() %>%
                    mutate(
                      test_num = i,
                      logit_prev = log(prev / (1 - prev))
                    )
                }
                
                return(test_prev)
            
          }
          
          # Step 2: Update covariate data with test-specific prevalences
          Klaus_fn_update_covariates <- function(cov_data, 
                                                 test_data_list) {
            
                  cov_per_test <- list()
                  
                  for (i in 1:length(test_data_list)) {
                    # Get studies in this test with their prevalences
                    studies_in_test <- test_data_list[[i]] %>%
                      group_by(Study) %>%
                      summarise(
                        prev_GAD = first(prev),
                        .groups = 'drop'
                      ) %>%
                      mutate(
                        logit_prev_GAD = log(prev_GAD / (1 - prev_GAD))
                      )
                    
                    # Filter cov_data to ONLY studies in this test
                    cov_per_test[[i]] <- cov_data %>%
                      filter(Study %in% studies_in_test$Study) %>%
                      select(-prev_GAD, -logit_prev_GAD) %>%
                      left_join(studies_in_test, by = "Study")
                  }
                  
                  return(cov_per_test)
            
          }
          
          # Step 3: Create indicator matrix
          Klaus_fn_create_indicator <- function( cov_data, 
                                                 test_data_list) {
            
                n_studies <- nrow(cov_data)
                n_tests <- length(test_data_list)
                indicator <- matrix(0, n_studies, n_tests)
                
                for (i in 1:n_tests) {
                  # Which rows in cov_data match studies in test i
                  studies_in_test <- cov_data %>%
                    mutate(row_id = row_number()) %>%
                    inner_join(
                      test_data_list[[i]] %>% select(Author, Year) %>% distinct(),
                      by = c("Author", "Year")
                    ) %>%
                    pull(row_id)
                  
                  indicator[studies_in_test, i] <- 1
                }
                
                return(indicator)
            
          }
          
          # # REPLACE Klaus_fn_create_indicator with this:
          Klaus_fn_create_indicator_from_study_ids <- function(cov_data,
                                                               test_data_list) {
            
            n_studies <- nrow(cov_data)
            n_tests <- length(test_data_list)
            indicator <- matrix(0, n_studies, n_tests)
            
            for (i in 1:n_tests) {
              # Get Study IDs that appear in this test's data
              studies_in_test <- test_data_list[[i]] %>%
                filter(!is.na(Study)) %>%
                pull(Study) %>%
                unique()
              
              # Mark those rows in indicator
              indicator[cov_data$Study %in% studies_in_test, i] <- 1
            }
            
            return(indicator)
          }
          
          ##
          # Robust function to match studies between datasets
          Klaus_fn_create_study_ids_fuzzy <- function( test_data_list, 
                                                       cov_data) {
                
                ##
                ## Clean names for matching
                ##
                clean_for_matching <- function(x) {
                  x %>%
                    trimws() %>%                    
                    gsub("\\s+", " ", .) %>%        
                    gsub("[[:punct:]]", " ", .) %>% 
                    gsub("\\s+", " ", .) %>%        
                    tolower() %>%                   
                    trimws()                         
                }
                ##
                ## Create clean version for matching in cov_data
                ##
                cov_data_match <- cov_data %>%
                                  mutate(Author_clean = clean_for_matching(Author))
                ##
                ## Process each test
                ##
                for (i in 1:length(test_data_list)) {
                  
                        # Check if Study column exists, if not create it
                        if (!"Study" %in% names(test_data_list[[i]])) {
                          test_data_list[[i]]$Study <- NA_integer_
                        }
                        ##
                        ## Get unique authors in this test
                        ##
                        test_authors <- test_data_list[[i]] %>%
                                        distinct(Author, Year) %>%
                                        mutate(Author_clean = clean_for_matching(Author))
                        ##
                        ## Match to cov_data - try Author + Year first
                        ##
                        matched <-  test_authors %>%
                                    left_join(
                                      cov_data_match %>% select( Study_match = Study,
                                                                 Author_clean,
                                                                 Year),
                                      by = c("Author_clean", "Year")
                                    )
                        
                        # For unmatched, try year-agnostic matching
                        unmatched <- matched %>% filter(is.na(Study_match))
                        
                        if (nrow(unmatched) > 0) {
                          
                                year_agnostic <-  unmatched %>%
                                                  select(Author, Author_clean) %>%
                                                  left_join(
                                                    cov_data_match %>% 
                                                      group_by(Author_clean) %>%
                                                      slice(1) %>%  # Take first match if multiple years
                                                      ungroup() %>%
                                                      select(Study_match2 = Study, Author_clean),
                                                    by = "Author_clean"
                                                  )
                                
                                # Update matched with year-agnostic matches
                                matched <-  matched %>%
                                            left_join(year_agnostic %>% select(Author, Study_match2), by = "Author") %>%
                                            mutate(Study_match = coalesce(Study_match, Study_match2)) %>%
                                            select(-Study_match2)
                          
                        }
                        
                        # Apply Study IDs to test data
                        test_data_list[[i]] <-  test_data_list[[i]] %>%
                                                left_join(
                                                  matched %>% select(Author, Year, Study_match), 
                                                  by = c("Author", "Year")
                                                ) %>%
                                                mutate(Study = coalesce(Study, Study_match)) %>%
                                                select(-Study_match)
                        
                        # Report any still unmatched
                        still_unmatched <-  test_data_list[[i]] %>%
                                            filter(is.na(Study)) %>%
                                            distinct(Author)
                        
                        if (nrow(still_unmatched) > 0) {
                            cat("\nTest", i, "- Truly unmatched authors (no covariate data):\n")
                            print(still_unmatched)
                        }
                  
                }
                
                return(test_data_list)
            
          }
    
    }
 
    
    
    
    {
 
        test_data_list <- Klaus_fn_fix_prevalences(test_data_list)
        
        # CRITICAL: Remove old (wrong) Study IDs before fuzzy matching
        for (i in 1:length(test_data_list)) {
          if ("Study" %in% names(test_data_list[[i]])) {
            test_data_list[[i]] <- test_data_list[[i]] %>% select(-Study)
          }
        }
        
        # 2. Match Study IDs (now using corrected cov_data_Klaus_3)
        test_data_list <- Klaus_fn_create_study_ids_fuzzy(test_data_list, cov_data_Klaus_3)
        
        # 3. Create derived objects
        test_prev <- Klaus_fn_get_test_prevalences(test_data_list)
        cov_data_per_test <- Klaus_fn_update_covariates(cov_data_Klaus_3, test_data_list)
        
        # Rebuild cov_data_per_test using Study IDs (not Author/Year)
        cov_data_per_test <- list()
        for (i in 1:length(test_data_list)) {
          study_ids <- test_data_list[[i]] %>%
            filter(!is.na(Study)) %>%
            pull(Study) %>%
            unique()
          
          cov_data_per_test[[i]] <- cov_data_Klaus_3 %>%
            filter(Study %in% study_ids)
        }
        
        # 4. Create indicator FROM STUDY IDs (not Author/Year matching)
        indicator_index_test_in_study <- Klaus_fn_create_indicator_from_study_ids(cov_data_Klaus_3, test_data_list)
    }
    
    # Check
    cat("Studies per test:", colSums(indicator_index_test_in_study), "\n")
    cat("Total studies:", nrow(cov_data_Klaus_3), "\n")
 
    # Look at the original raw data file for Study 32 (before any processing)
    data_gad2 %>% 
      filter(Study == 32) %>%
      select(Study, Author, Cutoff, TP, FP, TN, FN, prev)
    
    
    # ============================================================
    # CHECK: Are there STILL duplicate Study IDs in cov_data_Klaus_3?
    # ============================================================
    
    cov_data_Klaus_3 %>%
      group_by(Study) %>%
      filter(n() > 1) %>%
      arrange(Study) %>%
      select(Study, Author, Year) %>%
      print(n = 100)

    
    # ============================================================
    # STEP 4: Verify the fix
    # ============================================================
    
    for (i in 1:8) {
      expected <- sum(indicator_index_test_in_study[, i])
      actual <- nrow(cov_data_per_test[[i]])
      match <- ifelse(expected == actual, "✓", "❌")
      cat("Test", i, ":", match, "Indicator =", expected, ", cov_data =", actual, "\n")
    }
 
    
    
    # Verify
    for (i in 1:length(test_data_list)) {
      na_count <- sum(is.na(test_data_list[[i]]$Study))
      total_studies <- n_distinct(test_data_list[[i]]$Study, na.rm = TRUE)
      cat("Test", i, ": ", total_studies, "studies,", na_count, "NA rows\n")
    }
    
    # Check prevalence consistency within each study
    for (i in 1:length(test_data_list)) {
      cat("\n=== Test", i, "===\n")
      
      prev_check <- test_data_list[[i]] %>%
        group_by(Study) %>%
        summarise(
          n_cutoffs = n(),
          n_unique_prev = n_distinct(prev),
          min_prev = min(prev),
          max_prev = max(prev),
          prev_diff = max(prev) - min(prev)
        ) %>%
        filter(n_unique_prev > 1)  # Only show problems
      
      if (nrow(prev_check) > 0) {
        cat("❌ PROBLEM STUDIES with multiple prevalences:\n")
        print(prev_check)
        
        # Show details for first problem study
        problem_study <- prev_check$Study[1]
        cat("\nDetails for Study", problem_study, ":\n")
        test_data_list[[i]] %>%
          filter(Study == problem_study) %>%
          select(Study, Author, Cutoff, TP, FP, TN, FN, prev) %>%
          print()
      } else {
        cat("✅ All studies have consistent prevalences\n")
      }
    }
    
    # Also check Study 32 specifically (the one that had issues before)
    if (32 %in% test_data_list[[1]]$Study) {
      cat("\n=== Checking Study 32 specifically ===\n")
      test_data_list[[1]] %>%
        filter(Study == 32) %>%
        select(Study, Author, Cutoff, prev) %>%
        distinct() %>%
        print()
    }
    
    
    
    


outs_covs$MR_model_5


# }


print(categorical_covs)


cov_data_Klaus_3 %>% 
  select(-Year, -Country, -Female_pct,
         -low_RoB_QUADAS, -low_RoB_liberal,
         -index_test_admin_type,
         -DSM_type, -Ref_test,
         -logit_prev_AAD)  %>% 
  print(n = 20, width = 50000)


 
# Check the actual columns
names(cov_data_per_test[[1]])

# Check Ref_test_clean values
for (i in 1:8) {
  cat("\nTest", i, ":\n")
  cat("  Unique Ref_test values: ", unique(cov_data_per_test[[i]]$Ref_test_clean), "\n")
  cat("  NAs in Ref_test_clean: ", sum(is.na(cov_data_per_test[[i]]$Ref_test_clean)), "\n")
}



##
outs_covs_subset <- outs_covs$MR_model_5
##
X <- list(outs_covs_subset$X_nd, outs_covs_subset$X_d)



# CHECK 1: Structure - should be list[2][4]
cat("Output structure: List of", length(X), "x", 
    length(X[[1]]), "\n")
cat("Should be: List of 2 (ND/D) x 8 (tests)\n\n")

# CHECK 2: Dimensions match indicator matrix
for (i in 1:8) {
  expected_studies <- sum(indicator_index_test_in_study[,i])
  actual_nd <- nrow(X[[1]][[i]])
  actual_d <- nrow(X[[2]][[i]])
  
  cat("Test", i, ": Expected", expected_studies, "studies,",
      "Got ND=", actual_nd, ", D=", actual_d, "\n")
}

# CHECK 3: Test-specific prevalences are correct (WITH SCALING)
for (i in 1:2) {  # Bernstein is in tests 1 and 2
  # Get Bernstein's original values
  bernstein_orig_1 <- cov_data_per_test[[1]] %>% 
    filter(Author == "Bernstein") %>% 
    pull(logit_prev_GAD)
  bernstein_orig_2 <- cov_data_per_test[[2]] %>% 
    filter(Author == "Bernstein") %>% 
    pull(logit_prev_GAD)
  
  # Get scaled values from matrices (assuming Bernstein is row 3)
  bernstein_scaled_1 <- X[[1]][[1]][3, "logit_prev_GAD"]
  bernstein_scaled_2 <- X[[1]][[2]][3, "logit_prev_GAD"]
  
  cat("\nBernstein prevalences:\n")
  cat("  Original: Test 1 =", bernstein_orig_1, ", Test 2 =", bernstein_orig_2, "\n")
  cat("  Scaled:   Test 1 =", bernstein_scaled_1, ", Test 2 =", bernstein_scaled_2, "\n")
  cat("  Difference preserved?", 
      sign(bernstein_orig_1 - bernstein_orig_2) == sign(bernstein_scaled_1 - bernstein_scaled_2), "\n")
}


# CHECK 4: Verify column names are consistent
cat("\nColumns per test:\n")
for (i in 1:8) {
  cat("Test", i, ": ND has", ncol(X[[1]][[i]]), "cols, D has", 
      ncol(X[[2]][[i]]), "cols\n")
  cat("  Columns:", colnames(X[[1]][[i]]), "\n")
}

# CHECK 5: ND and D should be identical (since you used same data for both)
all_match <- TRUE
for (i in 1:8) {
  if (!all(X[[1]][[i]] == X[[2]][[i]])) {
    cat("WARNING: Test", i, "ND and D matrices differ!\n")
    all_match <- FALSE
  }
}
if (all_match) cat("\nND and D matrices are identical (as expected)\n")








# # ============================================================
# # DIAGNOSTIC: Find the mismatch between indicator and cov_data_per_test
# # ============================================================
# 
 cat("\n========== DIAGNOSTIC ==========\n\n")


cat("\n========== ACTUAL VERIFICATION ==========\n")

for (i in 1:8) {
  indicator_n <- sum(indicator_index_test_in_study[, i])
  test_data_n <- n_distinct(test_data_list[[i]]$Study, na.rm = TRUE)
  cov_per_test_n <- nrow(cov_data_per_test[[i]])
  
  all_match <- (indicator_n == test_data_n) && (indicator_n == cov_per_test_n)
  symbol <- ifelse(all_match, "✓", "❌")
  
  cat("Test", i, ":", symbol, 
      "Indicator =", indicator_n, 
      ", test_data =", test_data_n,
      ", cov_per_test =", cov_per_test_n, "\n")
}


cat("\n========== X MATRIX CHECK ==========\n")

for (i in 1:8) {
  indicator_n <- sum(indicator_index_test_in_study[, i])
  x_nd_n <- nrow(outs_covs$MR_model_5$X_nd[[i]])
  x_d_n <- nrow(outs_covs$MR_model_5$X_d[[i]])
  
  match <- (indicator_n == x_nd_n) && (indicator_n == x_d_n)
  symbol <- ifelse(match, "✓", "❌")
  
  cat("Test", i, ":", symbol, 
      "Indicator =", indicator_n, 
      ", X_nd =", x_nd_n,
      ", X_d =", x_d_n, "\n")
}

 
##
## ---------------------------------------------------------------------------------------------------------------------------------------------------
##
# MR_model_for_real_data <- "MR_model_6"
##
# # Model 1 (baseline meta-reg): 2 covariates
# intercept + prevalence
##
# # Model 2: 3 covariates  
# intercept +  SCID + Structured
##
# # Model 3: 3 covariates
# intercept +  setting_2 + setting_3
##
MR_model_for_real_data <- "MR_model_1"
# MR_model_for_real_data <- "MR_model_2"
# MR_model_for_real_data <- "MR_model_3"
# MR_model_for_real_data <- "MR_model_4"
##
## FULL MR model:
##
# MR_model_for_real_data <- "MR_model_5"
# MR_model_for_real_data <- "MR_model_5_excl_RoB"  
# ##
# ## 2-way MR models:
# ##
# MR_model_for_real_data <- "MR_model_6"  ## (intercept + prev_GAD + ref_test)    ### Best model using K-fold, on GAD data w/ 4 tests
# MR_model_for_real_data <- "MR_model_7"  ## (intercept + prev_GAD + study_setting)
# MR_model_for_real_data <- "MR_model_8"  ## (intercept + prev_GAD + RoB)
# MR_model_for_real_data <- "MR_model_9"  ## (intercept + ref_test + study_setting)
# MR_model_for_real_data <- "MR_model_10" ## (intercept + ref_test + RoB)
# MR_model_for_real_data <- "MR_model_11" ## (intercept + study_setting + RoB)
# # ##
# # ## 3-way MR models:
# # ##
# MR_model_for_real_data <- "MR_model_12" ## (intercept + prev_GAD + ref_test + study_setting)
# MR_model_for_real_data <- "MR_model_13" ## (intercept + prev_GAD + ref_test + RoB)
# MR_model_for_real_data <- "MR_model_14" ## (intercept + prev_GAD + study_setting + RoB)
# MR_model_for_real_data <- "MR_model_15" ## (intercept + ref_test + study_setting + RoB)
# ##
# ##
#
##
cov_data <- NULL
##
{
        ##
        if (MR_model_for_real_data == "intercept_only") {
          
                intercept_only <- TRUE
                X <- subset_X_covariates( X = X,
                                          covariates_to_keep = c("intercept"))
                baseline_case <- c(1)
              
          } else if (MR_model_for_real_data == "MR_model_1") { 
            
              
                intercept_only <- FALSE
                X <- subset_X_covariates( X = X, 
                                          covariates_to_keep = c("intercept", 
                                                                 "logit_prev_GAD"))
                baseline_case <- c(1, 
                                   0.0)
            
          } else if (MR_model_for_real_data == "MR_model_2") {
            
                intercept_only <- FALSE
                X <- subset_X_covariates( X = X, 
                                          covariates_to_keep = c("intercept", 
                                                                 "Ref_test_clean_SCID", "Ref_test_clean_Structured"))
                baseline_case <- c(1, 
                                   0, 0)
            
          } else if (MR_model_for_real_data == "MR_model_3") {
            
                intercept_only <- FALSE
                X <- subset_X_covariates( X = X, 
                                          covariates_to_keep = c("intercept", 
                                                                 "study_setting_2", "study_setting_3"))
                baseline_case <- c(1, 
                                   0, 1)
            
          } else if (MR_model_for_real_data == "MR_model_4") {
            
                intercept_only <- FALSE
                X <- subset_X_covariates( X = X,
                                          covariates_to_keep = c("intercept",
                                                                 "low_RoB_QUADAS_clean"))
                baseline_case <- c(1, 
                                   0)
            
          } else if (MR_model_for_real_data == "MR_model_5") {
            
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
            
          } else if (MR_model_for_real_data == "MR_model_5_excl_RoB") {
          
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
              
        } else if (MR_model_for_real_data == "MR_model_6") { ## (intercept + logit_prev_GAD + ref_test)
          
              intercept_only <- FALSE
              X <- subset_X_covariates( X = X,
                                        covariates_to_keep = c("intercept",
                                                               "logit_prev_GAD",
                                                               "Ref_test_clean_SCID", "Ref_test_clean_Structured"))
                                                          
              baseline_case <- c(1, 
                                 0.0,
                                 0, 0)
              
        } else if (MR_model_for_real_data == "MR_model_7") { ## (intercept + logit_prev_GAD + study_setting)
          
              intercept_only <- FALSE
              X <- subset_X_covariates( X = X,
                                        covariates_to_keep = c("intercept",
                                                               "logit_prev_GAD",
                                                               "study_setting_2", "study_setting_3"))
              baseline_case <- c(1,
                                 0.0, 
                                 0, 1)
          
        } else if (MR_model_for_real_data == "MR_model_8") {  ## (intercept + logit_prev_GAD + RoB)
              
              intercept_only <- FALSE
              X <- subset_X_covariates( X = X,
                                        covariates_to_keep = c("intercept",
                                                               "logit_prev_GAD",
                                                               "low_RoB_QUADAS_clean"))
              baseline_case <- c(1, 
                                 0.0, 
                                 0)
          
        } else if (MR_model_for_real_data == "MR_model_9") { ## (intercept + ref_test + study_setting)
          
              intercept_only <- FALSE
              X <- subset_X_covariates( X = X,
                                        covariates_to_keep = c("intercept",
                                                               "study_setting_2", "study_setting_3",
                                                               "Ref_test_clean_SCID", "Ref_test_clean_Structured"))
              baseline_case <- c(1, 
                                 0, 0,
                                 0, 1)
          
        } else if (MR_model_for_real_data == "MR_model_10") { ## (intercept + ref_test + RoB)
          
              intercept_only <- FALSE
              X <- subset_X_covariates( X = X,
                                        covariates_to_keep = c("intercept",
                                                               "low_RoB_QUADAS_clean",
                                                               "Ref_test_clean_SCID", "Ref_test_clean_Structured"))
              baseline_case <- c(1,
                                 0, 
                                 0, 0)
          
        } else if (MR_model_for_real_data == "MR_model_11") { ## (intercept + study_setting + RoB)
          
              intercept_only <- FALSE
              X <- subset_X_covariates( X = X,
                                        covariates_to_keep = c("intercept",
                                                               "low_RoB_QUADAS_clean",
                                                               "study_setting_2", "study_setting_3"))
              baseline_case <- c(1,
                                 0,
                                 0, 1)
          
        } else if (MR_model_for_real_data == "MR_model_12") { ## (intercept + logit_prev_GAD + ref_test + study_setting)
              
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
          
        } else if (MR_model_for_real_data == "MR_model_13") { ## (intercept + logit_prev_GAD + ref_test + RoB)
          
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
          
        } else if (MR_model_for_real_data == "MR_model_14") { ## (intercept + logit_prev_GAD + study_setting + RoB)
          
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
          
        } else if (MR_model_for_real_data == "MR_model_15") { ## (intercept + ref_test + study_setting + RoB)
              
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
        # X <- MetaOrdDTA:::R_fn_expand_covariates_to_all_studies(
        #   X,
        #   indicator_index_test_in_study = indicator_index_test_in_study)
        ##
        cov_data <- list()
        ##
        X_not_expanded <- X
        cov_data$X <- X
        ##
        cov_data$baseline_case_nd <- rep(list(baseline_case), n_index_tests)
        cov_data$baseline_case_d  <- rep(list(baseline_case), n_index_tests)



}
##
source(file.path(getwd(), "R", "R_fn_covariates.R"))
##
# Check if prevalences are properly calculated in test_data_list
for (i in 1:length(test_data_list)) {
  
      cat("\nTest", i, ":\n")
      missing_prev <- test_data_list[[i]] %>%
        filter(is.na(prev) | is.infinite(prev))
      cat("Rows with missing prev:", nrow(missing_prev), "\n")
      
      # Check unique studies and their prevalences
      study_prev <- test_data_list[[i]] %>%
        group_by(Study) %>%
        summarise(
          prev = first(prev),
          has_prev = !is.na(prev) & !is.infinite(prev)
        )
      cat("Studies with valid prev:", sum(study_prev$has_prev), "out of", nrow(study_prev), "\n")
  
}
##
# Check cov_data_per_test
cov_data_per_test[[1]] %>%
  filter(is.na(logit_prev_GAD) | logit_prev_GAD == 999) %>%
  select(Study, Author, prev_GAD, logit_prev_GAD)
##
# Are the Study IDs matching between datasets?
test_studies <- test_data_list[[1]] %>% distinct(Study) %>% pull(Study)
cov_studies <- cov_data_Klaus_3$Study
##
# Which studies in test 1 are NOT in cov_data?
setdiff(test_studies, cov_studies)
##
# Which studies in cov_data are NOT in test 1?
setdiff(cov_studies, test_studies)
##
cov_data_per_test[[1]] %>% print(width = 9000)
## 
# cov_data_per_test[[1]]$prev_AAD - 
# round(cov_data_per_test[[1]]$prev_GAD, 2)

## Process all tests for NMA
# # Score ranges for each test
# if (n_index_tests == 4) {
#       min_scores <- c(0, 0, 0, 0)
#       max_scores <- c(6, 21, 21, 63)
# } else {
#       min_scores <- c(0, 0, 0, 0,
#                       0, 0, 0, 0)
#       max_scores <- c(6, 21, 21, 63,
#                       20,
#                       80 - 20,
#                       80 - 20,
#                       40 - 8)
# }
#   
#   
##
# YOUR PERSONAL CODE
if (n_index_tests == 4) {
  test_names <- c("GAD-2", "GAD-7", "HADS", "BAI")
  max_scores <- c(6, 21, 21, 63)
  # test_data_list <- list(data_gad2, data_gad7, data_HADS, data_BAI)
} else {
  test_names <- c("GAD-2", "GAD-7", "HADS", "BAI", 
                  "OASIS", "STAI_S", "STAI_T", "PROMIS")
  max_scores <- c(6, 21, 21, 63, 
                  20, 60, 60, 32)
  # test_data_list <- list(data_gad2, data_gad7, data_HADS, data_BAI, 
  #                        data_OASIS, data_STAI_S, data_STAI_T, data_PROMIS)
}
##

 
# ============================================================
# CRITICAL: Filter test_data_list to only studies with covariate data
# ============================================================
cat("\n=== Filtering to studies with covariate data ===\n")
for (i in 1:length(test_data_list)) {
  before_n <- n_distinct(test_data_list[[i]]$Study)
  valid_studies <- cov_data_per_test[[i]]$Study
  test_data_list[[i]] <- test_data_list[[i]] %>%
    filter(Study %in% valid_studies)
  after_n <- n_distinct(test_data_list[[i]]$Study)
  
  if (before_n != after_n) {
    cat("Test", i, ": FILTERED", before_n, "→", after_n, "studies\n")
  } else {
    cat("Test", i, ": OK (", after_n, "studies)\n")
  }
}

# NOW prepare the NMA data
NMA_data_list <- prepare_NMA_data(
  test_data_list = test_data_list,
  test_names = test_names,
  max_scores = max_scores
)
 



##
## correct data entry errors (PERSONAL CODE):
if (target_condition == "GAD") {
      {
          ##
          ## ---- Correction for data entry error:
          ##
          for (t in 1:4) { 
            for (c in 1:2) {
              for (s in 1:76) {
                if (any(NMA_data_list$x[[t]][[c]][s, ] == 33)) {
                    print(paste("t = ", t))
                    print(paste("c = ", c))
                    print(paste("s = ", s))
                    ##
                    print(NMA_data_list$x[[t]][[c]][s, ])
                }
              }
            }
          }
          ## at: t = 4, c = 2, s = 1 + position/threshold 20:
          NMA_data_list$x[[4]][[2]][1, 20]
          ## Replace with 23:
          NMA_data_list$x[[4]][[2]][1, 20] <- 23
          ##
          NMA_data_list$x[[4]][[2]][1, 20]
      }
} else { 
  
  ##
  ## ---- Correction for data entry error:
  ##
  for (t in 1:8) { 
    for (c in 1:2) {
      for (s in 1:80) {
        if (any(NMA_data_list$x[[t]][[c]][s, ] == 87)) {
              print(paste("t = ", t))
              print(paste("c = ", c))
              print(paste("s = ", s))
              ##
              print(NMA_data_list$x[[t]][[c]][s, ])
        }
      }
    }
  }
}
##
# Build the real_data structure:
##
real_data <- list()
real_data$x <- NMA_data_list$x
real_data$indicator_index_test_in_study <- NMA_data_list$indicator_index_test_in_study
##
## Add covariates:
##
{
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
}
##
if (n_index_tests == 4) {
  saveRDS(object = real_data,
          file = file.path(local_pkg_dir, "temp_old", "Klaus_et_al_data",
                           paste0("real_data_4_tests", "_", "GAD")))
} else { 
  saveRDS(object = real_data,
          file = file.path(local_pkg_dir, "temp_old", "Klaus_et_al_data",
                           paste0("real_data_8_tests", "_", "GAD")))
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
##
## Study 7 appears in tests 1, 2, 3:
##
cov_data$X[[1]][[1]][7, ]  # Study 7 in test 1
cov_data$X[[1]][[2]][7, ]  # Study 7 in test 2  
cov_data$X[[1]][[3]][7, ]  # Study 7 in test 3
##
## Original data for Study 7:
##
cov_data_Klaus_3 %>% filter(Study == 7)
##
sum(X[[2]][[1]][, 3]) + 
sum(X[[2]][[2]][, 3]) + 
sum(X[[2]][[3]][, 3]) + 
sum(X[[2]][[4]][, 3]) 
##
sum(X[[2]][[1]][, 4]) + 
sum(X[[2]][[2]][, 4]) + 
sum(X[[2]][[3]][, 4]) + 
sum(X[[2]][[4]][, 4]) 
##
30/76
##
 

 

# ============================================================
# DIAGNOSTIC: Compare indicator matrices
# ============================================================

# The Klaus indicator (full - 138 rows)
cat("Klaus indicator dimensions:", dim(indicator_index_test_in_study), "\n")
cat("Klaus indicator studies per test:", colSums(indicator_index_test_in_study), "\n")

# Model indicator (reduced - 81 rows)
model_indicator <- create_indicator_index_test_in_study(real_data$x)
cat("\nModel indicator dimensions:", dim(model_indicator), "\n")
cat("Model indicator studies per test:", colSums(model_indicator), "\n")

# Do they match per-test?
cat("\nPer-test comparison:\n")
for (i in 1:8) {
  klaus_n <- sum(indicator_index_test_in_study[, i])
  model_n <- sum(model_indicator[, i])
  match <- ifelse(klaus_n == model_n, "✓", "❌")
  cat("Test", i, ":", match, "Klaus =", klaus_n, ", Model =", model_n, "\n")
}


# ============================================================
# CRITICAL: Do X matrices match the MODEL indicator?
# ============================================================

cat("\nX matrix rows vs Model indicator:\n")
for (i in 1:8) {
  model_n <- sum(model_indicator[, i])
  x_nd_n <- nrow(outs_covs$MR_model_5$X_nd[[i]])
  x_d_n <- nrow(outs_covs$MR_model_5$X_d[[i]])
  
  match <- ifelse(model_n == x_nd_n, "✓", "❌")
  cat("Test", i, ":", match, "Model =", model_n, ", X_nd =", x_nd_n, ", X_d =", x_d_n, "\n")
}



# ============================================================
# Find which studies are in covariate data but NOT in model data
# ============================================================

# Get studies from test_data_list (covariate source)
cov_studies_test2 <- test_data_list[[2]] %>% 
  select(Author, Year) %>% 
  distinct()

cov_studies_test4 <- test_data_list[[4]] %>% 
  select(Author, Year) %>% 
  distinct()

# Get studies from NMA_data_list (model source) 
# We need to check what prepare_NMA_data actually kept
cat("Test 2: Covariate has", nrow(cov_studies_test2), "studies, Model has", 
    sum(model_indicator[, 2]), "\n")
cat("Test 4: Covariate has", nrow(cov_studies_test4), "studies, Model has", 
    sum(model_indicator[, 4]), "\n")

# Check for studies with missing/invalid data that might get filtered
cat("\n=== Checking Test 2 for problematic data ===\n")
test_data_list[[2]] %>%
  group_by(Author, Year) %>%
  summarise(
    n_rows = n(),
    any_na_TP = any(is.na(TP)),
    any_na_FN = any(is.na(FN)),
    any_na_FP = any(is.na(FP)),
    any_na_TN = any(is.na(TN)),
    any_zero_total = any((TP + FN + FP + TN) == 0),
    .groups = "drop"
  ) %>%
  filter(any_na_TP | any_na_FN | any_na_FP | any_na_TN | any_zero_total) %>%
  print()

cat("\n=== Checking Test 4 for problematic data ===\n")
test_data_list[[4]] %>%
  group_by(Author, Year) %>%
  summarise(
    n_rows = n(),
    any_na_TP = any(is.na(TP)),
    any_na_FN = any(is.na(FN)),
    any_na_FP = any(is.na(FP)),
    any_na_TN = any(is.na(TN)),
    any_zero_total = any((TP + FN + FP + TN) == 0),
    .groups = "drop"
  ) %>%
  filter(any_na_TP | any_na_FN | any_na_FP | any_na_TN | any_zero_total) %>%
  print()






# Which studies are in model but NOT in covariates for Test 4?
model_studies_test4 <- test_data_list[[4]] %>%
  select(Author, Year) %>%
  distinct()

cov_studies_test4 <- cov_data_per_test[[4]] %>%
  select(Author, Year) %>%
  distinct()

# Studies in model but not in covariate data
model_studies_test4 %>%
  anti_join(cov_studies_test4, by = c("Author", "Year"))



# Check the RAW data before any processing
data_BAI %>%
  select(Author, Year) %>%
  distinct() %>%
  print()






# Rebuild NMA_data_list with current (filtered) test_data_list
NMA_data_list <- prepare_NMA_data(
  test_data_list = test_data_list,
  test_names = test_names,
  max_scores = max_scores
)

# Rebuild real_data
real_data$x <- NMA_data_list$x
##
indicator_index_test_in_study <- NMA_data_list$indicator_index_test_in_study
real_data$indicator_index_test_in_study <- NMA_data_list$indicator_index_test_in_study

# Now check again
model_indicator <- create_indicator_index_test_in_study(real_data$x)
cat("Model indicator studies per test:", colSums(model_indicator), "\n")









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







# Which study is in model but NOT in covariates for Test 2?
model_studies_test2 <- test_data_list[[2]] %>%
  select(Author, Year) %>%
  distinct()

cov_studies_test2 <- cov_data_per_test[[2]] %>%
  select(Author, Year) %>%
  distinct()

# This will show the problem study
model_studies_test2 %>%
  anti_join(cov_studies_test2, by = c("Author", "Year"))








# Check: Do these 6 studies have Study IDs in test_data_list?
test_data_list[[2]] %>%
  filter(Author %in% c("Ahn (+ Oh 2018)", "Ayers", "Micoulaud-Franchi ", 
                       "Pranckeviciene ", "Shih ", "Stockton")) %>%
  select(Author, Year, Study) %>%
  distinct()

# Check: Do they exist in cov_data_Klaus_3?
cov_data_Klaus_3 %>%
  filter(grepl("Ahn|Ayers|Micoulaud|Pranck|Shih|Stockton", Author, ignore.case = TRUE)) %>%
  select(Study, Author, Year)





# Rebuild cov_data_per_test using Study IDs (not Author/Year)
cov_data_per_test <- list()
for (i in 1:length(test_data_list)) {
  study_ids <- test_data_list[[i]] %>%
    filter(!is.na(Study)) %>%
    pull(Study) %>%
    unique()
  
  cov_data_per_test[[i]] <- cov_data_Klaus_3 %>%
    filter(Study %in% study_ids)
  
  cat("Test", i, ":", length(study_ids), "studies\n")
}



cat("\n=== FINAL VERIFICATION ===\n")
model_indicator <- create_indicator_index_test_in_study(real_data$x)

for (i in 1:8) {
    model_n <- sum(model_indicator[, i])
    cov_n <- nrow(cov_data_per_test[[i]])
    match <- ifelse(model_n == cov_n, "✓", "❌")
    cat("Test", i, ":", match, "Model =", model_n, ", cov_data =", cov_n, "\n")
}

 








##
source(file.path(getwd(), "R", "R_fn_covariates.R"))
validate_X_matrices(cov_data$X)
##
X
##
check_X_column_order(X)
check_X_column_order(cov_data$X)
##
## --------------------------------------------------------------------------------------------------------------------------------------------
##
# intercept_only <- TRUE
intercept_only <- FALSE
##
# cts <- box_cox <- softplus <- FALSE
# custom_file_name <- NULL ##  "DTA_NMA_Nyaga_Xu_RANDthr_kappa.stan"
# # ##
model_parameterisation <- "Xu"
##
source(file.path(getwd(), "R", "R_fns_misc.R"))
source(file.path(getwd(), "R", "R_fn_K_fold_CV.R"))
##
##
random_thresholds <- FALSE ; Dirichlet_random_effects_type <- "fixed"  ; compound_symmetry <- TRUE  ## Model A
# random_thresholds <- FALSE ; Dirichlet_random_effects_type <- "fixed"  ; compound_symmetry <- FALSE  ## Model B
# random_thresholds <- TRUE  ; Dirichlet_random_effects_type <- "kappa"  ; compound_symmetry <- TRUE   ## Model C
# random_thresholds <- TRUE  ; Dirichlet_random_effects_type <- "kappa"  ; compound_symmetry <- FALSE  ## Model D
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
                  init_lists_per_chain = NULL)
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
## ---- Sample:  ---------------------------------------------------------------------
##
{
    tictoc::tic("timer")
    ##
    seed <- 2
    ##
    n_chains <- 16
    n_burnin <- 500
    n_iter <- 500
    ##
    adapt_delta <- 0.80
    max_treedepth <- 10
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
      init_lists_per_chain = init_lists_per_chain
      ##
      # cmdstanr_args = final_stan_args
    )
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
##
# model_summary_and_trace_obj$extract_params(params = c("kappa")) %>% print(n = 20)
# ##
# intercept_only
# 
# min(model_summary_and_trace_obj$extract_params(params = c("log_lik_study"))$n_eff)
# 
# model_summary_and_trace_obj$get_efficiency_metrics()
# 
# model_summary_and_trace_obj$extract_params(params = c("beta_mu")) %>% print(n = 20)
# 
# ## w/ log(200/50), df = 15, log_sd = 0.75 priors:
# # A tibble: 8 × 10
# # param_group parameter   mean     sd `2.5%` `50%` `97.5%` n_eff  Rhat n_Rhat
# # <chr>       <chr>      <dbl>  <dbl>  <dbl> <dbl>   <dbl> <dbl> <dbl> <lgl> 
# #   1 kappa       kappa[1,1] 105.   24.8    65.9 102.    160.   1179  1.01 NA    
# # 2 kappa       kappa[2,1]  25.0   5.51   16.2  24.2    37.3  1379  1.01 NA    
# # 3 kappa       kappa[1,2] 202.   31.7   149.  199.    273.    169  1.07 NA    
# # 4 kappa       kappa[2,2]  97.5  20.3    66.1  94.5   144.    291  1.04 NA    
# # 5 kappa       kappa[1,3] 374.   71.1   262.  364.    540.    147  1.08 NA    
# # 6 kappa       kappa[2,3] 260.   97.6   131.  237.    502.    101  1.11 NA    
# # 7 kappa       kappa[1,4] 688.  240.    376.  636.   1293.    180  1.07 NA    
# # 8 kappa       kappa[2,4] 974.  560.    347.  814.   2412.     35  1.38 NA 
#
#
##
## ---- Save full model output:
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
    file_name_string <- paste0( "dummy_data_", dummy_data,
                                ##
                                "n_tests_", n_index_tests, ## new
                                "tar_cond_",  target_condition, ## new
                                ##
                                "intercept_only_", intercept_only,
                                "N_covs_", stan_data_list$n_covariates_max,
                                "param_", model_parameterisation,
                                "rand_thr_", random_thresholds,
                                "CS_", compound_symmetry,
                                "het_sigma_", 0)
    if (intercept_only == FALSE) {
      file_name_string <- paste0(file_name_string, 
                                 "_", MR_model_for_real_data)
    }
    ##
    # Save with overwrite check
    results_list_file <- file.path(output_dir, target_condition, paste0("seed_", seed, "_", file_name_string, "_applied_results.RDS"))
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
      results_list_file <- file.path(getwd(), output_dir, target_condition, paste0("seed_", seed, "_", file_name_string, "_applied_results.RDS"))
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

 
model_summary_and_trace_obj$get_HMC_diagnostic_info()


model_summary_and_trace_obj$get_efficiency_metrics()$time_total


compute_study_overlap_prop(indicator_index_test_in_study)
# Total studies: 76 
# Studies with 1 test: 38 
# Studies with 2+ tests: 26 
# Study overlap proportion: 0.342 
# [1] 0.3421053

priors$prior_beta_mu_SD






model_summary_and_trace_obj$extract_params(params = c("kappa")) %>% print(n = 20)
beta_mu <- model_summary_and_trace_obj$extract_params(params = c("beta_mu")) %>% print(n = 100)
## Min ESS w/ new priors after 500 iter + 16 chains for coeffs = 156 (~ 0)
## Min ESS w/ new priors after 500 iter + 16 chains for coeffs = 91 (but >1% divs)


R_fn_map_beta_coefficients(beta_tibble = beta_mu, 
                           covariate_names = dimnames(real_data$cov_data$X[[1]][[1]])[[2]], 
                           test_names = test_names, 
                           alpha = 0.05)



min(model_summary_and_trace_obj$extract_params(params = c("log_lik_study"))$n_eff)
##
## seed_1_dummy_data_0n_tests_8tar_cond_AADintercept_only_1N_covs_1param_Xurand_thr_FALSECS_1het_sigma_0_applied_results.RDS
## seed_2_dummy_data_0n_tests_8tar_cond_GADintercept_only_1N_covs_1param_Xurand_thr_FALSECS_TRUEhet_sigma_0_applied_results.RDS
Model_A <- readRDS(file.path(
  MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
  "old",
  # target_condition,
  "seed_1_dummy_data_0n_tests_8tar_cond_AADintercept_only_1N_covs_1param_Xurand_thr_FALSECS_1het_sigma_0_applied_results.RDS"))
##
Model_B <- readRDS(file.path(
  MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
  target_condition,
  "seed_1_dummy_data_0n_tests_8tar_cond_AADintercept_only_1N_covs_1param_Xurand_thr_FALSECS_0het_sigma_0_applied_results.RDS"))
##
Model_C <- readRDS(file.path(
  MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
  target_condition,
  "seed_1_dummy_data_0n_tests_8tar_cond_AADintercept_only_1N_covs_1param_Xurand_thr_TRUECS_1het_sigma_0_applied_results.RDS"))
##
Model_D <- readRDS(file.path(
  MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
  target_condition,
  "seed_1_dummy_data_0n_tests_8tar_cond_AADintercept_only_1N_covs_1param_Xurand_thr_TRUECS_0het_sigma_0_applied_results.RDS"))
##
model_list <- list(Model_A)




##
## For all AUCs:
##
for (i in 1:length(test_names)) {
  cat(test_names[i], "&", extract_param_row(model_list, paste0("AUC[", i, "]")), "\\\\\n")
}


##
## ---- Se/Sp table 1: -----------------------------------------------------------------------------------------------------------
##
##
## Get all Se/Sp sections at once:
##
all_se_sp <- create_all_se_sp_sections( model_list, 
                                        n_tests = 8, 
                                        max_n_thr = 63, 
                                        digits = 3)
cat(paste(all_se_sp, collapse = "\n"))
##
## Define thresholds
##
thr_screening_reference_tests <- list(
  "GAD-2"  = c(1, 2, 3),
  "GAD-7"  = c(8, 9, 10, 11),
  "HADS"   = c(8, 9, 10, 11),
  "BAI"    = c(14, 15, 16, 17, 18, 20),
  "OASIS"  = c(6, 7, 8, 9, 10, 11),
  "STAI_S" = c(64, 68, 72, 74, 75, 76) - 20,
  "STAI_T" = c(64, 68, 72, 74, 75, 76) - 20,
  "PROMIS" = c(18, 20, 21, 22, 23, 24)
)
##
## Optional: Define labels for key thresholds
##
thr_labels <- list(
  "GAD-2"  = c("(standard screening thr.)"),
  "GAD-7"  = c(" ", " ", "(moderate; standard)", " "),
  "HADS"   = c("(possible)", " ", " ", "(probable)"),
  "BAI"    = c(" ", " ", "(moderate)", " ", " ", " "),
  "OASIS"  = c(" ", " ", "(standard)", " ", " ", " "),
  "STAI_S" = c(" ", " ", " ", " ", " ", " "),
  "STAI_T" = c(" ", " ", " ", " ", " ", " "),
  "PROMIS" = c(" ", " ", " ", " ", " ", " ")
)
##
## Define test number mapping (if different from list order)
##
test_numbers <- c(
  "GAD-2"  = 1,
  "GAD-7"  = 2,
  "HADS"   = 3,
  "BAI"    = 4,
  "OASIS"  = 5,
  "STAI_S" = 6,
  "STAI_T" = 7,
  "PROMIS" = 8
)
##
## Generate LaTeX lines
##
latex_lines <- create_all_se_sp_sections(
  model_list = model_list,
  thresholds_list = thr_screening_reference_tests,
  threshold_labels = thr_labels,  # Optional - can be NULL
  digits = 2,
  n_tests = 8,
  max_n_thr = 63,
  test_numbers = test_numbers,  # Optional - can be NULL if order matches
  add_midrule_between_tests = TRUE
)
##
cat(latex_lines, sep = "\n")


##
## ---- Se/Sp table 2 (single/best model): ---------------------------------------------------------------------------------------------------
##
# Define thresholds and labels
# thr_screening_reference_tests <- list(
#   "GAD-2"  = c(1, 2, 3),
#   "GAD-7"  = c(8, 9, 10, 11),
#   "HADS"   = c(8, 9, 10, 11),
#   "BAI"    = c(14, 15, 16, 17, 18, 20),
#   "OASIS"  = c(6, 7, 8, 9, 10, 11),
#   "STAI_S" = c(44, 48, 52, 54, 55, 56),
#   "STAI_T" = c(44, 48, 52, 54, 55, 56),
#   "PROMIS" = c(18, 20, 21, 22, 23, 24)
# )
##
thr_screening_reference_tests <- list(
  "GAD-2"  = c(1, 2, 3),
  "GAD-7"  = c(4, 5, 6, 7),
  "HADS"   = c(5, 6, 7, 8),
  "BAI"    = c(2, 4, 6, 8, 10, 12),
  "OASIS"  = c(1, 2, 3, 4, 5, 6),
  "STAI_S" = c(52, 54, 56, 58, 60, 62, 64, 66, 68) - 20,
  "STAI_T" = c(52, 54, 56, 58, 60, 62, 64, 66, 68) - 20,
  "PROMIS" = c(12, 14, 16, 18, 20, 22)
)
##
# thr_labels <- list(
#   "GAD-2"  = c("", "", "(screening)"),
#   "GAD-7"  = c("", "", "(moderate)", ""),
#   "HADS"   = c("(possible)", "", "", "(probable)"),
#   "BAI"    = c("", "", "(moderate)", "", "", ""),
#   "OASIS"  = c("", "", "(standard)", "", "", ""),
#   "STAI_S" = c("", "", "", "", "", ""),
#   "STAI_T" = c("", "", "", "", "", ""),
#   "PROMIS" = c("", "", "", "", "", "")
# )
thr_labels <- list(
  "GAD-2"  = c("", "", ""),
  "GAD-7"  = c("", "", "", ""),
  "HADS"   = c("", "", "", ""),
  "BAI"    = c("", "", "", "", "", ""),
  "OASIS"  = c("", "", "", "", "", ""),
  "STAI_S" = c("", "", "", "", "", ""),
  "STAI_T" = c("", "", "", "", "", ""),
  "PROMIS" = c("", "", "", "", "", "")
)
##
# Generate table
latex_table <-  create_se_sp_table_compact(
                model_list = model_list,
                thresholds_list = thr_screening_reference_tests,
                threshold_labels = thr_labels,
                digits = 2,
                n_tests = 8,
                max_n_thr = 63,
                caption = "Sensitivity and specificity estimates from baseline NMA-DTA model",
                label = "tab:se_sp_baseline")
##
writeLines(latex_table)

 
 


##
## Get all estimates for table 2:
##
table2_lines <- create_table2_latex(model_list, digits = 3)
cat(paste(table2_lines, collapse = "\n"))




Model_A$model_summary_and_trace_obj$get_HMC_info()

 
Model_real <-  Model_A

# model_summary_and_trace_obj_REAL <- Model_real$model_summary_and_trace_obj

source(file.path("R", "R_fn_NMA_computations.R"))
##
AUC_outs_REAL <- R_fn_extract_NMA_AUC_summary( tibble_gq = Model_real$model_summary_and_trace_obj$get_summary_generated_quantities(),
                                               n_thr = n_thr,
                                               test_names = test_names,
                                                network = TRUE
                                               )


AUC_outs_REAL$auc
AUC_outs_REAL$auc_diff
AUC_outs_REAL$auc_pred

# AUC_outs_REAL$auc
# # A tibble: 8 × 7
# test test_name AUC_mean AUC_median AUC_lower AUC_upper AUC_sd
# <dbl> <chr>        <dbl>      <dbl>     <dbl>     <dbl>  <dbl>
#   1     1 GAD-2        0.827      0.828     0.785     0.864 0.0202
# 2     2 GAD-7        0.857      0.858     0.823     0.887 0.0161
# 3     3 HADS         0.859      0.860     0.824     0.891 0.0170
# 4     4 BAI          0.847      0.850     0.768     0.909 0.0362
# 5     5 OASIS        0.790      0.796     0.659     0.891 0.0595
# 6     6 STAI_S       0.781      0.787     0.636     0.888 0.0639
# 7     7 STAI_T       0.773      0.780     0.612     0.889 0.0713
# 8     8 PROMIS       0.822      0.827     0.701     0.910 0.0538
# > AUC_outs_REAL$auc_diff
# # A tibble: 28 × 10
# test1 test2 test1_name test2_name AUC_diff_mean AUC_diff_median AUC_diff_lower AUC_diff_upper AUC_diff_sd prob_test1_better
# <int> <int> <chr>      <chr>              <dbl>           <dbl>          <dbl>          <dbl>       <dbl>             <dbl>
#   1     1     2 GAD-2      GAD-7           -0.0302        -0.0300          -0.0709        0.00823      0.0201                NA
# 2     1     3 GAD-2      HADS            -0.0326        -0.0325          -0.0828        0.0165       0.0252                NA
# 3     2     3 GAD-7      HADS            -0.0200        -0.0221          -0.0925        0.0633       0.0400                NA
# 4     1     4 GAD-2      BAI              0.0367         0.0314          -0.0682        0.169        0.0606                NA
# 5     2     4 GAD-7      BAI              0.0459         0.0404          -0.0693        0.193        0.0670                NA
# 6     3     4 HADS       BAI              0.0542         0.0470          -0.0702        0.218        0.0740                NA
# 7     1     5 GAD-2      OASIS            0.00494       -0.000336        -0.0909        0.126        0.0556                NA
# 8     2     5 GAD-7      OASIS           -0.00245       -0.00264         -0.0456        0.0413       0.0220                NA
# 9     3     5 HADS       OASIS            0.0102         0.00784         -0.0586        0.0916       0.0384                NA
# 10     4     5 BAI        OASIS            0.0669         0.0618          -0.0356        0.198        0.0597                NA
# # ℹ 18 more rows
# # ℹ Use `print(n = ...)` to see more rows
# > AUC_outs_REAL$auc_pred
# # A tibble: 8 × 7
# test test_name AUC_mean AUC_median AUC_lower AUC_upper AUC_sd
# <dbl> <chr>        <dbl>      <dbl>     <dbl>     <dbl>  <dbl>
#   1     1 GAD-2        0.804      0.823     0.544     0.959  0.109
# 2     2 GAD-7        0.831      0.851     0.573     0.970  0.104
# 3     3 HADS         0.837      0.858     0.585     0.972  0.102
# 4     4 BAI          0.824      0.845     0.550     0.972  0.110
# 5     5 OASIS        0.768      0.792     0.446     0.960  0.135
# 6     6 STAI_S       0.754      0.780     0.411     0.955  0.143
# 7     7 STAI_T       0.746      0.771     0.405     0.950  0.144
# 8     8 PROMIS       0.801      0.826     0.496     0.970  0.125
# > 


##
## ---- Make LaTeX  AUC table:
##
AUC_latex_table_REAL <- create_AUC_latex_table( auc_outs = AUC_outs_REAL, 
                                                digits = 3)
cat(AUC_latex_table_REAL)












## sROC plot:
relevant_thresholds = list( "GAD-2" = c(1:6), 
                            "GAD-7" =  c(3:18), 
                            "HADS" =  c(3:18),
                            "BAI" = c(3:38),
                            "OASIS" = c(1:20),
                            "STAI_S" = c(3:40),
                            "STAI_I" = c(3:40),
                            "PROMIS" = c(1:32))


# relevant_thresholds <- NULL
tibbles_provided = FALSE
tibbles = NULL
##
base_size = 16
point_size = 3
linewidth = 0.5
##
##
conf_region_colour = "blue"
pred_region_colour = "green"
##
conf_region_alpha = 0.15
pred_region_alpha = 0.3

stan_model_file_name = model_summary_and_trace_obj_REAL$internal_obj$outs_stan_model_name$stan_model_file_name
stan_mod_samples     =  model_summary_and_trace_obj_REAL$internal_obj$outs_stan_sampling$stan_mod_samples

##
source(file.path("R", "R_fn_test_accuracy_plots_NMA.R"))
##
plots_REAL <- R_fn_sROC_plot_NMA(     stan_model_file_name = model_summary_and_trace_obj_REAL$internal_obj$outs_stan_model_name$stan_model_file_name,
                                      stan_mod_samples     =  model_summary_and_trace_obj_REAL$internal_obj$outs_stan_sampling$stan_mod_samples,
                                      ##
                                      test_names = test_names,
                                      ##
                                      relevant_thresholds = relevant_thresholds,,
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






# Load JUST the compiled library without attaching the package
library.dynam("MetaOrdDTA", package = "MetaOrdDTA")

# Verify it's loaded
"MetaOrdDTA" %in% names(getLoadedDLLs())

# Now check what's available
getDLLRegisteredRoutines("MetaOrdDTA")$.Call












# quick_kappa_prior(log_mean = log(200),
#                   log_sd = 1,
#                   df = 5)
# # outs <- induced_Dirichlet_ppc_plot(N = 5000,
# #                                    n_cat = 63,
# #                                    log_mean = log(200),
# #                                    log_sd = 1,
# #                                    df = 5)
# ## New prior; D- group:
# print(paste("New prior; d- group:"))
# quick_kappa_prior(log_mean = log(200),
#                   log_sd = 0.75,
#                   df = 15)
# # outs <- induced_Dirichlet_ppc_plot(N = 5000,
# #                                    n_cat = 63,
# #                                    log_mean = log(500),
# #                                    log_sd = 0.5,
# #                                    df = 15)
# ##
# ## Old prior; D+ group:
# print(paste("Old prior; d+ group:"))
# quick_kappa_prior(log_mean = log(50),
#                   log_sd = 1,
#                   df = 5)
# # outs <- induced_Dirichlet_ppc_plot(N = 5000,
# #                                    n_cat = 63,
# #                                    log_mean = log(50),
# #                                    log_sd = 1,
# #                                    df = 5)
# ## New prior; D+ group:
# print(paste("New prior; d+ group:"))
# quick_kappa_prior(log_mean = log(50),
#                   log_sd = 0.75,
#                   df = 15)
# # outs <- induced_Dirichlet_ppc_plot(N = 5000,
# #                                    n_cat = 63,
# #                                    log_mean = log(125),
# #                                    log_sd = 0.5,
# #                                    df = 15)
#
##
## ---- K-fold ---------------------------------------------------------------------------------------------------------------------------
##
source(file.path(getwd(), "R", "R_fns_misc.R"))
source(file.path(getwd(), "R", "R_fn_K_fold_CV.R"))
##
# intercept_only <- TRUE
intercept_only <- FALSE ; MR_model_for_real_data
##
random_thresholds <- FALSE ; Dirichlet_random_effects_type <- "fixed"  ; compound_symmetry <- TRUE  ## Model A
# random_thresholds <- FALSE ; Dirichlet_random_effects_type <- "fixed"  ; compound_symmetry <- FALSE  ## Model B
# random_thresholds <- TRUE ; Dirichlet_random_effects_type <- "kappa" ; compound_symmetry <- TRUE  ## Model C
# random_thresholds <- TRUE ; Dirichlet_random_effects_type <- "kappa"  ; compound_symmetry <- FALSE  ## Model D
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
                  init_lists_per_chain = NULL)
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
# priors$prior_beta_mu_SD
# for (kk in 1:8) {
#    priors$prior_beta_mu_SD[[kk]] <- array(1, dim = c(n_class, stan_data_list$n_index_tests))
# }
##
## ALways use same seed for every model:
set.seed(123, kind = "Mersenne-Twister") # set seed kind to ensure consistency
##
K <- 10
##
fold_assignments <- MetaOrdDTA::create_folds(
                    K = K,
                    indicator_index_test_in_study = real_data$indicator_index_test_in_study,
                    seed = 123,
                    network = TRUE)
##
# Inspect pattern to make sure its the same:
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
sum(fold_assignments ==  c(1L, 5L, 9L, 10L, 3L, 7L, 3L, 9L, 1L, 1L, 4L, 1L, 2L, 3L, 2L, 
                           7L, 1L, 2L, 10L, 10L, 8L, 2L, 8L, 1L, 3L, 1L, 1L, 5L, 9L, 1L, 
                           1L, 1L, 1L, 6L, 10L, 2L, 5L, 9L, 1L, 6L, 3L, 5L, 4L, 2L, 8L, 
                           9L, 3L, 5L, 7L, 4L, 3L, 2L, 6L, 7L, 1L, 6L, 7L, 8L, 4L, 1L, 4L, 
                           1L, 2L, 5L, 1L, 2L, 6L, 1L, 1L, 3L, 4L, 2L, 4L, 2L, 6L, 8L, 7L, 
                           3L, 8L, 1L, 10L))/n_studies
##
dput(fold_assignments)
##
## ---- And/or: run k-fold:
##
source(file.path(getwd(), "R", "R_fns_misc.R"))
source(file.path(getwd(), "R", "R_fn_K_fold_CV.R"))
##
n_chains <- 4
##
if (random_thresholds == FALSE) {
  n_burnin <- 500
  n_iter   <- 500
} else {
  # n_burnin <- 1000
  n_burnin <- 500
  n_iter   <- 1000
}
## Make sure have 4000 total iter for random-effects, and 1000-2000 total iter for fixed:
##
random_thresholds
n_iter*n_chains
n_chains*K
##
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
                                           cmdstanr_args = NULL,
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
                                           use_BayesMVP_for_faster_summaries = TRUE)
# #
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
                                  "het_sigma_", 0)
      ##
      if (intercept_only == FALSE) {
        file_name_string <- paste0(file_name_string, "_", MR_model_for_real_data)
      }
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
##
outs_kfold$min_ESS_vec
##
stan_data_list$n_covariates_max
##
compound_symmetry
random_thresholds
##
##
# model_names <- c("Intercept_only", "Prevalence", "Ref_test", "Setting")
##
## Note: accidently saved these w/ hetero_sigma = 1, eventhough it's 0 (+ redundant now anyway as removed from Stan models).
## seed_2_K_fold_dummy_data_0n_tests_8tar_cond_AADintercept_only_1N_covs_1param_Xurand_thr_TRUECS_1het_sigma_0_applied_results.RDS
##
kfold_Model_A <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                   target_condition,
                                   "seed_1_K_fold_dummy_data_0n_tests_8tar_cond_GADintercept_only_1N_covs_1param_Xurand_thr_FALSECS_TRUEhet_sigma_0_applied_results.RDS"))
kfold_Model_B <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                   target_condition,
                                   "seed_1_K_fold_dummy_data_0n_tests_8tar_cond_GADintercept_only_1N_covs_1param_Xurand_thr_FALSECS_FALSEhet_sigma_0_applied_results.RDS"))
kfold_Model_C <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                   target_condition,
                                   "seed_2_K_fold_dummy_data_0n_tests_8tar_cond_AADintercept_only_1N_covs_1param_Xurand_thr_TRUECS_TRUEhet_sigma_0_applied_results.RDS"))
kfold_Model_D <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                   target_condition,
                                   "seed_1_K_fold_dummy_data_0n_tests_8tar_cond_AADintercept_only_1N_covs_1param_Xurand_thr_TRUECS_FALSEhet_sigma_0_applied_results.RDS"))
##
kfold_results_list <- list(
  kfold_Model_A,
  kfold_Model_B
  # kfold_Model_C,
  # kfold_Model_D
)
##
# kfold_Model_A$outs_kfold$K
# kfold_Model_B$outs_kfold$K
# kfold_Model_C$outs_kfold$K
# kfold_Model_D$outs_kfold$K
##
model_names <- c(
  "Model_A (FIXED-C + CS)",
  "Model_B (FIXED-C + UN)"
  # "Model_C (RAND-C + CS)",
  # "Model_D (RAND-C + UN)"
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
                                        min_ess_threshold = 25)
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
## ---- Meta-regression k-fold:
##
kfold_Model_A <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                   target_condition,
                                   "seed_1_K_fold_dummy_data_0n_tests_8tar_cond_GADintercept_only_1N_covs_1param_Xurand_thr_FALSECS_TRUEhet_sigma_0_applied_results.RDS"))
##
kfold_MR_Model_1 <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                      target_condition,
                                      "seed_1_K_fold_dummy_data_0n_tests_8tar_cond_GADintercept_only_0N_covs_2param_Xurand_thr_FALSECS_TRUEhet_sigma_0_MR_model_1_applied_results.RDS"))
kfold_MR_Model_2 <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                      target_condition,
                                      "seed_1_K_fold_dummy_data_0n_tests_8tar_cond_AADintercept_only_0N_covs_3param_Xurand_thr_FALSECS_TRUEhet_sigma_0_MR_model_2_applied_results.RDS"))
kfold_MR_Model_3 <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                      target_condition,
                                      "seed_1_K_fold_dummy_data_0n_tests_8tar_cond_AADintercept_only_0N_covs_3param_Xurand_thr_FALSECS_TRUEhet_sigma_0_MR_model_3_applied_results.RDS"))
kfold_MR_Model_4 <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                      target_condition,
                                      "seed_1_K_fold_dummy_data_0n_tests_8tar_cond_AADintercept_only_0N_covs_2param_Xurand_thr_FALSECS_TRUEhet_sigma_0_MR_model_4_applied_results.RDS"))
##
kfold_MR_Model_5 <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                      target_condition,
                                      "seed_1_K_fold_dummy_data_0n_tests_8tar_cond_AADintercept_only_0N_covs_7param_Xurand_thr_FALSECS_TRUEhet_sigma_0_MR_model_5_applied_results.RDS"))
##
kfold_MR_Model_6 <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                      target_condition,
                                      "seed_1_K_fold_dummy_data_0n_tests_8tar_cond_AADintercept_only_0N_covs_4param_Xurand_thr_FALSECS_TRUEhet_sigma_0_MR_model_6_applied_results.RDS"))
kfold_MR_Model_7 <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                      target_condition,
                                      "seed_1_K_fold_dummy_data_0n_tests_8tar_cond_AADintercept_only_0N_covs_4param_Xurand_thr_FALSECS_TRUEhet_sigma_0_MR_model_7_applied_results.RDS"))
kfold_MR_Model_8 <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                      target_condition,
                                      "seed_1_K_fold_dummy_data_0n_tests_8tar_cond_AADintercept_only_0N_covs_4param_Xurand_thr_FALSECS_TRUEhet_sigma_0_MR_model_8_applied_results.RDS"))
kfold_MR_Model_9 <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                      target_condition,
                                      "seed_1_K_fold_dummy_data_0n_tests_8tar_cond_AADintercept_only_0N_covs_5param_Xurand_thr_FALSECS_TRUEhet_sigma_0_MR_model_9_applied_results.RDS"))
kfold_MR_Model_10 <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                       target_condition,
                                       "seed_1_K_fold_dummy_data_0n_tests_8tar_cond_AADintercept_only_0N_covs_4param_Xurand_thr_FALSECS_TRUEhet_sigma_0_MR_model_10_applied_results.RDS"))
kfold_MR_Model_11 <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                       target_condition,
                                       "seed_1_K_fold_dummy_data_0n_tests_8tar_cond_AADintercept_only_0N_covs_4param_Xurand_thr_FALSECS_TRUEhet_sigma_0_MR_model_11_applied_results.RDS"))
##
kfold_MR_Model_12 <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                       target_condition,
                                       "seed_1_K_fold_dummy_data_0n_tests_8tar_cond_AADintercept_only_0N_covs_6param_Xurand_thr_FALSECS_TRUEhet_sigma_0_MR_model_12_applied_results.RDS"))
kfold_MR_Model_13 <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                       target_condition,
                                       "seed_1_K_fold_dummy_data_0n_tests_8tar_cond_AADintercept_only_0N_covs_5param_Xurand_thr_FALSECS_TRUEhet_sigma_0_MR_model_13_applied_results.RDS"))
kfold_MR_Model_14 <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                       target_condition,
                                       "seed_1_K_fold_dummy_data_0n_tests_8tar_cond_AADintercept_only_0N_covs_5param_Xurand_thr_FALSECS_TRUEhet_sigma_0_MR_model_14_applied_results.RDS"))
kfold_MR_Model_15 <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                       target_condition,
                                       "seed_1_K_fold_dummy_data_0n_tests_8tar_cond_AADintercept_only_0N_covs_6param_Xurand_thr_FALSECS_TRUEhet_sigma_0_MR_model_15_applied_results.RDS"))
##
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
  kfold_MR_Model_1
  # kfold_MR_Model_2,
  # kfold_MR_Model_3,
  # kfold_MR_Model_4,
  # kfold_MR_Model_5,
  # ## 2-way:
  # kfold_MR_Model_6,
  # kfold_MR_Model_7,
  # kfold_MR_Model_8,
  # kfold_MR_Model_9,
  # kfold_MR_Model_10,
  # kfold_MR_Model_11,
  # ## 3-way:
  # kfold_MR_Model_12,
  # kfold_MR_Model_13,
  # kfold_MR_Model_14,
  # kfold_MR_Model_15
)
##
model_names <- c("MR_Model_0 (intercept-only - all MR models using FIXED-C + CS)",
                 "MR_Model_1 (intercept + prev_GAD)"
                 # "MR_Model_2 (intercept + ref_test)",
                 # "MR_Model_3 (intercept + study_setting)",
                 # "MR_Model_4 (intercept + RoB)",
                 # ##
                 # "MR_Model_FULL (intercept + prev_GAD + ref_test + study_setting + RoB)",
                 # ##
                 # ## First 2-way MR model as prev_GAD and ref_test gave best K-fold's out of univariate MR's:
                 # ##
                 # "MR_Model_6  (intercept + prev_GAD + ref_test)",
                 # "MR_Model_7  (intercept + prev_GAD + study_setting)",
                 # "MR_Model_8  (intercept + prev_GAD + RoB)",
                 # "MR_Model_9  (intercept + ref_test + study_setting)",
                 # "MR_Model_10 (intercept + ref_test + RoB)",
                 # "MR_Model_11 (intercept + study_setting + RoB)",
                 # ##
                 # ## 3-way models:
                 # ##
                 # "MR_Model_12 (intercept + prev_GAD + ref_test + study_setting)",
                 # "MR_Model_13 (intercept + prev_GAD + ref_test + RoB)",
                 # "MR_Model_14 (intercept + prev_GAD + study_setting + RoB)",
                 # "MR_Model_15 (intercept + ref_test + study_setting + RoB)"
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
comparison <- MetaOrdDTA:::compare_models_k_fold_CV( kfold_results_list = kfold_results_list,
                                        model_names = model_names,
                                        min_ess_threshold = 10)
##
summarise_model_comparison_k_fold_CV(comparison)
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













 
MR_Model_5_FULL <- readRDS(file.path(MetaOrdDTA::metaorddta_output_path("application_results_real_data"),
                                     "seed_1_dummy_data_0intercept_only_0N_covs_7param_Xurand_thr_TRUECS_0het_sigma_0_MR_model_5_applied_results.RDS"))
##
model_summary_and_trace_obj <- MR_Model_5_FULL$model_summary_and_trace_obj
n_covariates <- 7

# 
# n_covariates <- stan_data_list$n_covariates_max
# n_covariates
##
 

{
    n_thr_vec <- c(6, 21, 21, 63)
    max_thr <- max(n_thr_vec)
    max_cat <- max_thr + 1
    ##
    try({ 
      dir_cat_SDs_sigma <- model_summary_and_trace_obj$extract_params(params = c("dir_cat_SDs_sigma")) %>% print(n = 20)
    })
    try({ 
      C       <- model_summary_and_trace_obj$extract_params(params = c("C_mu")) %>% print(n = 20)
    })
    try({
      C       <- model_summary_and_trace_obj$extract_params(params = c("C_array")) %>% print(n = 20)
    })
    try({ 
      alpha   <- model_summary_and_trace_obj$extract_params(params = c("alpha")) %>% print(n = 20)
    })
    beta_mu <- model_summary_and_trace_obj$extract_params(params = c("beta_mu")) %>% print(n = 20)
    ##
    beta_sigma <- model_summary_and_trace_obj$extract_params(params = c("beta_sigma")) %>% print(20)
    log_beta_sigma_MU <- model_summary_and_trace_obj$extract_params(params = c("log_beta_sigma_MU")) %>% print(20)
    log_beta_sigma_SD <- model_summary_and_trace_obj$extract_params(params = c("log_beta_sigma_SD")) %>% print(20)
    ##
    beta_tau <- model_summary_and_trace_obj$extract_params(params = c("beta_tau")) %>% print(20)
    beta_corr <- model_summary_and_trace_obj$extract_params(params = c("beta_corr")) %>% print(20)
    beta_corr$`50%`
    ##
    length_1 <- 2*n_index_tests*max_cat
    dir_cat_SDs_sigma_array <- array(-999, dim = c(2, n_index_tests, max_cat))
    alpha_array <- array(-999, dim = c(2, n_index_tests, max_cat))
    ##
    length_2 <- 2*n_index_tests*max_thr
    C_array     <- array(-999, dim = c(2, n_index_tests, max_thr))
    ##
    length_3 <- n_index_tests*2*  n_covariates
    beta_mu_array <- array(-999, dim = c(n_index_tests, 2,   n_covariates))
    ##
    length_4 <- n_index_tests*2
    beta_sigma_array <- array(-999, dim = c(n_index_tests, 2))
    beta_tau_array   <- array(-999, dim = c(n_index_tests, 2))
    beta_corr_val <- beta_corr$`50%` ; beta_corr_val
    log_beta_sigma_MU_vals <- log_beta_sigma_MU$`50%` ; log_beta_sigma_MU_vals
    log_beta_sigma_SD_vals <- log_beta_sigma_SD$`50%` ; log_beta_sigma_SD_vals
    ##
}
for (i in 1:length_1) { 
    # Extract indices from parameter name using regex
    param_name <- alpha$parameter[i]
    ##
    # Extract the three indices [c,t,k]
    indices <- str_match(param_name, "alpha\\[(\\d+),(\\d+),(\\d+)\\]")[,2:4]
    ##
    if (!is.na(indices[1])) {
        c_idx <- as.numeric(indices[1])  # class index
        t_idx <- as.numeric(indices[2])  # test index
        k_idx <- as.numeric(indices[3])  # threshold index
        ##
        # Fill the array with the mean value
        alpha_array[c_idx, t_idx, k_idx] <- alpha$`50%`[i]
        dir_cat_SDs_sigma_array[c_idx, t_idx, k_idx] <- dir_cat_SDs_sigma$`50%`[i]
    }
}
for (i in 1:length_2) { 
    param_name <- C$parameter[i]
    indices <- str_match(param_name, "C_array\\[(\\d+),(\\d+),(\\d+)\\]")[,2:4]
    ##
    if (!is.na(indices[1])) {
      c_idx <- as.numeric(indices[1])  # class index
      t_idx <- as.numeric(indices[2])  # test index
      k_idx <- as.numeric(indices[3])  # threshold index
      ##
      C_array[c_idx, t_idx, k_idx] <- C$`50%`[i]
    }
}
for (i in 1:length_3) { 
  param_name <- beta_mu$parameter[i]
  indices <- str_match(param_name, "beta_mu\\[(\\d+),(\\d+),(\\d+)\\]")[,2:4]
  ##
  if (!is.na(indices[1])) {
    t_idx <- as.numeric(indices[1])  # test index
    c_idx <- as.numeric(indices[2])  # class index
    k_idx <- as.numeric(indices[3])  # cov index
    ##
    beta_mu_array[t_idx, c_idx, k_idx] <- beta_mu$`50%`[i]
  }
}
for (i in 1:length_4) { 
  param_name <- beta_tau$parameter[i]
  indices <- str_match(param_name, "beta_tau\\[(\\d+),(\\d+)\\]")[,2:3]
  ##
  if (!is.na(indices[1])) {
    t_idx <- as.numeric(indices[1])  # test index
    c_idx <- as.numeric(indices[2])  # class index
    ##
  ##  beta_sigma_array[t_idx, c_idx] <- beta_sigma$`50%`[i]
    beta_tau_array[t_idx, c_idx]   <- beta_tau$`50%`[i]
  }
}
  for (i in 1:2) { 
    param_name <- beta_sigma$parameter[i]
    indices <- str_match(param_name, "beta_sigma\\[(\\d+),(\\d+)\\]")[,2:3]
    ##
    if (!is.na(indices[1])) {
      t_idx <- as.numeric(indices[1])  # test index
      c_idx <- as.numeric(indices[2])  # class index
      ##
      beta_sigma_array[t_idx, c_idx] <- beta_sigma$`50%`[i]
      beta_tau_array[t_idx, c_idx]   <- beta_tau$`50%`[i]
    }
}
beta_sigma_array
beta_tau_array




beta_mu <- model_summary_and_trace_obj$extract_params(params = c("beta_mu")) %>% print(n = 100)
##
outs <- R_fn_map_beta_coefficients(beta_tibble = beta_mu, 
                                   covariate_names = dimnames(X[[1]][[1]])[[2]], 
                                   test_names = test_names, 
                                   alpha = 0.05)
outs$full_table %>% print(n = 100)




dput(round(C_array, 3))
dput(round(alpha_array, 3))
dput(round(beta_mu_array, 3))
dput(round(beta_sigma_array, 3))
dput(round(beta_tau_array, 3))
dput(round(beta_corr_val, 3))
# log_beta_sigma_MU_vals
# log_beta_sigma_SD_vals
##
beta_sigma
##
alpha %>% print(n = 100)
beta_mu %>% print(n = 100)
beta_sigma %>% print(n = 100)
beta_tau %>% print(n = 100)
beta_corr %>% print(n = 100)
# log_beta_sigma_MU %>% print(n = 100)
# log_beta_sigma_SD %>% print(n = 100)
##
model_summary_and_trace_obj$extract_params(params = c("rho")) %>% print(n = 20)
model_summary_and_trace_obj$extract_params(params = c("rho12")) %>% print(n = 20)
##
##

model_summary_and_trace_obj$extract_params(params = c("raw_dirichlet_phi_vec")) %>% print(n = 40)



# Non-diseased class means (beta_sigma[·,1])
non_diseased_means <- beta_sigma$`50%`[1:4]

# Diseased class means (beta_sigma[·,2])
diseased_means <-  beta_sigma$`50%`[5:8]

# Calculate CV for non-diseased
mean_nd <- mean(non_diseased_means)
sd_nd <- sd(non_diseased_means)
cv_nd <- (sd_nd / mean_nd) * 100

# Calculate CV for diseased
mean_d <- mean(diseased_means)
sd_d <- sd(diseased_means)
cv_d <- (sd_d / mean_d) * 100

cat("Non-diseased class (beta_sigma[·,1]):\n")
cat(sprintf("Values: %s\n", paste(non_diseased_means, collapse = ", ")))
cat(sprintf("Mean: %.3f\n", mean_nd))
cat(sprintf("SD: %.3f\n", sd_nd))
cat(sprintf("CV: %.1f%%\n", cv_nd))

cat("\nDiseased class (beta_sigma[·,2]):\n")
cat(sprintf("Values: %s\n", paste(diseased_means, collapse = ", ")))
cat(sprintf("Mean: %.3f\n", mean_d))
cat(sprintf("SD: %.3f\n", sd_d))
cat(sprintf("CV: %.1f%%\n", cv_d))



# ##
# ## ---- Plot sROC curve(s) (NMA - if WITH covariates):
# ##
# new_cov_data <- list()
# ##
# if (intercept_only) { 
#   baseline_case <- c(1)
# } else { 
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
                     0)   ## study_setting_3
# }
# new_cov_data$baseline_case_nd <- baseline_case
# new_cov_data$baseline_case_d  <- baseline_case
# ##
# ##
# if (model_parameterisation == "R&G") { 
#   new_cov_data$baseline_case    <- rep(list(cov_data_list$baseline_case), n_index_tests)
# } else { 
#   new_cov_data$baseline_case_nd <- rep(list(cov_data_list$baseline_case_nd), n_index_tests)
#   new_cov_data$baseline_case_d  <- rep(list(cov_data_list$baseline_case_d), n_index_tests)
# }
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


  model_summary_and_trace_obj$extract_params(params = c("log_lik_study"))  %>%  print(n = 100)


stan_mod_samples <- model_summary_and_trace_obj$internal_obj$outs_stan_sampling$stan_mod_samples
str(stan_mod_samples)
stan_mod_samples

  
  
  
  # Simple diagnostic for BAI without modifying Stan code
  # This extracts existing parameters to diagnose the AUC issue
  
# Extract key parameters
params_to_extract <- c("Se_baseline", "Sp_baseline", "AUC", "C_mu")
draws_array <- stan_mod_samples$draws(variables = params_to_extract, format = "array")

# Get dimensions
n_iter <- dim(draws_array)[1]
n_chains <- dim(draws_array)[2]

# Extract BAI (test 4) values from first iteration
iter <- 1
chain <- 1

# Find all Se_baseline parameters for test 4
all_params <- dimnames(draws_array)[[3]]
se_params <- grep("^Se_baseline\\[4,", all_params, value = TRUE)
n_thr_bai <- length(se_params)

cat("Number of thresholds for BAI:", n_thr_bai, "\n\n")

# Extract Se and Sp for BAI
se_bai <- numeric(n_thr_bai)
sp_bai <- numeric(n_thr_bai)

for (k in 1:n_thr_bai) {
  se_bai[k] <- draws_array[iter, chain, paste0("Se_baseline[4,", k, "]")]
  sp_bai[k] <- draws_array[iter, chain, paste0("Sp_baseline[4,", k, "]")]
}

# Find valid values
valid_idx <- which(se_bai > -0.5)
se_valid <- se_bai[valid_idx]
sp_valid <- sp_bai[valid_idx]

# Create simple diagnostic plot
par(mfrow = c(2, 2))

# 1. Se and Sp by threshold
plot(valid_idx, se_valid, type = "b", col = "red", 
     main = "BAI: Se and Sp by Threshold",
     xlab = "Threshold Index", ylab = "Value",
     ylim = c(0, 1))
lines(valid_idx, sp_valid, type = "b", col = "blue")
legend("right", c("Sensitivity", "Specificity"), col = c("red", "blue"), lty = 1)

# 2. ROC curve
fpr <- 1 - sp_valid
plot(fpr, se_valid, type = "b", col = "purple",
     main = "BAI: ROC Curve", 
     xlab = "False Positive Rate", ylab = "Sensitivity",
     xlim = c(0, 1), ylim = c(0, 1))
abline(0, 1, lty = 2)

# Add (0,0) and (1,1) points
points(c(0, 1), c(0, 1), col = "purple", pch = 16)

# 3. Check monotonicity
se_diffs <- diff(se_valid)
sp_diffs <- diff(sp_valid)

plot(valid_idx[-1], se_diffs, type = "h", 
     main = "BAI: Se Differences (should be ≤ 0)",
     xlab = "Threshold", ylab = "Se[k] - Se[k-1]")
abline(h = 0, col = "red", lty = 2)

plot(valid_idx[-1], sp_diffs, type = "h",
     main = "BAI: Sp Differences (should be ≥ 0)", 
     xlab = "Threshold", ylab = "Sp[k] - Sp[k-1]")
abline(h = 0, col = "red", lty = 2)

par(mfrow = c(1, 1))

# Calculate AUC manually
cat("\n=== AUC Calculation Diagnostic ===\n")

# Method 1: Simple trapezoidal rule with original order
auc_simple <- 0
fpr_with_bounds <- c(0, fpr, 1)
tpr_with_bounds <- c(0, se_valid, 1)

# Sort by FPR
ord <- order(fpr_with_bounds)
fpr_sorted <- fpr_with_bounds[ord]
tpr_sorted <- tpr_with_bounds[ord]

for (i in 1:(length(fpr_sorted)-1)) {
  width <- fpr_sorted[i+1] - fpr_sorted[i]
  height <- (tpr_sorted[i] + tpr_sorted[i+1]) / 2
  auc_simple <- auc_simple + width * height
}

# Get Stan's AUC
stan_auc <- draws_array[iter, chain, "AUC[4]"]

cat("Manual AUC calculation:", round(auc_simple, 4), "\n")
cat("Stan AUC for BAI:", round(stan_auc, 4), "\n")
cat("Difference:", round(auc_simple - stan_auc, 4), "\n\n")

# Diagnostic checks
cat("=== Diagnostic Checks ===\n")
cat("Se range:", range(se_valid), "\n")
cat("Sp range:", range(sp_valid), "\n")
cat("Se monotonic decreasing?", all(se_diffs <= 0.001), "\n")
cat("Sp monotonic increasing?", all(sp_diffs >= -0.001), "\n")
cat("Any NaN values?", any(is.nan(c(se_valid, sp_valid))), "\n")

# Check if ROC curve is inverted
if (auc_simple < 0.5 || stan_auc < 0.5) {
  cat("\n*** WARNING: AUC < 0.5 indicates inverted ROC curve! ***\n")
  cat("Possible causes:\n")
  cat("1. Sign error in Se/Sp calculation formulas\n") 
  cat("2. BAI scores interpreted opposite to other tests\n")
  cat("3. Cutpoints ordered in reverse\n")
  
  # Try calculating AUC with inverted Se
  se_inverted <- 1 - se_valid
  auc_inverted <- 0
  tpr_inv_with_bounds <- c(0, se_inverted, 1)
  
  ord_inv <- order(fpr_with_bounds)
  fpr_sorted_inv <- fpr_with_bounds[ord_inv]
  tpr_sorted_inv <- tpr_inv_with_bounds[ord_inv]
  
  for (i in 1:(length(fpr_sorted_inv)-1)) {
    width <- fpr_sorted_inv[i+1] - fpr_sorted_inv[i]
    height <- (tpr_sorted_inv[i] + tpr_sorted_inv[i+1]) / 2
    auc_inverted <- auc_inverted + width * height
  }
  
  cat("\nAUC with inverted sensitivity:", round(auc_inverted, 4), "\n")
  cat("This would match the expected pattern from your SROC plot.\n")
}

# Print cutpoint values to check ordering
cat("\n=== Cutpoint Values (first 10) ===\n")
c_nd_params <- grep("^C_mu\\[1,4,", all_params, value = TRUE)
c_d_params <- grep("^C_mu\\[2,4,", all_params, value = TRUE)

if (length(c_nd_params) >= 10) {
  cat("Non-diseased cutpoints:\n")
  for (k in 1:10) {
    c_val <- draws_array[iter, chain, paste0("C_mu[1,4,", k, "]")]
    cat(sprintf("  C[%d] = %.3f\n", k, c_val))
  }
  
  # Check if monotonic
  c_vals <- numeric(min(10, length(c_nd_params)))
  for (k in 1:length(c_vals)) {
    c_vals[k] <- draws_array[iter, chain, paste0("C_mu[1,4,", k, "]")]
  }
  cat("Cutpoints monotonic increasing?", all(diff(c_vals) >= -0.001), "\n")
}



###########


# Let's check how many valid points BAI actually has
# and look for the smoking gun

# First, let's see how many thresholds each test has
cat("Number of thresholds per test:\n")
for (t in 1:4) {
  se_params_t <- grep(paste0("^Se_baseline\\[", t, ","), 
                      dimnames(draws_array)[[3]], value = TRUE)
  cat(sprintf("Test %d: %d thresholds\n", t, length(se_params_t)))
}

# Now let's check if BAI has any peculiarities
# Extract all BAI Se/Sp values for first iteration
iter <- 1
chain <- 1

# Count how many are actually valid (> -0.5)
se_bai_all <- numeric(63)
sp_bai_all <- numeric(63)
for (k in 1:63) {
  se_bai_all[k] <- draws_array[iter, chain, paste0("Se_baseline[4,", k, "]")]
  sp_bai_all[k] <- draws_array[iter, chain, paste0("Sp_baseline[4,", k, "]")]
}

n_valid <- sum(se_bai_all > -0.5)
cat("\nBAI valid thresholds:", n_valid, "out of 63\n")

# Here's my theory: The AUC calculation in Stan might be including 
# the invalid (-1) values in some way

# Let's manually calculate what Stan should be getting
# Using ONLY the valid points

valid_idx <- which(se_bai_all > -0.5)
se_valid <- se_bai_all[valid_idx]
sp_valid <- sp_bai_all[valid_idx]
fpr_valid <- 1 - sp_valid

# Add boundary points
fpr_full <- c(0, fpr_valid, 1)
tpr_full <- c(0, se_valid, 1)

# Sort
ord <- order(fpr_full)
fpr_sorted <- fpr_full[ord]
tpr_sorted <- tpr_full[ord]

# Calculate AUC
auc_manual <- 0
for (i in 1:(length(fpr_sorted)-1)) {
  auc_manual <- auc_manual + 0.5 * (tpr_sorted[i] + tpr_sorted[i+1]) * 
    (fpr_sorted[i+1] - fpr_sorted[i])
}

cat("\nDiagnostic summary:\n")
cat("Valid thresholds:", n_valid, "\n")
cat("Total points for AUC (including 0,0 and 1,1):", length(fpr_full), "\n")
cat("Manual AUC:", round(auc_manual, 4), "\n")
cat("Stan AUC:", round(draws_array[iter, chain, "AUC[4]"], 4), "\n")
cat("Difference:", round(auc_manual - draws_array[iter, chain, "AUC[4]"], 4), "\n")

# Check for duplicates in FPR which might cause sorting issues
cat("\nChecking for duplicate FPR values:\n")
fpr_table <- table(round(fpr_valid, 6))
duplicates <- fpr_table[fpr_table > 1]
if (length(duplicates) > 0) {
  cat("Found", length(duplicates), "duplicate FPR values\n")
  print(duplicates)
} else {
  cat("No duplicate FPR values found\n")
}

# Theory: Maybe Stan's bubble sort has an issue with the array bounds
# when there are exactly 63 thresholds?
cat("\nArray indices check:\n")
cat("n_thr[4] =", 63, "\n")
cat("n_points = n_thr[4] + 2 =", 65, "\n")
cat("Bubble sort outer loop: 1 to", 64, "\n")
cat("Bubble sort inner loop: 1 to (65-i)\n")

# Another theory: Check if there's something special about the middle point
middle_idx <- round(n_valid / 2)
cat("\nMiddle threshold check:\n")
cat(sprintf("Threshold %d: Se = %.4f, Sp = %.4f, FPR = %.4f\n",
            middle_idx, se_valid[middle_idx], sp_valid[middle_idx], 
            fpr_valid[middle_idx]))






# ##
# (dput(round(C_nd, 7)))
# (dput(round(C_d, 7)))
# ##
# (dput(round(dirichlet_cat_SDs_sigma_nd, 7)))
# (dput(round(dirichlet_cat_SDs_sigma_d, 7)))
##
(dput(round(alpha_nd, 7)))
(dput(round(alpha_d, 7)))
##
##
##
# (dput(round(100*Se_vec, 3)))
# (dput(round(100*Fp_vec, 3)))
# (dput(round(100*Sp_vec, 3)))
##
##
##
(dput(signif(beta_mu_nd, 3)))
(dput(signif(beta_mu_d, 3)))
##

test
# Klaus_data_list$X_nd

##
##
Stan_mod_sample$summary(c("beta_mu")) %>% print(n = 100)
Stan_mod_sample$summary(c("beta_SD")) %>% print(n = 100)
Stan_mod_sample$summary(c("beta_Omega")) %>% print(n = 100)
##
test
##

####
#### ---------------------- REAL data + WITH covariates:
####
alpha %>% print(n = 100)
# A tibble: 512 × 10
param_group parameter       mean      sd  `2.5%`  `50%` `97.5%` n_eff   Rhat n_Rhat
<chr>       <chr>          <dbl>   <dbl>   <dbl>  <dbl>   <dbl> <dbl>  <dbl> <lgl> 
1 alpha       alpha[1,1,1]   20.4    8.41    8.50   18.8    40.9    129   1.05 NA    
2 alpha       alpha[2,1,1]    2.74   1.36    0.961   2.45    6.05   293   1.03 NA    
3 alpha       alpha[1,2,1]    2.94   0.621   1.90    2.87    4.27    87   1.08 NA    
4 alpha       alpha[2,2,1]    2.33   0.756   1.13    2.24    4.05   240   1.02 NA    
5 alpha       alpha[1,3,1]    7.08   2.37    3.57    6.72   12.6    173   1.04 NA    
6 alpha       alpha[2,3,1]   17.6   13.8     4.44   13.3    57.8     40   1.15 NA    
7 alpha       alpha[1,4,1]    2.51   1.16    0.956   2.27    5.27   343   1.02 NA    
8 alpha       alpha[2,4,1]    7.51   5.84    1.40    5.82   23.2    200   1.04 NA    
9 alpha       alpha[1,1,2]   24.2    6.68   13.4    23.4    39.2    374   1.02 NA    
10 alpha       alpha[2,1,2]    3.29   1.18    1.50    3.12    6.01   365   1.03 NA    
11 alpha       alpha[1,2,2]    4.41   0.815   3.02    4.35    6.10   160   1.04 NA    
12 alpha       alpha[2,2,2]    1.28   0.385   0.704   1.21    2.17   433   1.01 NA    
13 alpha       alpha[1,3,2]   11.8    3.35    6.38   11.3    19.2    217   1.03 NA    
14 alpha       alpha[2,3,2]    5.45   3.91    1.51    4.20   16.0     34   1.17 NA    
15 alpha       alpha[1,4,2]    4.04   1.82    1.54    3.69    8.44   272   1.03 NA    
16 alpha       alpha[2,4,2]    3.64   3.05    0.698   2.72   11.4    267   1.03 NA    
17 alpha       alpha[1,1,3]   25.8    6.07   15.9    25.2    39.7    833   1.01 NA    
18 alpha       alpha[2,1,3]    6.26   1.73    3.45    6.06   10.0    505   1.01 NA    
19 alpha       alpha[1,2,3]    6.97   1.13    4.98    6.87    9.36   355   1.03 NA    
20 alpha       alpha[2,2,3]    1.58   0.474   0.847   1.52    2.68   310   1.02 NA    
21 alpha       alpha[1,3,3]   17.8    4.58   10.3    17.3    28.0    264   1.03 NA    
22 alpha       alpha[2,3,3]   15.5   10.7     4.43   12.5    45.8     34   1.17 NA    
23 alpha       alpha[1,4,3]    5.55   2.53    2.18    5.04   11.9    219   1.04 NA    
24 alpha       alpha[2,4,3]    5.97   5.07    1.12    4.58   18.5    265   1.03 NA    
25 alpha       alpha[1,1,4]   11.2    3.00    6.49   10.8    18.0    388   1.01 NA    

# ℹ 412 more rows
# ℹ Use `print(n = ...)` to see more rows
> beta_mu %>% print(n = 100)
# A tibble: 56 × 10
param_group parameter           mean     sd  `2.5%`     `50%` `97.5%` n_eff  Rhat n_Rhat
<chr>       <chr>              <dbl>  <dbl>   <dbl>     <dbl>   <dbl> <dbl> <dbl> <lgl> 
1 beta_mu     beta_mu[1,1,1] -0.759    0.258  -1.28   -0.752    -0.271    184 1.05  NA    
2 beta_mu     beta_mu[2,1,1] -1.53     0.236  -1.97   -1.54     -1.07     212 1.04  NA    
3 beta_mu     beta_mu[3,1,1] -1.07     0.236  -1.54   -1.08     -0.592    312 1.02  NA    
4 beta_mu     beta_mu[4,1,1] -0.880    0.482  -1.81   -0.897     0.0945  1191 1.00  NA    
5 beta_mu     beta_mu[1,2,1]  0.504    0.337  -0.180   0.510     1.14     442 1.02  NA    
6 beta_mu     beta_mu[2,2,1] -0.186    0.297  -0.764  -0.185     0.413    707 1.01  NA    
7 beta_mu     beta_mu[3,2,1]  0.116    0.304  -0.487   0.119     0.703    287 1.03  NA    
8 beta_mu     beta_mu[4,2,1] -0.170    0.529  -1.21   -0.168     0.848   1985 1.00  NA    
9 beta_mu     beta_mu[1,1,2]  0.146    0.0650  0.0206  0.146     0.272   1949 1.00  NA    
10 beta_mu     beta_mu[2,1,2]  0.229    0.0911  0.0474  0.228     0.413   1614 1.00  NA    
11 beta_mu     beta_mu[3,1,2] -0.000176 0.0637 -0.121   0.000213  0.124   1723 1.00  NA    
12 beta_mu     beta_mu[4,1,2]  0.203    0.327  -0.479   0.214     0.810   1227 1.01  NA    
13 beta_mu     beta_mu[1,2,2]  0.347    0.118   0.114   0.348     0.570   2030 1.00  NA    
14 beta_mu     beta_mu[2,2,2]  0.346    0.126   0.0957  0.347     0.591   1558 1.01  NA    
15 beta_mu     beta_mu[3,2,2] -0.186    0.0783 -0.333  -0.188    -0.0289  2130 1.00  NA    
16 beta_mu     beta_mu[4,2,2]  0.170    0.346  -0.594   0.191     0.792   1748 1.01  NA    
17 beta_mu     beta_mu[1,1,3] -0.0175   0.187  -0.375  -0.0224    0.341   1622 1.00  NA    
18 beta_mu     beta_mu[2,1,3] -0.128    0.259  -0.638  -0.126     0.372   1606 1.00  NA    
19 beta_mu     beta_mu[3,1,3] -0.107    0.126  -0.357  -0.109     0.145   1573 1.00  NA    
20 beta_mu     beta_mu[4,1,3] -0.00861  0.999  -1.94    0.00414   1.93    5328 1.00  NA    
21 beta_mu     beta_mu[1,2,3] -0.329    0.331  -0.977  -0.329     0.314   1851 1.00  NA    
22 beta_mu     beta_mu[2,2,3] -0.477    0.337  -1.14   -0.479     0.168   1431 1.01  NA    
23 beta_mu     beta_mu[3,2,3] -0.113    0.152  -0.403  -0.116     0.197   2548 1.00  NA    
24 beta_mu     beta_mu[4,2,3]  0.0152   1.00   -2.00    0.0165    1.96    5543 0.999 NA    
25 beta_mu     beta_mu[1,1,4]  0.0823   0.197  -0.327   0.0849    0.458   1605 1.00  NA    

> beta_sigma %>% print(n = 100)
# A tibble: 8 × 10
param_group parameter        mean     sd `2.5%` `50%` `97.5%` n_eff  Rhat n_Rhat
<chr>       <chr>           <dbl>  <dbl>  <dbl> <dbl>   <dbl> <dbl> <dbl> <lgl> 
1 beta_sigma  beta_sigma[1,1] 0.220 0.0564 0.108  0.221   0.332   602  1.01 NA    
2 beta_sigma  beta_sigma[2,1] 0.236 0.0705 0.113  0.230   0.389   466  1.01 NA    
3 beta_sigma  beta_sigma[3,1] 0.215 0.0531 0.110  0.216   0.319   558  1.01 NA    
4 beta_sigma  beta_sigma[4,1] 0.224 0.0684 0.107  0.219   0.379   822  1.01 NA    
5 beta_sigma  beta_sigma[1,2] 0.228 0.0983 0.0812 0.210   0.460   436  1.02 NA    
6 beta_sigma  beta_sigma[2,2] 0.237 0.114  0.0798 0.212   0.518   342  1.02 NA    
7 beta_sigma  beta_sigma[3,2] 0.178 0.0678 0.0664 0.171   0.327  1129  1.01 NA    
8 beta_sigma  beta_sigma[4,2] 0.216 0.104  0.0763 0.197   0.447   977  1.01 NA    
> beta_tau %>% print(n = 100)
# A tibble: 8 × 10
param_group parameter       mean     sd  `2.5%`  `50%` `97.5%` n_eff  Rhat n_Rhat
<chr>       <chr>          <dbl>  <dbl>   <dbl>  <dbl>   <dbl> <dbl> <dbl> <lgl> 
1 beta_tau    beta_tau[1,1] 0.153  0.0887 0.00886 0.149    0.329   524  1.01 NA    
2 beta_tau    beta_tau[2,1] 0.368  0.109  0.124   0.373    0.567   443  1.01 NA    
3 beta_tau    beta_tau[3,1] 0.135  0.0786 0.00757 0.132    0.292   424  1.01 NA    
4 beta_tau    beta_tau[4,1] 0.297  0.284  0.0112  0.213    1.06   1342  1.00 NA    
5 beta_tau    beta_tau[1,2] 0.386  0.151  0.0830  0.390    0.670   449  1.02 NA    
6 beta_tau    beta_tau[2,2] 0.539  0.146  0.213   0.543    0.809   414  1.02 NA    
7 beta_tau    beta_tau[3,2] 0.0937 0.0703 0.00448 0.0802   0.259  1322  1.00 NA    
8 beta_tau    beta_tau[4,2] 0.319  0.298  0.0123  0.238    1.09   1673  1.00 NA    
> beta_corr %>% print(n = 100)
# A tibble: 1 × 10
param_group parameter  mean    sd `2.5%` `50%` `97.5%` n_eff  Rhat n_Rhat
<chr>       <chr>     <dbl> <dbl>  <dbl> <dbl>   <dbl> <dbl> <dbl> <lgl> 
1 beta_corr   beta_corr 0.356 0.318 -0.354 0.392   0.858   920  1.00 NA    
> log_beta_sigma_MU %>% print(n = 100)
# A tibble: 2 × 10
param_group       parameter             mean    sd `2.5%` `50%` `97.5%` n_eff  Rhat n_Rhat
<chr>             <chr>                <dbl> <dbl>  <dbl> <dbl>   <dbl> <dbl> <dbl> <lgl> 
1 log_beta_sigma_MU log_beta_sigma_MU[1] -1.51 0.251  -2.06 -1.49  -1.08    479  1.01 NA    
2 log_beta_sigma_MU log_beta_sigma_MU[2] -1.59 0.367  -2.38 -1.58  -0.928   575  1.01 NA    
> log_beta_sigma_SD %>% print(n = 100)
# A tibble: 2 × 10
param_group       parameter             mean    sd  `2.5%` `50%` `97.5%` n_eff  Rhat n_Rhat
<chr>             <chr>                <dbl> <dbl>   <dbl> <dbl>   <dbl> <dbl> <dbl> <lgl> 
1 log_beta_sigma_SD log_beta_sigma_SD[1] 0.177 0.136 0.00746 0.149   0.502  2075  1.00 NA    
2 log_beta_sigma_SD log_beta_sigma_SD[2] 0.243 0.167 0.0107  0.219   0.612  1688  1.00 NA    
> 






####
#### ---------------------- SIMULATED data + WITH covariates:
####
  
  ##
  > alpha %>% print(n = 100)
# A tibble: 512 × 10
param_group parameter       mean      sd  `2.5%`  `50%` `97.5%` n_eff   Rhat n_Rhat
<chr>       <chr>          <dbl>   <dbl>   <dbl>  <dbl>   <dbl> <dbl>  <dbl> <lgl> 
  1 alpha       alpha[1,1,1]   22.2    5.09   14.2    21.7    33.1    595   1.01 NA    
2 alpha       alpha[2,1,1]    2.85   0.702   1.73    2.79    4.40   751   1.01 NA    
3 alpha       alpha[1,2,1]    3.33   0.597   2.29    3.29    4.56   357   1.02 NA    
4 alpha       alpha[2,2,1]    3.10   0.634   2.04    3.04    4.44   375   1.02 NA    
5 alpha       alpha[1,3,1]    8.43   1.89    5.44    8.17   12.5    357   1.02 NA    
6 alpha       alpha[2,3,1]   22.4    7.89   11.4    21.0    39.2     20   1.32 NA    
7 alpha       alpha[1,4,1]    4.60   1.93    1.87    4.26    9.10   113   1.05 NA    
8 alpha       alpha[2,4,1]   11.0    4.80    4.45   10.3    21.8     33   1.18 NA    
9 alpha       alpha[1,1,2]   28.8    5.94   18.9    28.3    41.1    714   1.01 NA    
10 alpha       alpha[2,1,2]    4.02   0.833   2.61    3.93    5.78  1088   1.00 NA    
11 alpha       alpha[1,2,2]    4.88   0.853   3.39    4.81    6.63   418   1.02 NA    
12 alpha       alpha[2,2,2]    1.61   0.388   0.952   1.57    2.44   559   1.01 NA    
13 alpha       alpha[1,3,2]   11.3    2.51    7.09   11.0    16.8    421   1.01 NA    
14 alpha       alpha[2,3,2]    5.03   1.92    2.25    4.66    9.26    26   1.24 NA    
15 alpha       alpha[1,4,2]    4.85   1.97    1.99    4.51    9.13   131   1.05 NA    
16 alpha       alpha[2,4,2]    4.15   2.05    1.38    3.77    8.97    49   1.11 NA    
17 alpha       alpha[1,1,3]   34.5    6.79   23.0    34.0    48.1    683   1.01 NA    
18 alpha       alpha[2,1,3]    6.88   1.27    4.71    6.78    9.48  1003   1.01 NA    
19 alpha       alpha[1,2,3]    7.40   1.11    5.39    7.33    9.67   449   1.01 NA    
20 alpha       alpha[2,2,3]    1.63   0.358   1.03    1.59    2.40   508   1.01 NA    
21 alpha       alpha[1,3,3]   18.1    3.58   12.1    17.7    25.6    376   1.02 NA    
22 alpha       alpha[2,3,3]   18.7    6.72    9.19   17.5    33.6     20   1.33 NA    
23 alpha       alpha[1,4,3]   13.0    4.13    6.49   12.5    22.0     62   1.09 NA    
24 alpha       alpha[2,4,3]    5.22   2.56    1.73    4.74   10.6     52   1.11 NA    
25 alpha       alpha[1,1,4]   15.1    3.02   10.1    14.8    21.3    703   1.01 NA    

# ℹ 412 more rows
# ℹ Use `print(n = ...)` to see more rows
> beta_mu %>% print(n = 100)
# A tibble: 56 × 10
param_group parameter           mean     sd   `2.5%`     `50%` `97.5%` n_eff  Rhat n_Rhat
<chr>       <chr>              <dbl>  <dbl>    <dbl>     <dbl>   <dbl> <dbl> <dbl> <lgl> 
  1 beta_mu     beta_mu[1,1,1] -0.776    0.235  -1.21    -0.783    -0.325    974 1.00  NA    
2 beta_mu     beta_mu[2,1,1] -1.45     0.226  -1.88    -1.45     -1.03    1045 1.00  NA    
3 beta_mu     beta_mu[3,1,1] -0.230    0.649  -1.50    -0.223     0.979   2327 0.999 NA    
4 beta_mu     beta_mu[4,1,1] -0.0301   1.02   -1.94    -0.0744    1.89    4244 0.998 NA    
5 beta_mu     beta_mu[1,2,1]  0.436    0.340  -0.251    0.431     1.08    1314 0.998 NA    
6 beta_mu     beta_mu[2,2,1]  0.217    0.330  -0.401    0.203     0.848   1359 1.00  NA    
7 beta_mu     beta_mu[3,2,1]  0.340    0.324  -0.286    0.338     0.946   1269 1.01  NA    
8 beta_mu     beta_mu[4,2,1] -0.00413  0.875  -1.81    -0.0150    1.62    4303 0.999 NA    
9 beta_mu     beta_mu[1,1,2]  0.224    0.128  -0.0175   0.227     0.461   1633 1.00  NA    
10 beta_mu     beta_mu[2,1,2]  0.166    0.0818  0.00792  0.166     0.316   1989 1.00  NA    
11 beta_mu     beta_mu[3,1,2]  0.0303   0.339  -0.611    0.0261    0.686   1663 1.00  NA    
12 beta_mu     beta_mu[4,1,2]  0.0118   0.985  -1.95    -0.00983   1.89    4358 0.999 NA    
13 beta_mu     beta_mu[1,2,2]  0.162    0.128  -0.0898   0.163     0.409   1713 1.00  NA    
14 beta_mu     beta_mu[2,2,2]  0.323    0.106   0.118    0.323     0.529   1914 1.00  NA    
15 beta_mu     beta_mu[3,2,2] -0.470    0.195  -0.845   -0.475    -0.0764  1445 1.00  NA    
16 beta_mu     beta_mu[4,2,2]  0.00746  0.660  -1.33     0.0192    1.28    3603 0.998 NA    
17 beta_mu     beta_mu[1,1,3] -0.391    0.302  -0.989   -0.397     0.159   2728 1.00  NA    
18 beta_mu     beta_mu[2,1,3]  0.00654  0.246  -0.463   -0.000654  0.481   1908 1.00  NA    
19 beta_mu     beta_mu[3,1,3] -0.112    0.546  -1.17    -0.107     0.920   1898 1.00  NA    
20 beta_mu     beta_mu[4,1,3] -0.0125   1.01   -1.98    -0.00924   1.95    5164 0.997 NA    
21 beta_mu     beta_mu[1,2,3] -0.236    0.501  -1.22    -0.233     0.716   2367 0.999 NA    
22 beta_mu     beta_mu[2,2,3] -0.0461   0.465  -0.967   -0.0512    0.816   1768 1.00  NA    
23 beta_mu     beta_mu[3,2,3] -0.281    0.250  -0.767   -0.283     0.188   2087 1.00  NA    
24 beta_mu     beta_mu[4,2,3]  0.0117   1.02   -1.99     0.000640  1.96    4570 0.998 NA    
25 beta_mu     beta_mu[1,1,4]  0.0108   0.965  -1.90    -0.00527   1.82    6292 0.997 NA    

> beta_sigma %>% print(n = 100)
# A tibble: 8 × 10
param_group parameter        mean     sd `2.5%` `50%` `97.5%` n_eff  Rhat n_Rhat
<chr>       <chr>           <dbl>  <dbl>  <dbl> <dbl>   <dbl> <dbl> <dbl> <lgl> 
  1 beta_sigma  beta_sigma[1,1] 0.216 0.0574  0.104 0.216   0.326   322  1.03 NA    
2 beta_sigma  beta_sigma[2,1] 0.233 0.0730  0.107 0.227   0.381    95  1.06 NA    
3 beta_sigma  beta_sigma[3,1] 0.234 0.0722  0.111 0.227   0.383   190  1.04 NA    
4 beta_sigma  beta_sigma[4,1] 0.206 0.0581  0.101 0.202   0.321   382  1.02 NA    
5 beta_sigma  beta_sigma[1,2] 0.257 0.0840  0.114 0.249   0.440   296  1.02 NA    
6 beta_sigma  beta_sigma[2,2] 0.261 0.0819  0.122 0.251   0.429   267  1.03 NA    
7 beta_sigma  beta_sigma[3,2] 0.231 0.0551  0.117 0.231   0.336   395  1.02 NA    
8 beta_sigma  beta_sigma[4,2] 0.243 0.0704  0.111 0.242   0.390   458  1.01 NA    
> beta_tau %>% print(n = 100)
# A tibble: 8 × 10
param_group parameter      mean     sd  `2.5%` `50%` `97.5%` n_eff  Rhat n_Rhat
<chr>       <chr>         <dbl>  <dbl>   <dbl> <dbl>   <dbl> <dbl> <dbl> <lgl> 
  1 beta_tau    beta_tau[1,1] 0.153 0.0848 0.0152  0.147   0.324   313  1.03 NA    
2 beta_tau    beta_tau[2,1] 0.290 0.0987 0.0613  0.305   0.448    90  1.07 NA    
3 beta_tau    beta_tau[3,1] 0.288 0.102  0.0595  0.297   0.462   239  1.03 NA    
4 beta_tau    beta_tau[4,1] 0.126 0.0885 0.00732 0.113   0.326   593  1.01 NA    
5 beta_tau    beta_tau[1,2] 0.355 0.129  0.0947  0.359   0.602   306  1.02 NA    
6 beta_tau    beta_tau[2,2] 0.350 0.114  0.0771  0.364   0.538   269  1.03 NA    
7 beta_tau    beta_tau[3,2] 0.113 0.0753 0.00536 0.103   0.270   455  1.02 NA    
8 beta_tau    beta_tau[4,2] 0.224 0.135  0.0172  0.209   0.509   606  1.01 NA    
> beta_corr %>% print(n = 100)
# A tibble: 1 × 10
param_group parameter  mean    sd `2.5%` `50%` `97.5%` n_eff  Rhat n_Rhat
<chr>       <chr>     <dbl> <dbl>  <dbl> <dbl>   <dbl> <dbl> <dbl> <lgl> 
  1 beta_corr   beta_corr 0.319 0.281 -0.278 0.339   0.802   291  1.02 NA    
> log_beta_sigma_MU %>% print(n = 100)
# A tibble: 2 × 10
param_group       parameter             mean    sd `2.5%` `50%` `97.5%` n_eff  Rhat n_Rhat
<chr>             <chr>                <dbl> <dbl>  <dbl> <dbl>   <dbl> <dbl> <dbl> <lgl> 
  1 log_beta_sigma_MU log_beta_sigma_MU[1] -1.52 0.268  -2.12 -1.50  -1.09    211  1.04 NA    
2 log_beta_sigma_MU log_beta_sigma_MU[2] -1.42 0.266  -2.02 -1.41  -0.962   293  1.02 NA    
> log_beta_sigma_SD %>% print(n = 100)
# A tibble: 2 × 10
param_group       parameter             mean    sd  `2.5%` `50%` `97.5%` n_eff  Rhat n_Rhat
<chr>             <chr>                <dbl> <dbl>   <dbl> <dbl>   <dbl> <dbl> <dbl> <lgl> 
  1 log_beta_sigma_SD log_beta_sigma_SD[1] 0.188 0.137 0.00719 0.164   0.495  1290  1.00 NA    
2 log_beta_sigma_SD log_beta_sigma_SD[2] 0.183 0.136 0.0103  0.155   0.499  1078  1.01 NA    

# 










#### ------------- REAL data + NO covariates -----------------------------------------------------------------------

alpha %>% print(n = 50)
beta_mu %>% print(n = 100)
beta_sigma %>% print(n = 100)
beta_tau %>% print(n = 100)
beta_corr %>% print(n = 100)
 



> alpha %>% print(n = 50)
# A tibble: 512 × 10
param_group parameter     mean     sd `2.5%` `50%` `97.5%` n_eff  Rhat n_Rhat
<chr>       <chr>        <dbl>  <dbl>  <dbl> <dbl>   <dbl> <dbl> <dbl> <lgl> 
1 alpha       alpha[1,1,1] 20.6   9.38   8.23  18.7    43.4    228  1.03 NA    
2 alpha       alpha[2,1,1]  2.86  1.48   1.07   2.51    6.49   234  1.04 NA    
3 alpha       alpha[1,2,1]  2.92  0.614  1.85   2.88    4.19   206  1.04 NA    
4 alpha       alpha[2,2,1]  2.37  0.759  1.22   2.24    4.11   379  1.02 NA    
5 alpha       alpha[1,3,1]  6.66  2.24   3.33   6.29   11.8    148  1.06 NA    
6 alpha       alpha[2,3,1] 16.3   9.87   4.72  13.6    41.4     30  1.19 NA    
7 alpha       alpha[1,4,1]  2.63  1.29   0.951  2.36    5.86   319  1.03 NA    
8 alpha       alpha[2,4,1]  8.05  6.64   1.55   6.28   24.5    125  1.05 NA    
9 alpha       alpha[1,1,2] 23.9   7.17  12.7   23.0    40.1    341  1.02 NA    
10 alpha       alpha[2,1,2]  3.37  1.26   1.61   3.13    6.27   400  1.02 NA    
11 alpha       alpha[1,2,2]  4.38  0.800  2.94   4.32    6.03   480  1.03 NA    
12 alpha       alpha[2,2,2]  1.26  0.370  0.695  1.21    2.10   486  1.01 NA    
13 alpha       alpha[1,3,2] 11.3   3.39   6.10  10.8    19.1    174  1.05 NA    
14 alpha       alpha[2,3,2]  5.26  2.83   1.69   4.61   12.2     31  1.18 NA    
15 alpha       alpha[1,4,2]  4.27  2.02   1.53   3.84    9.09   238  1.04 NA    
16 alpha       alpha[2,4,2]  3.80  3.45   0.655  2.78   13.5    240  1.03 NA    
17 alpha       alpha[1,1,3] 25.2   6.14  14.9   24.5    38.9    810  1.01 NA    
18 alpha       alpha[2,1,3]  6.39  1.80   3.60   6.17   10.6    537  1.02 NA    
19 alpha       alpha[1,2,3]  6.93  1.14   4.95   6.85    9.29   650  1.01 NA    
20 alpha       alpha[2,2,3]  1.60  0.453  0.891  1.55    2.60   477  1.01 NA    
21 alpha       alpha[1,3,3] 17.4   4.70  10.1   16.8    28.1    185  1.05 NA    
22 alpha       alpha[2,3,3] 15.4   8.93   4.39  13.2    38.8     31  1.19 NA    
23 alpha       alpha[1,4,3]  5.77  2.76   2.14   5.15   12.6    211  1.05 NA    
24 alpha       alpha[2,4,3]  6.06  5.08   0.982  4.47   18.7    133  1.05 NA    
25 alpha       alpha[1,1,4] 10.8   2.80   6.44  10.4    17.2    503  1.01 NA    


beta_mu %>% print(n = 100)
# A tibble: 8 × 10
param_group parameter        mean    sd `2.5%`  `50%` `97.5%` n_eff  Rhat n_Rhat
<chr>       <chr>           <dbl> <dbl>  <dbl>  <dbl>   <dbl> <dbl> <dbl> <lgl> 
1 beta_mu     beta_mu[1,1,1] -0.614 0.230 -1.05  -0.617 -0.162    174  1.04 NA    
2 beta_mu     beta_mu[2,1,1] -1.45  0.138 -1.72  -1.45  -1.19      96  1.06 NA    
3 beta_mu     beta_mu[3,1,1] -0.768 0.126 -1.01  -0.767 -0.541     65  1.09 NA    
4 beta_mu     beta_mu[4,1,1] -1.53  0.249 -1.99  -1.55  -0.999   1006  1.01 NA    
5 beta_mu     beta_mu[1,2,1]  0.347 0.245 -0.136  0.341  0.818    281  1.04 NA    
6 beta_mu     beta_mu[2,2,1] -0.256 0.171 -0.582 -0.251  0.0730   304  1.03 NA    
7 beta_mu     beta_mu[3,2,1]  0.106 0.175 -0.224  0.109  0.439    143  1.04 NA    
8 beta_mu     beta_mu[4,2,1] -0.371 0.270 -0.886 -0.377  0.175   1311  1.00 NA    



beta_sigma %>% print(n = 100)
# A tibble: 8 × 10
param_group parameter        mean     sd `2.5%` `50%` `97.5%` n_eff  Rhat n_Rhat
<chr>       <chr>           <dbl>  <dbl>  <dbl> <dbl>   <dbl> <dbl> <dbl> <lgl> 
1 beta_sigma  beta_sigma[1,1] 0.266 0.0672  0.130 0.267   0.395   433  1.02 NA    
2 beta_sigma  beta_sigma[2,1] 0.291 0.0840  0.143 0.284   0.472   236  1.05 NA    
3 beta_sigma  beta_sigma[3,1] 0.249 0.0545  0.135 0.251   0.355   484  1.02 NA    
4 beta_sigma  beta_sigma[4,1] 0.272 0.0807  0.136 0.266   0.454   692  1.01 NA    
5 beta_sigma  beta_sigma[1,2] 0.324 0.121   0.138 0.306   0.587   311  1.02 NA    
6 beta_sigma  beta_sigma[2,2] 0.346 0.138   0.136 0.320   0.642   238  1.02 NA    
7 beta_sigma  beta_sigma[3,2] 0.274 0.0773  0.135 0.272   0.431   693  1.01 NA    
8 beta_sigma  beta_sigma[4,2] 0.305 0.112   0.128 0.290   0.560   602  1.01 NA    



beta_tau %>% print(n = 100)
# A tibble: 8 × 10
param_group parameter      mean     sd  `2.5%` `50%` `97.5%` n_eff  Rhat n_Rhat
<chr>       <chr>         <dbl>  <dbl>   <dbl> <dbl>   <dbl> <dbl> <dbl> <lgl> 
1 beta_tau    beta_tau[1,1] 0.218 0.106  0.0202  0.224   0.415   382  1.03 NA    
2 beta_tau    beta_tau[2,1] 0.379 0.125  0.0980  0.390   0.599   134  1.06 NA    
3 beta_tau    beta_tau[3,1] 0.122 0.0774 0.00654 0.114   0.280   330  1.02 NA    
4 beta_tau    beta_tau[4,1] 0.341 0.253  0.0216  0.289   0.969  1742  1.01 NA    
5 beta_tau    beta_tau[1,2] 0.458 0.168  0.0877  0.470   0.768   359  1.02 NA    
6 beta_tau    beta_tau[2,2] 0.538 0.180  0.117   0.559   0.836   223  1.02 NA    
7 beta_tau    beta_tau[3,2] 0.141 0.0950 0.00605 0.129   0.345  1037  1.00 NA    
8 beta_tau    beta_tau[4,2] 0.320 0.269  0.0125  0.257   0.974  1645  1.00 NA    


beta_corr %>% print(n = 100)
# A tibble: 1 × 10
param_group parameter  mean    sd  `2.5%` `50%` `97.5%` n_eff  Rhat n_Rhat
<chr>       <chr>     <dbl> <dbl>   <dbl> <dbl>   <dbl> <dbl> <dbl> <lgl> 
1 beta_corr   beta_corr 0.475 0.231 -0.0153 0.490   0.864   943  1.00 NA    



#### ------------- SIMULATED data + NO covariates -----------------------------------------------------------------------

> alpha %>% print(n = 100)
# A tibble: 512 × 10
param_group parameter       mean      sd  `2.5%`  `50%` `97.5%` n_eff   Rhat n_Rhat
<chr>       <chr>          <dbl>   <dbl>   <dbl>  <dbl>   <dbl> <dbl>  <dbl> <lgl> 
1 alpha       alpha[1,1,1]   15.7    6.15    7.05   14.6    29.8     39   1.16 NA    
2 alpha       alpha[2,1,1]    3.58   1.49    1.61    3.27    7.23   135   1.05 NA    
3 alpha       alpha[1,2,1]    4.39   1.03    2.75    4.26    6.61    79   1.09 NA    
4 alpha       alpha[2,2,1]    3.44   0.921   2.00    3.32    5.53   123   1.08 NA    
5 alpha       alpha[1,3,1]    5.41   1.50    3.14    5.23    8.73    49   1.13 NA    
6 alpha       alpha[2,3,1]   14.8    6.66    5.77   13.9    29.6     45   1.14 NA    
7 alpha       alpha[1,4,1]    4.90   1.46    2.66    4.71    8.21   258   1.02 NA    
8 alpha       alpha[2,4,1]    9.16   3.97    3.54    8.50   18.7     41   1.15 NA    
9 alpha       alpha[1,1,2]   24.9    5.98   14.9    24.3    36.6     90   1.08 NA    
10 alpha       alpha[2,1,2]    3.70   1.08    2.06    3.54    6.18   230   1.03 NA    
11 alpha       alpha[1,2,2]    5.35   1.07    3.48    5.25    7.61   136   1.05 NA    
12 alpha       alpha[2,2,2]    2.05   0.561   1.18    1.96    3.30    99   1.07 NA    
13 alpha       alpha[1,3,2]    8.75   2.17    5.35    8.45   13.6    176   1.05 NA    
14 alpha       alpha[2,3,2]    5.23   2.62    2.04    4.53   11.8     45   1.13 NA    
15 alpha       alpha[1,4,2]    5.95   1.70    3.33    5.74    9.76   303   1.02 NA    
16 alpha       alpha[2,4,2]    3.77   1.81    1.25    3.46    7.96    44   1.13 NA    
17 alpha       alpha[1,1,3]   32.1    6.25   21.1    31.6    44.5    440   1.02 NA    
18 alpha       alpha[2,1,3]    7.37   1.42    5.05    7.24   10.3    671   1.02 NA    
19 alpha       alpha[1,2,3]    7.93   1.33    5.63    7.79   10.7    166   1.05 NA    
20 alpha       alpha[2,2,3]    3.23   0.765   1.96    3.15    4.83    95   1.08 NA    
21 alpha       alpha[1,3,3]   17.9    4.01   11.4    17.4    26.4    209   1.04 NA    
22 alpha       alpha[2,3,3]   13.3    5.71    5.28   12.3    26.5     49   1.12 NA    
23 alpha       alpha[1,4,3]    9.08   2.43    5.15    8.73   14.3    276   1.03 NA    
24 alpha       alpha[2,4,3]    7.65   3.89    2.50    6.84   17.1     37   1.16 NA    
25 alpha       alpha[1,1,4]   16.1    3.63   10.0    15.8    23.8    192   1.03 NA    


beta_mu %>% print(n = 100)
# A tibble: 8 × 10
param_group parameter         mean    sd `2.5%`   `50%` `97.5%` n_eff  Rhat n_Rhat
<chr>       <chr>            <dbl> <dbl>  <dbl>   <dbl>   <dbl> <dbl> <dbl> <lgl> 
1 beta_mu     beta_mu[1,1,1] -0.854  0.226 -1.28  -0.855  -0.426     35  1.16 NA    
2 beta_mu     beta_mu[2,1,1] -1.30   0.129 -1.53  -1.29   -1.05      72  1.09 NA    
3 beta_mu     beta_mu[3,1,1] -1.02   0.138 -1.28  -1.02   -0.762     22  1.29 NA    
4 beta_mu     beta_mu[4,1,1] -1.58   0.181 -1.91  -1.59   -1.23     198  1.04 NA    
5 beta_mu     beta_mu[1,2,1]  0.175  0.259 -0.333  0.179   0.652    137  1.05 NA    
6 beta_mu     beta_mu[2,2,1]  0.0827 0.132 -0.164  0.0817  0.330    115  1.07 NA    
7 beta_mu     beta_mu[3,2,1]  0.0726 0.141 -0.202  0.0754  0.330     60  1.11 NA    
8 beta_mu     beta_mu[4,2,1] -0.399  0.169 -0.728 -0.403  -0.0761   280  1.02 NA    


beta_sigma %>% print(n = 100)
# A tibble: 8 × 10
param_group parameter        mean     sd `2.5%` `50%` `97.5%` n_eff  Rhat n_Rhat
<chr>       <chr>           <dbl>  <dbl>  <dbl> <dbl>   <dbl> <dbl> <dbl> <lgl> 
1 beta_sigma  beta_sigma[1,1] 0.270 0.0602  0.148 0.273   0.378   213  1.03 NA    
2 beta_sigma  beta_sigma[2,1] 0.298 0.0789  0.151 0.297   0.454   163  1.04 NA    
3 beta_sigma  beta_sigma[3,1] 0.304 0.0668  0.177 0.306   0.426   130  1.06 NA    
4 beta_sigma  beta_sigma[4,1] 0.295 0.0821  0.149 0.293   0.454   208  1.03 NA    
5 beta_sigma  beta_sigma[1,2] 0.291 0.0912  0.135 0.280   0.496    97  1.07 NA    
6 beta_sigma  beta_sigma[2,2] 0.286 0.0853  0.130 0.281   0.460   113  1.06 NA    
7 beta_sigma  beta_sigma[3,2] 0.276 0.0568  0.160 0.278   0.378   240  1.03 NA    
8 beta_sigma  beta_sigma[4,2] 0.283 0.0773  0.136 0.281   0.436   241  1.04 NA  

> beta_tau %>% print(n = 100)
# A tibble: 8 × 10
param_group parameter      mean     sd  `2.5%` `50%` `97.5%` n_eff  Rhat n_Rhat
<chr>       <chr>         <dbl>  <dbl>   <dbl> <dbl>   <dbl> <dbl> <dbl> <lgl> 
1 beta_tau    beta_tau[1,1] 0.131 0.0818 0.00902 0.125   0.293   181  1.04 NA    
2 beta_tau    beta_tau[2,1] 0.341 0.110  0.103   0.349   0.524   128  1.05 NA    
3 beta_tau    beta_tau[3,1] 0.203 0.111  0.0149  0.208   0.405   141  1.05 NA    
4 beta_tau    beta_tau[4,1] 0.345 0.176  0.0402  0.331   0.733   328  1.02 NA    
5 beta_tau    beta_tau[1,2] 0.589 0.130  0.341   0.590   0.842    83  1.07 NA    
6 beta_tau    beta_tau[2,2] 0.482 0.102  0.263   0.485   0.667   200  1.04 NA    
7 beta_tau    beta_tau[3,2] 0.136 0.0873 0.00786 0.128   0.314   206  1.02 NA    
8 beta_tau    beta_tau[4,2] 0.267 0.157  0.0197  0.255   0.585   287  1.02 NA 

> beta_corr %>% print(n = 100)
# A tibble: 1 × 10
param_group parameter  mean    sd `2.5%` `50%` `97.5%` n_eff  Rhat n_Rhat
<chr>       <chr>     <dbl> <dbl>  <dbl> <dbl>   <dbl> <dbl> <dbl> <lgl> 
1 beta_corr   beta_corr 0.520 0.196  0.135 0.518   0.869   127  1.04 NA  

 




#### ---- SIMULATED data (seed = 1; from fixed-C model) + NO covariates + random-C (fitted model):

beta_mu %>% print(n = 100)
# A tibble: 8 × 10
param_group parameter        mean    sd  `2.5%`  `50%` `97.5%` n_eff  Rhat n_Rhat
<chr>       <chr>           <dbl> <dbl>   <dbl>  <dbl>   <dbl> <dbl> <dbl> <lgl> 
1 beta_mu     beta_mu[1,1,1] -0.872 0.277 -1.43   -0.872  -0.377    38  1.17 NA    
2 beta_mu     beta_mu[2,1,1] -0.929 0.157 -1.25   -0.925  -0.641    19  1.35 NA    
3 beta_mu     beta_mu[3,1,1] -0.778 0.176 -1.13   -0.780  -0.447    31  1.20 NA    
4 beta_mu     beta_mu[4,1,1] -1.38  0.198 -1.78   -1.39   -1.01     60  1.09 NA    
5 beta_mu     beta_mu[1,2,1]  0.489 0.294 -0.104   0.490   1.04    131  1.03 NA    
6 beta_mu     beta_mu[2,2,1]  0.147 0.167 -0.169   0.140   0.474    31  1.18 NA    
7 beta_mu     beta_mu[3,2,1]  0.208 0.151 -0.0961  0.213   0.470    58  1.11 NA    
8 beta_mu     beta_mu[4,2,1] -0.183 0.180 -0.531  -0.186   0.138   290  1.02 NA    




#### ---- SIMULATED data (seed = 1; from fixed-C model) + NO covariates + fixed-C (fitted model):
beta_mu %>% print(n = 100)
# A tibble: 8 × 10
param_group parameter         mean    sd `2.5%`   `50%` `97.5%` n_eff  Rhat n_Rhat
<chr>       <chr>            <dbl> <dbl>  <dbl>   <dbl>   <dbl> <dbl> <dbl> <lgl> 
  1 beta_mu     beta_mu[1,1,1] -0.960  0.391 -1.72  -0.965   -0.215  2292 0.999 NA    
2 beta_mu     beta_mu[2,1,1] -1.20   0.228 -1.64  -1.20    -0.768  2367 1.00  NA    
3 beta_mu     beta_mu[3,1,1] -1.05   0.229 -1.49  -1.06    -0.628  2074 1.00  NA    
4 beta_mu     beta_mu[4,1,1] -1.50   0.198 -1.87  -1.50    -1.10   1911 1.00  NA    
5 beta_mu     beta_mu[1,2,1]  0.353  0.390 -0.394  0.344    1.11   1928 1.00  NA    
6 beta_mu     beta_mu[2,2,1]  0.0296 0.234 -0.424  0.0278   0.466  1523 1.00  NA    
7 beta_mu     beta_mu[3,2,1]  0.137  0.221 -0.292  0.138    0.545  2613 1.00  NA    
8 beta_mu     beta_mu[4,2,1] -0.219  0.190 -0.592 -0.225    0.138  1854 1.00  NA    





remove_only_minus_one_rows <- function(x) {
  
  # Helper function to remove rows with only -1 from a single matrix
  remove_minus_one_from_matrix <- function(mat) {
    # Check which rows have all values equal to -1
    rows_to_keep <- apply(mat, 1, function(row) {
      !all(row == -1) && !(row[1] != -1 && all(row[-1] == -1))
    })
    
    # Return matrix with only the rows to keep
    # Preserve matrix structure even if only one row remains
    mat[rows_to_keep, , drop = FALSE]
  }
  
  # Handle different input types
  if (is.matrix(x)) {
    # If input is a single matrix
    return(remove_minus_one_from_matrix(x))
    
  } else if (is.list(x)) {
    # If input is a list (of matrices)
    return(lapply(x, function(element) {
      if (is.matrix(element)) {
        remove_minus_one_from_matrix(element)
      } else {
        # Return element unchanged if it's not a matrix
        element
      }
    }))
    
  } else {
    stop("Input must be a matrix or a list of matrices")
  }
}



str(real_data$x)

x_single_test <- x[[1]]

x_single_test <- remove_only_minus_one_rows(x_single_test)


model_prep_obj <- MetaOrdDTA::MetaOrd_model$new(  
  debugging = TRUE,
  ##
  x = x_single_test,
  indicator_index_test_in_study =  NULL,
  ##
  intercept_only = TRUE,
  cov_data = NULL,
  ##
  n_chains = n_chains,
  ##
  cts = cts,
  ##
  network = FALSE,
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
  init_lists_per_chain = NULL,
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

priors$K_fold_CV_indicator <- rep(1, nrow(x_single_test[[1]]))


stan_data_list$n_obs_cutpoints

priors

# 
# length(init_lists_per_chain[[1]]$C_raw_vec[[1]])
# 
# priors$n_total_C_if_random
# 
# priors
# 
# sum(n_thr * n_studies_per_test)
# 
# stan_data_list$n_studies
# stan_data_list$n_thr
# stan_data_list$n_studies_per_test
# stan_data_list$x

##
## ----  Sample model: ----------------------------------------------------------------
##
{
  max_treedepth <- 9
  adapt_delta <- 0.65
  n_burnin <- 500
  n_iter   <- 300
}
##
model_samples_obj <-  model_prep_obj$sample(   
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
  init_lists_per_chain = init_lists_per_chain
)



##
## ----  Summarise + output results: -------------------------------------------------
##
model_summary_and_trace_obj <- model_samples_obj$summary(
  compute_main_params = TRUE,
  compute_transformed_parameters = TRUE, 
  compute_generated_quantities = TRUE,
  ##
  save_log_lik_trace = TRUE,
  ##
  use_BayesMVP_for_faster_summaries = TRUE)
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

model_summary_and_trace_obj$extract_params(params = c("kappa")) %>% print(n = 20)

# Se_true_PHQ_9 <- Baron_data_list$Se_true_PHQ_9
# Sp_true_PHQ_9 <- Baron_data_list$Sp_true_PHQ_9
# 
# Se_medians <- Se$`50%` * 100









































































