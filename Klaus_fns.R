



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


