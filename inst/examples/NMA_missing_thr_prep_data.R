


## setwd("/home/enzocerullo/Documents/Work/PhD_work/DTA_MA_NMA_w_missing_thresholds")

 
require(dplyr)
require(cmdstanr)


n_tests <- 5
n_index_tests <- n_tests - 1
##
n_studies <- 25
N_per_study_mean <- 2500
N_per_study_SD <- round(N_per_study_mean/4)
assume_perfect_GS <- 1

seed <- 123
set.seed(seed)

missing_indicator <- -1

options(max.print = 100000)


# Run simulated data - this simulates data from FIVE (5) diagnostic tests (1 BINARY reference test + 4 ORDINAL index tests)
sim_results <- MetaOrdDTA::R_fn_sim_data_ordinal_MA(     n_studies = n_studies,
                                                          N_per_study_mean = N_per_study_mean,
                                                          N_per_study_SD = N_per_study_SD,
                                                          assume_perfect_GS = assume_perfect_GS,
                                                          seed = 123,
                                                          threshsolds_latent_between_study_SD = 0.10)
 index_test_chosen_index <- 1 + 2

{
  
    y_list <- sim_results$y_list
    str(y_list)
    
    sim_results$n_thr_per_test
    sim_results$N_per_study_vec

    # n_thr <- 10
    n_thr_vec <- sim_results$n_thr_per_test_excl_ref ; n_thr_vec
    # ## subset to remove ref test:
    # n_thr <- n_thr[-c(1)]
    n_cat <- n_thr_vec + 1 ; n_cat
    n_thr_max <- max(n_thr_vec) ; n_thr_max
    n_cat_max <- max(n_cat) ; n_cat_max
    
}


{
  true_DGM_Fp <- 100 * sim_results$true_DGM_Fp[[index_test_chosen_index]]
  true_DGM_Sp <- 100 * sim_results$true_DGM_Sp[[index_test_chosen_index]]
  true_DGM_Se <- 100 * sim_results$true_DGM_Se[[index_test_chosen_index]]
}

sim_results$Se_OVERALL_all_tests_all_thresholds

sim_results$Se_OVERALL_all_tests_all_thresholds
# sim_results$Se_OVERALL_all_tests_all_thresholds[1, 1:n_thr_vec[1]]
sim_results$Se_OVERALL_all_tests_all_thresholds[2, 1:n_thr_vec[2]]
 
 true_Se_OVERALL <- sim_results$Se_OVERALL_all_tests_all_thresholds[index_test_chosen_index, 1:n_thr_vec[index_test_chosen_index - 1]] ; true_Se_OVERALL
 true_Sp_OVERALL <- sim_results$Sp_OVERALL_all_tests_all_thresholds[index_test_chosen_index, 1:n_thr_vec[index_test_chosen_index - 1]] ; true_Sp_OVERALL

  # plot(y = true_Se_OVERALL, x = 1 - true_Sp_OVERALL)

# sim_results$Se_per_study_all_tests_all_thresholds_list

## Now for the first example we will only take one index test, as initially we are evaluating just a "simple" model
## where there's only a single index test with 12 thresholds (so 13 categories:
##

##
y_list_example_1 <- list()
for (s in 1:n_studies) { 
    N <- nrow(y_list[[s]])
    y <- array(NA, dim = c(N, n_index_tests + 1))
    for (t in 1:(n_index_tests + 1)) { 
      y[, t] <- y_list[[s]][, t]
    }
    y_list_example_1[[s]] <- y
}

min(y[,index_test_chosen_index])
max(y[,index_test_chosen_index])

# y_list_example_1


## Convert to AGGREGATE data as these more basic NMA/MA models don't use individual-level data:
agg_data_cumulative <- convert_to_aggregate_counts( y_list = y_list_example_1, 
                                                    n_studies = n_studies,
                                                    n_tests_inc_ref = n_index_tests + 1,
                                                    n_thr = n_thr_vec)

agg_data_cumulative$Se_per_study_list


## Get OBSERVED "true" values BEFORE you remove any thresholds and BEFORE you remove any tests for a given study s!!!
Se_per_study_list <- agg_data_cumulative$Se_per_study_list
Sp_per_study_list <- agg_data_cumulative$Sp_per_study_list

true_Se_OVERALL <- true_Sp_OVERALL <- list() 
true_Se_OVERALL_weighted <- true_Sp_OVERALL_weighted <- list() 

agg_data_cumulative$x_d_list

# 0.9897257 0.9566693 0.8985080 0.7359064 0.3932815 
# 0.2462610 0.5247035 0.7138989 0.9045126 0.9878030  # Sp

for (t in 1:n_index_tests) {
  
      ##
      ## this is the means of the OBSERVED (i.e. specific to the seed set for the given simulation) study-specific Se's at each threshold:
      ##
      true_Se_OVERALL[[t]] <- 100 * colMeans(Se_per_study_list[[t]]) ; true_Se_OVERALL[[t]]
      true_Sp_OVERALL[[t]] <- 100 * colMeans(Sp_per_study_list[[t]]) ; true_Sp_OVERALL[[t]]
      ##
      ## this is equal to: Se @ threshold k  = {number of TOTAL true-positives across all studies} / {TOTAL N across all studies}:
      ##
      true_Se_OVERALL_weighted[[t]] <- sim_results$Se_OVERALL_all_tests_all_thresholds[t + 1, 1:n_thr_vec[t]]*100 ; true_Se_OVERALL_weighted[[t]]
      true_Sp_OVERALL_weighted[[t]] <- sim_results$Sp_OVERALL_all_tests_all_thresholds[t + 1, 1:n_thr_vec[t]]*100 ; true_Sp_OVERALL_weighted[[t]]
      ##
      ##
      print(true_Se_OVERALL[[t]] - true_Se_OVERALL_weighted[[t]])
      print(true_Sp_OVERALL[[t]] - true_Sp_OVERALL_weighted[[t]])

}


true_Se_OVERALL[[3]] - true_Se_OVERALL_weighted[[3]]
true_Sp_OVERALL[[3]] - true_Sp_OVERALL_weighted[[3]]


true_Se_OVERALL_weighted
true_Sp_OVERALL


# #### ------------------- Apply ** MISSING TESTS ** (hence "NMA" - optional): --------------------------------------------------------------------

outs_NMA <- create_test_in_study_for_NMA( seed = seed,
                                          n_studies = n_studies,
                                          n_tests_inc_ref = n_index_tests + 1,
                                          prob_present = 0.40,
                                          min_index_tests_per_study = 1)

indicator_index_test_in_study <- outs_NMA$indicator_index_test_in_study ; indicator_index_test_in_study
n_index_tests_per_study      <- outs_NMA$n_index_tests_per_study ; n_index_tests_per_study

# 
# ## Then update the "true values" list:
# for (t in 1:n_index_tests) {
#   indicator_mat <- array(1.0, dim = c(n_studies, n_thr_vec[t]))
#   for (k in 1:n_thr_vec[t]) {
#      indicator_mat[,k] <- ifelse(indicator_index_test_in_study[,t] == 0, NA, 1)
#   }
#   Se_per_study_list[[t]] <- Se_per_study_list[[t]] * indicator_mat
# }
# # 
# 
# ## Then re-calculate:
# for (t in 1:n_index_tests) {
#   ##
#   ## this is the means of the OBSERVED (i.e. specific to the seed set for the given simulation) study-specific Se's at each threshold:
#   ##
#   true_Se_OVERALL[[t]] <- 100 * colMeans(Se_per_study_list[[t]], na.rm = TRUE) ; true_Se_OVERALL[[t]]
#   true_Sp_OVERALL[[t]] <- 100 * colMeans(Sp_per_study_list[[t]], na.rm = TRUE) ; true_Sp_OVERALL[[t]]
#   ##
#   ## this is equal to: Se @ threshold k  = {number of TOTAL true-positives across all studies} / {TOTAL N across all studies}:
#   ##
#   true_Se_OVERALL_weighted[[t]] <- sim_results$Se_OVERALL_all_tests_all_thresholds[t + 1, 1:n_thr_vec[t]]*100 ; true_Se_OVERALL_weighted[[t]]
#   true_Sp_OVERALL_weighted[[t]] <- sim_results$Sp_OVERALL_all_tests_all_thresholds[t + 1, 1:n_thr_vec[t]]*100 ; true_Sp_OVERALL_weighted[[t]]
#   ##
#   ##
#   print(true_Se_OVERALL[[t]] - true_Se_OVERALL_weighted[[t]])
#   print(true_Sp_OVERALL[[t]] - true_Sp_OVERALL_weighted[[t]])
#   
# }



agg_data_cumulative$x_nd_list
####
#### Non-diseased group:
####
x_nd <- agg_data_cumulative$x_nd_list
# n_nd <- agg_data_cumulative$n_nd_list
n_total_nd   <- agg_data_cumulative$n_total_nd
print(paste("x_non_diseased_cumulative = ")) ; print(x_nd)
####
#### Diseased group:
####
x_d <- agg_data_cumulative$x_d_list
# n_d <- agg_data_cumulative$n_d_list
n_total_d      <- agg_data_cumulative$n_total_d
print(paste("x_diseased_cumulative = ")) ; print(x_d)

###
### Prep / save final NMA simulated data outputs:
x_NMA <- list(x_nd, x_d)
x <- x_NMA

data_example_1_NMA_list <- list()
data_example_1_NMA_list$x_NMA <- x_NMA
data_example_1_NMA_list$indicator_index_test_in_study <- indicator_index_test_in_study
data_example_1_NMA_list$n_index_tests_per_study <- n_index_tests_per_study

saveRDS(object = data_example_1_NMA_list, file = file.path(getwd(), "inst", "examples", "data_example_1_NMA_list.RDS"))





#### ------------------- Apply ** MISSING THRESHOLDS ** : -----------------------------------------------------------------------------

agg_data_cumulative_with_missing_thresholds <- agg_data_cumulative

agg_data_cumulative_with_missing_thresholds$x_d_list

x_PHQ9 <- list()
PHQ9_index <- 4

for (c in 1:2) {
  x_PHQ9[[c]] <- x[[c]][[PHQ9_index]]
}
x_PHQ9



x_PHQ9[[1]]


## Apply PHQ-9 "missing thresholds":
## (1) 1st quintile (e.g., studies 1-5   if 25 total studies)  0% of studies report @ ALL threshodls, so do nothing for these.
## (2) 2nd quintile (e.g., studies 6-10  if 25 total studies) report 
## (3) 3rd quintile (e.g., studies 11-15 if 25 total studies) report ...
## (3) 4th quintile (e.g., studies 16-20 if 25 total studies) report ...
## (3) 5th quintile (e.g., studies 21-25 if 25 total studies) report ...

# n_groups <- 5
# n_studies_div_frac <- ceiling(n_studies/n_groups)
# n_groups * n_studies_div_frac
# n_studies
# last_group <- n_studies - frac*n_studies_div_frac
# last_group


study_thr_groups <- divide_studies(n_studies = n_studies, n_groups = 5)
##
set.seed(seed)
##
{
    ## When a study ONLY reports at the single cut-off of 10
    ## (i.e. the "standard" PHQ-9 screening cut-off):
    case_1 <- list()
    ##
    case_1$studies <- study_thr_groups[[1]]
    ##
    case_1$thr_combo_vec_1 <- c(10)
    ##
    print(case_1)
}
{
    ## This is based off of the thresholds reported in Levis et al, 2019 
    ## (BMJ PHQ-9 paper):
    case_2 <- list()
    ##
    case_2$studies <- study_thr_groups[[2]]
    ##
    case_2$thr_combo_vec_1 <- c(5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
    ##
    print(case_2)
}
{
    ## When a study reports at a SINGLE cut-off, cut NOT the standard cut-off of 10
    ## (can be anything between 8 and 15):
    case_3 <- list()
    ##
    case_3$studies <- study_thr_groups[[3]] 
    ##
    random_single_thr <- sample(8:15, 1)
    case_3$thr_combo_vec_1 <- c(random_single_thr)
    ##
    print(case_3)
}
{
    case_4 <- list()
    ##
    case_4$studies <- study_thr_groups[[4]]
    ##
    case_4$thr_combo_vec_1  <- c(8, 9, 10, 11, 12)
    case_4$thr_combo_vec_2  <- c(10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21)
    case_4$thr_combo_vec_3  <- c(9, 10, 11, 12, 13)
    case_4$thr_combo_vec_4  <- c(10, 11, 12, 13, 14)
    case_4$thr_combo_vec_5  <- c(9, 10, 11, 12, 13, 14, 15)
    case_4$thr_combo_vec_6  <- c(6, 7, 8)
    case_4$thr_combo_vec_7  <- c(6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
    case_4$thr_combo_vec_8  <- c(11, 12, 13)
    case_4$thr_combo_vec_9  <- c(5, 6, 7, 8, 9, 10, 11, 12)
    case_4$thr_combo_vec_10 <- c(4, 5, 6, 7, 8, 9, 10)
    ## e.g. a study did 10 for screening and 15 for DIAGNOSIS (not common w/ PHQ-9 
    ## but some studies do this or investigate it):
    case_4$thr_combo_vec_11 <- c(10, 15) 
    ## this is also based on the study that looked at PHQ-9 for both screening (i.e., using the
    ## standard/lower thresholds) AND DIAGNOSIS (higher thresholds):
    case_4$thr_combo_vec_12 <- c(9, 10, 11, 14, 15, 16) 
    ##
    print(case_4)
}
{
    case_5 <- list()
    ##
    case_5$studies <- study_thr_groups[[5]]
    ##
    case_5$thr_combo_vec_1 <- c(1:27) ## i.e. these studies report at ALL thresholds
    ##
    print(case_5)
}

# Combine all cases into a single list
all_cases <- list(
  case_1 = case_1,
  case_2 = case_2,
  case_3 = case_3,
  case_4 = case_4,
  case_5 = case_5
)


x_PHQ9

x_PHQ9_w_missing_thr <- MetaOrdDTA:::R_fn_apply_missingness_pattern_PHQ_9( x_PHQ = x_PHQ9,
                                                              case_list = all_cases)

x_PHQ9_w_missing_thr
 
# x_PHQ9_w_missing_thr_old <- x_PHQ9_w_missing_thr
# x_PHQ9_w_missing_thr[[1]] - x_PHQ9_w_missing_thr_old[[1]]

## Compute the % of missing thr data:
100 * sum(x_PHQ9_w_missing_thr[[1]] == -1) / (n_studies * ncol(x_PHQ9_w_missing_thr[[1]]))

save_list_obj <- list()
x_NMA <- x
save_list_obj$indicator_index_test_in_study <- indicator_index_test_in_study
save_list_obj$n_index_tests_per_study <- n_index_tests_per_study
save_list_obj$x_NMA <- x_NMA
save_list_obj$x_PHQ9 <- x_PHQ9
save_list_obj$x_PHQ9_w_missing_thr <- x_PHQ9_w_missing_thr
saveRDS(object = save_list_obj,
        file.path(getwd(), "inst", "examples", "data_example_1_NMA_list.RDS"))

# {
# 
#     ## First 2 studies only report data for the "middle 8" thresholds (starting from threshold 3, so thr = {3, 4, 5, 6, 7, 8, 9, 10})
#     studies_subset_vec <- c(1, 2)
#     missing_thr_subset_vec <- c(1:2, 11:12)
#     agg_data_cumulative_with_missing_thresholds <-  apply_thr_missingness(  agg_data_cumulative = agg_data_cumulative_with_missing_thresholds,
#                                                                             studies_subset_vec = studies_subset_vec,
#                                                                             missing_thr_subset_vec = missing_thr_subset_vec,
#                                                                             missing_indicator = missing_indicator)
#     ## The next 2 studies only report at the "middle 4" thresholds (starting at threshold 5): so thr = {5, 6, 7, 8}:
#     studies_subset_vec <- c(3, 4)
#     missing_thr_subset_vec <- c(1:4, 9:12)
#     agg_data_cumulative_with_missing_thresholds <-  apply_thr_missingness(  agg_data_cumulative = agg_data_cumulative_with_missing_thresholds,
#                                                                             studies_subset_vec = studies_subset_vec,
#                                                                             missing_thr_subset_vec = missing_thr_subset_vec,
#                                                                             missing_indicator = missing_indicator)
#     ## The next two studies only report data at only a single threshold - threshold #7:
#     studies_subset_vec <- c(5, 6)
#     missing_thr_subset_vec <- c(1:6, 8:12)
#     agg_data_cumulative_with_missing_thresholds <-  apply_thr_missingness(  agg_data_cumulative = agg_data_cumulative_with_missing_thresholds,
#                                                                             studies_subset_vec = studies_subset_vec,
#                                                                             missing_thr_subset_vec = missing_thr_subset_vec,
#                                                                             missing_indicator = missing_indicator)
#     ## The next two studies only report data at 4 non-adjacent thresholds - so thr = {3, 5, 7, 9}:
#     studies_subset_vec <- c(7, 8)
#     missing_thr_subset_vec <- c(1, 2, 4, 6, 8, 10, 11, 12)
#     agg_data_cumulative_with_missing_thresholds <-  apply_thr_missingness(  agg_data_cumulative = agg_data_cumulative_with_missing_thresholds,
#                                                                             studies_subset_vec = studies_subset_vec,
#                                                                             missing_thr_subset_vec = missing_thr_subset_vec,
#                                                                             missing_indicator = missing_indicator)
#     ## And finally, the last 2 studies report data at ALL thresholds:
#     #  -----  Do nothing here - as no missing thresholds for study #10!
# 
# }
# 
# ## Now let's look at the overall % of missing thresholds:
# total_missing_thr <- sum(agg_data_cumulative_with_missing_thresholds$x_diseased == 0.999) ; total_missing_thr
# prop_missing_thr <- total_missing_thr / (n_studies * n_thr) ; prop_missing_thr
# ## Just under half the threshold data is missing (46.67%).
# 
# #  
#  









