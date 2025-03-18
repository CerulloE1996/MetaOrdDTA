


## setwd("/home/enzocerullo/Documents/Work/PhD_work/DTA_MA_NMA_w_missing_thresholds")
# 
MetaOrdinal_admin_path_root <- setwd("/home/enzo/Documents/Work/PhD_work/R_packages/MetaOrdDTA")
##
MetaOrdinal_admin_path_R           <- file.path(MetaOrdinal_admin_path_root, "R")
MetaOrdinal_admin_path_inst        <- file.path(MetaOrdinal_admin_path_root, "inst")
MetaOrdinal_admin_path_examples    <- file.path(MetaOrdinal_admin_path_inst, "examples")
MetaOrdinal_admin_path_stan_models <- file.path(MetaOrdinal_admin_path_inst, "stan_models")


source(file.path(MetaOrdinal_admin_path_root, "R_fn_sim_data_ordinal_MA.R"))


####  -------------------------------------------
require(TruncatedNormal)
require(dplyr)
require(cmdstanr)
require(ggplot2)

####  -------------------------------------------
options(max.print = 100000)


####  -------------------------------------------
{
  outs <- set_pkg_eg_path_and_wd()
  outs$user_root_dir
  outs$user_MetaOrdDTA_dir
  outs$pkg_example_path
}


####  -------------------------------------------
seed <- 123
set.seed(seed, kind = "L'Ecuyer-CMRG")


####  ------------------------------------------- 
{
  n_studies <- 25
  N_per_study_mean <- 2500
  N_per_study_SD <- 250
  assume_perfect_GS <- 1
  missing_indicator <- -1
}
 

####  -------------------------------------------
## Run simulated data - this simulates data from FIVE (5) diagnostic tests (1 BINARY reference test + 4 ORDINAL index tests)
sim_results <- MetaOrdDTA:::R_fn_sim_data_ordinal_MA(  seed = seed,
                                          n_studies = n_studies,
                                          N_per_study_mean = N_per_study_mean,
                                          N_per_study_SD = N_per_study_SD,
                                          assume_perfect_GS = assume_perfect_GS)

y_list <- sim_results$y_list
str(y_list)


sim_results$n_total_obs_thr_per_test
sim_results$N_per_study_vec

index_test_chosen_index <- 5
n_thr <- 27 ## PHQ-9 has 9 Qs w/ score 0-3 each qs. 

sim_results$Se_OVERALL_all_tests_all_thresholds


sim_results$Se_OVERALL_all_tests_all_thresholds[index_test_chosen_index, 1:n_thr]

 

# true_Se_OVERALL <- sim_results$Se_OVERALL_all_tests_all_thresholds[index_test_chosen_index, 1:n_thr] ; true_Se_OVERALL
# true_Sp_OVERALL <- sim_results$Sp_OVERALL_all_tests_all_thresholds[index_test_chosen_index, 1:n_thr] ; true_Sp_OVERALL

### plot(y = true_Se_OVERALL, x = 1 - true_Sp_OVERALL)

# sim_results$Se_per_study_all_tests_all_thresholds_list

## Now for the first example we will only take one index test, as initially we are evaluating just a "simple" model
## where there's only a single index test with 12 thresholds (so 13 categories:
##

##
n_cat <- n_thr + 1
y_list_example_1 <- list()
for (s in 1:n_studies) { 
    N <- nrow(y_list[[s]])
    y <- array(NA, dim = c(N, 2))
    y[, 1] <- y_list[[s]][, 1]
    y[, 2] <- y_list[[s]][, index_test_chosen_index] ## take test #5 as the index test for this example 
    y_list_example_1[[s]] <- y
}


min(unlist(y_list_example_1))
max(unlist(y_list_example_1))
# y_list_example_1

agg_data_cumulative <-  MetaOrdDTA:::convert_to_aggregate_counts_single_test( y_list_example_1, 
                                                      n_studies,
                                                      n_thr)

# Convert to AGGREGATE data as these more basic NMA/MA models don't use individual-level data:
agg_data_cumulative <- convert_to_aggregate_counts_single_test( y_list_example_1, 
                                                                n_studies,
                                                                n_thr)

agg_data_cumulative

str(y_list_example_1)

agg_data_cumulative$Se_per_study
agg_data_cumulative$Sp_per_study
# this is the means of the OBSERVED (i.e. specific to the seed set for the given simulation) study-specific Se's at each threshold:
true_Se_OVERALL <- 100 * colMeans(agg_data_cumulative$Se_per_study) ; true_Se_OVERALL 
true_Sp_OVERALL <- 100 * colMeans(agg_data_cumulative$Sp_per_study) ; true_Sp_OVERALL
# this is equal to: Se @ threshold k  = {number of TOTAL true-positives across all studies} / {TOTAL N across all studies}:
true_Se_OVERALL_weighted <- sim_results$Se_OVERALL_all_tests_all_thresholds[index_test_chosen_index, 1:n_thr]*100 ; true_Se_OVERALL_weighted 
true_Sp_OVERALL_weighted <- sim_results$Sp_OVERALL_all_tests_all_thresholds[index_test_chosen_index, 1:n_thr]*100 ; true_Sp_OVERALL_weighted 
##

true_Se_OVERALL - true_Se_OVERALL_weighted
true_Sp_OVERALL - true_Sp_OVERALL_weighted

## sim_results$Se_per_study_all_tests_all_thresholds_list

## Apply missing thresholds:
agg_data_cumulative_with_missing_thresholds <- agg_data_cumulative



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
 
# ## Diseased group:
x_d <- agg_data_cumulative_with_missing_thresholds$x_d
print(paste("x_d = ")) ; print(x_d)
# n_d <- agg_data_cumulative_with_missing_thresholds$n_d
n_total_d <- agg_data_cumulative_with_missing_thresholds$n_total_d

## Non-diseased group:
x_nd <- agg_data_cumulative_with_missing_thresholds$x_nd
print(paste("x_nd = ")) ; print(x_nd)
# n_nd <- agg_data_cumulative_with_missing_thresholds$n_nd
n_total_nd <- agg_data_cumulative_with_missing_thresholds$n_total_nd
print(paste("n_total_nd = ")) ; print(n_total_nd)



# x_non_diseased <- x_nd
# x_diseased <- x_d

# 
# 
# Stan_data_list <- R_fn_prep_MA_data(  x_diseased = x_d, 
#                                        x_non_diseased = x_nd, 
#                                        n_thr = n_thr)
# 
# 
# 
# 
# 

   
   

## agg_data_cumulative_with_missing_thresholds$x_d
# x_nd_cat <- convert_cumulative_to_category( cumulative_matrix = n_nd)
# x_non_diseased_cumulative
# x_nd_cat
# 
# 
# x_d_cat <- convert_cumulative_to_category( cumulative_matrix = n_diseased_cumulative)
# x_d_cat
# 
# x_cat <- list(x_nd_cat, x_d_cat)


# y_nd <- categorical_to_individual(x_cat, missing_indicator = -1, binary_disease_indicator = 1)



# Stan_data_list

# 
# reorder_increasing <- function(matrix_data) {
#   
#       # Apply to each row
#       t(apply(matrix_data, 1, function(row) {
#         # Sort in ascending order
#         sort(row)
#       }))
#   
# }
# 
# 
# reorder_decreasing <- function(matrix_data) {
#   
#   # Apply to each row
#   t(apply(matrix_data, 1, function(row) {
#     # Sort in ascending order
#     sort(row, decreasing = TRUE)
#   }))
#   
# }
# 
# 
# {
#       Stan_data_list$x_diseased     <- reorder_increasing(Stan_data_list$x_diseased)
#       Stan_data_list$x_non_diseased <- reorder_increasing(Stan_data_list$x_non_diseased)
#       
#       for (c in 1:2) { 
#         Stan_data_list$n[[c]] <- reorder_increasing(Stan_data_list$n[[c]])
#         Stan_data_list$x[[c]] <- reorder_increasing(Stan_data_list$x[[c]])
#         Stan_data_list$x_with_missings[[c]] <- reorder_increasing(Stan_data_list$x_with_missings[[c]])
#       }
# }
# 
# 
# 
##

# 
# 
# 
# print(paste("true_Se_OVERALL = ")) ; print(true_Se_OVERALL)
# print(paste("true_Sp_OVERALL = ")) ; print(true_Sp_OVERALL)
# 


