


## setwd("/home/enzocerullo/Documents/Work/PhD_work/DTA_MA_NMA_w_missing_thresholds")

source("R_fn_load_data_ordinal_MA_LC_MVP_sim.R")

require(dplyr)
require(cmdstanr)



n_studies <- 25
N_per_study_mean <- 1000
N_per_study_SD <- round(N_per_study_mean/4)
assume_perfect_GS <- 1
seed <- 123

set.seed(seed)

missing_indicator <- -1

options(max.print = 100000)


n_studies = n_studies
N_per_study_mean = N_per_study_mean
N_per_study_SD = N_per_study_SD
assume_perfect_GS = assume_perfect_GS
seed = seed
true_Mean_prev = 0.20

# Run simulated data - this simulates data from FIVE (5) diagnostic tests (1 BINARY reference test + 4 ORDINAL index tests)
sim_results <- R_fn_sim_data_ordinal_MA(   n_studies = n_studies,
                                                          N_per_study_mean = N_per_study_mean,
                                                          N_per_study_SD = N_per_study_SD,
                                                          assume_perfect_GS = assume_perfect_GS,
                                                          seed = seed)

y_list <- sim_results$y_list
str(y_list)


sim_results$n_total_obs_thr_per_test
sim_results$N_per_study_vec

index_test_chosen_index <- 3
n_thr <- 10

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

y_list_example_1



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

true_Se_OVERALL

## sim_results$Se_per_study_all_tests_all_thresholds_list

## Apply missing thresholds:
agg_data_cumulative_with_missing_thresholds <- agg_data_cumulative

true_Fp_OVERALL <- 100 - true_Sp_OVERALL

plot(x = true_Fp_OVERALL, y = true_Se_OVERALL)

# true_Se_OVERALL
# true_Sp_OVERALL
# 
# plot(x = true_Sp_OVERALL, y = true_Se_OVERALL)
 
# plot(x = 1:27,  y = true_Se_OVERALL, col = "red", main = "Se (red line) and Fp (green line)", ylim = c(0, 100))
# lines(x = 1:27, y = 100 - true_Sp_OVERALL, col = "green")

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
n_d <- agg_data_cumulative_with_missing_thresholds$n_d
n_total_d <- agg_data_cumulative_with_missing_thresholds$n_total_d
print(paste("x_d = ")) ; print(x_d)
## Non-diseased group:
x_nd <- agg_data_cumulative_with_missing_thresholds$x_nd
n_nd <- agg_data_cumulative_with_missing_thresholds$n_nd
n_total_nd <- agg_data_cumulative_with_missing_thresholds$n_total_nd
print(paste("x_nd = ")) ; print(x_nd)

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

