##
## ---- The time-to-event model of Hoyer, Hirt and Kuss (2018) fitted to one data set of the sim study, simulated or real ---------
##
## Hoyer A, Hirt S, Kuss O. Meta-analysis of full ROC curves using bivariate time-to-event models for interval-censored
## data. Res Synth Methods 2018;9(1):62-72. Fitted by dtametaTMB::fitHoyer() (CRAN, version 0.1.4).
##
## A simulated data set is the Stan data list of a MetaOrdDTA fit of the sim:
##   x[group, study, position]             group 1 = non-diseased, group 2 = diseased; position 1 = the number of people of
##                                         the group in the study; position j + 1 = the number of them with a test result at
##                                         or above the study's j-th reported threshold (-1: not reported)
##   cutpoint_index[group, study, position] the threshold (1, ..., n_thr) of position j + 1 (-1: not reported)
##   n_obs_cutpoints[study]                the number of positions of the study after the first
##
## A real data set is one of the files of sim_study/real_data_for_the_DGM_derivation_cleaned (real_x_GAD2.RDS,
## real_x_HADS.RDS, real_x_BAI.RDS), the same files the sim runner reads for the real-data fits of the other models:
##   real_x$nondiseased[study, column], real_x$diseased[study, column]
##                                         column 1 = the number of people of the group in the study; column k + 1 = the
##                                         number of them with a test result at or above threshold k (-1: not reported)
##
## The test results are scores 0, 1, ..., n_thr; threshold k means "positive when the score is k or more". The model needs
## positive test values, a value below the lowest threshold (smallest) and a value above the highest (largest); all three
## are arguments of the functions here, given by the caller.
##

##
## ---- One row per study and reported threshold: TP, FP, FN, TN (a simulated data set) ---------------------------------------------
##
table_of_counts_of_each_study_and_threshold_from_a_Stan_data_list <- function(Stan_data_list) {

      rows <-  vector("list", Stan_data_list$n_studies)
      for (s in seq_len(Stan_data_list$n_studies)) {
            positions <-  1 + seq_len(Stan_data_list$n_obs_cutpoints[s])
            thresholds <-  Stan_data_list$cutpoint_index[1, s, positions]
            reported <-  thresholds != -1
            n_non_diseased <-  Stan_data_list$x[1, s, 1]
            n_diseased <-  Stan_data_list$x[2, s, 1]
            FP <-  Stan_data_list$x[1, s, positions][reported]
            TP <-  Stan_data_list$x[2, s, positions][reported]
            if (any(FP < 0) || any(TP < 0)) {
                  stop(paste("a reported threshold of study", s, "has a count that is not reported"))
            }
            rows[[s]] <-  data.frame(study      = s,
                                     threshold  = thresholds[reported],
                                     TP         = TP,
                                     FN         = n_diseased - TP,
                                     FP         = FP,
                                     TN         = n_non_diseased - FP)
      }
      return(do.call(rbind, rows))

}

##
## ---- One row per study and reported threshold: TP, FP, FN, TN (a real data set) --------------------------------------------------
##
table_of_counts_of_each_study_and_threshold_from_the_real_data <- function(real_x) {

      n_studies <-  nrow(real_x$nondiseased)
      n_thr <-  ncol(real_x$nondiseased) - 1
      rows <-  vector("list", n_studies)
      for (s in seq_len(n_studies)) {
            FP_of_every_threshold <-  real_x$nondiseased[s, 1 + seq_len(n_thr)]
            TP_of_every_threshold <-  real_x$diseased[s, 1 + seq_len(n_thr)]
            if (any((FP_of_every_threshold == -1) != (TP_of_every_threshold == -1))) {
                  stop(paste("study", s, "reports a threshold in one group only"))
            }
            reported <-  FP_of_every_threshold != -1
            rows[[s]] <-  data.frame(study      = s,
                                     threshold  = seq_len(n_thr)[reported],
                                     TP         = TP_of_every_threshold[reported],
                                     FN         = real_x$diseased[s, 1] - TP_of_every_threshold[reported],
                                     FP         = FP_of_every_threshold[reported],
                                     TN         = real_x$nondiseased[s, 1] - FP_of_every_threshold[reported])
      }
      return(do.call(rbind, rows))

}

##
## ---- The fit to a table of counts, and Se and Sp at every threshold with their confidence limits --------------------------------
##
fit_the_Hoyer_2018_model_to_a_table_of_counts <- function(counts,
                                                          n_thr,
                                                          distribution,
                                                          smallest,
                                                          largest,
                                                          confidence_level) {

      tictoc::tic()
      fit <-  dtametaTMB::fitHoyer(data            = counts,
                                   TP              = TP,
                                   FP              = FP,
                                   FN              = FN,
                                   TN              = TN,
                                   study           = study,
                                   threshold       = threshold,
                                   smallest        = smallest,
                                   largest         = largest,
                                   testdirection   = "greater",
                                   eval_threshold  = seq_len(n_thr),
                                   conflevel       = confidence_level,
                                   dist            = distribution,
                                   verbose         = TRUE)
      elapsed <-  tictoc::toc(quiet = TRUE)
      seconds <-  unname(elapsed$toc - elapsed$tic)
      return(list(fit           = fit,
                  counts        = counts,
                  Se_and_Sp     = fit$sensspec,
                  convergence   = fit$fit$convergence,
                  message       = fit$fit$message,
                  seconds       = seconds))

}

##
## ---- The fit to a simulated data set -------------------------------------------------------------------------------------------------
##
fit_the_Hoyer_2018_model_to_a_Stan_data_list <- function(Stan_data_list,
                                                         distribution,
                                                         smallest,
                                                         largest,
                                                         confidence_level) {

      counts <-  table_of_counts_of_each_study_and_threshold_from_a_Stan_data_list(Stan_data_list = Stan_data_list)
      return(fit_the_Hoyer_2018_model_to_a_table_of_counts(counts            = counts,
                                                           n_thr             = Stan_data_list$n_thr,
                                                           distribution      = distribution,
                                                           smallest          = smallest,
                                                           largest           = largest,
                                                           confidence_level  = confidence_level))

}
