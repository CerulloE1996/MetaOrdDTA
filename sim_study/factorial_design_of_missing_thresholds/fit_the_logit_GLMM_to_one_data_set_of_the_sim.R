##
## ---- The logit GLMM for multiple thresholds fitted to one data set of the sim study, simulated or real --------------------------
##
## The model called "logit GLMM" in Stoye, Raths, Hapfelmeier, Fomenko, Kuss and Hoyer (2026), Biometrical Journal 68(4):e70147,
## their equation (6), which they attribute to Hoyer and Kuss (2018), Stat Methods Med Res 27(5):1410-1421, and fitted with
## lme4. For study i and threshold t:
##       logit(Sp_it) = b01 + b02 * t + u1_i        (TN_it out of the non-diseased of study i)
##       logit(Se_it) = b11 + b12 * t + u2_i        (TP_it out of the diseased of study i)
##       (u1_i, u2_i) bivariate normal with SDs s1, s2 and correlation rho
## The counts of one study at its several thresholds enter as separate binomial rows, and the threshold is not transformed.
## Se and Sp at threshold k are those of a study whose random effects are 0: plogis(b11 + b12 * k) and plogis(b01 + b02 * k);
## their confidence limits are Wald limits on the logit scale from the covariance matrix of the fixed effects.
## The logit GLMM is a competitor of the sim only: it is fitted by lme4::glmer and no prior of any model is taken from it.
##
## The table of counts is the one of fit_the_Hoyer_2018_time_to_event_model_to_one_data_set_of_the_sim.R (one row per study
## and reported threshold: study, threshold, TP, FN, FP, TN); that file is sourced before this one.
##

##
## ---- The fit to a table of counts, and Se and Sp at every threshold with their confidence limits --------------------------------
##
fit_the_logit_GLMM_to_a_table_of_counts <- function(counts,
                                                    n_thr,
                                                    confidence_level) {

      ## (two rows per row of counts: the non-diseased row counts the test negatives, the diseased row the test positives)
      long <-  rbind(data.frame(study                       = counts$study,
                                threshold                   = counts$threshold,
                                is_non_diseased             = 1,
                                is_diseased                 = 0,
                                n_with_the_outcome          = counts$TN,
                                n_without_the_outcome       = counts$FP),
                     data.frame(study                       = counts$study,
                                threshold                   = counts$threshold,
                                is_non_diseased             = 0,
                                is_diseased                 = 1,
                                n_with_the_outcome          = counts$TP,
                                n_without_the_outcome       = counts$FN))
      long$study <-  factor(long$study)
      ## (the threshold of a row, in the column of its own group and 0 in the other: one slope per group)
      long$threshold_of_the_non_diseased <-  long$is_non_diseased * long$threshold
      long$threshold_of_the_diseased <-  long$is_diseased * long$threshold
      messages_of_the_fit <-  character(0)
      tictoc::tic()
      fit <-  withCallingHandlers(
                    lme4::glmer(cbind(n_with_the_outcome, n_without_the_outcome) ~
                                      0 + is_non_diseased + threshold_of_the_non_diseased + is_diseased + threshold_of_the_diseased +
                                      (0 + is_non_diseased + is_diseased | study),
                                data     = long,
                                family   = stats::binomial(link = "logit"),
                                verbose  = 1),
                    warning = function(w) {
                          messages_of_the_fit <<-  c(messages_of_the_fit, conditionMessage(w))
                          invokeRestart("muffleWarning")
                    })
      elapsed <-  tictoc::toc(quiet = TRUE)
      seconds <-  unname(elapsed$toc - elapsed$tic)
      ##
      b <-  lme4::fixef(fit)
      V <-  as.matrix(stats::vcov(fit))
      z <-  stats::qnorm(1 - (1 - confidence_level) / 2)
      k <-  seq_len(n_thr)
      logit_Sp <-  b[["is_non_diseased"]] + b[["threshold_of_the_non_diseased"]] * k
      logit_Se <-  b[["is_diseased"]] + b[["threshold_of_the_diseased"]] * k
      SE_logit_Sp <-  sqrt(V["is_non_diseased", "is_non_diseased"] + k^2 * V["threshold_of_the_non_diseased", "threshold_of_the_non_diseased"] +
                           2 * k * V["is_non_diseased", "threshold_of_the_non_diseased"])
      SE_logit_Se <-  sqrt(V["is_diseased", "is_diseased"] + k^2 * V["threshold_of_the_diseased", "threshold_of_the_diseased"] +
                           2 * k * V["is_diseased", "threshold_of_the_diseased"])
      Se_and_Sp <-  data.frame(threshold     = k,
                               conflevel     = confidence_level,
                               Sens          = stats::plogis(logit_Se),
                               SensCI_Lower  = stats::plogis(logit_Se - z * SE_logit_Se),
                               SensCI_Upper  = stats::plogis(logit_Se + z * SE_logit_Se),
                               Spec          = stats::plogis(logit_Sp),
                               SpecCI_Lower  = stats::plogis(logit_Sp - z * SE_logit_Sp),
                               SpecCI_Upper  = stats::plogis(logit_Sp + z * SE_logit_Sp))
      SDs_and_correlation <-  as.data.frame(lme4::VarCorr(fit))
      parameters <-  c(b01  = b[["is_non_diseased"]],
                       b02  = b[["threshold_of_the_non_diseased"]],
                       b11  = b[["is_diseased"]],
                       b12  = b[["threshold_of_the_diseased"]],
                       s1   = SDs_and_correlation$sdcor[which(SDs_and_correlation$var1 == "is_non_diseased" & is.na(SDs_and_correlation$var2))],
                       s2   = SDs_and_correlation$sdcor[which(SDs_and_correlation$var1 == "is_diseased" & is.na(SDs_and_correlation$var2))],
                       rho  = SDs_and_correlation$sdcor[which(!is.na(SDs_and_correlation$var2))])
      return(list(fit                    = fit,
                  counts                 = counts,
                  Se_and_Sp              = Se_and_Sp,
                  parameters             = parameters,
                  code_of_the_fit        = fit@optinfo$conv$opt,
                  messages_of_the_fit    = c(messages_of_the_fit, unlist(fit@optinfo$conv$lme4$messages)),
                  singular               = lme4::isSingular(fit),
                  seconds                = seconds))

}
