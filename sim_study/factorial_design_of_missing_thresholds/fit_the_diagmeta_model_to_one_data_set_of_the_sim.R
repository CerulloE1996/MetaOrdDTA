##
## ---- The multiple-cutoffs linear mixed model of Steinhauser, Schumacher and Rücker (2016) fitted to one data set ---------------
##
## Steinhauser S, Schumacher M, Rücker G. Modelling multiple thresholds in meta-analysis of diagnostic test accuracy studies.
## BMC Med Res Methodol 2016;16:97. Fitted by diagmeta::diagmeta() and evaluated at thresholds 1..n_thr via diagmeta::diagstats().
##
## Uses the same counts table as fit_the_Hoyer_2018_time_to_event_model_to_one_data_set_of_the_sim.R and
## fit_the_logit_GLMM_to_one_data_set_of_the_sim.R (columns: study, threshold, TP, FN, FP, TN).
##

fit_the_diagmeta_model_to_a_table_of_counts <- function(counts,
                                                        n_thr,
                                                        model_of_diagmeta,
                                                        distr_of_diagmeta,
                                                        confidence_level) {

      messages_of_the_fit <-  character(0)
      tictoc::tic()
      fit <-  withCallingHandlers(
                    diagmeta::diagmeta(TP         = counts$TP,
                                       FP         = counts$FP,
                                       TN         = counts$TN,
                                       FN         = counts$FN,
                                       cutoff     = counts$threshold,
                                       studlab    = counts$study,
                                       distr      = distr_of_diagmeta,
                                       model      = model_of_diagmeta,
                                       equalvar   = FALSE,
                                       log.cutoff = FALSE,
                                       level      = confidence_level,
                                       silent     = TRUE),
                    warning = function(w) {
                          messages_of_the_fit <<-  c(messages_of_the_fit, conditionMessage(w))
                          invokeRestart("muffleWarning")
                    })
      elapsed <-  tictoc::toc(quiet = TRUE)
      seconds <-  unname(elapsed$toc - elapsed$tic)
      ##
      k <-  seq_len(n_thr)
      st <-  diagmeta::diagstats(fit, cutoff = k, level = confidence_level)
      Se_and_Sp <-  data.frame(threshold     = k,
                               conflevel     = confidence_level,
                               Sens          = st$Sens,
                               SensCI_Lower  = st$lower.Sens,
                               SensCI_Upper  = st$upper.Sens,
                               Spec          = st$Spec,
                               SpecCI_Lower  = st$lower.Spec,
                               SpecCI_Upper  = st$upper.Spec)
      lmer_fit <-  fit$result.lmer
      parameters <-  c(alpha0 = fit$regr$alpha0,
                       beta0  = fit$regr$beta0,
                       alpha1 = fit$regr$alpha1,
                       beta1  = fit$regr$beta1)
      return(list(fit                    = fit,
                  counts                 = counts,
                  Se_and_Sp              = Se_and_Sp,
                  parameters             = parameters,
                  code_of_the_fit        = lmer_fit@optinfo$conv$opt,
                  messages_of_the_fit    = c(messages_of_the_fit, unlist(lmer_fit@optinfo$conv$lme4$messages)),
                  singular               = lme4::isSingular(lmer_fit),
                  seconds                = seconds))

}
