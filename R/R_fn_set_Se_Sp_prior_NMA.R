##
## ==== R_fn_set_Se_Sp_prior_NMA.R ==============================================================================
##
## ---- The direct Se / Sp prior of an NMA ordinal-bivariate model: priors$Se_Sp_prior is a list with one entry
##      per index test, each entry as the MA list of R_fn_set_Se_Sp_prior_MA.R (type, Sp, Se), all of one type.
##      Without priors$Se_Sp_prior the same fields are passed with use_Se_Sp_prior = 0: the model's own prior on
##      the cutpoints of every test.
##
#' R_fn_set_Se_Sp_prior_NMA
#' @keywords internal
#' @export
R_fn_set_Se_Sp_prior_NMA <-  function(priors,
                                      model_parameterisation,
                                      random_thresholds,
                                      n_index_tests,
                                      n_thr) {

      n_thr_max <-  max(n_thr)
      ##
      ## ---- the fields of every test with the prior off (the model's own prior on the cutpoints):
      ##
      Se_Sp_prior <-  priors[["Se_Sp_prior"]]
      priors[["Se_Sp_prior"]] <-  NULL
      priors$use_Se_Sp_prior <-  0
      priors$Se_Sp_prior_type <-  4
      priors$n_chosen_cutoffs_max <-  n_thr_max
      priors$n_chosen_cutoffs <-  cbind(n_thr, n_thr)
      priors$chosen_cutoffs <-  array(1L, dim = c(n_index_tests, 2, n_thr_max))
      priors$prior_Se_Sp_dirichlet_A <-  array(1, dim = c(n_index_tests, 2, n_thr_max + 1))
      priors$prior_Se_Sp_beta <-  array(0, dim = c(n_index_tests, 2, n_thr_max, 2))
      priors$prior_Se_Sp_link_mean_SD <-  array(0, dim = c(n_index_tests, 2, n_thr_max, 2))
      priors$prior_Se_Sp_link_mean_SD[, , , 2] <-  1
      for (t in 1:n_index_tests) {
            priors$chosen_cutoffs[t, 1, seq_len(n_thr[t])] <-  seq_len(n_thr[t])
            priors$chosen_cutoffs[t, 2, seq_len(n_thr[t])] <-  seq_len(n_thr[t])
      }
      if (is.null(Se_Sp_prior)) {
            return(priors)
      }
      ##
      ## ---- each test's own direct prior, built as the MA prior of that test and put in its row:
      ##
      stopifnot(is.list(Se_Sp_prior), length(Se_Sp_prior) == n_index_tests)
      type <-  NULL
      for (t in 1:n_index_tests) {
            of_the_test <-  R_fn_set_Se_Sp_prior_MA(priors = list(Se_Sp_prior = Se_Sp_prior[[t]],
                                                                 use_probit_link = priors$use_probit_link),
                                                   cts = FALSE,
                                                   model_parameterisation = model_parameterisation,
                                                   random_thresholds = random_thresholds,
                                                   HSROC_Dirichlet_prior = "shared_cutpoints",
                                                   n_thr = n_thr[t])
            if (is.null(type)) {
                  type <-  of_the_test$Se_Sp_prior_type
            }
            stopifnot(type == of_the_test$Se_Sp_prior_type)
            priors$n_chosen_cutoffs[t, ] <-  of_the_test$n_chosen_cutoffs
            for (c in 1:2) {
                  n_chosen_cutoffs_of_this_class <-  of_the_test$n_chosen_cutoffs[c]
                  priors$chosen_cutoffs[t, c, ] <-  1L
                  priors$chosen_cutoffs[t, c, seq_len(n_chosen_cutoffs_of_this_class)] <-
                        of_the_test$chosen_cutoffs[c, seq_len(n_chosen_cutoffs_of_this_class)]
                  priors$prior_Se_Sp_dirichlet_A[t, c, ] <-  1
                  priors$prior_Se_Sp_dirichlet_A[t, c, seq_len(n_chosen_cutoffs_of_this_class + 1)] <-
                        of_the_test$prior_Se_Sp_dirichlet_A[[c]][seq_len(n_chosen_cutoffs_of_this_class + 1)]
                  priors$prior_Se_Sp_beta[t, c, seq_len(n_chosen_cutoffs_of_this_class), ] <-
                        of_the_test$prior_Se_Sp_beta[[c]][seq_len(n_chosen_cutoffs_of_this_class), ]
                  priors$prior_Se_Sp_link_mean_SD[t, c, seq_len(n_chosen_cutoffs_of_this_class), ] <-
                        of_the_test$prior_Se_Sp_link_mean_SD[[c]][seq_len(n_chosen_cutoffs_of_this_class), ]
            }
      }
      priors$use_Se_Sp_prior <-  1
      priors$Se_Sp_prior_type <-  type
      ##
      return(priors)

}
