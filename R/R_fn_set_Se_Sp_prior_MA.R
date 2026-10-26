##
## ==== R_fn_set_Se_Sp_prior_MA.R ===============================================================================
##
## ---- The direct Se / Sp prior of an MA model from the user's medians and 95% limits at 1-2 cutoffs per class:
##      priors$Se_Sp_prior = list(Sp = list(cutoffs, median, lower, upper), Se = list(cutoffs, median, lower,
##      upper)); a Normal on the link scale at those cutoffs (mean = qnorm(median), SD from the 95% limits),
##      nothing typed at the other cutoffs. Class order: Sp then Se in the ordinal and Jones models, Se then Sp in
##      the stratified-bivariate model. Without priors$Se_Sp_prior the priors are returned as they are.
##
#' R_fn_set_Se_Sp_prior_MA
#' @keywords internal
#' @export
R_fn_set_Se_Sp_prior_MA <-  function(priors,
                                     cts,
                                     model_parameterisation,
                                     n_thr) {

      if (is.null(priors[["Se_Sp_prior"]])) {
            return(priors)
      }
      ##
      ## ---- the two classes in the order of the Stan file ([[ ]] throughout: $ matches Se_Sp_prior_type too):
      ##
      strat_biv <-  grepl("strat_biv|Bivariate", model_parameterisation, ignore.case = TRUE)
      classes <-  if (strat_biv) c("Se", "Sp") else c("Sp", "Se")
      for (class in classes) {
            stopifnot(!is.null(priors[["Se_Sp_prior"]][[class]]),
                      all(c("cutoffs", "median", "lower", "upper") %in% names(priors[["Se_Sp_prior"]][[class]])))
      }
      n_chosen <-  sapply(classes, function(class) length(priors[["Se_Sp_prior"]][[class]]$cutoffs))
      n_chosen_max <-  max(n_chosen)
      stopifnot(n_chosen_max >= 1, n_chosen_max <= n_thr)
      ##
      ## ---- the Normal on the link scale (probit or logit) of each class at its cutoffs:
      ##
      fn_link <-  if (priors$use_probit_link == 1) stats::qnorm else stats::qlogis
      fn_link_mean_and_SD <-  function(prior_of_the_class) {

            cbind(mean = fn_link(prior_of_the_class$median),
                  SD = (fn_link(prior_of_the_class$upper) - fn_link(prior_of_the_class$lower)) /
                       (2 * stats::qnorm(0.975)))

      }
      chosen_cutoffs <-  matrix(1L, 2, n_chosen_max)
      link_mean_SD <-  list(matrix(0, n_chosen_max, 2), matrix(0, n_chosen_max, 2))
      for (c in 1:2) {
            prior_of_the_class <-  priors[["Se_Sp_prior"]][[classes[c]]]
            stopifnot(all(prior_of_the_class$cutoffs %in% 1:n_thr),
                      all(prior_of_the_class$lower < prior_of_the_class$median),
                      all(prior_of_the_class$median < prior_of_the_class$upper))
            chosen_cutoffs[c, seq_len(n_chosen[c])] <-  as.integer(prior_of_the_class$cutoffs)
            link_mean_SD[[c]][seq_len(n_chosen[c]), ] <-  fn_link_mean_and_SD(prior_of_the_class)
      }
      ##
      ## ---- the fields of the _Se_Sp_prior Stan files (the user's list itself is not Stan data):
      ##
      priors$use_Se_Sp_prior <-  1
      priors$Se_Sp_prior_type <-  4
      priors$n_chosen_max <-  n_chosen_max
      priors$n_chosen <-  unname(n_chosen)
      priors$chosen_cutoffs <-  chosen_cutoffs
      priors$prior_Se_Sp_beta <-  list(matrix(0, n_chosen_max, 2), matrix(0, n_chosen_max, 2))
      priors$prior_Se_Sp_link_mean_SD <-  link_mean_SD
      if (!strat_biv && cts == FALSE) {
            priors$prior_Se_Sp_dirichlet_A <-  list(rep(0, n_chosen_max + 1), rep(0, n_chosen_max + 1))
      }
      ##
      cutoffs_of_each_class <-  sapply(classes, function(class) {
            paste(priors[["Se_Sp_prior"]][[class]]$cutoffs, collapse = ", ")
      })
      message(paste0("Direct Se/Sp prior: ",
                     paste0(classes, " at cutoffs ", cutoffs_of_each_class, collapse = "; "),
                     " (Normal on the link scale, type 4); the other cutoffs follow the model's own prior"))
      priors[["Se_Sp_prior"]] <-  NULL
      ##
      return(priors)

}






















