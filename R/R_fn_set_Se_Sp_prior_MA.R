##
## ==== R_fn_set_Se_Sp_prior_MA.R ===============================================================================
##
## ---- The direct Se / Sp prior of an MA model from the user's medians and 95% limits at 1-2 cutoffs per class:
##      priors$Se_Sp_prior = list(Sp = list(cutoffs, median, lower, upper), Se = list(cutoffs, median, lower,
##      upper)); a Normal on the link scale at those cutoffs (mean = qnorm(median), SD from the 95% limits),
##      nothing typed at the other cutoffs. Class order: Sp then Se in the ordinal and Jones models, Se then Sp in
##      the stratified-bivariate model. priors$Se_Sp_prior$type picks the prior (Stan's Se_Sp_prior_type):
##      4 (default) the Normal on the link scale from median / lower / upper; 1 a Dirichlet on the segment sums
##      (each class gives A, a vector of length n_chosen_cutoffs + 1); 2 a Beta on each Sp_k / Se_k and 3 a Beta
##      on each conditional step (each class gives beta, an n_chosen_cutoffs x 2 matrix of a and b). Without
##      priors$Se_Sp_prior the same fields are passed with use_Se_Sp_prior = 0: the model's own prior on the
##      cutpoints; except the ordinal-HSROC random-cutpoint model with HSROC_Dirichlet_prior = "both_classes" (an
##      option; the default "shared_ordinal_priors" keeps the Dirichlet on the one shared cutpoint set), whose
##      Dirichlet(prior_dirichlet_alpha) goes on the induced probabilities of both classes at every cutoff (type 1).
##
#' R_fn_set_Se_Sp_prior_MA
#' @keywords internal
#' @export
R_fn_set_Se_Sp_prior_MA <-  function(priors,
                                     cts,
                                     model_parameterisation,
                                     random_thresholds,
                                     HSROC_Dirichlet_prior,
                                     n_thr) {

      ##
      ## ---- the fields of the Se / Sp prior converted to the Jones model's parameters (off unless set):
      ##
      if (is.null(priors$use_converted_Se_Sp_prior)) priors$use_converted_Se_Sp_prior <-  0
      if (is.null(priors$prior_converted_mean)) priors$prior_converted_mean <-  list(rep(0, 3), rep(0, 3))
      if (is.null(priors$prior_converted_SD)) priors$prior_converted_SD <-  list(rep(1, 3), rep(1, 3))
      if (is.null(priors$prior_converted_L_Omega)) priors$prior_converted_L_Omega <-  list(diag(3), diag(3))
      ##
      ## ---- without priors$Se_Sp_prior: the fields with the prior off, the model's own prior on the cutpoints:
      ##
      if (is.null(priors[["Se_Sp_prior"]])) {
            is_the_ord_HSROC <-  grepl("ord_HSROC|ordinal_HSROC|R&G|HSROC|Gatsonis", model_parameterisation, ignore.case = TRUE)
            on_both_classes <-  is_the_ord_HSROC && isTRUE(random_thresholds) &&
                                identical(HSROC_Dirichlet_prior, "both_classes")
            priors$use_Se_Sp_prior <-  if (on_both_classes) 1 else 0
            priors$Se_Sp_prior_type <-  if (on_both_classes) 1 else 4
            priors$n_chosen_cutoffs_max <-  n_thr
            priors$n_chosen_cutoffs <-  c(n_thr, n_thr)
            priors$chosen_cutoffs <-  rbind(1:n_thr, 1:n_thr)
            alpha_of_both_classes <-  if (on_both_classes) priors$prior_dirichlet_alpha else rep(1, n_thr + 1)
            stopifnot(length(alpha_of_both_classes) == n_thr + 1)
            priors$prior_Se_Sp_dirichlet_A <-  list(alpha_of_both_classes, alpha_of_both_classes)
            priors$prior_Se_Sp_beta <-  list(matrix(0, n_thr, 2), matrix(0, n_thr, 2))
            priors$prior_Se_Sp_link_mean_SD <-  list(cbind(rep(0, n_thr), rep(1, n_thr)),
                                                    cbind(rep(0, n_thr), rep(1, n_thr)))
            return(priors)
      }
      ##
      ## ---- the two classes in the order of the Stan file ([[ ]] throughout: $ matches Se_Sp_prior_type too):
      ##
      strat_biv <-  grepl("strat_biv|Bivariate", model_parameterisation, ignore.case = TRUE)
      classes <-  if (strat_biv) c("Se", "Sp") else c("Sp", "Se")
      type <-  if (is.null(priors[["Se_Sp_prior"]][["type"]])) 4 else priors[["Se_Sp_prior"]][["type"]]
      stopifnot(type %in% 1:4)
      fields_of_the_type <-  switch(as.character(type),
                                    "1" = c("cutoffs", "A"),
                                    "2" = c("cutoffs", "beta"),
                                    "3" = c("cutoffs", "beta"),
                                    "4" = c("cutoffs", "median", "lower", "upper"))
      for (class in classes) {
            stopifnot(!is.null(priors[["Se_Sp_prior"]][[class]]),
                      all(fields_of_the_type %in% names(priors[["Se_Sp_prior"]][[class]])))
      }
      n_chosen_cutoffs <-  sapply(classes, function(class) length(priors[["Se_Sp_prior"]][[class]]$cutoffs))
      n_chosen_cutoffs_max <-  max(n_chosen_cutoffs)
      stopifnot(n_chosen_cutoffs_max >= 1, n_chosen_cutoffs_max <= n_thr)
      ##
      ## ---- the Normal on the link scale (probit or logit) of each class at its cutoffs:
      ##
      fn_link <-  if (priors$use_probit_link == 1) stats::qnorm else stats::qlogis
      fn_link_mean_and_SD <-  function(prior_of_the_class) {

            cbind(mean = fn_link(prior_of_the_class$median),
                  SD = (fn_link(prior_of_the_class$upper) - fn_link(prior_of_the_class$lower)) /
                       (2 * stats::qnorm(0.975)))

      }
      chosen_cutoffs <-  matrix(1L, 2, n_chosen_cutoffs_max)
      prior_Se_Sp_link_mean_SD <-  list(matrix(0, n_chosen_cutoffs_max, 2), matrix(0, n_chosen_cutoffs_max, 2))
      dirichlet_A <-  list(rep(1, n_chosen_cutoffs_max + 1), rep(1, n_chosen_cutoffs_max + 1))
      prior_Se_Sp_beta <-  list(matrix(0, n_chosen_cutoffs_max, 2), matrix(0, n_chosen_cutoffs_max, 2))
      for (c in 1:2) {
            prior_of_the_class <-  priors[["Se_Sp_prior"]][[classes[c]]]
            stopifnot(all(prior_of_the_class$cutoffs %in% 1:n_thr))
            chosen_cutoffs[c, seq_len(n_chosen_cutoffs[c])] <-  as.integer(prior_of_the_class$cutoffs)
            if (type == 4) {
                  stopifnot(all(prior_of_the_class$lower < prior_of_the_class$median),
                            all(prior_of_the_class$median < prior_of_the_class$upper))
                  prior_Se_Sp_link_mean_SD[[c]][seq_len(n_chosen_cutoffs[c]), ] <-
                        fn_link_mean_and_SD(prior_of_the_class)
            } else if (type == 1) {
                  stopifnot(length(prior_of_the_class$A) == n_chosen_cutoffs[c] + 1,
                            all(prior_of_the_class$A > 0))
                  dirichlet_A[[c]][seq_len(n_chosen_cutoffs[c] + 1)] <-  prior_of_the_class$A
            } else {
                  stopifnot(all(dim(prior_of_the_class$beta) == c(n_chosen_cutoffs[c], 2)),
                            all(prior_of_the_class$beta > 0))
                  prior_Se_Sp_beta[[c]][seq_len(n_chosen_cutoffs[c]), ] <-  prior_of_the_class$beta
            }
      }
      ##
      ## ---- the fields of the _Se_Sp_prior Stan files (the user's list itself is not Stan data):
      ##
      priors$use_Se_Sp_prior <-  1
      priors$Se_Sp_prior_type <-  type
      priors$n_chosen_cutoffs_max <-  n_chosen_cutoffs_max
      priors$n_chosen_cutoffs <-  unname(n_chosen_cutoffs)
      priors$chosen_cutoffs <-  chosen_cutoffs
      priors$prior_Se_Sp_dirichlet_A <-  dirichlet_A
      priors$prior_Se_Sp_beta <-  prior_Se_Sp_beta
      priors$prior_Se_Sp_link_mean_SD <-  prior_Se_Sp_link_mean_SD
      ##
      cutoffs_of_each_class <-  sapply(classes, function(class) {
            paste(priors[["Se_Sp_prior"]][[class]]$cutoffs, collapse = ", ")
      })
      name_of_the_type <-  c("Dirichlet on the segment sums", "Beta on each Se / Sp", "Beta on each step",
                             "Normal on the link scale")[type]
      message(paste0("Direct Se/Sp prior: ",
                     paste0(classes, " at cutoffs ", cutoffs_of_each_class, collapse = "; "),
                     " (", name_of_the_type, ", type ", type, "); the other cutoffs follow the model's own",
                     " prior"))
      priors[["Se_Sp_prior"]] <-  NULL
      ##
      return(priors)

}






















