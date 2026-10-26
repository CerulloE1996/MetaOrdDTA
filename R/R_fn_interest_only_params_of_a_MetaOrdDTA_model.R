##
## ==== R_fn_interest_only_params_of_a_MetaOrdDTA_model.R ======================================================
##
## ---- The parameter families of interest of each of MetaOrdDTA's own DTA models, chosen by the package, never by
##      the user: the MAIN-block parameter families that the summary accuracy (Se_baseline, Sp_baseline) is built
##      from. NicoStan's trajectory-length (tau) criterion monitors only these families during burn-in
##      (NicoStan's $sample() argument interest_only, with tau_adaptation_block "main"), and the min ESS of the
##      summary is computed over them. A nuisance block (z blocks; C_raw in the random-cutpoint models) is never
##      listed: in the random-cutpoint models the summary cutpoints are the median over the studies of C, which
##      is built from the nuisance block C_raw, so the main-block family that sets where those cutpoints sit,
##      dirichlet_phi, is listed in its place.
##
#' The parameter families of interest of one of MetaOrdDTA's own DTA models
#'
#' Returns the names of the main-block parameter families from which Se_baseline and Sp_baseline are computed,
#' for NicoStan's interest_only argument (the families its tau criterion monitors and its min ESS is computed
#' over).
#'
#' @param name_of_the_model The model's name without ".stan", e.g. "DTA_MA_ord_biv_FIXEDthr".
#' @param model_parameterisation "Xu", "ordinal_HSROC" (or the older "R&G"), "Jones", "Bivariate" or "NMA"
#'   (used in the message of a model whose families are not known).
#' @param random_thresholds TRUE for the random-cutpoint models, FALSE for the fixed-cutpoint ones.
#' @return A character vector of parameter family names, or NULL for a model whose families are not known (then
#'   NicoStan monitors every main-block row, as before).
#' @export
fn_interest_only_params_of_a_MetaOrdDTA_model <-  function( name_of_the_model,
                                                            model_parameterisation,
                                                            random_thresholds ) {

      stopifnot(length(name_of_the_model) == 1, nzchar(name_of_the_model))
      stopifnot(length(model_parameterisation) == 1, nzchar(model_parameterisation))
      stopifnot(length(random_thresholds) == 1, is.logical(random_thresholds))
      ##
      ## ---- the MA models: ord_biv, ord_HSROC, Jones and stratified bivariate
      ##
      ## (ord_biv fixed-cutpoint: Se / Sp from C, C_raw_vec and beta_mu)
      if (grepl("^DTA_MA_ord_biv_FIXEDthr", name_of_the_model)) {
            return(c("beta_mu", "C_raw_vec"))
      }
      ## (ord_biv random-cutpoint: Se / Sp from C_mu, median C and beta_mu; dirichlet_phi sets C)
      if (grepl("^DTA_MA_ord_biv_RANDthr", name_of_the_model)) {
            return(c("beta_mu", "dirichlet_phi"))
      }
      ## (ordinal-HSROC fixed-cutpoint: Se / Sp from C (C_raw_vec), beta_mu and the scale from raw_scale_mu)
      if (grepl("^DTA_MA_ord_HSROC_FIXEDthr", name_of_the_model)) {
            return(c("beta_mu", "raw_scale_mu", "C_raw_vec"))
      }
      ## (ordinal-HSROC random-cutpoint: Se / Sp from C_mu (nuisance C_raw), beta_mu and raw_scale_mu)
      if (grepl("^DTA_MA_ord_HSROC_RANDthr", name_of_the_model)) {
            return(c("beta_mu", "raw_scale_mu", "dirichlet_phi"))
      }
      ## (Jones: Se / Sp from beta_mu, raw_scale_mu and the two Box-Cox lambdas)
      if (grepl("^DTA_MA_Jones", name_of_the_model)) {
            return(c("beta_mu", "raw_scale_mu", "lambda", "lambda_upper"))
      }
      ## (stratified bivariate: Se / Sp from mu at each threshold)
      if (grepl("^DTA_MA_strat_bivariate", name_of_the_model)) {
            return(c("mu"))
      }
      ##
      ## ---- the NMA models
      ##
      ## (Nyaga ord_biv fixed-cutpoint: Se / Sp from C (C_raw_vec) and beta_mu)
      if (grepl("^DTA_NMA_Nyaga_ord_biv_FIXEDthr", name_of_the_model)) {
            return(c("beta_mu", "C_raw_vec"))
      }
      ## (Nyaga ord_biv random-cutpoint: Se / Sp from median C (nuisance C_raw_vec) and beta_mu;
      ##  raw_dirichlet_phi_vec sets C)
      if (grepl("^DTA_NMA_Nyaga_ord_biv_RANDthr", name_of_the_model)) {
            return(c("beta_mu", "raw_dirichlet_phi_vec"))
      }
      ##
      message(paste0("\033[33m", "fn_interest_only_params_of_a_MetaOrdDTA_model: no parameter families of interest ",
                     "are known for the model \"", name_of_the_model, "\" (parameterisation \"",
                     model_parameterisation, "\", random_thresholds = ", random_thresholds,
                     "): NicoStan monitors every main-block row.", "\033[0m"))
      return(NULL)

}






















