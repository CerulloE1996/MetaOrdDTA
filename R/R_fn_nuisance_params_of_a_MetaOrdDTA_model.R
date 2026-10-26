##
## ==== R_fn_nuisance_params_of_a_MetaOrdDTA_model.R ===========================================================
##
## ---- The nuisance blocks of each of MetaOrdDTA's own DTA models, named by the package: every z block,
##      everything that grows with N, and C_raw only in the random-cutpoint models (..._RANDthr_kappa).
##
#' The nuisance blocks of one of MetaOrdDTA's own DTA models
#'
#' Returns the names of the parameter blocks NicoStan holds as the nuisance block for the model: every z
#' block, everything that grows with N, and C_raw only in the random-cutpoint models.
#'
#' @param name_of_the_model The model's name without ".stan", e.g. "DTA_MA_ord_biv_FIXEDthr" or
#'   "DTA_NMA_Nyaga_ord_biv_FIXEDthr".
#' @param model_parameterisation "Xu", "ordinal_HSROC", "Jones", "Bivariate" or "NMA" (used in the message of a model whose
#'   blocks are not known).
#' @param random_thresholds TRUE for the random-cutpoint models, FALSE for the fixed-cutpoint ones.
#' @return A character vector of nuisance block names, never NULL. A model whose blocks are not known stops, with
#'   its own name and its own Stan file in the message.
#' @export
fn_nuisance_params_of_a_MetaOrdDTA_model <-  function( name_of_the_model,
                                                       model_parameterisation,
                                                       random_thresholds ) {

      stopifnot(length(name_of_the_model) == 1, nzchar(name_of_the_model))
      stopifnot(length(model_parameterisation) == 1, nzchar(model_parameterisation))
      stopifnot(length(random_thresholds) == 1, is.logical(random_thresholds))
      ##
      ## ---- the MA models: ord_biv, ord_HSROC and Jones
      if (grepl("^DTA_MA_ord_biv_FIXEDthr", name_of_the_model)) {
            return(c("beta_z"))
      }
      if (grepl("^DTA_MA_ord_biv_RANDthr", name_of_the_model) && grepl("_centred_C_", name_of_the_model)) {
            return(c("beta_z", "C_v"))
      }
      if (grepl("^DTA_MA_ord_biv_RANDthr", name_of_the_model)) {
            return(c("beta_z", "C_raw"))
      }
      if (grepl("^DTA_MA_ord_HSROC_FIXEDthr", name_of_the_model)) {
            return(c("beta_z", "raw_scale_z"))
      }
      if (grepl("^DTA_MA_ord_HSROC_RANDthr", name_of_the_model)) {
            return(c("beta_z", "raw_scale_z", "C_raw"))
      }
      if (grepl("^DTA_MA_Jones", name_of_the_model)) {
            return(c("beta_z", "raw_scale_z"))
      }
      if (grepl("^DTA_MA_strat_bivariate", name_of_the_model)) {
            return(c("z"))
      }
      ##
      ## ---- the NMA models
      if (grepl("^DTA_NMA_Nyaga_ord_biv_FIXEDthr", name_of_the_model)) {
            return(c("beta_eta_z", "beta_delta_z_vec"))
      }
      if (grepl("^DTA_NMA_Nyaga_ord_biv_RANDthr", name_of_the_model)) {
            return(c("beta_eta_z", "beta_delta_z_vec", "C_raw_vec"))
      }
      if (grepl("^DTA_NMA_Nyaga_Jones", name_of_the_model)) {
            return(c("log_beta_sigma_z", "beta_eta_z", "beta_delta_z", "log_raw_scale_sigma_z",
                     "raw_scale_eta_z", "raw_scale_delta_z"))
      }
      stop(paste0("fn_nuisance_params_of_a_MetaOrdDTA_model: no nuisance blocks are known for the model \"",
                  name_of_the_model, "\" (parameterisation \"", model_parameterisation,
                  "\", random_thresholds = ", random_thresholds,
                  "). Its own Stan file must name them first."))
}






















