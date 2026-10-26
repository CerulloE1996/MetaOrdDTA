##
## ==== R_fn_reason_of_a_nuisance_block_of_a_MetaOrdDTA_model.R ================================================
##
## ---- Why each block of a MetaOrdDTA model is a nuisance block, printed with the blocks before every fit
##      (the blocks themselves come from fn_nuisance_params_of_a_MetaOrdDTA_model()).
##
#' Why one block of a MetaOrdDTA model is a nuisance block
#'
#' @param name_of_a_block The block's name, as returned by fn_nuisance_params_of_a_MetaOrdDTA_model().
#' @param name_of_the_model The model's name without ".stan".
#' @return One sentence saying why the block is nuisance.
#' @export
fn_reason_of_a_nuisance_block_of_a_MetaOrdDTA_model <-  function( name_of_a_block,
                                                                  name_of_the_model ) {

      stopifnot(length(name_of_a_block) == 1, nzchar(name_of_a_block))
      stopifnot(length(name_of_the_model) == 1, nzchar(name_of_the_model))
      inside_an_NMA_model <-  grepl("^DTA_NMA_", name_of_the_model)
      reason <-  switch(name_of_a_block,
            beta_z                 = "the study-specific random effects (a z block, one per study)",
            raw_scale_z            = "the study-specific raw-scale random effects (a z block, one per study)",
            C_raw                  = "the study-specific cutpoints (a random-cutpoint model)",
            C_raw_vec              = "the study-specific cutpoints (a random-cutpoint model)",
            z                      = "the per-threshold, per-study random effects (a z block)",
            beta_eta_z             = "the test- and study-specific random effects (a z block, one per study)",
            beta_delta_z_vec       = "the test-specific effects (a z block, one per observation)",
            beta_delta_z           = "the test- and study-specific effects (a z block, one per study)",
            raw_scale_eta_z        = "the test- and study-specific raw-scale effects (a z block, one per study)",
            raw_scale_delta_z      = "the test- and study-specific raw-scale effects (a z block, one per study)",
            log_beta_sigma_z       = "the between-study log SDs of the tests (a z block, one per test)",
            log_raw_scale_sigma_z  = "the between-study raw-scale log SDs of the tests (a z block, one per test)",
            paste0("a nuisance block of ", if (inside_an_NMA_model) "the NMA model" else "the model",
                   " (its name ends in _z, or it is a random-cutpoint block)"))
      return(reason)
}






















