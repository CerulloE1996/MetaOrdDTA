##
## ---- Read a native prior bank from a sampled meta-analysis model:
##
## Use the full retained CmdStan fit and its recorded data file. The prepared model's mutable data
## list is not used as evidence of the data supplied to the sampler.
##


fn_accuracy_prior_bank_from_a_MetaOrdDTA_model <- function(model_samples_obj,
                                                           raw_cutoff_grid,
                                                           targets) {

    if (!inherits(model_samples_obj, "MetaOrd_model")) {
        stop("Supply the MetaOrd_model returned by directly sampling the native prior-only model.")
    }

    internal <- model_samples_obj$internal_obj
    native_fit <- internal$outs_stan_sampling$stan_mod_samples

    if (!inherits(native_fit, "CmdStanMCMC")) {
        stop(paste0("The model must retain its full CmdStanMCMC prior-only fit. ",
                    "A prepared model or the accuracy-checker's returned summaries cannot supply it."))
    }

    required_methods <- c("metadata", "return_codes", "data_file", "output_files", "code", "draws")

    for (method in required_methods) {

        if (!is.function(native_fit[[method]])) {
            stop(paste0("The retained native fit lacks the public ", method, "() method."))
        }

    }

    stages <- c("Check the retained native fit",
                "Read the actual sampling data",
                "Read all paired native prior draws",
                "Connect the requested accuracy targets")
    progress <- utils::txtProgressBar(min = 0, max = length(stages), initial = 0, style = 3)
    on.exit(close(progress), add = TRUE)
    message(stages[1])

    metadata <- native_fit$metadata()
    return_codes <- native_fit$return_codes()

    if (!is.list(metadata) || !identical(metadata$method, "sample")) {
        stop("The retained native fit must record a sampling run.")
    }

    if (!is.numeric(return_codes) || length(return_codes) == 0
        || any(!is.finite(return_codes)) || any(return_codes != 0)) {
        stop("Every process in the retained native fit must have a recorded successful return code.")
    }

    metadata_to_model_code <- c(
        DTA_MA_ord_biv_FIXEDthr_model = "DTA_MA_ord_biv_FIXEDthr",
        DTA_MA_Xu_FIXEDthr_model = "DTA_MA_ord_biv_FIXEDthr",
        DTA_MA_ord_HSROC_RANDthr_kappa_model = "DTA_MA_ord_HSROC_RANDthr_kappa",
        DTA_MA_Gat_RANDthr_kappa_model = "DTA_MA_ord_HSROC_RANDthr_kappa",
        DTA_MA_Jones_model = "DTA_MA_Jones"
    )
    native_model_name <- metadata$model_name

    if (!is.character(native_model_name) || length(native_model_name) != 1
        || is.na(native_model_name) || !(native_model_name %in% names(metadata_to_model_code))) {
        stop(paste0("The retained fit must identify an ordinal-bivariate fixed-cutpoint, ordinal-HSROC ",
                    "random-cutpoint or Jones model."))
    }

    model_code <- unname(metadata_to_model_code[native_model_name])
    recorded_model_file <- internal$outs_stan_compile$stan_model_file_path

    if (!is.character(recorded_model_file) || length(recorded_model_file) != 1
        || is.na(recorded_model_file) || !nzchar(recorded_model_file)) {
        stop("The sampled model must retain the matching original native Stan source file.")
    }

    legacy_source_to_current <- c(
        DTA_MA_Xu_FIXEDthr.stan = "DTA_MA_ord_biv_FIXEDthr.stan",
        DTA_MA_Gat_RANDthr_kappa.stan = "DTA_MA_ord_HSROC_RANDthr_kappa.stan"
    )
    recorded_basename <- basename(recorded_model_file)
    expected_basename <- paste0(model_code, ".stan")

    if (identical(recorded_basename, expected_basename)) {
        model_file <- recorded_model_file
    } else if (recorded_basename %in% names(legacy_source_to_current)
               && identical(unname(legacy_source_to_current[recorded_basename]), expected_basename)) {
        model_file <- file.path(dirname(recorded_model_file),
                                unname(legacy_source_to_current[recorded_basename]))
    } else {
        stop("The sampled model must retain the matching original native Stan source file.")
    }

    if (!file.exists(model_file)) {
        stop("The current renamed native Stan source file is unavailable.")
    }

    model_file <- normalizePath(path = model_file, winslash = "/", mustWork = TRUE)
    retained_code <- native_fit$code()

    if (!is.character(retained_code) || length(retained_code) == 0 || anyNA(retained_code)
        || !identical(unname(retained_code), readLines(con = model_file, warn = FALSE))) {
        stop("The retained fit's top-level Stan source must match the recorded native source file.")
    }

    source_limitations <- paste0("The retained fit identifies its top-level Stan source, ",
                                 "not the historical contents of recursively included files.")
    utils::setTxtProgressBar(pb = progress, value = 1)
    message(stages[2])

    data_file <- native_fit$data_file()

    if (!is.character(data_file) || length(data_file) != 1 || is.na(data_file)
        || !nzchar(data_file) || !file.exists(data_file)) {
        stop(paste0("The actual data file recorded by the native fit is unavailable. ",
                    "It cannot be replaced by the prepared model's current data list."))
    }

    data_file <- normalizePath(path = data_file, winslash = "/", mustWork = TRUE)

    if (!requireNamespace(package = "jsonlite", quietly = TRUE)) {
        stop("The jsonlite package is required to read the native fit's recorded data file.")
    }

    native_stan_data <- jsonlite::fromJSON(txt = data_file,
                                           simplifyVector = TRUE,
                                           simplifyDataFrame = FALSE,
                                           simplifyMatrix = TRUE,
                                           flatten = FALSE)

    if (!is.list(native_stan_data)
        || !all(c("n_studies", "K_fold_CV_indicator") %in% names(native_stan_data))
        || !is.numeric(native_stan_data$n_studies) || length(native_stan_data$n_studies) != 1
        || !is.finite(native_stan_data$n_studies) || native_stan_data$n_studies < 1
        || !is.numeric(native_stan_data$K_fold_CV_indicator)
        || length(native_stan_data$K_fold_CV_indicator) != native_stan_data$n_studies
        || anyNA(native_stan_data$K_fold_CV_indicator)
        || any(native_stan_data$K_fold_CV_indicator != 0)) {
        stop("The recorded data must identify a prior-only fit for every study.")
    }

    output_files <- native_fit$output_files(include_failed = TRUE)

    if (!is.character(output_files) || length(output_files) == 0 || anyNA(output_files)
        || any(!nzchar(output_files)) || any(!file.exists(output_files)) || anyDuplicated(output_files)) {
        stop("Every recorded native sampling CSV must remain available, without duplicated file paths.")
    }

    output_files <- normalizePath(path = output_files, winslash = "/", mustWork = TRUE)
    utils::setTxtProgressBar(pb = progress, value = 2)
    message(stages[3])

    native_prior_draws <- native_fit$draws(variables = NULL, inc_warmup = FALSE, format = "matrix")
    utils::setTxtProgressBar(pb = progress, value = 3)
    message(stages[4])

    bank <- fn_accuracy_prior_bank_from_MetaOrdDTA(native_prior_draws = native_prior_draws,
                                                   native_stan_data = native_stan_data,
                                                   model_code = model_code,
                                                   raw_cutoff_grid = raw_cutoff_grid,
                                                   targets = targets,
                                                   model_file = model_file,
                                                   draw_source = paste(output_files, collapse = "\n"))

    bank$provenance$native_fit <- list(method = metadata$method,
                                       model_name = native_model_name,
                                       recorded_model_file = recorded_model_file,
                                       data_file = data_file,
                                       output_files = output_files,
                                       return_codes = return_codes,
                                       top_level_source_matches = TRUE,
                                       source_limitations = source_limitations)

    utils::setTxtProgressBar(pb = progress, value = length(stages))
    message(paste0("Native prior bank: ", model_code, "; ", nrow(native_prior_draws), " paired draws."))
    message(source_limitations)

    bank

}






















