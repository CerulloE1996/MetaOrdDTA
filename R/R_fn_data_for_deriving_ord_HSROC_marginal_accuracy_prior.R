##
## ---- Retain the paired native random-cutpoint ord_HSROC integration bank:
##


fn_ord_HSROC_native_prior_integration_bank <- function(native_prior_draws,
                                                       native_stan_data,
                                                       raw_cutoff_grid) {

    if (!is.matrix(native_prior_draws) || !is.numeric(native_prior_draws) || nrow(native_prior_draws) < 2
        || is.null(colnames(native_prior_draws)) || anyDuplicated(colnames(native_prior_draws))) {
        stop("Supply complete paired native prior draws as a numeric matrix with unique column names.")
    }

    required <- c("n_studies", "n_thr", "n_covariates", "baseline_case", "prior_beta_mu_mean",
                  "prior_beta_mu_SD", "prior_raw_scale_mu_mean", "prior_raw_scale_mu_SD",
                  "prior_dirichlet_alpha", "softplus", "use_probit_link", "K_fold_CV_indicator")

    if (!is.list(native_stan_data) || !all(required %in% names(native_stan_data))) {
        stop("Supply the actual native ord_HSROC data, including every baseline and prior input.")
    }

    for (name in c("n_studies", "n_thr", "n_covariates")) {
        value <- native_stan_data[[name]]
        if (!is.numeric(value) || length(value) != 1 || !is.finite(value)
            || value < 1 || value != as.integer(value)) {
            stop(paste0(name, " must be the actual positive native integer."))
        }
    }

    n_thr <- native_stan_data$n_thr
    n_covariates <- native_stan_data$n_covariates
    indicators <- native_stan_data$K_fold_CV_indicator
    if (!is.numeric(indicators) || length(indicators) != native_stan_data$n_studies
        || anyNA(indicators) || any(indicators != 0)) {
        stop("Every recorded study likelihood indicator must be zero.")
    }

    if (!is.numeric(raw_cutoff_grid) || !is.null(dim(raw_cutoff_grid)) || length(raw_cutoff_grid) != n_thr
        || any(!is.finite(raw_cutoff_grid)) || any(diff(raw_cutoff_grid) <= 0)) {
        stop("Supply the complete strictly increasing raw-cutoff grid for the actual native model.")
    }

    for (name in c("softplus", "use_probit_link")) {
        value <- native_stan_data[[name]]
        if (!(is.numeric(value) || is.logical(value)) || length(value) != 1
            || !is.finite(value) || !(value %in% c(0, 1))) {
            stop(paste0(name, " must be the actual native zero-or-one flag."))
        }
    }

    C_names <- paste0("C_mu[", seq_len(n_thr), "]")
    phi_names <- paste0("dirichlet_phi[", seq_len(n_thr + 1), "]")
    required_draws <- c(C_names, phi_names, "Xraw_scale_baseline")
    if (!all(required_draws %in% colnames(native_prior_draws))) {
        stop("Retain actual C_mu medians, dirichlet_phi and Xraw_scale_baseline from the native fit.")
    }

    cutpoint_bank <- unname(native_prior_draws[, C_names, drop = FALSE])
    phi_bank <- unname(native_prior_draws[, phi_names, drop = FALSE])
    raw_scale_bank <- as.numeric(native_prior_draws[, "Xraw_scale_baseline"])
    if (any(!is.finite(cutpoint_bank)) || any(!is.finite(phi_bank))
        || any(phi_bank <= 0 | phi_bank > 1) || any(!is.finite(raw_scale_bank))) {
        stop("Native bank entries must be finite, with strictly positive simplex entries; no rows are dropped.")
    }
    for (j in seq_len(nrow(cutpoint_bank))) {
        if (any(diff(cutpoint_bank[j, ]) < 0)) stop("Native C_mu medians must retain their cutoff order.")
    }

    baseline <- native_stan_data$baseline_case
    if (!is.numeric(baseline) || !is.null(dim(baseline)) || length(baseline) != n_covariates
        || any(!is.finite(baseline))) {
        stop("Supply the actual native baseline vector, matching its covariate count.")
    }

    projections <- list()
    for (block in c("beta", "raw_scale")) {
        means <- native_stan_data[[paste0("prior_", block, "_mu_mean")]]
        SDs <- native_stan_data[[paste0("prior_", block, "_mu_SD")]]
        if (!is.numeric(means) || !is.null(dim(means)) || length(means) != n_covariates
            || !is.numeric(SDs) || !is.null(dim(SDs)) || length(SDs) != n_covariates
            || any(!is.finite(means)) || any(!is.finite(SDs)) || any(SDs <= 0)) {
            stop("Each active native normal coefficient-prior vector must match the baseline vector.")
        }
        projections[[paste0(block, "_mean")]] <- sum(baseline * means)
        projections[[paste0(block, "_SD")]] <- sqrt(sum((baseline * SDs)^2))
    }
    if (any(!is.finite(unlist(projections))) || projections$beta_SD <= 0 || projections$raw_scale_SD <= 0) {
        stop("The actual baseline normal projections must have finite means and strictly positive SDs.")
    }

    alpha <- native_stan_data$prior_dirichlet_alpha
    if (!is.numeric(alpha) || !is.null(dim(alpha)) || length(alpha) != n_thr + 1
        || any(!is.finite(alpha)) || any(alpha <= 0)) {
        stop("Supply every actual native prior_dirichlet_alpha concentration, without substitution.")
    }

    message(paste0("Native ord_HSROC bank: ", nrow(cutpoint_bank), " paired rows, ", n_thr,
                   " cutoffs and ", native_stan_data$n_studies, " studies. No bank values are reconstructed."))

    list(n_thr = as.integer(n_thr),
         raw_cutoff_grid = unname(raw_cutoff_grid),
         n_bank = nrow(cutpoint_bank),
         cutpoint_bank = cutpoint_bank,
         phi_bank = phi_bank,
         raw_scale_bank = raw_scale_bank,
         bank_raw_scale_mean = projections$raw_scale_mean,
         bank_raw_scale_SD = projections$raw_scale_SD,
         bank_dirichlet_alpha = unname(alpha),
         softplus = as.integer(native_stan_data$softplus),
         use_probit_link = as.integer(native_stan_data$use_probit_link),
         baseline_case = baseline,
         n_covariates = as.integer(n_covariates),
         n_studies = as.integer(native_stan_data$n_studies),
         reference_beta_baseline_mean = projections$beta_mean,
         reference_beta_baseline_SD = projections$beta_SD,
         native_stan_data = native_stan_data,
         scope = "native R&G baseline marginal accuracy; not predictive accuracy or joint curve matching")

}


##
## ---- Prepare the explicitly selected ord_HSROC marginal composite likelihood:
##


fn_data_for_deriving_ord_HSROC_marginal_accuracy_prior <- function(reference_Se_draws,
                                                                  reference_Sp_draws,
                                                                  raw_cutoffs,
                                                                  native_prior_bank,
                                                                  objective,
                                                                  target_grain_size,
                                                                  estimate_beta_baseline_mean,
                                                                  estimate_beta_baseline_SD,
                                                                  estimate_raw_scale_baseline_mean,
                                                                  estimate_raw_scale_baseline_SD,
                                                                  estimate_dirichlet_alpha,
                                                                  fixed_beta_baseline_mean,
                                                                  fixed_beta_baseline_SD,
                                                                  fixed_raw_scale_baseline_mean,
                                                                  fixed_raw_scale_baseline_SD,
                                                                  fixed_dirichlet_alpha) {

    if (!identical(objective, "marginal_composite")) {
        stop("Select objective = marginal_composite explicitly; no joint curve likelihood is implemented.")
    }
    if (!is.numeric(target_grain_size) || length(target_grain_size) != 1 || !is.finite(target_grain_size)
        || target_grain_size < 1 || target_grain_size != as.integer(target_grain_size)) {
        stop("Supply target_grain_size explicitly in target accuracy points per minimum task.")
    }

    required_bank <- c("n_thr", "raw_cutoff_grid", "n_bank", "cutpoint_bank", "phi_bank", "raw_scale_bank",
                       "bank_raw_scale_mean", "bank_raw_scale_SD", "bank_dirichlet_alpha", "softplus",
                       "use_probit_link", "baseline_case", "n_covariates", "n_studies", "scope")
    if (!is.list(native_prior_bank) || !all(required_bank %in% names(native_prior_bank))) {
        stop("Prepare the bank with fn_ord_HSROC_native_prior_integration_bank first.")
    }

    if (!is.matrix(reference_Se_draws) || !is.numeric(reference_Se_draws)
        || !is.matrix(reference_Sp_draws) || !is.numeric(reference_Sp_draws)
        || !identical(dim(reference_Se_draws), dim(reference_Sp_draws))
        || nrow(reference_Se_draws) < 2 || ncol(reference_Se_draws) < 1) {
        stop("Supply complete paired Se and Sp matrices, with draws in rows and raw cutoffs in columns.")
    }
    probabilities <- c(as.vector(reference_Sp_draws), as.vector(reference_Se_draws))
    if (any(!is.finite(probabilities)) || any(probabilities <= 0 | probabilities >= 1)) {
        stop("Every supplied accuracy must be strictly inside zero and one; none are clipped or removed.")
    }
    if (!is.numeric(raw_cutoffs) || !is.null(dim(raw_cutoffs))
        || length(raw_cutoffs) != ncol(reference_Se_draws)
        || any(!is.finite(raw_cutoffs)) || any(diff(raw_cutoffs) <= 0)) {
        stop("Supply the strictly increasing raw scores of the reference matrix columns.")
    }
    cutoff_index <- match(raw_cutoffs, native_prior_bank$raw_cutoff_grid)
    if (anyNA(cutoff_index)) stop("Each target cutoff must occur exactly in the complete native raw-cutoff grid.")

    if (length(raw_cutoffs) > 1) {
        following <- seq.int(2, length(raw_cutoffs))
        preceding <- following - 1
        if (any(reference_Se_draws[, following, drop = FALSE]
                >= reference_Se_draws[, preceding, drop = FALSE])
            || any(reference_Sp_draws[, following, drop = FALSE]
                   <= reference_Sp_draws[, preceding, drop = FALSE])) {
            stop("Each reference curve must have decreasing Se and increasing Sp across its raw cutoffs.")
        }
    }

    masks <- list(estimate_beta_baseline_mean = estimate_beta_baseline_mean,
                  estimate_beta_baseline_SD = estimate_beta_baseline_SD,
                  estimate_raw_scale_baseline_mean = estimate_raw_scale_baseline_mean,
                  estimate_raw_scale_baseline_SD = estimate_raw_scale_baseline_SD,
                  estimate_dirichlet_alpha = estimate_dirichlet_alpha)
    fixed <- list(fixed_beta_baseline_mean = fixed_beta_baseline_mean,
                  fixed_beta_baseline_SD = fixed_beta_baseline_SD,
                  fixed_raw_scale_baseline_mean = fixed_raw_scale_baseline_mean,
                  fixed_raw_scale_baseline_SD = fixed_raw_scale_baseline_SD,
                  fixed_dirichlet_alpha = fixed_dirichlet_alpha)
    initial <- vector(mode = "list", length = length(masks))
    names(initial) <- sub(pattern = "^estimate_", replacement = "estimated_", x = names(masks))

    for (i in seq_along(masks)) {
        size <- if (i == length(masks)) native_prior_bank$n_thr + 1 else 1
        mask <- masks[[i]]
        value <- fixed[[i]]
        if (!is.logical(mask) || !is.null(dim(mask)) || length(mask) != size || anyNA(mask)) {
            stop(paste0(names(masks)[i], " must explicitly select every applicable entry with TRUE or FALSE."))
        }
        if (!is.numeric(value) || !is.null(dim(value)) || length(value) != size
            || any(!is.finite(value)) || (i %in% c(2, 4, 5) && any(value <= 0))) {
            stop(paste0(names(fixed)[i], " requires finite values and strictly positive SDs/concentrations."))
        }
        masks[[i]] <- as.integer(mask)
        initial[[i]] <- unname(value[mask])
    }

    n_reference <- nrow(reference_Se_draws)
    n_selected <- ncol(reference_Se_draws)
    status <- rep(c(1, 2), each = n_reference * n_selected)
    cutoff <- rep(rep(cutoff_index, each = n_reference), times = 2)
    metadata <- c("baseline_case", "n_covariates", "n_studies", "scope")
    stan_data <- native_prior_bank[setdiff(required_bank, metadata)]
    stan_data$n_target <- length(probabilities)
    stan_data$target_status <- status
    stan_data$target_cutoff_index <- cutoff
    stan_data$target_probability <- probabilities
    stan_data$target_grain_size <- as.integer(target_grain_size)
    stan_data$marginal_composite_likelihood <- 1
    stan_data <- c(stan_data, masks, fixed)

    message(paste0("ord_HSROC target: ", length(probabilities),
                   " marginal observations; location and raw-scale prior inputs are shared by both statuses."))

    list(stan_data = stan_data,
         initial = initial,
         target_map = data.frame(draw_row = rep(seq_len(n_reference), times = 2 * n_selected),
                                 quantity = ifelse(status == 1, "Sp", "Se"),
                                 raw_cutoff = native_prior_bank$raw_cutoff_grid[cutoff],
                                 native_cutoff_index = cutoff),
         baseline_case = native_prior_bank$baseline_case,
         n_covariates = native_prior_bank$n_covariates,
         n_studies = native_prior_bank$n_studies,
         native_stan_data = native_prior_bank$native_stan_data,
         objective = objective,
         scope = native_prior_bank$scope,
         approximation = "self-normalised importance integration over the paired native R&G prior bank")

}






















