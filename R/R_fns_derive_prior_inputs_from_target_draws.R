##
## ---- Derive inputs for an explicitly selected observable-distribution likelihood:
##
## The supplied target draws and native integration bank remain unchanged. Initial values and
## optimiser or sampling controls are explicit inputs. These functions never compile a native model.
## An optimised likelihood is not a check of the realised prior in the receiving native model.
##


fn_optimise_a_prior_derivation_likelihood <- function(stan_data,
                                                     initial,
                                                     compiled_model,
                                                     expected_Stan_file,
                                                     optimisation_arguments) {

    if (!is.list(stan_data) || !is.list(initial) || is.null(names(initial))
        || anyNA(names(initial)) || any(!nzchar(names(initial))) || anyDuplicated(names(initial))) {
        stop("Supply the prepared Stan data and uniquely named initial parameter blocks.")
    }

    if (is.null(compiled_model) || !is.function(compiled_model$optimize)
        || !is.function(compiled_model$stan_file)
        || !identical(basename(compiled_model$stan_file()), expected_Stan_file)) {
        stop(paste0("Supply the compiled ", expected_Stan_file, " model, not a fitted native model."))
    }

    if (!is.list(optimisation_arguments) || is.null(names(optimisation_arguments))
        || anyNA(names(optimisation_arguments)) || any(!nzchar(names(optimisation_arguments)))
        || anyDuplicated(names(optimisation_arguments))
        || any(c("data", "init") %in% names(optimisation_arguments))) {
        stop("Supply uniquely named optimiser arguments; data and initial values are supplied separately.")
    }

    required <- c("seed", "algorithm", "iter", "jacobian", "refresh", "show_messages")

    if (!all(required %in% names(optimisation_arguments))
        || !identical(optimisation_arguments$jacobian, FALSE)
        || !identical(optimisation_arguments$show_messages, TRUE)
        || !is.numeric(optimisation_arguments$refresh) || length(optimisation_arguments$refresh) != 1
        || !is.finite(optimisation_arguments$refresh) || optimisation_arguments$refresh <= 0) {
        stop("Supply explicit optimiser controls, with jacobian FALSE, show_messages TRUE and positive refresh.")
    }

    message(paste0("Observable-likelihood model: ", compiled_model$stan_file(), "."))
    message("Initial parameter values supplied by the caller:")
    print(initial)
    message("Explicit optimiser controls; the parameter-transform Jacobian is excluded:")
    print(optimisation_arguments)
    progress <- utils::txtProgressBar(min = 0, max = 1, initial = 0, style = 3)
    on.exit(close(progress), add = TRUE)

    fit <- do.call(what = compiled_model$optimize,
                   args = c(list(data = stan_data, init = list(initial)), optimisation_arguments))
    return_codes <- fit$return_codes()

    if (!is.numeric(return_codes) || length(return_codes) == 0
        || any(!is.finite(return_codes)) || any(return_codes != 0)) {
        stop("The optimisation did not finish successfully; no derived inputs are cleared for use.")
    }

    utils::setTxtProgressBar(pb = progress, value = 1)
    message(paste0("The optimiser process exited successfully. Its termination output must also confirm ",
                   "convergence; a successful exit alone does not exclude the iteration limit."))
    fit

}


fn_read_prior_derivation_values <- function(fit,
                                           variable,
                                           expected_names) {

    values <- fit$summary(variables = variable)

    if (!is.data.frame(values) || !all(c("variable", "estimate") %in% names(values))
        || anyDuplicated(values$variable) || !setequal(values$variable, expected_names)) {
        stop(paste0("The optimisation output does not identify every ", variable, " entry exactly."))
    }

    result <- values$estimate[match(expected_names, values$variable)]

    if (!is.numeric(result) || any(!is.finite(result))) {
        stop(paste0("The optimisation output contains nonfinite ", variable, " entries."))
    }

    stats::setNames(result, expected_names)

}


fn_derive_the_cutpoint_prior <- function(request,
                                        initial_alpha,
                                        compiled_model,
                                        optimisation_arguments) {

    if (!is.list(request) || !is.list(request$stan_data)
        || !all(c("n_draws", "n_thr", "cutpoint_draws", "anchor", "use_probit_link")
                %in% names(request$stan_data))) {
        stop("Prepare the target with fn_data_for_deriving_the_cutpoint_prior first.")
    }

    n_alpha <- request$stan_data$n_thr + 1

    if (!is.numeric(initial_alpha) || !is.null(dim(initial_alpha)) || length(initial_alpha) != n_alpha
        || any(!is.finite(initial_alpha)) || any(initial_alpha <= 0)) {
        stop("Supply one finite positive initial alpha for each interval separated by the target cutpoints.")
    }

    message(request$scope)
    message(request$alpha_scope)
    message("All supplied cutpoint concentrations are estimated; no random location or scale is integrated.")

    fit <- fn_optimise_a_prior_derivation_likelihood(
                stan_data = request$stan_data,
                initial = list(alpha = unname(initial_alpha)),
                compiled_model = compiled_model,
                expected_Stan_file = "derive_cutpoint_prior_target_draws.stan",
                optimisation_arguments = optimisation_arguments)

    alpha <- fn_read_prior_derivation_values(fit = fit,
                                             variable = "alpha",
                                             expected_names = paste0("alpha[", seq_len(n_alpha), "]"))

    if (any(alpha <= 0)) stop("The returned concentrations must remain strictly positive.")

    print(alpha)
    message("These are cutpoint-law inputs, not cleared priors for a model with random location or scale.")

    list(alpha = alpha,
         optimisation_fit = fit,
         raw_cutoffs = request$raw_cutoffs,
         accuracy_measure = request$accuracy_measure,
         alpha_scope = request$alpha_scope,
         scope = request$scope,
         optimisation_converged = NA,
         native_model_check_passed = NA)

}


##
## ---- Sample the selected likelihood with a flat constrained-input measure:
##


fn_sample_a_prior_derivation_likelihood <- function(stan_data,
                                                    initial_specification,
                                                    positive_blocks,
                                                    initial_values_per_chain,
                                                    compiled_model,
                                                    expected_Stan_file,
                                                    sampling_arguments,
                                                    probs) {

    if (!is.list(stan_data) || !is.list(initial_specification)
        || is.null(names(initial_specification)) || anyDuplicated(names(initial_specification))) {
        stop("Supply prepared Stan data and the complete named parameter-block specification.")
    }

    if (is.null(compiled_model) || !is.function(compiled_model$sample)
        || !is.function(compiled_model$stan_file)
        || !identical(basename(compiled_model$stan_file()), expected_Stan_file)) {
        stop(paste0("Supply the compiled ", expected_Stan_file, " model."))
    }

    if (!is.list(sampling_arguments) || is.null(names(sampling_arguments))
        || anyNA(names(sampling_arguments)) || any(!nzchar(names(sampling_arguments)))
        || anyDuplicated(names(sampling_arguments))
        || any(c("data", "init", "jacobian") %in% names(sampling_arguments))) {
        stop(paste0("Supply unique named sampling arguments; data/init are separate and ",
                     "jacobian is not a sample option."))
    }

    available <- names(formals(compiled_model$sample))
    if (any(!names(sampling_arguments) %in% available)) {
        stop("A supplied argument is not an explicit argument of this compiled model's sample method.")
    }

    required <- c("seed", "chains", "parallel_chains", "threads_per_chain", "iter_warmup", "iter_sampling",
                  "adapt_delta", "metric", "max_treedepth", "refresh", "show_messages", "fixed_param",
                  "output_dir")
    if (!all(required %in% names(sampling_arguments))) {
        stop(paste0("Supply every applicable sampling control explicitly: ",
                     paste(required, collapse = ", "), "."))
    }

    integer_controls <- c("seed", "chains", "parallel_chains", "threads_per_chain", "iter_warmup",
                          "iter_sampling", "max_treedepth", "refresh")
    for (name in integer_controls) {

        value <- sampling_arguments[[name]]
        if (!is.numeric(value) || length(value) != 1 || !is.finite(value)
            || value < 0 || (name != "seed" && value == 0) || value != floor(value)) {
            stop(paste0(name, " must be a supplied positive integer; seed may also be zero."))
        }

    }

    if (!identical(sampling_arguments$fixed_param, FALSE)
        || !identical(sampling_arguments$show_messages, TRUE)
        || !is.numeric(sampling_arguments$adapt_delta) || length(sampling_arguments$adapt_delta) != 1
        || !is.finite(sampling_arguments$adapt_delta)
        || sampling_arguments$adapt_delta <= 0 || sampling_arguments$adapt_delta >= 1
        || !is.character(sampling_arguments$metric) || length(sampling_arguments$metric) != 1
        || !sampling_arguments$metric %in% c("unit_e", "diag_e", "dense_e")) {
        stop(paste0("Use HMC, live messages, an explicit metric and an explicit adapt_delta ",
                     "strictly inside zero and one."))
    }

    if (!is.character(sampling_arguments$output_dir) || length(sampling_arguments$output_dir) != 1
        || is.na(sampling_arguments$output_dir) || !dir.exists(sampling_arguments$output_dir)) {
        stop("Supply an existing explicit output directory to retain every sampling result.")
    }

    if (!is.numeric(probs) || !is.null(dim(probs)) || length(probs) < 1
        || any(!is.finite(probs)) || any(probs <= 0 | probs >= 1) || any(diff(probs) <= 0)) {
        stop("Supply the strictly increasing interior probabilities of the requested posterior intervals.")
    }

    if (!is.list(initial_values_per_chain)
        || length(initial_values_per_chain) != sampling_arguments$chains) {
        stop("Supply a complete initial-value list for each explicitly requested chain; no random starts.")
    }

    for (initial in initial_values_per_chain) {

        if (!is.list(initial) || !identical(names(initial), names(initial_specification))) {
            stop("Each chain must name every parameter block in the declared order, including empty blocks.")
        }
        for (block in names(initial_specification)) {

            value <- initial[[block]]
            if (!is.numeric(value) || !is.null(dim(value))
                || length(value) != length(initial_specification[[block]]) || any(!is.finite(value))
                || (block %in% positive_blocks && any(value <= 0))) {
                stop(paste0("Supply complete finite initial values on the declared support for ", block, "."))
            }

        }

    }

    active <- names(initial_specification)[lengths(initial_specification) > 0]
    if (length(active) == 0) stop("All inputs are fixed; return them without invoking the sampler.")

    message(paste0("Sampling derivation model: ", compiled_model$stan_file(), "."))
    message(paste0("The free baseline location inputs have a flat measure on their constrained space; any free ",
                   "raw-scale or Box-Cox input carries the hyperprior given in the Stan data."))
    message("HMC includes its parameter-transform Jacobian; posterior propriety is not established here.")
    message("Explicit sampling arguments and complete per-chain initial values:")
    print(sampling_arguments)
    print(initial_values_per_chain)
    message("Sampling stage progress; iteration progress is printed by the sampler:")
    progress <- utils::txtProgressBar(min = 0, max = 1, initial = 0, style = 3)
    on.exit(close(progress), add = TRUE)
    fit <- do.call(what = compiled_model$sample,
                   args = c(list(data = stan_data, init = initial_values_per_chain), sampling_arguments))
    codes <- fit$return_codes()

    if (!is.numeric(codes) || length(codes) != sampling_arguments$chains
        || any(!is.finite(codes)) || any(codes != 0)) {
        stop("At least one sampling process did not finish successfully; retained outputs require inspection.")
    }

    draws <- fit$draws(variables = active, inc_warmup = FALSE, format = "draws_array")
    draws <- unclass(draws)
    ## (NicoStan's summary: its rank-normalized split ESS (bulk, tail, theta-squared) and R-hat)
    summary_of_NicoStan <- generate_summary_tibble(n_threads = NULL,
                                                   input_array = draws,
                                                   iter_dim = 1,
                                                   chain_dim = 2,
                                                   param_names = dimnames(draws)[[3]],
                                                   compute_nested_rhat = FALSE,
                                                   n_superchains = NULL)
    summary <- tibble::tibble(variable = summary_of_NicoStan$parameter,
                              median = summary_of_NicoStan$`50%`,
                              rhat = summary_of_NicoStan$Rhat,
                              ess_bulk = summary_of_NicoStan$n_eff,
                              ess_tail = summary_of_NicoStan$n_eff_tail,
                              ess_sd = summary_of_NicoStan$n_eff_sd)
    summary$ess_minimum <- pmin(summary$ess_bulk, summary$ess_tail, summary$ess_sd)
    quantiles <- matrix(data = NA_real_, nrow = nrow(summary), ncol = length(probs),
                         dimnames = list(summary$variable, as.character(probs)))

    for (i in seq_len(nrow(summary))) {

        values <- as.vector(draws[, , summary$variable[i], drop = FALSE])
        if (any(!is.finite(values))) stop("Nonfinite parameter draws are retained but cannot be summarised.")
        quantiles[i, ] <- stats::quantile(values, probs = probs, names = FALSE)

    }

    diagnostics <- fit$sampler_diagnostics(inc_warmup = FALSE, format = "draws_array")
    if (!"divergent__" %in% dimnames(diagnostics)[[3]]) {
        stop("The actual post-warm-up divergence indicators are unavailable; do not report zero.")
    }
    divergent <- diagnostics[, , "divergent__", drop = FALSE]
    if (anyNA(divergent) || any(!divergent %in% c(0, 1))) stop("Invalid divergence indicators.")
    denominator <- length(divergent)
    if (denominator == 0) stop("No post-warm-up transitions are available for a divergence percentage.")
    divergence_percent <- 100 * sum(divergent) / denominator

    print(summary)
    print(quantiles)
    message(paste0("Divergent retained post-warm-up transitions: ", sum(divergent), " / ", denominator,
                   " (", divergence_percent, "%)."))
    message("Process completion does not establish posterior propriety, convergence or a native-prior match.")
    utils::setTxtProgressBar(pb = progress, value = 1)

    list(sampling_fit = fit,
         posterior_summary = summary,
         posterior_quantiles = quantiles,
         probs = probs,
         divergent_transitions = sum(divergent),
         post_warmup_transitions = denominator,
         divergence_percent = divergence_percent,
         divergence_denominator_scope = "retained post-warm-up transitions; thinning may omit transitions",
         sampling_arguments = sampling_arguments,
         actual_sampling_metadata = fit$metadata(),
         posterior_propriety_established = NA,
         sampling_converged = NA,
         native_model_check_passed = NA)

}


fn_derive_the_cutpoint_prior_by_sampling <- function(request,
                                                    initial_values_per_chain,
                                                    compiled_model,
                                                    sampling_arguments,
                                                    probs) {

    if (!is.list(request) || !is.list(request$stan_data)
        || !all(c("n_draws", "n_thr", "cutpoint_draws", "anchor", "use_probit_link")
                %in% names(request$stan_data))) {
        stop("Prepare the target with fn_data_for_deriving_the_cutpoint_prior first.")
    }

    n_alpha <- request$stan_data$n_thr + 1
    specification <- list(alpha = numeric(length = n_alpha))
    result <- fn_sample_a_prior_derivation_likelihood(
                  stan_data = request$stan_data,
                  initial_specification = specification,
                  positive_blocks = "alpha",
                  initial_values_per_chain = initial_values_per_chain,
                  compiled_model = compiled_model,
                  expected_Stan_file = "derive_cutpoint_prior_target_draws.stan",
                  sampling_arguments = sampling_arguments,
                  probs = probs)
    expected <- paste0("alpha[", seq_len(n_alpha), "]")
    index <- match(expected, result$posterior_summary$variable)
    if (anyNA(index)) stop("The sampling output does not identify every concentration exactly.")

    c(list(alpha = stats::setNames(result$posterior_summary$median[index], expected),
           raw_cutoffs = request$raw_cutoffs,
           accuracy_measure = request$accuracy_measure,
           alpha_scope = request$alpha_scope,
           scope = request$scope), result)

}


fn_check_the_beta_only_prior_measure <- function(stan_data,
                                                 model_kind) {

    if (!is.list(stan_data) || !is.character(model_kind) || length(model_kind) != 1
        || !model_kind %in% c("Jones", "R_and_G")) {
        stop("Supply the actual marginal-likelihood data and explicit receiving-model kind.")
    }

    data <- stan_data
    for (field in c("n_target", "n_bank", "n_thr")) {

        value <- data[[field]]
        if (!is.numeric(value) || length(value) != 1 || !is.finite(value)
            || value < 1 || value != floor(value)) {
            stop(paste0(field, " must be the actual positive integer of the prepared request."))
        }

    }
    for (field in c("softplus", "use_probit_link")) {

        value <- data[[field]]
        if (!is.numeric(value) || length(value) != 1 || !is.finite(value) || !value %in% c(0, 1)) {
            stop(paste0(field, " must be the actual native zero-or-one flag."))
        }

    }

    is_Jones <- model_kind == "Jones"
    n_blocks <- if (is_Jones) 2 else 1
    fields <- c("beta_baseline_mean", "beta_baseline_SD", "raw_scale_baseline_mean", "raw_scale_baseline_SD",
                if (is_Jones) "boxcox_prior" else "dirichlet_alpha")

    for (field in fields) {

        size <- if (field == "dirichlet_alpha") data$n_thr + 1 else n_blocks
        mask <- data[[paste0("estimate_", field)]]
        if (!is.numeric(mask) || !is.null(dim(mask)) || length(mask) != size
            || anyNA(mask) || any(!mask %in% c(0, 1))) {
            stop("Retain complete fixed/free masks for the actual finite-bank likelihood.")
        }
        if (!startsWith(field, "beta_") && any(mask == 1)) {
            stop("The beta-only propriety check requires all importance-weight inputs to remain fixed.")
        }

    }

    for (field in c("fixed_beta_baseline_mean", "fixed_beta_baseline_SD")) {

        value <- data[[field]]
        if (!is.numeric(value) || !is.null(dim(value)) || length(value) != n_blocks
            || any(!is.finite(value)) || (endsWith(field, "_SD") && any(value <= 0))) {
            stop("Supply finite baseline-normal inputs, with positive standard deviations.")
        }

    }

    for (field in c("target_status", "target_cutoff_index", "target_probability")) {

        value <- data[[field]]
        if (!is.numeric(value) || !is.null(dim(value)) || length(value) != data$n_target
            || any(!is.finite(value))) {
            stop("Every target must retain its finite probability, disease status and native cutoff index.")
        }

    }
    if (any(!data$target_status %in% c(1, 2)) || !all(c(1, 2) %in% data$target_status)
        || any(data$target_cutoff_index != floor(data$target_cutoff_index))
        || any(data$target_cutoff_index < 1 | data$target_cutoff_index > data$n_thr)
        || any(data$target_probability <= 0 | data$target_probability >= 1)) {
        stop("Retain both statuses, valid cutoff indices and strictly interior target probabilities.")
    }

    C <- data$cutpoint_bank
    raw <- data$raw_scale_bank
    expected_C <- if (is_Jones) c(2, data$n_bank, data$n_thr) else c(data$n_bank, data$n_thr)
    if (!is.numeric(C) || !identical(dim(C), as.integer(expected_C)) || any(!is.finite(C))) {
        stop("The native cutpoint-bank dimensions and finite entries must be retained exactly.")
    }
    if (is_Jones) {

        if (!is.numeric(raw) || !identical(dim(raw), as.integer(c(data$n_bank, 2)))
            || any(!is.finite(raw))) stop("The Jones raw-scale bank must be a finite bank-row by status matrix.")
        signed_raw <- raw

    } else {

        if (!is.numeric(raw) || !is.null(dim(raw)) || length(raw) != data$n_bank
            || any(!is.finite(raw))) stop("The ordinal-HSROC raw-scale bank must retain every finite native row.")
        signed_raw <- cbind(raw, -raw)

    }

    scale <- if (data$softplus == 1) {
        1.4426950408889634 * (pmax(signed_raw, 0) + log1p(exp(-abs(signed_raw))))
    } else exp(signed_raw)
    if (any(!is.finite(scale)) || any(scale <= 0)) stop("Every native scale must be finite and positive.")
    if (!is_Jones && any(!is.finite(1 / scale))) stop("Every ordinal-HSROC inverse scale must be finite.")

    z <- if (data$use_probit_link == 1) {
        stats::qnorm(p = data$target_probability, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
    } else {
        stats::qlogis(p = data$target_probability, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE)
    }
    native_argument <- ifelse(data$target_status == 1, -z, z)
    if (any(!is.finite(z)) || (data$use_probit_link == 1
        && any(native_argument < -37.5 | native_argument > 8.25))) {
        stop("Targets must remain inside the native interior-link branch without clipping.")
    }

    message("Checking the flat baseline-normal measure of the fixed finite-bank composite likelihood:")
    progress <- utils::txtProgressBar(min = 0, max = n_blocks, initial = 0, style = 3)
    on.exit(close(progress), add = TRUE)
    checks <- vector(mode = "list", length = n_blocks)

    for (block in seq_len(n_blocks)) {

        indices <- if (is_Jones) which(data$target_status == block) else seq_len(data$n_target)
        n <- length(indices)
        mean_free <- data$estimate_beta_baseline_mean[block] == 1
        SD_free <- data$estimate_beta_baseline_SD[block] == 1
        if (n < 1 || (SD_free && n <= 1 + as.integer(mean_free))) {
            stop("The flat baseline-normal measure has a nonintegrable large-SD tail or no target information.")
        }
        common_centres <- NULL
        all_sets_contain_fixed_mean <- TRUE

        if (SD_free) {

            for (position in seq_along(indices)) {

                i <- indices[position]
                status <- data$target_status[i]
                k <- data$target_cutoff_index[i]
                centres <- if (is_Jones) {
                    C[status, , k] + (if (status == 1) -1 else 1) * z[i] * scale[, status]
                } else {
                    (if (status == 1) C[, k] else -C[, k]) - z[i] / scale[, status]
                }
                if (any(!is.finite(centres))) stop("Every transformed native-bank centre must remain finite.")

                if (mean_free) {
                    common_centres <- if (position == 1) unique(centres) else intersect(common_centres, centres)
                } else {
                    all_sets_contain_fixed_mean <- all_sets_contain_fixed_mean &&
                        any(centres == data$fixed_beta_baseline_mean[block])
                }

            }

            zero_residual_assignment <- if (mean_free) length(common_centres) > 0 else all_sets_contain_fixed_mean
            if (zero_residual_assignment) {
                stop(paste0("The flat baseline-normal measure has a nonintegrable zero-SD limit: ",
                             "a positive-weight component assignment has zero residual sum of squares. ",
                             "No numerical tolerance or extra prior is selected."))
            }

        }

        checks[[block]] <- data.frame(baseline_block = block,
                                      n_target = n,
                                      mean_free = mean_free,
                                      SD_free = SD_free,
                                      finite_bank_measure_proper = TRUE)
        utils::setTxtProgressBar(pb = progress, value = block)

    }

    table <- do.call(what = rbind, args = checks)
    print(table)
    message(paste0("This uses exact equality of computed centres; ",
                   "near-degeneracy and numerical conditioning are not checked."))
    message(paste0("Propriety of this fixed finite-bank likelihood ",
                   "does not qualify its native approximation or HMC output."))

    list(finite_bank_measure_proper = TRUE,
         by_baseline_block = table,
         exact_equality_of_computed_centres = TRUE,
         near_degeneracy_checked = FALSE,
         scope = "flat baseline-normal inputs; fixed finite positive mixture weights and Jacobians only",
         native_model_check_passed = NA,
         sampling_converged = NA)

}


fn_derive_a_marginal_accuracy_prior_by_sampling <- function(request,
                                                          initial_values_per_chain,
                                                          compiled_model,
                                                          sampling_arguments,
                                                          probs,
                                                          model_kind) {

    if (!is.list(request) || !identical(request$objective, "marginal_composite")
        || !is.list(request$stan_data) || !model_kind %in% c("Jones", "R_and_G")) {
        stop("Prepare the explicitly selected native marginal-composite request first.")
    }

    data <- request$stan_data
    fields <- c("beta_baseline_mean", "beta_baseline_SD", "raw_scale_baseline_mean", "raw_scale_baseline_SD")
    fields <- c(fields, if (model_kind == "Jones") c("boxcox_mean", "boxcox_SD") else "dirichlet_alpha")
    supplied <- specification <- masks <- vector(mode = "list", length = length(fields))
    names(supplied) <- paste0("prior_", fields)
    names(specification) <- paste0("estimated_", fields)
    names(masks) <- fields

    for (field in fields) {

        is_boxcox <- startsWith(field, "boxcox_")
        size <- if (model_kind == "Jones") 2 else if (field == "dirichlet_alpha") {
            data$n_thr + 1
        } else 1
        mask <- data[[if (is_boxcox) "estimate_boxcox_prior" else paste0("estimate_", field)]]
        value <- data[[paste0("fixed_", field)]]
        positive <- endsWith(field, "_SD") || field == "dirichlet_alpha"

        if (!is.numeric(mask) || !is.null(dim(mask)) || length(mask) != size
            || anyNA(mask) || any(!mask %in% c(0, 1))
            || !is.numeric(value) || !is.null(dim(value)) || length(value) != size
            || any(!is.finite(value)) || (positive && any(value <= 0))) {
            stop(paste0("The request lacks the explicit fixed/free inputs for ", field, "."))
        }

        supplied[[paste0("prior_", field)]] <- value
        specification[[paste0("estimated_", field)]] <- unname(value[mask == 1])
        masks[[field]] <- mask

    }

    importance_fields <- setdiff(fields, c("beta_baseline_mean", "beta_baseline_SD"))
    any_importance_input_free <- any(vapply(X = masks[importance_fields],
                                            FUN = function(mask) {

                                                any(mask == 1)

                                            },
                                            FUN.VALUE = logical(length = 1)))
    ## (a free raw-scale, Box-Cox or concentration input has an improper posterior under a flat measure, so it
    ##  is allowed only under the proper hyperpriors the data carries: a normal on each mean, a half-normal on
    ##  each SD and concentration, every value supplied by the user)
    hyperprior_fields <- if (model_kind == "Jones") {
        c("hyperprior_raw_scale_baseline_mean_mean", "hyperprior_raw_scale_baseline_mean_SD",
          "hyperprior_raw_scale_baseline_SD_scale", "hyperprior_boxcox_mean_mean", "hyperprior_boxcox_mean_SD",
          "hyperprior_boxcox_SD_scale")
    } else {
        c("hyperprior_raw_scale_baseline_mean_mean", "hyperprior_raw_scale_baseline_mean_SD",
          "hyperprior_raw_scale_baseline_SD_scale", "hyperprior_dirichlet_alpha_scale")
    }
    hyperpriors_supplied <- all(hyperprior_fields %in% names(data)) &&
                            all(vapply(X = data[hyperprior_fields],
                                       FUN = function(value) {

                                           is.numeric(value) && length(value) == 1 && is.finite(value)

                                       },
                                       FUN.VALUE = logical(length = 1))) &&
                            all(unlist(data[grep("_SD$|_scale$", hyperprior_fields, value = TRUE)]) > 0)
    if (any_importance_input_free && !hyperpriors_supplied) {
        stop(paste0("A flat unbounded measure on free native-bank importance inputs gives an improper ",
                     "posterior: the finite-bank normalised mixture has positive limiting likelihood. ",
                     "Free raw-scale, Box-Cox and concentration inputs need the proper hyperpriors in the data: ",
                     paste(hyperprior_fields, collapse = ", "), " (finite; SDs and scales positive)."))
    }
    if (any_importance_input_free) {
        message(paste0("Free importance inputs under the supplied hyperpriors: ",
                       paste(paste0(hyperprior_fields, " = ", unlist(data[hyperprior_fields])), collapse = ", ")))
    }

    if (!is.numeric(data$n_target) || length(data$n_target) != 1 || !is.finite(data$n_target)
        || data$n_target < 1 || data$n_target != floor(data$n_target)
        || !is.numeric(data$target_status) || !is.null(dim(data$target_status))
        || length(data$target_status) != data$n_target || anyNA(data$target_status)
        || any(!data$target_status %in% c(1, 2))) {
        stop("Retain the actual target count and complete status map from the prepared marginal request.")
    }

    for (status in seq_len(if (model_kind == "Jones") 2 else 1)) {

        n_observations <- if (model_kind == "Jones") sum(data$target_status == status) else data$n_target
        mean_is_free <- masks$beta_baseline_mean[status] == 1
        SD_is_free <- masks$beta_baseline_SD[status] == 1
        if (SD_is_free && n_observations <= 1 + as.integer(mean_is_free)) {
            stop(paste0("The flat baseline-normal input measure has a nonintegrable large-SD tail: ",
                         n_observations, " target observations for baseline block ", status, ". ",
                         "A free SD needs more than one observation; freeing its mean too needs more than two. ",
                         "This necessary condition is mathematical, not a posterior qualification threshold."))
        }

    }

    information <- list(objective = request$objective,
                         scope = request$scope,
                         approximation = request$approximation,
                         baseline_case = request$baseline_case,
                         n_covariates = request$n_covariates,
                         n_studies = request$n_studies)
    message("This is a marginal composite likelihood, not a joint accuracy-curve likelihood.")
    message(request$scope)
    message(request$approximation)

    if (sum(lengths(specification)) == 0) {

        message("Every prior input is fixed; the exact supplied inputs are returned without sampling.")
        return(c(list(prior_values = supplied,
                      stan_data = request$stan_data,
                      sampling_fit = NULL,
                      posterior_summary = NULL,
                      posterior_quantiles = NULL,
                      native_model_check_passed = NA), information))

    }

    ## (the beta-only propriety preflight applies when every importance input is fixed; with a free one the
    ##  proper hyperpriors above make the posterior proper)
    measure_check <- if (any_importance_input_free) {
        list(proper = TRUE, note = "free importance inputs under the supplied proper hyperpriors")
    } else fn_check_the_beta_only_prior_measure(stan_data = data, model_kind = model_kind)
    result <- fn_sample_a_prior_derivation_likelihood(
                  stan_data = data,
                  initial_specification = specification,
                  positive_blocks = names(specification)[endsWith(names(specification), "_SD")
                                                         | names(specification) == "estimated_dirichlet_alpha"],
                  initial_values_per_chain = initial_values_per_chain,
                  compiled_model = compiled_model,
                  expected_Stan_file = if (model_kind == "Jones") {
                      "derive_Jones_marginal_accuracy_prior_target_draws.stan"
                  } else {
                      paste0("derive_the_", model_kind,
                             "_marginal_accuracy_prior_from_target_draws.stan")
                  },
                  sampling_arguments = sampling_arguments,
                  probs = probs)
    result$posterior_propriety_established <- measure_check$finite_bank_measure_proper
    result$posterior_propriety_scope <- measure_check$scope
    entry_map <- list()

    for (field in fields) {

        selected <- which(masks[[field]] == 1)
        if (length(selected) == 0) next
        variable <- paste0("estimated_", field, "[", seq_along(selected), "]")
        index <- match(variable, result$posterior_summary$variable)
        if (anyNA(index)) stop("The sampling output does not identify every selected input exactly.")
        supplied[[paste0("prior_", field)]][selected] <- result$posterior_summary$median[index]
        entry_map[[field]] <- data.frame(variable = variable, input_entry = selected,
                                         row.names = NULL, stringsAsFactors = FALSE)

    }

    print(supplied)
    message("Fixed inputs are unchanged; posterior medians replace only the selected free entries.")
    message("Baseline projection inputs require a separate native likelihood-off accuracy check before use.")

    c(list(prior_values = supplied,
           stan_data = request$stan_data,
           prior_measure_check = measure_check,
           free_input_map = entry_map), information, result)

}


fn_derive_the_Jones_marginal_accuracy_prior_by_sampling <- function(request,
                                                                   initial_values_per_chain,
                                                                   compiled_model,
                                                                   sampling_arguments,
                                                                   probs) {

    fn_derive_a_marginal_accuracy_prior_by_sampling(request = request,
                                                  initial_values_per_chain = initial_values_per_chain,
                                                  compiled_model = compiled_model,
                                                  sampling_arguments = sampling_arguments,
                                                  probs = probs,
                                                  model_kind = "Jones")

}


fn_derive_the_R_and_G_marginal_accuracy_prior_by_sampling <- function(request,
                                                                     initial_values_per_chain,
                                                                     compiled_model,
                                                                     sampling_arguments,
                                                                     probs) {

    fn_derive_a_marginal_accuracy_prior_by_sampling(request = request,
                                                  initial_values_per_chain = initial_values_per_chain,
                                                  compiled_model = compiled_model,
                                                  sampling_arguments = sampling_arguments,
                                                  probs = probs,
                                                  model_kind = "R_and_G")

}


fn_derive_the_Jones_marginal_accuracy_prior <- function(request,
                                                       compiled_model,
                                                       optimisation_arguments) {

    if (!is.list(request) || !identical(request$objective, "marginal_composite")
        || !is.list(request$stan_data)) {
        stop("Prepare the explicitly selected Jones marginal-composite request first.")
    }

    data <- request$stan_data
    fields <- c("beta_baseline_mean", "beta_baseline_SD", "raw_scale_baseline_mean",
                "raw_scale_baseline_SD", "boxcox_mean", "boxcox_SD")
    initial <- vector(mode = "list", length = length(fields))
    names(initial) <- paste0("estimated_", fields)
    supplied <- vector(mode = "list", length = length(fields))
    names(supplied) <- paste0("prior_", fields)

    for (field in fields) {

        is_boxcox <- startsWith(field, "boxcox_")
        size <- if (is_boxcox) 1 else 2
        mask <- data[[if (is_boxcox) "estimate_boxcox_prior" else paste0("estimate_", field)]]
        value <- data[[paste0("fixed_", field)]]

        if (!is.numeric(mask) || !is.null(dim(mask)) || length(mask) != size
            || anyNA(mask) || any(!(mask %in% c(0, 1)))
            || !is.numeric(value) || !is.null(dim(value)) || length(value) != size
            || any(!is.finite(value)) || (endsWith(field, "_SD") && any(value <= 0))) {
            stop(paste0("The Jones request lacks the explicit fixed/free inputs for ", field, "."))
        }

        initial[[paste0("estimated_", field)]] <- unname(value[mask == 1])
        supplied[[paste0("prior_", field)]] <- value

    }

    message(request$scope)
    message(request$approximation)
    message("This is a marginal composite likelihood, not a joint accuracy-curve likelihood.")

    if (sum(lengths(initial)) == 0) {

        message("Every prior input is fixed; no optimisation is performed.")
        evaluation <- fn_log_density_and_gradient_of_the_Jones_marginal_accuracy_prior(theta = initial,
                                                                                      stan_data = data)

        return(list(prior_values = supplied,
                    log_density = evaluation$log_density,
                    optimisation_fit = NULL,
                    objective = request$objective,
                    scope = request$scope,
                    approximation = request$approximation,
                    baseline_case = request$baseline_case,
                    n_covariates = request$n_covariates,
                    optimisation_converged = NA,
                    native_model_check_passed = NA))

    }

    fit <- fn_optimise_a_prior_derivation_likelihood(
                stan_data = data,
                initial = initial,
                compiled_model = compiled_model,
                expected_Stan_file = "derive_Jones_marginal_accuracy_prior_target_draws.stan",
                optimisation_arguments = optimisation_arguments)

    derived <- supplied

    for (field in fields) {

        mask <- data[[if (startsWith(field, "boxcox_")) {
            "estimate_boxcox_prior"
        } else paste0("estimate_", field)]]

        if (!any(mask == 1)) next

        variable <- if (startsWith(field, "boxcox_")) {
            sub(pattern = "boxcox_", replacement = "boxcox_lambda_", x = paste0("prior_", field))
        } else paste0("prior_", field)
        expected <- if (startsWith(field, "boxcox_")) variable else paste0(variable, "[", seq_len(2), "]")
        value <- fn_read_prior_derivation_values(fit = fit, variable = variable, expected_names = expected)

        if (endsWith(field, "_SD") && any(value <= 0)) {
            stop(paste0("The returned ", field, " values must be strictly positive."))
        }

        derived[[paste0("prior_", field)]][mask == 1] <- unname(value[mask == 1])

    }

    print(derived)
    message("Baseline projection inputs are not arbitrary regression-coefficient inputs.")
    message("Check the realised means and intervals with the native likelihood-off model before use.")

    list(prior_values = derived,
         optimisation_fit = fit,
         objective = request$objective,
         scope = request$scope,
         approximation = request$approximation,
         baseline_case = request$baseline_case,
         n_covariates = request$n_covariates,
         optimisation_converged = NA,
         native_model_check_passed = NA)

}






















