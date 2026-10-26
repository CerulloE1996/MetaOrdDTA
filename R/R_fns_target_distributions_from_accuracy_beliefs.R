##
## ---- Construct explicitly selected marginal target distributions from accuracy beliefs:
##
## These distributions describe the elicited targets, not the realised priors of a native model.
## The supplied constructors retain their own calculations. Independent marginal draws do not
## define a joint curve across cutoffs and must not be sorted into one.
##


fn_target_distributions_from_accuracy_beliefs <- function(targets,
                                                          family_of_each_target,
                                                          constructor_arguments,
                                                          constructors,
                                                          n_target_draws,
                                                          seed) {

    targets <- fn_accuracy_prior_targets(targets = targets)
    identifiers <- targets$target_id

    if (any(targets$lower <= 0 | targets$upper >= 1 | targets$mean <= targets$lower
            | targets$mean >= targets$upper)) {
        stop("These constructors require an interior mean strictly between two interior interval endpoints.")
    }

    if (!is.character(family_of_each_target) || !identical(names(family_of_each_target), identifiers)
        || anyNA(family_of_each_target)
        || any(!(family_of_each_target %in% c("probit_normal", "dirichlet_Beta")))) {
        stop("Name an explicit probit_normal or dirichlet_Beta family for every target, in target order.")
    }

    if (!is.list(constructor_arguments) || !identical(names(constructor_arguments), identifiers)
        || !is.list(constructors) || is.null(names(constructors)) || anyDuplicated(names(constructors))
        || !all(unique(family_of_each_target) %in% names(constructors))) {
        stop("Supply named argument lists for every target and the existing selected constructor functions.")
    }

    for (family in unique(family_of_each_target)) {
        if (!is.function(constructors[[family]])) stop("Each selected constructor must be a function.")
    }

    if (!is.numeric(n_target_draws) || length(n_target_draws) != 1 || !is.finite(n_target_draws)
        || n_target_draws < 2 || n_target_draws != as.integer(n_target_draws)
        || !is.numeric(seed) || length(seed) != 1 || !is.finite(seed)
        || seed < 0 || seed != as.integer(seed)) {
        stop("Supply an explicit integer target-draw count of at least two and a nonnegative integer seed.")
    }

    distributions <- stats::setNames(vector("list", length(identifiers)), identifiers)
    comparison <- targets
    comparison$family <- unname(family_of_each_target)
    comparison$constructed_mean <- NA_real_
    comparison$constructed_lower <- NA_real_
    comparison$constructed_upper <- NA_real_

    for (i in seq_along(identifiers)) {

        family <- family_of_each_target[i]
        arguments <- constructor_arguments[[i]]
        common <- list(target_mean = targets$mean[i],
                        target_interval = c(targets$lower[i], targets$upper[i]),
                        CI_level = targets$interval_level[i],
                        verbose = TRUE)

        if (!is.list(arguments) || (length(arguments) > 0
            && (is.null(names(arguments)) || anyNA(names(arguments)) || any(!nzchar(names(arguments)))
                || anyDuplicated(names(arguments)))) || any(names(arguments) %in% names(common))) {
            stop("Constructor arguments must be uniquely named and must not replace the target table.")
        }

        full_arguments <- c(common, arguments)
        formal_names <- names(formals(constructors[[family]]))

        if ("..." %in% formal_names || !setequal(names(full_arguments), formal_names)) {
            stop(paste0("Supply every argument of the existing ", family,
                        " constructor explicitly; unused arguments and default values are not accepted."))
        }

        if (family == "dirichlet_Beta") {
            expected_class <- if (targets$quantity[i] == "Se") "present" else "absent"
            if (!identical(full_arguments$latent_class, expected_class)) {
                stop("The existing Dirichlet constructor requires present for Se and absent for Sp.")
            }
        }

        value <- do.call(what = constructors[[family]], args = full_arguments)

        if (!is.list(value) || !is.numeric(value$achieved_mean) || length(value$achieved_mean) != 1
            || !is.numeric(value$achieved_interval) || length(value$achieved_interval) != 2
            || any(!is.finite(c(value$achieved_mean, value$achieved_interval)))) {
            stop("The existing constructor must return its actual achieved mean and interval endpoints.")
        }

        if (family == "probit_normal") {

            if (!is.numeric(value$prior_mean) || length(value$prior_mean) != 1
                || !is.numeric(value$prior_sd) || length(value$prior_sd) != 1
                || !is.finite(value$prior_mean) || !is.finite(value$prior_sd) || value$prior_sd <= 0) {
                stop("The probit-normal constructor must return a finite mean and positive SD.")
            }
            parameters <- list(mean = value$prior_mean, SD = value$prior_sd)

        } else {

            alpha <- value$alpha_vector
            cutoff <- full_arguments$default_cutoff

            if (!is.numeric(alpha) || !is.null(dim(alpha)) || any(!is.finite(alpha)) || any(alpha <= 0)
                || length(alpha) != full_arguments$n_categories || length(cutoff) != 1
                || !is.finite(cutoff) || cutoff != as.integer(cutoff)
                || cutoff < 2 || cutoff > length(alpha)) {
                stop("The Dirichlet constructor must return every positive category concentration.")
            }

            below <- sum(alpha[seq_len(cutoff - 1)])
            above <- sum(alpha[seq.int(cutoff, length(alpha))])
            parameters <- if (targets$quantity[i] == "Se") {
                list(shape1 = above, shape2 = below)
            } else list(shape1 = below, shape2 = above)

        }

        distributions[[i]] <- list(family = unname(family),
                                    parameters = parameters,
                                    constructor_arguments = full_arguments,
                                    constructor_result = value)
        comparison$constructed_mean[i] <- value$achieved_mean
        comparison$constructed_lower[i] <- value$achieved_interval[1]
        comparison$constructed_upper[i] <- value$achieved_interval[2]

    }

    ##
    ## ---- Retain the caller's random-number state whilst generating the explicit target sample:
    ##

    original_kind <- RNGkind()
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    original_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL

    on.exit({
        RNGkind(kind = original_kind[1], normal.kind = original_kind[2], sample.kind = original_kind[3])
        if (had_seed) {
            assign(".Random.seed", original_seed, envir = .GlobalEnv)
        } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            rm(list = ".Random.seed", envir = .GlobalEnv)
        }
    }, add = TRUE)

    set.seed(seed = seed, kind = "L'Ecuyer-CMRG", normal.kind = original_kind[2], sample.kind = original_kind[3])
    message(paste0("Target draws per marginal: ", n_target_draws, "; seed: ", seed,
                   "; RNG: ", paste0(RNGkind(), collapse = ", "), "."))
    progress <- utils::txtProgressBar(min = 0, max = length(identifiers), initial = 0, style = 3)
    on.exit(close(progress), add = TRUE)
    draws <- stats::setNames(vector("list", length(identifiers)), identifiers)

    for (i in seq_along(identifiers)) {

        distribution <- distributions[[i]]
        parameters <- distribution$parameters
        probability <- if (distribution$family == "probit_normal") {
            stats::pnorm(q = stats::rnorm(n = n_target_draws, mean = parameters$mean, sd = parameters$SD),
                          mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
        } else {
            stats::rbeta(n = n_target_draws, shape1 = parameters$shape1, shape2 = parameters$shape2, ncp = 0)
        }

        if (any(!is.finite(probability)) || any(probability <= 0 | probability >= 1)) {
            stop("The target sample reaches a numerical boundary; no probabilities are clipped or discarded.")
        }

        draws[[i]] <- data.frame(target_id = identifiers[i],
                                 draw_row = seq_len(n_target_draws),
                                 test = targets$test[i],
                                 quantity = targets$quantity[i],
                                 raw_cutoff = targets$raw_cutoff[i],
                                 estimand = targets$estimand[i],
                                 probability = probability)
        utils::setTxtProgressBar(pb = progress, value = i)

    }

    comparison$mean_error <- comparison$constructed_mean - comparison$mean
    comparison$lower_error <- comparison$constructed_lower - comparison$lower
    comparison$upper_error <- comparison$constructed_upper - comparison$upper
    print(comparison)
    message("These are separately constructed target marginals, not a native-model prior or a joint curve.")
    message("The target-draw count controls the artificial likelihood information used for derivation.")

    list(targets = targets,
         distributions = distributions,
         comparison = comparison,
         target_draws = do.call(what = rbind, args = unname(draws)),
         n_target_draws = n_target_draws,
         seed = seed,
         RNG = RNGkind(),
         scope = "separate elicited marginal distributions; no joint dependence specified")

}






















