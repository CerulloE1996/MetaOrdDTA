##
## ---- Direct Se/Sp priors: preparation, gradients and prior-only sampling:
##
fn_sample_call_Se_Sp_prior <- function(expression_to_search,
                                       name_of_the_argument) {

      if (is.call(expression_to_search)) {
            if (name_of_the_argument %in% names(expression_to_search)) return(expression_to_search)
            for (part in as.list(expression_to_search)) {
                  call_found <- fn_sample_call_Se_Sp_prior(part, name_of_the_argument)
                  if (!is.null(call_found)) return(call_found)
            }
      }
      return(NULL)

}
##
fn_Se_Sp_prior_lp_R <- function(u,
                                 parameter_names,
                                 Stan_data) {

      block <- function(name) {

            u[startsWith(parameter_names, paste0(name, ".")) | parameter_names == name]

      }
      stopifnot(Stan_data$n_covariates_max == 1, Stan_data$use_probit_link == 1,
                all(Stan_data$K_fold_CV_indicator == 0), length(Stan_data$use_Se_Sp_prior) == 1,
                Stan_data$use_Se_Sp_prior %in% c(0, 1))
      if (Stan_data$use_Se_Sp_prior == 1) {
            stopifnot(Stan_data$Se_Sp_prior_type %in% c(1, 4),
                      all(Stan_data$n_chosen == Stan_data$n_thr))
      }
      mu <- matrix(block("beta_mu"), nrow = 2)
      SD <- exp(block("beta_SD"))
      z <- block("beta_z")
      rho_u <- block("beta_corr")
      rho_p <- stats::plogis(rho_u)
      rho_width <- Stan_data$beta_corr_ub - Stan_data$beta_corr_lb
      rho <- Stan_data$beta_corr_lb + rho_width * rho_p
      raw_u <- matrix(block("C_raw_vec"), nrow = Stan_data$n_thr)
      raw_p <- stats::plogis(raw_u)
      raw_C <- -15 + 17.5 * raw_p
      lp <- sum(stats::dnorm(mu, Stan_data$prior_beta_mu_mean, Stan_data$prior_beta_mu_SD,
                             log = TRUE))
      lp <- lp + sum(stats::dnorm(SD, Stan_data$prior_beta_SD_mean, Stan_data$prior_beta_SD_SD,
                                   log = TRUE)) + sum(log(SD))
      lp <- lp + sum(stats::dnorm(z, log = TRUE))
      lp <- lp + (Stan_data$prior_beta_corr_LKJ - 1) * log1p(-rho^2) -
                 lbeta(0.5, Stan_data$prior_beta_corr_LKJ)
      lp <- lp + log(rho_width) + log(rho_p) + log1p(-rho_p)
      lp <- lp + sum(log(17.5) + log(raw_p) + log1p(-raw_p))
      for (c in 1:2) {

            raw <- raw_C[, c]
            increments <- if (Stan_data$softplus == 1) log1p(exp(raw[-1])) else exp(raw[-1])
            C <- cumsum(c(raw[1], increments))
            baseline_case <- if (c == 1) Stan_data$baseline_case_nd else Stan_data$baseline_case_d
            C_baseline <- C - sum(baseline_case * mu[c, ])
            if (Stan_data$use_Se_Sp_prior == 0) {
                  q <- diff(c(0, stats::pnorm(C), 1))
                  A <- Stan_data$prior_dirichlet_alpha[[c]]
                  lp <- lp + lgamma(sum(A)) - sum(lgamma(A)) + sum((A - 1) * log(q))
                  lp <- lp + sum(stats::dnorm(C, log = TRUE))
            } else if (Stan_data$Se_Sp_prior_type == 1) {
                  q <- diff(c(0, stats::pnorm(C_baseline), 1))
                  A <- Stan_data$prior_Se_Sp_dirichlet_A[[c]]
                  lp <- lp + lgamma(sum(A)) - sum(lgamma(A)) + sum((A - 1) * log(q))
                  lp <- lp + sum(stats::dnorm(C_baseline, log = TRUE))
            } else {
                  mean_SD <- Stan_data$prior_Se_Sp_link_mean_SD[[c]]
                  link_accuracy <- if (c == 1) C_baseline else - C_baseline
                  lp <- lp + sum(stats::dnorm(link_accuracy, mean_SD[, 1], mean_SD[, 2], log = TRUE))
            }
            lp <- lp + if (Stan_data$softplus == 1) sum(log(stats::plogis(raw[-1]))) else sum(raw[-1])

      }
      return(lp)

}
##
fn_Se_Sp_prior_chunk_gradient <- function(so_file,
                                           JSON_text,
                                           points,
                                           main_positions,
                                           nuisance_positions,
                                           seed) {

      library_name <- tools::file_path_sans_ext(basename(so_file))
      dyn.load(so_file)
      construction <- .C("bs_model_construct_R", as.character(JSON_text), as.integer(seed),
                          ptr_out = raw(8), err_msg = as.character(""), err_ptr = raw(8),
                          PACKAGE = library_name)
      if (all(construction$ptr_out == 0)) stop(construction$err_msg)
      on.exit(.C("bs_model_destruct_R", as.raw(construction$ptr_out), PACKAGE = library_name))
      main_names <- strsplit(.C("bs_param_unc_names_R", as.raw(construction$ptr_out),
                                 names_out = as.character(""), PACKAGE = library_name)$names_out, ",")[[1]]
      if (!identical(main_names, names(points[[1]])[main_positions])) {
            stop("The checked library has a different main-parameter ordering.")
      }
      return(lapply(points, function(point) {

            evaluation <- .C("nicostan_bs_log_density_gradient_with_the_nuisance_block_R",
                              as.raw(construction$ptr_out), as.integer(0), as.integer(1),
                              as.double(point[main_positions]), as.double(point[nuisance_positions]),
                              nuisance_gradient = double(length(nuisance_positions)),
                              as.integer(length(nuisance_positions)), val = double(1),
                              gradient = double(length(main_positions)), return_code = as.integer(0),
                              err_msg = as.character(""), err_ptr = raw(8), NAOK = TRUE, PACKAGE = library_name)
            if (evaluation$return_code != 0) stop(evaluation$err_msg)
            gradient <- numeric(length(point))
            gradient[main_positions] <- evaluation$gradient
            gradient[nuisance_positions] <- evaluation$nuisance_gradient
            return(list(val = evaluation$val, gradient = gradient))

      }))

}
##
fn_check_Se_Sp_prior_gradients <- function(prepared,
                                            Stan_file,
                                            init_lists_per_chain,
                                            seed,
                                            output_prefix,
                                            tolerance,
                                            finite_difference_step) {

      data_file <- paste0(output_prefix, "_data.json")
      cmdstanr::write_stan_json(data = prepared$Stan_data_list, file = data_file)
      JSON_text <- paste(readLines(data_file, warn = FALSE), collapse = "\n")
      functions_folder <- system.file("stan_functions", package = "MetaOrdDTA")
      message("[2/4] Building the independent full-parameter gradient reference.")
      reference_library <- bridgestan::compile_model(
                                  Stan_file,
                                  stanc_args = paste0("--include-paths=", functions_folder))
      reference_model <- bridgestan::StanModel$new(reference_library, JSON_text, seed)
      parameter_names <- reference_model$param_names(include_tp = FALSE, include_gq = FALSE)
      parameter_blocks <- sub("\\..*$", "", parameter_names)
      unconstrained_names <- reference_model$param_unc_names()
      points <- lapply(init_lists_per_chain, function(init) {

            constrained <- numeric(length(parameter_names))
            for (block in unique(parameter_blocks)) {
                  positions <- which(parameter_blocks == block)
                  block_init <- init[[block]]
                  if (block == "C_raw_vec") {
                        block_init <- do.call(rbind, block_init)
                  }
                  block_values <- as.vector(block_init)
                  stopifnot(length(block_values) == length(positions))
                  constrained[positions] <- as.numeric(block_values)
            }
            names(constrained) <- parameter_names
            for (c in 1:2) {
                  C_names <- paste0("C_raw_vec.", c, ".", seq_len(prepared$Stan_data_list$n_thr))
                  stopifnot(all(constrained[C_names] == init$C_raw_vec[[c]]))
            }
            point <- reference_model$param_unconstrain(theta = constrained)
            return(stats::setNames(as.numeric(point), unconstrained_names))

      })
      nuisance_positions <- which(startsWith(unconstrained_names, "beta_z."))
      main_positions <- setdiff(seq_along(unconstrained_names), nuisance_positions)
      checked <- callr::r(func = fn_Se_Sp_prior_chunk_gradient,
                           args = list(so_file = prepared$so_file, JSON_text = JSON_text, points = points,
                                       main_positions = main_positions, nuisance_positions = nuisance_positions,
                                       seed = seed))
      rows <- lapply(seq_along(points), function(i) {

            u <- points[[i]]
            reference <- reference_model$log_density_gradient(theta_unc = u, propto = FALSE, jacobian = TRUE)
            lp_R <- fn_Se_Sp_prior_lp_R(u, unconstrained_names, prepared$Stan_data_list)
            gradient_R <- vapply(seq_along(u), function(j) {

                  e <- replace(numeric(length(u)), j, finite_difference_step)
                  (fn_Se_Sp_prior_lp_R(u - 2 * e, unconstrained_names, prepared$Stan_data_list) -
                   8 * fn_Se_Sp_prior_lp_R(u - e, unconstrained_names, prepared$Stan_data_list) +
                   8 * fn_Se_Sp_prior_lp_R(u + e, unconstrained_names, prepared$Stan_data_list) -
                   fn_Se_Sp_prior_lp_R(u + 2 * e, unconstrained_names, prepared$Stan_data_list)) /
                   (12 * finite_difference_step)

            }, numeric(1))
            scale <- max(abs(reference$gradient))
            row <- data.frame(point = i,
                              lp_R_relative = abs(lp_R - reference$val) / max(1, abs(reference$val)),
                              lp_NicoStan_relative = abs(checked[[i]]$val - reference$val) /
                                                    max(1, abs(reference$val)),
                              gradient_R_main_relative = max(abs(gradient_R - reference$gradient)[
                                                               main_positions]) / scale,
                              gradient_R_nuisance_relative = max(abs(gradient_R - reference$gradient)[
                                                                   nuisance_positions]) / scale,
                              gradient_NicoStan_main_relative = max(abs(checked[[i]]$gradient -
                                                                         reference$gradient)[main_positions]) /
                                                                 scale,
                              gradient_NicoStan_nuisance_relative = max(abs(checked[[i]]$gradient -
                                                                             reference$gradient)[
                                                                             nuisance_positions]) / scale)
            print(row, digits = 8)
            return(row)

      })
      rows <- do.call(rbind, rows)
      passed <- all(is.finite(as.matrix(rows))) && all(as.matrix(rows[, -1]) <= tolerance)
      result <- list(passed = passed, rows = rows, points = points, main_positions = main_positions,
                     nuisance_positions = nuisance_positions, checked_library = prepared$so_file,
                     reference_library = reference_library, tolerance = tolerance,
                     finite_difference_step = finite_difference_step)
      saveRDS(result, paste0(output_prefix, "_gradients.rds"))
      if (!passed) stop("Gradient gate failed; no sampling started.")
      message("PASS: full-parameter density and main/nuisance gradients; sampling gate cleared.")
      return(result)

}
##
fn_validate_Se_Sp_prior_comparison_reference <- function(Stan_data,
                                                          comparison_reference) {

      analytic_supported <- FALSE
      if (length(Stan_data$use_Se_Sp_prior) != 1 ||
          !Stan_data$use_Se_Sp_prior %in% c(0, 1)) {
            stop("use_Se_Sp_prior must be 0 or 1.")
      }
      if (Stan_data$use_Se_Sp_prior == 1) {
            if (length(Stan_data$n_chosen) != 2 || any(!is.finite(Stan_data$n_chosen)) ||
                any(Stan_data$n_chosen != as.integer(Stan_data$n_chosen)) ||
                any(Stan_data$n_chosen < 0) || any(Stan_data$n_chosen > Stan_data$n_thr) ||
                sum(Stan_data$n_chosen) == 0) {
                  stop("Active Se/Sp comparisons need valid selected-cutoff counts.")
            }
            for (c in seq_len(2)) {
                  m <- Stan_data$n_chosen[c]
                  if (m > 0) {
                        chosen <- Stan_data$chosen_cutoffs[c, seq_len(m)]
                        if (length(chosen) != m || any(!is.finite(chosen)) ||
                            any(chosen != as.integer(chosen)) || any(chosen < 1) ||
                            any(chosen > Stan_data$n_thr) || any(diff(chosen) <= 0)) {
                              stop("Chosen cutoffs must be increasing integer indices within n_thr.")
                        }
                  }
            }
            if (is.null(comparison_reference) && Stan_data$Se_Sp_prior_type == 1) {
                  analytic_supported <- all(vapply(seq_len(2), function(c) {

                        m <- Stan_data$n_chosen[c]
                        if (m == 0) return(TRUE)
                        A <- Stan_data$prior_Se_Sp_dirichlet_A[[c]]
                        length(A) >= m + 1 && all(is.finite(A[seq_len(m + 1)])) &&
                              all(A[seq_len(m + 1)] > 0)

                  }, logical(1)))
            }
            if (is.null(comparison_reference) && Stan_data$Se_Sp_prior_type == 4) {
                  analytic_supported <- all(vapply(seq_len(2), function(c) {

                        m <- Stan_data$n_chosen[c]
                        if (m == 0) return(TRUE)
                        mean_SD <- Stan_data$prior_Se_Sp_link_mean_SD[[c]]
                        is.matrix(mean_SD) && nrow(mean_SD) >= m && ncol(mean_SD) == 2 &&
                              all(is.finite(mean_SD[seq_len(m), , drop = FALSE])) &&
                              all(mean_SD[seq_len(m), 2] > 0) &&
                              length(unique(mean_SD[seq_len(m), 1])) == 1 &&
                              length(unique(mean_SD[seq_len(m), 2])) == 1

                  }, logical(1)))
            }
      }
      if (is.null(comparison_reference)) {
            if (!analytic_supported) {
                  stop("A supplied reference is required for this prior comparison.")
            }
            return("typed-prior analytical target")
      }
      stopifnot(is.list(comparison_reference), all(c("Se", "Sp") %in% names(comparison_reference)))
      for (variable in c("Se", "Sp")) {
            reference <- comparison_reference[[variable]]
            if (!is.matrix(reference) || nrow(reference) < 1 || ncol(reference) != Stan_data$n_thr ||
                any(!is.finite(reference)) || any(reference < 0) || any(reference > 1)) {
                  stop("Reference draws must provide finite Se and Sp matrices by cutoff.")
            }
      }
      stopifnot(nrow(comparison_reference$Se) == nrow(comparison_reference$Sp))
      return("supplied reference draws")

}
##
fn_Se_Sp_prior_comparison <- function(draws_array,
                                       Stan_data,
                                       probabilities,
                                       gap_limit,
                                       comparison_reference) {

      rows <- list()
      reference_type <- fn_validate_Se_Sp_prior_comparison_reference(Stan_data, comparison_reference)
      stopifnot(Stan_data$use_probit_link == 1, length(probabilities) > 0,
                all(is.finite(probabilities)), all(probabilities >= 0), all(probabilities <= 1),
                length(dim(draws_array)) == 3, !is.null(dimnames(draws_array)[[3]]))
      for (c in 1:2) {
            variable <- if (c == 1) "Sp" else "Se"
            m <- if (Stan_data$use_Se_Sp_prior == 0) Stan_data$n_thr else Stan_data$n_chosen[c]
            if (m == 0) next
            cutoffs <- if (Stan_data$use_Se_Sp_prior == 0) seq_len(Stan_data$n_thr) else
                                                                  Stan_data$chosen_cutoffs[c, seq_len(m)]
            if (reference_type == "supplied reference draws") {
                  reference_values <- comparison_reference[[variable]][, cutoffs, drop = FALSE]
            } else if (Stan_data$Se_Sp_prior_type == 1) {
                  A <- Stan_data$prior_Se_Sp_dirichlet_A[[c]][seq_len(m + 1)]
            } else {
                  mean_SD <- Stan_data$prior_Se_Sp_link_mean_SD[[c]][seq_len(m), , drop = FALSE]
            }
            for (j in seq_len(m)) {
                  cutoff <- cutoffs[j]
                  name <- paste0(variable, "_baseline[", cutoff, "]")
                  if (!name %in% dimnames(draws_array)[[3]]) {
                        stop(paste0("Required draw is missing: ", name))
                  }
                  values <- as.numeric(draws_array[, , name])
                  if (any(!is.finite(values))) stop(paste0("Non-finite draws found for ", name))
                  if (reference_type == "supplied reference draws") {
                        expected <- as.numeric(stats::quantile(reference_values[, j], probabilities))
                  } else if (Stan_data$Se_Sp_prior_type == 1) {
                        first <- sum(A[seq_len(j)])
                        last <- sum(A) - first
                        expected <- if (c == 1) stats::qbeta(probabilities, first, last) else
                                                      stats::qbeta(probabilities, last, first)
                  } else {
                        first <- if (c == 1) j else m + 1 - j
                        last <- if (c == 1) m + 1 - j else j
                        expected <- stats::qbeta(probabilities, first, last)
                        expected <- stats::pnorm(mean_SD[j, 1] + mean_SD[j, 2] *
                                                   stats::qnorm(expected))
                  }
                  observed <- as.numeric(stats::quantile(values, probabilities))
                  rows[[length(rows) + 1]] <- data.frame(variable = variable, cutoff = cutoff,
                                                         probability = probabilities, target = expected,
                                                         observed = observed,
                                                         absolute_gap = abs(observed - expected))
            }
      }
      table <- if (length(rows) == 0) {

            data.frame(variable = character(), cutoff = integer(), probability = numeric(), target = numeric(),
                       observed = numeric(), absolute_gap = numeric())

      } else do.call(rbind, rows)
      summary <- do.call(rbind, lapply(c("Se", "Sp", "both"), function(variable) {

            gaps <- if (variable == "both") table$absolute_gap else
                                              table$absolute_gap[table$variable == variable]
            return(data.frame(variable = variable,
                              mean_gap = if (length(gaps) == 0) NA_real_ else mean(gaps),
                              worst_gap = if (length(gaps) == 0) NA_real_ else max(gaps),
                              n_over_limit = sum(gaps > gap_limit), n_cells = length(gaps),
                              gap_limit = gap_limit))

      }))
      return(list(table = table, summary = summary, reference_type = reference_type))

}
##
fn_run_Se_Sp_prior_check <- function(Stan_data,
                                       Stan_file,
                                       init_lists_per_chain,
                                       sampler_args,
                                       n_chains,
                                       n_burnin,
                                       n_iter,
                                       seed,
                                       output_prefix,
                                       gradient_tolerance,
                                       finite_difference_step,
                                       probabilities,
                                       gap_limit,
                                       comparison_reference) {

      result_file <- paste0(output_prefix, "_results.rds")
      inputs_file <- paste0(output_prefix, "_inputs.rds")
      stopifnot(length(Stan_data$use_Se_Sp_prior) == 1, Stan_data$use_Se_Sp_prior %in% c(0, 1))
      requested_inputs <- list(Stan_data = Stan_data, sampler_args = sampler_args,
                               init_lists_per_chain = init_lists_per_chain,
                               n_chains = n_chains, n_burnin = n_burnin, n_iter = n_iter, seed = seed,
                               comparison_reference = comparison_reference,
                               Stan_file = normalizePath(Stan_file, mustWork = FALSE))
      if (file.exists(result_file)) {
            if (!file.exists(inputs_file)) stop("Saved inputs are missing; use a different output prefix.")
            saved_inputs <- readRDS(inputs_file)
            saved_inputs$Stan_file <- normalizePath(saved_inputs$Stan_file, mustWork = FALSE)
            inputs_match <- vapply(names(requested_inputs), function(name) {

                  identical(saved_inputs[[name]], requested_inputs[[name]])

            }, logical(1))
            if (!all(inputs_match)) {
                  stop(paste0("Saved inputs differ: ", paste(names(inputs_match)[!inputs_match], collapse = ", "),
                              ". Use a different output prefix."))
            }
            result <- readRDS(result_file)
            stopifnot(identical(result$gradients$tolerance, gradient_tolerance),
                      identical(result$gradients$finite_difference_step, finite_difference_step),
                      identical(unique(result$comparison$table$probability), probabilities),
                      identical(unique(result$comparison$summary$gap_limit), gap_limit))
            result$reuse_status <- "archived: no current-source gradient check or fit performed"
            message(paste0("Archived result: ", result_file))
            message("Saved inputs match; no current-source gradient check or fit was performed.")
            print(result$comparison$summary)
            return(invisible(result))
      }
      stopifnot(all(Stan_data$K_fold_CV_indicator == 0), length(init_lists_per_chain) == n_chains)
      stopifnot(Stan_data$n_covariates_max == 1, Stan_data$use_probit_link == 1)
      stopifnot(length(probabilities) > 0, all(is.finite(probabilities)),
                all(probabilities >= 0), all(probabilities <= 1), length(gap_limit) == 1,
                is.finite(gap_limit), gap_limit >= 0)
      if (Stan_data$use_Se_Sp_prior == 1) {
            stopifnot(Stan_data$Se_Sp_prior_type %in% c(1, 4),
                      all(Stan_data$n_chosen == Stan_data$n_thr))
      }
      fn_validate_Se_Sp_prior_comparison_reference(Stan_data, comparison_reference)
      stopifnot(n_burnin <= 250, n_iter <= 250,
                identical(sampler_args$nuisance_params, "beta_z"))
      message(paste0("NicoStan: ", n_chains, " chains; ", n_iter, " sampling per chain (",
                     n_chains * n_iter, " retained draws); ", n_burnin, " burn-in per chain (",
                     n_chains * n_burnin, " burn-in transitions); total ", n_chains * (n_burnin + n_iter), "."))
      message("Nuisance: beta_z only. Fixed-cutpoint C_raw_vec remains main.")
      print(list(Stan_file = Stan_file, Stan_data = Stan_data, sampler_args = sampler_args,
                 n_chains = n_chains, n_burnin = n_burnin, n_iter = n_iter, seed = seed))
      files_of_Stan_fns <- list.files(system.file("stan_functions", package = "MetaOrdDTA"),
                                      pattern = "\\.stan$", full.names = TRUE)
      names(files_of_Stan_fns) <- basename(files_of_Stan_fns)
      source_text <- list(Stan_model = readLines(Stan_file, warn = FALSE),
                          Stan_functions = lapply(files_of_Stan_fns, readLines, warn = FALSE))
      saveRDS(c(requested_inputs, list(source_text = source_text)), inputs_file)
      preparation_args <- sampler_args
      preparation_args$Stan_file_of_the_fused_version <- Stan_file
      preparation_args$prepare_only <- TRUE
      message("[1/4] Preparing the NicoStan library, without sampling.")
      preparation_start <- Sys.time()
      prepared <- MetaOrdDTA:::R_fn_sample_model_with_NicoStan(
                                    stan_data_list = Stan_data, stan_model_file_name = basename(Stan_file),
                                    n_chains = n_chains, init_lists_per_chain = init_lists_per_chain, seed = seed,
                                    n_burnin = n_burnin, n_iter = n_iter, NicoStan_args = preparation_args)
      preparation_seconds <- as.numeric(difftime(Sys.time(), preparation_start, units = "secs"))
      gradients <- fn_check_Se_Sp_prior_gradients(prepared, Stan_file, init_lists_per_chain, seed, output_prefix,
                                                   gradient_tolerance, finite_difference_step)
      fit_args <- sampler_args
      fit_args$nuisance_params <- NULL
      fit_args$sample_nuisance <- NULL
      fit_args$Stan_file_of_the_fused_version <- NULL
      fit_args$prepare_only <- NULL
      fit_args$nuisance_to_track <- "beta_z"
      fit_args$init_lists_per_chain <- init_lists_per_chain
      fit_args$n_chains_burnin <- n_chains
      fit_args$n_chains_sampling <- n_chains
      fit_args$n_burnin <- n_burnin
      fit_args$n_iter <- n_iter
      fit_args$seed <- seed
      message("[3/4] Gradient gate passed. Starting the prior-only fit.")
      start <- Sys.time()
      do.call(prepared$NicoStan_model_obj$sample, fit_args)
      fit_seconds <- as.numeric(difftime(Sys.time(), start, units = "secs"))
      summary_obj <- prepared$NicoStan_model_obj$summary(compute_main_params = TRUE,
                                                          compute_transformed_parameters = FALSE,
                                                          compute_generated_quantities = TRUE,
                                                          save_log_lik_trace = FALSE)
      draws <- summary_obj$get_posterior_draws()
      divergences <- summary_obj$get_divergences()
      selected <- grep("^(Se|Sp)_baseline\\[", dimnames(draws)[[3]])
      diagnostics <- posterior::summarise_draws(posterior::as_draws_array(draws[, , selected, drop = FALSE]),
                                                 "mean", "sd", "rhat", "ess_bulk", "ess_tail", "ess_sd")
      comparison <- fn_Se_Sp_prior_comparison(draws, Stan_data, probabilities, gap_limit,
                                               comparison_reference)
      result <- list(comparison = comparison, diagnostics = diagnostics, draws_array = draws,
                     divergences = divergences, gradients = gradients, Stan_data = Stan_data,
                     sampler_args = fit_args, Stan_file = Stan_file, compiled_library = prepared$so_file,
                     source_text = source_text,
                     reuse_status = "fresh: gradient gate and prior-only fit performed",
                     fit_seconds = fit_seconds, preparation_seconds = preparation_seconds,
                     process_times = proc.time(), n_chains = n_chains, n_burnin = n_burnin, n_iter = n_iter,
                     seed = seed, bounded_raw_cutpoints = c(-15, 2.5))
      saveRDS(result, result_file)
      utils::write.csv(comparison$table, paste0(output_prefix, "_Se_Sp.csv"), row.names = FALSE)
      utils::write.csv(diagnostics, paste0(output_prefix, "_diagnostics.csv"), row.names = FALSE)
      message(paste0("[4/4] Actual native Se/Sp compared with ", comparison$reference_type, "."))
      print(comparison$table, digits = 6)
      print(comparison$summary, digits = 6)
      print(diagnostics)
      print(divergences)
      message(paste0("Fit seconds: ", fit_seconds, "; results: ", result_file))
      message("Native raw-cutpoint bounds are retained; the comparison tests their effect too.")
      return(invisible(result))

}





















