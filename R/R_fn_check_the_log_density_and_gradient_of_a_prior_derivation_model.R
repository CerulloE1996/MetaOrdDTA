##
## ---- Compare a Stan model with an independent reference at explicitly supplied saved points:
##
## Both callbacks return log_density and a named gradient_unconstrained vector. They evaluate the
## same likelihood without a parameter-transform Jacobian. Neither callback estimates prior values.
## Integration-bank rows are data; all parameters of these derivation models belong to the main block.
## Checks of a receiving model's own main and nuisance blocks remain separate.
##


fn_check_the_log_density_and_gradient_of_a_prior_derivation_model <- function(saved_parameter_points,
                                                                             fn_Stan_evaluation,
                                                                             fn_reference_evaluation,
                                                                             relative_tolerance) {

      if (!is.list(saved_parameter_points) || length(saved_parameter_points) < 1
            || !all(vapply(saved_parameter_points, is.list, logical(1)))) {
            stop("saved_parameter_points must be a non-empty list of explicitly supplied parameter lists.");
      }

      if (!is.function(fn_Stan_evaluation) || !is.function(fn_reference_evaluation)) {
            stop("Supply the actual Stan evaluation callback and an independent reference callback.");
      }

      if (!is.numeric(relative_tolerance) || length(relative_tolerance) != 1
            || !is.finite(relative_tolerance) || relative_tolerance <= 0) {
            stop("relative_tolerance must be supplied explicitly as a finite positive number.");
      }

      n_points <- length(saved_parameter_points)
      rows <- vector(mode = "list", length = n_points)
      required_fields <- c("log_density", "gradient_unconstrained")
      progress <- utils::txtProgressBar(min = 0, max = n_points, initial = 0, style = 3)
      on.exit(close(progress), add = TRUE)

      for (point_i in seq_len(n_points)) {

            point <- saved_parameter_points[[point_i]]
            Stan <- fn_Stan_evaluation(point)
            reference <- fn_reference_evaluation(point)

            for (evaluation in list(Stan, reference)) {

                  if (!is.list(evaluation) || !all(required_fields %in% names(evaluation))) {
                        stop("Each callback must return log_density and gradient_unconstrained.");
                  }

                  log_density <- evaluation$log_density
                  gradient <- evaluation$gradient_unconstrained

                  if (!is.numeric(log_density) || length(log_density) != 1 || !is.finite(log_density)
                        || !is.numeric(gradient) || length(gradient) < 1 || !is.null(dim(gradient))
                        || any(!is.finite(gradient)) || is.null(names(gradient))
                        || anyNA(names(gradient)) || any(names(gradient) == "")
                        || anyDuplicated(names(gradient))) {
                        stop(paste0("Each evaluation must have a finite scalar density and a finite ",
                                    "gradient with unique explicit parameter names."));
                  }

            }

            if (!identical(names(Stan$gradient_unconstrained), names(reference$gradient_unconstrained))) {
                  stop(paste0("Stan and reference gradient names or order differ at point ", point_i,
                              ". Parameter order is never inferred."));
            }

            log_density_absolute_error <- abs(Stan$log_density - reference$log_density)
            gradient_absolute_error <- max(abs(Stan$gradient_unconstrained
                                                - reference$gradient_unconstrained))
            log_density_scale <- abs(reference$log_density)
            gradient_scale <- max(abs(reference$gradient_unconstrained))

            ##
            ## ---- Zero reference scales have no relative error; retain their absolute errors explicitly:
            ##

            log_density_relative_error <- if (log_density_scale > 0) {
                  log_density_absolute_error / log_density_scale
            } else {
                  NA_real_
            }

            gradient_relative_error <- if (gradient_scale > 0) {
                  gradient_absolute_error / gradient_scale
            } else {
                  NA_real_
            }

            scales_are_nonzero <- log_density_scale > 0 && gradient_scale > 0
            passed <- if (scales_are_nonzero) {
                  log_density_relative_error <= relative_tolerance &&
                        gradient_relative_error <= relative_tolerance
            } else {
                  NA
            }

            status <- if (!scales_are_nonzero) {
                  "ZERO_REFERENCE_SCALE"
            } else if (passed) {
                  "PASS"
            } else {
                  "FAIL"
            }

            rows[[point_i]] <- data.frame(point = point_i,
                                          n_main_parameters = length(Stan$gradient_unconstrained),
                                          n_nuisance_parameters = 0,
                                          Stan_log_density = Stan$log_density,
                                          reference_log_density = reference$log_density,
                                          log_density_absolute_error = log_density_absolute_error,
                                          log_density_relative_error = log_density_relative_error,
                                          gradient_absolute_error = gradient_absolute_error,
                                          gradient_relative_error = gradient_relative_error,
                                          relative_tolerance = relative_tolerance,
                                          passed = passed,
                                          status = status)

            utils::setTxtProgressBar(pb = progress, value = point_i)
            message(paste0("Point ", point_i, "/", n_points, ": ", status,
                           "; relative density error ", log_density_relative_error,
                           "; relative gradient error ", gradient_relative_error, "."))

      }

      results <- do.call(what = rbind, args = rows)
      print(results, row.names = FALSE)

      list(results = results,
           all_points_passed = all(!is.na(results$passed) & results$passed),
           parameter_names = names(Stan$gradient_unconstrained),
           scope = "derivation-model density and all main gradients; no receiving-model nuisance pass implied")

}






















