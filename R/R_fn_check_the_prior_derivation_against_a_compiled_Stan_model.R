##
## ---- Check a compiled derivation model against an independent reference at supplied saved points:
##
## The compiled model already contains its explicitly supplied Stan data. No compilation, fit or
## optimisation is started here. Only vector parameter blocks are supported; their names and order
## must match the compiled model. The density includes constants and excludes the transform Jacobian.
##


fn_check_the_prior_derivation_against_a_compiled_Stan_model <- function(saved_parameter_points,
                                                                       compiled_model,
                                                                       fn_reference_evaluation,
                                                                       relative_tolerance) {

      if (!inherits(compiled_model, "StanModel")) {
            stop("Supply the actual compiled BridgeStan model with its matching Stan data already attached.")
      }

      required_methods <- c("param_names", "param_unc_names", "param_num", "param_unc_num",
                            "param_unconstrain", "log_density_gradient", "name", "model_info")

      for (method in required_methods) {
            if (!is.function(compiled_model[[method]])) {
                  stop(paste0("The compiled model does not expose the required method ", method, "."))
            }
      }

      actual_constrained_names <- compiled_model$param_names(include_tp = FALSE, include_gq = FALSE)
      actual_unconstrained_names <- compiled_model$param_unc_names()
      n_constrained <- compiled_model$param_num(include_tp = FALSE, include_gq = FALSE)
      n_unconstrained <- compiled_model$param_unc_num()

      if (n_constrained < 1 || n_constrained != n_unconstrained) {
            stop(paste0("This check requires non-empty vector blocks with only entrywise transforms. ",
                        "Use the fixed-input check-only path when there are no free parameters."))
      }

      if (length(actual_constrained_names) != n_constrained
            || length(actual_unconstrained_names) != n_unconstrained) {
            stop("The compiled parameter names do not match its declared parameter counts.")
      }

      fn_Stan_evaluation <- function(point) {

            if (!is.list(point) || is.null(names(point)) || anyNA(names(point))
                  || any(names(point) == "") || anyDuplicated(names(point))) {
                  stop("Every saved point must contain uniquely named vector parameter blocks.")
            }

            expected_Stan_names <- character(0)
            reference_names <- character(0)

            for (block in names(point)) {

                  values <- point[[block]]

                  if (!is.numeric(values) || !is.null(dim(values)) || any(!is.finite(values))) {
                        stop(paste0("The saved parameter block ", block, " must be a finite numeric vector."))
                  }

                  if (length(values) > 0) {
                        expected_Stan_names <- c(expected_Stan_names,
                                                 paste0(block, ".", seq_along(values)))
                        reference_names <- c(reference_names,
                                              paste0(block, "[", seq_along(values), "]"))
                  }

            }

            if (!identical(actual_constrained_names, expected_Stan_names)
                  || !identical(actual_unconstrained_names, expected_Stan_names)) {
                  stop(paste0("Saved-point block names, sizes or order differ from the compiled model. ",
                              "No positional renaming or parameter reordering is used."))
            }

            constrained <- unlist(point, recursive = TRUE, use.names = FALSE)
            unconstrained <- compiled_model$param_unconstrain(theta = constrained)

            if (length(unconstrained) != n_unconstrained || any(!is.finite(unconstrained))) {
                  stop("The compiled model could not transform the supplied saved point to finite coordinates.")
            }

            evaluation <- compiled_model$log_density_gradient(theta_unc = unconstrained,
                                                               propto = FALSE,
                                                               jacobian = FALSE)

            if (!is.list(evaluation) || !all(c("val", "gradient") %in% names(evaluation))
                  || length(evaluation$gradient) != n_unconstrained) {
                  stop("The compiled model did not return the full density and gradient.")
            }

            list(log_density = evaluation$val,
                 gradient_unconstrained = stats::setNames(evaluation$gradient, reference_names))

      }

      message(paste0("Compiled model: ", compiled_model$name(), "."))
      message(paste0("Compiled model information: ", compiled_model$model_info(), "."))
      message(paste0("Checking ", n_unconstrained, " main parameters; no nuisance parameters are declared."))
      message("Density constants included; parameter-transform Jacobian excluded; no prior values estimated.")

      result <- fn_check_the_log_density_and_gradient_of_a_prior_derivation_model(
                      saved_parameter_points = saved_parameter_points,
                      fn_Stan_evaluation = fn_Stan_evaluation,
                      fn_reference_evaluation = fn_reference_evaluation,
                      relative_tolerance = relative_tolerance)

      result$compiled_model_name <- compiled_model$name()
      result$compiled_model_information <- compiled_model$model_info()
      result$compiled_constrained_parameter_names <- actual_constrained_names
      result$compiled_unconstrained_parameter_names <- actual_unconstrained_names

      return(result)

}






















