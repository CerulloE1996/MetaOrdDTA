#' Locate generated MetaOrdDTA outputs
#'
#' Keeps generated analyses and simulation results outside the R package source
#' tree. The root directory is selected, in order, from the
#' `MetaOrdDTA.output_dir` option, the `METAORDDTA_OUTPUT_DIR` environment
#' variable, a `MetaOrdDTA_outputs` directory beside the current MetaOrdDTA
#' project, or `~/MetaOrdDTA_outputs`.
#'
#' @param ... One or more path components below the output root.
#' @param create If `TRUE`, create the resulting directory recursively.
#'
#' @return A normalized character path. The path need not exist unless
#'   `create = TRUE`.
#' @export
metaorddta_output_path <- function(..., create = FALSE) {
  root <- getOption("MetaOrdDTA.output_dir", NULL)

  if (is.null(root) || length(root) != 1L || !nzchar(root)) {
    root <- Sys.getenv("METAORDDTA_OUTPUT_DIR", unset = "")
  }

  if (!nzchar(root)) {
    project_root <- find_metaorddta_project_root()
    root <- if (is.null(project_root)) {
      file.path(path.expand("~"), "MetaOrdDTA_outputs")
    } else {
      file.path(dirname(project_root), "MetaOrdDTA_outputs")
    }
  }

  path <- file.path(root, ...)
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)

  if (isTRUE(create) && !dir.exists(path)) {
    if (!dir.create(path, recursive = TRUE, showWarnings = FALSE) &&
        !dir.exists(path)) {
      stop("Could not create MetaOrdDTA output directory: ", path)
    }
  }

  path
}


find_metaorddta_project_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = FALSE)

  repeat {
    description <- file.path(current, "DESCRIPTION")
    if (file.exists(description)) {
      fields <- tryCatch(
        read.dcf(description, fields = "Package"),
        error = function(e) NULL
      )
      if (!is.null(fields) && identical(unname(fields[1L]), "MetaOrdDTA")) {
        return(current)
      }
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      return(NULL)
    }
    current <- parent
  }
}
