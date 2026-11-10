##
## ---- Read the real-data fits of the sim runner and print the DGM values of each test in the sim data function's form -----------
##
## Usage: Rscript derive_the_DGM_values_from_the_real_data_fits.R <folder of GAD_2> <folder of HADS> <folder of BAI> <out file>
## Each folder holds the results lists of the runner's fits to the REAL data of one test: ord-biv fixed-C, ord-biv random-C
## pop-median and R&G random-C pop-median. The DGM values are the posterior medians (the "50%" column of the summary):
## ord_biv_random_C: beta_mu, beta_SD, beta_corr, C_nd / C_d (= C_mu of each class), alpha_nd / alpha_d (= 0.01 + phi * kappa,
## the model's own line), C_nd_SD_prob_scale / C_d_SD_prob_scale (compute_dirichlet_SDs of the alpha, the function of the
## sim data file, the first n_thr of the n_cat category SDs, as the old values); ord_biv_fixed_C: beta_mu, beta_SD, beta_corr, C_nd / C_d (the cutpoints from C_raw_vec: C[1] = C_raw[1],
## C[k] = C[k - 1] + exp(C_raw[k]), the fixed-cutpoint model's own transform with softplus = 0); R&G_random_C: beta_mu,
## beta_SD, raw_scale_mu, raw_scale_SD, C (= C_mu), alpha, C_SD_prob_scale.
args <- commandArgs(trailingOnly = TRUE)
folders <- args[1:3]
names(folders) <- c("GAD_2", "HADS", "BAI")
out_file <- args[4]
##
## ---- compute_dirichlet_SDs, copied from R_fn_sim_data_ord_MA_v2.R ------------------------------------------------------------------
##
compute_dirichlet_SDs <- function(alpha_vec) {

  alpha_0 <- sum(alpha_vec)
  n_cat <- length(alpha_vec)

  # SD for each category probability
  sds <- sqrt((alpha_vec * (alpha_0 - alpha_vec)) / (alpha_0^2 * (alpha_0 + 1)))

  return(sds)

}
##
## ---- Helpers --------------------------------------------------------------------------------------------------------------------
##
median_of <- function(t, name) {
      rows <- t[t$parameter == name | grepl(paste0("^", gsub("([\\[\\]])", "\\\\\\1", name), "\\["), t$parameter), ]
      if (nrow(rows) == 0) stop(paste("no rows for", name))
      ## order by the index inside the brackets so that C_mu[1,2] follows C_mu[1,1]
      idx <- sub("^[^\\[]*\\[?", "", rows$parameter)
      idx <- sub("\\]$", "", idx)
      if (any(nzchar(idx))) {
            key <- do.call(rbind, lapply(strsplit(idx, ","), function(v) as.integer(v)))
            rows <- rows[do.call(order, as.data.frame(key)), ]
      }
      as.numeric(rows[["50%"]])
}
fmt <- function(v) paste0("c(", paste(signif(v, 6), collapse = ", "), ")")
fmt1 <- function(v) as.character(signif(v, 4))
find_fit <- function(folder, model, random_thresholds) {
      f <- list.files(folder, pattern = "_results_list[.]RDS$", full.names = TRUE)
      f <- f[grepl(paste0("param_", model), f, fixed = TRUE) &
             grepl(paste0("rand_thr_", if (random_thresholds) "TRUE" else "FALSE"), f, fixed = TRUE)]
      if (length(f) == 0) return(NULL)
      if (length(f) != 1) stop(paste("expected one fit, found", length(f), "for", model, random_thresholds, "in", folder))
      readRDS(f)$tibble_all
}
## (the Jones model's cutpoints: fn_Stan_box_cox of the thresholds, then C = -2 + 4 (C_raw - min) / (max - min))
box_cox_cutpoints_scaled <- function(x, lambda) {
      C_raw <- if (lambda == 0) log(x) else (x^lambda - 1) / lambda
      -2 + 4 * (C_raw - min(C_raw)) / (max(C_raw) - min(C_raw))
}
cutpoints_from_C_raw <- function(C_raw) {
      C <- C_raw
      if (length(C_raw) > 1) for (k in 2:length(C_raw)) C[k] <- C[k - 1] + exp(C_raw[k])
      C
}
lines_out <- character(0)
emit <- function(...) lines_out <<- c(lines_out, paste0(...))
##
## ---- One test at a time ----------------------------------------------------------------------------------------------------------
##
for (test in names(folders)) {
      L <- paste0("true_mean_", test, "_params_list")
      ## ---- ord-biv random-C (pop-median):
      t <- find_fit(folders[[test]], "ord_bivariate", TRUE)
      if (is.null(t)) emit("## ---- ", test, ", DGM ord_biv_random_C: NO FIT SAVED") else {
      n_thr <- length(median_of(t, "C_mu")) / 2
      C_mu <- matrix(median_of(t, "C_mu"), nrow = 2, byrow = TRUE)
      phi <- matrix(median_of(t, "dirichlet_phi"), nrow = 2, byrow = TRUE)
      kappa <- median_of(t, "kappa")
      alpha_nd <- 0.01 + phi[1, ] * kappa[1]
      alpha_d  <- 0.01 + phi[2, ] * kappa[2]
      emit("## ---- ", test, ", DGM ord_biv_random_C (real-data fit of the ord-biv random-cutpoint pop-median model):")
      emit(L, "$Ord_bivariate$beta_mu   <- ", fmt(median_of(t, "beta_mu")))
      emit(L, "$Ord_bivariate$beta_SD   <- ", fmt(median_of(t, "beta_SD")), "*reduction_factor")
      emit(L, "$Ord_bivariate$beta_corr <- ", fmt1(median_of(t, "beta_corr")))
      emit(L, "$Ord_bivariate$C_nd <- ", fmt(C_mu[1, ]))
      emit(L, "$Ord_bivariate$C_d  <- ", fmt(C_mu[2, ]))
      emit(L, "$Ord_bivariate$C_nd_SD_prob_scale <- ", fmt(compute_dirichlet_SDs(alpha_nd)[1:n_thr]))
      emit(L, "$Ord_bivariate$C_d_SD_prob_scale  <- ", fmt(compute_dirichlet_SDs(alpha_d)[1:n_thr]))
      emit(L, "$Ord_bivariate$alpha_nd          <- ", fmt(alpha_nd), " / (reduction_factor^2)")
      emit(L, "$Ord_bivariate$alpha_d           <- ", fmt(alpha_d), " / (reduction_factor^2)")
      emit("## kappa (median) nd / d: ", fmt(kappa), "; sum(alpha) nd / d: ", fmt(c(sum(alpha_nd), sum(alpha_d))))
      }
      ## ---- ord-biv fixed-C:
      t <- find_fit(folders[[test]], "ord_bivariate", FALSE)
      if (is.null(t)) emit("## ---- ", test, ", DGM ord_biv_fixed_C: NO FIT SAVED") else {
      C_raw <- matrix(median_of(t, "C_raw_vec"), nrow = 2, byrow = TRUE)
      emit("## ---- ", test, ", DGM ord_biv_fixed_C (real-data fit of the ord-biv fixed-cutpoint model):")
      emit(L, "$Ord_bivariate$beta_mu <- ", fmt(median_of(t, "beta_mu")))
      emit(L, "$Ord_bivariate$beta_SD <- ", fmt(median_of(t, "beta_SD")), "*reduction_factor")
      emit(L, "$Ord_bivariate$beta_corr <- ", fmt1(median_of(t, "beta_corr")))
      emit(L, "$Ord_bivariate$C_nd <- ", fmt(cutpoints_from_C_raw(C_raw[1, ])))
      emit(L, "$Ord_bivariate$C_d  <- ", fmt(cutpoints_from_C_raw(C_raw[2, ])))
      }
      ## ---- R&G random-C (pop-median):
      t <- find_fit(folders[[test]], "R&G", TRUE)
      if (is.null(t)) emit("## ---- ", test, ", DGM R&G_random_C: NO FIT SAVED") else {
      alpha <- 0.01 + median_of(t, "dirichlet_phi") * median_of(t, "kappa")
      emit("## ---- ", test, ", DGM R&G_random_C (real-data fit of the R&G random-cutpoint pop-median model):")
      emit(L, "$Ord_HSROC$beta_mu <- ", fmt1(median_of(t, "beta_mu")))
      emit(L, "$Ord_HSROC$beta_SD <-  ", fmt1(median_of(t, "beta_SD")), "*reduction_factor")
      emit(L, "$Ord_HSROC$raw_scale_mu <-  ", fmt1(median_of(t, "raw_scale_mu")))
      emit(L, "$Ord_HSROC$raw_scale_SD <-   ", fmt1(median_of(t, "raw_scale_SD")), "*reduction_factor")
      emit(L, "$Ord_HSROC$C <- ", fmt(median_of(t, "C_mu")))
      emit(L, "$Ord_HSROC$C_SD_prob_scale <-  ", fmt(compute_dirichlet_SDs(alpha)[1:(length(alpha) - 1)]))
      emit(L, "$Ord_HSROC$alpha <- ", fmt(alpha), " / (reduction_factor^2)")
      emit("## kappa (median): ", fmt1(median_of(t, "kappa")), "; sum(alpha): ", fmt1(sum(alpha)))
      }
      ## ---- Jones (the Box-Cox cutpoints from the lambda medians, as the model makes them: the Box-Cox of the thresholds
      ##      1, ..., n_thr, then scaled to [-2, 2]):
      ##      (a fit under the one-hill Se / Sp prior, 5th argument "one_hill": its parameters beta_mu,
      ##      raw_scale_mu and lambda are the prior's 3 standard Normals, so the model's baseline location and
      ##      raw scale are read from Xbeta_baseline_* and Xraw_scale_baseline_*, and its lambda from
      ##      lambda_transformed, which the summary of the fit must hold)
      t <- find_fit(folders[[test]], "Jones", FALSE)
      one_hill <- identical(args[5], "one_hill")
      if (is.null(t)) emit("## ---- ", test, ", DGM Jones: NO FIT SAVED") else
      if (one_hill && !any(grepl("^lambda_transformed\\[", t$parameter))) {
            emit("## ---- ", test, ", DGM Jones: lambda_transformed is not in the summary of this fit: ",
                 "NO DGM VALUES")
      } else {
      lambda <- median_of(t, if (one_hill) "lambda_transformed" else "lambda")
      location <- if (one_hill) c(median_of(t, "Xbeta_baseline_nd"), median_of(t, "Xbeta_baseline_d")) else
                  median_of(t, "beta_mu")
      raw_scale <- if (one_hill) c(median_of(t, "Xraw_scale_baseline_nd"),
                                   median_of(t, "Xraw_scale_baseline_d")) else
                   median_of(t, "raw_scale_mu")
      n_thr_of_the_test <- c(GAD_2 = 6, HADS = 21, BAI = 63)[[test]]
      emit("## ---- ", test, ", DGM Jones (real-data fit of the Jones model):")
      emit(L, "$Jones$beta_mu <- ", fmt(location))
      emit(L, "$Jones$beta_SD <- ", fmt(median_of(t, "beta_SD")), "*reduction_factor")
      emit(L, "$Jones$beta_corr <- ", fmt1(median_of(t, "beta_corr")))
      emit(L, "$Jones$C_nd <- ", fmt(box_cox_cutpoints_scaled(1:n_thr_of_the_test, lambda[1])))
      emit(L, "$Jones$C_d  <- ", fmt(box_cox_cutpoints_scaled(1:n_thr_of_the_test, lambda[2])))
      emit(L, "$Jones$raw_scale_mu <- ", fmt(raw_scale))
      emit(L, "$Jones$raw_scale_SD <- ", fmt(median_of(t, "raw_scale_SD")), "*reduction_factor")
      emit(L, "$Jones$raw_scale_corr <- ", fmt1(median_of(t, "raw_scale_corr")))
      emit("## lambda (median) nd / d: ", fmt(lambda))
      }
      emit("")
}
writeLines(lines_out, out_file)
cat(lines_out, sep = "\n")
