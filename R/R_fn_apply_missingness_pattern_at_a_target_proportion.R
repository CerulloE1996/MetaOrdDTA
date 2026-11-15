#' apply_missingness_pattern_at_a_target_proportion
#'
#' Makes thresholds missing (-1) in complete simulated data so that a chosen proportion of all
#' threshold-by-study cells is missing, whilst keeping the shape of a real reporting pattern
#' (i.e., the thresholds that the real studies report most often stay the ones reported most often).
#'
#' The proportion of studies missing each threshold in the real pattern is moved up or down on the
#' odds scale by one common factor, chosen so that the mean over the thresholds equals the target.
#' The same cells are made missing in the non-diseased and the diseased data.
#'
#' @param x_complete List with 2 matrices (non-diseased, diseased) - complete simulated data;
#'   column 1 is kept as it is, columns 2 onwards are the thresholds.
#' @param x_pattern List with 2 matrices - the real data whose reporting pattern gives the shape
#'   (-1 = threshold not reported); the first matrix is used.
#' @param target_proportion_missing The proportion of threshold-by-study cells to make missing (e.g. 0.85).
#' @param min_thresholds_per_study The smallest number of thresholds every study must report (e.g. 1).
#' @param min_studies_per_threshold The smallest number of studies that must report every threshold (e.g. 1).
#' @param missing_indicator The value that marks a missing threshold (e.g. -1).
#'
#' @export
apply_missingness_pattern_at_a_target_proportion <- function( x_complete,
                                                              x_pattern,
                                                              target_proportion_missing,
                                                              min_thresholds_per_study,
                                                              min_studies_per_threshold,
                                                              missing_indicator) {

        if (!is.list(x_complete) || !is.list(x_pattern) || length(x_complete) != 2 || length(x_pattern) != 2) {
          stop("Both inputs must be lists with exactly 2 elements (non-diseased and diseased matrices)")
        }
        ##
        n_studies <- nrow(x_complete[[1]])
        pattern_matrix <- as.matrix(x_pattern[[1]])
        n_studies_pattern <- nrow(pattern_matrix)
        n_thr <- min(ncol(x_complete[[1]]), ncol(pattern_matrix)) - 1
        ##
        ## ---- The largest proportion that still leaves every study its minimum number of thresholds:
        ##
        max_proportion_missing <- 1 - min_thresholds_per_study / n_thr
        if (target_proportion_missing > max_proportion_missing + 1e-12) {
          stop(sprintf("target_proportion_missing = %.4f is not possible: with %d thresholds and at least %d reported per study the largest is %.4f",
                       target_proportion_missing, n_thr, min_thresholds_per_study, max_proportion_missing))
        }
        if (target_proportion_missing < 0) stop("target_proportion_missing must be at least 0")
        if (n_studies * n_thr * (1 - target_proportion_missing) + 1e-9 < n_thr * min_studies_per_threshold) {
          stop(sprintf("target_proportion_missing = %.4f is not possible: %d studies cannot report each of %d thresholds at least %d time(s)",
                       target_proportion_missing, n_studies, n_thr, min_studies_per_threshold))
        }
        ##
        ## ---- The shape: the proportion of the real studies missing each threshold
        ## ---- (kept half a study away from 0 and 1, so that every threshold can move):
        ##
        prop_missing_real <- colMeans(pattern_matrix[, 2:(n_thr + 1), drop = FALSE] == missing_indicator, na.rm = TRUE)
        prop_missing_real <- pmin(pmax(prop_missing_real, 0.5 / n_studies_pattern), 1 - 0.5 / n_studies_pattern)
        ##
        ## ---- One common move on the log-odds scale, so that the mean equals the target:
        ##
        prop_missing_moved <- function(shift) stats::plogis(stats::qlogis(prop_missing_real) + shift)
        if (target_proportion_missing == 0) {
          prop_missing <- rep(0, n_thr)
        } else {
          shift <- stats::uniroot(f = function(shift) mean(prop_missing_moved(shift)) - target_proportion_missing,
                                  lower = -50, upper = 50, tol = 1e-10)$root
          prop_missing <- prop_missing_moved(shift)
        }
        ##
        ## ---- The studies that miss each threshold (a whole number of studies per threshold, rounded at random
        ## ---- so that the expected number is exact):
        ##
        is_missing <- matrix(FALSE, nrow = n_studies, ncol = n_thr)
        for (k in 1:n_thr) {
          expected_n_missing <- n_studies * prop_missing[k]
          n_missing <- floor(expected_n_missing) + (stats::runif(1) < (expected_n_missing - floor(expected_n_missing)))
          if (n_missing > 0) is_missing[sample.int(n_studies, n_missing), k] <- TRUE
        }
        ##
        ## ---- Every study reports its minimum number of thresholds: a study below it is given a threshold
        ## ---- (the more often reported ones more likely), taken from a study that can spare it when there is one,
        ## ---- so that the number of studies per threshold stays the same:
        ##
        n_reported <- n_thr - rowSums(is_missing)
        for (s in which(n_reported < min_thresholds_per_study)) {
          while ((n_thr - sum(is_missing[s, ])) < min_thresholds_per_study) {
            candidates <- which(is_missing[s, ])
            weights <- pmax(1 - prop_missing[candidates], 1e-12)
            k <- if (length(candidates) == 1) candidates else sample(candidates, 1, prob = weights)
            is_missing[s, k] <- FALSE
            donors <- which(!is_missing[, k] & (n_thr - rowSums(is_missing)) > min_thresholds_per_study)
            donors <- setdiff(donors, s)
            if (length(donors) > 0) {
              donor <- if (length(donors) == 1) donors else sample(donors, 1)
              is_missing[donor, k] <- TRUE
            }
          }
        }
        ##
        ## ---- Every threshold is reported by its minimum number of studies: a threshold below it is given a study,
        ## ---- and a cell is taken from the threshold that is furthest above its expected number of studies
        ## ---- (from a study that can spare it), so that the proportion missing stays the same:
        ##
        expected_n_reporting <- n_studies * (1 - prop_missing)
        while (any((n_studies - colSums(is_missing)) < min_studies_per_threshold)) {
          n_reporting <- n_studies - colSums(is_missing)
          short <- which(n_reporting < min_studies_per_threshold)
          k <- if (length(short) == 1) short else sample(short, 1)
          candidates <- which(is_missing[, k])
          s <- if (length(candidates) == 1) candidates else sample(candidates, 1)
          is_missing[s, k] <- FALSE
          n_reported_by_study <- n_thr - rowSums(is_missing)
          n_reporting <- n_studies - colSums(is_missing)
          can_spare <- !is_missing & matrix(n_reported_by_study > min_thresholds_per_study, nrow = n_studies, ncol = n_thr) &
                       matrix(n_reporting > min_studies_per_threshold, nrow = n_studies, ncol = n_thr, byrow = TRUE)
          can_spare[s, k] <- FALSE
          if (any(can_spare)) {
            surplus <- ifelse(colSums(can_spare) > 0, n_reporting - expected_n_reporting, -Inf)
            best <- which(surplus == max(surplus))
            k_from <- if (length(best) == 1) best else sample(best, 1)
            donors <- which(can_spare[, k_from])
            donor <- if (length(donors) == 1) donors else sample(donors, 1)
            is_missing[donor, k_from] <- TRUE
          }
        }
        ##
        ## ---- The same cells in both populations:
        ##
        result <- list(non_diseased = x_complete[[1]],
                       diseased     = x_complete[[2]])
        for (pop_idx in 1:2) {
          thresholds <- result[[pop_idx]][, 2:(n_thr + 1), drop = FALSE]
          thresholds[is_missing] <- missing_indicator
          result[[pop_idx]][, 2:(n_thr + 1)] <- thresholds
        }

        return(result)

}




