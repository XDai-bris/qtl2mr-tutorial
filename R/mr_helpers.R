# Small helpers used in Task 04. The qtl2mr calls stay visible in tasks/.

# Two-sample locus block: X effects fitted in one cohort, Y effects in another.
# No mouse is in both cohorts, so the X-Y covariance of the estimates is zero.
two_sample_block <- function(fit_x, fit_y) {
  stopifnot(identical(fit_x$marker, fit_y$marker),
            !length(intersect(fit_x$sample_ids, fit_y$sample_ids)))
  k <- length(fit_x$beta_fc)
  list(marker = fit_x$marker, chr = fit_x$chr,
       beta_x = fit_x$beta_fc, beta_y = fit_y$beta_fc,
       Sigma_x = fit_x$vcov_fc, Sigma_y = fit_y$vcov_fc,
       Sigma_xy = matrix(0, k, k, dimnames = list(names(fit_x$beta_fc), names(fit_y$beta_fc))),
       n_overlap = 0L, fit_x = fit_x, fit_y = fit_y)
}

# One point per locus for the MR scatter: the X founder pattern (from the
# exposure cohort) and how much of that same pattern shows up in Y.
locus_point <- function(block) {
  u <- block$beta_x / sqrt(sum(block$beta_x^2))
  data.frame(marker = block$marker,
             x = sum(u * block$beta_x), se_x = sqrt(drop(t(u) %*% block$Sigma_x %*% u)),
             y = sum(u * block$beta_y), se_y = sqrt(drop(t(u) %*% block$Sigma_y %*% u)))
}
