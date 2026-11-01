# Small helpers used in the tasks. The qtl2mr calls stay visible in tasks/.

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

# Founder effects (centred, one per founder A-H) and their standard errors.
founder_effects <- function(fit) {
  H <- fit$founder_H
  data.frame(founder = rownames(H), effect = drop(H %*% fit$beta_fc),
             se = sqrt(diag(H %*% fit$vcov_fc %*% t(H))))
}

# One point per locus for the MR scatter: the X founder pattern (from the
# exposure cohort) and how much of that same pattern shows up in Y.
locus_point <- function(block) {
  u <- block$beta_x / sqrt(sum(block$beta_x^2))
  data.frame(marker = block$marker,
             x = sum(u * block$beta_x), se_x = sqrt(drop(t(u) %*% block$Sigma_x %*% u)),
             y = sum(u * block$beta_y), se_y = sqrt(drop(t(u) %*% block$Sigma_y %*% u)))
}

# A biallelic SNP splits the 8 founders into carriers (1) and non-carriers (0).
# Its effect is the founder fit read along that one split:
#   d = H's,  effect = d'V^-1 b / d'V^-1 d,  variance = 1 / d'V^-1 d
snp_effect <- function(fit, carriers) {
  s <- as.numeric(rownames(fit$founder_H) %in% carriers)
  d <- drop(crossprod(fit$founder_H, s))
  Vi <- solve(fit$vcov_fc)
  info <- drop(t(d) %*% Vi %*% d)
  c(effect = drop(t(d) %*% Vi %*% fit$beta_fc) / info, se = sqrt(1 / info))
}

# The SNP a SNP-based study would pick: the founder split (out of 127) with the
# strongest association with the trait in this fit.
best_snp <- function(fit) {
  founders <- rownames(fit$founder_H)
  splits <- as.matrix(expand.grid(rep(list(0:1), 7)))[-1, ]
  z <- apply(splits, 1, function(s) {
    e <- snp_effect(fit, founders[c(s, 0) == 1])
    abs(e[["effect"]] / e[["se"]])
  })
  founders[c(splits[which.max(z), ], 0) == 1]
}
