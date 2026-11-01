# PART A | One QTL: eight founder effects, or one SNP?
cat("\nTASK 05A | A QTL is an 8-founder pattern\n")
locus <- blocks[[match("chr11_027", markers)]]
founder_x <- founder_effects(locus$fit_x)
show_table(founder_x)
snp <- best_snp(locus$fit_x)
cat("Best single SNP: founders", paste(snp, collapse = ""), "against the rest\n")
plot_founder_vs_snp(founder_x, snp, snp_effect(locus$fit_x, snp), locus$marker)

# PART B | Extension: founder-based vs SNP-based MR on the same mice
cat("\nTASK 05B | Founder-based vs SNP-based two-sample MR\n")
snp_table <- do.call(rbind, lapply(blocks, function(b) {
  s <- best_snp(b$fit_x)  # picked with X in Cohort A, as a SNP study would
  ex <- snp_effect(b$fit_x, s)
  ey <- snp_effect(b$fit_y, s)
  founder_signal <- drop(t(b$beta_x) %*% solve(b$Sigma_x) %*% b$beta_x)
  data.frame(marker = b$marker, chr = b$chr, snp = paste(s, collapse = ""),
             bx = ex[["effect"]], se_x = ex[["se"]], by = ey[["effect"]], se_y = ey[["se"]],
             share = (ex[["effect"]] / ex[["se"]])^2 / founder_signal)
}))
show_table(snp_table[, c("marker", "snp", "share")])
mr_snp <- second_order_gls_ivw(snp_table$bx, snp_table$by,
                               diag(snp_table$se_x^2), diag(snp_table$se_y^2))
comparison <- data.frame(
  method = c("Founder", "SNP"),
  estimate = c(mr_2s$theta, mr_snp$theta),
  lower = c(mr_2s$profile_ci_lower, mr_snp$profile_ci_lower),
  upper = c(mr_2s$profile_ci_upper, mr_snp$profile_ci_upper)
)
comparison$width <- comparison$upper - comparison$lower
show_table(comparison)
plot_founder_vs_snp_precision(comparison, snp_table)
cat("Milestone 5: founder effects keep the whole pattern; a SNP keeps one split.\n")
