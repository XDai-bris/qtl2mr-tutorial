# PART A | Instruments from Cohort A, looked up in Cohort B
cat("\nTASK 04A | Find X instruments in Cohort A\n")
resourcesA <- qtl2mr_prepare_resources(
  allele_probs = cohortA$genoprobs, kinship = cohortA$kinship,
  pheno = cohortA$pheno, covar = NULL, marker_map = cohortA$map
)
instruments_a <- qtl2mr_select_instruments(
  resourcesA, cohortA$scan, exposures = "X",
  lod_thresh = 3, F_X_min = 4, progress = FALSE
)
markers <- instruments_a$marker[instruments_a$selected]
show_table(instruments_a[instruments_a$selected, c("marker", "chr", "pos", "lod", "F_X_block")])
cat("Mice with both X and Y:", length(intersect(rownames(cohortA$pheno), rownames(cohortB$pheno))), "\n")
plot_miami(cohortA, cohortB, markers)

# PART B | Two-sample MR
cat("\nTASK 04B | X effects from Cohort A, Y effects from Cohort B\n")
x_in_A <- setNames(cohortA$pheno$X, rownames(cohortA$pheno))
y_in_B <- setNames(cohortB$pheno$Y, rownames(cohortB$pheno))
blocks <- lapply(markers, function(m) {
  chr <- cohortA$map$chr[cohortA$map$marker == m]
  fit_x <- fit_founder_lmm_block(cohortA$genoprobs, x_in_A, kinship = cohortA$kinship[[chr]], marker = m)
  fit_y <- fit_founder_lmm_block(cohortB$genoprobs, y_in_B, kinship = cohortB$kinship[[chr]], marker = m)
  two_sample_block(fit_x, fit_y)
})
mr_2s <- fit_multilocus_founder_mr(blocks, allow_cross_cov = FALSE,
                                   compute_independence = FALSE, heterogeneity_method = "chisq")
show_table(data.frame(n_loci = length(blocks), estimate = mr_2s$theta,
                      lower = mr_2s$profile_ci_lower, upper = mr_2s$profile_ci_upper,
                      interval = mr_2s$confidence_set_type))

# PART C | The MR scatter: one point per instrument
cat("\nTASK 04C | Plot the instruments\n")
points_2s <- do.call(rbind, lapply(blocks, locus_point))
show_table(points_2s)
plot_mr_scatter(points_2s, theta = mr_2s$theta)
founders_2s <- do.call(rbind, lapply(blocks, function(b) {
  data.frame(marker = b$marker, founder = rownames(b$fit_x$founder_H),
             x = founder_effects(b$fit_x)$effect, y = founder_effects(b$fit_y)$effect)
}))
plot_founder_scatter(founders_2s, theta = mr_2s$theta)
cat("Milestone 4: no mouse has both traits, yet MR estimates the effect of X on Y.\n")
