# PART A | Instruments from Cohort A, read off in Cohort B
cat("\nTASK 04A | Find X instruments in Cohort A\n")

# Step 1. Collect Cohort A's qtl2 inputs (Cohort A measured only X).
resourcesA <- qtl2mr_prepare_resources(
  allele_probs = cohortA$genoprobs, kinship = cohortA$kinship,
  pheno = cohortA$pheno, covar = NULL, marker_map = cohortA$map
)

# Step 2. Choose instruments for X from Cohort A's scan, with the same rules
# as in Task 03: LOD >= 3, F statistic >= 4.
instruments_a <- qtl2mr_select_instruments(
  resourcesA, cohortA$scan, exposures = "X",
  lod_thresh = 3, F_X_min = 4, progress = FALSE
)

# Keep the names of the selected markers; Part B uses them in both cohorts.
markers <- instruments_a$marker[instruments_a$selected]

# Show the instruments.
show_table(instruments_a[instruments_a$selected, c("marker", "chr", "pos", "lod", "F_X_block")])

# Check: how many mice are in both Cohort A and Cohort B?
in_both <- intersect(rownames(cohortA$pheno), rownames(cohortB$pheno))
cat("Mice with both X and Y:", length(in_both), "\n")

# Miami plot: the X scan of Cohort A above, the Y scan of Cohort B below.
plot_miami(cohortA, cohortB, markers)

# PART B | Two-sample MR
cat("\nTASK 04B | X effects from Cohort A, Y effects from Cohort B\n")

# X measured in Cohort A and Y measured in Cohort B, named by mouse.
x_in_A <- setNames(cohortA$pheno$X, rownames(cohortA$pheno))
y_in_B <- setNames(cohortB$pheno$Y, rownames(cohortB$pheno))

# For each instrument: fit the founder effects on X in Cohort A, fit the
# founder effects on Y in Cohort B, then join the two fits into one block.
blocks <- lapply(markers, function(m) {
  # The chromosome of this marker, to pick the matching kinship matrix.
  chr <- cohortA$map$chr[cohortA$map$marker == m]

  # Founder effects on X, from the mice of Cohort A.
  fit_x <- fit_founder_lmm_block(cohortA$genoprobs, x_in_A,
                                 kinship = cohortA$kinship[[chr]], marker = m)

  # Founder effects on Y, from the mice of Cohort B.
  fit_y <- fit_founder_lmm_block(cohortB$genoprobs, y_in_B,
                                 kinship = cohortB$kinship[[chr]], marker = m)

  # Join the two fits; the helper records that their errors are independent.
  two_sample_block(fit_x, fit_y)
})

# Find the one slope that best links the Y patterns to the X patterns.
# allow_cross_cov = FALSE: no mouse is in both cohorts.
mr_2s <- fit_multilocus_founder_mr(blocks, allow_cross_cov = FALSE,
                                   compute_independence = FALSE, heterogeneity_method = "chisq")

# Show the estimate and its 95% interval.
show_table(data.frame(n_loci = length(blocks), estimate = mr_2s$theta,
                      lower = mr_2s$profile_ci_lower, upper = mr_2s$profile_ci_upper,
                      interval = mr_2s$confidence_set_type))

# PART C | The MR scatter: one point per instrument
cat("\nTASK 04C | Plot the instruments\n")

# One point per instrument: its effect on X (Cohort A) and the effect of the
# same founder pattern on Y (Cohort B), with standard errors.
points_2s <- do.call(rbind, lapply(blocks, locus_point))
show_table(points_2s)

# Draw the points with the MR slope and the truth.
plot_mr_scatter(points_2s, theta = mr_2s$theta)

cat("Milestone 4: no mouse has both traits, yet MR estimates the effect of X on Y.\n")
