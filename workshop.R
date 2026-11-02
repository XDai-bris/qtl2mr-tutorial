# All four workshop tasks in one script (the same code as tasks/01-04).
# Open qtl2mr-tutorial.Rproj, start a fresh R session, then run from the top.

# Task 01: get ready.
# Run this at home before the workshop, and again at the start of class.

# Install the bundled qtl2mr 0.6.0 into this project's own library/ folder.
# Nothing outside the project folder is changed.
source("R/install.R")

# Load the small workshop helpers (tables and figures) and start a new
# results folder, results/run_..., where every figure is saved.
source("R/setup.R")

# List the simulated data files that the next tasks will load.
cat("Data files:", paste(list.files("data", pattern = "[.]rds$"), collapse = ", "), "\n")

cat("Milestone 1: qtl2mr 0.6.0 is loaded and a results folder is ready.\n")

# Task 02: meet the simulated world. Load the three cohorts and look at them.
cat("\nTASK 02 | Meet the simulated world\n")

# ---- Load the three cohorts ------------------------------------------------
# Three groups of mice drawn from ONE simulated Diversity Outbred population.
# No mouse is in two groups; they differ in size and in which trait was measured.

# Cohort 1: 400 mice, both X and Y measured.
cohort1 <- readRDS("data/cohort1.rds")

# Cohort A: 500 mice, only X measured.
cohortA <- readRDS("data/cohortA.rds")

# Cohort B: 700 mice, only Y measured.
cohortB <- readRDS("data/cohortB.rds")

# A small table: how many mice each cohort has, and what was measured.
show_table(data.frame(cohort = c("Cohort 1", "Cohort A", "Cohort B"),
                      mice = c(nrow(cohort1$pheno), nrow(cohortA$pheno), nrow(cohortB$pheno)),
                      measured = c("X and Y", "X only", "Y only")))

# Each cohort file is a list holding the usual qtl2 objects.
cat("Each cohort file holds:", paste(names(cohort1), collapse = ", "), "\n")

# ---- The traits --------------------------------------------------------------
# The truth built into the simulation:
# X raises Y by 0.5; Y has no effect on X; a hidden factor U raises both.

# The first five mice of Cohort 1: one row per mouse, X and Y in simulated units.
show_table(data.frame(mouse = rownames(cohort1$pheno)[1:5], head(cohort1$pheno, 5)))

# ---- The genotypes -----------------------------------------------------------
# The marker map: how many markers, on how many chromosomes.
cat("Markers:", nrow(cohort1$map), "on", length(unique(cohort1$map$chr)), "chromosomes\n")

# The first three markers, with their positions in centimorgans (cM).
show_table(head(cohort1$map, 3))

# Founder mosaic: which founder each stretch of chromosome 2 comes from,
# in six mice; the dashed line marks chr2_067, an instrument in Tasks 03 and 04.
plot_founder_mosaic(cohort1, chr = "2", marker = "chr2_067")

# Kinship: how related the first 60 mice of Cohort 1 are to each other.
plot_kinship(cohort1)

cat("Milestone 2: founder probabilities, kinship, traits and scans are loaded.\n")

# PART A | Correlation: the observational slope
cat("\nTASK 03A | Cohort 1: how much does Y change with X?\n")

# Fit a straight line of Y against X in Cohort 1 (an ordinary regression).
obs <- lm(Y ~ X, data = cohort1$pheno)

# The 95% interval of the slope.
obs_ci <- confint(obs)["X", ]

# Show the slope and its interval.
show_table(data.frame(slope = coef(obs)[["X"]], lower = obs_ci[[1]], upper = obs_ci[[2]]))

cat("The true causal effect is 0.5. The hidden U also moves X and Y together.\n")

# PART B | MR from X to Y
cat("\nTASK 03B | MR from X to Y in Cohort 1\n")

# Step 1. Collect Cohort 1's qtl2 inputs in one object: founder probabilities,
# kinship (one matrix per left-out chromosome), traits and the marker map.
resources1 <- qtl2mr_prepare_resources(
  allele_probs = cohort1$genoprobs, kinship = cohort1$kinship,
  pheno = cohort1$pheno, covar = NULL, marker_map = cohort1$map
)

# Step 2. Choose instruments for X from the qtl2 scan of X: the top marker
# of each region with LOD >= 3, kept only if its F statistic is at least 4.
instruments_x <- qtl2mr_select_instruments(
  resources1, cohort1$scan[, "X", drop = FALSE], exposures = "X",
  lod_thresh = 3, F_X_min = 4, progress = FALSE
)

# Show the instruments that passed both rules (F_X_block is the F statistic).
show_table(instruments_x[instruments_x$selected, c("marker", "chr", "pos", "lod", "F_X_block")])

# Step 3. Run MR from X to Y with these instruments.
# standardize_phenotypes = FALSE keeps X and Y in their own units, so the
# answer can be compared with the truth, 0.5. The slower optional analyses
# (bootstrap, leave-one-out, sensitivity, reports) are switched off.
mr_xy <- run_qtl2mr_analysis(
  qtl2mr_pair_input(resources1, instruments_x, exposure = "X", outcome = "Y"),
  output_dir = file.path(run_dir, "mr_x_to_y"), standardize_phenotypes = FALSE,
  heterogeneity_method = "chisq", diagnostic_bootstrap_B = 0L,
  run_loo_mr = FALSE, run_instrument_sensitivity = FALSE,
  make_html = FALSE, make_figures = FALSE, verbose = FALSE, progress = FALSE
)

# Show the MR estimate and its 95% interval.
show_table(summary(mr_xy)[, c("exposure", "outcome", "n_loci", "estimate", "ci_lower", "ci_upper")])

# PART C | MR from Y to X, and the three answers side by side
cat("\nTASK 03C | MR from Y to X in Cohort 1\n")

# Step 1. Candidate instruments for Y, from the qtl2 scan of Y, with the
# same two rules (LOD >= 3, F statistic >= 4).
instruments_y <- qtl2mr_select_instruments(
  resources1, cohort1$scan[, "Y", drop = FALSE], exposures = "Y",
  lod_thresh = 3, F_X_min = 4, progress = FALSE
)

# Step 2. Keep only loci linked to Y but not to X:
# look up the X LOD at the same marker, and drop the marker if it is 3 or more.
instruments_y$x_lod <- cohort1$scan[instruments_y$marker, "X"]
instruments_y$selected <- instruments_y$selected & instruments_y$x_lod < 3

# Show the Y instruments that remain, with their Y and X evidence.
keep_y <- instruments_y[instruments_y$selected, ]
show_table(data.frame(marker = keep_y$marker, chr = keep_y$chr, lod_Y = keep_y$lod,
                      F_Y = keep_y$F_X_block, lod_X = keep_y$x_lod))

# Step 3. Run MR from Y to X, with the same settings as in Part B.
mr_yx <- run_qtl2mr_analysis(
  qtl2mr_pair_input(resources1, instruments_y, exposure = "Y", outcome = "X"),
  output_dir = file.path(run_dir, "mr_y_to_x"), standardize_phenotypes = FALSE,
  heterogeneity_method = "chisq", diagnostic_bootstrap_B = 0L,
  run_loo_mr = FALSE, run_instrument_sensitivity = FALSE,
  make_html = FALSE, make_figures = FALSE, verbose = FALSE, progress = FALSE
)

# Step 4. Put the three answers side by side: the association and the two
# MR estimates, each with its 95% interval and the truth.
s_xy <- summary(mr_xy)
s_yx <- summary(mr_yx)
answers <- data.frame(
  question = c("Observational (Y on X)", "MR X -> Y", "MR Y -> X"),
  estimate = c(coef(obs)[["X"]], s_xy$estimate, s_yx$estimate),
  lower = c(obs_ci[[1]], s_xy$ci_lower, s_yx$ci_lower),
  upper = c(obs_ci[[2]], s_xy$ci_upper, s_yx$ci_upper),
  truth = c(NA, 0.5, 0)
)
show_table(answers)

# Draw the three answers in one figure, with the truth as dashed lines.
answers$label <- c("Observational\n(Y on X)", paste("MR X", arrow, "Y"), paste("MR Y", arrow, "X"))
plot_summary(answers)

cat("Milestone 3: the association (slope well above 0.5) is not the causal effect (0.5); MR recovers the arrow.\n")

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

