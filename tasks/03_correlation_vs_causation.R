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

# PART C | MR from Y to X, and the answers in both directions
cat("\nTASK 03C | MR from Y to X in Cohort 1\n")

# Step 1. The association in this direction: fit X against Y.
obs_rev <- lm(X ~ Y, data = cohort1$pheno)
obs_rev_ci <- confint(obs_rev)["Y", ]

# Step 2. Candidate instruments for Y, from the qtl2 scan of Y, with the
# same two rules (LOD >= 3, F statistic >= 4).
instruments_y <- qtl2mr_select_instruments(
  resources1, cohort1$scan[, "Y", drop = FALSE], exposures = "Y",
  lod_thresh = 3, F_X_min = 4, progress = FALSE
)

# Step 3. Keep only loci linked to Y but not to X:
# look up the X LOD at the same marker, and drop the marker if it is 3 or more.
instruments_y$x_lod <- cohort1$scan[instruments_y$marker, "X"]
instruments_y$selected <- instruments_y$selected & instruments_y$x_lod < 3

# Show the Y instruments that remain, with their Y and X evidence.
keep_y <- instruments_y[instruments_y$selected, ]
show_table(data.frame(marker = keep_y$marker, chr = keep_y$chr, lod_Y = keep_y$lod,
                      F_Y = keep_y$F_X_block, lod_X = keep_y$x_lod))

# Step 4. Run MR from Y to X, with the same settings as in Part B.
mr_yx <- run_qtl2mr_analysis(
  qtl2mr_pair_input(resources1, instruments_y, exposure = "Y", outcome = "X"),
  output_dir = file.path(run_dir, "mr_y_to_x"), standardize_phenotypes = FALSE,
  heterogeneity_method = "chisq", diagnostic_bootstrap_B = 0L,
  run_loo_mr = FALSE, run_instrument_sensitivity = FALSE,
  make_html = FALSE, make_figures = FALSE, verbose = FALSE, progress = FALSE
)

# Step 5. Put the answers for both directions side by side: in each
# direction, the association and the MR estimate, with 95% intervals.
# The truth is the causal effect: 0.5 for X -> Y and 0 for Y -> X.
s_xy <- summary(mr_xy)
s_yx <- summary(mr_yx)
answers <- data.frame(
  question = c("Observational (Y on X)", "MR X -> Y", "Observational (X on Y)", "MR Y -> X"),
  estimate = c(coef(obs)[["X"]], s_xy$estimate, coef(obs_rev)[["Y"]], s_yx$estimate),
  lower = c(obs_ci[[1]], s_xy$ci_lower, obs_rev_ci[[1]], s_yx$ci_lower),
  upper = c(obs_ci[[2]], s_xy$ci_upper, obs_rev_ci[[2]], s_yx$ci_upper),
  truth = c(NA, 0.5, NA, 0)
)
show_table(answers)

# Draw the four answers in one figure, grouped by direction, with the truth
# of each direction as a dashed line.
answers$label <- c("Observational\n(Y on X)", paste("MR X", arrow, "Y"),
                   "Observational\n(X on Y)", paste("MR Y", arrow, "X"))
plot_summary(answers)

cat("Milestone 3: the association is not the causal effect and shows no direction; MR recovers the arrow.\n")
