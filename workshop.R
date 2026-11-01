# All five workshop tasks in one script (the same code as tasks/01-05).
# Open qtl2mr-tutorial.Rproj, start a fresh R session, then run from the top.

source("R/install.R")
source("R/setup.R")
cat("Data files:", paste(list.files("data", pattern = "[.]rds$"), collapse = ", "), "\n")
cat("Milestone 1: qtl2mr 0.6.0 is loaded and a results folder is ready.\n")

cat("\nTASK 02 | Meet the simulated world\n")
# Three groups of 400 mice from one simulated Diversity Outbred population.
cohort1 <- readRDS("data/cohort1.rds")  # X and Y measured
cohortA <- readRDS("data/cohortA.rds")  # only X measured
cohortB <- readRDS("data/cohortB.rds")  # only Y measured
show_table(data.frame(cohort = c("Cohort 1", "Cohort A", "Cohort B"),
                      mice = c(nrow(cohort1$pheno), nrow(cohortA$pheno), nrow(cohortB$pheno)),
                      measured = c("X and Y", "X only", "Y only")))
cat("Each cohort file holds:", paste(names(cohort1), collapse = ", "), "\n")
# The truth built into the simulation:
# X raises Y by 0.5; Y has no effect on X; a hidden factor U raises both.
show_table(data.frame(mouse = rownames(cohort1$pheno)[1:5], head(cohort1$pheno, 5)))
cat("Markers:", nrow(cohort1$map), "on", length(unique(cohort1$map$chr)), "chromosomes\n")
show_table(head(cohort1$map, 3))
plot_founder_bars(cohort1, marker = "chr2_050")
plot_kinship(cohort1)
cat("Milestone 2: founder probabilities, kinship, traits and scans are loaded.\n")

# PART A | Correlation: the observational slope
cat("\nTASK 03A | Cohort 1: how much does Y change with X?\n")
obs <- lm(Y ~ X, data = cohort1$pheno)
obs_ci <- confint(obs)["X", ]
show_table(data.frame(slope = coef(obs)[["X"]], lower = obs_ci[[1]], upper = obs_ci[[2]]))
cat("The true causal effect is 0.5. The hidden U also moves X and Y together.\n")

# PART B | MR from X to Y
cat("\nTASK 03B | MR from X to Y in Cohort 1\n")
resources1 <- qtl2mr_prepare_resources(
  allele_probs = cohort1$genoprobs, kinship = cohort1$kinship,
  pheno = cohort1$pheno, covar = NULL, marker_map = cohort1$map
)
# Instruments: one peak per region with LOD >= 3, kept if its founder F >= 4.
instruments_x <- qtl2mr_select_instruments(
  resources1, cohort1$scan[, "X", drop = FALSE], exposures = "X",
  lod_thresh = 3, F_X_min = 4, progress = FALSE
)
show_table(instruments_x[instruments_x$selected, c("marker", "chr", "pos", "lod", "F_X_block")])
mr_xy <- run_qtl2mr_analysis(
  qtl2mr_pair_input(resources1, instruments_x, exposure = "X", outcome = "Y"),
  output_dir = file.path(run_dir, "mr_x_to_y"), standardize_phenotypes = FALSE,
  heterogeneity_method = "chisq", diagnostic_bootstrap_B = 0L,
  run_loo_mr = FALSE, run_instrument_sensitivity = FALSE,
  make_html = FALSE, make_figures = FALSE, verbose = FALSE, progress = FALSE
)
show_table(summary(mr_xy)[, c("exposure", "outcome", "n_loci", "estimate", "ci_lower", "ci_upper")])

# PART C | MR from Y to X, and the three answers side by side
cat("\nTASK 03C | MR from Y to X in Cohort 1\n")
instruments_y <- qtl2mr_select_instruments(
  resources1, cohort1$scan[, "Y", drop = FALSE], exposures = "Y",
  lod_thresh = 3, F_X_min = 4, progress = FALSE
)
# Use loci linked to Y but not to X: X LOD below 3 at the same marker.
instruments_y$x_lod <- cohort1$scan[instruments_y$marker, "X"]
instruments_y$selected <- instruments_y$selected & instruments_y$x_lod < 3
keep_y <- instruments_y[instruments_y$selected, ]
show_table(data.frame(marker = keep_y$marker, chr = keep_y$chr, lod_Y = keep_y$lod,
                      F_Y = keep_y$F_X_block, lod_X = keep_y$x_lod))
mr_yx <- run_qtl2mr_analysis(
  qtl2mr_pair_input(resources1, instruments_y, exposure = "Y", outcome = "X"),
  output_dir = file.path(run_dir, "mr_y_to_x"), standardize_phenotypes = FALSE,
  heterogeneity_method = "chisq", diagnostic_bootstrap_B = 0L,
  run_loo_mr = FALSE, run_instrument_sensitivity = FALSE,
  make_html = FALSE, make_figures = FALSE, verbose = FALSE, progress = FALSE
)
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
answers$label <- c("Observational\n(Y on X)", paste("MR X", arrow, "Y"), paste("MR Y", arrow, "X"))
plot_summary(answers)
cat("Milestone 3: the association (slope about 1) is not the causal effect (0.5); MR recovers the arrow.\n")

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

