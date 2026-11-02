# =============================================================================
# simulate_world.R
#
# Re-creates the simulated Diversity Outbred (DO) mice used in the workshop:
# data/cohort1.rds, data/cohortA.rds and data/cohortB.rds.
#
# The truth built in:
#   - X raises Y:              1 unit of X adds 0.5 units of Y   (X -> Y = 0.5)
#   - Y has no effect on X:                                      (Y -> X = 0)
#   - a hidden factor U raises both X and Y, so the plain association of Y
#     with X is about 1, twice the causal effect.
#
# Who does what:
#   simcross        breeds the DO population and simulates meiosis, giving each
#                   mouse its mosaic of the eight founder genomes;
#   qtl2founderSim  builds the traits from that mosaic: founder effects at the
#                   QTL, the hidden U, the arrow X -> Y, noise, and which mice
#                   have X and/or Y measured;
#   qtl2 / qtl2mr   kinship and genome scans for each cohort, and the design
#                   checks.
#
# Run from the project folder (it takes a few minutes):
#   Rscript data/simulation/simulate_world.R
#
# See README.md in this folder for the reasoning behind every number.
# =============================================================================


# ---- 0. Packages and helper functions ---------------------------------------

suppressPackageStartupMessages({
  library(simcross)        # DO pedigree and meiosis
  library(qtl2founderSim)  # traits from founder genotypes
  library(qtl2)            # kinship and genome scans
  library(qtl2mr)          # instrument selection for the design checks
})

source(file.path("data", "simulation", "simulation_functions.R"))

check_package_versions(c(simcross = "0.10", qtl2founderSim = "0.0.0.9009",
                         qtl2mr = "0.6.0"))

# Two optional hooks, used only to test this script; a normal run ignores them.
#   output_dir          - if set before running, files go there instead of
#                         the project folder (default ".", so they land in data/);
#   quick_test_settings - if set, a list that shrinks SETTINGS for a fast test
#                         run (applied right after the settings, step 1).
if (!exists("output_dir")) output_dir <- "."


# ---- 1. Settings (all fixed before any data were made) ----------------------

SETTINGS <- list(

  # Seeds: attempt k uses seed first_seed + k (see step 10).
  first_seed   = 202611025,
  max_attempts = 5,

  # Genome: 1,500 markers on the 19 autosomes.
  founders  = LETTERS[1:8],
  chr       = as.character(1:19),
  n_markers = 1500,

  # DO breeding (simcross): 20 generations, 144 breeding pairs,
  # 12 pups per pair in the last generation (1,728 mice to choose from).
  ngen   = 20,
  npairs = 144,
  nkids  = 12,
  ccgen  = NULL,     # NULL = simcross's default pre-CC inbreeding generations

  # Cohorts: three groups of different mice, sizes chosen in advance.
  cohort_sizes = c(cohort1 = 400,    # X and Y measured
                   cohortA = 500,    # only X measured
                   cohortB = 700),   # only Y measured

  # The causal model.
  theta   = 0.5,                 # X -> Y
  U_to_X  = sqrt(0.10),          # a: U explains 10% of the variance of X
  U_to_Y  = 0.5 / sqrt(0.10),    # b: chosen so that a * b = 0.5
  qtl_var = c(X = 0.08,          # each of the 10 X-QTL explains 0.08
              Y = 0.45),         # each of the 4 Y-QTL adds 0.45 to Y

  # Polygenic background (shared through kinship) and noise, per trait.
  polygenic_var = c(U = 0, X = 0.05, Y = 0),
  noise_var     = c(U = 1, X = 0.05, Y = 0.25),

  # Instrument rule used by the design checks (the same rule as the tasks).
  instrument_rule = list(lod_min = 3, F_min = 4, x_lod_max = 3, window_cM = 10),

  # Design checks (see run_design_checks() for what each one means).
  checks = list(freq = c(0.09, 0.16), min_copies = 20, x_recovered = 8,
                slope = c(0.80, 1.20), y_recovered = 3),

  kinship_digits = 6
)

# Test hook (see step 0): a quick test run may shrink these settings.
if (exists("quick_test_settings")) SETTINGS <- modifyList(SETTINGS, quick_test_settings)

TRAITS <- c("U", "X", "Y")


# ---- 2. The QTL: where they are and their founder patterns ------------------

# 10 QTL for X and 4 QTL that act on Y directly, each on its own chromosome.
# `frac` is the position along the chromosome (0 = start, 1 = end).
QTL <- data.frame(
  id    = c(sprintf("QX%d", 1:10), sprintf("QY%d", 1:4)),
  trait = c(rep("X", 10), rep("Y", 4)),
  chr   = as.character(c(1, 2, 4, 5, 7, 9, 11, 13, 15, 17, 3, 8, 12, 18)),
  frac  = c(.35, .60, .25, .70, .45, .55, .30, .65, .40, .50, .40, .60, .35, .55),
  type  = c("biallelic 1:7", "biallelic 2:6", "three-level", "biallelic 3:5",
            "eight-level", "biallelic 4:4", "two private alleles", "five-level",
            "biallelic 3:5", "irregular", "biallelic 1:7", "three-level",
            "biallelic 4:4", "four-level"),
  stringsAsFactors = FALSE
)

# How the eight founders (A-H) differ at each QTL. Five X-QTL split the
# founders into two groups, like a SNP ("biallelic"); the other five have three
# or more levels, which only a founder-based analysis can see in full.
# sdp_effect(c("A", "D")) means "founders A and D differ from the rest".
PATTERNS <- list(
  QX1  = sdp_effect("F", founders = SETTINGS$founders),
  QX2  = sdp_effect(c("A", "D"), founders = SETTINGS$founders),
  QX3  = c(1, 0, -1, 1, 0, -1, 1, 0),
  QX4  = sdp_effect(c("B", "E", "G"), founders = SETTINGS$founders),
  QX5  = c(1.5, -0.5, 3.5, -2.5, 0.5, -3.5, 2.5, -1.5),
  QX6  = sdp_effect(c("A", "B", "C", "D"), founders = SETTINGS$founders),
  QX7  = c(0, 0, 2, 0, 0, 0, -2, 0),
  QX8  = c(-1, 2, -1, 0, 1, 0, -2, 1),
  QX9  = sdp_effect(c("C", "F", "H"), founders = SETTINGS$founders),
  QX10 = c(1.2, -0.3, -0.8, 0.2, -1.1, 0.4, 0.9, -0.5),
  QY1  = sdp_effect("E", founders = SETTINGS$founders),
  QY2  = c(1, -1, 0, 0, 1, -1, 0, 0),
  QY3  = sdp_effect(c("B", "C", "G", "H"), founders = SETTINGS$founders),
  QY4  = c(-2, 0, 1, 1, 0, 0, -1, 1)
)

PATTERNS <- lapply(PATTERNS, centre_pattern, founders = SETTINGS$founders)

# The patterns must be clearly different from each other.
pattern_cor <- abs(cor(do.call(cbind, PATTERNS)))
diag(pattern_cor) <- 0
stopifnot(max(pattern_cor) < 0.80)


# ---- 3. Marker map ----------------------------------------------------------

# Real mouse chromosome lengths in cM (Cox et al. map, shipped with simcross).
CHR_LENGTHS <- simcross::mouseL_cox[SETTINGS$chr]

MAP <- make_marker_map(CHR_LENGTHS, SETTINGS$n_markers)

QTL <- place_qtl(QTL, MAP, CHR_LENGTHS)

stopifnot(nrow(MAP$table) == SETTINGS$n_markers, !anyDuplicated(QTL$chr))


# ---- Steps 4-9: one simulated world from one seed ---------------------------

simulate_world_once <- function(attempt) {

  seed <- SETTINGS$first_seed + attempt
  cat(sprintf("\n== Attempt %d (seed %d)\n", attempt, seed))

  # qtl2founderSim's random-number streams for this seed (traits, choice of
  # mice, cohorts). simcross uses a separate stream (set in step 4).
  streams <- rng_streams(seed = seed)


  # ---- 4. The DO population (simcross) ----

  set_genotype_stream(seed)

  population <- simulate_do_population(SETTINGS, CHR_LENGTHS, MAP)
  probs      <- population$probs

  cat(sprintf("Bred %d generations; %d mice in the last generation.\n",
              SETTINGS$ngen, dim(probs[[1]])[1]))


  # ---- 5. Choose the mice of the study population (qtl2founderSim) ----

  # Kinship of all last-generation mice, then a random draw of 1,600 of them.
  kinship_all <- qtl2::calc_kinship(probs, type = "overall", cores = 1)

  n_population <- sum(SETTINGS$cohort_sizes)

  # The hand-off from simcross to qtl2founderSim: the founder probabilities of
  # step 4 and the kinship computed from them.
  #   kinship_scale = "qtl2"  the kinship is on qtl2's scale (shared founder
  #                           alleles; a mouse with itself is about 0.57 on
  #                           average). qtl2founderSim doubles it to the usual
  #                           relationship scale (self about 1.14) and draws the
  #                           polygenic part of X with covariance 0.05 times this
  #                           matrix. Each mouse's polygenic value then has
  #                           variance about 0.05 x 1.14 = 0.057; because the
  #                           mice also share alleles with each other, the
  #                           variance across these mice is about 0.043.
  #   provenance = "overall"  a label: genome-wide kinship, not leave-one-out.
  #   n, stream               draw n mice with the design random-number stream.
  #   markers                 the QTL markers, where founder effects will act.
  #   map                     marker positions, only recorded with the design.
  design <- make_empirical_design(
    genoprobs = unclass(probs)[SETTINGS$chr], K_gen = kinship_all,
    kinship_scale = "qtl2", provenance = "overall",
    n = n_population, stream = streams,
    markers = QTL$marker, map = MAP$by_chr
  )

  population_ids <- design$ids


  # ---- 6. Founder effects at the QTL ----

  # Scale each pattern so the QTL explains its target variance in these mice.
  effects <- make_founder_effects(probs, population_ids, QTL, PATTERNS,
                                  SETTINGS$qtl_var, TRAITS)


  # ---- 7. The causal model (qtl2founderSim) ----

  # Arrows between traits: B[to, from]. U -> X, U -> Y and X -> Y.
  # There is no arrow from Y to X.
  B <- matrix(0, 3, 3, dimnames = list(TRAITS, TRAITS))
  B["X", "U"] <- SETTINGS$U_to_X
  B["Y", "U"] <- SETTINGS$U_to_Y
  B["Y", "X"] <- SETTINGS$theta

  polygenic <- diag(SETTINGS$polygenic_var)
  noise     <- diag(SETTINGS$noise_var)
  dimnames(polygenic) <- dimnames(noise) <- list(TRAITS, TRAITS)

  # Who is measured: Cohort 1 has X and Y, Cohort A only X, Cohort B only Y.
  sizes <- SETTINGS$cohort_sizes
  measurement_plan <- two_sample(
    n_x = sizes[["cohort1"]] + sizes[["cohortA"]],
    n_y = sizes[["cohort1"]] + sizes[["cohortB"]],
    overlap = sizes[["cohort1"]],
    x = "X", y = "Y", method = "random"
  )

  model <- causal_spec(
    design, TRAITS, B = B,
    roles = c(U = "other", X = "exposure", Y = "outcome"),
    latent = "U",                    # U is simulated but never measured
    f = effects$for_model,           # founder effects at the QTL
    Sigma_g = polygenic, Sigma_e = noise,
    membership = measurement_plan
  )


  # ---- 8. Simulate the traits ----

  # All 1,600 mice get U, X and Y from one model; then each trait is kept
  # only for the mice in which it is "measured".
  sim <- simulate_multitrait(design, model, stream = streams)

  # The total causal effects must be exactly the truth we built in.
  stopifnot(abs(sim$truth$total_effects["Y", "X"] - SETTINGS$theta) < 1e-12,
            sim$truth$total_effects["X", "Y"] == 0)


  # ---- 9. Split into the three cohorts ----

  cohort_ids <- split_cohorts(sim$truth$membership$sets)

  stopifnot(lengths(cohort_ids) == sizes[names(cohort_ids)],
            !anyDuplicated(unlist(cohort_ids)),
            setequal(unlist(cohort_ids), population_ids))

  measured <- list(cohort1 = c("X", "Y"), cohortA = "X", cohortB = "Y")

  cohorts <- lapply(names(cohort_ids), function(name) {
    make_cohort_file(probs, cohort_ids[[name]], sim$observed$y, measured[[name]],
                     MAP$table, SETTINGS$kinship_digits)
  })
  names(cohorts) <- names(cohort_ids)


  # ---- Design checks (no MR estimate is computed) ----

  checks <- run_design_checks(probs, population_ids, cohorts, QTL, SETTINGS,
                              attempt, seed)

  cat(sprintf("Checks: C1 %s, C2 %s, C3 %s, C4 %s\n",
              checks$C1, checks$C2, checks$C3, checks$C4))

  # Everything about this world. The cohorts and founder effects are saved in
  # step 11. `details` (pedigree, hidden U, the model object) is not saved by
  # this script; it is kept in memory for anyone who wants to inspect it.
  list(
    checks = checks, cohorts = cohorts, founder_effects = effects$table,
    details = list(seed = seed, pedigree = population$pedigree,
                   population_ids = population_ids, cohort_ids = cohort_ids,
                   model = model, truth = sim$truth,
                   instruments = attr(checks, "instruments"))
  )
}


# ---- 10. Attempts: keep the first world that passes all checks ---------------

# If a world fails a design check, it is discarded and the next seed is tried.
# The settings never change.
attempt_log <- NULL
world <- NULL

for (attempt in seq_len(SETTINGS$max_attempts) - 1) {

  candidate   <- simulate_world_once(attempt)
  attempt_log <- rbind(attempt_log, candidate$checks)

  if (isTRUE(candidate$checks$pass)) {
    world <- candidate
    break
  }
}

if (is.null(world)) stop("No attempt passed the design checks.", call. = FALSE)

cat(sprintf("\nAccepted attempt %d (seed %d).\n",
            world$checks$attempt, world$checks$seed))


# ---- 11. Save -------------------------------------------------------------------

data_dir  <- file.path(output_dir, "data")
truth_dir <- file.path(output_dir, "data", "simulation", "truth")
dir.create(truth_dir, recursive = TRUE, showWarnings = FALSE)

# The three cohort files used by the tasks.
for (name in names(world$cohorts)) {
  saveRDS(world$cohorts[[name]], file.path(data_dir, paste0(name, ".rds")),
          compress = "xz")
}

# The truth, as plain tables.
write.csv(QTL, file.path(truth_dir, "true_qtl.csv"), row.names = FALSE)

write.csv(world$founder_effects, file.path(truth_dir, "true_founder_effects.csv"),
          row.names = FALSE)

write.csv(attempt_log, file.path(truth_dir, "attempt_log.csv"), row.names = FALSE)

cat("Saved data/cohort1.rds, data/cohortA.rds, data/cohortB.rds and",
    "data/simulation/truth/.\n")
