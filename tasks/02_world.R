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
