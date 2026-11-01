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
