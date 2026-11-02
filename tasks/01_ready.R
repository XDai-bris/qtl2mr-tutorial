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
