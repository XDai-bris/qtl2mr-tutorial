# Run the nine workshop tasks the way a learner does: from the repository root,
# in one fresh R session, one code block after another, in order.
# The Console text of each block is saved as <out_dir>/<task>.txt (or
# <task>_a.txt and <task>_b.txt for tasks in two parts), with the same names as
# the files in reference/console/. compare_console.R then compares them.
#
# Usage, from the repository root:
#   Rscript .github/scripts/replay.R [out_dir]      (default out_dir: replay-console)
local({
  args <- commandArgs(trailingOnly = TRUE)
  out_dir <- if (length(args)) args[[1]] else "replay-console"
  if (!file.exists("qtl2mr-tutorial.Rproj")) {
    stop("Run this script from the repository root (the folder with qtl2mr-tutorial.Rproj).")
  }
  if (requireNamespace("qtl2mr", quietly = TRUE)) {
    stop("qtl2mr is already installed outside this project; use an R library without it.")
  }
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  cat("Replay with", R.version.string, "on", R.version$platform, "\n")
  for (f in sort(list.files("tasks", pattern = "\\.R$", full.names = TRUE))) {
    txt <- readLines(f, warn = FALSE)
    # Tasks in two parts start each part with a "# PART A | ..." comment line;
    # each part is one code block in the guide.
    start <- grep("^# PART [A-Z]", txt)
    if (!length(start)) start <- 1L
    ends <- c(start[-1] - 1L, length(txt))
    for (j in seq_along(start)) {
      block <- paste0(tools::file_path_sans_ext(basename(f)),
                      if (length(start) > 1) paste0("_", letters[j]) else "")
      con <- file(file.path(out_dir, paste0(block, ".txt")), "wt")
      t0 <- proc.time()[["elapsed"]]
      sink(con, split = TRUE)
      tryCatch(eval(parse(text = txt[start[j]:ends[j]]), envir = .GlobalEnv),
               finally = {
                 sink()
                 close(con)
               })
      cat(sprintf("[replay] %s finished in %.1f s\n", block, proc.time()[["elapsed"]] - t0))
    }
  }
  # The results folder name is random; record it so the figures can be checked.
  writeLines(get("run_dir", envir = .GlobalEnv), file.path(out_dir, "RUN_DIR.txt"))
  cat("[replay] All nine tasks completed.\n")
})
