# Check that workshop.R holds exactly the code of tasks/01-04: after its
# three header lines, each task file in order, each followed by one blank line.
#
# Usage, from the repository root:
#   Rscript .github/scripts/check_workshop_script.R      (exit status 1 if they differ)
local({
  tasks <- sort(list.files("tasks", pattern = "\\.R$", full.names = TRUE))
  want <- unlist(lapply(tasks, function(f) c(readLines(f, warn = FALSE), "")))
  have <- readLines("workshop.R", warn = FALSE)[-(1:3)]
  if (!identical(have, want)) {
    n <- min(length(have), length(want))
    first <- which(have[seq_len(n)] != want[seq_len(n)])[1]
    cat("FAIL: workshop.R differs from tasks/ (first difference at line ",
        if (is.na(first)) n + 1 else first, " after the header).\n", sep = "")
    quit(status = 1)
  }
  cat("PASS: workshop.R matches tasks/ (", length(tasks), " task files).\n", sep = "")
})
