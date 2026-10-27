# Compare the Console text written by replay.R with reference/console/, and
# check that every figure in reference/figures/ was written by the replay.
#
# Rules for the Console text:
#  * The random results-folder name (results/run_...) is ignored.
#  * The line "Installing the fs dependency..." is ignored; it appears only when
#    fs is not installed yet.
#  * Each line is split into tokens at white space, so column alignment does not
#    matter. Lines must match one to one.
#  * Text tokens must match exactly.
#  * A number may differ only in its last printed decimal place:
#    |replay - reference| <= 1.5 * 10^-d, where d is the number of decimals.
#    The workshop rounds to 3 decimals and drops trailing zeros (0.880 is
#    printed as 0.88), so d is 3 for every number with a decimal point (or
#    more, if a token shows more decimals). Whole numbers (counts) must match
#    exactly. Numbers in scientific notation are below 0.001, so they are
#    compared on the same 3-decimal scale: such tiny values are floating-point
#    noise and may differ between computers.
#
# PNG bytes are not compared: drawing details differ between operating systems.
#
# Usage, from the repository root:
#   Rscript .github/scripts/compare_console.R [replay_dir] [reference_dir]
# Defaults: replay-console and reference. Exit status 1 if anything differs.
local({
  args <- commandArgs(trailingOnly = TRUE)
  replay_dir <- if (length(args) >= 1) args[[1]] else "replay-console"
  ref_dir <- if (length(args) >= 2) args[[2]] else "reference"
  problems <- character()
  note <- function(...) problems <<- c(problems, paste0(...))

  read_console <- function(path) {
    x <- readLines(path, warn = FALSE, encoding = "UTF-8")
    x <- sub("\r$", "", x)
    x <- gsub("results[\\\\/]+run_[0-9A-Fa-f]+", "results/run_<RUN>", x)
    x <- x[!grepl("^Installing the fs dependency", x)]
    # Trailing blank lines do not matter.
    while (length(x) && !nzchar(trimws(x[length(x)]))) x <- x[-length(x)]
    x
  }
  tokens <- function(line) {
    t <- strsplit(trimws(line), "[[:space:]]+")[[1]]
    t[nzchar(t)]
  }
  number_re <- "^[-+]?([0-9]+[.]?[0-9]*|[.][0-9]+)([eE][-+]?[0-9]+)?$"
  decimals <- function(tok) {
    mantissa <- sub("[eE].*$", "", tok)
    if (grepl(".", mantissa, fixed = TRUE)) nchar(sub("^[^.]*[.]", "", mantissa)) else 0L
  }
  same_token <- function(ref, new) {
    if (identical(ref, new)) return(TRUE)
    if (!grepl(number_re, ref) || !grepl(number_re, new)) return(FALSE)
    a <- as.numeric(ref)
    b <- as.numeric(new)
    scientific <- grepl("[eE]", ref) || grepl("[eE]", new)
    whole <- !scientific && decimals(ref) == 0L && decimals(new) == 0L
    if (whole) return(a == b)
    d <- if (scientific) 3L else max(3L, decimals(ref), decimals(new))
    abs(a - b) <= 1.5 * 10^-d + 1e-12
  }

  ref_files <- sort(list.files(file.path(ref_dir, "console"), pattern = "[.]txt$"))
  if (!length(ref_files)) note("No reference Console files found in ", file.path(ref_dir, "console"))
  extra <- setdiff(list.files(replay_dir, pattern = "[.]txt$"), c(ref_files, "RUN_DIR.txt"))
  if (length(extra)) note("Replay wrote unexpected Console files: ", paste(extra, collapse = ", "))
  n_lines <- 0L
  for (f in ref_files) {
    new_path <- file.path(replay_dir, f)
    if (!file.exists(new_path)) {
      note(f, ": missing from the replay")
      next
    }
    ref <- read_console(file.path(ref_dir, "console", f))
    new <- read_console(new_path)
    if (length(ref) != length(new)) {
      note(f, ": ", length(new), " lines in the replay, ", length(ref), " in the reference")
    }
    for (i in seq_len(min(length(ref), length(new)))) {
      n_lines <- n_lines + 1L
      rt <- tokens(ref[i])
      nt <- tokens(new[i])
      ok <- length(rt) == length(nt) &&
        all(vapply(seq_along(rt), function(k) same_token(rt[k], nt[k]), logical(1)))
      if (!ok) {
        note(f, " line ", i, " differs\n    reference: ", ref[i], "\n    replay:    ", new[i])
      }
    }
    if (length(ref) != length(new)) {
      longer <- if (length(new) > length(ref)) new else ref
      who <- if (length(new) > length(ref)) "only in replay" else "only in reference"
      for (i in seq(min(length(ref), length(new)) + 1L, length(longer))) {
        note(f, " line ", i, " ", who, ": ", longer[i])
      }
    }
  }

  figures <- sort(list.files(file.path(ref_dir, "figures"), pattern = "[.]png$"))
  run_file <- file.path(replay_dir, "RUN_DIR.txt")
  if (!file.exists(run_file)) {
    note("RUN_DIR.txt is missing, so the figures cannot be checked")
  } else {
    fig_dir <- file.path(readLines(run_file, warn = FALSE)[1], "figures")
    for (png in figures) {
      p <- file.path(fig_dir, png)
      if (!file.exists(p) || file.size(p) == 0) {
        note("Figure not written: ", p)
      } else {
        con <- file(p, "rb")
        sig <- readBin(con, "raw", 8L)
        close(con)
        if (!identical(sig, as.raw(c(0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a)))) {
          note("Not a PNG file: ", p)
        }
      }
    }
  }

  cat("Compared", length(ref_files), "Console files (", n_lines, "lines ) and",
      length(figures), "figure files.\n")
  if (length(problems)) {
    cat("FAIL:", length(problems), "difference(s)\n")
    cat(paste0("- ", problems), sep = "\n")
    quit(status = 1)
  }
  cat("PASS: the Console output matches reference/console and all figures were written.\n")
})
