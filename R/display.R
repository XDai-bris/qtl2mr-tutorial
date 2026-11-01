# Printing only: stored numbers keep full precision.
format_number <- function(x) {
  vapply(x, function(z) {
    if (is.na(z)) return("NA")
    if (!is.finite(z)) return(as.character(z))
    if (z != 0 && abs(z) < 0.001) return(formatC(z, format = "e", digits = 3))
    sub("\\.?0+$", "", sprintf("%.3f", z))
  }, character(1))
}
show_table <- function(x) {
  if (is.matrix(x)) x <- as.data.frame(x)
  out <- x
  for (nm in names(out)) {
    if (is.numeric(out[[nm]])) out[[nm]] <- format_number(out[[nm]])
  }
  print(out, row.names = FALSE)
  invisible(out)
}
