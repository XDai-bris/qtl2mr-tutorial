# Presentation only: stored numbers and fitting inputs keep full precision.
format_number <- function(x, p_value = FALSE) {
  vapply(x, function(z) {
    if (is.na(z)) return("NA")
    if (!is.finite(z)) return(as.character(z))
    if (p_value && z < 0.001) return("<0.001")
    if (z != 0 && abs(z) < 0.001) return(formatC(z, format = "e", digits = 3))
    sub("\\.?0+$", "", sprintf("%.3f", z))
  }, character(1))
}
show_table <- function(x) {
  if (is.matrix(x)) x <- as.data.frame(x)
  if (is.numeric(x)) x <- data.frame(value = x)
  out <- x
  for (nm in names(out)) {
    if (is.numeric(out[[nm]])) {
      out[[nm]] <- format_number(out[[nm]], grepl("(^p$|_p$|^Q_p|p_value)", nm))
    }
  }
  print(out, row.names = FALSE)
  invisible(out)
}
