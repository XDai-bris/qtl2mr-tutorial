# Figure helpers. Colour-blind-safe Okabe-Ito colours; large text for a projector.
founder_colours <- function() {
  c(A = "#E69F00", B = "#56B4E9", C = "#009E73", D = "#F0E442",
    E = "#0072B2", F = "#D55E00", G = "#CC79A7", H = "#333333")
}
col_obs <- "#D55E00"; col_xy <- "#0072B2"; col_yx <- "#009E73"; col_truth <- "#333333"
arrow <- "\u2192"  # a right arrow; the escape keeps this file plain ASCII

big_text <- function(...) {
  par(cex.lab = 1.35, cex.axis = 1.15, cex.main = 1.45, mgp = c(2.9, 0.8, 0), ...)
}

# Draw in the Plots pane, then save the same figure as a PNG in the results folder.
draw_workshop_figure <- function(filename, draw, width = 1600, height = 900) {
  if (interactive()) tryCatch(draw(), error = function(e) {
    cat("Could not draw this figure in the Plots pane (", conditionMessage(e), ").\n",
        "It is still saved as ", file.path(run_dir, "figures", filename),
        "; enlarge the Plots pane and run the plot line again.\n", sep = "")
  })
  png(file.path(run_dir, "figures", filename), width = width, height = height, res = 150)
  tryCatch(draw(), finally = dev.off())
  cat("Figure saved:", file.path("figures", filename), "\n")
  invisible(filename)
}

plot_founder_bars <- function(cohort, marker, n = 6) {
  chr <- cohort$map$chr[cohort$map$marker == marker]
  p <- cohort$genoprobs[[chr]][seq_len(n), , marker]
  draw_workshop_figure("01_founder_probabilities.png", function() {
    big_text(mar = c(5, 5, 4, 1))
    barplot(t(p), col = founder_colours(), border = "white", ylim = c(0, 1.3), axes = FALSE,
            names.arg = rownames(p), xlab = "Mouse", ylab = "Founder probability",
            main = paste0("Each mouse carries two founder haplotypes (marker ", marker, ")"))
    axis(2, at = seq(0, 1, by = 0.2), las = 1)  # probabilities: the axis stops at 1
    legend("top", legend = names(founder_colours()), fill = founder_colours(),
           ncol = 8, bty = "n", cex = 1.2, border = "white")
  })
}

plot_kinship <- function(cohort, n = 60) {
  K <- cohort$kinship_overall[seq_len(n), seq_len(n)]
  cols <- hcl.colors(64, "YlGnBu", rev = TRUE)
  lim <- c(0, max(K))
  draw_workshop_figure("02_kinship.png", function() {
    layout(matrix(1:2, nrow = 1), widths = c(1, 0.25))
    big_text(mar = c(5, 5, 4, 1))
    image(seq_len(n), seq_len(n), K[, n:1], col = cols, zlim = lim, axes = FALSE,
          xlab = "Mouse", ylab = "Mouse",
          main = paste("Kinship among the first", n, "mice of Cohort 1"))
    axis(1, at = c(1, n), labels = rownames(K)[c(1, n)])
    axis(2, at = c(n, 1), labels = rownames(K)[c(1, n)], las = 1)
    box()
    par(mar = c(5, 0.5, 4, 3.5))
    v <- seq(lim[1], lim[2], length.out = 64)
    image(1, v, matrix(v, nrow = 1), col = cols, zlim = lim, axes = FALSE, xlab = "", ylab = "")
    axis(4, las = 1)
    box()
    mtext("kinship", side = 3, line = 1, cex = 1.2)
  })
}

plot_summary <- function(est) {
  draw_workshop_figure("03_summary_estimates.png", function() {
    big_text(mar = c(5, 15, 4, 2))
    y <- rev(seq_len(nrow(est)))
    cols <- c(col_obs, col_xy, col_yx)
    plot(NA, xlim = range(c(-0.2, 1.3, est$lower, est$upper)), ylim = c(0.4, nrow(est) + 0.6),
         yaxt = "n", ylab = "", xlab = "Slope (outcome units per exposure unit)",
         main = "Cohort 1: association vs causal effect")
    abline(v = 0, col = "grey80")
    axis(2, at = y, labels = est$label, las = 1, tick = FALSE)
    for (i in seq_len(nrow(est))) {
      if (is.finite(est$truth[i])) {
        segments(est$truth[i], y[i] - 0.38, est$truth[i], y[i] + 0.38,
                 lty = 2, lwd = 3, col = col_truth)
        text(est$truth[i], y[i] + 0.45, paste("truth", est$truth[i]), cex = 1.1, col = col_truth)
      }
      segments(est$lower[i], y[i], est$upper[i], y[i], lwd = 5, col = cols[i])
      points(est$estimate[i], y[i], pch = 19, cex = 2.2, col = cols[i])
    }
  })
}

genome_positions <- function(map, gap = 15) {
  chrs <- unique(map$chr)
  len <- tapply(map$pos, map$chr, max)[chrs]
  start <- c(0, cumsum(len + gap))[seq_along(chrs)]
  names(start) <- chrs
  list(x = map$pos + start[map$chr], start = start, len = len, chrs = chrs)
}

plot_miami <- function(cohort_x, cohort_y, markers) {
  map <- cohort_x$map
  g <- genome_positions(map)
  lx <- cohort_x$scan[map$marker, "X"]
  ly <- cohort_y$scan[map$marker, "Y"]
  top <- ceiling(max(lx, ly)) + 1
  draw_workshop_figure("04_miami.png", function() {
    big_text(mar = c(5, 6, 4, 2))
    plot(NA, xlim = range(g$x) + c(-5, 5), ylim = c(-top, top), xaxt = "n", yaxt = "n", xaxs = "i",
         xlab = "Chromosome", ylab = "LOD",
         main = "X instruments found in Cohort A, read off in Cohort B")
    ticks <- pretty(c(0, top))
    axis(2, at = c(-ticks, ticks), labels = c(ticks, ticks), las = 1)
    text(min(g$x), top * 0.92, "X in Cohort A", adj = 0, cex = 1.3, font = 2, col = "#0072B2")
    text(min(g$x), -top * 0.92, "Y in Cohort B", adj = 0, cex = 1.3, font = 2, col = "#D55E00")
    for (i in seq_along(g$chrs)) {
      ix <- map$chr == g$chrs[i]
      odd <- i %% 2 == 1
      lines(g$x[ix], lx[ix], col = if (odd) "#0072B2" else "#56B4E9", lwd = 2)
      lines(g$x[ix], -ly[ix], col = if (odd) "#D55E00" else "#E69F00", lwd = 2)
    }
    abline(h = 0, col = "grey60")
    abline(h = c(-3, 3), lty = 3, col = "grey50")
    abline(v = g$x[match(markers, map$marker)], lty = 2, lwd = 1.6, col = col_truth)
    axis(1, at = g$start + g$len / 2, labels = g$chrs, tick = FALSE, cex.axis = 1)
    legend("topright", c("X instrument (from Cohort A)", "LOD = 3"), lty = c(2, 3),
           lwd = c(1.6, 1), col = c(col_truth, "grey50"), bty = "n", cex = 1.1)
  }, width = 2000)
}

plot_mr_scatter <- function(pts, theta, truth = 0.5) {
  draw_workshop_figure("05_mr_scatter.png", function() {
    big_text(mar = c(5, 5.5, 4, 2))
    xr <- c(0, max(pts$x + 1.96 * pts$se_x) * 1.05)
    yr <- range(c(0, pts$y - 1.96 * pts$se_y, pts$y + 1.96 * pts$se_y))
    plot(NA, xlim = xr, ylim = yr, las = 1,
         xlab = "Instrument effect on X (Cohort A)",
         ylab = "Effect on Y of the same pattern (Cohort B)",
         main = "Two-sample MR: one point per instrument")
    abline(h = 0, col = "grey80")
    segments(pts$x - 1.96 * pts$se_x, pts$y, pts$x + 1.96 * pts$se_x, pts$y, col = "grey55", lwd = 2)
    segments(pts$x, pts$y - 1.96 * pts$se_y, pts$x, pts$y + 1.96 * pts$se_y, col = "grey55", lwd = 2)
    abline(0, truth, lty = 2, lwd = 3, col = col_truth)
    abline(0, theta, lwd = 4, col = col_xy)
    points(pts$x, pts$y, pch = 21, bg = col_xy, col = "white", cex = 2.4)
    legend("topleft", c(sprintf("MR slope %.2f", theta), sprintf("truth %.1f", truth)),
           lty = c(1, 2), lwd = c(4, 3), col = c(col_xy, col_truth), bty = "n", cex = 1.2)
  })
}

plot_founder_scatter <- function(fe, theta, truth = 0.5) {
  draw_workshop_figure("06_founder_scatter.png", function() {
    big_text(mar = c(5, 5.5, 4, 2))
    plot(fe$x, fe$y, pch = 21, bg = founder_colours()[fe$founder], col = "white", cex = 2,
         las = 1, xlab = "Founder effect on X (Cohort A)", ylab = "Founder effect on Y (Cohort B)",
         main = "Every founder at every instrument")
    abline(h = 0, v = 0, col = "grey80")
    abline(0, truth, lty = 2, lwd = 3, col = col_truth)
    abline(0, theta, lwd = 4, col = col_xy)
    legend("topleft", legend = names(founder_colours()), pt.bg = founder_colours(), pch = 21,
           col = "white", ncol = 4, bty = "n", cex = 1.1, pt.cex = 1.8)
  })
}

plot_founder_vs_snp <- function(fe, carriers, snp, marker) {
  s <- as.numeric(fe$founder %in% carriers)
  snp_line <- snp[["effect"]] * (s - mean(s))
  draw_workshop_figure("07_founders_vs_snp.png", function() {
    big_text(mar = c(5, 5.5, 4, 2))
    x <- seq_len(nrow(fe))
    lim <- range(c(fe$effect - 1.96 * fe$se, fe$effect + 1.96 * fe$se, snp_line))
    plot(NA, xlim = c(0.5, nrow(fe) + 0.5), ylim = lim + c(0, diff(lim) * 0.15), xaxt = "n", las = 1,
         xlab = "Founder", ylab = "Effect on X (Cohort A)",
         main = paste("One QTL (", marker, "): 8 founder effects vs one SNP split", sep = ""))
    axis(1, at = x, labels = fe$founder)
    abline(h = 0, col = "grey80")
    segments(x, fe$effect - 1.96 * fe$se, x, fe$effect + 1.96 * fe$se, lwd = 3,
             col = founder_colours()[fe$founder])
    points(x, fe$effect, pch = 21, bg = founder_colours()[fe$founder], col = "white", cex = 2.6)
    segments(x - 0.35, snp_line, x + 0.35, snp_line, lwd = 5, col = "grey55")
    legend("topright", c("Founder effect (95% interval)",
                         paste0("SNP: carriers ", paste(carriers, collapse = ""), " vs the rest")),
           pch = c(21, NA), pt.bg = c("grey30", NA), lty = c(NA, 1), lwd = c(NA, 5),
           col = c("white", "grey55"), pt.cex = 2, bty = "n", cex = 1.1)
  })
}

plot_founder_vs_snp_precision <- function(tab, share, truth = 0.5) {
  draw_workshop_figure("08_founder_vs_snp_precision.png", function() {
    layout(matrix(1:2, nrow = 1), widths = c(1, 1.15))
    big_text(mar = c(5, 6.5, 4, 1))
    y <- c(2, 1)
    plot(NA, xlim = range(c(0, 1, tab$lower, tab$upper)), ylim = c(0.5, 2.6), yaxt = "n", ylab = "",
         xlab = "Two-sample MR estimate", main = "Same mice, same instruments")
    axis(2, at = y, labels = tab$method, las = 1, tick = FALSE)
    segments(truth, 0.5, truth, 2.4, lty = 2, lwd = 3, col = col_truth)
    text(truth, 2.52, paste("truth", truth), cex = 1.1, col = col_truth)
    segments(tab$lower, y, tab$upper, y, lwd = 5, col = c(col_xy, "#E69F00"))
    points(tab$estimate, y, pch = 19, cex = 2.2, col = c(col_xy, "#E69F00"))
    text(tab$upper, y - 0.2, sprintf("width %.2f", tab$upper - tab$lower), adj = 1, cex = 1.1)
    big_text(mar = c(5, 5, 4, 1))
    barplot(share$share, names.arg = share$chr, ylim = c(0, 1.1), las = 1, border = "white",
            col = "#E69F00", xlab = "Instrument (chromosome)",
            ylab = "Share of the X association kept", main = "What the best single SNP keeps")
    abline(h = 1, lty = 2, col = col_truth)
  }, width = 2000)
}
