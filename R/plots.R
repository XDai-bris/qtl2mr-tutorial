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

# For one mouse, split the founder probabilities along a chromosome (founders x
# markers, exactly 0, 1/2 or 1 here) into two chromosome copies. At each marker
# the mouse carries two founders; each copy keeps the founder it had at the
# previous marker whenever it can, so that long stretches stay together.
founder_tracks <- function(p) {
  founders <- rownames(p)
  tracks <- matrix(NA_character_, 2, ncol(p))
  for (j in seq_len(ncol(p))) {
    o <- order(p[, j], decreasing = TRUE)
    pair <- if (p[o[1], j] > 0.75) founders[o[c(1, 1)]] else founders[o[1:2]]
    if (j > 1 && sum(rev(pair) == tracks[, j - 1]) > sum(pair == tracks[, j - 1])) {
      pair <- rev(pair)
    }
    tracks[, j] <- pair
  }
  tracks
}

# Founder mosaic: two strips per mouse (its two copies of one chromosome),
# coloured by the founder each stretch comes from; a dashed line marks a marker.
plot_founder_mosaic <- function(cohort, chr, marker, n = 6) {
  map <- cohort$map[cohort$map$chr == chr, ]
  stopifnot(marker %in% map$marker)
  pos <- map$pos
  step <- diff(pos)
  edges <- c(max(0, pos[1] - step[1] / 2), pos[-1] - step / 2, pos[length(pos)] + step[length(step)] / 2)
  # n mice spread across the cohort, so that close relatives are unlikely.
  rows <- round(seq(1, dim(cohort$genoprobs[[chr]])[1], length.out = n))
  probs <- cohort$genoprobs[[chr]][rows, , map$marker, drop = FALSE]
  at <- pos[map$marker == marker]
  cols <- founder_colours()
  draw_workshop_figure("01_founder_mosaic.png", function() {
    big_text(mar = c(5, 7, 4, 1))
    plot(NA, xlim = range(edges), ylim = c(0.45, n + 1.5), xaxs = "i", yaxt = "n", bty = "n",
         xlab = paste0("Position on chromosome ", chr, " (cM)"), ylab = "",
         main = "Each mouse is a mosaic of the eight founders")
    for (i in seq_len(n)) {
      tracks <- founder_tracks(probs[i, , ])
      y <- n + 1 - i
      for (k in 1:2) {
        runs <- rle(tracks[k, ])
        last <- cumsum(runs$lengths)
        first <- c(1, last[-length(last)] + 1)
        bottom <- if (k == 1) y + 0.02 else y - 0.34
        rect(edges[first], bottom, edges[last + 1], bottom + 0.32, col = cols[runs$values], border = NA)
      }
    }
    axis(2, at = n + 1 - seq_len(n), labels = dimnames(probs)[[1]], las = 1, tick = FALSE)
    segments(at, 0.55, at, n + 0.45, lty = 2, lwd = 3, col = col_truth)
    text(at, n + 0.75, marker, cex = 1.1, font = 2, col = col_truth)
    legend("top", legend = names(cols), fill = cols, border = "white", horiz = TRUE,
           bty = "n", cex = 1.1, x.intersp = 0.6)
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
  # Rows 1-2 answer "does X affect Y?" and rows 3-4 "does Y affect X?". In each
  # group the truth (the causal effect) is one dashed line across both rows.
  draw_workshop_figure("03_summary_estimates.png", function() {
    big_text(mar = c(5, 15, 4, 2))
    y <- c(5.2, 4.2, 2.2, 1.2)
    group <- c(1, 1, 2, 2)
    titles <- c("Does X affect Y?", "Does Y affect X?")
    cols <- c(col_obs, col_xy, col_obs, col_yx)
    xlim <- range(c(-0.2, 1.3, est$lower, est$upper))
    plot(NA, xlim = xlim, ylim = c(0.6, 6.2), yaxt = "n", ylab = "",
         xlab = "Slope (outcome units per exposure unit)",
         main = "Cohort 1: association vs causal effect")
    abline(v = 0, col = "grey80")
    abline(h = 3.2, col = "grey60", lty = 3, lwd = 2)
    axis(2, at = y, labels = est$label, las = 1, tick = FALSE)
    for (g in 1:2) {
      rows <- which(group == g)
      top <- max(y[rows])
      text(xlim[2], top + 0.72, titles[g], adj = c(1, 0.5), font = 2, cex = 1.2)
      truth <- est$truth[rows][is.finite(est$truth[rows])]
      if (length(truth)) {
        segments(truth, min(y[rows]) - 0.38, truth, top + 0.38, lty = 2, lwd = 3, col = col_truth)
        text(truth, top + 0.6, paste("truth", truth), cex = 1.1, col = col_truth)
      }
    }
    for (i in seq_len(nrow(est))) {
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
    big_text(mar = c(6.8, 6, 4, 2))
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
    # The legend sits below the plot, so that it covers no instrument line.
    legend("bottom", inset = c(0, -0.3), xpd = NA, horiz = TRUE,
           c("X instrument (from Cohort A)", "LOD = 3"), lty = c(2, 3),
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
