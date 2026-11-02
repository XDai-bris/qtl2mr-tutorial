# =============================================================================
# simulation_functions.R
#
# Helper functions for simulate_world.R. Each function has a header that says
# what it does, what goes in, what comes out, and how it works.
#
# The functions fall into five groups:
#   A. Set-up            check_package_versions, set_genotype_stream
#   B. Genome            make_marker_map, place_qtl
#   C. Mice              simulate_do_population, as_qtl2_founder_probs,
#                        subset_founder_probs
#   D. Traits & cohorts  centre_pattern, make_founder_effects, split_cohorts,
#                        make_cohort_file
#   E. Design checks     founder_balance, select_instruments, count_recovered,
#                        run_design_checks
# =============================================================================


# -----------------------------------------------------------------------------
# A. Set-up
# -----------------------------------------------------------------------------

# check_package_versions()
#
# What:   Stops with a clear message unless the package versions used to make
#         the workshop data are loaded.
# In:     versions - named character vector, e.g. c(simcross = "0.10").
# Out:    nothing (invisible TRUE); stops on a mismatch.
# Why:    Simulated data depend on the exact code that drew the random
#         numbers. Another version may give different mice.
check_package_versions <- function(versions) {

  for (pkg in names(versions)) {

    found <- as.character(utils::packageVersion(pkg))

    if (found != versions[[pkg]]) {
      stop(sprintf("This script needs %s %s, but %s is loaded.",
                   pkg, versions[[pkg]], found), call. = FALSE)
    }
  }

  invisible(TRUE)
}


# set_genotype_stream()
#
# What:   Sets R's random-number generator to the stream used for genotypes.
# In:     seed - the master seed of this attempt.
# Out:    nothing; changes R's global random-number state.
# How:    One master seed gives several independent random-number streams
#         (R's "L'Ecuyer-CMRG" generator). qtl2founderSim uses streams 1-4
#         for the traits, the choice of mice and the cohorts; simcross uses
#         stream 5 for the genotypes. Because every step has its own stream,
#         changing one step never changes the random numbers of another.
set_genotype_stream <- function(seed) {

  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  set.seed(seed, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion",
           sample.kind = "Rejection")

  # Move forward five streams: stream 5 is the genotype stream.
  state <- .Random.seed
  for (i in 1:5) state <- parallel::nextRNGStream(state)

  assign(".Random.seed", state, envir = .GlobalEnv)
  invisible(NULL)
}


# -----------------------------------------------------------------------------
# B. Genome
# -----------------------------------------------------------------------------

# make_marker_map()
#
# What:   Places markers evenly along the 19 mouse autosomes.
# In:     chr_lengths - named vector of chromosome lengths in cM.
#         n_markers   - total number of markers.
# Out:    list(by_chr = list of named position vectors, one per chromosome,
#              table  = data.frame(marker, chr, pos)).
# How:    Each chromosome gets markers in proportion to its length (largest
#         remainders get the leftover markers). Within a chromosome, marker i
#         of n sits at (i - 0.5) * length / n. Names look like "chr2_050".
make_marker_map <- function(chr_lengths, n_markers) {

  # How many markers each chromosome gets.
  share     <- n_markers * chr_lengths / sum(chr_lengths)
  n_per_chr <- floor(share)

  leftover <- n_markers - sum(n_per_chr)
  if (leftover > 0) {
    largest_remainders <- order(share - n_per_chr, decreasing = TRUE)[seq_len(leftover)]
    n_per_chr[largest_remainders] <- n_per_chr[largest_remainders] + 1
  }

  # Evenly spaced positions on each chromosome.
  by_chr <- lapply(names(chr_lengths), function(chr) {
    n   <- n_per_chr[[chr]]
    pos <- (seq_len(n) - 0.5) * chr_lengths[[chr]] / n
    names(pos) <- sprintf("chr%s_%03d", chr, seq_len(n))
    pos
  })
  names(by_chr) <- names(chr_lengths)

  # The same map as one table, the format qtl2mr expects.
  table <- data.frame(
    marker = unlist(lapply(by_chr, names), use.names = FALSE),
    chr    = rep(names(by_chr), lengths(by_chr)),
    pos    = round(unlist(by_chr, use.names = FALSE), 4),
    stringsAsFactors = FALSE
  )

  list(by_chr = by_chr, table = table)
}


# place_qtl()
#
# What:   Puts each QTL on the marker nearest to its chosen position.
# In:     qtl         - data.frame with columns chr and frac (relative position,
#                       0 = start of the chromosome, 1 = end).
#         map         - output of make_marker_map().
#         chr_lengths - chromosome lengths in cM.
# Out:    qtl with two new columns: marker and pos (cM).
place_qtl <- function(qtl, map, chr_lengths) {

  qtl$marker <- vapply(seq_len(nrow(qtl)), function(i) {
    positions <- map$by_chr[[qtl$chr[i]]]
    target    <- qtl$frac[i] * chr_lengths[[qtl$chr[i]]]
    names(positions)[which.min(abs(positions - target))]
  }, character(1))

  qtl$pos <- map$table$pos[match(qtl$marker, map$table$marker)]

  qtl
}


# -----------------------------------------------------------------------------
# C. Mice
# -----------------------------------------------------------------------------

# simulate_do_population()
#
# What:   Breeds a Diversity Outbred (DO) population with simcross and returns
#         the founder make-up of every mouse in the last generation.
# In:     settings - the SETTINGS list of simulate_world.R (uses ngen, npairs,
#                    nkids, ccgen, founders).
#         chr_lengths, map - chromosome lengths (cM) and marker map.
# Out:    list(probs = founder probabilities in qtl2 format (see
#                      as_qtl2_founder_probs), pedigree = the simcross pedigree).
# How:    1. sim_do_pedigree() writes down who mates with whom: the eight
#            founder strains are crossed into pre-CC lines (ccgen = NULL uses
#            simcross's default inbreeding generations), then the DO is bred
#            for `ngen` generations of random mating avoiding siblings
#            (npairs breeding pairs; nkids pups per pair in the last generation).
#         2. sim_from_pedigree() simulates meiosis along that pedigree: every
#            chromosome a pup inherits is a recombined copy of its parents', so
#            after many generations each chromosome is a mosaic of short
#            pieces from the founders (typically all eight).
#         3. For each marker we read which two founder pieces a mouse carries
#            (one from each parent) and turn that into probabilities: 1 for a
#            founder carried twice, 0.5 + 0.5 for two different founders.
#            These are exact, because we know the truth; in real data qtl2
#            estimates them from genotypes.
#         Call set_genotype_stream() first so the population is reproducible.
simulate_do_population <- function(settings, chr_lengths, map) {

  # 1. The breeding plan. ccgen is passed only when it is set; otherwise
  #    simcross uses its own default.
  pedigree_args <- list(ngen = settings$ngen, npairs = settings$npairs,
                        nkids_per = settings$nkids, design = "nosib")
  if (!is.null(settings$ccgen)) pedigree_args$ccgen <- settings$ccgen

  pedigree <- do.call(simcross::sim_do_pedigree, pedigree_args)

  # 2. Meiosis along the pedigree (autosomes only).
  crossovers <- simcross::sim_from_pedigree(
    pedigree, L = chr_lengths, xchr = rep(FALSE, length(chr_lengths))
  )

  # Keep the mice of the last DO generation.
  last_gen <- as.character(pedigree$id[pedigree$do & pedigree$gen == settings$ngen])

  crossovers <- simcross::collapse_do_alleles(
    lapply(crossovers, function(chr) chr[last_gen])
  )

  # 3. Founder pair at each marker -> founder probabilities.
  genotypes <- simcross::convert2geno(crossovers, map$by_chr)

  mouse_ids <- sprintf("M%04d", seq_along(last_gen))
  founders  <- settings$founders

  probs <- lapply(names(chr_lengths), function(chr) {

    g <- genotypes[[chr]]          # mice x markers x 2 (the two founder copies)

    p <- array(0, dim = c(length(last_gen), length(founders), dim(g)[2]),
               dimnames = list(mouse_ids, founders, dimnames(g)[[2]]))

    for (f in seq_along(founders)) {
      p[, f, ] <- ((g[, , 1] == f) + (g[, , 2] == f)) / 2
    }

    p
  })
  names(probs) <- names(chr_lengths)

  list(probs = as_qtl2_founder_probs(probs, founders), pedigree = pedigree)
}


# as_qtl2_founder_probs()
#
# What:   Labels a list of probability arrays as a qtl2 "allele probabilities"
#         object, so qtl2 and qtl2mr functions accept it.
# In:     probs_list - one array per chromosome: mice x founders x markers.
#         founders   - founder names (A-H).
# Out:    the same list with the attributes qtl2 expects.
as_qtl2_founder_probs <- function(probs_list, founders) {

  attr(probs_list, "is_x_chr")    <- setNames(rep(FALSE, length(probs_list)),
                                              names(probs_list))
  attr(probs_list, "crosstype")   <- "do"
  attr(probs_list, "alleles")     <- founders
  attr(probs_list, "alleleprobs") <- TRUE
  class(probs_list) <- c("calc_genoprob", "list")

  probs_list
}


# subset_founder_probs()
#
# What:   Keeps only some mice from a founder-probability object.
# In:     probs - qtl2 founder probabilities; ids - mouse IDs to keep.
# Out:    founder probabilities for those mice, still in qtl2 format.
subset_founder_probs <- function(probs, ids) {

  kept <- lapply(unclass(probs), function(a) a[ids, , , drop = FALSE])

  as_qtl2_founder_probs(kept, attr(probs, "alleles"))
}


# -----------------------------------------------------------------------------
# D. Traits and cohorts
# -----------------------------------------------------------------------------

# centre_pattern()
#
# What:   Centres a founder pattern so its eight values average 0.
# In:     pattern - eight numbers, one per founder A-H.
# Out:    the centred pattern, named A-H.
# Why:    A QTL is a pattern of differences between founders. Centring removes
#         the overall level, which belongs to the trait mean, not to the QTL.
centre_pattern <- function(pattern, founders = LETTERS[1:8]) {

  centred <- unname(pattern) - mean(pattern)
  names(centred) <- founders

  centred
}


# make_founder_effects()
#
# What:   Turns each QTL's founder pattern into founder effects of the size we
#         want.
# In:     probs    - founder probabilities of all mice.
#         ids      - the mice of the simulated population (the 1,600).
#         qtl      - the QTL table (id, trait, marker).
#         patterns - list of centred founder patterns, named by QTL id.
#         variance - named vector: the share of variance each QTL should explain,
#                    by trait, e.g. c(X = 0.08, Y = 0.45).
#         traits   - all trait names, e.g. c("U", "X", "Y").
# Out:    list(for_model = named list, one 8 x 3 matrix per QTL marker
#                          (founders x traits), the format qtl2founderSim wants,
#              table     = long data.frame of the true founder effects).
# How:    For each QTL, a mouse's genetic value is its founder probabilities
#         at the QTL marker times the pattern (for example 0.5 x pattern[B] +
#         0.5 x pattern[F]). The pattern is then multiplied by one number so
#         that these genetic values have exactly the target variance in the
#         simulated population. Only the QTL's own trait gets the effect; the
#         other trait columns stay 0.
make_founder_effects <- function(probs, ids, qtl, patterns, variance, traits) {

  founders  <- attr(probs, "alleles")
  for_model <- list()
  table     <- list()

  for (i in seq_len(nrow(qtl))) {

    trait   <- qtl$trait[i]
    marker  <- qtl$marker[i]
    pattern <- patterns[[qtl$id[i]]]

    # Founder probabilities of every mouse at this QTL marker.
    p_at_qtl <- unclass(probs)[[qtl$chr[i]]][ids, , marker]

    # Genetic value of each mouse before scaling, then the scale factor.
    value <- drop(p_at_qtl %*% pattern)
    scale <- sqrt(variance[[trait]] / mean((value - mean(value))^2))

    effects <- matrix(0, nrow = length(founders), ncol = length(traits),
                      dimnames = list(founders, traits))
    effects[, trait] <- pattern * scale

    for_model[[marker]] <- effects

    table[[i]] <- data.frame(qtl = qtl$id[i], trait = trait, marker = marker,
                             founder = founders, effect = effects[, trait])
  }

  list(for_model = for_model, table = do.call(rbind, table))
}


# split_cohorts()
#
# What:   Reads the three cohorts from the simulated measurement plan.
# In:     sets - the mice with X measured and the mice with Y measured
#                (sim$truth$membership$sets from qtl2founderSim).
# Out:    list(cohort1 = X and Y measured, cohortA = X only, cohortB = Y only),
#         each a sorted vector of mouse IDs.
# How:    qtl2founderSim's two_sample() shuffles the population once and deals
#         out the mice: the first group gets both traits measured, the next
#         only X, the last only Y. Cohort 1 is the overlap of the X and Y sets.
split_cohorts <- function(sets) {

  list(
    cohort1 = sort(intersect(sets$X, sets$Y)),
    cohortA = sort(setdiff(sets$X, sets$Y)),
    cohortB = sort(setdiff(sets$Y, sets$X))
  )
}


# make_cohort_file()
#
# What:   Builds everything the workshop needs for one cohort.
# In:     probs    - founder probabilities of all mice.
#         ids      - the cohort's mouse IDs.
#         pheno    - matrix of simulated traits (all mice; unmeasured = NA).
#         measured - the traits measured in this cohort, e.g. "X".
#         map      - marker map table.
#         kinship_digits - kinship values are rounded to this many digits.
# Out:    list(genoprobs, kinship, kinship_overall, pheno, map, scan), the
#         structure of data/cohort1.rds, cohortA.rds and cohortB.rds.
# How:    Kinship is computed with qtl2 from this cohort's own mice, both
#         "leave one chromosome out" (LOCO, used by qtl2mr) and overall (for the
#         heatmap). The genome scan is qtl2::scan1 with the LOCO kinship, the
#         same scan a qtl2 user would run.
make_cohort_file <- function(probs, ids, pheno, measured, map, kinship_digits = 6) {

  genoprobs <- subset_founder_probs(probs, ids)

  kinship_loco <- qtl2::calc_kinship(genoprobs, type = "loco", cores = 1)
  kinship_loco <- lapply(kinship_loco, round, kinship_digits)

  kinship_overall <- qtl2::calc_kinship(genoprobs, type = "overall", cores = 1)
  kinship_overall <- round(kinship_overall, kinship_digits)

  cohort_pheno <- as.data.frame(pheno[ids, measured, drop = FALSE])
  stopifnot(all(is.finite(as.matrix(cohort_pheno))))

  scan <- qtl2::scan1(genoprobs, cohort_pheno, kinship = kinship_loco, cores = 1)

  list(genoprobs = genoprobs, kinship = kinship_loco,
       kinship_overall = kinship_overall, pheno = cohort_pheno,
       map = map, scan = scan)
}


# -----------------------------------------------------------------------------
# E. Design checks (no MR estimate is computed anywhere in this folder)
# -----------------------------------------------------------------------------

# founder_balance()
#
# What:   Two numbers that show the eight founders are well represented.
# In:     probs - founder probabilities of all mice; ids - the population;
#         cohorts - list of cohort files; qtl_markers - the QTL markers.
# Out:    list(freq = each founder's share of the genome in the population,
#              min_copies = the fewest copies of any founder at any QTL marker
#                           in any cohort).
founder_balance <- function(probs, ids, cohorts, qtl_markers) {

  arrays    <- unclass(probs)
  n_markers <- sum(vapply(arrays, function(a) dim(a)[3], numeric(1)))

  # Share of the genome from each founder, over all mice and markers.
  per_chr <- vapply(arrays, function(a) apply(a[ids, , , drop = FALSE], 2, sum),
                    numeric(length(attr(probs, "alleles"))))
  freq <- rowSums(per_chr) / (length(ids) * n_markers)

  # Copies of each founder at each QTL marker, in each cohort (2 per mouse).
  copies <- sapply(cohorts, function(cohort) {
    sapply(qtl_markers, function(marker) {
      chr <- cohort$map$chr[match(marker, cohort$map$marker)]
      p   <- unclass(cohort$genoprobs)[[chr]][, , marker]
      min(2 * colSums(p))
    })
  })

  list(freq = freq, min_copies = min(copies))
}


# select_instruments()
#
# What:   Picks instruments for one trait in one cohort with qtl2mr, using the
#         workshop rule: the top marker of each region with LOD >= 3, kept only
#         if its F statistic is at least 4. (F is the joint test of the eight
#         founder effects at that marker, divided by its 7 degrees of freedom.)
# In:     cohort - a cohort file; trait - "X" or "Y"; lod_min, F_min - the rule.
# Out:    qtl2mr's instrument table (one row per region; column `selected`).
select_instruments <- function(cohort, trait, lod_min, F_min) {

  resources <- qtl2mr::qtl2mr_prepare_resources(
    cohort$genoprobs, cohort$kinship, cohort$pheno,
    covar = NULL, marker_map = cohort$map
  )

  qtl2mr::qtl2mr_select_instruments(
    resources, cohort$scan[, trait, drop = FALSE], exposures = trait,
    lod_thresh = lod_min, F_X_min = F_min, progress = FALSE
  )
}


# count_recovered()
#
# What:   Compares selected instruments with the true QTL.
# In:     instruments - output of select_instruments(); truth - true QTL rows;
#         window_cM - how close counts as "found".
# Out:    list(n_recovered = true QTL with a selected instrument within the
#              window, recovered = their ids, n_selected, n_chance = selected
#              instruments with no true QTL within the window).
count_recovered <- function(instruments, truth, window_cM) {

  chosen <- instruments[instruments$selected, , drop = FALSE]

  found <- vapply(seq_len(nrow(truth)), function(i) {
    any(chosen$chr == truth$chr[i] & abs(chosen$pos - truth$pos[i]) <= window_cM)
  }, logical(1))

  chance <- vapply(seq_len(nrow(chosen)), function(j) {
    !any(truth$chr == chosen$chr[j] & abs(truth$pos - chosen$pos[j]) <= window_cM)
  }, logical(1))

  list(n_recovered = sum(found), recovered = truth$id[found],
       n_selected = nrow(chosen), n_chance = sum(chance))
}


# run_design_checks()
#
# What:   The four checks a simulated world must pass before it is used. They
#         look only at design properties, never at an MR result.
#           C1 founders  - every founder makes up 9-16% of the genome, and
#                          has at least 20 copies at every QTL in every cohort;
#           C2 X-QTL     - at least 8 of the 10 X-QTL are found as instruments
#                          in Cohort 1 and in Cohort A;
#           C3 bias      - the plain regression slope of Y on X in Cohort 1 is
#                          between 0.8 and 1.2 (the bias we want to show);
#           C4 Y-QTL     - at least 3 of the 4 Y-QTL are found as Y-but-not-X
#                          instruments in Cohort 1 (Y LOD >= 3, F >= 4, and
#                          X LOD < 3 at that marker).
# In:     world pieces (probs, population ids, cohort files, QTL table) and
#         the settings.
# Out:    one-row data.frame with every checked value and pass = TRUE/FALSE,
#         plus the instrument tables (attribute "instruments").
run_design_checks <- function(probs, ids, cohorts, qtl, settings, attempt, seed) {

  rule  <- settings$instrument_rule
  limit <- settings$checks

  # C1: founder representation.
  balance <- founder_balance(probs, ids, cohorts, qtl$marker)

  # C2: X instruments in Cohort 1 and in Cohort A.
  x_in_cohort1 <- select_instruments(cohorts$cohort1, "X", rule$lod_min, rule$F_min)
  x_in_cohortA <- select_instruments(cohorts$cohortA, "X", rule$lod_min, rule$F_min)

  x_truth  <- qtl[qtl$trait == "X", ]
  found_x1 <- count_recovered(x_in_cohort1, x_truth, rule$window_cM)
  found_xA <- count_recovered(x_in_cohortA, x_truth, rule$window_cM)

  # C3: the plain association of Y with X in Cohort 1.
  slope <- unname(coef(lm(Y ~ X, data = cohorts$cohort1$pheno))[2])

  # C4: Y instruments that are not X regions, in Cohort 1.
  y_in_cohort1 <- select_instruments(cohorts$cohort1, "Y", rule$lod_min, rule$F_min)
  x_lod_there  <- unname(cohorts$cohort1$scan[y_in_cohort1$marker, "X"])
  y_in_cohort1$selected <- y_in_cohort1$selected & x_lod_there < rule$x_lod_max

  found_y1 <- count_recovered(y_in_cohort1, qtl[qtl$trait == "Y", ], rule$window_cM)

  # Pass or fail.
  C1 <- all(balance$freq >= limit$freq[1] & balance$freq <= limit$freq[2]) &&
    balance$min_copies >= limit$min_copies
  C2 <- found_x1$n_recovered >= limit$x_recovered &&
    found_xA$n_recovered >= limit$x_recovered
  C3 <- slope >= limit$slope[1] && slope <= limit$slope[2]
  C4 <- found_y1$n_recovered >= limit$y_recovered

  result <- data.frame(
    attempt = attempt, seed = seed,
    C1_founder_freq_min = min(balance$freq), C1_founder_freq_max = max(balance$freq),
    C1_min_copies = balance$min_copies,
    C2_xqtl_found_cohort1 = found_x1$n_recovered,
    C2_xqtl_found_cohortA = found_xA$n_recovered,
    C3_slope_cohort1 = slope,
    C4_yqtl_found_cohort1 = found_y1$n_recovered,
    info_x_selected_cohort1 = found_x1$n_selected,
    info_x_chance_cohort1 = found_x1$n_chance,
    info_x_selected_cohortA = found_xA$n_selected,
    info_x_chance_cohortA = found_xA$n_chance,
    info_y_selected_cohort1 = found_y1$n_selected,
    info_y_chance_cohort1 = found_y1$n_chance,
    C1 = C1, C2 = C2, C3 = C3, C4 = C4,
    pass = C1 && C2 && C3 && C4
  )

  attr(result, "instruments") <- list(x_cohort1 = x_in_cohort1, x_cohortA = x_in_cohortA,
                                      y_cohort1 = y_in_cohort1)
  result
}
