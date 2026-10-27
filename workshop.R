source("R/install.R")
source("R/setup.R")
data("qtl2mr_example", package="qtl2mr")
validate_qtl2mr_input(qtl2mr_example)
cat("Built-in example validated: TRUE\n")
plot_install_check(qtl2mr_example)
cat("Milestone 1: qtl2mr 0.6.0 is loaded and the example plot is visible.\n")


cat("\nTASK 02 | One genetic base, two possible directions\n")
demo <- readRDS("data/single_pair.rds")
cat("Mice:", nrow(demo$Forward), "| Markers:", nrow(demo$marker_map), "\n")
show_table(data.frame(mouse_id=rownames(demo$Forward)[1:4], head(demo$Forward, 4)))
show_table(data.frame(mouse_id=rownames(demo$covar)[1:4], head(demo$covar, 4)))
show_table(head(demo$marker_map, 4))
resources_forward <- qtl2mr_prepare_resources(
  allele_probs=demo$allele_probs, kinship=demo$kinship,
  pheno=demo$Forward, covar=demo$covar, marker_map=demo$marker_map
)
resources_reverse <- qtl2mr_prepare_resources(
  allele_probs=demo$allele_probs, kinship=demo$kinship,
  pheno=demo$Reverse, covar=demo$covar, marker_map=demo$marker_map
)
plot_single_pair_inputs(demo)
cat("Milestone 2: 400 mice loaded; both possible directions use the same genetic base.\n")


# PART A | From a scan to regional peaks
cat("\nTASK 03A | Many markers, one peak per region\n")
scan_forward <- readRDS("data/forward_scan.rds")
peaks_forward <- qtl2mr_select_instruments(
  resources_forward, scan_forward, exposures="x",
  lod_thresh=3, F_X_min=NULL, max_loci=NULL, progress=FALSE
)
cat("Markers scanned:", nrow(scan_forward), "| Candidate peaks:", nrow(peaks_forward), "\n")
show_table(head(peaks_forward[,c("marker","chr","pos","lod")],6))
plot_regional_selection(demo$marker_map, scan_forward, peaks_forward)
cat("Checkpoint: identify the highest marker in the shaded region.\n")

# PART B | Keep only peaks with adequate exposure strength
cat("\nTASK 03B | Apply the prespecified marginal F threshold\n")
instruments_forward <- qtl2mr_select_instruments(
  resources_forward, scan_forward, exposures="x",
  lod_thresh=3, F_X_min=4, max_loci=NULL, progress=FALSE
)
show_table(head(instruments_forward[,c("marker","lod","F_X_block","selected","drop_reason")],8))
cat("Retained:", sum(instruments_forward$selected),
    "| Rejected by F:", sum(instruments_forward$drop_reason=="F"), "\n")
plot_marginal_strength(instruments_forward)
cat("Milestone 3: regional peaks and the marginal F gate have produced the IV set.\n")


# PART A | Fit Forward once
cat("\nTASK 04A | Run Forward MR\n")
input_forward <- qtl2mr_pair_input(
  resources_forward, instruments_forward, exposure="x", outcome="y"
)
result_forward <- run_qtl2mr_analysis(
  input_forward, output_dir=file.path(run_dir,"forward"),
  standardize_phenotypes=FALSE, reported_estimator_policy="diagnostic_only",
  cross_cov_model="star", heterogeneity_method="chisq", diagnostic_bootstrap_B=0L,
  run_loo_mr=TRUE, run_instrument_sensitivity=FALSE,
  make_html=TRUE, make_figures=TRUE, verbose=FALSE, progress=FALSE
)
show_table(summary(result_forward)[,c("exposure","outcome","n_loci","estimate","ci_lower","ci_upper")])
plot_forward_result(result_forward)

# PART B | Read compatibility and influence
cat("\nTASK 04B | Read diagnostics with the estimate\n")
show_table(diagnostic_overview(result_forward))
loo_forward <- read_pair_table(result_forward,"loo_mr")
show_table(loo_forward[,c("marker","theta_minus_m","profile_ci_lower","profile_ci_upper")])
plot_forward_influence(result_forward,loo_forward)
cat("Milestone 4: I can read the estimate, interval, Q and leave-one-locus-out results.\n")


cat("\nTASK 05 | Repeat the workflow under Reverse\n")
scan_reverse <- readRDS("data/reverse_scan.rds")
instruments_reverse <- qtl2mr_select_instruments(
  resources_reverse, scan_reverse, exposures="x",
  lod_thresh=3, F_X_min=4, max_loci=NULL, progress=FALSE
)
cat("Candidate peaks:", nrow(instruments_reverse),
    "| Retained:", sum(instruments_reverse$selected), "\n")
input_reverse <- qtl2mr_pair_input(
  resources_reverse, instruments_reverse, exposure="x", outcome="y"
)
result_reverse <- run_qtl2mr_analysis(
  input_reverse, output_dir=file.path(run_dir,"reverse"),
  standardize_phenotypes=FALSE, reported_estimator_policy="diagnostic_only",
  cross_cov_model="star", heterogeneity_method="chisq", diagnostic_bootstrap_B=0L,
  run_loo_mr=TRUE, run_instrument_sensitivity=FALSE,
  make_html=TRUE, make_figures=TRUE, verbose=FALSE, progress=FALSE
)
show_table(diagnostic_overview(result_reverse))
plot_pair_comparison(result_forward,result_reverse)
cat("Milestone 5: Reverse also gives an MR estimate; its high Q needs explanation.\n")


# PART A | Read the SF report and calculate one partial R-squared
cat("\nTASK 06A | SF compares associations; it does not remove IVs\n")
sf_forward <- read_pair_table(result_forward,"block_first_stage")
sf_reverse <- read_pair_table(result_reverse,"block_first_stage")
direction <- direction_table(sf_forward,sf_reverse)
show_table(direction)

# Here: 400 mice, 5 covariate columns including the intercept,
# and 7 independent founder contrasts: residual df = 400 - 5 - 7.
k <- 7
residual_df <- 400 - 5 - k
F_example <- sf_forward$F_X_block[sf_forward$marker == "chr01_m040"]
R2_example <- k * F_example / (k * F_example + residual_df)
show_table(data.frame(marker="chr01_m040",F_X=F_example,
                      k=k,residual_df=residual_df,R2_X=R2_example))
plot_direction(sf_forward,sf_reverse)

# PART B | Explicitly choose a separate sensitivity IV set; still analyse X -> Y
cat("\nTASK 06B | Compare the original five IVs with a manual two-IV set\n")
keep_markers <- c("chr10_m040","chr11_m040")
iv_manual <- subset(instruments_reverse, selected & marker %in% keep_markers)
show_table(iv_manual[,c("marker","chr","F_X_block")])
input_manual <- qtl2mr_pair_input(
  resources_reverse, iv_manual, exposure="x", outcome="y"
)
result_manual <- run_qtl2mr_analysis(
  input_manual, output_dir=file.path(run_dir,"reverse_manual"),
  standardize_phenotypes=FALSE, reported_estimator_policy="diagnostic_only",
  cross_cov_model="star", heterogeneity_method="chisq", diagnostic_bootstrap_B=0L,
  run_loo_mr=TRUE, run_instrument_sensitivity=FALSE,
  make_html=TRUE, make_figures=TRUE, verbose=FALSE, progress=FALSE
)
manual_comparison <- compare_manual_iv_sets(result_reverse,result_manual)
show_table(manual_comparison[,c("analysis","n_IV","MR","CI_lower","CI_upper",
                                "Q","Q_df","Q_per_df","Q_p")])
show_table(manual_comparison[,c("analysis","warning","inference_status")])
plot_manual_iv_comparison(manual_comparison)
check_report_only(result_reverse,instruments_reverse,"Original five-IV analysis")
check_report_only(result_manual,iv_manual,"Manual two-IV sensitivity")
cat("Milestone 6: the primary fit keeps 5 IVs; the separate manual fit uses 2.\n")


cat("\nTASK 07 | Case A: three exposures, one outcome\n")
case_a <- readRDS("data/case_a_input.rds")
resources_a <- qtl2mr_prepare_resources(
  allele_probs=case_a$allele_probs, kinship=case_a$kinship,
  pheno=case_a$pheno, covar=case_a$covar, marker_map=case_a$marker_map
)
instruments_a <- read.csv("data/instruments.csv")
pairs_a <- data.frame(exposure=c("x1","x2","x3"),outcome="y_A")
show_table(pairs_a)
batch_a <- run_qtl2mr_batch(
  resources_a, pairs_a, instruments_a, output_dir=file.path(run_dir,"case_a_batch"),
  standardize_phenotypes=FALSE, heterogeneity_method="chisq",
  diagnostic_bootstrap_B=0L, n_workers=1L, progress=FALSE
)
show_table(batch_a$main_results[,c("exposure","status","n_loci","theta_reported","reported_ci_lower","reported_ci_upper")])
comparison_a <- compare_three_exposures(case_a,batch_a)
show_table(comparison_a)
plot_three_exposures(comparison_a)
cat("Milestone 7: three separate univariable MR results, alongside observational estimates.\n")


cat("\nTASK 08 | Case B: why founder patterns matter\n")
patterns <- readRDS("data/case_b_patterns.rds")
fit_proportional <- fit_multilocus_founder_mr(
  list(patterns$Proportional), allow_cross_cov=TRUE,
  compute_independence=FALSE, heterogeneity_method="chisq",q_sim_n=0L
)
fit_discordant <- fit_multilocus_founder_mr(
  list(patterns$Discordant), allow_cross_cov=TRUE,
  compute_independence=FALSE, heterogeneity_method="chisq",q_sim_n=0L
)
pattern_results <- compare_patterns(fit_proportional,fit_discordant)
show_table(pattern_results)
plot_patterns(patterns,fit_proportional,fit_discordant)
cat("Milestone 8: similar slopes can hide very different founder-pattern compatibility.\n")


cat("\nTASK 09 | Save and find your results\n")
write.csv(direction,file.path(run_dir,"direction_report.csv"),row.names=FALSE)
write.csv(manual_comparison,file.path(run_dir,"manual_iv_comparison.csv"),row.names=FALSE)
write.csv(comparison_a,file.path(run_dir,"three_exposures.csv"),row.names=FALSE)
write.csv(pattern_results,file.path(run_dir,"founder_patterns.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(run_dir,"sessionInfo.txt"))
saveRDS(list(forward=result_forward,reverse=result_reverse,manual=result_manual,
             manual_comparison=manual_comparison,direction=direction,
             batch=batch_a,comparison=comparison_a,patterns=pattern_results),
        file.path(run_dir,"workshop_results.rds"))
show_saved_files(run_dir)
plot_saved_summary(result_forward,result_reverse,batch_a,pattern_results)
cat("Milestone 9: tables, figures, session information and three native reports are saved.\n")
cat("Finished. Output folder:",run_dir,"\n")
