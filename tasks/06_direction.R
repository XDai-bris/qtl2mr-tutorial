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
