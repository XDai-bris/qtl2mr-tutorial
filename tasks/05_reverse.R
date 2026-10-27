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
