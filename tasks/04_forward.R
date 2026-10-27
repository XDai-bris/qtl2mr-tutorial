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
