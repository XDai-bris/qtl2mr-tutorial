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
