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
