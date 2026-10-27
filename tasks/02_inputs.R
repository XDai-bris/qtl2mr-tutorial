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
