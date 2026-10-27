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
