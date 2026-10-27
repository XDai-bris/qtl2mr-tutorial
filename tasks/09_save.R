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
