# Workshop display helpers. Scientific calls remain visible in tasks/.
founder_colours <- function() {
  c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#333333")
}
plot_install_check <- function(example) {
  draw_workshop_figure("00_install_check.png",function() {
    p <- example$allele_probs[[1]][1:6,,1]
    par(mar=c(4,4,3,1))
    barplot(t(p),col=founder_colours(),border=NA,ylim=c(0,1.2),
            xlab="Mouse in the built-in package example",ylab="Founder allele probability",
            main="Installation check: the bundled example is readable")
    legend("top",legend=colnames(p),fill=founder_colours(),ncol=8,bty="n",cex=.8)
  })
}
plot_single_pair_inputs <- function(demo) {
  ids <- rownames(demo$Forward)
  stopifnot(length(ids)==400,identical(ids,rownames(demo$Reverse)),
            identical(ids,rownames(demo$covar)),identical(ids,rownames(demo$kinship)))
  draw_workshop_figure("01_single_pair_inputs.png",function() {
    layout(matrix(1:3,nrow=1),widths=c(1,1,.18))
    par(cex=1,mar=c(5,4,3,1))
    p <- demo$allele_probs[[1]][1:10,,40]
    barplot(t(p),col=founder_colours(),border=NA,ylim=c(0,1.22),
            xlab="Mouse (first ten)",ylab="Founder allele probability",main="One locus, eight founders")
    legend("top",legend=colnames(p),fill=founder_colours(),ncol=4,bty="n",cex=.75)
    par(mar=c(5,4,3,1))
    colours<-hcl.colors(64,"YlGnBu");limits<-c(-.1,.5)
    image(1:20,1:20,demo$kinship[1:20,1:20],col=colours,zlim=limits,
          xlab="Mouse",ylab="Mouse",main="Kinship K (first 20 mice)")
    par(mar=c(5,0,3,3))
    values<-seq(limits[1],limits[2],length.out=64)
    image(1,values,matrix(values,nrow=1),col=colours,zlim=limits,
          axes=FALSE,xlab="",ylab="")
    axis(4,at=c(-.1,0,.25,.5),las=1,cex.axis=.8);box()
    mtext("K",side=3,line=1)
  },width=1800)
}
regional_islands <- function(marker_map,scan,threshold=3) {
  z <- marker_map;z$lod <- as.numeric(scan[z$marker,1]);out<-list()
  for(ch in unique(z$chr)) {
    a<-z[z$chr==ch,];r<-rle(is.finite(a$lod)&a$lod>=threshold)
    ends<-cumsum(r$lengths);starts<-ends-r$lengths+1
    for(i in which(r$values)) {
      ix<-starts[i]:ends[i];k<-ix[which.max(a$lod[ix])]
      out[[length(out)+1]]<-data.frame(chr=ch,start=a$pos[starts[i]],end=a$pos[ends[i]],
        n_markers=length(ix),marker=a$marker[k],pos=a$pos[k],lod=a$lod[k])
    }
  }
  do.call(rbind,out)
}
plot_regional_selection <- function(marker_map,scan,peaks) {
  regions<-regional_islands(marker_map,scan)
  stopifnot(setequal(regions$marker,peaks$marker))
  first_chr<-unique(marker_map$chr)[unique(marker_map$chr)%in%regions$chr][1]
  zoom<-regions[regions$chr==first_chr,];zoom<-zoom[order(-zoom$n_markers,zoom$start),][1,]
  cat("Zoom:",zoom$n_markers,"markers in the shaded island ->",zoom$marker,"\n")
  draw_workshop_figure("02_regional_peaks.png",function() {
    par(mfrow=c(2,1),mar=c(4,4,3,1))
    z<-marker_map;z$lod<-as.numeric(scan[z$marker,1]);xx<-(as.integer(z$chr)-1)*85+z$pos
    plot(xx,z$lod,type="n",xaxt="n",xlab="Chromosome",ylab="Exposure LOD",main="1,520 scan markers -> 37 candidate regions")
    for(ch in unique(z$chr)) {ix<-which(z$chr==ch);lines(xx[ix],z$lod[ix],col="#486575")}
    axis(1,at=(0:18)*85+40,labels=1:19);abline(h=3,lty=2,col="#bd693a")
    ii<-match(peaks$marker,z$marker);points(xx[ii],z$lod[ii],col="#087f8c",pch=1)
    a<-z[z$chr==first_chr,]
    plot(a$pos,a$lod,type="n",xlab="Chromosome 1 position (synthetic cM)",ylab="Exposure LOD",
         main=paste(zoom$n_markers,"markers -> one regional peak"),ylim=c(0,max(a$lod)*1.15))
    rect(zoom$start-.5,0,zoom$end+.5,max(a$lod)*1.15,col="#e4f2ef",border=NA)
    lines(a$pos,a$lod,col="#486575");points(a$pos,a$lod,pch=16,cex=.5,col="#486575")
    abline(h=3,lty=2,col="#bd693a");points(zoom$pos,zoom$lod,pch=18,col="#087f8c",cex=1.5)
    text(zoom$pos,zoom$lod,zoom$marker,pos=3,cex=.8)
    extra<-regions[regions$chr==first_chr & regions$marker!=zoom$marker,]
    for(i in seq_len(nrow(extra))) {
      rect(extra$start[i]-.5,0,extra$end[i]+.5,3,
           col=adjustcolor("#c87936",alpha.f=.2),border=NA)
    }
    points(extra$pos,extra$lod,pch=18,col="#ad6329",cex=1.3)
    label_y<-max(a$lod)*c(.35,.47)[seq_len(nrow(extra))]
    segments(extra$pos,extra$lod+.4,extra$pos,label_y-.6,col="#ad6329")
    text(extra$pos,label_y,extra$marker,cex=.72,col="#ad6329")
    legend("topright",c("LOD = 3","Peak in highlighted island","Peaks in other islands"),
           lty=c(2,NA,NA),pch=c(NA,18,18),
           col=c("#bd693a","#087f8c","#ad6329"),bty="n",cex=.75)
  },height=1100)
}
plot_marginal_strength <- function(instruments) {
  draw_workshop_figure("03_marginal_F.png",function() {
    par(mar=c(5,4,3,1));z<-instruments
    plot(seq_len(nrow(z)),z$F_X_block,pch=ifelse(z$selected,19,1),
         col=ifelse(z$selected,"#087f8c","#c87936"),cex=1.2,
         xlab="Candidate peak (ordered by exposure LOD)",ylab="Marginal founder F",
         main=paste(sum(z$selected),"retained from",nrow(z),"candidate peaks"),ylim=c(0,max(z$F_X_block)*1.12))
    abline(h=4,lty=2,col="#52697a")
    legend("topright",c("Retained","Rejected by F","F = 4"),pch=c(19,1,NA),lty=c(NA,NA,2),
           col=c("#087f8c","#c87936","#52697a"),bty="n")
  })
}
diagnostic_overview <- function(result) {
  s<-summary(result)
  data.frame(estimate=s$estimate,lower=s$ci_lower,upper=s$ci_upper,Q=s$Q,Q_p=s$Q_p,
             warning=ifelse(nzchar(s$warnings),s$warnings,"none"))
}
plot_forward_result <- function(result) {
  draw_workshop_figure("04_forward_mr.png",function() {
    s<-summary(result);par(mar=c(5,4,3,1))
    plot(NA,xlim=range(c(0,.9,s$ci_lower,s$ci_upper))+c(-.05,.1),ylim=c(.6,1.5),
         yaxt="n",ylab="",xlab="Y units per X unit",main="Forward: MR estimate and profile interval")
    abline(v=c(0,.9),lty=c(3,2),col=c("grey70","#c87936"))
    segments(s$ci_lower,1,s$ci_upper,1,lwd=3,col="#087f8c");points(s$estimate,1,pch=19,cex=1.6,col="#087f8c")
    text(s$estimate,1.13,sprintf("%.3f [%.3f, %.3f]",s$estimate,s$ci_lower,s$ci_upper))
    legend("topright","Dashed: structural effect 0.900",lty=2,col="#c87936",bty="n",cex=.85)
  },height=650)
}
plot_forward_influence <- function(result,loo) {
  draw_workshop_figure("05_forward_LOO.png",function() {
    par(mar=c(5,8,3,1));n<-nrow(loo);yy<-rev(seq_len(n));s<-summary(result)
    lim<-range(loo$profile_ci_lower,loo$profile_ci_upper,s$estimate)
    plot(NA,xlim=lim,ylim=c(.5,n+.7),yaxt="n",ylab="",xlab="MR estimate after omitting one locus",
         main="Leave one locus out: influence is not a deletion rule")
    axis(2,at=yy,labels=loo$marker,las=1,cex.axis=.85)
    abline(v=s$estimate,lty=2,col="grey50")
    segments(loo$profile_ci_lower,yy,loo$profile_ci_upper,yy,col="#087f8c",lwd=2)
    points(loo$theta_minus_m,yy,col="#087f8c",pch=19)
  },height=650)
}
plot_pair_comparison <- function(forward,reverse) {
  draw_workshop_figure("06_forward_reverse_mr.png",function() {
    s<-rbind(summary(forward),summary(reverse));par(mar=c(5,7,3,1))
    plot(NA,xlim=range(c(0,.9,s$ci_lower,s$ci_upper))+c(-.06,.1),ylim=c(.5,2.7),yaxt="n",
         ylab="",xlab="X-to-Y slope (original units)",main="Similar pooled slopes, very different compatibility")
    axis(2,at=2:1,labels=c("Forward","Reverse"),las=1);abline(v=0,lty=3,col="grey60")
    segments(s$ci_lower,2:1,s$ci_upper,2:1,lwd=3,col=c("#087f8c","#c87936"))
    points(s$estimate,2:1,pch=19,cex=1.4,col=c("#087f8c","#c87936"))
    points(c(.9,0),2:1,pch=4,cex=1.4)
    legend("topright","Cross: specified structural X-to-Y effect",pch=4,bty="n",cex=.8)
  },height=700)
}
direction_table <- function(forward,reverse) {
  do.call(rbind,lapply(c("Forward","Reverse"),function(s) {
    z<-if(s=="Forward") forward else reverse
    data.frame(scenario=s,marker=z$marker,F_X=z$F_X_block,R2_X=z$partial_R2_X_block,R2_Y=z$partial_R2_Y_block,
               direction_report=ifelse(is.na(z$steiger_pass),"Missing",ifelse(z$steiger_pass,"Exposure-dominant","Outcome-dominant")))
  }))
}
check_report_only <- function(result,instruments,label) {
  keep<-as.character(instruments$marker[instruments$selected])
  actual<-as.character(result$pilot$selected_loci$marker)
  stopifnot(identical(keep,actual),nrow(result$pilot$dropped_loci)==0,
            is.null(result$steiger_sensitivity),is.null(result$instrument_sensitivity))
  cat(label,": supplied IV set unchanged = TRUE; retained =",length(keep),"; automatic SF subset refit = FALSE\n")
}
plot_direction <- function(forward,reverse) {
  draw_workshop_figure("07a_SF_report.png",function() {
    par(mfrow=c(1,2),mar=c(4.5,4.5,3,1),oma=c(2,0,0,0),pty="s")
    lim<-max(forward$partial_R2_X_block,forward$partial_R2_Y_block,reverse$partial_R2_X_block,reverse$partial_R2_Y_block)*1.2
    for(s in c("Forward","Reverse")) {
      z<-if(s=="Forward")forward else reverse
      plot(NA,xlim=c(0,lim),ylim=c(0,lim),xaxs="i",yaxs="i",asp=1,
           xlab="Association with X (partial R-squared)",ylab="Association with Y (partial R-squared)",
           main=paste0(s," (",nrow(z)," IVs)"))
      polygon(c(0,0,lim),c(0,lim,lim),col="#edf0f5",border=NA);abline(0,1,lty=2,col="grey55")
      cc<-ifelse(as.integer(z$chr)<=9,"#3b67b2","#087f8c")
      points(z$partial_R2_X_block,z$partial_R2_Y_block,col=cc,pch=19,cex=1.2)
      pp<-ifelse(z$steiger_pass,4,2);pp[!z$steiger_pass&z$chr==2]<-3;pp[!z$steiger_pass&z$chr==3]<-4
      text(z$partial_R2_X_block,z$partial_R2_Y_block,paste("Chr",z$chr),pos=pp,cex=.72)
      labs<-if(s=="Forward")"Source region" else c("Source region","X-local region")
      legend("bottomright",labs,col=c("#3b67b2","#087f8c")[seq_along(labs)],pch=19,bty="n",cex=.7)
    }
    mtext("Above diagonal: Y association is larger. Colours show the known simulation mechanism.",side=1,outer=TRUE,line=.5,cex=.85)
  },height=780)
}
compare_manual_iv_sets <- function(full,manual) {
  one_row <- function(result,label) {
    s<-summary(result);fit<-result$pilot$fit
    data.frame(analysis=label,n_IV=s$n_loci,MR=s$estimate,
      CI_lower=s$ci_lower,CI_upper=s$ci_upper,Q=s$Q,
      Q_df=fit$Q_chisq_df,Q_per_df=s$Q/fit$Q_chisq_df,Q_p=s$Q_p,
      warning=ifelse(nzchar(s$warnings),s$warnings,"none"),
      confidence_set=s$confidence_set_type,primary_reliable=s$primary_reliable,
      inference_status=fit$inference_status)
  }
  rbind(one_row(full,"All five IVs"),one_row(manual,"Manual Chr10 + Chr11"))
}
plot_manual_iv_comparison <- function(z) {
  draw_workshop_figure("07b_manual_iv_comparison.png",function() {
    cc<-c("#b36c36","#087f8c")
    layout(matrix(c(1,2),nrow=2),heights=c(1.4,1))
    par(mar=c(4,10,3,2))
    xlim<-range(c(0,z$CI_lower,z$CI_upper))+c(-.07,.07)
    plot(NA,xlim=xlim,ylim=c(.55,2.65),yaxt="n",ylab="",
         xlab="X-to-Y MR estimate (95% profile interval, original units)",
         main="Frozen Reverse data: explicitly changing the IV set")
    axis(2,at=2:1,labels=c("All five IVs","Manual Chr10 + Chr11"),las=1,cex.axis=.85)
    abline(v=0,lty=2,col="grey55")
    segments(z$CI_lower,2:1,z$CI_upper,2:1,lwd=3,col=cc)
    points(z$MR,2:1,pch=19,cex=1.4,col=cc)
    text(z$MR,c(2.25,1.25),sprintf("%.3f [%.3f, %.3f]",z$MR,z$CI_lower,z$CI_upper),cex=.85)
    legend("topright",legend="Dashed: specified X-to-Y effect = 0",lty=2,col="grey55",bty="n",cex=.8)
    par(mar=c(1,1,1,1));plot.new();plot.window(xlim=c(0,1),ylim=c(0,1))
    xx<-c(.06,.46,.61,.77,.90)
    text(xx,.88,c("Analysis","Q","df","Q / df","Q p-value"),adj=c(0,.5),font=2,cex=.9)
    segments(.05,.76,.98,.76,col="#d6e2e4")
    for(i in 1:2) {
      yy<-c(.62,.42)[i];p<-if(z$Q_p[i]<.001)"<0.001" else sprintf("%.3f",z$Q_p[i])
      vals<-c(z$analysis[i],sprintf("%.3f",z$Q[i]),as.character(z$Q_df[i]),sprintf("%.3f",z$Q_per_df[i]),p)
      text(xx,yy,vals,adj=c(0,.5),cex=.9,col=cc[i])
    }
    text(.05,.17,"Same data and fitting settings. SF reports only; the two-IV set is specified by the user.",adj=0,cex=.78)
    text(.05,.05,"Both fits: hsq_boundary; nominal profile intervals; inference_status = not_calibrated.",adj=0,cex=.78)
  },width=1800,height=1050)
}
show_saved_files <- function(run_dir) {
  paths<-c("direction_report.csv","manual_iv_comparison.csv","three_exposures.csv","founder_patterns.csv","sessionInfo.txt",
           "workshop_results.rds","forward/report.html","reverse/report.html","reverse_manual/report.html")
  stopifnot(all(file.exists(file.path(run_dir,paths))))
  show_table(data.frame(file=paths,saved=TRUE))
}
plot_saved_summary <- function(forward,reverse,batch,patterns) {
  draw_workshop_figure("10_take_home.png",function() {
    par(mar=c(2,2,3,2));plot.new();title("Your take-home analysis record")
    text(.05,.82,"1. Forward + Reverse",adj=0,font=2);text(.05,.72,"Read SF; compare the primary and manually specified IV sets.",adj=0)
    text(.05,.55,"2. Case A: three exposures",adj=0,font=2);text(.05,.45,"Separate univariable MR fits; compare with observational results.",adj=0)
    text(.05,.28,"3. Case B: founder patterns",adj=0,font=2);text(.05,.18,"Similar slopes can conceal incompatible founder effects.",adj=0)
  },height=650)
}

read_pair_table <- function(result, name) {
  read.csv(file.path(result$output_dir, "data", paste0(name, ".csv")))
}

draw_workshop_figure <- function(filename, draw, width = 1500, height = 900) {
  # Draw in the Plots pane first. If that fails (for example "figure margins
  # too large" in a small Plots pane), carry on: the PNG file below is still saved.
  if (interactive()) tryCatch(draw(), error = function(e) {
    cat("Could not draw this figure in the Plots pane (", conditionMessage(e), ").\n",
        "It will still be saved as ", file.path(run_dir, "figures", filename),
        "; enlarge the Plots pane and run that plot line again to see it there.\n", sep = "")
  })
  png(file.path(run_dir, "figures", filename), width = width, height = height, res = 150)
  tryCatch(draw(), finally = dev.off())
  cat("Figure saved:", file.path("figures", filename), "\n")
  invisible(filename)
}

compare_three_exposures <- function(d, batch) {
  exposures <- c("x1", "x2", "x3")
  obs <- do.call(rbind, lapply(exposures, function(x) {
    fit <- lm(reformulate(c(x, "batch", "age"), response = "y_A"),
              data = cbind(d$pheno, d$covar))
    ci <- confint(fit, parm = x)
    data.frame(exposure = x, observational = unname(coef(fit)[x]),
               obs_lower = ci[1], obs_upper = ci[2])
  }))
  main <- batch$main_results[match(exposures, batch$main_results$exposure), ]
  data.frame(obs, structural_theta = c(-0.25, 0.25, 0), MR = main$theta_reported,
             MR_lower = main$reported_ci_lower, MR_upper = main$reported_ci_upper)
}

plot_three_exposures <- function(x) {
  draw_workshop_figure("08_case_a_batch.png", function() {
    par(mar = c(5, 5, 3, 1))
    lim <- range(x[, c("obs_lower", "obs_upper", "MR_lower", "MR_upper")])
    plot(NA, xlim = lim + c(-0.03, 0.03), ylim = c(0.5, 3.9), yaxt = "n",
         xlab = "Y units per X unit", ylab = "", main = "Three exposures, one outcome")
    axis(2, at = 3:1, labels = toupper(x$exposure), las = 1)
    abline(v = 0, lty = 3, col = "grey60")
    for (i in 1:3) {
      yy <- 4 - i
      segments(x$obs_lower[i], yy + 0.12, x$obs_upper[i], yy + 0.12,
               col = "#c87936", lwd = 2)
      points(x$observational[i], yy + 0.12, col = "#c87936", pch = 1, cex = 1.2)
      segments(x$MR_lower[i], yy - 0.12, x$MR_upper[i], yy - 0.12,
               col = "#087f8c", lwd = 2)
      points(x$MR[i], yy - 0.12, col = "#087f8c", pch = 16, cex = 1.2)
      points(x$structural_theta[i], yy, pch = 3, cex = 1.3)
    }
    legend("topright", c("Observational", "qtl2mr", "Specified structural effect"),
           col = c("#c87936", "#087f8c", "black"), pch = c(1, 16, 3),
           bty = "n", cex = 0.85)
  })
}

compare_patterns <- function(a, b) {
  data.frame(pattern = c("Proportional", "Discordant"),
             estimate = c(a$theta, b$theta),
             lower = c(a$profile_ci_lower, b$profile_ci_lower),
             upper = c(a$profile_ci_upper, b$profile_ci_upper),
             Q = c(a$Q, b$Q), Q_p = c(a$Q_p_chisq, b$Q_p_chisq))
}

plot_patterns <- function(patterns, a, b) {
  fits <- list(Proportional = a, Discordant = b)
  draw_workshop_figure("09_founder_patterns.png", function() {
    par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
    for (nm in names(patterns)) {
      block <- patterns[[nm]]
      H <- block$founder_H
      ax <- drop(H %*% block$beta_x); ay <- drop(H %*% block$beta_y)
      plot(ax, ay, xlim = c(-2.4, 2.4), ylim = c(-1.6, 1.8), pch = 16,
           col = "#087f8c", xlab = "Founder effect on X", ylab = "Founder effect on Y",
           main = nm)
      abline(a = 0, b = fits[[nm]]$theta, col = "#b8484a", lwd = 2)
      text(ax, ay, rownames(H), pos = 3, cex = 0.8)
    }
    for (nm in names(patterns)) {
      block <- patterns[[nm]]; theta <- fits[[nm]]$theta; H <- block$founder_H
      r <- block$beta_y - theta * block$beta_x
      V <- block$Sigma_y + theta^2 * block$Sigma_x -
        theta * (block$Sigma_xy + t(block$Sigma_xy))
      e <- drop(H %*% r)
      score <- e / sqrt(diag(H %*% V %*% t(H)))
      barplot(score, names.arg = rownames(H), ylim = c(-7, 7), border = NA,
              col = ifelse(score >= 0, "#b8484a", "#087f8c"),
              main = paste(nm, "residuals"), xlab = "Founder",
              ylab = "Residual / marginal SD")
      abline(h = 0, col = "grey50")
    }
  }, width = 1500, height = 1200)
}
