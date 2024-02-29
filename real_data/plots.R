
setwd("C:/Users/ayoung/Desktop/Thesis/real_data")
getwd()
draw_plot <- function(data,name){
  # 
  # # # labs compared to mean
  # pdf(file = "figure_3.pdf", width = 7, height = 7)
  # par(mar = c(4.5, 4.5, 0.5, 0.5), yaxs = "i")
  # plot.new()
  # plot.window(xlim = c(-2, 2), c(n.lab + 1, 0))
  # abline(h = 1:n.lab, col = 8, lty = 3)
  # box()
  # axis(1, at = seq(-2, 2, 1), labels = 2^seq(-2, 2, 1))
  # axis(2, at = 1:n.lab, labels = lab.data$Strain_no)
  # abline(v = 0, lty = 2)
  # title(xlab = "Fold difference in MIC compared to the expected mean", ylab = "Laboratory number")
  # with(lab.data, {
  #   points(diff.log.MIC, 1:n.lab, pch = 15)
  #   segments(lower.diff.log.MIC, 1:n.lab, upper.diff.log.MIC, 1:n.lab)
  # })
  # dev.off()
  # system(paste("open", "figure_3.pdf"))
  # 
  
  n <- ncol(X.samplepred)
  d <- 0.25
  
  
  pdf(file = paste0("figure_",name,".pdf"), width = 7, height = 7 * sqrt(2))
  par(mar = c(8,5.5,2, 5), yaxs = "i");
  plot.new();
  plot.window(xlim = c(-10, 8), c(n+0.5, 0.5))
  
  abline(h = seq(0.5, n+0.5, 1), col = 8);
  abline(h = c(10.5, 20.5), lwd = 2); box();
  axis(1, at = seq(-9, 7, 2), labels = signif(2^seq(-9, 7, 2), 3), cex.axis = 0.7)
  axis(2, at = 1:n, labels = with(data.newdata, levels(interaction(Strain_no, sep = " "))), las = 1, cex.axis = 0.7)
  title(xlab = "MIC")
  
  with(data.newdata, {
    points(mode.log.MIC, 1:n, pch = 0, cex = 0.7) #open squares - mod MICs
    points(E.log.MIC, 1:n-d, pch = 15, cex = 0.7) #solid square - Predicted mean MICs and 
    segments(lower.log.MIC, 1:n-d, upper.log.MIC, 1:n-d) 
    points(lower.log.MIC.ref, 1:n+d, col = 1, cex = 0.7) #open circles - lower boundary accounting for interval censoring
    points(upper.log.MIC.ref, 1:n+d, col = 1, pch = 16, cex = 0.7) # solid circles - reference MICs
    segments(lower.log.MIC.ref, 1:n+d, upper.log.MIC.ref, 1:n+d, col = 1)
  })
  par(xpd=TRUE)
  par(new=T)
  par(fig=c(0, 1, 0, 0.2), mar=c(2,2,2,2))
  #lower boundary accounting for 
  legend("bottomright", legend = c("Mode MICs", "mean MICs", "interval censoring", " reference MICs"),
         pch = c(0, 15, 1, 16, 1), col = c("black", "black", "black", "black"),
         cex = 0.7, bty = "n")
  dev.off()
  system(paste("open", "figure_4.pdf"))
  
}


