
setwd("C:/Users/ayoung/Desktop/Thesis/real_data")
getwd()
draw_plot <- function(data,name,len_col){

  
  n <- len_col

  species_labels <- ifelse(substr(data$Species, 1, 1) == "1", paste(levels(data$Strain_no), "faecalis", sep = " "), paste(levels(data$Strain_no), "faecium", sep = " "))
  d <- 0.25

  
  pdf(file = paste0("figure_", name, ".pdf"), width = 7, height = 7 * sqrt(2))
  par(mar = c(8, 5.5, 2, 5), yaxs = "i")
  plot.new()
  plot.window(xlim = c(0, 8), c(n + 0.5, 0.5))

  abline(h = seq(0.5, n + 0.5, 1))
  abline(h = c(10.5, 20.5), lwd = 2)
  box()

  # Generate the sequence with log2 steps
  at_values <- 2^(0:8) 
  at_labels <- 2^(-1:8) 

  axis(1, at = seq(-1,8), labels = at_labels, cex.axis = 0.7) # Labels with log2 steps
  axis(2, at = 1:n, labels = with(data, species_labels), las = 1, cex.axis = 0.7)
  title(xlab = "MIC")

  with(data, {

    points(mode.log.MIC, 1:n, pch = 0, cex = 0.7) #open squares - mod MICs
    points(E.log.MIC, 1:n-d, pch = 15, cex = 0.7) #solid square - Predicted mean MICs and
    segments(lower.log.MIC, 1:n-d, upper.log.MIC, 1:n-d)
    points(lower.log.MIC.ref, 1:n+d, col = 1, cex = 0.7) #open circles - lower boundary accounting for interval censoring
    points(upper.log.MIC.ref, 1:n+d, col = 1, pch = 16, cex = 0.7) # solid circles - reference MICs
 
    segments(lower.log.MIC.ref, 1:n+d, upper.log.MIC.ref, 1:n+d, col = 1)
  })

  par(xpd = TRUE)
  par(new = TRUE)
  par(fig = c(0, 1, 0, 0.2), mar = c(2, 2, 2, 2))
  legend("bottomright", legend = c("Mode MICs", "mean MICs", "interval censoring", " reference MICs"),
         pch = c(0, 15, 1, 16, 1), col = c("black", "black", "black", "black"),
         cex = 0.7, bty = "n")
  dev.off()

# 
#   pdf(file = paste0("figure_",name,".pdf"), width = 7, height = 7 * sqrt(2))
#   par(mar = c(8,5.5,2, 5), yaxs = "i");
#   plot.new();
#   plot.window(xlim = c(0,8), c(n+0.5, 0.5))
# 
#   abline(h = seq(0.5, n+0.5, 1), col = 8);
#   abline(h = c(10.5, 20.5), lwd = 2); box();
#   axis(1, at = seq(0,8, 2), labels = signif(2^seq(0, 8, 2), 3), cex.axis = 0.7)
#   axis(2, at = 1:n, labels = with(data, species_labels), las = 1, cex.axis = 0.7)
#   title(xlab = "MIC")
# 
#   with(data, {
#     points(mode.log.MIC, 1:n, pch = 0, cex = 0.7) #open squares - mod MICs
#     points(E.log.MIC, 1:n-d, pch = 15, cex = 0.7) #solid square - Predicted mean MICs and
#     segments(lower.log.MIC, 1:n-d, upper.log.MIC, 1:n-d)
#     points(lower.log.MIC.ref, 1:n+d, col = 1, cex = 0.7) #open circles - lower boundary accounting for interval censoring
#     points(upper.log.MIC.ref, 1:n+d, col = 1, pch = 16, cex = 0.7) # solid circles - reference MICs
#     segments(lower.log.MIC.ref, 1:n+d, upper.log.MIC.ref, 1:n+d, col = 1)
#   })
#   par(xpd=TRUE)
#   par(new=T)
#   par(fig=c(0, 1, 0, 0.2), mar=c(2,2,2,2))
#   #lower boundary accounting for
#   legend("bottomright", legend = c("Mode MICs", "mean MICs", "interval censoring", " reference MICs"),
#          pch = c(0, 15, 1, 16, 1), col = c("black", "black", "black", "black"),
#          cex = 0.7, bty = "n")
#   dev.off()


  system(paste("open",paste0("figure_",name,".pdf")))
  
}


