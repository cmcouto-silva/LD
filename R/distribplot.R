distribplot <- function(ihs.data, out, file.type = 'png', 
         lty = 1, 
         lwd = 1.5, 
         col = c("blue", "red"),
         main = "Genome-wide distribution",
         xlab = "iHS value",
         cex.main = 1.5,
         cex.lab = 1.25,
         qqplot = TRUE) {
  
  substrRight <- function(x, n) {
    substr(x, nchar(x)-n+1, nchar(x))
    }
  
  if(file.type == 'png') {
    out <- ifelse(substrRight(out, 4L) == ".png", paste0(out, ".ihsplot"), paste0(out, ".ihsplot.png"))
    png(out, width = 800, height = 600, pointsize = 16)
  }
  
  if(file.type == 'svg') {
    out <- ifelse(substrRight(out, 4L) == ".svg", paste0(out, ".ihsplot"), paste0(out, ".ihsplot.svg"))
    svg(out, width = 800, height = 600, pointsize = 16)
  }
  
  if(file.type == 'pdf') {
    out <- ifelse(substrRight(out, 4L) == ".pdf", paste0(out, ".ihsplot"), paste0(out, ".ihsplot.pdf"))
    pdf(out, width = 800, height = 600, pointsize = 16)
  }
  
  par(mar = c(5, 5, 4, 2) + 0.1)
  
  plot(density(ihs.data, na.rm = TRUE), main = main, xlab = xlab, col = col[1], 
       lty = lty, lwd = lwd, cex.main = cex.main, cex.lab = cex.lab)
  
  curve(dnorm, col = col[2], add = TRUE, bty='n')
  
  legend("topright", c("Observed", "Gaussian"), bty = "n", 
         col = col, lty = lty, lwd = lwd)
  
  dev.off()
  
  if(qqplot) {
    
    par(mar = c(5, 5, 4, 2) + 0.1)
    
    if(file.type == 'png') {
      out <- ifelse(substrRight(out, 4L) == ".png", paste0(out, ".ihsplot"), paste0(out, ".ihsplot.png"))
      png(out, width = 800, height = 600, pointsize = 16)
    }
    
    if(file.type == 'svg') {
      out <- ifelse(substrRight(out, 4L) == ".svg", paste0(out, ".ihsplot"), paste0(out, ".ihsplot.svg"))
      svg(out, width = 800, height = 600, pointsize = 16)
    }
    
    if(file.type == 'pdf') {
      out <- ifelse(substrRight(out, 4L) == ".pdf", paste0(out, ".ihsplot"), paste0(out, ".ihsplot.pdf"))
      pdf(out, width = 800, height = 600, pointsize = 16)
    }
    
    qqnorm(ihs.data[!is.na(ihs.data)], cex.main = cex.main, cex.lab = cex.lab, pch = 16, cex = 0.75)
    abline(a = 0, b = 1, lty = 2)
    dev.off()
  }
  
}
