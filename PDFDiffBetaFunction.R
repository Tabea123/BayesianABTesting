################################################################################
##                                                                            ##
##   Analytical computation of the difference between two beta functions      ##
##                                                                            ##
################################################################################

# this function is based on a formula published in
# Pham-Gia, T., Turkkan, N., & Eng, P. (1993). Bayesian analysis of the 
# difference of two proportions. Communications in Statistics-Theory and Methods, 
# 22, 1755-1771.

# p = p1 - p2

pdf.diff <- function(a1, b1, a2, b2){
  
  logA <- lbeta(a1, b1) + lbeta(a2, b2)
  
  # for 0 < p =< 1
  smaller1 <- function(p){
    exp(lbeta(a2, b1) + (b1+b2-1)*log(p) + (a2+b1-1)*log(1-p) +
          log(appell.F1(b1, a1+b1+a2+b2-2, 1-a1, b1+a2, (1-p), (1-p^2))) - logA)
    
  }
  
  ps   <- seq(0.003, 1, 0.001)

  par(cex.main = 1.5,  mar = c(5, 5, 3, 3) + 0.1, mgp = c(3.5, 1, 0), 
      cex.lab = 1.5, font.lab = 2, cex.axis = 1.6, bty = "n", las = 1)
  
  plotting <- function(ps){
    
    plot(x = ps, y = mapply(smaller1, ps), xlim = c(-1, 1), ylim = c(0, 6),
         type = "l", bty = "n", lwd = 2, xlab = "", ylab = "", axes = F,
         cex.lab = 1.5, cex.axis = 1.8, cex.main = 1.5, las = 1, bg = "grey")
    
    axis(1, at = seq(-1, 1, 0.25), 
         labels = seq(-1, 1, 0.25), 
         lwd = 2, lwd.ticks = 2, line = -0.1)
    axis(2, at = seq(0, 6, 1),  
         lwd = 2, lwd.ticks = 2, line = -0.2)
    
    mtext(expression(paste("Difference", ~delta)), 
          side = 1, line = 3, las = 1, cex = 2, font = 0.2, adj = 0.5)
    mtext("Density", side = 2, line = 2.5, cex = 2, font = 2, las = 0)
    
  }
  
  plot.smaller1 <- plotting(ps)
  
  # for -1 <= p < 0
  smaller0 <- function(p){
    exp(lbeta(a1, b2) + (b1+b2-1)*log(-p) + (a1+b2-1)*log(1+p) +
          log(appell.F1(b2, 1-a2, a1+a2+b1+b2-2, a1+b2, (1-p^2), (1+p))) - logA)
  }
  
  ps <- seq(-1, -0.005, 0.001)
  
  plot.smaller0 <- lines(x = ps, y = mapply(smaller0, ps), 
                         xlim = c(-1, 0), type = "l", lwd = 2)
  
  # p = 0
  point.zero <- points(x = 0,
                       y = exp(lbeta(a1+a2-1, b1 + b2 -1) - logA),
                       pch = 20, cex = 1.5)
  
  return(list(plot.smaller1, plot.smaller0, point.zero))
}



# combined function for exact difference
pdf.diff2 <- function(p){
  
  results <- numeric(length(p))
  
  logA <- lbeta(a1, b1) + lbeta(a2, b2)
  
  for(i in 1:length(p)){
    if(p[i] > 0 && p[i] <= 1){
      results[i] <- (exp(lbeta(a2, b1) + (b1+b2-1)*log(p[i]) + 
                           (a2+b1-1)*log(1-p[i]) + log(appell.F1(b1, a1+b1+a2+b2-2, 1-a1, b1+a2, 
                                                                 (1-p[i]), (1-p[i]^2))) - logA))
    } else if(p[i] <= -1 && p[i] < 0){
      results[i] <- (exp(lbeta(a1, b2) + (b1+b2-1)*log(-p[i]) + (a1+b2-1)*log(1+p[i]) +
                           log(appell.F1(b2, 1-a2, a1+a2+b1+b2-2, a1+b2, (1-p[i]^2), (1+p[i]))) - logA))
    } else {
      results[i] <- (exp(lbeta(a1+a2-1, b1 + b2 -1) - logA))
    }
  }
  
  return(results)
  
}