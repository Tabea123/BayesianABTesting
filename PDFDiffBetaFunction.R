################################################################################
##                                                                            ##
##   Analytical computation of the difference between two beta functions      ##
##                                                                            ##
################################################################################

# this function is based on a formula published in
# Pham-Gia, T., Turkkan, N., & Eng, P. (1993). Bayesian analysis of the 
# difference of two proportions. Communications in Statistics-Theory and Methods, 
# 22, 1755-1771.


pdf.diff <- function(a1, b1, a2, b2){
  
  logA <- lbeta(a1, b1) + lbeta(a2, b2)
  
  # for 0 < p =< 1
  smaller1 <- function(p){
    exp(lbeta(a2, b1) + (b1+b2-1)*log(p) + (a2+b1-1)*log(1-p) +
          log(appell.F1(b1, a1+b1+a2+b2-2, 1-a1, b1+a2, (1-p), (1-p^2))) - logA)
    
  }
  
  ps   <- seq(0.003, 1, 0.001)
  
  plotting <- function(ps){
    plot(x = ps,
         y = mapply(smaller1, ps), 
         xlim = c(-1, 1), ylim = c(0, 6),
         type = "l", bty = "n", lwd = 2,
         xlab = "", ylab = "", axes = F)
    
    axis1 <- axis(1)
    axis2 <- axis(2)
    
    mtext(expression(paste("Difference", ~delta)), side = 1, line = 2.5, cex = 1.5)
    mtext("Density", side = 2, line = 2.7, cex = 1.5)
  }
  
  plot.smaller1 <- plotting(ps)
  
  # for -1 <= p < 0
  smaller0 <- function(p){
    exp(lbeta(a1, b2) + (b1+b2-1)*log(-p) + (a1+b2-1)*log(1+p) +
          log(appell.F1(b2, 1-a2, a1+a2+b1+b2-2, a1+b2, (1-p^2), (1+p))) - logA)
  }
  
  ps <- seq(-1, -0.005, 0.001)
  
  plot.smaller0 <- lines(x = ps,
                         y = mapply(smaller0, ps), 
                         xlim = c(-1, 0),
                         type = "l", lwd = 2)
  
  # p = 0
  point.zero <- points(x = 0,
                       y = exp(lbeta(a1+a2-1, b1 + b2 -1) - logA),
                       pch = 20, cex = 1.5)
  
  return(list(plot.smaller1, plot.smaller0, point.zero))
}
