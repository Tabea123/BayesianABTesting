################################################################################
##                                                                            ##
##                      Comparison P(diff > 0)                                ##
##                                                                            ##
################################################################################

# note: run AppellF1Function.R, and PDFDiffBetaFunction.R before running this script

a1 <- 1
b1 <- 1
a2 <- 1
b2 <- 1

### exact diff
delta.exact  <- numeric(length(data$y1))
for(i in 1:length(data$y1)){
  
  a2 <- a2 + data$y1[i]
  b2 <- b2 + ifelse(data$y1 == 1, 0, 1)[i]
  a1 <- a1 + data$y2[i]
  b1 <- b1 + ifelse(data$y2 == 1, 0, 1)[i]
  
  delta.exact[i] <- integrate(pdf.diff2, 0, 1, stop.on.error = F)$value
}

a1 <- 1
b1 <- 1
a2 <- 1
b2 <- 1

### Normal Approximation
delta.approx <- numeric(length(data$y1))
for(i in 1:length(data$y1)){
  
a1 <- a1 + data$y1[i]
b1 <- b1 + ifelse(data$y1 == 1, 0, 1)[i]
a2 <- a2 + data$y2[i]
b2 <- b2 + ifelse(data$y2 == 1, 0, 1)[i]

mu1       <- a1/(a1+b1)
sigma1    <- sqrt((a1*b1)/((a1+b1)^2*(a1+b1+1)))
mu2       <- a2/(a2+b2)
sigma2    <- sqrt((a2*b2)/((a2+b2)^2*(a2+b2+1)))
mu.ges    <- mu2 - mu1
sigma.ges <- sigma1 + sigma2

norm.approx <- function(x){dnorm(x, mean = mu.ges, sd = sigma.ges)}
delta.approx[i] <- integrate(norm.approx, 0, 1)$value
}

#### H+ vs H-
posprob_Hplus <- numeric(length(data$y1))

for(i in 1:length(data$y1)){
conversion <- list(y1 = cumsum(data$y1)[1:i], y2 = cumsum(data$y2)[1:i], 
                   n1 = 1:length(data$y1[1:i]), n2 = 1:length(data$y2[1:i]))


plus.minus <- c(0, 1/2, 1/2, 0) # H+ vs H0
names(plus.minus) <- c("H1", "H+", "H-", "H0")

AB2 <- ab_test(conversion, prior_prob = plus.minus)  # uses default normal prior
posprob_Hplus[i] <- AB2$post_prob[2]

print(paste("Can I have", i, "scoops of ice cream?"))
}


### Plotting
plot(delta.approx, type = "l", ylim = c(0, 1), xlim = c(0, 100),
     bty = "n", xlab = "n", ylab = "P(delta > 0)", las = 1, col = "lightblue")
lines(delta.exact, lty = 2, col = "pink")
lines(posprob_Hplus, lty = 3, col = "orange")
legend(50, 1, lty = c(1, 2, 3), bty = "n",
       legend = c("normal approx","exact diff", "post prob H+"), 
       col = c("lightblue", "pink", "orange"))






