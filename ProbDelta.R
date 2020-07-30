################################################################################
##                                                                            ##
##                      Comparison P(diff > 0)                                ##
##                                                                            ##
################################################################################

# note: run AppellF1Function.R, and PDFDiffBetaFunction.R before running this script

conversion <- as.list(read.csv2("example_data2.csv")[,-1])

x <- c(conversion$y1[-1],35)
A <- x - conversion$y1
hitsA <- c(0, A[-201])

y <- c(conversion$y2[-1],50)
B <- y - conversion$y2
hitsB <- c(0, B[-201])

data2 <- data.frame(y1 = hitsA, n1 = conversion$n1, y2 = hitsB ,n2 = conversion$n2)
data2 <- data2[-1,]

a1 <- 1
b1 <- 1
a2 <- 1
b2 <- 1

### exact diff
delta.exact  <- numeric(length(data2$y1))
for(i in 1:length(data2$y1)){
  
  a2 <- a2 + data2$y1[i]
  b2 <- b2 + ifelse(data2$y1 == 1, 0, 1)[i]
  a1 <- a1 + data2$y2[i]
  b1 <- b1 + ifelse(data2$y2 == 1, 0, 1)[i]
  
  delta.exact[i] <- integrate(pdf.diff2, 0, 1, stop.on.error = F)$value
}

a1 <- 1
b1 <- 1
a2 <- 1
b2 <- 1

### Normal Approximation
delta.approx <- numeric(length(data2$y1))
for(i in 1:length(data2$y1)){
  
a1 <- a1 + data2$y1[i]
b1 <- b1 + ifelse(data2$y1 == 1, 0, 1)[i]
a2 <- a2 + data2$y2[i]
b2 <- b2 + ifelse(data2$y2 == 1, 0, 1)[i]

mu1       <- a1/(a1+b1)
sigma1    <- sqrt((a1*b1)/((a1+b1)^2*(a1+b1+1)))
mu2       <- a2/(a2+b2)
sigma2    <- sqrt((a2*b2)/((a2+b2)^2*(a2+b2+1)))
mu.ges    <- mu2 - mu1
sigma.ges <- sigma1 + sigma2

norm.approx <- function(x){dnorm(x, mean = mu.ges, sd = sigma.ges)}
delta.approx[i] <- integrate(norm.approx, 0, 1)$value
}

#### H+ vs H

posprob_Hplus <- numeric(length(data2$y1))

for(i in 1:length(data2$y1)){

  conversion1 <- list(y1 = cumsum(data2$y1)[1:i], y2 = cumsum(data2$y2)[1:i], 
                     n1 = 1:length(data2$y1[1:i]), n2 = 1:length(data2$y2[1:i]))
  
  
  
plus.minus <- c(0, 1/2, 1/2, 0) # H+ vs H0
names(plus.minus) <- c("H1", "H+", "H-", "H0")

AB2 <- ab_test(conversion1, prior_prob = plus.minus)  # uses default normal prior
posprob_Hplus[i] <- AB2$post_prob[2]

print(paste("Can I have", i, "scoops of ice cream?"))
}


### Plotting
png("example_pdeltazero.png", width = 18, height = 18, units = "cm", res = 600)
par(cex.main = 1.5, mar = c(5, 5, 3, 3) + 0.1, mgp = c(3.5, 1, 0), 
    cex.lab = 2, font.lab = 2, cex.axis = 1.2, bty = "n", las = 1)

plot(delta.approx, type = "l", ylim = c(0, 1), xlim = c(0, 200),
     bty = "n", xlab = "", ylab = "", axes = F, las = 1, col = "lightblue", 
     lwd = 2, cex.lab = 1.5, cex.axis = 1.8, cex.main = 1.5)

axis(1, at = seq(0, 200, 20), labels = seq(0, 200, 20), 
     lwd = 2, lwd.ticks = 2, line = -0.1)
axis(2, at = seq(0, 1, 0.1), lwd = 2, lwd.ticks = 2, line = -0.2)

mtext(expression(paste("P(", ~delta, " > 0)")), side = 2, 
      line = 3, las = 0, cex = 2, font = 0.2, adj = 0.5)
mtext("n", side = 1, line = 2.5, cex = 2, font = 2, las = 1)

lines(delta.exact, lty = 2, lwd = 2, col = "pink")
lines(posprob_Hplus, lty = 3, lwd = 2, col = "orange")


legend(110, 0.2, lty = c(1, 2, 3), bty = "n",
       legend = c("Normal Approximation","Exact Difference", 
                  "Posterior Probability H+"), 
       col = c("lightblue", "pink", "orange"))


dev.off()



