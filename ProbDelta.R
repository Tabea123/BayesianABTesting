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

png("sequential_diff_n80.png")
plot(delta.approx, type = "l", ylim = c(0, 1), xlim = c(0, 80),
     bty = "n", xlab = "n", ylab = "P(delta > 0)", las = 1, col = "lightblue")
lines(delta.exact, lty = 2, col = "pink")
lines(posprob_Hplus, lty = 3, col = "orange")
legend(50, 1, lty = c(1, 2, 3), bty = "n",
       legend = c("normal approx","exact diff", "post prob H+"), 
       col = c("lightblue", "pink", "orange"))
dev.off()


png("bothdifferencedistributions.png")
pdf.diff(a2, b2, a1, b1)
lines(density(pdf.delta.approx), lwd = 4, lty = 2, col = "pink")
dev.off()



#### H+ vs H-
posprob_Hplus <- numeric(length(data$y1))

CI.upper <- numeric(length(data$y1))
CI.lower <- numeric(length(data$y1))

for(i in 1:length(data$y1)){
conversion <- list(y1 = cumsum(data$y1)[1:i], y2 = cumsum(data$y2)[1:i], 
                   n1 = 1:length(data$y1[1:i]), n2 = 1:length(data$y2[1:i]))


plus.minus <- c(0, 1/2, 1/2, 0) # H+ vs H0
names(plus.minus) <- c("H1", "H+", "H-", "H0")

AB2 <- ab_test(conversion, prior_prob = plus.minus)  # uses default normal prior
posprob_Hplus[i] <- AB2$post_prob[2]

AB3 <-ab_test(conversion, prior_prob = plus.minus, 
              posterior = T, nsamples = 10000)
CI.upper[i] <- as.numeric(hdi(AB3$post$Hplus$arisk)[1])
CI.lower[i] <- as.numeric(hdi(AB3$post$Hplus$arisk)[2])

print(paste("Can I have", i, "scoops of ice cream?"))
}



CI.df <- data.frame(CI.lower, CI.upper)
y.upper <- CI.df[,1]
y.lower <- CI.df[,2]

# plot
png("tryout1.png", 
    width = 18, height = 18,
    units = "cm", res = 600,
    pointsize = 10)

par(cex.main = 1.5, mar = c(5, 5, 3, 3) + 0.1, mgp = c(3.5, 1, 0), 
    cex.lab = 1.5, font.lab = 2, cex.axis = 1.6, bty = "n", las = 1)

plot(1:length(conversion$n1), ylim = c(0, 1), 
     type = "n", xlab = "", ylab = "", yaxt = "n", xaxt = "n", bty = "n", 
     bg = "grey")

axis(1, at = seq(0, 100, by = 10), labels = seq(0, 100, 10), 
     lwd = 2, lwd.ticks = 2, line = -0.1)
axis(2, at = seq(0, 1, 0.1), labels = seq(0, 1, 0.1),
     lwd = 2, lwd.ticks = 2, line = -0.2)

mtext("n", side = 1, line = 3, las = 1, cex = 2, font = 2, adj = 0.5)
mtext("Width of CI", side = 2, line = 3.25, cex = 2, font = 2, las = 0)

polygon(c(1:length(conversion$n1),length(conversion$n1):1), c(y.upper, rev(y.lower)), 
        col = "lightgrey", border = NA)

lines(posprob_Hplus, lwd = 2, col = "lightblue")

legend(65, 0.8, "P(delta > 0)", col = "lightblue", lty = 1, bty = "n")
dev.off()


# plot
png("tryout2.png", 
    width = 18, height = 18,
    units = "cm", res = 600,
    pointsize = 10)

par(cex.main = 1.5, mar = c(5, 5, 3, 3) + 0.1, mgp = c(3.5, 1, 0), 
    cex.lab = 1.5, font.lab = 2, cex.axis = 1.6, bty = "n", las = 1)

plot(1:length(conversion$n1), ylim = c(0, 1), 
     type = "n", xlab = "", ylab = "", yaxt = "n", xaxt = "n", bty = "n", 
     bg = "grey")

axis(1, at = seq(0, 100, by = 10), labels = seq(0, 100, 10), 
     lwd = 2, lwd.ticks = 2, line = -0.1)
axis(2, at = seq(0, 1, 0.1), labels = seq(0, 1, 0.1),
     lwd = 2, lwd.ticks = 2, line = -0.2)

mtext("n", side = 1, line = 3, las = 1, cex = 2, font = 2, adj = 0.5)
mtext("Width of CI", side = 2, line = 3.25, cex = 2, font = 2, las = 0)

polygon(c(1:length(conversion$n1),length(conversion$n1):1), c(y.upper, rev(y.lower)), 
        col = "lightgrey", border = NA)

lines(posprob_Hplus, lwd = 2, col = "lightblue")

legend(65, 0.8, "P(delta > 0)", col = "lightblue", lty = 1, bty = "n")
dev.off()

