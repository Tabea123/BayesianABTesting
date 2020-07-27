################################################################################
##                                                                            ##
##                     Toy Example for Preregistration                        ##
##                                                                            ##
################################################################################

# input:  -
# output: -
# note: run AppellF1Function.R, PDFDiffBetaFunction.R and ProbABFunction.R 
# before running this script

#-------------------------------------------------------------------------------
#                                                               
# 0. Set Working Directory, Load Packages, Create Fictional Data
#                                                               
#-------------------------------------------------------------------------------

# rm(list = ls())

# packages needed
library(bayesAB)
library(hypergeo)
library(tolerance)
library(abtest)
library(HDInterval)
library(grDevices)

# ### Data Fabrication
# y1 <- sample(c(rep(0, 65),rep(1, 35)))
# n1 <- length(y1)
# y2 <- sample(c(rep(0, 50),rep(1, 50)))
# n2 <- length(y1)
# 
# A.successprob <- sum(y1)/max(n1)
# B.successprob <- sum(y2)/max(n2)
# 
# data <- data.frame(y1, n1 = 1:n1, y2,n2 = 1:n2)
# write.csv2(data, file = "example_data.csv")

# read in data
data <- read.csv2("example_data.csv")

# greate list of data for later calculations
conversion <- list(y1 = cumsum(data$y1), y2 = cumsum(data$y2), 
                   n1 = 1:length(data$y1), n2 = 1:length(data$y2))

#-------------------------------------------------------------------------------
#                                                               
# 2. A/B test with R package 'bayesAB' (Portman, 2019)                 
#                                                               
#-------------------------------------------------------------------------------

# specify uniform prior for both success probabilities + fit bayesTest object
AB1 <- bayesTest(data$y1, 
                 data$y2, 
                 priors = c('alpha' = 1, 'beta' = 1), 
                 n_samples = 1e5, 
                 distribution = 'bernoulli')

## Visualize the results

# plot prior distribution
plot(AB1, posteriors = FALSE, samples = FALSE) 

# plot posterior distributions,
png("example_posterior.png",
    width = 18, height = 18,
    units = "cm", res = 600,
    pointsize = 10)
plot(AB1, priors = FALSE, samples = FALSE)
dev.off()
dev.off()


# plot probability that version A is better than version B

# first switch A and B 
AB11 <- bayesTest(data$y2, 
                  data$y1, 
                  priors = c('alpha' = 1, 'beta' = 1), 
                  n_samples = 1e5, 
                  distribution = 'bernoulli')

png("example_montecarlo.png",
    width = 18, height = 18,
    units = "cm", res = 600,
    pointsize = 10)
par(mar=c(5, 5, 5, 5) + 0.5)
plot(AB11, priors = FALSE, posteriors = FALSE)
dev.off()
dev.off()

# check results
summary(AB11) # A correpsonds to version B


#-------------------------------------------------------------------------------
#                                                               
# 3. Analytically Compute P(B > A)
#                                                               
#-------------------------------------------------------------------------------

# get parameter of posterior distributions
success.A  <- sum(data$y1)
failures.A <- length(data$y1) - success.A
success.B  <- sum(data$y2)
failures.B <- length(data$y2) - success.B

a1 <- 1 + success.A
b1 <- 1 + failures.A
a2 <- 1 + success.B
b2 <- 1 + failures.B

# probability theta2 > theta1
prob.ab(a1, b1, a2, b2)

#-------------------------------------------------------------------------------
#                                                               
# 4. Analytically Compute P(B) - P(A)
#                                                               
#-------------------------------------------------------------------------------

png("example_pdfdiff.png",
    width = 18, height = 18,
    units = "cm", res = 600,
    pointsize = 10)

pdf.diff(a2, b2, a1, b1)

dev.off()

## alternatively normal approximation for large a1, b1, a2, b2
# normal approximation of beta 1
mu1    <- a1/(a1+b1)
sigma1 <- sqrt((a1*b1)/((a1+b1)^2*(a1+b1+1)))
# normal approximation of beta 2
mu2    <- a2/(a2+b2)
sigma2 <- sqrt((a2*b2)/((a2+b2)^2*(a2+b2+1)))
# parameter values for difference
mu.ges    <- mu2 - mu1
sigma.ges <- sigma2 + sigma1

delta <- rnorm(1000000, mu.ges, sigma.ges)

# plot 
png("example_approximation.png",  width = 18, height = 18, 
    units = "cm", res = 600, pointsize = 10)

par(cex.main = 1.5,  mar = c(5, 5, 3, 3) + 0.1, mgp = c(3.5, 1, 0), 
    cex.lab = 1.5, font.lab = 2, cex.axis = 1.6, bty = "n", las = 1)

hist(delta, 
     freq = F, main = "", xlab = "", ylab = " ", 
     xlim = c(-1, 1), ylim = c(0, 6),
     axes = FALSE, breaks = 17, yaxt = "n", xaxt = "n", col = "grey")

axis(1, at = c(-1, -0.75, -0.50, -0.25, 0, 0.25, 0.50, 0.75, 1), 
     labels = c(-1, -0.75, -0.50, -0.25, 0, 0.25, 0.50, 0.75, 1),
     lwd = 2, lwd.ticks = 2, line = -0.1)
axis(2, at = seq(0, 6, 1),  
     lwd = 2, lwd.ticks = 2, line = -0.2)

mtext(expression(paste("Difference", ~delta)), 
      side = 1, line = 3, las = 1, cex = 2, font = 2, adj = 0.5)
mtext("Density", side = 2, line = 2.5, cex = 2, font = 2, las = 0)

lines(density(delta), lwd = 4)

HDI <- hdi(delta)
arrows(x0 = HDI[1], y0 = 4.5, x1 = HDI[2], y1 = 4.5, angle = 90, 
       length = 0.1, code = 3, lwd = 2.2)

upper <- round(HDI[1],3)
lower <- round(HDI[2],3)
text(1, 5.5, paste0("95% HDI: [", upper, ";", lower, "]"), cex = 1.5, pos = 2)
text(1, 6, paste0("median= ", round(median(delta),3)), cex = 1.5, pos = 2)

dev.off()

## alternatively compute posterior probability of H+
plus.minus <- c(0, 1/2, 1/2, 0) # H+ vs H0
names(plus.minus) <- c("H1", "H+", "H-", "H0")
AB2 <- ab_test(conversion, prior_prob = plus.minus, 
               posterior = T, nsamples = 1000000)

png("example_hpluspost.png", width = 18, height = 18, units = "cm", res = 600)

par(cex.main = 1.5,  mar = c(5, 5, 3, 3) + 0.1, mgp = c(3.5, 1, 0), 
    cex.lab = 1.5, font.lab = 2, cex.axis = 1.6, bty = "n", las = 1)

plot(density(AB2$post$Hplus$arisk), type = "n", main = "", bty = "n",
     yaxt = "n", xaxt = "n", axes = F, xlim = c(-1, 1), ylim = c(0, 6))

lines(density(AB2$post$Hplus$arisk), lwd = 2)

axis(1, at = seq(-1, 1, 0.25), 
     labels = seq(-1, 1, 0.25), 
     lwd = 2, lwd.ticks = 2, line = -0.1)
axis(2, at = seq(0, 6, 1),  
     lwd = 2, lwd.ticks = 2, line = -0.2)

dev.off()


#-------------------------------------------------------------------------------
#                                                               
# #### 5. A/B test with R package 'abtest' (Gronau & Wagenmakers, 2019) ####
#                                                               
#-------------------------------------------------------------------------------

## Hypothesis Testing

# prior model probabilities: H+ = 0.5; H0 = 0.5
plus.null <- c(0, 1/2, 0, 1/2) # H+ vs H0
names(plus.null) <- c("H1", "H+", "H-", "H0")
AB3 <- ab_test(conversion, prior_prob = plus.null, posterior = T)  # uses default normal prior

# print results of ab_test
print(AB3) 

# visualize prior and posterior probabilities of the hypotheses 
# as probability wheels
png("example_probwheel.png", 
    width = 15, height = 8,
    units = "cm", res = 600,
    pointsize = 10)
par(mai = c(0.1, 0, 0.1, 0),mar = c(0.1,0,0.1,0))
prob_wheel(AB3)
dev.off()

# plot posterior probabilities of the hypotheses sequentially
# 530 / 72 (width) by 400 / 72 (height); in pixels, 530 (width) by 400 (height).
cairo_pdf("example_sequential.pdf", width = 530/72, height = 400/72)
plot_sequential(AB3)
dev.off()

## Parameter Estimation

# plot posterior distribution of log odds ratio
png("example_priorposterior.png",
    width = 18, height = 18,
    units = "cm", res = 600,
    pointsize = 10)
plot_posterior(AB3, what = "logor")
dev.off()

plot_posterior(AB3, what = "p1p2")


#-------------------------------------------------------------------------------
#                                                               
# #### 6. Comparison P(diff > 0)   
#                                                               
#-------------------------------------------------------------------------------

### Increase N to 200
y1 <- sample(c(rep(0, 130),rep(1, 70)))
n1 <- length(y1)
y2 <- sample(c(rep(0, 100),rep(1, 100)))
n2 <- length(y1)

data2 <- data.frame(y1, n1 = 1:n1, y2,n2 = 1:n2)

### Exact Difference
a1 <- 1
b1 <- 1
a2 <- 1
b2 <- 1

delta.exact  <- numeric(length(data2$y1))
for(i in 1:length(data2$y1)){
  
  a2 <- a2 + data2$y1[i]
  b2 <- b2 + ifelse(data2$y1 == 1, 0, 1)[i]
  a1 <- a1 + data2$y2[i]
  b1 <- b1 + ifelse(data2$y2 == 1, 0, 1)[i]
  
  delta.exact[i] <- integrate(pdf.diff2, 0, 1, stop.on.error = F)$value
}


### Normal Approximation
a1 <- 1
b1 <- 1
a2 <- 1
b2 <- 1

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


#### Posterior Probability of H+
posprob_Hplus <- numeric(length(data2$y1))

for(i in 1:length(data2$y1)){
  conversion <- list(y1 = cumsum(data2$y1)[1:i], y2 = cumsum(data2$y2)[1:i], 
                     n1 = 1:length(data2$y1[1:i]), n2 = 1:length(data2$y2[1:i]))
  
  
  plus.minus <- c(0, 1/2, 1/2, 0) # H+ vs H0
  names(plus.minus) <- c("H1", "H+", "H-", "H0")
  
  AB4 <- ab_test(conversion, prior_prob = plus.minus)  # uses default normal prior
  posprob_Hplus[i] <- AB4$post_prob[2]
}


### Plot the three options
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


#-------------------------------------------------------------------------------
#                                                               
# #### 7. Sequential Analysis of Delta
#                                                               
#-------------------------------------------------------------------------------

seq_posprob <- numeric(length(data$y1))

CI.upper <- numeric(length(data$y1))
CI.lower <- numeric(length(data$y1))

mean.ar <- numeric(length(data$y1))

for(i in 1:length(data$y1)){
  conversion <- list(y1 = cumsum(data$y1)[1:i], y2 = cumsum(data$y2)[1:i], 
                     n1 = 1:length(data$y1[1:i]), n2 = 1:length(data$y2[1:i]))
  
  
  plus.minus <- c(0, 1/2, 1/2, 0) # H+ vs H0
  names(plus.minus) <- c("H1", "H+", "H-", "H0")
  
  AB4 <- ab_test(conversion, prior_prob = plus.minus)  # uses default normal prior
  seq_posprob[i] <- AB4$post_prob[2]
  
  AB5 <- ab_test(conversion, prior_prob = plus.minus, 
                posterior = T, nsamples = 10000)
  CI.upper[i] <- as.numeric(hdi(AB5$post$H1$arisk)[1])
  CI.lower[i] <- as.numeric(hdi(AB5$post$H1$arisk)[2])
  
  mean.ar[i] <- mean(AB5$post$H1$arisk)
  
  print(paste("Can I have", i, "scoops of ice cream?"))
}


# plot
png("example_sequentialdiff.png", 
    width = 18, height = 18,
    units = "cm", res = 600,
    pointsize = 10)

par(cex.main = 1.5, mar = c(5, 5, 3, 3) + 0.1, mgp = c(3.5, 1, 0), 
    cex.lab = 1.5, font.lab = 2, cex.axis = 1.6, bty = "n", las = 1)

plot(1:length(conversion$n1), ylim = c(-1, 1), 
     type = "n", xlab = "", ylab = "", yaxt = "n", xaxt = "n", bty = "n", 
     bg = "grey")

axis(1, at = seq(0, 100, by = 10), labels = seq(0, 100, 10), 
     lwd = 2, lwd.ticks = 2, line = -0.1)
axis(2, at = seq(-1, 1, 0.2), labels = seq(-1, 1, 0.2),
     lwd = 2, lwd.ticks = 2, line = -0.2)

mtext("n", side = 1, line = 3, las = 1, cex = 2, font = 2, adj = 0.5)
mtext(expression(paste("Difference", ~delta)), side = 2, line = 3.25, cex = 2, font = 2, las = 0)

polygon(c(1:length(conversion$n1),length(conversion$n1):1), c(CI.upper, rev(CI.lower)), 
        col = "lightgrey", border = NA)

lines(mean.ar, lwd = 2, col = "orange")
abline(h = 0, lwd = 2, col = "darkgrey", lty = 2)

dev.off()
