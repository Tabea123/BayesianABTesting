################################################################################
##                                                                            ##
##                     Toy Example for Preregistration                        ##
##                                                                            ##
################################################################################

# input:  -
# output: -
# note: run appellF1.R, pfdiff.R, and probab.R before running this script

#-------------------------------------------------------------------------------
#                                                               
# 0. Set Working Directory, Load Packages, Create Fictional Data
#                                                               
#-------------------------------------------------------------------------------

rm(list = ls())

# packages needed
library(dplyr)
library(bayesAB)
library(hypergeo)
library(tolerance)
library(abtest)
library(HDInterval)
library(grDevices)

### Data Fabrication  
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

data <- read.csv2("example_data.csv")
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

# posterior distributions,
png("example_posterior.png",
    width = 18, height = 18,
    units = "cm", res = 600,
    pointsize = 10)
plot(AB1, priors = FALSE, samples = FALSE)
dev.off()
dev.off()

# and monte carlo 'integrated' samples (probability that version A is better
# than version B) 

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

#-------------------------------------------------------------------------------
#                                                               
# #### 5. A/B test with R package 'abtest' (Gronau & Wagenmakers, 2019) ####
#                                                               
#-------------------------------------------------------------------------------

conversion <- list(y1 = cumsum(data$y1), y2 = cumsum(data$y2), 
                   n1 = 1:length(data$y1), n2 = 1:length(data$y2))

## Hypothesis Testing

# prior model probabilities: H+ = 0.5; H0 = 0.5
plus.null <- c(0, 1/2, 0, 1/2) # H+ vs H0
names(plus.null) <- c("H1", "H+", "H-", "H0")
AB2 <- ab_test(conversion, prior_prob = plus.null)  # uses default normal prior

# print results of ab_test
print(AB2) 

# visualize prior and posterior probabilities of the hypotheses 
# as probability wheels
png("example_probwheel.png", 
    width = 15, height = 8,
    units = "cm", res = 600,
    pointsize = 10)
par(mai = c(0.1, 0, 0.1, 0),mar = c(0.1,0,0.1,0))
prob_wheel(AB2)
dev.off()

# plot posterior probabilities of the hypotheses sequentially
# 530 / 72 (width) by 400 / 72 (height); in pixels, 530 (width) by 400 (height).
png("example_sequential.png", width = 530/72, height = 400/72, units = "in", res = 600)
plot_sequential(AB2)
dev.off()

## Parameter Estimation

# plot posterior distribution of log odds ratio
png("example_priorposterior.png",
    width = 18, height = 18,
    units = "cm", res = 600,
    pointsize = 10)
plot_posterior(AB2, what = "logor")
dev.off()

plot_posterior(AB2, what = "p1p2")

# sequential analysis of CI 
CI.upper <- NULL
CI.lower <- NULL

# monitor CI for every data point
for(i in seq(1, length(conversion$n1), 1)){
  ab1 <- ab_test(lapply(conversion, '[', 1:i), prior_prob = plus.null, posterior = T)
  # store values in two separate vectors
  CI.lower <- c(CI.lower, hdi(ab1$post$Hplus$logor)[1])
  CI.upper <- c(CI.upper, hdi(ab1$post$Hplus$logor)[2])
}

CI.df <- data.frame(CI.lower, CI.upper)
y.upper <- CI.df[,1]
y.lower <- CI.df[,2]

# plot
png("example_sequentialCI.png", 
    width = 18, height = 18,
    units = "cm", res = 600,
    pointsize = 10)

par(cex.main = 1.5, mar = c(5, 5, 3, 3) + 0.1, mgp = c(3.5, 1, 0), 
    cex.lab = 1.5, font.lab = 2, cex.axis = 1.6, bty = "n", las = 1)

plot(1:length(conversion$n1), ylim = c(0,3), 
     type = "n", xlab = "", ylab = "", yaxt = "n", xaxt = "n", bty = "n", 
     bg = "grey")

axis(1, at = seq(0, 100, by = 10), labels = seq(0, 100, 10), 
     lwd = 2, lwd.ticks = 2, line = -0.1)
axis(2, at = seq(0, 3, 0.5), labels = seq(0, 3, 0.5),
     lwd = 2, lwd.ticks = 2, line = -0.2)

mtext("n", side = 1, line = 3, las = 1, cex = 2, font = 2, adj = 0.5)
mtext("Width of CI", side = 2, line = 3.25, cex = 2, font = 2, las = 0)

polygon(c(1:length(conversion$n1),length(conversion$n1):1), c(y.upper, rev(y.lower)), 
        col = "lightgrey", border = NA)
dev.off()
