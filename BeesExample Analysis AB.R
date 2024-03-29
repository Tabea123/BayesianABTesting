################################################################################
##                                                                            ##
##                         Bees Example Analyses                             ##
##                                                                            ##
################################################################################

# input: bees_data1.csv, bees_data2.csv
# output: example_prior.png, example_posterior.png, example_pdfdiff.png,
#        example_approximation.png, example_hpluspost.png, example_probwheel.png,
#        example_sequential.pdf, example_robustness.png, example_priorposterior.png,
#         example_sequentialdiff.png
#
# note: This script depends on functions that are defined in the scripts appellF1.R,
# pfdiff.R, and probab.R. Make sure to run these scripts before running this script.


#-------------------------------------------------------------------------------
#                                                               
#  #### 1. Set Working Directory, Load Packages, Create Fictional Data ####
#                                                               
#-------------------------------------------------------------------------------

# rm(list = ls())

# packages needed
library(bayesAB) # version 1.1.3
library(abtest) # version 1.0.1
library(hypergeo)
library(tolerance)
library(HDInterval)
library(grDevices)

# read in data
bees1 <- read.csv2("bees_data1.csv")
bees2 <- as.list(read.csv2("bees_data2.csv")[-1,-1])


#-------------------------------------------------------------------------------
#                                                               
#       #### 2. A/B test with R package 'bayesAB' (Portman, 2019) ####
#                                                               
#-------------------------------------------------------------------------------

# specify uniform prior for both success probabilities + fit bayesTest object
AB1 <- bayesTest(bees1$y1, 
                 bees1$y2, 
                 priors = c('alpha' = 1, 'beta' = 1), 
                 n_samples = 1e5, 
                 distribution = 'bernoulli')

## Visualize the results

# plot prior distribution
png("example_prior.png",
    width = 18, height = 18,
    units = "cm", res = 600,
    pointsize = 10)
plot(AB1, posteriors = FALSE, samples = FALSE) 
dev.off()


# plot posterior distributions,
png("example_posterior.png",
    width = 18, height = 18,
    units = "cm", res = 600,
    pointsize = 10)
plot(AB1, priors = FALSE, samples = FALSE)
dev.off()


# plot probability that version A is better than version B

# first switch A and B 
AB11 <- bayesTest(bees1$y2, 
                  bees1$y1, 
                  priors = c('alpha' = 1, 'beta' = 1), 
                  n_samples = 1e5, 
                  distribution = 'bernoulli')

png("example_uplift.png",
    width = 18, height = 18,
    units = "cm", res = 600,
    pointsize = 10)
par(mar=c(5, 5, 5, 5) + 0.5)
plot(AB11, priors = FALSE, posteriors = FALSE)
dev.off()


# check results
summary(AB11) # A correpsonds to version B


#-------------------------------------------------------------------------------
#                                                               
#               #### 3. Analytically Compute P(B > A) ####
#                                                               
#-------------------------------------------------------------------------------

# get parameter of posterior distributions
success.A  <- sum(bees1$y1)
failures.A <- length(bees1$y1) - success.A
success.B  <- sum(bees1$y2)
failures.B <- length(bees1$y2) - success.B

a1 <- 1 + success.A
b1 <- 1 + failures.A
a2 <- 1 + success.B
b2 <- 1 + failures.B

# probability theta2 > theta1
prob.ab(a1, b1, a2, b2)

#-------------------------------------------------------------------------------
#                                                               
#              #### 4. Analytically Compute P(B) - P(A) ####
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
mu1  <- a1/(a1+b1)
var1 <- ((a1*b1)/((a1+b1)^2*(a1+b1+1)))
# normal approximation of beta 2
mu2  <- a2/(a2+b2)
var2 <- ((a2*b2)/((a2+b2)^2*(a2+b2+1)))
# parameter values for difference
mu.ges  <- mu2 - mu1
var.ges <- (var1 + var2)

delta <- rnorm(1000000, mu.ges, sqrt(var.ges))

# plot 
png("normal_approximation.png",  width = 18, height = 18, 
    units = "cm", res = 600, pointsize = 10)

par(cex.main = 1.5,  mar = c(5, 5, 3, 3) + 0.1, mgp = c(3.5, 1, 0), 
    cex.lab = 1.5, font.lab = 2, cex.axis = 1.6, bty = "n", las = 1, xpd = T)

plot(density(delta), lwd = 4, xlab = "", ylab = " ", 
     xlim = c(-1, 1), ylim = c(0, 7), axes = FALSE, yaxt = "n", xaxt = "n")

axis(1, at = c(-1, -0.75, -0.50, -0.25, 0, 0.25, 0.50, 0.75, 1), 
     labels = c(-1, -0.75, -0.50, -0.25, 0, 0.25, 0.50, 0.75, 1),
     lwd = 2, lwd.ticks = 2, line = -0.1)
axis(2, at = seq(0, 7, 1),  
     lwd = 2, lwd.ticks = 2, line = -0.2)

mtext(expression(paste("Difference", ~delta)), 
      side = 1, line = 3, las = 1, cex = 2, font = 2, adj = 0.5)
mtext("Density", side = 2, line = 2.5, cex = 2, font = 2, las = 0)

HDI <- hdi(delta)
arrows(x0 = HDI[1], y0 = 6.2, x1 = HDI[2], y1 = 6.2, angle = 90, 
       length = 0.1, code = 3, lwd = 2.2)

upper <- round(HDI[1],3)
lower <- round(HDI[2],3)
text(1, 7, paste0("95% HDI: [", upper, ";", lower, "]"), cex = 1.5, pos = 2)
text(1, 8, paste0("median= ", round(median(delta),3)), cex = 1.5, pos = 2)

dev.off()



#-------------------------------------------------------------------------------
#                                                               
#                 #### 5. Sequential Analysis of Delta ####
#                                                               
#-------------------------------------------------------------------------------

CI.upper <- numeric(length(bees2$y1))
CI.lower <- numeric(length(bees2$y1))
mean.diff <- numeric(length(bees2$y1))

conversion <- lapply(bees2, function(x) c(0, x))

a1 <- conversion$y1
b1 <- conversion$n1 - conversion$y1
a2 <- conversion$y2
b2 <- conversion$n2 - conversion$y2

for(i in 1:length(conversion$y1)){
  
  diff.values <- rbeta(1000000, 1 + a2[i], 1 + b2[i]) - rbeta(1000000, 1 + a1[i], 1 + b1[i])
   
  CI.upper[i] <- as.numeric(hdi(diff.values)[1])
  CI.lower[i] <- as.numeric(hdi(diff.values)[2])
   
  mean.diff[i] <- mean(diff.values)
  
  print(paste("I am counting", i, "snowflakes"))
}

png("example_sequentialdiff.png",
    width = 18, height = 18,
    units = "cm", res = 600,
    pointsize = 10)

par(cex.main = 1.5, mar = c(5, 5, 3, 3) + 0.1, mgp = c(3.5, 1, 0), 
    cex.lab = 1.5, font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)

plot(1:length(conversion$y1), ylim = c(-1, 1), 
     type = "n", xlab = "", ylab = "", yaxt = "n", xaxt = "n", bty = "n", 
     bg = "grey")

axis(1, at = seq(0, 200, by = 10), labels = seq(0, 200, 10), 
     lwd = 2, lwd.ticks = 2, line = -0.1)
axis(2, at = seq(-1, 1, 0.2), labels = seq(-1, 1, 0.2),
     lwd = 2, lwd.ticks = 2, line = -0.2)

mtext("n", side = 1, line = 3, las = 1, cex = 2, font = 2, adj = 0.5)
mtext(expression(paste("Difference", ~delta)), side = 2, line = 3.25, cex = 2, font = 2, las = 0)

polygon(c(1:length(CI.upper),length(CI.upper):1), c((CI.upper), rev(CI.lower)), 
        col = "lightgrey", border = NA)

abline(h = 0,  lwd = 2, col = "darkgrey", lty = 2)
lines(mean.diff, lwd = 2, col = "orange")

dev.off()


#-------------------------------------------------------------------------------
#                                                               
# #### 6. A/B test with R package 'abtest' (Gronau & Wagenmakers, 2019) ####
#                                                               
#-------------------------------------------------------------------------------

## Hypothesis Testing

# prior model probabilities: H+ = 0.5; H0 = 0.5
plus.null <- c(0, 1/2, 0, 1/2) # H+ vs H0
names(plus.null) <- c("H1", "H+", "H-", "H0")
AB2 <- ab_test(bees2, prior_prob = plus.null, posterior = T)  # uses default normal prior

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
cairo_pdf("example_sequential.pdf", width = 530/72, height = 400/72)
plot_sequential(AB2)
dev.off()

# plot BF robustness plot
png("example_robustness.png", 
    width = 18, height = 18,
    units = "cm", res = 600,
    pointsize = 10)
plot_robustness(AB2, mu_range = c(0, 2), sigma_range = c(0.1, 1), bftype = "BF+0")
dev.off()

## Parameter Estimation

# plot posterior distribution of log odds ratio
png("example_priorposterior.png",
    width = 18, height = 18,
    units = "cm", res = 600,
    pointsize = 10)
plot_posterior(AB2, what = 'logor')
dev.off()

# independent posterior distributions
plot_posterior(AB2, what = "p1p2")


