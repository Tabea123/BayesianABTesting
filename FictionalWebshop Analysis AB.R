################################################################################
##                                                                            ##
##                   Analysis Webshop Data A/B Test                           ## 
##                                                                            ##
################################################################################

# input:  -
# output: SimulatedWebshopData.csv

# note: run appellF1.R, pfdiff.R, and probab.R before running this script

#-------------------------------------------------------------------------------
#                                                               
# 0. Set Working Directory                
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

#-------------------------------------------------------------------------------
#                                                               
#      ####  1. A/B test with R package 'bayesAB' (Portman, 2019) ####                
#                                                               
#-------------------------------------------------------------------------------

# read in data
data <- read.csv2("SimulatedWebshopData1.csv")

A_success <- ifelse(subset(data, version == "A")$successes == 1, 1, 0)
B_success <- ifelse(subset(data, version == "B")$successes == 1, 1, 0)

# specify uniform prior for both success probabilities + fit bayesTest object
AB1 <- bayesTest(A_success, 
                 B_success, 
                 priors = c('alpha' = 1, 'beta' = 1), 
                 n_samples = 1e5, 
                 distribution = 'bernoulli')

# plot prior distribution
plot(AB1, posteriors = FALSE, samples = FALSE) 

# posterior distributions,
plot(AB1, priors = FALSE, samples = FALSE)

# first switch A and B because of our hypothesis that B > A
AB11 <- bayesTest(B_success, 
                  A_success, 
                  priors = c('alpha' = 1, 'beta' = 1), 
                  n_samples = 1e5, 
                  distribution = 'bernoulli')

# check results
summary(AB11) # A correpsonds to version B

# monte carlo 'integrated' samples (probability that version A is 
# better than version B) 
plot(AB11, priors = FALSE, posteriors = FALSE)

## Analytically compute P(A > B)
# get parameter of posterior distributions
success.A  <- sum(A_success)
failures.A <- length(A_success) - success.A
success.B  <- sum(B_success)
failures.B <- length(B_success) - success.B

a1 <- 1 + success.A
b1 <- 1 + failures.A
a2 <- 1 + success.B
b2 <- 1 + failures.B

# probability theta2 > theta1
prob.ab(a1, b1, a2, b2)


## Analytically compute P(B) - P(A)
pdf.diff(a2, b2, a1, b1)


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


# plot 
delta <- rnorm(1000000, mu.ges, sigma.ges)

par(cex.main = 1.5,  mar = c(5, 5, 3, 3) + 0.1, mgp = c(3.5, 1, 0), 
    cex.lab = 1.5, font.lab = 2, cex.axis = 1.6, bty = "n", las = 1)

hist(delta, 
     freq = F, main = "", xlab = "", ylab = " ", 
     xlim = c(0, 0.05), ylim = c(0, 70),
     axes = FALSE, breaks = 17, yaxt = "n", xaxt = "n", col = "grey")

axis(1, at = c(0, 0.01, 0.02, 0.03, 0.04, 0.05), 
     labels = c(0, 0.01, 0.02, 0.03, 0.04, 0.05), 
     lwd = 2, lwd.ticks = 2, line = -0.1)
axis(2, at = seq(0, 70, 10),  
     lwd = 2, lwd.ticks = 2, line = -0.2)

mtext(expression(paste("Difference", ~delta)), 
      side = 1, line = 3, las = 1, cex = 2, font = 2, adj = 0.5)
mtext("Density", side = 2, line = 3, cex = 2, font = 2, las = 0)

lines(density(delta), lwd = 4)

HDI <- hdi(delta)
arrows(x0 = HDI[1], y0 = 65, x1 = HDI[2], y1 = 65, angle = 90, 
       length = 0.1, code = 3, lwd = 2.2)

upper <- round(HDI[1],3)
lower <- round(HDI[2],3)
text(0.05, 55, paste0("95% HDI: [", upper, ";", lower, "]"), cex = 1.5, pos = 2)
text(0.05, 60, paste0("median = ", round(median(delta),3)), cex = 1.5, pos = 2)

#-------------------------------------------------------------------------------
#                                                               
#   #### 2. A/B test with R package 'abtest' (Gronau & Wagenmakers, 2019) ####
#                                                               
#-------------------------------------------------------------------------------

### Data Preprocessing

data <- read.csv2("SimulatedWebshopData2.csv") 
conversion <- as.list(data[-1,])

### Analyze the results

## Specify Three Priors
# indifferent person
indifferent <- elicit_prior(q = c(-1.25, 0, 1.25), 
                            prob = c(.025, .5, .975), 
                            what = "logor")

# conservative person; OR: 1.05
conservative <- elicit_prior(q = c(0, 0.05, 0.1), prob = c(.025, .5, .975), 
                             what = "logor")

# optimistic person; OR: 1.2
optimistic <- elicit_prior(q = c(0.13, 0.18, 0.23), prob = c(.025, .5, .975), 
                           what = "logor")

# plot three priors for log odds ratio
# png("priors.png", width = 500, height = 600)
# par(mfrow = c(3, 1))
# plot_prior_mod(indifferent, what = "logor", hypothesis = "H+", main = "Indifferent") # truncated normal
# plot_prior_mod(conservative, what = "logor", hypothesis = "H+", main = "Conservative")
# plot_prior_mod(optimistic, what = "logor", hypothesis = "H+", main = "Optimistic")
# dev.off()

## Hypothesis Testing

# success probabilities
round(tail(conversion$y1,1)/tail(conversion$n1, 1), 3)
round(tail(conversion$y2,1)/tail(conversion$n2, 1), 3)

# prior model probabilities: H+ = 0.5; H0 = 0.5
plus.null <- c(0, 1/2, 0, 1/2) # H+ vs H0
names(plus.null) <- c("H1", "H+", "H-", "H0")
AB.indifferent   <- ab_test(conversion, indifferent, prior_prob = plus.null) 
AB.conservative  <- ab_test(conversion, conservative, prior_prob = plus.null) 
AB.optimistic    <- ab_test(conversion, optimistic, prior_prob = plus.null) 

AB   <- ab_test(conversion, prior_prob = plus.null) 

# print results of ab_test
print(AB.indifferent) 
print(AB.conservative)
print(AB.optimistic)

# visualize prior and posterior probabilities of the hypotheses 
# as probability wheels
png("probwheel.png", width = 600, height = 600)
par(mfrow = c(3, 2), cex.main = 1.5, mar = c(0.25, 1, 1.5, 1)) 
prob_wheel(AB.indifferent, type = "prior")
title("Prior Probabilities")
prob_wheel(AB.indifferent)
title("Posterior Probabilities")
prob_wheel(AB.conservative, type = "prior")
prob_wheel(AB.conservative)
prob_wheel(AB.optimistic, type = "prior")
prob_wheel(AB.optimistic)
dev.off()

# plot posterior probabilities of the hypotheses sequentially
# plot_sequential(AB.indifferent)
# plot_sequential(AB.conservative)
# plot_sequential(AB.optimistic)


## Parameter Estimation

# plot posterior distribution of log odds ratio
png("posterior_ind.png", width = 600, height = 600)
plot_posterior(AB.indifferent, what = "logor")
dev.off()
png("posterior_con.png", width = 600, height = 600)
plot_posterior(AB.conservative, what = "logor")
dev.off()
png("posterior_opt.png", width = 600, height = 600)
plot_posterior(AB.optimistic, what = "logor")
dev.off()

# sequential analysis of CI 
CI.upper <- NULL
CI.lower <- NULL

# monitor CI 
for(i in seq(7, length(conversion$n1), 1)){
  ab1 <- ab_test(lapply(conversion, '[', 1:i), prior_prob = plus.null, posterior = T)
  # store values in two separate vectors
  CI.lower <- c(CI.lower, hdi(ab1$post$Hplus$logor)[1])
  CI.upper <- c(CI.upper, hdi(ab1$post$Hplus$logor)[2])
}

CI.df <- data.frame(CI.lower, CI.upper)
y.upper <- CI.df[,1]
y.lower <- CI.df[,2]

plot(1:length(conversion$n1), ylim = c(0,3), 
     type = "n", xlab = "", ylab = "", yaxt = "n", bty = "n", axes = FALSE, 
     cex.lab = 1.3, cex.axis = 1.3, cex.main = 2, las = 1, pch = 21, lwd = 4,   
     bg = "grey")

axis(1, at = seq(0, 20000, by = 1000), labels = seq(0, 20000, by = 1000))
axis(2, at = seq(0, 3, 0.5), labels = seq(0, 3, 0.5))

mtext("n", side = 1, line = 3, las = 1, cex = 2, font = 0.2, adj = 0.5)
mtext("Width of CI", side = 2, line = 3, cex = 2, font = 2, las = 0)


polygon(c(7:length(conversion$n1),length(conversion$n1):7), c(y.upper, rev(y.lower)), 
        col = "lightgrey", border = NA)


## Compare Marginal Likelihoods
round(exp(AB.conservative$logml$logmlplus - AB.indifferent$logml$logmlplus) ,1)
round(exp(AB.optimistic$logml$logmlplus - AB.indifferent$logml$logmlplus) ,1)
round(exp(AB.optimistic$logml$logmlplus - AB.conservative$logml$logmlplus) ,1)


