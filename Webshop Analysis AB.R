################################################################################
##                                                                            ##
##                   Analysis Webshop Data A/B Test                           ## 
##                                                                            ##
################################################################################

# input:  -
# output: SimulatedWebshopData.csv


#-------------------------------------------------------------------------------
#                                                               
# 0. Set Working Directory                
#                                                               
#-------------------------------------------------------------------------------

rm(list = ls())

# packages needed
library(dplyr)
library(bayesAB)
library(abtest)

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

# check results
summary(AB1) # A correpsonds to version B

## Visualize the results

# plot prior distribution
plot(AB1, posteriors = FALSE, samples = FALSE) # A correpsonds to version B

# posterior distributions,
plot(AB1, priors = FALSE, samples = FALSE)

# and monte carlo 'integrated' samples (probability that version A is better
# than version B) 

# first switch A and B 
AB11 <- bayesTest(B_success, 
                  A_success, 
                  priors = c('alpha' = 1, 'beta' = 1), 
                  n_samples = 1e5, 
                  distribution = 'bernoulli')
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
pdf.diff(a1, b1, a2, b2)


## for large a1, b1, a2, b2: normal approximation of difference between beta distributions
# normal approximation of beta 1
mu1    <- a1/(a1+b1)
sigma1 <- sqrt(a1*b1)/((a1+b1)^2*(a1+b1+1))
# normal approximation of beta 2
mu2    <- a2/(a2+b2)
sigma2 <- sqrt(a2*b2)/((a2+b2)^2*(a2+b2+1))
# parameter values for difference
mu.ges    <- mu2 - mu1
sigma.ges <- sigma2 - sigma1
# plot density
plot(density(rnorm(100, mu.ges, sigma.ges)))


#-------------------------------------------------------------------------------
#                                                               
#   #### 2. A/B test with R package 'abtest' (Gronau & Wagenmakers, 2019) ####
#                                                               
#-------------------------------------------------------------------------------

### Data Preprocessing

data <- read.csv2("SimulatedWebshopData.csv") 
conversion <- as.list(data)

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

## Compare Marginal Likelihoods
round(exp(AB.conservative$logml$logmlplus - AB.indifferent$logml$logmlplus) ,1)
round(exp(AB.optimistic$logml$logmlplus - AB.indifferent$logml$logmlplus) ,1)
round(exp(AB.optimistic$logml$logmlplus - AB.conservative$logml$logmlplus) ,1)


