################################################################################
##                                                                            ##
##                     Toy Example for Preregistration                        ##
##                                                                            ##
################################################################################

# input:  -
# output: -

#-------------------------------------------------------------------------------
#                                                               
# 0. Set Working Directory, Load Packages, Create Fictional Data
#                                                               
#-------------------------------------------------------------------------------

# rm(list = ls())

setwd("C:/Users/Tabea Hoffmann/Documents/Master/Thesis/Preregistration/")

# packages needed
library(dplyr)
library(bayesAB)
library(hypergeo)
library(tolerance)
library(abtest)

y1 <- sample(c(rep(0, 65),rep(1, 35)))
n1 <- length(y1)
y2 <- sample(c(rep(0, 50),rep(1, 50)))
n2 <- length(y1)

A.successprob <- sum(y1)/max(n1)
B.successprob <- sum(y2)/max(n2)

#-------------------------------------------------------------------------------
#                                                               
# 2. A/B test with R package 'bayesAB' (Portman, 2019)                 
#                                                               
#-------------------------------------------------------------------------------

# specify uniform prior for both success probabilities + fit bayesTest object
AB1 <- bayesTest(y2, 
                 y1, 
                 priors = c('alpha' = 1, 'beta' = 1), 
                 n_samples = 1e5, 
                 distribution = 'bernoulli')

# check results
summary(AB1) # A correpsonds to version B

## Visualize the results

# plot prior distribution
plot(AB1, posteriors = FALSE, samples = FALSE) # A correpsonds to version B

# posterior distributions,
# png("example_posterior.png")
plot(AB1, priors = FALSE, samples = FALSE)
# dev.off()

# and monte carlo 'integrated' samples (probability that version A is better
# than version B) 
# png("example_montecarlo.png")
plot(AB1, priors = FALSE, posteriors = FALSE)
# dev.off()


#-------------------------------------------------------------------------------
#                                                               
# 3. Analytically Compute P(B > A)
#                                                               
#-------------------------------------------------------------------------------

# get parameter of posterior distributions
success.A  <- sum(y1)
failures.A <- length(y1) - success.A
success.B  <- sum(y2)
failures.B <- length(y2) - success.B

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

# png("example_pdfdiff.png")
pdf.diff(a1, b1, a2, b2)
# dev.off()


#-------------------------------------------------------------------------------
#                                                               
# #### 5. A/B test with R package 'abtest' (Gronau & Wagenmakers, 2019) ####
#                                                               
#-------------------------------------------------------------------------------

conversion <- list(y1 = cumsum(y1), y2 = cumsum(y2), 
                   n1 = 1:length(y1), n2 = 1:length(y2))

## Hypothesis Testing

# prior model probabilities: H+ = 0.5; H0 = 0.5
plus.null <- c(0, 1/2, 0, 1/2) # H+ vs H0
names(plus.null) <- c("H1", "H+", "H-", "H0")
AB2 <- ab_test(conversion, prior_prob = plus.null)  # uses default normal prior

# print results of ab_test
print(AB2) 

# visualize prior and posterior probabilities of the hypotheses 
# as probability wheels
png("example_probwheel.png")
prob_wheel(AB2)
dev.off()

# plot posterior probabilities of the hypotheses sequentially
plot_sequential(AB2)


## Parameter Estimation

# plot posterior distribution of log odds ratio
plot_posterior(AB2, what = "logor")

