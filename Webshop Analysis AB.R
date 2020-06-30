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
setwd("C:/Users/Tabea Hoffmann/Documents/Master/Thesis/Preregistration")

# packages needed
library(dplyr)
library(abtest)

#-------------------------------------------------------------------------------
#                                                               
# 1. A/B test with R package 'abtest' (Gronau & Wagenmakers, 2019)
#                                                               
#-------------------------------------------------------------------------------

### Data Preprocessing

data <- read.csv2("SimulatedWebshopData.csv") 
conversion <- as.list(data)

### Analyze the results

## Specify Three Priors
# indifferent person
indifferent <- elicit_prior(q = c(-1.25, 0, 1.25), prob = c(.025, .5, .975), 
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


