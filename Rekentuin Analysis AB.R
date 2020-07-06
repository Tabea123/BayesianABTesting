################################################################################
##                                                                            ##
##               Preregistration Analysis Bayesian A/B Testing                ##
##                                                                            ##
################################################################################

# input:  SimluatedRekentuinData.csv, Rekentuin_ABTestData.csv
# output: -

#-------------------------------------------------------------------------------
#                                                               
# 0. Set Working Directory and Load/Preprocess Data                    
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

setwd("C:/Users/Tabea Hoffmann/Documents/Master/Thesis/Preregistration/")

# load simulated data
data <- read.csv2("SimluatedRekentuinData2.csv")

#-------------------------------------------------------------------------------
#                                                               
#             #### 1. Exclude Players + Select Last Game  ####
#                                                               
#-------------------------------------------------------------------------------

# exclude players that did not play any crown games
data2 <- data %>% group_by(user_id) %>% filter(any(game_type == "crown", na.rm = T))

# select last game + exclude players who's only crown game was their last game
data3 <- data2 %>% group_by(user_id)  %>% 
  
  # count frequency of games
  add_tally(game_type == "crown") %>% 
  
  # select last game for each child
  slice(which.max(as.numeric(time))) %>%   
  
  # exclude players 
  filter(!any(n == 1 && game_type == "crown"))  


# store information on last game for each version in separate vector
A_lastgame <- subset(data3, variant_id == "con")
B_lastgame <- subset(data3, variant_id == "exp")

# success has to correspond to not playing crown games
A_success_ncrown <- ifelse(A_lastgame$clicked_crown_game == 1, 0, 1) 
B_success_ncrown <- ifelse(B_lastgame$clicked_crown_game == 1, 0, 1) 
# 1 now corresponds to playing a non-crown game


#-------------------------------------------------------------------------------
#                                                               
#      ####  2. A/B test with R package 'bayesAB' (Portman, 2019) ####                
#                                                               
#-------------------------------------------------------------------------------

## Analyze the results 

# specify uniform prior for both success probabilities + fit bayesTest object
AB1 <- bayesTest(B_success_ncrown, 
                 A_success_ncrown, 
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
plot(AB1, priors = FALSE, posteriors = FALSE)


#-------------------------------------------------------------------------------
#                                                               
#                   #### 3. Analytically Compute P(B > A) ####
#                                                               
#-------------------------------------------------------------------------------

# get parameter of posterior distributions
success.A  <- sum(A_success_ncrown)
failures.A <- length((A_success_ncrown)) - success.A
success.B  <- sum(B_success_ncrown)
failures.B <- length((B_success_ncrown)) - success.B

a1 <- 1 + success.A
b1 <- 1 + failures.A
a2 <- 1 + success.B
b2 <- 1 + failures.B


if (a1*b2 < a2*b1){ # why is this true for our case 
  
  # normal series
  logz <- lgamma(a1 + a2) + lgamma(b1 + b2) - lgamma(a1 + b1 + a2 + b2 - 1) +
    log(genhypergeo(U = c(1, 1-a1, 1-b2), L = c(b1 + 1, a2 + 1), z = 1, 
                    tol = 1e-15,  series = T)) - log(b1*a2)
  
} else {
  
  # flipped series
  logz <- lgamma(a1 + a2) + lgamma(b1 + b2) - lgamma(a1 + b1 + a2 + b2 - 1) +
    log(genhypergeo(U = c(1, 1-b1, 1-a2), L = c(a1 + 1, b2 + 1), z = 1, 
                    tol = 1e-15,  series = T)) - log(a1*b2)
}

# posterior probability of the event theta1 > theta2 
exp(logz - lbeta(a1, b1) - lbeta(a2, b2))


#-------------------------------------------------------------------------------
#                                                               
#                 #### 4. Analytically Compute P(B) - P(A) ####
#                                                               
#-------------------------------------------------------------------------------

# joint posterior 
A <- beta(a1, b1)*beta(a2, b2)

# for 0 < p =< 1
smaller1 <- function(p){
  beta(a2, b1)*(p)^(b1+b2-1)*(1-p)^(a2+b1-1)*
    F1(b1, a1+b1+a2+b2-2, 1-a1, b1+a2, (1-p), (1-p^2)) / A
}

ps <- seq(0.01, 1, 0.01)

ymax <- beta(a1+a2-1, b1 + b2 -1) / A
plot(x = ps,y = mapply(smaller1, ps), xlim = c(-1, 1), type = "l", ylab = "pdf")

# for -1 <= p < 0
smaller0 <- function(p){
  beta(a1, b2)*(-p)^(b1+b2-1)*(1+p)^(a1+b2-1)*
    F1(b2, 1-a2, a1+a2+b1+b2-2, a1+b2, (1-p^2), (1+p)) / A
}

ps <- seq(-1, -0.01, 0.01)
lines(x = ps,y = mapply(smaller0, ps), xlim = c(-1, 0),type = "l")

# p = 0
points(x = 0, y = beta(a1+a2-1, b1 + b2 -1) / A, pch = 20)


#-------------------------------------------------------------------------------
#                                                               
# #### 5. A/B test with R package 'abtest' (Gronau & Wagenmakers, 2019) ####
#                                                               
#-------------------------------------------------------------------------------

#### Load Preprocessed Data ####

data4 <- read.csv2("Rekentuin_ABTestData2.csv") 
conversion <- as.list(data4)

#### Analyze the results ####

# success probabilities
tail(conversion$y1,1)/tail(conversion$n1, 1)
tail(conversion$y2,1)/tail(conversion$n2, 1)


## Hypothesis Testing

# prior model probabilities: H+ = 0.5; H0 = 0.5
plus.null <- c(0, 1/2, 0, 1/2) # H+ vs H0
names(plus.null) <- c("H1", "H+", "H-", "H0")
AB2 <- ab_test(conversion, prior_prob = plus.null)  # uses default normal prior

# print results of ab_test
print(AB2) 

# visualize prior and posterior probabilities of the hypotheses 
# as probability wheels
prob_wheel(AB2)

# plot posterior probabilities of the hypotheses sequentially
plot_sequential(AB2)


## Parameter Estimation

# plot posterior distribution of log odds ratio
plot_posterior(AB2, what = "logor")


# sequential analysis of CI 
CI.upper <- NULL
CI.lower <- NULL

# monitor CI in tranches of 10 data points
for(i in seq(10, 100, 10)){
  ab1 <- ab_test(lapply(conversion, '[', 1:i), prior_prob = plus.null, posterior = T)
  # store values in two separate vectors
  CI.lower <- c(CI.lower, hdi(ab1$post$Hplus$logor[1:10000])[1])
  CI.upper <- c(CI.upper, hdi(ab1$post$Hplus$logor[1:10000])[2])
}

CI.df <- data.frame(CI.lower, CI.upper)