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

# rm(list = ls())

# packages needed
library(dplyr)
library(bayesAB)
library(hypergeo)
library(tolerance)
library(abtest)
library(HDInterval)

# load simulated data
data <- read.csv2("SimluatedRekentuinData.csv")

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
AB1 <- bayesTest(A_success_ncrown, 
                 B_success_ncrown, 
                 priors = c('alpha' = 1, 'beta' = 1), 
                 n_samples = 1e5, 
                 distribution = 'bernoulli')

# plot prior distribution
plot(AB1, posteriors = FALSE, samples = FALSE) # A correpsonds to version B

# posterior distributions,
plot(AB1, priors = FALSE, samples = FALSE)

# switch A and B because of hypothesis B > A
AB11 <- bayesTest(B_success_ncrown, 
                 A_success_ncrown, 
                 priors = c('alpha' = 1, 'beta' = 1), 
                 n_samples = 1e5, 
                 distribution = 'bernoulli')

# check results
summary(AB1) # A correpsonds to version B

# and monte carlo 'integrated' samples (probability that version A is better
# than version B) 
plot(AB11, priors = FALSE, posteriors = FALSE)


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

# probability theta2 > theta1
prob.ab(a1, b1, a2, b2)


#-------------------------------------------------------------------------------
#                                                               
#                 #### 4. Analytically Compute P(B) - P(A) ####
#                                                               
#-------------------------------------------------------------------------------

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
HDI <- hdi(delta)

par(cex.main = 1.5, mar = c(5.5, 5.5, 5.9, 3) + 0.1, mgp = c(3.5, 1, 0), 
    cex.lab = 1.5, font.lab = 2, cex.axis = 1.8, bty = "n", las = 1)

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
      side = 1, line = 3, cex = 2.4, font = 2, adj = 0.5)
mtext("Density", side = 2, line = 3, cex = 2.4, font = 2, las = 0)

lines(density(delta), lwd = 4)

arrows(x0 = HDI[1], y0 = 4, x1 = HDI[2], y1 = 4, angle = 90, 
       length = 0.1, code = 3, lwd = 2.2)
text("95% HDI", x = mean(HDI), y = 4.5, cex = 1.8)


#-------------------------------------------------------------------------------
#                                                               
#   #### 5. A/B test with R package 'abtest' (Gronau & Wagenmakers, 2019) ####
#                                                               
#-------------------------------------------------------------------------------

## Load Preprocessed Data 

data4 <- read.csv2("Rekentuin_ABTestData.csv") 
conversion <- as.list(data4)

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
for(i in seq(5, length(conversion$n1), 1)){
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

axis(1, at = seq(0, 140, by = 10), labels = seq(0, 140, 10))
axis(2, at = seq(0, 3, 0.5), labels = seq(0, 3, 0.5))

mtext("n", side = 1, line = 3, las = 1, cex = 2, font = 0.2, adj = 0.5)
mtext("Width of CI", side = 2, line = 3, cex = 2, font = 2, las = 0)


polygon(c(5:length(conversion$n1),length(conversion$n1):5), c(y.upper, rev(y.lower)), 
        col = "lightgrey", border = NA)
