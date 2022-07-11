################################################################################
##                                                                            ##
##                   Analysis Webshop Data A/B Test                           ## 
##                                                                            ##
################################################################################

# input:  SimulatedWebshopData1.csv, SimulatedWebshopData2.csv
#
# output: webshop_posterior.png, webshop_uplift.png, webshop_delta.png,
#         webshop_prior1.png, webshop_prior2.png, webshop_prior3.png,
#         webshop_probwheel1.png, webshop_probwheel2.png,
#         webshop_probwheel3.png, indifferent_sequential.pdf, 
#         optimistic_sequential.pdf, conservative_sequential.pdf, 
#         indifferent_robustness.pdf, posterior_ind.png, posterior_con.png,
#         posterior_opt.png, webshop_sequentialdiff.png
#
# note: This script depends on functions that are defined in the scripts appellF1.R,
# pfdiff.R, and probab.R. Make sure to run these scripts before running this script.


rm(list = ls())

#-------------------------------------------------------------------------------
#                                                               
#           #### 0. Set Working Directory and Load/Preprocess Data ####                    
#                                                               
#-------------------------------------------------------------------------------

# packages needed
library(bayesAB) # version 1.1.3
library(abtest) # version 1.0.1
library(HDInterval)

data <- read.csv2("SimulatedWebshopData1.csv")
data3 <- read.csv2("SimulatedWebshopData2.csv") 
conversion <- as.list(data3[-1,])

#-------------------------------------------------------------------------------
#                                                               
#      ####  1. A/B test with R package 'bayesAB' (Portman, 2019) ####                
#                                                               
#-------------------------------------------------------------------------------

A_success <- ifelse(subset(data, version == "A")$successes == 1, 1, 0)
B_success <- ifelse(subset(data, version == "B")$successes == 1, 1, 0)

# success probabilities
sum(A_success)/length(A_success)
sum(B_success)/length(B_success)

# specify uniform prior for both success probabilities + fit bayesTest object
AB1 <- bayesTest(A_success, 
                 B_success, 
                 priors = c('alpha' = 1, 'beta' = 1), 
                 n_samples = 1e5, 
                 distribution = 'bernoulli')

# plot prior distribution
plot(AB1, posteriors = FALSE, samples = FALSE) 

# posterior distributions,
png("webshop_posterior.png", width = 18, height = 18, units = "cm", res = 600)
plot(AB1, priors = FALSE, samples = FALSE)
dev.off()

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
png("webshop_uplift.png", width = 18, height = 18, units = "cm", res = 600)
plot(AB11, priors = FALSE, posteriors = FALSE)
dev.off()


#-------------------------------------------------------------------------------
#                                                               
#                #### 2. Analytically Compute P(B > A) ####
#       
#-------------------------------------------------------------------------------

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


#-------------------------------------------------------------------------------
#                                                               
#                #### 3. Analytically Compute P(B) - P(A)  ####
#       
#-------------------------------------------------------------------------------

# exact calculation breaks down
# pdf.diff(a2, b2, a1, b1)

## normal approximation of difference distribution
a1 <- 1 + success.A
b1 <- 1 + failures.A
a2 <- 1 + success.B
b2 <- 1 + failures.B

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
png("webshop_delta.png", width = 18, height = 18, units = "cm", res = 600)

par(cex.main = 1.5,  mar = c(5, 5, 3, 3) + 0.1, mgp = c(3.5, 1, 0), 
    cex.lab = 1.5, font.lab = 2, cex.axis = 1.6, bty = "n", las = 1)

plot(density(delta), lwd = 4,
     main = "", xlab = "", ylab = " ",  
     xlim = c(-0.1, 0.1), ylim = c(0, 100),
     axes = FALSE, yaxt = "n", xaxt = "n")

axis(1, at = seq(-0.1, 0.1, 0.025), 
     labels = seq(-0.1, 0.1, 0.025), 
     lwd = 2, lwd.ticks = 2, line = -0.1)
axis(2, at = seq(0, 100, 10),  
     lwd = 2, lwd.ticks = 2, line = -0.2)

mtext(expression(paste("Difference", ~delta)), 
      side = 1, line = 3, las = 1, cex = 2, font = 2, adj = 0.5)
mtext("Density", side = 2, line = 3, cex = 2, font = 2, las = 0)

lines(density(delta), lwd = 4)

HDI <- hdi(delta)
arrows(x0 = HDI[1], y0 = 89, x1 = HDI[2], y1 = 89, angle = 90, 
       length = 0.1, code = 3, lwd = 2.2)

upper <- round(HDI[1],3)
lower <- round(HDI[2],3)
text(0.1, 95, paste0("95% HDI: [", upper, ";", lower, "]"), cex = 1.5, pos = 2)
text(0.1, 100, paste0("median = ", round(median(delta),3)),  cex = 1.5, pos = 2)

dev.off()



#-------------------------------------------------------------------------------
#                                                               
#   #### 4. A/B test with R package 'abtest' (Gronau & Wagenmakers, 2019) ####
#                                                               
#-------------------------------------------------------------------------------

### Analyze the results

## Specify Three Priors
# conservative person; OR: 1.05
conservative <- elicit_prior(q = c(0, 0.05, 0.1), prob = c(.025, .5, .975), 
                             what = "logor")

# optimistic person; OR: 1.2
optimistic <- elicit_prior(q = c(0.13, 0.18, 0.23), prob = c(.025, .5, .975), 
                           what = "logor")

# plot three priors for log odds ratio
png("webshop_prior1.png", width = 18, height = 18, units = "cm", res = 600)
plot_prior(what = "logor", hypothesis = "H+") # truncated normal
dev.off()
png("webshop_prior2.png", width = 18, height = 18, units = "cm", res = 600)
plot_prior(conservative, what = "logor", hypothesis = "H+")
dev.off()
png("webshop_prior3.png", width = 18, height = 18, units = "cm", res = 600)
plot_prior(optimistic, what = "logor", hypothesis = "H+")
dev.off()

## Hypothesis Testing

# success probabilities
round(tail(conversion$y1,1)/tail(conversion$n1, 1), 3)
round(tail(conversion$y2,1)/tail(conversion$n2, 1), 3)

# prior model probabilities: H+ = 0.5; H0 = 0.5
plus.null <- c(0, 1/2, 0, 1/2) # H+ vs H0
names(plus.null) <- c("H1", "H+", "H-", "H0")
AB.indifferent   <- ab_test(conversion, prior_prob = plus.null) 
AB.conservative  <- ab_test(conversion, conservative, prior_prob = plus.null) 
AB.optimistic    <- ab_test(conversion, optimistic, prior_prob = plus.null) 

# print results of ab_test
print(AB.indifferent) 
round(AB.indifferent$bf$bfplus0,3)
print(AB.conservative)
round(AB.conservative$bf$bfplus0,3)
print(AB.optimistic)
round(AB.optimistic$bf$bfplus0,3)

# visualize prior and posterior probabilities of the hypotheses 
# as probability wheels
png("webshop_probwheel1.png", width = 18, height = 12, units = "cm", res = 600)
prob_wheel(AB.indifferent)
dev.off()
png("webshop_probwheel2.png", width = 18, height = 12, units = "cm", res = 600)
prob_wheel(AB.conservative)
dev.off()
png("webshop_probwheel3.png", width = 18, height = 12, units = "cm", res = 600)
prob_wheel(AB.optimistic)
dev.off()

# plot posterior probabilities of the hypotheses sequentially
cairo_pdf("indifferent_sequential.pdf", width = 530/72, height = 400/72)
plot_sequential(AB.indifferent)
dev.off()

cairo_pdf("optimistic_sequential.pdf", width = 530/72, height = 400/72)
plot_sequential(AB.optimistic)
dev.off()

cairo_pdf("conservative_sequential.pdf", width = 530/72, height = 400/72)
plot_sequential(AB.conservative)
dev.off()

# plot BF robustness
cairo_pdf("indifferent_robustness.pdf", width = 530/72, height = 400/72)
plot_robustness(AB.indifferent)
dev.off()

## Parameter Estimation

# plot posterior distribution of log odds ratio
png("posterior_ind.png", width = 18, height = 18, units = "cm", res = 600)
plot_posterior(AB.indifferent, what = "logor")
dev.off()
png("posterior_con.png", width = 18, height = 18, units = "cm", res = 600)
plot_posterior(AB.conservative, what = "logor")
dev.off()
png("posterior_opt.png", width = 18, height = 18, units = "cm", res = 600)
plot_posterior(AB.optimistic, what = "logor")
dev.off()


## Compare Marginal Likelihoods
round(exp(AB.conservative$logml$logmlplus - AB.indifferent$logml$logmlplus) ,1)
round(exp(AB.optimistic$logml$logmlplus - AB.indifferent$logml$logmlplus) ,1)
round(exp(AB.optimistic$logml$logmlplus - AB.conservative$logml$logmlplus) ,1)


#-------------------------------------------------------------------------------
#                                                               
#                 #### 5. Sequential Analysis of Delta ####
#                                                               
#-------------------------------------------------------------------------------

CI.upper <- numeric(length(conversion$y1))
CI.lower <- numeric(length(conversion$y1))
mean.diff <- numeric(length(conversion$y1))

conversion <- lapply(conversion, function(x) c(0, x))

a1 <- conversion$y1
b1 <- conversion$n1 - conversion$y1
a2 <- conversion$y2
b2 <- conversion$n2 - conversion$y2

for(i in 1:length(conversion$y1)){
  
  diff.values <- rbeta(10000, 1 + a2[i], 1 + b2[i]) - rbeta(10000, 1 + a1[i], 1 + b1[i])
  
  CI.upper[i] <- as.numeric(hdi(diff.values)[1])
  CI.lower[i] <- as.numeric(hdi(diff.values)[2])
  
  mean.diff[i] <- mean(diff.values)
  
  print(paste("Can I have", i, "scoops of ice cream?"))
}

# plot
png("webshop_sequentialdiff.png", 
    width = 18, height = 18,
    units = "cm", res = 600,
    pointsize = 10)

par(cex.main = 1.5, mar = c(5, 5, 3, 3) + 0.1, mgp = c(3.5, 1, 0), 
    cex.lab = 1.5, font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)

plot(1:length(conversion$y1), ylim = c(-1, 1), 
     type = "n", xlab = "", ylab = "", yaxt = "n", xaxt = "n", bty = "n", 
     bg = "grey")

axis(1, at = seq(0, 20000, by = 1000), labels = seq(0, 20000, 1000), 
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



