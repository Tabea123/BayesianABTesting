################################################################################
##                                                                            ##
##                      Probability P(B) > P(A)                               ##
##                                                                            ##
################################################################################

# this function is based on a formula published in
# Schmidt, M. N., & Morup, M. (2019). Efficient computation for Bayesian 
# comparison of two proportions. Statistics & Probability Letters, 145, 57-62.

prob.ab <- function(a1, b1, a2 ,b2){
  
  if (a1*b2 > a2*b1){ 
    
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
  return(exp(logz - lbeta(a1, b1) - lbeta(a2, b2)))
}
