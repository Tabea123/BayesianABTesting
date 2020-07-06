################################################################################
##                                                                            ##
##                            Appell series F1                                ##
##                                                                            ##
################################################################################

# this function is based on the F1 function of the 'tolerance' package
# written by Derek S. Young (derek.young@uky.edu)
# see https://github.com/cran/tolerance/blob/master/R/F1.R


appell.F1 <- function(a, b, b.prime, c, x, y, ...)
{
  A1.simple <- function(u, a, b, b.prime, c, x, y)
  {
    one <- (a-1)*log(u) 
    two <- (c-a-1)*log(1-u)
    three <- (-b)*log(1-u*x)
    four <- (-b.prime)*log((1-u*y))
    return(exp(one + two + three + four))
  }
  gamma(c)/(gamma(a)*gamma(c-a))*as.numeric(integrate(A1.simple, 0, 1, a = a, b = b, b.prime = b.prime, c = c, x = x, y = y, ...)$value)
}
