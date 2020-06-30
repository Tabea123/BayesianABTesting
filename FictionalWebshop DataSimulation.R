################################################################################
##                                                                            ##
##                         A/B Test Data Simulation                           ## 
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
library(dplyr)

#-------------------------------------------------------------------------------
#                                                               
# 1. Simulate Data 
#                                                               
#-------------------------------------------------------------------------------

## create users
user <- 1:20000

## simulate successes
# probability of a success in version A (baseline)
successA <- rbinom(max(user)/2, 1, seq(.10, .11, by = 0.001))
# probability of a success in version B ('test it out')
successB <- rbinom(max(user)/2, 1, seq(.12, .13, by = 0.001))

## add version information
version <- rep(c("A", "B"), each = max(user)/2)

## combine in a data frame
data <- data.frame(user.id = user, successes = c(successA, successB), version)

## put users in a random order to simulate sequential data
random.order <- sample(user, length(user))
data1 <- data %>% slice(match(random.order, user.id))

## count conversions
conversion.A <- numeric(nrow(data1))
conversion.B <- numeric(nrow(data1))
n.A <- numeric(nrow(data1))
n.B <- numeric(nrow(data1))

for(i in 1:nrow(data1))
  if(data1[i,]$successes == 1 && 
     data1[i,]$version == "A"){
    
    conversion.A[i+1] <- conversion.A[i] + 1
    conversion.B[i+1] <- conversion.B[i] 
    n.A[i+1] <- n.A[i] + 1 
    n.B[i+1] <- n.B[i] 
    
  } else if(data1[i,]$successes == 0 && 
            data1[i,]$version == "A"){
    
    conversion.A[i+1] <- conversion.A[i] 
    conversion.B[i+1] <- conversion.B[i] 
    n.A[i+1] <- n.A[i] + 1
    n.B[i+1] <- n.B[i] 
    
  } else if(data1[i,]$successes == 1 && 
            data1[i,]$version == "B"){
    
    conversion.A[i+1] <- conversion.A[i] 
    conversion.B[i+1] <- conversion.B[i] + 1
    n.A[i+1] <- n.A[i] 
    n.B[i+1] <- n.B[i] + 1
    
  } else if(data1[i,]$successes == 0 && 
            data1[i,]$version == "B"){
    
    conversion.A[i+1] <- conversion.A[i] 
    conversion.B[i+1] <- conversion.B[i] 
    n.A[i+1] <- n.A[i] 
    n.B[i+1] <- n.B[i] + 1
    
  }


#-------------------------------------------------------------------------------
#                                                               
# 3. Export Data
#                                                               
#-------------------------------------------------------------------------------

## combine everything in one data set
conversion <- list(n1 = n.A, y1 = conversion.A, n2 = n.B, y2 = conversion.B)

## export data
write.csv2(as.data.frame(conversion), "SimulatedWebshopData.csv", 
           row.names = FALSE)


