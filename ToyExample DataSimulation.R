################################################################################
##                                                                            ##
##                      Toy Example Data Simulation                           ## 
##                                                                            ##
################################################################################

# input:  -
# output: example_data1.csv, example_data2.csv


user <- 1:200

## simulate successes
A <- sample(c(rep(0, 65),rep(1, 35)))
B <- sample(c(rep(0, 50),rep(1, 50)))

## add version information
version <- rep(c("A", "B"), each = 100)

## combine in a data frame
data <- data.frame(user, successes = c(A, B), version)

## put users in a random order to simulate sequential data
random.order <- sample(user, length(user))
data1 <- data %>% slice(match(random.order, user))

y1 <- subset(data1, version == "A")$successes
y2 <- subset(data1, version == "B")$successes
data2 <- data.frame(y1, n1 = 1:100, y2, n2 = 1:100)
write.csv2(data2, file = "example_data1.csv")

## count conversions
conversion.A <- numeric(nrow(data1))
conversion.B <- numeric(nrow(data1))
n.A <- numeric(nrow(data1))
n.B <- numeric(nrow(data1))

for(i in 1:nrow(data1)){
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
}

conversion <- list(n1 = n.A, y1 = conversion.A, n2 = n.B, y2 = conversion.B)
write.csv2(conversion, file = "example_data2.csv")
