################################################################################
##                                                                            ##
##                 Rekentuin Data Pre-Processing for A/B Test                 ## 
##                                                                            ##
################################################################################

# input:  SimluatedRekentuinData2.csv
# output: Rekentuin_ABTestData.csv

#-------------------------------------------------------------------------------
#                                                               
# 0. Set Working Directory, Load Data + Load Packages Needed             
#                                                               
#-------------------------------------------------------------------------------

rm(list = ls())

library(dplyr)

setwd("C:/Users/Tabea Hoffmann/Documents/Master/Thesis/Preregistration/")

data <- read.csv2("SimluatedRekentuinData2.csv")

#-------------------------------------------------------------------------------
#                                                               
# 1. Exclude Players + Select Last Game
#                                                               
#-------------------------------------------------------------------------------

# exclude players that did not play any crown games
data2 <- data %>% group_by(user_id) %>% filter(any(game_type == "crown", na.rm = T))

# select last game and exclude players who's only crown game was their last game
data3 <- data2 %>% group_by(user_id)  %>% 
  
  # count frequency of games
  add_tally(game_type == "crown", name = "freq_crown") %>% 
  
  # select last game for each child
  slice(which.max(as.numeric(time))) %>%   
  
  # exclude players 
  filter(!any(freq_crown == 1 && game_type == "crown"))  


#-------------------------------------------------------------------------------
#                                                               
# 2. Create Vectors with the Conversions and N for Each Version
#                                                               
#-------------------------------------------------------------------------------  

# order data according to time
data3.reorder <- data3 %>% arrange(time)

conversion.A <- numeric(nrow(data3))
conversion.B <- numeric(nrow(data3))
n.A <- numeric(nrow(data3))
n.B <- numeric(nrow(data3))

# clicked_crown_game = 1 = crown
# conversion: success = non-crown = 1

for(i in 1:nrow(data3)){
  if(data3.reorder[i,]$clicked_crown_game == 0 && 
     data3.reorder[i,]$variant_id == "con"){
    
    conversion.A[i+1] <- conversion.A[i] + 1
    conversion.B[i+1] <- conversion.B[i] 
    n.A[i+1] <- n.A[i] + 1 
    n.B[i+1] <- n.B[i] 
    
  } else if(data3.reorder[i,]$clicked_crown_game == 1 && 
            data3.reorder[i,]$variant_id == "con"){
    
    conversion.A[i+1] <- conversion.A[i] 
    conversion.B[i+1] <- conversion.B[i] 
    n.A[i+1] <- n.A[i] + 1
    n.B[i+1] <- n.B[i] 
    
  } else if(data3.reorder[i,]$clicked_crown_game == 0 && 
            data3.reorder[i,]$variant_id == "exp"){
    
    conversion.A[i+1] <- conversion.A[i] 
    conversion.B[i+1] <- conversion.B[i] + 1
    n.A[i+1] <- n.A[i] 
    n.B[i+1] <- n.B[i] + 1
    
  } else if(data3.reorder[i,]$clicked_crown_game == 1 && 
            data3.reorder[i,]$variant_id == "exp"){
    
    conversion.A[i+1] <- conversion.A[i] 
    conversion.B[i+1] <- conversion.B[i] 
    n.A[i+1] <- n.A[i] 
    n.B[i+1] <- n.B[i] + 1
    
  }
}

conversion <- list(n1 = n.A, y1 = conversion.A, n2 = n.B, y2 = conversion.B)

#-------------------------------------------------------------------------------
#                                                               
# 3. Export Data for The Analysis in JASP
#                                                               
#-------------------------------------------------------------------------------

write.csv2(as.data.frame(conversion), "Rekentuin_ABTestData2.csv", 
           row.names = FALSE)



