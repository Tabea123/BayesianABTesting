################################################################################
##                                                                            ##
##                 Rekentuin Data Pre-Processing for T-Test                   ## 
##                                                                            ##
################################################################################

# input:  SimluatedRekentuinData2.csv
# output: Rekentuin_TTestData2.csv

#-------------------------------------------------------------------------------
#                                                               
# 0. Set Working Directory, Load Data + Load Packages Needed                              
#                                                               
#-------------------------------------------------------------------------------

rm(list = ls())

setwd("C:/Users/Tabea Hoffmann/Documents/Master/Thesis/Preregistration/")

library(dplyr)
library(tibble)

data <- read.csv2("SimluatedRekentuinData2.csv")

#-------------------------------------------------------------------------------
#                                                               
# 1. Exclude Players 
#                                                               
#-------------------------------------------------------------------------------

# exclude players that did not play one crown game during the time of testing
data2 <- data %>% group_by(user_id) %>% filter(any(game_type == "crown", na.rm = T))

# identify players who's only crown game was their last game
rm <- data2 %>% group_by(user_id)  %>% 
  
  # count frequency of games
  add_tally(game_type == "crown", name = "freq_crown") %>% 
  
  # select last game for each child
  slice(which.max(as.numeric(time))) %>%   
  
  # select players who's only crown game was their last game
  filter(any(freq_crown == 1 && game_type == "crown"))  


# exclude players who's only crown game was their last game
data3 <- data2 %>% filter(!user_id %in% rm$user_id)


#-------------------------------------------------------------------------------
#                                                               
# 2. Calculate Playtime 
#                                                               
#-------------------------------------------------------------------------------

# substract previous game start from current game start per user 
playtime1 <- data3 %>% group_by(user_id) %>% 
  transmute(playtime1 = time - lag(time))
# delete first row to move every value one up to match game with it's playtime
playtime2 <- playtime1 %>% ungroup() %>% select(playtime1) %>% slice(-1) 
# add NA to last 
playtime3 <- c(as.numeric(unlist(playtime2)), NA)

# add playtime as column to data frame
data4 <- data3 %>% add_column(playtime = playtime3)

# indicate NA when playtime is too long 
# because it likely means that they left the website and returned for a new session
data4 %>% filter(playtime > 600)
data5 <- data4 %>% mutate(playtime = replace(playtime, playtime > 600, NA))

# note: this is a crude approximation of a user's playtime as we cannot
# know when a user ended his/her last game

#-------------------------------------------------------------------------------
#                                                               
# 3. Sort Data According to the Last Game Played
#                                                               
#-------------------------------------------------------------------------------

# select last game and reverse order
# first finished user comes first, last finished comes last
sorted <- data3 %>% 
  group_by(user_id) %>%  
  slice(which.max(as.numeric(time))) %>% 
  arrange(time)

users_sorted <- sorted$user_id

#-------------------------------------------------------------------------------
#                                                               
# 4. Sum up Playtime for Crown Games and All Games
#                                                               
#-------------------------------------------------------------------------------

# sum up time spent on crown games per user
onlycrown <- data5 %>% 
  group_by(user_id) %>% 
  filter(game_type == "crown") %>% 
  summarise(sum_crown = sum(playtime, na.rm = T))

# sort summed up time sequentially
crowntime <- onlycrown %>% slice(match(users_sorted, user_id))

# sum up time spent on all games
allgames <- data5 %>% 
  group_by(user_id) %>% 
  summarise(sum_all = sum(playtime, na.rm = T))

# sort summed up time sequentially
alltime  <- allgames %>% slice(match(users_sorted, user_id))

# divide time spent on crown games by the time spent on all games for 
# each user
proportion.crown <- as.numeric(crowntime$sum_crown)/as.numeric(alltime$sum_all)

#-------------------------------------------------------------------------------
#                                                               
# 5. Create Dataframe 
#                                                               
#-------------------------------------------------------------------------------

df <- data.frame(user_id = users_sorted, 
                 crown_playtime = as.numeric(crowntime$sum_crown),
                 all_playtime = as.numeric(alltime$sum_all),
                 proportion_crown = proportion.crown)

# include version information for analysis in JASP
df$variant_id <- ifelse(df$user_id > 100, "exp", "con")


#-------------------------------------------------------------------------------
#                                                               
# 6. Export Data for the Analysis in JASP
#                                                               
#-------------------------------------------------------------------------------

write.csv2(df, "Rekentuin_TTestData2.csv")

