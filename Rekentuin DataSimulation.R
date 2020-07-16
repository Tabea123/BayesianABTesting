################################################################################
##                                                                            ##
##                        Rekentuin Data Simulation                           ## 
##                                                                            ##
################################################################################

# input:  -
# output: SimluatedRekentuinData.csv

#-------------------------------------------------------------------------------
#                                                               
# 0. Set Working Directory + Load All Packages Needed             
#                                                               
#-------------------------------------------------------------------------------

rm(list = ls())

library(truncnorm)

#-------------------------------------------------------------------------------
#                                                               
# 1. Simulate Data for Version A
#                                                               
#-------------------------------------------------------------------------------

## create 100 players that play 4 games ##
user_id <- rep(1:100, each = 4)

## simulate which game (crown vs. non-crown) was played ##

# probability of crown game ranges from 0.5 - 0.8
clicked_crown_game <- rbinom(length(user_id), 1, seq(.5, .8, by = 0.001))

# clicked_crown_game is the variable name in the real data set; it is not
# logical because successes need to correspond to not playing a crown game
game_type <- ifelse(clicked_crown_game == 1, "crown", "non-crown")


## add time stamp; 0 = first game ever played, 90 = 90 secs later ##

# sampling time of first game for each person in a range of 8 hours
first_game  <- rtruncnorm(n = 100, a = 0, b = 28800, mean = 14400, sd = 3600)

# sampling seconds to add for second game from truncated normal distribution
second_game <- first_game + rtruncnorm(n = 100, a = 60, b = 240, mean = 120, sd = 30)

# sampling seconds to add for third game
third_game  <- second_game + rtruncnorm(n = 100, a = 60, b = 240, mean = 120, sd = 30)
# add a minimum of a day to some games
x <- sample(user_id, 10)
third_game[x] <- third_game[x] + 
                 rtruncnorm(n = length(x), a = 86400, b = 100800, 
                            mean = 86520, sd = 28800)

# sampling seconds to add for fourth game
fourth_game <- third_game  + 
              rtruncnorm(n = 100, a = 60, b = 240, mean = 120, sd = 30)

# combine all games
all_games   <- cbind(first_game, second_game, third_game, fourth_game)

# all other games
time <- NULL
for(i in 1:100){
  time <- c(time, as.numeric(all_games[i,]))
}

# first game ever played
time[1] <- 0

# version information
variant_id <- rep("con", length(user_id))

# create data frame with user id, start and end time of play session, and game type 
VersionA <- data.frame(user_id, clicked_crown_game, game_type, time, variant_id)

# some players did not play at some days 
delete <- sample(2:100, 10)
A <- VersionA[-delete,]


#-------------------------------------------------------------------------------
#                                                               
# 2. Simulate Data for Version B
#                                                               
#-------------------------------------------------------------------------------

## create 100 players that play 4 games ##
user_id <- rep(101:200, each = 4)

## simulate which game (crown vs. non-crown) was played ##

# probability of crown game ranges from 0.1-0.4
clicked_crown_game <- rbinom(length(user_id), 1, seq(.1, .4, by = 0.001))

# clicked_crown_game is the variable name in the real data set; it is not
# logical because successes need to correspond to not playing a crown game
game_type <- ifelse(clicked_crown_game == 1, "crown", "non-crown")


## add time stamp; 0 = first game ever played, 90 = 90 secs later ##

# sampling time of first game for each person in a range of 8 hours
first_game <- rtruncnorm(n = 100, a = 1, b = 28800, mean = 14400, sd = 3600)

# sampling seconds to add for second game from truncated normal distribution
# mean is lower, sd higher in this group
second_game <- first_game + 
               rtruncnorm(n = 100, a = 30, b = 180, mean = 80, sd = 120)

# sampling seconds to add for third game
third_game  <- second_game + 
              rtruncnorm(n = 100, a = 30, b = 180, mean = 80, sd = 60)
# additionally add a minimum of a day to some games
x <- sample(1:100, 10)
third_game[x] <- third_game[x] + 
  rtruncnorm(n = length(x), a = 86400, b = 100800, 
             mean = 86520, sd = 28800)

# sampling seconds to add for fourth game
fourth_game <- third_game  + rtruncnorm(n = 100, a = 30, b = 180, mean = 80, sd = 60)

# combine all games
all_games <- cbind(first_game, second_game, third_game, fourth_game)

# all other games
time <- NULL
for(i in 1:100){
  time <- c(time, as.numeric(all_games[i,]))
}

# version information
variant_id <- rep("exp", length(user_id))

# create data frame with user id, start and end time of play session, and game type 
VersionB <- data.frame(user_id, clicked_crown_game, game_type, time, variant_id)

# some players did not play at some days; more in this version
delete <- sample(2:100, 30)
B <- VersionB[-delete,]


#-------------------------------------------------------------------------------
#                                                               
# 3. Create Dataframe With Both Versions
#                                                               
#-------------------------------------------------------------------------------

## merging of the two dataframes
all.data <- rbind(A, B)

#-------------------------------------------------------------------------------
#                                                               
# 4. Export Data
#                                                               
#-------------------------------------------------------------------------------

## export data
write.csv2(all.data, "SimluatedRekentuinData.csv")

