# Living Planet Index -  Calculating CV over multiple moving windows for each .csv

#Author(s): Reilly O'Connor
#Version: 2024-01-15

#Pkgs
library(tidyverse)
library(RColorBrewer)
library(beepr)
library(zoo)
library(RcppRoll)
library(easystats)
library(reshape2)
library(tseries)
library(urca)
library(wql)
library(trend)

#load data
#All TS
#df_lpi <- read.csv("../Mammalian-Fast-Slow-Stability//Data/Living Planet Index/TS_lpi_all.csv", header = T)

#All TS without significant trends
#df_lpi <- read.csv("../Mammalian-Fast-Slow-Stability//Data/Living Planet Index/TS_lpi_notrend.csv", header = T)

#All TS without trends + detrended
df_lpi <- read.csv("../Mammalian-Fast-Slow-Stability//Data/Living Planet Index/TS_lpi_detrended.csv", header = T)

##### CV Calculation #####
#Filter by Maximum number of Years in TS
#df_lpi <- df_lpi %>% filter(unique_years > 19)

unique(df_lpi$Binomial)
unique(df_lpi$ID)

#List of Main IDs
main_ids <- unique(df_lpi$ID)

#define function for calculating cv
cv <- function(x) {
  sd(x) / mean(x)
}

#Initialize the dataframe to store the results
df_lpi_cv <- data.frame(ID = integer(), cv_window = integer(), cv_value = numeric())

for (id in main_ids) {
  
  df_cv_id <- df_lpi %>%
    filter(ID == id)
  
  df_id <- df_cv_id %>% group_by(Year) %>%
    reframe(Population = sum(Value)) %>% 
    arrange(Year) %>%
    dplyr::select(Population)
  
  time_series_length <- nrow(df_id)
  
  #Calculate and store CVs for window sizes from 3 to the minimum of 5 or the time series length
  max_window_size <- min(5, time_series_length)
  
  for (window_size in 3:max_window_size) {
    #Calculate rolling CV for the current window size
    rolling_cv <- rollapply(data = df_id$Population, width = window_size, FUN = cv, by = 1, align = 'center', partial = TRUE)
    mean_cv_for_window <- mean(rolling_cv, na.rm = TRUE)
    
    #Store the results in the dataframe
    df_lpi_cv <- rbind(df_lpi_cv, data.frame(ID = id, cv_window = window_size, cv_value = mean_cv_for_window))
  }
}

ggplot(df_lpi_cv, aes(y = cv_value, x = as.factor(cv_window))) +
  geom_boxplot(outliers = F) +
  theme_classic() +
  ylim(0,1.25)


df_lpi_info <- df_lpi %>% dplyr::select(-X, -Year, -Value)
df_lpi_info_unique <- as.data.frame(unique(df_lpi_info))

df_lpi_cv_fin <- merge(df_lpi_cv, df_lpi_info_unique, by = "ID")

#write.csv(df_lpi_cv_fin, "../Mammalian-Fast-Slow-Stability//Data/lpi mammal population cv all.csv")
#write.csv(df_lpi_cv_fin, "../Mammalian-Fast-Slow-Stability//Data/lpi mammal population cv no trends.csv")
#write.csv(df_lpi_cv_fin, "../Mammalian-Fast-Slow-Stability//Data/lpi mammal population cv detrended.csv")
