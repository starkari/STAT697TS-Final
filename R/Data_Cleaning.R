
# Cleans Up and Analysis Raw data

library(tidyverse)
library(readr)



# raw data for NOAA
dat_full <- read_csv("Data/weather_data.csv")

# data frame of missing valuse
missing_data <- dat_full %>% 
  filter(is.na(TOBS))

# last missing value for variable of interest is in 2009 so data is from 2010 till
# end of 2021
# add year and month variables
# also add day in year and number of day in year and ratio of day/# days in year
# this assist in sinusoidal regression components
data <- dat_full %>% 
  filter(DATE>=as.Date.character("2010-01-01")) %>% 
  select(c(DATE,TOBS)) %>% 
  mutate(year = format(as.Date(DATE, format="%Y/%m/%d"),"%Y")) %>% 
  mutate(month = format(as.Date(DATE, format="%Y/%m/%d"),"%m")) %>% 
  mutate(day.in.year = lubridate::yday(DATE)) %>% 
  group_by(year) %>% 
  mutate(num.days.in.year=max(day.in.year)) %>% 
  ungroup() %>% 
  mutate(omega = day.in.year/num.days.in.year)

# gets index of the row that is at the 90% mark of my data
ninty.percent.cutoff <- floor(0.9*nrow(data))

data.90.percent <- data[1:ninty.percent.cutoff,]


save(data,file="Data/data.Rdata")
save(data.90.percent,file="Data/data.90.percent.Rdata")
save(missing_data,file="Data/missing_data.Rdata")



