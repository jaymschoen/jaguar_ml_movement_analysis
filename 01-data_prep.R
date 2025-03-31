# Jaguar Machine Learning Movement Model
## 01 - Data import and cleaning

#### Setup ####
library(sp)
library(sf)
library(tidyverse)
library(terra)
library(moveHMM)
library(lubridate)

setwd("~/Columbia E3B/PhD/Dissertation/Ch3-jaguar_movement/scripts/git_push/")

data1 <- read_csv("jaguar_movement_data.csv")
data2 <- read_csv("Jaguar_additional_information_fixed.csv")

data_raw <- left_join(data1, data2, 'ID')

long_lat <- "epsg:4326"

prop.table(table(data_raw$country)) # Mostly in South America

# Limiting analysis to Brazil, Paraguay, and Argentina
data <- data_raw %>%
  filter(country %in% c("Brazil", "Paraguay", "Argentina"))

IDs <- unique(data$ID)            # 108/117 individuals in these three countries

#### Date ranges for individuals ####
data <- data %>%
  mutate(timestamp = as.POSIXct(timestamp, format = "%m/%d/%Y %H:%M")) %>%
  mutate(date = format(.$timestamp, format = "%m/%d/%Y"),
         year = format(.$timestamp, format = "%Y")) %>%
  mutate(year = as.numeric(year))

# filling in 16 pts with NA in yr (appears to be rare error in lubridate functions)
for(i in 1:nrow(data)){                  
  if(is.na(data$year[i])) data$year[i] = data$year[i-1]
}

ind_date_ranges <- data %>%
  group_by(ID) %>%
  dplyr::reframe(date_range = range(timestamp, na.rm = T)) %>%
  mutate(start_end = rep(c("start", "end"), nrow(.)/2))


#### Regularly sampled subset ####

# Checking regularity of sampling (needs to be mostly regular for HMMs)
intervals <- data %>%
  group_by(ID) %>%
  mutate(interval = round(difftime(timestamp, lag(timestamp), units = "hours"), 0),
         int_n = as.numeric(interval)) %>%
  drop_na(int_n)

# Checking long intervals
View(intervals%>%filter(int_n>1000))

# long_int <- c(25, 27, 36, 53, 89)

# jag25 has two unique time periods, both in same area*
# jag27 has 3 time periods, low proportion of regular intervals, looks to use two areas in all time periods
# jag35 and jag36 have very few fixes nad no regular interval
# jag53 has two unique time periods, looks to use two areas in both
# jag89 has two unique time periods, both in same area*


## Based on below analysis, jags 25 and 89 are the only individuals in jags_reg_int subset
## Because they are in the same area in both periods, will leave all data together (effect should be minimal)
## Minimal land cover change between years of their data (visualized in QGIS), so either year for covariates will work

# Function to get proportion of intervals for individuals
## based on id = individual and b1/2 = buffer times around mode interval

interval_prop <- function(id, b1, b2) {
  j <- intervals %>% 
    filter(ID == id) %>%
    dplyr::select(ID, `Planned Schedule (h)`, int_n)
  mn <- mean(j$int_n, na.rm = T)
  med <- median(j$int_n, na.rm = T)
  mod <- as.numeric(names(sort(-table(j$int_n)))[1])
  buf1 <- c(seq(mod-b1,mod+b1,1))       # allowing fixes +/- "b" hours of mode interval time
  buf2 <- c(seq(mod-b2,mod+b2,1))       # larger buffer
  s0 <- filter(j, int_n == mod)       # fixes == mod
  s1 <- filter(j, int_n %in% buf1)     # fixes == mod + buffer 
  s2 <- filter(j, int_n %in% buf2)     # fixes == mod + buffer 
  p0 <- nrow(s0)/nrow(j)
  p1 <- nrow(s1)/nrow(j)
  p2 <- nrow(s2)/nrow(j)
  df <- data.frame(ID = id, 
                   planned_int = as.numeric(first(j$`Planned Schedule (h)`)),
                   mean_int = mn,
                   med_int = med, 
                   mod_int = mod,
                   prop_reg0 = p0,           # == mod
                   prop_reg1 = p1,           # == mod + buffer1
                   prop_reg2 = p2)           # == mod + buffer2
  return(df)
}


# applying to all individuals
ind_interval_props <- data.frame(ID = NA, 
                                 planned_int = NA,
                                 mean_int = NA,
                                 med_int = NA,
                                 mod_int = NA,
                                 prop_reg0 = NA,     # proportion regular with no buffer
                                 prop_reg2 = NA,     # prop reg with 2 hr buffer
                                 prop_reg4 = NA)     # prop reg with 4 hr buffer
for (i in IDs) {
  ind_interval_props[i,] <- interval_prop(i,2,4)
}
ind_interval_props <- drop_na(ind_interval_props, ID)

# Regular intervals (>80% based off of Yeh et al. 2012):
ind_reg_int <- ind_interval_props %>%
  filter(prop_reg4 >0.799,
         ID != 114) %>%               # very irregular movement pattern indicative of GPS error
  mutate(ID = as.character(ID),
         dif2 = prop_reg2 - prop_reg0,
         dif4 = prop_reg4 - prop_reg0)

# Checking sensitivity based on interval buffer allowed for "regular"
median(ind_reg_int$dif2); sd(ind_reg_int$dif2)
median(ind_reg_int$dif4); sd(ind_reg_int$dif4)

# 4 hour buffer strikes balance between strict regularity and including individuals with >80% regularity

# median
median(ind_reg_int$mod_int)
mean(ind_reg_int$mod_int)

# IDs of jaguars with regularly sampled data as defined above
## will be used in movement analysis
jags_reg_int <- unique(ind_reg_int$ID)
