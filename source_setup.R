# Jaguar Movement Database Analysis

#### Setup ####
library(tidyverse)
library(terra)
library(lubridate)


setwd("~/Columbia E3B/PhD/Dissertation/Ch3-jaguar_movement/scripts/git_push_R/")

data1 <- read_csv("data/jaguar_movement_data.csv")
data2 <- read_csv("data/Jaguar_additional_information_fixed.csv")

data_raw <- left_join(data1, data2, 'ID')

# Limiting initial analysis to Brazil, Paraguay, and Argentina
data <- data_raw %>%
  filter(country %in% c("Brazil", "Paraguay", "Argentina"))

IDs <- unique(data$ID)            # 108/117 jaguars in these three countries

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
  dplyr::reframe(date_range = range(timestamp, na.rm = T))

#### Regularly sampled subset ####

# Checking regularity of sampling (needs to be mostly regular for HMMs)
intervals <- data %>%
  group_by(ID) %>%
  mutate(interval = round(difftime(timestamp, lag(timestamp), units = "hours"), 0),
         int_n = as.numeric(interval)) %>%
  drop_na(int_n)

# Function to get proportion of intervals for individuals

interval_prop <- function(id) {
  j <- intervals %>% 
    filter(ID == id) %>%
    dplyr::select(ID, `Planned Schedule (h)`, int_n)
  mn <- mean(j$int_n, na.rm = T)
  med <- median(j$int_n, na.rm = T)
  mod <- as.numeric(names(sort(-table(j$int_n)))[1])
  r <- c(seq(med-4,med+4,1))                      # allowing fixes +/- 4 hour of mode interval time
  s <- filter(j, int_n %in% r)
  p <- nrow(s)/nrow(j)
  df <- data.frame(ID = id, 
                   planned_int = as.numeric(first(j$`Planned Schedule (h)`)),
                   mean_int = mn,
                   med_int = med, 
                   mod_int = mod,
                   prop_reg = p)
  return(df)
}


# applying to all individuals
ind_interval_props <- data.frame(ID = NA, 
                                 planned_int = NA,
                                 mean_int = NA,
                                 med_int = NA,
                                 mod_int = NA,
                                 prop_reg = NA)
for (i in IDs) {
  ind_interval_props[i,] <- interval_prop(i)
}

# Regular intervals (>80% based off of Yeh et al. 2012):
ind_reg_int <- ind_interval_props %>%
  filter(prop_reg >0.799,
         ID != 114) %>%               # very irregular movement pattern indicative of GPS error
  mutate(ID = as.character(ID))
jags_reg_int <- c(ind_reg_int$ID) 
