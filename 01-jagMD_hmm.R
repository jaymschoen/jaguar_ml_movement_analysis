# Jaguar Machine Learning Movement Model
## Hidden Markov Models for movement state classification

#### Setup ####
library(sp)
library(sf)
# library(rgdal)
library(tidyverse)
# library(raster)
library(terra)
library(moveHMM)
library(lubridate)

setwd("~/Columbia E3B/Paraguay Jaguar Project/Jaguar Movement Database/")
path <- "C:/Users/jayms/Documents/Columbia E3B/PhD/Dissertation/Ch3-jaguar_movement/"

# save.image("~/Columbia E3B/PhD/Dissertation/Ch3-jaguar_movement/scripts/HMM_workspace.RData")
load("~/Columbia E3B/PhD/Dissertation/Ch3-jaguar_movement/scripts/HMM_workspace.RData")

# library(spatialEco) 

data1 <- read_csv("jaguar_movement_data.csv")
data2 <- read_csv("Jaguar_additional_information_fixed.csv")

data_raw <- left_join(data1, data2, 'ID')

long_lat <- "epsg:4326"

table(data_raw$country) # Mostly in South America

# Limiting analysis to Brazil, Paraguay, and Argentina
data <- data_raw %>%
  filter(country %in% c("Brazil", "Paraguay", "Argentina"))

IDs <- unique(data$ID)            # 108/117 individuals in these three countries

# Visualizing
data_samp_sp <- sample_n(data, 1e4) %>%
  vect(geom = c("location.long", "location.lat"))
# writeVector(data_samp_sp, "~/Columbia E3B/PhD/Dissertation/Ch3-jaguar_movement/data_sample_10k.shp")

SA <- vect("~/Columbia E3B/Paraguay Jaguar Project/Paraguay Jaguar/Paraguay Shape Files/Admin areas/SouthAmerica.shp")

plot(SA)
plot(data_samp_sp, add = T, col = data_samp_sp$ID, pch = 1);box();axis(1);axis(2)


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
  dplyr::summarize(date_range = range(timestamp, na.rm = T))


#### Regularly sampled subset ####

# Checking regularity of sampling (needs to be mostly regular for HMMs)
intervals <- data %>%
  group_by(ID) %>%
  mutate(interval = round(difftime(timestamp, lag(timestamp), units = "hours"), 0),
         int_n = as.numeric(interval)) %>%
  drop_na(int_n)

# Checking long intervals
# View(intervals%>%filter(int_n>1000)) 

# long_int <- c(25, 27, 36, 53, 89)

# jag25 has two unique time periods, both in same area*
# jag27 has 3 time periods, low proportion of regular intervals, looks to use two areas in all time periods
# jag35 and jag36 have very few fixes nad no regular interval
# jag53 has two unique time periods, looks to use two areas in both
# jag89 has two unique time periods, both in same area*

# jags 25 and 89 may need to be separated into two time periods
## Because they are in the same area in both periods, will leave all data together (effect should be minimal)


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

median(ind_interval_props$mod_int)

# Regular intervals (>80% based off of Yeh et al. 2012):
ind_reg_int <- ind_interval_props %>%
  filter(prop_reg0 >0.799,
         ID != 114) %>%               # very irregular movement pattern indicative of GPS error
  mutate(ID = as.character(ID),
         dif2 = prop_reg2 - prop_reg0,
         dif4 = prop_reg4 - prop_reg0)

median(ind_reg_int$dif2); sd(ind_reg_int$dif2)
median(ind_reg_int$dif4); sd(ind_reg_int$dif4)

ind_reg_int

jags_reg_int0 <- c(ind_reg_int$ID)
# jags_reg_int2 <- c(ind_reg_int$ID)
# jags_reg_int4 <- c(ind_reg_int$ID)
# 
setdiff(jags_reg_int4,jags_reg_int2)

# Polygons without buffer: 1,2,3,7,10,12,13

# Lost IDs without buffer on reg_int
# poly1: 80
# poly2: 7,9,10
# poly4: 76,77
# poly5: 73
# poly6: 14,15,19,25,68,69,79,84,86,87,
# poly7: 74,75
# poly8: 105,108,114,115
# poly9: 65
# poly10: 81,88
# poly11: 89
# poly13: 99



# median
median(ind_reg_int$mod_int)
mean(ind_reg_int$mod_int)

# Visualizing
jags_reg_sp <- data%>%
  filter(ID %in% jags_reg_int) %>%
  vect(geom = c("location.long", "location.lat"), crs = long_lat)
plot(SA);plot(jags_reg_sp, add = T, col = jags_reg_sp$ID, pch = 1)

# writeVector(jags_reg_sp, "jags_SA_reg_int.shp", filetype = "ESRI Shapefile")
jags_reg_sp <- vect("ind_reg_samp_st.shp")
as_tibble(jags_reg_sp) %>% filter(year == 1999)

#### Hidden Markov Models ####
# using moveHMM to classify behavioral states

library(moveHMM)

# prepData for entire data set to visualize outputs
full_data <- prepData(data %>%
                        dplyr::select(location.long, 
                                      location.lat,
                                      ID,
                                      year) %>%
                        rename(x = location.long, y = location.lat) %>%
                        mutate(Event_ID = data$Event_ID) %>%
                        as.data.frame()
                      )

str(full_data)
summary(full_data)
hist(full_data$step)
hist(full_data$angle)

# See example
plot(full_data[full_data$ID == "25",])



# function to prep data by individual ID
hmm_prep <- function(data, id) {
  xy <- data %>%
    filter(ID == id) %>%
    dplyr::select(location.long, location.lat) %>%
    rename(x = location.long, y = location.lat)
  s <- data %>%
    filter(ID == id) %>%
    dplyr::select(Sex) %>%
    first()
  a <- data %>%
    filter(ID == id) %>%
    dplyr::select(`Estimated Age`) %>%
    first()
  yr <- data %>%
    filter(ID == id) %>%
    dplyr::select(year)
  d <- prepData(as.data.frame(xy)) %>%
    mutate(ID = str_glue("jag{id}"),
           sex = s,
           age = a,
           year = yr) %>%
    mutate(year = year$year)
  return(d)
}

# applying hmm_prep to all individuals
for (i in unique(data$ID)) {
  assign(str_glue("jag{i}_data"), hmm_prep(data, i))
}

# # applying hmm_prep to all individuals with >=80% regular sampling
# for (i in jags_reg_int) {
#   assign(str_glue("jag{i}_data"), hmm_prep(data, i))
# }



#### Extents/Background ####

##### Buffered convex hulls ####
# of individuals for covariate generation in GEE

# Buffer function for convex hull
ind_poly <- function(id, k) {
  d <- as.matrix(get(str_glue("jag{id}_data"))[4:5])
  v <- vect(d, crs = long_lat)
  h <- convHull(v)
  b <- buffer(h, k) 
  b$ID = id
  # s <- as(b, "Spatial")     # Converting to spatial polygons for easy export (writeVect wasn't working)
  return(b)
}

t <- ind_poly(10, 1e4)
plot(t)

ind_polys_10k_buf <- lapply(jags_reg_int, function(x) ind_poly(x, 1e4))
all_ind_polys_10k_buf <- do.call(rbind, ind_polys_10k_buf) # %>%
  # aggregate()
plot(all_ind_polys_10k_buf, col = all_ind_polys_10k_buf$ID)

## To create individual polygons
# for (i in jags_reg_int) {
#   assign(str_glue("jag{i}_poly_10k_buf"), ind_polys_10k_buf[[i]])
# }

# writeVector(all_ind_polys_10k_buf, str_glue("{path}ind_polys_convex_hull/all_ind_polys_10k_buf.shp"))

# for (i in jags_reg_int){
#   export <- get(str_glue("jag{i}_poly_10k_buf"))
#   filename <- str_glue("~/Columbia E3B/PhD/Dissertation/Ch3-jaguar_movement/ind_polygons_10k_buf/jag{i}_poly_10k_buf.shp")
#   shapefile(export, filename, overwrite = T)
# }




#### HMMs ####

# Visualizing step/turning angle histograms (random sample of 20)
par(mfrow = c(4,5))

set.seed(1234)
for(i in sample(jags_reg_int, 20)) {
  hist(full_data[full_data$ID==i,][,2])
}

for(i in sample(jags_reg_int, 20)) {
  hist(full_data[full_data$ID==i,][,3])
}

par(mfrow=c(1,1))


# 2-state HMMs seem to be missing rare events (exploratory)
## See old code for details/trials



##### 3-state HMMS ####

# function to apply fit 3-state hmm to individuals based on summary metrics
hmm3_ind <- function(data, id) {
  p <- data %>%
    filter(ID == id)
  
  step_quant <- quantile(p$step, probs = seq(0,1,0.05), na.rm = T)
  mu0 <- as.numeric(c(step_quant[5],step_quant[13],step_quant[20])) # step means based in 20/60/95 %iles
  sigma0 <- c(0.05,1,3) # step SD (higher states given more variation)
  zeromass0 <- c(0.1,0.01,1e-3) # step zero-mass (generally rare and unlikely to have any 0 steps in second state)
  stepPar1 <- c(mu0,sigma0,zeromass0) # if 0 distance steps
  stepPar0 <- c(mu0,sigma0)   # if no 0 distance steps
  
  stepPar <- NA
  if (length(which(p$step ==0)) >0) {   # need to remove zero-mass argument when no 0-distance steps
    stepPar = stepPar1
  } else {
    stepPar = stepPar0
  }
  
  angleMean0 <- c(0,pi,pi/2) # angle mean (more turning in state 1)
  kappa0 <- c(1,1,1) # angle concentration
  anglePar0 <- c(angleMean0,kappa0)
  
  m <- fitHMM(data = p, nbStates=3, stepPar0=stepPar, anglePar0=anglePar0)
  return(m)
}

jag7_hmm3 <- hmm3_ind(full_data, 7)
plot(jag7_hmm3)

# applying it to full dataset
# for (i in jags_reg_int) {
#   assign(str_glue("jag{i}_hmm3"),hmm3_ind(full_data, i))
# }

## 3-state HMMs doing good job for most individuals, but several models are perfrom poorly
## However, unsure if I am generating the BEST model


###### Optimized 3-state HMM ####

# Building in model optimization for each individual (based off of package tutorial)
# (Improving chances of finding global maximums over local maximums for parameters)
# Based on animal ID and # iterations used for optimization
# Best models appear to be the most commonly found by hmm algorithm, NOT highest likelihood score

hmm3_ind_opt <- function(id, n_iter) {
  mod_list <- list()
  
  for(i in 1:n_iter){
    d <- get(str_glue("jag{id}_data")) 
    q <- quantile(d$step, probs = seq(0,1,0.05), na.rm = T) # Quantiles by 5% increments
    
    # Step means from random numbers within bounds based on step length quantiles
    stepMean0 <- runif(3,
                       min = c(q[3], q[9], q[14]),
                       max = c(q[9], q[14], q[20])
                       )
    
    # Step length standard deviation  (good models showing SD close to step mean)
    stepSD0 <- stepMean0/1.5
    # Step zero-mass (generally rare and unlikely to have any 0 steps in second/third states)
    zeromass0 <- c(0.001, 1e-4, 1e-5)
    
    # Final step parameters with flexibility for 0 distance steps
    stepPar1 <- c(stepMean0,stepSD0,zeromass0) # if 0 distance steps
    stepPar0 <- c(stepMean0,stepSD0)   # if no 0 distance steps
    stepPar <- NA
      if (length(which(d$step ==0)) >0) {   # need to remove zero-mass argument when no 0-distance steps
        stepPar = stepPar1
      } else {
        stepPar = stepPar0
      }
    
    # Turning angle mean (higher in lower movement states)
    angleMean0 <- c(pi,pi/2,0) 
    # Turning angle concentration (relatively random within these bounds based on initial analysis)
    angleCon0 <- runif(3, 0.1, 2)
    
    # Final angle parameters
    anglePar0 <- c(angleMean0,angleCon0)
    
    # Fitting model (allowing for errors in parameter search space with try())
    mod_list[[i]] <- try(fitHMM(data = get(str_glue("jag{id}_data")), 
                                nbStates = 3, 
                                stepPar0 = stepPar,
                                anglePar0 = anglePar0),
                         silent = TRUE)
  }
  
  # Using likelihood values from list of models to pull best model (most consistently chosen)
  # Trials showing consistent models look better visually than lowest log likelihood models
  
  # First, only allowing models with no errors and 3 states
  mod_list <- mod_list[lapply(mod_list, length)==7]   # If no error, output is list of length = 7 
  n_states <- mod_list %>%
    lapply(viterbi) %>%       # built-in function to classify all points into states based on model
    lapply(unique) %>%
    lapply(length)
  mod_list <- mod_list[which(n_states==3)]    # Only keeping true 3-state models
  
  # Selecting most chosen model
  l <- round(unlist(lapply(mod_list, function(m) m$mod$minimum)), 2)        # log-likelihood scores
  best <- first(which(l == as.numeric(names(which.max(table(l))))))         # most chosen 

  # Optimized model
  opt_m <- mod_list[[best]]
  return(opt_m)
}

# Tests
# system.time(jag7_hmm3_test31 <- hmm3_ind_opt(id = 7, n_iter = 31))
# 
# system.time(jag41_hmm3_test5 <- hmm3_ind_opt(id = 41, n_iter = 5))
# plot(jag41_data)
# plot(jag41_hmm3_test5)
# 
# states <- viterbi(jag41_hmm3_test5)
# prop.table(table(states))


# Applying it to full dataset
# 21 attempts to ensure strong model is found
set.seed(234)   #for replicability

system.time(
  for (i in jags_reg_int) {
  assign(str_glue("jag{i}_hmm3_opt"),hmm3_ind_opt(i, 21))
  }
)
# ~ 4 hours to complete; 100% successful in finding optimum model


# Saving workspace
# save.image(str_glue("{path}scripts/HMM_workspace.RData"))

# Import workspace with all models (all code run up to here):
# load(str_glue("{path}scripts/HMM_workspace.RData"))


##### Proportion of states for all individuals ####

prop_states <- function(id) {
  s <- viterbi(get(str_glue("jag{id}_hmm3_opt")))
  t <- prop.table(table(s))
  return(t)
}

state_prop <- do.call(bind_rows, lapply(jags_reg_int, prop_states)) %>%
  mutate(ID = jags_reg_int) %>%
  dplyr::select(ID,1,2,3) %>%
  rename(prop_1 = '1', prop_2 = '2', prop_3 = '3')

# mean/median state proportions
apply(state_prop[2:4], 2, mean)
apply(state_prop[2:4], 2, median)

# SD/IQR
apply(state_prop[2:4], 2, sd)
apply(state_prop[2:4], 2, IQR)
apply(state_prop[2:4], 2, quantile)


##### Summary stats for step/angle ####

## Mean step per state
# step_m <- full_data %>%
#   filter(ID %in% jags_reg_int) %>%
#   group_by(ID, state) %>%
#   summarize(mean_step = mean(step, na.rm=T)) %>%
#   right_join(ind_reg_int) %>%
#   mutate(m_step_hr = ifelse(mod_int != 0, 
#                              mean_step/mod_int,
#                              mean_step/mean_int))

add_summary <- function(id) {
  s <- viterbi(get(str_glue("jag{id}_hmm3_opt")))
  d <- get(str_glue("jag{id}_data")) %>%
    mutate(ID = as.numeric(str_sub(ID, 4))) %>%
    left_join(ind_interval_props, by = "ID") %>%
    mutate(state = s) 
  mean_step <- d %>%
    group_by(state) %>%
    summarize(mean_step = mean(step, na.rm=T))
  summary <- d %>%
    left_join(mean_step, by = "state") %>%
    mutate(step_hrly = ifelse(mod_int != 0, 
                              mean_step/mod_int,
                              mean_step/mean_int)) 
  return(summary)
}

data_summary <- do.call(bind_rows, lapply(jags_reg_int, add_summary))

# Medians
data_summary %>%
  # group_by(state) %>%
  summarise(across(c(step, step_hrly, angle), ~median(.x, na.rm = T)))

# SDs
data_summary %>%
  # group_by(state) %>%
  summarise(across(c(step, step_hrly, angle), ~sd(.x, na.rm = T)))

#IQR
data_summary %>%
  # group_by(state) %>%
  reframe(across(c(step, step_hrly, angle), ~quantile(.x, probs = seq(0,1,0.1), na.rm = T)))


## Check one interesting data set (jag114): 

plot(jag114_hmm3_opt)
jag114_int_table <- intervals %>% 
  filter(ID==114) %>%
  dplyr::select(int_n) %>%
  table %>%
  prop.table()

# proportion long intervals
intervals %>% 
  filter(ID==114, int_n > 12) %>%
  nrow()/nrow(jag114_data)

# Appears to be ok; ~17% of intervals >7 hours, ~8.6% over 12 hours
# Still a very strange movement pattern; will inspect further in QGIS

### Removing jag114 from analysis ###
