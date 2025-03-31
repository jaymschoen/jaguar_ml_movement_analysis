# Jaguar Machine Learning Movement Model
## 02 - Building Hidden Markov Models for movement state classification

#### Setup ####

setwd("~/Columbia E3B/PhD/Dissertation/Ch3-jaguar_movement/scripts/git_push/")

source("01-data_prep.R")


#### Hidden Markov Models ####
# using moveHMM to classify behavioral states

# using only regularly sampled jaguars (as defined in 01-data_prep)
data_reg <- filter(data, ID %in% jags_reg_int)

# prepData for entire data set to visualize outputs
full_data <- prepData(data_reg %>%
                        dplyr::select(location.long, 
                                      location.lat,
                                      ID,
                                      year) %>%
                        rename(x = location.long, y = location.lat) %>%
                        mutate(Event_ID = data_reg$Event_ID) %>%
                        as.data.frame()
                      )

str(full_data)
summary(full_data)
hist(full_data$step)
hist(full_data$angle)

# See example
plot(full_data[full_data$ID == "25",])



# Function to prep data by individual ID
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
for (i in unique(data_reg$ID)) {
  assign(str_glue("jag{i}_data"), hmm_prep(data_reg, i))
}


#### Extents/Background ####

##### Buffered convex hulls ####
# Areas of individuals' movement for covariate generation in Google Earth Engine

# Buffer function for convex hull
ind_poly <- function(id, k) {
  d <- as.matrix(get(str_glue("jag{id}_data"))[4:5])
  v <- vect(d, crs = long_lat)
  h <- convHull(v)
  b <- buffer(h, k) 
  b$ID = id
  return(b)
}

# Creating 10km buffers for each individual

ind_polys_10k_buf <- lapply(jags_reg_int, function(x) ind_poly(x, 1e4))

all_ind_polys_10k_buf <- do.call(rbind, ind_polys_10k_buf)

# To create individual polygons
for (i in jags_reg_int) {
  assign(str_glue("jag{i}_poly_10k_buf"), ind_polys_10k_buf[[i]])
}

# writeVector(all_ind_polys_10k_buf, str_glue("{path}ind_polys_convex_hull/all_ind_polys_10k_buf.shp"))

# for (i in jags_reg_int){
#   export <- get(str_glue("jag{i}_poly_10k_buf"))
#   filename <- str_glue("{path}ind_polygons_10k_buf/jag{i}_poly_10k_buf.shp")
#   shapefile(export, filename, overwrite = T)
# }


#### Building HMMs ####

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


# 2-state HMMs missed rarer movement events (i.e., "exploratory" movements)
## See old code for details/trials

##### 3-state HMMS ####

# function to fit 3-state hmm to individuals based on summary metrics
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

# Applying it to full dataset (will take time to run)
for (i in jags_reg_int) {
  assign(str_glue("jag{i}_hmm3"), hmm3_ind(full_data, i))
}

## 3-state HMMs doing good job for most individuals, but several models perform poorly
## Poor performance indicated by irregularly fit gamma/von mises curves
## Need to generate the BEST model possible


###### Optimized 3-state HMM ####

# Building in model optimization for each individual (based off of moveHMM package tutorial)
# (Improving chances of finding global maximums over local maximums for parameters)
# Based on animal ID and # iterations used for optimization
# Best models appear to be the most commonly found by hmm algorithm, NOT highest likelihood score (which often returns irregular curves)

hmm3_ind_opt <- function(id, n_iter) {
  mod_list <- list()
  
  # Building n_iter number of models to select from
  for(i in 1:n_iter){
    d <- get(str_glue("jag{id}_data")) 
    q <- quantile(d$step, probs = seq(0,1,0.05), na.rm = T) # Quantiles by 5% increments
    
    # Step means from random numbers within bounds based on step length quantiles
    stepMean0 <- runif(n = 3,
                       min = c(q[3], q[9], q[14]),
                       max = c(q[9], q[14], q[20]))
    
    # Step length standard deviation  (good models showing SD close to step mean)
    stepSD0 <- stepMean0/1.5
    # Step zero-mass (generally rare and unlikely to have any 0 steps in second/third states)
    zeromass0 <- c(0.001, 1e-4, 1e-5)
    
    # Final step parameters with flexibility for 0 distance steps
    stepPar1 <- c(stepMean0,stepSD0,zeromass0) # if 0 distance steps
    stepPar0 <- c(stepMean0,stepSD0)   # if no 0 distance steps
    stepPar <- NA
      if (length(which(d$step ==0)) >0) {   # need to remove zero-mass argument when no 0-distance steps
        stepPar = stepPar1}
      else {
        stepPar = stepPar0}
    
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
  
  # Trials showing consistent models look better visually (better curve fits) than lowest log-likelihood models
  
  # Only allowing models with no errors and 3 states
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
system.time(jag41_hmm3_test5 <- hmm3_ind_opt(id = 41, n_iter = 5))
plot(jag41_data)
plot(jag41_hmm3_test5)

states <- viterbi(jag41_hmm3_test5)
prop.table(table(states))


# Applying it to full dataset (this will take a while)
set.seed(234)   # for replication

system.time(
  for (i in jags_reg_int) {
  assign(str_glue("jag{i}_hmm3_opt"),hmm3_ind_opt(i, 21)) # 21 attempts to ensure strong model is found
    } 
)

# ~ 4 hours to complete on my device (Windows Surface 3 Laptop, 16GB RAM)
# 100% successful in finding optimum model


# Saving optimized model results as R workspace
# save.image(str_glue("{path}scripts/HMM_workspace.RData"))

# Import workspace with all models (all code run up to here):
load(str_glue("../HMM_workspace.RData"))


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
add_state <- function(id) {
  s <- viterbi(get(str_glue("jag{id}_hmm3_opt")))
  d <- get(str_glue("jag{id}_data")) %>%
    mutate(state = s)
  return(d)
}

data_state <- do.call(bind_rows, lapply(jags_reg_int, add_state))

# Medians
data_state %>%
  group_by(state) %>%
  summarise(across(c(step, angle), ~median(.x, na.rm = T)))

# SDs
data_state %>%
  group_by(state) %>%
  summarise(across(c(step, angle), ~sd(.x, na.rm = T)))

#IQR
data_state %>%
  group_by(state) %>%
  reframe(across(c(step, angle), ~quantile(.x, na.rm = T)))