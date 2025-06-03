# Jaguar Machine Learning Movement Model
## 03 - Background Areas

### Using optimized HMM states to create background areas for models

#### Setup ####

# library(sp)
# library(sf)
library(tidyverse)
# library(raster)
library(terra)
library(moveHMM)

# Set to directory with cloned repo
setwd("~/Columbia E3B/PhD/Dissertation/Ch3-jaguar_movement/scripts/git_push_R/")

# Loading data from 02-hmm_building
load("~/Columbia E3B/PhD/Dissertation/Ch3-jaguar_movement/scripts/HMM_workspace.RData")

#### Adding movement states ####
# Appending data with movement states 

# Updating each individual's data with predicted movement state
for(i in jags_reg_int) {
  assign(str_glue("jag{i}_data"),
         get(str_glue("jag{i}_data")) %>%
           as_tibble() %>%
           mutate(s = viterbi(get(str_glue("jag{i}_hmm3_opt"))),
                  ID = as.factor(str_sub(as.character(ID), 4)))
  )
}

# Adding to full dataset
states <- data %>%
  filter(ID %in% jags_reg_int) %>%
  rename(x = location.long,
         y = location.lat,
         sex = Sex,
         age = `Estimated Age`) %>% 
  dplyr::select(c(1,3:4,7,10:11,22)) %>%
  mutate(s = NA,
         ID = as.factor(ID))

for(i in jags_reg_int) {
  d <- get(str_glue("jag{i}_data")) %>%
    as.data.frame() %>%
    # mutate(s = viterbi(get(str_glue("jag{i}_hmm3_opt"))),
    #        ID = as.factor(str_sub(as.character(ID), 4))) %>%
    dplyr::select(-c(step, angle))
  states[states$ID == i,] <- right_join(states[states$ID == i,], d)
}

# re-adding Event_ID for alignment with original data 
states <- mutate(states, 
                 Event_ID = as.data.frame(data)[data$ID %in% jags_reg_int, "Event_ID"]) %>%
  left_join(full_data)
str(states) 

# Exporting data with all regularly sampled individuals with behavioral states from HMMs
# write_csv(states, str_glue("data/ind_reg_samp_st.csv"))

# Spatial
states_sp <- states %>%
  vect(geom = c("x", "y"), crs = long_lat)

# writeVector(states_sp, str_glue("data/ind_reg_samp_st.shp"))

#### Presence points ####

# Functions to make spatial vectors for each individual
all_pts <- function(id) {
  m <- get(str_glue("jag{id}_hmm3_opt"))    # getting ID specific 3-state model
  states <- viterbi(m) 
  pts <- m$data %>%                         # selecting all points (all movement states)
    as_tibble() %>%
    mutate(Event_ID = as.data.frame(data)[data$ID==id, "Event_ID"],
           year = as.data.frame(data)[data$ID==id, "year"],
           state = states)  %>%
    dplyr::select(5,1:length(.)) %>%
    vect(geom = c("x", "y"), crs = long_lat) 
  return(pts)
}

## Longest step length paths considered "exploratory"
# Function to pull points from maximum mean step length state
expl_pts <- function(id) {
  m <- get(str_glue("jag{id}_hmm3_opt"))    # getting ID specific 3-state model
  e_state <- names(which.max(m$mle$stepPar[1,])) %>%
    str_sub("-1")                           # pulling the number of state with largest mean step length
  states <- viterbi(m) 
  pts <- m$data %>%                         # selecting points in that movement state (exploratory)
    as_tibble() %>%
    mutate(Event_ID = as.data.frame(data)[data$ID==id, "Event_ID"],
           year = as.data.frame(data)[data$ID==id, "year"],
           state = states) %>%
    filter(state == as.integer(e_state)) %>%
    dplyr::select(5,1:length(.)) %>%
    vect(geom = c("x", "y"), crs = long_lat) 
  return(pts)
}

# testing

jag6_all_pts <- all_pts(6)
summary(jag6_all_pts)
plot(jag6_all_pts)
plot(as.lines(jag6_all_pts), add=T)

jag6_expl_pts <- expl_pts(6)
summary(jag6_expl_pts)
prop_states(6)
plot(jag6_expl_pts, col = 'red', add = T)
plot(as.lines(jag6_all_pts), col = 'darkred', add=T)

# Non-exploratory pts
jag6_ne_pts<- jag6_all_pts %>%
  subset(!(.$Event_ID %in% jag6_expl_pts$Event_ID))
plot(jag6_ne_pts, col = 'blue', add = T)

# applying to full dataset
## all/expl points
for (i in jags_reg_int) {
  assign(str_glue("jag{i}_all_pts"), all_pts(i))
  assign(str_glue("jag{i}_all_lines"), all_pts(i) %>% as.lines())
  assign(str_glue("jag{i}_expl_pts"), expl_pts(i))
  assign(str_glue("jag{i}_expl_lines"), all_pts(i) %>% as.lines())
}

## non-expl points
for (i in jags_reg_int) {
  a <- get(str_glue("jag{i}_all_pts"))
  e <- get(str_glue("jag{i}_expl_pts"))
  ne <- a %>%
    subset(!(.$Event_ID %in% e$Event_ID))
  
  assign(str_glue("jag{i}_ne_pts"), ne)
  assign(str_glue("jag{i}_ne_lines"), get(str_glue("jag{i}_ne_pts")) %>% as.lines())
}

# test
plot(jag6_all_pts)
plot(jag6_all_lines, add = T)
plot(jag6_expl_pts, col = "red", add = T)
plot(jag6_expl_lines, col = "darkred", add = T)
plot(jag6_ne_pts, col = "blue", add = T)
plot(jag6_ne_lines, col = "darkblue", add = T)

## Lines used just for visualization
## Due to behavioral segmentation/long intervals, only points will be used in models (point selection rather than path)


# Bind all individual exploratory/non-exploratory points
expl_list <- list()
for (i in jags_reg_int) {
  expl_list[[i]] <- get(str_glue("jag{i}_expl_pts"))
}

ne_list <- list()
for (i in jags_reg_int) {
  ne_list[[i]] <- get(str_glue("jag{i}_ne_pts"))
}

expl_pts_sp <- vect(expl_list)
ne_pts_sp <- vect(ne_list)


# Export as .csv to use for point data extraction
expl_pts_df <- as_tibble(expl_pts_sp) %>%
  mutate(ID = as.integer(str_sub(as.character(ID),4)),
         x = geom(expl_pts_sp)[,3],
         y = geom(expl_pts_sp)[,4]) %>%
  dplyr::select(Event_ID, ID, x, y, step, angle, state)
# write_csv(expl_pts_df, str_glue("data/ind_reg_samp_expl.csv"))

ne_pts_df <- as_tibble(ne_pts_sp) %>%
  mutate(ID = as.integer(str_sub(as.character(ID),4)),
         x = geom(ne_pts_sp)[,3],
         y = geom(ne_pts_sp)[,4]) %>%
  dplyr::select(Event_ID, ID, x, y, step, angle, state)
# write_csv(ne_pts_df, str_glue("data/ind_reg_samp_ne.csv"))



#### Background pts ####

# Buffering each pt for each individual based on mean step length
## Buffering at 6 * mean hourly step length, then merging for individual background area
### Because mean hourly step is close to resolution (250m) of spatial data,
### Using mean 6 hr step length to connect movement buffers (avoid islands)
### Based on visual inspection of several buffer sizes


# Mean hourly step lengths for regularly sampled jags
step_m <- full_data %>%
  filter(ID %in% jags_reg_int) %>%
  group_by(ID) %>%
  summarize(mean_step = mean(step, na.rm=T)) %>%
  right_join(ind_reg_int) %>%
  mutate(m_step_dly = ifelse(mod_int != 0, 
                             mean_step * 24/mod_int,
                             mean_step * 24/mean_int),
         m_step_hr = m_step_dly/24)

# Summary stats
str(step_m)
summary(step_m)


# Buffer equal to 6 * mean hourly step length (buffering each each point)
## Using "id" = individual id and "st" = movement state
ind_poly_pt <- function(id, st) {
  data <- get(str_glue("jag{id}_{st}_pts"))
  buf <- as.numeric(step_m[step_m$ID == id, "m_step_hr"] *6000)     # to km * 6 hr
  
  pbufs <- list()
  for(i in 1:length(data)) {
    pbufs[[i]] <- buffer(data[i], buf)
  }
  
  merged <- do.call(rbind, pbufs)
  dissolved <- aggregate(merged)
  dissolved$ID <- id
  
  return(dissolved)
}

# Testing
t86_all <- ind_poly_pt(86, "all")
t86_expl <- ind_poly_pt(86, "expl")
t86_ne <- ind_poly_pt(86, "ne")
plot(t86_ne, lwd = 1.5)
plot(t86_expl, border = 'red', lwd = 1.5, add = T)
plot(jag86_ne_pts, add = T)
plot(jag86_expl_pts, col = "red", add = T)

# Applying to all individuals
polys_6hr_buf_all_list <- lapply(jags_reg_int, function(x) ind_poly_pt(x, "all"))
polys_6hr_buf_expl_list <- lapply(jags_reg_int, function(x) ind_poly_pt(x, "expl"))
polys_6hr_buf_ne_list <- lapply(jags_reg_int, function(x) ind_poly_pt(x, "ne"))

polys_6hr_buf_all <- do.call(rbind, polys_6hr_buf_all_list)
polys_6hr_buf_expl <- do.call(rbind, polys_6hr_buf_expl_list)
polys_6hr_buf_ne <- do.call(rbind, polys_6hr_buf_ne_list)

# Creating individual shapes
for (i in 1:length(jags_reg_int)) {
  id <- jags_reg_int[i]
  assign(str_glue("jag{id}_polys_6hr_buf_all"), 
         polys_6hr_buf_all_list[[i]])
  assign(str_glue("jag{id}_polys_6hr_buf_expl"), 
         polys_6hr_buf_expl_list[[i]])
  assign(str_glue("jag{id}_polys_6hr_buf_ne"), 
         polys_6hr_buf_ne_list[[i]])
}


# Generating random (background) points within buffer = 6*mean_step_hr for each point

## 1:1 random sample equally distributed throughout buffer shapes (for each individual)
# calculating area (to distribute points evenly)
# polys_6hr_buf_all_agg <- terra::aggregate(polys_6hr_buf_all) %>% disagg()
polys_6hr_buf_all$id_p <- 1:length(polys_6hr_buf_all)
polys_6hr_buf_all$area <- expanse(polys_6hr_buf_all, unit = "km")
sum(polys_6hr_buf_all$area)

### Using polys from individual 6 hr bufs (id), polygon vector (p_vect), total points (t_pt), sampling method (m)
ind_bg_prop <- function(id, p_vect, t_pt, m) {
  l <- get(p_vect)                # state-specific polygon vector (all, expl, non-expl)
  p <- subset(l, l$ID == id)      # target polygon (based on jag id)
  
  # polygon area proportional to total
  t_area <- sum(l$area)           # total area
  p_area <- p$area/t_area         # proportion of individuals bg area/total area
  n <- (t_pt*p_area) %>% round(0) # number of points proportional to area of individual's bg area
  
  pts <- get(str_glue("jag{id}_all_pts"))
  y <- mean(pts$year)%>%round(0)
  
  set.seed(234)
  ### terra::spatSample showed North bias at time of analysis; used sf::st_sample instead
  # samp <- spatSample(p, n, method = m)
  samp_sf <- sf::st_sample(sf::st_as_sf(p), size = n, type = m)
  samp <- vect(samp_sf)
  samp$year <- y
  samp$ID <- id
  
  return(samp)
}

t <- ind_bg_prop(6, "polys_6hr_buf_all", length(states_sp), "random")
plot(jag6_all_pts); plot(t, col = "red", cex = 0.4, add = T); plot(polys_6hr_buf_all[2], add = T)

# Applying to all individuals
bg_rand_prop_6hr_buf_all_list <- lapply(jags_reg_int, function(x) ind_bg_prop(x, "polys_6hr_buf_all", length(states_sp), "random"))
bg_rand_prop_6hr_buf_all <- do.call(rbind, bg_rand_prop_6hr_buf_all_list)

# Creating individual vectors
for (i in 1:length(jags_reg_int)) {
  id <- jags_reg_int[i]
  assign(str_glue("jag{id}_bg_rand_prop_6hr_buf_all"), 
         bg_rand_prop_6hr_buf_all_list[[i]])
}

# tests
plot(subset(polys_6hr_buf_all, polys_6hr_buf_all$ID==14))
plot(jag14_all_pts, add = T)
plot(jag14_bg_rand_prop_6hr_buf_all, col = "red", cex = 0.4, add = T)

plot(subset(polys_6hr_buf_all, polys_6hr_buf_all$ID==10))
plot(jag10_all_pts, add = T)
plot(jag10_bg_rand_prop_6hr_buf_all, col = "red", cex = 0.4, add = T)

# # Writing to disk for next steps
# writeVector(bg_rand_prop_6hr_buf_all, "data/bg_rand_prop_step_6hr_mean_buf_all.shp", overwrite = T)
