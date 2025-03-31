# 02-states_bg_prep
## Using optimized HMM states to create spatial layers and background points for model

library(sp)
library(sf)
library(tidyverse)
library(raster)
library(terra)
library(moveHMM)

setwd("~/Columbia E3B/PhD/Dissertation/Ch3-jaguar_movement/")
path <- "C:/Users/jayms/Documents/Columbia E3B/PhD/Dissertation/Ch3-jaguar_movement/"

# Saving/loading data from 01-jagMD-hmm
# save.image("~/Columbia E3B/PhD/Dissertation/Ch3-jaguar_movement/scripts/HMM_workspace.RData")
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

# write_csv(states, str_glue("{path}ind_reg_samp_st.csv"))

# Spatial
states_sp <- states %>%
  vect(geom = c("x", "y"), crs = long_lat)
# writeVector(states_sp, str_glue("{path}ind_reg_samp_st.shp"), overwrite = T)

states_sp <- vect("ind_reg_samp_st.shp")

#### Presence points ####

# Functions to make spatial vectors for each individual
all_pts <- function(id) {
  m <- get(str_glue("jag{id}_hmm3_opt"))    # getting ID specific 3-state model
  states <- viterbi(m) 
  pts <- m$data %>%                     # selecting all points (all movement states)
    as_tibble() %>%
    mutate(Event_ID = as.data.frame(data)[data$ID==id, "Event_ID"],
           year = as.data.frame(data)[data$ID==id, "year"],
           state = states)  %>%
    dplyr::select(5,1:length(.)) %>%
    vect(geom = c("x", "y"), crs = long_lat) 
  return(pts)
}

# Longest step length paths considered "exploratory"
# function to pull points from maximum mean step length state
expl_pts <- function(id) {
  m <- get(str_glue("jag{id}_hmm3_opt"))    # getting ID specific 3-state model
  e_state <- names(which.max(m$mle$stepPar[1,])) %>%
    str_sub("-1")                       # pulling the number of state with largest mean step length
  states <- viterbi(m) 
  pts <- m$data %>%                     # selecting points in that movement state (exploratory)
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
plot(jag6_expl_pts)
plot(as.lines(jag6_expl_pts), add=T)

# Non-exploratory pts
t_ne <- jag6_all_pts %>%
  subset(!(.$Event_ID %in% jag6_expl_pts$Event_ID))

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
  assign(str_glue("jag{i}_ne_lines"),
         get(str_glue("jag{i}_ne_pts")) %>% as.lines())
}

# test
plot(jag5_all_pts)
plot(jag5_all_lines, add = T)
plot(jag5_expl_pts, col = "red", add = T)
plot(jag5_expl_lines, col = "darkred", add = T)
plot(jag5_ne_pts, col = "blue", add = T)
plot(jag5_ne_lines, col = "darkblue", add = T)



# Exporting
# for (i in jags_reg_int){
#   export <- get(str_glue("jag{i}_all_pts"))
#   filename <- str_glue("{path}ind_all_pts_states/jag{i}_all_pts.shp")
#   writeVector(export, filename, overwrite = T)
# }
# 
# for (i in jags_reg_int){
#   export <- get(str_glue("jag{i}_all_lines"))
#   filename <- str_glue("{path}ind_all_pts_states/jag{i}_all_lines.shp")
#   writeVector(export, filename, overwrite = T)
# }
# 
# for (i in jags_reg_int){
#   export <- get(str_glue("jag{i}_expl_pts"))
#   filename <- str_glue("{path}ind_expl_pts/jag{i}_expl_pts.shp")
#   writeVector(export, filename, overwrite = T)
# }
# 
# for (i in jags_reg_int){
#   export <- get(str_glue("jag{i}_expl_lines"))
#   filename <- str_glue("{path}ind_expl_pts/jag{i}_expl_lines.shp")
#   writeVector(export, filename, overwrite = T)
# }
#
# for (i in jags_reg_int){
#   export <- get(str_glue("jag{i}_ne_pts"))
#   filename <- str_glue("{path}ind_ne_pts/jag{i}_expl_pts.shp")
#   writeVector(export, filename, overwrite = T)
# }
# 
# for (i in jags_reg_int){
#   export <- get(str_glue("jag{i}_ne_lines"))
#   filename <- str_glue("{path}ind_ne_pts/jag{i}_expl_lines.shp")
#   writeVector(export, filename, overwrite = T)
# }



## Importing

# Year not in exported shapes for some reason, adding on import
ind_y <- data %>%
  group_by(ID) %>%
  dplyr::summarise(year = round(mean(year), 0))

for (i in jags_reg_int){
  file <- vect(str_glue("ind_all_pts_states/jag{i}_all_pts.shp"))
  file$year <- as.data.frame(data)[data$ID==i, "year"]
  name <- str_glue("jag{i}_all_pts")
  assign(name, file)
}

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

# plot(SA); plot(expl_pts_sp, add = TRUE)


# Export as .csv to use for point data extraction
expl_pts_df <- as_tibble(expl_pts_sp) %>%
  # dplyr::select(Event_ID, ID, state) %>%
  mutate(ID = as.integer(str_sub(as.character(ID),4)),
         x = geom(expl_pts_sp)[,3],
         y = geom(expl_pts_sp)[,4]) %>%
  dplyr::select(Event_ID, ID, x, y, step, angle, state)
# write_csv(expl_pts_df, str_glue("{path}ind_reg_samp_expl.csv"))

ne_pts_df <- as_tibble(ne_pts_sp) %>%
  # dplyr::select(Event_ID, ID, state) %>%
  mutate(ID = as.integer(str_sub(as.character(ID),4)),
         x = geom(ne_pts_sp)[,3],
         y = geom(ne_pts_sp)[,4]) %>%
  dplyr::select(Event_ID, ID, x, y, step, angle, state)
# write_csv(ne_pts_df, str_glue("{path}ind_reg_samp_ne.csv"))



#### Background pts ####

# Buffering each pt for each individual based on mean step length
## Buffering at 6* mean hourly step length
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

str(step_m)
summary(step_m)


# Buffer equal to 6 * mean hourly step length (function for each point)
## Using individual id and st = movement state
ind_poly_pt <- function(id, st) {
  d <- get(str_glue("jag{id}_{st}_pts"))
  b <- as.numeric(step_m[step_m$ID == id, "m_step_hr"] *6000)     # to km * 6 hr
  
  pbufs <- list()
  for(i in 1:length(d)) {
    pbufs[[i]] <- buffer(d[i], b)
  }
  
  m <- do.call(rbind, pbufs)
  f <- aggregate(m)
  f$ID <- id
  return(f)
}


t86_all <- ind_poly_pt(86, "all")
t86_expl <- ind_poly_pt(86, "expl")
plot(t86_all)
plot(t86_expl, add = T)
plot(jag86_all_pts, add = T)
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

# writeVector(polys_6hr_buf_all, str_glue("{path}ind_polys_pts/polys_step_6hr_mean_buf_all.shp"), overwrite = T)
# writeVector(polys_6hr_buf_expl, str_glue("{path}ind_polys_pts/polys_step_6hr_mean_buf_expl.shp"), overwrite = T)
# writeVector(polys_6hr_buf_ne, str_glue("{path}ind_polys_pts/polys_step_6hr_mean_buf_ne.shp"), overwrite = T)

polys_6hr_buf_all <- vect("ind_polys_pts/polys_step_6hr_mean_buf_all.shp")

# Generating random (background) points within buffer = 6*mean_step_hr for each point
## 1:1 regular sample within 6 hr buffer polygons
### Using individual id, st = movement state
ind_bg_reg <- function(id, st) {
  n <- which(jags_reg_int==id)
  l <- get(str_glue("polys_6hr_buf_{st}_list"))
  b <- l[[n]]
  
  pts <- get(str_glue("jag{id}_{st}_pts"))
  y <- mean(pts$year)%>%round(0)
  
  set.seed(234)
  p <- spatSample(b, nrow(get(str_glue("jag{id}_{st}_pts"))), method = "regular")
  p$year <- y
  return(p)
}
polys_6hr_buf_all_list[[1]]
jag5_all_pts
ind_bg_reg(5, "all")
plot(jag5_all_pts); plot(ind_bg_reg(5, "all"), col = "red", add = T); plot(jag5_polys_6hr_buf_all, add = T)

# Applying to all individuals
bg_reg_6hr_buf_all_list <- lapply(jags_reg_int, function(x) ind_bg_reg(x, "all"))
bg_reg_6hr_buf_expl_list <- lapply(jags_reg_int, function(x) ind_bg_reg(x, "expl"))
bg_reg_6hr_buf_ne_list <- lapply(jags_reg_int, function(x) ind_bg_reg(x, "ne"))

bg_reg_6hr_buf_all <- do.call(rbind, bg_reg_6hr_buf_all_list)
bg_reg_6hr_buf_expl <- do.call(rbind, bg_reg_6hr_buf_expl_list)
bg_reg_6hr_buf_ne <- do.call(rbind, bg_reg_6hr_buf_ne_list)

# Creating individual shapes
for (i in 1:length(jags_reg_int)) {
  id <- jags_reg_int[i]
  assign(str_glue("jag{id}_bg_reg_6hr_buf_all"), 
         bg_reg_6hr_buf_all_list[[i]])
  assign(str_glue("jag{id}_bg_reg_6hr_buf_expl"), 
         bg_reg_6hr_buf_expl_list[[i]])
  assign(str_glue("jag{id}_bg_reg_6hr_buf_ne"), 
         bg_reg_6hr_buf_ne_list[[i]])
}


plot(jag14_polys_6hr_buf_all)
plot(jag14_all_pts, add = T)
plot(jag14_bg_reg_6hr_buf_all, col = "red", cex = 0.4, add = T)

plot(jag105_polys_6hr_buf_all)
plot(jag105_all_pts, add = T)
plot(jag105_bg_reg_6hr_buf_all, col = "red", cex = 0.4, add = T)

# writeVector(bg_reg_6hr_buf_all, str_glue("{path}ind_polys_pts/bg_reg_step_6hr_mean_buf_all.shp"), overwrite = T)
# writeVector(bg_reg_6hr_buf_expl, str_glue("{path}ind_polys_pts/bg_reg_step_6hr_mean_buf_expl.shp"), overwrite = T)
# writeVector(bg_reg_6hr_buf_ne, str_glue("{path}ind_polys_pts/bg_reg_step_6hr_mean_buf_ne.shp"), overwrite = T)


## 1:1 regular sample equally distributed throughout buffer shapes (for each individual)
# calculating area (to distribute points evenly)
# polys_6hr_buf_all_agg <- terra::aggregate(polys_6hr_buf_all) %>% disagg()
polys_6hr_buf_all$id_p <- 1:length(polys_6hr_buf_all)
polys_6hr_buf_all$area <- expanse(polys_6hr_buf_all, unit = "km")
sum(polys_6hr_buf_all$area)

### Using polys from individual 6 hr bufs (poly_i), polygon vector (p_vect), total points (t_pt)
ind_bg_reg4 <- function(id, p_vect, t_pt) {
  l <- get(p_vect)
  p <- subset(l, l$ID == id) # target polygon
  
  # polygon area proportional to total
  t_area <- sum(l$area)
  p_area <- p$area/t_area
  n <- (t_pt*p_area) %>% round(0)
  
  pts <- get(str_glue("jag{id}_all_pts"))
  y <- mean(pts$year)%>%round(0)
  
  set.seed(234)
  samp <- spatSample(p, n, method = "regular")
  samp$year <- y
  return(samp)
}
polys_6hr_buf_all_list[[1]]
jag5_all_pts
t <- ind_bg_reg4(5, "polys_6hr_buf_all", length(states_sp))
plot(jag5_all_pts); plot(t, col = "red", add = T); plot(jag5_polys_6hr_buf_all, add = T)

# Applying to all individuals
bg_reg4_6hr_buf_all_list <- lapply(jags_reg_int, function(x) ind_bg_reg4(x, "polys_6hr_buf_all", length(states_sp)))

bg_reg4_6hr_buf_all <- do.call(rbind, bg_reg4_6hr_buf_all_list)

# Creating individual shapes
for (i in 1:length(jags_reg_int)) {
  id <- jags_reg_int[i]
  assign(str_glue("jag{id}_bg_reg4_6hr_buf_all"), 
         bg_reg4_6hr_buf_all_list[[i]])
}


plot(jag14_polys_6hr_buf_all)
plot(jag14_all_pts, add = T)
plot(jag14_bg_reg4_6hr_buf_all, col = "red", cex = 0.4, add = T)

plot(jag105_polys_6hr_buf_all)
plot(jag105_all_pts, add = T)
plot(jag105_bg_reg4_6hr_buf_all, col = "red", cex = 0.4, add = T)

# writeVector(bg_reg4_6hr_buf_all, str_glue("{path}ind_polys_pts/bg_reg4_step_6hr_mean_buf_all.shp"))


## 1:1 random sample equally distributed throughout buffer shapes (for each individual)
# calculating area (to distribute points evenly)
# polys_6hr_buf_all_agg <- terra::aggregate(polys_6hr_buf_all) %>% disagg()
polys_6hr_buf_all$id_p <- 1:length(polys_6hr_buf_all)
polys_6hr_buf_all$area <- expanse(polys_6hr_buf_all, unit = "km")
sum(polys_6hr_buf_all$area)

### Using polys from individual 6 hr bufs (poly_i), polygon vector (p_vect), total points (t_pt), sampling method (m)
ind_bg_prop <- function(id, p_vect, t_pt, m) {
  l <- get(p_vect)
  p <- subset(l, l$ID == id) # target polygon
  
  # polygon area proportional to total
  t_area <- sum(l$area)
  p_area <- p$area/t_area
  n <- (t_pt*p_area) %>% round(0)
  
  pts <- get(str_glue("jag{id}_all_pts"))
  y <- mean(pts$year)%>%round(0)
  
  set.seed(234)
  # terra::spatSample showing heavy North bias; using sf:st_sample instead
  # samp <- spatSample(p, n, method = m)
  samp_sf <- sf::st_sample(st_as_sf(p), size = n, type = m)
  samp <- vect(samp_sf)
  samp$year <- y
  samp$ID <- id
  
  return(samp)
}

t <- ind_bg_prop(5, "polys_6hr_buf_all", length(states_sp), "random")
plot(jag5_all_pts); plot(t, col = "red", add = T); plot(polys_6hr_buf_all[1], add = T)

# Applying to all individuals
bg_rand_prop_6hr_buf_all_list <- lapply(jags_reg_int, function(x) ind_bg_prop(x, "polys_6hr_buf_all", length(states_sp), "random"))
bg_rand_prop_6hr_buf_all <- do.call(rbind, bg_rand_prop_6hr_buf_all_list)

# Creating individual shapes
for (i in 1:length(jags_reg_int)) {
  id <- jags_reg_int[i]
  assign(str_glue("jag{id}_bg_rand_prop_6hr_buf_all"), 
         bg_rand_prop_6hr_buf_all_list[[i]])
}


plot(subset(polys_6hr_buf_all, polys_6hr_buf_all$ID==14))
plot(jag14_all_pts, add = T)
plot(jag14_bg_rand_prop_6hr_buf_all, col = "red", cex = 0.4, add = T)

plot(subset(polys_6hr_buf_all, polys_6hr_buf_all$ID==105))
plot(jag105_all_pts, add = T)
plot(jag105_bg_rand_prop_6hr_buf_all, col = "red", cex = 0.4, add = T)

plot(subset(polys_6hr_buf_all, polys_6hr_buf_all$ID==10))
plot(jag10_all_pts, add = T)
plot(jag10_bg_rand_prop_6hr_buf_all, col = "red", cex = 0.4, add = T)

######### test #####
# noticing a North bias for random sampling, even when stratifying as below, and when removing crs...
ptest <- polys_6hr_buf_all
crs(ptest)<- NULL
t <- spatSample(polys_6hr_buf_all, 1e3, strata = "ID")
plot(polys_6hr_buf_all[1:5]); plot(t, col = "red", cex = 0.4, add = T)

# will post and ask about it

t <- st_sample(sf::st_as_sf(polys_6hr_buf_all[1:5]), 500, by_polygon = T)


writeVector(bg_rand_prop_6hr_buf_all, "ind_polys_pts/bg_rand_prop_step_6hr_mean_buf_all.shp", overwrite = T)


# 4 random points within buffer of each movement point (4:1 ratio)
### Using individual id, n = number of background pts/per pt, st = movement state
### adding ID column to link bg pts to correct temporal data
ind_bg <- function(id, n, st) {
  d <- get(str_glue("jag{id}_{st}_pts"))
  b <- as.numeric(step_m[step_m$ID == id, "m_step_hr"] *6000)     # to km * 6 hr
  y <- first(d$year)
  
  pbufs <- list()
  for(i in 1:length(d)) {
    pbufs[[i]] <- buffer(d[i], b)
  }
  
  m <- do.call(rbind, pbufs)
  set.seed(234)
  p_bg <- spatSample(m, rep(n, length(m)), method = "random")   # reg sampling not working (2 points showing up instead of 4)
  p_bg$year <- y
  return(p_bg)
}

t14_all_bg4 <- ind_bg(14, 4, "all")
t14_all_bg4
plot(jag14_all_pts)
plot(t14_all_bg4, col = "blue", cex = 0.4, add = T)

# Applying to all individuals
bg4_6hr_buf_all_list <- lapply(jags_reg_int, function(x) ind_bg(x, 4, "all"))
bg4_6hr_buf_expl_list <- lapply(jags_reg_int, function(x) ind_bg(x, 4, "expl"))
bg4_6hr_buf_ne_list <- lapply(jags_reg_int, function(x) ind_bg(x, 4, "ne"))

bg4_6hr_buf_all <- do.call(rbind, bg4_6hr_buf_all_list) 
bg4_6hr_buf_expl <- do.call(rbind, bg4_6hr_buf_expl_list) 
bg4_6hr_buf_ne <- do.call(rbind, bg4_6hr_buf_ne_list) 


# Creating individual shapes
for (i in 1:length(jags_reg_int)) {
  id <- jags_reg_int[i]
  assign(str_glue("jag{id}_bg4_6hr_buf_all"), 
         bg4_6hr_buf_all_list[[i]])
  assign(str_glue("jag{id}_bg4_6hr_buf_expl"), 
         bg4_6hr_buf_expl_list[[i]])
  assign(str_glue("jag{id}_bg4_6hr_buf_ne"), 
         bg4_6hr_buf_ne_list[[i]])
}

plot(jag65_polys_6hr_buf_expl)
plot(jag65_expl_pts, add = T)
plot(jag65_bg4_6hr_buf_expl, col = "blue", cex = 0.3, add = T)
plot(jag65_bg_reg_6hr_buf_expl, col = "red", cex = 0.5, add = T)

# writeVector(bg4_6hr_buf_all, str_glue("{path}ind_polys_pts/bg4_step_6hr_mean_buf_all.shp"), overwrite=T)
# writeVector(bg4_6hr_buf_expl, str_glue("{path}ind_polys_pts/bg4_step_6hr_mean_buf_expl.shp"), overwrite=T)
# writeVector(bg4_6hr_buf_ne, str_glue("{path}ind_polys_pts/bg4_step_6hr_mean_buf_ne.shp"), overwrite=T)

# ##### All pts ####
# ### Polys
# polys_6hr_buf_all_list <- lapply(jags_reg_int, function(x) ind_poly_pt(x, "all"))
# 
# plot(polys_6hr_buf_all_list[[21]], main = polys_6hr_buf_all_list[[21]]$ID)
# plot(jag65_all_pts, add = T)
# 
# polys_6hr_buf_all <- do.call(rbind, polys_6hr_buf_all_list)
# 
# 
# 
# ### BG points
# set.seed(234)
# 
# bg4_6hr_buf_all_list <- lapply(jags_reg_int, function(x) ind_bg(x, 4, "all"))
# 
# bg4_6hr_buf_all <- do.call(rbind, bg4_6hr_buf_all_list) 
# 
# plot(bg4_6hr_buf_all, col = "red", cex = 0.4, add = T)
# plot(jag7_all_pts, add = T)
# 
# 
# # Pulling individuals' data
# for (i in 1:length(jags_reg_int)) {
#   id <- jags_reg_int[i]
#   assign(str_glue("jag{id}_polys_6hr_buf_all"), 
#          polys_6hr_buf_all_list[[i]])
#   assign(str_glue("jag{id}_bg4_6hr_buf_all"), 
#          bg4_6hr_buf_all_list[[i]])
# }
# 
# # testing
# plot(jag89_polys_6hr_buf_all)
# plot(jag89_all_pts, add = T)
# plot(jag89_bg4_6hr_buf, cex = 0.3, col = "red", add = T)
# 
# # Exporting
# # writeVector(polys_6hr_buf_all, str_glue("{path}ind_polys_pts/polys_step_6hr_mean_buf_all.shp"), overwrite = T)
# # writeVector(bg4_6hr_buf_all, str_glue("{path}ind_polys_pts/bg4_step_6hr_mean_buf_all.shp"))
# 
# 
# ##### Expl pts ####
# ### Polys
# polys_6hr_buf_expl_list <- lapply(jags_reg_int, function(x) ind_poly_pt(x, "expl"))
# 
# plot(polys_6hr_buf_expl_list[[21]], main = polys_6hr_buf_expl_list[[21]]$ID)
# plot(jag65_expl_pts, add = T)
# 
# polys_6hr_buf_expl <- do.call(rbind, polys_6hr_buf_expl_list)
# 
# 
# ### BG points
# set.seed(234)
# 
# bg4_6hr_buf_expl_list <- lapply(jags_reg_int, function(x) ind_bg(x, 4, "expl"))
# 
# bg4_6hr_buf_expl <- do.call(rbind, bg4_6hr_buf_expl_list) 
# # plot(jag20_expl_pts)
# # plot(bg4_pts_st_m_buf_expl, col = "red", cex = 0.4, add = T)
# 
# # writeVector(polys_6hr_buf_expl, str_glue("{path}ind_polys_pts/polys_step_6hr_mean_buf_expl.shp"), overwrite = T)
# # writeVector(bg4_6hr_buf_expl, str_glue("{path}ind_polys_pts/bg4_step_6hr_mean_buf_expl.shp"))
# 
# 
# # Pulling individuals' data
# for (i in 1:length(jags_reg_int)) {
#   id <- jags_reg_int[i]
#   assign(str_glue("jag{id}_polys_6hr_buf_expl"),
#          polys_6hr_buf_expl_list[[i]])
#   assign(str_glue("jag{id}_bg4_6hr_buf_expl"),
#          bg4_6hr_buf_expl_list[[i]])
# }
# 
# 
# ##### Non-expl points ####
# ### Polys
# polys_6hr_buf_ne_list <- lapply(jags_reg_int, function(x) ind_poly_pt(x, "ne"))
# 
# plot(polys_6hr_buf_ne_list[[21]], main = polys_6hr_buf_ne_list[[21]]$ID)
# plot(jag65_ne_pts, add = T)
# 
# polys_6hr_buf_ne <- do.call(rbind, polys_6hr_buf_ne_list)
# 
# ### BG points
# set.seed(234)
# 
# bg4_6hr_buf_ne_list <- lapply(jags_reg_int, function(x) ind_bg(x, 4, "ne"))
# 
# bg4_6hr_buf_ne <- do.call(rbind, bg4_6hr_buf_ne_list) 
# plot(jag20_all_pts)
# plot(jag20_ne_pts, col = "blue", add = T)
# plot(bg4_pts_st_m_buf_ne, col = "red", cex = 0.4, add = T)
# 
# # Pulling individuals' data
# for (i in 1:length(jags_reg_int)) {
#   id <- jags_reg_int[i]
#   assign(str_glue("jag{id}_polys_6hr_buf_ne"), 
#          polys_6hr_buf_ne_list[[i]])
#   assign(str_glue("jag{id}_bg4_6hr_buf_ne"), 
#          bg4_6hr_buf_ne_list[[i]])
# }
# 
# # writeVector(bg4_6hr_buf_ne, str_glue("{path}ind_polys_pts/bg4_step_6hr_mean_buf_ne.shp"))
