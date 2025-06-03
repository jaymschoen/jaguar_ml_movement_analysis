# Jaguar Machine Learning Movement Analysis
## 08 - Data thinning

### Thinning data from short-interval jaguars to reduce oversampling of these individuals/areas


#### Setup ####
library(terra)
library(tidyverse)

# Set to directory with cloned repo
setwd("~/Columbia E3B/PhD/Dissertation/Ch3-jaguar_movement/scripts/git_push_R/")

# Make this path where you stored environmental covariate data
path <- "C:/Users/jayms/Documents/Columbia E3B/PhD/Dissertation/Ch3-jaguar_movement/covariate_data/"

#### Importing data ####
# Using polygon-standardized values with 1:1 random bg points (distributed proportionally within individual bg area) for each individual within 6x mean hrly step buffer
data_all <- read_csv(str_glue("{path}/cov_data_all_bg_rand_prop_ind_pt_6hr_buf.csv")) %>%
  mutate(land_cover = as.factor(land_cover),
         pres = as.factor(pres)) %>%
  unique()                                    # removing erroneous repeats

# Using previously exported presences to re-add EventID and state data
pres_st <- read_csv("../../ind_reg_samp_st.csv")

data_all_pres <- data_all %>%
  filter(pres==1) %>%
  mutate(Event_ID = pres_st$Event_ID) %>%
  dplyr::select(Event_ID, 2:length(.))

covars <- names(data_all)[4:20]

covars_short <- c(str_sub(covars[1:16], 1,-6), covars[17])

#### Sampling within data (thinning) ####

## To reduce oversampling of short interval individuals, sampling <2 hr interval individuals at 1 obs per 2 hours

# Adding individual ID's to match with EventIDs
ev_ids <- read_csv("../../event_id_ind_id.csv")

data_all_pres <- data_all_pres %>%
  left_join(ev_ids)

data_all_abs <- data_all %>%
  filter(pres==0)


### Summary stats on interval data

# Loading interval data from source_setup.R script
## source_setup.R is a shortened 01-data_prep script used to import interval data

source("source_setup.R")

## Resetting WD
setwd("C:/Users/jayms/Documents/Columbia E3B/PhD/Dissertation/Ch3-jaguar_movement/")

ind_reg_int

# Adjusting one individual with <1 hr sampling (ind42 was 0.5 hrs but interval function rounded to 0)
ind_reg_int[ind_reg_int$ID==42,"mod_int"] <- 0.5

# ind_state <- read_csv("ind_reg_samp_st.csv") %>%
ind_state <- pres_st %>%
  mutate(ID = as.character(ID)) %>%
  left_join(ind_reg_int %>% dplyr::select(ID, mod_int),
            by = "ID")

summary(ind_state$mod_int)
unique(ind_state$mod_int)

ind_state %>%
  group_by(mod_int) %>%
  summarise(prop = n()/nrow(.))

##### Presence points ####

# Function to sample individuals based on mod_int
## Thinning high resolution animals to limit model bias from their higher sampling

## based on id = individual ID, and t = thinning parameter (baseline interval)

smp_thin <- function(data, id, t) {
  d <- data %>% filter(ID==id)
  int <- ind_reg_int %>% filter(ID==id) %>% dplyr::select(mod_int) %>% as.numeric()
  
  p <- round((t/int),0)    # thinning parameter (1 in every "p" steps will be kept from data)
  
  s <- tibble()
  if(p>1){
    s <- d[seq(1, nrow(d), p),]} else{s<-d}
  # regularly sampling data based on thinning parameter (per id)
  # if individual's mod_int >t, keeping all points the same
  
  return(s)
}


# Checking for lowest interval individual
smp_thin(data_all_pres, 42, 2)
data_all_pres%>%filter(ID==42)
ind_state%>%filter(ID==42)%>%dplyr::select(mod_int)

# Re-sampling to 2 hrs
## result is 30/43 ind sampled at 2 hrs, 6 at 3 hrs, 5 at 4 hrs, 1 at 5/6 hrs each
data_all_pres_smp_2hr_l <- lapply(jags_reg_int, 
                                  function(x) smp_thin(data_all_pres, x, 2)) %>%
  `names<-`(jags_reg_int)
data_all_pres_smp_2hr <- do.call(rbind, data_all_pres_smp_2hr_l) %>%
  mutate(ID = as.character(ID))

for(i in jags_reg_int) {
  assign(str_glue("jag{i}_smp_2hr_thin_all"), data_all_pres_smp_2hr_l[[i]])
}

##### Background points ####
# For BG_rand_prop
## adding ID to data_all_abs
abs_v <- vect("ind_polys_pts/bg_rand_prop_step_6hr_mean_buf_all.shp")

abs_df <- as_tibble(abs_v) %>%
  mutate(u_id = data_all_abs$u_id) %>%
  dplyr::select(u_id, 1:length(.)) %>%
  right_join(data_all_abs, by = "u_id") %>%
  dplyr::select(1, 4:length(.), 2:3)

# thinning parameter =2 hrs for each individual
th_p <- data_all_pres %>%
  mutate(ID = as.character(ID)) %>%
  group_by(ID) %>%
  summarise(n_pts= n()) %>%
  left_join(ind_reg_int[,c("ID", "mod_int")]) %>%
  mutate(p = round((2/mod_int),0))   # thinning parameter (1 in every "p" steps will be kept from data)

# data_all_abs <- data_all %>%
#   filter(pres==0)


# Checking:
identical(abs_df[,-c(23:24)], data_all_abs)


# Check function for background data
smp_thin(abs_df, 13, 2)
abs_df%>%filter(ID==13)
ind_reg_int%>%filter(ID==13)


data_all_abs_smp_2hr_l <- lapply(jags_reg_int, function(x) smp_thin(abs_df, x, 2)) %>%
  `names<-`(jags_reg_int)
data_all_abs_smp_2hr <- do.call(rbind, data_all_abs_smp_2hr_l)

# writeVector(vect(data_all_abs_smp_2hr, geom = c("x", "y")),
#             "ind_polys_pts/bg_rand_prop_2hr_thin_6hr_buf_all.shp")

data_all_smp_2hr <- bind_rows(data_all_pres_smp_2hr[,-1], data_all_abs_smp_2hr[,-1]) %>%
  mutate(u_id = 1:nrow(.)) %>%
  dplyr::select(u_id, 1:length(.))


#### Adding bg data to movement state subsets ####

## Will model "all" (states 1-3), "exploratory" (state 3), and "non-exploratory" (states 1-2)

# Exploratory data 
data_expl_smp_2hr <- data_all_pres %>%
  left_join(ind_state %>% dplyr::select(Event_ID, s, year)) %>%
  filter(Event_ID %in% data_all_pres_smp_2hr$Event_ID,
         s==3) %>%
  dplyr::select(-c(Event_ID, s)) %>%
  mutate(ID = as.character(ID)) %>%
  bind_rows(data_all_abs_smp_2hr[,-1])


# Non-expl
data_ne_smp_2hr <- data_all_pres %>%
  left_join(ind_state %>% dplyr::select(Event_ID, s, year)) %>%
  filter(Event_ID %in% data_all_pres_smp_2hr$Event_ID,
         s!=3) %>%
  dplyr::select(-c(Event_ID, s)) %>%
  mutate(ID = as.character(ID)) %>%
  bind_rows(data_all_abs_smp_2hr[,-1])

# # Exporting
# write_csv(data_all_smp_2hr,
#           "covariate_data/cov_data_all_bg_rand_prop_2hr_thin_ind_pt_6hr_buf.csv")
# write_csv(data_expl_smp_2hr,
#           "covariate_data/cov_data_expl_bg_rand_prop_2hr_thin_ind_pt_6hr_buf.csv")
# write_csv(data_ne_smp_2hr,
#           "covariate_data/cov_data_ne_bg_rand_prop_2hr_thin_ind_pt_6hr_buf.csv")