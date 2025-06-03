# Jaguar Machine Learning Movement Analysis
## 07 - Temporally aligned covariate extraction


### Extracting data to match year of observed data
### Presence and background points matched with temporally dynamic covariates

library(terra)
library(tidyverse)
library(ggthemes)

# Set to directory with cloned repo
setwd("~/Columbia E3B/PhD/Dissertation/Ch3-jaguar_movement/scripts/git_push_R/")

# Make this path where you stored environmental covariate data
path <- "C:/Users/jayms/Documents/Columbia E3B/PhD/Dissertation/Ch3-jaguar_movement/covariate_data/"

long_lat <- "epsg:4326"

#### Importing movement data ####
data_all <- read_csv("~/Columbia E3B/PhD/Dissertation/Ch3-jaguar_movement/ind_reg_samp_st.csv") %>%
  mutate(ID = as.character(ID),
         state = as.character(s)) %>%
  dplyr::select(-s)

data_expl <- read_csv("~/Columbia E3B/PhD/Dissertation/Ch3-jaguar_movement/ind_reg_samp_expl.csv") %>%
  mutate(ID = as.character(ID),
         state = as.character(state)) %>%
  left_join(data_all)

data_ne <- read_csv("~/Columbia E3B/PhD/Dissertation/Ch3-jaguar_movement/ind_reg_samp_ne.csv") %>%
  mutate(ID = as.character(ID),
         state = as.character(state)) %>%
  left_join(data_all)


data_all_sp <- vect(data_all, geom = c("x", "y"), crs = long_lat)
data_expl_sp <- vect(data_expl, geom = c("x", "y"), crs = long_lat)
data_ne_sp <- vect(data_ne, geom = c("x", "y"), crs = long_lat)


#### Temporal alignment ####

# Function to extract temporally aligned data from stack of annual rasters
t_align <- function(covar, pts){
  e <- terra::extract(covar, pts) %>%
    `colnames<-`(c("ID", 2000:2016)) %>%
    mutate(year = pts$year,
           var = NA)
  
  for(i in 1:nrow(e)){
    for(j in 2:length(e)){
      if(e$year[i] == colnames(e)[j]) {
        e$var[i] = e[i,j]
      }
    }
  }
  return(e$var %>%
           as.data.frame() %>%
           `names<-`(substr(names(covar[[1]]), 
                            1, nchar(names(covar[[1]])) -10)) # -10 for _p_sc covariates
  )
}


#### Extractions ####

##### Presence points ####

# Importing layers

## Polygon scaled rasters
stack_covar_p_sc <- list.files(str_glue("{path}cov_stack_polygon_scaled/"), 
                               pattern = ".tif",
                               all.files = TRUE,
                               full.names = TRUE) %>%
  rast()

covars <- c("dist_roads", "dens_roads", "elev", "slope", "pop_dens", "hii", "gpp", "pct_tc",
            "dist_tc_p20", "dist_tc_p40", "dist_non_tc_p20", "dist_non_tc_p40", 
            "dist_nat", "dist_non_nat", "dist_wat_perm", "dist_wat_seas", "land_cover")
for(i in covars[1:16]) {
  assign(str_glue("{i}_p_sc"), c(stack_covar_p_sc[i]))
}

# Land cover separate (not normalized by polygon)
lc <- list.files(str_glue("{path}landcover/"),
                 pattern = "lc_2",
                 all.files = TRUE, 
                 full.names = TRUE) %>%
  rast() %>%
  `names<-`(str_glue("lc_{2000:2016}"))

# Binding temporally aligned covariate extractions to others
cov1_list <- list(pop_dens_p_sc, hii_p_sc,gpp_p_sc, pct_tc_p_sc, dist_tc_p20_p_sc,
                  dist_tc_p40_p_sc, dist_non_tc_p20_p_sc, dist_non_tc_p40_p_sc,
                  dist_nat_p_sc, dist_non_nat_p_sc, dist_wat_perm_p_sc, dist_wat_seas_p_sc)

covar1 <- lapply(cov1_list, function(x) (t_align(x, data_all_sp))) %>%
  bind_cols()

covar2 <- terra::extract(c(dist_roads_p_sc, dens_roads_p_sc, elev_p_sc, slope_p_sc), 
                         data_all_sp) %>%
  rename(u_id = ID)

# Extracting land cover at higher resolution for accuracy
covar_lc <- t_align(lc, data_all_sp) %>%
  `names<-`("land_cover") %>%
  mutate(land_cover = as.factor(land_cover))

# Binding all pres data together ("p_sc" = polygon scaled)
pres_data <- bind_cols(data_all, covar1, covar2, covar_lc) %>%
  dplyr::select(23, 1:length(.)) %>%
  mutate(pres = 1) 
names(pres_data)[c(12:23)] <- str_glue("{names(pres_data)[c(12:23)]}_p_sc")

# Exporting
# write_csv(pres_data, str_glue("{path}pres_data_covar.csv"))

# Importing
# pres_data <- read_csv(str_glue("{path}pres_data_covar.csv")) %>%
#   mutate(land_cover = as.factor(land_cover),
#          ID = as.factor(ID))

##### Background points ####
## Extracting

# Importing bg points created from 6*mean_hrly_step_buffer for each individual's pts
## 1:1 bg pts proportional to each individual's 6 hr poly buffer

bg_prop_all <- vect("~/Columbia E3B/PhD/Dissertation/Ch3-jaguar_movement/ind_polys_pts/bg_rand_prop_step_6hr_mean_buf_all.shp")

bg_all_covar1 <- lapply(cov1_list, function(x) (t_align(x, bg_prop_all))) %>%
  bind_cols()
bg_all_covar2 <- terra::extract(c(dist_roads_p_sc, dens_roads_p_sc, elev_p_sc, slope_p_sc),
                                bg_prop_all) %>%
  rename(u_id = ID)

# Adding landcover
bg_all_lc <- t_align(lc, bg_prop_all) %>%
  `names<-`("land_cover") %>%
  mutate(land_cover = as.factor(land_cover))

# Binding all bg data together ("p_sc" = polygon scaled)
bg_all_data <- bind_cols(bg_all_covar1, bg_all_covar2, bg_all_lc) %>%
  mutate(
    x = geom(bg_prop_all)[,3],
    y = geom(bg_prop_all)[,4],
    year = bg_prop_all$year,
    pres = 0) %>%
  as_tibble() %>%
  dplyr::select(13, 19:20, 1:length(.))

names(bg_all_data)[c(4:15)] <- str_glue("{names(bg_all_data)[c(4:15)]}_p_sc")
names(bg_all_data)

# Exporting
# write_csv(bg_all_data, "covariate_data/bg_rand_prop_cov_data.csv")

# Importing
# bg_all_data <- read_csv(str_glue("{path}/bg_rand_prop_cov_data.csv")) %>%
#   mutate(land_cover = as.factor(land_cover))


##### Full Data ####

# Adding bg to pres for full data set (presence + background) and extracted covariate values
## All three subsets ("all", "expl", and "non-expl") will be modeled with same background data

cov_data_all <- bind_rows(pres_data[,-c(2,5:11)],
                          bg_all_data[,-21]) %>%
  mutate(u_id = 1:nrow(.))

cov_data_expl <- pres_data %>%
  filter(state == 3) %>%
  dplyr::select(-c(2,5:11)) %>%
  bind_rows(bg_all_data[,-21])

cov_data_ne <- pres_data %>%
  filter(state != 3) %>%
  dplyr::select(-c(2, 5:11)) %>%
  bind_rows(bg_all_data[,-21])

###### Adding polygon ID ####
# Adding polygon column to data

polys_10k <- vect("~/Columbia E3B/PhD/Dissertation/Ch3-jaguar_movement/jags_reg_int_ind_polygons_10k_buf.shp") %>%
  disagg() %>%
  `values<-`(NULL)
polys_10k$poly_ID <- 1:length(polys_10k)

for(i in 1:length(polys_10k)) {
  assign(str_glue("poly{i}"), subset(polys_10k, polys_10k$poly_ID==i))
}


cov_data_all_sp <- vect(cov_data_all, geom = c("x","y"), crs = long_lat)
cov_data_expl_sp <- vect(cov_data_expl, geom = c("x","y"), crs = long_lat)
cov_data_ne_sp <- vect(cov_data_ne, geom = c("x","y"), crs = long_lat)

## All data
poly_extract_all <- terra::extract(polys_10k, cov_data_all_sp) %>%
  mutate(x = geom(cov_data_all_sp)[,3],
         y = geom(cov_data_all_sp)[,4]) %>%#,
  dplyr::select(-id.y)
head(poly_extract_all)

## Expl data
poly_extract_expl <- terra::extract(polys_10k, cov_data_expl_sp) %>%
  mutate(x = geom(cov_data_expl_sp)[,3],
         y = geom(cov_data_expl_sp)[,4]) %>%#,
  dplyr::select(-id.y)
head(poly_extract_expl)

## Non-expl data
poly_extract_ne <- terra::extract(polys_10k, cov_data_ne_sp) %>%
  mutate(x = geom(cov_data_ne_sp)[,3],
         y = geom(cov_data_ne_sp)[,4]) %>%#,
  dplyr::select(-id.y)
head(poly_extract_ne)

cov_data_all <- left_join(cov_data_all, poly_extract_all)
cov_data_expl <- left_join(cov_data_expl, poly_extract_expl)
cov_data_ne <- left_join(cov_data_ne, poly_extract_ne)

names(cov_data_ne)

##### Quality check ####
# Check for NAs
summary(cov_data_all)
summary(cov_data_expl)
summary(cov_data_ne)

# NAs are all in dens_roads and hii
# Locating NAs

ids <- cov_data_all %>% 
  filter(is.na(dens_roads_p_sc)) %>%
  # filter(if_any(everything(), is.na)) %>%
  dplyr::select(u_id) %>%
  unique()
cov_data_all %>% filter(u_id %in% ids$u_id) %>% dplyr::select(poly_ID) %>% unique()

### All in polygon 10 for dens_roads

# Since polygon 10 (Pantanal) has 0 road density, scaling turned it to NA
## Setting road_dens here to 0

### HII is NA in water as well; setting NAs in HII to zero

cov_data_all <- mutate(cov_data_all, 
                       dens_roads_p_sc = replace_na(dens_roads_p_sc, 0),
                       hii_p_sc = replace_na(hii_p_sc, 0))
cov_data_expl <- mutate(cov_data_expl, 
                        dens_roads_p_sc = replace_na(dens_roads_p_sc, 0),
                        hii_p_sc = replace_na(hii_p_sc, 0))
cov_data_ne <- mutate(cov_data_ne, 
                      dens_roads_p_sc = replace_na(dens_roads_p_sc, 0),
                      hii_p_sc = replace_na(hii_p_sc, 0))
summary(cov_data_all)
summary(cov_data_expl)
summary(cov_data_ne)

# No NAs


# # Writing to disk for next steps of analysis
# write_csv(cov_data_all, str_glue("{path}cov_data_all_bg_rand_prop_ind_pt_6hr_buf.csv"))
# write_csv(cov_data_expl, str_glue("{path}cov_data_expl_bg_rand_prop_ind_pt_6hr_buf.csv"))
# write_csv(cov_data_ne, str_glue("{path}cov_data_ne_bg_rand_prop_ind_pt_6hr_buf.csv"))
