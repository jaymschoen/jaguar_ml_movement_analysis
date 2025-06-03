# Jaguar Machine Learning Movement Analysis
## 05 - Covariate imports

### Script to import covariate data (loaded by 06 - Covariate scaling script)


library(terra)
library(tidyverse)

# Make this path where you stored environmental covariate data
path <- "C:/Users/jayms/Documents/Columbia E3B/PhD/Dissertation/Ch3-jaguar_movement/covariate_data/"

#### Importing Covariate Data ####

##### Roads ####
###### Distance to major roads ####
dist_roads <- rast(str_glue("{path}roads/for_analysis/dist_roads_mjr.tif")) %>%
  `names<-`("dist_roads_mjr")

# plot(dist_roads)

###### Road Density (all roads) ####
dens_roads <- rast(str_glue("{path}roads/for_analysis/dens_roads_min.tif")) %>%
  `names<-`("road_density")

##### Elevation/Slope ####
elev <- rast(str_glue("{path}elevation_slope/elevation.tif"))

slope <- rast(str_glue("{path}elevation_slope/slope.tif"))

##### Population Density ####
pop_dens <- list.files(str_glue("{path}population/"),
                       pattern = ".tif",
                       all.files = TRUE, 
                       full.names = TRUE) %>%
  rast() %>%
  `names<-`(str_glue("pop_dens_{2000:2016}"))

##### Human Impact Index (WCS Human Footprint) ####
hii <- list.files(str_glue("{path}human_footprint/"),
                  pattern = ".tif",
                  all.files = TRUE, 
                  full.names = TRUE) %>%
  rast() %>%
  `names<-`(str_glue("hii_{2001:2016}"))

# Data start in 2001; copying 2001 layer for year 2000
hii_2000 <- hii[[1]] %>%
  `names<-`(str_glue("hii_2000"))

hii <- c(hii_2000, hii)

##### Gross Primary Production ####
gpp <- list.files(str_glue("{path}gpp/"),
                  pattern = ".tif$",
                  all.files = TRUE, 
                  full.names = TRUE) %>%
  rast() %>%
  `names<-`(str_glue("gpp_{2000:2016}"))

##### Tree Cover ####
###### Percent Tree Cover ####
pct_tc <- list.files(str_glue("{path}percent_tree_cover/"), 
                     pattern = ".tif$",
                     all.files = TRUE, 
                     full.names = TRUE) %>%
  rast() %>%
  `names<-`(str_glue("pct_tc_{2000:2016}"))

###### Distance to Tree Cover (distance to forest) ####
# >20th %ile tree cover
dist_tc_p20 <- list.files(str_glue("{path}distance_tree_cover_p20/"),
                          pattern = ".tif$",
                          all.files = TRUE, 
                          full.names = TRUE) %>%
  rast() %>%
  `names<-`(str_glue("dist_tc_p20_{2000:2016}"))

# >40th %ile tree cover
dist_tc_p40 <- list.files(str_glue("{path}distance_tree_cover_p40/"),
                          pattern = ".tif$",
                          all.files = TRUE, 
                          full.names = TRUE) %>%
  rast() %>%
  `names<-`(str_glue("dist_tc_p40_{2000:2016}"))

###### Distance to Non-Tree Cover (distance to forest edge) ####
# <20th %ile tree cover
dist_non_tc_p20 <- list.files(str_glue("{path}distance_non_tree_cover_p20/"),
                              pattern = ".tif$",
                              all.files = TRUE, 
                              full.names = TRUE) %>%
  rast()  %>%
  `names<-`(str_glue("dist_non_tc_p20_{2000:2016}"))

# <40th %ile tree cover
dist_non_tc_p40 <- list.files(str_glue("{path}distance_non_tree_cover_p40/"),
                              pattern = ".tif$",
                              all.files = TRUE, 
                              full.names = TRUE) %>%
  rast()  %>%
  `names<-`(str_glue("dist_non_tc_p40_{2000:2016}"))

##### MapBiomas Land Cover ####
# Keeping land cover separate to keep finer resolution (100m) for extraction. Will manually append data later.

lc <- list.files(str_glue("{path}landcover/"),
                 pattern = "lc_2",
                 all.files = TRUE, 
                 full.names = TRUE) %>%
  rast() %>%
  `names<-`(str_glue("lc_{2000:2016}"))

###### Distance to natural areas ####
dist_nat <- list.files(str_glue("{path}distance_natural/"),
                       pattern = ".tif$",
                       all.files = TRUE, 
                       full.names = TRUE) %>%
  rast() #%>%

###### Distance to non-natural areas (distance to matrix) ####
dist_non_nat <- list.files(str_glue("{path}distance_non_natural/"),
                           pattern = ".tif$",
                           all.files = TRUE, 
                           full.names = TRUE) %>%
  rast() #%>%      


##### JRC Water ####
# - Layers created from JRC annual water product in GEE library
# - Seasonal water = permanent + seasonal (larger extent)

##### Distance to permanent water ####
dist_wat_perm <- list.files(str_glue("{path}distance_water_permanent/"),
                            pattern = ".tif$",
                            all.files = TRUE,
                            full.names = TRUE) %>%
  rast() %>%
  `names<-`(str_glue("dist_wat_perm_{2000:2016}"))


##### Distance to seasonal water ####
dist_wat_seas <- list.files(str_glue("{path}distance_water_seasonal/"),
                            pattern = ".tif$",
                            all.files = TRUE,
                            full.names = TRUE) %>%
  rast() %>%
  `names<-`(str_glue("dist_wat_seas_{2000:2016}"))

#### Stacking covariates ####
## Keeping land cover separate for now (different resolution)
stack_covar <- c(dist_roads, dens_roads, elev, slope, pop_dens, hii, gpp, pct_tc, 
                 dist_tc_p20, dist_tc_p40, dist_non_tc_p20, dist_non_tc_p40,
                 dist_nat, dist_non_nat, dist_wat_perm, dist_wat_seas)

