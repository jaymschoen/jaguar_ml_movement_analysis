# Jaguar Machine Learning Movement Analysis
## 06 - Covariate scaling

### Scaling covariates by polygon for global model
### Assumption: jaguars will respond to RELATIVE effect of environmental variables (see manuscript for discussion)

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

## List of regularly sampled jags used for analysis
# jags_reg_int <- c(5,6,7,9,10,12,13,14,15,16,18,19,20,22,23,25,
#                   41,42,50,52,65,68,69,73,74,75,76,77,79,80,81,
#                   84,86,87,88,89,96,99,105,108,114,115,116,117)

data_all_sp <- vect(data_all, geom = c("x", "y"), crs = long_lat)
data_expl_sp <- vect(data_expl, geom = c("x", "y"), crs = long_lat)
data_ne_sp <- vect(data_ne, geom = c("x", "y"), crs = long_lat)

polys_10k <- vect("~/Columbia E3B/PhD/Dissertation/Ch3-jaguar_movement/jags_reg_int_ind_polygons_10k_buf.shp") %>%
  disagg() 

SA <- vect("~/Columbia E3B/Paraguay Jaguar Project/Paraguay Jaguar/Paraguay Shape Files/Admin areas/SouthAmerica.shp")

# Visualize
plot(polys_10k, col = "orange4"); plot(SA, add = T)
plot(sample(data_all_sp, 1e3), col = 'red', cex = 0.3, add = T)


# Load covariate data from 05-covariate_import script
source("05-covariate_import.R")


#### Covariate visualizations ####
# Visualizing covariate distributions
covars <- c("dist_roads", "dens_roads", "elev", "slope", "pop_dens", "hii", "gpp", "pct_tc",
                 "dist_tc_p20", "dist_tc_p40", "dist_non_tc_p20", "dist_non_tc_p40", 
                 "dist_nat", "dist_non_nat", "dist_wat_perm", "dist_wat_seas")

covars_full <- c("Distance to roads", "Road density", "Elevation", "Slope", 
                 "Population density", "Human impact index", "Gross primary production",
                 "Percent tree cover", "Distance to medium tree cover", 
                 "Distance to high tree cover", "Distance to very low tree cover",
                 "Distance to low tree cover", "Distance to natural areas", 
                 "Distance to matrix", "Distance to permanent water", 
                 "Distance to seasonal water")

colors16 <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "gold1", "skyblue2", 
              "palegreen2", "#FDBF6F", "gray70", "maroon", "orchid1", "darkturquoise", "darkorange4", "brown")

colors13 <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "gold1", "skyblue2", 
              "#FDBF6F", "gray70", "maroon", "darkturquoise", "steelblue4", "brown")

##### Density plots ####

## 2014 data as example

## Combining static and dynamic covariate data (for 2014)
c1 <- c(dist_roads, dens_roads, elev, slope)
c2 <- rast(lapply(covars[5:length(covars)], function(x) get(x)[[15]]))

covar_2014 <- c(c1, c2)

# Visualizing distributions of covariate data between polygons
poly_list <- lapply(1:length(polys_10k), function(x) polys_10k[x]) %>%
  `names<-`(1:length(.))
for(j in 1:length(poly_list)) {poly_list[[j]]$poly = j}

dens_list <- list()
for(i in 1:nlyr(covar_2014)) {
  dens_list[[i]] <- do.call(bind_rows, lapply(poly_list, function(p) {
    v = covars_full[i]
    values = covar_2014[[i]] %>% crop(p) %>% values(mat = F, na.rm = T)
    poly = p$poly %>% first()
    
    return(tibble(variable = v,
                  values = values,
                  polygon = poly)
    )
  }))
}

dens_data <- do.call(bind_rows, dens_list)

dens_plot <- ggplot(dens_data, aes(x = values, color = as.factor(polygon))) +
  stat_density(geom = "line", linewidth = 0.5, trim = F, position = "identity") +
  facet_wrap(~variable, scales = "free") +
  scale_color_manual(values = colors13) +
  theme_bw() +
  xlab("Value") +
  ylab("Density") +
  labs(color = "Polygon (region)") +
  theme(strip.text = element_text(face = "bold", color = "white"),
        strip.background = element_rect(fill = "black"),
        legend.key.width = unit(1, "cm")) +
  guides(color = guide_legend(override.aes = list(linewidth = 1)))

plot(dens_plot)

# Large variation in densities for several covariates
## Will normalize data by polygon (region) for relative selection for each polygon (region)


#### Scaling rasters by polygon ####
polys_10k <- polys_10k %>%
  disagg() %>%
  `values<-`(NULL)
polys_10k$poly_ID <- 1:length(polys_10k)

for(i in 1:length(polys_10k)) {
  assign(str_glue("poly{i}"), subset(polys_10k, polys_10k$poly_ID==i))
}

covars_p_sc <- c(str_glue("{covars}_p_sc"))

# Function to scale each covariate within the polygon
poly_scale <- function(poly, covar) {
  p <- get(str_glue("poly{poly}"))
  scaled <- get(covar) %>%
    crop(p, mask = T) %>%
    scale(center = FALSE) %>%
    `names<-`(str_glue("{names(.)}_p{poly}_sc"))
  return(scaled)
}

# Applying to each covariate within each polygon (region)

terraOptions(steps = 5)  # helping with memory allocation (breaking into more steps)

list_p_sc <- list()
for(i in covars) {
  l1 <- lapply(polys_10k$poly_ID, function(x) poly_scale(x, i))
  list_p_sc[[i]] <- lapply(1:nlyr(l1[[1]]), function(y) do.call(merge, sapply(l1, "[[", y))) %>%
    rast() 
  if(nlyr(list_p_sc[[i]]) >1) {
    names(list_p_sc[[i]]) <- str_glue("{i}_{2000:2016}_p_sc")
  }
  else{names(list_p_sc[[i]]) <- str_glue("{i}_p_sc")}
}

for(i in covars) {
  assign(str_glue("{i}_p_sc"), list_p_sc[[i]])
}

# Writing polygon scaled raster stacks to disk

# for(i in covars_p_sc) {
#   c <- get(i)
#   for(j in 1:nlyr(c)) {
#     export <- c[[j]]
#     filename <- str_glue("{path}cov_stack_polygon_scaled/{names(c[[j]])}.tif")
#     writeRaster(export, filename)
#   }
# }