# Jaguar Machine Learning Movement Analysis
## 12 - Predictions

### Model predictions
### Generating rasterized predictions of 3 models for visualization and interpretation


#### Setup ####
library(terra)
library(ranger)
library(tidyverse)
library(tidymodels)
# library(pdp)
# library(plotly)


# Set to directory with cloned repo
setwd("~/Columbia E3B/PhD/Dissertation/Ch3-jaguar_movement/scripts/git_push_R/")

# Make this path where you stored environmental covariate data
path <- "C:/Users/jayms/Documents/Columbia E3B/PhD/Dissertation/Ch3-jaguar_movement/covariate_data/"

#### Imports ####
# Importing 2-hr interval thinned data

data_all_smp_2hr <- read_csv(str_glue("{path}/cov_data_all_bg_rand_prop_2hr_thin_ind_pt_6hr_buf.csv"))
data_expl_smp_2hr <- read_csv(str_glue("{path}/cov_data_expl_bg_rand_prop_2hr_thin_ind_pt_6hr_buf.csv"))
data_ne_smp_2hr <- read_csv(str_glue("{path}/cov_data_ne_bg_rand_prop_2hr_thin_ind_pt_6hr_buf.csv"))

# Filtering for presence points
data_all <- data_all_smp_2hr %>%
  mutate(pres = as.factor(pres),
         land_cover = as.factor(land_cover)) %>%
  filter(pres==1)
data_expl <- data_expl_smp_2hr %>%
  mutate(pres = as.factor(pres),
         land_cover = as.factor(land_cover)) %>%
  filter(pres==1)
data_ne <- data_ne_smp_2hr %>%
  mutate(pres = as.factor(pres),
         land_cover = as.factor(land_cover)) %>%
  filter(pres==1)

covars <- names(data_all)[4:20]

covars_short <- c(str_sub(covars[1:16], 1,-6), covars[17])

data_all_sp <- vect(data_all, geom = c("x", "y"), crs = 'epsg:4326')
data_expl_sp <- vect(data_expl, geom = c("x", "y"), crs = 'epsg:4326')
data_ne_sp <- vect(data_ne, geom = c("x", "y"), crs = 'epsg:4326')


#### Load models ####
# Direct to where models from 10-rf_models were saved
load(str_glue("{path}/../model_global/bg_rand_prop_2hr_thin_mtry3_tr100_minn_800_final/tidy_global_all.RData"))
load(str_glue("{path}/../model_global/bg_rand_prop_2hr_thin_mtry3_tr100_minn_800_final/tidy_global_expl.RData"))
load(str_glue("{path}/../model_global/bg_rand_prop_2hr_thin_mtry3_tr100_minn_800_final/tidy_global_ne.RData"))

ranger_a_obj <- extract_fit_parsnip(final_a_model)$fit
ranger_e_obj <- extract_fit_parsnip(final_e_model)$fit
ranger_ne_obj <- extract_fit_parsnip(final_ne_model)$fit


#### Predictions ####

## Splitting polygons for predictions
# 10km buffered polygons w ID
polys_10k <- vect(str_glue("{path}../jags_reg_int_ind_polygons_10k_buf.shp")) %>%
  disagg() %>%
  `values<-`(NULL)
polys_10k$poly_ID <- 1:length(polys_10k)

for(i in 1:length(polys_10k)) {
  assign(str_glue("poly{i}"), subset(polys_10k, polys_10k$poly_ID==i))
}

## Polygon scaled rasters
stack_covar_p_sc <- list.files(str_glue("{path}../covariate_data/cov_stack_polygon_scaled/"), 
                               pattern = ".tif",
                               all.files = TRUE,
                               full.names = TRUE) %>%
  rast()

for(i in covars_short[1:16]) {
  assign(str_glue("{i}_p_sc"), c(stack_covar_p_sc[i]))
}

lc <- list.files(str_glue("{path}../covariate_data/landcover/250m/"),
                 pattern = "lc_2",
                 all.files = TRUE, 
                 full.names = TRUE) %>%
  rast() %>%
  `names<-`(str_glue("lc_{2000:2016}"))

## Using 2014 (for dynamic covariates) to generate predictions

predict_stack <- c(dist_roads_p_sc, dens_roads_p_sc, slope_p_sc, elev_p_sc,
                   pop_dens_p_sc[[15]], hii_p_sc[[15]], gpp_p_sc[[15]],
                   pct_tc_p_sc[[15]], dist_tc_p20_p_sc[[15]], dist_tc_p40_p_sc[[15]],
                   dist_non_tc_p20_p_sc[[15]], dist_non_tc_p40_p_sc[[15]],
                   dist_nat_p_sc[[15]], dist_non_nat_p_sc[[15]],
                   dist_wat_perm_p_sc[[15]],dist_wat_seas_p_sc[[15]],
                   lc[[15]])

names(predict_stack)[5:16] <- str_sub(names(predict_stack)[5:16],1,-11)
names(predict_stack)[5:16] <- str_glue("{names(predict_stack)[5:16]}_p_sc")
names(predict_stack)[17] <- "land_cover"
names(predict_stack)


### Function for poly prediction
pr_poly <- function(poly, mod) {
  p <- get(str_glue("poly{poly}"))
  s <- predict_stack %>% crop(p) %>% terra::mask(p)
  lc <- s[[17]]
  lc[lc==0] <- NA
  f <- c(s[[1:16]], lc)
  f[is.na(f)] <- 7                          # set NA areas to rare lc/value (can't have NAs)
  
  prd <- terra::predict(f, mod, factors = "land_cover", na.rm = T)
  prd_c <- prd[[2]] %>% crop(p) %>% terra::mask(p)      # 2nd layer is prob of pres (1)
  
  return(prd_c)
}


## Poly 3 visualiztion

par(mfrow = c(1,1), mar = c(3,4,4,1))
plot(predict_stack$land_cover%>%crop(poly3))
pr_poly3_a <- pr_poly(3, ranger_a_obj)
plot(pr_poly3_a, main = 'All points'); plot(data_all_sp%>%crop(poly3), cex = 0.3, add=T)
pr_poly3_e <- pr_poly(3, ranger_e_obj)
plot(pr_poly3_e, main = 'Exploratory points'); plot(data_expl_sp%>%crop(poly3), cex = 0.3, add=T)
pr_poly3_ne <- pr_poly(3, ranger_ne_obj)
plot(pr_poly3_ne, main = 'Non-exploratory'); plot(data_ne_sp%>%crop(poly3), cex = 0.3, add=T)


### Predicting for all polys

pr_polys_all <- lapply(1:length(polys_10k), function(x) pr_poly(x, ranger_a_obj))
pr_polys_expl <- lapply(1:length(polys_10k), function(x) pr_poly(x, ranger_e_obj))
pr_polys_ne <- lapply(1:length(polys_10k), function(x) pr_poly(x, ranger_ne_obj))

# Exporting "all" projection tests
# for (i in 1:length(pr_polys_all)) {
#   export <- pr_polys_all[[i]]
#   filename <- str_glue("polys_yr_projections/tests/smp_2hr_bg1_rand_group_id/poly{i}_2014_samp_2hr_bg1_rand_6hr_buf_id_group.tif")
#   writeRaster(export, filename)
# }

# Visualizing several target polygons
par(mfrow = c(4,3), mar = c(3,4,4,1))

for(i in c(1,3,6,9)) {
  plot(pr_polys_all[[i]], main = "All points")
  plot(pr_polys_expl[[i]], main = "Exploratory points")
  plot(pr_polys_ne[[i]], main = "Non-exploratory points")
}

# Combining for export as raster stack
pr_all <- do.call(merge, pr_polys_all)
pr_expl <- do.call(merge, pr_polys_expl)
pr_ne <- do.call(merge, pr_polys_ne)

# Comparing means
lapply(c(pr_all, pr_expl, pr_ne), values) %>% lapply(function(x) mean(x, na.rm=T))
## "All points model predicts generally higher probabilities (likely due to having the most data input)

pr_stack <- c(pr_all, pr_expl, pr_ne) %>%
  `names<-`(c("all", "expl", "ne"))


#### Exports ####
# # Exporting stack
# writeRaster(pr_stack, c("model_global/pr_all_6hr_buf_inv.tif",
#                         "model_global/pr_expl_6hr_buf_inv.tif",
#                         "model_global/pr_ne_6hr_buf_inv.tif"))

# # Exporting separately as well
# for (i in 1:length(pr_polys_all)) {
#   export <- pr_polys_all[[i]]
#   filename <- str_glue("polys_yr_projections/bg_rand_prop_2hr_thin_mtry3_tr100_minn_800_all/poly{i}_2014_samp_2hr_bg1_rand_6hr_buf_id_group.tif")
#   writeRaster(export, filename)
# }
# 
# for (i in 1:length(pr_polys_expl)) {
#   export <- pr_polys_expl[[i]]
#   filename <- str_glue("polys_yr_projections/bg_rand_prop_2hr_thin_mtry3_tr100_minn_800_expl/poly{i}_2014_samp_2hr_bg1_rand_6hr_buf_id_group.tif")
#   writeRaster(export, filename)
# }
# 
# for (i in 1:length(pr_polys_ne)) {
#   export <- pr_polys_ne[[i]]
#   filename <- str_glue("polys_yr_projections/bg_rand_prop_2hr_thin_mtry3_tr100_minn_800_ne/poly{i}_2014_samp_2hr_bg1_rand_6hr_buf_id_group.tif")
#   writeRaster(export, filename)
# }
# 
