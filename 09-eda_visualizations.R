# Jaguar Machine Learning Movement Analysis
## 09 - Exploratory Data Analysis and Visualizations

### Visualizing covariate data prior to building models


#### Setup ####
library(terra)
library(tidyverse)

# Set to directory with cloned repo
setwd("~/Columbia E3B/PhD/Dissertation/Ch3-jaguar_movement/scripts/git_push_R/")

# Make this path where you stored environmental covariate data
path <- "C:/Users/jayms/Documents/Columbia E3B/PhD/Dissertation/Ch3-jaguar_movement/covariate_data/"

#### Importing data ####
# Importing 2-hr interval thinned data

data_all_smp_2hr <- read_csv(str_glue("{path}/cov_data_all_bg_rand_prop_2hr_thin_ind_pt_6hr_buf.csv"))
data_expl_smp_2hr <- read_csv(str_glue("{path}/cov_data_expl_bg_rand_prop_2hr_thin_ind_pt_6hr_buf.csv"))
data_ne_smp_2hr <- read_csv(str_glue("{path}/cov_data_ne_bg_rand_prop_2hr_thin_ind_pt_6hr_buf.csv"))

data_all <- data_all_smp_2hr %>%
  mutate(pres = as.factor(pres),
         land_cover = as.factor(land_cover)) 
data_expl <- data_expl_smp_2hr %>%
  mutate(pres = as.factor(pres),
         land_cover = as.factor(land_cover))
data_ne <- data_ne_smp_2hr %>%
  mutate(pres = as.factor(pres),
         land_cover = as.factor(land_cover))

# Export final pres/bg covariate extractions as vectors for visualization

final_vects <- lapply(list(data_all_smp_2hr, data_expl_smp_2hr, data_ne_smp_2hr),
                      function(x) vect(x, geom = c('x', 'y'))) %>%
  `names<-`(c('data_all_smp_2hr', 'data_expl_smp_2hr', 'data_ne_smp_2hr'))

# for(i in 1:length(final_vects)){
#   lyr = final_vects[[i]]
#   name = names(final_vects[i])
#   writeVector(lyr, str_glue("{path}/{name}_bg_rand_prop_6_hr_buf.shp"))
# }



## For sensitivity to reg_interval (Supplementary sensitivity analysis)

# Filtering jags_reg_int for 0 hr/2 hr/4 hr buffer to test sensitivity of models to this cutoff
# Use ind_reg_int dataframe from 01-data_prep to select subset of jags within cutoff

# data_all <- data_all %>%
#   filter(ID %in% jags_reg_int)   # manipulate jags_reg_int to select different subsets
# 
# data_expl <- data_expl %>%
#   filter(ID %in% jags_reg_int)   # manipulate jags_reg_int to select different subsets
# 
# data_ne <- data_ne %>%
#   filter(ID %in% jags_reg_int)   # manipulate jags_reg_int to select different subsets



#### EDA ####

vars_grph <- c("Pop. dens.", "HII", "GPP", "%TC",
               "Dist - mod. TC", "Dist - high TC",
               "Dist - very low TC", "Dist - low TC",
               "Dist - natural", "Dist - matrix",
               "Dist - perm. water", "Dist - seas. water",
               "Dist - maj. roads", "Dens - min roads",
               "Elevation", "Slope", "Land cover")

### Correlation
names(data_all)
c <- cor(data_all[,-c(1:3,20:24)])
c[c>0.6]
View(c)
## no correlations >0.6 (only diagonals ==1)


##### Densities ####

# Plotting everything together
par(mfrow = c(4,4), mar = c(4,3,3,1))

for(i in 4:19) {
  d <- as.data.frame(data_all)[,i]
  
  plot(density(d, width = 0.25),
       xlim = c(quantile(d, probs = seq(0,1,0.01))[2],
                quantile(d, probs = seq(0,1,0.01))[97]),
       main = vars_grph[i-3],
       xlab = " ",
       ylab = "density",
       lwd = 2)}

# Land cover
par(mfrow = c(1,1), mar = c(4,3,3,1))

ggplot(data_all, aes(x = as.factor(land_cover))) +
  geom_bar() +
  theme_bw() +
  xlab("Land Cover") +
  ylab("# Points") +
  scale_x_discrete(labels = c("Forested", "Natural non-forested", "Pasture",
                              "Agriculture", "Mixed Ag/Pasture", "Forest plantation",
                              "Non-vegetated", "Water", "Floodplain")) +
  theme(axis.text.x = element_text(angle = 20, vjust = 0.4))

# For density inserts (matching x axis of partial dependency plots produced later)
par(mfrow = c(4,4), mar = c(4,3,3,1))

plot(density(data_all$pop_dens_p_sc, width = 0.25), 
     xlim = c(0,1.5),
     xlab = " ",
     ylab = "density",
     lwd = 2)

plot(density(data_all$hii_p_sc, width = 0.25), 
     xlim = c(0,1.75),
     xlab = " ",
     ylab = "density",
     lwd = 2)

plot(density(data_all$gpp_p_sc, width = 0.25), 
     xlim = c(0,1.5),
     xlab = " ",
     ylab = "density",
     lwd = 2)

plot(density(data_all$pct_tc_p_sc, width = 0.25), 
     xlim = c(0,2.25),
     xlab = " ",
     ylab = "density",
     lwd = 2)

plot(density(data_all$dist_tc_p40_p_sc, width = 0.25), 
     xlim = c(0,0.6),
     xlab = " ",
     ylab = "density",
     lwd = 2)

plot(density(data_all$dist_non_tc_p20_p_sc, width = 0.25), 
     xlim = c(0,2.5),
     xlab = " ",
     ylab = "density",
     lwd = 2)

plot(density(data_all$dist_non_tc_p40_p_sc, width = 0.25), 
     xlim = c(0,2.75),
     xlab = " ",
     ylab = "density",
     lwd = 2)

plot(density(data_all$dist_non_nat_p_sc, width = 0.25), 
     xlim = c(0,2),
     xlab = " ",
     ylab = "density",
     lwd = 2)

plot(density(data_all$dist_wat_perm_p_sc, width = 0.25), 
     xlim = c(0,2.25),
     xlab = " ",
     ylab = "density",
     lwd = 2)

plot(density(data_all$dist_wat_seas_p_sc, width = 0.25), 
     xlim = c(0,2.25),
     xlab = " ",
     ylab = "density",
     lwd = 2)

plot(density(data_all$dist_roads_p_sc, width = 0.25), 
     xlim = c(0,2),
     xlab = " ",
     ylab = "density",
     lwd = 2)

plot(density(data_all$elev_p_sc, width = 0.25), 
     xlim = c(0.85,1.15),
     xlab = " ",
     ylab = "density",
     lwd = 2)

plot(density(data_all$slope_p_sc, width = 0.25), 
     xlim = c(0,2.25),
     xlab = " ",
     ylab = "density",
     lwd = 2)

hist(as.numeric(data_all$land_cover), 
     breaks = seq(0,9, 0.5), 
     col = "black",
     xaxt = "n")
axis(1, at = c(1:9))


##### Presence/background distributions ####
# Visualizing covariate distributions for presence vs. background points for all data + subsets
par(mfrow = c(1,1), mar = c(4,3,3,1))

###### Continuous variables ####
# All pts
var_plot_all <- data_all %>%
  dplyr::select(-c(u_id, x, y, land_cover, poly_ID, ID, year)) %>%
  slice_sample(n = 1e5) %>%
  pivot_longer(cols = -pres) %>%
    ggplot(., aes(pres, value)) +
    theme_bw() +
    geom_violin(col = "grey80") +
    geom_boxplot(alpha = 0.3, col = "grey20") +
    facet_wrap(~name, scales = "free_y",
               nrow = 3, ncol = 6) + 
    labs(title = "All pts")
plot(var_plot_all)

# Exploratory pts
var_plot_expl <- data_expl %>%
  dplyr::select(-c(x, y, land_cover, poly_ID, ID, year)) %>%
  pivot_longer(cols = -pres) %>%
    ggplot(., aes(pres, value)) +
    theme_bw() +
    geom_violin(col = "grey80") +
    geom_boxplot(alpha = 0.3, col = "grey20") +
    facet_wrap(~name, scales = "free_y",
               nrow = 2, ncol = 8) +
    labs(title = "Exploratory pts")
plot(var_plot_expl)

# Non-exploratory pts
var_plot_ne <- data_ne %>%
  dplyr::select(-c(x, y, land_cover, poly_ID, ID, year)) %>%
  pivot_longer(cols = -pres) %>%
    ggplot(., aes(pres, value)) +
    theme_bw() +
    geom_violin(col = "grey80") +
    geom_boxplot(alpha = 0.3, col = "grey20") +
    facet_wrap(~name, scales = "free_y",
               nrow = 2, ncol = 8) +
    labs(title = "Non-exploratory pts")
plot(var_plot_ne)

###### Discrete variable (land cover) ####

# Land cover classes proportional to presence (0/1)

prop_lc_all <- data_all %>%
  select(land_cover, pres) %>%
  group_by(pres) %>%
  mutate(total = n()) %>%
  ungroup() %>%
  group_by(land_cover, pres) %>%
  dplyr::reframe(sum = n(),
                 prop = sum/total) %>%
  unique() %>% mutate(model = "all")

prop_lc_exp <- data_expl %>%
  select(land_cover, pres) %>%
  group_by(pres) %>%
  mutate(total = n()) %>%
  ungroup() %>%
  group_by(land_cover, pres) %>%
  dplyr::reframe(sum = n(),
                 prop = sum/total) %>%
  unique() %>% mutate(model = "expl")

prop_lc_ne <- data_ne %>%
  select(land_cover, pres) %>%
  group_by(pres) %>%
  mutate(total = n()) %>%
  ungroup() %>%
  group_by(land_cover, pres) %>%
  dplyr::reframe(sum = n(),
                 prop = sum/total) %>%
  unique() %>% mutate(model = "ne")

prop_lc_comp <- bind_rows(prop_lc_all, prop_lc_exp, prop_lc_ne) %>%
  ggplot(aes(x = land_cover, y = prop, fill = model, alpha = pres)) + 
  theme_bw() +
  geom_col(position = "dodge") +
  scale_fill_manual(values =c("black", "blue", "grey70")) + 
  scale_alpha_manual(values = c(0.5, 1)) +
  scale_x_discrete(labels = c("Forest",
                              "Natural non-forest",
                              "Pasture",
                              "Agriculture",
                              "Mixed pasture/ag",
                              "Forest plantation",
                              "Non-vegetated",
                              "Water",
                              "Floodplain")) +
  xlab("Land cover category") +
  ylab("Proportion") + 
  labs(fill = "Data", alpha = "Pres",
       title = "Land cover comparison")
plot(prop_lc_comp)


##### Water visualization ####

# Land cover analysis shows increased presence in water
## Looking at points classified as lc = water in permanent/seasonal water

wat_plot <- data_all %>%
  filter(pres==1 & land_cover %in% 8) %>%
  select(land_cover, dist_wat_perm_p_sc, dist_wat_seas_p_sc) %>%
  pivot_longer(-land_cover, values_to = "distance") %>%
    ggplot(.,aes(x = land_cover, y = distance)) +
    theme_bw() +
    geom_violin(col = "grey80", fill = "grey90") +
    geom_boxplot(alpha = 0.3, col = "grey20") + 
    facet_wrap(~name)
plot(wat_plot)

# Majority of these presences are in/near seasonal water, not permanent
## Remember when interpreting variable importance