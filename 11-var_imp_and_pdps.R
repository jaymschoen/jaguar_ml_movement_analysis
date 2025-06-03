# Jaguar Machine Learning Movement Analysis
## 11 - Variable Importance and Partial Dependence Plots

### Model comparisons and predictions
### Generating:
      # - variable importance comparison plot
      # - partial dependence plots


#### Setup ####
library(terra)
library(ranger)
library(tidyverse)
library(tidymodels)
library(pdp)
library(plotly)

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

covars <- names(data_all)[4:20]

covars_short <- c(str_sub(covars[1:16], 1,-6), covars[17])

#### Load models ####
# Direct to where models from 10-rf_models were saved
load(str_glue("{path}/../model_global/bg_rand_prop_2hr_thin_mtry3_tr100_minn_800_final/tidy_global_all.RData"))
load(str_glue("{path}/../model_global/bg_rand_prop_2hr_thin_mtry3_tr100_minn_800_final/tidy_global_expl.RData"))
load(str_glue("{path}/../model_global/bg_rand_prop_2hr_thin_mtry3_tr100_minn_800_final/tidy_global_ne.RData"))

ranger_a_obj <- extract_fit_parsnip(final_a_model)$fit
ranger_e_obj <- extract_fit_parsnip(final_e_model)$fit
ranger_ne_obj <- extract_fit_parsnip(final_ne_model)$fit


#### Variable Importance comparison ####
var_imp <- tibble(var = covars,
                  all = ranger_a_obj$variable.importance,
                  expl = ranger_e_obj$variable.importance,
                  ne = ranger_ne_obj$variable.importance) %>%
  # need to standardize (add to 1) for comparability
  mutate(all = all/sum(all),
         expl = expl/sum(expl),
         ne = ne/sum(ne))

vars_grph <- c("Pop. dens.", "HII", "GPP", "%TC",
               "Dist - mod. TC", "Dist - high TC",
               "Dist - very low TC", "Dist - low TC",
               "Dist - natural", "Dist - matrix",
               "Dist - perm. water", "Dist - seas. water",
               "Dist - maj. roads", "Dens - min roads",
               "Elevation", "Slope", "Land cover")

var_imp_plot <- var_imp %>%
  mutate(var = factor(var, levels = covars)) %>%
  pivot_longer(-var, values_to = "perm_imp", names_to = "Model") %>%
  ggplot(., aes(x = var, y = perm_imp)) +
  theme_bw() +
  geom_col(aes(fill = Model), position = "dodge") +
  scale_fill_manual(values = c("black", "blue", "grey"),
                    labels = c("All", "Exploratory", "Non-exploratory")) +
  labs(title = "Variable Importance Comparison") +
  xlab("Variable") +
  ylab("Permutation importance (relative)") +
  theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust = 1)) +
  scale_x_discrete(label = vars_grph)
var_imp_plot

# ggsave("model_global/bg_rand_prop_2hr_thin_mtry3_tr100_minn_800-final/vi_comparison.png",
#        width = 9,
#        height = 6,
#        dpi = 350)


#### PDPs ####
# Normalized partial dependence plots (normalized for inter-model comparability)

# function to create partial dependency plots for each predictor (imported)
pd_vis_all_tidy <- lapply(covars,
                          function(x) {partial(object = ranger_a_obj,
                                               train = data_all,
                                               # train = data_a_train,
                                               pred.var = x,
                                               which.class = 2,  # pres=1 is response
                                               chull = TRUE,
                                               prob = TRUE,
                                               progress = TRUE,
                                               trim.outliers = TRUE)
                          })

pd_vis_expl_tidy <- lapply(covars,
                           function(x) {partial(object = ranger_e_obj,
                                                train = data_expl,
                                                # train = data_e_train,
                                                pred.var = x,
                                                which.class = 2,  # pres=1 is response
                                                chull = TRUE,
                                                prob = TRUE,
                                                progress = TRUE,
                                                trim.outliers = TRUE)
                           })

pd_vis_ne_tidy <- lapply(covars,
                         function(x) {partial(object = ranger_ne_obj,
                                              train = data_ne,
                                              # train = data_ne_train,
                                              pred.var = x,
                                              which.class = 2,  # pres=1 is response
                                              chull = TRUE,
                                              prob = TRUE,
                                              progress = TRUE,
                                              trim.outliers = TRUE)
                         })

# Removing unimportant variables (no pdp generated)
pd_vis_all_tidy <- pd_vis_all_tidy[-c(5,9,14)]
pd_vis_expl_tidy <- pd_vis_expl_tidy[-c(5,9,14)]
pd_vis_ne_tidy <- pd_vis_ne_tidy[-c(5,9,14)]


# # To save pdp data
# save(pd_vis_all_tidy, file = str_glue("{path}../model_global/bg_rand_prop_2hr_thin_mtry3_tr100_minn_800_final/pd_vis_all_tidy.RData"))
# save(pd_vis_expl_tidy, file = str_glue("{path}../model_global/bg_rand_prop_2hr_thin_mtry3_tr100_minn_800_final/pd_vis_expl_tidy.RData"))
# save(pd_vis_ne_tidy, file = str_glue("{path}../model_global/bg_rand_prop_2hr_thin_mtry3_tr100_minn_800_final/pd_vis_ne_tidy.RData"))

# # To import pdp data
load(str_glue("{path}../model_global/bg_rand_prop_2hr_thin_mtry3_tr100_minn_800_final/pd_vis_all_tidy.RData"))
load(str_glue("{path}../model_global/bg_rand_prop_2hr_thin_mtry3_tr100_minn_800_final/pd_vis_expl_tidy.RData"))
load(str_glue("{path}../model_global/bg_rand_prop_2hr_thin_mtry3_tr100_minn_800_final/pd_vis_ne_tidy.RData"))

# Comparing response curves

## adding full variable names for pdp graphics
vars_pdp <- c("Population density", "HII", "GPP", "%TC",
              "Distance - high TC", "Distance - very low TC", "Distance - low TC",
              "Distance - matrix", "Distance - perm. water", "Distance - seas. water",
              "Distance - maj. roads", "Elevation", "Slope", "Land cover")
# vars_pdp <- c("Population density", "HII", "GPP", "%TC",
#                "Distance - very low TC", "Distance - low TC",
#                "Distance - matrix", "Distance - perm. water", "Distance - seas. water",
#                "Distance - maj. roads", "Elevation", "Slope", "Land cover")

# Normalizing responses to plot together
# Scaled and centered so each curve is relative to others based on sd of pdp values
pd_vis_all_norm <- pd_vis_all_tidy
for(i in 1:length(pd_vis_all_norm)) {
  pd_vis_all_norm[[i]]$yhat <- scale(pd_vis_all_norm[[i]]$yhat) 
  colnames(pd_vis_all_norm[[i]]) <- c(vars_pdp[i], "yhat") 
}

pd_vis_expl_norm <- pd_vis_expl_tidy
for(i in 1:length(pd_vis_expl_norm)) {
  pd_vis_expl_norm[[i]]$yhat<- scale(pd_vis_expl_norm[[i]]$yhat)
  colnames(pd_vis_expl_norm[[i]]) <- c(vars_pdp[i], "yhat") 
}

pd_vis_ne_norm <- pd_vis_ne_tidy
for(i in 1:length(pd_vis_ne_norm)) {
  pd_vis_ne_norm[[i]]$yhat <- scale(pd_vis_ne_norm[[i]]$yhat)
  colnames(pd_vis_ne_norm[[i]]) <- c(vars_pdp[i], "yhat") 
}


# Plotting all covariate pdps together

## to save
# png("model_global/bg_rand_prop_2hr_thin_mtry3_tr100_minn_800_final/pdp_cov_plots_norm.png",
#      height = 1000, width = 1500, pointsize = 6, res = 220)

# svg("model_global/bg_rand_prop_2hr_thin_mtry3_tr100_minn_800_final/pdp_cov_plots_norm2.svg",
#     height = 10, width = 15, pointsize = 12)

par(mfrow = c(4,4), mar = c(4,3,3,1))
for(i in 1:length(pd_vis_all_norm)) {
  all_l <- pd_vis_all_norm
  expl_l <- pd_vis_expl_norm
  ne_l <- pd_vis_ne_norm
  
  # expl_l[[5]]$yhat <- pd_vis_ne_norm[[5]]$yhat    #covering up NAs
  
  if(isFALSE(is.factor(all_l[[i]][,1]))) {
    plot(all_l[[i]],
         type = "l",
         lwd = 3,
         col = "black",
         ylim = c(min(rbind(all_l[[i]]$yhat,
                            expl_l[[i]]$yhat,
                            ne_l[[i]]$yhat)),
                  max(rbind(all_l[[i]]$yhat,
                            expl_l[[i]]$yhat,
                            ne_l[[i]]$yhat))))}
  else {
    plot(pd_vis_all_norm[[i]],
         type = "p",
         pch = 95,
         cex = 4,
         col = "black",
         ylim = c(min(rbind(all_l[[i]]$yhat,
                            expl_l[[i]]$yhat,
                            ne_l[[i]]$yhat), na.rm=T),
                  max(rbind(all_l[[i]]$yhat,
                            expl_l[[i]]$yhat,
                            ne_l[[i]]$yhat), na.rm=T)))};
  if(isFALSE(is.factor(pd_vis_expl_norm[[i]][,1]))) {
    lines(pd_vis_expl_norm[[i]],
          lwd = 3,
          col = "blue")}
  else {
    points(pd_vis_expl_norm[[i]],
           pch = 95,
           cex = 4,
           col = "blue")};
  if(isFALSE(is.factor(pd_vis_ne_norm[[i]][,1]))) {
    lines(pd_vis_ne_norm[[i]],
          lwd = 3,
          col = "grey")}
  else {
    points(pd_vis_ne_norm[[i]],
           pch = 95,
           cex = 4,
           col = "grey")}
}

# dev.off()

par(mfrow = c(1,1))


#### Variable Interactions ####

##### 3d plot function
# Interactive 3D partial dependence plot with coloring scale

int_plot3d <- function(pd, var1, var2) {
  # Interpolate the partial dependence values
  dens <- akima::interp(x = get(pd)[,var1],
                        y = get(pd)[,var2],
                        z = get(pd)[,"yhat"])
  # 3D partial dependence plot with a coloring scale
  p <- plot_ly(x = dens$x, 
               y = dens$y, 
               z = dens$z,
               colors = c("red2", "yellow3", "green4", "blue3"),
               type = "surface")
  # Add axis labels added individually
  # This default plot seems to be switching x and y axes, so I'm correcting the labels below
}



##### Distance to non-natural areas (matrix) and % tree cover ####

## All pts
pd_all_int_dist_mat_tc <- partial(object = ranger_a_obj,
                                  # train = data_all,
                                  train = data_a_train, # results the same for subset
                                  pred.var = c("dist_non_nat_p_sc", "pct_tc_p_sc"),
                                  which.class = 2,  # pres=1 is response
                                  chull = TRUE,
                                  prob = TRUE,
                                  progress = TRUE,
                                  trim.outliers = TRUE)

## Expl pts
pd_expl_int_dist_mat_tc <- partial(object = ranger_e_obj,
                                   # train = data_expl,
                                   train = data_e_train, # results the same for subset
                                   pred.var = c("dist_non_nat_p_sc", "pct_tc_p_sc"),
                                   which.class = 2,  # pres=1 is response
                                   chull = TRUE,
                                   prob = TRUE,
                                   progress = TRUE,
                                   trim.outliers = TRUE)

## Non-expl pts
pd_ne_int_dist_mat_tc <- partial(object = ranger_ne_obj,
                                 # train = data_ne,
                                 train = data_ne_train, # results the same for subset
                                 pred.var = c("dist_non_nat_p_sc", "pct_tc_p_sc"),
                                 which.class = 2,  # pres=1 is response
                                 chull = TRUE,
                                 prob = TRUE,
                                 progress = TRUE,
                                 trim.outliers = TRUE)

# Save/load exported pdp interation results
# save(pd_all_int_dist_mat_tc,
#      file = str_glue("{path}/../model_global/bg_rand_prop_2hr_thin_mtry3_tr100_minn_800_final/pd_all_int_dist_mat_tc.RData"))
# save(pd_expl_int_dist_mat_tc,
#      file = str_glue("{path}/../model_global/bg_rand_prop_2hr_thin_mtry3_tr100_minn_800_final/pd_expl_int_dist_mat_tc.RData"))
# save(pd_ne_int_dist_mat_tc,
#      file = str_glue("{path}/../model_global/bg_rand_prop_2hr_thin_mtry3_tr100_minn_800_final/pd_ne_int_dist_mat_tc.RData"))

load(str_glue("{path}../model_global/bg_rand_prop_2hr_thin_mtry3_tr100_minn_800_final/pd_all_int_dist_mat_tc.RData"))
load(str_glue("{path}../model_global/bg_rand_prop_2hr_thin_mtry3_tr100_minn_800_final/pd_expl_int_dist_mat_tc.RData"))
load(str_glue("{path}../model_global/bg_rand_prop_2hr_thin_mtry3_tr100_minn_800_final/pd_ne_int_dist_mat_tc.RData"))

plot_all_int_dist_mat_tc <- int_plot3d("pd_all_int_dist_mat_tc", "dist_non_nat_p_sc", "pct_tc_p_sc") %>%
  layout(scene = list(xaxis = list(title = "% Tree Cover"),
                      yaxis = list(title = "Distance to Matrix"),
                      zaxis = list(title = "Partial Dependence")))
plot_all_int_dist_mat_tc

plot_expl_int_dist_mat_tc <- int_plot3d("pd_expl_int_dist_mat_tc", "dist_non_nat_p_sc", "pct_tc_p_sc")%>%
  layout(scene = list(xaxis = list(title = "% Tree Cover"),
                      yaxis = list(title = "Distance to Matrix"),
                      zaxis = list(title = "Partial Dependence"))) %>%
  colorbar(limits = c(0.21,0.27))
plot_expl_int_dist_mat_tc

plot_ne_int_dist_mat_tc <- int_plot3d("pd_ne_int_dist_mat_tc", "dist_non_nat_p_sc", "pct_tc_p_sc") %>% 
  layout(scene = list(xaxis = list(title = "% Tree Cover"),
                      yaxis = list(title = "Distance to Matrix"),
                      zaxis = list(title = "Partial Dependence")))
plot_ne_int_dist_mat_tc



##### % tree cover & distance to tree cover ####
## All pts
pd_all_int_pct_tc_dist_tc_p40 <- partial(object = ranger_a_obj,
                                         train = data_a_train,
                                         pred.var = c("pct_tc_p_sc", "dist_tc_p40_p_sc"),
                                         which.class = 2,  # pres=1 is response
                                         chull = TRUE,
                                         prob = TRUE,
                                         progress = TRUE,
                                         trim.outliers = TRUE)

## Expl pts
pd_expl_int_pct_tc_dist_tc_p40 <- partial(object = ranger_e_obj,
                                          train = data_e_train,
                                          pred.var = c("pct_tc_p_sc", "dist_tc_p40_p_sc"),
                                          which.class = 2,  # pres=1 is response
                                          chull = TRUE,
                                          prob = TRUE,
                                          progress = TRUE,
                                          trim.outliers = TRUE)

## Non-expl pts
pd_ne_int_pct_tc_dist_tc_p40 <- partial(object = ranger_ne_obj,
                                        train = data_ne_train,
                                        pred.var = c("pct_tc_p_sc", "dist_tc_p40_p_sc"),
                                        which.class = 2,  # pres=1 is response
                                        chull = TRUE,
                                        prob = TRUE,
                                        progress = TRUE,
                                        trim.outliers = TRUE)


plot_all_int_pct_tc_dist_tc_p40 <- int_plot3d("pd_all_int_pct_tc_dist_tc_p40", "pct_tc_p_sc", "dist_tc_p40_p_sc") %>% 
  layout(scene = list(xaxis = list(title = "Distance to Tree Cover (40th %ile)"),
                      yaxis = list(title = "% Tree Cover"),
                      zaxis = list(title = "Partial Dependence")))
plot_all_int_pct_tc_dist_tc_p40

plot_expl_int_pct_tc_dist_tc_p40 <- int_plot3d("pd_expl_int_pct_tc_dist_tc_p40", "pct_tc_p_sc", "dist_tc_p40_p_sc")%>% 
  layout(scene = list(xaxis = list(title = "Distance to Tree Cover (40th %ile)"),
                      yaxis = list(title = "% Tree Cover"),
                      zaxis = list(title = "Partial Dependence")))
plot_expl_int_pct_tc_dist_tc_p40

plot_ne_int_pct_tc_dist_tc_p40 <- int_plot3d("pd_ne_int_pct_tc_dist_tc_p40", "pct_tc_p_sc", "dist_tc_p40_p_sc") %>% 
  layout(scene = list(xaxis = list(title = "Distance to Tree Cover (40th %ile)"),
                      yaxis = list(title = "% Tree Cover"),
                      zaxis = list(title = "Partial Dependence")))
plot_ne_int_pct_tc_dist_tc_p40

##### Distance to roads and % tree cover ####
## All pts
pd_all_int_dist_rd_tc <- partial(object = ranger_a_obj,
                                 train = data_a_train,
                                 pred.var = c("dist_roads_p_sc", "pct_tc_p_sc"),
                                 which.class = 2,  # pres=1 is response
                                 chull = TRUE,
                                 prob = TRUE,
                                 progress = TRUE,
                                 trim.outliers = TRUE)

## Expl pts
pd_expl_int_dist_rd_tc <- partial(object = ranger_e_obj,
                                  train = data_e_train,
                                  pred.var = c("dist_roads_p_sc", "pct_tc_p_sc"),
                                  which.class = 2,  # pres=1 is response
                                  chull = TRUE,
                                  prob = TRUE,
                                  progress = TRUE,
                                  trim.outliers = TRUE)

## Non-expl pts
pd_ne_int_dist_rd_tc <- partial(object = ranger_ne_obj,
                                train = data_ne_train,
                                pred.var = c("dist_roads_p_sc", "pct_tc_p_sc"),
                                which.class = 2,  # pres=1 is response
                                chull = TRUE,
                                prob = TRUE,
                                progress = TRUE,
                                trim.outliers = TRUE)

plot_all_int_dist_rd_tc <- int_plot3d("pd_all_int_dist_rd_tc", "dist_roads_p_sc", "pct_tc_p_sc") %>% 
  layout(scene = list(xaxis = list(title = "% Tree Cover"),
                      yaxis = list(title = "Distance to Roads"),
                      zaxis = list(title = "Partial Dependence")))
plot_all_int_dist_rd_tc

plot_expl_int_dist_rd_tc <- int_plot3d("pd_expl_int_dist_rd_tc", "dist_roads_p_sc", "pct_tc_p_sc") %>%
  layout(scene = list(xaxis = list(title = "% Tree Cover"),
                      yaxis = list(title = "Distance to Roads"),
                      zaxis = list(title = "Partial Dependence")))
plot_expl_int_dist_rd_tc

plot_ne_int_dist_rd_tc <- int_plot3d("pd_ne_int_dist_rd_tc", "dist_roads_p_sc", "pct_tc_p_sc") %>% 
  layout(scene = list(xaxis = list(title = "% Tree Cover"),
                      yaxis = list(title = "Distance to Roads"),
                      zaxis = list(title = "Partial Dependence")))
plot_ne_int_dist_rd_tc


##### Distance to non-tc p20 & seasonal water
## All pts
pd_all_int_dist_non_tc_p20_wat_seas <- partial(object = ranger_a_obj,
                                               train = data_a_train,
                                               pred.var = c("dist_non_tc_p20_p_sc", "dist_wat_seas_p_sc"),
                                               which.class = 2,  # pres=1 is response
                                               chull = TRUE,
                                               prob = TRUE,
                                               progress = TRUE,
                                               trim.outliers = TRUE)

## Expl pts
pd_expl_int_dist_non_tc_p20_wat_seas <- partial(object = ranger_e_obj,
                                                train = data_e_train,
                                                pred.var = c("dist_non_tc_p20_p_sc", "dist_wat_seas_p_sc"),
                                                which.class = 2,  # pres=1 is response
                                                chull = TRUE,
                                                prob = TRUE,
                                                progress = TRUE,
                                                trim.outliers = TRUE)

## Non-expl pts
pd_ne_int_dist_non_tc_p20_wat_seas <- partial(object = ranger_ne_obj,
                                              train = data_ne_train,
                                              pred.var = c("dist_non_tc_p20_p_sc", "dist_wat_seas_p_sc"),
                                              which.class = 2,  # pres=1 is response
                                              chull = TRUE,
                                              prob = TRUE,
                                              progress = TRUE,
                                              trim.outliers = TRUE)

plot_all_int_dist_non_tc_p20_wat_seas <- int_plot3d("pd_all_int_dist_non_tc_p20_wat_seas", "dist_non_tc_p20_p_sc", "dist_wat_seas_p_sc") %>% 
  layout(scene = list(xaxis = list(title = "Distance to Water (seasonal)"),
                      yaxis = list(title = "Distance to Low Tree Cover (20th %ile)"),
                      zaxis = list(title = "Partial Dependence")))
plot_all_int_dist_non_tc_p20_wat_seas

plot_expl_int_dist_non_tc_p20_wat_seas <- int_plot3d("pd_expl_int_dist_non_tc_p20_wat_seas", "dist_non_tc_p20_p_sc", "dist_wat_seas_p_sc") %>% 
  layout(scene = list(xaxis = list(title = "Distance to Water (seasonal)"),
                      yaxis = list(title = "Distance to Low Tree Cover (20th %ile)"),
                      zaxis = list(title = "Partial Dependence")))
plot_expl_int_dist_non_tc_p20_wat_seas

plot_ne_int_dist_non_tc_p20_wat_seas <- int_plot3d("pd_ne_int_dist_non_tc_p20_wat_seas", "dist_non_tc_p20_p_sc", "dist_wat_seas_p_sc") %>% 
  layout(scene = list(xaxis = list(title = "Distance to Water (seasonal)"),
                      yaxis = list(title = "Distance to Low Tree Cover (20th %ile)"),
                      zaxis = list(title = "Partial Dependence")))
plot_ne_int_dist_non_tc_p20_wat_seas
