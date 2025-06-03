# Jaguar Machine Learning Movement Analysis
## 10 - Random Forest Models

### Building Random Forest models for "all", "exploratory", and "non-exploratory" data sets


#### Setup ####
library(terra)
library(ranger)
library(tidyverse)
library(tidymodels)
library(modEvA)
library(vip)

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

#### Univariate performance ####
# Using fast fitting ranger models to assess performance of each variable independently (for feature selection)

# All data
rf_univ_a <- list()
univ_oob_a <- tibble(var = covars, oob_a = NA) 
for(i in 1:length(covars)) {
  rf_univ_a[[i]] <- ranger(
    formula = as.formula(str_glue("pres ~ {covars[i]}")), 
    data = data_all[,-c(1:3, 22:24)],
    num.trees = 100,
    importance = "permutation",
    probability = TRUE)
  univ_oob_a[i,2] <- rf_univ_a[[i]]$prediction.error
}
names(rf_univ_a) <- covars

univ_oob_a

# Exploratory data
rf_univ_e <- list()
univ_oob_e <- tibble(var = covars, oob_e = NA) 
for(i in 1:length(covars)) {
  rf_univ_e[[i]] <- ranger(
    formula = as.formula(str_glue("pres ~ {covars[i]}")), 
    data = data_expl[,-c(1:2, 21:23)],
    num.trees = 100,
    importance = "permutation",
    probability = TRUE)
  univ_oob_e[i,2] <- rf_univ_e[[i]]$prediction.error
}
names(rf_univ_e) <- covars

univ_oob_e

# Non-exploratory data

rf_univ_ne <- list()
univ_oob_ne <- tibble(var = covars, oob_ne = NA) 
for(i in 1:length(covars)) {
  rf_univ_ne[[i]] <- ranger(
    formula = as.formula(str_glue("pres ~ {covars[i]}")), 
    data = data_ne[,-c(1:2, 21:23)],
    num.trees = 100,
    importance = "permutation",
    probability = TRUE)
  univ_oob_ne[i,2] <- rf_univ_ne[[i]]$prediction.error
}
names(rf_univ_ne) <- covars

univ_oob_ne

## GPP, dist_roads, HII consistently lowest OOB error
## Other variables varied by data subset

# Will keep all variables for modeling due to generally low univariate OOB error


#### Fitting models ####
# Following tidymodels workflow for final models

##### All pts ####

## Grouping by ID on CV and split
# Ensures independence of training/testing data

## Cross-validating using "all" data
# Will apply optimum parameters for all three models (for comparability)

set.seed(234)

# Splitting and setting up for cross-validation
data_a_split <- group_initial_split(data_all, group = ID, prop = 3/4)
data_a_train <- training(data_a_split)
data_a_test <- testing(data_a_split)
data_a_cv <- group_vfold_cv(data_a_train, group = ID, v = 10, prop = 3/4)

# Specifying formula and pre-processing
data_a_recipe <- recipe(pres ~ ., data = data_all[,-c(1:3,22:24)]) 

data_a_prep <- data_a_recipe %>%
  prep() %>%
  juice()

# Specify model (random forest using "ranger" algorithm)
rf_a_mod <- rand_forest() %>%
  set_args(mtry = tune(), trees = tune(), min_n = tune()) %>%   # when testing all parameters
  # set_args(mtry = 3, trees = 100, min_n = tune()) %>%         # if further tuning min_n
  set_engine("ranger", importance = "permutation") %>%
  set_mode("classification")

# Workflow (put it together)
rf_a_workflow <- workflow() %>%
  add_recipe(data_a_recipe) %>%
  add_model(rf_a_mod)

# Tune parameters
rf_a_grid <- expand.grid(mtry = c(3, 4, 5),
                         trees = c(100, 300),
                         min_n = c(100,400,800))


###### Accuracy + ROC/AUC ####

# Extract cross-validation results
## Using overall accuracy and ROC/AUC for first step
rf_a_tune_results <- rf_a_workflow %>%
  tune_grid(resamples = data_a_cv, # CV object
            grid = rf_a_grid, # grid of all parameters' values to try
            metrics = metric_set(accuracy, roc_auc), # target metrics
            control = control_resamples(save_pred = T))

collect_metrics(rf_a_tune_results) %>%
  arrange(desc(mean)) %>%
  print(n = nrow(.))

### Little separation between CV fold results (all within 1 std error of "top" roc/auc)
### Will defer to Boyce index results

# To view full results
cv_metrics <- collect_metrics(rf_a_tune_results, summarize = F)
View(cv_metrics)


###### Boyce index ####

# Generating predictions for Boyce index analysis (manually)
cv_pred <- collect_predictions(rf_a_tune_results) %>%
  # filter(mtry == 3 & trees == 100) %>%
  as.data.frame()

samp_grid <- expand.grid(id = c(str_glue("Resample0{1:9}"), "Resample10"),
                         min_n = unique(cv_metrics$min_n),
                         .metric = "boyce",
                         .estimator = "binary",
                         .estimate = NA,
                         .config = "Preprocessor1_Model1") %>%
  as_tibble()

for(i in str_pad(1:10, 2, "left", "0")) {
  for(j in unique(cv_metrics$min_n)){
    b <- get(str_glue("boyce_fold_{i}_min_n_{j}"))[[2]]
    samp_grid[samp_grid$id == str_glue("Resample{i}") & samp_grid$min_n == j, 
              ".estimate"] = b
  }
}

cv_metrics_boyce <- cv_metrics %>%
  bind_rows(samp_grid)

# Compare metrics for parameter values
cv_metrics_boyce %>%
  group_by(min_n, .metric) %>%
  # group_by(mtry, .metric) %>%
  # group_by(trees, .metric) %>%
  dplyr::summarize(mn = mean(.estimate),
                   sd = sd(.estimate)) %>%
  print(n=nrow(.))

# Still minimal separation betwen parameter values (within 1 std dev)
## mtry and trees increase computational load; will keep lowest parameters and further tune min_n

### min_n parameter will be particularly important for controlling overfitting
### Will tune min_n parameter on full training data (high variance in CV fold models)


# Testing on full data set 
set.seed(234)
param_a_final <- tibble(mtry=3, trees=100, min_n=800)  # change min_n to test further
param_a_final

# Finalize parameters
rf_a_workflow1 <- rf_a_workflow %>%
  finalize_workflow(param_a_final)
rf_a_workflow1

# Fit model
rf_a_fit <- rf_a_workflow1 %>%
  # fit on the training set and evaluate on test set
  last_fit(data_a_split)

test_a_performance <- rf_a_fit %>% collect_metrics()
test_a_performance

## Further tweaking of min_n parameter results

# min_n = 400: acc = 0.725, roc = 0.642
# min_n = 600: acc = 0.715, roc = 0.648
# min_n = 800: acc = 0.706, roc = 0.655
# min_n = 1000: acc = 0.687, roc = 0.648
# min_n = 1200: acc = 0.699, roc = 0.657
# min_n = 1400: acc = 0.699, roc = 0.655
# min_n = 2000: acc = 0.683, roc = 0.654

# Note general decrease in "accuracy" but increase in ROC/AUC when increasing min_n

### Trusting ROC/AUC more for these data (overall accuracy can be biased/overfit)
## min_n = 800, 1200, 1400, 2000 have highest ROC/AUC

## Final min_n parameter value decided using Boyce index
### Boyce index shown to perform well on presence only data

# generate predictions from the test set
test_a_predictions <- rf_a_fit %>% collect_predictions()

# Boyce index
test_a_boyce <- Boyce(obs = test_a_predictions$pres, pred = test_a_predictions$.pred_1,
                      main = str_glue("All pts_minn_{param_a_final$min_n}"))

## Boyce index results
# min_n = 400: 0.919
# min_n = 600: 0.965
# min_n = 800: 0.990
# min_n = 1000: 0.924
# min_n = 1200: 0.957
# min_n = 1400: 0.960
# min_n = 2000: 0.937


# generate a confusion matrix
test_a_predictions %>%
  conf_mat(truth = pres, estimate = .pred_class) %>%
  summary(event_level = "second") 

# ROC curve
test_a_predictions %>%
  roc_curve(pres, .pred_0) %>%
  ggplot(aes(1 - specificity, y = sensitivity)) +
  geom_path() +
  geom_abline(lty = 3) +
  coord_equal() +
  theme_bw()

# Final model
final_a_model <- fit(rf_a_workflow1, data_all)

# Variable importance
final_a_model %>%
  extract_fit_parsnip() %>%
  vip(num_features = 17)

ranger_a_obj <- extract_fit_parsnip(final_a_model)$fit
ranger_a_obj$variable.importance


# min_n = 800 showed highest Boyce index, also in top models for ROC/AUC


######### FINAL PARAMETERS: mtry = 3, trees = 100, min_n = 800 ###
# Will apply these parameters for all 3 models ("all", "expl", "non-expl")


##### Expl pts ####
set.seed(234)

data_e_split <- group_initial_split(data_expl, group = ID, prop = 3/4)
data_e_recipe <- recipe(pres ~ ., data = data_expl[,-c(1:2,21:23)]) 
data_e_prep <- data_e_recipe %>%
  prep() %>%
  juice()

rf_e_mod <- rand_forest() %>%
  set_args(mtry=3, trees=100, min_n=800) %>%
  set_engine("ranger", importance = "permutation") %>%
  set_mode("classification")

rf_e_workflow <- workflow() %>%
  add_recipe(data_e_recipe) %>%
  add_model(rf_e_mod)

rf_e_fit <- rf_e_workflow %>%
  last_fit(data_e_split)

test_e_performance <- rf_e_fit %>% collect_metrics()
test_e_performance

test_e_predictions <- rf_e_fit %>% collect_predictions()

test_e_boyce <- Boyce(obs = test_e_predictions$pres, pred = test_e_predictions$.pred_1,
                      main = str_glue("Expl pts_minn_{param_a_final$min_n}"))

# Final model
final_e_model <- fit(rf_e_workflow, data_expl)

# Variable importance
final_e_model %>%
  extract_fit_parsnip() %>%
  vip(num_features = 17)

ranger_e_obj <- extract_fit_parsnip(final_e_model)$fit
ranger_e_obj$variable.importance


##### Non-expl pts ####
set.seed(234)

data_ne_split <- group_initial_split(data_ne, group = ID, prop = 3/4)
data_ne_recipe <- recipe(pres ~ ., data = data_ne[,-c(1:2,21:23)]) 
data_ne_prep <- data_ne_recipe %>%
  prep() %>%
  juice()

rf_ne_mod <- rand_forest() %>%
  set_args(mtry = 3, trees = 100, min_n = 800) %>%
  set_engine("ranger", importance = "permutation") %>%
  set_mode("classification")

rf_ne_workflow <- workflow() %>%
  add_recipe(data_ne_recipe) %>%
  add_model(rf_ne_mod)

rf_ne_fit <- rf_ne_workflow %>%
  last_fit(data_ne_split)

test_ne_performance <- rf_ne_fit %>% collect_metrics()
test_ne_performance

test_ne_predictions <- rf_ne_fit %>% collect_predictions()

test_ne_boyce <- Boyce(obs = test_ne_predictions$pres, pred = test_ne_predictions$.pred_1,
                       main = str_glue("Non-expl pts_minn_{param_a_final$min_n}"))

# Final model
final_ne_model <- fit(rf_ne_workflow, data_ne)

# Variable importance
final_ne_model %>%
  extract_fit_parsnip() %>%
  vip(num_features = 17)

ranger_ne_obj <- extract_fit_parsnip(final_ne_model)$fit
ranger_ne_obj$variable.importance


# ### Saving models
# save(final_a_model,  file = str_glue("{path}../model_global/bg_rand_prop_2hr_thin_mtry3_tr100_minn_800_final/tidy_global_all.RData"))
# save(final_e_model,  file = str_glue("{path}../model_global/bg_rand_prop_2hr_thin_mtry3_tr100_minn_800_final/tidy_global_expl.RData"))
# save(final_ne_model, file = str_glue("{path}../model_global/bg_rand_prop_2hr_thin_mtry3_tr100_minn_800_final/tidy_global_ne.RData"))