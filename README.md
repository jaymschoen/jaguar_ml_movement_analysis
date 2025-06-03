**Code for behavior-specific machine learning jaguar movement analysis, published [here](https://doi.org/10.1016/j.biocon.2025.110978).**

First two R scripts and Google Earth Engine scripts shared here; remaining code will be committed shortly.

Note on order of operations:
- R script 01 ingests and pepares movement data
- R script 02 builds Hidden Markov Models
- R script 03 prepares areas for background sampling
  - background areas used to generate covariate data in Google Earth Engine
- GEE scripts (in no order) use an exported vector file to extract data for focal areas
- R script 04 builds, aligns, and exports covariate data
- R script 05 imports final covariate data
- R script 06 scales covariate data by polygon
- R script 07 extracts covariate data for presence and background points
- R script 08 thins data to balance sampling between individuals
- R script 09 visualizes data prior to modeling
- R script 10 builds RandomForest models
- R script 11 creates variable importance and partial dependence plots
- R script 12 generates rasterized predictions for the models