**Code for behavior-specific machine learning jaguar movement analysis, published [here](https://doi.org/10.1016/j.biocon.2025.110978).**

First two R scripts and Google Earth Engine scripts shared here; remaining code will be committed shortly.

Note on order of operations:
    - R script 01 ingests and pepares movement data 
    - R script 02 builds Hidden Markov Models
    - R script 03 prepares areas for background sampling 
        - background areas use to generate covariate data in Google Earth Engine
    - GEE scripts (in no order) use an exported vector file to extract data for focal areas
    - R script 04 ingests and aligns covariate data