04-covariate_prep
================
Jay Schoen

Jaguar Machine Learning Movement Model

04 - Covariate Prep

Creating several covariate data layers here; others were created in
Google Earth Engine and are imported and prepared for analysis.

To save and access future steps, create directories (folders) within
“data” folder for each covariate (e.g., “gpp”, “population”, etc.)

# Setup

``` r
library(sp)
```

    ## The legacy packages maptools, rgdal, and rgeos, underpinning the sp package,
    ## which was just loaded, will retire in October 2023.
    ## Please refer to R-spatial evolution reports for details, especially
    ## https://r-spatial.org/r/2023/05/15/evolution4.html.
    ## It may be desirable to make the sf package available;
    ## package maintainers should consider adding sf to Suggests:.
    ## The sp package is now running under evolution status 2
    ##      (status 2 uses the sf package in place of rgdal)

``` r
library(sf)
```

    ## Warning: package 'sf' was built under R version 4.3.3

    ## Linking to GEOS 3.11.2, GDAL 3.8.2, PROJ 9.3.1; sf_use_s2() is TRUE

``` r
library(terra)
```

    ## terra 1.7.39

``` r
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.2     ✔ readr     2.1.4
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.0
    ## ✔ ggplot2   3.4.2     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.2     ✔ tidyr     1.3.0
    ## ✔ purrr     1.0.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ tidyr::extract() masks terra::extract()
    ## ✖ dplyr::filter()  masks stats::filter()
    ## ✖ dplyr::lag()     masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(rgeos)
```

    ## rgeos version: 0.6-4, (SVN revision 699)
    ##  GEOS runtime version: 3.11.2-CAPI-1.17.2 
    ##  Please note that rgeos will be retired during October 2023,
    ## plan transition to sf or terra functions using GEOS at your earliest convenience.
    ## See https://r-spatial.org/r/2023/05/15/evolution4.html for details.
    ##  GEOS using OverlayNG
    ##  Linking to sp version: 2.0-0 
    ##  Polygon checking: TRUE 
    ## 
    ## 
    ## Attaching package: 'rgeos'
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     symdiff

``` r
long_lat <- "epsg:4326"
merc <- "epsg:3395"
utm21s <- "epsg:32721"

# Set to directory with cloned repo
setwd("~/Columbia E3B/PhD/Dissertation/Ch3-jaguar_movement/scripts/git_push_R/")

# Make this path where you stored environmental covariate data
path <- "C:/Users/jayms/Documents/Columbia E3B/PhD/Dissertation/Ch3-jaguar_movement/covariate_data/"


# Importing polygon shapefile for regularly sampled jaguars (buffered by 10km)
polys <- vect(str_glue("data/jags_reg_int_ind_polygons_10k_buf.shp")) %>%
  disagg() %>%
  # project(merc) %>%
  project(utm21s)     # projecting for distance calculations


# 60km buffer around individual polygons for distance radius
polys_60k_buf <- buffer(polys, 5e4) %>%   # "polys" already buffered by 10km
  aggregate() %>%
  disagg()
plot(polys_60k_buf)

## Stock raster with all polygon areas set to 0 (for aligning and filling NAs in polys)
stock_r <- rast(str_glue("{path}gpp/gpp_2000.tif")) %>%
  `names<-`("stock_raster")
stock_r[!is.na(stock_r)] <- 0 # clearing all values but leaving NAs for masking
```

    ## |---------|---------|---------|---------|=========================================                                          |---------|---------|---------|---------|=========================================                                          |---------|---------|---------|---------|=========================================                                          

``` r
# writeVector(polys_60k_buf, "data/polys_60k_buf.shp", overwrite = T)

plot(polys_60k_buf); plot(polys, col = "blue", add = T)
```

![](04-covariate_prep_files/figure-gfm/Setup-1.png)<!-- -->

# Static Covariate Data

## Population Density

``` r
# Setting water masked areas to 0
pop_dens <- list.files(str_glue("{path}population/raw/"),
                  pattern = ".tif",
                  all.files = TRUE, 
                  full.names = TRUE) %>%
  rast() %>%
  `names<-`(str_glue("pop_dens_{2000:2016}"))
  
pop_dens <- cover(pop_dens[[1:17]], stock_r)
```

    ## |---------|---------|---------|---------|=========================================                                          

``` r
plot(pop_dens[[1]])
```

![](04-covariate_prep_files/figure-gfm/Population%20Density-1.png)<!-- -->

``` r
# Exporting
# for (i in 1:nlyr(pop_dens)){
#   export <- pop_dens[[i]]
#   filename <- str_glue("data/population/{names(pop_dens[[i]])}.tif")
#   writeRaster(export, filename)
# }
```

## Roads

``` r
# Roads data from https://data.humdata.org/ (OSM data)
## cropped to polys_60k_buf shape in QGIS then imported here for analysis

roads_polys_60k_buf <- vect(str_glue("{path}roads/roads_rgn_clip_60k_buf.shp")) %>%
  project(utm21s)

## Subsetting with major/paved roads
### Note: most of these "major" roads are unpaved 1-2 lane roads, likely use by jaguars when moving
roads_mjr <- c("trunk", "trunk_link", "motorway", "motorway_link", "primary",
               "primary_link", "secondary", "secondary_link", "tertiary", "tertiary_link")

roads_mjr_polys_60k_buf <- subset(roads_polys_60k_buf, 
                                  roads_polys_60k_buf$highway %in% roads_mjr)
roads_min_polys_60k_buf <- subset(roads_polys_60k_buf, 
                                  !(roads_polys_60k_buf$highway %in% roads_mjr))
```

### Distance to roads (“major” subset)

``` r
# Creating blank rasters with extent of each polygon
for(i in 1:length(polys_60k_buf)) {
  assign(str_glue("poly{i}_r"), rast(polys_60k_buf[i], resolution = 250, vals = 1e5))
}

# Creating distance rasters for each polygon extent
polys_60k_dist_roads_mjr <- lapply(1:length(polys_60k_buf), function(x)
  distance(get(str_glue("poly{x}_r")),
           crop(roads_mjr_polys_60k_buf, get(str_glue("poly{x}_r")))))
```

    ## |---------|---------|---------|---------|=========================================                                          |---------|---------|---------|---------|=========================================                                          |---------|---------|---------|---------|=========================================                                          |---------|---------|---------|---------|=========================================                                          |---------|---------|---------|---------|=========================================                                          |---------|---------|---------|---------|=========================================                                          |---------|---------|---------|---------|=========================================                                          |---------|---------|---------|---------|=========================================                                          |---------|---------|---------|---------|=========================================                                          |---------|---------|---------|---------|=========================================                                          

``` r
# plot(polys_60k_dist_roads_mjr[[3]])

# Clipping to original shape to avoid overlap before merge
polys_60k_dist_roads_mjr_list <- lapply(polys_60k_dist_roads_mjr, 
                                        function(x){crop(x, polys_60k_buf, mask = T)})
                                               # mask(x, polys_60k_buf)})

# Merging into one large raster
polys_60k_dist_roads_mjr <- do.call(merge, polys_60k_dist_roads_mjr_list) %>%
  crop(polys, mask = T) #%>%
  # mask(polys)

# plot(polys_60k_dist_roads_mjr)

## Reproject/align
dist_roads_mjr <- polys_60k_dist_roads_mjr %>%
  `names<-`("dist_roads_mjr") %>%
  project(stock_r)    # re-aligning extent
```

    ## |---------|---------|---------|---------|=========================================                                          

``` r
# Exporting
writeRaster(polys_60k_dist_roads_mjr, str_glue("data/roads/dist_roads_mjr.tif"), overwrite = T)
```

### Road Density (minor roads)

``` r
# Minor road density in poly rasters created above (dist_roads)
## Modified from Robert Hijmans' response to https://stackoverflow.com/questions/69035297/calculating-road-density-raster-from-road-shapefile

# Building function to apply to all polygons
dens_roads <- function(n, l) {
  p <- as.polygons(get(str_glue("poly{n}_r")))
  roads <- l %>% crop(p)
  
  rs <- rast(get(str_glue("poly{n}_r"))) %>%
    `res<-`(250) %>%                           # 250m resolution
    `values<-`(1:ncell(.)) %>%
    `names<-`("rast")
  rsp <- as.polygons(rs)
  
  rp <- intersect(roads, rsp)
  rp$length <- perim(rp) / 1000                # length of roads/1 km perimeter
  x <- tapply(rp$length, rp$rast, sum)
  
  r <- rast(rs)
  r[as.integer(names(x))] <- as.vector(x)
  
  return(r)
}

# Minor Roads
## Applying it to all polygons
polys_dens_roads_min_list <- lapply(1:length(polys_60k_buf),
                                    function(x) dens_roads(x, roads_min_polys_60k_buf))

polys_dens_roads_min <- do.call(merge, polys_dens_roads_min_list)
```

    ## |---------|---------|---------|---------|=========================================                                          

``` r
plot(polys_dens_roads_min)
```

![](04-covariate_prep_files/figure-gfm/Road%20Density%20(minor%20roads)-1.png)<!-- -->

``` r
## Individual polygons
for(i in 1:length(polys_60k_buf)) {
  assign(str_glue("poly{i}_dens_roads_all"), 
         polys_dens_roads_min_list[[i]])  
}

## Reproject/align
# dens_roads_min <- polys_dens_roads_min %>%
#   project(stock_r) %>%    
#   cover(stock_r) %>%               # setting background in polygons to 0
#   crop(polys %>% project(long_lat), mask = T) %>%
#   project(stock_r)                 # re-aligning extent
dens_roads_min <- polys_dens_roads_min %>%
  crop(polys, mask = T) %>%
  project(stock_r) %>%    
  cover(stock_r)                # setting background in polygons to 0
```

    ## |---------|---------|---------|---------|=========================================                                          

``` r
# Exporting
writeRaster(polys_dens_roads_min, "data/roads/dens_roads_min.tif", overwrite = T)
```

# Importing Annual Covariate Data

## Primary production

### GPP

``` r
# Importing from GEE exports
gpp <- list.files(str_glue("{path}gpp/raw/"),
                  pattern = ".tif$",
                  all.files = TRUE, 
                  full.names = TRUE) %>%
  rast() %>%
  `names<-`(str_glue("gpp_{2000:2016}"))

## setting NA areas in polygons (water mask) to 0
gpp <- cover(gpp[[1:17]], stock_r)
```

    ## |---------|---------|---------|---------|=========================================                                          

``` r
# # Exporting
# for (i in 1:nlyr(gpp)){
#   export <- gpp[[i]]
#   filename <- str_glue("data/gpp/{names(gpp[[i]])}.tif")
#   writeRaster(export, filename)
# }
```

## Tree Cover

### Percent Tree Cover

``` r
pct_tc <- list.files(str_glue("{path}percent_tree_cover/raw"), 
                     pattern = ".tif$",
                     all.files = TRUE, 
                     full.names = TRUE) %>%
  rast() %>%
  `names<-`(str_glue("pct_tc_{2000:2016}"))

# Filling water masked areas with 0
pct_tc <- cover(pct_tc[[1:17]], stock_r)
```

    ## |---------|---------|---------|---------|=========================================                                          

``` r
# Exporting
# for (i in 1:nlyr(pct_tc)){
#   export <- pct_tc[[i]]
#   filename <- str_glue("data/percent_tree_cover/{names(pct_tc[[i]])}.tif")
#   writeRaster(export, filename)
# }
```

### Distance to Non-Tree Cover (distance to forest edge)

*These are imported at 150m resolution then resampled here bc GEE’s
resampling masked out small non-tc areas within high tree cover during
export*

#### \<20th %ile tree cover

``` r
dist_non_tc_p20 <- list.files(str_glue("{path}distance_non_tree_cover_p20/raw/res_150m/"),
                          pattern = ".tif$",
                          all.files = TRUE, 
                          full.names = TRUE) %>%
  rast() %>%
  project(stock_r) %>%      
  `names<-`(str_glue("dist_non_tc_p20_{2000:2016}"))
```

    ## |---------|---------|---------|---------|=========================================                                          

``` r
# Some polygons in annual images show up NA bc no non-TC areas in 10k buffered extent. Assigning max distance value of respective annual image to those polygons

for (i in 1:nlyr(dist_non_tc_p20)) {
  img <- dist_non_tc_p20[[i]]         # pulling annual image
  max <- stock_r                      # stock_r has values within all polygon areas 
  max[!is.na(max)] <- minmax(img)[2]  # setting all areas within polygons to max
  f <- cover(img, max)                # filling NA areas of annual image with max value
  dist_non_tc_p20[[i]] <- f           # appending raster stack
}
```

    ## |---------|---------|---------|---------|=========================================                                          |---------|---------|---------|---------|=========================================                                          |---------|---------|---------|---------|=========================================                                          |---------|---------|---------|---------|=========================================                                          |---------|---------|---------|---------|=========================================                                          |---------|---------|---------|---------|=========================================                                          |---------|---------|---------|---------|=========================================                                          |---------|---------|---------|---------|=========================================                                          |---------|---------|---------|---------|=========================================                                          |---------|---------|---------|---------|=========================================                                          |---------|---------|---------|---------|=========================================                                          |---------|---------|---------|---------|=========================================                                          |---------|---------|---------|---------|=========================================                                          |---------|---------|---------|---------|=========================================                                          |---------|---------|---------|---------|=========================================                                          |---------|---------|---------|---------|=========================================                                          |---------|---------|---------|---------|=========================================                                          |---------|---------|---------|---------|=========================================                                          |---------|---------|---------|---------|=========================================                                          |---------|---------|---------|---------|=========================================                                          |---------|---------|---------|---------|=========================================                                          |---------|---------|---------|---------|=========================================                                          |---------|---------|---------|---------|=========================================                                          |---------|---------|---------|---------|=========================================                                          |---------|---------|---------|---------|=========================================                                          |---------|---------|---------|---------|=========================================                                          |---------|---------|---------|---------|=========================================                                          

``` r
# Exporting
# for (i in 1:nlyr(dist_non_tc_p20)){
#   export <- dist_non_tc_p20[[i]]
#   filename <- str_glue("data/distance_non_tree_cover_p20/{names(dist_non_tc_p20[[i]])}.tif")
#   writeRaster(export, filename, overwrite = TRUE)
# }
```

#### \<40th %ile tree cover

``` r
dist_non_tc_p40 <- list.files(str_glue("{path}distance_non_tree_cover_p40/raw/res_150m/"),
                          pattern = ".tif$",
                          all.files = TRUE, 
                          full.names = TRUE) %>%
  rast() %>%
  project(stock_r) %>%      
  `names<-`(str_glue("dist_non_tc_p40_{2000:2016}"))
```

    ## |---------|---------|---------|---------|=========================================                                          

``` r
# Exporting
# for (i in 1:nlyr(dist_non_tc_p40)){
#   export <- dist_non_tc_p40[[i]]
#   filename <- str_glue("data/distance_non_tree_cover_p40/{names(dist_non_tc_p40[[i]])}.tif")
#   writeRaster(export, filename)
# }
```

## Land Cover

``` r
# Importing Earth Engine exports (in several regions)
## subsetting to 2000-2016
mb_chaco <- rast(str_glue("{path}landcover/mb_chaco_polys.tif")) %>%
  subset(1:17)

mb_af <- rast(str_glue("{path}landcover/mb_af_polys.tif")) %>%
  subset(1:17)

mb_bra1 <- rast(str_glue("{path}landcover/mb_bra1_polys.tif")) %>%
  subset(16:32)

mb_bra2 <- rast(str_glue("{path}landcover/mb_bra2_polys.tif")) %>%
  subset(16:32)

mb_bra3 <- rast(str_glue("{path}landcover/mb_bra3_polys.tif")) %>%
  subset(16:32)

# Merging entire region by year
merge_yrs <- function(x) {
  merge(mb_bra1[[x]], 
        mb_bra2[[x]], 
        mb_bra3[[x]], 
        mb_chaco[[x]], 
        mb_af[[x]])
}

mb_rgn <- lapply(1:17, merge_yrs) %>%       
  # sds()   # to SpatRasterDataset (brick/stack analog)
  rast()  
```

    ## |---------|---------|---------|---------|=====================                                          |---------|---------|---------|---------|==================                                          |---------|---------|---------|---------|==================                                          |---------|---------|---------|---------|======================                                          |---------|---------|---------|---------|===================                                          |---------|---------|---------|---------|=====================                                          |---------|---------|---------|---------|==================                                          |---------|---------|---------|---------|=====================                                          |---------|---------|---------|---------|==================                                          |---------|---------|---------|---------|==================                                          |---------|---------|---------|---------|=====================                                          |---------|---------|---------|---------|==================                                          |---------|---------|---------|---------|=====================                                          |---------|---------|---------|---------|==================                                          |---------|---------|---------|---------|=====================                                          |---------|---------|---------|---------|==================                                          |---------|---------|---------|---------|==================                                          

``` r
# test
plot(mb_rgn[[17]])
```

![](04-covariate_prep_files/figure-gfm/Land%20Cover-1.png)<!-- -->

``` r
## Consolidating land cover classes (see MapBiomas legends for details):

# classification matrix for classify function in terra
m <- c(1,1,3,1,4,1,5,1,49,1,                                 # 1: Forested
       10,2,11,2,12,2,32,2,29,2,13,2,43,2,                   # 2: Natural non-forested
       15,3,                                                 # 3: Pasture
       18,4,19,4,39,4,20,4,40,4,41,4,36,4,46,4,47,4,48,4,    # 4: Agriculture
       21,5,                                                 # 5: Ag/pasture mixed
       9,6,                                                  # 6: Forest Plantation
       22,7,23,7,24,7,30,7,25,7,                             # 7: Non-vegetated
       26,8,33,8,31,8,                                       # 8: Water
       6,9)                                                  # 9: Floodplain
rcl_m <- matrix(m, ncol=2, byrow=TRUE)

mb_rgn_cons <- classify(mb_rgn, rcl_m) 
```

    ## |---------|---------|---------|---------|=========================================                                          

    ## Warning: [classify] Estimated disk space needed without compression: 41GB.
    ## Available: 24 GB.

``` r
mb_rgn_cons <- mb_rgn_cons %>%
  crop(polys%>%project(long_lat), mask = T)
```

    ## |---------|---------|---------|---------|=========================================                                          |---------|---------|---------|---------|=========================================                                          

    ## Warning: [crop] Estimated disk space needed without compression: 41GB.
    ## Available: 28 GB.

``` r
# test
plot(mb_rgn_cons[[1]])
```

![](04-covariate_prep_files/figure-gfm/Land%20Cover-2.png)<!-- -->

``` r
# Exporting annual consolidated rasters
# for (i in 1:length(mb_rgn_cons)){
#   export <- mb_rgn_cons[i]
#   filename <- str_glue("data/landcover/mb_lc_cons_{names(mb_rgn_cons[i])}.tif")
#   writeRaster(export, filename, overwrite = T)
# }
```

### Distance to natural areas

``` r
dist_nat <- list.files(str_glue("{path}distance_natural/raw/"),
                       pattern = ".tif$",
                       all.files = TRUE, 
                       full.names = TRUE) %>%
  rast() %>%
  project(stock_r) %>%
  `names<-`(str_glue("dist_nat_{2000:2016}"))
```

    ## |---------|---------|---------|---------|=========================================                                          

``` r
# These are at 150m resolution originally then resampled here (to 250m) bc GEE's resampling masked out some small areas 

# Exporting
# for (i in 1:nlyr(dist_nat)){
#   export <- dist_nat[[i]]
#   filename <- str_glue("data/distance_natural/{names(dist_nat[[i]])}.tif")
#   writeRaster(export, filename)
# }
```

### Distance to non-natural areas

*These are at 150m resolution originally then resampled here (to 250m)
bc GEE’s resampling masked out some small areas*

``` r
dist_non_nat <- list.files(str_glue("{path}distance_non_natural/raw/"),
                           pattern = ".tif$",
                           all.files = TRUE, 
                           full.names = TRUE) %>%
  rast() %>%
  project(stock_r) %>%     
  `names<-`(str_glue("dist_non_nat_{2000:2016}"))
```

    ## |---------|---------|---------|---------|=========================================                                          

``` r
# One polygon in 2000 image (in the dry chaco) shows up NA bc no non-natural areas in 10k buffered extent. Assigning max value of 2000 image to this polygon's area

img2000 <- dist_non_nat[[1]]
plot(img2000)
```

![](04-covariate_prep_files/figure-gfm/Distance%20to%20non-natural%20areas-1.png)<!-- -->

``` r
max2000 <- stock_r
max2000[!is.na(max2000)] <- minmax(img2000)[2]
plot(max2000)
```

![](04-covariate_prep_files/figure-gfm/Distance%20to%20non-natural%20areas-2.png)<!-- -->

``` r
img2000 <- cover(img2000, max2000)

# img2000[is.na(img2000)] <- minmax(img2000)[2] 

# img2000 <- img2000 %>%            # if extents are slightly off
#   crop(polys) %>%
#   mask(polys)

plot(img2000)
```

![](04-covariate_prep_files/figure-gfm/Distance%20to%20non-natural%20areas-3.png)<!-- -->

``` r
dist_non_nat[[1]] <- img2000

# Exporting
# for (i in 1:nlyr(dist_non_nat)){
#   export <- dist_non_nat[[i]]
#   filename <- str_glue("data/distance_non_natural/{names(dist_non_nat[[i]])}.tif")
#   writeRaster(export, filename)
# }
```
