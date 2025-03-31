/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var pop = ee.ImageCollection("WorldPop/GP/100m/pop"),
    polys = ee.FeatureCollection("users/jaymschoen/jags_reg_int_ind_polygons_10k_buf");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
//// ===== Population Density ===== ////

var roi = polys.geometry();
var years = ee.List.sequence(2000,2016,1);

// Selecting images for focal countries/areas
var popFocal = pop.filter(
  ee.Filter.inList('country', ee.List(['PRY','BRA','BOL','ARG'])))
  .map(function(x){return x.clip(polys)});

print(popFocal);

Map.addLayer(popFocal, {min:0, max:0.1}, 'population');

// Mosaicking countries by year

function annMos (coll, yr) {
  return coll.filter(ee.Filter.eq('year', yr)).mosaic()
    .set('id', yr);
}

// test
var annMos2007 = annMos(popFocal, 2007);
print(annMos2007);
Map.addLayer(annMos2007, {}, '2007 Mosaic');

// Apply to IC
var annMosaics = ee.ImageCollection(years.map(function(x){
  return annMos(popFocal, x);
}));

print('Annual Mosaics', annMosaics);

// Batch Download
var lonlat = "EPSG:4326";

var batch = require('users/fitoprincipe/geetools:batch');

// batch.Download.ImageCollection.toDrive(annMosaics, 'dissertation_columbia//ch3_covariates', 
//                 {scale: 250, 
//                 region: polys, 
//                 type: 'float',
//                 crs: lonlat});






