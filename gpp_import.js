/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var polys = ee.FeatureCollection("users/jaymschoen/jags_reg_int_ind_polygons_10k_buf");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
// ===== Median annual GPP composites ===== //
Map.setOptions('SATELLITE');

// Original dataset is 8 day cumualtive composite of daily product (MODIS)

var roi = polys.geometry();
var years = ee.List.sequence(2000,2016,1);

var medGpp = ee.ImageCollection.fromImages(years.map(function(y) {
  var filtered = ee.ImageCollection('MODIS/006/MOD17A2H')
    .select('Gpp')
    .filter(ee.Filter.calendarRange({
      start: y,
      field: 'year'}));
  var composites = filtered
    .map(function(image){return image.clip(roi)})
    .median();
  return composites
    .set('year', y);
}));

print(medGpp); 

var gppVis = {
  min: 0.0,
  max: 600.0,
  palette: ['bbe029', '0a9501', '074b03'],
};

// Testing
Map.addLayer(medGpp.first(), gppVis, 'Median GPP 2000');


// Adding layers with forloop

// Get the size of the image list to use evaluate() (client-side of map() to allow forloop)

var medGpp_list = medGpp.toList(17);
var len = medGpp_list.size();

len.evaluate(function(l) {
  for (var i=0; i < l; i++) {
    var img = ee.Image(medGpp_list.get(i));
    var yr = img.get('year').getInfo();
    Map.addLayer(img, gppVis, 'Median GPP '+yr, false);
  } 
});



// Batch Download
var longlat = "EPSG:4326";

var batch = require('users/fitoprincipe/geetools:batch');

// batch.Download.ImageCollection.toDrive(medGpp, 'gpp/', 
//                 {scale: 250, 
//                 region: polys, 
//                 type: 'float',
//                 crs: longlat});






