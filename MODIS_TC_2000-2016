/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var pgy = ee.FeatureCollection("users/jaymschoen/pgy_adm1_utm21s"),
    polys = ee.FeatureCollection("users/jaymschoen/jags_reg_int_ind_polygons_10k_buf");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
// ===== Clipping/exporting MODIS Tree Cover 2000-2016 to target polygons ===== //

Map.setCenter(-55.4597, -24.0721, 6);
Map.setOptions('SATELLITE');

var roi = polys.geometry();

// --- MODIS Tree Cover Images 2000-2016 --- //
var tc_coll = ee.ImageCollection('MODIS/006/MOD44B')
        .select('Percent_Tree_Cover')
        .filterDate('2000-01-01', '2016-12-31')
        .map(function(image){return image.clip(roi)});

// print('tree cover collection', tc_coll);


// Adding layers to map with forloop
var tcList = tc_coll.toList(17);
var len = tcList.size();

len.evaluate(function(l) {
  for (var i=0; i < l; i++) {
    var img = ee.Image(tcList.get(i));
    var yr_str = ee.String(img.get('system:index')).slice(0,4).getInfo();
    Map.addLayer(img, {palette: ['lightgrey', 'green', 'darkgreen']}, 'Tree Cover '+yr_str.toString(), false);
  } 
});


var utm21s = "EPSG:32721";
var longlat = "EPSG:4326";

var batch = require('users/fitoprincipe/geetools:batch');

// batch.Download.ImageCollection.toDrive(tc_coll, 'percent_tree_cover', 
//                 {scale: 250, 
//                 region: polys, 
//                 type: 'float',
//                 crs: longlat});