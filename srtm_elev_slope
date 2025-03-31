/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var polys = ee.FeatureCollection("users/jaymschoen/jags_reg_int_ind_polygons_10k_buf"),
    srtm = ee.Image("CGIAR/SRTM90_V4");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
//// ===== Elevation and slope images ===== ////

var roi = polys.geometry();

var elevation = srtm.clip(roi).select('elevation');
var slope = ee.Terrain.slope(elevation);
Map.addLayer(slope, {min: 0, max: 30}, 'slope');
Map.addLayer(elevation, {min: 0, max: 500}, 'elevation');

var longlat = "EPSG:4326";

// Export.image.toDrive({
//   image: elevation,
//   description: 'elevation',
//   folder: 'elevation_slope',
//   region: roi,
//   scale: 250,
//   crs: longlat,
//   maxPixels:1e12});

// Export.image.toDrive({
//   image: slope,
//   description: 'slope',
//   folder: 'elevation_slope',
//   region: roi,
//   scale: 250,
//   crs: longlat,
//   maxPixels:1e12});