/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var mb_af = ee.Image("projects/mapbiomas_af_trinacional/public/collection1/mapbiomas_atlantic_forest_collection1_integration_v1"),
    mb_brazil = ee.Image("projects/mapbiomas-public/assets/brazil/lulc/collection9/mapbiomas_collection90_integration_v1"),
    mb_chaco = ee.Image("projects/mapbiomas-chaco/public/collection2/mapbiomas_chaco_collection2_integration_v1"),
    mb_brazil_old = ee.Image("projects/mapbiomas-workspace/public/collection6/mapbiomas_collection60_integration_v1"),
    polys = ee.FeatureCollection("users/jaymschoen/jags_reg_int_ind_polygons_10k_buf"),
    chaco = 
    /* color: #d63000 */
    /* shown: false */
    /* displayProperties: [
      {
        "type": "rectangle"
      }
    ] */
    ee.Geometry.Polygon(
        [[[-60.68969166083577, -19.76576268149378],
          [-60.68969166083577, -23.59327812422116],
          [-57.78380787177327, -23.59327812422116],
          [-57.78380787177327, -19.76576268149378]]], null, false),
    af = 
    /* color: #98ff00 */
    /* shown: false */
    /* displayProperties: [
      {
        "type": "rectangle"
      }
    ] */
    ee.Geometry.Polygon(
        [[[-54.65907082612495, -25.471310674534053],
          [-54.65907082612495, -25.916810007846568],
          [-54.0795420175312, -25.916810007846568],
          [-54.0795420175312, -25.471310674534053]]], null, false),
    bra1 = 
    /* color: #00ff00 */
    /* shown: false */
    /* displayProperties: [
      {
        "type": "rectangle"
      }
    ] */
    ee.Geometry.Polygon(
        [[[-57.76163837889257, -16.691766922922813],
          [-57.76163837889257, -20.29600197624174],
          [-52.58707783201757, -20.29600197624174],
          [-52.58707783201757, -16.691766922922813]]], null, false),
    bra2 = 
    /* color: #0000ff */
    /* shown: false */
    /* displayProperties: [
      {
        "type": "rectangle"
      }
    ] */
    ee.Geometry.Polygon(
        [[[-48.96137427501609, -8.236214050886549],
          [-48.96137427501609, -16.301938056303484],
          [-42.01801490001609, -16.301938056303484],
          [-42.01801490001609, -8.236214050886549]]], null, false),
    bra3 = 
    /* color: #999900 */
    /* shown: false */
    /* displayProperties: [
      {
        "type": "rectangle"
      }
    ] */
    ee.Geometry.Polygon(
        [[[-65.03708047678923, -2.828319528886852],
          [-65.03708047678923, -3.2287596245959356],
          [-64.64157266428923, -3.2287596245959356],
          [-64.64157266428923, -2.828319528886852]]], null, false);
/***** End of imports. If edited, may not auto-convert in the playground. *****/
var roi = polys.geometry();


Map.addLayer(roi, {}, 'polys');

var mb_af_polys = mb_af.clip(af).clip(roi);
var mb_chaco_polys = mb_chaco.clip(chaco).clip(roi);
var mb_bra1_polys = mb_brazil.clip(bra1).clip(roi);
var mb_bra2_polys = mb_brazil.clip(bra2).clip(roi);
var mb_bra3_polys = mb_brazil.clip(bra3).clip(roi);

// print(mb_chaco_polys);

Map.addLayer(mb_chaco_polys.randomVisualizer(), {}, 'MapBiomas Chaco');
Map.addLayer(mb_af_polys.randomVisualizer(), {}, 'MapBiomas Atlantic Forest');
Map.addLayer(mb_bra1_polys.randomVisualizer(), {}, 'MapBiomas Brazil1');
Map.addLayer(mb_bra2_polys.randomVisualizer(), {}, 'MapBiomas Brazil2');
Map.addLayer(mb_bra3_polys.randomVisualizer(), {}, 'MapBiomas Brazil3');

var lonlat = "epsg:4326";
// Export.image.toDrive({
//         image: mb_chaco_polys,
//         folder: 'land_cover_mb',
//         description: 'mb_chaco_polys',
//         scale: 100, 
//         maxPixels: 1e12,
//         region: chaco, 
//         crs: lonlat
// });

// Export.image.toDrive({
//         image: mb_bra1_polys,
//         folder: 'land_cover_mb',
//         description: 'mb_bra1_polys',
//         scale: 100, 
//         maxPixels: 1e12,
//         region: bra1, 
//         crs: lonlat
// });

// Export.image.toDrive({
//         image: mb_af_polys,
//         folder: 'land_cover_mb',
//         description: 'mb_af_polys',
//         scale: 100, 
//         maxPixels: 1e12,
//         region: af,
//         crs: lonlat
// });
