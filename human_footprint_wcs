/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var human_impact_index = ee.ImageCollection("projects/HII/v1/hii"),
    infrastructure = ee.ImageCollection("projects/HII/v1/driver/infrastructure"),
    land_use = ee.ImageCollection("projects/HII/v1/driver/land_use"),
    population_density = ee.ImageCollection("projects/HII/v1/driver/population_density"),
    power = ee.ImageCollection("projects/HII/v1/driver/power"),
    railways = ee.ImageCollection("projects/HII/v1/driver/railways"),
    roads = ee.ImageCollection("projects/HII/v1/driver/roads"),
    water = ee.ImageCollection("projects/HII/v1/driver/water"),
    polys = ee.FeatureCollection("users/jaymschoen/jags_reg_int_ind_polygons_10k_buf");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
// Asset specifications for Human Impact Index and drivers
// Use Inspector tool to inspect ImageCollection details and pixel series charts

Map.setOptions('SATELLITE');

// // Driver Visualization
var hiiviz = {min: 5, max: 5000, palette: ["224f1a","a3ff76","feff6f","a09568","ffa802","f7797c","fb0102","d87136","a90086","7a1ca5","421137","000000"]};
// var infviz = {min: 0, max: 1000, palette: ['f1edec', 'dfbcb0', 'd08b73', 'c0583b', 'a62225', '730e27', '3c0912']};
// var luviz = {min: 0, max: 1000, palette: ["172313", "144b2a", "187328", "5f920c", "aaac20", "e1cd73", "fffdcd"]};
// var popviz = {min: 0, max: 1000, palette: ['151d44', '156c72', '7eb390', 'fdf5f4', 'db8d77', '9c3060', '340d35']};
// var powerviz = {min: 0, max: 1000, palette: ["000000", "800026", "bd0026", "e31a1c", "fc4e2a", "fd8d3c", "feb24c", "fed976", "ffeda0", "ffffcc"]};
// var railviz = {min: 0, max: 1000, palette: ['feedb0', 'f7b37c', 'eb7858', 'ce4356', '9f2462', '66185c', '2f0f3e']};
// var roadviz = {min: 0, max: 1000, palette: ['feedb0', 'f7b37c', 'eb7858', 'ce4356', '9f2462', '66185c', '2f0f3e']};
// var waterviz = {min: 0, max: 1000, palette: ['e6f1f1', 'a2cee2', '76a4e5', '7871d5', '7642a5', '621d62', '360e24']};

// Map.addLayer(infrastructure, infviz, "Infrastructure driver", false);
// Map.addLayer(land_use, luviz, "Land cover driver", false);
// Map.addLayer(population_density, popviz, "Population driver", false);
// Map.addLayer(power, powerviz, "Power driver", false);
// Map.addLayer(railways, railviz, "Railways driver", false);
// Map.addLayer(roads, roadviz, "Roads driver", false);
// Map.addLayer(water, waterviz, "Navigable water driver", false);
// Map.addLayer(human_impact_index, hiiviz, "Human Impact Index", false);

var roi = polys.geometry();
var years = ee.List.sequence(2000,2016,1);

// Selecting images for focal countries/areas
var hii = human_impact_index
      .filterDate('2000-01-01', '2016-12-31')
      .map(function(x){return x.clip(roi)});

var hii_l = hii.toList(16);
print(hii);
// ------- NOTE: images start in 2001; will need to use 2001 HII image for 2000

Map.addLayer(hii.first(), hiiviz, "Human Impact Index 2001");

// Adding layers with forloop

// Get the size of the image list (this is a server side ee.Number object).
var len = hii_l.size();

len.evaluate(function(l) {
  for (var i=0; i < l; i++) {
    var img = ee.Image(hii_l.get(i));
    var yr = img.date().get('year').getInfo();
    Map.addLayer(img, hiiviz, 'Human Footprint '+yr, false);
  } 
});



//  Batch Download
var longlat = "EPSG:4326";

var batch = require('users/fitoprincipe/geetools:batch');

// batch.Download.ImageCollection.toDrive(hii, 'human_footprint', 
//                 {scale: 250, 
//                 region: polys, 
//                 type: 'float',
//                 crs: longlat});
