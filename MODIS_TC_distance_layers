/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var geometry = 
    /* color: #d63000 */
    /* shown: false */
    /* displayProperties: [
      {
        "type": "rectangle"
      }
    ] */
    ee.Geometry.Polygon(
        [[[-60.21877837636518, -22.31334564969466],
          [-60.21877837636518, -22.792116176569426],
          [-59.431882624412054, -22.792116176569426],
          [-59.431882624412054, -22.31334564969466]]], null, false),
    polys = ee.FeatureCollection("users/jaymschoen/jags_reg_int_ind_polygons_10k_buf");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
//// ===== Distance to tree cover (TC) and non-tree cover (nonTC) images ===== ////

// See "Functions" section at end of script for functions used

/**/ // Setup
Map.setCenter(-54.888, -19.168, 6);
Map.setOptions('SATELLITE');

var roi = polys.geometry();

// MODIS Tree Cover Images 2000-2016
var tc_coll = ee.ImageCollection('MODIS/006/MOD44B')
        .select('Percent_Tree_Cover')
        .filterDate('2000-01-01', '2016-12-31')
        .map(function(image){return image.clip(roi)});

/**/ // Splitting TC into Percentiles

// Visually inspecting percentile cutoffs for "high" and "low" tree cover

// Apply pctl_cut function for precentile cutoffs
var treeCover_p10 = tc_coll.map(function(x){return(pctl_cut(x, "10"))});
var treeCover_p20 = tc_coll.map(function(x){return(pctl_cut(x, "20"))});
var treeCover_p30 = tc_coll.map(function(x){return(pctl_cut(x, "30"))});
var treeCover_p40 = tc_coll.map(function(x){return(pctl_cut(x, "40"))});
var treeCover_p50 = tc_coll.map(function(x){return(pctl_cut(x, "50"))});

// // Visualizing 
Map.addLayer(treeCover_p10.first(), {}, 'tree cover p10 2000');
Map.addLayer(treeCover_p20.first(), {}, 'tree cover p20 2000');
Map.addLayer(treeCover_p40.first(), {}, 'tree cover p40 2000');

// Using 40th %ile for high TC based on visual analysis of results (see paper & supplement for details)
// Using 20th %ile for low TC based on visual analysis of results (see paper & supplement for details)


/**/ // Sieving 

// Filtering out (sieving) small TC/nonTC patches (<2.5 sq km or 40 pixels) within target patches

// Sieving TC images with sieve() function
var TCp40_sieved40 = treeCover_p40.map(function(x){return sieve(x, 40)});
var TCp20_sieved40 = treeCover_p20.map(function(x){return sieve(x, 40)});

// Sieving nonTC images 
var nonTCp40_sieved40 = treeCover_p40.map(function(x){return sieve(ee.Image(x).not(), 40)});  
var nonTCp20_sieved40 = treeCover_p20.map(function(x){return sieve(ee.Image(x).not(), 40)});

// Visual tests (reprojecting to force 250m scale)
Map.addLayer(ee.Image(treeCover_p40.first().reproject({crs:'epsg:4326', scale: 250})), {}, 'p40 unsieved');
Map.addLayer(ee.Image(TCp40_sieved40.first().reproject({crs:'epsg:4326', scale: 250})), {}, 'p40 sieve test');

Map.addLayer(ee.Image(treeCover_p20.first().reproject({crs:'epsg:4326', scale: 250})), {}, 'p20 unsieved');
Map.addLayer(ee.Image(TCp20_sieved40.first().reproject({crs:'epsg:4326', scale: 250})), {}, 'p20 sieve test');

// sieve() function working as expected


/**/ // Distance rasters 
var imageNumbers = ee.List.sequence(0, treeCover_p40.size().subtract(1));

// Distance to TC/nonTC images using sieve40 (removing patches <10 km2 within target area [TC or nonTC])
var dist_TCp40 = imageNumbers.map(function(x){
  return distance(TCp40_sieved40, x);
});
var dist_TCp20 = imageNumbers.map(function(x){
  return distance(TCp20_sieved40, x);
});
var dist_nonTCp40 = imageNumbers.map(function(x){
  return distance(nonTCp40_sieved40, x);
});
var dist_nonTCp20 = imageNumbers.map(function(x){
  return distance(nonTCp20_sieved40, x);
});

/**/ // Adding layers with forloop

// Change percentile (p20/p40) and tc/nonTC to visualize desired collection

// Get the size of the image list to use evaluate() (client-side of map() to allow forloop)
var len = imageNumbers.size();

len.evaluate(function(l) {
  for (var i=0; i < l; i++) {
    var img = ee.Image(nonTCp40_sieved40.toList(nonTCp40_sieved40.size()).get(i));
    var img1 = ee.Image(dist_nonTCp40.get(i));
    var yr_str = ee.String(img.get('system:index')).slice(0,4).getInfo();
    Map.addLayer(img.reproject({crs:'epsg:4326', scale: 250}), 
                  {}, 'nonTCp40_sieved40_'+yr_str, false);
    Map.addLayer(img1.reproject({crs:'epsg:4326', scale: 250}),
                  {max:3000}, 'nonTCp40_dist_'+yr_str, false);
  } 
});

/**/ // Batch Download
var toExport = ee.ImageCollection.fromImages(dist_TCp20); // change collection for p20/p40 and dist_TC/dist_nonTC
var utm21s = "EPSG:32721";
var longlat = "EPSG:4326";

var batch = require('users/fitoprincipe/geetools:batch');

// batch.Download.ImageCollection.toDrive(toExport, 'distance_tree_cover_p20', 
//                 {scale: 250, 
//                  region: polys, 
//                  type: 'float',
//                  crs: longlat});

/**/ // Functions

function pctl_cut(image, thresh) {
  var quant10_list = ee.List.sequence(10,100,10);
  var cuts = image.reduceRegion({reducer: ee.Reducer.percentile(quant10_list), 
                                 maxPixels:1e10});
  var threshold = ee.Number(cuts.get("Percent_Tree_Cover_p"+thresh));
  
  return image.gte(threshold);
}

function sieve(img, nPix) {
  var con = ee.Image(img).connectedPixelCount(100, true); // 8 neighbor connected
  var mask = con.gt(nPix);
  var sieved = ee.Image(img).updateMask(mask).unmask();
  
  return(sieved);
}

function distance(coll, y) {
  var list = coll.toList(imageNumbers.size()); // converting to list to call by indexed number
  var tar = ee.Image(list.get(y));             // target year image
  var yr = ee.String(tar.get('system:index')).slice(0,4);
  var non_tc = tar.lte(0);   // creating "non TC" mask
  var tc = tar.gt(0);        // creating "TC" mask
  var dist = non_tc.cumulativeCost({           // distance via cost function
                      source: tc,
                      maxDistance: 50000
    });
  
  return dist
    .set('year', yr);
}
