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
    polys = ee.FeatureCollection("users/jaymschoen/jags_reg_int_ind_polygons_10k_buf"),
    lc_coll = ee.ImageCollection("users/jaymschoen/MapBiomas_jagMD_rgn");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
// ===== Distance to natural/non-natural (matrix) land covers ===== //

/**/ //  Setup 

// IMPORTANT NOTE: FOR DISCRETE LAYERS LIKE LAND COVER, BE SURE TO SET "PYRAMID TYPE" TO "MODE" UPON INGESTION

// Map.centerObject(polys, 5);
Map.setOptions('SATELLITE'); 
var utm21s = "EPSG:32721";
var longlat = "EPSG:4326";

var roi = polys.geometry();
// Map.addLayer(roi, {}, 'polys');

// Standard long/lat projection
var proj = "EPSG:4326";


Map.addLayer(ee.Image(lc_coll.toList(17).get(16)).clip(roi).randomVisualizer(), {}, 'lc_2016');
// Map.addLayer(ee.Image(lc_coll.toList(17).get(0)).clip(roi).randomVisualizer(), {}, 'lc_2000');
// Map.addLayer(ee.Image(lc_coll.toList(17).get(14)).clip(roi).randomVisualizer(), {}, 'lc_2014');


// Creating natural/non-natural land cover mask
// Natural categories: 1,2,8,9
// Non-natural: 3-7

// Testing on sample image (2016 land cover)
var img = ee.Image(lc_coll.toList(17).get(16)).clip(roi);

var nat = img.lt(3).or(img.gt(7));             // 1-2, 8-9 natural; 3-7 non-natural
var non_nat = nat.eq(0);                       // see MapBiomas R script for consolidation details

var wat = img.eq(8);
var non_wat = wat.eq(0);

Map.addLayer(force100m(nat), {}, "natural", false);
Map.addLayer(force100m(non_nat), {}, "non-natural", false);
Map.addLayer(wat, {}, "water", false);
Map.addLayer(non_wat, {}, "non-water", false);



/**/ // Sieving 

//// Removing small patches of landcover within other areas
// Testing on sample image (2016 land cover)

// 0.5 sq km (5 pixels)
var nat_sieve5 = sieve(nat, 5);
var non_nat_sieve5 = sieve(non_nat, 5);

// 1 sq km (10 pixels)
var nat_sieve10 = sieve(nat, 10);
var non_nat_sieve10 = sieve(non_nat, 10);

// 5 sq km (50 pixels)
var nat_sieve50 = sieve(nat, 50);
var non_nat_sieve50 = sieve(non_nat, 50);

// Using force100m function to visualize in resolution layers will be analyzed in
Map.addLayer(force100m(nat_sieve5), {}, "nat_sieve5", false);
Map.addLayer(force100m(non_nat_sieve5), {}, "non_nat_sieve5", false);

Map.addLayer(force100m(nat_sieve10), {}, "nat_sieve10", false);
Map.addLayer(force100m(non_nat_sieve10), {}, "non_nat_sieve10", false);

Map.addLayer(force100m(nat_sieve50), {}, "nat_sieve50", false);
Map.addLayer(force100m(non_nat_sieve50), {}, "non_nat_sieve50", false);


// Visual analysis shows that 0.5km cutoff is best (losing too many patches with higher cutoffs)


/**/ //  Distance rasters

// Distance via cumulativeCost function

// Testing on sample image (2016 land cover)
var dist_non_nat_500m = nat_sieve5.cumulativeCost({
                      source: non_nat_sieve5,
                      maxDistance: 5e5
}).clip(roi);

var dist_nat_500m = non_nat_sieve5.cumulativeCost({           
                      source: nat_sieve5,
                      maxDistance: 5e5
}).clip(roi);

var dist_wat = non_wat.cumulativeCost({
                      source: wat,
                      maxDistance: 5e5
}).clip(roi);

Map.addLayer(dist_nat_500m, {max: 1000}, "distance to nat (500m)", false);
Map.addLayer(dist_non_nat_500m, {max: 1000}, "distance to non-nat (500m)", false);
Map.addLayer(dist_wat, {max: 5000}, "distance to water", false);



// Apply distance() function to MapBiomas land cover ImageCollection
var img_n = ee.List.sequence(0, lc_coll.size().subtract(1));

var dist_imgs = img_n.map(function(x){
  return distance(lc_coll, x);
});

// Use to_coll function to split into separate lists for "dist_nat" and "dist_non_nat"

var dist_nat500m_list = img_n.map(function(x){
  return to_coll(dist_imgs, x, 0);
});

var dist_non_nat500m_list = img_n.map(function(x){
  return to_coll(dist_imgs, x, 1);
});

/**/ // Adding layers with forloop

// Get the size of the image list to use evaluate() (client-side of map() to allow forloop)
var len = img_n.size();

len.evaluate(function(l) {
  for (var i=0; i < l; i++) {
    var img = ee.Image(dist_nat500m_list.get(i));
    var img1 = ee.Image(dist_non_nat500m_list.get(i));
    var yr_str = ee.String(img.get('year')).getInfo();
    Map.addLayer(img, {max:1000}, 'dist_nat500m_'+yr_str, false);
    Map.addLayer(img1, {max:1000}, 'dist_non_nat500m_' +yr_str, false);
  } 
});


/**/ // Batch Download
var natColl = ee.ImageCollection.fromImages(dist_nat500m_list);
var non_natColl = ee.ImageCollection.fromImages(dist_non_nat500m_list);

var utm21s = "EPSG:32721";
var longlat = "EPSG:4326";

var batch = require('users/fitoprincipe/geetools:batch');


// batch.Download.ImageCollection.toDrive(natColl, 'distance_natural', 
//                 {scale: 250, 
//                 region: polys, 
//                 type: 'float',
//                 crs: longlat});

// batch.Download.ImageCollection.toDrive(non_natColl, 'distance_non_natural', 
//                 {scale: 250, 
//                 region: polys, 
//                 type: 'float',
//                 crs: longlat});


/**/ // Functions

function sieve(img, nPix) {
  var con = ee.Image(img).connectedPixelCount(100, true); // 8 neighbor connected
  var mask = con.gt(nPix);
  var sieved = ee.Image(img).updateMask(mask).unmask();
  
  return(sieved);
}

function force100m(img) {return img.reproject('epsg:4326', null, 100)}

function distance(coll, y) {
  var list = coll.toList(img_n.size());
  var targ = ee.Image(list.get(y));         //target year image
  var nat = ee.Image(targ).lt(3).or(
            ee.Image(targ).gt(7));
  var non_nat = nat.eq(0);
  var nat_500m = sieve(nat, 5);//.gt(0);
  var non_nat_500m = sieve(non_nat, 5);//.gt(0);
  
  var dist_nat = non_nat_500m.cumulativeCost({
                            source: nat_500m,
                            maxDistance: 5e5
    }).set('id', 'dist_nat','id_number', ee.Number(y), 
           'year', ee.String(targ.get('system:index')).slice(-4));
  
  var dist_non_nat = nat_500m.cumulativeCost({
                            source: non_nat_500m,
                            maxDistance: 5e5
    }).set('id', 'dist_non_nat','id_number', ee.Number(y), 
           'year', ee.String(targ.get('system:index')).slice(-4));
 
  return [dist_nat, dist_non_nat];
}

function to_coll(list,n,y) {
  return ee.List(list.get(n)).get(y);
}
