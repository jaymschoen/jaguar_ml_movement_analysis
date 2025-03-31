/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var polys = ee.FeatureCollection("users/jaymschoen/jags_reg_int_ind_polygons_10k_buf"),
    wat = ee.ImageCollection("JRC/GSW1_3/YearlyHistory");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
// ===== Distance to annual and seasonal water layers ===== //

/**/ // Setup
Map.setOptions('satellite');
var longlat = "EPSG:4326";

var roi = polys.geometry();
var buf60k = polys.geometry().buffer(5e4); // 60k buffer (for distance calculations)
var years = ee.List.sequence(2000,2016,1);

// Selecting images for focal countries/areas
var watFocal = wat
      .select('waterClass')
      .filterDate('2000-01-01', '2016-12-31')
      .map(function(x){return x.clip(buf60k)});

print('Water Focal Collection', watFocal);

var visualization = {
  bands: ['waterClass'],
  min: 0.0,
  max: 3.0,
  palette: ['cccccc', 'ffffff', '99d9ea', '0000ff']
};

Map.addLayer(polys, {}, '10k polygons', false);
Map.addLayer(watFocal, visualization, 'water', false);


/**/ // Distance Calculations 
// Distance to permanent/seasonal water

// testing on 2016 image
var img = ee.Image(wat.toList(17).get(16)).clip(buf60k);
var bg = ee.Image(1).clip(buf60k);
var watPerm = img.eq(3);
var watSeas = img.gte(2);

Map.addLayer(img, visualization, '2016', false);
Map.addLayer(watPerm, {}, 'permanent water', false);
Map.addLayer(watSeas, {}, 'seasonal water', false);
Map.addLayer(bg, {}, 'background', false);

var dist_watPerm = bg.cumulativeCost({           // distance via cost function
                      source: watPerm,
                      maxDistance: 1e5
});

var dist_watSeas = bg.cumulativeCost({           // distance via cost function
                      source: watSeas,
                      maxDistance: 1e5
});

Map.addLayer(dist_watPerm, {min:100, max: 10000}, "distance to permanent water", false);
Map.addLayer(dist_watSeas, {min:100, max: 10000}, "distance to seasonal water", false);



// Function for ImageCollection
var img_n = ee.List.sequence(0, watFocal.size().subtract(1));

function distance(coll, y) {
  var list = coll.toList(img_n.size());
  var watPerm = ee.Image(list.get(y)).eq(3)
    .set('year', ee.Image(list.get(y)).get('year'));
  var watSeas = ee.Image(list.get(y)).gte(2)
    .set('year', ee.Image(list.get(y)).get('year'));

  var dist_watPerm = bg.cumulativeCost({
                            source: watPerm,
                            maxDistance: 1e5
  }).set({id: 'dist_watPerm',
          id_number: ee.Number(y),
          year: watPerm.get('year')});
  
  var dist_watSeas = bg.cumulativeCost({
                            source: watSeas,
                            maxDistance: 1e5
  }).set('id', 'dist_watSeas','id_number', ee.Number(y), 'year', watSeas.get('year'));
 
  return [dist_watPerm.clip(roi), dist_watSeas.clip(roi)];
}

var dist_imgs = img_n.map(function(x){
  return distance(watFocal, x);
});

print('Distance images', dist_imgs);

// Separate lists for permanent and seasonal water
var dist_watPerm_list = img_n.map(function(x){
  return ee.List(dist_imgs.get(x)).get(0);
});

var dist_watSeas_list = img_n.map(function(x){
  return ee.List(dist_imgs.get(x)).get(1);
});

// print(dist_watPerm_list);
// print(dist_watSeas_list);

/**/ // Adding layers with forloop

// Get the size of the image list to use evaluate() (client-side of map() to allow forloop)
var len = img_n.size();

len.evaluate(function(l) {
  for (var i=0; i < l; i++) {
    var img = ee.Image(dist_watPerm_list.get(i));
    var img1 = ee.Image(dist_watSeas_list.get(i));
    var yr = img.get('year').getInfo();
    Map.addLayer(img, {min:100, max:10000}, 'dist_watPerm_'+yr, false);
    Map.addLayer(img1, {min:100, max:10000}, 'dist_watSeas_' +yr, false);
  } 
});


/**/ // Batch Download
var permColl = ee.ImageCollection.fromImages(dist_watPerm_list);
var seasColl = ee.ImageCollection.fromImages(dist_watSeas_list);

var utm21s = "EPSG:32721";
var longlat = "EPSG:4326";

var batch = require('users/fitoprincipe/geetools:batch');

// batch.Download.ImageCollection.toDrive(permColl, 'distance_water_permanent', 
//                 {scale: 250, 
//                 region: polys, 
//                 type: 'float',
//                 crs: longlat});