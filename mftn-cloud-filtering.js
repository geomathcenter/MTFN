/****
 * Autor: Mg. Rolando Renee Badaracco Meza, Mg. Jhon Brayan Guerrero Salinas, Mg. Michel Aristóteles Choccña Rejas 
 * Institución: GEOMATH CENTER ---> https://www.facebook.com/GeoMathCenter          
 * Contacto: geomathcenter@gmail.com                                                
 * Fecha de creación: 2026-01-01                                                    
 * Última modificación: 2026-05-01                                                   
 * Versión: 1.0                                                                      
 *                                                                                    
 * Descripción:                                                                       
 * Aplicación en Google Earth Engine para generar compuestos Landsat.                   
 *                                                                                             
 * /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 * APP GEE — Landsat 5/7/8/9 por PATH/ROW
 * Filtros seleccionables: nubes, sombras, cirrus, nieve, fill, saturación y MFTN simple.
 * Mascara de Filtrado Temporal de Nubes (MFTN), es un modelo simplificado del modelo ---> Multi-Temporal Cloud Detection (MTCD) ---> https://www.sciencedirect.com/science/article/abs/pii/S0034425710000908

                            MFTN simple	                                                         MTCD original
  ----------------------------------------------------------------||--------------------------------------------------------------------------------------------------------------------------------------
  Usa una mediana temporal de toda la colección como referencia.  || Compara imágenes en una serie temporal procesada en orden cronológico.
  Detecta aumento en la banda azul.	                              || También detecta aumento en la banda azul.
  Usa la banda roja como filtro simple.	                          || Además evalúa si el píxel tiene un espectro más “blanco” que antes.
  No analiza vecindarios espaciales.	                            || El MTCD original combina el aumento en azul con una prueba de correlación lineal de vecindarios de píxeles entre imágenes sucesivas.
  No procesa imagen por imagen de forma dependiente.	            || El método original requiere procesar la serie en orden, porque compara la imagen actual con imágenes anteriores o referencias previas.
  Es más simple y práctico para Google Earth Engine.	            || Es un algoritmo más completo de detección multitemporal.
  ----------------------------------------------------------------||--------------------------------------------------------------------------------------------------------------------------------------
 
 * Expansión de máscara: 3x3, 5x5 o 7x7 para eliminar bordes alrededor de píxeles enmascarados.
 * Correcciones opcionales: BRDF - ADPATADO DE (https://github.com/ndminhhus/geeguide/blob/master/04.topo_correction.md)
 * Compuestos disponibles: mediana, medoide, media, mínimo, máximo, desviación estándar y percentil.
 * Exportación: crea tarea Export.image.toDrive o Export.image.toAsset en proyección geográfica EPSG:4326.
 ****/

// ==================== PALETA DE COLORES GLOBAL ====================
var PALETTE = {
  primary:    '#2563eb',
  secondary:  '#059669',
  accent:     '#dc2626',
  warning:    '#d97706',
  text:       '#1f2937',
  muted:      '#6b7280',
  bg:         '#f9fafb',
  border:     '#e5e7eb',
  dark:       '#111827',
  light:      '#ffffff'
};

// ==================== MAPA BASE DARK THEME ====================
var DARKMAP = [
  { stylers: [{ invert_lightness: true }, { saturation: -100 }] },
  { elementType: 'labels', stylers: [{ visibility: 'off' }] },
  { featureType: 'water', stylers: [{ color: '#1a1a2e' }] },
  { featureType: 'landscape', stylers: [{ color: '#16213e' }] },
  { featureType: 'road', stylers: [{ visibility: 'off' }] },
  { featureType: 'poi', stylers: [{ visibility: 'off' }] },
  { featureType: 'transit', stylers: [{ visibility: 'off' }] }
];

Map.setOptions('Dark', { 'Dark': DARKMAP });
Map.style().set('cursor', 'crosshair');

// =====================================================================
// CONFIGURACIÓN GENERAL
// =====================================================================
var PI = Math.PI;

var OPTICAL_BANDS = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2'];

var VIS_FALSE_COLOR = {
  bands: ['swir2', 'nir', 'red'],
  min: [0.01, 0.05, 0.03],
  max: [0.2, 0.6, 0.25]
};

var EXPORT_CRS = 'EPSG:4326';

var COEF = {
  blue:  {fiso: 0.0774, fgeo: 0.0079, fvol: 0.0372},
  green: {fiso: 0.1306, fgeo: 0.0178, fvol: 0.0580},
  red:   {fiso: 0.1690, fgeo: 0.0227, fvol: 0.0574},
  nir:   {fiso: 0.3093, fgeo: 0.0330, fvol: 0.1535},
  swir1: {fiso: 0.3430, fgeo: 0.0453, fvol: 0.1154},
  swir2: {fiso: 0.2658, fgeo: 0.0387, fvol: 0.0639}
};

var lastComposite = null;
var lastRawCollection = null;
var lastProcessedCollection = null;
var lastRegion = null;
var lastParams = null;

// =====================================================================
// PREPROCESAMIENTO LANDSAT C2 L2
// =====================================================================
function maskLandsatC2L2(img, params) {
  params = params || {};

  var qa = img.select('QA_PIXEL');
  var radsat = img.select('QA_RADSAT');

  var mask = ee.Image.constant(1);

  if (params.maskFill) {
    mask = mask.and(qa.bitwiseAnd(1 << 0).eq(0));
  }

  if (params.maskDilatedCloud) {
    mask = mask.and(qa.bitwiseAnd(1 << 1).eq(0));
  }

  if (params.maskCirrus) {
    mask = mask.and(qa.bitwiseAnd(1 << 2).eq(0));
  }

  if (params.maskCloud) {
    mask = mask.and(qa.bitwiseAnd(1 << 3).eq(0));
  }

  if (params.maskCloudShadow) {
    mask = mask.and(qa.bitwiseAnd(1 << 4).eq(0));
  }

  if (params.maskSnow) {
    mask = mask.and(qa.bitwiseAnd(1 << 5).eq(0));
  }

  if (params.maskSaturation) {
    mask = mask.and(radsat.eq(0));
  }

  return mask;
}

function scaleOptical(img, inputBands) {
  return img
    .select(inputBands, OPTICAL_BANDS)
    .multiply(0.0000275)
    .add(-0.2)
    .toFloat();
}

function applyMaskExpansion(image, windowSize) {
  windowSize = parseInt(windowSize, 10);

  if (!windowSize || windowSize < 3) {
    return ee.Image(image).set('MASK_EXPANSION', 'none');
  }

  var radius = (windowSize - 1) / 2;
  var validMask = ee.Image(image).select('blue').mask();

  var expandedValidMask = validMask.focal_min({
    radius: radius,
    kernelType: 'square',
    units: 'pixels'
  });

  return ee.Image(image)
    .updateMask(expandedValidMask)
    .set('MASK_EXPANSION', windowSize + 'x' + windowSize);
}

function prepL457(img, params) {
  var mask = maskLandsatC2L2(img, params);

  var optical = scaleOptical(
    img,
    ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7']
  ).updateMask(mask);

  return optical
    .copyProperties(img, img.propertyNames())
    .set('sensor_group', 'L5_L7')
    .set('date_ymd', ee.Date(img.get('system:time_start')).format('YYYYMMdd'));
}

function prepL89(img, params) {
  var mask = maskLandsatC2L2(img, params);

  var optical = scaleOptical(
    img,
    ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7']
  ).updateMask(mask);

  return optical
    .copyProperties(img, img.propertyNames())
    .set('sensor_group', 'L8_L9')
    .set('date_ymd', ee.Date(img.get('system:time_start')).format('YYYYMMdd'));
}

function baseFilters(path, row, startDate, endDate) {
  return ee.Filter.and(
    ee.Filter.eq('WRS_PATH', path),
    ee.Filter.eq('WRS_ROW', row),
    ee.Filter.date(startDate, endDate)
  );
}

function buildLandsatCollection(params) {
  var filter = baseFilters(
    params.path,
    params.row,
    params.startDate,
    params.endDate
  );

  var col = ee.ImageCollection([]);

  if (params.useL5) {
    col = col.merge(
      ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
        .filter(filter)
        .map(function(img) {
          return prepL457(img, params);
        })
    );
  }

  if (params.useL7) {
    col = col.merge(
      ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
        .filter(filter)
        .map(function(img) {
          return prepL457(img, params);
        })
    );
  }

  if (params.useL8) {
    col = col.merge(
      ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
        .filter(filter)
        .map(function(img) {
          return prepL89(img, params);
        })
    );
  }

  if (params.useL9) {
    col = col.merge(
      ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
        .filter(filter)
        .map(function(img) {
          return prepL89(img, params);
        })
    );
  }

  return col.sort('system:time_start');
}

// =====================================================================
// MFTN — MÁSCARA DE FILTRADO TEMPORAL DE NUBES POR GRUPO ESPECTRAL
// =====================================================================
function applyMTCDSimpleToGroup(subCollection, params, groupName) {
  subCollection = ee.ImageCollection(subCollection);

  var reference = subCollection.select(OPTICAL_BANDS).median();
  var threshold = ee.Number(params.mtcdThreshold);

  return subCollection.map(function(img) {
    var blueIncrease = img.select('blue')
      .subtract(reference.select('blue'));

    var redIncrease = img.select('red')
      .subtract(reference.select('red'));

    var blueTest = blueIncrease.gt(threshold);
    var redChangeDominates = redIncrease.gt(blueIncrease.multiply(1.5));
    var possibleCloud = blueTest.and(redChangeDominates.not());

    return img
      .updateMask(possibleCloud.not())
      .set('MTCD_simple_applied', 1)
      .set('MTCD_blue_threshold', params.mtcdThreshold)
      .set('MTCD_reference_group', groupName)
      .set('MTCD_reference_type', 'median_by_sensor_group');
  });
}

function applyMTCDSimple(collection, params) {
  if (!params.applyMTCD) {
    return collection;
  }

  var collectionL457 = collection.filter(ee.Filter.eq('sensor_group', 'L5_L7'));
  var collectionL89 = collection.filter(ee.Filter.eq('sensor_group', 'L8_L9'));

  var filteredL457 = applyMTCDSimpleToGroup(collectionL457, params, 'L5_L7');
  var filteredL89 = applyMTCDSimpleToGroup(collectionL89, params, 'L8_L9');

  return ee.ImageCollection(filteredL457.merge(filteredL89))
    .sort('system:time_start');
}

function applyMaskExpansionToCollection(collection, params) {
  return collection.map(function(img) {
    return applyMaskExpansion(img, params.maskExpansionWindow)
      .copyProperties(img, img.propertyNames());
  });
}

// =====================================================================
// BRDF
// =====================================================================
function getSunAngles(image) {
  var sunElev = ee.Number(image.get('SUN_ELEVATION'));
  var sunAz = ee.Number(image.get('SUN_AZIMUTH'));

  sunElev = ee.Number(ee.Algorithms.If(sunElev, sunElev, 45.0));
  sunAz = ee.Number(ee.Algorithms.If(sunAz, sunAz, 180.0));

  var sunZen = ee.Image(90).subtract(sunElev).multiply(PI).divide(180);
  var sunAzRad = ee.Image.constant(sunAz).multiply(PI).divide(180);

  return {
    sunZen: sunZen,
    sunAzRad: sunAzRad,
    sunElev: sunElev,
    sunAz: sunAz
  };
}

function getViewAngles() {
  return {
    viewZen: ee.Image.constant(0.0),
    viewAz: ee.Image.constant(0.0)
  };
}

function getSunZenOut(geometry) {
  var centerLat = ee.Number(geometry.centroid(30).coordinates().get(1));

  var lat2 = centerLat.multiply(centerLat);
  var lat3 = lat2.multiply(centerLat);
  var lat4 = lat3.multiply(centerLat);
  var lat5 = lat4.multiply(centerLat);
  var lat6 = lat5.multiply(centerLat);

  var sunZenOutDeg = ee.Image(0).expression(
    '31.0076 - 0.1272 * lat + 0.01187 * lat2 + 0.000024 * lat3 - 0.000000948 * lat4 - 0.00000000195 * lat5 + 0.0000000000615 * lat6',
    {
      lat: centerLat,
      lat2: lat2,
      lat3: lat3,
      lat4: lat4,
      lat5: lat5,
      lat6: lat6
    }
  );

  return sunZenOutDeg.multiply(PI).divide(180);
}

function rossThick(sunZen, viewZen, relAz) {
  var cosPhase = ee.Image(0).expression(
    'cos(sz) * cos(vz) + sin(sz) * sin(vz) * cos(raz)',
    {
      sz: sunZen,
      vz: viewZen,
      raz: relAz
    }
  ).clamp(-1, 1);

  var phase = cosPhase.acos();

  return ee.Image(0).expression(
    '((PI / 2 - phase) * cosPhase + sin(phase)) / (cos(sz) + cos(vz)) - PI / 4',
    {
      phase: phase,
      cosPhase: cosPhase,
      sz: sunZen,
      vz: viewZen,
      PI: PI
    }
  );
}

function liThin(sunZen, viewZen, relAz) {
  var b_r = 1;
  var h_b = 2;

  var tanSzPrime = ee.Image(0).expression(
    'b_r * tan(sz)',
    {
      b_r: b_r,
      sz: sunZen
    }
  ).max(0);

  var tanVzPrime = ee.Image(0).expression(
    'b_r * tan(vz)',
    {
      b_r: b_r,
      vz: viewZen
    }
  ).max(0);

  var szPrime = tanSzPrime.atan();
  var vzPrime = tanVzPrime.atan();

  var distance = ee.Image(0).expression(
    'sqrt(pow(tanSz, 2) + pow(tanVz, 2) - 2 * tanSz * tanVz * cos(raz))',
    {
      tanSz: tanSzPrime,
      tanVz: tanVzPrime,
      raz: relAz
    }
  );

  var temp = ee.Image(0).expression(
    '1 / cos(sz) + 1 / cos(vz)',
    {
      sz: szPrime,
      vz: vzPrime
    }
  );

  var cosT = ee.Image(0).expression(
    'h_b * sqrt(pow(d, 2) + pow(tanSz * tanVz * sin(raz), 2)) / temp',
    {
      h_b: h_b,
      d: distance,
      tanSz: tanSzPrime,
      tanVz: tanVzPrime,
      raz: relAz,
      temp: temp
    }
  ).clamp(-1, 1);

  var t = cosT.acos();

  var overlap = ee.Image(0).expression(
    '(1 / PI) * (t - sin(t) * cosT) * temp',
    {
      t: t,
      cosT: cosT,
      temp: temp,
      PI: PI
    }
  ).max(0);

  var cosPhasePrime = ee.Image(0).expression(
    'cos(sz) * cos(vz) + sin(sz) * sin(vz) * cos(raz)',
    {
      sz: szPrime,
      vz: vzPrime,
      raz: relAz
    }
  );

  return ee.Image(0).expression(
    'overlap - temp + 0.5 * (1 + cosPhasePrime) / (cos(sz) * cos(vz))',
    {
      overlap: overlap,
      temp: temp,
      cosPhasePrime: cosPhasePrime,
      sz: szPrime,
      vz: vzPrime
    }
  );
}

function applyBRDF(image, brdfMultiplier) {
  var sun = getSunAngles(image);
  var view = getViewAngles();

  var sunZenOut = getSunZenOut(image.geometry());
  var viewZenOut = ee.Image.constant(0.0);
  var relAz = sun.sunAzRad.subtract(view.viewAz);
  var relAzOut = ee.Image.constant(0.0);

  var kvol = rossThick(sun.sunZen, view.viewZen, relAz);
  var kgeo = liThin(sun.sunZen, view.viewZen, relAz);
  var kvol0 = rossThick(sunZenOut, viewZenOut, relAzOut);
  var kgeo0 = liThin(sunZenOut, viewZenOut, relAzOut);

  var correctedImages = [];

  for (var i = 0; i < OPTICAL_BANDS.length; i++) {
    var bandName = OPTICAL_BANDS[i];
    var coef = COEF[bandName];

    var fiso = ee.Image.constant(coef.fiso);
    var fgeo = ee.Image.constant(coef.fgeo);
    var fvol = ee.Image.constant(coef.fvol);

    var brdf = fiso
      .add(kvol.multiply(fvol).multiply(brdfMultiplier))
      .add(kgeo.multiply(fgeo));

    var brdf0 = fiso
      .add(kvol0.multiply(fvol).multiply(brdfMultiplier))
      .add(kgeo0.multiply(fgeo));

    var cFactor = brdf0.divide(brdf).clamp(0.5, 2.0);

    correctedImages.push(
      image.select(bandName).multiply(cFactor).rename(bandName)
    );
  }

  return ee.Image.cat(correctedImages)
    .copyProperties(image, image.propertyNames())
    .set('BRDF_applied', 1)
    .set('BRDF_multiplier', brdfMultiplier)
    .set('BRDF_sun_elevation_deg', sun.sunElev)
    .set('BRDF_sun_azimuth_deg', sun.sunAz);
}

// =====================================================================
// COMPUESTOS Y MÉTRICAS
// =====================================================================
function medianComposite(collection) {
  return collection
    .select(OPTICAL_BANDS)
    .median()
    .rename(OPTICAL_BANDS)
    .toFloat()
    .set('COMPOSITE_METHOD', 'median');
}

function medoidComposite(collection) {
  var median = collection.select(OPTICAL_BANDS).median();

  var withQuality = collection.map(function(img) {
    var distance = img.select(OPTICAL_BANDS)
      .subtract(median)
      .pow(2)
      .reduce(ee.Reducer.sum())
      .sqrt()
      .multiply(-1)
      .rename('quality');

    return img.addBands(distance);
  });

  return withQuality
    .qualityMosaic('quality')
    .select(OPTICAL_BANDS)
    .rename(OPTICAL_BANDS)
    .toFloat()
    .set('COMPOSITE_METHOD', 'medoid');
}

function meanComposite(collection) {
  return collection
    .select(OPTICAL_BANDS)
    .mean()
    .rename(OPTICAL_BANDS)
    .toFloat()
    .set('COMPOSITE_METHOD', 'mean');
}

function minComposite(collection) {
  return collection
    .select(OPTICAL_BANDS)
    .min()
    .rename(OPTICAL_BANDS)
    .toFloat()
    .set('COMPOSITE_METHOD', 'min');
}

function maxComposite(collection) {
  return collection
    .select(OPTICAL_BANDS)
    .max()
    .rename(OPTICAL_BANDS)
    .toFloat()
    .set('COMPOSITE_METHOD', 'max');
}

function stdDevComposite(collection) {
  return collection
    .select(OPTICAL_BANDS)
    .reduce(ee.Reducer.stdDev())
    .rename(OPTICAL_BANDS)
    .toFloat()
    .set('COMPOSITE_METHOD', 'stdDev');
}

function percentileComposite(collection, percentileValue) {
  var p = parseFloat(percentileValue);

  if (isNaN(p)) {
    p = 90;
  }

  p = Math.max(0, Math.min(100, p));

  return collection
    .select(OPTICAL_BANDS)
    .reduce(ee.Reducer.percentile([p]))
    .rename(OPTICAL_BANDS)
    .toFloat()
    .set('COMPOSITE_METHOD', 'percentile')
    .set('PERCENTILE_VALUE', p);
}

function applyOptionalCorrections(collection, params) {
  var out = collection;

  if (params.applyBRDF) {
    out = out.map(function(img) {
      return applyBRDF(img, params.brdfMultiplier);
    });
  }

  return out;
}

function makeComposite(collection, params) {
  var composite;

  if (params.compositeMethod === 'Medoid') {
    composite = medoidComposite(collection);

  } else if (params.compositeMethod === 'Median') {
    composite = medianComposite(collection);

  } else if (params.compositeMethod === 'Mean') {
    composite = meanComposite(collection);

  } else if (params.compositeMethod === 'Minimum') {
    composite = minComposite(collection);

  } else if (params.compositeMethod === 'Maximum') {
    composite = maxComposite(collection);

  } else if (params.compositeMethod === 'Standard Deviation') {
    composite = stdDevComposite(collection);

  } else if (params.compositeMethod.indexOf('Percentile') === 0) {
    composite = percentileComposite(collection, params.percentileValue);

  } else {
    composite = medianComposite(collection);
  }

  var output = ee.Image(composite)
    .set('PATH', params.path)
    .set('ROW', params.row)
    .set('START_DATE', params.startDate)
    .set('END_DATE', params.endDate)
    .set('BRDF_requested', params.applyBRDF ? 1 : 0)
    .set('LANDSAT_5_used', params.useL5 ? 1 : 0)
    .set('LANDSAT_7_used', params.useL7 ? 1 : 0)
    .set('LANDSAT_8_used', params.useL8 ? 1 : 0)
    .set('LANDSAT_9_used', params.useL9 ? 1 : 0)
    .set('MASK_cloud', params.maskCloud ? 1 : 0)
    .set('MASK_cloud_shadow', params.maskCloudShadow ? 1 : 0)
    .set('MASK_cirrus', params.maskCirrus ? 1 : 0)
    .set('MASK_dilated_cloud', params.maskDilatedCloud ? 1 : 0)
    .set('MASK_snow', params.maskSnow ? 1 : 0)
    .set('MASK_fill', params.maskFill ? 1 : 0)
    .set('MASK_saturation', params.maskSaturation ? 1 : 0)
    .set('MTCD_simple_enabled', params.applyMTCD ? 1 : 0)
    .set('MTCD_blue_threshold', params.mtcdThreshold)
    .set(
      'MASK_expansion_kernel',
      params.maskExpansionWindow && params.maskExpansionWindow >= 3
        ? params.maskExpansionWindow + 'x' + params.maskExpansionWindow
        : 'none'
    )
    .set('PERCENTILE_VALUE', params.percentileValue);

  return ee.Image(output);
}

// =====================================================================
// INTERFAZ PROFESIONAL
// =====================================================================
function styledLabel(text, options) {
  options = options || {};

  return ui.Label({
    value: text,
    style: {
      fontWeight: options.bold ? 'bold' : 'normal',
      fontSize: options.size || '12px',
      color: options.color || PALETTE.text,
      margin: options.margin || '0 6px 0 0',
      padding: options.padding || '4px 0 0 0',
      backgroundColor: options.bg || 'rgba(0,0,0,0)'
    }
  });
}

function styledCheck(label, value, color) {
  return ui.Checkbox({
    label: label,
    value: value,
    style: {
      fontWeight: 'bold',
      color: color || PALETTE.primary,
      margin: '2px 16px 5px 0'
    }
  });
}

var panel = ui.Panel({
  style: {
    width: '400px',
    padding: '12px',
    backgroundColor: PALETTE.bg
  }
});

// HEADER
var headerPanel = ui.Panel({
  widgets: [
    ui.Label({
      value: 'GEOMATH CENTER',
      style: {
        fontSize: '10px',
        fontWeight: 'bold',
        color: PALETTE.primary,
        margin: '0 0 2px 0'
      }
    }),
    ui.Label({
      value: 'Landsat Composite Studio',
      style: {
        fontSize: '18px',
        fontWeight: 'bold',
        color: PALETTE.dark,
        margin: '0 0 4px 0'
      }
    }),
    ui.Label({
      value: 'Multi-sensor · BRDF · MFTN Cloud Filtering',
      style: {
        fontSize: '10px',
        color: PALETTE.muted,
        margin: '0 0 8px 0'
      }
    }),
    ui.Label({
      value: 'v2.3 — Global Release',
      style: {
        fontSize: '9px',
        color: PALETTE.muted,
        margin: '0',
        fontStyle: 'italic'
      }
    })
  ],
  layout: ui.Panel.Layout.flow('vertical'),
  style: {
    backgroundColor: PALETTE.light,
    border: '1px solid ' + PALETTE.border,
    padding: '14px',
    margin: '0 0 10px 0',
    borderRadius: '8px',
    textAlign: 'center'
  }
});

panel.add(headerPanel);

function createSection(title, description, widgets) {
  var section = ui.Panel({
    style: {
      backgroundColor: PALETTE.light,
      border: '1px solid ' + PALETTE.border,
      padding: '10px',
      margin: '0 0 8px 0',
      borderRadius: '6px'
    }
  });

  section.add(ui.Label({
    value: title,
    style: {
      fontWeight: 'bold',
      fontSize: '13px',
      color: PALETTE.primary,
      margin: '0 0 2px 0'
    }
  }));

  if (description) {
    section.add(ui.Label({
      value: description,
      style: {
        fontSize: '10px',
        color: PALETTE.muted,
        margin: '0 0 8px 0',
        fontStyle: 'italic'
      }
    }));
  }

  for (var i = 0; i < widgets.length; i++) {
    section.add(widgets[i]);
  }

  return section;
}

// SECCIÓN 1
var pathBox = ui.Textbox({
  placeholder: 'Ex: 6',
  value: '6',
  style: {
    width: '70px',
    margin: '0',
    padding: '2px'
  }
});

var rowBox = ui.Textbox({
  placeholder: 'Ex: 66',
  value: '66',
  style: {
    width: '70px',
    margin: '0',
    padding: '2px'
  }
});

var wrsPanel = ui.Panel({
  widgets: [
    ui.Label({
      value: 'PATH:',
      style: {
        fontWeight: 'bold',
        fontSize: '12px',
        color: PALETTE.primary,
        margin: '0 6px 0 0',
        padding: '4px 0 0 0'
      }
    }),
    pathBox,
    ui.Label({
      value: 'ROW:',
      style: {
        fontWeight: 'bold',
        fontSize: '12px',
        color: PALETTE.primary,
        margin: '0 6px 0 20px',
        padding: '4px 0 0 0'
      }
    }),
    rowBox
  ],
  layout: ui.Panel.Layout.flow('horizontal'),
  style: {
    stretch: 'horizontal',
    margin: '0',
    padding: '0'
  }
});

panel.add(createSection(
  '1. Orbit Path / Row (WRS-2)',
  'Define the Landsat orbit coordinates for your area of interest',
  [wrsPanel]
));

// SECCIÓN 2
var startBox = ui.Textbox({
  placeholder: 'YYYY-MM-DD',
  value: '2020-01-01',
  style: {
    width: '100px',
    margin: '0',
    padding: '2px'
  }
});

var endBox = ui.Textbox({
  placeholder: 'YYYY-MM-DD',
  value: '2020-12-31',
  style: {
    width: '100px',
    margin: '0',
    padding: '2px'
  }
});

var datePanel = ui.Panel({
  widgets: [
    ui.Label({
      value: 'Start:',
      style: {
        fontWeight: 'bold',
        fontSize: '12px',
        color: PALETTE.primary,
        margin: '0 6px 0 0',
        padding: '4px 0 0 0'
      }
    }),
    startBox,
    ui.Label({
      value: 'End:',
      style: {
        fontWeight: 'bold',
        fontSize: '12px',
        color: PALETTE.primary,
        margin: '0 6px 0 20px',
        padding: '4px 0 0 0'
      }
    }),
    endBox
  ],
  layout: ui.Panel.Layout.flow('horizontal'),
  style: {
    stretch: 'horizontal',
    margin: '0',
    padding: '0'
  }
});

panel.add(createSection(
  '2. Time Period',
  'Select the date range for image acquisition',
  [datePanel]
));

// SECCIÓN 3
var l5Check = styledCheck('Landsat 5', true, PALETTE.primary);
var l7Check = styledCheck('Landsat 7', true, PALETTE.primary);
var l8Check = styledCheck('Landsat 8', true, PALETTE.primary);
var l9Check = styledCheck('Landsat 9', true, PALETTE.primary);

var sensorPanel = ui.Panel({
  widgets: [l5Check, l7Check, l8Check, l9Check],
  layout: ui.Panel.Layout.flow('horizontal'),
  style: {
    stretch: 'horizontal'
  }
});

panel.add(createSection(
  '3. Satellite Sensors',
  'Select one or more Landsat missions to include',
  [sensorPanel]
));

// SECCIÓN 4
var maskCloudCheck = styledCheck('Clouds', true, PALETTE.accent);
var maskCloudShadowCheck = styledCheck('Cloud Shadows', true, PALETTE.accent);
var maskCirrusCheck = styledCheck('Cirrus', true, PALETTE.accent);
var maskDilatedCloudCheck = styledCheck('Dilated Cloud', true, PALETTE.accent);
var maskSnowCheck = styledCheck('Snow / Ice', true, PALETTE.accent);
var maskFillCheck = styledCheck('Fill / No Data', true, PALETTE.accent);
var maskSaturationCheck = styledCheck('Radiometric Saturation', true, PALETTE.accent);

var qualityPanel = ui.Panel({
  widgets: [
    maskCloudCheck,
    maskCloudShadowCheck,
    maskCirrusCheck,
    maskDilatedCloudCheck,
    maskSnowCheck,
    maskFillCheck,
    maskSaturationCheck
  ],
  layout: ui.Panel.Layout.flow('vertical'),
  style: {
    margin: '0 0 8px 0'
  }
});

var mtcdCheck = ui.Checkbox({
  label: 'MFTN (Mascara de Filtrado Temporal de Nubes)',
  value: false,
  style: {
    fontWeight: 'bold',
    color: PALETTE.secondary,
    margin: '8px 0 4px 0'
  }
});

var mtcdThresholdSelect = ui.Textbox({
  placeholder: 'Ex: 0.03',
  value: '0.03',
  style: {
    width: '70px',
    margin: '0',
    padding: '2px'
  }
});

var mtcdPanel = ui.Panel({
  widgets: [
    ui.Label({
      value: 'Blue band delta threshold >',
      style: {
        fontSize: '11px',
        color: PALETTE.muted,
        margin: '0 6px 0 0',
        padding: '4px 0 0 0'
      }
    }),
    mtcdThresholdSelect
  ],
  layout: ui.Panel.Layout.flow('horizontal'),
  style: {
    stretch: 'horizontal',
    margin: '0 0 0 12px'
  }
});

panel.add(createSection(
  '4. Quality Masking & Cloud Filtering',
  'Apply QA_PIXEL masks and optional multi-temporal cloud detection',
  [qualityPanel, mtcdCheck, mtcdPanel]
));

// SECCIÓN 5
var maskExpansionSelect = ui.Select({
  items: ['No expansion', '3x3', '5x5', '7x7'],
  value: 'No expansion',
  style: {
    width: '130px'
  }
});

var expansionPanel = ui.Panel({
  widgets: [
    styledLabel('Kernel size:', {
      bold: true,
      color: PALETTE.primary,
      margin: '0 6px 0 0',
      padding: '4px 0 0 0'
    }),
    maskExpansionSelect
  ],
  layout: ui.Panel.Layout.flow('horizontal'),
  style: {
    stretch: 'horizontal'
  }
});

panel.add(createSection(
  '5. Mask Edge Expansion',
  'Expand cloud/shadow mask boundaries to clean pixel edges',
  [expansionPanel]
));

// SECCIÓN 6
var brdfCheck = styledCheck('BRDF Correction (RossThick-LiThin)', false, PALETTE.secondary);

panel.add(createSection(
  '6. Radiometric Corrections',
  'Optional BRDF correction',
  [brdfCheck]
));

// SECCIÓN 7
var compositeSelect = ui.Select({
  items: [
    'Median',
    'Medoid',
    'Mean',
    'Minimum',
    'Maximum',
    'Standard Deviation',
    'Percentile 10',
    'Percentile 25',
    'Percentile 50',
    'Percentile 75',
    'Percentile 90',
    'Percentile 95'
  ],
  value: 'Median',
  style: {
    width: '180px'
  }
});

var compositePanel = ui.Panel({
  widgets: [
    styledLabel('Metric:', {
      bold: true,
      color: PALETTE.primary,
      margin: '0 6px 0 0',
      padding: '4px 0 0 0'
    }),
    compositeSelect
  ],
  layout: ui.Panel.Layout.flow('horizontal'),
  style: {
    stretch: 'horizontal'
  }
});

panel.add(createSection(
  '7. Composite Metric',
  'Statistical method for combining the time-series',
  [compositePanel]
));

// SECCIÓN 8
var exportFolderBox = ui.Textbox({
  placeholder: 'GEE_Results',
  value: 'GEE_Results',
  style: {
    width: '160px',
    margin: '0 8px 0 0',
    padding: '2px'
  }
});

var exportNameBox = ui.Textbox({
  placeholder: 'landsat_composite',
  value: 'landsat_composite',
  style: {
    width: '180px',
    margin: '0',
    padding: '2px'
  }
});

var assetIdBox = ui.Textbox({
  placeholder: 'projects/your-project/assets/composite',
  value: 'projects/your-project/assets/composite',
  style: {
    width: '260px',
    margin: '0',
    padding: '2px'
  }
});

var exportPanel = ui.Panel({
  widgets: [
    ui.Panel({
      widgets: [
        ui.Label({
          value: 'Drive folder:',
          style: {
            fontWeight: 'bold',
            fontSize: '12px',
            color: PALETTE.primary,
            margin: '0 6px 0 0',
            padding: '4px 0 0 0',
            width: '85px'
          }
        }),
        exportFolderBox
      ],
      layout: ui.Panel.Layout.flow('horizontal'),
      style: {
        stretch: 'horizontal',
        margin: '0 0 4px 0'
      }
    }),
    ui.Panel({
      widgets: [
        ui.Label({
          value: 'File name:',
          style: {
            fontWeight: 'bold',
            fontSize: '12px',
            color: PALETTE.primary,
            margin: '0 6px 0 0',
            padding: '4px 0 0 0',
            width: '85px'
          }
        }),
        exportNameBox
      ],
      layout: ui.Panel.Layout.flow('horizontal'),
      style: {
        stretch: 'horizontal',
        margin: '0 0 4px 0'
      }
    }),
    ui.Panel({
      widgets: [
        ui.Label({
          value: 'Asset ID:',
          style: {
            fontWeight: 'bold',
            fontSize: '12px',
            color: PALETTE.primary,
            margin: '0 6px 0 0',
            padding: '4px 0 0 0',
            width: '85px'
          }
        }),
        assetIdBox
      ],
      layout: ui.Panel.Layout.flow('horizontal'),
      style: {
        stretch: 'horizontal'
      }
    })
  ],
  layout: ui.Panel.Layout.flow('vertical'),
  style: {
    stretch: 'horizontal'
  }
});

panel.add(createSection(
  '8. Export Configuration',
  'Set Google Drive folder and/or Earth Engine Asset destination',
  [exportPanel]
));

// BOTONES
var processButton = ui.Button({
  label: '▶  Process & Visualize',
  style: {
    color: '#1f2937',
    fontWeight: 'bold',
    fontSize: '13px',
    stretch: 'horizontal',
    margin: '12px 0 4px 0',
    border: '2px solid ' + PALETTE.primary,
    backgroundColor: '#dbeafe'
  }
});

var exportButton = ui.Button({
  label: '⬇  Export to Google Drive',
  style: {
    color: '#1f2937',
    fontWeight: 'bold',
    fontSize: '13px',
    stretch: 'horizontal',
    margin: '4px 0',
    border: '2px solid ' + PALETTE.secondary,
    backgroundColor: '#d1fae5'
  }
});

var exportAssetButton = ui.Button({
  label: '☁  Export to Earth Engine Asset',
  style: {
    color: '#1f2937',
    fontWeight: 'bold',
    fontSize: '13px',
    stretch: 'horizontal',
    margin: '4px 0',
    border: '2px solid #7c3aed',
    backgroundColor: '#ede9fe'
  }
});

panel.add(processButton);
panel.add(exportButton);
panel.add(exportAssetButton);

// STATUS BAR
var statusBar = ui.Panel({
  widgets: [
    ui.Label({
      value: '● Ready — Configure parameters and click Process',
      style: {
        color: PALETTE.secondary,
        fontSize: '11px',
        fontWeight: 'bold'
      }
    })
  ],
  style: {
    backgroundColor: '#f0fdf4',
    border: '1px solid #bbf7d0',
    padding: '8px 10px',
    margin: '10px 0',
    borderRadius: '6px'
  }
});

panel.add(statusBar);

var countLabel = ui.Label({
  value: '',
  style: {
    color: PALETTE.muted,
    fontSize: '11px',
    margin: '4px 0 0 0',
    textAlign: 'center'
  }
});

panel.add(countLabel);

// FOOTER
var footerPanel = ui.Panel({
  widgets: [
    ui.Label({
      value: 'Powered by GEOMATH CENTER',
      style: {
        fontSize: '9px',
        color: '#9ca3af',
        margin: '12px 0 2px 0'
      }
    }),
    ui.Label({
      value: 'https://www.facebook.com/GeoMathCenter',
      style: {
        fontSize: '9px',
        color: PALETTE.primary,
        margin: '0 0 2px 0'
      }
    }),
    ui.Label({
      value: 'geomathcenter@gmail.com',
      style: {
        fontSize: '9px',
        color: '#9ca3af',
        margin: '0'
      }
    })
  ],
  style: {
    textAlign: 'center'
  }
});

panel.add(footerPanel);

ui.root.insert(0, panel);

// =====================================================================
// FUNCIONES DE LA APP
// =====================================================================
function parseNumberFromBox(box, fallback) {
  var value = parseFloat(box.getValue());
  return isNaN(value) ? fallback : value;
}

function getMaskExpansionWindowValue() {
  var selected = maskExpansionSelect.getValue();

  if (selected === '3x3') return 3;
  if (selected === '5x5') return 5;
  if (selected === '7x7') return 7;

  return 0;
}

function getPercentileFromCompositeName(name) {
  if (name === 'Percentile 10') return 10;
  if (name === 'Percentile 25') return 25;
  if (name === 'Percentile 50') return 50;
  if (name === 'Percentile 75') return 75;
  if (name === 'Percentile 90') return 90;
  if (name === 'Percentile 95') return 95;

  return 50;
}

function getParamsFromUI() {
  var selectedComposite = compositeSelect.getValue();

  return {
    path: parseInt(pathBox.getValue(), 10),
    row: parseInt(rowBox.getValue(), 10),
    startDate: startBox.getValue(),
    endDate: endBox.getValue(),

    useL5: l5Check.getValue(),
    useL7: l7Check.getValue(),
    useL8: l8Check.getValue(),
    useL9: l9Check.getValue(),

    compositeMethod: selectedComposite,
    percentileValue: getPercentileFromCompositeName(selectedComposite),

    maskExpansionWindow: getMaskExpansionWindowValue(),

    applyBRDF: brdfCheck.getValue(),
    brdfMultiplier: 1.0,

    maskCloud: maskCloudCheck.getValue(),
    maskCloudShadow: maskCloudShadowCheck.getValue(),
    maskCirrus: maskCirrusCheck.getValue(),
    maskDilatedCloud: maskDilatedCloudCheck.getValue(),
    maskSnow: maskSnowCheck.getValue(),
    maskFill: maskFillCheck.getValue(),
    maskSaturation: maskSaturationCheck.getValue(),

    applyMTCD: mtcdCheck.getValue(),
    mtcdThreshold: parseNumberFromBox(mtcdThresholdSelect, 0.03)
  };
}

function validateParams(params) {
  if (isNaN(params.path) || isNaN(params.row)) {
    return 'PATH and ROW must be valid integers.';
  }

  if (!params.useL5 && !params.useL7 && !params.useL8 && !params.useL9) {
    return 'Select at least one Landsat sensor.';
  }

  if (!params.startDate || !params.endDate) {
    return 'Start and end dates are required.';
  }

  if (params.applyMTCD && (params.mtcdThreshold <= 0 || params.mtcdThreshold > 0.3)) {
    return 'MFTN threshold out of range. Typical values: 0.03, 0.05, 0.07, 0.10.';
  }

  if (params.compositeMethod.indexOf('Percentile') === 0) {
    if (params.percentileValue < 0 || params.percentileValue > 100) {
      return 'Percentile must be between 0 and 100.';
    }
  }

  return null;
}

function updateStatus(message, type) {
  var color = PALETTE.secondary;
  var bg = '#f0fdf4';
  var border = '#bbf7d0';

  if (type === 'error') {
    color = PALETTE.accent;
    bg = '#fef2f2';
    border = '#fecaca';

  } else if (type === 'warning') {
    color = PALETTE.warning;
    bg = '#fffbeb';
    border = '#fde68a';

  } else if (type === 'processing') {
    color = PALETTE.primary;
    bg = '#eff6ff';
    border = '#bfdbfe';
  }

  statusBar.style().set('backgroundColor', bg);
  statusBar.style().set('border', '1px solid ' + border);

  statusBar.widgets().set(0, ui.Label({
    value: message,
    style: {
      color: color,
      fontSize: '11px',
      fontWeight: 'bold'
    }
  }));
}

function processAndVisualize() {
  var params = getParamsFromUI();
  var validationError = validateParams(params);

  if (validationError) {
    updateStatus('✕ ' + validationError, 'error');
    return;
  }

  updateStatus('⏳ Building collection...', 'processing');
  countLabel.setValue('');

  var rawCollection = buildLandsatCollection(params);

  rawCollection.size().evaluate(function(n) {
    if (!n || n === 0) {
      updateStatus('✕ No images found. Try expanding the date range or check PATH/ROW.', 'error');
      countLabel.setValue('');
      return;
    }

    updateStatus('⏳ Processing ' + n + ' images — applying masks & corrections...', 'processing');

    var mtcdCollection = applyMTCDSimple(rawCollection, params);
    var maskedCollection = applyMaskExpansionToCollection(mtcdCollection, params);
    var processedCollection = applyOptionalCorrections(maskedCollection, params);
    var composite = makeComposite(processedCollection, params);

    var firstImage = ee.Image(rawCollection.first());
    var region = firstImage.geometry();

    lastRawCollection = rawCollection;
    lastProcessedCollection = processedCollection;
    lastComposite = ee.Image(composite).clip(region);
    lastRegion = region;
    lastParams = params;

    Map.layers().reset();
    Map.centerObject(region, 8);

    var labelBase = 'Landsat ' + params.compositeMethod +
      ' P' + params.path + ' R' + params.row;

    Map.addLayer(lastComposite, VIS_FALSE_COLOR, labelBase, true);

    updateStatus('✓ Success — ' + n + ' images processed | Composite ready', 'success');

    countLabel.setValue(
      n + ' images · ' +
      params.compositeMethod + ' composite · ' +
      params.startDate + ' to ' +
      params.endDate
    );

    print('Parameters:', params);
    print('Raw collection:', rawCollection);
    print('MFTN collection:', mtcdCollection);
    print('Masked collection:', maskedCollection);
    print('Processed collection:', processedCollection);
    print('Final composite:', lastComposite);
  });
}

function sanitizeName(name) {
  return String(name || 'landsat_composite')
    .replace(/[^a-zA-Z0-9_-]/g, '_')
    .replace(/_+/g, '_');
}

function sanitizeAssetId(assetId) {
  return String(assetId || '')
    .replace(/\s+/g, '')
    .replace(/\/+/g, '/');
}

function exportLastComposite() {
  if (!lastComposite || !lastRegion || !lastParams) {
    updateStatus('✕ First click "Process & Visualize" to generate a composite.', 'error');
    return;
  }

  var exportName = sanitizeName(exportNameBox.getValue());
  var folder = sanitizeName(exportFolderBox.getValue());

  Export.image.toDrive({
    image: lastComposite.select(OPTICAL_BANDS).toFloat(),
    description: exportName,
    folder: folder,
    fileNamePrefix: exportName,
    region: lastRegion,
    scale: 30,
    crs: EXPORT_CRS,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  });

  updateStatus('✓ Drive export task created — Check Tasks tab', 'success');
}

function exportLastCompositeToAsset() {
  if (!lastComposite || !lastRegion || !lastParams) {
    updateStatus('✕ First click "Process & Visualize" to generate a composite.', 'error');
    return;
  }

  var exportName = sanitizeName(exportNameBox.getValue());
  var assetId = sanitizeAssetId(assetIdBox.getValue());

  if (!assetId) {
    updateStatus('✕ Enter a valid Earth Engine Asset ID.', 'error');
    return;
  }

  Export.image.toAsset({
    image: lastComposite.select(OPTICAL_BANDS).toFloat(),
    description: exportName + '_asset',
    assetId: assetId,
    region: lastRegion,
    scale: 30,
    crs: EXPORT_CRS,
    maxPixels: 1e13
  });

  updateStatus('✓ Asset export task created — Check Tasks tab', 'success');
}

processButton.onClick(processAndVisualize);
exportButton.onClick(exportLastComposite);
exportAssetButton.onClick(exportLastCompositeToAsset);

// Vista inicial global
Map.setCenter(-59.88, -7.39, 4);