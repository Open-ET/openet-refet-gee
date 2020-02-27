import daily
import ee
import math
import refet

ee.Initialize()

# random image from collection
maca_img = ee.ImageCollection('IDAHO_EPSCOR/MACAv2_METDATA')\
  .filterMetadata('ensemble', 'equals', "r1i1p1")\
  .filterMetadata('model', 'equals', "BNU-ESM").first()

# Set the output crs and crsTransform to match the maca image
crs = maca_img.projection().getInfo()['crs']
geo = maca_img.projection().getInfo()['transform']

# test point for reducer
lat = 39.52
lon = -119.30

test_pnt = ee.Geometry.Point(lon, lat)

# gee_refet eto
gee_eto = daily.Daily.maca(maca_img).eto()\
      .rename(['eto'])\
      .reduceRegion(ee.Reducer.mean(), test_pnt, crs=crs, crsTransform=geo)\
      .getInfo()['eto']

# gee_refet etr
gee_etr = daily.Daily.maca(maca_img).etr()\
      .rename(['etr'])\
      .reduceRegion(ee.Reducer.mean(), test_pnt, crs=crs, crsTransform=geo)\
      .getInfo()['etr']

# get point data for RefET calc/comparison
pt_data = maca_img.reduceRegion(ee.Reducer.mean(), test_pnt, crs=crs,
                                crsTransform=geo).getInfo()
# image doy
doy = ee.Number.parse(maca_img.date().format("dd")).getInfo()

# elevation pt extraction (gridmet elevation asset)
elev_img = ee.Image('projects/earthengine-legacy/assets/'
                    'projects/climate-engine/gridmet/elevation')
elev_crs = elev_img.projection().getInfo()['crs']
elev_geo = elev_img.projection().getInfo()['transform']
elev_pt = elev_img.reduceRegion(ee.Reducer.mean(), test_pnt, crs=elev_crs,
                                crsTransform=elev_geo).getInfo()

# gridcell lat from elevation asset
lat_img =ee.Image('projects/earthengine-legacy/assets/'
               'projects/climate-engine/gridmet/elevation') \
    .multiply(0).add(ee.Image.pixelLonLat().select('latitude'))
cell_lat = lat_img.reduceRegion(ee.Reducer.mean(), test_pnt, crs=elev_crs,
                                crsTransform=elev_geo).getInfo()['b1']

# air pressure from gridmet elevation using refet module
pair_kpa = refet.calcs._air_pressure(elev_pt['b1'], method='asce')
# actual vapor pressure (kg/kg) using refet module
ea_kpa = refet.calcs._actual_vapor_pressure(pt_data['huss'], pair_kpa)

# resultant wind speed
windspeed = math.sqrt(pt_data['uas']**2 + pt_data['vas']**2)

# refet etr calc
etr = refet.Daily(
    tmin=pt_data['tasmin'], tmax=pt_data['tasmax'], ea=ea_kpa,
    rs=pt_data['rsds']*0.0864,
    uz=windspeed, zw=10, elev=elev_pt['b1'],
    lat=cell_lat, doy=doy, method='asce',
    input_units={'tmin': 'K', 'tmax': 'K', 'uz': 'm/s',
                 'lat': 'deg', 'ea': 'kpa', 'rs': 'mj m-2 day-1'}
    ).etr()

# refet eto calc
eto = refet.Daily(
    tmin=pt_data['tasmin'], tmax=pt_data['tasmax'], ea=ea_kpa,
    rs=pt_data['rsds']*0.0864,
    uz=windspeed, zw=10, elev=elev_pt['b1'],
    lat=cell_lat, doy=doy, method='asce',
    input_units={'tmin': 'K', 'tmax': 'K', 'uz': 'm/s',
                 'lat': 'deg', 'ea': 'kpa', 'rs': 'mj m-2 day-1'}
    ).eto()

# Compare output
print('\nRefET-GEE ETo: {} mm'.format(gee_eto))
print('RefET ETo: {} mm'.format(float(eto)))
print('\nRefET-GEE ETr: {} mm'.format(gee_etr))
print('RefET ETr: {} mm'.format(float(etr)))

