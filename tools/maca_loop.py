import geerefet.daily
import ee
# import geopandas as gpd

ee.Initialize()

# test range
# start_date = '2020-01-01'
# end_date = '2020-01-31'
# Full rcp45 and rcp85 date range (2006,2099)
start_date = '1980-01-01'
end_date = '1980-12-31'

scenario_list = ['rcp45', 'rcp85', 'historical']
# future_scenario_list = ['rcp45', 'rcp85']

model_list = ['bcc-csm1-1',
    'bcc-csm1-1-m',
    'BNU-ESM',
    'CanESM2',
    'CCSM4',
    'CNRM-CM5',
    'CSIRO-Mk3-6-0',
    'GFDL-ESM2G',
    'GFDL-ESM2M',
    'HadGEM2-CC365',
    'HadGEM2-ES365',
    'inmcm4',
    'IPSL-CM5A-MR',
    'IPSL-CM5A-LR',
    'IPSL-CM5B-LR',
    'MIROC5',
    'MIROC-ESM',
    'MIROC-ESM-CHEM',
    'MRI-CGCM3',
    'NorESM1-M']

# TESTING LISTS
future_scenario_list = ['historical']
model_list = ['bcc-csm1-1', 'bcc-csm1-1-m']

# ftr from ee asset upload (add .shp to ftr option)
ftr = ee.FeatureCollection('users/ChrisCPearson86/Urban_Reproject')
# import itertools
# for i, j in itertools.product(range(x), range(y)):
for scenario in future_scenario_list:
    for model in model_list:
        print('Running {}, {}'.format(scenario, model))
        maca_coll = ee.ImageCollection('IDAHO_EPSCOR/MACAv2_METDATA')\
            .filterMetadata('scenario', 'equals', scenario)\
            .filterMetadata('model', 'equals', model)
            # .filterDate(start_date, end_date)


        # print(maca_coll.first())
        crs = maca_coll.first().projection().getInfo()['crs']
        geo = maca_coll.first().projection().getInfo()['transform']

        # // Extract the values by running reduceRegions over each image
        # // in the image collection.
        def maca_reduceregions(i):
            dateStr = i.date()
            datenum = ee.Image.constant(ee.Number.parse(dateStr.format("YYYYMMdd")))
            # Add dateNum Band to Image
            dateYear= ee.Number.parse(dateStr.format("YYYY"))
            dateMonth= ee.Number.parse(dateStr.format("MM"))
            dateDay= ee.Number.parse(dateStr.format("dd"))
            i = i.addBands(ee.Image(datenum).rename('datenum'))
            i = i.addBands(ee.Image(dateYear).rename('year'))
            i = i.addBands(ee.Image(dateMonth).rename('month'))
            i = i.addBands(ee.Image(dateDay).rename('day'))
            # geerefet eto calc
            eto = geerefet.daily.Daily.maca(i).eto()
            i = i.addBands(ee.Image(eto).rename('eto_mm'))
            # geerefet etr calc
            etr = geerefet.daily.Daily.maca(i).etr()
            i = i.addBands(ee.Image(etr).rename('etr_mm'))
            # resultant wind speed from vectors
            ws = ee.Image(i.select(['uas'])).pow(2) \
                    .add(ee.Image(i.select(['vas'])).pow(2)) \
                    .sqrt().rename(['ws'])
            i = i.addBands(ee.Image(ws).rename('uz'))
            zw = 10 # maca wind vector measurement height (meters)
            # Wind speed(Eqn 67) adjust to 2m
            u2 = ws.expression(
                'uz * 4.87 / log(67.8 * zw - 5.42)', {'uz': ws, 'zw': zw})
            i = i.addBands(ee.Image(u2).rename('ws_2m'))

            # Unit changes
            tmmx = i.select(['tasmax'], ['tmax_c']).subtract(273.15) # K to C
            i = i.addBands(ee.Image(tmmx), ['tmax_c'])
            tmmn = i.select(['tasmin'], ['tmin_c']).subtract(273.15) # K to C
            i = i.addBands(ee.Image(tmmn), ['tmin_c'])

            # reduce regions (select desired met variables and rename)
            fc = i.select(['datenum', 'year', 'month', 'day', 'tmax_c', 'tmin_c', 'rsds', 'uz', 'huss', 'pr',
                 'ws_2m', 'eto_mm', 'etr_mm'],
                ['datenum', 'year', 'month', 'day', 'tmax_c', 'tmin_c', 'srad_wm2', 'uz', 'q_kgkg', 'prcp_mm',
                 'ws_2m', 'eto_mm', 'etr_mm']).reduceRegions(collection=ftr,
                                 reducer=ee.Reducer.mean(),
                                 crs=crs,
                                 crsTransform=geo)
            return ee.FeatureCollection(fc).set('date', ee.String(dateStr.format('YYYY-MM-dd')))

        output = maca_coll.map(maca_reduceregions)
        # print(output.getInfo())

        # Export Summary Table
        ee.batch.Export.table.toDrive(
            collection=ee.FeatureCollection(output.flatten()),
            selectors=['datenum', 'year', 'month', 'day', 'UACE10', 'tmax_c', 'tmin_c', 'srad_wm2', 'uz', 'q_kgkg',
                       'prcp_mm', 'ws_2m', 'eto_mm', 'etr_mm'],
            description='{}_{}'.format(scenario, model),
            fileFormat='CSV',
            folder='ee_exports').start()





# DEV ZONE

# var maca = ee.ImageCollection('IDAHO_EPSCOR/MACAv2_METDATA')
#   .select("tasmax")
#   // .filterMetadata('ensemble', 'equals', "r1i1p1")
#   .filterMetadata('model', 'equals', "BNU-ESM").first()
#   // .filterMetadata('month', 'greater_than', 5)
#   // .filterMetadata('month', 'less_than', 11).first()
#   // .filterBounds(ee.Geometry.Point(25.8544, -18.08874))
#   // .filterDate('2020-01-01', '2029-12-31')
#   // .mean().subtract(273.15).multiply(9/5).add(32);
# print(maca);


#
#
# test = fcToGdf(output)
# print(test.head(5))
# output_df = output.flatten()

# shapefile = gpd.read_file("C:\Users\cpearson\Desktop\urban\Census2010_UrbanAreas_WUS_200kAmmend.shp")
#
# features = []
# for i in range(shapefile.shape[0]):
#     geom = shapefile.iloc[i:i+1,:]
#     jsonDict = eval(geom.to_json())
#     geojsonDict = jsonDict['features'][0]
#     features.append(ee.Feature(geojsonDict))
#
# fc = ee.FeatureCollection(features)

# def fcToGdf(fc, crs={'init': 'epsg:4326'}):
#     """converts a featurecollection to a geoPandas GeoDataFrame. Use this function only if all features have a geometry.
#
#     caveats:
#     Currently only supports non-geodesic (planar) geometries because of limitations in geoPandas and Leaflet. Geodesic geometries are simply interpreted as planar geometries.
#     FeatureCollections larger than memory are currently not supported. Consider splitting data and merging (server side) afterwards.
#
#     Args:
#         fc (ee.FeatureCollection) : the earth engine feature collection to convert.
#         crs (dictionary, optional) : the coordinate reference system in geopandas format. Defaults to {'init' :'epsg:4326'}
#
#     Returns:
#         gdf (geoPandas.GeoDataFrame or pandas.DataFrame) : the corresponding (geo)dataframe.
#
#     """
#     crs = {'init': 'epsg:4326'}
#
#     features = fc.getInfo()['features']
#     dictarr = []
#
#     for f in features:
#         # geodesic = ee.Feature(f).geometry().edgesAreGeodesics()
#         # if geodesic:
#         attr = f['properties']
#         attr['geometry'] = f['geometry']
#         attr['geometry']
#         dictarr.append(attr)
#
#     gdf = gpd.GeoDataFrame(dictarr)
#     gdf['geometry'] = list(map(lambda s: shapely.geometry.shape(s), gdf.geometry))
#     gdf.crs = crs
#     return gdf