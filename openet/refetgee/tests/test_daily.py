import ee
import pytest

from openet.refetgee import Daily
import openet.refetgee.units as units

# Eventually move to conftest.py or a separate file
# Fallon AgriMet site parameters
s_args = {
    'elev': 1208.5,
    'lat': 39.4575,
    'lon': -118.77388,
    'pair': 87.81876435813037,
    'pair_asce': 87.80710537212929,
    'zw': 3.0,
}

# Daily test parameters for 2015-07-01
# Ea was computed from q for both asce and refet methods
d_args = {
    'doy': 182,
    'ea': 1.2206674169951346,
    'eto_asce': 7.9422320475712835,
    'eto_refet': 7.9422320475712835,
    'etr_asce': 10.626087665395694,
    'etr_refet': 10.571314344056955,
    'etr_rso_simple': 10.628137858930051,
    'q': 0.008691370735727117,          # Computed from Ea from Tdew
    'q_asce': 0.008692530868140688,     # Computed from Ea from Tdew
    'rs': 674.07 * 0.041868,            # Conversion from Langleys to MJ m-2
    'rso': 31.565939444861765,
    'tdew': units._f2c(49.84),
    'tmin': units._f2c(66.65),
    'tmax': units._f2c(102.80),
    # 'tmean': f2c(84.725),
    'uz': 4.80 * 0.44704,               # Conversion from mph to m s-1
    'u2': 1.976111757722194,
}

constant_geom = ee.Geometry.Rectangle([0, 0, 10, 10], 'EPSG:32613', False)


# # Test full daily functions with positional inputs
# def test_refet_daily_input_positions():
#     refet = Daily(
#         ee.Image.constant(d_args['tmax']), ee.Image.constant(d_args['tmin']),
#         ee.Image.constant(d_args['ea']), ee.Image.constant(d_args['rs']),
#         ee.Image.constant(d_args['uz']), ee.Number(s_args['zw']),
#         ee.Number(s_args['elev']), ee.Number(s_args['lat']),
#         ee.Number(d_args['doy']), method='refet')
#     output = refet.etr\
#         .reduceRegion(ee.Reducer.first(), geometry=constant_geom, scale=1)\
#         .getInfo()
#     assert float(output['etr']) == pytest.approx(d_args['etr'])


# Test full daily calculations with keyword inputs
# Test surface, rso_type, and rso inputs
def test_refet_daily_surface_etr():
    refet = Daily(
        tmax=ee.Image.constant(d_args['tmax']),
        tmin=ee.Image.constant(d_args['tmin']),
        ea=ee.Image.constant(d_args['ea']),
        rs=ee.Image.constant(d_args['rs']),
        uz=ee.Image.constant(d_args['uz']), zw=ee.Number(s_args['zw']),
        elev=ee.Number(s_args['elev']), lat=ee.Number(s_args['lat']),
        doy=ee.Number(d_args['doy']), method='refet')
    output = refet.etr\
        .reduceRegion(ee.Reducer.first(), geometry=constant_geom, scale=1)\
        .getInfo()
    assert float(output['etr']) == pytest.approx(d_args['etr_refet'])


def test_refet_daily_surface_eto():
    refet = Daily(
        tmax=ee.Image.constant(d_args['tmax']),
        tmin=ee.Image.constant(d_args['tmin']),
        ea=ee.Image.constant(d_args['ea']),
        rs=ee.Image.constant(d_args['rs']),
        uz=ee.Image.constant(d_args['uz']), zw=ee.Number(s_args['zw']),
        elev=ee.Number(s_args['elev']), lat=ee.Number(s_args['lat']),
        doy=ee.Number(d_args['doy']), method='refet')
    output = refet.eto\
        .reduceRegion(ee.Reducer.first(), geometry=constant_geom, scale=1)\
        .getInfo()
    assert float(output['eto']) == pytest.approx(d_args['eto_refet'])


def test_refet_daily_rso_type_simple():
    refet = Daily(
        tmax=ee.Image.constant(d_args['tmax']),
        tmin=ee.Image.constant(d_args['tmin']),
        ea=ee.Image.constant(d_args['ea']),
        rs=ee.Image.constant(d_args['rs']),
        uz=ee.Image.constant(d_args['uz']), zw=ee.Number(s_args['zw']),
        elev=ee.Number(s_args['elev']), lat=ee.Number(s_args['lat']),
        doy=ee.Number(d_args['doy']), method='refet', rso_type='simple')
    output = refet.etr\
        .reduceRegion(ee.Reducer.first(), geometry=constant_geom, scale=1)\
        .getInfo()
    assert float(output['etr']) == pytest.approx(d_args['etr_rso_simple'])


def test_refet_daily_rso_type_array():
    refet = Daily(
        tmax=ee.Image.constant(d_args['tmax']),
        tmin=ee.Image.constant(d_args['tmin']),
        ea=ee.Image.constant(d_args['ea']),
        rs=ee.Image.constant(d_args['rs']),
        uz=ee.Image.constant(d_args['uz']), zw=ee.Number(s_args['zw']),
        elev=ee.Number(s_args['elev']), lat=ee.Number(s_args['lat']),
        doy=ee.Number(d_args['doy']), method='refet',
        rso_type='array', rso=ee.Number(d_args['rso']))
    output = refet.etr\
        .reduceRegion(ee.Reducer.first(), geometry=constant_geom, scale=1)\
        .getInfo()
    assert float(output['etr']) == pytest.approx(d_args['etr_refet'])


def test_refet_daily_rso_type_exception():
    with pytest.raises(ValueError):
        refet = Daily(
            tmax=ee.Image.constant(d_args['tmax']),
            tmin=ee.Image.constant(d_args['tmin']),
            ea=ee.Image.constant(d_args['ea']),
            rs=ee.Image.constant(d_args['rs']),
            uz=ee.Image.constant(d_args['uz']), zw=ee.Number(s_args['zw']),
            elev=ee.Number(s_args['elev']), lat=ee.Number(s_args['lat']),
            doy=ee.Number(d_args['doy']), rso_type='nonsense', method='refet')


def test_refet_daily_asce():
    refet = Daily(
        tmax=ee.Image.constant(d_args['tmax']),
        tmin=ee.Image.constant(d_args['tmin']),
        ea=ee.Image.constant(d_args['ea']),
        rs=ee.Image.constant(d_args['rs']),
        uz=ee.Image.constant(d_args['uz']), zw=ee.Number(s_args['zw']),
        elev=ee.Number(s_args['elev']), lat=ee.Number(s_args['lat']),
        doy=ee.Number(d_args['doy']), method='asce')
    output = refet.etr.reduceRegion(
        reducer=ee.Reducer.first(),
        geometry=ee.Geometry.Rectangle([0, 0, 10, 10], 'EPSG:32613', False),
        scale=1).getInfo()
    assert float(output['etr']) == pytest.approx(d_args['etr_asce'])


@pytest.mark.parametrize(
    'surface, expected',
    [['etr', d_args['etr_refet']],
     ['alfalfa', d_args['etr_refet']],
     ['tall', d_args['etr_refet']],
     ['eto', d_args['eto_refet']],
     ['grass', d_args['eto_refet']],
     ['short', d_args['eto_refet']]])
def test_refet_daily_etsz(surface, expected):
    refet = Daily(
        tmax=ee.Image.constant(d_args['tmax']),
        tmin=ee.Image.constant(d_args['tmin']),
        ea=ee.Image.constant(d_args['ea']),
        rs=ee.Image.constant(d_args['rs']),
        uz=ee.Image.constant(d_args['uz']), zw=ee.Number(s_args['zw']),
        elev=ee.Number(s_args['elev']), lat=ee.Number(s_args['lat']),
        doy=ee.Number(d_args['doy']), method='refet')
    output = refet.etsz(surface).rename(['etsz']).reduceRegion(
        reducer=ee.Reducer.first(),
        geometry=ee.Geometry.Rectangle([0, 0, 10, 10], 'EPSG:32613', False),
        scale=1).getInfo()
    assert float(output['etsz']) == pytest.approx(expected)


def test_refet_daily_gridmet_etr():
    """Generate a fake GRIDMET image from the test values"""
    gridmet_img = ee.Image.constant([
            d_args['tmax'] + 273.15, d_args['tmin'] + 273.15,
            d_args['q_asce'], d_args['rs'] / 0.0864, d_args['uz']])\
        .rename(['tmmx', 'tmmn', 'sph', 'srad', 'vs'])\
        .set('system:time_start', ee.Date('2015-07-01').millis())
    refet = Daily.gridmet(
        ee.Image(gridmet_img), elev=ee.Number(s_args['elev']),
        lat=ee.Number(s_args['lat']), zw=ee.Number(s_args['zw']),
        method='asce')
    output = refet.etr\
        .reduceRegion(ee.Reducer.first(), geometry=constant_geom, scale=1)\
        .getInfo()
    assert float(output['etr']) == pytest.approx(d_args['etr_asce'])


def test_refet_daily_maca_etr():
    """Generate a fake MACA image from the test values"""
    maca_img = ee.Image.constant([
            d_args['tmax'] + 273.15, d_args['tmin'] + 273.15,
            d_args['q_asce'], d_args['rs'] / 0.0864,
            d_args['uz'] / (2 ** 0.5), d_args['uz'] / (2 ** 0.5)])\
        .rename(['tasmax', 'tasmin', 'huss', 'rsds', 'uas', 'vas'])\
        .set('system:time_start', ee.Date('2015-07-01').millis())
    refet = Daily.maca(
        ee.Image(maca_img), elev=ee.Number(s_args['elev']),
        lat=ee.Number(s_args['lat']), zw=ee.Number(s_args['zw']),
        method='asce')
    output = refet.etr\
        .reduceRegion(ee.Reducer.first(), geometry=constant_geom, scale=1)\
        .getInfo()
    assert float(output['etr']) == pytest.approx(d_args['etr_asce'])


def test_refet_daily_nldas_etr():
    """Generate a fake NLDAS image from the test values

    Convert the test Rs from MJ m-2 d-1 to W m-2, then allocate half to each image

    """
    band_names = ['temperature', 'specific_humidity',
                  'shortwave_radiation', 'wind_u', 'wind_v']

    wind_u = d_args['uz'] / (2 ** 0.5)
    nldas_coll = ee.ImageCollection.fromImages([
        ee.Image.constant([d_args['tmin'], d_args['q_asce'],
                           0.0, wind_u, wind_u]) \
            .double().rename(band_names) \
            .set({'system:time_start': ee.Date('2015-07-01T00:00:00', 'UTC').millis()}),
        ee.Image.constant([d_args['tmax'], d_args['q_asce'],
                           d_args['rs'] / 0.0036, wind_u, wind_u]) \
            .double().rename(band_names) \
            .set({'system:time_start': ee.Date('2015-07-01T12:00:00', 'UTC').millis()})
    ])

    refet = Daily.nldas(
        nldas_coll, elev=ee.Number(s_args['elev']),
        lat=ee.Number(s_args['lat']), zw=ee.Number(s_args['zw']), method='asce')
    output = refet.etr\
        .reduceRegion(ee.Reducer.first(), geometry=constant_geom, scale=1)\
        .getInfo()
    assert float(output['etr']) == pytest.approx(d_args['etr_asce'])


def test_refet_daily_nldas_eto():
    """Compare values to previous calculations"""
    test_point = ee.Geometry.Point(-120.113, 36.336)

    nldas_coll = ee.ImageCollection('NASA/NLDAS/FORA0125_H002')\
        .filterDate('2017-07-01T06:00:00', '2017-07-02T06:00:00')
    output = Daily.nldas(nldas_coll, method='asce').eto \
        .reduceRegion(ee.Reducer.first(), geometry=test_point, scale=1)\
        .getInfo()['eto']

    expected = ee.Image('projects/eddi-noaa/nldas/daily/20170701')\
        .select(['ETo'])\
        .reduceRegion(ee.Reducer.first(), geometry=test_point, scale=1)\
        .getInfo()['ETo']

    assert output == pytest.approx(expected, rel=0.001)


def test_refet_daily_cfsv2_etr():
    """Generate a fake CFSv2 image from the test values

    Convert the test Rs from MJ m-2 d-1 to W m-2

    """
    band_names = [
        'Maximum_temperature_height_above_ground_6_Hour_Interval',
        'Minimum_temperature_height_above_ground_6_Hour_Interval',
        'Specific_humidity_height_above_ground',
        'Downward_Short-Wave_Radiation_Flux_surface_6_Hour_Average',
        'u-component_of_wind_height_above_ground',
        'v-component_of_wind_height_above_ground',
    ]

    wind_u = d_args['uz'] / (2 ** 0.5)

    cfsv2_coll = ee.ImageCollection.fromImages([
        ee.Image.constant([d_args['tmax'] + 273.15, d_args['tmin'] + 273.15,
                           d_args['q_asce'], d_args['rs'] / 0.0864,
                           wind_u, wind_u]) \
            .double().rename(band_names) \
            .set({'system:time_start': ee.Date('2015-07-01T00:00:00', 'UTC').millis()}),
        ee.Image.constant([d_args['tmax'] + 273.15, d_args['tmin'] + 273.15,
                           d_args['q_asce'], d_args['rs'] / 0.0864,
                           wind_u, wind_u]) \
            .double().rename(band_names) \
            .set({'system:time_start': ee.Date('2015-07-01T06:00:00', 'UTC').millis()}),
        ee.Image.constant([d_args['tmax'] + 273.15, d_args['tmin'] + 273.15,
                           d_args['q_asce'], d_args['rs'] / 0.0864,
                           wind_u, wind_u]) \
             .double().rename(band_names) \
             .set({'system:time_start': ee.Date('2015-07-01T12:00:00', 'UTC').millis()}),
        ee.Image.constant([d_args['tmax'] + 273.15, d_args['tmin'] + 273.15,
                           d_args['q_asce'], d_args['rs'] / 0.0864,
                           wind_u, wind_u]) \
             .double().rename(band_names) \
             .set({'system:time_start': ee.Date('2015-07-01T18:00:00', 'UTC').millis()}),
    ])

    refet = Daily.cfsv2(
        cfsv2_coll, elev=ee.Number(s_args['elev']),
        lat=ee.Number(s_args['lat']), zw=ee.Number(s_args['zw']), method='asce')
    output = refet.etr\
        .reduceRegion(ee.Reducer.first(), geometry=constant_geom, scale=1)\
        .getInfo()
    assert float(output['etr']) == pytest.approx(d_args['etr_asce'])


def test_refet_daily_rtma_etr():
    """Generate a fake RTMA image from the test values"""
    band_names = ['TMP', 'SPFH', 'WIND']

    rtma_coll = ee.ImageCollection.fromImages([
        ee.Image.constant([d_args['tmin'], d_args['q_asce'], d_args['uz']]) \
            .double().rename(band_names) \
            .set({'system:time_start': ee.Date('2015-07-01T00:00:00', 'UTC').millis()}),
        ee.Image.constant([d_args['tmax'], d_args['q_asce'], d_args['uz']]) \
            .double().rename(band_names) \
            .set({'system:time_start': ee.Date('2015-07-01T12:00:00', 'UTC').millis()})
    ])

    refet = Daily.rtma(
        rtma_coll,
        rs=ee.Image.constant(d_args['rs']),
        elev=ee.Number(s_args['elev']), lat=ee.Number(s_args['lat']),
        zw=ee.Number(s_args['zw']), method='asce')
    output = refet.etr\
        .reduceRegion(ee.Reducer.first(), geometry=constant_geom, scale=1)\
        .getInfo()
    assert float(output['etr']) == pytest.approx(d_args['etr_asce'])


# TODO: Add a test for using the default Rs when one is not provided
