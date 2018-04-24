import ee
import pytest

from geerefet.daily import Daily
import geerefet.units as units

ee.Initialize()


# Eventually move to conftest.py or a separate file
# Fallon AgriMet site parameters
s_args = {
    'elev': 1208.5,
    'lat': 39.4575,
    'lon': -118.77388,
    # DEADBEEF
    # 'lat': units._deg2rad(39.4575),
    # 'lon': units._deg2rad(-118.77388),
    'zw': 3.0,
}

# Daily test parameters for 2015-07-01
d_args = {
    'doy': 182,
    'ea': 1.2206674169951346,
    'ea_asce': 1.2205053588697359,
    'eto': 7.9422320475712835,
    'etr': 10.571314344056955,
    'etr_asce': 10.62628103954038,
    'etr_rso_simple': 10.628137858930051,
    'q': 0.008691370735727117,
    'rs': 674.07 * 0.041868,  # Conversion from Langleys to MJ m-2
    'rso': 31.565939444861765,
    'tdew': units._f2c(49.84),
    'tmin': units._f2c(66.65),
    'tmax': units._f2c(102.80),
    # 'tmean': f2c(84.725),
    'uz': 4.80 * 0.44704,  # Conversion from mph to m s-1
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
#     output = refet.etr()\
#         .reduceRegion(ee.Reducer.first(), geometry=constant_geom, scale=1)\
#         .getInfo()
#     assert float(output['etr']) == pytest.approx(d_args['etr'])


# Test full daily calculations with keyword inputs
# Test surface, rso_type, and rso inputs
def test_refet_daily_surface_etr():
    refet = Daily(
        tmax=ee.Image.constant(d_args['tmax']),
        tmin=ee.Image.constant(d_args['tmin']),
        ea=ee.Image.constant(d_args['ea']), rs=ee.Image.constant(d_args['rs']),
        uz=ee.Image.constant(d_args['uz']), zw=ee.Number(s_args['zw']),
        elev=ee.Number(s_args['elev']), lat=ee.Number(s_args['lat']),
        doy=ee.Number(d_args['doy']), method='refet')
    output = refet.etr()\
        .reduceRegion(ee.Reducer.first(), geometry=constant_geom, scale=1)\
        .getInfo()
    assert float(output['etr']) == pytest.approx(d_args['etr'])


def test_refet_daily_surface_eto():
    refet = Daily(
        tmax=ee.Image.constant(d_args['tmax']),
        tmin=ee.Image.constant(d_args['tmin']),
        ea=ee.Image.constant(d_args['ea']), rs=ee.Image.constant(d_args['rs']),
        uz=ee.Image.constant(d_args['uz']), zw=ee.Number(s_args['zw']),
        elev=ee.Number(s_args['elev']), lat=ee.Number(s_args['lat']),
        doy=ee.Number(d_args['doy']), method='refet')
    output = refet.eto()\
        .reduceRegion(ee.Reducer.first(), geometry=constant_geom, scale=1)\
        .getInfo()
    assert float(output['eto']) == pytest.approx(d_args['eto'])


def test_refet_daily_rso_type_simple():
    refet = Daily(
        tmax=ee.Image.constant(d_args['tmax']),
        tmin=ee.Image.constant(d_args['tmin']),
        ea=ee.Image.constant(d_args['ea']), rs=ee.Image.constant(d_args['rs']),
        uz=ee.Image.constant(d_args['uz']), zw=ee.Number(s_args['zw']),
        elev=ee.Number(s_args['elev']), lat=ee.Number(s_args['lat']),
        doy=ee.Number(d_args['doy']), method='refet', rso_type='simple')
    output = refet.etr()\
        .reduceRegion(ee.Reducer.first(), geometry=constant_geom, scale=1)\
        .getInfo()
    assert float(output['etr']) == pytest.approx(d_args['etr_rso_simple'])


def test_refet_daily_rso_type_array():
    refet = Daily(
        tmax=ee.Image.constant(d_args['tmax']),
        tmin=ee.Image.constant(d_args['tmin']),
        ea=ee.Image.constant(d_args['ea']), rs=ee.Image.constant(d_args['rs']),
        uz=ee.Image.constant(d_args['uz']), zw=ee.Number(s_args['zw']),
        elev=ee.Number(s_args['elev']), lat=ee.Number(s_args['lat']),
        doy=ee.Number(d_args['doy']), method='refet',
        rso_type='array', rso=ee.Number(d_args['rso']))
    output = refet.etr()\
        .reduceRegion(ee.Reducer.first(), geometry=constant_geom, scale=1)\
        .getInfo()
    assert float(output['etr']) == pytest.approx(d_args['etr'])


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
        ea=ee.Image.constant(d_args['ea_asce']),
        rs=ee.Image.constant(d_args['rs']),
        uz=ee.Image.constant(d_args['uz']), zw=ee.Number(s_args['zw']),
        elev=ee.Number(s_args['elev']), lat=ee.Number(s_args['lat']),
        doy=ee.Number(d_args['doy']), method='asce')
    output = refet.etr().reduceRegion(
        reducer=ee.Reducer.first(),
        geometry=ee.Geometry.Rectangle([0, 0, 10, 10], 'EPSG:32613', False),
        scale=1).getInfo()
    assert float(output['etr']) == pytest.approx(d_args['etr_asce'])


def test_refet_daily_gridmet():
    # Convert input units to GRIDMET units
    # Overriding GRIDMET windspeed height from 10m to 3m
    gridmet_img = ee.Image.constant([
            d_args['tmax'] + 273.15, d_args['tmin'] + 273.15,
            d_args['q'], d_args['rs'] / 0.0864, d_args['uz']])\
        .rename(['tmmx', 'tmmn', 'sph', 'srad', 'vs'])\
        .set('system:time_start', ee.Date('2015-07-01').millis())
    refet = Daily.gridmet(
        ee.Image(gridmet_img), elev=ee.Number(s_args['elev']),
        lat=ee.Number(s_args['lat']), zw=ee.Number(s_args['zw']), method='asce')
    output = refet.etr()\
        .reduceRegion(ee.Reducer.first(), geometry=constant_geom, scale=1)\
        .getInfo()
    assert float(output['etr']) == pytest.approx(d_args['etr_asce'])
