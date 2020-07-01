import ee
import pytest

from openet.refetgee import Hourly
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

# Hourly test parameters for 2015-07-01 18:00 UTC (11:00 AM PDT)
# DEADBEEF - The eto_refet value changes...
h_args = {
    'doy': 182,
    'ea': 1.1990099614301906,
    'es': 5.09318785259078,
    'eto_refet': 0.6068613650177562,
    'eto_asce': 0.6063515410076268,
    'etr_refet': 0.7201865213918281,
    'etr_asce': 0.7196369609713682,
    'q': 0.008536365803069757,          # Computed from Ea from Tdew
    'q_asce': 0.00853750513849305,      # Computed from Ea from Tdew
    'ra': 4.30824147948541,
    'rnl': 0.22897874401150786,
    'rs': 61.16 * 0.041868,             # Conversion from Langleys to MJ m-2
    'tdew': units._f2c(49.36),
    'time': 18.0,
    'time_mid': 18.5,
    'tmean': units._f2c(91.80),
    'uz': 3.33 * 0.44704,               # Conversion from mph to m s-1
    'u2': 1.3709275319197722,
}

constant_geom = ee.Geometry.Rectangle([0, 0, 10, 10], 'EPSG:32613', False)


# Test full hourly functions with positional inputs
def test_refet_hourly_input_positions():
    refet = Hourly(
        ee.Image.constant(h_args['tmean']), ee.Image.constant(h_args['ea']),
        ee.Image.constant(h_args['rs']), ee.Image.constant(h_args['uz']),
        ee.Image.constant(s_args['zw']), ee.Number(s_args['elev']),
        ee.Number(s_args['lat']), ee.Number(s_args['lon']),
        ee.Number(h_args['doy']), ee.Number(h_args['time']), method='refet')
    output = refet.etr\
        .reduceRegion(ee.Reducer.first(), geometry=constant_geom, scale=1)\
        .getInfo()
    assert float(output['etr']) == pytest.approx(h_args['etr_refet'])


# Test full hourly calculations with keyword inputs
def test_refet_hourly_default_method_etr():
    refet = Hourly(
        tmean=ee.Image.constant(h_args['tmean']),
        ea=ee.Image.constant(h_args['ea']),
        rs=ee.Image.constant(h_args['rs']), uz=ee.Image.constant(h_args['uz']),
        zw=ee.Image.constant(s_args['zw']), elev=ee.Number(s_args['elev']),
        lat=ee.Number(s_args['lat']), lon=ee.Number(s_args['lon']),
        doy=ee.Number(h_args['doy']), time=ee.Number(h_args['time']))
    output = refet.etr\
        .reduceRegion(ee.Reducer.first(), geometry=constant_geom, scale=1)\
        .getInfo()
    assert float(output['etr']) == pytest.approx(h_args['etr_asce'])


def test_refet_hourly_asce_method_etr():
    refet = Hourly(
        tmean=ee.Image.constant(h_args['tmean']),
        ea=ee.Image.constant(h_args['ea']),
        rs=ee.Image.constant(h_args['rs']), uz=ee.Image.constant(h_args['uz']),
        zw=ee.Number(s_args['zw']), elev=ee.Number(s_args['elev']),
        lat=ee.Number(s_args['lat']), lon=ee.Number(s_args['lon']),
        doy=ee.Number(h_args['doy']), time=ee.Number(h_args['time']),
        method='asce')
    output = refet.etr\
        .reduceRegion(ee.Reducer.first(), geometry=constant_geom, scale=1)\
        .getInfo()
    assert float(output['etr']) == pytest.approx(h_args['etr_asce'])


def test_refet_hourly_refet_method_etr():
    refet = Hourly(
        tmean=ee.Image.constant(h_args['tmean']),
        ea=ee.Image.constant(h_args['ea']),
        rs=ee.Image.constant(h_args['rs']), uz=ee.Image.constant(h_args['uz']),
        zw=ee.Image.constant(s_args['zw']), elev=ee.Number(s_args['elev']),
        lat=ee.Number(s_args['lat']), lon=ee.Number(s_args['lon']),
        doy=ee.Number(h_args['doy']), time=ee.Number(h_args['time']),
        method='refet')
    output = refet.etr\
        .reduceRegion(ee.Reducer.first(), geometry=constant_geom, scale=1)\
        .getInfo()
    assert float(output['etr']) == pytest.approx(h_args['etr_refet'])


def test_refet_hourly_default_method_eto():
    refet = Hourly(
        tmean=ee.Image.constant(h_args['tmean']),
        ea=ee.Image.constant(h_args['ea']),
        rs=ee.Image.constant(h_args['rs']), uz=ee.Image.constant(h_args['uz']),
        zw=ee.Image.constant(s_args['zw']), elev=ee.Number(s_args['elev']),
        lat=ee.Number(s_args['lat']), lon=ee.Number(s_args['lon']),
        doy=ee.Number(h_args['doy']),
        time=ee.Number(h_args['time']))
    output = refet.eto\
        .reduceRegion(ee.Reducer.first(), geometry=constant_geom, scale=1)\
        .getInfo()
    assert float(output['eto']) == pytest.approx(h_args['eto_asce'])


def test_refet_hourly_asce_method_eto():
    refet = Hourly(
        tmean=ee.Image.constant(h_args['tmean']),
        ea=ee.Image.constant(h_args['ea']),
        rs=ee.Image.constant(h_args['rs']), uz=ee.Image.constant(h_args['uz']),
        zw=ee.Image.constant(s_args['zw']), elev=ee.Number(s_args['elev']),
        lat=ee.Number(s_args['lat']), lon=ee.Number(s_args['lon']),
        doy=ee.Number(h_args['doy']), time=ee.Number(h_args['time']),
        method='asce')
    output = refet.eto\
        .reduceRegion(ee.Reducer.first(), geometry=constant_geom, scale=1)\
        .getInfo()
    assert float(output['eto']) == pytest.approx(h_args['eto_asce'])


def test_refet_hourly_asce_method_eto():
    refet = Hourly(
        tmean=ee.Image.constant(h_args['tmean']),
        ea=ee.Image.constant(h_args['ea']),
        rs=ee.Image.constant(h_args['rs']), uz=ee.Image.constant(h_args['uz']),
        zw=ee.Image.constant(s_args['zw']), elev=ee.Number(s_args['elev']),
        lat=ee.Number(s_args['lat']), lon=ee.Number(s_args['lon']),
        doy=ee.Number(h_args['doy']), time=ee.Number(h_args['time']),
        method='refet')
    output = refet.eto\
        .reduceRegion(ee.Reducer.first(), geometry=constant_geom, scale=1)\
        .getInfo()
    assert float(output['eto']) == pytest.approx(h_args['eto_refet'])


@pytest.mark.parametrize(
    'surface, expected',
    [['etr', h_args['etr_refet']],
     ['alfalfa', h_args['etr_refet']],
     ['tall', h_args['etr_refet']],
     ['eto', h_args['eto_refet']],
     ['grass', h_args['eto_refet']],
     ['short', h_args['eto_refet']]])
def test_refet_daily_etsz(surface, expected):
    refet = Hourly(
        tmean=ee.Image.constant(h_args['tmean']),
        ea=ee.Image.constant(h_args['ea']),
        rs=ee.Image.constant(h_args['rs']), uz=ee.Image.constant(h_args['uz']),
        zw=ee.Image.constant(s_args['zw']), elev=ee.Number(s_args['elev']),
        lat=ee.Number(s_args['lat']), lon=ee.Number(s_args['lon']),
        doy=ee.Number(h_args['doy']), time=ee.Number(h_args['time']),
        method='refet')
    output = refet.etsz(surface).rename(['etsz']).reduceRegion(
        reducer=ee.Reducer.first(),
        geometry=ee.Geometry.Rectangle([0, 0, 10, 10], 'EPSG:32613', False),
        scale=1).getInfo()
    assert float(output['etsz']) == pytest.approx(expected)


def test_refet_hourly_nldas_etr():
    """Generate a fake NLDAS image from the test values"""
    nldas_time = ee.Date(
        '2015-07-01T{}:00:00'.format(int(h_args['time'])), 'UTC').millis()
    wind_u = h_args['uz'] / (2 ** 0.5)
    nldas_img = ee.Image.constant([
            h_args['tmean'], h_args['q_asce'], h_args['rs'] / 0.0036,
            wind_u, wind_u])\
        .rename(['temperature', 'specific_humidity', 'shortwave_radiation',
                 'wind_u', 'wind_v'])\
        .set('system:time_start', nldas_time)
    refet = Hourly.nldas(
        ee.Image(nldas_img), elev=ee.Number(s_args['elev']),
        lat=ee.Number(s_args['lat']), lon=ee.Number(s_args['lon']),
        zw=ee.Number(s_args['zw']), method='asce')
    output = refet.etr\
        .reduceRegion(ee.Reducer.first(), geometry=constant_geom, scale=1)\
        .getInfo()
    assert float(output['etr']) == pytest.approx(h_args['etr_asce'])
