import ee
import pytest

from openet.refetgee import Hourly
import openet.refetgee.units as units
import utils

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
    'tdew': units.f2c(49.36),
    'time': 18.0,
    'time_mid': 18.5,
    'tmean': units.f2c(91.80),
    'uz': 3.33 * 0.44704,               # Conversion from mph to m s-1
    'u2': 1.3709275319197722,
}

# Compute wind x and y components
h_args['uz_x'] = h_args['uz'] / (2 ** 0.5)
h_args['uz_y'] = h_args['uz'] / (2 ** 0.5)
# h_args['uz_x'] = h_args['uz']
# h_args['uz_y'] = 0.0


# Test full hourly functions with positional inputs
def test_refet_hourly_input_positions():
    refet = Hourly(
        ee.Image.constant(h_args['tmean']), ee.Image.constant(h_args['ea']),
        ee.Image.constant(h_args['rs']), ee.Image.constant(h_args['uz']),
        ee.Image.constant(s_args['zw']), ee.Number(s_args['elev']),
        ee.Number(s_args['lat']), ee.Number(s_args['lon']),
        ee.Number(h_args['doy']), ee.Number(h_args['time']), method='refet'
    )

    output = utils.constant_image_value(refet.etr)

    assert float(output['etr']) == pytest.approx(h_args['etr_refet'])


# Test full hourly calculations with keyword inputs
def test_refet_hourly_default_method_etr():
    refet = Hourly(
        tmean=ee.Image.constant(h_args['tmean']), ea=ee.Image.constant(h_args['ea']),
        rs=ee.Image.constant(h_args['rs']), uz=ee.Image.constant(h_args['uz']),
        zw=ee.Image.constant(s_args['zw']), elev=ee.Number(s_args['elev']),
        lat=ee.Number(s_args['lat']), lon=ee.Number(s_args['lon']),
        doy=ee.Number(h_args['doy']), time=ee.Number(h_args['time']),
    )

    output = utils.constant_image_value(refet.etr)

    assert float(output['etr']) == pytest.approx(h_args['etr_asce'])


def test_refet_hourly_asce_method_etr():
    refet = Hourly(
        tmean=ee.Image.constant(h_args['tmean']), ea=ee.Image.constant(h_args['ea']),
        rs=ee.Image.constant(h_args['rs']), uz=ee.Image.constant(h_args['uz']),
        zw=ee.Number(s_args['zw']), elev=ee.Number(s_args['elev']),
        lat=ee.Number(s_args['lat']), lon=ee.Number(s_args['lon']),
        doy=ee.Number(h_args['doy']), time=ee.Number(h_args['time']), method='asce',
    )

    output = utils.constant_image_value(refet.etr)

    assert float(output['etr']) == pytest.approx(h_args['etr_asce'])


def test_refet_hourly_refet_method_etr():
    refet = Hourly(
        tmean=ee.Image.constant(h_args['tmean']), ea=ee.Image.constant(h_args['ea']),
        rs=ee.Image.constant(h_args['rs']), uz=ee.Image.constant(h_args['uz']),
        zw=ee.Image.constant(s_args['zw']), elev=ee.Number(s_args['elev']),
        lat=ee.Number(s_args['lat']), lon=ee.Number(s_args['lon']),
        doy=ee.Number(h_args['doy']), time=ee.Number(h_args['time']), method='refet',
    )

    output = utils.constant_image_value(refet.etr)

    assert float(output['etr']) == pytest.approx(h_args['etr_refet'])


def test_refet_hourly_default_method_eto():
    refet = Hourly(
        tmean=ee.Image.constant(h_args['tmean']), ea=ee.Image.constant(h_args['ea']),
        rs=ee.Image.constant(h_args['rs']), uz=ee.Image.constant(h_args['uz']),
        zw=ee.Image.constant(s_args['zw']), elev=ee.Number(s_args['elev']),
        lat=ee.Number(s_args['lat']), lon=ee.Number(s_args['lon']),
        doy=ee.Number(h_args['doy']), time=ee.Number(h_args['time']),
    )

    output = utils.constant_image_value(refet.eto)

    assert float(output['eto']) == pytest.approx(h_args['eto_asce'])


def test_refet_hourly_asce_method_eto():
    refet = Hourly(
        tmean=ee.Image.constant(h_args['tmean']), ea=ee.Image.constant(h_args['ea']),
        rs=ee.Image.constant(h_args['rs']), uz=ee.Image.constant(h_args['uz']),
        zw=ee.Image.constant(s_args['zw']), elev=ee.Number(s_args['elev']),
        lat=ee.Number(s_args['lat']), lon=ee.Number(s_args['lon']),
        doy=ee.Number(h_args['doy']), time=ee.Number(h_args['time']), method='asce',
    )

    output = utils.constant_image_value(refet.eto)

    assert float(output['eto']) == pytest.approx(h_args['eto_asce'])


def test_refet_hourly_asce_method_eto():
    refet = Hourly(
        tmean=ee.Image.constant(h_args['tmean']), ea=ee.Image.constant(h_args['ea']),
        rs=ee.Image.constant(h_args['rs']), uz=ee.Image.constant(h_args['uz']),
        zw=ee.Image.constant(s_args['zw']), elev=ee.Number(s_args['elev']),
        lat=ee.Number(s_args['lat']), lon=ee.Number(s_args['lon']),
        doy=ee.Number(h_args['doy']), time=ee.Number(h_args['time']), method='refet',
    )

    output = utils.constant_image_value(refet.eto)

    assert float(output['eto']) == pytest.approx(h_args['eto_refet'])


@pytest.mark.parametrize(
    'surface, expected',
    [
        ['etr', h_args['etr_refet']],
        ['alfalfa', h_args['etr_refet']],
        ['tall', h_args['etr_refet']],
        ['eto', h_args['eto_refet']],
        ['grass', h_args['eto_refet']],
        ['short', h_args['eto_refet']],
    ]
)
def test_refet_daily_etsz(surface, expected):
    refet = Hourly(
        tmean=ee.Image.constant(h_args['tmean']), ea=ee.Image.constant(h_args['ea']),
        rs=ee.Image.constant(h_args['rs']), uz=ee.Image.constant(h_args['uz']),
        zw=ee.Image.constant(s_args['zw']), elev=ee.Number(s_args['elev']),
        lat=ee.Number(s_args['lat']), lon=ee.Number(s_args['lon']),
        doy=ee.Number(h_args['doy']), time=ee.Number(h_args['time']), method='refet',
    )

    output = utils.constant_image_value(refet.etsz(surface).rename(['etsz']))

    assert float(output['etsz']) == pytest.approx(expected)


def test_refet_hourly_nldas_etr():
    """Generate a mock NLDAS image from the test values"""
    time_start = ee.Date(f'2015-07-01T{int(h_args["time"])}:00:00', 'UTC').millis()

    input_img = (
        ee.Image.constant([
            h_args['tmean'], h_args['q_asce'], h_args['rs'] / 0.0036,
            h_args['uz_x'], h_args['uz_y']
        ])
        .rename(['temperature', 'specific_humidity', 'shortwave_radiation', 'wind_u', 'wind_v'])
        .set('system:time_start', time_start)
    )

    refet = Hourly.nldas(
        ee.Image(input_img), elev=ee.Number(s_args['elev']),
        lat=ee.Number(s_args['lat']), lon=ee.Number(s_args['lon']),
        zw=ee.Number(s_args['zw']), method='asce',
    )

    output = utils.constant_image_value(refet.etr)

    assert float(output['etr']) == pytest.approx(h_args['etr_asce'])


def test_refet_hourly_rtma_etr():
    """Generate a mock RTMA image from the test values"""
    time_start = ee.Date(f'2015-07-01T{int(h_args["time"])}:00:00', 'UTC').millis()

    input_img = (
        ee.Image.constant([h_args['tmean'], h_args['q_asce'], h_args['uz']])
        .rename(['TMP', 'SPFH', 'WIND'])
        .set('system:time_start', time_start)
    )

    refet = Hourly.rtma(
        ee.Image(input_img), rs=ee.Image.constant(h_args['rs']),
        elev=ee.Number(s_args['elev']),
        lat=ee.Number(s_args['lat']), lon=ee.Number(s_args['lon']),
        zw=ee.Number(s_args['zw']), method='asce',
    )

    output = utils.constant_image_value(refet.etr)

    assert float(output['etr']) == pytest.approx(h_args['etr_asce'])


def test_refet_hourly_era5_etr():
    """Generate a mock ERA5 image from the test values"""
    time_start = ee.Date(f'2015-07-01T{int(h_args["time"])}:00:00', 'UTC').millis()

    input_img = (
        ee.Image.constant([
            h_args['tmean'] + 273.15, h_args['tdew'] + 273.15,
            h_args['rs'] * 1000000, h_args['uz_x'], h_args['uz_y']
        ])
        .rename(['temperature_2m', 'dewpoint_temperature_2m',
                 'surface_solar_radiation_downwards',
                 'u_component_of_wind_10m', 'v_component_of_wind_10m'])
        .set('system:time_start', time_start)
    )

    refet = Hourly.era5(
        ee.Image(input_img), elev=ee.Number(s_args['elev']),
        lat=ee.Number(s_args['lat']), lon=ee.Number(s_args['lon']),
        zw=ee.Number(s_args['zw']), method='asce',
    )

    output = utils.constant_image_value(refet.etr)

    assert float(output['etr']) == pytest.approx(h_args['etr_asce'])


def test_refet_hourly_era5_land_etr():
    """Generate a mock ERA5-Land image from the test values"""
    time_start = ee.Date(f'2015-07-01T{int(h_args["time"])}:00:00', 'UTC').millis()

    input_img = (
        ee.Image.constant([
            h_args['tmean'] + 273.15, h_args['tdew'] + 273.15,
            h_args['rs'] * 1000000, h_args['uz_x'], h_args['uz_y']
        ])
        .rename(['temperature_2m', 'dewpoint_temperature_2m',
                 'surface_solar_radiation_downwards_hourly',
                 'u_component_of_wind_10m', 'v_component_of_wind_10m'])
        .set('system:time_start', time_start)
    )

    refet = Hourly.era5_land(
        ee.Image(input_img), elev=ee.Number(s_args['elev']),
        lat=ee.Number(s_args['lat']), lon=ee.Number(s_args['lon']),
        zw=ee.Number(s_args['zw']), method='asce',
    )

    output = utils.constant_image_value(refet.etr)

    assert float(output['etr']) == pytest.approx(h_args['etr_asce'])


def test_refet_hourly_era5_land_fill_edge_cells():
    """Check that the fill_edge_cells flag works for an edge cell along the coast of England"""
    input_img = ee.Image('ECMWF/ERA5_LAND/HOURLY/20150701T20')
    output = utils.point_image_value(Hourly.era5_land(input_img, fill_edge_cells=False).etr, xy=[0.0, 50.7])
    assert output['etr'] is None
    output = utils.point_image_value(Hourly.era5_land(input_img, fill_edge_cells=True).etr, xy=[0.0, 50.7])
    assert output['etr'] is not None


# TODO: Add a test for using the default Rs when one is not provided
