import ee
import pytest

from openet.refetgee import Daily
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

# Daily test parameters for 2015-07-01
# Ea was computed from q for both asce and refet methods
d_args = {
    'doy': 182,
    'ea': 1.2206674169951346,
    'eto_asce': 7.9422320475712835,
    'eto_refet': 7.9422320475712835,
    'etr_asce': 10.626087665395694,
    'etr_refet': 10.571314344056955,
    'etw_refet': 6.242411580074248,
    'eto_fs1': 4.373303512213786,
    'eto_fs2': 3.5689285353574958,
    'etr_rso_simple': 10.628137858930051,
    'q': 0.008691370735727117,          # Computed from Ea from Tdew
    'q_asce': 0.008692530868140688,     # Computed from Ea from Tdew
    'pet_hargreaves': 8.247962376780558,
    'rs': 674.07 * 0.041868,            # Conversion from Langleys to MJ m-2
    'rso': 31.565939444861765,
    'tdew': units.f2c(49.84),
    'tmin': units.f2c(66.65),
    'tmax': units.f2c(102.80),
    # 'tmean': f2c(84.725),
    'uz': 4.80 * 0.44704,               # Conversion from mph to m s-1
    'u2': 1.976111757722194,
}

# Compute wind x and y components
d_args['uz_x'] = d_args['uz'] / (2 ** 0.5)
d_args['uz_y'] = d_args['uz'] / (2 ** 0.5)
# d_args['uz_x'] = d_args['uz']
# d_args['uz_y'] = 0.0


# # Test full daily functions with positional inputs
# def test_refet_daily_input_positions():
#     refet = Daily(
#         ee.Image.constant(d_args['tmax']), ee.Image.constant(d_args['tmin']),
#         ee.Image.constant(d_args['ea']), ee.Image.constant(d_args['rs']),
#         ee.Image.constant(d_args['uz']), ee.Number(s_args['zw']),
#         ee.Number(s_args['elev']), ee.Number(s_args['lat']),
#         ee.Number(d_args['doy']), method='refet',
#     )
#     output = utils.constant_image_value(refet.etr)
#     assert float(output['etr']) == pytest.approx(d_args['etr'])


# Test full daily calculations with keyword inputs
# Test surface, rso_type, and rso inputs
def test_refet_daily_etr():
    refet = Daily(
        tmax=ee.Image.constant(d_args['tmax']), tmin=ee.Image.constant(d_args['tmin']),
        ea=ee.Image.constant(d_args['ea']), rs=ee.Image.constant(d_args['rs']),
        uz=ee.Image.constant(d_args['uz']), zw=ee.Number(s_args['zw']),
        elev=ee.Number(s_args['elev']), lat=ee.Number(s_args['lat']),
        doy=ee.Number(d_args['doy']), method='refet',
    )

    output = utils.constant_image_value(refet.etr)

    assert float(output['etr']) == pytest.approx(d_args['etr_refet'])


def test_refet_daily_eto():
    refet = Daily(
        tmax=ee.Image.constant(d_args['tmax']), tmin=ee.Image.constant(d_args['tmin']),
        ea=ee.Image.constant(d_args['ea']), rs=ee.Image.constant(d_args['rs']),
        uz=ee.Image.constant(d_args['uz']), zw=ee.Number(s_args['zw']),
        elev=ee.Number(s_args['elev']), lat=ee.Number(s_args['lat']),
        doy=ee.Number(d_args['doy']), method='refet',
    )

    output = utils.constant_image_value(refet.eto)

    assert float(output['eto']) == pytest.approx(d_args['eto_refet'])


def test_refet_daily_etw():
    refet = Daily(
        tmax=ee.Image.constant(d_args['tmax']), tmin=ee.Image.constant(d_args['tmin']),
        ea=ee.Image.constant(d_args['ea']), rs=ee.Image.constant(d_args['rs']),
        uz=ee.Image.constant(d_args['uz']), zw=ee.Number(s_args['zw']),
        elev=ee.Number(s_args['elev']), lat=ee.Number(s_args['lat']),
        doy=ee.Number(d_args['doy']), method='refet',
    )

    output = utils.constant_image_value(refet.etw)

    assert float(output['etw']) == pytest.approx(d_args['etw_refet'])


def test_refet_daily_eto_fs1():
    refet = Daily(
        tmax=ee.Image.constant(d_args['tmax']), tmin=ee.Image.constant(d_args['tmin']),
        ea=ee.Image.constant(d_args['ea']), rs=ee.Image.constant(d_args['rs']),
        uz=ee.Image.constant(d_args['uz']), zw=ee.Number(s_args['zw']),
        elev=ee.Number(s_args['elev']), lat=ee.Number(s_args['lat']),
        doy=ee.Number(d_args['doy']), method='refet',
    )

    output = utils.constant_image_value(refet.eto_fs1)

    assert float(output['eto_fs1']) == pytest.approx(d_args['eto_fs1'])


def test_refet_daily_eto_fs2():
    refet = Daily(
        tmax=ee.Image.constant(d_args['tmax']), tmin=ee.Image.constant(d_args['tmin']),
        ea=ee.Image.constant(d_args['ea']), rs=ee.Image.constant(d_args['rs']),
        uz=ee.Image.constant(d_args['uz']), zw=ee.Number(s_args['zw']),
        elev=ee.Number(s_args['elev']), lat=ee.Number(s_args['lat']),
        doy=ee.Number(d_args['doy']), method='refet',
    )

    output = utils.constant_image_value(refet.eto_fs2)

    assert float(output['eto_fs2']) == pytest.approx(d_args['eto_fs2'])


def test_refet_daily_pet_hargreaves():
    refet = Daily(
        tmax=ee.Image.constant(d_args['tmax']), tmin=ee.Image.constant(d_args['tmin']),
        ea=ee.Image.constant(0), rs=ee.Image.constant(0),
        uz=ee.Image.constant(0), zw=ee.Number(s_args['zw']),
        elev=ee.Number(0), lat=ee.Number(s_args['lat']),
        doy=ee.Number(d_args['doy']), method='asce',
    )

    output = utils.constant_image_value(refet.pet_hargreaves)

    assert float(output['pet_hargreaves']) == pytest.approx(d_args['pet_hargreaves'])


def test_refet_daily_rso_type_simple():
    refet = Daily(
        tmax=ee.Image.constant(d_args['tmax']), tmin=ee.Image.constant(d_args['tmin']),
        ea=ee.Image.constant(d_args['ea']), rs=ee.Image.constant(d_args['rs']),
        uz=ee.Image.constant(d_args['uz']), zw=ee.Number(s_args['zw']),
        elev=ee.Number(s_args['elev']), lat=ee.Number(s_args['lat']),
        doy=ee.Number(d_args['doy']), method='refet', rso_type='simple',
    )

    output = utils.constant_image_value(refet.etr)

    assert float(output['etr']) == pytest.approx(d_args['etr_rso_simple'])


def test_refet_daily_rso_type_array():
    refet = Daily(
        tmax=ee.Image.constant(d_args['tmax']), tmin=ee.Image.constant(d_args['tmin']),
        ea=ee.Image.constant(d_args['ea']), rs=ee.Image.constant(d_args['rs']),
        uz=ee.Image.constant(d_args['uz']), zw=ee.Number(s_args['zw']),
        elev=ee.Number(s_args['elev']), lat=ee.Number(s_args['lat']),
        doy=ee.Number(d_args['doy']), method='refet',
        rso_type='array', rso=ee.Number(d_args['rso']),
    )

    output = utils.constant_image_value(refet.etr)

    assert float(output['etr']) == pytest.approx(d_args['etr_refet'])


def test_refet_daily_rso_type_exception():
    with pytest.raises(ValueError):
        Daily(
            tmax=ee.Image.constant(d_args['tmax']),
            tmin=ee.Image.constant(d_args['tmin']),
            ea=ee.Image.constant(d_args['ea']),
            rs=ee.Image.constant(d_args['rs']),
            uz=ee.Image.constant(d_args['uz']), zw=ee.Number(s_args['zw']),
            elev=ee.Number(s_args['elev']), lat=ee.Number(s_args['lat']),
            doy=ee.Number(d_args['doy']), rso_type='nonsense', method='refet',
        )


def test_refet_daily_etr_asce():
    refet = Daily(
        tmax=ee.Image.constant(d_args['tmax']), tmin=ee.Image.constant(d_args['tmin']),
        ea=ee.Image.constant(d_args['ea']), rs=ee.Image.constant(d_args['rs']),
        uz=ee.Image.constant(d_args['uz']), zw=ee.Number(s_args['zw']),
        elev=ee.Number(s_args['elev']), lat=ee.Number(s_args['lat']),
        doy=ee.Number(d_args['doy']), method='asce',
    )

    output = utils.constant_image_value(refet.etr)

    assert float(output['etr']) == pytest.approx(d_args['etr_asce'])


@pytest.mark.parametrize(
    'surface, expected',
    [
        ['etr', d_args['etr_refet']],
        ['alfalfa', d_args['etr_refet']],
        ['tall', d_args['etr_refet']],
        ['eto', d_args['eto_refet']],
        ['grass', d_args['eto_refet']],
        ['short', d_args['eto_refet']],
    ]
)
def test_refet_daily_etsz(surface, expected):
    refet = Daily(
        tmax=ee.Image.constant(d_args['tmax']), tmin=ee.Image.constant(d_args['tmin']),
        ea=ee.Image.constant(d_args['ea']), rs=ee.Image.constant(d_args['rs']),
        uz=ee.Image.constant(d_args['uz']), zw=ee.Number(s_args['zw']),
        elev=ee.Number(s_args['elev']), lat=ee.Number(s_args['lat']),
        doy=ee.Number(d_args['doy']), method='refet',
    )

    output = utils.constant_image_value(refet.etsz(surface).rename(['etsz']))

    assert float(output['etsz']) == pytest.approx(expected)


def test_refet_daily_gridmet_etr():
    """Generate a mock GRIDMET image from the test values"""
    input_img = (
        ee.Image.constant([
            d_args['tmax'] + 273.15, d_args['tmin'] + 273.15,
            d_args['q_asce'], d_args['rs'] / 0.0864, d_args['uz']
        ])
        .rename(['tmmx', 'tmmn', 'sph', 'srad', 'vs'])
        .set('system:time_start', ee.Date('2015-07-01').millis())
    )

    refet = Daily.gridmet(
        ee.Image(input_img), elev=ee.Number(s_args['elev']),
        lat=ee.Number(s_args['lat']), zw=ee.Number(s_args['zw']), method='asce',
    )

    output = utils.constant_image_value(refet.etr)

    assert float(output['etr']) == pytest.approx(d_args['etr_asce'])


def test_refet_daily_maca_etr():
    """Generate a mock MACA image from the test values"""
    input_img = (
        ee.Image.constant([
            d_args['tmax'] + 273.15, d_args['tmin'] + 273.15,
            d_args['q_asce'], d_args['rs'] / 0.0864,
            d_args['uz_x'], d_args['uz_y']])
        .rename(['tasmax', 'tasmin', 'huss', 'rsds', 'uas', 'vas'])
        .set('system:time_start', ee.Date('2015-07-01').millis())
    )

    refet = Daily.maca(
        ee.Image(input_img), elev=ee.Number(s_args['elev']),
        lat=ee.Number(s_args['lat']), zw=ee.Number(s_args['zw']), method='asce',
    )

    output = utils.constant_image_value(refet.etr)

    assert float(output['etr']) == pytest.approx(d_args['etr_asce'])


def test_refet_daily_nldas_etr():
    """Generate a mock NLDAS image from the test values

    Convert the test Rs from MJ m-2 d-1 to W m-2, then allocate half to each image

    """
    band_names = [
        'temperature', 'specific_humidity', 'shortwave_radiation', 'wind_u', 'wind_v'
    ]

    input_coll = ee.ImageCollection.fromImages([
        ee.Image.constant([d_args['tmin'], d_args['q_asce'], 0.0, d_args['uz_x'], d_args['uz_y']])
            .double().rename(band_names)
            .set({'system:time_start': ee.Date('2015-07-01T00:00:00', 'UTC').millis()}),
        ee.Image.constant([d_args['tmax'], d_args['q_asce'], d_args['rs'] / 0.0036,
                           d_args['uz_x'], d_args['uz_y']])
            .double().rename(band_names)
            .set({'system:time_start': ee.Date('2015-07-01T12:00:00', 'UTC').millis()})
    ])

    refet = Daily.nldas(
        input_coll, elev=ee.Number(s_args['elev']),
        lat=ee.Number(s_args['lat']), zw=ee.Number(s_args['zw']), method='asce',
    )

    output = utils.constant_image_value(refet.etr)

    assert float(output['etr']) == pytest.approx(d_args['etr_asce'])


def test_refet_daily_nldas_eto():
    """Compare values to previous calculations"""
    test_point = ee.Geometry.Point(-120.113, 36.336)

    input_coll = (
        ee.ImageCollection('NASA/NLDAS/FORA0125_H002')
        .filterDate('2017-07-01T06:00:00', '2017-07-02T06:00:00')
    )
    output = utils.get_info(
        Daily.nldas(input_coll, method='asce').eto
        .reduceRegion(ee.Reducer.first(), geometry=test_point, scale=1)
    )
    expected = utils.get_info(
        ee.Image('projects/eddi-noaa/nldas/daily/20170701')
        .select(['eto_asce'])
        .reduceRegion(ee.Reducer.first(), geometry=test_point, scale=1)
    )

    assert output['eto'] == pytest.approx(expected['eto_asce'], rel=0.001)


def test_refet_daily_cfsv2_etr():
    """Generate a mock CFSv2 image from the test values

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

    input_coll = ee.ImageCollection.fromImages([
        ee.Image.constant([d_args['tmax'] + 273.15, d_args['tmin'] + 273.15,
                           d_args['q_asce'], d_args['rs'] / 0.0864,
                           d_args['uz_x'], d_args['uz_y']])
            .double().rename(band_names)
            .set({'system:time_start': ee.Date('2015-07-01T00:00:00', 'UTC').millis()}),
        ee.Image.constant([d_args['tmax'] + 273.15, d_args['tmin'] + 273.15,
                           d_args['q_asce'], d_args['rs'] / 0.0864,
                           d_args['uz_x'], d_args['uz_y']])
            .double().rename(band_names)
            .set({'system:time_start': ee.Date('2015-07-01T06:00:00', 'UTC').millis()}),
        ee.Image.constant([d_args['tmax'] + 273.15, d_args['tmin'] + 273.15,
                           d_args['q_asce'], d_args['rs'] / 0.0864,
                           d_args['uz_x'], d_args['uz_y']])
            .double().rename(band_names)
            .set({'system:time_start': ee.Date('2015-07-01T12:00:00', 'UTC').millis()}),
        ee.Image.constant([d_args['tmax'] + 273.15, d_args['tmin'] + 273.15,
                           d_args['q_asce'], d_args['rs'] / 0.0864,
                           d_args['uz_x'], d_args['uz_y']])
            .double().rename(band_names)
            .set({'system:time_start': ee.Date('2015-07-01T18:00:00', 'UTC').millis()}),
    ])

    refet = Daily.cfsv2(
        input_coll, elev=ee.Number(s_args['elev']), lat=ee.Number(s_args['lat']),
        zw=ee.Number(s_args['zw']), method='asce',
    )

    output = utils.constant_image_value(refet.etr)

    assert float(output['etr']) == pytest.approx(d_args['etr_asce'])


def test_refet_daily_rtma_etr():
    """Generate a mock RTMA image from the test values"""
    band_names = ['TMP', 'SPFH', 'WIND']

    input_coll = ee.ImageCollection.fromImages([
        ee.Image.constant([d_args['tmin'], d_args['q_asce'], d_args['uz']])
            .double().rename(band_names)
            .set({'system:time_start': ee.Date('2015-07-01T00:00:00', 'UTC').millis()}),
        ee.Image.constant([d_args['tmax'], d_args['q_asce'], d_args['uz']])
            .double().rename(band_names)
            .set({'system:time_start': ee.Date('2015-07-01T12:00:00', 'UTC').millis()})
    ])

    refet = Daily.rtma(
        input_coll, rs=ee.Image.constant(d_args['rs']),
        elev=ee.Number(s_args['elev']), lat=ee.Number(s_args['lat']),
        zw=ee.Number(s_args['zw']), method='asce',
    )

    output = utils.constant_image_value(refet.etr)

    assert float(output['etr']) == pytest.approx(d_args['etr_asce'])


def test_refet_daily_era5_etr():
    """Generate a mock ERA5 image from the test values"""
    band_names = [
        'temperature_2m', 'dewpoint_temperature_2m',
        'surface_solar_radiation_downwards',
        'u_component_of_wind_10m', 'v_component_of_wind_10m',
    ]

    # Each image needs half of the total Rs since hourly images are summed
    input_coll = ee.ImageCollection.fromImages([
        ee.Image.constant([d_args['tmin'] + 273.15, d_args['tdew'] + 273.15,
                           0.0, d_args['uz_x'], d_args['uz_y']])
            .double().rename(band_names)
            .set({'system:time_start': ee.Date('2015-07-01T00:00:00', 'UTC').millis()}),
        ee.Image.constant([d_args['tmax'] + 273.15, d_args['tdew'] + 273.15,
                           d_args['rs'] * 1000000, d_args['uz_x'], d_args['uz_y']])
            .double().rename(band_names)
            .set({'system:time_start': ee.Date('2015-07-01T12:00:00', 'UTC').millis()})
    ])

    refet = Daily.era5(
        input_coll, elev=ee.Number(s_args['elev']), lat=ee.Number(s_args['lat']),
        zw=ee.Number(s_args['zw']), method='asce',
    )

    output = utils.constant_image_value(refet.etr)

    assert float(output['etr']) == pytest.approx(d_args['etr_asce'])


def test_refet_daily_era5_land_etr():
    """Generate a mock ERA5-Land image from the test values"""
    band_names = [
        'temperature_2m', 'dewpoint_temperature_2m',
        'surface_solar_radiation_downwards_hourly',
        'u_component_of_wind_10m', 'v_component_of_wind_10m',
    ]

    # Each image needs half of the total Rs since hourly images are summed
    input_coll = ee.ImageCollection.fromImages([
        ee.Image.constant([d_args['tmin'] + 273.15, d_args['tdew'] + 273.15,
                           0.0, d_args['uz_x'], d_args['uz_y']])
            .double().rename(band_names)
            .set({'system:time_start': ee.Date('2015-07-01T00:00:00', 'UTC').millis()}),
        ee.Image.constant([d_args['tmax'] + 273.15, d_args['tdew'] + 273.15,
                           d_args['rs'] * 1000000, d_args['uz_x'], d_args['uz_y']])
            .double().rename(band_names)
            .set({'system:time_start': ee.Date('2015-07-01T12:00:00', 'UTC').millis()})
    ])

    refet = Daily.era5_land(
        input_coll, elev=ee.Number(s_args['elev']), lat=ee.Number(s_args['lat']),
        zw=ee.Number(s_args['zw']), method='asce',
    )

    output = utils.constant_image_value(refet.etr)

    assert float(output['etr']) == pytest.approx(d_args['etr_asce'])


def test_refet_daily_era5_land_fill_edge_cells():
    """Check that the fill_edge_cells flag works for an edge cell along the coast of England"""
    input_coll = ee.ImageCollection('ECMWF/ERA5_LAND/HOURLY').filterDate('2015-07-01', '2015-07-02')
    output = utils.point_image_value(Daily.era5_land(input_coll, fill_edge_cells=False).etr, xy=[0.0, 50.7])
    assert output['etr'] is None
    output = utils.point_image_value(Daily.era5_land(input_coll, fill_edge_cells=True).etr, xy=[0.0, 50.7])
    assert output['etr'] is not None


# TODO: Add a test for using the default Rs when one is not provided
