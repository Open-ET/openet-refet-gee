import ee
import pytest

from openet.refetgee import Monthly
import openet.refetgee.units as units


# Eventually move to conftest.py or a separate file
# Fallon AgriMet site parameters
s_args = {
    'elev': 1208.5,
    'lat': 39.4575,
    'lon': -118.77388,
}

# Monthy test parameters for 2015-07-01
m_args = {
    'month': 7,
    'pet_thornthwaite': 0,
    'tmean': units._f2c(84.725),
    'tmean_monthly': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
}

constant_geom = ee.Geometry.Rectangle([0, 0, 10, 10], 'EPSG:32613', False)


def test_refet_monthly_pet_hargreaves():
    time_start = ee.Date.fromYMD(2015, m_args['month'], 1).millis()
    refet = Monthly(
        tmean=ee.Image.constant(m_args['tmean']).rename(['tmean'])
            .set({'system:time_start': time_start}),
        lat=ee.Image.constant(ee.Number(s_args['lat'])),
        tmean_monthly=ee.Image.constant(m_args['tmean_monthly']),
    )
    output = refet.pet_thornthwaite\
        .reduceRegion(ee.Reducer.first(), geometry=constant_geom, scale=1)\
        .getInfo()
    assert float(output['pet_thornthwaite']) == pytest.approx(m_args['pet_thornthwaite'])


# def test_refet_monthly_prism_pet():
#     """Generate a fake PRISM monthly image from the test values"""
#     prism_img = ee.Image.constant([m_args['tmean'] + 273.15])\
#         .rename(['tmean'])\
#         .set('system:time_start', ee.Date('2015-07-01').millis())
#     refet = Monthly.prism(ee.Image(prism_img), lat=ee.Number(s_args['lat']))
#     output = refet.pet_thornthwaite\
#         .reduceRegion(ee.Reducer.first(), geometry=constant_geom, scale=1)\
#         .getInfo()
#     assert float(output['pet_thornthwaite']) == pytest.approx(m_args['pet_thornthwaite'])
