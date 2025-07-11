import ee
import pytest

import openet.refetgee.calcs as calcs
import openet.refetgee.units as units
import utils

# Eventually move to conftest.py or a separate file
# Fallon AgriMet site parameters
s_args = {
    'elev': 1208.5,
    'lat': units.deg2rad(39.4575),
    'lon': units.deg2rad(-118.77388),
    'pair': 87.81876435813037,
    'pair_asce': 87.80710537212929,
    'zw': 3.0,
}
# Daily test parameters for 2015-07-01
d_args = {
    'doy': 182,
    'delta': 0.40352881013673136,
    'delta_asce': 0.4029517192078854,
    'doy_frac': 3.132985550429273,
    'dr': 0.9670012223491632,
    'ea': 1.2206674169951346,
    'es': 4.6747236227258835,
    'es_slope': 0.23489129849801055,
    'es_slope_asce': 0.23488581814172638,
    # 'eto': 7.942481120179387,
    # 'etr': 10.571560006380153,
    'fcd': 0.8569860867772078,
    'omega_s': 1.9298904620748385,
    'q': 0.008691370735727117,          # Computed from Ea from Tdew
    'q_asce': 0.008692530868140688,     # Computed from Ea from Tdew
    'ra': 41.67610845067083,
    'ra_asce': 41.64824567735701,
    'rn': 15.174377350374279,
    'rnl': 6.556533974825727,
    'rs': 674.07 * 0.041868,            # Conversion from Langleys to MJ m-2
    'rso': 31.565939444861765,
    'rso_simple': 32.26439287925584,
    'sc': -0.05874166519510547,
    'tdew': units.f2c(49.84),
    'tmin': units.f2c(66.65),
    'tmax': units.f2c(102.80),
    # 'tmean': f2c(84.725),
    'w': 17.107650595384076,
    'uz': 4.80 * 0.44704,               # Conversion from mph to m s-1
    'u2': 1.976111757722194,
}
# Hourly test parameters for 2015-07-01 18:00 UTC (11:00 AM PDT)
h_args = {
    'doy': 182,
    'ea': 1.1990099614301906,           # Computed from Tdew
    'es': 5.09318785259078,
    # 'eto': 0.6065255163817055,
    # 'etr': 0.7201865213918281,
    'fcd': 0.6816142001345745,
    'fcd_asce': 0.6816142001345745,
    # 'omega': -0.3866777826605525,     # Computed from time_mid
    'omega': -0.5175774765601271,       # Computed from time to match IN2
    'q': 0.008536365803069757,          # Computed from Ea from Tdew
    'q_asce': 0.00853750513849305,      # Computed from Ea from Tdew
    'ra': 4.30824147948541,
    'ra_asce': 4.30635461533285,
    'rnl': 0.22897874401150786,
    'rs': 61.16 * 0.041868,             # Conversion from Langleys to MJ m-2
    'rso': 3.350936122776373,
    'rso_asce': 3.350851392549374,
    'rso_simple': 3.33531130617322,
    # 'solar_time': -1.4770003318617722,  # Computed from time_mid
    'solar_time': -1.9770003318617722,  # Computed from time to match IN2
    'tdew': units.f2c(49.36),
    'time': 18.0,
    'time_mid': 18.5,
    'tmean': units.f2c(91.80),
    'uz': 3.33 * 0.44704,  # Conversion from mph to m s-1
    'u2': 1.3709275319197722,
}

# # Playing around with passing a dictionary of function keyword arguments
# daily_args = {
#     k: v for k, v in daily.items()
#     if k in ['tmin', 'tmax', 'ea', 'rs', 'uz', 'zw', 'doy']}
# daily_args.update({'ref_type':'etr', 'rso_type': 'full', 'rso': None})
# hourly_args = {
#     k: v for k, v in daily.items()
#     if k in ['tmean', 'ea', 'rs', 'uz', 'zw', 'doy', 'time']}
# hourly_args.update({'ref_type':'etr'})


# Test ancillary functions with positional inputs
def test_air_pressure_default(elev=s_args['elev'], expected=s_args['pair_asce']):
    output = utils.get_info(calcs.air_pressure(elev=ee.Number(elev)))
    assert float(output) == pytest.approx(expected)

def test_air_pressure_asce(elev=s_args['elev'], expected=s_args['pair_asce']):
    output = utils.get_info(calcs.air_pressure(elev=ee.Number(elev), method='asce'))
    assert float(output) == pytest.approx(expected)

def test_air_pressure_refet(elev=s_args['elev'], expected=s_args['pair']):
    output = utils.get_info(calcs.air_pressure(elev=ee.Number(elev), method='refet'))
    assert float(output) == pytest.approx(expected)

def test_air_pressure_position(elev=s_args['elev'], expected=s_args['pair']):
    output = utils.get_info(calcs.air_pressure(ee.Number(elev), 'refet'))
    assert float(output) == pytest.approx(expected)

def test_air_pressure_image(elev=s_args['elev'], expected=s_args['pair_asce']):
    output = utils.constant_image_value(calcs.air_pressure(elev=ee.Image.constant(elev)))
    assert float(output['constant']) == pytest.approx(expected)


@pytest.mark.parametrize(
    'tdew, expected',
    [
        [d_args['tdew'], d_args['ea']],
        [h_args['tmean'], h_args['es']],
        [h_args['tdew'], h_args['ea']],
    ]
)
def test_sat_vapor_pressure_number(tdew, expected):
    output = utils.get_info(calcs.sat_vapor_pressure(ee.Number(tdew)))
    assert float(output) == pytest.approx(expected)

def test_sat_vapor_pressure_keyword(tdew=d_args['tdew'], expected=d_args['ea']):
    output = utils.get_info(calcs.sat_vapor_pressure(temperature=ee.Number(tdew)))
    assert float(output) == pytest.approx(expected)

def test_sat_vapor_pressure_image(tdew=d_args['tdew'], expected=d_args['ea']):
    output = utils.constant_image_value(calcs.sat_vapor_pressure(ee.Image.constant(tdew)))
    assert float(output['constant']) == pytest.approx(expected)


@pytest.mark.parametrize(
    'ea, pair, expected',
    [
        [d_args['ea'], s_args['pair'], d_args['q']],
        [d_args['ea'], s_args['pair_asce'], d_args['q_asce']],
        [h_args['ea'], s_args['pair'], h_args['q']],
        [h_args['ea'], s_args['pair_asce'], h_args['q_asce']],
    ]
)
def test_specific_humidity_number(ea, pair, expected):
    output = utils.get_info(calcs.specific_humidity(
        ea=ee.Number(ea), pair=ee.Number(pair)
    ))
    assert float(output) == pytest.approx(expected)

def test_specific_humidity_position(ea=d_args['ea'], pair=s_args['pair'], expected=d_args['q']):
    output = utils.get_info(calcs.specific_humidity(ee.Number(ea), ee.Number(pair)))
    assert float(output) == pytest.approx(expected)

def test_specific_humidity_image(ea=d_args['ea'], pair=s_args['pair'], expected=d_args['q']):
    output = utils.constant_image_value(calcs.specific_humidity(
        ea=ee.Image.constant(ea), pair=ee.Number(pair)
    ))
    assert float(output['constant']) == pytest.approx(expected)


def test_actual_vapor_pressure_number(q=d_args['q'], pair=s_args['pair'], expected=d_args['ea']):
    output = utils.get_info(calcs.actual_vapor_pressure(q=ee.Number(q), pair=ee.Number(pair)))
    assert float(output) == pytest.approx(expected)

def test_actual_vapor_pressure_asce(
        q=d_args['q_asce'], pair=s_args['pair_asce'], expected=d_args['ea']
):
    output = utils.get_info(calcs.actual_vapor_pressure(q=ee.Number(q), pair=ee.Number(pair)))
    assert float(output) == pytest.approx(expected)

def test_actual_vapor_pressure_position(q=d_args['q'], pair=s_args['pair'], expected=d_args['ea']):
    output = utils.get_info(calcs.actual_vapor_pressure(ee.Number(q), ee.Number(pair)))
    assert float(output) == pytest.approx(expected)

def test_actual_vapor_pressure_image(q=d_args['q'], pair=s_args['pair'], ea=d_args['ea']):
    output = utils.constant_image_value(calcs.actual_vapor_pressure(
        q=ee.Image.constant(q), pair=ee.Number(pair)
    ))
    assert float(output['constant']) == pytest.approx(ea)


def test_vpd_number(es=d_args['es'], ea=d_args['ea']):
    output = utils.get_info(calcs.vpd(es=ee.Number(es), ea=ee.Number(ea)))
    assert float(output) == pytest.approx(float(es - ea))

    # Check that negative VPD's are set to 0
    output = utils.get_info(calcs.vpd(es=ee.Number(es), ea=ee.Number(es + 1)))
    assert float(output) == pytest.approx(0)

def test_vpd_image(es=d_args['es'], ea=d_args['ea']):
    output = utils.constant_image_value(calcs.vpd(es=ee.Image.constant(es), ea=ee.Number(ea)))
    assert float(output['constant']) == pytest.approx(float(es-ea))


def test_es_slope_default(tmin=d_args['tmin'], tmax=d_args['tmax'], expected=d_args['es_slope_asce']):
    output = utils.get_info(calcs.es_slope(tmean=ee.Number(0.5 * (tmin + tmax))))
    assert float(output) == pytest.approx(expected)

def test_es_slope_asce(tmin=d_args['tmin'], tmax=d_args['tmax'], expected=d_args['es_slope_asce']):
    output = utils.get_info(calcs.es_slope(
        tmean=ee.Number(0.5 * (tmin + tmax)), method='asce'
    ))
    assert float(output) == pytest.approx(expected)

def test_es_slope_refet(tmin=d_args['tmin'], tmax=d_args['tmax'], expected=d_args['es_slope']):
    output = utils.get_info(calcs.es_slope(tmean=ee.Number(0.5 * (tmin + tmax)), method='refet'))
    assert float(output) == pytest.approx(expected)

def test_es_slope_position(tmin=d_args['tmin'], tmax=d_args['tmax'], expected=d_args['es_slope']):
    output = utils.get_info(calcs.es_slope(ee.Number(0.5 * (tmin + tmax)), 'refet'))
    assert float(output) == pytest.approx(expected)

def test_es_slope_image(tmin=d_args['tmin'], tmax=d_args['tmax'], expected=d_args['es_slope_asce']):
    output = utils.constant_image_value(calcs.es_slope(
        tmean=ee.Image.constant(0.5 * (tmin + tmax))
    ))
    assert float(output['constant']) == pytest.approx(expected)


def test_precipitable_water_number(ea=d_args['ea'], pair=s_args['pair'], expected=d_args['w']):
    output = utils.get_info(calcs.precipitable_water(ee.Number(ea), ee.Number(pair)))
    assert float(output) == pytest.approx(expected)

def test_precipitable_water_image_1(ea=d_args['ea'], pair=s_args['pair'], expected=d_args['w']):
    output = utils.constant_image_value(calcs.precipitable_water(
        ee.Image.constant(ea), ee.Number(pair)
    ))
    assert float(output['constant']) == pytest.approx(expected)

def test_precipitable_water_image_2(ea=d_args['ea'], pair=s_args['pair'], expected=d_args['w']):
    output = utils.constant_image_value(calcs.precipitable_water(
        ee.Image.constant(ea), ee.Image.constant(pair)
    ))
    assert float(output['constant']) == pytest.approx(expected)


def test_doy_fraction_number(doy=d_args['doy'], expected=d_args['doy_frac']):
    assert float(utils.get_info(calcs.doy_fraction(ee.Number(doy)))) == pytest.approx(expected)

def test_doy_fraction_keyword(doy=d_args['doy'], expected=d_args['doy_frac']):
    output = utils.get_info(calcs.doy_fraction(doy=ee.Number(doy)))
    assert float(output) == pytest.approx(expected)

def test_doy_fraction_image(doy=d_args['doy'], expected=d_args['doy_frac']):
    output = utils.constant_image_value(calcs.doy_fraction(ee.Image.constant(doy)))
    assert float(output['constant']) == pytest.approx(expected)


def test_declination_default(doy=d_args['doy'], expected=d_args['delta_asce']):
    assert float(utils.get_info(calcs.declination(ee.Number(doy)))) == pytest.approx(expected)

def test_declination_asce(doy=d_args['doy'], expected=d_args['delta_asce']):
    output = utils.get_info(calcs.declination(ee.Number(doy), method='asce'))
    assert float(output) == pytest.approx(expected)

def test_declination_refet(doy=d_args['doy'], expected=d_args['delta']):
    output = utils.get_info(calcs.declination(ee.Number(doy), method='refet'))
    assert float(output) == pytest.approx(expected)

def test_declination_keyword(doy=d_args['doy'], expected=d_args['delta_asce']):
    assert float(utils.get_info(calcs.declination(doy=ee.Number(doy)))) == pytest.approx(expected)

def test_declination_image(doy=d_args['doy'], expected=d_args['delta_asce']):
    output = utils.constant_image_value(calcs.declination(ee.Image.constant(doy)))
    assert float(output['constant']) == pytest.approx(expected)


def test_dr_number(doy=d_args['doy'], expected=d_args['dr']):
    assert float(utils.get_info(calcs.dr(ee.Number(doy)))) == pytest.approx(expected)

def test_dr_keyword(doy=d_args['doy'], expected=d_args['dr']):
    assert float(utils.get_info(calcs.dr(doy=ee.Number(doy)))) == pytest.approx(expected)

def test_dr_image(doy=d_args['doy'], expected=d_args['dr']):
    output = utils.constant_image_value(calcs.dr(ee.Image.constant(doy)))
    assert float(output['constant']) == pytest.approx(expected)


def test_seasonal_correction_number(doy=d_args['doy'], expected=d_args['sc']):
    output = utils.get_info(calcs.seasonal_correction(ee.Number(doy)))
    assert float(output) == pytest.approx(expected)

def test_seasonal_correction_keyword(doy=d_args['doy'], expected=d_args['sc']):
    output = utils.get_info(calcs.seasonal_correction(doy=ee.Number(doy)))
    assert float(output) == pytest.approx(expected)

def test_seasonal_correction_image(doy=d_args['doy'], expected=d_args['sc']):
    output = utils.constant_image_value(calcs.seasonal_correction(ee.Image.constant(doy)))
    assert float(output['constant']) == pytest.approx(expected)


def test_solar_time_rad_number(
        lon=s_args['lon'], time_mid=h_args['time'], sc=d_args['sc'], expected=h_args['solar_time']
):
    output = utils.get_info(calcs.solar_time_rad(
        ee.Number(lon), ee.Number(time_mid), ee.Number(sc)
    ))
    assert float(output) == pytest.approx(expected)


def test_solar_time_rad_image(
        lon=s_args['lon'], time_mid=h_args['time'], sc=d_args['sc'], expected=h_args['solar_time']
):
    output = utils.constant_image_value(calcs.solar_time_rad(
        ee.Image.constant(lon), ee.Number(time_mid), ee.Number(sc)
    ))
    assert float(output['constant']) == pytest.approx(expected)


# # The code will fail if the lon is a number but the time_mid is an image
# def test_solar_time_rad_image_fail(lon=s_args['lon'], time_mid=h_args['time'],
#                                    sc=d_args['sc'], expected=h_args['solar_time']):
#     output = utils.constant_image_value(calcs.solar_time_rad(
#         ee.Image.number(lon), ee.Image.constant(time_mid), ee.Number(sc)
#     ))
#     assert float(output['constant']) == pytest.approx(expected)


def test_sha_number(solar_time=h_args['solar_time'], expected=h_args['omega']):
    output = utils.get_info(calcs.sha(ee.Number(solar_time)))
    assert float(output) == pytest.approx(expected)

def test_sha_keyword(solar_time=h_args['solar_time'], expected=h_args['omega']):
    output = utils.get_info(calcs.sha(solar_time=ee.Number(solar_time)))
    assert float(output) == pytest.approx(expected)

def test_sha_image(solar_time=h_args['solar_time'], expected=h_args['omega']):
    output = utils.constant_image_value(calcs.sha(ee.Image.constant(solar_time)))
    assert float(output['constant']) == pytest.approx(expected)


@pytest.mark.parametrize(
    'x, x_min, x_max, expected',
    [
        [1.2, 1.2, 1.5, 1.2],
        [1.1, 1.2, 1.5, 1.4],
        [1.6, 1.2, 1.5, 1.3],
        [2.0, 1.2, 1.5, 1.4],
    ]
)
def test_wrap(x, x_min, x_max, expected):
    assert float(utils.get_info(calcs.wrap(ee.Number(x), x_min, x_max))) == pytest.approx(expected)


def test_sha_sunset_number(lat=s_args['lat'], delta=d_args['delta'], expected=d_args['omega_s']):
    output = utils.get_info(calcs.sha_sunset(lat=ee.Number(lat), delta=ee.Number(delta)))
    assert float(output) == pytest.approx(expected)

def test_sha_sunset_position(lat=s_args['lat'], delta=d_args['delta'], expected=d_args['omega_s']):
    output = utils.get_info(calcs.sha_sunset(ee.Number(lat), ee.Number(delta)))
    assert float(output) == pytest.approx(expected)

def test_sha_sunset_image(lat=s_args['lat'], delta=d_args['delta'], expected=d_args['omega_s']):
    output = utils.constant_image_value(calcs.sha_sunset(
        lat=ee.Image.constant(lat), delta=ee.Number(delta)
    ))
    assert float(output['constant']) == pytest.approx(expected)


# Ra Daily
def test_ra_daily_default(lat=s_args['lat'], doy=d_args['doy'], expected=d_args['ra_asce']):
    output = utils.get_info(calcs.ra_daily(lat=ee.Number(lat), doy=ee.Number(doy)))
    assert float(output) == pytest.approx(expected)

def test_ra_daily_asce(lat=s_args['lat'], doy=d_args['doy'], expected=d_args['ra_asce']):
    output = utils.get_info(calcs.ra_daily(
        lat=ee.Number(lat), doy=ee.Number(doy), method='asce'
    ))
    assert float(output) == pytest.approx(expected)

def test_ra_daily_refet(lat=s_args['lat'], doy=d_args['doy'], expected=d_args['ra']):
    output = utils.get_info(calcs.ra_daily(
        lat=ee.Number(lat), doy=ee.Number(doy), method='refet'
    ))
    assert float(output) == pytest.approx(expected)

def test_ra_daily_position(lat=s_args['lat'], doy=d_args['doy'], expected=d_args['ra_asce']):
    output = utils.get_info(calcs.ra_daily(ee.Number(lat), ee.Number(doy)))
    assert float(output) == pytest.approx(expected)

def test_ra_daily_image(lat=s_args['lat'], doy=d_args['doy'], expected=d_args['ra_asce']):
    output = utils.constant_image_value(calcs.ra_daily(
        lat=ee.Image.constant(lat), doy=ee.Number(doy)
    ))
    assert float(output['constant']) == pytest.approx(expected)


# Ra Hourly
def test_ra_hourly_default(lat=s_args['lat'], lon=s_args['lon'],
                           doy=h_args['doy'], time_mid=h_args['time_mid'],
                           expected=h_args['ra_asce']):
    output = utils.get_info(calcs.ra_hourly(
        lat=ee.Number(lat), lon=ee.Number(lon), doy=ee.Number(doy), time_mid=ee.Number(time_mid)
    ))
    assert float(output) == pytest.approx(expected)

def test_ra_hourly_asce(lat=s_args['lat'], lon=s_args['lon'],
                        doy=h_args['doy'], time_mid=h_args['time_mid'],
                        expected=h_args['ra_asce']):
    output = utils.get_info(calcs.ra_hourly(
        lat=ee.Number(lat), lon=ee.Number(lon), doy=ee.Number(doy),
        time_mid=ee.Number(time_mid), method='asce',
    ))
    assert float(output) == pytest.approx(expected)

def test_ra_hourly_refet(lat=s_args['lat'], lon=s_args['lon'],
                         doy=h_args['doy'], time_mid=h_args['time_mid'],
                         expected=h_args['ra']):
    output = utils.get_info(calcs.ra_hourly(
        lat=ee.Number(lat), lon=ee.Number(lon), doy=ee.Number(doy),
        time_mid=ee.Number(time_mid), method='refet',
    ))
    assert float(output) == pytest.approx(expected)

def test_ra_hourly_position(lat=s_args['lat'], lon=s_args['lon'],
                            doy=h_args['doy'], time_mid=h_args['time_mid'],
                            expected=h_args['ra_asce']):
    output = utils.get_info(calcs.ra_hourly(
        ee.Number(lat), ee.Number(lon), ee.Number(doy), ee.Number(time_mid)
    ))
    assert float(output) == pytest.approx(expected)

# def test_ra_hourly_image(lat=s_args['lat'], lon=s_args['lon'],
#                          doy=h_args['doy'], time_mid=h_args['time_mid'],
#                          expected=h_args['ra_asce']):
#     output = utils.constant_image_value(calcs.ra_hourly(
#         lat=ee.Number(lat), lon=ee.Number(lon), doy=ee.Number(doy),
#         time_mid=ee.Number(time_mid)
#     ))
#     assert float(output['constant']) == pytest.approx(expected)


# Rso Full Daily
def test_rso_daily_number(
        ea=d_args['ea'], pair=s_args['pair'], ra=d_args['ra'], doy=d_args['doy'],
        lat=s_args['lat'], expected=d_args['rso']
):
    output = utils.get_info(calcs.rso_daily(
        ea=ee.Number(ea), pair=ee.Number(pair), ra=ee.Number(ra), doy=ee.Number(doy),
        lat=ee.Number(lat),
    ))
    assert float(output) == pytest.approx(expected)

def test_rso_daily_position(
        ea=d_args['ea'], pair=s_args['pair'], ra=d_args['ra'], doy=d_args['doy'],
        lat=s_args['lat'], expected=d_args['rso']
):
    output = utils.get_info(calcs.rso_daily(
        ee.Number(ea), ee.Number(pair), ee.Number(ra), ee.Number(doy), ee.Number(lat)
    ))
    assert float(output) == pytest.approx(expected)

def test_rso_daily_image(
        ea=d_args['ea'], pair=s_args['pair'], ra=d_args['ra'], doy=d_args['doy'],
        lat=s_args['lat'], expected=d_args['rso']
):
    output = utils.constant_image_value(
        calcs.rso_daily(
            ea=ee.Image.constant(ea), pair=ee.Number(pair), ra=ee.Number(ra),
            doy=ee.Number(doy), lat=ee.Number(lat),
        )
    )
    assert float(output['constant']) == pytest.approx(expected)

def test_rso_daily_image_doy_number(
        ea=d_args['ea'], ra=d_args['ra'], pair=s_args['pair'], doy=d_args['doy'],
        lat=s_args['lat'], expected=d_args['rso']
):
    # Check that the function works correctly when DOY is a number but other inputs are images
    # The current formulation requires that pair be an image if lat is
    output = utils.constant_image_value(calcs.rso_daily(
        ea=ee.Image.constant(ea), pair=ee.Image.constant(pair), ra=ee.Number(ra),
        doy=ee.Number(doy), lat=ee.Image.constant(lat),
    ))
    assert float(output['constant']) == pytest.approx(expected)


# Rso Full Hourly
def test_rso_hourly_default(ea=h_args['ea'], ra=h_args['ra'],
                            pair=s_args['pair'], doy=h_args['doy'],
                            time_mid=h_args['time_mid'], lat=s_args['lat'],
                            lon=s_args['lon'], expected=h_args['rso_asce']):
    output = utils.get_info(calcs.rso_hourly(
        ea=ee.Number(ea), ra=ee.Number(ra), pair=ee.Number(pair),
        doy=ee.Number(doy), time_mid=ee.Number(time_mid), lat=ee.Number(lat),
        lon=ee.Number(lon),
    ))
    assert float(output) == pytest.approx(expected)

def test_rso_hourly_asce(ea=h_args['ea'], ra=h_args['ra'], pair=s_args['pair'],
                         doy=h_args['doy'], time_mid=h_args['time_mid'],
                         lat=s_args['lat'], lon=s_args['lon'],
                         expected=h_args['rso_asce']):
    output = utils.get_info(calcs.rso_hourly(
        ea=ee.Number(ea), ra=ee.Number(ra), pair=ee.Number(pair),
        doy=ee.Number(doy), time_mid=ee.Number(time_mid), lat=ee.Number(lat),
        lon=ee.Number(lon), method='asce',
    ))
    assert float(output) == pytest.approx(expected)

def test_rso_hourly_refet(ea=h_args['ea'], ra=h_args['ra'],
                          pair=s_args['pair'], doy=h_args['doy'],
                          time_mid=h_args['time_mid'], lat=s_args['lat'],
                          lon=s_args['lon'], expected=h_args['rso']):
    output = utils.get_info(calcs.rso_hourly(
        ea=ee.Number(ea), ra=ee.Number(ra), pair=ee.Number(pair),
        doy=ee.Number(doy), time_mid=ee.Number(time_mid),
        lat=ee.Number(lat), lon=ee.Number(lon), method='refet',
    ))
    assert float(output) == pytest.approx(expected)

def test_rso_hourly_position(ea=h_args['ea'], ra=h_args['ra'],
                             pair=s_args['pair'], doy=h_args['doy'],
                             time_mid=h_args['time_mid'], lat=s_args['lat'],
                             lon=s_args['lon'], expected=h_args['rso_asce']):
    output = utils.get_info(calcs.rso_hourly(
        ee.Number(ea), ee.Number(ra), ee.Number(pair), ee.Number(doy),
        ee.Number(time_mid), ee.Number(lat), ee.Number(lon)
    ))
    assert float(output) == pytest.approx(expected)

def test_rso_hourly_image(ea=h_args['ea'], ra=h_args['ra'],
                          pair=s_args['pair'],
                          doy=h_args['doy'], time_mid=h_args['time_mid'],
                          lat=s_args['lat'], lon=s_args['lon'],
                          expected=h_args['rso_asce']):
    output = utils.constant_image_value(calcs.rso_hourly(
        ea=ee.Image.constant(ea), ra=ee.Image.constant(ra), pair=ee.Number(pair),
        doy=ee.Number(doy), time_mid=ee.Number(time_mid), lat=ee.Number(lat),
        lon=ee.Number(lon),
    ))
    assert float(output['constant']) == pytest.approx(expected)


# Rso Simple
@pytest.mark.parametrize(
    'ra, elev, expected',
    [
        [d_args['ra'], s_args['elev'], d_args['rso_simple']],
        [h_args['ra'], s_args['elev'], h_args['rso_simple']]
    ]
)
def test_rso_simple_number(ra, elev, expected):
    output = utils.get_info(calcs.rso_simple(ra=ee.Number(ra), elev=ee.Number(elev)))
    assert float(output) == pytest.approx(expected)

def test_rso_simple_position(ra=d_args['ra'], elev=s_args['elev'], expected=d_args['rso_simple']):
    output = utils.get_info(calcs.rso_simple(ee.Number(ra), ee.Number(elev)))
    assert float(output) == pytest.approx(expected)

def test_rso_simple_image(ra=d_args['ra'], elev=s_args['elev'], expected=d_args['rso_simple']):
    output = utils.constant_image_value(calcs.rso_simple(
        ra=ee.Image.constant(ra), elev=ee.Number(elev)
    ))
    assert float(output['constant']) == pytest.approx(expected)


# Fcd daily
def test_fcd_daily_number(rs=d_args['rs'], rso=d_args['rso'], expected=d_args['fcd']):
    output = utils.get_info(calcs.fcd_daily(rs=ee.Number(rs), rso=ee.Number(rso)))
    assert float(output) == pytest.approx(expected)

def test_fcd_daily_position(rs=d_args['rs'], rso=d_args['rso'], expected=d_args['fcd']):
    output = utils.get_info(calcs.fcd_daily(ee.Number(rs), ee.Number(rso)))
    assert float(output) == pytest.approx(expected)

def test_fcd_daily_image(rs=d_args['rs'], rso=d_args['rso'], expected=d_args['fcd']):
    output = utils.constant_image_value(calcs.fcd_daily(
        rs=ee.Image.constant(rs), rso=ee.Number(rso)
    ))
    assert float(output['constant']) == pytest.approx(expected)


# Fcd hourly
def test_fcd_hourly_default(rs=h_args['rs'], rso=h_args['rso'],
                            doy=h_args['doy'], time_mid=h_args['time_mid'],
                            lat=s_args['lat'], lon=s_args['lon'],
                            expected=h_args['fcd_asce']):
    output = utils.get_info(calcs.fcd_hourly(
        rs=ee.Number(rs), rso=ee.Number(rso), doy=ee.Number(doy),
        time_mid=ee.Number(time_mid), lat=ee.Number(lat), lon=ee.Number(lon),
    ))
    assert float(output) == pytest.approx(expected)

def test_fcd_hourly_asce(rs=h_args['rs'], rso=h_args['rso'],
                         doy=h_args['doy'], time_mid=h_args['time_mid'],
                         lat=s_args['lat'], lon=s_args['lon'],
                         expected=h_args['fcd_asce']):
    output = utils.get_info(calcs.fcd_hourly(
        rs=ee.Number(rs), rso=ee.Number(rso), doy=ee.Number(doy),
        time_mid=ee.Number(time_mid), lat=ee.Number(lat), lon=ee.Number(lon),
        method='asce',
    ))
    assert float(output) == pytest.approx(expected)

def test_fcd_hourly_refet(rs=h_args['rs'], rso=h_args['rso'],
                          doy=h_args['doy'], time_mid=h_args['time_mid'],
                          lat=s_args['lat'], lon=s_args['lon'],
                          expected=h_args['fcd']):
    output = utils.get_info(calcs.fcd_hourly(
        rs=ee.Number(rs), rso=ee.Number(rso), doy=ee.Number(doy),
        time_mid=ee.Number(time_mid), lat=ee.Number(lat), lon=ee.Number(lon),
        method='refet',
    ))
    assert float(output) == pytest.approx(expected)

def test_fcd_hourly_night(rs=h_args['rs'], rso=h_args['rso'], doy=h_args['doy'],
                          time_mid=6.5, lat=s_args['lat'], lon=s_args['lon'],
                          expected=1):
    # For now, check that nighttime fcd values are set to 1
    output = utils.get_info(calcs.fcd_hourly(
        rs=ee.Number(rs), rso=ee.Number(rso), doy=ee.Number(doy),
        time_mid=ee.Number(time_mid), lat=ee.Number(lat), lon=ee.Number(lon),
        method='refet',
    ))
    assert float(output) == pytest.approx(expected)

def test_fcd_hourly_position(rs=h_args['rs'], rso=h_args['rso'],
                             doy=h_args['doy'], time_mid=h_args['time_mid'],
                             lat=s_args['lat'], lon=s_args['lon'],
                             expected=h_args['fcd_asce']):
    output = utils.get_info(calcs.fcd_hourly(
        ee.Number(rs), ee.Number(rso), ee.Number(doy), ee.Number(time_mid),
        ee.Number(lat), ee.Number(lon),
    ))
    assert float(output) == pytest.approx(expected)

def test_fcd_hourly_image(rs=h_args['rs'], rso=h_args['rso'],
                          doy=h_args['doy'], time_mid=h_args['time_mid'],
                          lat=s_args['lat'], lon=s_args['lon'],
                          expected=h_args['fcd_asce']):
    output = utils.constant_image_value(calcs.fcd_hourly(
        rs=ee.Image.constant(rs), rso=ee.Number(rso), doy=ee.Number(doy),
        time_mid=ee.Number(time_mid), lat=ee.Number(lat), lon=ee.Number(lon),
    ))
    assert float(output['constant']) == pytest.approx(expected)


# Rnl Daily
def test_rnl_daily_number(tmin=d_args['tmin'], tmax=d_args['tmax'],
                          ea=d_args['ea'], fcd=d_args['fcd'], expected=d_args['rnl']):
    output = utils.get_info(calcs.rnl_daily(
        tmax=ee.Number(tmax), tmin=ee.Number(tmin), ea=ee.Number(ea),
        fcd=ee.Number(fcd),
    ))
    assert float(output) == pytest.approx(expected)

def test_rnl_daily_position(tmin=d_args['tmin'], tmax=d_args['tmax'],
                            ea=d_args['ea'], fcd=d_args['fcd'],
                            expected=d_args['rnl']):
    output = utils.get_info(calcs.rnl_daily(
        ee.Number(tmax), ee.Number(tmin), ee.Number(ea), ee.Number(fcd)
    ))
    assert float(output) == pytest.approx(expected)

def test_rnl_daily_image(tmin=d_args['tmin'], tmax=d_args['tmax'],
                         ea=d_args['ea'], fcd=d_args['fcd'], expected=d_args['rnl']):
    output = utils.constant_image_value(calcs.rnl_daily(
        tmax=ee.Image.constant(tmax), tmin=ee.Number(tmin),
        ea=ee.Number(ea), fcd=ee.Number(fcd),
    ))
    assert float(output['constant']) == pytest.approx(expected)


# Rnl Hourly
def test_rnl_hourly_number(tmean=h_args['tmean'], ea=h_args['ea'],
                           fcd=h_args['fcd'], expected=h_args['rnl']):
    output = utils.get_info(calcs.rnl_hourly(
        tmean=ee.Number(tmean), ea=ee.Number(ea), fcd=ee.Number(fcd)
    ))
    assert float(output) == pytest.approx(expected)

def test_rnl_hourly_position(tmean=h_args['tmean'], ea=h_args['ea'],
                             fcd=h_args['fcd'], expected=h_args['rnl']):
    output = utils.get_info(calcs.rnl_hourly(
        ee.Number(tmean), ee.Number(ea), ee.Number(fcd)
    ))
    assert float(output) == pytest.approx(expected)

def test_rnl_hourly_image(tmean=h_args['tmean'], ea=h_args['ea'],
                          fcd=h_args['fcd'], expected=h_args['rnl']):
    output = utils.constant_image_value(calcs.rnl_hourly(
        tmean=ee.Image.constant(tmean), ea=ee.Number(ea), fcd=ee.Number(fcd)
    ))
    assert float(output['constant']) == pytest.approx(expected)


# Rn
def test_rn_number(rs=d_args['rs'], rnl=d_args['rnl'], expected=d_args['rn']):
    output = utils.get_info(calcs.rn(rs=ee.Number(rs), rnl=ee.Number(rnl)))
    assert float(output) == pytest.approx(expected)

def test_rn_keywods(rs=d_args['rs'], rnl=d_args['rnl'], expected=d_args['rn']):
    output = utils.get_info(calcs.rn(ee.Number(rs), ee.Number(rnl)))
    assert float(output) == pytest.approx(expected)

def test_rn_image(rs=d_args['rs'], rnl=d_args['rnl'], expected=d_args['rn']):
    output = utils.constant_image_value(calcs.rn(
        rs=ee.Number(rs), rnl=ee.Image.constant(rnl)
    ))
    assert float(output['constant']) == pytest.approx(expected)


# Wind
def test_wind_height_adjust_2m(uz=2.5, zw=2.0, expected=2.5):
    output = utils.get_info(calcs.wind_height_adjust(uz=ee.Number(uz), zw=ee.Number(zw)))
    assert float(output) == pytest.approx(expected, abs=0.001)

def test_wind_height_adjust_position(uz=d_args['uz'], zw=s_args['zw'], expected=d_args['u2']):
    output = utils.get_info(calcs.wind_height_adjust(ee.Number(uz), ee.Number(zw)))
    assert float(output) == pytest.approx(expected)

@pytest.mark.parametrize(
    'uz, zw, expected',
    [
        [d_args['uz'], s_args['zw'], d_args['u2']],
        [h_args['uz'], s_args['zw'], h_args['u2']]
    ]
)
def test_wind_height_adjust_number(uz, zw, expected):
    output = utils.get_info(calcs.wind_height_adjust(uz=ee.Number(uz), zw=ee.Number(zw)))
    assert float(output) == pytest.approx(expected)

@pytest.mark.parametrize(
    'uz, zw, expected',
    [
        [d_args['uz'], s_args['zw'], d_args['u2']],
        [h_args['uz'], s_args['zw'], h_args['u2']]
    ]
)
def test_wind_height_adjust_image(uz, zw, expected):
    output = utils.constant_image_value(calcs.wind_height_adjust(
        uz=ee.Image.constant(uz), zw=ee.Number(zw)
    ))
    assert float(output['constant']) == pytest.approx(expected, abs=0.001)
