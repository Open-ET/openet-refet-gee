import math

import ee


def _air_pressure(elev, method='asce'):
    """Mean atmospheric pressure at station elevation (Eqs. 3 & 34)

    Parameters
    ----------
    elev : ee.Image or ee.Number
        Elevation [m].
    method : {'asce' (default), 'refet'}, optional
        Calculation method:
        * 'asce' -- Calculations will follow ASCE-EWRI 2005 [1] equations.
        * 'refet' -- Calculations will follow RefET software.

    Returns
    -------
    pair : ee.Image or ee.Number
        Air pressure [kPa].

    Notes
    -----
    The current calculation in Ref-ET:
        101.3 * (((293 - 0.0065 * elev) / 293) ** (9.8 / (0.0065 * 286.9)))
    Equation 3 in ASCE-EWRI 2005:
        101.3 * (((293 - 0.0065 * elev) / 293) ** 5.26)
    Per Dr. Allen, the calculation with full precision:
        101.3 * (((293.15 - 0.0065 * elev) / 293.15) ** (9.80665 / (0.0065 * 286.9)))

    """
    if method == 'asce':
        return elev.multiply(-0.0065).add(293).divide(293)\
            .pow(5.26).multiply(101.3)
    elif method == 'refet':
        return elev.multiply(-0.0065).add(293).divide(293)\
            .pow(9.8 / (0.0065 * 286.9)).multiply(101.3)


def _sat_vapor_pressure(temperature):
    """Saturation vapor pressure from temperature (Eq. 7)

    Parameters
    ----------
    temperature : ee.Image
        Air temperature [C].

    Returns
    -------
    e : ee.Image
        Saturation vapor pressure [kPa].

    Notes
    -----
    es = 0.6108 * exp(17.27 * temperature / (temperature + 237.3))

    """
    return temperature.add(237.3).pow(-1).multiply(temperature)\
        .multiply(17.27).exp().multiply(0.6108)
    # return temperature.multiply(17.27).divide(temperature.add(237.3)).exp()\
    #     .multiply(0.6108)


def _es_slope(tmean, method='asce'):
    """Slope of the saturation vapor pressure-temperature curve (Eq. 5)

    Parameters
    ----------
    tmean : ee.Image or ee.Number
        Mean air temperature [C].
    method : {'asce' (default), 'refet'}, optional
        Calculation method:
        * 'asce' -- Calculations will follow ASCE-EWRI 2005 [1] equations.
        * 'refet' -- Calculations will follow RefET software.

    Returns
    -------
    ee.Image or ee.Number

    Notes
    -----
    4098 * 0.6108 * exp(17.27 * T / (T + 237.3)) / ((T + 237.3) ** 2))

    """
    if method == 'refet':
        return _sat_vapor_pressure(tmean)\
            .multiply(4098.0).divide(tmean.add(237.3).pow(2))
    elif method == 'asce':
        return tmean.add(237.3).pow(-1).multiply(tmean).multiply(17.27).exp()\
            .multiply(2503.0).divide(tmean.add(237.3).pow(2))


def _actual_vapor_pressure(q, pair):
    """"Actual vapor pressure from specific humidity

    Parameters
    ----------
    q : ee.Image or ee.Number
        Specific humidity [kg/kg].
    pair : ee.Image or ee.Number
        Air pressure [kPa].

    Returns
    -------
    ea : ee.Image or ee.Number
        Actual vapor pressure [kPa].

    Notes
    -----
    ea = q * pair / (0.622 + 0.378 * q)

    """
    return q.multiply(0.378).add(0.622).pow(-1).multiply(q).multiply(pair)
    # return q.multiply(pair).divide(q.multiply(0.378).add(0.622))


def _specific_humidity(ea, pair):
    """"Specific humidity from actual vapor pressure

    Parameters
    ----------
    ea : ee.Image or ee.Number
        Specific humidity [kPa].
    pair : ee.Image or ee.Number
        Air pressure [kPa].

    Returns
    -------
    q : ee.Image or ee.Number
        Specific humidity [kg/kg].

    Notes
    -----
    q = 0.622 * ea / (pair - 0.378 * ea)

    """
    return ea.multiply(-0.378).add(pair).pow(-1).multiply(ea).multiply(0.622)
    # return ea.multiply(0.622).divide(ea.multiply(-0.378).add(pair))


def _vpd(es, ea):
    """Vapor pressure deficit

    Parameters
    ----------
    es : ee.Image or ee.Number
        Saturated vapor pressure [kPa].
    ea : ee.Image or ee.Number
        Actual vapor pressure [kPa].

    Returns
    -------
    ee.Image or ee.Number
        Vapor pressure deficit [kPa].

    """

    return es.subtract(ea).max(0)


def _precipitable_water(ea, pair):
    """Precipitable water in the atmosphere (Eq. D.3)

    Parameters
    ----------
    ea : ee.Image
        Vapor pressure [kPa].
    pair : ee.Image or ee.Number
        Air pressure [kPa].

    Returns
    -------
    ee.Image or ee.Number
        Precipitable water [mm].


    Notes
    -----
    w = pair * 0.14 * ea + 2.1

    """
    return ea.multiply(pair).multiply(0.14).add(2.1)


def _doy_fraction(doy):
    """Fraction of the DOY in the year (Eq. 50)

    Parameters
    ----------
    doy : ee.Image or ee.Number
        Day of year.

    Returns
    -------
    ee.Image or ee.Number
        DOY fraction [radians].

    """
    return doy.multiply(2.0 * math.pi / 365)


def _delta(doy, method='asce'):
    """Earth declination (Eq. 51)

    Parameters
    ----------
    doy : ee.Image or ee.Number
        Day of year.
    method : {'asce' (default), 'refet'}, optional
        Calculation method:
        * 'asce' -- Calculations will follow ASCE-EWRI 2005 [1] equations.
        * 'refet' -- Calculations will follow RefET software.

    Returns
    -------
    ee.Image or ee.Number
        Earth declination [radians].

    Notes
    -----
    Original equation in Duffie & Beckman (1980) (with conversions to radians):
        23.45 * (pi / 180) * sin(2 * pi * (doy + 284) / 365)
    Equation 24 in ASCE-EWRI (2005):
        0.409 * sin((2 * pi * doy / 365) - 1.39)

    """
    if method == 'asce':
        return _doy_fraction(doy).subtract(1.39).sin().multiply(0.409)
    else:
        return doy.add(284).multiply(2 * math.pi / 365).sin()\
            .multiply(23.45 * (math.pi / 180))


def _dr(doy):
    """Inverse square of the Earth-Sun Distance (Eq. 50)

    Parameters
    ----------
    doy : ee.Image or ee.Number
        Day of year.

    Returns
    -------
    ee.Image or ee.Number

    Notes
    -----
    This function returns 1 / d^2, not d, for direct use in radiance to
      TOA reflectance calculation
    pi * L * d^2 / (ESUN * cos(theta)) -> pi * L / (ESUN * cos(theta) * d)

    """
    return _doy_fraction(doy).cos().multiply(0.033).add(1.0)


def _seasonal_correction(doy):
    """Seasonal correction for solar time (Eqs. 57 & 58)

    Parameters
    ----------
    doy : ee.Image or ee.Number
        Day of year.

    Returns
    ------
    ee.Image or ee.Number
        Seasonal correction [hour]

    Notes
    -----
    sc = 0.1645 * sin(2 * b) - 0.1255 * cos(b) - 0.0250 * sin(b)

    """
    b = doy.subtract(81).divide(364.0).multiply(2 * math.pi)
    return b.multiply(2).sin().multiply(0.1645)\
        .subtract(b.cos().multiply(0.1255)).subtract(b.sin().multiply(0.0250))


def _solar_time_rad(lon, time_mid, sc):
    """Solar time (i.e. noon is 0) (Eq. 55)

    Parameters
    ----------
    lon : ee.Image or ee.Number
        Longitude [radians].
    time_mid : ee.Image or ee.Number
        UTC time at midpoint of period [hours].
    sc : ee.Image or ee.Number
        Seasonal correction [hours].

    Returns
    -------
    ee.Image or ee.Number
        Solar time [hours].

    Notes
    -----
    This function could be integrated into the _omega() function since they are
    always called together (i.e. _omega(_solar_time_rad()).  It was built
    independently from _omega to eventually support having a separate
    solar_time functions for longitude in degrees.
    time = time_mid + (lon * 24 / (2 * math.pi)) + sc - 12

    """
    return time_mid.add(lon.multiply(24 / (2 * math.pi))).add(sc).subtract(12)


def _omega(solar_time):
    """Solar hour angle (Eq. 55)

    Parameters
    ----------
    solar_time : ee.Image or ee.Number
        Solar time (i.e. noon is 0) [hours].

    Returns
    -------
    omega : ee.Image or ee.Number
        Hour angle [radians].

    """
    omega = solar_time.multiply(2 * math.pi / 24.0)

    # Need to adjust omega so that the values go from -pi to pi
    # Values outside this range are wrapped (i.e. -3*pi/2 -> pi/2)
    omega = _wrap(omega, -math.pi, math.pi)
    return omega


def _wrap(x, x_min, x_max):
    """Wrap floating point values into range

    Parameters
    ----------
    x : ee.Image or ee.Number
        Values to wrap.
    x_min : float
        Minimum value in output range.
    x_max : float
        Maximum value in output range.

    Returns
    -------
    ee.Image or ee.Number

    Notes
    -----
    This formula is used to mimic the Python modulo operator.
    Javascript/EE mod operator has the same sign as the dividend,
        so negative values stay negative.
    Python mod operator has the same sign as the divisor,
        so negative values wrap to positive.

    """
    x_range = x_max - x_min
    return x.subtract(x_min).mod(x_range).add(x_range).mod(x_range).add(x_min)
    # return np.mod((x - x_min), (x_max - x_min)) + x_min


def _omega_sunset(lat, delta):
    """Sunset hour angle (Eq. 59)

    Parameters
    ----------
    lat : ee.Image or ee.Number
        Latitude [radians].
    delta : ee.Image or ee.Number
        Earth declination [radians].

    Returns
    -------
    ee.Image or ee.Number
        Sunset hour angle [radians].

    Notes
    -----
    acos(-tan(lat) * tan(delta))

    """
    return lat.tan().multiply(-1).multiply(delta.tan()).acos()


def _ra_daily(lat, doy, method='asce'):
    """Daily extraterrestrial radiation (Eq. 21)

    Parameters
    ----------
    lat : ee.Image or ee.Number
        latitude [radians].
    doy : ee.Image or ee.Number
        Day of year.
    method : {'asce' (default), 'refet'}, optional
        Calculation method:
        * 'asce' -- Calculations will follow ASCE-EWRI 2005 [1] equations.
        * 'refet' -- Calculations will follow RefET software.

    Returns
    -------
    ra : ee.Image or ee.Number
        Daily extraterrestrial radiation [MJ m-2 d-1].

    Notes
    -----
    Equation in ASCE-EWRI 2005 uses a solar constant of ~1366.666... W m-2
    Equation in Duffie & Beckman (?) uses a solar constant of 1367 W m-2

    """
    delta = _delta(doy, method)
    omegas = _omega_sunset(lat, delta)
    theta = omegas.multiply(lat.sin()).multiply(delta.sin())\
        .add(lat.cos().multiply(delta.cos()).multiply(omegas.sin()))

    if method == 'asce':
        # (24. / math.pi) * 4.92 * _dr(doy) * theta
        ra = theta.multiply(_dr(doy)).multiply((24. / math.pi) * 4.92)
    else:
        # ra = (24. / math.pi) * (1367 * 0.0036) * _dr(doy) * theta
        ra = theta.multiply(_dr(doy)).multiply((24. / math.pi) * (1367 * 0.0036))
    return ra


def _ra_hourly(lat, lon, doy, time_mid, method='asce'):
    """Hourly extraterrestrial radiation (Eq. 48)

    Parameters
    ----------
    lat : ee.Image or ee.Number
        Latitude [radians].
    lon : ee.Image or ee.Number
        Longitude [radians].
    doy : ee.Image or ee.Number
        Day of year.
    time_mid : ee.Image or ee.Number
        UTC time at midpoint of period [hours].
    method : {'asce' (default), 'refet'}, optional
        Calculation method:
        * 'asce' -- Calculations will follow ASCE-EWRI 2005 [1] equations.
        * 'refet' -- Calculations will follow RefET software.

    Returns
    -------
    ra : ee.Image or ee.Number
        Hourly extraterrestrial radiation [MJ m-2 h-1].

    Notes
    -----
    Equation in ASCE-EWRI 2005 uses a solar constant of ~1366.666... W m-2
    Equation in Duffie & Beckman (?) uses a solar constant of 1367 W m-2

    """
    omega = _omega(_solar_time_rad(lon, time_mid, _seasonal_correction(doy)))
    delta = _delta(doy, method)
    omegas = _omega_sunset(lat, delta)

    # Solar time as start and end of period (Eqs. 53 & 54)
    # Modify omega1 and omega2 at sunrise and sunset (Eq. 56)
    omega1 = omega.subtract(math.pi / 24).max(omegas.multiply(-1)).min(omegas)
    omega2 = omega.add(math.pi / 24).max(omegas.multiply(-1)).min(omegas)
    omega1 = omega1.min(omega2)

    # Extraterrestrial radiation (Eq. 48)
    theta = omega2.subtract(omega1).multiply(lat.sin()).multiply(delta.sin())\
        .add(lat.cos().multiply(delta.cos()).multiply(omega2.sin().subtract(omega1.sin())))
    if method == 'asce':
        # ra = (12. / math.pi) * 4.92 * _dr(doy) * theta
        ra = theta.multiply(_dr(doy)).multiply((12. / math.pi) * 4.92)
    else:
        # ra = (12. / math.pi) * (1367 * 0.0036) * _dr(doy) * theta
        ra = theta.multiply(_dr(doy)).multiply((12. / math.pi) * (1367 * 0.0036))
    return ra


def _rso_daily(ea, ra, pair, doy, lat):
    """Full daily clear sky solar radiation formulation (Appendix D)

    Parameters
    ----------
    ea : ee.Image or ee.Number
        Actual vapor pressure [kPa].
    ra : ee.Image or ee.Number
        Extraterrestrial radiation [MJ m-2 d-1].
    pair : ee.Image or ee.Number
        Air pressure [kPa].
    doy : ee.Image or ee.Number
        Day of year.
    lat : ee.Image or ee.Number
        Latitude [rad].

    Returns
    -------
    rso : ee.Image or ee.Number
        Daily clear sky solar radiation [MJ m-2 d-1].
        Output data type will match "ea" data type.

    """
    # sin of the angle of the sun above the horizon (D.5 and Eq. 62)
    sin_beta_24 = _doy_fraction(doy).subtract(1.39).sin().multiply(lat).multiply(0.3)\
        .add(0.85).subtract(lat.pow(2).multiply(0.42)).sin().max(0.1)

    # Precipitable water
    w = _precipitable_water(ea, pair)

    # Clearness index for direct beam radiation (Eq. D.2)
    # Limit sin_beta >= 0.01 so that KB does not go undefined
    kb = w.divide(sin_beta_24).pow(0.4).multiply(-0.075)\
        .add(pair.multiply(-0.00146).divide(sin_beta_24))\
        .exp().multiply(0.98)
    # kb = ea.expression(
    #     '0.98 * exp((-0.00146 * pair) / sin_beta_24 - '
    #     '           0.075 * (w / sin_beta_24) ** 0.4)',
    #     {'pair': pair, 'sin_beta_24': sin_beta_24, 'w': w})

    # Transmissivity index for diffuse radiation (Eq. D.4)
    kd = kb.multiply(-0.36).add(0.35).min(kb.multiply(0.82).add(0.18))

    # print('{:>10s}: {:>8.3f}'.format('sin_beta_24', float(sin_beta_24)))
    # print('{:>10s}: {:>8.3f}'.format('w', float(w)))
    # print('{:>10s}: {:>8.3f}'.format('kb', float(kb)))
    # print('{:>10s}: {:>8.3f}'.format('kd', float(kd)))

    rso = kb.add(kd).multiply(ra)
    return rso


def _rso_hourly(ea, ra, pair, doy, time_mid, lat, lon, method='asce'):
    """Full hourly clear sky solar radiation formulation (Appendix D)

    Parameters
    ----------
    ea : ee.Image or ee.Number
        Actual vapor pressure [kPa].
    ra : ee.Image or ee.Number
        Extraterrestrial radiation [MJ m-2 h-1].
    pair : ee.Image or ee.Number
        Air pressure [kPa].
    doy : ee.Image or ee.Number
        Day of year.
    time_mid : ee.Image or ee.Number
        UTC time at midpoint of period [hours].
    lat : ee.Image or ee.Number
        Latitude [rad].
    lon : ee.Image or ee.Number
        Longitude [rad].
    method : {'asce' (default), 'refet'}, optional
        Calculation method:
        * 'asce' -- Calculations will follow ASCE-EWRI 2005 [1] equations.
        * 'refet' -- Calculations will follow RefET software.
        Passed through to declination calculation (_delta()).

    Returns
    -------
    rso : ee.Image or ee.Number
        Hourly clear sky solar radiation [MJ m-2 h-1].
        Output data type will match "ra" data type.

    """
    sc = _seasonal_correction(doy)
    omega = _omega(_solar_time_rad(lon, time_mid, sc))

    # sin of the angle of the sun above the horizon (D.6 and Eq. 62)
    delta = _delta(doy, method)
    sin_beta = lat.sin().multiply(delta.sin())\
        .add(lat.cos().multiply(delta.cos()).multiply(omega.cos()))

    # Precipitable water
    w = _precipitable_water(ea, pair)

    # Clearness index for direct beam radiation (Eq. D.2)
    # Limit sin_beta >= 0.01 so that KB does not go undefined
    kt = 1.0
    kb = w.divide(sin_beta.max(0.01)).pow(0.4).multiply(-0.075)\
        .add(pair.multiply(-0.00146).divide(sin_beta.max(0.01).multiply(kt)))\
        .exp().multiply(0.98)
    # kb = ea.expression(
    #     '0.98 * exp((-0.00146 * pair) / (kt * sin_beta) - '
    #     '           0.075 * (w / sin_beta) ** 0.4))',
    #     {'pair': pair, 'kt': kt, 'sin_beta': sin_beta.max(0.01), 'w': w})

    # Transmissivity index for diffuse radiation (Eq. D.4)
    kd = kb.multiply(-0.36).add(0.35).min(kb.multiply(0.82).add(0.18))

    rso = kb.add(kd).multiply(ra)
    return rso


def _rso_simple(ra, elev):
    """Simplified daily/hourly clear sky solar formulation (Eqs. 19 & 45)

    Parameters
    ----------
    ra : ee.Image or ee.Number
        Extraterrestrial radiation [MJ m-2 d-1 or MJ m-2 h-1].
    elev : ee.Image or ee.Number
        Elevation [m].

    Returns
    -------
    rso : ee.Image or ee.Number
        Clear sky solar radiation [MJ m-2 d-1 or MJ m-2 h-1].
        Output data type will match "ra" data type.

    Notes
    -----
    rso = (0.75 + 2E-5 * elev) * ra

    """
    return ra.multiply(elev.multiply(2E-5).add(0.75))


def _fcd_daily(rs, rso):
    """Daytime cloudiness fraction (Eq. 18)

    Parameters
    ----------
    rs : ee.Image or ee.Number
        Measured solar radiation [MJ m-2 d-1].
    rso : ee.Image or ee.Number
        Clear sky solar radiation [MJ m-2 d-1].

    Returns
    -------
    fcd : ee.Image or ee.Number
        Output data type will match "rs" data type.

    Notes
    -----
    fcd = 1.35 * min(max(rs / rso, 0.3), 1.0) - 0.35

    """
    return rs.divide(rso).max(0.3).min(1.0).multiply(1.35).subtract(0.35)


def _fcd_hourly(rs, rso, doy, time_mid, lat, lon, method='asce'):
    """Cloudiness fraction (Eq. 45)

    Parameters
    ----------
    rs : ee.Image or ee.Number
        Measured solar radiation [MJ m-2 h-1].
    rso : ee.Image or ee.Number
        Clear sky solar radiation [MJ m-2 h-1].
    doy : ee.Image or ee.Number
        Day of year.
    time_mid : ee.Image or ee.Number
        UTC time at midpoint of period [hours].
    lat : ee.Image or ee.Number
        Latitude [rad].
    lon : ee.Image or ee.Number
        Longitude [rad].
    method : {'asce' (default), 'refet'}, optional
        Calculation method:
        * 'asce' -- Calculations will follow ASCE-EWRI 2005 [1] equations.
        * 'refet' -- Calculations will follow RefET software.
        Passed through to declination calculation (_delta()).

    Returns
    -------
    fcd : ee.Image or ee.Number
        Output data type will match "rs" data type.

    """
    # DEADBEEF - These values are only needed for identifying low sun angles
    sc = _seasonal_correction(doy)
    delta = _delta(doy, method)
    omega = _omega(_solar_time_rad(lon, time_mid, sc))
    beta = lat.sin().multiply(delta.sin())\
        .add(lat.cos().multiply(delta.cos()).multiply(omega.cos())).asin()

    fcd = rs.divide(rso).max(0.3).min(1).multiply(1.35).subtract(0.35)

    # Intentionally not using where() so that function will work with ee.Number
    fcd = fcd.max(beta.lt(0.3))
    # fcd = fcd.where(beta.lt(0.3), 1)

    # # DEADBEEF - Code from NumPy functions
    # fcd = np.ones(rso.shape)
    # fcd[rso > 0] = 1.35 * np.clip(rs[rso > 0] / rso[rso > 0], 0.3, 1) - 0.35
    #
    # # For now just set fcd to 1 for low sun angles
    # # DEADBEEF - Still need to get daytime value of fcd when beta > 0.3
    # #   Get closest value in time (array space) when beta > 0.3
    # fcd[beta < 0.3] = 1

    return fcd


def _rnl_daily(tmax, tmin, ea, fcd):
    """Daily net long-wave radiation  (Eq. 17)

    Parameters
    ----------
    tmax : ee.Image or ee.Number
        Daily maximum air temperature [C].
    tmin : ee.Image or ee.Number
        Daily minimum air temperature [C].
    ea : ee.Image or ee.Number
        Actual vapor pressure [kPa].
    fcd : ee.Image or ee.Number
        cloudiness fraction.

    Returns
    -------
    ee.Image or ee.Number
        Daily net long-wave radiation [MJ m-2 d-1].
        Output data type will match "tmax" data type.

    Notes
    -----
    rnl = 4.901E-9 * fcd * (0.34 - 0.14 * sqrt(ea)) *
          0.5 * ((tmax + 273.16) ** 4 + (tmin + 273.16) ** 4))

    """
    return tmax.add(273.16).pow(4).add(tmin.add(273.16).pow(4)).multiply(0.5)\
        .multiply(ea.sqrt().multiply(-0.14).add(0.34))\
        .multiply(fcd).multiply(4.901E-9)


def _rnl_hourly(tmean, ea, fcd):
    """Hourly net long-wave radiation  (Eq. 44)

    Parameters
    ----------
    tmean : ee.Image or ee.Number
        Mean hourly air temperature [C].
    ea : ee.Image or ee.Number
        Actual vapor pressure [kPa].
    fcd : ee.Image or ee.Number
        Cloudiness fraction.

    Returns
    -------
    ee.Image or ee.Number
        Hourly net long-wave radiation [MJ m-2 h-1].
        Output data type will match "tmean" data type.

    Notes
    -----
    rnl = 2.042E-10 * fcd * (0.34 - 0.14 * sqrt(ea)) * ((tmean + 273.16) ** 4)

    """
    return tmean.add(273.16).pow(4)\
        .multiply(ea.sqrt().multiply(-0.14).add(0.34))\
        .multiply(fcd).multiply(2.042E-10)


def _rn(rs, rnl):
    """Net daily/hourly radiation (Eqs. 15 & 16)

    Parameters
    ----------
    rs : ee.Image or ee.Number
        Measured solar radiation [MJ m-2 d-1 or MJ m-2 h-1].
    rnl : ee.Image or ee.Number
        Hourly net long-wave radiation [MJ m-2 d-1 or MJ m-2 h-1].

    Returns
    -------
    ee.Image or ee.Number
        Hourly net long-wave radiation [MJ m-2 d-1 or MJ m-2 h-1].
        Output data type will match "rnl" data type.

    Notes
    -----
    Switching calculation to work from rnl (which is computed from temperature)
    rnl = 0.77 * rs - rnl

    """
    return rnl.multiply(-1).add(rs.multiply(0.77))
    # return rs.multiply(0.77).subtract(rnl)


def _wind_height_adjust(uz, zw):
    """Wind speed at 2 m height based on full logarithmic profile (Eq. 33)

    Parameters
    ----------
    uz : ee.Image or ee.Number
        Wind speed at measurement height [m s-1].
    zw : ee.Image or ee.Number
        Wind measurement height [m].

    Returns
    -------
    ee.Image or ee.Number
        Wind speed at 2 m height [m s-1].

    Notes
    -----
    u2 = uz * 4.87 / log(67.8 * zw - 5.42)

    """
    return uz.multiply(4.87).divide(zw.multiply(67.8).subtract(5.42).log())
