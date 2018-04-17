import math

import ee
# import numpy as np

from . import calcs


def daily(tmin, tmax, ea, rs, uz, zw, elev, lat, doy, surface,
          method='refet', rso_type=None, rso=None):
    """ASCE Daily Standardized Reference Evapotranspiration (ET)

    Arguments
    ---------
    tmin : ee.Image
        Minimum daily temperature [C].
    tmax : ee.Image
        Maximum daily temperature [C].
    ea : ee.Image
        Actual vapor pressure [kPa].
    rs : ee.Image
        Incoming shortwave solar radiation [MJ m-2 day-1].
    uz : ndarray
        Wind speed [m/s].
    zw : float
        Wind speed height [m].
    elev : ee.Image or ee.Number
        Elevation [m].
    lat : ee.Image or ee.Number
        Latitude [radians].
    doy : ndarray
        Day of year.
    surface : {'eto', 'etr', 'grass', 'alfalfa', 'short', 'tall'}
        Specifies which reference crop surface.
        * 'etr', 'alfalfa', 'tall' -- Tall reference crop
        * 'eto', 'grass', 'short' -- Short reference crop
    method : {'refet', 'asce'}, optional
        Specifies which calculation method to use.
        * 'refet' -- Calculations will follow RefET software (default).
        * 'asce' -- Calculations will follow ASCE-EWRI 2005 [1] equations exactly.
    rso_type : {None (default), 'full' , 'simple', 'array'}, optional
        Specifies which clear sky solar radiation (Rso) model to use.
        * None -- Rso type will be determined from "method" parameter
        * 'full' -- Full clear sky solar formulation
        * 'simple' -- Simplified clear sky solar formulation
        * 'array' -- Read Rso values from "rso" function parameter
    rso : array_like or None, optional
        Clear sky solar radiation [MJ m-2 day-1] (the default is None).
        Only used if rso_type == 'array'.

    Returns
    -------
    etsz : ee.Image or ee.Number
        Standardized reference ET [mm].

    Raises
    ------
    ValueError
        If 'surface', 'method' or 'rso_type' parameter is invalid.
        If latitude values are outside the range [-pi/2, pi/2].

    Notes
    -----
    cn: 900 for ETo, 1600 for ETr
    cd: 0.34 for ETo, 0.38 for ETr
    Divide solar radiation values by 0.0864 to convert MJ m-2 day-1 to W m-2

    References
    ----------
    .. [1] ASCE-EWRI (2005). The ASCE standardized reference evapotranspiration
        equation. ASCE-EWRI Standardization of Reference Evapotranspiration
        Task Committee Rep., ASCE Reston, Va.
        http://www.kimberly.uidaho.edu/water/asceewri/ascestzdetmain2005.pdf
        http://www.kimberly.uidaho.edu/water/asceewri/appendix.pdf
    """

    # # Convert all inputs to NumPy arrays
    # tmin = np.array(tmin, copy=True, ndmin=1)
    # tmax = np.array(tmax, copy=True, ndmin=1)
    # ea = np.array(ea, copy=True, ndmin=1)
    # rs = np.array(rs, copy=True, ndmin=1)
    # uz = np.array(uz, copy=True, ndmin=1)
    # elev = np.array(elev, copy=True, ndmin=1)
    # lat = np.array(lat, copy=True, ndmin=1)

    # # Check that latitudes are in radians
    # if np.any(np.fabs(lat) > (0.5 * math.pi)):
    #     raise ValueError('latitudes must be in radians [-pi/2, pi/2]')

    if method.lower() not in ['asce', 'refet']:
        raise ValueError('method must be "asce" or "refet"')

    if surface.lower() in ['eto', 'grass', 'short']:
        # Tall reference crop parameters
        cn, cd = 900, 0.34
    elif surface.lower() in ['etr', 'alfalfa', 'tall']:
        # Short reference crop parameters
        cn, cd = 1600, 0.38
    else:
        raise ValueError('surface must be "etr" or "eto"')

    if rso_type is None:
        pass
    elif rso_type.lower() not in ['simple', 'full', 'array']:
        raise ValueError('rso_type must be None, "simple", "full", or "array')
    elif rso_type.lower() in 'array':
        # Check that rso is an array
        pass

    # To match standardized form, pair is calculated from elevation
    pair = calcs._air_pressure(elev, method)

    psy = pair.multiply(0.000665)

    tmean = tmax.add(tmin).multiply(0.5)
    es_slope = calcs._es_slope(tmean, method)

    # Saturated vapor pressure
    es = calcs._sat_vapor_pressure(tmax).add(calcs._sat_vapor_pressure(tmin))\
        .multiply(0.5)

    # Vapor pressure deficit
    vpd = calcs._vpd(es, ea)

    # Extraterrestrial radiation
    ra = calcs._ra_daily(lat, doy, method)

    # Clear sky solar radiation
    # If rso_type is not set, use the method
    # If rso_type is set, use rso_type directly
    if rso_type is None :
        if method.lower() == 'asce':
            rso = calcs._rso_simple(ra, elev)
        elif method.lower() == 'refet':
            rso = calcs._rso_daily(ra, ea, pair, doy, lat)
    elif rso_type.lower() == 'simple':
        rso = calcs._rso_simple(ra, elev)
    elif rso_type.lower() == 'full':
        rso = calcs._rso_daily(ra, ea, pair, doy, lat)
    elif rso_type.lower() == 'array':
        # Use rso array passed to function
        pass

    # Cloudiness fraction
    fcd = calcs._fcd_daily(rs, rso)

    # Net long-wave radiation
    rnl = calcs._rnl_daily(tmax, tmin, ea, fcd)

    # Net radiation (Eqs. 15 and 16)
    rn = rs.multiply(0.77).subtract(rnl)

    # Wind speed
    u2 = calcs._wind_height_adjust(uz, zw)

    # Daily reference ET (Eq. 1)
    etsz = tmean.expression(
        '(0.408 * es_slope * rn + (psy * cn * u2 * vpd / (tmean + 273))) / '
        '(es_slope + psy * (cd * u2 + 1)))',
        {'cd': cd, 'cn': cn, 'es_slope': es_slope, 'psy': psy, 'rn': rn,
         'tmean': tmean, 'u2': u2, 'vpd': vpd})
    # etsz = (
    #     (0.408 * es_slope * rn + (psy * cn * u2 * vpd / (tmean + 273))) /
    #     (es_slope + psy * (cd * u2 + 1)))

    return etsz
