import math

import ee
# import numpy as np

from . import calcs


def hourly(tmean, ea, rs, uz, zw, elev, lat, lon, doy, time, surface,
           method='refet'):
    """ASCE Hourly Standardized Reference Evapotranspiration (ET)

    .. warning:: Cloudiness fraction at night is not being computed correctly

    Arguments
    ---------
    tmean : ee.Image
        Average hourly temperature [C].
    ea : ee.Image
        Actual vapor pressure [kPa].
    rs : ee.Image
        Shortwave solar radiation [MJ m-2 hr-1].
    uz : ee.Image
        Wind speed [m/s].
    zw : float
        Wind speed measurement/estimated height [m].
    elev : ee.Image or ee.Number
        Elevation [m]
    lat : ee.Image
        Latitude [radians]
    lon : ndarray
        Longitude [radians].
    doy : ndarray
        Day of year.
    time : ndarray
        UTC hour at start of time period.
    surface : {'eto', 'etr', 'grass', 'alfalfa', 'short', 'tall'}
        Specifies which reference crop surface to use.
        * 'etr', 'alfalfa', 'tall' -- Tall reference crop
        * 'eto', 'grass', 'short' -- Short reference crop
    method : {'refet' (default), 'asce'}, optional
        Specifies which calculation method to use.
        * 'refet' -- Calculations will follow RefET software.
        * 'asce' -- Calculations will follow ASCE-EWRI 2005 [1] equations exactly.

    Returns
    -------
    etsz : ee.Image or ee.Number
        Standardized reference ET [mm].

    Raises
    ------
    ValueError
        If 'surface' or 'method' parameter is invalid.
        If latitude values are outside the range [-pi/2, pi/2].
        If longitude values are outside the range [-pi, pi].

    Notes
    -----
    Divide solar radiation values by 0.0036 to convert MJ m-2 hr-1 to W m-2

    References
    ----------
    .. [1] ASCE-EWRI (2005). The ASCE standardized reference evapotranspiration
        equation. ASCE-EWRI Standardization of Reference Evapotranspiration
        Task Committee Rep., ASCE Reston, Va.
        http://www.kimberly.uidaho.edu/water/asceewri/ascestzdetmain2005.pdf
        http://www.kimberly.uidaho.edu/water/asceewri/appendix.pdf
    """

    # # Convert all inputs to NumPy arrays
    # tmean = np.array(tmean, copy=True, ndmin=1)
    # ea = np.array(ea, copy=True, ndmin=1)
    # rs = np.array(rs, copy=True, ndmin=1)
    # uz = np.array(uz, copy=True, ndmin=1)
    # elev = np.array(elev, copy=True, ndmin=1)
    # lat = np.array(lat, copy=True, ndmin=1)
    # lon = np.array(lon, copy=True, ndmin=1)
    # doy = np.array(doy, copy=True, ndmin=1)
    # time = np.array(time, copy=True, ndmin=1)
    time_mid = time + 0.5

    # # Check that latitude & longitude are in radians
    # if np.any(np.fabs(lat) > (0.5 * math.pi)):
    #     raise ValueError('latitudes must be in radians [-pi/2, pi/2]')
    # elif np.any(np.fabs(lon) > math.pi):
    #     raise ValueError('longitudes must be in radians [-pi, pi]')

    if method.lower() not in ['asce', 'refet']:
        raise ValueError('method must be "asce" or "refet"')

    if surface.lower() in ['eto', 'grass', 'short']:
        # Short reference crop parameters
        cn_day = 37.0
        cd_day = 0.24
        g_rn_day = 0.1
        cn_night = 37.0
        cd_night = 0.96
        g_rn_night = 0.5
    elif surface.lower() in ['etr', 'alfalfa', 'tall']:
        # Tall reference crop parameters
        cn_day = 66.0
        cd_day = 0.25
        g_rn_day = 0.04
        cn_night = 66.0
        cd_night = 1.7
        g_rn_night = 0.2
    else:
        raise ValueError('surface must be "etr" or "eto"')

    # To match standardized form, psy is calculated from elevation based pair
    pair = calcs._air_pressure(elev, method)
    psy = 0.000665 * pair
    es = calcs._sat_vapor_pressure(tmean)
    es_slope = calcs._es_slope(tmean, method)

    # Extraterrestrial radiation
    ra = calcs._ra_hourly(lat, lon, doy, time_mid, method)

    # Clear sky solar radiation
    if method == 'asce':
        rso = calcs._rso_simple(ra, elev)
    elif method == 'refet':
        rso = calcs._rso_hourly(ra, ea, pair, doy, time_mid, lat, lon, method)

    # Cloudiness fraction
    # Intentionally not using time_mid to match Beta value in IN2 file
    # In IN2, "Beta" is computed for the start of the time period,
    #   but "SinBeta" is computed for the midpoint.
    # Beta (not SinBeta) is used for clamping fcd.
    fcd = calcs._fcd_hourly(rs, rso, doy, time, lat, lon, method)

    # Net long-wave radiation
    rnl = calcs._rnl_hourly(tmean, ea, fcd)

    # Net radiation (Eqs. 42 and 43)
    rn = rs.multiply(0.77).subtract(rnl)

    # Adjust coefficients for daytime/nighttime
    # Nighttime is defined as when Rn < 0 (pg 44)
    cn = np.zeros(rn.shape)
    cd = np.zeros(rn.shape)
    g_rn = np.zeros(rn.shape)
    cn[:] = cn_day
    cd[:] = cd_day
    g_rn[:] = g_rn_day
    rn_mask = rn < 0
    cn[rn_mask] = cn_night
    cd[rn_mask] = cd_night
    g_rn[rn_mask] = g_rn_night

    # Soil heat flux (Eqs. 65 and 66)
    g = rn.multiply(g_rn)

    # Wind speed
    u2 = calcs._wind_height_adjust(uz, zw)

    # Hourly reference ET (Eq. 1)
    etsz = tmean.expression(
        '(0.408 * es_slope * (rn - g) + (psy * cn * u2 * (es - ea) / (tmean + 273))) / '
        '(es_slope + psy * (cd * u2 + 1)))',
        {'cd': cd, 'cn': cn, 'ea': ea, 'es': es, 'es_slope': es_slope, 'g': g,
         'psy': psy, 'rn': rn, 'tmean': tmean, 'u2': u2})
    # etsz = (
    #     (0.408 * es_slope * (rn - g) + (psy * cn * u2 * (es - ea) / (tmean + 273))) /
    #     (es_slope + psy * (cd * u2 + 1)))

    return etsz
