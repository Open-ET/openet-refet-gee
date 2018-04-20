import math

import ee

from . import calcs

ee.Initialize()


class Hourly():
    def __init__(self, tmean, ea, rs, uz, zw, elev, lat, lon, doy, time,
                 method='asce'):
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
        zw : ee.Number
            Wind speed measurement/estimated height [m].
        elev : ee.Image or ee.Number
            Elevation [m]
        lat : ee.Image or ee.Number
            Latitude [degrees]
        lon : ee.Image or ee.Number
            Longitude [degrees].
        doy : ee.Number
            Day of year.
        time : ee.Number
            UTC hour at start of time period.
        method : {'asce' (default), 'refet'}, optional
            Specifies which calculation method to use.
            * 'asce' -- Calculations will follow ASCE-EWRI 2005 [1].
            * 'refet' -- Calculations will follow RefET software.

        Raises
        ------
        ValueError
            If 'method' parameter is invalid.

        Notes
        -----
        Divide solar radiation values by 0.0036 to convert MJ m-2 hr-1 to W m-2.
        Latitude & longitude units are degress, not radians.

        References
        ----------
        .. [1] ASCE-EWRI (2005). The ASCE standardized reference evapotranspiration
            equation. ASCE-EWRI Standardization of Reference Evapotranspiration
            Task Committee Rep., ASCE Reston, Va.
            http://www.kimberly.uidaho.edu/water/asceewri/ascestzdetmain2005.pdf
            http://www.kimberly.uidaho.edu/water/asceewri/appendix.pdf

        """
        if method.lower() not in ['asce', 'refet']:
            raise ValueError('method must be "asce" or "refet"')

        # Get time_start from tmin
        # Shoudl time_start be set in init?
        self.time_start = ee.Image(tmean).get('system:time_start')
        self.date = ee.Date(self.time_start)

        # Do these all need to be set onto self?
        self.tmean = tmean
        self.ea = ea
        self.rs = rs
        self.uz = uz
        self.zw = zw
        self.elev = elev
        self.lat = lat.multiply(math.pi / 180)
        self.lon = lon.multiply(math.pi / 180)
        self.doy = doy
        self.time = time

        # To match standardized form, psy is calculated from elevation based pair
        self.pair = calcs._air_pressure(self.elev, method)
        self.psy = self.pair.multiply(0.000665)

        self.es = calcs._sat_vapor_pressure(self.tmean)
        self.es_slope = calcs._es_slope(self.tmean, method)

        # Vapor pressure deficit
        self.vpd = self.es.subtract(self.ea)
        # self.vpd = calcs._vpd(self.es, self.ea)

        # Extraterrestrial radiation
        time_mid = self.time.add(0.5)
        self.ra = calcs._ra_hourly(
            self.lat, self.lon, self.doy, time_mid, method)

        # Clear sky solar radiation
        if method == 'asce':
            self.rso = calcs._rso_simple(self.ra, self.elev)
        elif method == 'refet':
            self.rso = calcs._rso_hourly(
                self.ra, self.ea, self.pair, self.doy, time_mid,
                self.lat, self.lon, method)

        # Cloudiness fraction
        # Intentionally not using time_mid to match Beta value in IN2 file
        # In IN2, "Beta" is computed for the start of the time period,
        #   but "SinBeta" is computed for the midpoint.
        # Beta (not SinBeta) is used for clamping fcd.
        self.fcd = calcs._fcd_hourly(
            self.rs, self.rso, self.doy, self.time, self.lat, self.lon, method)

        # Net long-wave radiation
        self.rnl = calcs._rnl_hourly(self.tmean, self.ea, self.fcd)

        # Net radiation (Eqs. 42 and 43)
        self.rn = self.rs.multiply(0.77).subtract(self.rnl)

        # Wind speed
        self.u2 = calcs._wind_height_adjust(self.uz, self.zw)

    def _etsz(self):
        """Hourly reference ET (Eq. 1)


        Returns
        -------
        etsz : ee.Image
            Standardized reference ET [mm].

        """
        # return self.u2.multiply(self.vpd).multiply(self.cn).multiply(self.psy)\
        #     .divide(self.tmean.add(273))\
        #     .add(self.es_slope.multiply(self.rn.subtract(self.g)).multiply(0.408))\
        #     .divide(self.u2.multiply(self.cd).add(1)\
        #                 .multiply(self.psy).add(self.es_slope))

        return self.tmean.expression(
            '(0.408 * es_slope * (rn - g) + (psy * cn * u2 * vpd / (tmean + 273))) / '
            '(es_slope + psy * (cd * u2 + 1))',
            {'cd': self.cd, 'cn': self.cn, 'es_slope': self.es_slope,
             'g': self.g, 'psy': self.psy, 'rn': self.rn, 'tmean': self.tmean,
             'u2': self.u2, 'vpd': self.vpd})

    def eto(self):
        """Short (grass) reference surface"""
        self.cn = ee.Number(37.0)
        cd_day = 0.24
        g_rn_day = 0.1
        cd_night = 0.96
        g_rn_night = 0.5

        # Adjust coefficients for daytime/nighttime
        # Nighttime is defined as when Rn < 0 (pg 44)
        self.cd = self.rn.multiply(0).add(cd_day).where(self.rn.lt(0), cd_night)
        self.g_rn = self.rn.multiply(0).add(g_rn_day)\
            .where(self.rn.lt(0), g_rn_night)
        # self.cd = ee.Image.constant(cd_day).where(self.rn.lt(0), cd_night)
        # self.g_rn = ee.Image.constant(g_rn_day).where(self.rn.lt(0), g_rn_night)

        # Soil heat flux (Eqs. 65 and 66)
        self.g = self.rn.multiply(self.g_rn)

        return ee.Image(self._etsz().rename(['eto'])
            .set('system:time_start', self.time_start))

    def etr(self):
        """Tall (alfalfa) reference surface"""
        self.cn = ee.Number(66.0)
        cd_day = 0.25
        g_rn_day = 0.04
        cd_night = 1.7
        g_rn_night = 0.2

        # Adjust coefficients for daytime/nighttime
        # Nighttime is defined as when Rn < 0 (pg 44)
        self.cd = self.rn.multiply(0).add(cd_day).where(self.rn.lt(0), cd_night)
        self.g_rn = self.rn.multiply(0).add(g_rn_day)\
            .where(self.rn.lt(0), g_rn_night)
        # self.cd = ee.Image.constant(cd_day).where(self.rn.lt(0), cd_night)
        # self.g_rn = ee.Image.constant(g_rn_day).where(self.rn.lt(0), g_rn_night)

        # Soil heat flux (Eqs. 65 and 66)
        self.g = self.rn.multiply(self.g_rn)

        return ee.Image(self._etsz().rename(['etr'])
            .set('system:time_start', self.time_start))

