import math

import ee

from . import calcs


class Hourly():
    def __init__(self, tmean, ea, rs, uz, zw, elev, lat, lon, doy, time,
                 method='asce'):
        """ASCE Hourly Standardized Reference Evapotranspiration (ET)

        .. warning:: Cloudiness fraction at night is not being computed per [1]_

        Arguments
        ---------
        tmean : ee.Image
            Average hourly temperature [C].
        ea : ee.Image
            Actual vapor pressure [kPa].
        rs : ee.Image
            Shortwave solar radiation [MJ m-2 hr-1].
        uz : ee.Image
            Wind speed [m s-1].
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
        self.pair = calcs._air_pressure(self.elev, method=method)

        # Psychrometric constant (Eq. 35)
        self.psy = self.pair.multiply(0.000665)

        self.es = calcs._sat_vapor_pressure(self.tmean)
        self.es_slope = calcs._es_slope(self.tmean, method)

        # Vapor pressure deficit
        self.vpd = self.es.subtract(self.ea)
        # self.vpd = calcs._vpd(es=self.es, ea=self.ea)

        # Extraterrestrial radiation
        time_mid = self.time.add(0.5)
        self.ra = calcs._ra_hourly(
            lat=self.lat, lon=self.lon, doy=self.doy, time_mid=time_mid,
            method=method)

        # Clear sky solar radiation
        if method == 'asce':
            self.rso = calcs._rso_simple(ra=self.ra, elev=self.elev)
        elif method == 'refet':
            self.rso = calcs._rso_hourly(
                ea=self.ea, ra=self.ra, pair=self.pair, doy=self.doy,
                time_mid=time_mid, lat=self.lat, lon=self.lon, method=method)

        # Cloudiness fraction
        # Intentionally not using time_mid to match Beta value in IN2 file
        # In IN2, "Beta" is computed for the start of the time period,
        #   but "SinBeta" is computed for the midpoint.
        # Beta (not SinBeta) is used for clamping fcd.
        self.fcd = calcs._fcd_hourly(
            rs=self.rs, rso=self.rso, doy=self.doy, time_mid=self.time,
            lat=self.lat, lon=self.lon, method=method)

        # Net long-wave radiation
        self.rnl = calcs._rnl_hourly(tmean=self.tmean, ea=self.ea, fcd=self.fcd)

        # Net radiation
        self.rn = calcs._rn(self.rs, self.rnl)

        # Wind speed
        self.u2 = calcs._wind_height_adjust(uz=self.uz, zw=self.zw)

    def etsz(self, surface):
        """Standardized reference ET

        Parameters
        ----------
        surface : {'alfalfa', 'etr', 'tall', 'grass', 'eto', 'short'}
            Reference surface type.

        Returns
        -------
        ee.Image

        """
        if surface.lower() in ['alfalfa', 'etr', 'tall']:
            return self.etr()
        elif surface.lower() in ['grass', 'eto', 'short']:
            return self.eto()
        else:
            raise ValueError('unsupported surface type: {}'.format(surface))

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

    @classmethod
    def nldas(cls, nldas_img, zw=None, elev=None, lat=None, lon=None,
              method='asce'):
        """Initialize daily RefET from an NLDAS image collection

        Parameters
        ----------
        nldas_img : ee.Image
            NLDAS hourly images from the collection NASA/NLDAS/FORA0125_H002.
        zw : ee.Number, optional
            Wind speed height [m] (the default is 10).
        elev : ee.Image or ee.Number, optional
            Elevation image [m].  The standard GRIDMET elevation image
            (projects/climate-engine/gridmet/elevtion) will be used if not set.
        lat : ee.Image or ee.Number
            Latitude image [degrees].  The latitude will be computed
            dynamically using ee.Image.pixelLonLat() if not set.
        lon : ee.Image or ee.Number
            Longitude image [degrees].  The longitude will be computed
            dynamically using ee.Image.pixelLonLat() if not set.
        method : {'asce' (default), 'refet'}, optional
            Specifies which calculation method to use.
            * 'asce' -- Calculations will follow ASCE-EWRI 2005.
            * 'refet' -- Calculations will follow RefET software.

        Notes
        -----
        Solar radiation is converted from W m-2 to MJ m-2 day-1.
        Actual vapor pressure is computed from specific humidity and air
            pressure (from elevation).

        """

        if zw is None:
            zw = ee.Number(10)
        if elev is None:
            # Reproject to the NLDAS grid
            elev = ee.Image("CGIAR/SRTM90_V4") \
                .reproject('EPSG:4326', [0.125, 0, -125, 0, -0.125, 53])
            elev = ee.Image('CGIAR/SRTM90_V4')
        if lat is None:
            # Reproject to the NLDAS grid
            lat = ee.Image.pixelLonLat().select('latitude') \
                .reproject('EPSG:4326', [0.125, 0, -125, 0, -0.125, 53])
            # lat = ee.Image.pixelLonLat().select('latitude')
        if lon is None:
            # Reproject to the NLDAS grid
            lon = ee.Image.pixelLonLat().select('longitude')\
                .reproject('EPSG:4326', [0.125, 0, -125, 0, -0.125, 53])
            # lon = ee.Image.pixelLonLat().select('longitude')
        image_date = ee.Date(nldas_img.get('system:time_start'))

        return cls(
            tmean=nldas_img.select(['temperature']),
            ea=calcs._actual_vapor_pressure(
                pair=calcs._air_pressure(elev, method),
                q=nldas_img.select(['specific_humidity'])),
            rs=nldas_img.select(['shortwave_radiation']).multiply(0.0036),
            uz=nldas_img.select(["wind_u"]).pow(2) \
                .add(nldas_img.select(["wind_v"]).pow(2)) \
                .sqrt().rename(['uz']),
            zw=zw,
            elev=elev,
            lat=lat,
            lon=lon,
            doy=ee.Number(image_date.getRelative('day', 'year')).add(1).double(),
            # time=ee.Number(image_date.getRelative('hour', 'day')),
            time=ee.Number(image_date.get('hour')),
            method=method,
        )
