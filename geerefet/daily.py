import ee

from . import calcs

ee.Initialize()


class Daily():
    """"""

    def __init__(self, tmax, tmin, ea, rs, uz, zw, elev, lat, doy,
                 method='asce', rso_type=None, rso=None):
        """ASCE Daily Standardized Reference Evapotranspiration (ET)

        Arguments
        ---------
        tmax : ee.Image or ee.Number
            Maximum daily temperature [C].
        tmin : ee.Image or ee.Number
            Minimum daily temperature [C].
        ea : ee.Image or ee.Number
            Actual vapor pressure [kPa].
        rs : ee.Image or ee.Number
            Incoming shortwave solar radiation [MJ m-2 day-1].
        uz : ee.Image or ee.Number
            Wind speed [m/s].
        zw : ee.Number or float
            Wind speed height [m].
        elev : ee.Image or ee.Number
            Elevation [m].
        lat : ee.Image or ee.Number
            Latitude [radians].
        doy : ee.Number
            Day of year.
        method : {'asce' (default), 'refet'}, optional
            Specifies which calculation method to use.
            * 'asce' -- Calculations will follow ASCE-EWRI 2005 [1].
            * 'refet' -- Calculations will follow RefET software.
        rso_type : {None (default), 'full' , 'simple', 'array'}, optional
            Specifies which clear sky solar radiation (Rso) model to use.
            * None -- Rso type will be determined from "method" parameter
            * 'full' -- Full clear sky solar formulation
            * 'simple' -- Simplified clear sky solar formulation
            * 'array' -- Read Rso values from "rso" function parameter
        rso : ee.Image, ee.Number, or None, optional
            Clear sky solar radiation [MJ m-2 day-1] (the default is None).
            Only used if rso_type == 'array'.

        Raises
        ------
        ValueError
            If 'method' or 'rso_type' parameter is invalid.

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

        if rso_type is None:
            pass
        elif rso_type.lower() not in ['simple', 'full', 'array']:
            raise ValueError(
                'rso_type must be None, "simple", "full", or "array')
        elif rso_type.lower() in 'array':
            # Check that rso is an ee.Image or ee.Number?
            pass

        # Do these all need to be set onto self?
        self.tmin = tmin
        self.tmax = tmax
        self.ea = ea
        self.rs = rs
        self.uz = uz
        self.zw = zw
        self.elev = elev
        self.lat = lat
        self.doy = doy

        # To match standardized form, pair is calculated from elevation
        self.pair = calcs._air_pressure(self.elev, method)
        self.psy = self.pair.multiply(0.000665)

        self.tmean = self.tmax.add(self.tmin).multiply(0.5)
        self.es_slope = calcs._es_slope(self.tmean, method)

        # Saturated vapor pressure
        self.es = calcs._sat_vapor_pressure(self.tmax).add(
            calcs._sat_vapor_pressure(self.tmin)) \
            .multiply(0.5)

        # Vapor pressure deficit
        self.vpd = calcs._vpd(self.es, self.ea)

        # Extraterrestrial radiation
        self.ra = calcs._ra_daily(self.lat, self.doy, method)

        # Clear sky solar radiation
        # If rso_type is not set, use the method
        # If rso_type is set, use rso_type directly
        if rso_type is None:
            if method.lower() == 'asce':
                self.rso = calcs._rso_simple(self.ra, self.elev)
            elif method.lower() == 'refet':
                self.rso = calcs._rso_daily(
                    self.ra, self.ea, self.pair, self.doy, self.lat)
        elif rso_type.lower() == 'simple':
            self.rso = calcs._rso_simple(self.ra, self.elev)
        elif rso_type.lower() == 'full':
            self.rso = calcs._rso_daily(
                self.ra, self.ea, self.pair, self.doy, self.lat)
        elif rso_type.lower() == 'array':
            # Use rso array passed to function
            self.rso = rso

        # Cloudiness fraction
        self.fcd = calcs._fcd_daily(self.rs, self.rso)

        # Net long-wave radiation
        self.rnl = calcs._rnl_daily(self.tmax, self.tmin, self.ea, self.fcd)

        # Net radiation (Eqs. 15 and 16)
        self.rn = self.rs.multiply(0.77).subtract(self.rnl)

        # Wind speed
        self.u2 = calcs._wind_height_adjust(self.uz, self.zw)

    def eto(self):
        """Grass reference surface"""
        self.cn = 900
        self.cd = 0.34
        return self._etsz()

    def etr(self):
        """Alfalfa reference surface"""
        self.cn = 1600
        self.cd = 0.38
        return self._etsz()

    def _etsz(self):
        """Daily reference ET (Eq. 1)

        Returns
        -------
        etsz : ee.Image or ee.Number
            Standardized reference ET [mm].

        """
        return self.u2.multiply(self.vpd).multiply(self.cn).multiply(self.psy)\
            .divide(self.tmean.add(273))\
            .add(self.es_slope.multiply(self.rn).multiply(0.408))\
            .divide(self.u2.multiply(self.cd).add(1)\
                        .multiply(self.psy).add(self.es_slope))
        # return self.tmean.expression(
        #     '(0.408 * es_slope * rn + (psy * cn * u2 * vpd / (tmean + 273))) / '
        #     '(es_slope + psy * (cd * u2 + 1)))',
        #     {'cd': self.cd, 'cn': self.cn, 'es_slope': self.es_slope,
        #      'psy': self.psy, 'rn': self.rn, 'tmean': self.tmean,
        #      'u2': self.u2, 'vpd': self.vpd})
        # return (
        #     (0.408 * es_slope * rn + (psy * cn * u2 * vpd / (tmean + 273))) /
        #     (es_slope + psy * (cd * u2 + 1)))

    @classmethod
    def gridmet(cls, gridmet_img, zw=10, elev=None, lat=None, method='asce',
                rso_type=None):
        """Initialize daily RefET from a GRIDMET image

        Parameters
        ----------
        gridmet_img : ee.Image
            GRIDMET image from the collection IDAHO_EPSCOR/GRIDMET.
        zw : ee.Number or float
            Wind speed height [m] (the default is 10).
        elev : ee.Image or ee.Number, optional
            Elevation image [m].  The standard GRIDMET elevation image
            (projects/climate-engine/gridmet/elevtion) will be used if not set.
        lat : ee.Image or ee.Number
            Latitude image [degrees].  The latitude will be computed
            dynamically using ee.Image.pixelLonLat() if not set.
        method : {'asce' (default), 'refet'}, optional
            Specifies which calculation method to use.
            * 'asce' -- Calculations will follow ASCE-EWRI 2005.
            * 'refet' -- Calculations will follow RefET software.
        rso_type : {None (default), 'full' , 'simple'}, optional
            Specifies which clear sky solar radiation (Rso) model to use.
            * None -- Rso type will be determined from "method" parameter
            * 'full' -- Full clear sky solar formulation
            * 'simple' -- Simplified clear sky solar formulation

        Notes
        -----
        Temperatures are converted from K to C.
        Solar radiation is converted from W m-2 to MJ m-2 day-1.
        Actual vapor pressure is computed from specific humidity (GRIDMET sph)
            and air pressure (from elevation).

        """

        if elev is None:
            elev = ee.Image('projects/climate-engine/gridmet/elevation')
        if lat is None:
            lat = ee.Image.pixelLonLat().select('latitude')

        return cls(
            tmin=gridmet_img.select(['tmmn']).subtract(273.15),
            tmax=gridmet_img.select(['tmmx']).subtract(273.15),
            ea=calcs._actual_vapor_pressure(
                pair=calcs._air_pressure(elev, method),
                q=gridmet_img.select(['sph'])),
            rs=gridmet_img.select(['srad']).multiply(0.0864),
            uz=gridmet_img.select(['vs']),
            zw=ee.Number(zw),
            elev=elev,
            lat=lat,
            doy=ee.Number(ee.Date(gridmet_img.get('system:time_start'))
                          .getRelative('day', 'year')).add(1).double(),
            method=method,
            rso_type=rso_type)