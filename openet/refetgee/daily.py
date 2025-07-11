import math

import ee

from . import calcs


def lazy_property(fn):
    """Decorator that makes a property lazy-evaluated

    https://stevenloria.com/lazy-properties/
    """
    attr_name = '_lazy_' + fn.__name__

    @property
    def _lazy_property(self):
        if not hasattr(self, attr_name):
            setattr(self, attr_name, fn(self))
        return getattr(self, attr_name)

    return _lazy_property


class Daily():
    def __init__(self, tmax, tmin, ea, rs, uz, zw, elev, lat, doy, method='asce', rso_type=None, rso=None):
        """ASCE Daily Standardized Reference Evapotranspiration (ET)

        Arguments
        ---------
        tmax : ee.Image
            Maximum daily temperature [C].
        tmin : ee.Image
            Minimum daily temperature [C].
        ea : ee.Image
            Actual vapor pressure [kPa].
        rs : ee.Image
            Incoming shortwave solar radiation [MJ m-2 day-1].
        uz : ee.Image
            Wind speed [m s-1].
        zw : ee.Number
            Wind speed height [m].
        elev : ee.Image or ee.Number
            Elevation [m].
        lat : ee.Image or ee.Number
            Latitude [degrees].
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

        Notes
        -----
        Latitude units are degrees, not radians.

        References
        ----------
        .. [1] ASCE-EWRI (2005). The ASCE standardized reference evapotranspiration
            equation. ASCE-EWRI Standardization of Reference Evapotranspiration
            Task Committee Rep., ASCE Reston, Va.

        """

        if method.lower() not in ['asce', 'refet']:
            raise ValueError('method must be "asce" or "refet"')

        if rso_type is None:
            pass
        elif rso_type.lower() not in ['simple', 'full', 'array']:
            raise ValueError('rso_type must be None, "simple", "full", or "array')
        elif rso_type.lower() in 'array':
            # Check that rso is an ee.Image or ee.Number?
            pass

        # Get time_start from tmin
        # Should time_start be set in init?
        self.time_start = ee.Image(tmin).get('system:time_start')
        self.date = ee.Date(self.time_start)

        # Do these all need to be set onto self?
        self.tmax = tmax
        self.tmin = tmin
        self.ea = ea
        self.rs = rs
        self.uz = uz
        self.zw = zw
        self.elev = elev
        self.lat = lat
        self.doy = doy

        # Convert latitude to radians
        self.lat = self.lat.multiply(math.pi / 180)

        # To match standardized form, pair is calculated from elevation
        self.pair = calcs.air_pressure(self.elev, method)

        # Psychrometric constant (Eq. 4)
        self.psy = self.pair.multiply(0.000665)

        self.tmean = self.tmax.add(self.tmin).multiply(0.5)
        self.es_slope = calcs.es_slope(self.tmean, method)

        # Saturated vapor pressure
        self.es = calcs.sat_vapor_pressure(self.tmax).add(calcs.sat_vapor_pressure(self.tmin)).multiply(0.5)

        # Vapor pressure deficit
        self.vpd = calcs.vpd(es=self.es, ea=self.ea)

        # Extraterrestrial radiation
        self.ra = calcs.ra_daily(lat=self.lat, doy=self.doy, method=method)

        # Clear sky solar radiation
        # If rso_type is not set, use the method
        # If rso_type is set, use rso_type directly
        if rso_type is None:
            if method.lower() == 'asce':
                self.rso = calcs.rso_simple(ra=self.ra, elev=self.elev)
            elif method.lower() == 'refet':
                self.rso = calcs.rso_daily(ea=self.ea, pair=self.pair, ra=self.ra, doy=self.doy, lat=self.lat)
        elif rso_type.lower() == 'simple':
            self.rso = calcs.rso_simple(ra=self.ra, elev=self.elev)
        elif rso_type.lower() == 'full':
            self.rso = calcs.rso_daily(ea=self.ea, pair=self.pair, ra=self.ra, doy=self.doy, lat=self.lat)
        elif rso_type.lower() == 'array':
            # Use rso array passed to function
            self.rso = rso

        # Cloudiness fraction
        self.fcd = calcs.fcd_daily(rs=self.rs, rso=self.rso)

        # Net long-wave radiation
        self.rnl = calcs.rnl_daily(tmax=self.tmax, tmin=self.tmin, ea=self.ea, fcd=self.fcd)

        # Net radiation
        self.rn = calcs.rn(self.rs, self.rnl)

        # Wind speed
        self.u2 = calcs.wind_height_adjust(uz=self.uz, zw=self.zw)

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
            return self.etr
        elif surface.lower() in ['grass', 'eto', 'short']:
            return self.eto
        else:
            raise ValueError('unsupported surface type: {}'.format(surface))

    @lazy_property
    def eto(self):
        """Short (grass) reference surface"""
        self.cn = 900
        self.cd = 0.34
        return ee.Image(self._etsz().rename(['eto']).set('system:time_start', self.time_start))

    @lazy_property
    def etr(self):
        """Tall (alfalfa) reference surface"""
        self.cn = 1600
        self.cd = 0.38
        return ee.Image(self._etsz().rename(['etr']).set('system:time_start', self.time_start))

    def _etsz(self):
        """Daily reference ET (Eq. 1)

        Returns
        -------
        etsz : ee.Image
            Standardized reference ET [mm].

        """
        return (
            self.u2.multiply(self.vpd).multiply(self.cn).multiply(self.psy).divide(self.tmean.add(273))
            .add(self.es_slope.multiply(self.rn).multiply(0.408))
            .divide(self.u2.multiply(self.cd).add(1).multiply(self.psy).add(self.es_slope))
        )
        # # CGM - Intentionally leaving expression form of calculation here for reference
        # return self.tmin.expression(
        #     '((0.408 * es_slope * rn) + (psy * cn * u2 * vpd / (tmean + 273))) / '
        #     '(es_slope + psy * (cd * u2 + 1))',
        #     {
        #         'cd': self.cd, 'cn': self.cn, 'es_slope': self.es_slope,
        #         'psy': self.psy, 'rn': self.rn, 'tmean': self.tmean,
        #         'u2': self.u2, 'vpd': self.vpd,
        #     }
        #  )

    @lazy_property
    def etw(self):
        """Priestley-Taylor evaporation (alpha = 1.26)

        Returns
        -------
        etw : ee.Image
            Priestley-Taylor ET [mm].

        References
        ----------
        https://wetlandscapes.github.io/blog/blog/penman-monteith-and-priestley-taylor/

        """
        return (
            self.es_slope
            .expression(
                '(alpha * es_slope * rn * 1000 / (2453 * (es_slope + psy)))',
                {'es_slope': self.es_slope, 'rn': self.rn, 'psy': self.psy, 'alpha': 1.26}
            )
            .rename(['etw'])
            .set('system:time_start', self.time_start)
        )

    @lazy_property
    def eto_fs1(self):
        """UF-IFAS Extension FS1 Radiation Term (ETrad)

        Returns
        -------
        eto_fs1 : ee.Image
            FS1 ETrad [mm].

        References
        ----------
        https://edis.ifas.ufl.edu/pdffiles/ae/ae45900.pdf

        """
        return (
            self.u2
            .expression(
                '(delta / (delta + psy * (1 + 0.34 * u2))) * (0.408 * rn)',
                {'delta': self.es_slope, 'psy': self.psy, 'u2': self.u2, 'rn': self.rn}
            )
            .rename(['eto_fs1'])
            .set('system:time_start', self.time_start)
        )

    @lazy_property
    def eto_fs2(self):
        """UF-IFAS Extension FS2 Wind Term (ETwind)

        Returns
        -------
        eto_fs2 : ee.Image
            FS2 ETwind [mm].

        References
        ----------
        https://edis.ifas.ufl.edu/pdffiles/ae/ae45900.pdf

        """
        # Temperature Term (Eq. 14)
        tt = self.u2.expression('(900 / (t + 273)) * u2', {'t': self.tmean, 'u2': self.u2})

        # Psi Term (Eq. 13)
        pt = self.u2.expression(
            'psy / (slope + psy * (1 + 0.34 * u2))',
            {'slope': self.es_slope, 'psy': self.psy, 'u2': self.u2}
        )

        return (
            self.es.subtract(self.ea).multiply(pt).multiply(tt)
            .rename(['eto_fs2'])
            .set('system:time_start', self.time_start)
        )

    @lazy_property
    def pet_hargreaves(self):
        """Hargreaves potential ET

        Returns
        -------
        hargreaves_pet : ee.Image
            Hargreaves ET [mm].

        References
        ----------

        """
        return (
            self.tmax
            .expression(
                '0.0023 * (tmean + 17.8) * ((tmax - tmin) ** 0.5) * 0.408 * ra',
                {'tmean': self.tmean, 'tmax': self.tmax, 'tmin': self.tmin, 'ra': self.ra}
            )
            .rename(['pet_hargreaves'])
            .set('system:time_start', self.time_start)
        )

    @classmethod
    def gridmet(cls, input_img, zw=None, elev=None, lat=None, method='asce', rso_type=None):
        """Initialize daily RefET from a GRIDMET image

        Parameters
        ----------
        input_img : ee.Image
            GRIDMET image from the collection IDAHO_EPSCOR/GRIDMET.
        zw : ee.Number, optional
            Wind speed height [m] (the default is 10).
        elev : ee.Image or ee.Number, optional
            Elevation image [m].  The standard GRIDMET elevation image
            (projects/openet/assets/meteorology/gridmet/elevation) will be used if not set.
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
        image_date = ee.Date(input_img.get('system:time_start'))

        if zw is None:
            zw = ee.Number(10)
        if elev is None:
            elev = ee.Image('projects/openet/assets/meteorology/gridmet/elevation')
        if lat is None:
            lat = ee.Image('projects/openet/assets/meteorology/gridmet/latitude')

        return cls(
            tmax=input_img.select(['tmmx']).subtract(273.15),
            tmin=input_img.select(['tmmn']).subtract(273.15),
            ea=calcs.actual_vapor_pressure(
                pair=calcs.air_pressure(elev, method),
                q=input_img.select(['sph']),
            ),
            rs=input_img.select(['srad']).multiply(0.0864),
            uz=input_img.select(['vs']),
            zw=zw,
            elev=elev,
            lat=lat,
            doy=ee.Number(image_date.getRelative('day', 'year')).add(1).double(),
            method=method,
            rso_type=rso_type,
        )

    @classmethod
    def maca(cls, input_img, zw=None, elev=None, lat=None, method='asce', rso_type=None):
        """Initialize daily RefET from a MACA image

        Parameters
        ----------
        input_img : ee.Image
            MACA image from the collection IDAHO_EPSCOR/MACAv2_METDATA.
        zw : ee.Number, optional
            Wind speed height [m] (the default is 10).
        elev : ee.Image or ee.Number, optional
            Elevation image [m].  The standard GRIDMET elevation image
            (projects/openet/assets/meteorology/gridmet/elevation) will be used if not set.
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
        Actual vapor pressure is computed from specific humidity (MACA huss)
            and air pressure (from elevation).
        Wind speed from east and north vector components m s-1.

        """
        image_date = ee.Date(input_img.get('system:time_start'))

        # TODO: CHECK MACA transform (Why the difference?)
        # maca_transform = [0.04163593073593073, 0, -124.7722,
        #                   0, 0.04159470085470086, 25.0631]
        # gridmet_transform = [0.041666666666666664, 0, -124.78749996666667,
        #                      0, -0.041666666666666664, 49.42083333333334]

        if zw is None:
            zw = ee.Number(10)
        if elev is None:
            # TODO: should we use gridmet elev? maca grid is shifted
            elev = ee.Image('projects/openet/assets/meteorology/gridmet/elevation')
        if lat is None:
            elev = ee.Image('projects/openet/assets/meteorology/gridmet/latitude')

        def wind_magnitude(input_img):
            """Compute daily wind magnitude from vectors"""
            return (
                ee.Image(input_img.select(['uas'])).pow(2)
                .add(ee.Image(input_img.select(['vas'])).pow(2))
                .sqrt().rename(['uz'])
            )

        return cls(
            tmax=input_img.select(['tasmax']).subtract(273.15),
            tmin=input_img.select(['tasmin']).subtract(273.15),
            ea=calcs.actual_vapor_pressure(
                pair=calcs.air_pressure(elev, method),
                q=input_img.select(['huss']),
            ),
            rs=input_img.select(['rsds']).multiply(0.0864),
            uz=ee.Image(wind_magnitude(input_img)).select(['uz']),
            zw=zw,
            elev=elev,
            lat=lat,
            doy=ee.Number(image_date.getRelative('day', 'year')).add(1).double(),
            method=method,
            rso_type=rso_type,
        )

    @classmethod
    def nldas(cls, input_coll, zw=None, elev=None, lat=None, method='asce', rso_type=None):
        """Initialize daily RefET from a hourly NLDAS image collection

        Parameters
        ----------
        input_coll : ee.ImageCollection
            Collection of NLDAS hourly images for a single day from the
            collection NASA/NLDAS/FORA0125_H002.
        zw : ee.Number, optional
            Wind speed height [m] (the default is 10).
        elev : ee.Image or ee.Number, optional
            Elevation image [m].  A custom NLDAS elevation image
            (projects/eddi-noaa/nldas/elevation) will be used if not set.
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
        Solar radiation is converted from W m-2 to MJ m-2 day-1.
        Actual vapor pressure is computed from specific humidity and air
            pressure (from elevation).

        """
        input_coll = ee.ImageCollection(input_coll)
        image_date = ee.Date(ee.Image(input_coll.first()).get('system:time_start'))

        if zw is None:
            zw = ee.Number(10)
        if elev is None:
            elev = ee.Image('projects/openet/assets/meteorology/nldas/ancillary/elevation')
        if lat is None:
            lat = ee.Image('projects/openet/assets/meteorology/nldas/ancillary/latitude')

        def wind_magnitude(input_img):
            """Compute hourly wind magnitude from vectors"""
            return (
                ee.Image(input_img.select(['wind_u'])).pow(2)
                .add(ee.Image(input_img.select(['wind_v'])).pow(2))
                .sqrt().rename(['uz'])
            )

        return cls(
            tmax=input_coll.select(['temperature']).max(),
            tmin=input_coll.select(['temperature']).min(),
            ea=calcs.actual_vapor_pressure(
                pair=calcs.air_pressure(elev, method),
                q=input_coll.select(['specific_humidity']).mean(),
            ),
            rs=input_coll.select(['shortwave_radiation']).sum().multiply(0.0036),
            uz=ee.Image(ee.ImageCollection(input_coll.map(wind_magnitude)).mean()),
            zw=zw,
            elev=elev,
            lat=lat,
            doy=ee.Number(image_date.getRelative('day', 'year')).add(1).double(),
            method=method,
            rso_type=rso_type,
        )

    @classmethod
    def cfsv2(cls, input_coll, zw=None, elev=None, lat=None, method='asce', rso_type=None):
        """Initialize daily RefET from a 6-hourly CFSv2 image collection

        Parameters
        ----------
        input_coll : ee.ImageCollection
            Collection of CFSv2 6 hourly images for a single day from the
            collection NOAA/CFSV2/FOR6H.
        zw : ee.Number, optional
            Wind speed height [m] (the default is 10).
        elev : ee.Image or ee.Number, optional
            Elevation image [m].  The SRTM elevation image (CGIAR/SRTM90_V4)
            will be reprojected to the CFSv2 grid if not set.
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
        Actual vapor pressure is computed from specific humidity and air
            pressure (from elevation).

        """
        input_coll = ee.ImageCollection(input_coll)
        image_date = ee.Date(ee.Image(input_coll.first()).get('system:time_start'))

        cfsv2_crs = 'EPSG:4326'
        # Transform for 2011 to Present
        cfsv2_transform = [0.20454520376789903, 0., -180, 0., -0.20442210122586923, 90]
        # Transform for 1979 to 2010
        # cfsv2_transform = [0.3124995757764987, 0, -180, 0, -0.31221217943274454, 90]

        if zw is None:
            zw = ee.Number(10)
        if elev is None:
            # TODO: Build a CFSv2 elevation asset
            elev = ee.Image('CGIAR/SRTM90_V4').reproject(cfsv2_crs, cfsv2_transform)
        if lat is None:
            lat = ee.Image.pixelLonLat().select('latitude').reproject(cfsv2_crs, cfsv2_transform)

        def wind_magnitude(input_img):
            """Compute hourly wind magnitude from vectors"""
            u_img = ee.Image(input_img).select(['u-component_of_wind_height_above_ground'])
            v_img = ee.Image(input_img).select(['v-component_of_wind_height_above_ground'])
            return u_img.pow(2).add(v_img.pow(2)).sqrt()

        return cls(
            tmax=input_coll.select(['Maximum_temperature_height_above_ground_6_Hour_Interval'])
                .max().subtract(273.15),
            tmin=input_coll.select(['Minimum_temperature_height_above_ground_6_Hour_Interval'])
                .min().subtract(273.15),
            ea=calcs.actual_vapor_pressure(
                pair=calcs.air_pressure(elev, method),
                q=input_coll.select(['Specific_humidity_height_above_ground']).mean()
            ),
            # TODO: Check the conversion on solar
            rs=input_coll.select(['Downward_Short-Wave_Radiation_Flux_surface_6_Hour_Average'])
                .mean().multiply(0.0864),
            uz=ee.Image(ee.ImageCollection(input_coll.map(wind_magnitude)).mean()),
            zw=zw,
            elev=elev,
            lat=lat,
            doy=ee.Number(image_date.getRelative('day', 'year')).add(1).double(),
            method=method,
            rso_type=rso_type,
        )

    @classmethod
    def rtma(cls, input_coll, rs=None, zw=None, elev=None, lat=None, method='asce', rso_type=None):
        """Initialize daily RefET from a hourly RTMA image collection

        Parameters
        ----------
        input_coll : ee.ImageCollection
            Collection of RTMA hourly images for a single day from the
            collection NOAA/NWS/RTMA.
        rs : ee.Image, str, optional
            Incoming solar radiation [MJ m-2 day-1].  The GRIDMET image for the
            concurrent day will be used if not set.
        zw : ee.Number, optional
            Wind speed height [m] (the default is 10).
        elev : ee.Image or ee.Number, optional
            Elevation image [m].  The RTMA elevation image asset will be used if not set
            (projects/openet/assets/meteorology/rtma/ancillary/elevation).
        lat : ee.Image or ee.Number
            Latitude image [degrees].  The RTMA latitude image asset will be used if not set
            (projects/openet/assets/meteorology/rtma/ancillary/latitude).
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
        Solar radiation is converted from W m-2 to MJ m-2 day-1.
        Actual vapor pressure is computed from specific humidity and air
            pressure (from elevation).

        """
        input_coll = ee.ImageCollection(input_coll)
        start_date = ee.Date(ee.Image(input_coll.first()).get('system:time_start'))

        # Parse the solar radiation input
        if isinstance(rs, ee.Image):
            pass
        elif isinstance(rs, ee.Number) or isinstance(rs, float) or isinstance(rs, int):
            rs = ee.Image.constant(rs)
        elif (rs is None) or (rs.upper() == 'GRIDMET'):
            rs_coll = (
                ee.ImageCollection('IDAHO_EPSCOR/GRIDMET')
                .filterDate(start_date, start_date.advance(1, 'day'))
                .select(['srad'])
            )
            rs = ee.Image(rs_coll.first()).multiply(0.0864)
        elif rs.upper() == 'NLDAS':
            rs_coll = (
                ee.ImageCollection('NASA/NLDAS/FORA0125_H002')
                .filterDate(start_date, start_date.advance(1, 'day'))
                .select(['shortwave_radiation'])
            )
            # TODO: Check this unit conversion
            rs = ee.Image(rs_coll.sum()).multiply(0.0036)
        else:
            raise ValueError('Unsupported Rs input')

        if zw is None:
            zw = ee.Number(10)
        if elev is None:
            elev = ee.Image('projects/openet/assets/meteorology/rtma/ancillary/elevation')
        if lat is None:
            lat = ee.Image('projects/openet/assets/meteorology/rtma/ancillary/latitude')

        # # Not computing wind speed from the components since wind speed is available
        # def wind_magnitude(input_img):
        #     """Compute hourly wind magnitude from vectors"""
        #     return (
        #         ee.Image(input_img.select(['UGRD'])).pow(2)
        #         .add(ee.Image(input_img.select(['VGRD'])).pow(2))
        #         .sqrt().rename(['WIND'])
        #     )

        return cls(
            tmax=input_coll.select(['TMP']).max(),
            tmin=input_coll.select(['TMP']).min(),
            ea=calcs.actual_vapor_pressure(
                pair=calcs.air_pressure(elev, method),
                q=input_coll.select(['SPFH']).mean()
            ),
            rs=rs,
            # Using wind speed directly instead of computing from components
            uz=input_coll.select(['WIND']).mean(),
            # uz=ee.Image(ee.ImageCollection(input_coll.map(wind_magnitude)).mean())
            zw=zw,
            elev=elev,
            lat=lat,
            doy=ee.Number(start_date.getRelative('day', 'year')).add(1).double(),
            method=method,
            rso_type=rso_type,
        )

    @classmethod
    def era5(cls, input_coll, zw=None, elev=None, lat=None, method='asce', rso_type=None):
        """Initialize daily RefET from an hourly ERA5 image collection

        Parameters
        ----------
        input_coll : ee.ImageCollection
            Collection of ERA5 hourly images for a single day from the
            collection ECMWF/ERA5/HOURLY.
        zw : ee.Number, optional
            Wind speed height [m] (the default is 10).
        elev : ee.Image or ee.Number, optional
            Elevation image [m].  The OpenET ERA5 elevation image
            (projects/openet/assets/meteorology/era5/ancillary/elevation)
            will be used if not set.
        lat : ee.Image or ee.Number
            Latitude image [degrees].  The OpenET ERA5 latitude image
            (projects/openet/assets/meteorology/era5/ancillary/latitude)
            will be used if not set.
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
        Solar radiation is summed and converted from J m-2 to MJ m-2 day-1.
        Actual vapor pressure is computed from dew point temperature.

        """
        input_coll = ee.ImageCollection(input_coll)
        start_date = ee.Date(ee.Image(input_coll.first()).get('system:time_start'))

        if zw is None:
            zw = ee.Number(10)
        if elev is None:
            elev = ee.Image('projects/openet/assets/meteorology/era5/ancillary/elevation')
        if lat is None:
            lat = ee.Image('projects/openet/assets/meteorology/era5/ancillary/latitude')

        def wind_magnitude(input_img):
            """Compute hourly wind magnitude from vectors"""
            return (
                ee.Image(input_img.select(['u_component_of_wind_10m'])).pow(2)
                .add(ee.Image(input_img.select(['v_component_of_wind_10m'])).pow(2))
                .sqrt().rename(['wind_10m'])
            )

        return cls(
            tmax=input_coll.select(['temperature_2m']).max().subtract(273.15),
            tmin=input_coll.select(['temperature_2m']).min().subtract(273.15),
            ea=calcs.sat_vapor_pressure(input_coll.select(['dewpoint_temperature_2m']).mean().subtract(273.15)),
            rs=input_coll.select(['surface_solar_radiation_downwards']).sum().divide(1000000),
            uz=ee.Image(ee.ImageCollection(input_coll.map(wind_magnitude)).mean()),
            zw=zw,
            elev=elev,
            lat=lat,
            doy=ee.Number(start_date.getRelative('day', 'year')).add(1).double(),
            method=method,
            rso_type=rso_type,
        )

    @classmethod
    def era5_land(
            cls,
            input_coll,
            zw=None,
            elev=None,
            lat=None,
            method='asce',
            rso_type=None,
            fill_edge_cells=False,
    ):
        """Initialize daily RefET from an hourly ERA5-Land image collection

        Parameters
        ----------
        input_coll : ee.ImageCollection
            Collection of ERA5-Land hourly images for a single day from the
            collection ECMWF/ERA5_LAND/HOURLY.
        zw : ee.Number, optional
            Wind speed height [m] (the default is 10).
        elev : ee.Image or ee.Number, optional
            Elevation image [m].  The OpenET ERA5-Land elevation image
            (projects/openet/assets/meteorology/era5land/ancillary/elevation)
            will be used if not set.
        lat : ee.Image or ee.Number
            Latitude image [degrees].  The OpenET ERA5-Land latitude image
            (projects/openet/assets/meteorology/era5land/ancillary/latitude)
            will be used if not set.
        method : {'asce' (default), 'refet'}, optional
            Specifies which calculation method to use.
            * 'asce' -- Calculations will follow ASCE-EWRI 2005.
            * 'refet' -- Calculations will follow RefET software.
        rso_type : {None (default), 'full' , 'simple'}, optional
            Specifies which clear sky solar radiation (Rso) model to use.
            * None -- Rso type will be determined from "method" parameter
            * 'full' -- Full clear sky solar formulation
            * 'simple' -- Simplified clear sky solar formulation
        fill_edge_cells : bool, optional
            If True, fill any masked pixels within a one cell buffer of the input image.
            This will fill in small holes and cells along the coasts.

        Notes
        -----
        Temperatures are converted from K to C.
        Solar radiation is summed and converted from J m-2 to MJ m-2 day-1.
        Actual vapor pressure is computed from dew point temperature.

        """
        input_coll = ee.ImageCollection(input_coll)
        start_date = ee.Date(ee.Image(input_coll.first()).get('system:time_start'))

        if zw is None:
            zw = ee.Number(10)
        if elev is None:
            elev = ee.Image('projects/openet/assets/meteorology/era5land/ancillary/elevation')
        if lat is None:
            lat = ee.Image('projects/openet/assets/meteorology/era5land/ancillary/latitude')

        def wind_magnitude(input_img):
            """Compute hourly wind magnitude from vectors"""
            return (
                ee.Image(input_img.select(['u_component_of_wind_10m'])).pow(2)
                .add(ee.Image(input_img.select(['v_component_of_wind_10m'])).pow(2))
                .sqrt().rename(['wind_10m'])
            )

        # Compute the input variables and convert to standard units
        tmax = input_coll.select(['temperature_2m']).max().subtract(273.15)
        tmin = input_coll.select(['temperature_2m']).min().subtract(273.15)
        ea = calcs.sat_vapor_pressure(input_coll.select(['dewpoint_temperature_2m']).mean().subtract(273.15))
        rs = input_coll.select(['surface_solar_radiation_downwards_hourly']).sum().divide(1000000)
        uz = ee.Image(ee.ImageCollection(input_coll.map(wind_magnitude)).mean())

        # Fill any masked pixels along the edge of the land mask
        # This will fill most coastal pixels and any small holes
        if fill_edge_cells:
            def fill(image):
                img = ee.Image(image)
                return img.unmask(
                    img.reduceNeighborhood('mean', ee.Kernel.square(1), 'kernel', False)
                    .reproject(img.projection())
                )
            tmax = fill(tmax)
            tmin = fill(tmin)
            ea = fill(ea)
            rs = fill(rs)
            uz = fill(uz)

        return cls(
            tmax=tmax,
            tmin=tmin,
            ea=ea,
            rs=rs,
            uz=uz,
            zw=zw,
            elev=elev,
            lat=lat,
            doy=ee.Number(start_date.getRelative('day', 'year')).add(1).double(),
            method=method,
            rso_type=rso_type,
        )

    # @classmethod
    # def prism(cls, input_img, lat=None):
    #     """Initialize daily RefET from a PRISM image
    #
    #     Parameters
    #     ----------
    #     input_img : ee.Image
    #         PPRISM image from the collection OREGONSTATE/PRISM/AN81d.
    #     lat : ee.Image or ee.Number
    #         Latitude image [degrees].  The latitude will be computed
    #         dynamically using ee.Image.pixelLonLat() if not set.
    #
    #     Notes
    #     -----
    #     PRISM can only be used to compute Hargreaves PET.
    #
    #     """
    #     image_date = ee.Date(input_img.get('system:time_start'))
    #
    #     if lat is None:
    #         lat = ee.Image('projects/earthengine-legacy/assets/'
    #                        'projects/climate-engine/prism/elevation')\
    #             .multiply(0).add(ee.Image.pixelLonLat().select('latitude'))\
    #             .rename(['latitude'])
    #         # prism_transform = [0.0416666666667, 0, -125.02083333333336,
    #         #                    0, -0.0416666666667, 49.93749999999975]
    #         # lat = ee.Image.pixelLonLat().select('latitude')\
    #         #     .reproject('EPSG:4326', prism_transform)
    #         # lat = input_img.select([0]).multiply(0)\
    #         #     .add(ee.Image.pixelLonLat().select('latitude'))
    #
    #     return cls(
    #         tmax=input_img.select(['tmmx']),
    #         tmin=input_img.select(['tmmn']),
    #         ea=0,
    #         rs=0,
    #         uz=0,
    #         zw=2,
    #         elev=0,
    #         lat=lat,
    #         doy=ee.Number(image_date.getRelative('day', 'year')).add(1).double(),
    #         method='asce',
    #     )
