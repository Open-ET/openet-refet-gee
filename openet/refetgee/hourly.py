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


class Hourly():
    def __init__(self, tmean, ea, rs, uz, zw, elev, lat, lon, doy, time, method='asce'):
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
        Latitude & longitude units are degrees, not radians.

        References
        ----------
        .. [1] ASCE-EWRI (2005). The ASCE standardized reference evapotranspiration
            equation. ASCE-EWRI Standardization of Reference Evapotranspiration
            Task Committee Rep., ASCE Reston, Va.

        """

        if method.lower() not in ['asce', 'refet']:
            raise ValueError('method must be "asce" or "refet"')

        # Get time_start from tmin
        # Should time_start be set in init?
        self.time_start = ee.Image(tmean).get('system:time_start')
        self.date = ee.Date(self.time_start)

        # Do these all need to be set onto self?
        self.tmean = tmean
        self.ea = ea
        self.rs = rs
        self.uz = uz
        self.zw = zw
        self.elev = elev
        self.lat = lat
        self.lon = lon
        self.doy = doy
        self.time = time

        # Convert latitude & longitude to radians
        self.lat = self.lat.multiply(math.pi / 180)
        self.lon = self.lon.multiply(math.pi / 180)

        # To match standardized form, psy is calculated from elevation based pair
        self.pair = calcs.air_pressure(self.elev, method=method)

        # Psychrometric constant (Eq. 35)
        self.psy = self.pair.multiply(0.000665)

        self.es = calcs.sat_vapor_pressure(self.tmean)
        self.es_slope = calcs.es_slope(self.tmean, method)

        # Vapor pressure deficit
        self.vpd = self.es.subtract(self.ea)
        # self.vpd = calcs.vpd(es=self.es, ea=self.ea)

        # Extraterrestrial radiation
        time_mid = self.time.add(0.5)
        self.ra = calcs.ra_hourly(lat=self.lat, lon=self.lon, doy=self.doy, time_mid=time_mid, method=method)

        # Clear sky solar radiation
        if method == 'asce':
            self.rso = calcs.rso_simple(ra=self.ra, elev=self.elev)
        elif method == 'refet':
            self.rso = calcs.rso_hourly(
                ea=self.ea, ra=self.ra, pair=self.pair, doy=self.doy,
                time_mid=time_mid, lat=self.lat, lon=self.lon, method=method,
            )

        # Cloudiness fraction
        # Intentionally not using time_mid to match Beta value in IN2 file
        # In IN2, "Beta" is computed for the start of the time period,
        #   but "SinBeta" is computed for the midpoint.
        # Beta (not SinBeta) is used for clamping fcd.
        self.fcd = calcs.fcd_hourly(
            rs=self.rs, rso=self.rso, doy=self.doy, time_mid=self.time,
            lat=self.lat, lon=self.lon, method=method,
        )

        # Net long-wave radiation
        self.rnl = calcs.rnl_hourly(tmean=self.tmean, ea=self.ea, fcd=self.fcd)

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
        self.cn = ee.Number(37.0)
        cd_day = 0.24
        g_rn_day = 0.1
        cd_night = 0.96
        g_rn_night = 0.5

        # Adjust coefficients for daytime/nighttime
        # Nighttime is defined as when Rn < 0 (pg 44)
        self.cd = self.rn.multiply(0).add(cd_day).where(self.rn.lt(0), cd_night)
        self.g_rn = self.rn.multiply(0).add(g_rn_day).where(self.rn.lt(0), g_rn_night)

        # Soil heat flux (Eqs. 65 and 66)
        self.g = self.rn.multiply(self.g_rn)

        return ee.Image(self._etsz().rename(['eto']).set('system:time_start', self.time_start))

    @lazy_property
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
        self.g_rn = self.rn.multiply(0).add(g_rn_day).where(self.rn.lt(0), g_rn_night)

        # Soil heat flux (Eqs. 65 and 66)
        self.g = self.rn.multiply(self.g_rn)

        return ee.Image(self._etsz().rename(['etr']).set('system:time_start', self.time_start))

    def _etsz(self):
        """Hourly reference ET (Eq. 1)

        Returns
        -------
        etsz : ee.Image
            Standardized reference ET [mm].

        """
        return (
            self.u2.multiply(self.vpd).multiply(self.cn).multiply(self.psy).divide(self.tmean.add(273))
            .add(self.es_slope.multiply(self.rn.subtract(self.g)).multiply(0.408))
            .divide(self.u2.multiply(self.cd).add(1).multiply(self.psy).add(self.es_slope))
        )
        # # CGM - Intentionally leaving expression form of calculation here for reference
        # return self.tmean.expression(
        #     '((0.408 * es_slope * (rn - g)) + (psy * cn * u2 * vpd / (tmean + 273))) / '
        #     '(es_slope + psy * (cd * u2 + 1))',
        #     {
        #         'cd': self.cd, 'cn': self.cn, 'es_slope': self.es_slope,
        #         'g': self.g, 'psy': self.psy, 'rn': self.rn, 'tmean': self.tmean,
        #         'u2': self.u2, 'vpd': self.vpd,
        #     }
        # )

    @classmethod
    def nldas(cls, input_img, zw=None, elev=None, lat=None, lon=None, method='asce'):
        """Initialize hourly RefET from an NLDAS image

        Parameters
        ----------
        input_img : ee.Image
            NLDAS hourly image from the collection NASA/NLDAS/FORA0125_H002.
        zw : ee.Number, optional
            Wind speed height [m] (the default is 10).
        elev : ee.Image or ee.Number, optional
            Elevation image [m].  A precomputed NLDAS elevation image
            (projects/openet/assets/meteorology/nldas/ancillary/elevation)
            will be used if not set.
        lat : ee.Image or ee.Number
            Latitude image [degrees].  A precomputed NLDAS latitude image
            (projects/openet/assets/meteorology/nldas/ancillary/latitude)
            will be used if not set.
        lon : ee.Image or ee.Number
            Longitude image [degrees].  A precomputed NLDAS longitude image
            (projects/openet/assets/meteorology/nldas/ancillary/longitude)
            will be used if not set.
        method : {'asce' (default), 'refet'}, optional
            Specifies which calculation method to use.
            * 'asce' -- Calculations will follow ASCE-EWRI 2005.
            * 'refet' -- Calculations will follow RefET software.

        Notes
        -----
        Solar radiation is converted from W m-2 to MJ m-2 hour-1.
        Actual vapor pressure is computed from specific humidity and air
            pressure (from elevation).

        """
        image_date = ee.Date(input_img.get('system:time_start'))

        if zw is None:
            zw = ee.Number(10)
        if elev is None:
            elev = ee.Image('projects/openet/assets/meteorology/nldas/ancillary/elevation')
        if lat is None:
            lat = ee.Image('projects/openet/assets/meteorology/nldas/ancillary/latitude')
        if lon is None:
            lon = ee.Image('projects/openet/assets/meteorology/nldas/ancillary/longitude')

        return cls(
            tmean=input_img.select(['temperature']),
            ea=calcs.actual_vapor_pressure(
                pair=calcs.air_pressure(elev, method),
                q=input_img.select(['specific_humidity']),
            ),
            rs=input_img.select(['shortwave_radiation']).multiply(0.0036),
            uz=input_img.select(['wind_u']).pow(2).add(input_img.select(['wind_v']).pow(2)).sqrt().rename(['uz']),
            zw=zw,
            elev=elev,
            lat=lat,
            lon=lon,
            doy=ee.Number(image_date.getRelative('day', 'year')).add(1).double(),
            time=ee.Number(image_date.get('hour')),
            method=method,
        )

    @classmethod
    def rtma(cls, input_img, rs=None, zw=None, elev=None, lat=None, lon=None, method='asce'):
        """Initialize hourly RefET from an RTMA image

        Parameters
        ----------
        input_img : ee.Image
            RTMA hourly image from the collection NOAA/NWS/RTMA.
        rs : ee.Image, ee.Number, str, optional
            Incoming solar radiation [MJ m-2 hr-1].
            The NLDAS2 image for the concurrent hour will be used if not set.
        zw : ee.Number, optional
            Wind speed height [m] (the default is 10).
        elev : ee.Image or ee.Number, optional
            Elevation image [m].  The RTMA elevation image asset will be used if not set
            (projects/openet/assets/meteorology/rtma/ancillary/elevation).
        lat : ee.Image or ee.Number
            Latitude image [degrees].  The RTMA latitude image asset will be used if not set
            (projects/openet/assets/meteorology/rtma/ancillary/latitude).
        lon : ee.Image or ee.Number
            Longitude image [degrees].  The RTMA longitude image asset will be used if not set
            (projects/openet/assets/meteorology/rtma/ancillary/longitude).
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
        start_date = ee.Date(input_img.get('system:time_start'))

        # Parse the solar radiation input
        if isinstance(rs, ee.Image):
            pass
        elif isinstance(rs, ee.Number) or isinstance(rs, float) or isinstance(rs, int):
            rs = ee.Image.constant(rs)
        elif (rs is None) or ((rs is str) and (rs.upper() in ['NLDAS', 'NLDAS2'])):
            rs = (
                ee.ImageCollection('NASA/NLDAS/FORA0125_H002')
                .filterDate(start_date, start_date.advance(30, 'minute'))
                .select(['shortwave_radiation'])
            )
            rs = ee.Image(rs.first()).multiply(0.0036)
        else:
            raise ValueError('Unsupported Rs input')

        if zw is None:
            zw = ee.Number(10)
        if elev is None:
            elev = ee.Image('projects/openet/assets/meteorology/rtma/ancillary/elevation')
        if lat is None:
            lat = ee.Image('projects/openet/assets/meteorology/rtma/ancillary/latitude')
        if lon is None:
            lon = ee.Image('projects/openet/assets/meteorology/rtma/ancillary/longitude')

        return cls(
            tmean=input_img.select(['TMP']),
            # Computing vapor pressure from specific humidity instead of dew point
            ea=calcs.actual_vapor_pressure(pair=calcs.air_pressure(elev, method), q=input_img.select(['SPFH'])),
            # ea=calcs.sat_vapor_pressure(input_img.select(['DPT'])),
            rs=rs,
            # Using wind speed band directly instead of computing from components
            uz=input_img.select(['WIND']),
            # uz=input_img.select(['UGRD']).pow(2).add(input_img.select(['VGRD']).pow(2)).sqrt(),
            zw=zw,
            elev=elev,
            lat=lat,
            lon=lon,
            doy=ee.Number(start_date.getRelative('day', 'year')).add(1).double(),
            time=ee.Number(start_date.get('hour')),
            method=method,
        )

    @classmethod
    def era5(cls, input_img, zw=None, elev=None, lat=None, lon=None, method='asce'):
        """Initialize hourly RefET from an ERA5 image

        Parameters
        ----------
        input_img : ee.Image
            ERA5 hourly image from the collection ECMWF/ERA5/HOURLY.
        zw : ee.Number, optional
            Wind speed height [m] (the default is 10).
        elev : ee.Image or ee.Number, optional
            Elevation image [m].  The OpenET ERA5-Land elevation image
            (projects/openet/assets/meteorology/era5/ancillary/elevation)
            will be used if not set.
        lat : ee.Image or ee.Number
            Latitude image [degrees].  The OpenET ERA5 latitude image
            (projects/openet/assets/meteorology/era5/ancillary/latitude)
            will be used if not set.
        lon : ee.Image or ee.Number
            Longitude image [degrees].  The OpenET ERA5 longitude image
            (projects/openet/assets/meteorology/era5/ancillary/longitude)
            will be used if not set.
        method : {'asce' (default), 'refet'}, optional
            Specifies which calculation method to use.
            * 'asce' -- Calculations will follow ASCE-EWRI 2005.
            * 'refet' -- Calculations will follow RefET software.

        Notes
        -----
        Temperatures are converted from K to C.
        Solar radiation is summed and converted from J m-2 to MJ m-2 hr-1.
        Actual vapor pressure is computed from dew point temperature.

        """
        start_date = ee.Date(input_img.get('system:time_start'))

        if zw is None:
            zw = ee.Number(10)
        if elev is None:
            elev = ee.Image('projects/openet/assets/meteorology/era5/ancillary/elevation')
        if lat is None:
            lat = ee.Image('projects/openet/assets/meteorology/era5/ancillary/latitude')
        if lon is None:
            lon = ee.Image('projects/openet/assets/meteorology/era5/ancillary/longitude')

        return cls(
            tmean=input_img.select(['temperature_2m']).subtract(273.15),
            ea=calcs.sat_vapor_pressure(input_img.select(['dewpoint_temperature_2m']).subtract(273.15)),
            rs=input_img.select(['surface_solar_radiation_downwards']).divide(1000000),
            uz=input_img.select(['u_component_of_wind_10m']).pow(2)
                .add(input_img.select(['v_component_of_wind_10m']).pow(2))
                .sqrt().rename(['wind_10m']),
            zw=zw,
            elev=elev,
            lat=lat,
            lon=lon,
            doy=ee.Number(start_date.getRelative('day', 'year')).add(1).double(),
            time=ee.Number(start_date.get('hour')),
            method=method,
        )

    @classmethod
    def era5_land(
            cls,
            input_img,
            zw=None,
            elev=None,
            lat=None,
            lon=None,
            method='asce',
            fill_edge_cells=False,
    ):
        """Initialize hourly RefET from an ERA5-Land image

        Parameters
        ----------
        input_img : ee.Image
            ERA5-Land hourly image from the collection ECMWF/ERA5_LAND/HOURLY.
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
        lon : ee.Image or ee.Number
            Longitude image [degrees].  The OpenET ERA5-Land longitude image
            (projects/openet/assets/meteorology/era5land/ancillary/longitude)
            will be used if not set.
        method : {'asce' (default), 'refet'}, optional
            Specifies which calculation method to use.
            * 'asce' -- Calculations will follow ASCE-EWRI 2005.
            * 'refet' -- Calculations will follow RefET software.
        fill_edge_cells : bool, optional
            If True, fill any masked pixels within a one cell buffer of the input image.
            This will fill in small holes and cells along the coasts.

        Notes
        -----
        Temperatures are converted from K to C.
        Solar radiation is summed and converted from J m-2 to MJ m-2 hr-1.
        Actual vapor pressure is computed from dew point temperature.

        """
        start_date = ee.Date(input_img.get('system:time_start'))

        if zw is None:
            zw = ee.Number(10)
        if elev is None:
            elev = ee.Image('projects/openet/assets/meteorology/era5land/ancillary/elevation')
        if lat is None:
            lat = ee.Image('projects/openet/assets/meteorology/era5land/ancillary/latitude')
        if lon is None:
            lon = ee.Image('projects/openet/assets/meteorology/era5land/ancillary/longitude')

        # Compute the input variables and convert to standard units
        tmean = input_img.select(['temperature_2m']).subtract(273.15)
        ea = calcs.sat_vapor_pressure(input_img.select(['dewpoint_temperature_2m']).subtract(273.15))
        rs = input_img.select(['surface_solar_radiation_downwards_hourly']).divide(1000000)
        uz = (
            input_img.select(['u_component_of_wind_10m']).pow(2)
            .add(input_img.select(['v_component_of_wind_10m']).pow(2))
            .sqrt().rename(['wind_10m'])
        )

        # Fill any masked pixels along the edge of the land mask
        # This will fill most coastal pixels and any small holes
        if fill_edge_cells:
            def fill(image):
                img = ee.Image(image)
                return img.unmask(
                    img.reduceNeighborhood('mean', ee.Kernel.square(1), 'kernel', False)
                    .reproject(img.projection())
                )
            tmean = fill(tmean)
            ea = fill(ea)
            rs = fill(rs)
            uz = fill(uz)

        return cls(
            tmean=tmean,
            ea=ea,
            rs=rs,
            uz=uz,
            zw=zw,
            elev=elev,
            lat=lat,
            lon=lon,
            doy=ee.Number(start_date.getRelative('day', 'year')).add(1).double(),
            time=ee.Number(start_date.get('hour')),
            method=method,
        )
