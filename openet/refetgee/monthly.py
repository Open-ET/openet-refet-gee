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


class Monthly():
    """"""

    def __init__(self, tmean, lat, tmean_monthly):
        """Monthly Potential Evapotranspiration (ET)

        Arguments
        ---------
        tmean : ee.Image
            Mean daily temperature [C].
        lat : ee.Image or ee.Number
            Latitude [degrees].
        tmean_monthly : ee.Image
            Mean daily temperature [C] for all 12 months in the year.

        Raises
        ------

        Notes
        -----
        Latitude units are degress, not radians.

        References
        ----------

        """

        # Get time_start from tmin
        # Should time_start be set in init?
        self.time_start = ee.Image(tmean).get('system:time_start')
        self.date = ee.Date(self.time_start)

        # Do these all need to be set onto self?
        self.tmean = tmean
        self.lat = lat

        # Convert latitude to radians
        self.lat = self.lat.multiply(math.pi / 180)

        # CGM - I'm not sure the best place to set these values
        #   We could change the method from a lazy_property and pass them in?
        self.tmean_monthly = tmean_monthly

    @lazy_property
    def pet_thornthwaite(self):
        """Thornthwaite potential ET

        Returns
        -------
        thornthwaite_pet : ee.Image
            Thornthwaite ET [mm].

        References
        ----------
        Thornthwaite, C. W. (1948). "An approach toward a rational classification
        of climate". Geographical Review. 38 (1): 55–94. doi:10.2307/210739

        """

        # Number of days in month
        # The image start date is probably already the month start date,
        #   but build it explicitly just to be sure
        month_start = ee.Date.fromYMD(self.date.get("year"), self.date.get("month"), 1)
        month_end = month_start.advance(1, "month")
        days_in_month = month_end.difference(month_start, "day")

        # TODO: Check if there is a day length equation in calcs, it not add it
        def day_length_func(doy):
            """
            Compute the day length in hours for the target day of year (Eq. 34, FAO 56)
            :param img: Earth Engine Image
            :return: Earth Engine Image
            """
            omegas = calcs._omega_sunset(self.lat, calcs._delta(ee.Number(doy)))
            return omegas.multiply(24.0 / math.pi)

        # Compute the average day length (hours) for the month
        doy_start = ee.Number(month_start.getRelative("day", "year")).add(1).double()
        doy_list = ee.List.sequence(doy_start, doy_start.add(days_in_month).subtract(1))
        day_length = ee.ImageCollection(doy_list.map(day_length_func)).mean()
        # For testing, use the day length for the midpoint of the month
        # test_doy = ee.Number(start_date.getRelative("day", "year")).add(1).double().add(14)
        # day_length = day_length_func(test_doy)

        # Average daily air temperature
        # Clamp temperature to >= 0
        ta = self.tmean.select("tmean").max(0)

        # CGM - Check if divide and power are being applied separately to each band
        heat_index = self.tmean_monthly.divide(5).pow(1.514).reduce(ee.Reducer.sum())

        a = heat_index.expression(
            "(6.75e-07 * hi ** 3) - (7.71e-05 * hi ** 2) + (1.792e-02 * hi) + 0.49239",
            {"hi": heat_index})
        print(a.getInfo())

        return (
            self.tmean.expression(
                "1.6 * (L / 12) * (N / 30) * ((10 * Ta  / I) ** a)",
                {"L": day_length, "N": days_in_month, "Ta": ta, "I": heat_index, "a": a})
            .rename(["pet_thornthwaite"])
            .set({"system:time_start": self.time_start})
        )

    @classmethod
    def prism(cls, input_img, lat=None):
        """Initialize monthly RefET from a PRISM image

        Parameters
        ----------
        input_img : ee.Image
            PPRISM image from the collection OREGONSTATE/PRISM/AN81m.
        lat : ee.Image or ee.Number
            Latitude image [degrees].  The latitude will be computed
            dynamically using ee.Image.pixelLonLat() if not set.

        Notes
        -----
        PRISM monthly can only be used to compute Thornthwaite PET.

        """
        image_date = ee.Date(input_img.get('system:time_start'))

        if lat is None:
            lat = ee.Image('projects/earthengine-legacy/assets/'
                           'projects/climate-engine/prism/elevation')\
                .multiply(0).add(ee.Image.pixelLonLat().select('latitude'))\
                .rename(['latitude'])
            # prism_transform = [0.0416666666667, 0, -125.02083333333336,
            #                    0, -0.0416666666667, 49.93749999999975]
            # lat = ee.Image.pixelLonLat().select('latitude')\
            #     .reproject('EPSG:4326', prism_transform)
            # lat = input_img.select([0]).multiply(0)\
            #     .add(ee.Image.pixelLonLat().select('latitude'))

        # CGM - I'm not sure where to set the monthly mean image
        # Are the monthly means always for the same year as the target image?
        start_date = ee.Date.fromYMD(image_date.get("year"), 1, 1)
        tmean_monthly = (
            ee.ImageCollection("OREGONSTATE/PRISM/AN81m")
            .filterDate(start_date, start_date.advance(1, "year"))
            .select(["tmean"])
            .toBands()
        )

        return cls(
            tmean=input_img.select(['tmean']),
            lat=lat,
            tmean_monthly=tmean_monthly,
        )



