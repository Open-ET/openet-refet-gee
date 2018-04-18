.. image:: https://travis-ci.org/cgmorton/RefET-GEE.svg?branch=master
   :target: https://travis-ci.org/cgmorton/RefET-GEE

.. image:: https://badge.fury.io/py/RefET-GEE.svg
   :target: https://badge.fury.io/py/RefET-GEE

=======================================================================
Google Earth Engine ASCE Standardized Reference Evapotranspiration (ET)
=======================================================================

Google Earth Engine (GEE) functions for computing daily and hourly reference ET.

Usage
=====

Daily
-----

The following demonstrates how to compute a single daily ETr value using weather data for 2015-07-01 from the `Fallon, NV AgriMet station <https://www.usbr.gov/pn/agrimet/agrimetmap/falnda.html>`__.
The necessary unit conversions are shown on the input values.
The raw input data is available `here <https://www.usbr.gov/pn-bin/daily.pl?station=FALN&year=2015&month=7&day=1&year=2015&month=7&day=1&pcode=ETRS&pcode=MN&pcode=MX&pcode=SR&pcode=YM&pcode=UA>`__.

.. code-block:: console

    import math
    import ee
    import eerefet

    # Unit conversions
    tmin_c = (66.65 - 32) * (5.0 / 9)                          # F -> C
    tmax_c = (102.80 - 32) * (5.0 / 9)                         # F -> C
    tdew_c = (57.26 - 32) * (5.0 / 9)                          # F -> C
    ea = 0.6108 * math.exp(17.27 * tdew_c / (tdew_c + 237.3))  # kPa
    rs = (674.07 * 0.041868)                                   # Langleys -> MJ m-2 d-1
    uz = 4.80 * 0.44704                                        # mpg -> m s-1
    lat_radians = (39.4575 * math.pi / 180)                    # degrees -> radians

    etr = eerefet.Daily(
        tmin=tmin_c, tmax=tmax_c, ea=ea, rs=rs, uz=uz, zw=3, elev=1208.5,
        lat=lat_radians, doy=182).etr().getInfo()

    print('ETr: {:.2f} mm'.format(float(etr)))

Hourly
------

The following demonstrates how to compute a single hourly ETr value using weather data for 18:00 UTC (11:00 AM PDT) on 2015-07-01 from the `Fallon, NV AgriMet station <https://www.usbr.gov/pn/agrimet/agrimetmap/falnda.html>`__.
The necessary unit conversions are shown on the input values.
The raw input data is available `here <https://www.usbr.gov/pn-bin/instant.pl?station=FALN&year=2015&month=7&day=1&year=2015&month=7&day=1&pcode=OB&pcode=EA&pcode=WS&pcode=SI&print_hourly=1>`__

.. code-block:: console

    import math
    import ee
    import eerefet

    # Unit conversions
    tmean_c = (91.80 - 32) * (5.0 / 9)           # F -> C
    ea = 1.20                                    # kPa
    rs = (61.16 * 0.041868)                      # Langleys -> MJ m-2 h-1
    uz = 3.33 * 0.44704                          # mph -> m s-1
    lat_radians = (39.4575 * math.pi / 180)      # degrees -> radians
    lon_radians = (-118.77388 * math.pi / 180)   # degrees -> radians

    etr = eerefet.Hourly(
        tmean=tmean_c, ea=ea, rs=rs, uz=uz, zw=3, elev=1208.5,
        lat=lat_radians, lon=lon_radians, doy=182, time=18).etr().getInfo()

    print('ETr: {:.2f} mm'.format(float(etr)))

GRIDMET
-------

.. code-block:: console
    import ee
    import eerefet

    gridmet_img = ee.Image(ee.ImageCollection('IDAHO_EPSCOR/GRIDMET').first())
    etr = eerefet.Daily.gridmet(gridmet_img).etr().getInfo()

    print('ETr: {:.2f} mm'.format(float(etr)))

Input Parameters
================

Required Parameters (hourly & daily)
-----------------------------------

========  ===================  =================================================
Variable  Type                 Description [units]
========  ===================  =================================================
ea        ee.Image, ee.Number  Actual vapor pressure [kPa]
rs        ee.Image, ee.Number  Incoming shortwave solar radiation [MJ m-2 day-1]
uz        ee.Image, ee.Number  Wind speed [m/s]
zw        ee.Number              Wind speed height [m]
elev      ee.Image, ee.Number  Elevation [m]
lat       ee.Image, ee.Number  Latitude [radians]
doy       ee.Image, ee.Number  Day of year
========  ===================  =================================================

Required Daily Parameters
-------------------------

========  ===================  =================================================
Variable  Type                 Description [units]
========  ===================  =================================================
tmin      ee.Image, ee.Number  Minimum daily temperature [C]
tmax      ee.Image, ee.Number  Maximum daily temperature [C]
========  ===================  =================================================

Required Hourly Parameters
--------------------------

========  ===================  =================================================
Variable  Type                 Description [units]
========  ===================  =================================================
tmean     ee.Image, ee.Number  Average hourly temperature [C]
lon       ee.Image, ee.Number  Longitude [radians]
time      ee.Image, ee.Number  UTC hour at start of time period
========  ===================  =================================================

Optional Parameters
-------------------

========  =========  ====================================================
Variable  Type       Description [units]
========  =========  ====================================================
method    str        | Calculation method
                       * 'refet' -- Calculations will follow RefET software (default)
                       * 'asce' -- Calculations will follow ASCE-EWRI 2005 equations
rso_type  str        | Clear sky solar radiation (Rso) model
                       * 'full' -- Full clear sky solar formulation (default)
                       * 'simple' -- Simplified clear sky solar formulation (Eq. 19)
                       * 'array' -- Read Rso values from "rso" function parameter
rso       ee.Image   | Clear sky solar radiation [MJ m-2 day-1]
          ee.Number    * Only needed if rso_type is 'array'
                       * Defaults to None if not set
========  =========  ====================================================


Limitations
===========

Currently the user must handle all of the file I/O and unit conversions.

Cloudiness Fraction (hourly)
----------------------------

The hourly reference ET calculation is currently performed independently for each time step.  The cloudiness fraction (fcd) for very low sun angles (i.e. at night) is hard coded to 1 for very low sun angles instead of being derived from the .

Installation
============

To install the RefET-GEE python module:

.. code-block:: console

   pip install eerefet

Validation
==========

Please see the `validation document <VALIDATION.md>`__ for additional details on the source of the test values and the comparison of the functions to the Ref-ET software.

Dependencies
============

 * `earthengine-api <https://github.com/google/earthengine-api>`__

Modules needed to run the test suite:

 * `pandas <http://pandas.pydata.org>`__
 * `pytest <https://docs.pytest.org/en/latest/>`__
 * `pytz <http://pythonhosted.org/pytz/>`__

References
==========

ASCE-EWRI Standardized Reference Evapotranspiration Equation (2005)

 * `Report <http://www.kimberly.uidaho.edu/water/asceewri/ascestzdetmain2005.pdf>`__
 * `Appendix <http://www.kimberly.uidaho.edu/water/asceewri/appendix.pdf>`__