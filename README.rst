=======================================================================
Google Earth Engine ASCE Standardized Reference Evapotranspiration (ET)
=======================================================================

|version| |build|

Google Earth Engine (GEE) functions for computing daily and hourly reference ET following the ASCE Standardized Reference Evapotranspiration Equations (ASCE2005_).

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
    import geerefet

    # Unit conversions
    tmin_c = (66.65 - 32) * (5.0 / 9)                          # F -> C
    tmax_c = (102.80 - 32) * (5.0 / 9)                         # F -> C
    tdew_c = (57.26 - 32) * (5.0 / 9)                          # F -> C
    ea = 0.6108 * math.exp(17.27 * tdew_c / (tdew_c + 237.3))  # kPa
    rs = (674.07 * 0.041868)                                   # Langleys -> MJ m-2 d-1
    uz = 4.80 * 0.44704                                        # mpg -> m s-1
    lat = 39.4575                                              # degrees

    etr = geerefet.Daily(
        tmin=tmin_c, tmax=tmax_c, ea=ea, rs=rs, uz=uz, zw=3, elev=1208.5,
        lat=lat, doy=182).etr().getInfo()

    print('ETr: {:.2f} mm'.format(float(etr)))

Hourly
------

The following demonstrates how to compute a single hourly ETr value using weather data for 18:00 UTC (11:00 AM PDT) on 2015-07-01 from the `Fallon, NV AgriMet station <https://www.usbr.gov/pn/agrimet/agrimetmap/falnda.html>`__.
The necessary unit conversions are shown on the input values.
The raw input data is available `here <https://www.usbr.gov/pn-bin/instant.pl?station=FALN&year=2015&month=7&day=1&year=2015&month=7&day=1&pcode=OB&pcode=EA&pcode=WS&pcode=SI&print_hourly=1>`__

.. code-block:: console

    import math
    import ee
    import geerefet

    # Unit conversions
    tmean_c = (91.80 - 32) * (5.0 / 9)           # F -> C
    ea = 1.20                                    # kPa
    rs = (61.16 * 0.041868)                      # Langleys -> MJ m-2 h-1
    uz = 3.33 * 0.44704                          # mph -> m s-1
    lat = 39.4575                                # degrees
    lon = -118.77388                             # degrees

    etr = geerefet.Hourly(
        tmean=tmean_c, ea=ea, rs=rs, uz=uz, zw=3, elev=1208.5,
        lat=lat, lon=lon_radians, doy=182, time=18).etr().getInfo()

    print('ETr: {:.2f} mm'.format(float(etr)))

GRIDMET
-------

A helper function for computing daily ETo and ETr for `GRIDMET <http://www.climatologylab.org/gridmet.html>`__ images is available.

.. code-block:: console

    import ee
    import geerefet

    gridmet_img = ee.Image(ee.ImageCollection('IDAHO_EPSCOR/GRIDMET').first())
    etr = geerefet.Daily.gridmet(gridmet_img).etr().getInfo()

    print('ETr: {:.2f} mm'.format(float(etr)))

NLDAS
-----

Helper functions for computing daily/hourly ETo/ETr for `NLDAS <https://ldas.gsfc.nasa.gov/nldas/NLDAS2forcing.php>`__ images are available.

For the daily function, the NLDAS collection must be filtered to a single 24 hour period.

.. code-block:: console

    import ee
    import geerefet

    nldas_coll = ee.ImageCollection('NASA/NLDAS/FORA0125_H002')\
        .filterDate('2017-07-01', '2017-07-02)
    etr = geerefet.Daily.nldas(nldas_coll).etr().getInfo()

    print('ETr: {:.2f} mm'.format(float(etr)))

.. code-block:: console

    import ee
    import geerefet

    nldas_img = ee.Image(ee.ImageCollection('NASA/NLDAS/FORA0125_H002').first())
    etr = geerefet.Hourly.nldas(nldas_img).etr().getInfo()

    print('ETr: {:.2f} mm'.format(float(etr)))

Input Parameters
================

Required Parameters (hourly & daily)
------------------------------------

========  ===================  =================================================
Variable  Type                 Description [units]
========  ===================  =================================================
ea        ee.Image             Actual vapor pressure [kPa]
rs        ee.Image             Incoming shortwave solar radiation [MJ m-2 day-1]
uz        ee.Image             Wind speed [m s-1]
zw        ee.Number            Wind speed height [m]
elev      ee.Image, ee.Number  Elevation [m]
lat       ee.Image, ee.Number  Latitude [degrees]
doy       ee.Image, ee.Number  Day of year
========  ===================  =================================================

Required Daily Parameters
-------------------------

========  ===================  =================================================
Variable  Type                 Description [units]
========  ===================  =================================================
tmin      ee.Image             Minimum daily temperature [C]
tmax      ee.Image             Maximum daily temperature [C]
========  ===================  =================================================

Required Hourly Parameters
--------------------------

========  ===================  =================================================
Variable  Type                 Description [units]
========  ===================  =================================================
tmean     ee.Image             Average hourly temperature [C]
lon       ee.Image, ee.Number  Longitude [degrees]
time      ee.Number            UTC hour at start of time period
========  ===================  =================================================

Optional Parameters
-------------------

========  ===================  ====================================================
Variable  Type                 Description [units]
========  ===================  ====================================================
method    str                  | Calculation method

                               * 'asce' -- Calculations will follow ASCE-EWRI 2005 (default)
                               * 'refet' -- Calculations will follow RefET software

rso_type  str                  | Override default clear sky solar radiation (Rso) calculation
                               | Defaults to None if not set

                               * 'full' -- Full clear sky solar formulation (default)
                               * 'simple' -- Simplified clear sky solar formulation (Eq. 19)
                               * 'array' -- Read Rso values from "rso" function parameter

rso       ee.Image, ee.Number  | Clear sky solar radiation [MJ m-2 day-1]

                               * Only needed if rso_type is 'array'
                               * Defaults to None if not set

========  ===================  ====================================================

Issues
======

Currently the user must handle all of the file I/O and unit conversions.

Cloudiness Fraction (hourly)
----------------------------

The cloudiness fraction (fcd) is computed as the ratio of the measured solar radiation (Rs) to the theoretical clear sky solar radiation (Rso).  This ratio cannot be computed directly at night since Rso is 0.  ASCE-EWRI 2005 suggests computing a representative nighttime fcd based on the fcd at sunset and/or sunrise.

In the RefET module fcd is hard coded to 1 for all time steps with very low sun angles since the hourly reference ET is computed independently for each time step.

Calculation Method - ASCE vs. RefET
===================================

The main difference between the two "methods" is that the "asce" method attempts to follow the equations in ASCE2005_, whereas the "refet" method attempts to follow the calculations of the `RefET Software <https://www.uidaho.edu/cals/kimberly-research-and-extension-center/research/water-resources/ref-et-software>`__ as closely as possible.  The difference in output between these methods is generally negligible (if not identical for realistic numbers of significant digits).  Note that the default is set to "asce" to best match the calculations a user would expect to have happen. The "refet" method was added in order to help validate this code to the RefET Software.

Installation
============

To install the RefET-GEE python module:

.. code-block:: console

    pip install geerefet

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

.. [ASCE2005]
 | ASCE-EWRI (2005). The ASCE standardized reference evapotranspiration equation.
 | `https://ascelibrary.org/doi/book/10.1061/9780784408056 <https://ascelibrary.org/doi/book/10.1061/9780784408056>`__

.. |build| image:: https://github.com/cgmorton/RefET-GEE/workflows/build/badge.svg
   :alt: Build status
   :target: https://github.com/cgmorton/RefET-GEE
.. |version| image:: https://badge.fury.io/py/geerefet.svg
   :alt: Latest version on PyPI
   :target: https://badge.fury.io/py/geerefet
