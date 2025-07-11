================================================================================
OpenET - Google Earth Engine ASCE Standardized Reference Evapotranspiration (ET)
================================================================================

|version| |build|

This repository provides `Google Earth Engine <https://earthengine.google.com/>`__ Python API based implementation of the ASCE Standardized Reference Evapotranspiration Equations (ASCE2005_) for computing daily and hourly reference ET.

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
    import openet.refetgee

    # Unit conversions
    tmin_c = (66.65 - 32) * (5.0 / 9)                          # F -> C
    tmax_c = (102.80 - 32) * (5.0 / 9)                         # F -> C
    tdew_c = (57.26 - 32) * (5.0 / 9)                          # F -> C
    ea = 0.6108 * math.exp(17.27 * tdew_c / (tdew_c + 237.3))  # kPa
    rs = (674.07 * 0.041868)                                   # Langleys -> MJ m-2 d-1
    uz = 4.80 * 0.44704                                        # mpg -> m s-1
    lat = 39.4575                                              # degrees

    etr = openet.refetgee.Daily(
        tmin=tmin_c, tmax=tmax_c, ea=ea, rs=rs, uz=uz, zw=3, elev=1208.5,
        lat=lat, doy=182).etr.getInfo()

    print('ETr: {:.2f} mm'.format(float(etr)))

Hourly
------

The following demonstrates how to compute a single hourly ETr value using weather data for 18:00 UTC (11:00 AM PDT) on 2015-07-01 from the `Fallon, NV AgriMet station <https://www.usbr.gov/pn/agrimet/agrimetmap/falnda.html>`__.
The necessary unit conversions are shown on the input values.
The raw input data is available `here <https://www.usbr.gov/pn-bin/instant.pl?station=FALN&year=2015&month=7&day=1&year=2015&month=7&day=1&pcode=OB&pcode=EA&pcode=WS&pcode=SI&print_hourly=1>`__

.. code-block:: console

    import math
    import ee
    import openet.refetgee

    # Unit conversions
    tmean_c = (91.80 - 32) * (5.0 / 9)           # F -> C
    ea = 1.20                                    # kPa
    rs = (61.16 * 0.041868)                      # Langleys -> MJ m-2 h-1
    uz = 3.33 * 0.44704                          # mph -> m s-1
    lat = 39.4575                                # degrees
    lon = -118.77388                             # degrees

    etr = openet.refetgee.Hourly(
        tmean=tmean_c, ea=ea, rs=rs, uz=uz, zw=3, elev=1208.5,
        lat=lat, lon=lon, doy=182, time=18
    ).etr.getInfo()

    print('ETr: {:.2f} mm'.format(float(etr)))

GRIDMET
-------

A helper function for computing daily ETo and ETr for `GRIDMET <http://www.climatologylab.org/gridmet.html>`__ images is available.

.. code-block:: console

    import ee
    import openet.refetgee

    source_img = ee.Image(ee.ImageCollection('IDAHO_EPSCOR/GRIDMET').first())
    etr = (
        openet.refetgee.Daily.gridmet(source_img).etr
        .reduceRegion(reducer=ee.Reducer.first(),
                      geometry=ee.Geometry.Point(-118.77388, 39.4575),
                      scale=1000)
        .getInfo()
    )

    print('ETr: {:.2f} mm'.format(float(etr['etr'])))

NLDAS
-----

Helper functions for computing daily/hourly ETo/ETr for `NLDAS <https://ldas.gsfc.nasa.gov/nldas/NLDAS2forcing.php>`__ images are available.

For the daily function, the NLDAS collection must be filtered to a single 24 hour period.

.. code-block:: console

    import ee
    import openet.refetgee

    source_coll = ee.ImageCollection('NASA/NLDAS/FORA0125_H002')\
        .filterDate('2015-07-01', '2015-07-02')
    etr = (
        openet.refetgee.Daily.nldas(source_coll).etr
        .reduceRegion(reducer=ee.Reducer.first(),
                      geometry=ee.Geometry.Point(-118.77388, 39.4575),
                      scale=1000)
        .getInfo()
    )

    print('ETr: {:.2f} mm'.format(float(etr['etr'])))

.. code-block:: console

    import ee
    import openet.refetgee

    source_img = ee.Image('NASA/NLDAS/FORA0125_H002/A20150701_2000')
    etr = (
        openet.refetgee.Hourly.nldas(source_img).etr
        .reduceRegion(reducer=ee.Reducer.first(),
                      geometry=ee.Geometry.Point(-118.77388, 39.4575),
                      scale=1000)
        .getInfo()
    )

    print('ETr: {:.2f} mm'.format(float(etr['etr'])))

CFSv2
-----

A helper function for computing daily ETo and ETr for `CFSv2 <http://>`__ images is available.

For the daily function, the CFSv2 collection must be filtered to a single 24 hour period.

.. code-block:: console

    import ee
    import openet.refetgee

    source_coll = ee.ImageCollection('NOAA/CFSV2/FOR6H').filterDate('2015-07-01', '2015-07-02')
    etr = (
        openet.refetgee.Daily.cfsv2(source_coll).etr
        .reduceRegion(reducer=ee.Reducer.first(),
                      geometry=ee.Geometry.Point(-118.77388, 39.4575),
                      scale=1000)
        .getInfo()
    )

    print('ETr: {:.2f} mm'.format(float(etr['etr'])))

RTMA
-----

Helper functions for computing daily/hourly ETo/ETr for `RTMA <https://>`__ images are available.

For the daily function, the RTMA collection must be filtered to a single 24 hour period.

.. code-block:: console

    import ee
    import openet.refetgee

    source_coll = ee.ImageCollection('NOAA/NWS/RTMA').filterDate('2015-07-01', '2015-07-02')
    etr = (
        openet.refetgee.Daily.rtma(source_coll).etr
        .reduceRegion(reducer=ee.Reducer.first(),
                      geometry=ee.Geometry.Point(-118.77388, 39.4575),
                      scale=1000)
        .getInfo()
    )

    print('ETr: {:.2f} mm'.format(float(etr['etr'])))

.. code-block:: console

    import ee
    import openet.refetgee

    source_img = ee.Image('NOAA/NWS/RTMA/2015070120')
    etr = (
        openet.refetgee.Hourly.nldas(source_img).etr
        .reduceRegion(reducer=ee.Reducer.first(),
                      geometry=ee.Geometry.Point(-118.77388, 39.4575),
                      scale=1000)
        .getInfo()
    )

    print('ETr: {:.2f} mm'.format(float(etr['etr'])))

ERA5-Land
---------

Helper functions for computing daily/hourly ETo/ETr for `ERA5-Land <https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land>`__ images are available.

For the daily function, the ERA5-Land collection must be filtered to a single 24 hour period.

.. code-block:: console

    import ee
    import openet.refetgee

    source_coll = ee.ImageCollection('ECMWF/ERA5_LAND/HOURLY')\
        .filterDate('2015-07-01', '2015-07-02')
    etr = (
        openet.refetgee.Daily.era5_land(source_coll).etr
        .reduceRegion(reducer=ee.Reducer.first(),
                      geometry=ee.Geometry.Point(-118.77388, 39.4575),
                      scale=1000)
        .getInfo()
    )

    print('ETr: {:.2f} mm'.format(float(etr['etr'])))

.. code-block:: console

    import ee
    import openet.refetgee

    source_img = ee.Image('ECMWF/ERA5_LAND/HOURLY/20150701T20')
    etr = (
        openet.refetgee.Hourly.era5_land(source_img).etr
        .reduceRegion(reducer=ee.Reducer.first(),
                      geometry=ee.Geometry.Point(-118.77388, 39.4575),
                      scale=1000)
        .getInfo()
    )

    print('ETr: {:.2f} mm'.format(float(etr['etr'])))

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

The OpenET RefET GEE python module can be installed via pip:

.. code-block:: console

    pip install openet-refet-gee

OpenET Namespace Package
========================

Each OpenET model is stored in the "openet" folder (namespace).  The model can then be imported as a "dot" submodule of the main openet module.

.. code-block:: console

    import openet.refetgee as refetgee

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

.. |build| image:: https://github.com/Open-ET/openet-refet-gee/workflows/build/badge.svg
   :alt: Build status
   :target: https://github.com/Open-ET/openet-refet-gee
.. |version| image:: https://badge.fury.io/py/openet-refet-gee.svg
   :alt: Latest version on PyPI
   :target: https://badge.fury.io/py/openet-refet-gee
