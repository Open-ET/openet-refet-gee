# Validation

Validation of the RefET functions was done using [AgriMet](https://www.usbr.gov/pn/agrimet/) data for all of 2015 from the [Fallon, NV station](https://www.usbr.gov/pn/agrimet/agrimetmap/falnda.html).  The data were downloaded and prepared using the setup_test_data.py Python script in in the tests/data folder.

## Daily

Daily AgriMet data query for 2015:
```
https://www.usbr.gov/pn-bin/daily.pl?station=FALN&year=2015&month=1&day=1&year=2015&month=12&day=31&pcode=MN&pcode=MX&pcode=SR&pcode=YM&pcode=UA
```

Daily AgriMet parameters and units:
* MN = Minimum Daily Air Temperature (F)
* MX = Maximum Daily Air Temperature (F)
* SR = Daily Global Solar Radiation (langleys)
* YM = Mean Daily Dewpoint Temperature (F)
* UA = Daily Average Wind Speed (mph)

Using Tdew instead of Rh for validation (Table 3, page 12) since Ea is not available at the daily time step.

## Hourly

Hourly Agrimet data query for 2015:
```
https://www.usbr.gov/pn-bin/instant.pl?station=FALN&year=2015&month=1&day=1&year=2015&month=12&day=31&pcode=OB&pcode=EA&pcode=WS&pcode=SI&print_hourly=1
```

Hourly AgriMet parameters and units:
* OB = Air Temperature, 15 Minute Instantaneous (degrees F)
* TP = Dew Point Temperature, 15 Minute Average (degrees F)
* WS = Wind Speed, Hourly Average (mph)
* SI = 15 minute Solar Radiation (langleys/hour)

Using Tdew for the validation, but the preferred method in the documentation (Table 4, page 30) is Ea, Tdew, then Rh.
* EA = Actual Vapor Pressure, 15 Minute Average (KPa)
* TU = Relative Humidity, 15 Minute Average (percent)

## Ref-ET

Ref-ET v4.2 (v4.1.9.8.2017) was used to generate the validation dataset.
