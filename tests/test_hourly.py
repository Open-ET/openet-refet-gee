import datetime as dt
import math
import os

import ee
import numpy as np
import pandas as pd
import pytest
import pytz

from geerefet.hourly import Hourly
import geerefet.units as units

ee.Initialize()


# Eventually move to conftest.py or a separate file
# Fallon AgriMet site parameters
s_args = {
    'elev': 1208.5,
    'lat': 39.4575,
    'lon': -118.77388,
    # DEADBEEF
    # 'lat': units._deg2rad(39.4575),
    # 'lon': units._deg2rad(-118.77388),
    'zw': 3.0,
}

# Hourly test parameters for 2015-07-01 18:00 UTC (11:00 AM PDT)
h_args = {
    'doy': 182,
    'ea': 1.1990099614301906,
    'es': 5.09318785259078,
    'eto': 0.6065255163817055,
    'etr': 0.7201865213918281,
    'etr_asce': 0.7196369609713682,
    'ra': 4.30824147948541,
    'rnl': 0.22897874401150786,
    'rs': 61.16 * 0.041868,  # Conversion from Langleys to MJ m-2
    'tdew': units._f2c(49.36),
    'time': 18.0,
    'time_mid': 18.5,
    'tmean': units._f2c(91.80),
    'uz': 3.33 * 0.44704,  # Conversion from mph to m s-1
    'u2': 1.3709275319197722,
}

# # Playing around with passing a dictionary of function keyword arguments
# hourly_args = {
#     k: v for k, v in hourly.items()
#     if k in ['tmean', 'ea', 'rs', 'uz', 'zw', 'doy', 'time']}
# hourly_args.update({'surface':'etr'})


# Test full hourly functions with positional inputs
def test_refet_hourly_input_positions():
    etr = Hourly(
        ee.Number(h_args['tmean']), ee.Number(h_args['ea']),
        ee.Number(h_args['rs']), ee.Number(h_args['uz']),
        ee.Number(s_args['zw']), ee.Number(s_args['elev']),
        ee.Number(s_args['lat']), ee.Number(s_args['lon']),
        ee.Number(h_args['doy']), ee.Number(h_args['time']),
        method='refet').etr().getInfo()
    assert float(etr) == pytest.approx(h_args['etr'])


# Test full hourly calculations with keyword inputs
# Test surface, rso_type, and rso inputs
def test_refet_hourly_asce():
    etr = Hourly(
        tmean=ee.Number(h_args['tmean']), ea=ee.Number(h_args['ea']),
        rs=ee.Number(h_args['rs']), uz=ee.Number(h_args['uz']),
        zw=ee.Number(s_args['zw']), elev=ee.Number(s_args['elev']),
        lat=ee.Number(s_args['lat']), lon=ee.Number(s_args['lon']),
        doy=ee.Number(h_args['doy']), time=ee.Number(h_args['time']),
        method='asce').etr().getInfo()
    assert float(etr) == pytest.approx(h_args['etr_asce'])


# Test hourly functions using actual RefET input/output files
# DEADBEEF - This doesn't work if I move it to conftest.py
class HourlyData():
    """Setup hourly validation data from Fallon AgriMet station"""
    val_ws = os.path.join(os.getcwd(), 'tests', 'data')
    # val_ws = os.path.join(os.path.dirname(os.getcwd()), 'tests', 'data')

    csv_path = os.path.join(val_ws, 'FALN_Agrimet_hourly_raw_2015.csv')
    # in2_path = os.path.join(val_ws, 'FALN_Agrimet_hourly_raw_2015.in2')
    out_path = os.path.join(val_ws, 'FALN_Agrimet_hourly_raw_2015.out')

    # Read in the inputs CSV file using pandas
    csv_df = pd.read_csv(csv_path, engine='python', na_values='NO RECORD')
    csv_df.rename(
        columns={'OB': 'TEMP', 'TP': 'TDEW', 'WS': 'WIND', 'SI': 'RS'},
        inplace=True)
    # AgriMet times are local with DST (this will drop one hour)
    # DEADBEEF - Can't set timezone as variable in class?
    csv_df['DATETIME'] = csv_df[['YEAR', 'MONTH', 'DAY', 'HOUR']].apply(
        lambda x: pytz.timezone('US/Pacific').localize(dt.datetime(*x)), axis=1)
    # To match RefET IN2 values exactly, compute DOY using localtime (not UTC)
    csv_df['DOY'] = csv_df['DATETIME'].apply(lambda x: int(x.strftime('%j')))
    csv_df['DATETIME'] = csv_df['DATETIME'].apply(
        lambda x: x.tz_convert('UTC').strftime('%Y-%m-%d %H:00'))
    csv_df.set_index('DATETIME', inplace=True, drop=True)

    # Convert inputs units
    csv_df['TEMP'] = units._f2c(csv_df['TEMP'])
    csv_df['TDEW'] = units._f2c(csv_df['TDEW'])
    csv_df['WIND'] *= 0.44704
    csv_df['RS'] *= 0.041868  # Conversion from Langleys to MJ m-2 to match RefET
    # csv_df['RS'] *= 0.041840  # Alternate conversion from Langleys to MJ m-2
    csv_df['EA'] = 0.6108 * np.exp(17.27 * csv_df['TDEW'] / (csv_df['TDEW'] + 237.3))

    # Eventually compare ancillary functions directly to IN2 values
    # # Identify the row number of the IN2 data
    # with open(in2_path) as in2_f:
    #     in2_data = in2_f.readlines()
    # in2_start = [i for i, x in enumerate(in2_data)
    #              if x.startswith(' Mo Da Year ')][0]
    # # Read in the IN2 file using pandas
    # in2_df = pd.read_table(
    #     in2_path, delim_whitespace=True, skiprows=in2_start, header=[0, 1, 2])
    # # Flatten multi-row header
    # in2_df.columns = [
    #     ' '.join(col).replace('-', '').strip()
    #     for col in in2_df.columns.values]
    # in2_df.rename(
    #     columns={'Year': 'YEAR', 'Mo': 'MONTH', 'Da': 'DAY', 'DoY': 'DOY',
    #              'HrMn': 'HOUR'},
    #     inplace=True)
    # in2_df['HOUR'] = (in2_df['HOUR'] / 100).astype(int)
    # # AgriMet times are local with DST (this will drop one hour)
    # in2_df['DATETIME'] = in2_df[['YEAR', 'MONTH', 'DAY', 'HOUR']].apply(
    #     lambda x: pytz.timezone('US/Pacific').localize(dt.datetime(*x)),
    #     axis=1)
    # in2_df['DATETIME'] = in2_df['DATETIME'].apply(
    #     lambda x: x.tz_convert('UTC').strftime('%Y-%m-%d %H:00'))
    # in2_df.set_index('DATETIME', inplace=True, drop=True)

    # Identify the row number of the OUT data
    with open(out_path) as out_f:
        out_data = out_f.readlines()
    out_start = [
        i for i, x in enumerate(out_data) if x.startswith(' Mo Day Yr')][0]
    # Read in the OUT file using pandas (skip header and units)
    out_df = pd.read_table(
        out_path, delim_whitespace=True, index_col=False,
        skiprows=list(range(out_start)) + [out_start + 1])
    out_df.rename(
        columns={'Yr': 'YEAR', 'Mo': 'MONTH', 'Day': 'DAY', 'HrMn': 'HOUR',
                 'Tmax': 'TMAX', 'Tmin': 'TMIN', 'DewP': 'TDEW',
                 'Wind': 'WIND', 'Rs': 'RS'},
        inplace=True)
    out_df['HOUR'] = (out_df['HOUR'] / 100).astype(int)
    # AgriMet times are local with DST (this will drop one hour)
    out_df['DATETIME'] = out_df[['YEAR', 'MONTH', 'DAY', 'HOUR']].apply(
        lambda x: pytz.timezone('US/Pacific').localize(dt.datetime(*x)),
        axis=1)
    out_df['DATETIME'] = out_df['DATETIME'].apply(
        lambda x: x.tz_convert('UTC').strftime('%Y-%m-%d %H:00'))
    # out_df['DATETIME'] = out_df[['YEAR', 'MONTH', 'DAY', 'HOUR']].apply(
    #     lambda x: dt.datetime(*x).strftime('%Y-%m-%d %H:00'), axis=1)
    out_df.set_index('DATETIME', inplace=True, drop=True)

    # Read the station properties from the IN2 file for now
    # The values should probably be read using a regular expression
    for line in out_data:
        if line.strip().startswith('The anemometer height is'):
            zw = float(line.split(':')[1].split()[0])
        elif line.strip().startswith('The weather station elevation is'):
            elev = float(line.split(':')[1].split()[0])
        elif line.strip().startswith('The weather station latitude is'):
            lat = float(line.split(':')[1].split()[0])
        elif line.strip().startswith('The weather station longitude is'):
            lon = float(line.split(':')[1].split()[0])
            if 'West' in line:
                lon *= -1

    # Get list of date strings from the input CSV file
    values, ids = [], []
    for test_date in list(csv_df.index):
        # Datetime that has issues with fcd calculation
        if not test_date.startswith('2015-07-01'):
            continue

        # Only check day time values for now
        if float(csv_df.loc[test_date, 'RS']) <= 1.0:
            continue

        test_dt = dt.datetime.strptime(test_date, '%Y-%m-%d %H:%M')
        # Can the surface type be parameterized inside pytest_generate_tests?
        for surface in ['ETr', 'ETo']:
            date_values = csv_df \
                .loc[test_date, ['TEMP', 'EA', 'RS', 'WIND', 'DOY']] \
                .rename({
                    'DOY': 'doy', 'TEMP': 'tmean', 'EA': 'ea', 'RS': 'rs',
                    'WIND': 'uz'}) \
                .to_dict()
            date_values.update({
                'surface': surface.lower(),
                'expected': out_df.loc[test_date, surface],
                # 'doy': int(test_dt.strftime('%j')),
                'time': test_dt.hour,
                'zw': zw,
                'elev': elev,
                'lat': lat,
                'lon': lon
                # DEADBEEF
                # 'lat': lat * math.pi / 180,
                # 'lon': lon * math.pi / 180
            })
            values.append(date_values)
            ids.append('{}-{}'.format(test_date, surface))


@pytest.fixture(scope='module')
def hourly_data():
    _hourly = HourlyData()
    return _hourly


def pytest_generate_tests(metafunc):
    # Read in inputs for each hourly timestep
    # Set dictionary keys to input variable names
    hourly = hourly_data()

    if 'hourly_params' in metafunc.fixturenames:
        metafunc.parametrize('hourly_params', hourly.values, ids=hourly.ids)


def test_refet_hourly_func_values(hourly_params):
    """Test hourly RefET calculation at a single point and time"""
    # If I don't copy, the pop changes the test values in hourly_data()
    inputs = hourly_params.copy()
    surface = inputs.pop('surface')
    expected = inputs.pop('expected')
    # print('ETr: {}'.format(expected))

    # Cast all numeric inputs to ee.Number type
    inputs = {
        k: ee.Number(v) if k != 'rso_type' else v
        for k, v in inputs.items()}

    if surface.lower() == 'etr':
        assert float(Hourly(
            **inputs).etr().getInfo()) == pytest.approx(expected, abs=0.01)
    elif surface.lower() == 'eto':
        assert float(Hourly(
            **inputs).eto().getInfo()) == pytest.approx(expected, abs=0.01)
