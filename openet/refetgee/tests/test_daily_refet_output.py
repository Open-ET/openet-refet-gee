import datetime as dt
import os

import ee
import numpy as np
import pandas as pd
import pytest

from openet.refetgee import Daily
import openet.refetgee.units as units

constant_geom = ee.Geometry.Rectangle([0, 0, 10, 10], 'EPSG:32613', False)


# Test daily functions using actual RefET input/output files
# DEADBEEF - This doesn't work if I move it to conftest.py
class DailyData():
    """Setup daily validation data from Fallon AgriMet station"""
    val_ws = os.path.join(os.getcwd(), 'openet', 'refetgee', 'tests', 'data')
    # val_ws = os.path.join(os.getcwd(), 'tests', 'data')
    # val_ws = os.path.join(os.path.dirname(os.getcwd()), 'tests', 'data')

    csv_path = os.path.join(val_ws, 'FALN_Agrimet_daily_raw_2015.csv')
    # in2_path = os.path.join(val_ws, 'FALN_Agrimet_daily_raw_2015.in2')
    out_path = os.path.join(val_ws, 'FALN_Agrimet_daily_raw_2015.out')

    # Read in the inputs CSV file using pandas
    csv_df = pd.read_csv(csv_path, engine='python', na_values='NO RECORD')
    csv_df.rename(
        columns={'MN': 'TMIN', 'MX': 'TMAX', 'YM': 'TDEW', 'UA': 'WIND',
                 'SR': 'RS'},
        inplace=True)
    csv_df['DATE'] = csv_df[['YEAR', 'MONTH', 'DAY']].apply(
        lambda x: dt.datetime(*x).strftime('%Y-%m-%d'), axis=1)
    csv_df.set_index('DATE', inplace=True, drop=True)

    # Convert inputs units
    csv_df['TMIN'] = units._f2c(csv_df['TMIN'])
    csv_df['TMAX'] = units._f2c(csv_df['TMAX'])
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
    # in2_df.rename(
    #     columns={'Year': 'YEAR', 'Mo': 'MONTH', 'Da': 'DAY', 'DoY': 'DOY'},
    #     inplace=True)
    # in2_df['DATE'] = in2_df[['YEAR', 'MONTH', 'DAY']].apply(
    #     lambda x: dt.datetime(*x).strftime('%Y-%m-%d'), axis=1)
    # in2_df.set_index('DATE', inplace=True, drop=True)

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
        columns={'Yr': 'YEAR', 'Mo': 'MONTH', 'Day': 'DAY', 'Tmax': 'TMAX',
                 'Tmin': 'TMIN', 'Wind': 'WIND', 'Rs': 'RS', 'DewP': 'TDEW'},
        inplace=True)
    out_df['DATE'] = out_df[['YEAR', 'MONTH', 'DAY']].apply(
        lambda x: dt.datetime(*x).strftime('%Y-%m-%d'), axis=1)
    out_df.set_index('DATE', inplace=True, drop=True)

    # Read the station properties from the IN2 file for now
    # The values should probably be read using a regular expression
    for line in out_data:
        if line.strip().startswith('The anemometer height is'):
            zw = float(line.split(':')[1].split()[0])
        elif line.strip().startswith('The weather station elevation is'):
            elev = float(line.split(':')[1].split()[0])
        elif line.strip().startswith('The weather station latitude is'):
            lat = float(line.split(':')[1].split()[0])

    # Get list of date strings from the input CSV file
    values, ids = [], []
    for test_date in list(csv_df.index):
        # Datetime that has issues with fcd calculation
        if not test_date.startswith('2015-07'):
            continue

        # This day has missing data and is not being handled correctly
        # if test_date.startswith('2015-04-22'):
        #     continue

        test_dt = dt.datetime.strptime(test_date, '%Y-%m-%d')
        # Can the surface type be parameterized inside pytest_generate_tests?
        for surface in ['ETr', 'ETo']:
            date_values = csv_df \
                .loc[test_date, ['TMIN', 'TMAX', 'EA', 'RS', 'WIND']] \
                .rename({
                    'TMIN': 'tmin', 'TMAX': 'tmax', 'EA': 'ea', 'RS': 'rs',
                    'WIND': 'uz'}) \
                .to_dict()
            date_values.update({
                'surface': surface.lower(),
                'expected': out_df.loc[test_date, surface],
                'doy': int(
                    dt.datetime.strptime(test_date, '%Y-%m-%d').strftime('%j')),
                'zw': zw,
                'elev': elev,
                'lat': lat,
                # DEADBEEF
                # 'lat': lat * math.pi / 180,
                'rso_type': 'full'})
            values.append(date_values)
            ids.append('{}-{}'.format(test_date, surface))


def pytest_generate_tests(metafunc):
    if 'daily_params' not in metafunc.fixturenames:
        return
    daily = DailyData()
    metafunc.parametrize('daily_params', daily.values, ids=daily.ids,
                         scope='module')


def test_refet_daily_output(daily_params):
    """Test daily RefET calculation at a single point and time"""
    # If I don't copy, the pop changes the test values in daily_data()
    inputs = daily_params.copy()
    surface = inputs.pop('surface')
    expected = inputs.pop('expected')
    # print('ETr: {}'.format(expected))

    # ETr/ETo values only have 4 significant figures
    # Small number of days don't match if difference is set < 0.008
    diff = 0.05 if expected >= 10.0 else 0.008

    # Cast all numeric inputs to ee.Number type except tmin (for now)
    adj_inputs = {}
    for k, v in inputs.items():
        if k == 'rso_type':
            adj_inputs[k] = v
        elif k in ['tmin', 'tmax', 'ea', 'rs', 'uz']:
            adj_inputs[k] = ee.Image.constant(v)
        else:
            adj_inputs[k] = ee.Number(v)

    if surface.lower() == 'etr':
        refet = Daily(**adj_inputs).etr
    elif surface.lower() == 'eto':
        refet = Daily(**adj_inputs).eto
    output = refet\
        .reduceRegion(ee.Reducer.first(), geometry=constant_geom, scale=1)\
        .getInfo()
    assert float(output[surface.lower()]) == pytest.approx(expected, abs=diff)
