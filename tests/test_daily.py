import datetime as dt
import math
import os

import ee
import numpy as np
import pandas as pd
import pytest

from geerefet.daily import Daily
import geerefet.calcs as calcs
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
    'pair': 87.81876435813037,
    'pair_asce': 87.80710537212929,
    'zw': 3.0,
}
# Daily test parameters for 2015-07-01
d_args = {
    'doy': 182,
    'ea': 1.2206674169951346,
    'ea_asce': 1.2205053588697359,
    'eto': 7.9422320475712835,
    'etr': 10.571314344056955,
    'etr_asce': 10.626087665395694,
    'etr_rso_simple': 10.628137858930051,
    'q': 0.008691370735727117,
    'rs': 674.07 * 0.041868,  # Conversion from Langleys to MJ m-2
    'rso': 31.565939444861765,
    'tdew': units._f2c(49.84),
    'tmin': units._f2c(66.65),
    'tmax': units._f2c(102.80),
    # 'tmean': f2c(84.725),
    'uz': 4.80 * 0.44704,  # Conversion from mph to m s-1
    'u2': 1.976111757722194,
}

# # Playing around with passing a dictionary of function keyword arguments
# daily_args = {
#     k: v for k, v in daily.items()
#     if k in ['tmin', 'tmax', 'ea', 'rs', 'uz', 'zw', 'doy']}
# daily_args.update({'surface':'etr', 'rso_type': 'full', 'rso': None})


# Test full daily functions with positional inputs
def test_refet_daily_input_positions():
    etr = Daily(
        ee.Number(d_args['tmax']), ee.Number(d_args['tmin']),
        ee.Number(d_args['ea']), ee.Number(d_args['rs']),
        ee.Number(d_args['uz']), ee.Number(s_args['zw']),
        ee.Number(s_args['elev']), ee.Number(s_args['lat']),
        ee.Number(d_args['doy']), method='refet').etr().getInfo()
    assert float(etr) == pytest.approx(d_args['etr'])


# Test full daily calculations with keyword inputs
# Test surface, rso_type, and rso inputs
def test_refet_daily_surface_etr():
    etr = Daily(
        tmax=ee.Number(d_args['tmax']), tmin=ee.Number(d_args['tmin']),
        ea=ee.Number(d_args['ea']), rs=ee.Number(d_args['rs']),
        uz=ee.Number(d_args['uz']), zw=ee.Number(s_args['zw']),
        elev=ee.Number(s_args['elev']), lat=ee.Number(s_args['lat']),
        doy=ee.Number(d_args['doy']), method='refet').etr().getInfo()
    assert float(etr) == pytest.approx(d_args['etr'])


def test_refet_daily_surface_eto():
    eto = Daily(
        tmax=ee.Number(d_args['tmax']), tmin=ee.Number(d_args['tmin']),
        ea=ee.Number(d_args['ea']), rs=ee.Number(d_args['rs']),
        uz=ee.Number(d_args['uz']), zw=ee.Number(s_args['zw']),
        elev=ee.Number(s_args['elev']), lat=ee.Number(s_args['lat']),
        doy=ee.Number(d_args['doy']), method='refet').eto().getInfo()
    assert float(eto) == pytest.approx(d_args['eto'])


def test_refet_daily_rso_type_simple():
    etr = Daily(
        tmax=ee.Number(d_args['tmax']), tmin=ee.Number(d_args['tmin']),
        ea=ee.Number(d_args['ea']), rs=ee.Number(d_args['rs']),
        uz=ee.Number(d_args['uz']), zw=ee.Number(s_args['zw']),
        elev=ee.Number(s_args['elev']), lat=ee.Number(s_args['lat']),
        doy=ee.Number(d_args['doy']), method='refet',
        rso_type='simple').etr().getInfo()
    assert float(etr) == pytest.approx(d_args['etr_rso_simple'])


def test_refet_daily_rso_type_array():
    etr = Daily(
        tmax=ee.Number(d_args['tmax']), tmin=ee.Number(d_args['tmin']),
        ea=ee.Number(d_args['ea']), rs=ee.Number(d_args['rs']),
        uz=ee.Number(d_args['uz']), zw=ee.Number(s_args['zw']),
        elev=ee.Number(s_args['elev']), lat=ee.Number(s_args['lat']),
        doy=ee.Number(d_args['doy']), method='refet',
        rso_type='array', rso=ee.Number(d_args['rso'])).etr().getInfo()
    assert float(etr) == pytest.approx(d_args['etr'])


def test_refet_daily_rso_type_exception():
    with pytest.raises(ValueError):
        etr = Daily(
            tmax=ee.Number(d_args['tmax']), tmin=ee.Number(d_args['tmin']),
            ea=ee.Number(d_args['ea']), rs=ee.Number(d_args['rs']),
            uz=ee.Number(d_args['uz']), zw=ee.Number(s_args['zw']),
            elev=ee.Number(s_args['elev']), lat=ee.Number(s_args['lat']),
            doy=ee.Number(d_args['doy']), rso_type='nonsense',
            method='refet')
        # assert float(etr) == pytest.approx(d_args['etr'])


def test_refet_daily_asce():
    etr = Daily(
        tmax=ee.Number(d_args['tmax']), tmin=ee.Number(d_args['tmin']),
        ea=ee.Number(d_args['ea']), rs=ee.Number(d_args['rs']),
        uz=ee.Number(d_args['uz']), zw=ee.Number(s_args['zw']),
        elev=ee.Number(s_args['elev']), lat=ee.Number(s_args['lat']),
        doy=ee.Number(d_args['doy']), method='asce').etr().getInfo()
    assert float(etr) == pytest.approx(d_args['etr_asce'])


def test_refet_daily_gridmet():
    # Convert test values to GRIDMET units
    # MJ m-2 d-1 -> W m-2 and C -> K
    #
    gridmet_img = ee.Image.constant([
            d_args['tmax'] + 273.15, d_args['tmin'] + 273.15, d_args['q'],
            d_args['rs'] / 0.0864, d_args['uz']])\
        .rename(['tmmx', 'tmmn', 'sph', 'srad', 'vs'])\
        .set('system:time_start', ee.Date('2015-07-01').millis())
    refet = Daily.gridmet(
        ee.Image(gridmet_img), elev=ee.Number(s_args['elev']),
        lat=ee.Number(s_args['lat']), method='asce')
    # Output bands are not currently being renamed, defaults to 'vs'
    etr = refet.etr().reduceRegion(
        reducer=ee.Reducer.first(),
        geometry=ee.Geometry.Rectangle([0, 0, 10, 10], 'EPSG:32613', False),
        scale=1).getInfo()['vs']
    print(etr)
    assert float(etr) == pytest.approx(d_args['etr_asce'])


# Test daily functions using actual RefET input/output files
# DEADBEEF - This doesn't work if I move it to conftest.py
class DailyData():
    """Setup daily validation data from Fallon AgriMet station"""
    val_ws = os.path.join(os.getcwd(), 'tests', 'data')
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


@pytest.fixture(scope='module')
def daily_data():
    _daily = DailyData()
    return _daily


def pytest_generate_tests(metafunc):
    # Read in inputs for each daily timestep
    # Set dictionary keys to input variable names
    daily = daily_data()

    if 'daily_params' in metafunc.fixturenames:
        metafunc.parametrize('daily_params', daily.values, ids=daily.ids)


def test_refet_daily_values(daily_params):
    """Test daily RefET calculation at a single point and time"""
    # If I don't copy, the pop changes the test values in daily_data()
    inputs = daily_params.copy()
    surface = inputs.pop('surface')
    expected = inputs.pop('expected')
    # print('ETr: {}'.format(expected))

    # ETr/ETo values only have 4 significant figures
    # Small number of days don't match if difference is set < 0.008
    diff = 0.05 if expected >= 10.0 else 0.008

    # Cast all numeric inputs to ee.Number type
    inputs = {
        k: ee.Number(v) if k != 'rso_type' else v
        for k, v in inputs.items()}

    if surface.lower() == 'etr':
        assert float(Daily(
            **inputs).etr().getInfo()) == pytest.approx(expected, abs=diff)
    elif surface.lower() == 'eto':
        assert float(Daily(
            **inputs).eto().getInfo()) == pytest.approx(expected, abs=diff)
