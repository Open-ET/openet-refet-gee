import datetime as dt
from future.standard_library import install_aliases
install_aliases()

from urllib.request import urlopen


# Input parameters
station = 'FALN'
year = 2015
hourly_csv = '{}_Agrimet_hourly_raw_{}.csv'.format(station, year)
daily_csv = '{}_Agrimet_daily_raw_{}.csv'.format(station, year)



# Retrieve the data from the Agrimet site
print('Retrieving raw hourly data')
hourly_url = (
    "https://www.usbr.gov/pn-bin/instant.pl?"
    "station={0}&year={1}&month=1&day=1&year={1}&month=12&day=31&"
    "pcode=OB&pcode=TP&pcode=WS&pcode=SI&"
    # "pcode=OB&pcode=TU&pcode=EA&pcode=TP&pcode=WS&pcode=SI&"
    "print_hourly=1").format(station, year)
hourly_txt = urlopen(hourly_url).read().decode('utf-8')

# Select the data (i.e. strip off header and footer)
hourly_list = [x.strip() for x in hourly_txt.split('\n')]
start_i = hourly_list.index('BEGIN DATA')
end_i = hourly_list.index('END DATA')
hourly_list = hourly_list[start_i+1:end_i]

# Remove station name from column names
hourly_list[0] = hourly_list[0].replace(station, '').replace('DATE', 'DATE,')

# Clean up formatting
# This will probably fail if there are missing values
hourly_list = [
    ','.join([x.strip() for x in line.split(',')])
    for line in hourly_list]

# Change date to separate year, month, day columns
# Change time to hour
hourly_list[0] = hourly_list[0].replace('DATE,TIME', 'YEAR,MONTH,DAY,HOUR')
for i, line in enumerate(hourly_list[1:]):
    line_dt = line.split(',')[0]
    hourly_list[i+1] = line.replace(
        line_dt,
        dt.datetime.strptime(line_dt, '%m/%d/%Y %H:00').strftime('%Y,%m,%d,%H'))

# Write the raw data to a CSV file
print('Writing raw hourly data')
with open(hourly_csv, 'w') as f:
    for line in hourly_list:
        f.write(line + '\n')



# Retrieve the data from the Agrimet site
print('Retrieving raw daily data')
daily_url = (
    "https://www.usbr.gov/pn-bin/daily.pl?"
    "station={0}&year={1}&month=1&day=1&year={1}&month=12&day=31&"
    "pcode=MN&pcode=MX&pcode=SR&pcode=YM&pcode=UA&"
    "pcode=ETRS&pcode=ETOS").format(station, year)
daily_txt = urlopen(daily_url).read().decode('utf-8')

# Select the data (i.e. strip off header and footer)
daily_list = [x.strip() for x in daily_txt.split('\n')]
start_i = daily_list.index('BEGIN DATA')
end_i = daily_list.index('END DATA')
daily_list = daily_list[start_i+1:end_i]

# Remove station name from column names
daily_list[0] = daily_list[0].replace(station, '')

# Clean up formatting
# This will probably fail if there are missing values
daily_list = [
    ','.join([x.strip() for x in line.split(',')])
    for line in daily_list]

# Change date to separate year, month, day columns
# Change time to hour
daily_list[0] = daily_list[0].replace('DATE', 'YEAR,MONTH,DAY')
for i, line in enumerate(daily_list[1:]):
    line_dt = line.split(',')[0]
    daily_list[i+1] = line.replace(
        line_dt,
        dt.datetime.strptime(line_dt, '%m/%d/%Y').strftime('%Y,%m,%d'))

# Write the raw data to a CSV file
print('Writing raw daily data')
with open(daily_csv, 'w') as f:
    for line in daily_list:
        f.write(line + '\n')
