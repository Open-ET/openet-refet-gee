import base64
import logging
import os

import ee
import google.oauth2.credentials
import pytest


@pytest.fixture(scope="session", autouse=True)
def test_init():
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    logging.getLogger('googleapiclient').setLevel(logging.ERROR)
    logging.debug('Test Setup')

    # For GitHub Actions authenticate using workload identify federation
    if "ACTION_EE_TOKEN" in os.environ:
        creds = google.oauth2.credentials.Credentials(os.getenv('ACTION_EE_TOKEN'))
        ee.Initialize(creds, project=os.getenv('PROJECT_ID'))
    elif 'EE_PRIVATE_KEY_B64' in os.environ:
        print('Writing privatekey.json from environmental variable ...')
        content = base64.b64decode(os.environ['EE_PRIVATE_KEY_B64']).decode('ascii')
        GEE_KEY_FILE = 'privatekey.json'
        with open(GEE_KEY_FILE, 'w') as f:
            f.write(content)
        ee.Initialize(ee.ServiceAccountCredentials('', key_file=GEE_KEY_FILE))
    else:
        ee.Initialize()
