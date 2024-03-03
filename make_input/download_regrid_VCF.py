import requests
import numpy as np
import xarray as xr
from pyhdf.SD import SD, SDC
from osgeo import gdal, gdalconst, osr

from modis_tools.auth import ModisSession
from modis_tools.resources import CollectionApi, GranuleApi
from modis_tools.granule_handler import GranuleHandler
from earthdata_login_details import username, password

import sys
sys.path.append('../libs/')
sys.path.append('libs/')

from combine_path_and_make_dir import *

from pdb import set_trace

# Define latitude and longitude ranges for your region
lat_min, lat_max = -90.0, 90.0
lon_min, lon_max = -180.0, 180.0

temp_dir = '../temp'
new_download = True

# Authenticate a session
session = ModisSession(username=username, password=password)

# Query the MODIS catalog for collections
collection_client = CollectionApi(session=session)
collections = collection_client.query(short_name="MOD44B", version="061")

# Query the selected collection for granules
granule_client = GranuleApi.from_collection(collections[0], session=session)

# Filter the selected granules via spatial and temporal parameters
bbox = [lon_min, lat_min, lon_max, lat_max]
granules = granule_client.query(start_date="2000-01-01", end_date="2023-12-31", bounding_box=bbox)

path = combine_path_and_make_dir(temp_dir, 'glob-MODIS')

# Download the granules
GranuleHandler.download_from_granules(granules, session, path = path)

