import requests
import numpy as np
import xarray as xr
from pyhdf.SD import SD, SDC
from osgeo import gdal, gdalconst, osr

from modis_tools.auth import ModisSession
from modis_tools.resources import CollectionApi, GranuleApi
from modis_tools.granule_handler import GranuleHandler

import sys
sys.path.append('../libs/')
sys.path.append('libs/')

from combine_path_and_make_dir import *

from pdb import set_trace

# Define latitude and longitude ranges for your region
lat_min, lat_max = 50.0, 60.0
lon_min, lon_max = -11.0, 3.0

temp_dir = '../temp'

username = ""  # Update this line
password = ''  # Update this line

# Authenticate a session
session = ModisSession(username=username, password=password)

# Query the MODIS catalog for collections
collection_client = CollectionApi(session=session)
collections = collection_client.query(short_name="MOD44B", version="061")

# Query the selected collection for granules
granule_client = GranuleApi.from_collection(collections[0], session=session)

# Filter the selected granules via spatial and temporal parameters
nigeria_bbox = [lon_min, lon_max, lat_min, lat_max]
nigeria_granules = granule_client.query(start_date="2016-01-01", end_date="2018-12-31", bounding_box=nigeria_bbox)

# Download the granules
path = combine_path_and_make_dir(temp_dir, 'MODIS')
GranuleHandler.download_from_granules(nigeria_granules, session, path = "../temp/MODIS")

'''
# Define MODIS product URL and download the data
url = "https://e4ftl01.cr.usgs.gov/MOLT/MOD44B.061/2017.03.06/MOD44B.A2017065.h17v03.061.2022284051332.hdf"
response = requests.get(url)
with open("MODIS_VCF.hdf", "wb") as f:
    f.write(response.content)

# Open MODIS HDF file
modis_file = SD("MODIS_VCF.hdf", SDC.READ)

# Read the necessary dataset
dataset = modis_file.select("Majority_Land_Cover_Type_1")
data = dataset.get()
'''
