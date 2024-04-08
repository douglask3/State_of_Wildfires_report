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


def dowload_vcf(temp_dir2, bbox):
    # Authenticate a session
    session = ModisSession(username=username, password=password)

    # Query the MODIS catalog for collections
    collection_client = CollectionApi(session=session)
    collections = collection_client.query(short_name="MOD44B", version="061")

    # Query the selected collection for granules
    granule_client = GranuleApi.from_collection(collections[0], session=session)

    # Filter the selected granules via spatial and temporal parameters
     
    granules = granule_client.query(start_date="2000-01-01", end_date="2023-12-31", bounding_box=bbox)

    path = combine_path_and_make_dir(temp_dir1, temp_dir2)

    # Download the granules
    GranuleHandler.download_from_granules(granules, session, path = path)


temp_dir1 = '../temp/'
# Define latitude and longitude ranges for your region
#bbox = [17.50, 33.0, 30.0, 43.0]
#temp_dir2 = 'Greece-MODIS/'
#dowload_vcf(temp_dir2, bbox)

#bbox = [-82.0, -55.0, -55.0, -5.0]
#temp_dir2 = 'BoliviaChile-MODIS/'
#dowload_vcf(temp_dir2, bbox)


#bbox = [-20.0, 50.0, -55.0, 55.0]
#temp_dir2 = 'MED-MODIS/'
#dowload_vcf(temp_dir2, bbox)

#bbox = [-85.0, -60.0, -30.0, 15.0]
#temp_dir2 = 'SouthAmerica-MODIS/'
#dowload_vcf(temp_dir2, bbox)


bbox = [-12.0, 49.0, 3.0, 63.0]
temp_dir2 = 'UK-MODIS/'
dowload_vcf(temp_dir2, bbox)

