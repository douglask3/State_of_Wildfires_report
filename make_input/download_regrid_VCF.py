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
if new_download:
    # Download the granules
    GranuleHandler.download_from_granules(granules, session, path = path)




import os
import numpy as np
from osgeo import gdal
from pyresample import geometry, kd_tree, load_area, create_area_def

def regrid_modis(input_file, output_file, target_resolution=0.5):
    # Open MODIS HDF file
    
    dataset = gdal.Open(input_file, gdal.GA_ReadOnly)
    if dataset is None:
        print("Could not open file:", input_file)
        return

    # Read data and geotransform
    data = dataset.ReadAsArray()
    geotransform = dataset.GetGeoTransform()
    projection = dataset.GetProjection()

    # Define MODIS projection information
    modis_area_id = 'MODIS_Swath_Type_GEO'
    modis_proj_dict = {'a': geotransform[1], 'b': 0.0, 'c': geotransform[0],
                       'd': 0.0, 'e': geotransform[5], 'f': geotransform[3]}

    # Define area definition for MODIS
    modis_width = dataset.RasterXSize
    modis_height = dataset.RasterYSize
    modis_area_extent = (geotransform[0], geotransform[3] + geotransform[5]*modis_height, 
                         geotransform[0] + geotransform[1]*modis_width, geotransform[3])
    modis_area_def = geometry.AreaDefinition(modis_area_id, 'MODIS swath', 'MODIS swath', modis_proj_dict,
                                         modis_width, modis_height, modis_area_extent)

    # Define target area definition (0.5 degrees)
    target_area_id = 'Equidistant_Cylindrical_0.5deg'
    target_area_extent = (-180, -90, 180, 90)
    target_area_def = geometry.AreaDefinition(target_area_id, 'Global Grid 0.5x0.5 degree', 'Global Grid 0.5x0.5 degree',
                                              {'proj': 'eqc', 'lon_0': 0, 'lat_0': 0}, target_resolution, target_resolution,
                                              target_area_extent)

    # Resample MODIS data to target resolution
    modis_swath_data = geometry.SwathDefinition(lons=modis_area_def.get_lon(), lats=modis_area_def.get_lat())
    target_swath_data = geometry.SwathDefinition(lons=target_area_def.get_lon(), lats=target_area_def.get_lat())
    modis_to_target = kd_tree.resample_nearest(modis_area_def, data, target_area_def, radius_of_influence=5000, fill_value=None)
    modis_regridded = kd_tree.get_sample_from_neighbour_info('nn', target_swath_data, modis_to_target)

    # Write regridded data to GeoTIFF
    driver = gdal.GetDriverByName("GTiff")
    outdata = driver.Create(output_file, modis_regridded.shape[1], modis_regridded.shape[0], 1, gdal.GDT_Float32)
    outdata.SetGeoTransform((target_area_extent[0], target_resolution, 0, target_area_extent[3], 0, -target_resolution))
    outdata.SetProjection(projection)
    outdata.GetRasterBand(1).WriteArray(modis_regridded)
    outdata.FlushCache()

    print("Regridded file saved as:", output_file)

# Example usage
input_directory = path
output_directory = path = combine_path_and_make_dir('../data/data/', 'MODIS')

# Loop through MODIS HDF files in the input directory
for file_name in os.listdir(input_directory):
    if file_name.endswith(".hdf"):
        input_file = os.path.join(input_directory, file_name)
        output_file = os.path.join(output_directory, file_name.replace(".hdf", "_regridded.tif"))
        
        regrid_modis(input_file, output_file)






import iris
import numpy as np
from iris.util import equalise_attributes

# Define the target grid (0.5-degree resolution)
target_lons = np.arange(lon_min, lon_max, 0.5)
target_lats = np.arange(lat_min, lat_max, 0.5)

#cube_list = iris.load_raw(path + '*.hdf')

hdf_files = [f for f in os.listdir(path) if f.endswith('.hdf')]
def open_dataset(file):
    set_trace()
    xr.open_dataset(file, engine = 'h5netcdf')

ds_list = [open_dataset(path + file) for file in hdf_files]
set_trace()

