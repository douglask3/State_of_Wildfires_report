
#Region mask notes:
#https://regionmask.readthedocs.io/en/v0.9.0/
#https://regionmask.readthedocs.io/en/v0.9.0/notebooks/mask_2D.html
#https://regionmask.readthedocs.io/en/v0.9.0/defined_scientific.html

import iris
import iris.coord_categorisation as icc
import cartopy.io.shapereader as shpreader
try:
    from ascend import shape
except:
    print("WARNING: can't load shape from ascent. Some constraints calls might fail")
import numpy as np
import cartopy.crs as ccrs
import geopandas as gp
import regionmask
from pdb import set_trace
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import datetime
import geopandas as gpd

def constrain_cube_to_lonlat_range(cube, lon_min, lon_max, lat_min, lat_max):
    
    constraint = iris.Constraint(latitude=lambda cell: (lat_min) <= cell <= (lat_max), 
                                 longitude=lambda cell: (lon_min) <= cell <= (lon_max))
    constrained_cube = cube.extract(constraint)

    return(constrained_cube)

def constrain_to_data(cube):
    """constrain cube to the lats and lons that contain data that isn't 'nan'   
    Arguments:
        cube -- iris cube with 'latitude' and 'londitude' coordinates. 
                Can have extra dimentions.
    Returns:
        cube of same or smaller extent restircited to range of lats and lons where data.
        cube has useable data.
    """
    lats = cube.coord('latitude').points
    lons = cube.coord('longitude').points

    lat_valid = lats[~np.isnan(cube.data).all(axis=(0,2))]
    lon_valid = lons[~np.isnan(cube.data).all(axis=(0,1))]

    lat_min, lat_max = np.min(lat_valid), np.max(lat_valid)
    lon_min, lon_max = np.min(lon_valid), np.max(lon_valid)

    constraint = iris.Constraint(latitude=lambda cell: (lat_min) <= cell <= (lat_max), 
                                 longitude=lambda cell: (lon_min) <= cell <= (lon_max))
    constrained_cube = cube.extract(constraint)

    return constrain_cube_to_lonlat_range(cube, lon_min, lon_max, lat_min, lat_max)
    
def ar6_region(cube, region_code):
    """constrain cube an ar6 region
    Arguments:
        cube -- iris cube with 'latitude' and 'londitude' coordinates. 
                Can have extra dimentions.
	region_code -- AR6 region, which can be defined either as a number, 
                region code or region name. Can be a list of multiple regions if you want to 
                constrain to more than one
    Returns:
        cube of same or smaller extent restircited to range of lats and lons with areas outside 
	of AR6 region masked out.
    """
    lats = cube.coord('latitude').points
    lons = cube.coord('longitude').points
    
    if not isinstance(region_code, list): region_code = [region_code]
    region_code = [regionmask.defined_regions.ar6.all.region_ids[rc] for rc in region_code] 
    
    mask = regionmask.defined_regions.ar6.all.mask(lons, lats)
    region_mask = mask.isin(region_code)

    true_indices = region_mask.where(region_mask, drop=True).stack(index=('lat', 'lon')).index
    
    min_lat = true_indices['lat'].min().values
    max_lat = true_indices['lat'].max().values
    min_lon = true_indices['lon'].min().values
    max_lon = true_indices['lon'].max().values

    
    cube = constrain_cube_to_lonlat_range(cube, min_lon, max_lon, min_lat, max_lat)
    constrained_mask = region_mask.sel(lat=slice(max_lat, min_lat), lon=slice(min_lon, max_lon))
    
    masked_data = np.where(constrained_mask, cube.data, np.nan)
    masked_cube = cube.copy()   
    masked_cube.data = masked_data
    masked_cube.data[ masked_cube.data>1e15] = np.nan
    
    cube_out = constrain_to_data(masked_cube)
    
    return cube_out

def add_lat_lon_bounds(cube):
    """add latitude and longitude bounds to a cube. Attempts each corrdinate in turn and add 
        bounds if the coordinate exists but the bound doesn't
    Arguments:
        cube -- iris cube
    Returns:
        cube with latitude and longitude bounds
    """
    try:
        cube.coord('latitude').guess_bounds()
    except:
        pass
    try:
        cube.coord('longitude').guess_bounds()
    except:
        pass
    return cube

def annual_average(cube, annual_aggregate = None):
    """calculates cube annual average
    Arguments:
        cube -- iris cube
        annual_aggregate -- If None, does nothing. Otherwise, iris.analysis function for
            aggratating years. If you want annual summed fractions, for example (i.e Burnt area)
            use iris.analysis.SUM. If you want average flux, for example, use iris.analysis.MEAN
    Returns:
        cube of averaged over time, noting annual_agregrate. And 2 item list of year 
        range of original cube
    """
    cube = add_lat_lon_bounds(cube)
    if annual_aggregate is not None:
        cube = cube.copy().aggregated_by('year', annual_aggregate)
    
    year_range = [np.min(cube.coord('year').points), np.max(cube.coord('year').points)]
    cube = cube.collapsed('time',  iris.analysis.MEAN)
    
    return cube, year_range

def make_time_series(cube, annual_aggregate = None, year_range = None):
    """converts cube into collapse time series cube.
    Arguments:
        cube -- iris cube
        annual_aggregate -- Boolean. If True, time series is annual average.
        year_range -- range of time series
    Returns:
        cube of summed over time
    """
    if year_range is not None:
        cube = sub_year_range(cube, year_range)
    
    cube = add_lat_lon_bounds(cube)

    ## in km2
    weights = iris.analysis.cartography.area_weights(cube)/1000000000000.0 

    cube.data[np.isnan(cube.data)] = 0.0
    collapsed_cube = cube.collapsed(['latitude', 'longitude'], iris.analysis.SUM, 
                                    weights=weights)
    
    if annual_aggregate is not None:
        collapsed_cube = collapsed_cube.aggregated_by('year', annual_aggregate)    
        
    return collapsed_cube

def sub_year_range(cube, year_range):
    """Selects months of a year from data   
    Arguments:
        cube -- iris cube with time array with year information.
        year_range -- numeric list of first to last year to cut
    Returns:
        cube of just years between to years provided.
    """
    
    try:
        icc.add_year(cube, 'time')
    except:
        pass
    
    constraint = iris.Constraint(year=lambda cell: (year_range[0]-0.95) <= cell <= (year_range[1]+0.95))
    
    return cube.extract(constraint)
    
    
def sub_year_months(cube, months_of_year):
    """Selects months of a year from data   
    Arguments:
        data -- iris cube with time array we can add add_month_number too.
        months_of_year -- numeric, month of the year you are interested in
                from 0 (Jan) to 11 (Dec)
    Returns:
        cube of just months we are interested in.
    """
    try: 
        icc.add_month_number(cube, 'time')
    except:
        pass  
           
    months_of_year = np.array(months_of_year)+1
    season = iris.Constraint(month_number = lambda cell, mnths = months_of_year: \
                             np.any(np.abs(mnths - cell[0])<0.5))
    return cube.extract(season)

def constrain_cube_by_cube_and_numericIDs(cube, regions, region):
    """constrains a cube to region identifies in 'mask'
        Assumes that the cubes aere iris and on the same grid
    Arguments:
        cube -- an iris cube 
        mask -- an iris cube at same resultion as 'cube' and where regions area identified
                by number
        region -- numeric list of regions to pick in 'mask'
            You can pick more than one:
            
    Returns:
    Input cube with areas outside of selected regions masked out, 
    constrained to that region
    """
    
    cube = add_lat_lon_bounds(cube)
    regions = add_lat_lon_bounds(regions)
    regions = regions.regrid(cube, iris.analysis.Linear())        
    
    mask = regions.copy()
    mask.data[:] = 0.0
    for reg in region: mask.data += regions.data == reg
    
    mask.data[mask.data.mask] = 0.0
    mask = mask.data == 0

    if not cube.data.mask.shape == cube.data.shape:
         cube.data.mask = np.isnan(cube.data)
    for layer in cube.data:
        
        layer.mask[mask] = False
        layer[mask] = np.nan
    
    cube_out = constrain_to_data(cube)
    return cube_out

def constrain_GFED(cube, region, *args, **kw):
    """constrains a cube to GFED region
        Assumes that the cube is iris and on a 0.5 degree grid
    Arguments:
        cube -- an iris cube at 0.5 degrees
        region -- numeric list (i.e [3, 7, 8]) where numbers pick GFED region.
            You can pick more than one:
            1 BONA
            2 TENA
            3 CEAM
            4 NHSA
            5 SHSA
            6 EURO
            7 MIDE
            8 NHAF
            9 SHAF
            10 BOAS
            11 CEAS
            12 SEAS
            13 EQAS
            14 AUST
    Returns:
    Input cube with areas outside of selected GFEDregions masked out, 
    constrained to that region
    """
    
    regions = iris.load_cube('data/GFEDregions.nc')
    return constrain_cube_by_cube_and_numericIDs(cube, regions, region)

def constrain_olson(cube, ecoregions):
    """constrains a cube to Olson ecoregion
    Assumes that the cube is iris and on a 0.5 degree grid

    Arguments:

    cube -- an iris cube at 0.5 degrees
    ecoregions -- numeric list (i.e [3, 7, 8]) where numbers pick Olson biomes.
        You can pick more than one:
            1 Tropical and subtropical moist broadleaf forests
            2 Tropical and subtropical dry broadleaf forests
            3 Tropical and suptropical coniferous forests
            4 Temperate broadleaf and mixed forests
            5 Temperate Coniferous Forest
            6 Boreal forests / Taiga
            7 Tropical and subtropical grasslands, savannas and shrublands
            8 Temperate grasslands, savannas and shrublands
            9 Flooded grasslands and savannas
            10 Montane grasslands and shrublands
            11 Tundra
            12 Mediterranean Forests, woodlands and scrubs
            13 Deserts and xeric shrublands
            14 Mangroves

    Returns:
    Input cube with areas outside of selected Olson biomes masked out, 
    constrained to that region
    """
    biomes = iris.load_cube('data/wwf_terr_ecos_0p5.nc')
    return constrain_cube_by_cube_and_numericIDs(cube, biomes, ecoregions)

def constrain_natural_earth(cube, Country, Continent = None, shpfilename = None, 
                            constrain = True, *args, **kw):
    
    """constrains a cube to Natural Earth Country or continent.
    Assumes that the cube is iris
    If Country is defined, it wont select Continent.

    Arguments:

    cube -- an iris cube
    Continent -- name of continent. Options are:
        'South America'
        'Oceania'
        'Europe'
        'Afria'
        'North America'
        'Asia' 
    shpfilename -- path and filename of natural earth shapefile.
                   If set to None, it will look for temp. file version and 
                   download if it does not exist.
    Returns:
    Input cube constrained to the extent to that country or continent, 
    and mask areas outside of it. Uses Natural Earth
    """
    if shpfilename is None:
        shpfilename = shpreader.natural_earth(resolution='110m', 
                                              category='cultural', name='admin_0_countries')
    natural_earth_file = shape.load_shp(shpfilename)
    if Country is not None:
        NAMES = [i.attributes.get('NAME') for i in natural_earth_file]
        NAME = [s for s in NAMES if Country in s][0]
        CountrySelect = shape.load_shp(shpfilename, NAME=NAME)
    elif Continent is not None:
        CountrySelect = shape.load_shp(shpfilename, Continent='South America')
        CountrySelect = Country.unary_union()
    
    if constrain: cube = CountrySelect[0].constrain_cube(cube)
    
    cube = CountrySelect[0].mask_cube(cube)
    return cube

def constrain_region(cube, ecoregions = None, Country = None, Continent = None, *args, **kw):
    """ checks if any spatial constrains are set and contrains according to comments in functions above

    Arguments:

    cube -- an iris cube
    ecoregions -- see constrain_olson help
    Country, Continent -- see constrain_natural_earth help
      
    Returns:
    Input cube constrained to the extent to defined country or continent, 
    and mask areas outside of it and defined ecoregions.
    """
    if ecoregions is not None:
        cube = constrain_olson(cube, ecoregions)

    if Continent is not None or Country is not None:
        cube = constrain_natural_earth(cube, Country, Continent, *args, **kw)

    return(cube)


def contrain_coords(cube, extent):
    
    longitude_constraint = iris.Constraint(longitude=lambda cell: extent[0] <= cell.point <= extent[1])
    latitude_constraint = iris.Constraint(latitude=lambda cell: extent[2] <= cell.point <= extent[3])
    
    return cube.extract(longitude_constraint & latitude_constraint)

    
def constrain_BR_biomes(cube, biome_ID):
    if len(biome_ID) == 1 and biome_ID[0] == 0: return cube
    mask = iris.load_cube('data/BR_Biomes.nc') 
    #mask = iris.load_cube('data/Pantanal_basin.nc') 
    return constrain_cube_by_cube_and_numericIDs (cube, mask, biome_ID)



