import numpy as np
from pdb import set_trace

import iris
from libs.iris_plus import *
from libs.constrain_cubes_standard import *
from libs.read_variable_from_netcdf import *

def read_variable_from_netcdf(filename, dir = '', subset_function = None, 
                              make_flat = False, units = None, 
                              subset_function_args = None,
                              time_series = None, 
                              time_points = None, return_time_points = False,
                              *args, **kw):
    """Read data from a netCDF file 
        Assumes that the variables in the netcdf file all have the name "variable"
        Assunes that values < -9E9, you dont want. This could be different in some circumstances
    Arguments:
        filename -- a string with filename or two element python list containing the name of 
            the file and the target variable name. If just the sting of "filename" 
            assumes variable name is "variable"
        dir -- directory file is in. Path can be in "filename" and None means no 
            additional directory path needed.
        subset_function -- a function or list of functions to be applied to each data set
        subset_function_args -- If subset_function is a function, dict arguments for that function. 
                            If subset_function is a list, a  list or dict constaining arguments 
                                for those functions in turn.
        make_flat -- should the output variable to flattened or remain cube
        time_series -- list comtaining range of years. If making flat and returned a time series, 
            checks if that time series contains year.
    Returns:
        Y - if make_flat, a numpy vector of the target variable, otherwise returns iris cube
    """

    print("Opening:")
    print(filename)
    if filename[0] == '~' or filename[0] == '/' or filename[0] == '.': dir = ''
    try:
        if isinstance(filename, str):        
            dataset = iris.load_cube(dir + filename, callback=sort_time)
        else:
            dataset = iris.load_cube(dir + filename[0], filename[1], callback=sort_time)
    except:
        try:
            dataset = iris.load_cube(dir + filename)
        except:
            print("==============\nERROR!")
            print("can't open data.")
            print("Check directory (''" + dir + "''), filename (''" + filename + \
              "'') or file format")
            print("==============")
            set_trace()
    coord_names = [coord.name() for coord in dataset.coords()]
    if dataset is None: return None
    if time_points is not None:     
        if 'time' in coord_names:
            dataset = dataset.interpolate([('time', time_points)], iris.analysis.Linear())
        else:   
            def addTime(time_point):
                time = iris.coords.DimCoord(np.array([time_point]), standard_name='time',
                                            units = 'days since 1661-01-01 00:00:00')
                dataset_cp = dataset.copy()
                dataset_cp.add_aux_coord(time)
                return dataset_cp

            dataset_time = [addTime(time_point) for time_point in time_points]
            dataset = iris.cube.CubeList(dataset_time).merge_cube()
            dataset0 = dataset.copy()
    if units is not None: dataset.units = units
    if subset_function is not None:
        if isinstance(subset_function, list):
            for FUN, args in zip(subset_function, subset_function_args):
                try:    
                    dataset = FUN(dataset, **args)
                except:
                    print("Warning! function: " + FUN.__name__ + " not applied to file: " + \
                          dir + filename)
        else:      
            dataset = subset_function(dataset, **subset_function_args) 
    if return_time_points: time_points = dataset.coord('time').points 
    
    
    if make_flat: 
        if time_series is not None: years = dataset.coord('year').points
        
        try:    
            dataset = dataset.data.flatten()
        except:
            set_trace()
            
        if time_series is not None:
            if not years[ 0] == time_series[0]:
                dataset = np.append(np.repeat(np.nan, years[ 0]-time_series[0]), dataset)
            if not years[-1] == time_series[1]:
                dataset = np.append(dataset, np.repeat(np.nan, time_series[1]-years[-1]))
            if return_time_points: set_trace()
        
    if return_time_points: dataset = (dataset, time_points)
    
    return dataset

def read_all_data_from_netcdf(y_filename, x_filename_list, CA_filename = None, add_1s_columne = False, 
                              y_threshold = None, x_normalise01 = False, scalers = None,
                              check_mask = True, frac_random_sample = 1.0, 
                              min_data_points_for_sample = None, *args, **kw):
                              
    """Read data from netCDF files 
        
    Arguments:
        y_filename -- a two element python list containing the name of the file and the target 
            variable name
        x_filename_list -- a python list of filename containing the feature variables
        CA_filename -- a python list of filename containing the area of the cover type
        y_threshold -- if converting y into boolean, the threshold we use to spit into 
            0's and 1's
        add_1s_columne -- useful for if using for regressions. Adds a variable of 
            just 1's t rperesent y = SUM(a_i * x_i) + c
        x_normalise01 -- Boolean. If True, then X's are normalised between 0 and 1.
        scalers -- None or numpy array of shape 2 by n. columns of X.
            Defines what scalers (min and max) to apply to each X column. 
            If None, doesn't apply anything.
        check_mask -- Boolean. If True, simple checks if there are any large negtaive numbers 
            and makes them out. Assunes that values < -9E9, you dont want. 
            This could be different in some circumstances
        frac_random_sample -- fraction of data to be returned
        see read_variable_from_netcdf comments for *arg and **kw.
    Returns:
        Y - a numpy array of the target variable
        X - an n-D numpy array of the feature variables 
    """
    Y, time_points = read_variable_from_netcdf(y_filename, make_flat = True, *args, return_time_points = True, **kw)
    
    if CA_filename is not None:
        CA = read_variable_from_netcdf(CA_filename, make_flat = True, 
                                       time_points = time_points, *args, **kw)
   
    # Create a new categorical variable based on the threshold
    if y_threshold is not None:
        Y = np.where(Y >= y_threshold, 0, 1)
        #count number of 0 and 1 
        counts = np.bincount(Y)
        #print(f"Number of 0's: {counts[0]}, Number of 1's: {counts[1]}")
    
    n=len(Y)
    m=len(x_filename_list)
    
    X = np.zeros([n,m])
    
    for i, filename in enumerate(x_filename_list):
        X[:, i]=read_variable_from_netcdf(filename, make_flat = True, time_points = time_points,
                                          *args, **kw)
    
    if add_1s_columne: 
        X = np.column_stack((X, np.ones(len(X)))) # add a column of ones to X 
    
    if check_mask:
        if CA_filename is not None:
            cells_we_want = np.array([np.all(rw > -9e9) and np.all(rw < 9e9) for rw in np.column_stack((X, Y, CA))])
            CA = CA[cells_we_want]
        else:
            cells_we_want = np.array([np.all(rw > -9e9) and np.all(rw < 9e9) for rw in np.column_stack((X,Y))])
        Y = Y[cells_we_want]
        X = X[cells_we_want, :]
        
    if x_normalise01: 
        try:
            scalers = np.array([np.min(X, axis=0), np.max(X, axis=0)])
        except:
            set_trace()
        squidge = (scalers[1,:]-scalers[0,:])/(X.shape[0])
        scalers[0,:] = scalers[0,:] - squidge
        scalers[1,:] = scalers[1,:] + squidge
        
        test = scalers[1,:] == scalers[0,:]
        scalers[0,test] = 0.0
        scalers[1,test] = 1.0
    

    if frac_random_sample is None: 
        frac_random_sample = 1000
    else:
        if min_data_points_for_sample is not None:
            min_data_frac = min_data_points_for_sample/len(Y)
            if min_data_frac > frac_random_sample: frac_random_sample = min_data_frac
    
    if frac_random_sample < 1:
        M = X.shape[0]
        selected_rows = np.random.choice(M, size = int(M * frac_random_sample), replace=False)
        Y = Y[selected_rows]
        X = X[selected_rows, :]
        if CA_filename is not None:
            CA = CA[selected_rows]
    
    if scalers is not None:
        X = (X-scalers[0, :]) / (scalers[1, :] - scalers[0, :])
        if check_mask: 
            if CA_filename is not None: return Y, X, CA, cells_we_want, scalers
        return Y, X, cells_we_want, scalers
        
    if check_mask or frac_random_sample: 
        if CA_filename is not None: return Y, X, CA, cells_we_want
    return Y, X, cells_we_want
    
    if CA_filename is not None: return Y, X, CA
    return Y, X 
