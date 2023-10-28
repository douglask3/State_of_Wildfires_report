import arviz as az
import numpy as np
import iris

import os
import sys
sys.path.append('fire_model/')
sys.path.append('libs/')

from combine_path_and_make_dir import * 
from MaxEntFire import MaxEntFire
from iris_plus import insert_data_into_cube

from pdb import set_trace

def select_post_param(trace):
    def select_post_param_name(name): 
        out = trace.posterior[name].values
        A = out.shape[0]
        B = out.shape[1]
        new_shape = ((A * B), *out.shape[2:])
        return np.reshape(out, new_shape)

    params = trace.to_dict()['posterior']
    params_names = params.keys()
    params = [select_post_param_name(var) for var in params_names]
    return params, params_names

def runSim_MaxEntFire(trace, sample_for_plot, X, eg_cube, lmask, run_name, 
                      dir_samples, grab_old_trace, 
                      class_object = MaxEntFire, method = 'burnt_area'):  
    def sample_model(i, run_name = 'control'):   
        dir_sample =  combine_path_and_make_dir(dir_samples, run_name)
        file_sample = dir_sample + '/sample' + str(i) + '.nc'
        
        if os.path.isfile(file_sample) and grab_old_trace:
            return iris.load_cube(file_sample)
        
        print("Generating Sample:" + file_sample)
        param_in = [param[i] if param.ndim == 1 else param[i,:] for param in params]
        param_in = dict(zip(params_names, param_in))
        obj = class_object(param_in)
        out = getattr(obj, method)(X)
        out = insert_data_into_cube(out, eg_cube, lmask)
        coord = iris.coords.DimCoord(i, "realization")
        out.add_aux_coord(coord)
        iris.save(out, file_sample)
        
        return out
    params, params_names = select_post_param(trace)  
    nits = len(trace.posterior.chain)*len(trace.posterior.draw)
    idx = range(0, nits, int(np.floor(nits/sample_for_plot)))
    out = np.array(list(map(lambda id: sample_model(id, run_name), idx)))

    return iris.cube.CubeList(out).merge_cube()

