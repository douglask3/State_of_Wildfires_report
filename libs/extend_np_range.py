import numpy as np
from pdb import set_trace
from itertools import combinations
import itertools
from itertools import repeat

#def rep(elem, each, length_out):
#    return list(itertools.chain.from_iterable(repeat(elem, each) for _ in range(length_out)))
#
def rep_list(input_list, each, length_out):
    repeated_list = [elem for elem in input_list for _ in range(each)]
    num_repeats = length_out // len(repeated_list)
    remainder = length_out % len(repeated_list)
    result = repeated_list * num_repeats + repeated_list[:remainder]
    return result

def generate_combinations(M):
    """Generate all combinations of columns."""
    all_combinations = []
    for r in range(1, M + 1):  # Generate combinations of all lengths from 1 to M
        all_combinations.extend(combinations(range(M), r))
    return all_combinations

def create_Xi(Xo, selected_columns):
    """Create Xi from Xo with selected columns."""    
    N, M = Xo.shape
    num_combinations = len(selected_columns)
    Xi = np.zeros((N, num_combinations))

    for i, cols in enumerate(selected_columns):
        set_trace()
        Xi[:, i] = Xo[:, cols]

    return Xi

def extend_np_range(X, axis = 0, indicies = None, N = 100000, range_extent_scale = 2.0):
    """ Takes 2D np array, changes range of number along an axis and converts to regular gaps.
    Arguments:
        X -- 2d numpy array
        axis -- either 0 or 1, axis over which to perform opperation. 
             If 0, range and spacing changed over col, if 1, over rows.
        N - size of final array with equal spacing.
        range_extent_scale -- numeric. How much to extend range by. 
            2.0 will double. 1.0 maintains current range.
   Returns:
        numpy array of (if axis 0) same no. cols or (if axis 1) same no. rows by N.
    """
    if axis == 1: X = np.transpose(X)

    # Get the range of values for each column
    min_vals = np.min(X, axis = 0)        
    max_vals = np.max(X, axis = 0)
    ranges = max_vals - min_vals

    # Create Xi array
    M = X.shape[1]
    N = int(np.floor(N**(1/len(indicies))))
    
    Xi = np.zeros((N, M))
    
    for i in range(M):
        # Generate 1000 equally spaced values with twice the range
        col_range = ranges[i] * (range_extent_scale - 1.0)
        min_val = min_vals[i] - col_range  # Extend the range by moving min_val left
        max_val = max_vals[i] + col_range  # Extend the range by moving max_val right
        xi_values = np.linspace(min_val, max_val, N)
        
        # Assign values to Xi
        Xi[:, i] = np.sort(xi_values)

    
    #all_combinations = generate_combinations(len(indicies))
    #Xi_list = [create_Xi(Xi, cols) for cols in all_combinations]
    #out = np.stack(Xi_list, axis=2)
    #Xi = X[i:, indicies]
    #out = np.concatenate([Xi[:, i].reshape(-1, 1) for i in range(M)], axis=1)

    
    

    N, M = Xi.shape
    out = np.zeros((N**len(indicies), M))
    index = range(0, N)
    for i in range(len(indicies)):
        id = rep_list(list(index), N**i, N**len(indicies))
        out[:, indicies[i]] = Xi[id, indicies[i]]
        if i == 0:
            for j in np.setdiff1d(np.arange(M), indicies):
                out[:, j] = Xi[id, i]
    

    
    Xi = out
    if axis == 1: Xi = np.transpose(Xi)
    return Xi

