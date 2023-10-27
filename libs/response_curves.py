
def non_masked_data(cube):
        return cube.data[cube.data.mask == False].data

def standard_response_curve(Sim, x, eg_cube, id_name, dir_samples,
                            run_name, grab_old_trace):  
    for col_to_keep in range(X.shape[1]-1):
        other_cols = np.arange(X.shape[1]-1)  # Create an array of all columns
        other_cols = other_cols[other_cols != col_to_keep]  # Exclude col_to_keep
        original_X = X[:, other_cols].copy()
        X[:, other_cols] = 0.0  
        
        Sim2 = runSim("_to_zero")
        
        fcol = math.floor(math.sqrt(X.shape[1]))
        frw = math.ceil(X.shape[1]/fcol)
        
        ax = plt.subplot(frw,fcol, col_to_keep + 1)  # Select the corresponding subplot
        
        variable_name = x_filen_list[col_to_keep].replace('.nc', '')
        ax.set_title(variable_name)
        
        num_bins = 10
        hist, bin_edges = np.histogram(X[:, col_to_keep], bins=num_bins)
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        median_values = []
        percentile_10 = []
        percentile_90 = []
        
        for i in range(num_bins):
        
            mask = (X[:, col_to_keep] >= bin_edges[i]) & (X[:, col_to_keep] < bin_edges[i + 1])
            values_in_bin = []
        
            for rw in range(Sim.shape[0]):
                sim_final = non_masked_data(Sim[rw]) - non_masked_data(Sim2[rw])                      
                values_in_bin.append(sim_final[mask])
            values_in_bin = np.array(values_in_bin).flatten()    
            #set_trace()    
            median_values.append(np.median(values_in_bin))
            percentile_10.append(np.percentile(values_in_bin, 10))
            percentile_90.append(np.percentile(values_in_bin, 90))
      
              
        set_trace()   
        ax.plot(bin_centers, median_values, marker='.', label='Median')
        ax.fill_between(bin_centers, percentile_10, percentile_90, alpha=0.3, label='10th-90th Percentiles')                           
                
        X[:, other_cols] = original_X
    fig_dir = combine_path_and_make_dir(dir_outputs, '/figs/')
    
    plt.savefig(fig_dir + 'control-response-curves.png')   


def sensitivity_reponse_curve(Sim, X, eg_cube, id_name, dir_samples, 
                            run_name, grab_old_trace):
    #plot sensitivity response curves
    plt.figure(figsize=(14, 12))
    for col in range(X.shape[1]-1):
        x_copy = X[:, col].copy()  # Copy the values of the current column
        
        print(col)

        dx = 0.001
        X[:, col] -= dx/2.0  # Subtract 0.1 of all values for the current column

        Sim3 = runSim("subtract_01")    
        
        X[:, col] = x_copy #restore values
        
        X[:, col] += dx/2.0 #add 0.1 to all values for the current column
        
        Sim4 = runSim("add_01") 
        
        fcol = math.floor(math.sqrt(X.shape[1]))
        frw = math.ceil(X.shape[1]/fcol)
        
        ax = plt.subplot(frw,fcol, col + 1)
        variable_name = x_filen_list[col].replace('.nc', '')
        ax.set_title(variable_name)

         num_bins = 20
        hist, bin_edges = np.histogram(X[:, col], bins=num_bins)
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        median_values = []
        percentile_10 = []
        percentile_90 = []
        
        for i in range(num_bins):
        
            mask = (X[:, col] >= bin_edges[i]) & (X[:, col] < bin_edges[i + 1])
            values_in_bin = []
        
            for rw in range(Sim.shape[0]):
                sim_final = (non_masked_data(Sim3[rw]) - non_masked_data(Sim4[rw]))/dx                      
                values_in_bin.append(sim_final[mask])
            values_in_bin = np.array(values_in_bin).flatten()   
            median_values.append(np.median(values_in_bin))
            try:
                percentile_10.append(np.percentile(values_in_bin, 10))
                percentile_90.append(np.percentile(values_in_bin, 90))
            except:
                percentile_10.append(np.nan)
                percentile_90.append(np.nan)
            
              
        ax.plot(bin_centers, median_values, marker='.', label='Median')
        ax.fill_between(bin_centers, percentile_10, percentile_90, alpha=0.3, label='10th-90th Percentiles')      
        
        #for rw in range(Sim.shape[0]):
            #ax.plot(X[:, col], non_masked_data(Sim3[rw]) - non_masked_data(Sim4[rw]), '.', color = "darkred", markersize = 0.5, linewidth=0.5)
    
        X[:, col] = x_copy 
        #set_trace() 
