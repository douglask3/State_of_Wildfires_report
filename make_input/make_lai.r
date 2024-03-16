library(terra)

files = list.files("../data/data/VIIRS-Land/", full.names=TRUE)

variable = 'LAI'

dat = rast(files, variable)
example = rast('../data/data/raw_bioclim/cvl_C.nc')
extent = c(-141.0, -53.0, 42, 83)
area_name = 'Canada'
extent = c(-13.0, 3.0, 49, 62)
area_name = 'UK'


example = crop(example, extent)

date  = sapply(files, function(i) strsplit(i, '-NPP_')[[1]][2])
date  = substr(date,1, 8)
year  = as.numeric(substr(date, 1, 4))
month = as.numeric(substr(date, 5, 6))
day   = as.numeric(substr(date, 7, 8))


replace_nans_with_local_mean <- function(rast) {
  # Define a function to calculate local mean of non-NaN values
  local_mean <- function(x) {
    # Calculate mean of non-NaN values
    mean_val <- mean(x, na.rm = TRUE)
    # Replace NaNs with the mean value
    
    x[is.nan(x)] <- mean_val
    return(x)
  }
  
  # Apply the local_mean function with a sliding window of size 3x3x3
  res <- focal(rast, w = matrix(1, nrow = 3, ncol = 3, byrow = TRUE), fun = local_mean)
  
  return(res)
}




#crop_resample <- function(i) 

cdat = crop(dat, extent)
rdat = resample(cdat, example, method = "bilinear")

spread_data <- function(layer, name = '') {
    tfile = paste(c('../temp/make_lai/_spread_latlon', extent, name, date[layer], '.nc'),
                  collapse = '-')
    if (file.exists(tfile)) return(rast(tfile))
    print(tfile)
    dat = rdat[[layer]]
    vdat = t(as.array(dat)[,,1])
    mask = is.na(vdat)
    vdat[mask] = 0.0
    nr = nrow(vdat)
    nc = ncol(vdat)

    mask = !mask
    smudge = smask = matrix(0, nr, nc)
    
    for (i in 0:2) { for (j in 0:2) {
        smudge[(1+i):nr, (1+j):nc] = smudge[(1+i):nr, (1+j):nc] + vdat[1:(nr-i), 1:(nc-j)]
        smask[(1+i):nr, (1+j):nc] = smask[(1+i):nr, (1+j):nc] + mask[1:(nr-i), 1:(nc-j)]

        
        smudge[1:(nr-i), 1:(nc-j)] = smudge[1:(nr-i), 1:(nc-j)] + vdat[(1+i):nr, (1+j):nc]
        smask[1:(nr-i), 1:(nc-j)] = smask[1:(nr-i), 1:(nc-j)] + mask[(1+i):nr, (1+j):nc]
    }}
    
    dat[as.vector((!mask))] = as.vector(((smudge[]/smask[])[!mask]))
    writeCDF(dat, tfile)
}

sdat = rast(lapply(1:nlyr(rdat), spread_data))
#rdat_no_nans <- replace_nans_with_local_mean(rdat)

minDat = round(3+ nlyr(sdat) *0.1)
gapfill <- function(i) {
    tfile = paste(c('../temp/make_lai/_gfill', xyFromCell(dat, i), date[1], tail(date, 1),
                    '.csv'), collapse = '-')

    if (file.exists(tfile)) read.csv(tfile, stringsAsFactors=FALSE)[,1]
    if (i == 1000*round(i/1000)) print(tfile)
    vs = svals[i,]
    if (all(is.na(vs))) return(vs)
    if (sum(!is.na(vs)) < minDat) return(vs)
    vs0 = vs
    vs = log((vs)^2 + 0.000001)
    # Generate indices of non-NaN values
    non_nan_indices <- !is.nan(vs)
    
    # Create a sequence of indices for the entire vector
    all_indices <- seq_along(vs)

    # Interpolate using spline function
    #interpolated_values = splinefun(x = all_indices[non_nan_indices], y = vs[non_nan_indices], method = "natural")(all_indices)
    interpolated_values <- approx(x = all_indices[non_nan_indices], y = vs[non_nan_indices], 
                                  xout = all_indices, method = "linear", rule = 2)$y
    

    
    # Replace NaN values with interpolated values
    vs[is.nan(vs)] <- interpolated_values[is.nan(vs)]
    vs = sqrt(exp(vs))
    write.csv(vs, tfile, row.names = FALSE) 
    
    return(vs)
}

svals = values(sdat)
values(sdat) = sapply(1:(nrow(sdat)*ncol(sdat)), gapfill)
sdat = rast(lapply(1:nlyr(rdat), spread_data), 'double')
browser()

