library(terra)

files = list.files("../data/data/VIIRS-Land/", full.names=TRUE)

variable = 'LAI'

#dat = rast(files, variable)
example = rast('../data/data/raw_bioclim/cvl_C.nc')
extent = c(-141.0, -53.0, 42, 83)
area_name = 'Canada'
#extent = c(-13.0, 3.0, 49, 62)
#area_name = 'UK'
variable = 'LAI'
output_dir = '../data/data/'

example = crop(example, extent)

date  = sapply(files, function(i) strsplit(i, '-NPP_')[[1]][2])
date  = substr(date,1, 8)
dates <- as.Date(date, format = "%Y%m%d")

# Generate a sequence of dates covering the range
min_date <- min(dates)
max_date <- max(dates)
all_dates <- seq.Date(min_date, max_date, by = "day")#[1:(365*2)]
all_dates = gsub('-', '', all_dates)

years  = as.numeric(substr(all_dates, 1, 4))
months = as.numeric(substr(all_dates, 5, 6))
days   = as.numeric(substr(all_dates, 7, 8))

print(area_name)
print(extent)

# Identify missing dates
spread_data_year <- function(year, name = '') {
    tfile = paste(c('../temp/make_lai/_spread_latlon', variable, extent, name, year, '.nc'),
                  collapse = '-')
    
    if (file.exists(tfile)) {
        #out = rast(tfile)
        #if (nlyr(out) != length(year_dt) || year == "2016") browser()

        return(rast(tfile))
    }
    print(tfile)
    year_dt = all_dates[substr(all_dates, 1, 4) == year]
    dat = lapply(year_dt, spread_data)
    dat = rast(dat)
    writeCDF(dat, tfile, overwrite=TRUE)
    gc()
    return(rast(tfile))
}

spread_data <- function(date) {
    print(date)
    fileID = grep(date, files)
    if (length(fileID) > 1) {
        id = which(substr(sapply(files[fileID], function(file) 
                            strsplit(file, 'NPP_')[[1]][2]), 1, 8) == date)[1]
        fileID = fileID[id]
    }
    
    if (length(fileID) == 0) {
        dat = example
        dat[] = NaN
    } else {
        file = files[fileID]
        dat = rast(file, variable)
        dat = crop(dat, extent)
        dat = resample(dat, example, method = "bilinear")
        
        #dat = rdat[[layer]]
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
        
    }
    if (nlyr(dat)>1) browser()
    return(dat)
}

sdat = rast(lapply(unique(years), spread_data_year))

minDat = round(3+ nlyr(sdat) *0.1)

fit_loss <- function(x, y, nx) {
    x = c(x, x + 1, x + 2)
    y = c(y, y, y)
    
    #nx = seq(1, 2, length.out = 365.24)
    predictions = predict(loess(y ~ x, span = 1/6), newdata = data.frame(x = nx))
  return(predictions)
}

gapfill <- function(i) {
    tfile = paste(c('../temp/make_lai/_gfill/xy-', xyFromCell(sdat, i), date[1], tail(date, 1),
                    '.csv'), collapse = '-')
    
    if (file.exists(tfile)) return(read.csv(tfile, stringsAsFactors=FALSE)[,1])
    #if (i == 100*round(i/100)) print(tfile)
    vs = svals[i,]
    if (all(is.na(vs))) {
        out = vs
    } else if (sum(!is.na(vs)) < minDat) {
        out = vs
    } else {
        vs0 = vs
        vs = log((vs)^2 + 0.000001)
        
        
    
        # Interpolate climatology
        x = ((months-1) + days/31)/12
        vs_clim = vs[sort.int(x, index.return=T)[[2]]]
        x_clim = sort(x)
        mask_clim = !is.na(vs_clim)
        clim = fit_loss(x = x_clim[mask_clim], y = vs_clim[mask_clim], x + 1)
        
        avs = vs-clim
    
        # Generate indices of non-NaN values
        non_nan_indices <- !is.nan(vs)
        
        # Create a sequence of indices for the entire vector
        all_indices <- seq_along(vs)
        
        
        #
        interpolated_values <- approx(x = all_indices[non_nan_indices], 
                                      y = avs[non_nan_indices], 
                                      xout = all_indices, method = "linear", rule = 2)$y
    
        # Replace NaN values with interpolated values
        vs[is.nan(vs)] <- interpolated_values[is.nan(vs)] + clim[is.nan(vs)]
        vs = sqrt(exp(vs))
        out = vs
    }   
    write.csv(out, tfile, row.names = FALSE) 
    return(vs)
}

gapfill_row <- function(row) {
    tfile = paste(c('../temp/make_lai/_gfill/row-', row, date[1], tail(date, 1), extent,
                    '.csv'), collapse = '-')
    print(tfile)
    print(row/nrow(sdat))
    if (file.exists(tfile)) return(as.matrix(read.csv(tfile, stringsAsFactors=FALSE)))
    index = seq(1 + (row-1) * ncol(sdat), length.out = ncol(sdat))
    out = sapply(index, gapfill)

    write.csv(out, tfile, row.names = FALSE)
    
    return(out)
}
tfile = paste(c('../temp/make_lai/_gfill/row-', date[1], tail(date, 1), extent,
                    '.csv'), collapse = '-')

if (file.exists(tfile)) {
    outs = read.csv(tfile, stringsAsFactors=FALSE)
} else {
    svals = values(sdat)
    outs = lapply(1:nrow(sdat), gapfill_row)
    outs = t(do.call(cbind, outs))
    write.csv(outs, tfile, row.names = FALSE)
}
print("all the hard stuff done")
#values(sdat) = t(outs)
#sdat = rast(lapply(1:nlyr(rdat), spread_data), 'double')

#out_file = paste0(output_dir, '/', area_name, '-', variable, '.nc')
#writeCDF(sdat, out_file, overwrite = TRUE)
browser()

