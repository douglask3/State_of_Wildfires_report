library(terra)

dir = "/net/data/users/dkelley/state_of_fires_report_20YY/data/raw_bioclim"

variables_info = list('Fuel-Moisture-Live' = c('Fuel-Moisture/Live/', 'LFMC'),
                      'd2m' = c('data.nc', 'd2m'),
                      't2m' = c('data.nc', 't2m'),
                      'tp'  = c('data.nc', 'tp'),
                      'Fuel-Moisture-Dead-Foilage' = c('Fuel-Moisture/Dead/', 'DFMC_Foliage'),
                      'Fuel-Moisture-Dead-Wood' = c('Fuel-Moisture/Dead/', 'DFMC_Wood'))

out_dir1 = '../data/data/driving_data/Canada/'
out_dir2 = 'nrt'

eg_rast = rast('../data/data/driving_data/Canada/isimp3a/obsclim/GSWP3-W5E5/period_2002_2019/consec_dry_mean.nc')[[1]]

year_range = c(2003, 2024)
extent = c(-141, -53, 42, 83)
out_dir = paste0(out_dir1, '/', out_dir2)
try(dir.create(out_dir))
open_and_process <- function(var_info, name) {
    file = paste0(dir, '/', var_info[1])
    fname_out = paste0(out_dir, '/', name, '.nc')
    if (file.exists(file)  && !dir.exists(file)) {
       
        if (file.exists(fname_out)) return(rast(fname_out))
        dat = rast(file, var_info[2])
        dat = dat[[sapply(unique(time(dat)), function(i) which(time(dat) == i)[1])]]
        yrs = as.numeric(substr(time(dat), 1, 4))
        dat = dat[[(yrs >= year_range[1]) & (yrs <= year_range[2]) ]]
        dat = crop(dat, extent)
        dat = resample(dat, eg_rast)
        dat[is.na(eg_rast)] = NaN
        
        writeCDF(dat, fname_out, overwrite=TRUE)
        
    } else {
        files = list.files(file, recursive=TRUE, full.names=TRUE)
        files = files[substr(files, nchar(files)-2, nchar(files)) == '.nc']
        openFile <- function(file) {
            print(file)
            tfile = paste0(c('../temp/nrt/', var_info[2], extent, 
                           gsub('/', '', strsplit(file, dir)[[1]][2])), collapse = '')
            
            if (file.exists(tfile)) return(rast(tfile))
            out = rast(file, var_info[2])
            if (nlyr(out) > 1) browser()
            out = crop(out, extent)
            out = resample(out, eg_rast)
            out[is.na(eg_rast)] = NaN
            out = writeCDF(out, tfile, overwrite = TRUE)
        }
        dat = rast(lapply(files, openFile))
        dat = writeCDF(dat, '../temp/nrt/temp_no_time_file.nc', overwrite = TRUE)

        dates = sapply(files, function(file) tail(strsplit(file, '_')[[1]], 2))
        
        years = as.numeric(dates[1,])
        mnths = as.numeric(substr(dates[2,], 1, 2))
        days = rep(15, length(years))
        add_date_to_file('../temp/nrt/temp_no_time_file.nc', fname_out, years, mnths, days, 
                         name)
        
    }
}

mapply(open_and_process, variables_info, names(variables_info))



                 
