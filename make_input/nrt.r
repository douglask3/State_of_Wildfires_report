library(terra)
source("libs/add_date_time.r")
dir = "/net/data/users/dkelley/state_of_fires_report_20YY/data/"
region = 'Canada'
extent = c(-180, 180, -90, 70)
data_file = paste0('raw_bioclim/data-', region, '.nc')
variables_info = list('d2m' = c(data_file, 'd2m'),
                      'snowCover' = c(paste0('raw_bioclim/data-', region, '2.nc'), "snowc"),
                      'burnt_area' = c('driving_data/burnt_area.nc', 'burnt_area'),
                      'cropland' = c('HYDE/cropland.nc', 'cropland'),
                      'grazing_land' = c('HYDE/grazing_land.nc', 'grazing_land'),
                      'pasture' = c('HYDE/pasture.nc', 'pasture'),
                      'population_density' = c('HYDE/population_density.nc', 'population_density'),
                      'rangeland' = c('HYDE/rangeland.nc', 'rangeland'),
                      'rural_population' = c('HYDE/rural_population.nc', 'rural_population'),
                      'total_irrigated' = c('HYDE/total_irrigated.nc', 'total_irrigated'),
                      'urban_population' = c('HYDE/urban_population.nc', 'urban_population'),
                      'VOD' = c('raw_bioclim/VOD/', 'var255'),
                      'Fuel-Moisture-Live' = c('raw_bioclim/Fuel-Moisture/Live/', 'LFMC'),
                      't2m' = c(data_file, 't2m'),
                      'tp'  = c(data_file, 'tp'),
                      'Fuel-Moisture-Dead-Foilage' = c('raw_bioclim/Fuel-Moisture/Dead/', 'DFMC_Foliage'),
                      'Fuel-Moisture-Dead-Wood' = c('raw_bioclim/Fuel-Moisture/Dead/', 'DFMC_Wood'))

out_dir1 = paste0('../data/data/driving_data/', region, '/')
out_dir2 = 'nrt'

eg_rast = rast(paste0('../data/data/driving_data/', region, '/isimp3a/obsclim/GSWP3-W5E5/period_2000_2019/burnt_area-2000-2023.nc'))[[1]]
eg_rast = crop(eg_rast, extent)
year_range = c(2012, 2023)

out_dir = paste0(out_dir1, '/', out_dir2, '/period_', year_range[1], '_', year_range[2], '/')
try(dir.create(out_dir, recursive = TRUE))

resample_not_scale <- function(dat, eg_rast) {
    r0 = range(dat[[1]][], na.rm = TRUE)
    r1 = range(resample(dat[[1]], dat[[1]])[], na.rm = TRUE)
    dat = resample(dat, eg_rast)
    dat = (dat - r1[1]) * diff(r0)/diff(r1) + r0[1]
    
    return(dat)
}
open_and_process <- function(var_info, name) {
    file = paste0(dir, '/', var_info[1])
    fname_out = paste0(out_dir, '/', name, '.nc')
    #if (file.exists(fname_out)) return(fname_out)
    if (name == 'VOD') year_range[1] = year_range[1] - 1
    if (file.exists(file)  && !dir.exists(file)) {
        dat = rast(file, var_info[2])
        dat = dat[[sapply(unique(time(dat)), function(i) which(time(dat) == i)[1])]]
        years = as.numeric(substr(time(dat), 1, 4))
        mnths = as.numeric(substr(time(dat), 6, 7))
        days  = as.numeric(substr(time(dat), 9, 10))
        test_dates = (years >= year_range[1]) & (years <= year_range[2])  & !is.na(years)
        dat = dat[[test_dates]]
        dat = resample_not_scale(dat, eg_rast)
        
    } else {
        files = list.files(file, recursive=TRUE, full.names=TRUE)
        files = files[substr(files, nchar(files)-2, nchar(files)) == '.nc']
        openFile <- function(file) {
            tfile = paste0(c('../temp/nrt/', var_info[2], region, 
                           gsub('/', '', strsplit(file, dir)[[1]][2])), collapse = '')
            
            #if (file.exists(tfile)) return(rast(tfile))
            print(file)
            out = rast(file, var_info[2])
            
            if (nlyr(out) > 1) browser()
            
            out = resample_not_scale(out, eg_rast)
            out[is.na(eg_rast)] = NaN
            if (name == 'VOD') out[is.na(out) & !is.na(eg_rast)] = 0
            
            out = writeCDF(out, tfile, overwrite = TRUE)
        }

        dates = sapply(files, function(file) tail(strsplit(file, '_')[[1]], 2))
        
        if (name == 'VOD') {
            years = as.numeric(substr(dates[2,], 1,4))
            mnths = as.numeric(substr(dates[2,], 5,6))
        } else {
            years = as.numeric(dates[1,])
            mnths = as.numeric(substr(dates[2,], 1, 2))
        }
        days = rep(15, length(years))

        test_dates = (years >= year_range[1]) & (years <= year_range[2]) 
        files = files[test_dates]
        dat = rast(lapply(files, openFile))
    }   
    
    years = years[test_dates]; mnths = mnths[test_dates]; days = days[test_dates]

    dat[is.na(dat)] = mean(dat[], na.rm = TRUE)
    dat[is.na(eg_rast)] = NaN
    #yay = sum(is.na(dat))[]
    #if (length(unique(yay)) !=2) browser()
    #if (max(yay) != 156 && max(yay) !=13) browser()
    #print(fname_out)
    #print(sum(!is.na(dat[])))

    if (name == "VOD") {
        index = 13:nlyr(dat)
        writeCDF(dat[[index]], '../temp/nrt/temp_no_time_file.nc', overwrite = TRUE)
        years = years[index]
        mnths = mnths[index]
        days = days[index]
        add_date_to_file('../temp/nrt/temp_no_time_file.nc', 
                         fname_out, years, mnths, days, name,
                         overwrite_date = TRUE)
        dat12 = rast(lapply( 13:nlyr(dat), function(mn) max(dat[[(mn-11):mn]])))
    
        dat12 = writeCDF(dat12, '../temp/nrt/temp_no_time_file.nc', overwrite = TRUE)
        add_date_to_file('../temp/nrt/temp_no_time_file.nc', 
                         paste0(substr(fname_out, 1, nchar(fname_out) - 3), '-12monthMax.nc'), 
                         years, mnths, days, name,
                         overwrite_date = TRUE)
    #for ( i in 1:nlyr(dat)) {
    #    dat[[i]][is.na(eg_rast)] = NaN
    #    mask = ((is.na(dat[[i]]) & (!is.na(eg_rast))))
    #    if (sum(mask[]) > 0) browser()
    #}
    } else {
        
        writeCDF(dat, '../temp/nrt/temp_no_time_file.nc', overwrite = TRUE)
        add_date_to_file('../temp/nrt/temp_no_time_file.nc', 
                         fname_out, years, mnths, days, name,
                         overwrite_date = TRUE)
    }
    plot(rast(fname_out)[[6]])
    
    return(fname_out)
}

fname_outs = mapply(open_and_process, variables_info, names(variables_info))


spread_clim <- function(name, file, vname = NULL) {
    file_out = paste0(out_dir, '/', name, '.nc')
    file.copy(fname_outs[1], file_out, overwrite = TRUE)
    

    if (file.exists(file) && !dir.exists(file)) {
        dat = rast(file, vname)
    } else {
        dat = rast(list.files(file, full.names=TRUE), vname)
    }
    
    dat_out = rast(file_out)

    dat = resample_not_scale(dat, dat_out[[1]])
    dat[is.na(dat)] = 0
    dat[is.na(eg_rast)] = NaN  
  
    for (mn in 1:12) 
        for (mn_out in seq(mn, nlyr(dat_out), by = 12)) dat_out[[mn_out]] = dat[[mn]]

    time(dat_out) = time(rast(file_out))
    dat_out = writeCDF(dat_out,file_out, overwrite = TRUE)
}

spread_clim('lightn',"/hpc//data/d00/hadea/isimip2b/ancils/lightning/LISOTD_HRMC_V2.3.2015.0p5ancil.nc")

mapply(spread_clim, name = c('LiveFuelFoilage', "LiveFuelWood", 
                             "DeadFuelFoilage", "DeadFuelWood"), 
       file = '../data/data/raw_bioclim/C_FUEL/',
       vname = c("Live_Leaf", "Live_Wood", "Dead_Foliage", "Dead_Wood"))
