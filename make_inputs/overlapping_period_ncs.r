library(raster)

inDir = 'isimip3a/driving_data/GSWP3-W5E5-20yrs/Brazil/historic_TS_2000_2019/obsclim/'
outDir = 'isimip3a/driving_data/GSWP3-W5E5-20yrs/Brazil/AllConFire_2000_2010/'
files = list.files(inDir)

layers_default = 1:120

fire_file = "../fireMIPbenchmarking/data/benchmarkData/GFED4s_v2.nc"

processFile <- function(file, dir = inDir, layers = layers_default, template = NULL) {
    if(substr(file, nchar(file)-2, nchar(file)) != '.nc') return()
    r = brick(paste0(dir, file))[[layers_default]]
   
   #nms =  as.Date(paste(rep(years, each = 12), 1:12, 15, sep = '-'))
    if (is.null(template)) {
        nms = substr(names(r), 2, nchar(names(r)[1]))
        nms = sapply(nms, function(i)
                            paste(strsplit(i, '.', fixed = TRUE)[[1]], collapse = '-'))
        nms = as.Date(nms)
    } else {
        nms = template[[2]]
        r = raster::crop(r, template[[1]][[1]])
        r = raster::resample(r, template[[1]][[1]])
    }
    names(r) = nms
    r = setZ(r, nms, 'Date')
    
    writeRaster(r, file = paste0(outDir, tail(strsplit(file, '/')[[1]], 1)), overwrite = TRUE)
    return(list(r, nms))
}

template = lapply(files, processFile)[[1]]
processFile(fire_file, '', layers_default + 25, template = template)
#fireObs = brick()[[25:204]]


