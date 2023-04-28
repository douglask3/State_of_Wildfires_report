library(raster)
source("../gitProjectExtras/gitBasedProjects/R/makeDir.r")
source("../gitProjectExtras/gitBasedProjects/R/sourceAllLibs.r")
sourceAllLibs("../gitProjectExtras/gitBasedProjects/R/")

inDir = 'isimip3a/driving_data/GSWP3-W5E5-20yrs/Brazil/historic_TS_2000_2019/obsclim/'
outDirs = paste0('isimip3a/driving_data/GSWP3-W5E5-20yrs/Brazil/', 
                 c('AllConFire_2000_2009/', 'AllConFire_2010_2019/', 'AllConFire_2010_2016/'))
files = list.files(inDir)

layers_defaults = list(1:120, 121:240, 121:204)
layers_fires = list(37:156, NULL, 157:240)

fire_file = "/home/h01/cburton/GitHub/ISIMIP3a/Observations/GFED4.1s_Burned_Fraction.nc"

forPeriod <- function(outDir, layers_default, layers_fires) {
    makeDir(outDir)
    processFile <- function(file, dir = inDir, 
                            layers = layers_default, template = NULL) {
        print(file)
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
        #browser()
        r = setZ(r, as.Date(nms), 'Date')
        
        if (exists("maskRaster")) r[maskRaster] = NaN
        else maskRaster <<- is.na(r[[1]])
        r = setZ(r, as.Date(nms), 'Date')
        fout = paste0(outDir, tail(strsplit(file, '/')[[1]], 1))
        print(fout)
        
        writeRaster.Standard(r, file = fout)
        return(list(r, nms))
    }
    
    template = lapply(files, processFile)[[1]]
    processFile(fire_file, '', layers_fire, template = template)
    #fireObs = brick()[[25:204]]
}
mapply(forPeriod, outDirs, layers_defaults, layers_fires)
