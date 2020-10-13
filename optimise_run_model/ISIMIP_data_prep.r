library(raster)
library(rasterExtras)
library(gitBasedProjects)
library(ncdf4)

data_in = 'data/ISIMIP_data/'

models = list.files('data/ISIMIP_data/')

#popDens = '../LimFIRE/data/hyde_land/'
#
#popDens = list.files(popDens, full.names = TRUE)
#popDens = popDens[ grepl('_pop', popDens)]
#popDens = popDens[!grepl('.zip', popDens)]
#popDens = list.files(popDens, full.names = TRUE)
#popDens = popDens[grepl('popd_', popDens)]
#popYrs  = as.numeric(sapply(popDens, function(i) substring(strsplit(i, 'popd_')[[1]][2], 1, 4)))
#
#convert2Raster <- function(file) {
#    dat = read.csv(file)
#    browser()
#}
#popDens = lapply(popDens, convert2Raster)
#popDenYs <- function(yr) {    
#    yri = which(yr == popYrs)
#    if (length(yri) == 0) {
#        diff = abs(yr - popYrs)
#        y1 = which.min(diff)#
#        y2 = which.min(diff[-y1])
#        if (y2 == y1) y2 = y2 + 1
#        y1 = popDenYs(popYrs[y1])
#        y2 = popDenYs(popYrs[y2])     
#    } else {
#        browser()
#    }
#}
#popDens = layer.apply(1995:2004, popDenYs)
fireFile = "../fireMIPbenchmarking/data/GFED4.fBA.r0d5.1995.2013.nc"
fireIndex = 8:115
datIndex =  13:120

burnt_area = brick(fireFile)[[fireIndex]]      
rs_ba = NULL
processModel <- function(model) {
    dir = paste0(data_in, model, '/historic/')
    files =  list.files(dir, full.names = TRUE)
    names =  list.files(dir)
    
    names = unlist(strsplit(names, '.nc', fixed = TRUE))
    dat = sapply(files, function(i) brick(i)[[datIndex]])
    mask = all(layer.apply(dat, function(i) !is.na(i[[1]])))
    if (is.null(rs_ba)) {
        rs_ba = raster::resample(burnt_area, mask)
        rs_ba[!mask] = NaN
        writeRaster.gitInfo(rs_ba, file = "data/ISIMIP_data/burnt_area_GFED4sObs.nc", zname = 'time', overwrite = TRUE)
    }
    rs_ba <<- rs_ba
    dat = c(dat, rs_ba)
    dat = sapply(dat,  function(i) as.vector(i[mask]))     
    colnames(dat) = c(names, "burnt_area")
    
    out_file = paste0("inputs/ISIMIP_inference/", model, ".csv")
    write.csv(dat, file = out_file)
}

lapply(models, processModel)
