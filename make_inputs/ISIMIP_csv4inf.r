library(raster)
source("../gitProjectExtras/gitBasedProjects/R/sourceAllLibs.r")
sourceAllLibs("../rasterextrafuns/rasterExtras/R/")

dir = paste0('isimip3a/driving_data//GSWP3-W5E5/Global/historic_TS_', 
             c('2000_2009', '2010_2019'), '/obsclim/')
temp_dir = 'isimip3a/temp/ISIMIP_csv4ins'
out_file = "isimip3a/Global/inference_data/GSWP5.csv"

files = list.files(dir[1])
files = files[grepl('.nc', files)]

fireObs = brick("../fireMIPbenchmarking/data/benchmarkData/GFED4s_v2.nc")
nl = nlayers(fireObs)

openDat <- function(file) {
    print(file)
    out = layer.apply(paste(dir, file, sep = '/'), brick)
    out[[1:nl]]
}

dats = lapply(files, openDat)
dats = c(raster::crop(fireObs, dats[[1]]), dats)
    
mask = !is.na(sum(layer.apply(dats, function(i) i[[1]])))

extract <- function(dat, name) {

    extractLayer <- function(i) dat[[i]][mask]
    
    print(name)
    tfile = paste(temp_dir, name, '.Rd', sep = '-')
    
    if (file.exists(tfile)) {load(tfile); return(mat)}
    mat = sapply(1:nlayers(dat), extractLayer)
    mat = as.vector(mat)
    
    save(mat, file = tfile)
    return(mat)
}
names = c('fireObs', sapply(files, function(i) substr(i, 1, nchar(i)-3)))

out = mapply(extract, dats, names)
colnames(out) = names
write.csv(out, file = out_file)
