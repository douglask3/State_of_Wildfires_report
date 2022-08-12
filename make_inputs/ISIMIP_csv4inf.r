library(raster)
source("../gitProjectExtras/gitBasedProjects/R/sourceAllLibs.r")
sourceAllLibs("../rasterextrafuns/rasterExtras/R/")

dir = '/data/users/dkelley/ConFIRE_ISIMIP/isimip3_inputs/Global/historic_TS_short/obsclim/'
temp_dir = '/data/users/dkelley/ConFIRE_ISIMIP_temp/ISIMIP_csv4ins'
out_file = "/data/users/dkelley/ConFIRE_ISIMIP/isimip3_inputs/Global/inference_data/GSWP5.csv"

files = list.files(dir)
files = files[grepl('.nc', files)]

fireObs = brick("../fireMIPbenchmarking/data/benchmarkData/GFED4s_v2.nc")
nl = nlayers(fireObs)

openDat <- function(file) {
    print(file)
    out = brick(paste(dir, file, sep = '/'))
    out[[1:nl]]
}

dats = lapply(files, openDat)
dats = c(raster::crop(fireObs, dats[[1]]), dats)
    
mask = !is.na(sum(layer.apply(dats, function(i) i[[1]])))

extract <- function(dat, name) {

    extractLayer <- function(i) dat[[i]][mask]
    
    print(name)
    tfile = paste(temp_dir, '.Rd', sep = '-')
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
