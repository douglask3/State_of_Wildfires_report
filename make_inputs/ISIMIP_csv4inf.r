graphics.off()
library(raster)
source("../gitProjectExtras/gitBasedProjects/R/sourceAllLibs.r")
sourceAllLibs("../rasterextrafuns/rasterExtras/R/")
sourceAllLibs("libs/")
dir = paste0('isimip3a/driving_data//GSWP3-W5E5/Global/historic_TS_', 
             c('2000_2009', '2010_2019'), '/obsclim/')
temp_dir = 'isimip3a/temp/ISIMIP_csv4ins'
out_file = "isimip3a/driving_data/GSWP3-W5E5/Global/inference_data/GSWP5.csv"

reds9   = c("#FFF7FB", "#F7DEEB", "#EFC6DB", "#E19ECA", "#D66BAE", 
            "#C64292", "#B52171", "#9C0851", "#6B0830")

files = list.files(dir[1])
files = files[grepl('.nc', files)]

fireObs = brick("../fireMIPbenchmarking/data/benchmarkData/GFED4s_v2.nc")[[25:204]]
nl = nlayers(fireObs)

grab_cache = FALSE

openDat <- function(file) {
    print(file)
    out = layer.apply(paste(dir, file, sep = '/'), brick)
    out[[1:nl]]
}

if (!exists("dats") || !exists("isFire") && grab_cache) {
    dats = lapply(files, openDat)
    dats = c(raster::crop(fireObs, dats[[1]]), dats)
    dats[[1]] = raster::resample(dats[[1]], dats[[2]])
    mask = !is.na(sum(layer.apply(dats, function(i) i[[1]])))
    notMask = !mask
    isFire = mask & (dats[[1]]>0)
}

names = c('fireObs', sapply(files, function(i) substr(i, 1, nchar(i)-3)))
convert2Clim <- function(dat) {
    monthMeans <- function(mn) mean(dat[[seq(mn, nlayers(dat), by = 12)]])
    layer.apply(1:12, monthMeans)
}

extractBurntCells <- function(r) 
    unlist(lapply(1:nlayers(r), function(i) r[[i]][isFire[[i]]]))
yFire = logit(extractBurntCells(dats[[1]]))
plotVar <- function(dat, name) {
    datC = convert2Clim(dat)
    mask0 = mean(datC) == 0
    datC[notMask] = NaN
    whichX <- function(FUN) {
        out = FUN(datC)
        out[mask0] = NaN
        out
    }
    wmx = whichX(which.max)
    wmn = whichX(which.min)
    dmx = max(datC)
    dmn = min(datC)
    
    cols = c("white", "red", "black")
    limits = unique(signif(quantile(c(dmn[dmn != 0], dmx[dmx != 0]) , 
                                    seq(0, 1, 0.1), na.rm = TRUE), 1))
    plotD <- function(datM) {
        plotStandardMap(datM, cols, limits)
        StandardLegend(cols, limits, datM, add = TRUE, oneSideLabels= FALSE)
    }
    plotD(dmx); mtext(name, side = 2, line = -1)
    plotD(dmn)
    cols = c("#2520e3", "#0bee27", "#e91f02", "#2520e3")
    limits = seq(0.5, 11.5)
    plotC <- function(dat) {
        plotStandardMap(dat, cols, limits)
        SeasonLegend(limits, cols, add = TRUE, plot_loc = c(0.85, 2, 0.0, 0.8))
    }
    plotC(wmx); plotC(wmn)
    if (name == names[1]) {
        plot.new()
        return()
    } 
    x = extractBurntCells(dat)
    if (max(abs(range(x) - c(0, 1))) < 0.1) x = logit(x)
    cols = densCols(x,yFire, colramp = colorRampPalette(reds9), bandwidth = 0.1)
                #cols = densCols(x,y, nbin = 128* 4, bandwidth = 0.1)
    par(mar = c(2.5, 2.5, 0.2, 0.2))
    plot(yFire~x, pch = 20, col = cols, cex = 0.5, xlab = '', ylab = '')
    par(mar = rep(0,4))
} 

png("figs/ISIMIP3a_vars2.png", height = length(dats) * 2, width = 14, units = 'in', res = 300)
layout(t(matrix(1:(length(dats)*5), nrow = 5)), widths = c(1, 1, 1, 1, 0.67))
par(mar = rep(0, 4)) #
mapply(plotVar, dats, names)
dev.off()

  

extract <- function(dat, name) {

    extractLayer <- function(i) dat[[i]][mask]
    
    print(name)
    tfile = paste(temp_dir, name, '.Rd', sep = '-')
    
    if (file.exists(tfile) && grab_cache) {load(tfile); return(mat)}
    mat = sapply(1:nlayers(dat), extractLayer)
    mat = as.vector(mat)
    
    save(mat, file = tfile)
    return(mat)
}

out = mapply(extract, dats, names)
colnames(out) = names
write.csv(out, file = out_file)
