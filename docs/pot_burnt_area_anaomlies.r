library(raster)
source("../gitProjectExtras/gitBasedProjects/R/sourceAllLibs.r")
sourceAllLibs("../rasterextrafuns/rasterPlotFunctions/R/")
sourceAllLibs("../rasterextrafuns/rasterExtras/R/")
sourceAllLibs("../gitProjectExtras/gitBasedProjects/R/")
library(ncdf4)
sourceAllLibs("../ConFIRE_attribute/libs/")
source("../ConFIRE_attribute/libs/plotStandardMap.r")
source("../LPX_equil/libs/legendColBar.r")
source("../Bayesian_fire_models/libs/find_levels.r")
graphics.off()

cols = rev(c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061'))
dir  = 'data/data/burnt_area/'
regions = list('NW_Amazon', 'Greece', 'Canada')
boxess = list(list(list(9, c(-70, -55, -8, 0)),list(10, c(-70, -55, -8, 0))),
             list(list(8, c(24, 27, 40, 42)),list(8, c(22, 25, 37, 39.5))),
             list(list(5, c(-128, -105, 53, 62)), list(6, c(-80, -70, 47, 58)),
                  list(7, c(-80, -70, 47, 58)), list(9, c(-125, -115, 57, 64))))
speedy = FALSE
plot_FUN <- function(region, boxes) {
    eg_extent_file = paste0("data/data/driving_data/", region, 
                            "/nrt/period_2013_2023/burnt_area.nc")
        
    openDat <- function(file) {
        print(file)
        dat = raster(file)
        extent(dat) = c(-180, 180, -90, 90)
        if (file == files[1]) {
            eg_extent_rast = raster(eg_extent_file)
            dat = crop(dat, eg_extent_rast)        
            eg_extent_rast = resample(eg_extent_rast, dat)
            eg_extent_rast = is.na(eg_extent_rast)
            eg_extent_rast <<- eg_extent_rast
        }
        dat[eg_extent_rast] = NaN
        return(dat)
    }
    
    out_file = paste0('data/data/driving_data/', region, '/raw_burnt_area.nc')
    if (file.exists(out_file)) {
        dat = brick(out_file)
    } else {
        files = list.files(dir, full.names = TRUE)
        dat = layer.apply(files, openDat)
        writeRaster(dat, out_file, overwrite = TRUE)
    }

    last_year = dat[[(nlayers(dat)-11):nlayers(dat)]]
    clim = layer.apply(1:12, function(mn) mean(dat[[seq(mn, nlayers(dat), by = 12)]]))
    anom = last_year - clim
    anom = anom / (4*area(anom))

    levels = find_levels_n(anom[[8]], 10, TRUE)
  
    plot_month <- function(r, mn, addX1 = FALSE, addY2 = FALSE, addX3 = FALSE, addY4 = FALSE,
                           boxes, ...) {
        print(month.abb[mn])
        if (is.list(r)) r = r[[1]]
        plotStandardMap(r, cols = cols, limits = levels, speedy = speedy, ...)
        mtext(month.abb[mn], side = 1, line = -1.5, adj = 0.03)
        
        if (addX1) axis(1)
        if (addX3) axis(3)
        if (addY2) axis(2)
        if (addY4) axis(4)
    
        addBox <- function(bx) 
            polygon(bx[[2]][c(1, 2, 2, 1, 1)], bx[[2]][c(3, 3, 4, 4, 3)])
        lapply(boxes[which(sapply(boxes, function(i) i[[1]]) == mn)], addBox)
        
    }

    anom = layer.apply(anom, function(i) list(i))
    hght = nrow(anom[[1]][[1]])/ncol(anom[[1]][[1]])
    lmat = matrix(1:12, ncol = 4)
    lmat = cbind(0, rbind(0, t(lmat), 13, 0), 0)
    heights = c(0.15, rep(hght, 4), 0.5, 0.15)
    widths = c(0.15, rep(1, 3), 0.15)
    
    png(paste0("outputs/figs/mnthly_BA_anaom", region, '.png'), 
        height = 2.5*sum(heights), width = 2.5*sum(widths), res = 300, units = 'in')
        layout(lmat, heights = heights, widths = widths)
        par(mar = rep(0.5, 4), oma = c(0, 0, 0, 0))
        mapply(plot_month, anom, 1:12,
               c(rep(F, 9), rep(T, 3)), rep(c(T, F, F), 4),
               c(rep(T, 3), rep(F, 9)), rep(c(F, F, T), 4), 
               MoreArgs = list(boxes = boxes))
        legendColBar(c(0.5, 0.7), c(0.1, 0.8), cols = cols, limits = levels, 
                     extend_min = T, extend_max = T, transpose = TRUE, oneSideLabels = NA)
    dev.off()
}

mapply(plot_FUN, regions, boxess)
