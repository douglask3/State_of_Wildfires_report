library(raster)
source("../gitProjectExtras/gitBasedProjects/R/sourceAllLibs.r")
sourceAllLibs("../rasterextrafuns/rasterPlotFunctions/R/")
sourceAllLibs("../rasterextrafuns/rasterExtras/R/")
sourceAllLibs('ConFire/libs/')
#sourceAllLibs('make_inputs/libs/')
source("../ConFIRE_attribute/libs/legendColBar.r")

source("../ConFIRE_attribute/libs/plotStandardMap.r")
source("../Bayesian_fire_models/libs/find_levels.r")
graphics.off()
source("../LPX_equil/libs/legendColBar.r")
dir = 'data/data/driving_data/Canada/nrt/period_2012_2023/'

files = list.files(dir, full.names = TRUE)

mn = 6
cols = c('#f7f7f7','#d9d9d9','#bdbdbd','#969696','#737373','#525252','#252525')
plot_file <- function(file) {  
    dat = brick(file) 
    dat = dat[[seq(mn, nlayers(dat), by = 12)]]
    aa = mean(dat)
    #mx = max(dat)

    levels = find_levels_n(aa, 9, not0 = TRUE)
    
    plotStandardMap(aa, cols =cols, limits = levels)
    mtext(side = 3, line = -1, gsub('.nc', '', tail(strsplit(file, '/')[[1]], 1)))
    legendColBar(add = TRUE, yy = c(-140, -60), xx = c(35, 40), cols, levels, 
                 oneSideLabels = FALSE, transpose = TRUE)
}
png("outputs/figs/nrt_input_maps.png", height = 15, width = 15, res = 300, units = 'in')
par(mfrow = c(7, 6), mar = c(0, 0, 1, 0))
lapply(files, plot_file)
dev.off()
