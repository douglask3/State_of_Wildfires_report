library(raster)
source("../gitProjectExtras/gitBasedProjects/R/sourceAllLibs.r")
sourceAllLibs("../rasterextrafuns/rasterPlotFunctions/R/")
sourceAllLibs("../rasterextrafuns/rasterExtras/R/")
sourceAllLibs("../gitProjectExtras/gitBasedProjects/R/")
library(ncdf4)
sourceAllLibs("../ConFIRE_attribute/libs/")
source("../ConFIRE_attribute/libs/plotStandardMap.r")
source("../LPX_equil/libs/legendColBar.r")
source("libs/find_levels.r")
graphics.off()

#########
## cfg ##
#########

colss = list(c('#ffffe5','#f7fcb9','#d9f0a3','#addd8e','#78c679','#41ab5d','#238443','#006837','#004529'), 
            c('#ffffd9','#edf8b1','#c7e9b4','#7fcdbb','#41b6c4','#1d91c0','#225ea8','#253494','#081d58'),
            c('#fff7f3','#fde0dd','#fcc5c0','#fa9fb5','#f768a1','#dd3497','#ae017e','#7a0177','#49006a'),
            c('#ffffff','#f0f0f0','#d9d9d9','#bdbdbd','#969696','#737373','#525252','#252525','#000000'))

runs = c('factual-', 'counterfactual-', 'early_industrial-')
dir = 'outputs/ConFire_Greece-tuning2/samples/_13-frac_points_0.8/'

controls = c('Fuel', 'Moisture', 'Ignitions', 'Suppression')

plot_control <- function(i, name, cols) {
    plot_run <- function(run, levels = NULL) {
        
        files = list.files(paste0(dir, '/', run, '/Standard_', i-1), full.names = TRUE)
        files = files[grepl('sample-pred', files)]
    
        dats = layer.apply(files, function(i) mean(brick(i)))
        dats[dats>9E9] = NaN     
        quants = apply(dats[], 1, quantile, c(0.1, 0.9), na.rm = TRUE)
        mn = mx = dats[[1]]
        mn[] = quants[1,]
        mx[] = quants[2,]
        
        if (is.null(levels))
            levels = find_levels(c(mn[!is.na(mn)], mx[!is.na(mx)]), seq(10, 90, 10))
        plotStandardMap(mn, cols = cols, limits = levels)
        if (i == 1) mtext(side = 2, '10%', adj = 0.15, line = 0)
        if(run == runs[[1]]) mtext(name, adj = -0.4, xpd = TRUE, line = 0.5)
        plotStandardMap(mx, cols = cols, limits = levels)
        if (i == 1) mtext(side = 2, '90%', adj = 0.85, line = 0)
        if (run == tail(runs, 1)) 
            legendColBar(c(0.1, 0.7), c(0.1, 0.8), cols = cols, limits = levels, 
                           extend_min = F, minLab = 0, transpose = TRUE)
        return(list(addLayer(mn, mx), levels))
    }
    quants = plot_run(runs[1])
    lapply(runs[-1], plot_run, quants[[2]])
}
#png("r-code/Canada_ConFire_histMaps.png", height = 8, width = 7.2, res = 300, units = 'in')
#layout(rbind(c(1:3, 7:9), c(4:6, 10:12)), widths = c(1, 1, 0.3, 1, 1, 0.3))

layout(matrix(1:(4*(length(runs)*2 + 1)), ncol = 4))
par(mar = c(0, 0.0, 0.5, 0), oma = c(0, 1, 2, 1))
out = mapply(plot_control, 1:length(controls), controls, colss)
#dev.off()
