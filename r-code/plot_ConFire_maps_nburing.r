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

cols = c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026')
dcols = rev(c('#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac'))


levels = c(0, 0.001, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)#find_levels(c(mn[!is.na(mn)], mx[!is.na(mx)]), seq(10, 90, 10))

dlevels = c(-0.1, -0.05, -0.02, -0.01, 0.01, 0.02, 0.05, 0.1)

dir = 'outputs/ConFire_UK/samples/crop_lightn_soilM_trees_csoil_pas_vpd_cveg_precip_tas_rhumid_totalVeg-frac_points_0.5/'

controls = c('factual', 'counterfactual', 'ss126_GFDL', 'ss585_GFDL')

plot_control <- function(name, cols, levels,dats0 = NULL, ...) {
    files = list.files(paste0(dir, name), full.names = TRUE)
    
    dats = layer.apply(files[1:20], function(i) mean(brick(i)[[nlayers(brick(i)) + (-31:0)]]))
    if (!is.null(dats0)) dats = dats - dats0
    dats[dats>9E9] = NaN
    mn = mx = dats[[1]]
    qu = apply(dats[], 1, quantile, c(0.25, 0.75), na.rm = TRUE)
    mn[] = qu[1,]
    mx[] = qu[2,]
    
    plotStandardMap(mn*100, cols = cols, limits = levels)
    if (i == 1) mtext(side = 3, '10%', adj = 0.15, line = -2)
    plotStandardMap(mx*100, cols = cols, limits = levels)
    if (i == 1) mtext(side = 3, '90%', adj = 0.85, line = -2)
    mtext(name, adj = -0.4, xpd = TRUE, line = -1.5)
    legendColBar(c(0.1, 0.7), c(0.1, 0.9), cols = cols, limits = levels, ...)
    return(dats)
}

plot_controls <- function(controls, name) {
    png(paste0("r-code/UK_ConFire_histMaps_fires-", name, ".png"), height = 4, width = 7, 
            res = 300, units = 'in')
        layout(rbind(c(1:3, 7:9), c(4:6, 10:12)), widths = c(1, 1, 0.5, 1, 1, 0.5))
        par(mar = c(0, 1, 0, 0), oma = c(0, 0, 0, 2))
        dats = plot_control(controls[1], cols, levels, extend_min = F, minLab = 0)
        lapply(controls[-1], plot_control, dcols, dlevels,dats, extend_min = T)
    dev.off()
}

plot_controls(c('factual', 'counterfactual', 'ss126_GFDL', 'ss585_GFDL'), "GFDL")
plot_controls(c('factual', 'counterfactual', 'ss126_IPSL', 'ss585_IPSL'), "IPSL")
plot_controls(c('factual', 'counterfactual', 'ss126_MPI' , 'ss585_MPI' ), "MPI")
plot_controls(c('factual', 'counterfactual', 'ss126_MRI' , 'ss585_MRI' ), "MRI")
