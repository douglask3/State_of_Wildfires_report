library(raster)
library(rasterExtras)
library(gitBasedProjects)
library(ncdf4)
source("libs/plotStandardMap.r")
source("libs/sd.raster.r")
source("libs/filename.noPath.r")
graphics.off()

dir = 'outputs/sampled_posterior_ConFire_ISIMIP_solutions/attempt2/'

models = c("GFDL-ESM2M", "HADGEM2-ES", "IPSL-CM5A-LR", "MIROC5")

periods = c("historic", "RCP2.6", "RCP6.0")

variable = "burnt_area_mean"

limits =  c(0,1, 2, 5, 10, 20, 40)
dlimits = c(-10, -5, -2, -1, -0.1, 0.1,  1, 2, 5, 10)
cols = c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026')

dcols = rev(c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695'))

sc = 12 * 100

obs = mean(brick("data/ISIMIP_data/burnt_area_GFED4sObs.nc"))


logit <- function(r) {
    r[r < 0.0000001] = 0.0000001
    log(r/(1-r))
}

logistic <- function(r) 
    1/(1+exp(r*(-1)))


obsL = logit(obs)

plotMap <- function(model, period, cols, limits, dcols = NULL, dlimits = NULL, dat0 = NULL,
                    anomolise = NULL) {
    file = paste(dir, model, period, "model_summary.nc", sep = '/')
    dat = brick(file, varname = variable)[[c(5, 50, 95)]]

    if (!is.null(anomolise)) 
        if (is.null(dat0)) datP = obs else {        
            dat = logistic(logit(dat) - anomolise)
            datP = dat
    } else datP = dat
    plotStandardMap(datP * sc, cols = cols,limits = limits)
    if (model == models[1]) mtext(side = 3, period)
    if (is.null(dat0)) {
        mtext(side = 2, adj = 0.2, model)
        plot.new()
    } else {
        plotStandardMap((dat-dat0) * sc, cols = dcols,limits = dlimits)
    }
    return(dat)
}

legendFun <- function(cols, limits, dat) 
    add_raster_legend2(cols, limits, dat = dat,
                           transpose = FALSE, srt = 0, oneSideLabels= TRUE,
                           plot_loc = c(0.1, 0.9, 0.73, 0.8), ylabposScling=0.8, add = FALSE)

plot <- function(fname, anomolise = FALSE) {
    png(paste0("figs/", fname, "anomIs_", anomolise, "_maps.png"),
        height = 10, width = 7.2, units = 'in',res  = 300)
        layout(rbind(matrix(1:(6*length(models)), ncol = 3), (6*length(models)) + c(1, 2, 2)))
        par(mar = rep(0, 4), oma = c(0, 2, 2, 0))
        
        if (anomolise) datA = obs else datA = NULL
        dat0 = lapply(models, plotMap, periods[1], cols, limits, anomolise = datA)
        if (anomolise) {
            datA = lapply(dat0, function(i) logit(i) - logit(obs))
            dat0 = list(obs)
        } else datA = list(NULL)
        for (p in periods[2:3]) 
            datP = mapply(plotMap, models, dat0 = dat0, anomolise = datA,
                    MoreArgs = list(p, cols, limits, dcols, dlimits))
        legendFun( cols,  limits, dat0[[1]])
        legendFun(dcols, dlimits, datP[[1]])
    dev.off.gitWatermark()
}
plot("BA")

limits =  c(0, 0.1, 0.5, 1, 5, 10, 50)
dlimits = c(-10, -5, -1, -0.5, -0.1, -0.05, -0.01, 0.01, 0.05, 0.1, 0.5, 1, 5, 10)

plot("BA", TRUE)

sc = sc * 12 * 100# * 10
cols = c('#ffffe5','#f7fcb9','#d9f0a3','#addd8e','#78c679','#41ab5d','#238443','#006837','#004529')
limits = seq(0, 0.9, 0.2)*100

dcols = c('#40004b','#762a83','#9970ab','#c2a5cf','#e7d4e8','#f7f7f7','#d9f0d3','#a6dba0','#5aae61','#1b7837','#00441b')
dlimits = c(-0.5, -0.2, -0.1, -0.01, -0.001, 0.001, 0.01, 0.1, 0.2, 0.5)*100
variable = "potential_fuel"
plot("fuel")


cols = c('#ffffe5','#fff7bc','#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#993404','#662506')

dcols = rev(c('#543005','#8c510a','#bf812d','#dfc27d','#f6e8c3','#f5f5f5','#c7eae5','#80cdc1','#35978f','#01665e','#003c30'))
variable = "potential_moisture"
plot("moisture")
