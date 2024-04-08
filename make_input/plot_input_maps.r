library(terra)
library(raster)
source("../../gitProjectExtras/gitBasedProjects/R/sourceAllLibs.r")
sourceAllLibs('../ConFire/libs/')
sourceAllLibs('libs/')
source("../../ConFIRE_attribute/libs/legendColBar.r")

graphics.off()
dir = '../data/data/driving_data/'

regions = c('Greece', 'Greece')
mnth = 8

isimip3b = c('historical', 'ssp126', 'ssp370', 'ssp585')
models = c('GFDL-ESM4', 'IPSL-CM6A-LR', 'MPI-ESM1-2-HR', 'MRI-ESM2-0', 'UKESM1-0-LL')
isimp3a = c('obsclim', 'counterclim')


isimip3a_period = c('period_2000_2019')
isimip3b_period = c('period_1994_2014', rep('period_2015_2099', length(isimip3b)-1))

runs = c(paste0('isimp3a/', isimp3a, '/GSWP3-W5E5/', isimip3a_period, '/'),
         paste0('isimp3b/', mapply(function(a, b) paste0(a, '/',  models, '/', b, '/'),
                                        isimip3b, isimip3b_period)))
spread_info4plot <- function(a, b) 
    c(a, rep(b, each = length(models)))


names(runs) = spread_info4plot(isimp3a, isimip3b)

variables = c('Consecutive Dry Days' = "consec_dry_mean.nc", "Dry Days" = "dry_days.nc",
              "Max. Temp" = "tas_max.nc", "Precip" = "pr_mean.nc", "Max. VPD" = "vpd_max.nc")
cols = list('Consecutive Dry Days' = c('#ffffe5','#fff7bc','#fee391','#fec44f',
                                       '#fe9929','#ec7014','#cc4c02','#993404','#662506'),
            'Dry Days' = c('#ffffe5','#fff7bc','#fee391','#fec44f',
                                       '#fe9929','#ec7014','#cc4c02','#993404','#662506'),
             "Max. Temp" = c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c',
                             '#fc4e2a','#e31a1c','#bd0026','#800026'),
             "Precip" = c('#ffffd9','#edf8b1','#c7e9b4','#7fcdbb','#41b6c4',
                          '#1d91c0','#225ea8','#253494','#081d58'),
             "Max. VPD" = c('#fff7f3','#fde0dd','#fcc5c0','#fa9fb5','#f768a1',
                            '#dd3497','#ae017e','#7a0177','#49006a'))

dcols = list('Consecutive Dry Days' = rev(c('#8c510a','#bf812d','#dfc27d','#f6e8c3','#f5f5f5',
                                       '#c7eae5','#80cdc1','#35978f','#01665e')),
             'Dry Days' = rev(c('#8c510a','#bf812d','#dfc27d','#f6e8c3','#f5f5f5',
                                       '#c7eae5','#80cdc1','#35978f','#01665e')),
              "Max. Temp" = rev(c('#d73027','#f46d43','#fdae61','#fee090','#ffffbf',
                                  '#e0f3f8','#abd9e9','#74add1','#4575b4')),
              "Precip" = c('#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7',
                           '#d1e5f0','#92c5de','#4393c3','#2166ac'),
              "Max. VPD" = rev(c('#c51b7d','#de77ae','#f1b6da','#fde0ef','#f7f7f7',
                                  '#e6f5d0','#b8e186','#7fbc41','#4d9221')))
region = regions[1]

find_lims <- function(pdat, nsig = 1) {
    lims = quantile(pdat[], na.rm = TRUE, seq(0, 1, 0.1))
    if (lims[1] < 0 && lims[2] < 0) {
        lims = sort(abs(lims))
        neg = TRUE
    } else { neg = FALSE }
    lims = unique(signif(lims, nsig))
    if (neg) lims = c(-rev(lims), lims)
    if (length(lims) < 5) lims = find_lims(pdat, nsig + 1)
    return(lims)
}

forVariable <- function(variable, title, col, dcol) { 
    openDat <- function(run, name_run, col, dat0 = NULL, lims = NULL) {
        dat = rast(paste0(dir, region, '/',run, variable))
        dat = dat[[seq(mnth, nlyr(dat), by = 12)]]
        if (!is.null(dat0)) dat = dat - dat0
        
        mean_dat = mean(dat)
        range_dat = quantile(dat, c(0.1, 0.9))
        pdat = c(range_dat[[1]], mean_dat, range_dat[[2]])

        if (is.null(lims)) lims = find_lims(pdat)
        for (i in 1:nlyr(pdat)) {
            plotStandardMap(raster(pdat[[i]]), cols = col, limits = lims)
            if (i == 2 && is.null(dat0)) mtext(title, side = 2, line = 0)
            if (i == 1 && variable == variables[1]) mtext(name_run, side = 3, line = 0)
        }
        return(list(dat, lims))
    }
    out = openDat(runs[1], names(runs)[1], col)
    legendColBar(cols = col, limits = out[[2]], oneSideLabels=FALSE)
    dat0 = out[[1]]
    print("yay")
    if (length(runs)>1) dlims = openDat(runs[2], names(runs)[2], dcol, dat0)[[2]]
    if (length(runs)>2) browser()
    legendColBar(cols = dcol, limits = dlims, oneSideLabels=FALSE)
}
runs = runs[grepl('isimp3a', runs)]

lmat = cbind(1:3, 4, 5:(4+(3*(length(runs)-1))))
lmat0 = lmat = cbind(lmat, max(lmat)+1)
for (i in 2:length(variables)) lmat = rbind(lmat, lmat0 + max(lmat))
lmat = cbind(0, lmat, 0)
lmat = rbind(0, lmat, 0)

widths = c(0.3, 1, 0.3, rep(1, ncol(lmat)-5), 0.3, 0.01)
heights = c(0.3, rep(1, nrow(lmat)-2), 0.01)
png(paste0("../outputs/figs/input_maps_vars_for_", region, '.png'), height = ncol(lmat) * 4, 
    width = sum(widths) *2, units = 'in', res = 300)
layout(lmat, widths = widths)
par(mar = c(0, 0,0, 0))

mapply(forVariable, variables, names(variables), cols, dcols)

dev.off()

