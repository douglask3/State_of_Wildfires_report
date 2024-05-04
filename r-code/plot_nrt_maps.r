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

colss = list(c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026'),
             c('#ffffe5','#f7fcb9','#d9f0a3','#addd8e','#78c679','#41ab5d','#238443','#006837','#004529'), 
             c('#ffffe5','#fff7bc','#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#993404','#662506'),
             c('#fff7f3','#fde0dd','#fcc5c0','#fa9fb5','#f768a1','#dd3497','#ae017e','#7a0177','#49006a'),
             c('#ffffff','#f0f0f0','#d9d9d9','#bdbdbd','#969696','#737373','#525252','#252525','#000000'))

dcolss = list(rev(c('#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac')),
              c('#c51b7d','#de77ae','#f1b6da','#fde0ef','#f7f7f7','#e6f5d0','#b8e186','#7fbc41','#4d9221'), 
              rev(c('#8c510a','#bf812d','#dfc27d','#f6e8c3','#f5f5f5','#c7eae5','#80cdc1','#35978f','#01665e')),
              c('#b35806','#e08214','#fdb863','#fee0b6','#f7f7f7','#d8daeb','#b2abd2','#8073ac','#542788'),
              c('#762a83','#9970ab','#c2a5cf','#e7d4e8','#f7f7f7','#d9f0d3','#a6dba0','#5aae61','#1b7837'))

cols_r = c('#f7f4f9','#e7e1ef','#d4b9da','#c994c7','#df65b0','#e7298a','#ce1256','#980043','#67001f')

dir = 'outputs/ConFire_Canada-nrt7/samples/_26-frac_points_0.2/baseline-/'

runs = c('Burnt area' = 'control', 'Fuel' = 'Standard_0', 'Moisture' = 'Standard_1',    
         'Weather' = 'Standard_2',
         'Ignitions' = 'Standard_3', 'Suppression' = 'Standard_4', 'Snow' = 'Standard_5')

mnths = c(6)


png("figs/nrt-Canada-map-obs.png", height = 3, width = 4, 
        res = 300, unit = 'in')
layout(rbind(1:3, 4:6, 7:9), widths = c(0.1, 1, 0.3))
par(mar = rep(0.2, 4))

obs = brick('data/data/driving_data/Canada/nrt/period_2011_2023/burnt_area.nc')

obs_clim = mean(layer.apply(mnths, function(mn) mean(obs[[seq(mn, nlayers(obs), by = 12)]])))
plot.new()
mtext('annual average', side = 4)
levels = find_levels(obs_clim*100, seq(10, 90, 10), not0 = TRUE)
plotStandardMap(obs_clim * 100, cols = colss[[1]], limits = levels)
legendColBar(c(0.1, 0.7), c(0.1, 0.8), cols = colss[[1]], limits = levels, 
                           extend_min = F, minLab = 0)

obs_anom = mean(layer.apply(mnths, function(mn) obs[[tail(seq(mn, nlayers(obs), 12), 1)]])) - obs_clim
plot.new()
mtext('anomoly', side = 4)
levels = find_levels(obs_anom*100, not0 = TRUE)
plotStandardMap(obs_anom * 100, cols = dcolss[[1]], limits = levels)
legendColBar(c(0.1, 0.7), c(0.1, 0.8), cols = dcolss[[1]], limits = levels, 
                           extend_min = T)

obs_rank = layer.apply(mnths, function(mn) sum(obs[[tail(seq(mn, nlayers(obs), 12), 1)]] > obs[[seq(mn, nlayers(obs), by = 12)]]))
plot.new()
mtext('anomoly', side = 4)
levels = find_levels(obs_rank, seq(10, 90, 10), not0 = TRUE)
plotStandardMap(obs_rank, cols = cols_r, limits = levels)
legendColBar(c(0.1, 0.7), c(0.1, 0.8), cols = cols_r, limits = levels, 
                           extend_min = F, minLab = 0)

dev.off()

plot_run  <- function(run, cols, dcols) {
    files = list.files(paste0(dir, run, '/'), full.names = TRUE)
    files = files[grepl('pred', files)]
    
    openFile <- function(file) {
        tfile = paste0('temp/plot_nrt_map/', gsub('/', '-', file))
       
        if (file.exists(tfile)) return(brick(tfile))
        print(file)
        dat0 =  brick(file)
        clim = layer.apply(1:12, function(mn) mean(dat0[[seq(mn, nlayers(dat0), by = 12)]]))
        dat = dat0[[(nlayers(dat0)-11):(nlayers(dat0))]]
        is_an = dat > clim
        rank =  layer.apply(1:12, function(mn) 
                            sum(dat[[mn]] > dat0[[seq(mn, nlayers(dat0), by = 12)]]))
        
        out = addLayer(dat, clim, dat- clim, is_an, rank)

        writeRaster(out, tfile, overwrite = TRUE)
        return(out)
    }
    dats = lapply(files, openFile)
    find_percentile <- function(i, experiment) {
        
        index = ((i-1)*12+1):(i*12)
        for_quantile <- function(qu) {
            for_mn <- function(mn) {
                print(mn)
                out = calc(layer.apply(dats, function(r) r[[mn]]), 
                           function(x) quantile(x, qu,na.rm = TRUE))
            }
            tfile = paste('temp/plot_nrt_map/', gsub('/', '-', dir), 
                          run, experiment, qu, '.nc', sep = '-')
            if (file.exists(tfile)) return(brick(tfile))
            out = layer.apply(index, for_mn)
            out = writeRaster(out, tfile, overwrite = TRUE)
            return(out)
        }
        out = lapply(c(0.1, 0.5, 0.9), for_quantile)
    }
    mn = mnths[1]
    
    dat = mapply(find_percentile,1:4, c('dat', 'cim', 'anom', 'is_an', 'rank'), 
                 SIMPLIFY = FALSE)

    
    #return(dat)
   
    dat = lapply(dat, function(r) layer.apply(r, function(x) x[[mn]]))
    
    lmat = t(matrix(1:16, nrow = 4))
    lmat = cbind(17:20, lmat, 0)
    lmat = rbind(21:26, lmat, 0)
    png(paste0("figs/nrt-Canada-map-", run, '.png'), height = 7.2, width = 7.3, 
        res = 300, unit = 'in')
    layout(lmat, heights = c(0.1, rep(1, 4), 0.1), widths = c(0.1, 1, 1, 1, 0.3, 0.1))
    par(mar = rep(0.2, 4))
    
    
    levels = find_levels(dat[[2]]*100, seq(10, 90, 10))
    layer.apply(dat[[2]] * 100, plotStandardMap, cols = cols, limits = levels)
    legendColBar(c(0.1, 0.7), c(0.1, 0.8), cols = cols, limits = levels, 
                           extend_min = F, minLab = 0)
    
    levels = find_levels(dat[[3]]*100)
    layer.apply(dat[[3]] * 100, plotStandardMap, cols = dcols, limits = levels)
    legendColBar(c(0.1, 0.7), c(0.1, 0.8), cols = dcols, limits = levels, 
                           extend_min = T)

    
    levels = find_levels(dat[[5]], seq(10, 90, 10), not0 = TRUE)
    layer.apply(dat[[5]], plotStandardMap, cols = cols_r, limits = levels)
    legendColBar(c(0.1, 0.7), c(0.1, 0.8), cols = cols_r, limits = levels, 
                           extend_min = F, minLab = 0)

    plot.new()
    plot.new()

    liki = mean(layer.apply(dats, function(r) r[[37:48]][[mn]]))
    levels = find_levels(liki, seq(10, 90, 10), not0 = TRUE)
    cols_l = c('#ffffff','#f0f0f0','#d9d9d9','#bdbdbd','#969696','#737373','#525252','#252525','#000000')

    plotStandardMap(liki, cols = cols_l, limits =  levels)
    legendColBar(c(0.1, 0.7), c(0.1, 0.8), cols = cols_l, limits = levels, 
                           extend_min = F, minLab = 0)
    
    mtext_pnew <- function(txt) {
        plot.new()        
        mtext(txt, side = 4, line = -1, srt = -90)
    }
    mtext_pnew("annual average burnt area")
    mtext_pnew("2023 anomoly")
    mtext_pnew("Year rank")
    mtext_pnew("likelihood of anomoly")
    mtext_pnew("10%")
    mtext_pnew("50%")
    mtext_pnew("90%")
    
    dev.off()
}

mapply(plot_run, runs, colss, dcolss)
