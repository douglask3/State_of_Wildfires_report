graphics.off()
source("../gitProjectExtras/gitBasedProjects/R/sourceAllLibs.r")
sourceAllLibs("../gitProjectExtras/gitBasedProjects/R/")
sourceAllLibs("libs/")
library(raster)
sourceAllLibs("../rasterextrafuns/rasterExtras/R/")
dir = "/data/users/dkelley/ConFIRE_ISIMIP/inputs_JULES_UKESM/"
countries = c("Global", "Brazil", "Indonesia")

vars = c("burnt_area", "trees", "totalVeg", "crop", "pas", "cveg", "cs_gb", "npp", "gpp", "soilMtot", "tas", "precip", "humid")
areaAv = rep(TRUE, 14)
annualAv = areaAv#c(TRUE, TRUE, TRUE)
monthAv = areaAv#c(TRUE, TRUE, TRUE)
TSannualAv = areaAv

cols = list(burnt_area = c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c',
                            '#fc4e2a','#e31a1c','#bd0026','#800026'),
            trees      = c('#f7fcf5','#e5f5e0','#c7e9c0','#a1d99b','#74c476',
                           '#41ab5d','#238b45','#006d2c','#00441b'),
            totalVeg   = c('#f7fcfd','#e5f5f9','#ccece6','#99d8c9','#66c2a4',
                           '#41ae76','#238b45','#006d2c','#00441b'),
            crop       = c('#fcfbfd','#efedf5','#dadaeb','#bcbddc','#9e9ac8',
                           '#807dba','#6a51a3','#54278f','#3f007d'),
            pas        = c('#fcfbfd','#efedf5','#dadaeb','#bcbddc','#9e9ac8',
                           '#807dba','#6a51a3','#54278f','#3f007d'),
            cveg       = c('#ffffe5','#f7fcb9','#d9f0a3','#addd8e','#78c679',
                           '#41ab5d','#238443','#006837','#004529'),
            cs_gb      = c('#ffffe5','#fff7bc','#fee391','#fec44f','#fe9929',
                           '#ec7014','#cc4c02','#993404','#662506'),
            npp        = c('#ffffe5','#fff7bc','#fee391','#fec44f',
                           '#fe9929','#ec7014','#cc4c02','#993404','#662506'),
            gpp        = c('#ffffe5','#fff7bc','#fee391','#fec44f',
                           '#fe9929','#ec7014','#cc4c02','#993404','#662506'),
            soilMtot   = c('#ffffd9','#edf8b1','#c7e9b4','#7fcdbb',
                           '#41b6c4','#1d91c0','#225ea8','#253494','#081d58'),
            tas        = c('#d53e4f','#f46d43','#fdae61','#fee08b','#e6f598',
                           '#abdda4','#66c2a5','#3288bd'),
            precip     = c('#f7fbff','#deebf7','#c6dbef','#9ecae1','#6baed6',
                           '#4292c6','#2171b5','#08519c','#08306b'),
            humid      = c('#f7fcf0','#e0f3db','#ccebc5','#a8ddb5','#7bccc4',
                           '#4eb3d3','#2b8cbe','#0868ac','#084081'))

limits = list(burnt_area = c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50), 
              trees      = c(0, 1, 2, 5, 10, 20, 40, 60, 80),
              totalVeg   = c(0, 1, 2, 5, 10, 20, 40, 60, 80),
              crop       = c(0, 10, 20, 30, 40, 50),
              pas        = c(0, 10, 20, 30, 40, 50),
              cveg       = c(0, 0.01, 0.1, 1, 10, 100, 1000),
              cs_gb      = c(0, 0.01, 0.1, 1, 10, 100, 1000),
              npp        = c(0, 0.01, 0.1, 1, 10, 100, 1000),
              gpp        = c(0, 0.01, 0.1, 1, 10, 100, 1000),
              soilMtot   = c(0, 1, 2, 5, 10, 15, 20, 25, 20),
              tas        = c(18, 20, 22, 24, 26, 28, 30, 32),
              precip     = c(0, 20, 50, 100, 200, 400, 600, 1000, 1500, 2000),
              humid      = c(0, 1, 2, 5, 10, 20))

scale = c(100*60*60*365*24, 100, 100, 100, 100, 100, 
          1, 1, 60*60*365*24, 60*60*365*24, 1, 
          60*60*365*24, 1, 1000)
shift = c(rep(0, 10), 273.15, 0, 0)

periods = c("historic_TS_short" = 2000, "ssp1_TS" = 2015, "ssp3_TS" = 2015, "ssp5_TS" = 2015)
col = c("black", "blue", "purple", "red")
mapYrs = list(c(2005, 2014), c(2015, 2024), c(2040, 2050), c(2089, 2099))

extend_mins = c(rep(F, 10), T, F, F)
extend_maxs = c(rep(F, 5), rep(T, 8))
maxLabs     = list(100, 100, 100, 100, 100, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)

lty = c(1)

forCountry <- function(country, variable, areaA, annualA, monthA, TSannualA,      
                       scale, shift, cols, limits,
                       extend_min, extend_max, maxLab) {

    #png(paste0("figs/", country, variable, ".png"), height = 10, width = 7, res = 300, units = 'in')
    #dev.new()
    mapID = 1+(1:length(mapYrs))
    par(oma = c(0, 1.5, 0, 1.5))
    layout(rbind(1,
                 c( 2,  4,  6,  8),
                 c( 3,  5,  7,  9),
                 c( 0,  10, 12, 14),
                 c(0, 11, 13, 15),
                    16),
#                 mapID, 
#                 mapID + 1*length(mapYrs), 
#                 mapID + 2*length(mapYrs), 
#                 mapID + 3*length(mapYrs), 
#                 mapID + 4*length(mapYrs), 
#                 mapID + 5*length(mapYrs), 
#                 mapID + 6*length(mapYrs), 7*length(mapYrs) + 2),
                 heights = c(2, 1, 1, 1, 1, 0.5))
    
    openPeriod <- function(period, yrs) {
        tfile = paste("temp/plot_JULES_UKESM", period, yrs, country, variable, 
                      areaA, annualA, monthA, ".Rd", sep = '-')
        if (file.exists(tfile)) {load(tfile); return(list(yrs, tss, maps))}
        dir = paste(dir, country, period, sep = '/')
        files = list.files(dir, recursive = TRUE)       
        
        openMod <- function(file) {
            dat = brick(paste0(dir, '/', file))
            areaR = raster::area(dat[[1]], na.rm = TRUE)
            ts = dat * areaR
            ts = apply(ts[], 2, sum, na.rm = TRUE)
            ts = ts/sum(areaR[], na.rm = TRUE)

            yrs = yrs + (1:nlayers(dat))/12-1/24
            findMap <- function(mapY) {
                index = which(yrs > mapY[1] & yrs < mapY[2])
                if (length(index) == 0) return(NULL)
                map = sum(dat[[index]])
                if (annualA) map = map /length(index)
                if (!monthA) map = map * 12
                map
            }
            maps = lapply(mapYrs, findMap)
            
            return(list(yrs, ts, maps, strsplit(file, '/', fixed = TRUE)[[1]][1]))
        }
        files = files[grepl(variable, files)]
        
        dat = lapply(files, openMod)
        
        yrs = dat[[1]][[1]]
        tss = lapply(dat, function(i)i[[2]])
        names(tss) = sapply(dat, function(i) i[[4]])
        
        rngMap <- function(i) {
            maps = lapply(dat, function(j) j[[3]][[i]])
            if (all(sapply(maps, is.null))) return(NULL)
            range(layer.apply(maps, function(i) i))
        }
        maps = lapply(1:length(mapYrs), rngMap)
        save(yrs, tss, maps, file = tfile)
        return(list(yrs, tss, maps))
    }
    
    tsMap = mapply( openPeriod, names(periods), periods)
    
    joinStarts <- function(tss) {
        joinStart <- function(h, f) c(tail(h, 12), f)
        mapply(joinStart, tsMap[[2,1]], tss, SIMPLIFY = FALSE)
    }
    
    if (TSannualA) {   
        tsMap[2,-1] = lapply(tsMap[2,-1], joinStarts)     
        runnuingMn <- function(ts, nm = 12)
            sapply(nm:length(ts), function(mn) mean(ts[(mn-nm+1):mn]))
        tsMap[2, ] = lapply(tsMap[2,], lapply, runnuingMn)
    }
    
    plot(c(min(unlist(tsMap[1,])), 2099), range(unlist(tsMap[2,]))*scale-shift, type = 'n', 
         xlab = '', ylab = '', xaxs = 'i')
    
    legend('topleft', names(tsMap[[2,1]]), lwd = 2, lty = lty, bty = 'n')
    legend('top', colnames(tsMap), col = col, lwd = 2, bty = 'n')
    title(paste(country, '-', variable))
    addTS <- function(tss) {
        yrs = tss[[1]]
        #if (TSannualA) browser()#yrs = yrs[12:length(yrs)]
        addMod <- function(ts, name) {
            #if (TSannualA) {
            #    yrs = yrs[12:length(yrs)]
            #    ts = sapply(12:length(ts), function(mn) mean(ts[(mn-11):mn]))
            #}
            if (length(ts) == (length(yrs)+1)) yrs = c(yrs[1] - diff(yrs[1:2]), yrs)

            if (length(yrs) != length(ts)) 
                yrs = yrs[(1+(length(yrs) - length(ts))):length(yrs)]
            lines(yrs, ts*scale-shift, lty = lty[names(lty) == name], col = tss[[4]])
        }
        mapply(addMod, tss[[2]], names(tss[[2]]))
    }

    apply(rbind(tsMap, col), 2, addTS)
   
    addMaps <- function(maps) {
        
        addMap <- function(map, mY)   {              
            if (!is.null(map)) for (id in 1:2) {
                if (id == 2 && nlayers(map) ==1) plot.new()
                else plotStandardMap(convert_pacific_centric_2_regular(map[[id]]), 
                                     cols = cols, limits = (limits + shift)/scale)
                
                if (id == 1 && addYrs) mtext(side = 3, paste(mY, collapse = '-'))
            } 
            
        }
        mapply(addMap, maps[[3]], mapYrs)
        #addYrs <<- FALSE
    }   
    mar = par("mar")
    par(mar = rep(0.5,4))

    addYrs <<- TRUE
    #browser()
    apply( tsMap[,c(1, 2, ncol(tsMap))], 2, addMaps)
        StandardLegend(cols, limits, tsMap[[3,1]][[1]][[1]], extend_max = extend_max, 
                       extend_min = extend_min, maxLab = maxLab)
    mtext(c("max", "min", "max", "min"), side = 4, outer = TRUE, adj = c(0.08, 0.25, 0.42, 0.59), line = -1)
    mtext(side = 2, outer = TRUE, "historic", adj = 0.57, line = -1)
    mtext(side = 2, outer = TRUE, line = -26, adj = c(0.17,0.57), c("RCP2.6", "RCP6.0"))
    par(mar = mar)
    
}

forVar <- function(...) {
    
    lapply(countries, forCountry, ...)

}

pdf("figs/CountryVarTSMap-withC.pdf", height = 10, width = 7)
mapply(forVar, vars, areaAv, annualAv, monthAv, TSannualAv,      
                       scale, shift, cols, limits, 
                       extend_mins, extend_maxs, maxLabs)
dev.off()
