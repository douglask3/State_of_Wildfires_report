graphics.off()
source("../gitProjectExtras/gitBasedProjects/R/sourceAllLibs.r")
sourceAllLibs("../gitProjectExtras/gitBasedProjects/R/")
sourceAllLibs("libs/")
library(raster)
sourceAllLibs("../rasterextrafuns/rasterExtras/R/")
dir = "../ConFIRE_ISIMIP/inputs2/"
countries = c("Brazil", "Brazil")#, "Malaysia", "Brazil", Paraguay = 'Paraguay')
extents = list(NULL, c(-52.5, -42.5, -17.5, -2.5))

vars = c("trees", "tas", "soilM_top", "soilM_bottom", "crop", "pas", "totalVeg", "precip", "humid", "cveg", "csoil")
areaAv = rep(TRUE, 11)
annualAv = areaAv#c(TRUE, TRUE, TRUE)
monthAv = areaAv#c(TRUE, TRUE, TRUE)
TSannualAv = areaAv

cols = list(c('#f7fcf5','#e5f5e0','#c7e9c0','#a1d99b','#74c476',
              '#41ab5d','#238b45','#006d2c','#00441b'),
            rev(c('#d73027','#f46d43','#fdae61','#fee090','#ffffbf',
              '#e0f3f8','#abd9e9','#74add1','#4575b4')),
           c('#fff7fb','#ece2f0','#d0d1e6','#a6bddb','#67a9cf','#3690c0','#02818a','#016c59','#014636'),
           c('#fff7fb','#ece2f0','#d0d1e6','#a6bddb','#67a9cf','#3690c0','#02818a','#016c59','#014636'), 
          c('#fff7f3','#fde0dd','#fcc5c0','#fa9fb5','#f768a1','#dd3497','#ae017e','#7a0177','#49006a'),
            c('#fff7f3','#fde0dd','#fcc5c0','#fa9fb5','#f768a1','#dd3497','#ae017e','#7a0177','#49006a'), c('#ffffe5','#f7fcb9','#d9f0a3','#addd8e','#78c679','#41ab5d','#238443','#006837','#004529'), 
         c('#f7fbff','#deebf7','#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#08519c','#08306b'), c("white", "cyan", "black"))
limits = list(c(0, 1, 2, 5, 10, 20, 40, 60, 80),
              c(18, 20, 22, 24, 26, 28, 30, 32),
              c(0, 1, 2, 5, 10, 15, 20, 25, 20),
              c(0, 1, 2, 5, 10, 15, 20, 25, 20),
              c(0, 10, 20, 30, 40, 50),
              c(0, 10, 20, 30, 40, 50),
              c(0, 1, 2, 5, 10, 20, 40, 60, 80),
             c(0, 20, 50, 100, 200, 400, 600, 1000, 1500, 2000),
             c(0, 1, 2, 5, 10, 20),
              c(0, 0.01, 0.1, 1, 10, 100, 1000),
              c(0, 0.01, 0.1, 1, 10, 100, 1000))
scale = c(100, 1, 1, 1, 100, 100, 100, 60*60*365*24, 1000, 1, 1)
shift = c(0, 273.15, 0, 0, 0, 0, 0, 0, 0, 0, 0)

periods = c("historic_TS_short" = 1960, "RCP2.6_TS" = 2006, "RCP6.0_TS" = 2006)
col = c("black", "blue", "red")
mapYrs = list(c(1960, 1970), c(1995, 2005), c(2040, 2050), c(2089, 2099))

lty = c("GFDL-ESM2M" = 1, "HADGEM2-ES" = 2, "IPSL-CM5A-LR" = 3, "MIROC5" = 4)

extend_mins = c(FALSE, TRUE, rep(FALSE, 9))
extend_maxs = c(FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE)
maxLabs     = list(100, NULL, NULL, NULL, 100, 100, 100, NULL, NULL, NULL, NULL)

forCountry <- function(country, extent, variable, areaA, annualA, monthA, TSannualA,      
                       scale, shift, cols, limits,
                       extend_min, extend_max, maxLab) {

    #png(paste0("figs/", country, variable, ".png"), height = 10, width = 7, res = 300, units = 'in')
    mapID = 1+(1:length(mapYrs))
    par(oma = c(0, 1.5, 0, 1.5))
    layout(rbind(1,
                 c( 2,  4,  6,  8),
                 c( 3,  5,  7,  9),
                 c( 0,  0, 10, 12),
                 c(14, 14, 11, 13),
                 c( 0,  0, 11, 13)),
#                 mapID, 
#                 mapID + 1*length(mapYrs), 
#                 mapID + 2*length(mapYrs), 
#                 mapID + 3*length(mapYrs), 
#                 mapID + 4*length(mapYrs), 
#                 mapID + 5*length(mapYrs), 
#                 mapID + 6*length(mapYrs), 7*length(mapYrs) + 2),
                 heights = c(2, 1, 1, 1, 0.5, 0.5))
    
    openPeriod <- function(period, yrs) {
        tfile = paste("temp/plot_ISIMIP-x", period, yrs, country, variable, 
                      areaA, annualA, monthA, sep = '-')
        if (is.null(extent)) tfile = paste(tfile, ".Rd", sep = '-')
            else  tfile = paste(c(tfile, extent, ".Rd"), collapse  = '-')        
        if (file.exists(tfile)) {load(tfile); return(list(yrs, tss, maps))}
        dir = paste(dir, country, period, sep = '/')
        files = list.files(dir, recursive = TRUE)
        openMod <- function(file) {
            dat = brick(paste0(dir, '/', file))
            print(nlayers(dat))
            #if (grepl("MIR", file)) browser()
            if (!is.null(extent)) dat = raster::crop(dat, extent)
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
    if (TSannualA) {        
        runnuingMn <- function(ts, nm = 12)
            sapply(nm:length(ts), function(mn) mean(ts[(mn-nm+1):mn]))#ts[nm:length(ts)]#
        
        tsMap = rbind(tsMap, lapply(tsMap[2,], lapply, runnuingMn))
    }
    
    plot(range(unlist(tsMap[1,])), range(unlist(tsMap[4,]))*scale-shift, type = 'n', 
         xlab = '', ylab = '')
    
    legend('topleft', names(tsMap[[2,1]]), lwd = 2, lty = lty, bty = 'n')
    legend('top', colnames(tsMap), col = col, lwd = 2, bty = 'n')
    title(paste(country, '-', variable))
    addTS <- function(tss) {
        yrs = tss[[1]]
        if (TSannualA) yrs12 = yrs[12:length(yrs)] else yrs12 = yrs
        tss[[2]] = lapply(tss[[2]], function(i) i * scale + shift)
        findAxis <- function(nsig = 1) {
            out = seq(min(unlist(tss[[2]])), max(unlist(tss[[2]])) , length.out = 6)
            out = unique(signif(out, nsig))
            if (length(out) < 4) out = findAxis(nsig+1)
            return(out)
        }
        ylabels = findAxis()
        maxT = max(c(ylabels, unlist(tss[[2]]))) 
        minT = min(c(ylabels, unlist(tss[[2]]))) 
        addMod <- function(ts12, ts, name) {
            #if (TSannualA) {
            #    yrs = yrs[12:length(yrs)]
            #    ts = sapply(12:length(ts), function(mn) mean(ts[(mn-11):mn]))
            #}
            test = yrs >=1996 & yrs <= 2005
            if (sum(test) >=12) {
                ts = ts[test]
                ts = sapply(1:12, function(mn) mean(ts[seq(mn, length(ts), by = 12)]))
            
                yrange =  par("usr")[3:4]
                yrange = yrange[1] + diff(yrange) * c(0.67, 0.97)
            
                xclim = seq(1990, 2011, length.out = 12)
                yclim = yrange[1] + (ts-minT) * diff(yrange) / (maxT-minT)
                lines(xclim, yclim, lty = lty[names(lty) == name])  
                
                lines(range(xclim), rep(yrange[1], 2)) 
                text(xclim, yrange[1] - 0.05*diff(yrange), sapply(month.abb, substr, 1, 1))

                lines(rep(xclim[1] - 0.1 * diff(xclim[1:2]), 2), yrange)
                yat = yrange[1] + (ylabels-minT) * diff(yrange) / (maxT-minT)
                text(xclim[1] - 0.1 * diff(xclim[1:2]), yat , ylabels, adj = 1.5, col = tss[[5]])
            }
            
            lines(yrs12, ts12*scale-shift, lty = lty[names(lty) == name], col = tss[[5]])
        }
        mapply(addMod, tss[[4]], tss[[2]], names(tss[[2]]))
    }

    apply(rbind(tsMap, col), 2, addTS)

    addMaps <- function(maps) {
        
        addMap <- function(map, mY)   {              
            if (!is.null(map)) for (id in 1:2) {
                plotStandardMap(map[[id]], cols = cols, limits = (limits + shift)/scale)
                if (id == 1 && addYrs) mtext(side = 3, paste(mY, collapse = '-'))
            } 
        }
        mapply(addMap, maps[[3]], mapYrs)
        #addYrs <<- FALSE
    }   
    mar = par("mar")
    par(mar = rep(0.5,4))

    addYrs <<- TRUE
    apply(tsMap, 2, addMaps)
        StandardLegend(cols, limits, tsMap[[3,1]][[1]][[1]], extend_max = extend_max, 
                       extend_min = extend_min, maxLab = maxLab)
    mtext(c("max", "min", "max", "min"), side = 4, outer = TRUE, adj = c(0.08, 0.25, 0.42, 0.59), line = -1)
    mtext(side = 2, outer = TRUE, "historic", adj = 0.57, line = -1)
    mtext(side = 2, outer = TRUE, line = -26, adj = c(0.17,0.57), c("RCP2.6", "RCP6.0"))
    par(mar = mar)
    
}

forVar <- function(...) {
    
    mapply(forCountry, countries, extents, MoreArgs = list(...))

}

pdf("figs/CountryVarTSMap-withC.pdf", height = 10*1.5, width = 7*1.5)
mapply(forVar, vars, areaAv, annualAv, monthAv, TSannualAv,      
                       scale, shift, cols, limits, 
                       extend_mins, extend_maxs, maxLabs)
dev.off()
