graphics.off()
library(raster)
source("libs/plotStandardMap.r")
source("../gitProjectExtras/gitBasedProjects/R/sourceAllLibs.r")
sourceAllLibs("../rasterextrafuns/rasterPlotFunctions/R/")
sourceAllLibs("../rasterextrafuns/rasterExtras/R/")
sourceAllLibs("libs/")
cols = c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a',
            '#e31a1c','#bd0026','#800026')

likiDir = "../ConFIRE_ISIMIP/outputs/sampled_posterior_ConFire_ISIMIP_solutions/attempt4-full/"
likiFiles = paste0(c("historic", "RCP2.6_2010s", "RCP6.0_2010s"), "/fullPost-10.nc")

ObsFile = "../ConFIRE_ISIMIP/longTermRecord_modis--Max_monthly_burnt_areaJan_2001-Dec_2019.nc"

seamaskFile = "../savanna_fire_feedback_test/data//seamask.nc"

pnts = list("Artic fires 2019" = c(120, 68), "Amazon fires 2019" = c(-53, -10.5),
            "SE Aus 2019/2020" = c(147.5, -35))
pnts = list("Artic fires 2019" = c(120, 68), "Amazon fires 2019" = c(-53, -10.5),
            "SE Aus 2019/2020" = c(151, -30))#, "Coastal Argentina 2011" = c( -60,  -37.5))

obs = raster(ObsFile)/4

openLiki <- function(file) 
    lapply(paste0(likiDirs, '/', file), brick)

likiDirs = list.dirs(likiDir, recursive = FALSE)[c(1, 2)]

likis = lapply(likiFiles, openLiki)

seamask = raster(seamaskFile)

obs = raster::resample(obs, likis[[1]][[1]])
obs[seamask == 0] = NaN
plotStandardMap(obs, limits = c(0, 1, 2, 5, 10, 20, 50, 100), col = cols, y_range = c(-60, 90))
lapply(pnts, function(i) points(i[1], i[2]))

px = names(likis[[1]][[1]])

convertLayer2Num <- function(x) {
    x0 = x
    x = substr(x, 2, nchar(x))
    
    if (substr(x, 1, 1) == '.') x = - as.numeric(substr(x, 2, nchar(x)))
    else x = as.numeric(x)
    #logistic(x)
}
px = sapply(px, convertLayer2Num)
forPnt <- function(pnt, nm) {
    xy = c(colFromX(obs, pnt[1]), rowFromY(obs, pnt[2]))
    
    addPoly <- function(liki, col = NULL, ncall) {
        py = lapply(liki, function(i) i[xy[2], xy[1]])
        
        py = py0 = lapply(py, function(i) i / sum(i))
        py = do.call('*', py)^(1/length(py))
        #lapply(py, function(y) polygon(px, y, col = make.transparent(col, 0.9), border = NA))
        
        c(px, py) := spline(px, py, 1000)#log(py+0.000000001), 1000)
        py[py<0] = 0
        #py = exp(py)
        py = py/max(py)
        if (is.null(col))     
            return(unlist(py))
        else {
            polygon(px, py, col = make.transparent(col, 0.95), border = NA)
            if (is.na(bah)) {
                bah = which((cumsum(py)/sum(py))>0.99)[1]
                bah <<- bah
                lines(px[c(bah, bah)], c(0, py[bah]*1.5), col = col, twd = 2, lty = 2)
                
                text(-5, 0.5, '1% likelihood 2010-2020', adj = 0.0,  length = 0.1, angle = 20)
                #arrows(-5.1, 0.5,px[bah], py[bah])
            } else if(ncall ==1 ) {
                lk = round(100*((sum(py[bah:length(py)])-1)/sum(py)), 1)
                text(-4, 0.4, "2090-2010 likelihood", adj = 0)
                if (col == "blue") {
                    text(-4, 0.33, paste0(lk, '% - RCP2.6'), 
                          adj = 0.0, col = 'blue')
                    #arrows(-4.1, 0.33, -5, 0.05, col = 'blue', length = 0.1, angle = 20)
                } else  {
                    text(-4, 0.26, paste0(lk, '% - RCP6.0'), 
                         adj = 0.0, col = 'red')
                    #arrows(-4.1, 0.26, -4.6, 0.08, col = 'red',  length = 0.1, angle = 20)
                }
                
            }
            index = bah:length(px)
            polygon(px[c(bah, index)], c(0, py[index]), col = make.transparent(col, 0.5), density = 50)
              #browser()  
           # }
        }
        return(py)
    }
    py = mapply(addPoly, likis)
    test = which(apply(py, 1, sum)>0)
    #test = c(max(1, test[1]-1), test, min(length(py), tail(test, 1)+1))
    xr = range(px[test])
    xr = range(px)
    xr[2] = 1.4
    #xr = logit(c(0.1, 30)/100)
    #dev.new()
    plot(xr, c(0, max(unlist(py))),
         type = 'n', axes = FALSE, xlab = '', ylab = '')
    bah <<- NaN
    for (i in 1:8) mapply(addPoly, likis, c("black", "red", "blue"), i)
    #polygon(px, py[,2], col = make.transparent('blue', 0.9), border = 'blue')
    #polygon(px, py[,3], col = make.transparent('red', 0.9), border = 'red')
    #polygon(px, py[,1], col = make.transparent(, 0.9))
    
   
    labels = c(0, signif(logistic(seq(xr[1], xr[2], length.out = 5)), 1), 1)
    at = logit(labels)
    labels = labels *100
    at[1] = min(px); at[length(at)] = max(px)
    axis(1, at = at, labels = labels)
    
    ii = xy[1] + -1:1
    jj = xy[2] + -1:1
    obsBased = FALSE
    if (!obsBased) {
        nmS = sapply(1:nchar(nm), function(i) substr(nm,i,i))
        nmS = paste0(nmS[is.na(as.numeric(nmS))], collapse = '')
       
        mtext(side = 3, line = -1.5, nmS)
        mtext.units(sise = 3, line = -2.5,  paste0(pnt[1], '~DEG~E ', pnt[2], '~DEG~N'))
        return()
    }
    pO = p0 = obs[jj, ii]/12
    pO = logit(mean(pO)/100) - 0.0

    lines(c(pO, pO), c(0, max(py)*0.67), lty = 2)
    pc =  apply(py, 2, function(i) sum(i[px > pO])/sum(i))
    if (pc[1] < 0.01) pc[1] = 0.01
    pc[2:3] = pc[2:3]/pc[1]
    pc[1] = pc[1] * 100
    pc = round(pc, 2)
    if (pc[1] <= 1) {
        pc[2:3] =  round((100/pc[2:3])-1)
        pc =  as.character(pc)       
        pc[1] = '< 1'   
        
        mtext(side = 3, line = -3, adj = 1, col = 'blue', 
              paste0('occurs one in ', pc[2], '  years by 2100 under RCP2.6'))
        mtext(side = 3, line = -4.5, adj = 1, col = 'red',  
              paste0('occurs one in ', pc[3], '  years by 2100 under RCP6.0s'))
    
    } else {
        mtext(side = 3, line = -3, adj = 1, col = 'blue', 
              paste0(pc[2], ' times more likely by 2100 under RCP2.6'))
        mtext(side = 3, line = -4.5, adj = 1, col = 'red',  
              paste0(pc[3], ' times more likely by 2100 under RCP6.0'))
    }
        
    #browser()
    #pc =
    #pc[is.infinite(pc)] = '100+'
    mtext(side = 3, line = -1.5, adj = 1,  paste0(nm, ' at ', pc[1], '% likihood'))

    #browser()
    
}

png("figs/likihoodEventCurves-noArrows.png", height = 6, width = 5, units = 'in', res = 300)
    par(mfrow = c(3, 1), mar = rep(1, 4), oma = c(2, 2, 0, 0))
    mapply(forPnt,pnts, names(pnts))
    mtext('Burnt area (%)', side = 1, line = 2, font = 2)
    mtext('Prob.', side = 2, line = 0, outer = TRUE, font = 2)
dev.off()
