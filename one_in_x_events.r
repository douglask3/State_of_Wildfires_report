graphics.off()
library(raster)
source("libs/plotStandardMap.r")
source("../gitProjectExtras/gitBasedProjects/R/sourceAllLibs.r")
sourceAllLibs("../rasterextrafuns/rasterPlotFunctions/R/")
sourceAllLibs("../rasterextrafuns/rasterExtras/R/")
sourceAllLibs("libs/")
cols = c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a',
            '#e31a1c','#bd0026','#800026')
dcols = rev(c('#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4'))
likiDir = "../ConFIRE_ISIMIP/outputs/sampled_posterior_ConFire_ISIMIP_solutions/attempt4-full/"#

RCPs = c("RCP2.6", "RCP6.0")
histName = 'historic'
years = c("2010s", "2020s", "2040s", "2090s")

Xevent = 100
Ninter = 1001

models = list.dirs(likiDir, recursive = FALSE, full.name = FALSE)

#postFile = "fullPost-10.nc"
postFile = "model_summary-10.nc"
percentile = T#F

openDat <- function(model, rcp, yrs)
    brick(paste0(likiDir, model, '/', rcp, '_', yrs, '/', postFile))

openMods <- function(...) 
    lapply(models, openDat, ...)

openYrs <- function(...)
    lapply(RCPs, openMods, ...)

futrs = lapply(years, openYrs)

hist = lapply(paste0(likiDir, models, '/', histName, '/', postFile), brick)

ba = substr(names(hist[[1]]), 2, nchar(names(hist[[1]])))
test = substr(ba, 1, 1) == '.'
ba[test] = paste0('-', substr(ba[test], 2, nchar(ba[test])))
ba = as.numeric(ba)
if (percentile) {
    ba = ba/100
    baApprox = seq(0, 1, length.out = Ninter)
} else {
    baApprox = logistic(approx(1:length(ba), ba, n = Ninter)[[2]])
}
one_in_x_event <- function(r, x = Xevent) {
    x = (1- 1/x)
    FUN <- function(i) {
        if (percentile) {
            id = which(ba == x)
            if (length(id) == 1) out = c(id, ba[id], i[id])
            else browser()
        } else {
            y = approx(1:length(i),i,  n = Ninter)[[2]]
            y = cumsum(y/sum(y))        
            xevent = which(y > x)[1]
            out = c(xevent, y[xevent], baApprox[xevent])
        }
        return(out)
    }    
    mask = sum(r) > 0 & r[[1]] < 9E9
    events = apply(r[mask], 1, FUN)
    
    out = r[[1:3]]
    out[[1]][mask] = events[1,]
    out[[2]][mask] = events[2,]
    out[[3]][mask] = events[3,]
    out[!mask] = NaN
    
    return(out)
}

calHistEvents <- function(r, mod) {
    #mod = tail(strsplit(dir, '/')[[1]], 1)
    temp_file = paste0('temp/oneIn100-', mod, '-', Xevent, '-', Ninter, '-hist-', postFile)
    print(temp_file)
    if (file.exists(temp_file)) return(brick(temp_file))
    event = one_in_x_event(r)
    writeRaster(event, file = temp_file, overwrite = TRUE)
    return(event)
}
events = mapply(calHistEvents, hist, models)

matrix2list <- function(x, side = 1) 
    unlist(apply(x, side, list), recursive=FALSE)


x_event_in <- function(r, event) {
    event0 = event
    x = (1-1/Xevent)
    cummOfEvent <- function(i, j, nf, k) {
        if (percentile) {
            y = approx(c(0, ba, 1),c(0, i, 1), xout = baApprox)[[2]]
           

            id = which.min(abs(y-k))
            index = (id+1):length(y)
            yi = y[index]
            yi[yi> 0.99999999] = 0.99999999
            yi = logit(yi)
            pr = baApprox[index]

            yapprox = seq(min(yi), max(yi),   length.out = Ninter)
            papprox = approx(yi, pr, xout = yapprox)[[2]]

            yapprox = head(yapprox, -1)
            papprox = diff(papprox)

            area = sum(logistic(yapprox)*papprox)/sum(papprox)

            sum(y[index]*baApprox[index])/sum(baApprox[index])
            prob = baApprox[which.min(abs(y-k))] * x/nf
            return(c(prob, area))
        } else {
            y = approx(1:length(i),i,  n = Ninter)[[2]]
            out = sum(y[1:j])/sum(y)        
            
            area = sum(y[j:length(y)]*baApprox[j:length(y)])/sum(y[j:length(y)])
            return(c(x*out/nf, area))
        }
    }
    mask = event[[1]] > 0 & event[[1]] < 9E9 & !is.na(event[[1]])
    out = event[[1:2]]
    vals = mapply(cummOfEvent, matrix2list(r[mask]), 
                       event[[1]][mask], event[[2]][mask], event[[3]][mask])
    
    for (i in 1:2) out[[i]][mask] = vals[i,]
    out[!mask] = NaN
    return(out)
}

calFutureOcc <- function(futr, yrs) {
    RCPOcc <- function(rs, rcp) {
        modOcc <- function(r, event, mod) {
            tfile = paste("temp/event100_prob",Xevent, Ninter,
                           yrs, rcp, mod, postFile, sep = '-')
            print(tfile)
            #if (file.exists(tfile)) return(brick(tfile))
            liki = x_event_in(r, event)
              
            writeRaster(liki, file = tfile, overwrite = TRUE)
        }
        mapply(modOcc, rs, events, models, SIMPLIFY = FALSE)
    }
    mapply(RCPOcc, futr, RCPs, SIMPLIFY = FALSE)
}

histOcc <- function(r, event, mod) {
    tfile = paste("temp/event100_prob",Xevent, Ninter,'hist', mod, postFile, sep = '-')
    if (file.exists(tfile)) return(raster(tfile))
    liki = x_event_in(r, event)
    writeRaster(liki, file = tfile, overwrite = TRUE)
}
base = mapply(histOcc, hist, events, models, SIMPLIFY = FALSE)
out = mapply(calFutureOcc, futrs, years, SIMPLIFY = FALSE)
global4Mod <- function(r) {
    r = (1-r[[1]])*r[[2]]
    ar = raster::area(r, na.rm = TRUE)
    sum.raster(r*ar, na.rm = TRUE)/sum.raster(ar, na.rm = TRUE)
}

#base = sapply(base, global4Mod)

#global4Mods <- function(rs) {
#    out = sapply(rs, global4Mod)
#    paste0(round(range(out)*100, 2), collapse = '-')
#}

glb = sapply(out, lapply, lapply, global4Mod)

normaliseGlb <- function(tots) {
    FUN <- function(x, base) {
        out = (unlist(x)/unlist(base))*(100/Xevent)
        paste0(round(range(out), 2), collapse = '-')
    }
    out = mapply(FUN, tots, glb[,1])
}
normGlb = apply(glb, 2, normaliseGlb)



