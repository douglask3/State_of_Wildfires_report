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

#countriesSubL = sort(c("Brazil", "Botswana", "Madaga", "Portugal", "Guinea", "Paraguay", "Ghana", "Russia", "Thailand", "Ivory Coast", "Israel", "Cambodia", "Australia", "Canada", "United States of America"))
countriesSubL = NULL
countriesSubL = c("Netherlands", "France", "Belgium", "Germany", "Luxe")
findCountry <- function(i) {
    id = which (substr(ckey, 1, nchar(i)) == i)
    if (length(id) > 1) id = id[nchar(as.character(ckey[id])) == nchar(i)]
    id
}

openDat <- function(model, rcp, yrs) {
    tfile = paste0('temp/', 'masked-', model, '_', rcp, '_', yrs, '_', postFile)
    
    if (file.exists(tfile)) return(brick(tfile))
    out = brick(paste0(likiDir, model, '/', rcp, '_', yrs, '/', postFile))
    
    out[out > 9E9] = NaN
    out = writeRaster(out, file = tfile, overwrite = TRUE)
    return(out)
}

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
            if (file.exists(tfile)) return(brick(tfile))
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

global4Mod <- function(r, areaNorm = TRUE) {
    if (nlayers(r) == 2) r = (1-r[[1]])*r[[2]]
    r[r==r[1]] = NaN
    ar = raster::area(r, na.rm = TRUE)
    out = out1 = sum.raster(r*ar, na.rm = TRUE)
    if(areaNorm) out = out/sum.raster(ar, na.rm = TRUE)
    return(out)
}

#base = sapply(base, global4Mod)

#global4Mods <- function(rs) {
#    out = sapply(rs, global4Mod)
#    paste0(round(range(out)*100, 2), collapse = '-')
#}

#glb = sapply(out, lapply, lapply, global4Mod)

normaliseGlb <- function(tots) {
    FUN <- function(x, base) {
        out = (unlist(x)/unlist(base))*(100/Xevent)
        paste0(round(range(out), 2), collapse = '-')
    }
    out = mapply(FUN, tots, glb[,1])
}
#normGlb = apply(glb, 2, normaliseGlb)

percentileArea <- function(rs, base = NULL, countries, ctid = NaN, ...) {
    rs = rs[[c(10, 90)]]
    rs[rs >9E9] = NaN
    if (! is.na(ctid)) rs[countries != ctid] = NaN
    
    out = layer.apply(rs, global4Mod, ...)
    out = unlist(out)#paste0(round(100*unlist(out), 2), collapse = '-')
    
    if (!is.null(base)) out = out/base
    return(out)   
}

RCPGlobaAnom <- function(rs, histBA, obsBA, ...) {
    out = mapply(percentileArea, rs, histBA, ...)* obsBA   
    out = paste(round(range(out), 12), collapse='-')
    return(out)   
}

obs = mean(brick('../ConFIRE_ISIMIP/burnt_area_GFED4sObs.nc'))*100*12 

forCountry <- function(i, name = '', areaNorm = FALSE, countries, ...) {
    tfile = paste0('temp/conutryBA-unity-', name,  i, '.Rd') 
    print(tfile)
    if (file.exists(tfile))  
        load(tfile)
    else {
        
        if (!is.na(i) && i > 0) obs[countries != i] = NaN
        else if (grepl('Forest', name)) obs[countries == 0] = NaN
        
        obsBA = global4Mod(obs, areaNorm = areaNorm)
        
        if (!areaNorm) obsBA = obsBA/1000000
        histBA = lapply(hist, percentileArea, ctid = i, 
                        areaNorm = areaNorm, countries = countries, ...)
        futrBA = sapply(futrs, lapply, RCPGlobaAnom, histBA, obsBA, ctid = i, 
                        areaNorm = areaNorm, MoreArgs = list(countries = countries,...))
        
        out = cbind(obsBA, futrBA)
        if (!is.na(i)) RCPs = paste0(ckey[i], '-', RCPs)
        rownames(out) = RCPs
        colnames(out) = c("Hist", years)   
        save(out, file = tfile)
    }
    
    return(out)
}


forMask <- function(countries, maskName) {
    #countries[is.na(countries)] = 0
    id = unique(countries[])#[1:10]
    if (!is.null(countriesSubL)) id = sapply(countriesSubL, findCountry)
    
    BAs = lapply(id, forCountry, paste0('byArea', maskName), areaNorm = TRUE, 
                 countries = countries)       
    BAkms = lapply(id, forCountry, maskName, countries = countries)

    BAsout = do.call(rbind, BAs)        
    scenNotSplit <- function(y) {
        y = strsplit(y, '-')[[1]]
        test = which(grepl('e', y))
        if (length(test)>0) 
            y = sapply(test, function(i) as.numeric(paste0(y[i], '-', y[i+1])))
        else 
            y = as.numeric(y)
        return(y)
    }

    #### find zeros
    allZeros <- function(i, j) {
        allZero <- function(x) {
            isZero <- function(y) {
                if (is.numeric(y)) {
                    return(round(y, 6) == 0 )
                } else {
                    return(all(round(scenNotSplit(y), 6)==0))
                }
            }
            #browser()
            out =  sapply(x, isZero)
            return(any(is.na(out)) | all(out))
        }
        #browser()
        out = c(allZero(i), allZero(j))
        if (any(is.na(out))) out = TRUE
        else out =  all(out)
        out
    }
    test = !mapply(allZeros, BAkms, BAs)
    
    BAkms = BAkms[test]
    BAs = BAs[test]

    ##### find biggest 
    biggest <- function(x) {
        big <- function(y) 
            abs(as.numeric(scenNotSplit(y)[2])/x[[1,1]])
    
        large = sapply(x[,-1], big)
        max(large)
    }

    largest = sapply(BAs, biggest)
    largest[is.na(largest)] = 0
    index = sort.int(largest, index.return=TRUE, decreasing=TRUE)[[2]]
    
    #if (is.null(countriesSubL)) {BAkms = BAkms[index]; BAs = BAs[index]}
    #BAkms = BAkms[1:7]; BAs = BAs[1:7]
    
    cols = c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a',
             '#e31a1c','#bd0026','#800026')
    limits = c(0, 0.01, 0.1, 1, 2, 5, 10, 20, 100)
    cols = make_col_vector(cols, limits = limits)

    dcols = rev(c('#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac'))
    dlimits = c(-100, -5, -4, -2, -1, -0.5, -0.2, -0.1, -0.05, -0.02, -0.01, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 100)
    dcols = make_col_vector(dcols, limits = dlimits)

    thsndConvert <- function(nl) {
        txt1 = floor(nl*1000)
        txt2 = round(nl*1000000 - txt1*1000)
    
        zeroise <- function(i) {
            if (i < 10 ) i = paste0('00', i)
            else if (i < 100) i = paste0( '0', i)
            i
        }
        txt2 = sapply(txt2, zeroise)
    
        out = paste0(paste0(paste0(txt1, ',', txt2),collapse = '-'), 'k~m2~')
        #if (grepl("452", out)) browser()
        return(out)
    }
    addLine <- function(y, BA, KM) {
        y = y*2
        font = 1; lwd = 1
        addCell <- function(x, v1, v2, offset, ht, col, lim, diff = 0) {
            if (x == 1) return()
            else if (x > 1) x = x-1
            x = x + 1
            if (is.numeric(v1)) {
                col = col[cut_results(v1, lim)]
                v1 = round(v1, 2)
                v2 = v2 / 100
                if (v2 > 1)  v2 = paste0(round(v2, 3), 'Mk~m2~')
                else if (v2 * 1000 > 1) v2  = thsndConvert(v2)
                else v2 = paste0(round(v2*1000000, 3), 'k~m2~')            
            } else {            
                nl1 = scenNotSplit(v1)
                nl2 = scenNotSplit(v2)/100
                if (max(nl2) > 1) {
                    nl2 = round(nl2, 3)
                    v2 = paste0(nl2[1], '-', nl2[2], 'Mk~m2~')
                } else if (max(nl2) * 1000 > 1)  {
                    v2 = thsndConvert(nl2)                
                } else {
                    nl2 = round(nl2*1000000, 3)
                    v2 = paste0(nl2[1], '-', nl2[2], 'k~m2~')
                }
                 
                v1 = paste(round(nl1, 2), collapse = '-')
            
                    nl = nl1 - diff
            
                if (!(nl[1] <= 0 && nl[2] >= 0)) {font = 2; lwd = 2}
                nl = cut_results(nl, lim)
                nl = nl[1]:nl[2]
                col = col[nl]
            }
            y = y + offset

            addPoly <- function(sx, col, dx, border = NA) {
                polygon(x + (sx + dx * 0.5*c(-1, -1, 1, 1, -1))*0.95, 
                        y + ht * 0.95*c(-1, 1, 1, -1, -1), 
                        col = make.transparent(col, 0.5), border = border, lwd = lwd)
            }
            if (length(col) == 1) addPoly(0, col, 1)
            else { 
                sx = seq(-0.5+0.5/length(col), 0.5-0.5/length(col), length.out = length(col))   
                mapply(addPoly, sx, col, 1/(length(col)))
            }
            addPoly(0, make.transparent('white', 1), border = tail(col, 1), 1)
        
            text.units(x = x-0.42, y = y - ht * 0.33, v2, adj = 0, font = font)
            text.units(x = x-0.42, y = y + ht * 0.33, paste0('(', v1, ')'), adj = 0, 
                       font = font)
            v2 = gsub('~', '', v2)
            v2 = gsub(',', '', v2)
            
            return(list(v1, v2))
        }
        if (is.na(BA[[1,1]])) return()
        base = addCell(0, BA[[1,1]], KM[[1,1]], 0, 1, cols, limits) 
        addCells <- function(v1, v2, offset) {
            out = mapply(addCell, 1:length(v1), v1, v2, 
                   MoreArgs = list(offset, 0.5, dcols, dlimits, diff = base[[1]]))
            cbind(base, matrix(unlist(out), nrow = 2))
        }
        ft1 = addCells(BA[1,-1], KM[1,-1], -0.5)
        ft2 = addCells(BA[2,-1], KM[2,-1], 0.5)
       
        nms = rownames(BA)
        rcps = substr(nms, nchar(nms)-5, nchar(nms))
        nm = substr(nms[1], 1, nchar(nms[1]) - 7)
        nm = strsplit(nm, ' ')[[1]]
        if (length(nm)> 1)
            nm = paste(c(paste(nm[1:floor(length(nm)/2)], collapse = ' '), 
                         paste(nm[(1+floor(length(nm)/2)):length(nm)], collapse = ' ')), 
                       collapse = '\n')
        if (length(nm) == 0) nm = "Global"
        #if (grepl("Armenia", nm)) browser()
        text(-1.25, y, nm, adj = 0)
        text( -0.1, y + c(-0.33, 0.33), rcps, adj = 0)
        out = rbind(ft1, ft2)
        
        rownames(out) = paste0(rep(paste(nm, '-', rcps), each = 2), ' (', c( '%', 'area'), ')')
        colnames(out) = c("Historic", years[-1])
        return(out)
    }

    pdf(paste0("figs/countryBAchange", maskName, ".pdf"), 
        height = 1.1*length(BAs) + 1, width = 14)
    par(mar = rep(0, 4))
    plot(0, 0, xlim = c(-1.75,5.5), ylim = c(length(BAs)*2+2, 0), type = 'n', 
         axes = FALSE, xlab = '', ylab = '')
    
    out = mapply(addLine, 1:length(BAs), BAs, BAkms, SIMPLIFY = FALSE)
    out = do.call(rbind, out)
    
    write.csv(out, file = paste0("outputs/countryBAchange-", maskName, ".csv"))
    
    text(-1.25, 0, 'Country', font = 2, adj = 0)
    text(-0.1, 0, 'RCP', font = 2, adj = 0)
    text(-0.45+(1:4), 0, c("Historic", years[-1]), font = 2, adj = 0)
    text.units(0.55, 0.45, 'k~m2~', adj = 0)
    text.units(0.55, 0.9, '(%)', adj = 0)
    ### km2 (%) coloured by increase/decrease
    
    dev.off()
}

countriesAll = raster('../ConFIRE_ISIMIP/data/countries.nc')
countriesForest = raster( '../ConFIRE_ISIMIP/data/countries-forest.nc')
ckey = read.csv("../ConFIRE_ISIMIP/data/countries_key.csv")[,2]

forMask(countriesAll, "All")
forMask(countriesForest, "Forest")
