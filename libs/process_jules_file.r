
source("libs/warningFuns.r")

process.jules.file <- function(file, level = NULL, varName, fracWeight = FALSE, mask = NULL) {
    print(file)
    nc = nc_open(file)
    #if (class(nc) == "try-error") return(NULL)
    vars = names(nc$var)
    
    if (all(vars != varName)) {
        noFileWarning(c(), varName)
        browser()
	return(NULL)
    }
	
    getVar <- function(var) {
	var = nc$var[[which(vars == var)]]
	dat = ncvar_get( nc, var)
	return(dat)
    }
    
    dat = getVar(varName)
    lat = getVar("latitude")
    lon = getVar("longitude")
    
    tim = getVar("time_bounds")
    if (!is.null(mask))     {
        extent = extent(mask)
        test = lon >= (extent[1]-0.5) & lon <= (extent[2]+0.5) & lat > (extent[3]-0.25) & lat < (extent[4]+0.25)
        if (length(dim(dat)) == 4) dat = dat[test,,,]
        else if (length(dim(dat)) == 3) dat = dat[test,,]
        else if (length(dim(dat)) == 2) dat = dat[test,]
        else if (length(dim(dat)) == 1) dat = dat[test]
        lat = lat[test]
        lon = lon[test]
    }
    l = length(lat)
    
    multiLayer <- function(mn, leveli = level) {
	mdat = dat[, leveli, mn]
	if (!is.null(dim(mdat)))
            mdat = apply(mdat,1 , sum)
	return(mdat)
    }


    cropMaskRasterisze <- function(lon, lat, r) {
        r = rasterFromXYZ(cbind(lon, lat, r))
        if (!is.null(extent) && class(extent) != "standardGeneric")
            r = raster::crop(r, extent)
        if (!is.null(mask)) r[mask] = NaN
        r
    }
	
    singleLayer <- function(mn) dat[, mn]
	
    monthizeData <- function(mn, FUN, ..., rasterize = TRUE) {
    	r = FUN(mn, ...)
	if (rasterize) r = cropMaskRasterisze(lon, lat, r)
	return(r)
    }
    
    if (length(dim(dat)) == 2 && ncol(dat) == 12)
        r = layer.apply(1:12, monthizeData, singleLayer)
    else if (length(dim(dat)) == 3) {
    	if (varName != "landCoverFrac" && fracWeight == TRUE) {
	    openWeightLayer <- function(fracLevel, varLevel) {
		frac = process.jules.file(file, fracLevel, "landCoverFrac")
		r = layer.apply(1:12, monthizeData, multiLayer, varLevel)
		return(frac * r / 100)
	    }
	    ri = mapply(openWeightLayer, list(1, 2, 3:5, 6:8, 9), 1:5)
	    r = ri[[1]] + ri[[2]] + ri[[3]] + ri[[4]] + ri[[5]]				
	} else  {	
            monAll <- function(...) sapply(1:12, monthizeData, multiLayer, 
                                                ..., rasterize = FALSE)  
            if (is.null(level))  {
                
                r = lapply(1:(dim(dat)[2]), monAll)
                r = lapply(r, function(i) cropMaskRasterisze(lon, lat, i))
            } else r = layer.apply(1:12, monthizeData, multiLayer)
        }
    } else {
        r = cropMaskRasterisze(lon, lat, dat)
        if (!is.null(level)) r = sum(r[[level]])
    }   
    return(r)
}
