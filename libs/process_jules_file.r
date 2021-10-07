process.jules.file <- function(file, level = NULL, varName, fracWeight = FALSE) {
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
	
    l = length(lat)
	
    multiLayer <- function(mn, leveli = level) {
	mdat = dat[, leveli, mn]
	if (!is.null(dim(mdat)))
            mdat = apply(mdat,1 , sum)
	return(mdat)
    }
	
    singleLayer <- function(mn) dat[, mn]
	
    monthizeData <- function(mn, FUN, ...) {
    	mdat = FUN(mn, ...)
	r = rasterFromXYZ(cbind(lon, lat, mdat))
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
            monAll <- function(...) layer.apply(1:12, monthizeData, multiLayer, ...)
            if (is.null(level))  r = lapply(1:(dim(dat)[2]), monAll)
            else r = layer.apply(1:12, monthizeData, multiLayer)
        }
    } else r = sum(rasterFromXYZ(cbind(lon, lat, dat))[[level]])
    
    return(r)
}
