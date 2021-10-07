openJulesVar <- function(dir, levels = c(1:5, 12:13), varname = 'landCoverFrac', splitPFTs = FALSE, annual_average = TRUE,
                          yrs = c(2000:2014)) {

	files = list.files(dir, full.names=TRUE)
	yrsf = sapply(files, function(i) strsplit(i, 'Monthly.')[[1]][2])
	if (all(is.na(yrsf))) yrsf = sapply(files, function(i) strsplit(i, 'Annual.')[[1]][2])
	yrsf = as.numeric(unlist(strsplit(yrsf, '.nc')))
    
	index =  apply(sapply(yrs, '==', yrsf), 1, any)
	
	files = files[index]
    openFun <- function(level) {
        TC = layer.apply(files, process.jules.file, level, varname)
        if (annual_average) TC = mean(TC)
        return(TC)
    }
    if (splitPFTs)  TC = layer.apply(levels, openFun)
        else TC = openFun(levels)
        
    TC = convert_pacific_centric_2_regular(TC)
    TC = raster::crop(TC, extent(c(-180, 180, -30, 30)))    
    return(TC)
}
