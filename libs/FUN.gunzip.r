FUN.gunzip <- function(FUN,filename,lon_range=NULL,lat_range=NULL,...) {
	if (exists('ptm')){
		cat("\ntime since start of last call:\n")
		print(proc.time() - ptm)
	}
	ptm <<- proc.time()

	print(filename)
	ext=substr(filename,nchar(filename)-2,nchar(filename))

	if (ext=='.gz' ||  ext=='zip') {
		if(ext=='.gz') filename_unzip=substr(filename,1,nchar(filename)-3)  else
			filename_unzip=substr(filename,1,nchar(filename)-4)

		try(gunzip(filename,destname=filename_unzip))
		out=FUN(filename_unzip,...)
	} else out=FUN(filename,...)

	if (nlayers(out)==1) out=out[[1]]

	if (!is.null(lon_range) && !is.null(lat_range)) {
		nl=nlayers(out)
		if (nl==1) out=out[[1]]
		iby=findidDenominator(nl)
		out=layer.apply(1:(nl/iby),crop.mem.safe,iby,out,extent(lon_range[1],lon_range[2],lat_range[1],lat_range[2]))
	} #else #=out*1

	if (ext=='.gz' ||  ext=='zip') {
	    out=out*1
	    try(gzip(filename_unzip))
	}

	cat("time to open file:\n")
	print(proc.time() - ptm)
	return(out)
}

crop.mem.safe <- function (i,iby,r,...) {
	i=((i-1)*iby+1):(i*iby)
	return(mem.safe.crop(r[[i]],...))
}



raster.gunzip <- function(...) FUN.gunzip (raster,...)
stack.gunzip  <- function(...) FUN.gunzip (stack,...)
brick.gunzip  <- function(...) FUN.gunzip (brick,...)
