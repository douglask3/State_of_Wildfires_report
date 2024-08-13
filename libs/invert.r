invert.raster <- function(a,vals=NULL) {
	if(is.null(vals)) vals=unique.raster(a)
	return(invert_gubbins(a,vals))
}

invert <- function(a,vals=NULL) {
	if(is.null(vals)) vals=unique(a)
	vals=vals[!is.na(vals)]
	return(invert_gubbins(a,vals))
}

invert_gubbins <- function(a,vals) {
	rvals=rev(vals)
	b=a
	for (i in 1:length(vals)) b[a==vals[i]]=rvals[i]
	return(b)
}