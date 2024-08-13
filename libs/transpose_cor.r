transpose_cor <- function(transpose=FALSE,x,y,z=NULL) {
	
	xx=x;yy=y
	if (transpose) {
		xx=y
		yy=x
		if(!is.null(z)) z=t(z)
	}
	if (is.null(z)) return(list(xx,yy))
	
	return(list(xx,yy,z))
}
