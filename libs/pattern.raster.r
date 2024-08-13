pattern.raster <- function(r,pattern="forward-diagonal",res=4,thick=0.7,mult=TRUE) {

	round.p1.d2 <- function(r,invert=FALSE) {
		r=(r+1)/2
		thick=1-thick
		r[r< thick]=0
		r[r>=thick]=1
		if (invert) r=1-r
		return(r)
	}
	
	a=r
	
	if (pattern=="check.board")
		values(a)=round.p1.d2(cos(res*xFromCell(a,1:length(a)))*cos(res*yFromCell(a,1:length(a))))
	if (pattern=="Diamond")
		values(a)=round.p1.d2(cos(res*xFromCell(a,1:length(a)))+cos(res*yFromCell(a,1:length(a))))
	if (pattern=="Circle") {
		thick=1-0.3*thick
		values(a)=round.p1.d2(cos(res*xFromCell(a,1:length(a)))+cos(res*yFromCell(a,1:length(a))),TRUE)	
	}	
	if (pattern=="romulan")	
		values(a)=round.p1.d2(cos(res*yFromCell(a,1:length(a))*sin(res*xFromCell(a,1:length(a)))))
	if (pattern=="circles-spirals")	
		values(a)=round.p1.d2(cos(res*yFromCell(a,1:length(a))*res*xFromCell(a,1:length(a))))
	if (pattern=="forward-diagonal")	
		values(a)=round.p1.d2(cos(res*yFromCell(a,1:length(a))+res*xFromCell(a,rev(1:length(a)))))
	if (pattern=="backward-diagonal")	
		values(a)=round.p1.d2(cos(res*yFromCell(a,1:length(a))-	res*xFromCell(a,rev(1:length(a)))))
	if (pattern=="horizontal")	
		values(a)=round.p1.d2(cos(res*xFromCell(a,1:length(a))))
	if (pattern=="vertical")	
		values(a)=round.p1.d2(cos(res*yFromCell(a,1:length(a))))
		
	if (mult) a=r*a
	
	return(a)
	
	
}
