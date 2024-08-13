##########################################################################################
## Gubbins: Putting everything together													##
##########################################################################################
plot_3way <- function(x,y,A,B,C,x_range=c(-180,180),y_range=c(-90,90),
	limits=c(0.33,0.5,0.67),cols=c("FF","CC","99","55","11"),
	add_legend=TRUE,smooth_image=FALSE,smooth_factor=5, add = FALSE,...) {


	## If the legend is to be included, it goes on bfore the main plot as there are
	## some white bits which might overlap otherwise


	ncols=length(cols)

	# normalise variables
	mag=A+B+C
	A=A/mag
	B=B/mag
	C=C/mag

	A[mag==0]=0.33
	B[mag==0]=0.33
	C[mag==0]=0.33

	out=rasterFromXYZ(cbind(x,y,A))
	out=addLayer(out,rasterFromXYZ(cbind(x,y,B)),rasterFromXYZ(cbind(x,y,C)))

	## Interpolates the data to be plotted so that low-res stuff looks nicer on the eye
	if (smooth_image) {
		c(x,y,A,B,C):=disagg_xyabc(x,y,A,B,C,smooth_factor)
	}

	c(x,y,A,B,C):=remove_nans(x,y,A,B,C)

	## Coverts data into numbers which match up to the colour in cols for that colour
	## channel based on the limits
	Az=cut_results(A,limits)
	Bz=cut_results(B,limits)
	Cz=cut_results(C,limits)

	## z is the index for the colour in zcols
	z=1:length(Az)
	zcols=paste("#",cols[Az],cols[Bz],cols[Cz],sep="")

	## convert z to a map
	z=rasterFromXYZ(cbind(x,y,z))
	lims = (min.raster(z, na.rm = TRUE):max.raster(z, na.rm = TRUE)-0.5)[-1]

	## plot
	plotFun <- function(add)
	plot_raster_from_raster(z, cols = zcols[sort(unique(z))], limits = lims,
							x_range = x_range, y_range = y_range,
							smooth_image = FALSE, smooth_factor = NULL,
		                    readyCut = TRUE, add_legend = FALSE, add =add, ...)
    plotFun(add)
	if (add_legend) add_raster_3way_legend(cols,limits,...)
	plotFun(TRUE)
	return(out)

}

##########################################################################################
## Some little functions to help out with the general effort							##
##########################################################################################


disagg_xyabc <- function(x,y,A,B,C,smooth_factor) {
	c(nn,nn,A):=disagg_xyz(x,y,A,smooth_factor=smooth_factor)
	c(nn,nn,B):=disagg_xyz(x,y,B,smooth_factor=smooth_factor)
	c(x,y,C):=disagg_xyz(x,y,C,smooth_factor=smooth_factor)

	return(list(x,y,A,B,C))
}

remove_nans <- function (x,y,A,B,C) {
	test=is.na(A+B+C)==FALSE
	A=A[test]
	B=B[test]
	C=C[test]
	x=x[test]
	y=y[test]

	return(list(x,y,A,B,C))
}

##########################################################################################
## The big scary legend function which needs tidying....								##
##########################################################################################

add_raster_3way_legend <- function(cols,limits,x_shift=0.25,y_shift=0,add=TRUE,plot_frac=0.17,
								   text_offset=0.1,rvar='',gvar='',bvar='',
								   text.cex=par("cex"),...) {

	ncols=length(cols)
	nsquares=100

	cal_frac_of_plot_covered <- function(parc,up=0) {
		xlims=par(parc)[1:2]
		xlims=xlims+(xlims[2]-xlims[1])*up
		xlims[2]=xlims[1]+(xlims[2]-xlims[1])*plot_frac
		return(xlims)
	}

	if (add) {
		xlims=cal_frac_of_plot_covered("xaxp",x_shift)
		ylims=cal_frac_of_plot_covered("yaxp",y_shift+plot_frac/3)
	} else {
		xlims=c(0,1)
		ylims=c(0,1)
	}


	x=y=z=matrix(0,nsquares,nsquares)
	zcols=matrix("transparent",nsquares,nsquares)
	p=0
	ii=0
	nsquares=nsquares-1
	for (i in (0:nsquares)/nsquares) {
		ii=ii+1
		y[,ii]=i
		jj=0
		for (j in 2*(0:nsquares)/nsquares) {
			jj=jj+1
			x[jj,ii]=j
			if (j>i && j<(2-i)) {
				p=p+1
				z[jj,ii]=p
				v=i
				u=(j-i)/2
				w=1-v-u

				convert_to_ncols <- function(x) x=round(x*(ncols-1)+1)
				v=convert_to_ncols(v)
				u=convert_to_ncols(u)
				w=convert_to_ncols(w)

				zcols[jj,ii]=paste("#",cols[u],cols[v],cols[w],sep="")
			}
		}
	}

	z=as.vector(z)+1

	y=ylims[1]+as.vector(y)*(ylims[2]-ylims[1])*sqrt(2)
	x=xlims[1]+as.vector(x)*(xlims[2]-xlims[1])
	points(x,y,col=zcols,pch=17,xpd=FALSE)

	xp=c(0,0.5,1,1 ,2,2,-1,-1,0	)*(max(x)-min(x))+min(x)
	yp=c(0,1  ,0,-1,-1,2,2,-1,-1)*(max(y)-min(y))+min(y)
	polygon(xp,yp,col="white",border="white")

	xp=c(-1,2,2,-1)*(max(x)-min(x))+min(x)
	yp=c(0,0,-1,-1)*(max(y)-min(y))+min(y)
	polygon(xp,yp,col="white",border="white")

	xp=c(0,0.5,1,0)*(max(x)-min(x))+min(x)
	yp=c(0,1  ,0,0)*(max(y)-min(y))+min(y)
	lines(xp,yp)

	xp=seq(min(x),max(x),length.out=length(limits)+2)[c(-1,-length(limits)-2)]-min(x)
	yp=seq(min(y),max(y),length.out=length(limits)+2)[c(-1,-length(limits)-2)]-min(y)
	xoffset=(max(x)-min(x))*text_offset
	yoffset=(max(y)-min(y))*text_offset
	text(x=min(x)+xp,y=min(y)-1.3*yoffset,rev(limits),cex=text.cex,srt=90)

	text(x=min(x)+xp/2-0.9*xoffset,y=min(y)+yp+0.9*yoffset,(limits),srt=60-90,cex=text.cex)

	text(x=min(x)+max(xp)/2+xp/2+2*xoffset*0.95,y=min(y)+rev(yp)+2*yoffset*0.95,(limits),
		srt=-60+90,
		cex=text.cex)

	xp=c(0.5,0.25,0.75)*(max(x)-min(x))+min(x)+c(0,-2*xoffset,2*xoffset)
	yp=c(0,0.5,0.5)*(max(y)-min(y))+min(y)+c(-3.2*yoffset,1.5*yoffset,2*yoffset)

	apply_labs <- function(i,xp,yp,var,str)
			text(x=xp[i],y=yp[i],var[i],srt=str[i],xpd=TRUE,cex=text.cex*1.33)

	lapply(1:3,apply_labs,xp,yp,c(bvar,gvar,rvar),c(0,50,-50))
	#text(x=xp,y=yp,c(bvar,gvar,rvar),srt=0,xpd=TRUE,cex=text.cex*1.33)


}
