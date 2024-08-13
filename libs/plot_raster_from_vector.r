plot_raster_from_vector <-function(x,y,z,x_range=c(-180,180),y_range=c(-90,90),
							   	   limits=seq(min(z),max(z),10),cols=rainbow(length(limits)),
								   coastline=NULL,coast.lwd=par("lwd"),
								   add_legend=TRUE,legend.pos='bottomleft',
								   smooth_image=FALSE,smooth_factor=5,...) {



	plot(x=range(x_range,na.rm=TRUE), y=range(y_range,na.rm=TRUE), type='n', axes=FALSE, ann=FALSE)

	if (smooth_image) c(x,y,z):=disagg_xyz(x,y,z,TRUE,smooth_factor)

	z=rasterFromXYZ(cbind(x,y,z))

	plot_raster_map(z,limits=limits,cols=cols,
					coastline=coastline,coast.lwd=coast.lwd,
					add_legend=add_legend,legend.pos=legend.pos,add=TRUE,...)
}

disagg_xyz <- function(x,y,z,rm.na=FALSE,smooth_factor=5) {
	z=rasterFromXYZ(cbind(x,y,z))
	z=disaggregate(z,smooth_factor,method="bilinear")
	xy=xyFromCell(z,1:length(values(z)))
	z=values(z)

	if (rm.na) {
		test=is.na(z)==FALSE
		xy=xy[test,]
		z=z[test]
	}

	return(list(xy[,1],xy[,2],z))
}
