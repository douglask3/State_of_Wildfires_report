FUN.raster 		<- function(FUN,a,...) 	FUN(values(a),...)
sum.raster 		<- function(...) 		FUN.raster(sum  ,...)
mean.raster 	<- function(...) 		FUN.raster(mean  ,...)
max.raster 		<- function(...) 		FUN.raster(max  ,...)
min.raster 		<- function(...) 		FUN.raster(min  ,...)
range.raster 	<- function(...) 		FUN.raster(range,...)