writeRaster.gitInfo <- function(x, filename,
	                            VersionNumber = TRUE, URL = TRUE,
								comment = NULL, ...) {

	x = writeRaster(x, filename, format= "CDF",...)

	filename = filename(x)

	if (filename != "")
		addGitInfo2RasterFile(filename, VersionNumber, URL, comment)
	else {
		filename = c()
		for (i in 1:nlayers(x)) filename = c(filename, filename(x[[i]]))
		filename = unique(filename)

		sapply(filename, addGitInfo2RasterFile, VersionNumber, URL, comment)
	}
	return(x)
}

addGitInfo2RasterFile <- function(filename, VersionNumber, URL, comment) {

	attPutStandard <- function(test, name, val)
		if (test) ncatt_put(nc, 0, name, val)

	nc = nc_open(filename, write = TRUE)
		attPutStandard(VersionNumber, "git Revision Number", gitVersionNumber())
		attPutStandard(URL          , "git Repository URL" , gitRemoteURL())
		#attPutStandard(TRUE         , "git info by"        , paste("gitBasedProjects",  packageVersion('gitBasedProjects')))

		if (!is.null(comment)) {
			if (is.null(names(comment)))
				names(comment) = paste(comment, 1:length(comment))
			mapply(attPutStandard, TRUE, names(comment), comment)
		}
	nc_close(nc)
}
