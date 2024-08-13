makeDir <- function(fname, ...) if(!file.exists(fname)) dir.create(fname, ...)

makeGlobDir <- function(obj, fname, ...) {
	makeDir(fname, ...)
	assignGlobal(obj, fname)
}
