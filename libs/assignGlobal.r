assignGlobal <- function (obj, fname) {
	assign(obj, fname)
	assign(obj, fname, envir = as.environment(1))
}
