setupProjectStructure <- function(names = paste(dirn,'_dir',sep=""),
					                dirn = c('outputs', 'data', "libs",
									         'temp', 'figs', 'docs'))
	mapply(makeGlobDir, names, paste(dirn, '/', sep = ''))


setupProject <- function(dir = '.', initaliseGit = TRUE, remoteURL = NULL, setupStructure = TRUE, ...) {
    if (dir != '.') {
        makeDir(dir)
        setwd(dir)
    }
    if (initaliseGit) git.init()
    if (!is.null(remoteURL)) git.remote(paste('add origin', remoteURL))
    if (setupStructure) setupProjectStructure(...)
}
