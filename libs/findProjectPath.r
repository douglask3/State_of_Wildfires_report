findProjectPath <- function(path, ncallMax = 10, warnOnly = TRUE) {
    path0 = path
    if (substr(path, nchar(path) - 3, nchar(path)) != '.git') {
        if (substr(path, nchar(path)  , nchar(path)) == '.' &&
            substr(path, nchar(path)-1, nchar(path)) != '..')
            path = paste(path, 'git', sep = '')
        else if (substr(path, nchar(path), nchar(path))=='/')
            path = paste(path, '.git', sep = '')
        else path = paste(path, '/.git', sep = '')
    }
    
	nc = 0
	while (!file.exists(path) && nc < (ncallMax + 1)) {
		path = file.path("..", path)
		nc   = nc + 1
	}
	if (nc > ncallMax) {
		if (warnOnly) {
			warning(paste('\nno git repo found between', path0, 'and', path, '\n'))
			return('Not_git_repo')
		} else {
			stop(paste('\nno git repo found between', path0, 'and', path, '\n'))
		}
	}
	return(path)
}
