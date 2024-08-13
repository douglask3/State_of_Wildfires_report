gitRemoteURL <- function(gitLoc = ".") {
	gitLoc   = findProjectPath(gitLoc)
	if (gitLoc == 'Not_git_repo') return('Not_git_repo')
	
	fileName = file.path(gitLoc,'config')
	URL = readChar(fileName, file.info(fileName)$size)
	URL = strsplit(URL,'url = ')[[1]][2]

	return(strsplit(URL,'\n\t')[[1]][1])
}
