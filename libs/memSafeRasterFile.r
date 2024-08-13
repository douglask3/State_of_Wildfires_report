memSafeFile.initialise <- function(filename='',ext='nc') {
    filename=paste(filename,'memSafe',sep="")

    memSafeFile.filename <<- filename
    memSafeFile.Count    <<- 0
    memSafeFile.ext      <<- ext
}

memSafeFile <- function() {
    if (!exists("memSafeFile.Count"))  memSafeFile.initialise()
    memSafeFile.Count <<- memSafeFile.Count+1
    return(paste(memSafeFile.filename,memSafeFile.Count,memSafeFile.ext,sep="."))
}

memSafeFile.remove <- function() {
    list.files.patternPath <- function(path,...) {
        splittd=strsplit(path,'/')[[1]]
        pattern=tail(splittd,1)
        path=paste(head(splittd,-1),collapse='/')
        return(list.files(path=path,pattern=pattern,...))
    }
    files=list.files.patternPath(memSafeFile.filename,full.name=TRUE)
    files=files[grepl(memSafeFile.ext,files)]
    unlink(files)
}

memSafeFile.crop <- function(...)
    crop(...,filename=memSafeFile())
