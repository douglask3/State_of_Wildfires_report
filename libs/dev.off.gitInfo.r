dev.off.gitWatermark <- function(...) {
    gitWatermark(...)
    dev.off()
}

gitFullInfo <- function(VersionNumber = TRUE, URL = TRUE,
                         timeAndDate = TRUE, comment = '') {
    txt = c('git Information:')
    if (VersionNumber) txt = c(txt, gitVersionNumber())
    if (URL          ) txt = c(txt, gitRemoteURL())
    if (timeAndDate  ) txt = c(txt, as.character(Sys.time()))
    txt = c(txt, comment)
    txt = paste(txt, collapse = '\n')

    return(txt)
}


gitWatermark <- function(VersionNumber = TRUE, URL = TRUE,
                         timeAndDate = TRUE, comment = '',
                         x = 0.5, y = 0.5,
                         col = 'black',transparency = 0.8, srt = 45, ...) {

    txt = gitFullInfo(VersionNumber, URL, timeAndDate, comment)

    par(fig = c(0, 1, 0, 1), mar = rep(0, 4))
    usr = par('usr')

    x = usr[1] + diff(usr[1:2]) * x
    y = usr[3] + diff(usr[3:4]) * y

    if (names(dev.cur()) == 'X11') col = make.lighter(col, transparency)
    else col = make.transparent(col, transparency)

    try(text(x, y, txt, srt = srt, col = col, xpd = TRUE,...), silent = TRUE)
    test = tail(names(warnings()), 1)
}

make.transparent <- function(col, transparency = 0.5) {
     ## Semitransparent colours
     tmp <- col2rgb(col)/255
     rgb(tmp[1,], tmp[2,], tmp[3,], alpha=1-transparency)
}

make.lighter <- function(color, factor = 0.5) {
    col    = col2rgb(color)
    factor = (255 - col) * factor
    col    = col + factor
    col    = rgb(t(col), maxColorValue=255)
    return(col)
}
