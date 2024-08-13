interpolateAnnual2Monthly <- function(dat) {
    interpolate <- function(d1, d2) {
        interpolateMonth <- function(mn)
            (d1 * (12-mn) + d2 * mn)/12
        
        dat = layer.apply(1:12, interpolateMonth)
    }
    dat = mapply(interpolate, head(dat, -1), dat[-1])
    dat = layer.apply(dat, function(i) i)
}
