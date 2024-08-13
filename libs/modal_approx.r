testModal <- function(r) {
    out0 = rep(NaN, 12)
    cellModel <- function(y) {
        if (any(is.na(y))) return(out0)
        if (sum(y) == 0) return(out0)
        
        out0[1] = 1 - min(y)/max(y)
        out0[7] = which.max(y)
        testy = y[c(which.min(y):12, 1:which.min(y))]
        dy = diff(testy)
        peaks = which(testy > c(testy[-1], testy[1]) & testy >= c(tail(testy, 1), head(testy,-1)))
        if (length(peaks) == 1) return(out0)
        
        peaks = peaks[sort.int(testy[peaks], index.return = TRUE)[[2]]]
        out0[8:(6+length(peaks))] = head(peaks - tail(peaks,1),-1)
        prominance <- function(p1, p2)
            (testy[p1] - min(testy[p1:p2]))/max(testy)
        
        out = mapply(prominance, head(peaks, -1), peaks[-1])
        out0[2:(length(out)+1)] = out
        if (length(out) > 5) browser()
        return(out0)
    }
    outr = r[[1:12]]
    outr[] = t(apply(r[], 1,cellModel))
    return(outr)
}
modal_approx <- function(dat) {
    modal = testModal(dat)
    modal_approx = layer.apply(2:6, function(i)
                               modal[[i]] * (0.5-cos(2*pi * modal[[i+6]]/12)/2))
    modal_approx = 1 + sum(modal_approx, na.rm = TRUE)/modal[[1]]
    return(modal_approx)
}
