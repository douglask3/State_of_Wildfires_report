source("libs/YearlySeason.r")


srank.raster <- function(r1, r2 = NULL, lab = '', name = '',season = NULL, plot = TRUE) {
    if(!is.null(season)) 
        r1 = YearlySeason(season, r1)

    if (is.null(r2)) 
        r2 = r1[[sample(1:nlayers(r1), nlayers(r1), replace = FALSE)]]
    else if(!is.null(season))    
        r2 = YearlySeason(season, r2)
    
    mask = !any(is.na(r1+r2))
    srank.cell <- function(v1, v2) 
         cor.test(v1, v2, method = "spearman")[[4]]
    
    out = r1[[1]]
    out[mask] = mapply(srank.cell, as.data.frame(t(r1[mask])), as.data.frame(t(r2[mask])))
    out[!mask] = NaN

    if (!plot) return(out)

    plotStandardMap(out, limits = limits_rank, cols = cols_rank)
    #mtext(name, side = 2, line = -2)
    addLetLab(lab, name, line = 0)
    StandardLegend(out, limits = limits_rank, cols = cols_rank,
                   extend_max = FALSE, maxLab = 1, add = TRUE, oneSideLabels = FALSE)
    return(out)
}
