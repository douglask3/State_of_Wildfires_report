source("../rasterextrafuns/rasterPlotFunctions/R/plot_raster_map.r")

file = 'outputs/ConFire_Canada-nrt4/figs/_17-frac_points_0.2-baseline-control_TS/points-Control.csv'

dat = read.csv(file, stringsAsFactors = TRUE)

cols = c('#ffffb2','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#b10026')
limits = c(0.001, 0.01, 0.02, 0.05, 0.1, 0.2)

dat = cut_results(dat, limits)

start = 2011 + 1/24
x = seq(start, length.out = ncol(dat), by = 1/12)

plot(range(x), c(-1, 1), type = 'n', axes = FALSE)

plot_col <- function(i) {
    xi = x[i] + c(-1, 1) * diff(x[1:2])/2
    zi = rev(sort(dat[, i]))
    for_row <- function(j) {
        lines(xi,  rep(j/nrow(dat), 2), col = cols[zi[j]], lwd = 5)
        lines(xi, -rep(j/nrow(dat), 2), col = cols[zi[j]], lwd = 5)
    }
    lapply(1:nrow(dat), for_row)
    lines(x[c(1, length(x))], c(0, 0))
}

lapply(1:length(x), plot_col)





