fname = "outputs/ConFire_UK/figs/crop_lightn_soilM_trees_csoil_pas_vpd_cveg_precip_tas_rhumid_totalVeg-frac_points_0.5-control_TS/points-"
fname = "outputs/ConFire_Canada-biasC3/figs/crop_lightn_soilM_tree_biascorrected_csoil_pas_vpd_cveg_precip_tas_rhumid_totalVeg_biascorrected-frac_points_0.02-"
fileEx = "-control_TS/points-Control.csv"

experiments = c("factual", "ss126_GFDL", "ss126_IPSL")
#experiments = c("counterfactual", "factual")
plot_experiments <- function(experiments, conPeriod = "today", expPeriod = "by 2090",
                             title = '') {
    cols = c('#0000FF', '#FF0000')
    openDat <- function(exp, ny) {
        dat = tail(t(log10(read.csv(paste0(fname, exp, fileEx)))), ny)
    }

    dats = lapply(experiments, openDat, ny = 20)
    if (length(dats)>2) dats = c(dats[1], list(unlist(dats[-1])))

    bins = range(unlist(dats)); bins = seq(bins[1], bins[2], length.out = 20)

    ys = lapply(dats, function(x) hist(unlist(x), bins, plot = FALSE)$density)
    x = bins[-1] - diff(bins)
    plot(range(x), range(unlist(ys)), type = 'n', axes = FALSE, xlab = '', ylab = '')
    axis(1, at = seq(-100, 100), 10^(seq(-100, 100)))

    mtext(title, line = -1, side = 3, font = 2)
    addPoly <- function(y, col) 
        polygon(c(x[1], x, tail(x, 1)), c(0, y, 0), col = paste0(col, '44'), border = NA)

    mapply(addPoly, c(ys, rev(ys)), c(cols, rev(cols)))

    extreme = quantile(unlist(dats[[1]]), 0.99)
    lines(rep(extreme, 2), c(0, 9E9), lty = 2)
    extreme_fut = round(1/(mean(dats[[2]] > extreme)))
    mtext(side = 3, adj = 0.5, line = -4, paste0("1-in-100 ", conPeriod,  " occurs\n1-in-", 
                                                 extreme_fut, " ", expPeriod))

}
png("UK_attribution.png", height = 3, width = 7, units = 'in', res = 300)
par(mfrow = c(1, 3), mar = c(3, 0.5, 3, 0.5))
plot_experiments(c("counterfactual", "factual"), "at PI", "today")
plot_experiments(c("factual", "ss126_GFDL", "ss126_IPSL", "ss126_MPI", "ss126_MRI"), title = "ssp126")
plot_experiments(c("factual", "ss585_GFDL", "ss585_IPSL", "ss585_MPI", "ss585_MRI"), title = "ssp585")
dev.off()               
