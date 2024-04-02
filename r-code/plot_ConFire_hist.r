graphics.off()
library(terra)


region = 'Greece'
date_test = '2023-08'

fname = "outputs/ConFire_UK/figs/crop_lightn_soilM_trees_csoil_pas_vpd_cveg_precip_tas_rhumid_totalVeg-frac_points_0.5-control_TS/points-"
fname = paste0("outputs/ConFire_", region, "-logitnormal12/figs/_13-frac_points_0.1-")
fileEx = "-control_TS/points-Control.csv"

burnt_area_data = paste0("data/data/driving_data/", region, "/isimp3a/obsclim/GSWP3-W5E5/period_2000_2019/burnt_area-2000-2023.nc")


experiments = c("factual", "ss126_GFDL", "ss126_IPSL")


burnt_area = rast(burnt_area_data)
date_test = substr(time(burnt_area), 1, 7) == date_test
#burnt_area_event = burnt_area[[date_test]]
gridArea = cellSize(burnt_area)

burnt_area_tot = sapply(1:nlyr(burnt_area), function(i) sum((burnt_area[[i]] * gridArea)[], na.rm = TRUE)) / sum(gridArea[], na.rm = T)

burnt_area_event = burnt_area_tot[date_test]

percentile = mean(burnt_area_tot <= burnt_area_event)^(1/12)
#burnt_area_event = burnt_area_event * cellSize(burnt_area_event) 
#burnt_area_event = sum(burnt_area_event[], na.rm = T)/ sum(gridArea[], na.rm = T)
#burnt_area_tot / sum(gridArea[], na.rm = T)
#experiments = c("counterfactual", "factual")

trans_fun <- function(x) log10(x)
itrans_fun <- function(x) 10^x


trans_fun <- function(x) x*100
itrans_fun <- function(x) x

extreme = trans_fun(burnt_area_event)
plot_experiments <- function(experiments, conPeriod = "today", expPeriod = "by 2090",
                            title = '') {
    cols = c('#0000FF', '#FF0000')
    openDat <- function(exp, ny) {
        dat = tail(t(trans_fun(0.000001+read.csv(paste0(fname, exp, fileEx)))), ny)
    }
    
    dats = lapply(experiments, openDat, ny = 20)
    if (length(dats)>2) dats = c(dats[1], list(unlist(dats[-1])))
    #dats = lapply(dats, function(i) i[i>0])

    scale = extreme/quantile(dats[[2]], c(percentile))
    dats = lapply(dats, '*', scale)
    bins = range(unlist(dats)); bins = seq(bins[1], bins[2], 
                 length.out = floor(sqrt(length(unlist(dats[[1]])))))

    ys = lapply(dats, function(x) hist(unlist(x), bins, plot = FALSE)$density)
    ys = lapply(ys, function(i) log(i + 1))
    x = bins[-1] - diff(bins)
    
    plot(range(x, extreme), range(unlist(ys)), 
         type = 'n', axes = FALSE, xlab = '', ylab = '')
    #axis(1, at = seq(-100, 100), itrans_fun(seq(-100, 100)))
    axis(1)
    mtext(title, line = -1, side = 3, font = 2)
    addPoly <- function(y, col) 
        polygon(c(x[1], x, tail(x, 1)), c(log(1), y, log(1)), col = paste0(col, '44'), border = NA)

    mapply(addPoly, c(ys, rev(ys)), c(cols, rev(cols)))
    
    #extreme = burnt_area#quantile(unlist(dats[[1]]), 0.99)
    lines(rep(extreme, 2), c(0, 9E9), lty = 2)
    
    extreme_fut = round(mean(dats[[2]] > extreme)/mean(dats[[1]] > extreme), 3)
    #mtext(side = 3, adj = 0.5, line = -4, paste0("1-in-100 ", conPeriod,  " occurs\n1-in-", 
    #                                             extreme_fut, " ", expPeriod))
    
    mtext(side = 3, adj = 0.5, line = -4, paste0("2023 risk ratio: ", extreme_fut))

}
png(paste0(region, "_attribution.png"), height = 3, width = 7, units = 'in', res = 300)
par(mfrow = c(1, 3), mar = c(3, 0.5, 3, 0.5))
plot_experiments(c("counterfactual", "factual"), "at PI", "today")
#plot_experiments(c("factual", "ss126_GFDL", "ss126_IPSL", "ss126_MPI", "ss126_MRI"), title = "ssp126")
#plot_experiments(c("factual", "ss585_GFDL", "ss585_IPSL", "ss585_MPI", "ss585_MRI"), title = "ssp585")
dev.off()               
