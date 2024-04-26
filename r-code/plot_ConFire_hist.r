graphics.off()
library(terra)
library(modi)

region = 'Greece'
date_test = '2023-08'

fname = paste0("outputs/ConFire_", region, "-tuning12/figs/_13-frac_points_0.5-")
fileEx = "-control_TS/points-Control.csv"

burnt_area_data = paste0("data/data/driving_data/", region, "/isimp3a/obsclim/GSWP3-W5E5/period_2000_2019/burnt_area-2000-2023.nc")


experiments = c("factual", "ss126_GFDL", "ss126_IPSL")


burnt_area = rast(burnt_area_data)
date_test = substr(time(burnt_area), 1, 7) == date_test
#burnt_area_event = burnt_area[[date_test]]
gridArea = cellSize(burnt_area[[1]])
vArea = values(gridArea)
mean_95 <- function(i) {
    #weighted.quantile(burnt_area[[i]][], gridArea[[i]], prob = 0.95, plot = FALSE)
    
    vr = values(burnt_area[[i]])
    val = weighted.quantile(vr, vArea, prob = 0.95, plot = FALSE)
    
    test = vr >= val & !is.na(vr)
    out = sum(vr[test] * vArea[test])/sum(vArea[test])
    return(out)
}

burnt_area_tot = sapply(1:nlyr(burnt_area), mean_95)#function(i) sum((burnt_area[[i]] * gridArea)[], na.rm = TRUE)) / sum(gridArea[], na.rm = T)

burnt_area_event = burnt_area_tot[date_test]
burnt_area_tot = burnt_area_tot[sort(unlist(lapply(7:8, function(i) seq(i, nlyr(burnt_area), by = 12))))]
percentile = mean(burnt_area_tot <= burnt_area_event)
#burnt_area_event = burnt_area_event * cellSize(burnt_area_event) 
#burnt_area_event = sum(burnt_area_event[], na.rm = T)/ sum(gridArea[], na.rm = T)
#burnt_area_tot / sum(gridArea[], na.rm = T)
#experiments = c("counterfactual", "factual")

trans_fun <- function(x) log10(x)
itrans_fun <- function(x) 10^x


trans_fun <- function(x) x*100
itrans_fun <- function(x) x

extreme = trans_fun(burnt_area_event)
plot_experiments <- function(experiments, cols = c('#FF0000', '#0000FF'), conPeriod = "today", expPeriod = "by 2090",
                            title = '') {
    
    openDat <- function(exp, ny) {
        dat = tail(t(trans_fun(0.000001+read.csv(paste0(fname, exp, fileEx)))), ny)
    }
    
    dats = lapply(experiments, openDat, ny = 20*3)
    if (length(dats)>2) dats = c(dats[1], list(unlist(dats[-1])))
    #dats = lapply(dats, function(i) i[i>0])
    
    scale = extreme/quantile(dats[[1]], c(percentile))
    dats = lapply(dats, '*', scale)
    
    bins = range(unlist(dats)); bins = seq(bins[1], bins[2], 
                 length.out = floor(sqrt(length(unlist(dats[[1]])))))
    bins = c(bins, tail(bins, 1) + diff(tail(bins, 2)) * 1:(round(0.1 * length(bins))))
    ys = lapply(dats, function(x) hist(unlist(x), bins, plot = FALSE)$density)
    ys = lapply(ys, function(i) log(i + 1))
    x = bins[-1] - diff(bins)
    
    plot(range(x, extreme), range(unlist(ys)), 
         type = 'n', axes = FALSE, xlab = '', ylab = '')
    #axis(1, at = seq(-100, 100), itrans_fun(seq(-100, 100)))
    axis(1)
    mtext(title, line = -1, side = 3, font = 2)
    addPoly <- function(y, col) {
        addLine <- function(xi, yi) 
            lines(rep(xi, 2), c(0, yi), col = paste0(col, '44'))
        mapply(addLine, x, y)
        smoothx = seq(min(x), max(x), length.out = 1000)
        ly = log(y + 0.0001)
        
        smoothy = exp(predict(smooth.spline(ly~x, spar = 0.75), x = smoothx)[[2]])
        
        polygon(c(smoothx[1], smoothx, tail(smoothx, 1)),
                c(log(1), smoothy, log(1)), col = paste0(col, '44'), border = NA)
    }

    mapply(addPoly, c(ys, rev(ys)), c(cols, rev(cols)))
    
    #extreme = burnt_area#quantile(unlist(dats[[1]]), 0.99)
    lines(rep(extreme, 2), c(0, 9E9), lty = 2)
    
    extreme_fut = round(mean(dats[[1]] > extreme)/mean(dats[[2]] > extreme), 2)
    pc = mean(dats[[1]] > extreme)
    
    
    pc_diff <- function(i) 100*(quantile(dats[[1]][,i], pc)/quantile(dats[[2]][,i], pc)-1)
    pcs = sapply(1:ncol(dats[[1]]), pc_diff)
    pcsq = round(quantile(pcs, c(0.1, 0.9)), 2)
    pval = round(mean(pcs<0), 2)
    
    mtext(side = 3, adj = 1, line = 1, paste(conPeriod, "vs", expPeriod), xpd = NA, font = 2)
    mtext(side = 3, adj = 1, line = -0.5, paste0("2023 risk ratio: ", extreme_fut))
    mtext(side = 3, adj = 1, line = -2, paste0("p-val: ", pval))
    mtext(side = 3, adj = 1, line = -3.5, paste0('Impact: ', pcsq[1], ' to ', pcsq[2], '%' ))

    for (pch in c(16, 19))
        legend('right', c(conPeriod, expPeriod, ''), pch = pch, 
               col = c(paste0(cols, '44'), '#FFFFFF00'),
               bty = 'n')
    return(quantile(dats[[2]], percentile))

}
png(paste0(region, "_attribution.png"), height = 5, width = 9, units = 'in', res = 300)
par(mfrow = c(1, 2), mar = c(3, 0.5, 3, 0.5))
extreme = plot_experiments(c("factual", "counterfactual"), c('#FF0000', '#0000FF'), "Factual", "Counter")
plot_experiments(c("counterfactual", "early_industrial"), c('#0000FF', '#000000'), "Counter", "Early Industrial")
#plot_experiments(c("factual", "ss126_GFDL", "ss126_IPSL", "ss126_MPI", "ss126_MRI"), title = "ssp126")
#plot_experiments(c("factual", "ss585_GFDL", "ss585_IPSL", "ss585_MPI", "ss585_MRI"), title = "ssp585")
dev.off()               
