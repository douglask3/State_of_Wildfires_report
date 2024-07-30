graphics.off()
library(terra)
library(modi)

spar = 0.75

runs = c('isimip-final', 'tuning12', 'final')
runs = c('tuning16-r2', 'tuning12', 'final')
mn_samples = c(1, NaN, 2)
date_tests = c('2023-06', '2023-08', '2023-10')
regions = c('Canada', 'Greece', 'NW_Amazon')
region_names = c('Canada - June', 'Greece - Aug', 'Western Amazonia - Sept/Oct')
pcs = c('pc-95', 'pc-95', 'pc-95')

mean_95 <- function(i, burnt_area) {
    #weighted.quantile(burnt_area[[i]][], gridArea[[i]], prob = 0.95, plot = FALSE)
    
    vr = values(burnt_area[[i]])
    val = weighted.quantile(vr, vArea, prob = 0.95, plot = FALSE)
    
    test = vr >= val & !is.na(vr)
    out = sum(vr[test] * vArea[test])/sum(vArea[test])
    return(out)
}


trans_fun <- function(x) log10(x)
itrans_fun <- function(x) 10^x


trans_fun <- function(x) x*100
itrans_fun <- function(x) x
reds9   = c("#FFF7FB", "#F7DEEB", "#EFC6DB", "#E19ECA", "#D66BAE", 
            "#C64292", "#B52171", "#9C0851", "#6B0830")
plot_region <- function(region , run, pc, date_test, region_name, mn_sample) {
    fname = paste0("outputs/ConFire_", region, "-", run, "/figs/_13-frac_points_0.5-")
    fileEx = paste0("-control_TS/", pc, "/points-Control.csv")
    load(paste0("outputs/obs_time_series/", region, "/outs.Rd"))
    if (percentile == 1) percentile = 1 - 1/60
    print(percentile)
    
    plot_experiments <- function(experiments, extreme, cols = c('#FF0000', '#0000FF'), 
                                 scols = reds9, conPeriod = "today", expPeriod = "by 2090",
                            title = '') {
    
        openDat <- function(exp, ny) {
            dat = tail(t(trans_fun(0.000001+read.csv(paste0(fname, exp, fileEx)))), ny)
            if (!is.nan(mn_sample)) dat = dat[seq(mn_sample, nrow(dat), by = 3),]
            return(dat)
        }
    
        dats = lapply(experiments, openDat, ny = 20*3)
        if (length(dats)>2) dats = c(dats[1], list(unlist(dats[-1])))
        #dats = lapply(dats, function(i) i[i>0])
        
        scale = extreme/quantile(dats[[1]], c(percentile))
        dats = lapply(dats, '*', scale)
        
        bins = range(unlist(dats))
        bins = seq(bins[1], bins[2], length.out = floor(sqrt(length(unlist(dats[[1]])))))
        bins0 = bins
        bins = c(sort(seq(bins[1], 0, by = -diff(bins[1:2]))[-1]), 
                    bins, tail(bins, 1) + diff(tail(bins, 2)) * 1:(round(0.1 * length(bins))))
    
        ys = lapply(dats, function(x) hist(unlist(x), bins, plot = FALSE)$density)
        mys = max(unlist(ys))
        ys = lapply(ys, '/', mys)
        #ys = lapply(ys, function(i) log(i + 1))
        x = bins[-1] - diff(bins)
        
        #plot(range(bins0, extreme),range(unlist(ys)), 
        #    type = 'n', axes = FALSE, xlab = '', ylab = '')
        #axis(1, at = seq(-100, 100), itrans_fun(seq(-100, 100)))
        #axis(1)
        #mtext(title, line = -1, side = 3, font = 2)
        addPoly <- function(y, col) {
            addLine <- function(xi, yi) 
                lines(rep(xi, 2), c(0, yi), col = paste0(col, '44'))
            #mapply(addLine, x, y)
            smoothx = seq(min(x), max(x), length.out = 1000)
            ly = log(y + 0.0001)
            xs = log(x)
            #browser()
            smoothy = exp(predict(smooth.spline(ly~x, spar = spar), x = smoothx)[[2]])
            smoothy = smoothy / max(smoothy) 
            
            #polygon(c(smoothx[1], smoothx, tail(smoothx, 1)),
            #        c(log(1), smoothy, log(1)), col = paste0(col, '44'), border = NA)
            sample(smoothx, 10000, replace = TRUE, prob = smoothy)
        }
    
        datsS = mapply(addPoly, c(ys, rev(ys)), c(cols, rev(cols)))
    
        #extreme = burnt_area#quantile(unlist(dats[[1]]), 0.99)
        #lines(rep(extreme, 2), c(0, 9E9), lty = 2)
        
        extreme_fut = round(mean(dats[[1]] > extreme)/mean(dats[[2]] > extreme), 1)
        extreme_fut = round(mean(datsS[,1] > extreme)/mean(datsS[,2] > extreme), 1)
        pc = mean(dats[[1]] > extreme)
        
        #quantile(dats[[2]][,i], 1-pc)#
        pc_diff <- function(i) 100*(quantile(dats[[1]][,i], 1-pc)/quantile(dats[[2]][,i], 1-pc)-1)
        pcs = sapply(1:ncol(dats[[1]]), pc_diff)
        pcsq = round(quantile(pcs, c(0.1, 0.5, 0.9)), 1)
        #dev.new()
        x = dats[[1]]
        y = 100*(-1 + dats[[1]]/dats[[2]])
        #y = dats[[1]] - dats[[2]]
                #cols = densCols(x,y, colramp = colorRampPalette(reds9), bandwidth = 0.1)
        cols = densCols(x,y, colramp = colorRampPalette(scols), bandwidth = 0.5)
        #cols = densCols(x,y, nbin = 128* 4, bandwidth = 0.1)
        
        plot(x, y, type = 'n', xlab = '', ylab = '', axes = FALSE, ylim = quantile(y, c(0.01, 0.99)))
        axis(1)
        for (cex in seq(2, 0.2, by = -0.2)) points(y~x, pch = 20, col = cols, cex = 1)
        grid()
        axis(1)
        axis(2)
        lines(rep(extreme, 2), c(-9E9, 9E9), lty = 2)
        lines( c(-9E9, 9E9), c(0, 0), lty = 2)
            #browser()
        #plot(dats[[1]], 100*(-1 + dats[[1]]/dats[[2]]))
        pval = round(100*mean(pcs<0), 1)
        if (pval < 50) pval = 100 - pval
        if (pval == 100) pval = '> 99.9'
        if (region == regions[1] && FALSE) {
            mtext(side = 2, adj = 0.5, line = 2.5, paste(conPeriod, "vs", expPeriod), 
                  xpd = NA, font = 2)
            arrows(0, 0, 0, 1, xpd = NA)
            mtext(side = 2, adj = 0.7, line = 0.5, 'Likelihood')        
        }
        if (experiments[1] == "factual") 
            mtext(side = 3, adj = 1, line = 0, region_name, xpd = NA, font = 2)
        #mtext(side = 3, adj = 1, line = -0.5, paste0("2023 risk ratio: ", extreme_fut))
        mtext(side = 3, adj = 1, line = -2, 
              paste0('Impact: ', pcsq[2], '% [', pcsq[1], '-',  pcsq[3], ']'))
        mtext(side = 3, adj = 1, line = -3.5, paste0("Likelihood: ", pval, '%'))
        if (region == regions[1]) {
            mtext(side = 2, paste(conPeriod, 'vs', expPeriod), line = 2)
        }
        #if (region == tail(regions,1)) {
        #    for (pch in c(16, 19))
        #        legend('right', c(conPeriod, expPeriod, ''), pch = pch, 
        #               col = c(paste0(cols, '44'), '#FFFFFF00'),
        #               bty = 'n')
        #}
        return(quantile(dats[[2]], percentile))

    }
    extreme = plot_experiments(c("factual", "counterfactual"),
                               extreme = trans_fun(burnt_area_event),
                                c('#FF0000', '#0000FF'), scols = reds9, "Factual", "Counter")
    
    plot_experiments(c("counterfactual", "early_industrial"), extreme, 
                     c('#0000FF', '#000000'), scols = blues9, "Counter", "Early Industrial")
}

plot_fun <- function() {
    par(mfcol = c(2, 3), oma = c(1, 3.25, 0, 0), mar = c(3, 1.5, 3, 0.5))
    mapply(plot_region, regions, runs, pcs, date_tests, region_names, mn_samples)
    mtext(side = 1, 'Mean monthly burned area fraction (%)', outer = TRUE)
    mtext(side = 2, 'Change in BA extent (%)', outer = TRUE, line = 2)
}
    

png("figs/Figure_14_factual_count_attribution.png", height = 6, width = 9, units = 'in', res = 300)
    plot_fun()
dev.off()

pdf("figs/Figure_14_factual_count_attribution.pdf", height = 6, width = 9)
    plot_fun()
dev.off()



               
