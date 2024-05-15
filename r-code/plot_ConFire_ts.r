graphics.off()

if (F) {
burnt_area_data = paste0("data/data/driving_data/", region, "/isimp3a/obsclim/GSWP3-W5E5/period_2000_2019/burnt_area-2000-2023.nc")
date_test = '2023-06'
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
}

rnning_mean <- function(r) 
    filter(r, rep(1 / 10, 10), sides = 1)

make_plot <- function(dir, file) {
    openDat <- function(exp) {
        print(exp)
        dirs = list.dirs(dir, full.name = TRUE, recursive=TRUE)
        dirs = dirs[grep(paste0('-', exp), dirs)]
        for (ptt in pattern) dirs = dirs[grep(paste0(ptt), dirs)]
        if (length(dirs) > 1) { 
            if (length(dirs) == 6) dirs = dirs[-1]
            else dirs = dirs[!(grepl('mean', dirs) | grepl('pc-', dirs))]
            if (length(dirs) > 1) if (length(dirs) == 6) dirs = dirs[-1]
        }

        openDir <- function(dir) 
            read.csv(paste0(dir, '/', file))
        
        dat = lapply(dirs, openDir)
        return(dat)
    }
    dats = lapply(experiments, openDat)
    scale = burnt_area_event/quantile(as.matrix(dats[[1]][[1]]), percentile)
    dats = lapply(dats, lapply, function(i)  i[,seq(2, dim(i)[2], by = 3)])
    dats = lapply(dats, lapply, function(i)  i * scale)

    join_dats <- function(dat, id) 
        mapply(cbind, dats[[id]], dat, SIMPLIFY = FALSE)
    dats[4:6] = lapply(dats[4:6], join_dats, 3)
    start[4:6] = start[3]

    years = mapply(function(x, st) st + (1:ncol(x[[1]])) - 1, dats, start)
    
    find_occurnaces <- function(pc) {
        print(pc)
        event = log(quantile(as.numeric(as.matrix(dats[[1]][[1]])), pc))
        find_occurnace <- function(dat) {
            for_run <- function(x) {

                find_pc <- function(r) {
                    r = log(unique(sort(r)))
                    xr = seq(0, 1, length.out = length(r)+1)
                    xr = xr[-1] - diff(xr)/2
                    xr = log(xr/(1-xr))
                    xp = seq(-30, 30, length.out = 2500)
                    rp = predict(smooth.spline(r ~ xr, spar = 0.5), xp)[[2]]
                    rp = sort(rp)
                    #xp = xr
                    #rp = r
                    if (event > tail(rp, 1)) {
                        #out = tail(xp, 1) + xp[which(rp >= 2*tail(rp, 1) - event)[1]]
                        out = -tail(xp, 1)
                    } else {
                        out = -xp[which(rp >= event)[1]]
                    }
                    
                    out = 1/(1+exp(-out))
                    if (is.na(out)) out = 0
                    
                    return(out)
                }
                out = apply(x, 2, find_pc)
               
                return(rnning_mean(out))
            }
            lapply(dat, for_run)
        }
        return(lapply(dats, find_occurnace))
    }
    freqs = lapply(percentiles, find_occurnaces)
    
    find_quatiles <- function(dat, year) { 
        for_run <- function(y) {                
            yq = apply(y, 1, rnning_mean)
            test = (year < 2019)
            #if (length(test) == 0) test = which(year == 2014)
            
            yq = apply(yq, 2, function(i) i / mean(i[test], na.rm = TRUE))
            apply(yq, 1, quantile, c(0.1, 0.9), na.rm = TRUE)
        }
        lapply(dat, for_run)
    }

    qdats = mapply(find_quatiles, dats, years, SIMPLIFY = FALSE)
    
    tfun <- function(x) return(100*x)
    add_Freq <- function(freq, pc, nme, ylab = '', axis1 = FALSE) {
        x_range = range(unlist(years, recursive = TRUE))
        y_range = range(tfun(unlist(freq, recursive = TRUE)), na.rm = TRUE)
        plot(x_range, y_range, type = 'n', xlab = '', ylab = '', axes = FALSE)
        mtext(nme, line = -1.5, side = 3)
        at = seq(1990, 2100, by = 10)
        axis(1, at = at, rep('', length(at)))
        if (axis1) axis(1, at = at)
        axis(2)
        mtext(ylab, side = 2, line = 2.5)
        add_experiment <- function(fq, year, col, alpha = '11') {
            add_run <- function(y) 
                lines(year, tfun(y), col = col, lty = 3, lwd = 0.5)
                
            if (length(fq) == 1) {
                lines(year, tfun(fq[[1]]), col = col, lwd = 2)
            } else {
                    
                pfq = apply(do.call(cbind, fq), 1, range)
                yearp = year[!is.na(pfq[1,])]
                pfq = pfq[,!is.na(pfq[1,])]
                if (nchar(col) == 7) colp = paste0(col, alpha) else colp = col
                polygon(c(yearp, rev(yearp)), tfun(c(pfq[1,], rev(pfq[2,]))), 
                        col = colp, border = NA)
                lines(yearp, tfun(pfq[1,]), col = col)
                lines(yearp, tfun(pfq[2,]), col = col)
                lapply(fq, add_run)
            }
            return(do.call(cbind, fq))
        }
        for (i in 1:4) mapply(add_experiment, freq, years, cols)
        fqs = mapply(add_experiment, freq, years, cols, '00')
    }
    png(paste0("futures_", region, ".png"), res = 300, height = 7, width = 7.2, units = 'in')
    par(mfrow = c(3, 2), mar = c(3, 3, 0, 0), oma = c(1, 1, 0, 0))
    fqss = mapply(add_Freq, freqs, percentiles, percentile_names, 
            c('', '', '% Likelihood of event', rep('', 3)), c(rep(F, 4), T, T))

    legend('topleft', legend = c('historic', 'ssp126', 'ssp370', 'ssp585'), col = cols[-(2:3)], pch = 19, pt.cex = 2)
    dev.off()
    return(fqss)
}

regions = c('Canada', 'Greece', 'NW_Amazon')



experiments = c("factual", "counterfactual", "historical", "ssp126", "ssp370", "ssp585")
experiments = c("factual", "counterfactual", "historical", "ssp126", "ssp370", "ssp585")
start = c(2000, 2000, 1994, 2015, 2015, 2015)
#cols = c('#000000', '#0000FF', '#0000FF', '#0000FF', '#FF0000', '#FF0000', '#FF0000')#, '#FF0000')
cols = c('#000000', '#648FFF', '#DC267F00', '#785EF0', '#FFB000', '#FE6100')
cols = c('#000000', '#999999', '#DC267F00', '#7570b3', '#1b9e77', '#d95f02')

percentiles = c(0.5, 0.9, 0.95, 0.99, 0.999)
percentile_names = 1/(1-percentiles)
percentiles = c(percentile, percentiles)
percentile_names = c('2023 event', paste0('1-in-', round(percentile_names, 0)))



start = c(2000, 2000, -999, 1994, 1994, 1994)

plot_region_fqi <- function(control_name, col_hint, region, fi = 1) {
    dir = paste0("outputs/ConFire_", region, "-tuning12/figs/")
    pattern = c("_13-frac_points_0.5-")
    file = paste0("points-", control_name, ".csv")

    load(paste0("outputs/obs_time_series/", region, "/outs.Rd"))
    fqss = make_plot(dir, file)
    fqs = fqss[,fi]

    years = seq(2019, 2099, by = 10)
    if(col_hint != "NULL") cols = sapply(cols, function(col) make_col_vector(c(col, col_hint), ncols = 3)[2])
    get_points_for_year <- function(year) {
        get_dat <- function(fq, st, id) {
            if (id ==2) return(NULL)
            if (is.null(dim(fq))) {
                yrs = st + 1:length(fq) - 1
                id = which(yrs == year)
                if (length(id) == 0) id = which(abs(yrs - year) == 1)
                if (length(id) == 0) return(NULL)
                return(fq[id] * 100)
            } else {
                yrs = st + 1:(dim(fq)[1])
                id = which(yrs == year)
                if (length(id) == 0) id = which(abs(yrs - year) == 1)
                if (length(id) == 0) return(NULL)
                return(fq[id,] * 100)
            }
        }
        mapply(get_dat, fqs, start, 1:length(fqs))
    }

    pnts = lapply(years, get_points_for_year)

    plot_for_year <- function(pnt, year) {
        plot_point <- function(y, col, off) {
            #col = make_col_vector(c(col, col_hint), ncols = 3)[2]
            if (is.null(y)) return()      
            x = year + (off-5)*2 + c(-1, 1)
            if (length(y) == 1) {  
                lines(range(years), c(y, y), col = col, lwd = 2, lty = 2)
            } else {
                yp = range(y)
                polygon(x[c(1, 2, 2, 1, 1)], yp[c(1, 1, 2, 2, 1)], 
                        col = paste0(col, '66'), border = paste0(col, 'BB'))
                lapply(y, function(i) lines(x, rep(i, 2), col = col))
            }
            return(col)
        }
        cols = mapply(plot_point, pnt, cols, 1:length(pnt))
    }

    yrange = range(unlist(pnts))


    plot(c(2010, 2110), yrange, type = 'n', axes = FALSE, xlab = '', ylab = '')
    axis(1, at = years + 1, labels = rep('', length(years)))
    if (control_name == tail(controls, 1)) axis(1, at = years + 1)
    colsp = mapply(plot_for_year, pnts, years)

    yrange4 = yrange/pnts[[1]][[1]]

    labels = round(seq(yrange4[1], yrange4[2], length.out = 6), 1)
    axis(2)
    
    axis(4, at = labels *pnts[[1]][[1]], labels = labels)
    if (control_name == controls[1]) mtext(side = 3, line = -1, region)
    if (region == regions[1]) {
        mtext(side = 2, line = 2, control_name)
        legend('topleft', experiments[-(2:3)], col = paste0(cols[-(2:3)], 'BB'), pt.cex = 2, pch = 15, bty = 'n')
        legend('topleft', experiments[-(2:3)], col = cols[-(2:3)], pt.cex = 2, pch = 1, bty = 'n')
    }
}
controls = c('Control', 'Fuel', 'Moisture', 'Ignition') 
cols_hint = c('NULL', '#00FF00', '#0000FF', '#FF0000')#, '#333333')
png(paste0("figs/box_futures-", region, "-", fi, ".png"), 
    res = 300, units = 'in', width = 8, height = 7)
par(mfcol = c(4, 3), oma = c(2, 4, 2, 4), mar = c(1, 2, 0, 2))
lapply(regions, function(region) mapply(plot_region_fqi, controls, cols_hint, MoreArgs = list(region)))
mtext(side = 2, 'Liklihood (%)', outer = TRUE, line = 2.5)
mtext(side = 4, line = 2, 'times more likely)', outer = TRUE, line = 2.5)
dev.off()
