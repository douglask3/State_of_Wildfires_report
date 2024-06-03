library(raster)
source("../gitProjectExtras/gitBasedProjects/R/sourceAllLibs.r")
sourceAllLibs("../rasterextrafuns/rasterPlotFunctions/R/")
sourceAllLibs("../rasterextrafuns/rasterExtras/R/")
sourceAllLibs("../gitProjectExtras/gitBasedProjects/R/")
library(ncdf4)
sourceAllLibs("../ConFIRE_attribute/libs/")
source("../ConFIRE_attribute/libs/plotStandardMap.r")
source("../LPX_equil/libs/legendColBar.r")
source("libs/find_levels.r")
graphics.off()

cols = c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026')
dcols = rev(c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061'))

dcols2 = rev(c('#7f3b08','#b35806','#e08214','#fdb863','#fee0b6','#f7f7f7','#d8daeb','#b2abd2','#8073ac','#542788','#2d004b'))


levels = c(0.001, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)#find_levels(c(mn[!is.na(mn)], mx[!is.na(mx)]), seq(10, 90, 10))

dlevels = c(0.5, 0.6, 0.8, 1, 1.25, 1.75,2)
dlevels2 = c(0.1, 0.2, 0.5, 0.75, 1, 1.5, 2, 5, 10)
nsample = 10


region = 'NW_Amazon'
dir = paste0("outputs/ConFire_", region, "-tuning12/samples/_13-frac_points_0.5/")
eg_rast = raster(paste0("outputs/ConFire_", region, "-tuning12/samples/_13-frac_points_0.5/baseline-/control/sample-pred0.nc"))
plt_width = c('Canada' = 1, 'NW_Amazon' = 1, 'Greece' = 0.67)[region]
plt_height = plt_width*nrow(eg_rast)/ncol(eg_rast)


cols = make_col_vector(cols, ncols = length(levels)+1)
dcols = make_col_vector(dcols, ncols = length(dlevels)+1)
dcols2 = make_col_vector(dcols2, ncols = length(dlevels2)+1)



plot_control <- function(run, control, name, years = c(2000, 2019), cols, levels, dats0 = NULL, addLab = FALSE, ...) {
    
    tfile = paste0(c('temp2/plot-ConFire_futureFire/', region, run, control, name, years), 
                   collapse = '-')
    get_dat <- function(id) {
        dirs = list.files(paste0(dir, id), full.names=TRUE)
        dirs = paste0(dirs, '/', control, '/')
        filess = lapply(dirs, list.files, pattern = 'pred', full.names = TRUE)
        filess = lapply(filess, function(x) x[seq(1, length(x), length.out = nsample)])
        open_model <- function(files) {
            open_file <- function(file) {
                print(file)
                dat = brick(file)
                yrs = as.numeric(substr(names(dat),2, 5))
                yrID = which((yrs >= years[1]) & (yrs < years[2]))
                if (length(yrID) == 0) return(NULL)
                dat[[yrID]]
            }
            lapply(files, open_file)
        }
        dat = lapply(filess, open_model)
    }
    combine_runs <- function(id, md) {
        tfile = paste(tfile, id, md, '.nc', sep = '-')
        print(tfile)
        if (file.exists(tfile)) return(raster(tfile))
        hist = dats[[1]][[md]][[id]]
        futr = dats[[2]][[md]][[id]]
        if (is.null(hist)) out = mean(futr)
            else out = mean(addLayer(hist, futr))
        out = writeRaster(out, file = tfile)
        return(out)
    }
    combine_mod <- function(md) 
        layer.apply(1:length(dats[[1]][[md]]), combine_runs, md)

    tfile_all= paste0(tfile, '-all.Rd')
    if (file.exists(tfile_all)) {
        load(tfile_all)
    } else {
        dats = lapply(c('historical', run), get_dat)
        dats = lapply(1:length(dats[[1]]), combine_mod)
        save(dats, file = tfile_all)
    }
    dats = do.call(addLayer, dats)
    dats_out = dats
    if (!is.null(dats0)) dats = dats/dats0 #dats = dats - dats0
    
    #if (ncol(dats) > 40) dats = raster::aggregate(dats, fact = ncol(dats)/40)

    tfile_cdat= paste0(tfile, ncol(dats), '-cut_dat.Rd')
    if (file.exists(tfile_cdat)) {
        load(tfile_cdat)
    } else {
        cdat = cut_results(dats, levels)
        save(cdat, file = tfile_cdat)
    }
    
    
    cdat = layer.apply(1:(length(levels) + 1), function(i) mean(cdat == i))
    mask = which(!is.na(cdat[[1]][]))
    vcdat = t(apply(cdat[mask], 1, function(i) rev(cumsum(rev(i)))))^0.5
    vcdat = 15 * vcdat / ncol(cdat)
    xydat = cbind(xyFromCell(cdat, mask)) 

    if (ncol(dats) > 40) {
        y = seq(min(xydat[,2]), max(xydat[,2]), length.out = 40)
        x = seq(min(xydat[,1]), max(xydat[,1]), by = diff(y[1:2]))
        xydat_new = cbind(rep(x, length(y)), rep(y, each = length(x)))

        assign_xy <- function(xy) 
            which.min((xy[1] - xydat_new[,1])**2 + (xy[2] - xydat_new[,2])**2)
       
        id = apply(xydat, 1, assign_xy)
        vcdat_new = matrix(NaN, ncol = dim(vcdat)[2], nrow = nrow(xydat_new))
        for (i in 1:nrow(xydat_new)) {
        
            test = which(id == i)
            if (length(test) == 1) vcdat_new[i,] = vcdat[test,]
            else if (length(test) > 1 ) vcdat_new[i,] = apply(vcdat[test,], 2, mean)
        }
        test = !is.na(vcdat_new[,1])
        xydat = xydat_new[test,]
        vcdat = vcdat_new[test,]
    }
    
    plot(range(xydat[,1]), range(xydat[,2]), type = 'n', xlab = '', ylab = '', 
              xaxt = 'n', yaxt = 'n', xaxs = 'i', yaxs ='i')
    addLevel <- function(i) points(xydat[,1], xydat[,2], pch = 19, 
                                   cex = plt_width*3*vcdat[,i], col = cols[i])
    lapply(1:length(cols), addLevel) 
    contour(is.na(eg_rast), levels = 0.5, drawlabels=FALSE, add = TRUE)   
    #legendColBar(c(0.1, 0.7), c(0.1, 0.9), cols = cols, limits = levels, ...)
    return(dats_out)
}

plot_controls <- function(years, runs, control, name) {
    lmat = t(matrix(1:9, nrow = 3))
    lmat = rbind(lmat, 10:12)

    lmat = rbind(c(1, 0, 2, 0, 3), 0, c(4, 0, 5, 0, 6), 0, c(7, 0, 8, 0, 9), c(10, 0, 11, 0, 12))
    lmat = cbind(0, rbind(0, lmat, 0), 0)
   
    
    heights = c(0.4, rep(c(0.02, plt_height), 3)[-1], 0.3, 0.1)
    widths = c(0.2, plt_width, 0.02, plt_width, 0.02, plt_width, 0.2)

    png(paste("figs/ConFire_histMaps_fires", region, years, control, name, ".png", sep = '-'), 
        height = 2.2*sum(heights), width = 2.2*sum(widths), 
            res = 300, units = 'in')
    layout(lmat, heights = heights, widths = widths)
    par(mar = c(0, 0, 0, 0))
    plot_ssp <- function(run) {
        dats = plot_control(run, control, name, c(2010, 2019), cols, levels, 
                            extend_min = F, minLab = 0, addLab = T)
        if (run == runs[1]) {
            axis(3)
            mtext(side = 3, '2010 - 2020', line = 2)
        }
        if (run == tail(runs, 1)) axis(1)
        mtext(side = 2, run, line = 2)
        axis(2)

        ddats = plot_control(run, control, name, years, dcols, dlevels, dats, extend_min = T)
        if (run == runs[1]) {
            axis(3)
            mtext(side = 3, paste(years[1], '-', years[2] + 1), line = 2)
        }
        if (run == tail(runs, 1)) axis(1)
        
        qx = logit(head(seq(0, 1, length.out = nlayers(dats)+2)[-1], -1))
        quantile_smooth <- function(x, qu = 0.99, na.rm = TRUE) {
            if (all(is.na(x))) return(NaN)
            y = sort(x)
            predict(smooth.spline(y~qx, spar = 0.5), logit(0.99))[[2]]
        } 
        gt_smooth <- function(x) {
            if (all(is.na(x))) return(NaN)
            target = x[1]
            y = sort(x[-1])
            
            1-logistic(predict(smooth.spline(qx~y, spar = 0.5), target)[[2]])
        } 
        q99 = calc(dats, function(x) quantile(x, 0.99, na.rm = TRUE))
        q99 = calc(dats, quantile_smooth)#function(x) quantile(x, 0.99, na.rm = TRUE))
        rchange = calc(addLayer(q99, ddats) , gt_smooth)*100
        #rchange = mean(ddats > q99) / mean(dats > q99)
        
        #if (ncol(rchange) > 40) rchange = raster::aggregate(rchange, fact = ncol(dats)/40)
        mask = which(!is.na(rchange[[1]][]))
        xydat = cbind(xyFromCell(rchange, mask), cut_results(rchange[mask], dlevels2))
        
        plot(range(xydat[,1]), range(xydat[,2]), type = 'n', xlab = '', ylab = '', 
              xaxt = 'n', yaxt = 'n', xaxs = 'i', yaxs ='i')
        add_poly <- function(vs) {
            polygon(vs[1] + 0.25 * c(-1, -1, 1, 1, -1), vs[2] + 0.25 *c(-1, 1, 1, -1, -1), 
                    col = dcols2[vs[3]], border = NA)
        }
        apply(xydat, 1, add_poly)
        contour(is.na(eg_rast), levels = 0.5, drawlabels=FALSE, add = TRUE)
        if (run == runs[1]) {
            mtext(side = 3, 'Change in\n1-in-100 event', line = 2)
            axis(3)
        }
        if (run == tail(runs, 1)) axis(1)
        axis(4)
    }
    lapply(runs, plot_ssp)
    
    add_leg <- function(cols, levels) {
        plot(c(0, 1), c(0, 1), type = 'n', axes = FALSE, xlab = '', ylab = '')
        xp = seq(0.05, 0.95, length.out = length(cols))
        points(xp, rep(0.05, length(cols)), pch = 19, col = cols)
        labs = c(0, levels, '+')
        xt = c(xp[1], xp[-1] - diff(xp)/2, tail(xp, 1))
        text(xt, rep(0.05, length(xt)), labs, adj = c(0.5, -1), srt = 30)
    }
    add_leg(cols,  levels)
    add_leg(dcols, dlevels) 
    legendColBar(c(0.0, 0.37), c(0, 1), dcols2, dlevels2, transpose = TRUE, oneSideLabels=FALSE)
    dev.off()
}

years = list(c(2030, 2039), c(2040, 2049), c(2090, 2099))


lapply(years, plot_controls, c('ssp126', 'ssp370', 'ssp585'), "control", "burnt_area")

