graphics.off()
library(stringr) 
source("libs/legendColBar.r")
source("libs/find_levels.r")
source("libs/make_col_vector.r")
source("libs/mtext.units.r")
regions = c('Canada', 'Greece', 'NW_Amazon')
pcs = c(95, 95, 95)
Anom_titles = c('Canada', 'Greece', 'Western Amazonia')
xlims = list(c(3, 9)*30, c(6, 9)*30, c(6, 12)*30)

dirs = paste0("outputs/ConFire_", regions, "-nrt-tuning12/figs/_12-frac_points_0.5-baseline-control_TS/pc-", pcs, "/")

files = list("points-Control.csv",
	     "points-standard-Fuel.csv",
	     "points-standard-Moisture.csv",
	     "points-standard-Weather.csv",
	     c("points-standard-Ignition.csv", "points-standard-Suppression.csv"),
             c("residual", "points-Evaluate.csv", "points-Control.csv"))

cols = list(c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c',
              '#fc4e2a','#e31a1c','#bd0026','#800026', '#400013'),
	    c('#ffffe5','#f7fcb9','#d9f0a3','#addd8e','#78c679',
              '#41ab5d','#238443','#006837','#004529', '#002315'),
	    rev(c('#ffffd9','#edf8b1','#c7e9b4','#7fcdbb','#41b6c4',
                  '#1d91c0','#225ea8','#253494','#081d58', '#041034')),
	    c('#f7f4f9','#e7e1ef','#d4b9da','#c994c7','#df65b0','#e7298a',
              '#ce1256','#980043','#67001f', '#340014'),
	    c('#fff5f0','#fee0d2','#fcbba1','#fc9272','#fb6a4a','#ef3b2c',
              '#cb181d','#a50f15','#67000d', '#340202'),
	    c('#fff5f0','#f0f0f0','#d9d9d9','#bdbdbd','#969696','#737373',
              '#525252','#252525','#131313', '#340202'))

dcols = list(rev(c('#35000f', '#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#f0f0f0',
                   '#d1e5f0','#92c5de','#4393c3','#2166ac','#053061', '#031531')),
	     c('#400131', '#8e0152','#c51b7d','#de77ae','#f1b6da','#fde0ef','#f7f7f7',
               '#e6f5d0','#b8e186','#7fbc41','#4d9221','#276419', '#153210'),
	     rev(c('#301503', '#543005','#8c510a','#bf812d','#dfc27d','#f6e8c3',
                   '#f5f5f5','#c7eae5','#80cdc1','#35978f','#01665e','#003c30', '#001c15')),
	     rev(c('#501a04', '#7f3b08','#b35806','#e08214','#fdb863','#fee0b6','#f7f7f7',
                   '#d8daeb','#b2abd2','#8073ac','#542788','#2d004b', '#19002b')),
	     rev(c('#730013', '#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf',
                   '#e0f3f8','#abd9e9','#74add1','#4575b4','#313695', '#100019')),
	     rev(c('#20002a', '#40004b','#762a83','#9970ab','#c2a5cf','#e7d4e8','#f7f7f7',
                   '#d9f0d3','#a6dba0','#5aae61','#1b7837','#00441b', '#002210')))

levels_controls_regions = list(Canada = list(levels = c(0, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50),
                                     dlevels = c(-2, -1, -0.5, -0.2,-0.1, 0,
                                                 0.1, 0.2, 0.5, 1, 2)),
                       Greece = list(levels = c(0, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50),
                                     dlevels = c(-2, -1, -0.5, -0.2, -0.1,0, 
                                                 0.1, 0.2, 0.5, 1, 2)),
                       NW_Amazon = list(levels = c(0, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50),
                                        dlevels = c(-30, -10, -3, -1, -0.3, -0.1, 0,
                                                    0.1, 0.3, 1, 3, 10, 30)))

levels_BA_regions = list(Canada = list(levels_BA = c(0, 0.001, 0.002, 0.004, 0.006, 
                                                    0.008, 0.01, 0.02, 0.04, 0.06, 0.08),
                                      dlevels_BA = c(-0.7, -0.3, -0.1, -0.05, -0.01, 0,
                                                     0.01, 0.05, 0.1, 0.3, 0.7)),
                        Greece = list(levels_BA = c(0, 0.001, 0.002, 0.005, 0.01, 0.02, 
                                                    0.04, 0.06, 0.08, 0.1),
                                      dlevels_BA = c(-0.3, -0.1, -0.05, -0.02, -0.01, 0,
                                                     0.01, 0.02, 0.05, 0.1, 0.3)),
                        NW_Amazon = list(levels_BA = c(0, 0.0001, 0.001, 0.002, 0.005, 0.01, 
                                                       0.02, 0.05, 0.1),
                                         dlevels_BA = c(-1, -0.5, -0.2, -0.1, -0.05, -0.02, 0,
                                                        0.02, 0.05, 0.1, 0.2, 0.5, 1)))


dlevels = c(-2, -1, -0.5, -0.2,-0.1, 0, 0.1, 0.2, 0.5, 1, 2)*100
levels_controls_regions = list(Canada = list(levels = c(0, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50),
                                     dlevels = dlevels),
                       Greece = list(levels = c(0, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50),
                                     dlevels =dlevels),
                       NW_Amazon = list(levels = c(0, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50),
                                        dlevels = dlevels))

dlevels_BA = c(-0.5, -0.2, -0.1, -0.05, -0.01, 0,  0.01, 0.05, 0.1, 0.2, 0.5) * 100
levels_BA_regions = list(Canada = list(levels_BA = c(0, 0.001, 0.002, 0.004, 0.006, 
                                                    0.008, 0.01, 0.02, 0.04, 0.06, 0.08),
                                      dlevels_BA = dlevels_BA),
                        Greece = list(levels_BA = c(0, 0.001, 0.002, 0.005, 0.01, 0.02, 
                                                    0.04, 0.06, 0.08, 0.1),
                                      dlevels_BA = dlevels_BA),
                        NW_Amazon = list(levels_BA = c(0, 0.0001, 0.001, 0.002, 0.005, 0.01, 
                                                       0.02, 0.05, 0.1),
                                         dlevels_BA = dlevels_BA))
levels_s2n = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

simple_control_cols = TRUE
control_cols = c('#004D40', '#FFC107', '#D81B60','#1E88E5')

bg_col = "white"
fg_col = "black"

cut_results <- function(x, breaks) {
	out = x
	out[] = 0
	for (i in 1:length(breaks)) 
		out[x > breaks[i]] = i
	return(out)
}

split_to_day <- function(x) 
	#approx(1:length(x), x, seq(1, length(x), by = 1/30))[[2]]
	predict(smooth.spline(1:length(x), x), 
		    seq(0.5, length(x)+0.5, by = 1/30))[[2]]

find_clim_av <- function(x) 
	 sapply(1:360, function(day) mean(x[seq(day, length(x), by = 360)]))

logit <- function(x) {
    x = 0.01/2.02+x/1.01
    log(x/(1-x))
}
logistic <- function(x) 1/(1+exp(-x))
 
calc_dat <- function(dir, file, do_last_year = FALSE, do_anom = FALSE, ...) {   
        print(dir)
        sort_clim = NULL
	if (length(file) > 1) {
            if (file[1] == "residual") {
                dat = lapply(file[-1], function(fl) (as.matrix(read.csv(paste0(dir, fl), 
                                                    stringsAsFactors = FALSE))))
                dat[[1]][is.na(dat[[1]])] = 0
                sort_clim = dat[[2]]
                sort_clim = apply(sort_clim, 1, split_to_day)
                sort_clim = apply(sort_clim, 2, find_clim_av)
                
                dat = 100-100*sqrt((dat[[1]] - dat[[2]])^2)/dat[[2]]
                
            } else {
                dat = lapply(file, function(fl)
                                    log(as.matrix(read.csv(paste0(dir, fl), 
                                                    stringsAsFactors = FALSE))))
                dat = exp(do.call('+', dat))*100
            }
        } else {
            dat = as.matrix(read.csv(paste0(dir, file), stringsAsFactors = FALSE))*100
        }
	idat = apply(dat, 1, split_to_day)
	clim = apply(idat, 2, find_clim_av)
        
	if (do_anom) {
            print(file)
	    if (length(file) == 1 && file == files[[1]]) {
                #clim0 <<- log(tail(idat, 360)) - log(clim)
                clim = 100*(tail(idat, 360) - clim)/clim
            } else if (file[1] == "residual") {
                clim = tail(idat, 360)
            } else {
                BA = as.matrix(read.csv(paste0(dir, files[1]), stringsAsFactors = FALSE))*100
                idat_BA = apply(BA, 1, split_to_day)
	        clim_BA = apply(idat_BA, 2, find_clim_av)

                clim = ((tail(idat, 360) * clim_BA/clim) - clim_BA)    
                clim = 100*clim /abs(tail(idat_BA, 360) - clim_BA)
                #clim[clim > 400] = 400
                #clim[clim < -400] = -400
                
                #clim = (tail(idat, 360) - clim)/clim
            }
        
	} else {
            if (do_last_year) clim = tail(idat, 360)
	}
        clim[is.na(clim)] = 0
        return(clim)
        
}

plot_for_region <- function(region, region_name, dir, xlim, levels_controls_region, levels_BA_region) {
    
    calc_datR <- function(...) {
        out = calc_dat(dir = dir, ...)        
        out = out[seq(xlim[1] + 15, xlim[2], 30),]
        return(out)
    }
    add_xaxis <- function(id = 3, side = 1, ...) {
        axis(at = seq(15, 360-15, 30), labels = month.abb, 
             side = side, pos = par("usr")[id], ...)
        axis(at = c(0, 360), labels = c('', ''), 
             side = side, pos = par("usr")[id], xpd = FALSE, ...)
    }

    
    #png(paste0("figs/control_stripes-3-", region, ".png"),
    #    height = 7, width = 20, units = 'in', res = 300)

    #par(mfcol = c(6, 3), mar = c(0, 3.5, 0, 0), oma = c(2,2, 3.5, 8))
    #cdats_clim = mapply(calc_datR, files, SIMPLIFY = FALSE)

    #cdats_last = mapply(calc_datR, files, MoreArgs = list(do_last_year = TRUE))
    cdats_anom = mapply(calc_datR, head(files, -1), head(dcols, -1), 
                        MoreArgs = list(do_anom = TRUE), SIMPLIFY = FALSE)
    conf = calc_datR(tail(files, 1)[[1]], do_last_year = TRUE)
    
    xlim[1] = xlim[1] #+ 10
    xlim[2] = xlim[2] #- 10
    plot(c(-15, 360), range(cdats_anom[[1]]), xlab = '', ylab = '', type = 'n', xlim = xlim, axes = FALSE)
    polygon(c(-9E9, 9E9, 9E9, -9E9), c(-9E9, -9E9, 9E9, 9E9), col = bg_col)
    add_vases <- function(dats, xi, cols, control_col, levels, cdats, scale = NaN, simple_col = FALSE, lwd = 2, wds = 1) {
        cdats = cut_results(cdats, levels)
        if (!is.na(scale)) dats = dats * scale
        add_vase <- function(x0) {
            dat = dats[x0,]; cdat = cdats[x0,]
            xy = hist(dat, 20, plot = FALSE)[c('mids', 'density')]
            xy = spline(xy[[1]], xy[[2]], 200)
            xy[[2]] = 1.5*0.98*xy[[2]]/max(xy[[2]])
            x = 1+(x0-1)*30 + (xi-3) *4.8
            if (xi == 1) x = x - 2
            x = x + xlim[1] + 15
            
            addLine <- function(y, z, col = NULL) {     
                if (is.null(col)) col = cols[z]
                id = which.min(abs(xy[[1]] - y))
                wdth = xy[[2]][id] +wds-1
                if (simple_col) col = control_col
                lines(x + c(-1, 1) * wdth, rep(y, 2), col = col, lwd = lwd)
                #lines(x + c(-1.05, 1.05) * wdth, rep(y, 2), col ='white', lwd = 3)
            }
            
            mapply(addLine, dat, cdat)
            return(x)
        }
        mapply(add_vase, 1:nrow(dats))
        
    }
    range1 = range(cdats_anom[[1]])
    range2 = sapply(cdats_anom[-1], quantile, c(0.1, 0.9))
    range2 = range(range2)
    
    scale = min(range1/range2)
    add_vases(cdats_anom[[1]], 1, fg_col, '', levels_s2n*100, conf, lwd = 3, wds = 1.1)
    mapply(add_vases, cdats_anom[-1], 2:length(cdats_anom), fg_col, 
               fg_col,
               MoreArgs = list(levels = dlevels_BA, cdat = cdats_anom[[1]], scale = scale, 
                               lwd = 3, wds = 1.1, simple_col = simple_control_cols))# 
    for (lwd in seq(3, 0.5, by = -0.5)) {
        add_vases(cdats_anom[[1]], 1, cols[[1]], '', levels_s2n*100, conf, lwd = lwd)
    
    
        mapply(add_vases, cdats_anom[-1], 2:length(cdats_anom), head(dcols[-1], -1), 
               control_cols,
               MoreArgs = list(levels = dlevels_BA, cdat = cdats_anom[[1]], scale = scale, 
                               lwd = lwd, simple_col = simple_control_cols))#  
    }
    axis(2)
    add_xaxis()

    labs4 = unique(signif(seq(range1[1]/scale, range1[2]/scale, length.out = 8), 1))
    axis(4, at = labs4 * scale, labels = labs4)
    lines(c(-100, 460), c(0, 0), col = fg_col, lty = 3)
    mtext(side = 3, line = -1.5, region_name, col = fg_col, adj = 0.01, font = 2)    
    lapply(seq(30, 360, 30), function(x) lines(c(x, x), c(-9E9, 9E9), col = fg_col, lwt = 2, lwd = 0.5))
    return(cdats_anom)
}
png("figs/control_vases.png", res = 300, height = 7.2, width = 7.2, units = 'in')
par(mfrow = c(4, 1), mar = c(2,4,0.5, 4))
nrt_clims = mapply(plot_for_region, regions, Anom_titles,  dirs, xlims, levels_controls_regions, levels_BA_regions, SIMPLIFY = FALSE)
save(nrt_clims, file = "temp2/nrt_clims.Rd")
mtext('Burned Area extent anomaly (%)', outer = TRUE, side = 2, adj = 1-3/8, line = -1.5)  
mtext('Control influnce (%)', outer = TRUE, side = 4, adj = 1-3/8, line = -1.5)  

text(y = -72, x = 185, 'BA', srt = 90, col = fg_col, adj = 0)
text(y = -72, x = 191, 'Fuel', srt = 90, col = fg_col, adj = 0)
text(y = -72, x = 191+4.8, 'Dryness', srt = 90, col = fg_col, adj = 0)
text(y = -72, x = 191+4.8*2, 'Weather', srt = 90, col = fg_col, adj = 0)
text(y = -72, x = 191+4.8*3, 'Human', srt = 90, col = fg_col, adj = 0)


plot(c(0, 1), c(0, 1), axes = FALSE, type = 'n', xlab = '', ylab = '')

legendColBar(yy = c(0.0, 0.5), xx = c(0.5, 0.75), cols = cols[[1]], limits = levels_s2n*100, add = TRUE, transpose = TRUE, extend_max = FALSE, extend_min = FALSE, maxLab = 100, minLab = 0)
mtext(side = 3, cex = 0.8, 'Variations in burned area explained by model (%)', line = -1, adj = 0.1)
legend(x = 0.55, y = 0.725, horiz = TRUE, col = control_cols, c('Fuel', 'Dryness', 'Weather', 'Human'), pch = 19)
dev.off()
# Joesy cols:
#dcols = c('#2c124c', '#66479a', '#aca2d0', '#e1e2ed', 
#          '#fae5c8', '#f6b15e', '#c56a28', '#813d1b')
#cols = c('#fff7ec', '#fee8c8', '#fdd49e', '#fdbb84', 
#         '#fc8d59', '#ef6548', '#d7301f', '#b30000', '#7f0000')

