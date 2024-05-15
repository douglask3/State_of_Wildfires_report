graphics.off()
library(stringr) 
source("../LPX_equil/libs/legendColBar.r")
source("libs/find_levels.r")
source("../rasterextrafuns/rasterPlotFunctions/R/make_col_vector.r")
source("../rasterextrafuns/rasterPlotFunctions/R/mtext.units.r")
regions = c('Canada', 'Greece', 'NW_Amazon')
pcs = c(95, 90, 95)
Anom_titles = c('Canada', 'Greece', 'South American Domain')
xlims = list(c(3, 9)*30, c(6, 9)*30, c(6, 12)*30)

dirs = paste0("outputs/ConFire_", regions, "-nrt-tuning10/figs/_12-frac_points_0.5-baseline-control_TS/pc-", pcs, "/")

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
                                     dlevels = c(-5, -2, -1, -0.5, -0.2,-0.1, 
                                                 0.1, 0.2, 0.5, 1, 2, 5)),
                       Greece = list(levels = c(0, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50),
                                     dlevels = c(-2, -1, -0.5, -0.2, -0.1,0, 
                                                 0.1, 0.2, 0.5, 1, 2)),
                       NW_Amazon = list(levels = c(0, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50),
                                        dlevels = c(-30, -10, -3, -1, -0.3, -0.1, 
                                                    0.1, 0.3, 1, 3, 10, 30)))

levels_BA_regions = list(Canada = list(levels_BA = c(0, 0.001, 0.002, 0.004, 0.006, 
                                                    0.008, 0.01, 0.02, 0.04, 0.06, 0.08),
                                      dlevels_BA = c(-0.7, -0.3, -0.1, -0.05, -0.01, 
                                                     0.01, 0.05, 0.1, 0.3, 0.7)),
                        Greece = list(levels_BA = c(0, 0.001, 0.002, 0.005, 0.01, 0.02, 
                                                    0.04, 0.06, 0.08, 0.1),
                                      dlevels_BA = c(-0.3, -0.1, -0.05, -0.02, -0.01, 
                                                     0.01, 0.02, 0.05, 0.1, 0.3)),
                        NW_Amazon = list(levels_BA = c(0, 0.0001, 0.001, 0.002, 0.005, 0.01, 
                                                       0.02, 0.05, 0.1),
                                         dlevels_BA = c(-1, -0.5, -0.2, -0.1, -0.05, -0.02, 
                                                        0.02, 0.05, 0.1, 0.2, 0.5, 1)))
levels_s2n = c(0, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4)
			 
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



plot_strpes <- function(cdat, xlim, levels, cols, addAxis = TRUE, sort_clim = NULL) {
        cdat0 = cdat
	cdat = cut_results(cdat, levels)
        cdat[cdat == 0] = 1
	plot(c(1, 360), ncol(cdat) * c(-1, 1), type = 'n', axes = FALSE, xlab = '', ylab = '', 
             yaxs = 'i', xlim = xlim)
        
	add_Day <- function(day) {
		add_enemble <- function(ens) {
			x = day + 0.5 * c(-1, 1)
			y = c(ens, ens)
			lines(x, y, col = cols[z[ens]], lwd = 2)
			lines(x, -y, col = cols[z[ens]], lwd = 2)
		}
                if (is.null(sort_clim))  {
                    z = rev(sort(cdat[day,]))
                } else {
                    z = rev(cdat[day, sort.int(sort_clim[day,], index.return=TRUE)[[2]]])
                }
		lapply(1:ncol(cdat), add_enemble)
	}

	sapply(xlim[1]:xlim[2], add_Day)
	if (addAxis) {
            axis(at = seq(15, 360-15, 30), labels = month.abb, side = 1, pos = 0)
	    axis(at = c(0, 360), labels = c('', ''), side = 1, pos = 0, xpd = FALSE)
        }
}

logit <- function(x) log(x/(1+x))
logistic <- function(x) 1/(1+exp(-x))

plot_cols <- function(Anom_title, dir, xlim,
                      file, cols, levels = NULL, do_last_year = FALSE, do_anom = FALSE, 
                      addLegend = TRUE, ...) {   
        
        sort_clim = NULL
	if (length(file) > 1) {
            if (file[1] == "residual") {
                dat = lapply(file[-1], function(fl) (as.matrix(read.csv(paste0(dir, fl), 
                                                    stringsAsFactors = FALSE))))
                dat[[1]][is.na(dat[[1]])] = 0
                sort_clim = dat[[2]]
                sort_clim = apply(sort_clim, 1, split_to_day)
                sort_clim = apply(sort_clim, 2, find_clim_av)
                dat = sqrt((dat[[1]] - dat[[2]])^2)/dat[[2]]
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
	    if (length(file) == 1 && file == files[[1]]) {
                #clim0 <<- log(tail(idat, 360)) - log(clim)
                clim = 100*(tail(idat, 360) - clim)
            } else if (file[1] == "residual") {
                clim = tail(idat, 360)
            } else {
                clim = 100*(tail(idat, 360) - clim)/clim
            }
	    #levels = quantile(abs(clim), seq(0, 1, length.out = ceiling(length(cols))/2))[-1]
	    #levels = c(rev(-levels), levels)
            if (file[1] == "residual")
                extend_min = FALSE
            else
                extend_min = TRUE
            minLab = ''
	} else {
            if (do_last_year) clim = tail(idat, 360)
	    #levels =  head(quantile(dat, seq(0, 1, length.out = length(cols) + 1))[-1], -1)
            extend_min = FALSE
            minLab = 0
            clim[clim < 0] = 0
	}
        clim[is.na(clim)] = 0
        if (is.null(levels)) levels = find_levels_n(clim, 9, TRUE)
        levels = unique(signif(levels), 2)
        cols =  make_col_vector(cols, ncols = length(levels) + 1)
	if (file == "points-standard-Moisture.csv") brows = TRUE else brows = FALSE
	plot_strpes(clim, xlim, levels, cols, ..., sort_clim = sort_clim)

        if (length(file) > 1) {
            if (!do_anom && !do_last_year) {
                if (file[1] == "residual") txt = 'Signal:Noise' else txt = 'Human'
                mtext(side = 2, line = 0, txt)
            }
        } else {
            if (file == files[1]) {
                if (do_anom) mtext(side = 3, line = 2, Anom_title, font = 2)
                    else if (do_last_year) mtext(side = 3, line = 2, '2023', font = 2)
                    else mtext(side = 3, line = 2, 'Climatology', font = 2)
            }
            if (!do_anom && !do_last_year) {
                if (grepl('Control', file)) {
                    txt = 'Burnt Area'
                } else {
                    txt = gsub('-', ' ', gsub('.csv', '', gsub("points-", "", file)))
                    
                    txt = str_to_title(txt)
                    txt = gsub('Moisture', 'Dryness', txt)
                    txt = gsub('Standard', '', txt)
                    #txt = gsub(' ', '\n', txt)
                }
                mtext(side = 2, line = 0, txt)
            }
        }
        if (addLegend) 
            legendColBar(xlim[2] + diff(xlim) * c(0.08333, .125), ncol(clim) * c(-0.9, 0.9), 
                         cols, levels, 
                         extend_max = TRUE,  extend_min = extend_min, minLab = minLab, 
                         add = TRUE, xtext_pos_scale= 1, oneSideLabels = FALSE)
        return(levels)
        
}
plot_for_region <- function(region, Anom_title, dir, xlim, levels_controls_region, levels_BA_region) {
    
    plt_colsR <- function(...)
        plot_cols(Anom_title = Anom_title, dir = dir, xlim = xlim, ...)
    
    add_xaxis <- function(id = 3, side = 1, ...) {
        axis(at = seq(15, 360-15, 30), labels = month.abb, 
             side = side, pos = par("usr")[id], ...)
        axis(at = c(0, 360), labels = c('', ''), 
             side = side, pos = par("usr")[id], xpd = FALSE, ...)
    }

    if (T) {
        png(paste0("figs/control_stripes-2-", region, ".png"),
            height = 7, width = 12, units = 'in', res = 300)

            par(mfcol = c(6, 3), mar = c(0, 3.5, 0, 0), oma = c(2,2, 3.5, 4))
            levels = mapply(plt_colsR, files, cols, SIMPLIFY = FALSE)
            mapply(plt_colsR, files, cols, levels = levels, 
                   MoreArgs = list(do_last_year = TRUE))
            mapply(plt_colsR, head(files, -1), head(dcols, -1), MoreArgs = list(do_anom = TRUE))
            plt_colsR(tail(files, 1)[[1]], tail(cols, 1)[[1]], do_last_year = TRUE)  
        dev.off()
    }
    levels_BA = levels_BA_region[[1]]
    dlevels_BA = levels_BA_region[[2]]
    levels = levels_controls_region[[1]]
    dlevels = levels_controls_region[[2]]
    if (T) {
        png(paste0("figs/control_stripes-2-JoeyStyle", region, ".png"), 
            height = 7, width = 12, units = 'in', res = 300)

            par(mfcol = c(6, 3), mar = c(0, 3.5, 0, 0), oma = c(2,2, 3.5, 4))
        
        for (do_last_year in c(FALSE, TRUE)) {
            plt_colsR(files[[1]], cols = cols[[1]], do_last_year = do_last_year,
                      levels = levels_BA, addAxis = FALSE)
            add_xaxis(4, side = 3)
        
            mapply(plt_colsR, head(files[-1], -1), cols = head(cols[-1], -1), 
                   MoreArgs = list(levels = levels, 
                                   do_last_year = do_last_year, addAxis = FALSE)) 
    
            plt_colsR(tail(files, 1)[[1]], cols = tail(cols, 1)[[1]],levels = levels_s2n, 
                      do_last_year = do_last_year, do_anom = FALSE , addAxis = FALSE)
            add_xaxis(3)
        }

        plt_colsR(files[[1]], levels = dlevels_BA, cols = dcols[[1]],
                  do_anom = TRUE , addAxis = FALSE)
        add_xaxis(4, side = 3)
        mapply(plt_colsR, head(files[-1], -1), cols = head(dcols[-1], -1), 
               MoreArgs = list(levels = dlevels, 
               do_anom = TRUE , addAxis = FALSE))
    
        plt_colsR(tail(files, 1)[[1]], cols = tail(cols, 1)[[1]],levels = levels_s2n, 
                 do_last_year = do_last_year, do_anom = FALSE , addAxis = FALSE)
        add_xaxis(3)
        
        dev.off()
    }
}

mapply(plot_for_region, regions,Anom_titles,  dirs, xlims, levels_controls_regions, levels_BA_regions)
# Joesy cols:
#dcols = c('#2c124c', '#66479a', '#aca2d0', '#e1e2ed', 
#          '#fae5c8', '#f6b15e', '#c56a28', '#813d1b')
#cols = c('#fff7ec', '#fee8c8', '#fdd49e', '#fdbb84', 
#         '#fc8d59', '#ef6548', '#d7301f', '#b30000', '#7f0000')

