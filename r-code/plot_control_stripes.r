graphics.off()
library(stringr) 
source("../LPX_equil/libs/legendColBar.r")
source("libs/find_levels.r")
source("../rasterextrafuns/rasterPlotFunctions/R/make_col_vector.r")
source("../rasterextrafuns/rasterPlotFunctions/R/mtext.units.r")
region = 'Canada'
region = 'Greece'
#region = 'NW_Amazon'

Anom_title = region
#Anom_title = 'South American Domain'
xlim = c(3, 9)*30
xlim = c(6, 9)*30
#xlim = c(6, 12)*30

dir = paste0("outputs/ConFire_", region, "-nrt-tuning10/figs/_12-frac_points_0.5-baseline-control_TS/pc-90/")

files = list("points-Control.csv",
	     "points-standard-Fuel.csv",
	     "points-standard-Moisture.csv",
	     "points-standard-Weather.csv",
	     c("points-standard-Ignition.csv", "points-standard-Suppression.csv"),
             c("residual", "points-Evaluate.csv", "points-Control.csv"))


levels = c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)
cols = list(c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026', '#400013'),
	    c('#ffffe5','#f7fcb9','#d9f0a3','#addd8e','#78c679','#41ab5d','#238443','#006837','#004529', '#002315'),
	    rev(c('#ffffd9','#edf8b1','#c7e9b4','#7fcdbb','#41b6c4','#1d91c0','#225ea8','#253494','#081d58', '#041034')),
	    c('#f7f4f9','#e7e1ef','#d4b9da','#c994c7','#df65b0','#e7298a','#ce1256','#980043','#67001f', '#340014'),
	    c('#fff5f0','#fee0d2','#fcbba1','#fc9272','#fb6a4a','#ef3b2c','#cb181d','#a50f15','#67000d', '#340202'),
	    c('#fff5f0','#f0f0f0','#d9d9d9','#bdbdbd','#969696','#737373','#525252','#252525','#131313', '#340202'))

dlevels = c(-0.05, -0.02, -0.01, -0.005, -0.002, -0.001, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05)
dcols = list(rev(c('#35000f', '#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#f0f0f0',
                   '#d1e5f0','#92c5de','#4393c3','#2166ac','#053061', '#031531')),
	     c('#400131', '#8e0152','#c51b7d','#de77ae','#f1b6da','#fde0ef','#f7f7f7',
               '#e6f5d0','#b8e186','#7fbc41','#4d9221','#276419', '#153210'),
	     rev(c('#301503', '#543005','#8c510a','#bf812d','#dfc27d','#f6e8c3','#f5f5f5','#c7eae5','#80cdc1','#35978f','#01665e','#003c30', '#001c15')),
	     rev(c('#501a04', '#7f3b08','#b35806','#e08214','#fdb863','#fee0b6','#f7f7f7','#d8daeb','#b2abd2','#8073ac','#542788','#2d004b', '#19002b')),
	     rev(c('#730013', '#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695', '#100019')),
	     rev(c('#20002a', '#40004b','#762a83','#9970ab','#c2a5cf','#e7d4e8','#f7f7f7','#d9f0d3','#a6dba0','#5aae61','#1b7837','#00441b', '#002210')))
			 
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



plot_strpes <- function(cdat, levels, cols, addAxis = TRUE, sort_clim = NULL) {
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

plot_cols <- function(file, cols, levels = NULL, do_last_year = FALSE, do_anom = FALSE, 
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
	plot_strpes(clim, levels, cols, ..., sort_clim = sort_clim)

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
            legendColBar(xlim[2] + diff(xlim) * c(0.08333, .125), ncol(clim) * c(-0.9, 0.9), cols, levels, 
                         extend_max = TRUE,  extend_min = extend_min, minLab = minLab, 
                         add = TRUE, xtext_pos_scale= 1, oneSideLabels = FALSE)
        return(levels)
        
}

if (T) {
png(paste0("figs/control_stripes-", region, ".png"), height = 7, width = 12, units = 'in', res = 300)
par(mfcol = c(6, 3), mar = c(0, 2, 0, 0), oma = c(0,2, 2, 3.5))

levels = mapply(plot_cols, files, cols, SIMPLIFY = FALSE)

mapply(plot_cols, files, cols, levels = levels, MoreArgs = list(do_last_year = TRUE))
mapply(plot_cols, files, dcols, MoreArgs = list(do_anom = TRUE))
dev.off()
}

if (T) {
png(paste0("figs/control_stripes-JoeyStyle", region, ".png"), height = 7, width = 12, units = 'in', res = 300)

#lmat = matrix(1:18, nrow = 6)
#lmat = rbind(lmat[1,], c(19, 19, 20), lmat[-1,], c(21, 21, 22))
#layout(lmat, height = c(1, 0.6, rep(1, 6), 0.6))
par(mfcol = c(6, 3), mar = c(0, 3.5, 0, 0), oma = c(2,2, 3.5, 4))
#dcols = c('#2c124c', '#66479a', '#aca2d0', '#e1e2ed', 
#          '#fae5c8', '#f6b15e', '#c56a28', '#813d1b')
#cols = c('#fff7ec', '#fee8c8', '#fdd49e', '#fdbb84', 
#         '#fc8d59', '#ef6548', '#d7301f', '#b30000', '#7f0000')

#levels = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50)
#levels = 1:50

#dlevels = seq(1, 6, by = 0.03)
#dlevels = 100*0.030*(((dlevels)^3.5)/(6^3.5))
#dlevels = unique(signif(dlevels, 1))
#dlevels = c(-rev(dlevels), dlevels)

#levels_BA = unique(signif(10^seq(-5, -1, length.out = 51), 2))

#dlevels_BA = unique(signif(10^seq(-5, -2, length.out = 25), 2))
#dlevels_BA = c(-rev(dlevels_BA), dlevels_BA)

levels = c(0, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50)
dlevels = c(-5, -2, -1, -0.5, -0.2,-0.1, 0.1, 0.2, 0.5, 1, 2, 5)

#levels = c(0, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50)
#dlevels = c(-2, -1, -0.5, -0.2, -0.1,0, 0.1, 0.2, 0.5, 1, 2)

#levels = c(0, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50)
#dlevels = c(-30, -10, -3, -1, -0.3, -0.1, 0.1, 0.3, 1, 3, 10, 30)

levels_BA = c(0, 0.001, 0.002, 0.004, 0.006, 0.008, 0.01, 0.02, 0.04, 0.06, 0.08)
dlevels_BA = c(-0.7, -0.3, -0.1, -0.05, -0.01, 0.01, 0.05, 0.1, 0.3, 0.7)

#levels_BA = c(0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1)
#dlevels_BA = c(-0.3, -0.1, -0.05, -0.02, -0.01, 0.01, 0.02, 0.05, 0.1, 0.3)

#levels_BA = c(0, 0.0001, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1)
#dlevels_BA = c(-0.01, -0.005, -0.002, -0.001, -0.0005, -0.0002, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01)*100

levels_s2n = c(0, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4)
add_xaxis <- function(id = 3, side = 1, ...) {
    axis(at = seq(15, 360-15, 30), labels = month.abb, side = side, pos = par("usr")[id], ...)
    axis(at = c(0, 360), labels = c('', ''), side = side, pos = par("usr")[id], xpd = FALSE, ...)
}

for (do_last_year in c(FALSE, TRUE)) {
    plot_cols(files[[1]], cols = cols[[1]], do_last_year = do_last_year,
              levels = levels_BA, addAxis = FALSE) #], addLegend = FALSE
    add_xaxis(4, side = 3)

    mapply(plot_cols, head(files[-1], -1), cols = head(cols[-1], -1), 
           MoreArgs = list(levels = levels, 
                           do_last_year = do_last_year, addAxis = FALSE)) #, addLegend = FALSE
    
    plot_cols(tail(files, 1)[[1]], cols = tail(cols, 1)[[1]],levels = levels_s2n, 
           do_last_year = do_last_year, do_anom = FALSE , addAxis = FALSE)
    add_xaxis(3)
}


#mapply(plot_cols, files, 
#       MoreArgs = list(levels = levels, cols = cols, 
#                       addLegend = FALSE, do_last_year = TRUE, addAxis = FALSE))



plot_cols(files[[1]], levels = dlevels_BA, cols = dcols[[1]],
          do_anom = TRUE , addAxis = FALSE) #addLegend = FALSE, 
add_xaxis(4, side = 3)
mapply(plot_cols, head(files[-1], -1), cols = head(dcols[-1], -1), MoreArgs = list(levels = dlevels, 
                                         do_anom = TRUE , addAxis = FALSE)) #addLegend = FALSE
    
plot_cols(tail(files, 1)[[1]], cols = tail(cols, 1)[[1]],levels = levels_s2n, 
           do_last_year = do_last_year, do_anom = FALSE , addAxis = FALSE)
add_xaxis(3)
#legendColBar(c(0.3, 0.7), c(0.0, 1), cols = cols, limits = levels_BA, extend_max = TRUE, extend_min = FALSE, minLab = 0, transpose = TRUE, oneSideLabels = FALSE)
#legendColBar(c(0.3, 0.7), c(0.0, 1), cols = dcols, limits = dlevels_BA, extend_max = TRUE, extend_min = TRUE, transpose = TRUE, oneSideLabels = FALSE)

#legendColBar(c(0.3, 0.7), c(0.0, 1), cols = cols, limits = levels, extend_max = TRUE, extend_min = FALSE, minLab = 0, transpose = TRUE, oneSideLabels = FALSE)
#legendColBar(c(0.3, 0.7), c(0.0, 1), cols = dcols, limits = dlevels, extend_max = TRUE, extend_min = TRUE, transpose = TRUE, oneSideLabels = FALSE)
dev.off()
}
