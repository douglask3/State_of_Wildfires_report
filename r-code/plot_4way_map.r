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
source("libs/number_to_hex.r")
graphics.off()

region = 'Greece'
mnths = c(8)
speedy = FALSE
pltHeight = 0.5*8/7

region = 'Canada'
mnths = c(6)
speedy = TRUE
pltHeight = 0.33

boxes = list(list(6, c(-80, -70, 47, 58)), list(5, c(-128, -105, 53, 62)),
                  list(7, c(-80, -70, 47, 58)), list(9, c(-125, -115, 57, 64)))

#region = 'NW_Amazon'
#mnths = c(10)
#speedy = FALSE
#pltHeight = 0.33

dir = paste0('outputs/ConFire_', region, '-nrt-tuning10/samples/_12-frac_points_0.5/baseline-/')

run_names = c('Burnt area' = 'control', 'Fuel' = 'Standard_0', 'Moisture' = 'Standard_1',    
         'Weather' = 'Standard_2',
         'Ignitions' = 'Standard_3', 'Suppression' = 'Standard_4', 'Snow' = 'Standard_5')

control_groups = list(1, 2, 3, 4, c(5, 6))

levels = c(0.01, 0.03, 0.1, 0.3)
nlevs = 3
cols = c('#00FF00', '#0000BB', '#FF0000')


make_col_map <- function(col1, col2) make_col_vector(c(col1, col2), ncols = 3)[2]

legLabs = list(c('Less Fuel', '', 'More Fuel'), 
                c('Drier Fuel', '', 'Wetter Fuel'),
                c('Low Fire Weather', '', 'High Fire Weather'),
                c('Supression', '', 'Human-caused increase'))

select_col <- function(i, j, k = 1, l = 1) {
    lab = mapply(function(x, y) x[y], head(legLabs, -1), c(i, j, k)) 
    lab = lab[lab != ""]
    c(paste0('#', number_to_hex(k, k), number_to_hex(i, i), number_to_hex(j, j)), l, paste0(lab, collapse = ', '))
    #c(paste0('#', number_to_hex(i, j), number_to_hex(j, k), number_to_hex(i,k)), l)
}
find_superLevel <- function(x) {
    if (any(is.na(x))) return(NaN)
    sum((x-1) * nlevs^(0:3)) + 1
}

open_runs  <- function(cntrl_ID, mnths) {
    runs = run_names[cntrl_ID]

    listFiles <- function(run) {
        files = list.files(paste0(dir, run, '/'), full.names = TRUE)
        files = files[grepl('pred', files)]
    }
        
    files = sapply(runs, listFiles)
    
    openFiles <- function(files) {      
        tfile = paste(c('temp2/plot_nrt_map_4ways/mn-', mnths, region, '-', gsub('/', '-', files)),
                      collapse = '-')
        
        if (file.exists(tfile)) return(brick(tfile))
        print(files)
        dats0 =  lapply(files, brick)
        dat0 = dats0[[1]]
        if (length(dats0) > 1)
            for (r in dats0[-1]) dat0 = dat0 * r
        
        clim = mean(layer.apply(mnths, function(mn) mean(dat0[[seq(mn, nlayers(dat0), by = 12)]])))
        dat2023 = mean(dat0[[((nlayers(dat0)-11):(nlayers(dat0)))[mnths]]])
        first_half = layer.apply(mnths, function(mn)
                    mean(dat0[[seq(mn, floor(nlayers(dat0)/2), by = 12)]]))
        second_half = layer.apply(mnths, function(mn)
                    mean(dat0[[seq(ceiling(nlayers(dat0)/2), nlayers(dat0), by = 12)]]))
        out = addLayer(clim, dat2023, first_half, second_half)
        out = writeRaster(out, file= tfile, overwrite = TRUE)
        return(out)
    }
    
    outs = apply(files, 1, openFiles) 
    return(outs)
}

run_for_box <- function(box) {
    mnths = box[[1]]
    runs = lapply(control_groups, open_runs, mnths)
    BA = runs[[1]]
    controls = runs[-1]

    root_square_diff <- function(rs) sqrt(sum(rs**2))
    diff_ <- function(rs) rs[[2]] - rs[[1]]

    get_control_range <- function(control, ctr,  metric, FUN) {
        tfile = paste(c('temp2/plot_nrt_map_4ways_qus/mn-', mnths, '-', region, run_names[ctr], 
                        metric, '.nc'),
                      collapse = '-')
        if (file.exists(tfile)) return(brick(tfile))
    
        dats = layer.apply(control, function(i) FUN(i[[metric]]))  
        for_quantile <- function(qu) {
            print(qu)
            calc(dats, function(x) quantile(x, qu,na.rm = TRUE))
        }
    
        out = layer.apply(c(0.3439, 0.5, 0.6561), for_quantile)
        out = writeRaster(out, file = tfile, overwrite = TRUE)
    }
    
    control_qu = mapply(get_control_range, controls, control_groups[-1], 
                        MoreArgs = list(1:2, diff_))

    plot_uncertain_corners <- function(control_i, txt, qus) {
        if (length(qus) == 2) {
            control_sample = rep(qus[2], 4)
            control_sample[control_i] = qus[1]
            map = do.call(addLayer, mapply(function(r, i) r[[i]], control_qu, control_sample))
        } else {
             map = do.call(addLayer, lapply(control_qu, function(r) r[[qus]]))
        }
        cmap = cut_results(map, levels)
        cmap = calc(cmap, find_superLevel)
    
        index = 1:nlevs
        pcols = dots = names = c()
        for (l in index) for (k in index) for (j in index) for (i in index) {
            pcols = c(pcols, select_col(i,j,k,l)[1])
            dots = c(dots, select_col(i,j,k,l)[2])
            names = c(names, select_col(i,j,k,l)[3])
        }
        
        plotStandardMap(cmap, pcols, limits = NULL, readyCut = TRUE, speedy = speedy)
        if (!is.null(txt)) mtext(paste('Max.', txt, 'limitation'), side = 3, font = 2)
        
        pnts = cbind(xyFromCell(cmap, which(!is.na(cmap[]))), cmap[!is.na(cmap)])
        
        dts = nlevs * (as.numeric(dots)/nlevs)^2
        points(pnts[,1], pnts[,2], cex = (25/ncol(cmap)) * (dts[pnts[,3]]-1)/2, 
            pch = 19, col = 'white')
        return(list(unique(cmap[!is.na(cmap)]), pcols, dots, names))
    }
    
    add_legend <- function(levels, cmap = NULL) {
        nlevs = length(levels) + 1
        if (!is.null(cmap)) cmap = lapply(cmap, cut_results, levels)
        
        
        leg_cube = 1:nlevs
        leg_cube = rep(leg_cube, length.out = nlevs^4)
        for (i in 2:4) 
            leg_cube = cbind(leg_cube, rep(leg_cube, each = 2^(i-1), , length.out = nlevs^4))
        
        plot_leg_square <- function(i, j, k = 1, l = 1) {
            print(i*1000 + j*100 + k*10 + l)
            x = i + (k - 1) * (nlevs + 1)
            y = j + (l - 1) * (nlevs + 1)

            col = select_col(i, j, k, l)
            if (is.null(cmap)) {
                polygon(x + c(-1, 0, 0, -1, -1), y + c(-1, -1, 0, 0, -1), col = col[1])
                points(x-0.5, y-0.5, col = 'white', pch = 19, cex = 1.3*(as.numeric(col[2])-1))
            } else {
                test = mapply(function(r, x) r == x, cmap, c(i, j, k, l))
                val = sapply(1:nlayers(test[[1]]), function(ly)
                                mean(all(layer.apply(test, function(r) 
                                    r[[ly]]))[], na.rm = TRUE))
                polygon(x + c(-1, 0, 0, -1, -1), y + c(-1, -1, 0, 0, -1), col = '#FFFFFF00', 
                        lty = 2, border = col)
                if (any(val > 0)) {
                    val = sort(val)
                    vcol = make_col_vector(c("white", "black"), ncol = length(val) + 1)[-1]
                    addShade_poly <- function(v, col)
                        polygon(x + c(-1, 0, 0, -1, -1), y -1 + c(0, 0, v, v, 0), col = col)
                    mapply(addShade_poly, val^(0.5), vcol)
                }        
            }
        
            ascale = diff(par("usr")[1:2])
            p1txt = -0.025 * ascale
            p1tck = -0.01 * ascale
            p1tit = -0.05 * ascale

            p2txt = -0.085 * ascale
            p2tck = -0.07 *ascale
            p2tit = -0.11 * ascale
            
            if (j == 1) {
                text(x, p1txt, levels[i], xpd = NA)
                lines(c(x, x),c(0, p1tck))
            }
            if (i == 1) {
                text(p1txt, y, levels[j], xpd = NA, srt = 90)
                lines(c(0, p1tck), c(y, y))
            }
            if (i == 1 && j == 1 && l == 1) {
                text(x+ nlevs, p2txt, levels[k], xpd = NA)
                
                lines(c(x-1, x + nlevs), rep(p2tck - p1tck, 2), xpd = NA)
                lines(c(x + nlevs, x + nlevs), c(p2tck - p1tck, p2tck), xpd = NA)
            }       
            if (i == 1 && j == 1 && k == 1) {
                text(p2txt, y, levels[l], xpd = TRUE, srt = 90)
                
                lines(rep(p2tck - p1tck, 2), c(y, y + nlevs), xpd = NA)
                lines(c(p2tck - p1tck, p2tck), c(y, y), xpd = NA)
            }    
            if (i == 1 && j == 1 && k == 1 && l == 1) {
                text(nlevs/2  , p1tit , 'Fuel', xpd = T)
                text(nlevs^2/2, p2tit, 'Weather', xpd = T)
           }   
            if (i == 1 && j == 1 && k == 1 && l == 1) {
                text(p1tit, nlevs/2  , 'Moisture', xpd = T, srt = 90)
                text(p2tit, nlevs^2/2, 'Human', xpd = T, srt = 90)
           }
        }
        par(mar = c(6, 6, 0, 0))
        plot(c(0, (nlevs+1)^2-nlevs), c(0, (nlevs+1)^2-nlevs), 
             type = 'n', xlab = '', ylab = '', axes = FALSE)
            
        lapply(1:nlevs, function(i) lapply(1:nlevs, function(j) 
              lapply(1:nlevs, function(k) 
                     lapply(1:nlevs, function(l) plot_leg_square(i,j, k, l)))))
    }
    if (!is.null(box[[2]])) control_qu = lapply(control_qu, crop, box[[2]])
    levels = find_levels_n(do.call(addLayer, control_qu),  nlevs-2, TRUE)
    pltHeight = nrow(control_qu[[1]])/ncol(control_qu[[1]])

    if (is.null(box[[2]])) {
        heights = c(pltHeight, pltHeight, 1, 1)    
        png(paste(c("figs/limitation_map", region, '-', box[[1]], box[[2]], speedy, ".png"), 
                  collape = '-'), 
            height = 3*sum(heights), width = 6, units = 'in', res = 300)
        layout(rbind(1:2, 3:4, 5, 6), heights = heights)
        par(mar = c(0, 1.5, 1.5, 0))
        mapply(plot_uncertain_corners, 1:4, c('Fuel', 'Moisture', 'Weather', 'Human'),
               MoreArgs = list(c(1, 3)))
        add_legend(levels, control_qu)
        add_legend(levels)
    } else {
        heights = c(pltHeight, 1)  
        
        png(paste(c("figs/limitation_map", region, '-', box[[1]], box[[2]], speedy, ".png"), 
                  collapse = '-'), 
            height = 6*sum(heights), width = 6, units = 'in', res = 300)
        layout(rbind(1, 2), heights = heights)
        par(mar = c(0, 1.5, 1.5, 0))
        cmap_pnts = plot_uncertain_corners(1, NULL, 2)

        obs = brick(paste0('data/data/driving_data/', region, '/nrt/period_2013_2023/burnt_area.nc'))
        obs = obs[[seq(box[[1]], nlayers(obs), by = 12)]]
        obs = crop(obs, box[[2]])
        cobs = obs[[nlayers(obs)]] > mean(obs[[1:(nlayers(obs)-1)]])
        cobs[is.na(cobs)] = 0
        cobs = focal(cobs, w=matrix(1, nrow=3, ncol=3), fun=change_cell, pad=TRUE, padValue=0)
        cobs = focal(cobs, w=matrix(1, nrow=3, ncol=3), fun=change_cell, pad=TRUE, padValue=0)
        contour(raster::disaggregate(cobs, fact = 5), add = TRUE)
        
        cmap_pnts[[1]] = sort(cmap_pnts[[1]])
        leg_cols = unique(cmap_pnts[[2]][cmap_pnts[[1]]])
        leg_labs = unique(cmap_pnts[[4]][cmap_pnts[[1]]])
        leg_labs[leg_labs == ''] = 'No anomolous control'
        plot.new()
        legend('top', pch = 19, col = leg_cols, leg_labs, bty = 'n', xpd = NA)
        
    }

    
    dev.off()
}
lapply(boxes, run_for_box)
