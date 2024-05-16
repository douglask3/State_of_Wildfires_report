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
speedy = TRUE

dir = paste0('outputs/ConFire_', region, '-nrt-tuning10/samples/_12-frac_points_0.5/baseline-/')

run_names = c('Burnt area' = 'control', 'Fuel' = 'Standard_0', 'Moisture' = 'Standard_1',    
         'Weather' = 'Standard_2',
         'Ignitions' = 'Standard_3', 'Suppression' = 'Standard_4', 'Snow' = 'Standard_5')

control_groups = list(1, 2, 3, 4, c(5, 6))

levels = c(0.01, 0.03, 0.1, 0.3)
nlevs = 4
cols = c('#00FF00', '#0000BB', '#FF0000')


make_col_map <- function(col1, col2) make_col_vector(c(col1, col2), ncols = 3)[2]

select_col <- function(i, j, k = 1, l = 1) 
    (c(paste0('#', number_to_hex(i, j), number_to_hex(j, k), number_to_hex(i,k)), l))

find_superLevel <- function(x) {
    if (any(is.na(x))) return(NaN)
    sum((x-1) * nlevs^(0:3)) + 1
}

open_runs  <- function(cntrl_ID) {
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

#runs = lapply(control_groups, open_runs)

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
    
    out = layer.apply(c(0.3439, 0.6561), for_quantile)
    out = writeRaster(out, file = tfile, overwrite = TRUE)
}

control_qu = mapply(get_control_range, controls, control_groups[-1], 
                    MoreArgs = list(1:2, diff_))
#control_qu = mapply(get_control_range, controls, control_groups[-1], 
#                    MoreArgs = list(3:4, root_square_diff))

#control_qu = mapply(get_control_range, controls, control_groups[-1], 
#                    MoreArgs = list(1, function(i) i))

levels = find_levels_n(do.call(addLayer, control_qu), nlevs-2, TRUE, zeroAt0 = TRUE)
levels = sort(c(levels, 0))
nlevs = length(levels) + 1

plot_uncertain_corners <- function(control_i, txt) {
    control_sample = rep(2, 4)
    control_sample[control_i] = 1
    map = do.call(addLayer, mapply(function(r, i) r[[i]], control_qu, control_sample))
    cmap = map#raster::disaggregate(map, fact = 10, method = 'bilinear')   
    
    cmap = cut_results(cmap, levels)
    #browser()
    cmap = calc(cmap, find_superLevel)

    index = 1:nlevs
    pcols = dots = c()
    for (l in index) for (k in index) for (j in index) for (i in index) {
        pcols = c(pcols, select_col(i,j,k,l)[1])
        dots = c(dots, select_col(i,j,k,l)[2])
    }
    
    plotStandardMap(cmap, pcols, limits = NULL, readyCut = TRUE, speedy = speedy)
    mtext(paste('Max.', txt, 'limitation'), side = 3, font = 2)
    cmap = map
    cmap = cut_results(cmap, levels)
    cmap = calc(cmap, find_superLevel)
    pnts = cbind(xyFromCell(cmap, which(!is.na(cmap[]))), cmap[!is.na(cmap)])
    #browser()
    points(pnts[,1], pnts[,2], cex = (as.numeric(dots)[pnts[,3]]-1)/2, pch = 19, col = 'white')
}

layout(rbind(1:2, 3:4, 5))
par(mar = c(0, 0, 1.5, 0))
mapply(plot_uncertain_corners, 1:4, c('Fuel', 'Moisture', 'Weather', 'Human'))

add_legend <- function(levels, cmap = NULL) {
    nlevs = length(levels) + 1
    if (!is.null(cmap)) cmap = lapply(cmap, cut_results, levels)


    leg_cube = 1:nlevs
    leg_cube = rep(leg_cube, length.out = nlevs^4)
    for (i in 2:4) 
        leg_cube = cbind(leg_cube, rep(leg_cube, each = 2^(i-1), , length.out = nlevs^4))
    
    plot_leg_square <- function(i, j, k = 1, l = 1) {
        x = i + (k - 1) * (nlevs + 1)
        y = j + (l - 1) * (nlevs + 1)

        col = select_col(i, j, k, l)
        if (is.null(cmap)) {
            polygon(x + c(-1, 0, 0, -1, -1), y + c(-1, -1, 0, 0, -1), col = col[1])
            points(x-0.5, y-0.5, col = 'white', pch = 19, cex = (as.numeric(col[2])-1)/4)
        } else {
            #in_bin <- function(r, x) 
                #quantile(unlist(layer.apply(r, function(rl) mean(rl[] == x, na.rm = TRUE))), c(0.1, 0.9))
                #range(unlist(layer.apply(r, function(rl) mean(rl[] == x, na.rm = TRUE))))
            
            test = mapply(function(r, x) r == x, cmap, c(i, j, k, l))
            val = sapply(1:nlayers(test[[1]]), function(ly) mean(all(layer.apply(test, function(r) r[[ly]]))[], na.rm = TRUE))
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
        
        if (j == 1) {
            text(x, -0.06 * ascale, levels[i], xpd = NA)
            lines(c(x, x),c(0, -0.02*ascale))
        }
        if (i == 1) {
            text(-0.03 * ascale, y, levels[j], xpd = NA, adj = 1)
            lines(c(0, -0.02*ascale), c(y, y))
        }
        if (i == 1 && j == 1 && l == 1) {
            text(x+ nlevs, -0.18 * ascale, levels[k], xpd = NA)
            
            lines(c(x-1, x + nlevs), rep(-0.12 *ascale, 2), xpd = NA)
            lines(c(x + nlevs, x + nlevs), c(-0.12 *ascale, -0.14*ascale), xpd = NA)
        }       
        if (i == 1 && j == 1 && k == 1) {
            text(-0.09 * ascale, y, levels[l], xpd = TRUE, adj = 1)
            
            lines(rep(-0.05 *ascale, 2), c(y, y + nlevs), xpd = NA)
            lines(c(-0.05 *ascale, -0.07*ascale), c(y, y), xpd = NA)
        }    
        if (i == 1 && j == 1 && k == 1 && l == 1) {
            text(nlevs/2, -0.09 * ascale, 'Fuel', xpd = T)
            text(nlevs^2/2, -0.21 * ascale, 'Weather', xpd = T)
       }
    }
    par(mar = c(6, 6, 0, 0))
    plot(c(0, (nlevs+1)^2), c(0, (nlevs+1)^2), type = 'n', xlab = '', ylab = '', axes = FALSE)
    
    lapply(1:nlevs, function(i) lapply(1:nlevs, function(j) 
        lapply(1:nlevs, function(k) lapply(1:nlevs, function(l) plot_leg_square(i,j, k, l)))))
}

add_legend(levels, control_qu)

