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
boxes = list(list(8, c(20, 27, 30.5, 42), 'Greece'))

#region = 'Canada'
#boxes = list(list(5, c(-123, -105, 53, 62), 'Western Shield'),
#             list(6, c(-80, -70, 47, 58), 'Quebec'), 
#             list(7, c(-80, -70, 47, 58), 'Quebec'), 
#             list(9, c(-124, -110, 57, 64), 'Western Shield'))

region = 'NW_Amazon'
boxes = list(list( 9, c(-70, -55, -8, 0), 'Western Amazon'),
             list(10, c(-70, -55, -8, 0), 'Western Amazon'))

dir = paste0('outputs/ConFire_', region, '-nrt-tuning10/samples/_12-frac_points_0.5/baseline-/')

run_names = c('Burnt area' = 'control', 'Fuel' = 'Standard_0', 'Moisture' = 'Standard_1',    
         'Weather' = 'Standard_2',
         'Ignitions' = 'Standard_3', 'Suppression' = 'Standard_4', 'Snow' = 'Standard_5')

control_groups = list(1, 2, 3, 4, c(5, 6))

levels = c(0.01, 0.03, 0.1, 0.3)
nlevs = 6
cols = c('#00FF00', '#0000BB', '#FF0000')


boxes_from_raster_bool <- function(pnts, dx = 0.25, dy = dx, ...)
                polygon(pnts[1] + dx * c(-1, 1, 1, -1, -1), 
                pnts[2] + dy * c(-1, -1, 1, 1, -1),...)

make_col_map <- function(col1, col2) make_col_vector(c(col1, col2), ncols = 3)[2]

legLabs = list(c('Less Fuel', '', 'More Fuel'), 
                c('Drier Fuel', '', 'Wetter Fuel'),
                c('Low Fire Weather', '', 'High Fire Weather'),
                c('Supression', '', 'Human-caused increase'))


conf_levels = c(0.5, 0.75, 0.9, 0.95, 0.99)

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
        tfile = paste(c('temp2/plot_nrt_map_4ways_agree/mn-', mnths, '-', 
                        region, run_names[ctr], metric, '.nc'), collapse = '-')
        if (file.exists(tfile)) return(brick(tfile))
    
        dats = layer.apply(control, function(i) FUN(i[[metric]]))  
        qus = layer.apply(c(0.05, 0.95), function(qu) calc(dats, function(x) quantile(x, qu ,na.rm = TRUE)))
        sig = mean(dats > 0)
        out = qus[[1]]
        out[sig > 0.5] = qus[[2]][sig > 0.5]
        out = addLayer(out, sig)
        out = writeRaster(out, file = tfile, overwrite = TRUE)
    }
    
    control_increase = mapply(get_control_range, controls, control_groups[-1], 
                        MoreArgs = list(1:2, diff_))

    plot_uncertain_corners <- function(cmap, dx = -0.25, dy = -0.25, pcols, direction = 1) {
        if (direction == 2) {
            cmap[[1]] = - cmap[[1]]
            cmap[[2]] = 1 - cmap[[2]]
        }
        test = cmap[[1]] < 0
        #cmap0 = cmap
        cmap[test] = NaN
        cmap[[1]] = cut_results(cmap[[1]], levels)
        cmap[[2]] = cut_results(cmap[[2]], conf_levels)
        #cmap[[2]] = length(conf_levels) + 1
        cmap[[2]] = 0.25 + 0.75 * cmap[[2]]/(length(conf_levels)+1)
        #cmap1 = cmap
        cmap[[2]] = 14* cmap[[2]] / ncol(cmap)
        
        #pcols = make_col_vector(cols, ncols = nlevs)
        
        dat = cbind(xyFromCell(cmap, 1:length(cmap[[1]][])), cmap[[1]][], cmap[[2]][])
        dat = dat[!is.na(dat[,3]), ]
        
        points(dat[,1] + dx/2, dat[,2] + dy/2, pch = 19, col = 'black', cex = dat[,4])
        points(dat[,1] + dx/2, dat[,2] + dy/2, pch = 19, col = pcols[dat[,3]], cex = dat[,4])
        return()
    }
    
    if (!is.null(box[[2]])) control_increase = lapply(control_increase, crop, box[[2]])
    pltHeight = nrow(control_increase[[1]])/ncol(control_increase[[1]])
    
    levels = find_levels_n(abs(do.call(addLayer, control_increase)),  nlevs-1, TRUE) 
    heights = c(0.25, pltHeight, 0.2, 0.15, pltHeight, 0.15, pltHeight, 0.6, 0.1)  
    lmat = cbind(0, rbind(0, 1, 2, 0, 3, 0, 4, 5, 0), 0)
    widths = c(0.25, 1, 0.25)
    png(paste(c("figs/limitation_map-2", region, '-', box[[1]], box[[2]], ".png"), 
              collapse = '-'), 
        height = 2.5*sum(heights), width = 2.5*sum(widths), units = 'in', res = 300)
        layout(lmat, heights = heights, widths = widths)
        par(mar = rep(0, 4))
        extent = extent(control_increase[[1]])
        
        
        obs = brick(paste0('data/data/driving_data/', region, '/nrt/period_2013_2023/burnt_area.nc'))
        obs = obs[[seq(box[[1]], nlayers(obs), by = 12)]]
        obs = crop(obs, box[[2]])
        obs[[1]][is.na(control_increase[[1]][[1]])] = NaN
        cobs = obs[[nlayers(obs)]] > mean(obs[[1:(nlayers(obs)-1)]])
        pobs = 100*(obs[[nlayers(obs)]] - mean(obs[[1:(nlayers(obs)-1)]]))

        newPlot <- function() {
            plot(extent[1:2], extent[3:4], xaxs = 'i', yaxs = 'i', 
                 type = 'n', xlab = '', ylab = '', axes = FALSE)
            contour(is.na(cobs), levels = 0.5, drawlabels=FALSE, add = TRUE)
        }
        newPlot()
        
            
        levels_BA = find_levels_n(pobs, 9, TRUE)
        cols_BA = rev(c('#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7',
                    '#d1e5f0','#92c5de','#4393c3','#2166ac'))
        cols_BA = make_col_vector(cols_BA, ncols = length(levels_BA) +1)

        addPoly <- function(pnts) {
            polygon(pnts[1] + 0.25 * c(-1, 1, 1, -1, -1), pnts[2] + 0.25 * c(-1, -1, 1, 1, -1), col = cols_BA[pnts[3]], border = 'grey', lty = 2)
        }
        vobs = cut_results(pobs, levels_BA)
        vobs = cbind(xyFromCell(vobs, 1:length(vobs)), vobs[])
        vobs = vobs[!is.na(vobs[,3]),] 
        apply(vobs, 1, addPoly)
        contour(is.na(cobs), levels = 0.5, drawlabels=FALSE, add = TRUE)
        axis(2)
        #axis(3)
        lines(extent[c(1, 2, 2, 1, 1)], extent[c(3, 3, 4, 4, 3)])
        mtext(side = 3, line = 0.5, 'Burnt area anomoly (%)')
        legendColBar(c(0.2, 0.7), c(0, 1), 
                     cols_BA, levels_BA, F, T, T, transpose = T, oneSideLabels = F)
                
        cols_controls = list(c('white', '#CCDDAA', '#225522'),
                     c('white', '#BBCCEE', '#222255'),
                     c('white', '#FFCCCC', '#663333'),
                     c('white', '#CCEEFF', '#225555'))
                
        cols_controls = lapply(cols_controls, make_col_vector, ncols = length(levels) + 1)
        
        plot_direction <- function(direction = 1) {
            newPlot()
            title = paste0(c(c('increase', 'decrease')[direction], ' from controls'), collapse = '')
            mtext(title, side = 3, line = 0.5)
            if (direction == 2) cobs = 1 - cobs
            cobs[is.na(cobs)] = -1
            coords = xyFromCell(cobs, (1:length(cobs)))

            apply(coords[cobs[] == 0,], 1, boxes_from_raster_bool, border = '#FFFFFF00', #density = 5 * ncol(cobs), 
                   lty = 2, col = c('#EEEEEE', '#EEEEEE')[direction])
            coords = coords[cobs[] == 1,] 
            
            apply(coords, 1, boxes_from_raster_bool, border = 'grey', #density = 5 * ncol(cobs), 
                   lty = 2)
            mapply(plot_uncertain_corners, control_increase, 
                   c(-0.25, -0.25, 0.25, 0.25), c(-0.25, 0.25, -0.25, 0.25),
                cols_controls, direction = direction)
            axis(2)
            if (direction == 2) axis(1)
            lines(extent[c(1, 2, 2, 1, 1)], extent[c(3, 3, 4, 4, 3)])
        }
        lapply(1:2, plot_direction)
        
        
        legX = seq(0.2, 0.7, length.out = length(cols_controls[[1]]))
        addLegend <- function(y, cols, lab) {
            for (col in list('black', cols)) 
                points(legX, rep(y, length(legX)), col = col, pch = 19, cex = 2, xpd = NA)
            text(x = 0.05, y = y, adj = 1, lab, xpd = NA) 
    
        }
        
        plot(c(0.05, 0.95), c(0, 1.1), axes = FALSE, type = 'n', xlab = '', ylab = '')
        mapply(addLegend, c(0.5, 0.35, 0.20, 0.05), cols_controls, 
               c("Fuel", "Moisture", "Weather", "Human"))

        text(c(legX[1], head(legX, -1) + diff(legX[1:2])/2, tail(legX, 1)), xpd = NA,
             rep(0.6, length(legX) + 1), c('0', levels*100, '+'), srt = 30, adj = 0)
        points(rep(0.85, length(conf_levels)), 
               y = seq(0.05, 0.5, length.out = length(conf_levels)), xpd = NA,
              pch = 19, cex = 2*(0.25 + 0.75*(1:length(conf_levels))/length(conf_levels)))
        text(paste('>', conf_levels), x = 0.95, y = seq(0.05, 0.5, length.out = length(conf_levels)), xpd = NA)
        text('Confidence level', x = 0.9, y = 0.6, xpd = NA)
        
        mtext(side = 3, outer = TRUE, font = 2,
              line = -2, paste(box[[3]], '-', month.name[box[[1]]]))

        mtext(side = 2, outer = TRUE, expression(paste( degree, " Latitude")), 
              xpd = NA, line = -2.75)#, adj = 1-(0.15 + hght*nrow/2)/sum(heights))
        mtext(side = 3, expression(paste( degree, " Longitude")), xpd = NA, line = -3.5)
   dev.off()
}
lapply(boxes, run_for_box)
