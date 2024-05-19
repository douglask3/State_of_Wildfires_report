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
nlevs = 6
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
        tfile = paste(c('temp2/plot_nrt_map_4ways_agree/mn-', mnths, '-', 
                        region, run_names[ctr], metric, '.nc'), collapse = '-')
        if (file.exists(tfile)) return(brick(tfile))
    
        dats = layer.apply(control, function(i) FUN(i[[metric]]))  
        out = addLayer(mean(dats), mean(dats > 0))
        out = writeRaster(out, file = tfile, overwrite = TRUE)
    }
    
    control_increase = mapply(get_control_range, controls, control_groups[-1], 
                        MoreArgs = list(1:2, diff_))

    plot_uncertain_corners <- function(cmap, dx = -0.25, dy = -0.25, cols, direction = 1) {
        test = cmap[[1]] < 0
        if (direction == 2) test = !test
        cmap[test] = NaN
        cmap[[1]] = cut_results(cmap[[1]], levels)
        cmap[[2]] = cut_results(cmap[[2]], c(0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99))/9
        
        pcols = make_col_vector(cols[[direction]], ncols = nlevs)
        
        dat = cbind(xyFromCell(cmap, 1:length(cmap[[1]][])), cmap[[1]][], cmap[[2]][])
        dat = dat[!is.na(dat[,3]), ]
        #cols = c(cols[1], 'white', cols[2])
        #dat[,3] = 3
        #plot_point <- function(pnts) {
            #polygon(pnts[1] + c(0, dx, dx, 0, 0), 
            #        pnts[2] + c(0, 0, dy, dy, 0), 
            #        col = cols[pnts[3]], border = NA)
        #}
        #apply(dat, 1, plot_point)
        
        points(dat[,1] + dx/2, dat[,2] + dy/2, pch = 19, col = pcols[dat[,3]], cex = dat[,4])
        return()
    }
    
    if (!is.null(box[[2]])) control_increase = lapply(control_increase, crop, box[[2]])
    pltHeight = nrow(control_increase[[1]])/ncol(control_increase[[1]])
    
    levels = find_levels_n(abs(do.call(addLayer, control_increase)),  nlevs-1, TRUE) 
    heights = c(pltHeight, 1, 1)  
        
    #png(paste(c("figs/limitation_map", region, '-', box[[1]], box[[2]], speedy, ".png"), 
    #          collapse = '-'), 
    #    height = 6*sum(heights), width = 6, units = 'in', res = 300)
        layout(rbind(1, 2, 3), heights = heights)
        par(mar = c(0, 1.5, 1.5, 0))
        extent = extent(control_increase[[1]])
        plot(extent[1:2], extent[3:4], xaxs = 'i', yaxs = 'i', type = 'n', xlab = '', ylab = '')
        
# list(c('black', '#1b9e77', '#66c2a5'),
#                    c('black', '#7570b3', '#8da0cb'),
#                    c('black', '#d95f02', '#fc8d62'),
#                    c('black', '#666666', '#DDDDDD')
#list(c('black', '#009900', '#00ff00'), 
#                   c('black', '#93A200', '#efff00'),
#                   c('black', '#c65a00', '#ff8900'),
#                   c('black', '#333333', '#AAAAAA'))
        
        obs = brick(paste0('data/data/driving_data/', region, '/nrt/period_2013_2023/burnt_area.nc'))
        obs = obs[[seq(box[[1]], nlayers(obs), by = 12)]]
        obs = crop(obs, box[[2]])
        cobs = obs[[nlayers(obs)]] > mean(obs[[1:(nlayers(obs)-1)]])
        cobs[is.na(cobs)] = 0
        coords = xyFromCell(cobs, (1:length(cobs)))
        coords = coords[cobs[] == 1,] 
    
        boxFire <- function(pnts) polygon(pnts[1] + 0.25 * c(-1, 1, 1, -1, -1), 
                                          pnts[2] + 0.25 * c(-1, -1, 1, 1, -1),
                                          border = 'grey', col = '#FFFFDD', lty = 2)
        apply(coords, 1, boxFire)
        

        mapply(plot_uncertain_corners, control_increase, c(-0.25, -0.25, 0.25, 0.25),
                                                         c(-0.25, 0.25, -0.25, 0.25),
              list(list(c('white', '#CCDDAA', '#225522'), c('white', '#EEEEBB', '#666633')),
                   rev(list(c('white', '#EEEEBB', '#666633'), c('white', '#BBCCEE', '#222255'))),
                   list(c('white', '#FFCCCC', '#663333'), c('white', '#CCEEFF', '#225555')),
                   list(c('white', '#DDDDDD', '#555555')), list(c('white', '#DDDDDD', '#555555'))))
        browser()
    dev.off()
}
lapply(boxes, run_for_box)
