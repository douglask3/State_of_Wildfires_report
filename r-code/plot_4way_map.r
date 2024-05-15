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

region = 'Greece'
mnths = c(8)


dir = paste0('outputs/ConFire_', region, '-nrt-tuning10/samples/_12-frac_points_0.5/baseline-/')

run_names = c('Burnt area' = 'control', 'Fuel' = 'Standard_0', 'Moisture' = 'Standard_1',    
         'Weather' = 'Standard_2',
         'Ignitions' = 'Standard_3', 'Suppression' = 'Standard_4', 'Snow' = 'Standard_5')

control_groups = list(1, 2, 3, 4, c(5, 6))

levels = c(0.01, 0.05, 0.1, 0.5)
nlevs = 5
cols = c('#1b9e77','#d95f02','#7570b3')
cols = c('#00FF00', '#0000BB', '#FF0000')
#cols = sapply(cols, function(col) make_col_vector(c(make_col_vector(c('black', col), ncols = 5)[2], col, 'white'), ncols = nlevs))
#cols = sapply(cols, function(col) (make_col_vector(c(col, 'white'), ncols = nlevs)))

make_col_map <- function(col1, col2) make_col_vector(c(col1, col2), ncols = 3)[2]
# Function to convert a number to hexadecimal
number_to_hex <- function(i, j = 1) {
    
    i = (i-1)/(nlevs-1)
    j = (j-1)/(nlevs-1)
    number =  (i+j)/2#sqrt(i^2 + j^2)
    #number = min(i, j)
    if (number > 1) number = 1
    number = round(number*255)

  if (number < 0 || number > 255) {
    stop("Number must be between 0 and 255")
  }
  hex <- sprintf("%02X", number)
  return(hex)
}


select_col <- function(i, j, k = 1, l = 1) 
    (c(paste0('#', number_to_hex(i, j), number_to_hex(j, k), number_to_hex(i,k)), l))
    #c(make_col_map(make_col_map(cols[i, 1], cols[j, 2]), cols[k,3]), l)


open_runs  <- function(cntrl_ID) {
    runs = run_names[cntrl_ID]

    listFiles <- function(run) {
        files = list.files(paste0(dir, run, '/'), full.names = TRUE)
        files = files[grepl('pred', files)]
    }
        
    files = sapply(runs, listFiles)
    
    openFiles <- function(files) {
        
        tfile = paste(c('temp2/plot_nrt_map_4ways/mn-', mnths, '-', gsub('/', '-', files)),
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

metric = 1
get_control_range <- function(control, ctr,  metric) {
    tfile = paste(c('temp2/plot_nrt_map_4ways/mn-', mnths, '-', region, run_names[ctr], 
                     metric, '.nc'),
                  collapse = '-')
    if (file.exists(tfile)) return(brick(tfile))
    dats = layer.apply(control, function(i) i[[metric]])
    for_quantile <- function(qu) {
        print(qu)
        calc(dats, function(x) quantile(x, qu,na.rm = TRUE))
    }
    
    out = layer.apply(c(0.3439, 0.6561), for_quantile)
    out = writeRaster(out, file = tfile, overwrite = TRUE)
}

control_qu = mapply(get_control_range, controls, control_groups[-1], 1)

map = addLayer(control_qu[[1]][[1]], control_qu[[2]][[2]], control_qu[[3]][[1]], control_qu[[4]][[1]])

levels = quantile(map[], head(seq(0, 1, length.out = nlevs+1)[-1], -1), na.rm = TRUE)

#map[[1]][] = 1
#map[[2]][] = 1
#map[[3]][] = 0
#map[[4]][] = 1

cmap = cmap0 = cut_results(map, levels)

find_superLevel <- function(x) {
    if (any(is.na(x))) return(NaN)
    #if (any(is.na(x))) browser()
    #x[1] + sum((x-1) * (nlevs^(0:3)-1))
    sum((x-1) * nlevs^(0:3)) + 1
}
cmap = calc(cmap, find_superLevel)

plevels = pcols = c()
index = 1:nlevs
for (l in index) for (k in index) for (j in index) for (i in index) {   
   plevels =  c(plevels, find_superLevel(c(i, j, k, l)))
   pcols = c(pcols, select_col(i,j,k,l)[1])
}
 plotStandardMap(cmap, pcols, limits = NULL, readyCut = TRUE)
browser()
#normalise_controls <- function(ens_id) {
#    browser()
#}

#controls = lapply(1:length(controls[[1]]), normalise_controls)



leg_cube = 1:nlevs
leg_cube = rep(leg_cube, length.out = nlevs^4)
for (i in 2:4) 
    leg_cube = cbind(leg_cube, rep(leg_cube, each = 2^(i-1), , length.out = nlevs^4))



plot_leg_square <- function(i, j, k = 1, l = 1) {
    col = select_col(i, j, k, l)
    x = i + (k - 1) * (nlevs + 1)
    y = j + (l - 1) * (nlevs + 1)
    polygon(x + c(-1, 0, 0, -1, -1), y + c(-1, -1, 0, 0, -1), col = col[1])
    points(x-0.5, y-0.5, col = 'white', pch = 19, cex = (as.numeric(col[2])-1)/4)
    #browser()
}
dev.new()
plot(c(0, (nlevs+1)^2), c(0, (nlevs+1)^2), type = 'n', xlab = '', ylab = '')

lapply(1:nlevs, function(i) lapply(1:nlevs, function(j) 
       lapply(1:nlevs, function(k) lapply(1:nlevs, function(l) plot_leg_square(i,j, k, l)))))

