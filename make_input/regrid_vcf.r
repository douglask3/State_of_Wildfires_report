library(terra)
library(gdalUtils)

source("libs/add_date_time.r")

path = '../temp/MODIS/'
temp_path = '../temp/regrid_vcf/'
output_path = '../data/data/driving_data/isimp3a/'
newproj = "+proj=longlat +datum=WGS84"
example_file = rast('../data/wwf_terr_ecos_0p5.nc')
extent = c(-170.0, -30.0, 30, 85)
area_name = 'Canada'
variables = rev(c("tree" = 1, "nontree" = 2, "nonveg" = 3))

eg_raster = rast(example_file)
eg_raster = crop(eg_raster, extent)

test_if_overlap <- function(r1, r2) {
    ext1 = ext(r1)
    ext2 = ext(r2)
    extc = intersect(ext1, ext2)
    return((extc[2] > extc[1]) & (extc[4] > extc[3]))
}

temp_path = paste0(temp_path, '/', area_name, '/') 
dir.create(temp_path, recursive = TRUE) 
files = rev(list.files(path, full.name = TRUE))


regrid_file <- function(file, band = 1, name = 'tree') {
    print(file)
    out_info0 = gsub('/', '', strsplit(file, path)[[1]][2], fixed = TRUE)
    temp_path = paste0(temp_path, '/', name, '/')
    dir.create(temp_path, recursive = TRUE) 
    
    out_info = paste0(temp_path, '-', band, '-', gsub('.hdf', '.txt', out_info0, fixed = TRUE))
    if (file.exists(out_info)) {
        info = read.table(out_info)[1,1]
        if (info == 'NoOverlap') return(NULL)
        return(rast(as.character(info)))
    }
    
    dat = rast(file, band)
    test = project(aggregate(dat, 100), newproj)

   
    overlap = test_if_overlap(test, eg_raster)
    if (!overlap) {
        writeLines('NoOverlap', out_info)
        return(NULL)
    }
   
    dat = terra::project(dat, newproj)

    find_area <- function(dat, ...) {
        #dat = aggregate(dat, 4)
        #dat = aggregate(dat, 0.5/rev(res(dat)), ...)
        #if ((ext(test)[2] - ext(test)[1])>180)  browser()
        dat = resample(dat, eg_raster, ...)
    }
    veg_cover = land_cover = dat
    veg_cover[veg_cover>150] = 0
    veg_cover = find_area(veg_cover)
    if (all(is.na(veg_cover[]))) {
        writeLines('NoOverlap', out_info)
        return(NULL)
    }
    land_cover[land_cover < 150] = 1
    land_cover[land_cover > 150] = 0
    land_cover = find_area(land_cover, 'sum')
    out = c(veg_cover, land_cover)
    
    nc_out = paste0(temp_path, sub('.hdf', '.nc',out_info0, fixed = TRUE))
    writeCDF(out, nc_out, overwrite = TRUE)
    writeLines(nc_out, out_info)
    
}

years = sapply(files, function(file) substr(strsplit(file, 'MOD44B.A')[[1]][2], 1, 4))
mn = 3
day = 6
forVegType <- function(band, name) {
    output_path = paste0(output_path, area_name, '_extended/')
    dir.create(output_path, recursive = TRUE) 
    output_fname = paste0(output_path, '/', name, '-raw.nc')
    #if (file.exists(output_fname)) return(rast(output_fname))
    
    dats = lapply(files, regrid_file, band, name)
    years = as.numeric(years)
    
    test = !sapply(dats, is.null)
    dats = dats[test]
    years = years[test]
    
    yearI = sort(unique(years))
    eg_raster[] = 0
    output = areaR = rep(eg_raster, max(2, length(yearI)))

    for (i in 1:length(dats)) {
        print(i)
        dat = dats[[i]]
        dat[is.na(dat)] = 0
        whichY = which(yearI == years[[i]])
    
        output[[whichY]] = output[[whichY]] + dat[[1]] * dat[[2]]
    
        areaR[[whichY]] = areaR[[whichY]] + dat[[2]]
    }
    
    cover = output/areaR
    writeCDF(cover, output_fname, overwrite=TRUE)
    mn = rep(mn, length(yearI))
    day = rep(day, length(yearI))
    day = day-(4*as.integer(yearI/4) == yearI)
    add_date_to_file(output_fname, yearI, mn, day, name, paste0(name, '-raw'), unit = '%')    
}

outs = mapply(forVegType, variables, names(variables))
