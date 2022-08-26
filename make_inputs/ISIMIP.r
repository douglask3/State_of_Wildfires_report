source("libs/filename.noPath.r")
library(ncdf4)
library(raster)
source("../gitProjectExtras/gitBasedProjects/R/sourceAllLibs.r")
sourceAllLibs("../gitProjectExtras/gitBasedProjects/R/")
sourceAllLibs("../rasterextrafuns/rasterExtras/R/")
source("libs/process_jules_file.r")
source("libs/writeRaster.Standard.r")
options(error=recover)
abcd <- function() source("make_inputs/ISIMIP.r")
countriesMap = raster('../ConFIRE_ISIMIP/data/countries.nc')
ckey = read.csv("../ConFIRE_ISIMIP/data/countries_key.csv")[,2]
genVarID = "genVar-C-cover2-FfsoilM-soilC-VDP"
#
overwrite_outputs = FALSE

histDir = "/hpc//data/d00/hadea/jules_output/u-cc669_isimip3a_es/GSWP3-W5E5/"
soilFile = '/hpc//data/d00/hadea/isimip3a/jules_ancils/qrparm.soil.latlon_fixed.nc'
popDir = "/hpc//data/d00/hadea/isimip3a/InputData/socioeconomic/pop/histsoc/"
lightFile = "/hpc/data/d00/hadea/isimip2b/ancils/lightning/lightning_cloud2ground.nc"
runYrLen = 10

Syears = seq(1880, 2010, by = runYrLen)
names = paste('historic_TS', Syears, Syears+runYrLen-1, sep = '_')
years = lapply(Syears, function(i) i:(i+runYrLen-1))
names(years) = names

dirs = rep(histDir, length(years))
names(dirs) = names

sm_sat = raster(soilFile, varname = "sm_sat")


countries = c(Global = NA)

fileIDs = c(cover = "ilamb", soilM = "ilamb", cveg = "ilamb", cs_gb = "ilamb", spres = "ilamb",
            precip = "ilamb", humid = "ilamb", tas = "ilamb")

varnames =  c(cover = "frac", soilM = "smc_tot", cveg = "cv", cs_gb = "cs_gb", spres = "pstar",
              precip = "precip", humid = "q1p5m_gb", tas = "t1p5m_gb")

models = c("obsclim", "counterclim")

temp_dir = '/data/users/dkelley/ConFIRE_ISIMIP_temp/-makeISIMIP3ins'
temp_dir_mem = '/data/users/dkelley/ConFIRE_ISIMIP_temp/memSafe/'
out_dir  = 'isimip3a/driving_data/GSWP3-W5E5/'

coverTypes = list(trees = c(1:7), totalVeg = c(1:13), crop = c(10, 12), pas = c(11, 13))

lightn = brick(lightFile)
popFiles = list.files(popDir, full.names = TRUE)
popFiles = popFiles[grepl('.nc', popFiles)]
popDens  = stack(popFiles)
extent(lightn) = extent(popDens)

makeDir(out_dir)
try(memSafeFile.remove())
memSafeFile.initialise(temp_dir_mem)
makeDat <- function(id, dir, years, out_dir, mask,  extent, country) {
    print(id)
    print(country)
    #years = c(years_out, tail(years_out, 1) + 1)
    forModel <- function(mod) {
        
        print(id)
        print(mod)
        out_dirM = paste0(out_dir , '/', id)
        
        makeDir(out_dirM)
        makeDir(paste0(out_dir, '/inference_data/'))
        out_dirM = paste(out_dir,  id, mod, '', sep = '/')
        makeDir(out_dirM)
        
        genVarFile = paste0(out_dirM, '/', genVarID, '.Rd')
        if(file.exists(genVarFile)) return()
        tfile0 = paste0(c(temp_dir, country, id, mod, range(years)), collapse = '-')
        
        #dir = paste0(dir, '/', mod, '/')
        files = list.files(dir, full.names = TRUE)
        files = files[grepl(mod, files)]
        ## select years
        files = files[apply(sapply(years, function(i) grepl(i, files)), 1, any)]
        files = files0 = files[substr(files, nchar(files)-2, nchar(files))=='.nc']
        if (length(files) == 0) browser()
        openVar <- function(fileID, vname) {
            tfileC = paste(tfile0 , fileID, vname, '-masked-corrected.Rd', sep = '-')
            print(tfileC)
            
            if (file.exists(tfileC)) {
                load(tfileC)
            } else {
                print(tfileC)
                
                files = files[grepl(fileID, files)] 
                
                if (substr(id, 1,3) == "RCP") 
                    files = files[grepl(paste0('rcp', substr(id, 4, 6)), files)]
                
                processSaveFile <- function(file, yr) {
                    dat =  process.jules.file(file, NULL, vname, mask = mask)
                    
                    if (!is.list(dat)) {
                        dat = writeRaster(dat, paste(tfile0 , fileID, vname, yr,
                                                      '.nc', sep = '-') ,overwrite=TRUE)
                    } else {
                        tfile = paste(tfile0 , fileID, vname, yr, 1:length(dat),
                                      '.nc', sep = '-')
                        dat = mapply(writeRaster, dat, tfile, overwrite = TRUE)
                    }
                    return(dat)
                }
                dat = mapply(processSaveFile, files, years, SIMPLIFY = FALSE)
                
                save(dat, file =  tfileC)
            }
            gc()
            return(dat)
        }
        dats = mapply(openVar, fileIDs, varnames)
        
        cover = dats[, 'cover']
       
        makeCover <- function(ty) {
            print(ty)
            group <- function(i) {
                ctfile = paste(c(tfile0, 'coverSummed', 
                                 strsplit(filename.noPath(i[[1]], TRUE), 'frac')[[1]][2],
                                 ty, '.nc'), collapse = '-')
                if (file.exists(ctfile)) return(brick(ctfile))
                cv = i[ty]
                out = cv[[1]]
                
                for (i in cv[-1]) 
                    out = out + i               
                out = writeRaster(out, ctfile, overwrite = TRUE)
                return(out)
            }
            coverTy = layer.apply(cover, group)
        }
        covers = lapply(coverTypes, makeCover)
        
        #soilMs = lapply(1:4, function(ly) layer.apply(dats[, 'soilM'], function(i) i[[ly]]))
        
        int2Month <- function(i, r) {
            out = layer.apply(r, function(j) j[[i]])
            out[[rep(1:nlayers(out), each = 12)]]
        }
        
        writeOut <- function(dat, name) {
            file = paste0(out_dirM,  name, '.nc')
            if (!overwrite_outputs && file.exists(file)) return(brick(file))
            if (name == "soilM")
                dat = layer.apply(dat, function(i) i/(3001*raster::resample(sm_sat, i)))
            print(file)
            date = as.Date(paste(rep(years, each = 12), 1:12, 15, sep = '-'))
            dat = setZ(dat, date, name = "Date")
            writeRaster.Standard(dat, file)
        }
          
        var2layerWrite <- function(name, file = name) {
            out = layer.apply(dats[, name], function(i) i)   
            out = writeOut(out, file)
        }
        out = list(covers = mapply(writeOut, covers, names(coverTypes)),
                   #soilMs = mapply(writeOut, soilMs, paste0('soilM', 1:4)),    
                   #soilCs = mapply(writeOut, soilCLayers, paste0('soilC', 1:4)),
                   precip = var2layerWrite("precip"),
                   humid  = var2layerWrite("humid"),
                   tas    = var2layerWrite("tas"),
                   cveg   = var2layerWrite("cveg"),
                   soilM  = var2layerWrite("soilM"),
                   spres  = var2layerWrite("spres"),
                   csoil  = var2layerWrite("cs_gb", "csoil"))       
        relativeHumidity <- function(humid, pres, temp) {
            #temp = temp - 273.15
            out = 0.00263 * pres * humid / exp((17.67 * (temp-273.15))/(temp - 29.65))
        }

        SVP <- function(temp) {
            temp = temp - 273.16
            610.78 * exp(temp / (temp +237.3) * 17.2694)
        }
        rh = relativeHumidity(out[['humid']], out[['spres']], out[['tas']])
        svp = SVP(out[['tas']])
        vpd = svp*(1-rh)
        vpd[vpd < 0] = 0
        writeOut(rh, "rhumid")
        writeOut(svp, "svp")
        writeOut(vpd, "vpd")

        popDensi = popDens[[c(min(years) - 1, years, max(years) + 1) - 1849]]
        popDensi = raster::crop(popDensi, vpd[[1]])
        interpolate_popdens <- function(mn) {
            mn = mn + 6
            mn = mn/12
            wd = mn - floor(mn)
            yr = ceiling(mn)
           return(popDensi[[yr]] * (1-wd) + popDensi[[yr + 1]] * wd)
        }
        
        popDens_int = layer.apply(1:nlayers(vpd), interpolate_popdens)
        writeOut(popDens_int, "popDens")

        lightn = raster::resample(lightn, vpd[[1]])
        lightn = do.call(addLayer, rep(c(lightn), length(years)))
        writeOut(lightn, "lightn")
        
        save(out, file = genVarFile)
        closeAllConnections()
        gc()
    }
    lapply(models, forModel)
    gc()
}

makeCountry <- function(country, cID) {
    if (is.na(cID)) {
        mask = NULL
        extent = NULL
    } else {
        mask = countriesMap!=which(ckey == cID)
        extent = as.vector(apply(xyFromCell(mask, which(mask[] == 0)), 2, range))
        mask = raster::crop(mask,extent)
    }
    out_dirC = paste0(out_dir, '/', country, '/')
    makeDir(out_dirC)
    #extent = as.vector(apply(xyFromCell(mask, mask[] ==1), 2, range))
    mapply(makeDat, names(dirs), dirs, years,
           MoreArgs = list(out_dir = out_dirC, mask = mask, extent = extent,
                           country = country))
}
countriesMap[is.na(countriesMap)] = 0
mapply(makeCountry, names(countries), countries)
memSafeFile.remove()
