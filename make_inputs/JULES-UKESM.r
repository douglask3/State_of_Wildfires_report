source("libs/filename.noPath.r")
library(ncdf4)
library(raster)
source("../gitProjectExtras/gitBasedProjects/R/sourceAllLibs.r")
sourceAllLibs("../gitProjectExtras/gitBasedProjects/R/")
sourceAllLibs("../rasterextrafuns/rasterExtras/R/")
source("libs/process_jules_file.r")
source("libs/writeRaster.Standard.r")
source("libs/convert_pacific_centric_2_regular.r")
options(error=recover)
abcd <- function() source("make_inputs/JULES-UKESM.r")
countriesMap = raster('../ConFIRE_ISIMIP/data/countries.nc')
ckey = read.csv("../ConFIRE_ISIMIP/data/countries_key.csv")[,2]
genVarID = "genVar"

histDir = "/hpc/data/d05/cburton/jules_output/u-bi607_HIST/"
futDir_ssp5  = "/hpc/data/d05/cburton/jules_output/u-cd136_SSP5/"
futDir_ssp3  = "/hpc/data/d05/cburton/jules_output/u-cd136_SSP3/"
futDir_ssp1  = "/hpc/data/d05/cburton/jules_output/u-cd136_SSP1/"
dirs = list(historic_TS_short = histDir,
            historic_TS  = histDir,
            ssp1_TS    = futDir_ssp1 ,
            ssp3_TS    = futDir_ssp3 ,
            ssp5_TS    = futDir_ssp5)

years = list(historic_TS_short = 2000:2014,
             historic_TS  = 1861:2014,
             ssp1_TS    = 2015:2100,
             ssp3_TS    = 2015:2100,
             ssp5_TS    = 2015:2100)

countries = c(#Global = NA, 
              Indonesia = 'Indonesia',  
              Brazil = 'Brazil', Madagascar = 'Madagascar', Russia = 'Russia', Australia = 'Australia',  USA = 'United States of America', UK = 'United Kingdom')

fileIDs = c(cover = "ilamb", soilMtot = "ilamb", burnt_area = "ilamb",
            cveg = "ilamb", cs_gb = "ilamb", gpp = "ilamb", npp = "ilamb",
            precip = "ilamb", humid = "ilamb", tas = "ilamb")

varnames =  c(cover = "frac", soilMtot = "smc_tot", burnt_area = "burnt_area_gb", 
              cveg = "cv", cs_gb = "cs_gb", gpp = "gpp_gb", npp = "npp_n_gb",
              precip = "precip", humid = "q1p5m_gb", tas = "t1p5m_gb")


temp_dir = '/data/users/dkelley/ConFIRE_ISIMIP_temp/-makeISIMIPins_JULES_UKESM'
temp_dir_mem = '/data/users/dkelley/ConFIRE_ISIMIP_temp/memSafe/'
out_dir  = '/data/users/dkelley/ConFIRE_ISIMIP/inputs_JULES_UKESM/'

coverTypes = list(trees = c(1:7), totalVeg = c(1:13), crop = c(10, 12), pas = c(11, 13))
makeDir(out_dir)
try(memSafeFile.remove())
memSafeFile.initialise(temp_dir_mem)

reg = raster('../jules_workshed/be397.nc')
countriesMap = raster::resample(convert_regular_2_pacific_centric(countriesMap), reg)

makeDat <- function(id, dir, years_out, out_dir, mask,  extent, country) {
    
    print(country)
    years = c(years_out[1] - 1, years_out)#, tail(years_out, 1) + 1)
    
    print(id)
    out_dirM = paste0(out_dir , '/', id, '/')
    makeDir(out_dirM)
        
    genVarFile = paste0(out_dirM, '/', genVarID, '.Rd')
    if(file.exists(genVarFile)) return()
    tfile0 = paste0(c(temp_dir, country, id, range(years)), collapse = '-')
    dir = paste0(dir, '/')
    files = list.files(dir, full.names = TRUE)
     
    ## select years
    files = files[apply(sapply(years, function(i) grepl(i, files)), 1, any)]
    files = files[substr(files, nchar(files)-2, nchar(files))=='.nc']
       
    openVar <- function(fileID, vname) {
        tfileC = paste(tfile0 , fileID, vname, '-masked-corrected.Rd', sep = '-')
        print(tfileC)
            
        if (file.exists(tfileC)) {
            load(tfileC)
        } else {
            print(tfileC)
            
            files = files[grepl(fileID, files)]  
            
            
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
        
    cover = dats[-1, 'cover']
        
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
        
    #tfile = paste(tfile0, 'Ycovers', '.Rd', sep = '-')
    #if (file.exists(tfile) && FALSE) load(tfile)
    #else {
        covers = lapply(coverTypes, makeCover)
          
    #    save(covers, file = tfile)
    #}
    #tfile = paste(tfile0, 'soils', '.Rd', sep = '-')
    #if (file.exists(tfile)) load(tfile)
    #else {
        #soilM_top    = layer.apply(dats[, 'soilM'], function(i) i[[1]])
        #soilM_bottom = layer.apply(dats[-1, 'soilM'], function(i) i[[1]])

        #st = 2:(nlayers(soilM_top)-11)
        #ed = 13:nlayers(soilM_top)
        #soil12 = mapply(function(s, e) sum(soilM_top[[s:e]]), st, ed)
        #soil12 = layer.apply(soil12, function(i) i)
        #soilM_top =  soilM_top[[-(1:12)]]
        #soil12 = soilM_top/soil12
       #     save(soilM_top, soilM_bottom, soil12, file = tfile)
       # }
    packVar <- function(vname) 
        layer.apply(dats[-1, vname], function(i) i)
    burnt_area = packVar("burnt_area")
    cveg       = packVar('cveg')
    cs_gb      = packVar('cs_gb')
    gpp        = packVar('gpp')
    npp        = packVar('npp')
    soilMtot   = packVar('soilMtot')
    precip     = packVar('precip')
    humid      = packVar('humid')
    tas        = packVar('tas')        
    
    writeOut <- function(dat, name) {
        file = paste0(out_dirM,  '/', name, '.nc')
        print(file)
        #dat = dat[[-1]]
        nl = 12*floor(nlayers(dat)/12)
        dat = dat[[1:nl]]
        
        writeRaster.Standard(dat, file)
        return(file)
    }
    out = c(mapply(writeOut, covers, names(coverTypes)),
           writeOut(burnt_area, 'burnt_area'),
           writeOut(soilMtot, 'soilMtot'),
           writeOut(cveg, 'cveg'),
           writeOut(gpp, 'gpp'),
           writeOut(npp, 'npp'),
           writeOut(cs_gb, 'cs_gb'),
           writeOut(precip, 'precip'),
           writeOut(humid, 'humid'),
           writeOut(tas, 'tas')) 
    save(out, file = genVarFile)
    closeAllConnections()
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
