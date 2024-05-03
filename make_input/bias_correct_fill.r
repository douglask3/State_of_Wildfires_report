library(terra)

region = 'Canada'
experiments = c('historical', 'ssp126', 'ssp370', 'ssp585')
periods = c('period_1994_2014', 'period_2015_2099/', 'period_2015_2099/', 'period_2015_2099/')
models = c('GFDL-ESM4', 'IPSL-CM6A-LR', 'MPI-ESM1-2-HR', 'MRI-ESM2-0', 'UKESM1-0-LL')

dir = paste0('../data/data/driving_data/', region, sep = '')
subDir = mapply(function(i, j) paste(i, models, j, '', sep = '/'), experiments, periods) 
dirs = c(paste0(dir, '/isimp3a/obsclim/GSWP3-W5E5/period_2000_2019/'),
         paste0(dir, '/isimp3a/obsclim/GSWP3-W5E5/period_1901_1920/'),
         paste0(dir, '/isimp3a/counterclim/GSWP3-W5E5/period_2000_2019/'),
         paste0(dir, '/isimp3a/counterclim/GSWP3-W5E5/period_1901_1920/'),
         paste0(dir, '/isimp3b/', subDir))

debiass = c('debiased_tree_cover_jules-es.nc', 'debiased_nonetree_cover_jules-es.nc')
orginals = c('tree_cover_jules-es.nc', 'tree_cover_jules-es.nc')
out_files = c('filled_debiased_tree_cover_jules-es.nc', 'filled_debiased_nonetree_cover_jules-es.nc')
out_all_file = 'filled_debiased_vegCover_jules-es.nc'

forDir <- function(dir) {
    print(dir)
    forVar <- function(debias, orginal, out_file) {
        deb = deb0 = rast(paste0(dir, debias))
        org = rast(paste0(dir, orginal))
        for (i in 1:nlyr(deb)) {
            print(i)
            mask = which(is.na(values(deb[[i]])) & !is.na(values(org[[i]])))
            deb[[i]][mask] = as.matrix(org[[i]][mask])
        }       
        
        time(deb) = time(org)
        out_file = paste0(dir, out_file)
        print(out_file)
        writeCDF(deb, out_file, overwrite = TRUE)
        return(deb)
    }
    out = mapply(forVar, debiass, orginals, out_files)
    outAll = out[[1]] + out[[2]]#100*(out[[2]]/100 + (out[[1]]/100) * (1-out[[2]]/100))
    
    writeCDF(outAll, paste0(dir, out_all_file), overwrite = TRUE)
}

lapply(dirs, forDir)
