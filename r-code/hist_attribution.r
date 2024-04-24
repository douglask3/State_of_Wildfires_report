library(terra)

factual_dir = 'outputs/ConFire_Greece-tuning6/samples/_15-frac_points_0.5/factual-/control/'
counter_dir = 'outputs/ConFire_Greece-tuning6/samples/_15-frac_points_0.5/counterfactual-/control/'

obs_file = 'data/data/driving_data/Greece/isimp3a/obsclim/GSWP3-W5E5/period_2000_2019/burnt_area-2000-2023.nc'
month = 8


obs = rast(obs_file)
obs = obs[[seq(month, nlyr(obs), by = 12)]]
tobs = substr(time(obs), 1, 7)

for_pc <- function(pc) {
    findXX <- function(i) {
        r = obs[[i]]
        cells = which(r[] >= quantile(r[], pc, na.rm = TRUE))
        return(list(cells, sum(cellSize(r)[cells] * r[cells])))
    }
    cells_pc = lapply(1:nlyr(obs), findXX)
    areas_pc = log(sapply(cells_pc, function(i) i[[2]]) + 0.00000000001)
    
    openDats <- function(name, dir, dir2 = NULL) {
        tfile = paste0('temp/hist_attribution/', name, '-', pc, '.csv')
        if (file.exists(tfile)) return(read.csv(tfile, stringsAsFactors = FALSE)[,-1])
        print(tfile)
        files = list.files(dir, pattern = 'pred', full.names = TRUE)
        if (!is.null(dir2)) files2 = list.files(dir2, pattern = 'pred', full.names = TRUE)
        openDat <- function(i) {
            openFile <- function(file) {
                dat = rast(file)
                
                index = lapply(tobs, function(i) which(substr(time(dat), 1, 7) == i))
                obs_index = which(sapply(index, length) > 0)
                mod_index = unlist(index)
                mod = dat[[mod_index]]
                modArea = cellSize(mod[[1]])
                grab_values <- function(i) {
                    sum(mod[[i]][cells_pc[[i]][[1]]] * modArea[cells_pc[[i]][[1]]], na.rm = T)
                }
                mod = log(sapply(1:nlyr(mod), grab_values) + 0.00000000001)
                return(list(mod, obs_index))
            }
    
            mod = openFile(files[i])
            obs_index = mod[[2]]
            mod = mod[[1]]
            if (is.null(dir2)) {
                mod0 = mod
            } else {
                mod0 = openFile(files2[i])[[1]]
            }
            
            #mod0 = mod0 - mean(mod)
            
            #mod0 = mod0 * sd(areapc[obs_index])/sd(mod)
            #mod0 = mod0 + mean(areapc[obs_index])
            return(mod0)
        }
        dat = sapply(1:length(files), openDat)
        write.csv(dat, tfile)
        return(dat)
    }
    
    bias_correct <- function(a, b, obs) {
        mnR = mean(a, na.rm = TRUE)/mean(b, na.rm = TRUE)
        a = a - mean(a)
        a = a * sd(obs)/sd(b)
        a = a + mean(event)*mnR
    }
    factual = unlist(openDats('Greece_factual2', factual_dir))
    counter = unlist(openDats('Greece_counter2', counter_dir))
    event = tail(areas_pc, 1)
    #counter = bias_correct(counter, factual, areas_pc)
    #factual = bias_correct(factual, factual, areas_pc)
    
    return(mean(factual > event)/mean(counter > event))
}

pcs = seq(0.5, 1, 0.05)
rrs = sapply(pcs, for_pc)

