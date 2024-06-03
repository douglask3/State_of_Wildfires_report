library(terra)
graphics.off()
factual_dir = 'outputs/ConFire_Greece-tuning9/samples/_13-frac_points_0.5/factual-/'
counter_dir = 'outputs/ConFire_Greece-tuning9/samples/_13-frac_points_0.5/counterfactual-/'

obs_file = 'data/data/driving_data/Greece/isimp3a/obsclim/GSWP3-W5E5/period_2000_2019/burnt_area-2000-2023.nc'
month = 8


obs = rast(obs_file)
obs = obs[[seq(month, nlyr(obs), by = 12)]]
tobs = substr(time(obs), 1, 7)

for_pc <- function(pc, plot = TRUE, subdir = 'control') {
    findXX <- function(i) {
        r = obs[[i]]
        cells = which(r[] >= quantile(r[], pc, na.rm = TRUE))
        return(list(cells, sum(cellSize(r)[cells] * r[cells])))
    }
    cells_pc = lapply(1:nlyr(obs), findXX)
    areas_pc = log(sapply(cells_pc, function(i) i[[2]]) + 0.00000000001)
    
    openDats <- function(name, dir) {
        tfile = paste0('temp/hist_attribution/', name, '-', subdir, '-', pc, '.csv')
        #if (file.exists(tfile)) return(read.csv(tfile, stringsAsFactors = FALSE)[,-1])
        print(tfile)
        files = list.files(paste0(dir, '/', subdir, '/'), pattern = 'pred', full.names = TRUE)
        
        openDat <- function(file) {
            tfile = paste0('temp/hist_attribution/-', pc, gsub('.nc', '.csv', gsub('/', '-', file)))
            if (file.exists(tfile)) return(read.csv(tfile, stringsAsFactors = FALSE)[,-1])
            dat = rast(file)
            print(file)
            index = lapply(tobs, function(i) which(substr(time(dat), 1, 7) == i))
            
            mod_index = unlist(index)
            mod = dat[[mod_index]]
           
            modArea = cellSize(mod[[1]])
            grab_values <- function(i) 
                sum(mod[[i]][cells_pc[[i]][[1]]] * modArea[cells_pc[[i]][[1]]], na.rm = T)
            out = log(sapply(1:nlyr(mod), grab_values) + 0.000001)
            if (any(is.na(out))) browser()
            write.csv(out, tfile)
            return(out)
        }
        dat = sapply(files, openDat)
        
        write.csv(dat, tfile)
        return(dat)
    }
    
    bias_correct <- function(a, b, obs) {
        mnR = mean(a, na.rm = TRUE)/mean(b, na.rm = TRUE)
        
        a = a - mean(a)
        a = a * sd(obs)/sd(b)
        a = a + mean(obs)*mnR
    }
    factual = unlist(openDats('Greece_factual2', factual_dir))
    counter = unlist(openDats('Greece_counter2', counter_dir))
    event = tail(areas_pc, 1)
    
    #counter = bias_correct(counter, factual, areas_pc)
    #factual = bias_correct(factual, factual, areas_pc)
    
    counter = exp(counter)
    factual = exp(factual)
    event = exp(event)
    if (plot) {
        bins = seq(min(counter, factual), max(counter, factual), length.out = sqrt(length(counter)))
        chist = hist(counter, breaks = bins, plot = FALSE)$density
        fhist = hist(factual, breaks = bins, plot = FALSE)$density  
        bins = bins[-1] - diff(bins)/2
        scale = max(chist, fhist)
        chist = chist / scale; fhist = fhist /scale
        plot(range(bins), c(0, 1), type = 'n', xlab = '', ylab = '', axes = FALSE, log = 'x')
        axis(1)
        for (i in 1:5) {
            polygon(c(bins[1], bins, tail(bins, 1)), c(0, chist, 0), 
                    col = '#0000FF11', border = '#0000FF')
            polygon(c(bins[1], bins, tail(bins, 1)), c(0, fhist, 0), 
                    col = '#FF000011', border = '#FF0000')
        }
        
    }
    browser()
    return(mean(factual > event)/mean(counter > event))
}

for_pc(0.0, plot = TRUE, 'control')
browser()
pcs = seq(0.5, 1, 0.05)
rrs = sapply(pcs, for_pc)

