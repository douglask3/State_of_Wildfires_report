library(terra)

dir = '../data/data/driving_data/'

regions = c('Greece', 'Greece')

isimip3b = c('historical', 'ssp126', 'ssp370', 'ssp585')
models = c('GFDL-ESM4', 'IPSL-CM6A-LR', 'MPI-ESM1-2-HR', 'MRI-ESM2-0', 'UKESM1-0-LL')
isimp3a = c('obsclim', 'counterclim')
isimip3a_cols = c('black', '#648FFF')
imimip3a_lwd = 2
imimip3a_lty = 1

isimip3b_cols = c('#DC267F', '#785EF0', '#FFB000', '#FE6100')
imimip3b_lwd = 1
imimip3b_lty = 2

isimip3a_period = c('period_2000_2019')
isimip3b_period = c('period_1994_2014', rep('period_2015_2099', length(isimip3b)-1))

runs = c(paste0('isimp3a/', isimp3a, '/GSWP3-W5E5/', isimip3a_period, '/'),
         paste0('isimp3b/', mapply(function(a, b) paste0(a, '/',  models, '/', b, '/'),
                                        isimip3b, isimip3b_period)))
spread_info4plot <- function(a, b) 
    c(a, rep(b, each = length(models)))


names(runs) = spread_info4plot(isimp3a, isimip3b)
cols = spread_info4plot(isimip3a_cols, isimip3b_cols)              
lwds = spread_info4plot(rep(imimip3a_lwd, 2), rep(imimip3b_lwd, each = length(isimip3b_period)))   
ltys = spread_info4plot(rep(imimip3a_lty, 2), rep(imimip3b_lty, each = length(isimip3b_period)))      

colL = c(isimip3a_cols, isimip3b_cols)
lwdL = c(imimip3a_lty, rep(imimip3b_lwd, length(isimip3b_period)))
ltyL = c(imimip3a_lty, rep(imimip3b_lty, length(isimip3b_period)))


#variables = c('% Tree cover - jules' = 'tree_cover_jules-es.nc', '% None-tree cover - jules' = 'nonetree_cover_jules-es.nc',
#              '% Tree cover - debiased' = 'debiased_tree_cover_jules-es.nc', '% None-tree cover - debiased' = 'debiased_nonetree_cover_jules-es.nc')
variables = c('Consecutive Dry Days' = "consec_dry_mean.nc", "Max. Temp" = "tas_max.nc")
variable = variables[1]
region = regions[1]
smooth = c(max, max)
forVariable <- function(variable, title) { 
    openDat <- function(run) {
        dat = rast(paste0(dir, region, '/',run, variable))
        tm = time(dat)
        tm = as.numeric(substr(tm, 1, 4)) + as.numeric(substr(tm, 6, 7))/12
        dat = dat * cellSize(dat)/expanse(dat)
    
        layer_sums <- rep(0, dim(dat)[3])
        for (i in 1:dim(dat)[3]) 
            layer_sums[i] <- sum(dat[[i]][], na.rm = TRUE)
        if (!is.null(smooth)) layer_sums = rollmean(layer_sums, k = 12, align = "right", fill = NA)
        return(cbind(tm, layer_sums))
    }

    dats = lapply(runs, openDat)

    xrange = range(sapply(dats, function(i) i[,1]), na.rm = TRUE)
    yrange = range(sapply(dats, function(i) i[,2]), na.rm = TRUE)
    plot(xrange, yrange, type = 'n', xlab = '', ylab = '')

    plotLine <- function(dat, ...) lines(dat[,1], dat[,2], ...)

    mapply(plotLine, dats, col = cols, lwd = lwds, lty = ltys)


    title(title)
}

nrows = ceiling(sqrt(length(variables)))
ncols = ceiling(length(variables)/nrows)
lmat = rbind(matrix(1 + 1:(nrows * ncols), nrow = nrows), 1)
png(paste0("../outputs/figs/input_vars_for_", region, '.png'), height = 3*(nrows + 0.1), width = 5 * ncols, units = 'in', res = 300)
layout(lmat, heights = c(rep(1, nrows), 0.2))
par(mar = c(2, 0,0, 0))
plot.new()
legend('center', legend = unique(names(runs)), col = colL, lty = ltyL, lwd = lwdL, ncol = 3)
par(mar = c(2, 2, 1, 1))
mapply(forVariable, variables, names(variables))

dev.off()
