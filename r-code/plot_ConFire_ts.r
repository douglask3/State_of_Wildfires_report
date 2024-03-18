
make_plot <- function() {
openDat <- function(exp, nav = 0) {
    dat = read.csv(paste0(fname, exp, fileEx))
    
    if (nav>0)  
        dat[,2:3] = apply(dat[,2:3], 2, function(x) filter(x, rep(1/nav, nav), sides = 1))
    #browser()
    return(dat)
}
dats = lapply(experiments, openDat, nav)


x_range = range(sapply(dats, function(i) i[,1]), na.rm = T)
y_range = range(sapply(dats, function(i) i[,-1]), na.rm = T)

plot(x_range, y_range, type = 'n', xlab = '', ylab = '')
grid()
mtext(side = 2, 'burnt area (%/month)', line = 2)
#polygon(c(-9E9, 9E9)[c(1, 2, 2, 1)], c(-9E9, 9E9)[c(1, 1, 2, 2)], col = 'black')

plot_exp <- function(dat, col) {
    
    polygon(c(dat[,1], rev(dat[,1])), c(dat[,2], rev(dat[,3])), 
            border = NA, col = paste0(col, '44'))
}


for (i in 1:5) mapply(plot_exp, dats, cols)
}

png("r-code/ConFire_TS.png", height = 6, width = 6, res = 300, units = 'in')
fname = "outputs/ConFire_Canada-biasC3/figs/crop_lightn_soilM_tree_biascorrected_csoil_pas_vpd_cveg_precip_tas_rhumid_totalVeg_biascorrected-frac_points_0.02-"
fileEx = "-control_TS/time_seriesControl.csv"

experiments = c("factual", "ss126_GFDL", "ss126_IPSL", "ss126_MPI",# "ss126_MRI",
                "ss585_GFDL", "ss585_IPSL", "ss585_MPI")#, "ss585_MRI")
cols = c('#000000', '#0000FF', '#0000FF', '#0000FF', '#FF0000', '#FF0000', '#FF0000')#, '#FF0000')
nav = 6
make_plot()

legend('topleft', legend = c('historic', 'ssp126', 'ssp585'), col = c('black', 'blue', 'red'), pch = 19, pt.cex = 2)
dev.off()

experiments = c("Fuel", "Moisture", "Ignition", "Suppression")
cols = c('#00FFFF', '#FFFF00', '#9999FF', '#9999FF', '#FF9999')
nav = 0
#make_plot()

