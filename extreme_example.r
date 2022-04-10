graphics.off()
library(raster)
source("libs/plotStandardMap.r")
source("../gitProjectExtras/gitBasedProjects/R/sourceAllLibs.r")
sourceAllLibs("../rasterextrafuns/rasterPlotFunctions/R/")
sourceAllLibs("../rasterextrafuns/rasterExtras/R/")
sourceAllLibs("libs/")


x = seq(0, 1, by = 0.001)

xl = logit(x)

muH = logit(0.3)
sdH = 0.6

muF = logit(0.32)
sdF = 0.7

distr <- function(mu, sd) {
    y = dnorm(xl, mu, sd)
    y = y/max(y)        
    cy = pnorm(xl, mu, sd)
    selectXatT <- function(thresh) 
        x[which(cy > thresh)[1]]
   
    list(y, selectXatT(0.05), selectXatT(0.95), selectXatT(0.99), selectXatT(0.995))
}

yH = distr(muH, sdH)
yF = distr(muF, sdF)

png("figs/BAchange_extreme_example.png", width = 7, height = 7, units = 'in', res = 300)
plot(c(0, 1), c(0, 1.3), type = 'n', axes = FALSE, xlab = '', ylab = '')

addPoly <- function(xp, y, col, border= col) 
    polygon(xp, y[[1]], col = make.transparent(col, 0.9), border = border)

addPolyLines <- function(y, col, pos, txt) {
    addPoly(x, y, col)
    lapply(y[2:3], function(x) lines(rep(x, 2), c(0, pos), lty = 2, col = col))
    lapply(y[2:3], function(x) lines(rep(x, 2), pos + c(-0.05, 0.05), col = col))
    lines(unlist(y[2:3]), rep(pos, 2), col = col)
    text(x = mean(unlist(y[2:3])), y = pos+0.05, txt, col = col)
}
#for (i in 1:3) {
addPolyLines(yH, "blue", 1, txt = 'Historic range')
addPolyLines(yF, "red", 1.1, txt = 'Future range')
#}

id = which(x == yH[[4]])
lines(rep(yH[[4]], 2), c(0, 0.3), lty = 3, col = 'blue', lwd = 2)
#lines(c(yH[[4]], 1), rep(yH[[1]][id], 2), lty = 3, col = 'black', lwd = 2)
#lines(c(yH[[4]], 1), rep(yF[[1]][id], 2), lty = 3, col = 'black', lwd = 2)
xp = x
xp[x < yH[[4]]] = yH[[4]]

for (i in 1:5) {
    addPoly(xp, yH, "blue", "transparent")
    addPoly(xp, yF, "red", "transparent")
}

arrows(0.8, 0.6, yH[[4]], 0.3, col = 'blue')
text(x = 0.82, y = 0.6, adj = 0, 'Historic extreme\n(1:100 event)', xpd = NA, col = 'blue')

arrows(0.8, 0.4, yH[[5]], 0.02, col = 'blue')
text(x = 0.82, y = 0.4, adj = 0, 'Prob. Historic\nextreme', col = 'blue')
arrows(0.8, 0.2, yF[[5]], 0.02, col = 'red')
text(x = 0.82, y = 0.2, adj = 0, 'Prob. Future\nextreme', col = 'red')

arrows(0, -0.05, 1, -0.05, xpd = NA, lwd = 2)
text(x = 0.5, y = -0.1, 'Burnt Area', xpd = NA, cex = 1.2)


arrows(-0.05, -0.0, -0.05, 1, xpd = NA, lwd = 2)
text(x = -0.1, y = 0.5, 'Probability', srt = 90, xpd = NA, cex = 1.2)

dev.off()
