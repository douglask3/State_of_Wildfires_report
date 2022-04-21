graphics.off()
library(raster)
source("libs/plotStandardMap.r")
source("../gitProjectExtras/gitBasedProjects/R/sourceAllLibs.r")
sourceAllLibs("../rasterextrafuns/rasterPlotFunctions/R/")
sourceAllLibs("../rasterextrafuns/rasterExtras/R/")
sourceAllLibs("libs/")


x = seq(0.1, 0.8, by = 0.0005)

xl = logit(x)
x = (x-0.1)*30 

muH = logit(0.3)
sdH = 0.6

muF = logit(0.32)
sdF = 0.7

distr <- function(mu, sd) {
    y = dnorm(xl, mu, sd)
    y = y + dnorm(xl, -2.19, 0.3) * 0.3
    y = y/max(y)        
    cy = pnorm(xl, mu, sd)
    selectXatT <- function(thresh) 
        x[which(cy > thresh)[1]]
   
    list(y, selectXatT(0.1), selectXatT(0.5), selectXatT(0.9), selectXatT(0.99))
}

yH = distr(muH, sdH)
yF = distr(muF, sdF)

png("figs/BAchange_extreme_example.png", width = 7.5*0.9, height = 5*0.9, units = 'in', res = 300)
plot(range(x), c(0, 1.05), type = 'n', axes = FALSE, xlab = '', ylab = '', xaxs = 'i', yaxs = 'i')
lines(rep(min(x), 2),c(-100, 100))
#axis(1)
mtext(side = 1, line = 0.5, "Burnt Area (%)", font = 2)
mtext(side = 2, line = 0.5, "Probability of occurance", font = 2)

addBox <- function(H, col, txt) {
    polygon(H[c(1, 2, 2, 1)], c(0, 0, 1.1, 1.1), border = NA, col = make.transparent(col, 0.0))
    polygon(H[c(1, 2, 2, 1)], c(0, 0, 1.1, 1.1), border = NA, col = make.transparent("white", 0.2))
    #lines(H[c(2,2)], c(0, 1), lty = 3, lwd = 2)
    text(x = mean(H), y = 0.9, txt)
}

cols = rev(c('#d7191c','#fdae61','#ffffbf','#abd9e9','#2c7bb6'))
test = x > yH[[5]]

addBox(c(min(x), yH[[2]]), cols[1], "No\nFire")
addBox(c(yH[[5]], max(x)), cols[5], "Wildfire")

pFUN <- function(y, ...)
    polygon(c(x[test][1], x[test], tail(x[test], 1)), c(0, y[[1]][test], 0), col = tail(cols, 1), border = NA, ...)
pFUN(yH)
pFUN(yF, density = 10, lwd = 4, lend = 1, angle = 60)

for (i in 2:4) 
    addBox(c(yH[[i]], yH[[i+1]]), cols[i], c("", "", "High Fire\nActivity")[i-1])
lines(x, yH[[1]], lwd = 4)
lines(x, yF[[1]], lwd = 4, lty = 2)

text(10, 0.53, srt = -50, 'Historic (2010s)')
text(12, 0.64, srt = -45, 'Future')
dev.off()
browser()
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
