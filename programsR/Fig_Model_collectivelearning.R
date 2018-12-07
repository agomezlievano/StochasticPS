

# Andres path
workingpath <- "~/../Dropbox/Harvard/LittleProjects/StochasticPS/"

# Oscar path
#workingpath <- "C:/Users/oskky/Dropbox/NHBcommentary/"


#  Declare other paths
inpath <- paste0(workingpath,"inputdata/")
outpath <- paste0(workingpath,"outputdata/")
figpath <- paste0(workingpath,"figures/")
tablespath <- paste0(workingpath,"tables/")
rcodepath <- paste0(workingpath,"programs/")


library(dplyr)
library(tidyr)
library(data.table)

library(ggplot2)
library(sfsmisc)
library(sfsmisc)
library(scales)

Mps <- seq(1,20,1)
Mp <- 8
deltaMp <- 1


rcs <- seq(0.05,0.95,0.05)
rc <- 0.8
drc <- deltaMp/Mp

sis <- seq(0.05,0.95,0.05)
si <- 0.2
dsi <- deltaMp/Mp

picp <- function(Mp, rc, si){
  exp(-Mp*(1-si)*(1-rc))
}

# pltpicp <- function(x, y1, y2){
#   
# }

textxlabelsi <- expression(s[i])
textxlabelrc <- expression(r[c])
textxlabelMp <- expression(M[p])
textylabel <- ""
textylabelall <- expression(Pr(X[paste(i, ",", c, ",", p)]==1))
numxlabelcex <- 0.8
textxlabelcex <- 2
seglength <- 0.5
posincrease <- "topleft"
posdecrease <- "topright"
colbottom <- alpha("#6baed6",0.8)
coltop <- alpha("#084594",0.8)
coltoparrow <- alpha("gray10",0.95)
lwdbottom <- 1


namefile <- "Fig_Model_param_comparison.png"
png(paste0(figpath,namefile), 
    height=4, width=7, units="in", res=300)


mult.fig(6, mfrow = c(2,3),
         marP = c(0, 0, 0, 0),
         mgp=c(2.2, .1, 0), oma=c(0, 1, 0, 0))

# varying Mp, fixing rc and si
plot(sis, picp(Mp, rc, sis),
     ylim = c(0,1),
     #main = "Tech. improvement", 
     xlab = textxlabelsi,
     ylab = textylabel,
     cex.lab = textxlabelcex,
     cex.axis = numxlabelcex,
     tck = -0.01,
     type = "l",
     lwd = lwdbottom,
     col = colbottom)
lines(sis, picp(Mp - deltaMp, rc, sis),
      lwd = 2,
      col = coltop)
legend(posincrease,
       bty = "n",
       bg="transparent",
       lty = 1,
       lwd = 2,
       col = coltop,
       legend = "Effect of a\ntech. improvement",
       seg.len = seglength)
arrows(0.5, 0.9*picp(Mp, rc, 0.5),
       0.5, 1.1*picp(Mp - deltaMp, rc, 0.5),
       col = coltoparrow,
       length = 0.1,
       lwd = 1)

# varying Mp, fixing rc and si
plot(rcs, picp(Mp, rcs, si),
     ylim = c(0,1),
     #main = "Individual learning", 
     xlab = textxlabelrc,
     ylab = textylabel,
     cex.lab = textxlabelcex,
     cex.axis = numxlabelcex,
     tck = -0.01,
     type = "l",
     lwd = lwdbottom,
     col = colbottom)
lines(rcs, picp(Mp, rcs, si + dsi),
      lwd = 2,
      col = coltop)
legend(posincrease,
       bty = "n",
       bg="transparent",
       lty = 1,
       lwd = 2,
       col = coltop,
       legend = "Effect of individual\nlearning",
       seg.len = seglength)
arrows(0.7, 0.8*picp(Mp, 0.7, si),
       0.7, 1.2*picp(Mp, 0.7, si + dsi),
       col = coltoparrow,
       length = 0.1,
       lwd = 1)


# varying Mp, fixing rc and si
plot(Mps, picp(Mps, rc, si),
     ylim = c(0,1),
     #main = "Collective learning", 
     xlab = textxlabelMp,
     ylab = textylabel,
     cex.lab = textxlabelcex,
     cex.axis = numxlabelcex,
     tck = -0.01,
     type = "l",
     lwd = lwdbottom,
     col = colbottom)
lines(Mps, picp(Mps, rc + drc, si),
      lwd = 2,
      col = coltop)
legend(posdecrease,
       bty = "n",
       bg="transparent",
       lty = 1,
       lwd = 2,
       col = coltop,
       legend = "Effect of collective\nlearning",
       seg.len = seglength)
arrows(Mp, 0.9*picp(Mp, rc, si),
       Mp, 1.1*picp(Mp, rc + drc, si),
       col = coltoparrow,
       length = 0.1,
       lwd = 1)


plot(rcs, picp(Mp, rcs, si),
     #main = "Tech. improvement", 
     xlab = textxlabelrc,
     ylab = textylabel,
     cex.lab = textxlabelcex,
     cex.axis = numxlabelcex,
     tck = -0.01,
     type = "l",
     lwd = lwdbottom,
     col = colbottom)
lines(rcs, picp(Mp - deltaMp, rcs, si),
      lwd = 2,
      col = coltop)
legend(posincrease,
       bty = "n",
       bg="transparent",
       lty = 1,
       lwd = 2,
       col = coltop,
       legend = "Effect of a\ntech. improvement",
       seg.len = seglength)
arrows(0.7, 0.8*picp(Mp, 0.7, si),
       0.7, 1.2*picp(Mp - deltaMp, 0.7, si),
       col = coltoparrow,
       length = 0.1,
       lwd = 1)


plot(Mps, picp(Mps, rc, si),
     ylim = c(0,1),
     #main = "Individual learning", 
     xlab = textxlabelMp,
     ylab = textylabel,
     cex.lab = textxlabelcex,
     cex.axis = numxlabelcex,
     tck = -0.01,
     type = "l",
     lwd = lwdbottom,
     col = colbottom)
lines(Mps, picp(Mps, rc, si + dsi),
      lwd = 2,
      col = coltop)
legend(posdecrease,
       bty = "n",
       bg="transparent",
       lty = 1,
       lwd = 2,
       col = coltop,
       legend = "Effect of individual\nlearning",
       seg.len = seglength)
arrows(Mp, 0.8*picp(Mp, rc, si),
       Mp, 1.2*picp(Mp, rc, si + dsi),
       col = coltoparrow,
       length = 0.1,
       lwd = 1)


plot(sis, picp(Mp, rc, sis),
     ylim = c(0,1),
     #main = "Collective learning", 
     xlab = textxlabelsi,
     ylab = textylabel,
     cex.lab = textxlabelcex,
     cex.axis = numxlabelcex,
     tck = -0.01,
     type = "l",
     lwd = lwdbottom,
     col = colbottom)
lines(sis, picp(Mp, rc+drc, sis),
      lwd = 2,
      col = coltop)
legend(posincrease,
       bty = "n",
       bg="transparent",
       lty = 1,
       lwd = 2,
       col = coltop,
       legend = "Effect of collective\nlearning",
       seg.len = seglength)
arrows(0.5, 0.95*picp(Mp, rc, 0.5),
       0.5, 1.05*picp(Mp, rc+drc, 0.5),
       col = coltoparrow,
       length = 0.1,
       lwd = 1)


mtext(textylabelall, side = 2, -3, outer=TRUE, las=0, cex=textxlabelcex)

dev.off()


