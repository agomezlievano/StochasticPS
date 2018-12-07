# EComm_0001_FigsToyWorlds.R

# #########################################################

rm(list = ls(all = TRUE))  # resets R to fresh
gc()

# ========================================================
# Libraries
# --------------------------------------------------------
library(sfsmisc)
library(scales)
library(pheatmap)
library(ggplot2)

library(ggpubr) # to make grids of ggplots, see: http://www.sthda.com/english/wiki/ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page
###############################
## try http:// if https:// URLs are not supported
# source("https://bioconductor.org/biocLite.R")
# biocLite("EBImage")
library("EBImage")

# - ########################################################################
put.fig.letter <- function(label, location="topleft", x=NULL, y=NULL, 
                           offset=c(0, 0), ...) {
  if(length(label) > 1) {
    warning("length(label) > 1, using label[1]")
  }
  if(is.null(x) | is.null(y)) {
    coords <- switch(location,
                     topleft = c(0.015,0.98),
                     topcenter = c(0.5525,0.98),
                     topright = c(0.985, 0.98),
                     bottomleft = c(0.015, 0.02), 
                     bottomcenter = c(0.5525, 0.02), 
                     bottomright = c(0.985, 0.02),
                     c(0.015, 0.98) )
  } else {
    coords <- c(x,y)
  }
  this.x <- grconvertX(coords[1] + offset[1], from="nfc", to="user")
  this.y <- grconvertY(coords[2] + offset[2], from="nfc", to="user")
  text(labels=label[1], x=this.x, y=this.y, xpd=T, ...)
}
# - ########################################################################


# ========================================================
# Some declarations
# --------------------------------------------------------
workingpath <- "~/../Dropbox/Harvard/LittleProjects/StochasticPS"

#  Declare other paths
folder.outputdata <- paste0(workingpath,"/outputdata")
folder.figures <- paste0(workingpath,"/figures")
folder.rprograms <- paste0(workingpath,"/programsR")

setwd(workingpath)


# --------------------------------------------------------
# useful functions
colfunc <- colorRampPalette(c("deepskyblue4", "deepskyblue", "cyan"))
col7 <- colorRampPalette(c("#f7f7f7", "#74add1", "#2166ac", "#053061")) 

# ========================================================
# Name of figures for UNIFORM
# --------------------------------------------------------

namefile1 <- "/EComm_0008_uniform_Mcp_matrix_5communities.png"
namefile2 <- "/EComm_0008_uniform_C_matrix_5communities.png"
namefile3 <- "/EComm_0008_uniform_Frequencies_5communities.png"
namefile4 <- "/EComm_0008_uniform_toyMcp_5communities.png"


# -------------------------------------------------------
# First, join the matrices
filenames <- c(namefile1, 
               namefile2)


foo <- list()
for(j in c(1:length(filenames))){
  foo[[j]] <- readImage(paste(folder.figures, filenames[j],sep=""))
} 

my.locations <- c("topleft", 
                  "topleft")
my.labels <- c("A", 
               "B")

namefile <- "/EComm_0001_uniform_matrices.png"
png(paste0(folder.figures, namefile), 
    height=3, width=10, units="in", res=300)

mult.fig(1, marP = c(0,1,1,0))
layout(matrix(c(1, 2), 
              nrow=1, ncol=2, byrow = TRUE), 
       heights=c(5), widths = c(6,4), respect=TRUE)


for(j in c(1:length(filenames))){
  #cat(j)
  display(foo[[j]], method="raster")
  #put.fig.letter(label=my.labels[j], location=my.locations[j], font=2, offset=c(0.010,-0.04))  
  put.fig.letter(label=my.labels[j], location=my.locations[j], font=2)  
}

dev.off()


# -------------------------------------------------------
# Second, join the previous with the density and the scatters
filenames <- c(namefile,
               namefile3, 
               namefile4)


foo <- list()
for(j in c(1:length(filenames))){
  foo[[j]] <- readImage(paste(folder.figures, filenames[j],sep=""))
} 

my.locations <- c("topleft", 
                  "topleft", 
                  "topleft")
my.labels <- c("AB",
               "C", 
               "D")

namefile_final <- "/EComm_0001_uniform_alltogether.png"
png(paste0(folder.figures, namefile_final), 
    height=5, width=5, units="in", res=300)

mult.fig(1, marP = c(0,1,1,0))
layout(matrix(c(1, 1, 2, 2, 3, 3), 
              nrow=3, ncol=2, byrow = TRUE), 
       heights=c(1, 2, 3), widths = c(5, 5), respect=TRUE)


for(j in c(1:length(filenames))){
  #cat(j)
  display(foo[[j]], method="raster")
  #put.fig.letter(label=my.labels[j], location=my.locations[j], font=2, offset=c(0.010,-0.04))  
  if(j>2){
    put.fig.letter(label=my.labels[j], location=my.locations[j], font=2, offset=c(-0.010,-0.04))  
  }
}

dev.off()




