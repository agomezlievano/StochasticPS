# EComm_0002_growthregressions.R

# ###########################################################
# Created on 2018-02-07
# Program for eigenspace regressions
# Author: Andres Gomez-Lievano
# ###########################################################


rm(list = ls(all = TRUE))  # resets R to fresh
gc()


library(ggplot2)
library(scales)
library(sfsmisc)
library(scales)
library(plyr)
library(dplyr)
library(tidyr)
library(pheatmap)
library(lfe)
library(data.table)
library(pryr)
library(stringr)
library(MASS)
#library(ineq)

#==========================================================================

# Andres path
workingfolder <- "~/../Dropbox/Harvard/LittleProjects/StochasticPS/"


#  Declare other paths
folder.data <- paste0(workingfolder,"data/")
folder.datain <- paste0(workingfolder,"inputdata/")
folder.dataout <- paste0(workingfolder,"outputdata/")
folder.figs <- paste0(workingfolder,"figures/")
folder.tables <- paste0(workingfolder,"tables/")
folder.rcodes <- paste0(workingfolder,"programsR/")

setwd(workingfolder)
#source("NHB_9999_functions.R")

#==========================================================================



#==========================================================================
# LOADING DATA

# Back From STATA
filefromstata <- "EComm_0002_2_fromStata2R.csv"
newdf <- as.data.table(rio::import(paste0(folder.dataout, filefromstata)))

ents <- newdf %>%
  group_by(year) %>%
  summarise(EntropyCommunities_ccp = first(EntropyCommunities_ccp), 
            EntropyCommunities_cp = first(EntropyCommunities_cp)) %>%
  ungroup() %>%
  mutate(EntropyCommunities_ccp = EntropyCommunities_ccp/mean(EntropyCommunities_ccp, na.rm=TRUE),
         EntropyCommunities_cp = EntropyCommunities_cp/mean(EntropyCommunities_cp, na.rm=TRUE)) %>%
  as.data.table()
  

#==========================================================================
# plotting the decline in the number of communities

# -----------------------------
namefile <- "EComm_0003_Worldcoherence.png"
png(paste0(folder.figs,namefile), 
    height=6, width=8, units="in", res=600)
# -----------------------------

mult.fig(1)
yrange <- range(c(ents$EntropyCommunities_cp, ents$EntropyCommunities_ccp), na.rm=TRUE)
plot(1, 1,
     type = "n",
     xlim = range(ents$year),
     ylim = yrange,
     xlab = "Year",
     ylab = "Decoherence")
grid()

abline(lm(ents$EntropyCommunities_ccp ~ ents$year), 
       col=alpha("darkred", 0.3), lwd=2)
points(ents$year, ents$EntropyCommunities_ccp, 
       pch=21, bg="darkred", col = "gray", cex=0.5)

abline(lm(ents$EntropyCommunities_cp ~ ents$year), 
       col=alpha("dodgerblue", 0.5), lwd=4)
points(ents$year, ents$EntropyCommunities_cp,
       pch=21, bg="dodgerblue", col = "gray", cex=1.5)

legend("topleft", bg="gray90",
       #bty = "n", 
       inset = c(0.03, 0.03),
       legend = c("From Exporter vs. Product", "From Exporter vs. Importer-Product"),
       pch = c(21, 21),
       pt.bg = c("dodgerblue", "darkred"),
       pt.cex = c(1.5, 0.5),
       col = c("gray", "gray"))

# -----------------------------
dev.off()
# -----------------------------

# =========================================================================
# Growth of distances?
mult.fig(4, mfrow = c(2,2))
plot(leftDistance2Origin_cp ~ year, data=newdf, type="n")
lines(newdf[exporter=="THA"]$year, newdf[exporter=="THA"]$leftDistance2Origin_cp, col = "red")
lines(newdf[exporter=="USA"]$year, newdf[exporter=="USA"]$leftDistance2Origin_cp, col = "green")
lines(newdf[exporter=="CHN"]$year, newdf[exporter=="CHN"]$leftDistance2Origin_cp, col = "blue")
lines(newdf[exporter=="IND"]$year, newdf[exporter=="IND"]$leftDistance2Origin_cp, col = "black")

plot(rightDistance2Origin_cp ~ year, data=newdf, type="n")
lines(newdf[exporter=="THA"]$year, newdf[exporter=="THA"]$rightDistance2Origin_cp, col = "red")
lines(newdf[exporter=="USA"]$year, newdf[exporter=="USA"]$rightDistance2Origin_cp, col = "green")
lines(newdf[exporter=="CHN"]$year, newdf[exporter=="CHN"]$rightDistance2Origin_cp, col = "blue")
lines(newdf[exporter=="IND"]$year, newdf[exporter=="IND"]$rightDistance2Origin_cp, col = "black")

plot(leftSimDirection2MostComplex_cp ~ year, data=newdf, type="n")
lines(newdf[exporter=="THA"]$year, newdf[exporter=="THA"]$leftSimDirection2MostComplex_cp, col = "red")
lines(newdf[exporter=="USA"]$year, newdf[exporter=="USA"]$leftSimDirection2MostComplex_cp, col = "green")
lines(newdf[exporter=="CHN"]$year, newdf[exporter=="CHN"]$leftSimDirection2MostComplex_cp, col = "blue")
lines(newdf[exporter=="IND"]$year, newdf[exporter=="IND"]$leftSimDirection2MostComplex_cp, col = "black")

plot(leftDistance2MostComplex_cp ~ year, data=newdf, type="n")
lines(newdf[exporter=="THA"]$year, newdf[exporter=="THA"]$leftDistance2MostComplex_cp, col = "red")
lines(newdf[exporter=="USA"]$year, newdf[exporter=="USA"]$leftDistance2MostComplex_cp, col = "green")
lines(newdf[exporter=="CHN"]$year, newdf[exporter=="CHN"]$leftDistance2MostComplex_cp, col = "blue")
lines(newdf[exporter=="IND"]$year, newdf[exporter=="IND"]$leftDistance2MostComplex_cp, col = "black")

# =========================================================================
# step AIC on levels
summary(newdf)
colnames(newdf)

# -------------------------
fullglm <- glm(loggdppc ~ ., data=subset(newdf, select = -c(exporter, year, population, 
                                                            #Distance2Origin_cp, SimDirection2MostComplex_ccp,
                                                            latitude, longitude, region, 
                                                            TotalExports, TotalImports, rX10yr, rY10yr,
                                                            GDPpcPPPconstant2011intUSD, GDPpcconst2010USD, incomeLevel, 
                                                            EntropyCommunities_ccp, EntropyCommunities_cp)))
step <- stepAIC(fullglm, trace=TRUE)
step$anova
print(step$formula)
summary(lm(as.formula(step$formula), data=newdf))
summary(lm(loggdppc ~ population + eci + 
             rightDistance2Origin_ccp + 
             rightSimDirection2MostComplex_cp + 
             leftDistance2MostComplex_cp + logpop, data=newdf))
# -------------------------
fullglm.left <- glm(loggdppc ~ ., data=subset(newdf, select = c("loggdppc", "logpop", "eci",
                                                           str_subset(colnames(newdf), fixed("left")))))
fullglm.right <- glm(loggdppc ~ ., data=subset(newdf, select = c("loggdppc", "logpop", "eci",
                                                                str_subset(colnames(newdf), fixed("right")))))

step.left <- stepAIC(fullglm.left, trace=TRUE)
step.right <- stepAIC(fullglm.right, trace=TRUE)
step.left$anova
step.right$anova
print(step.left$formula)
print(step.right$formula)
summary(lm(as.formula(step.left$formula), data=newdf))
summary(lm(as.formula(step.right$formula), data=newdf))

summary(lm(loggdppc ~ population + eci + 
             rightDistance2Origin_ccp + 
             rightSimDirection2MostComplex_cp + 
             leftDistance2MostComplex_cp + logpop, data=newdf))

step2 <- stepAIC(fullglm, ~.^2, trace=FALSE) # with pairwise interactions
step2$anova
summary(lm(as.formula(step2$formula), data=newdf))


# step AIC on growth
fullglm <- glm(rY10yr ~ ., data=subset(newdf, select = -c(exporter, year, population,
                                                          #Distance2Origin_cp, SimDirection2MostComplex_ccp, Distance2MostComplex_ccp,
                                                          latitude, longitude, region, 
                                                          TotalExports, TotalImports, rX10yr,
                                                          GDPpcPPPconstant2011intUSD, GDPpcconst2010USD, incomeLevel, 
                                                          EntropyCommunities_ccp, EntropyCommunities_cp)))
step <- stepAIC(fullglm, trace=TRUE)
step$anova
print(step$formula)
summary(lm(as.formula(step$formula), data=newdf))

step2 <- stepAIC(fullglm, ~.^2, trace=FALSE)
step2$anova
summary(lm(as.formula(step2$formula), data=newdf))
print(step2$formula)
summary(lm(rY10yr ~ Distance2Origin_ccp + SimDirection2MostComplex_cp + 
             logpop + loggdppc + Distance2Origin_ccp:loggdppc + 
             logpop:loggdppc + SimDirection2MostComplex_cp:logpop + SimDirection2MostComplex_cp:loggdppc + 
             Distance2Origin_ccp:logpop + eci:SimDirection2MostComplex_cp + 
             eci:Distance2MostComplex_cp + Distance2Origin_ccp:Distance2MostComplex_cp, data=newdf))





# ---------------
cols2use <- c("logpop", "loggdppc", "eci", 
              "leftDistance2Origin_cp", "leftDistance2Origin_ccp",
              "leftSimDirection2MostComplex_cp", 
              "leftDistance2MostComplex_cp", "leftDistance2MostComplex_ccp", 
              "rightDistance2Origin_cp", "rightDistance2Origin_ccp", 
              "rightSimDirection2MostComplex_cc", 
              "rightDistance2MostComplex_cp")

# Libraries
library(ellipse)
library(RColorBrewer)

# Use of the mtcars data proposed by R
data=cor(newdf[,.SD, .SDcols = cols2use])

library(corrplot)
corrplot(data, method = "color", outline = T, 
         addgrid.col = "darkgray", order="hclust", 
         addrect = 4, rect.col = "black", rect.lwd = 5,
         cl.pos = "b", tl.col = "indianred4", tl.cex = 0.8, tl.pos = "l",
         cl.cex = 1.5, addCoef.col = "white", number.digits = 2, 
         number.cex = 0.75, col = colorRampPalette(c("darkred","white","midnightblue"))(100))


