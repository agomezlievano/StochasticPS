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

datafilename <- "EComm_0003_constructingdataallyears.csv"
gdpfilename <- "WorldBank_GDPperCapita_1962_to_2015.csv"
cpyfilename <- "S2_final_cpy_all.dta"

# Trade Data with distances
tradedata <- as.data.table(rio::import(file=paste0(folder.dataout,datafilename)))

# GDP Data
gdpdata <- as.data.table(rio::import(file=paste0(folder.data,gdpfilename)))

# GDP Data
cpydata <- as.data.table(rio::import(file=paste0(folder.datain,cpyfilename)))
# cpydata <- fread(paste0(folder.datain,cpyfilename), 
#                  header = TRUE, data.table = TRUE)


#==========================================================================
# aggregating total exports
totalexportsbycountry <- cpydata %>%
  group_by(year, exporter) %>%
  summarise(TotalExports = sum(export_value, na.rm = TRUE)/1000000,
            TotalImports = sum(import_value, na.rm = TRUE)/1000000) %>%
  ungroup() %>% 
  mutate(year = as.integer(year)) %>%
  as.data.table()
rm(cpydata)
gc()
mem_used()


#==========================================================================
# Merging
cols2keep <- c("gdp_per_capita_PPP_constant_2011_international_dollar",
               "gdp_per_capita_constant2010USD",
               "incomeLevel",
               "iso3c",
               "latitude",
               "longitude",
               "region",
               "year")
gdpdata <- gdpdata[,.SD, .SDcols = c(cols2keep)]
newcols <- c("GDPpcPPPconstant2011intUSD",
             "GDPpcconst2010USD",
             "incomeLevel",
             "iso3c",
             "latitude",
             "longitude",
             "region",
             "year")
colnames(gdpdata) <- newcols
newdf <- as.data.table(merge(tradedata, gdpdata, 
                             by.x=c("exporter", "year"), by.y=c("iso3c", "year"), 
                             all.x = TRUE))

newdf <- as.data.table(merge(newdf, totalexportsbycountry, 
                             all.x = TRUE))


#==========================================================================
# Go to STATA
file2stata <- "EComm_0002_fromR2Stata.csv"
rio::export(newdf, paste0(folder.dataout, file2stata))


#==========================================================================
# Back From STATA
filefromstata <- "EComm_0002_2_fromStata2R.csv"
newdf <- as.data.table(rio::import(paste0(folder.dataout, filefromstata)))

ents <- newdf %>%
  group_by(year) %>%
  summarise(EntropyCommunities_ccp = first(EntropyCommunities_ccp), 
            EntropyCommunities_cp = first(EntropyCommunities_cp)) %>%
  as.data.table() %>%
  ungroup() %>%
  mutate(EntropyCommunities_ccp = EntropyCommunities_ccp/mean(EntropyCommunities_ccp, na.rm=TRUE),
         EntropyCommunities_cp = EntropyCommunities_cp/mean(EntropyCommunities_cp, na.rm=TRUE))

#==========================================================================
# plotting the decline in the number of communities

# -----------------------------
namefile <- "EComm_0002_Worldcoherence.png"
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
# step AIC on levels
library(MASS)
summary(newdf)
colnames(newdf)
fullglm <- glm(loggdppc ~ ., data=subset(newdf, select = -c(exporter, year, 
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


