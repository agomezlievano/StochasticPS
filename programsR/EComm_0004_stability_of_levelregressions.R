# EComm_0004_stability_of_levelregressions.R

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



# ============================================
# LOADING DATA
# --------------------------------------------

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

# --------------------------------------------


# ============================================
# plotting the decline in the number of communities
# --------------------------------------------

# --------------------------------------------




# ============================================
# plotting the decline in the number of communities
# --------------------------------------------

# --------------------------------------------




# ============================================
# plotting the decline in the number of communities
# --------------------------------------------

# --------------------------------------------



# ============================================
# plotting the decline in the number of communities
# --------------------------------------------

# --------------------------------------------
