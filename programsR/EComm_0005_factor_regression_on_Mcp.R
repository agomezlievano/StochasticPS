# EComm_0005_factor_regression_on_Mcp.R

# #########################################################
# Here we 
# Apply the idea that by looking at what places have,
# one can decompose the data into ``complexities''.
# 
# log(prob_cp) = M(1-rc)log(1-tp)
# #########################################################

rm(list = ls(all = TRUE))  # resets R to fresh
gc()

# ========================================================
# Libraries
# --------------------------------------------------------
library(sfsmisc)
library(scales)
library(data.table)
library(Matrix)
library(ggplot2)
library(stargazer)
library(fitdistrplus)
library(evd)
library(lfe)

# ========================================================
# Some declarations
# --------------------------------------------------------
workingpath <- "~/shared_space/cidgrowlab/AGL/StochasticPS"


#  Declare other paths
outpath <- paste0(workingpath,"/outputdata")
figpath <- paste0(workingpath,"/figures")
tablespath <- paste0(workingpath,"/tables")
inpath <- paste0(workingpath,"/inputdata")

setwd(workingpath)

# --------------------------------------------------------
# useful functions
colfunc <- colorRampPalette(c("deepskyblue4", "deepskyblue", "cyan"))
col7 <- colorRampPalette(c("#f7f7f7", "#74add1", "#2166ac", "#053061")) 



# ========================================================
# Load data
# --------------------------------------------------------

# edgelist of Mcp
inputfile.name <- "/EComm_0004_syntheticMcp_C200_P800.csv"
Mcp_longdf <- rio::import(file = paste0(inpath, inputfile.name))

Mcp_longdf$country_code <- as.factor(Mcp_longdf$country_code)
Mcp_longdf$product_code <- as.factor(paste0("P",Mcp_longdf$product_code))
Mcp_longdf$logmlogMcp <- 1.0 - 2.0*Mcp_longdf$P
head(Mcp_longdf)

cty_vec <- levels(Mcp_longdf$country_code)
prod_vec <- levels(Mcp_longdf$product_code)
length(cty_vec) + length(prod_vec)

# ========================================================
# Regressions
# --------------------------------------------------------

# with fixed effects
est2 <- felm(logmlogMcp  ~  1 | 
               country_code + product_code, 
             data=Mcp_longdf)
sumest2 <- summary(est2)
objects(est2)

est2fe <- est2$fe
names(est2fe)
length(est2fe$country_code)
length(est2fe$product_code)
head(est2fe$country_code)
head(est2fe$product_code)
est2$p

pred <- est2$fitted.values

estimatedfe2 <- getfe(est2)
head(estimatedfe2)
tail(estimatedfe2)



# with fixed effects and clustered errors
est3 <- felm(logmlogMcp  ~  1 | 
               country_code + product_code | 0 | country_code + product_code, 
             data=Mcp_longdf)
sumest3 <- summary(est3)

estimatedfe3 <- getfe(est3)
head(estimatedfe3)
tail(estimatedfe3)


# with fixed effects and clustered errors
est4 <- felm(P  ~  1 | 
               country_code + product_code | 0 | country_code + product_code, 
             data=Mcp_longdf)
sumest4 <- summary(est4)

estimatedfe4 <- getfe(est4)
head(estimatedfe4)
tail(estimatedfe4)


# ========================================================
# EXPORTING
# --------------------------------------------------------

outputfile.name <- "/EComm_0005_factors_estimated_C200_P800_P.csv"
rio::export(estimatedfe4, file = paste0(outpath, outputfile.name))


