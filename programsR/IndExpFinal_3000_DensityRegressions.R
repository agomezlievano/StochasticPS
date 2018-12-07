# IndExpFinal_3000_DensityRegressions.R

# #########################################################
# Here we 
# Do the density regressions and export the LaTeX tables.
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
library(entropy)
library(infotheo)
library(diverse)
library(treemap)
library(pheatmap)
library(ggplot2)
library(stargazer)
library(fitdistrplus)
library(evd)
library(lfe)

# ========================================================
# Some declarations
# --------------------------------------------------------
workingpath <- "~/shared_space/cidgrowlab/Colombia/IndustryExports/Andres"
sidworkingpath <- "~/shared_space/cidgrowlab/Colombia/IndustryExports/Sid"

#  Declare other paths
outpath <- paste0(workingpath,"/outputdata")
dopath <- paste0(workingpath,"/dofiles")
temppath <- paste0(workingpath,"/temps")
figpath <- paste0(workingpath,"/figures")
tablespath <- paste0(workingpath,"/tables")
rcodepath <- paste0(workingpath,"/Rprograms")
sidrcodepath <- paste0(sidworkingpath,"/Rprograms")
inpath <- paste0(workingpath,"/inputdata")

setwd(rcodepath)
source("EC_999_functions.R")

# Source the utils files
source(paste0(sidrcodepath,"/00_utils.R"))

# --------------------------------------------------------
# useful functions
colfunc <- colorRampPalette(c("deepskyblue4", "deepskyblue", "cyan"))
col7 <- colorRampPalette(c("#f7f7f7", "#74add1", "#2166ac", "#053061")) 



# ========================================================
# Load data
# --------------------------------------------------------

# City-product dataset
inputfile.name <- "/IndExpFinal_2000_Cities_with_densities_and_potentials_with_changes.dta"
cpy_df <- rio::import(file = paste0(outpath, inputfile.name))



# Definition of cities
locs_df <- rio::import(file=paste0(inpath,"/Colombia_city_key.csv"))


# ========================================================
# Getting some info of entities
# --------------------------------------------------------
years_vec <- unique(cpy_df$year)
citis_vec <- unique(cpy_df$city_code)
prods_vec <- unique(cpy_df$p)

years_vec <- years_vec[(years_vec!="" & years_vec!=".") & !is.na(years_vec)]
citis_vec <- citis_vec[(citis_vec!="" & citis_vec!=".") & !is.na(citis_vec)]
prods_vec <- prods_vec[(prods_vec!="" & prods_vec!=".") & !is.na(prods_vec)]


length(years_vec)
length(citis_vec)
length(prods_vec)

# ========================================================
# SETTING UP THE LOOPS FOR REGRESSIONS
# --------------------------------------------------------
# Creating factor variables
cpy_df$year <- factor(cpy_df$year)
cpy_df$city_code <- factor(cpy_df$city_code)
cpy_df$p <- factor(cpy_df$p)


# PAIRWISE CORRELATIONS BETWEEN DENSITIES
stargazer(cor(cpy_df[,substr(colnames(cpy_df),1,2)=="d_"], use="complete.obs"),
          out = paste0(tablespath,"/IndExpFinal_3000_DensityCors.tex"))


head(cpy_df[cpy_df$year==2010 & !is.na(cpy_df$d_worker_cp) & !is.na(cpy_df$pslogrca_mod),
               c("pslogrca_mod","employment","num_firms","d_worker_cp","d_city_cp","d_trad1_cp","d_trad2_cp")])

summary(data.frame(cpy_df[cpy_df$year==2010 & !is.na(cpy_df$d_worker_cp) & !is.na(cpy_df$pslogrca_mod),
           c("pslogrca_mod","employment","num_firms","d_worker_cp","d_city_cp","d_trad1_cp","d_trad2_cp")]))
cor(cpy_df[cpy_df$year==2010 & !is.na(cpy_df$d_worker_cp) & !is.na(cpy_df$pslogrca_mod),
           c("pslogrca_mod","employment","num_firms","d_worker_cp","d_city_cp","d_trad1_cp","d_trad2_cp")], 
    use = "complete.obs")
plot(cpy_df[cpy_df$year==2010 & !is.na(cpy_df$d_worker_cp) & !is.na(cpy_df$pslogrca_mod),
            c("pslogrca_mod","employment","num_firms","d_worker_cp","d_city_cp","d_trad1_cp","d_trad2_cp")])

pairs(cpy_df[cpy_df$year==2010 & !is.na(cpy_df$d_worker_cp) & !is.na(cpy_df$pslogrca_mod),
             c("pslogrca_mod","employment","num_firms","d_worker_cp","d_city_cp","d_trad1_cp","d_trad2_cp")], 
      col=rgb(1, 0.1, 0.5), pch=19, cex=0.5,
      diag.panel = panel.hist,
      lower.panel = panel.loess)


# ############### #
# What to report?
# For each of the three measures:
# 1) Include all four densities and do regressions for all 5 time windows
# 2) Include only the best density and do regressions for all 5 time windows
# 3) Pick the best density and the longest time window, but include the fixed effects sequentially.

# ------------------------
# # without fixed effects
# est <- felm(modrca_change_1yr  ~  pslogrca_mod + d_worker_cp + d_city_cp + d_trad1_cp + d_trad2_cp, 
#             data=cpy_df)
# summary(est)
# est <- felm(modrca_change_1yr  ~  pslogrca_mod + d_worker_cp + d_city_cp + d_trad1_cp + d_trad2_cp |
#               0, 
#             data=cpy_df)
# summary(est)
# 
# # with fixed effects
# est2 <- felm(modrca_change_1yr  ~  pslogrca_mod + d_worker_cp + d_city_cp + d_trad1_cp + d_trad2_cp | 
#                city_code + p, 
#              data=cpy_df)
# summary(est2)
# 
# with fixed effects and clustered errors
est3 <- felm(modrca_change_1yr  ~  pslogrca_mod + d_worker_cp + d_city_cp + d_trad1_cp + d_trad2_cp | 
               city_code + p | 0 | city_code + p, 
             data=cpy_df)
sumest3 <- summary(est3)
# ------------------------

# ############### #


measures <- c("modrca_change", "emp_change", "numfirms_change")
names(measures) <- c("pslogrca_mod", "employment", "num_firms")
windows <- c(1:5)
densities <- c("alldens", "d_worker_cp", "d_city_cp", "d_trad1_cp", "d_trad2_cp")
fixedeffects <- c("nofe", "city_code", "p", "city_code + p")

# -----
elasticities <- FALSE
standardize <- TRUE
if(elasticities){
  for(var in densities[2:5]){
    cpy_df[, var] <- cpy_df[, var] - min(cpy_df[, var], na.rm=TRUE)
    cpy_df[cpy_df[ , var]==0 & !(is.na(cpy_df[, var])), var] <- rep(NA, nrow(cpy_df[cpy_df[ , var]==0  & !(is.na(cpy_df[, var])), ])) 
    cpy_df[cpy_df[ , var]>0  & !(is.na(cpy_df[, var])), var] <- log(cpy_df[cpy_df[ , var]>0  & !(is.na(cpy_df[, var])), var])
  }
}
if(standardize){
  for(var in c(names(measures),
               densities[2:5],
               paste0(apply(expand.grid(measures,windows), 
                            MARGIN=1, 
                            FUN=function(x)paste(x, collapse="_")),
                      "yr"))){
    cpy_df[ , var] <- as.numeric(scale(cpy_df[ , var], center = TRUE, scale = TRUE))
  }
}
# -----

allformulas_df <- data.frame(Variable = character(),
                             DependentVar = character(),
                             TimeWindow = integer(),
                             Densities = character(),
                             WhichDensity = character(),
                             FixedEffect = character(),
                             Formula = character())
for(meas in names(measures)){
  for(w in windows){
    depvar <- paste0(measures[meas], "_", w, "yr")
    
    for(fe in fixedeffects){
      if(fe=="nofe"){
        fevars <- 0
      } else {
        fevars <- fe
      }
      
      for(dens in densities){
        if(dens=="alldens"){
          densvars <- paste(densities[2:5], collapse=" + ")
        } else {
          densvars <- dens
          dens <- "singledens"
        }
        
        
        # Adding regression to the mean term
        indepvars <- paste(meas, densvars, sep=" + ")
        
        # the 'lfe' package requires formulas like: y ~ x1 + x2 | f1 + f2 | (Q|W ~ x3+x4) | clu1 + clu2
        
        # Creating the full formula
        formiidf <- data.frame(Variable = meas,
                               DependentVar = depvar,
                               TimeWindow = w,
                               Densities = dens,
                               WhichDensity = densvars,
                               FixedEffect = fe,
                               Formula = paste(depvar, " ~ ", paste(indepvars, fevars, sep=" | ")))
        
        
        allformulas_df <- data.frame(rbind(allformulas_df, formiidf))
      }
    }
  }  
}


# ========================================================
# REGRESSIONS
# --------------------------------------------------------

# ----------------------
print("STARTING TO MEASURE TIME")
start.time <- Sys.time()
# ----------------------

# -------------------------
library(foreach)
library(doParallel)
library(doRNG)
numcores <- detectCores()
print(numcores)
cl <- makeCluster(10) #not to overload your computer
registerDoParallel(cl)
registerDoRNG(seed = 1984)
# -------------------------c(1:nrow(allformulas_df))
returnlist <- foreach(i=1:nrow(allformulas_df), 
                      .errorhandling='pass',
                      .export=c("cpy_df",  
                                "allformulas_df",
                                "densities"),
                      .packages = c("lfe", "data.table")) %dopar% {
                        
                        # -------------------------------------------------
                        # RUNNING THE REGRESSION
                        myfit_i <- felm(as.formula(as.character(allformulas_df[i,"Formula"])), data=cpy_df)
                        sum_myfit_i <- summary(myfit_i)
                        
                        # -------------------------------------------------
                        # CREATING DATA.FRAME
                        dfrowi <- data.frame(Variable = allformulas_df[i,"Variable"],
                                             DependentVar = allformulas_df[i,"DependentVar"],
                                             TimeWindow = allformulas_df[i,"TimeWindow"],
                                             Densities = allformulas_df[i,"Densities"],
                                             WhichDensity = allformulas_df[i,"WhichDensity"],
                                             FixedEffect = allformulas_df[i,"FixedEffect"],
                                             Formula = allformulas_df[i,"Formula"],
                                             N = sum_myfit_i$N,
                                             R2 = sum_myfit_i$r2adj,
                                             ReversionToMean = as.numeric(NA),
                                             se_ReversionToMean = as.numeric(NA),
                                             pval_ReversionToMean = as.numeric(NA),
                                             tstat_ReversionToMean = as.numeric(NA),
                                             d_worker_cp = as.numeric(NA),
                                             d_city_cp = as.numeric(NA),
                                             d_trad1_cp = as.numeric(NA),
                                             d_trad2_cp = as.numeric(NA),
                                             se_d_worker_cp = as.numeric(NA),
                                             se_d_city_cp = as.numeric(NA),
                                             se_d_trad1_cp = as.numeric(NA),
                                             se_d_trad2_cp = as.numeric(NA),
                                             upper_d_worker_cp = as.numeric(NA),
                                             upper_d_city_cp = as.numeric(NA),
                                             upper_d_trad1_cp = as.numeric(NA),
                                             upper_d_trad2_cp = as.numeric(NA),
                                             lower_d_worker_cp = as.numeric(NA),
                                             lower_d_city_cp = as.numeric(NA),
                                             lower_d_trad1_cp = as.numeric(NA),
                                             lower_d_trad2_cp = as.numeric(NA),
                                             pval_d_worker_cp = as.numeric(NA),
                                             pval_d_city_cp = as.numeric(NA),
                                             pval_d_trad1_cp = as.numeric(NA),
                                             pval_d_trad2_cp = as.numeric(NA),
                                             tstat_d_worker_cp = as.numeric(NA),
                                             tstat_d_city_cp = as.numeric(NA),
                                             tstat_d_trad1_cp = as.numeric(NA),
                                             tstat_d_trad2_cp = as.numeric(NA))
                        
                        if(allformulas_df[i,"Densities"]=="alldens"){
                          for(densvar in densities[2:5]){
                            dfrowi[1,"ReversionToMean"] <- sum_myfit_i$coefficients[1,"Estimate"]
                            dfrowi[1,"se_ReversionToMean"] <- sum_myfit_i$coefficients[1,"Std. Error"]
                            dfrowi[1,"pval_ReversionToMean"] <- sum_myfit_i$coefficients[1,"Pr(>|t|)"]
                            dfrowi[1,"tstat_ReversionToMean"] <- sum_myfit_i$coefficients[1,"t value"]
                            dfrowi[1,as.character(densvar)] <- sum_myfit_i$coefficients[as.character(densvar),"Estimate"]
                            dfrowi[1,paste0("se_",as.character(densvar))] <- sum_myfit_i$coefficients[as.character(densvar),"Std. Error"]
                            dfrowi[1,paste0("pval_",as.character(densvar))] <- sum_myfit_i$coefficients[as.character(densvar),"Pr(>|t|)"]
                            dfrowi[1,paste0("tstat_",as.character(densvar))] <- sum_myfit_i$coefficients[as.character(densvar),"t value"]
                            dfrowi[1,paste0("upper_",as.character(densvar))] <- dfrowi[1,as.character(densvar)] + 1.96*dfrowi[1,paste0("se_",as.character(densvar))]
                            dfrowi[1,paste0("lower_",as.character(densvar))] <- dfrowi[1,as.character(densvar)] - 1.96*dfrowi[1,paste0("se_",as.character(densvar))]
                          }
                        } else {
                          densvar <- as.character(allformulas_df[i,"WhichDensity"])
                          
                          dfrowi[1,"ReversionToMean"] <- sum_myfit_i$coefficients[1,"Estimate"]
                          dfrowi[1,"se_ReversionToMean"] <- sum_myfit_i$coefficients[1,"Std. Error"]
                          dfrowi[1,"pval_ReversionToMean"] <- sum_myfit_i$coefficients[1,"Pr(>|t|)"]
                          dfrowi[1,"tstat_ReversionToMean"] <- sum_myfit_i$coefficients[1,"t value"]
                          dfrowi[1,densvar] <- sum_myfit_i$coefficients[densvar,"Estimate"]
                          dfrowi[1,paste0("se_",densvar)] <- sum_myfit_i$coefficients[densvar,"Std. Error"]    
                          dfrowi[1,paste0("pval_",densvar)] <- sum_myfit_i$coefficients[densvar,"Pr(>|t|)"]    
                          dfrowi[1,paste0("tstat_",densvar)] <- sum_myfit_i$coefficients[densvar,"t value"]
                          dfrowi[1,paste0("upper_",densvar)] <- dfrowi[1,densvar] + 1.96*dfrowi[1,paste0("se_",densvar)]
                          dfrowi[1,paste0("lower_",densvar)] <- dfrowi[1,densvar] - 1.96*dfrowi[1,paste0("se_",densvar)]
                        }
                        
                        
                        # -------------------------------------------------
                        # APPENDING
                        dfrowi
                      }

newdf <- rbindlist(returnlist)
newdf$reg <- as.character(c(1:nrow(newdf)))
head(newdf)


# -----------------------
#stop cluster
stopCluster(cl)
# -----------------------

# ----------------------
print("END FOR LOOP")
end.time <- Sys.time()
time.taken <- end.time - start.time
print("TIME TAKEN:")
print(time.taken)
# ----------------------


numbydepvar <- nrow(newdf[newdf$DependentVar=="modrca_change_1yr",])
uniquedepvar <- as.character(unique(newdf$DependentVar))
newdf$regwithin <- numeric(nrow(newdf))
for(dvar in uniquedepvar){
    newdf[newdf$DependentVar==dvar, "regwithin"] <- seq(1,numbydepvar) - 0.5
}



# ########################################################
# ========================================================
# GGPLOTS
# --------------------------------------------------------
myshapes <- c(15,18,16,17)
# R2
ggplt.R2 <- ggplot(newdf, aes(x=regwithin, y=R2, colour=Variable, shape=FixedEffect, size=TimeWindow)) +
#   theme(panel.grid.major = element_line(colour="gray70", size=0.5, linetype = "dashed"),
#         panel.grid.minor = element_line(colour="gray80", size=0.5, linetype = "dotted")) +
  scale_x_continuous(name="Regression specification") + 
  scale_y_continuous(name=expression(paste("Adj.", R^2)), 
                     limits = c(0,0.4),
                     breaks = seq(0,1,0.2), 
                     minor_breaks=seq(0,1,0.05)) + 
  scale_color_brewer(palette="Set1") +
  geom_hline(yintercept = 0, alpha=0.5, size=1.5, color="black") + 
  geom_point(alpha=0.5) + 
  scale_shape_manual(values=myshapes) +
  theme_bw(base_size = 12, base_family = "Helvetica")

# Reversion to the mean
ggplt.mrtm <- ggplot(newdf, aes(x=regwithin, y=tstat_ReversionToMean, colour=Variable, shape=FixedEffect, size=TimeWindow)) +
  #   theme(panel.grid.major = element_line(colour="gray70", size=0.5, linetype = "dashed"),
  #         panel.grid.minor = element_line(colour="gray80", size=0.5, linetype = "dotted")) +
  scale_x_continuous(name="Regression specification") + 
  scale_y_continuous(name="Reversion to the Mean t-stat") + 
  scale_color_brewer(palette="Set1") +
  geom_hline(yintercept = 0, alpha=0.5, size=1.5, color="black") + 
  geom_point(alpha=0.5) + 
  scale_shape_manual(values=myshapes) +
  theme_bw(base_size = 12, base_family = "Helvetica")

# -------------------------
# estimate
# worker
ggplt.dworker.est <- ggplot(newdf[newdf$FixedEffect=="city_code + p",], aes(x=regwithin, y=d_worker_cp, colour=Variable, shape=FixedEffect)) +
  scale_x_continuous(name="Regression specification") + 
  #scale_y_continuous(name="Std. Coef. of density type 1", limits=c(-0.7, 0.7), breaks=seq(-0.7,0.7,0.1)) + 
  #scale_y_continuous(name="Coef. of density type 1") + 
  scale_y_continuous(name="Std. Coef. of density type 1", limits=c(-1, 2), breaks=seq(-1,2,0.2)) + 
  geom_hline(yintercept = 0, alpha=0.5, size=1.5, color="black") + 
  geom_errorbar(aes(ymin = lower_d_worker_cp, ymax=upper_d_worker_cp)) + 
  geom_point(aes(size=TimeWindow), alpha=0.5) + 
  scale_color_brewer(palette="Set1") +
  scale_shape_manual(values=17) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(legend.position="none")

# city
ggplt.dcity.est <- ggplot(newdf[newdf$FixedEffect=="city_code + p",], aes(x=regwithin, y=d_city_cp, colour=Variable, shape=FixedEffect)) +
#   theme(panel.grid.major = element_line(colour="gray70", size=0.5, linetype = "dashed"),
#         panel.grid.minor = element_line(colour="gray80", size=0.5, linetype = "dotted")) +
  scale_x_continuous(name="Regression specification") + 
  #scale_y_continuous(name="Std. Coef. of density type 2", limits=c(-0.7, 0.7), breaks=seq(-0.7,0.7,0.1)) + 
  #scale_y_continuous(name="Coef. of density type 2") + 
  scale_y_continuous(name="Std. Coef. of density type 2", limits=c(-1, 2), breaks=seq(-1,2,0.2)) + 
  geom_hline(yintercept = 0, alpha=0.5, size=1.5, color="black") + 
  geom_errorbar(aes(ymin = lower_d_city_cp, ymax=upper_d_city_cp)) + 
  geom_point(aes(size=TimeWindow), alpha=0.5) + 
  scale_color_brewer(palette="Set1") +
  scale_shape_manual(values=17) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(legend.position="none")

# trad1
ggplt.dtrad1.est <- ggplot(newdf[newdf$FixedEffect=="city_code + p",], aes(x=regwithin, y=d_trad1_cp, colour=Variable, shape=FixedEffect)) +
#   theme(panel.grid.major = element_line(colour="gray70", size=0.5, linetype = "dashed"),
#         panel.grid.minor = element_line(colour="gray80", size=0.5, linetype = "dotted")) +
  scale_x_continuous(name="Regression specification") + 
  #scale_y_continuous(name="Std. Coef. of density type 3", limits=c(-0.7, 0.7), breaks=seq(-0.7,0.7,0.1)) + 
  #scale_y_continuous(name="Coef. of density type 3") + 
  scale_y_continuous(name="Std. Coef. of density type 3", limits=c(-1, 2), breaks=seq(-1,2,0.2)) + 
  geom_hline(yintercept = 0, alpha=0.5, size=1.5, color="black") + 
  geom_errorbar(aes(ymin = lower_d_trad1_cp, ymax=upper_d_trad1_cp)) + 
  geom_point(aes(size=TimeWindow), alpha=0.5) + 
  scale_color_brewer(palette="Set1") +
  scale_shape_manual(values=17) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(legend.position="none")

# trad2
ggplt.dtrad2.est <- ggplot(newdf[newdf$FixedEffect=="city_code + p",], aes(x=regwithin, y=d_trad2_cp, colour=Variable, shape=FixedEffect)) +
#   theme(panel.grid.major = element_line(colour="gray70", size=0.5, linetype = "dashed"),
#         panel.grid.minor = element_line(colour="gray80", size=0.5, linetype = "dotted")) +
  scale_x_continuous(name="Regression specification") + 
  #scale_y_continuous(name="Std. Coef. of density type 4", limits=c(-0.7, 0.7), breaks=seq(-0.7,0.7,0.1)) + 
  #scale_y_continuous(name="Coef. of density type 4") + 
  scale_y_continuous(name="Std. Coef. of density type 4", limits=c(-1, 2), breaks=seq(-1,2,0.2)) + 
  geom_hline(yintercept = 0, alpha=0.5, size=1.5, color="black") + 
  geom_errorbar(aes(ymin = lower_d_trad2_cp, ymax=upper_d_trad2_cp)) + 
  geom_point(aes(size=TimeWindow), alpha=0.5) + 
  scale_color_brewer(palette="Set1") +
  scale_shape_manual(values=17) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(legend.position="none")

# -------------------------



gga.bottom <- ggarrange(ggplt.dworker.est, ggplt.dcity.est, ggplt.dtrad1.est, ggplt.dtrad2.est, 
                        #common.legend = TRUE,
                        #labels = c("B", "C", "D", "E"),
                        labels = c("C", "D", "E", "F"),
                        ncol = 4, nrow = 1,
                        widths = c(1,1,1,1), heights=1)

gga <- ggarrange(ggplt.R2,
                 ggplt.mrtm,
                 gga.bottom, 
                 common.legend = TRUE,
                 labels = c("A", "B", ""),
                 ncol = 1, nrow = 3,
                 widths = 2, heights=c(1,1,1))
gga

# ---------------------------------------------
namefile <- "/IndExpFinal_3000_DensityRegressionsPlots.png"
png(paste0(figpath,namefile), 
    height=12, 
    width=12, 
    units="in", 
    res=300)
# ---------------------------------------------

print(gga)
# ---------------------------------------------
dev.off()
# ---------------------------------------------




library(RColorBrewer)
brewer.pal(n=3, name="Set1")

# ########################################
# modRCA
# -------------------------
# Only modRCA
# estimate
# worker
ggplt.dworker.est.modRCA <- ggplot(newdf[newdf$FixedEffect=="city_code + p" & newdf$Variable=="pslogrca_mod",], 
                            aes(x=regwithin, y=d_worker_cp, colour=Variable, shape=FixedEffect)) +
  scale_x_continuous(name="Regression specification") + 
  #scale_y_continuous(name="Std. Coef. of density type 1", limits=c(-0.7, 0.7), breaks=seq(-0.7,0.7,0.1)) + 
  #scale_y_continuous(name="Coef. of density type 1") + 
  scale_y_continuous(name="Std. Coef. of density type 1", limits=c(-1, 2), breaks=seq(-1,2,0.2)) + 
  geom_hline(yintercept = 0, alpha=0.5, size=1.5, color="black") + 
  geom_errorbar(aes(ymin = lower_d_worker_cp, ymax=upper_d_worker_cp)) + 
  geom_point(aes(size=TimeWindow), alpha=0.5) + 
  scale_color_manual(values="#E41A1C") +
  scale_shape_manual(values=17) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(legend.position="none")

# city
ggplt.dcity.est.modRCA <- ggplot(newdf[newdf$FixedEffect=="city_code + p" & newdf$Variable=="pslogrca_mod",], 
                          aes(x=regwithin, y=d_city_cp, colour=Variable, shape=FixedEffect)) +
  #   theme(panel.grid.major = element_line(colour="gray70", size=0.5, linetype = "dashed"),
  #         panel.grid.minor = element_line(colour="gray80", size=0.5, linetype = "dotted")) +
  scale_x_continuous(name="Regression specification") + 
  #scale_y_continuous(name="Std. Coef. of density type 2", limits=c(-0.7, 0.7), breaks=seq(-0.7,0.7,0.1)) + 
  #scale_y_continuous(name="Coef. of density type 2") + 
  scale_y_continuous(name="Std. Coef. of density type 2", limits=c(-1, 2), breaks=seq(-1,2,0.2)) + 
  geom_hline(yintercept = 0, alpha=0.5, size=1.5, color="black") + 
  geom_errorbar(aes(ymin = lower_d_city_cp, ymax=upper_d_city_cp)) + 
  geom_point(aes(size=TimeWindow), alpha=0.5) + 
  scale_color_manual(values="#E41A1C") +
  scale_shape_manual(values=17) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(legend.position="none")

# trad1
ggplt.dtrad1.est.modRCA <- ggplot(newdf[newdf$FixedEffect=="city_code + p" & newdf$Variable=="pslogrca_mod",], 
                           aes(x=regwithin, y=d_trad1_cp, colour=Variable, shape=FixedEffect)) +
  #   theme(panel.grid.major = element_line(colour="gray70", size=0.5, linetype = "dashed"),
  #         panel.grid.minor = element_line(colour="gray80", size=0.5, linetype = "dotted")) +
  scale_x_continuous(name="Regression specification") + 
  #scale_y_continuous(name="Std. Coef. of density type 3", limits=c(-0.7, 0.7), breaks=seq(-0.7,0.7,0.1)) + 
  #scale_y_continuous(name="Coef. of density type 3") + 
  scale_y_continuous(name="Std. Coef. of density type 3", limits=c(-1, 2), breaks=seq(-1,2,0.2)) + 
  geom_hline(yintercept = 0, alpha=0.5, size=1.5, color="black") + 
  geom_errorbar(aes(ymin = lower_d_trad1_cp, ymax=upper_d_trad1_cp)) + 
  geom_point(aes(size=TimeWindow), alpha=0.5) + 
  scale_color_manual(values="#E41A1C") +
  scale_shape_manual(values=17) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(legend.position="none")

# trad2
ggplt.dtrad2.est.modRCA <- ggplot(newdf[newdf$FixedEffect=="city_code + p" & newdf$Variable=="pslogrca_mod",], 
                           aes(x=regwithin, y=d_trad2_cp, colour=Variable, shape=FixedEffect)) +
  #   theme(panel.grid.major = element_line(colour="gray70", size=0.5, linetype = "dashed"),
  #         panel.grid.minor = element_line(colour="gray80", size=0.5, linetype = "dotted")) +
  scale_x_continuous(name="Regression specification") + 
  #scale_y_continuous(name="Std. Coef. of density type 4", limits=c(-0.7, 0.7), breaks=seq(-0.7,0.7,0.1)) + 
  #scale_y_continuous(name="Coef. of density type 4") + 
  scale_y_continuous(name="Std. Coef. of density type 4", limits=c(-1, 2), breaks=seq(-1,2,0.2)) + 
  geom_hline(yintercept = 0, alpha=0.5, size=1.5, color="black") + 
  geom_errorbar(aes(ymin = lower_d_trad2_cp, ymax=upper_d_trad2_cp)) + 
  geom_point(aes(size=TimeWindow), alpha=0.5) + 
  scale_color_manual(values="#E41A1C") +
  scale_shape_manual(values=17) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(legend.position="none")
# -------------------------


# emp
# -------------------------
# Only employment
# estimate
# worker
ggplt.dworker.est.emp <- ggplot(newdf[newdf$FixedEffect=="city_code + p" & newdf$Variable=="employment",], 
                                   aes(x=regwithin, y=d_worker_cp, colour=Variable, shape=FixedEffect)) +
  scale_x_continuous(name="Regression specification") + 
  #scale_y_continuous(name="Std. Coef. of density type 1", limits=c(-0.7, 0.7), breaks=seq(-0.7,0.7,0.1)) + 
  #scale_y_continuous(name="Coef. of density type 1") + 
  scale_y_continuous(name="Std. Coef. of density type 1", limits=c(-1, 2), breaks=seq(-1,2,0.2)) + 
  geom_hline(yintercept = 0, alpha=0.5, size=1.5, color="black") + 
  geom_errorbar(aes(ymin = lower_d_worker_cp, ymax=upper_d_worker_cp)) + 
  geom_point(aes(size=TimeWindow), alpha=0.5) + 
  scale_color_manual(values="#377EB8") +
  scale_shape_manual(values=17) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(legend.position="none")

# city
ggplt.dcity.est.emp <- ggplot(newdf[newdf$FixedEffect=="city_code + p" & newdf$Variable=="employment",], 
                                 aes(x=regwithin, y=d_city_cp, colour=Variable, shape=FixedEffect)) +
  #   theme(panel.grid.major = element_line(colour="gray70", size=0.5, linetype = "dashed"),
  #         panel.grid.minor = element_line(colour="gray80", size=0.5, linetype = "dotted")) +
  scale_x_continuous(name="Regression specification") + 
  #scale_y_continuous(name="Std. Coef. of density type 2", limits=c(-0.7, 0.7), breaks=seq(-0.7,0.7,0.1)) + 
  #scale_y_continuous(name="Coef. of density type 2") + 
  scale_y_continuous(name="Std. Coef. of density type 2", limits=c(-1, 2), breaks=seq(-1,2,0.2)) + 
  geom_hline(yintercept = 0, alpha=0.5, size=1.5, color="black") + 
  geom_errorbar(aes(ymin = lower_d_city_cp, ymax=upper_d_city_cp)) + 
  geom_point(aes(size=TimeWindow), alpha=0.5) + 
  scale_color_manual(values="#377EB8") +
  scale_shape_manual(values=17) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(legend.position="none")

# trad1
ggplt.dtrad1.est.emp <- ggplot(newdf[newdf$FixedEffect=="city_code + p" & newdf$Variable=="employment",], 
                                  aes(x=regwithin, y=d_trad1_cp, colour=Variable, shape=FixedEffect)) +
  #   theme(panel.grid.major = element_line(colour="gray70", size=0.5, linetype = "dashed"),
  #         panel.grid.minor = element_line(colour="gray80", size=0.5, linetype = "dotted")) +
  scale_x_continuous(name="Regression specification") + 
  #scale_y_continuous(name="Std. Coef. of density type 3", limits=c(-0.7, 0.7), breaks=seq(-0.7,0.7,0.1)) + 
  #scale_y_continuous(name="Coef. of density type 3") + 
  scale_y_continuous(name="Std. Coef. of density type 3", limits=c(-1, 2), breaks=seq(-1,2,0.2)) + 
  geom_hline(yintercept = 0, alpha=0.5, size=1.5, color="black") + 
  geom_errorbar(aes(ymin = lower_d_trad1_cp, ymax=upper_d_trad1_cp)) + 
  geom_point(aes(size=TimeWindow), alpha=0.5) + 
  scale_color_manual(values="#377EB8") +
  scale_shape_manual(values=17) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(legend.position="none")

# trad2
ggplt.dtrad2.est.emp <- ggplot(newdf[newdf$FixedEffect=="city_code + p" & newdf$Variable=="employment",], 
                                  aes(x=regwithin, y=d_trad2_cp, colour=Variable, shape=FixedEffect)) +
  #   theme(panel.grid.major = element_line(colour="gray70", size=0.5, linetype = "dashed"),
  #         panel.grid.minor = element_line(colour="gray80", size=0.5, linetype = "dotted")) +
  scale_x_continuous(name="Regression specification") + 
  #scale_y_continuous(name="Std. Coef. of density type 4", limits=c(-0.7, 0.7), breaks=seq(-0.7,0.7,0.1)) + 
  #scale_y_continuous(name="Coef. of density type 4") + 
  scale_y_continuous(name="Std. Coef. of density type 4", limits=c(-1, 2), breaks=seq(-1,2,0.2)) + 
  geom_hline(yintercept = 0, alpha=0.5, size=1.5, color="black") + 
  geom_errorbar(aes(ymin = lower_d_trad2_cp, ymax=upper_d_trad2_cp)) + 
  geom_point(aes(size=TimeWindow), alpha=0.5) + 
  scale_color_manual(values="#377EB8") +
  scale_shape_manual(values=17) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(legend.position="none")
# -------------------------

# nf
# -------------------------
# Only num_firms
# estimate
# worker
ggplt.dworker.est.nf <- ggplot(newdf[newdf$FixedEffect=="city_code + p" & newdf$Variable=="num_firms",], 
                                aes(x=regwithin, y=d_worker_cp, colour=Variable, shape=FixedEffect)) +
  scale_x_continuous(name="Regression specification") + 
  #scale_y_continuous(name="Std. Coef. of density type 1", limits=c(-0.7, 0.7), breaks=seq(-0.7,0.7,0.1)) + 
  #scale_y_continuous(name="Coef. of density type 1") + 
  scale_y_continuous(name="Std. Coef. of density type 1", limits=c(-1, 2), breaks=seq(-1,2,0.2)) + 
  geom_hline(yintercept = 0, alpha=0.5, size=1.5, color="black") + 
  geom_errorbar(aes(ymin = lower_d_worker_cp, ymax=upper_d_worker_cp)) + 
  geom_point(aes(size=TimeWindow), alpha=0.5) + 
  scale_color_manual(values="#4DAF4A") +
  scale_shape_manual(values=17) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(legend.position="none")

# city
ggplt.dcity.est.nf <- ggplot(newdf[newdf$FixedEffect=="city_code + p" & newdf$Variable=="num_firms",], 
                              aes(x=regwithin, y=d_city_cp, colour=Variable, shape=FixedEffect)) +
  #   theme(panel.grid.major = element_line(colour="gray70", size=0.5, linetype = "dashed"),
  #         panel.grid.minor = element_line(colour="gray80", size=0.5, linetype = "dotted")) +
  scale_x_continuous(name="Regression specification") + 
  #scale_y_continuous(name="Std. Coef. of density type 2", limits=c(-0.7, 0.7), breaks=seq(-0.7,0.7,0.1)) + 
  #scale_y_continuous(name="Coef. of density type 2") + 
  scale_y_continuous(name="Std. Coef. of density type 2", limits=c(-1, 2), breaks=seq(-1,2,0.2)) + 
  geom_hline(yintercept = 0, alpha=0.5, size=1.5, color="black") + 
  geom_errorbar(aes(ymin = lower_d_city_cp, ymax=upper_d_city_cp)) + 
  geom_point(aes(size=TimeWindow), alpha=0.5) + 
  scale_color_manual(values="#4DAF4A") +
  scale_shape_manual(values=17) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(legend.position="none")

# trad1
ggplt.dtrad1.est.nf <- ggplot(newdf[newdf$FixedEffect=="city_code + p" & newdf$Variable=="num_firms",], 
                               aes(x=regwithin, y=d_trad1_cp, colour=Variable, shape=FixedEffect)) +
  #   theme(panel.grid.major = element_line(colour="gray70", size=0.5, linetype = "dashed"),
  #         panel.grid.minor = element_line(colour="gray80", size=0.5, linetype = "dotted")) +
  scale_x_continuous(name="Regression specification") + 
  #scale_y_continuous(name="Std. Coef. of density type 3", limits=c(-0.7, 0.7), breaks=seq(-0.7,0.7,0.1)) + 
  #scale_y_continuous(name="Coef. of density type 3") + 
  scale_y_continuous(name="Std. Coef. of density type 3", limits=c(-1, 2), breaks=seq(-1,2,0.2)) + 
  geom_hline(yintercept = 0, alpha=0.5, size=1.5, color="black") + 
  geom_errorbar(aes(ymin = lower_d_trad1_cp, ymax=upper_d_trad1_cp)) + 
  geom_point(aes(size=TimeWindow), alpha=0.5) + 
  scale_color_manual(values="#4DAF4A") +
  scale_shape_manual(values=17) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(legend.position="none")

# trad2
ggplt.dtrad2.est.nf <- ggplot(newdf[newdf$FixedEffect=="city_code + p" & newdf$Variable=="num_firms",], 
                               aes(x=regwithin, y=d_trad2_cp, colour=Variable, shape=FixedEffect)) +
  #   theme(panel.grid.major = element_line(colour="gray70", size=0.5, linetype = "dashed"),
  #         panel.grid.minor = element_line(colour="gray80", size=0.5, linetype = "dotted")) +
  scale_x_continuous(name="Regression specification") + 
  #scale_y_continuous(name="Std. Coef. of density type 4", limits=c(-0.7, 0.7), breaks=seq(-0.7,0.7,0.1)) + 
  #scale_y_continuous(name="Coef. of density type 4") + 
  scale_y_continuous(name="Std. Coef. of density type 4", limits=c(-1, 2), breaks=seq(-1,2,0.2)) + 
  geom_hline(yintercept = 0, alpha=0.5, size=1.5, color="black") + 
  geom_errorbar(aes(ymin = lower_d_trad2_cp, ymax=upper_d_trad2_cp)) + 
  geom_point(aes(size=TimeWindow), alpha=0.5) + 
  scale_color_manual(values="#4DAF4A") +
  scale_shape_manual(values=17) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(legend.position="none")
# -------------------------

# -------------------------


gga.modRCA <- ggarrange(ggplt.dworker.est.modRCA, ggplt.dcity.est.modRCA, ggplt.dtrad1.est.modRCA, ggplt.dtrad2.est.modRCA, 
                        #common.legend = TRUE,
                        #labels = c("B", "C", "D", "E"),
                        labels = c("C", "D", "E", "F"),
                        ncol = 4, nrow = 1,
                        widths = c(1,1,1,1), heights=1)

gga.emp <- ggarrange(ggplt.dworker.est.emp, ggplt.dcity.est.emp, ggplt.dtrad1.est.emp, ggplt.dtrad2.est.emp, 
                        #common.legend = TRUE,
                        #labels = c("B", "C", "D", "E"),
                        labels = c("G", "H", "I", "J"),
                        ncol = 4, nrow = 1,
                        widths = c(1,1,1,1), heights=1)

gga.nf <- ggarrange(ggplt.dworker.est.nf, ggplt.dcity.est.nf, ggplt.dtrad1.est.nf, ggplt.dtrad2.est.nf, 
                     #common.legend = TRUE,
                     #labels = c("B", "C", "D", "E"),
                     labels = c("K", "L", "M", "N"),
                     ncol = 4, nrow = 1,
                     widths = c(1,1,1,1), heights=1)

gga.all <- ggarrange(ggplt.R2,
                 ggplt.mrtm,
                 gga.modRCA, 
                 gga.emp,
                 gga.nf,
                 common.legend = TRUE,
                 labels = c("A", "B"),
                 ncol = 1, nrow = 5,
                 widths = 3, heights=c(1,1,1,1,1))
gga.all

# ---------------------------------------------
namefile <- "/IndExpFinal_3000_DensityRegressionsPlots_disaggr.png"
png(paste0(figpath,namefile), 
    height=18, 
    width=15, 
    units="in", 
    res=300)
# ---------------------------------------------

print(gga.all)
# ---------------------------------------------
dev.off()
# ---------------------------------------------




# ################################################
# Returning the table of the regression that makes most sense
rowi <- which(allformulas_df$FixedEffect=="city_code + p" & 
        allformulas_df$Variable=="num_firms" &
        allformulas_df$TimeWindow==5 & 
          allformulas_df$Densities=="alldens")
allformulas_df[rowi,"Formula"]

myfit_0 <- felm(numfirms_change_5yr  ~  num_firms | city_code + p, data=cpy_df)

#myfit_i <- felm(as.formula(as.character(allformulas_df[i,"Formula"])), data=cpy_df)
myfit_1 <- felm(numfirms_change_5yr  ~  num_firms + d_city_cp | city_code + p, data=cpy_df)
myfit_1lm <- lm(numfirms_change_5yr  ~  num_firms + d_city_cp + city_code + p, data=cpy_df)
head(summary(myfit_1lm)$coefficients)
summary(myfit_1lm)$coefficients[1,"Estimate"]

myfit_2 <- felm(numfirms_change_5yr  ~  num_firms + d_trad1_cp | city_code + p, data=cpy_df)
myfit_2lm <- lm(numfirms_change_5yr  ~  num_firms + d_trad1_cp + city_code + p, data=cpy_df)

myfit_3 <- felm(numfirms_change_5yr  ~  num_firms + d_city_cp + d_trad1_cp | city_code + p, data=cpy_df)
myfit_3lm <- lm(numfirms_change_5yr  ~  num_firms + d_city_cp + d_trad1_cp + city_code + p, data=cpy_df)
myfit_3p <- felm(numfirms_change_5yr  ~  num_firms + I(d_city_cp + d_worker_cp) + d_trad1_cp | city_code + p, data=cpy_df)

myfit_4 <- felm(numfirms_change_5yr  ~  num_firms + d_worker_cp + d_trad1_cp | city_code + p, data=cpy_df)
myfit_5 <- felm(numfirms_change_5yr  ~  num_firms + d_city_cp + d_trad1_cp, data=cpy_df)

myfit_6 <- felm(numfirms_change_1yr  ~  num_firms + d_city_cp + d_trad1_cp | city_code + p, data=cpy_df)

summary(myfit_3)
summary(myfit_4)
summary(myfit_5)
summary(myfit_6)

stargazer(myfit_0, myfit_1, myfit_2, myfit_3,
          dep.var.labels = "Change in number of firms in 5 years",
          covariate.labels = c("Number of firms initial year",
                               "$D^{(2)}$",
                               "$D^{(3)}$"),
          label = "tab:finaldensityresultsnumfirms",
          font.size = "scriptsize",
          out = paste0(tablespath,"/IndExpFinal_3000_DefinitiveDensityRegression_num_firms.tex"))



# ################################################
# Returning the table of the regression that makes most sense

myfit_0 <- felm(modrca_change_5yr  ~  pslogrca_mod | city_code + p, data=cpy_df)

#myfit_i <- felm(as.formula(as.character(allformulas_df[i,"Formula"])), data=cpy_df)
myfit_1 <- felm(modrca_change_5yr  ~  pslogrca_mod + d_city_cp | city_code + p, data=cpy_df)
myfit_1lm <- lm(modrca_change_5yr  ~  pslogrca_mod + d_city_cp + city_code + p, data=cpy_df)
head(summary(myfit_1lm)$coefficients)
summary(myfit_1lm)$coefficients[1,"Estimate"]

myfit_2 <- felm(modrca_change_5yr  ~  pslogrca_mod + d_trad1_cp | city_code + p, data=cpy_df)
myfit_2lm <- lm(modrca_change_5yr  ~  pslogrca_mod + d_trad1_cp + city_code + p, data=cpy_df)

myfit_3 <- felm(modrca_change_5yr  ~  pslogrca_mod + d_city_cp + d_trad1_cp | city_code + p, data=cpy_df)
myfit_3lm <- lm(modrca_change_5yr  ~  pslogrca_mod + d_city_cp + d_trad1_cp + city_code + p, data=cpy_df)
myfit_3p <- felm(modrca_change_5yr  ~  pslogrca_mod + I(d_city_cp+d_worker_cp) + d_trad1_cp | city_code + p, data=cpy_df)

myfit_4 <- felm(modrca_change_5yr  ~  pslogrca_mod + d_worker_cp + d_trad1_cp | city_code + p, data=cpy_df)
myfit_5 <- felm(modrca_change_5yr  ~  pslogrca_mod + d_city_cp + d_trad1_cp, data=cpy_df)

myfit_6 <- felm(modrca_change_5yr  ~  pslogrca_mod + d_city_cp + d_trad1_cp | city_code + p, data=cpy_df)

summary(myfit_3)
summary(myfit_4)
summary(myfit_5)
summary(myfit_6)

stargazer(myfit_0, myfit_1, myfit_2, myfit_3,
          dep.var.labels = "Change in modRCA in 5 years",
          covariate.labels = c("modRCA initial year",
                               "$D^{(2)}$",
                               "$D^{(3)}$"),
          label = "tab:finaldensityresultsmodrca",
          font.size = "scriptsize",
          out = paste0(tablespath,"/IndExpFinal_3000_DefinitiveDensityRegression_modrca.tex"))


# ################################################
# Returning the table of the regression that makes most sense

myfit_0 <- felm(emp_change_5yr  ~  employment | city_code + p, data=cpy_df)

#myfit_i <- felm(as.formula(as.character(allformulas_df[i,"Formula"])), data=cpy_df)
myfit_1 <- felm(emp_change_5yr  ~  employment + d_city_cp | city_code + p, data=cpy_df)
myfit_1lm <- lm(emp_change_5yr  ~  employment + d_city_cp + city_code + p, data=cpy_df)
head(summary(myfit_1lm)$coefficients)
summary(myfit_1lm)$coefficients[1,"Estimate"]

myfit_2 <- felm(emp_change_5yr  ~  employment + d_trad1_cp | city_code + p, data=cpy_df)
myfit_2lm <- lm(emp_change_5yr  ~  employment + d_trad1_cp + city_code + p, data=cpy_df)

myfit_3 <- felm(emp_change_5yr  ~  employment + d_city_cp + d_trad1_cp | city_code + p, data=cpy_df)
myfit_3lm <- lm(emp_change_5yr  ~  employment + d_city_cp + d_trad1_cp + city_code + p, data=cpy_df)

myfit_4 <- felm(emp_change_5yr  ~  employment + d_worker_cp + d_trad1_cp | city_code + p, data=cpy_df)
myfit_5 <- felm(emp_change_5yr  ~  employment + d_city_cp + d_trad1_cp, data=cpy_df)

myfit_6 <- felm(emp_change_5yr  ~  employment + d_city_cp + d_trad1_cp | city_code + p, data=cpy_df)

summary(myfit_3)
summary(myfit_4)
summary(myfit_5)
summary(myfit_6)

stargazer(myfit_0, myfit_1, myfit_2, myfit_3,
          dep.var.labels = "Change in employment in 5 years",
          covariate.labels = c("Employment initial year",
                               "$D^{(2)}$",
                               "$D^{(3)}$"),
          label = "tab:finaldensityresultsemp",
          font.size = "scriptsize",
          out = paste0(tablespath,"/IndExpFinal_3000_DefinitiveDensityRegression_emp.tex"))




#################################
# Making sure the regression to the mean has to expected behavior
# library(plotly)
# library(dplyr)

ggplot(newdf, aes(x=Variable, y=tstat_ReversionToMean)) + 
  geom_bar(stat="identity")


  


# p-values
ggplot(newdf, aes(x=regwithin, y=pval_d_worker_cp, colour=Variable, shape=FixedEffect, size=d_worker_cp)) +
  geom_point() + 
  scale_y_continuous(limits=c(0,0.001))

ggplot(newdf, aes(x=regwithin, y=pval_d_worker_cp, colour=Variable, shape=FixedEffect, size=d_worker_cp)) +
  geom_point() + 
  #scale_y_log10(limits=c(10^(-50),100)) +
  scale_y_continuous(limits=c(0,1)) +
  geom_hline(yintercept = 0.001, col="red") +
  annotate("text", 18, 0.002, label="p-value = 0.001", col="red") +
  geom_hline(yintercept = 0.01, col="red") +
  annotate("text", 18, 0.02, label="p-value = 0.01", col="red") +
  geom_hline(yintercept = 0.1, col="red") +
  annotate("text", 18, 0.2, label="p-value = 0.1", col="red")

# t-statistics
ggplot(newdf, aes(x=regwithin, y=tstat_d_worker_cp, colour=Variable, shape=FixedEffect, size=d_worker_cp)) +
  geom_point() + 
  scale_y_log10(limits=c(-10,100))

ggplot(newdf, aes(x=regwithin, y=tstat_d_worker_cp, colour=Variable, shape=FixedEffect, size=d_worker_cp, fill=Densities)) +
  geom_point() + 
  scale_y_log10(limits=c(-10,100))






# # ========================================================
# # GGPLOTS
# # --------------------------------------------------------
# numbydepvar <- nrow(newdf[newdf$DependentVar=="modrca_change_1yr",])
# uniquedepvar <- as.character(unique(newdf$DependentVar))
# newdf$regrank_d_worker_cp <- 0
# newdf$regrank_d_city_cp <- 0
# newdf$regrank_d_trad1_cp <- 0
# newdf$regrank_d_trad2_cp <- 0
# for(dvar in uniquedepvar){
#   for(dens in densities[2:5]){
#     newdf[newdf$DependentVar==dvar, paste0("regrank_", dens)]
#   }
# }
# 
# asdf <- c(1,2,3,9,8,7)
# order(asdf)
# 
# 
# 
# 
# # # =============================================================
# # values of density
# pvalsdensityvec <- newdf$pval_d_worker_cp
# summary(pvalsdensityvec)
# hist(pvalsdensityvec, 30)
# 
# var2plot <- "employment"
# tw2plot <- 1
# newdf <- newdf[order(newdf$d_worker_cp),]
# # tempdf <- newdf[newdf$Variable==var2plot & newdf$TimeWindow==tw2plot & !is.na(newdf$d_worker_cp),
# #                 c("reg", "d_worker_cp", "se_d_worker_cp")]
# tempdf <- newdf[newdf$Variable==var2plot & !is.na(newdf$d_worker_cp),
#                 c("reg", "d_worker_cp", "se_d_worker_cp")]
# 
# 
# 
# ######################################################
# mycolorfunc <- colorRampPalette(c("#9e0142",
#                                   "#d7191c",
#                                   "#fdae61",
#                                   "#abdda4",
#                                   "#2b83ba",
#                                   "#5e4fa2"))
# 
# groups <- as.character(c(1:nrow(tempdf)))
# bpColor <- mycolorfunc(length(groups))
# 
# 
# 
# # ----------------
# # For plotting
# namefile <- "/IndExpFinal_3000_DensityRegressions_test.png"
# png(paste(figpath,namefile,sep=""), 
#     height=6, 
#     width=7, 
#     units="in", 
#     res=600)
# # ----------------
# 
# mult.fig(1, marP=c(1,1,0,0))
# xlbl <- "Sample of workers from total"
# ylbl <- "Estimate and Std.Error"
# plot(NULL, type="n",
#      xlab = "",
#      ylab = "", 
#      xaxt="n",
#      yaxt="n",
#      xlim = c(0.5, nrow(tempdf) + 0.5),
#      ylim = c(-0.01, 1.12))
# 
# abline(h=0, lwd=3, col="black")
# abline(h=mean(tempdf$d_worker_cp), lty="dotted", lwd=2, col="gray")
# 
# # simulation estimates
# arrows(c(1:nrow(tempdf)), tempdf$d_worker_cp-tempdf$se_d_worker_cp, 
#        c(1:nrow(tempdf)), tempdf$d_worker_cp+tempdf$se_d_worker_cp, 
#        length=0.05, angle=90, code=3, lwd=4, col=alpha("#a63603",0.7))
# points(c(1:nrow(tempdf)), tempdf$d_worker_cp, 
#        bg="#7f2704", pch=21, col="gray70", cex=2.2)
# title(xlab=xlbl, line=2.8, cex.lab=1.9)
# title(ylab=ylbl, line=2.5, cex.lab=1.9)
# axis(1, at=c(1:nrow(tempdf)), labels=groups)
# axis(2, las=2)
# 
# 
# # real estimates
# arrows(c(1:nrow(tempdf)), realbetasdf$betas-realbetasdf$se, 
#        c(1:nrow(tempdf)), realbetasdf$betas+realbetasdf$se, 
#        length=0.05, angle=90, code=3, lwd=4, col=alpha("#08519c",0.7))
# points(realbetasdf$betas ~ realbetasdf$SampleSize, 
#        bg="#08306b", pch=21, col="gray70", cex=2.2)
# 
# 
# #Add an annotation within the graph using legend()
# # this adds to an existing plot
# legend("topleft", #where to put it; could be "top..." or "bottom...", "...left" or "....right"
#        legend = c(" ", " ",
#                   " ", " "), #what to write
#        lwd=c(2.1, -1,
#              2.1, -1),
#        col=c(alpha("#08519c",0.7),NA,
#              alpha("#a63603",0.7),NA),
#        inset = c(0.001,0.015),
#        x.intersp = 2,
#        xjust = 0.5,
#        bg="transparent",
#        seg.len = 1.2,
#        cex = 1.4,
#        bty = "n") #remove the box around the legend
# legend("topleft", #where to put it; could be "top..." or "bottom...", "...left" or "....right"
#        legend = c(" Real estimate", " with error bars = +/- 1 SE",
#                   " Mean estimate after spatial randomization", " with error bars = +/- 1 SD"), #what to write
#        pch=c(21,-1,
#              21,-1),
#        pt.bg=c("#08306b",-1,
#                "#7f2704",-1),
#        pt.cex=c(1.5,-1,
#                 1.5,-1),
#        col=c("gray70", NULL,
#              "gray70", NULL),
#        inset = c(0.025,0.015),
#        #xjust = 1,
#        bg="transparent",
#        cex = 1.4,
#        bty = "n") #remove the box around the legend
# 
# # ----------------
# dev.off()
# # ----------------
# 
# 
# 
# 
