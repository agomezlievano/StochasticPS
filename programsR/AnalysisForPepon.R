# Analysis for Pepon

library(data.table)
library(sfsmisc)
library(scales)
library(dplyr)
library(tidyr)

library(mdatools)
library(RColorBrewer)


mat2im <- function(M){
  dM <- as.data.frame(melt(M))
  dM <- as.matrix(dM[,"value"], ncol=1)
  attr(dM, "width") = ncol(M)
  attr(dM, "height") = nrow(M)
  attr(dM, "name") <- "Image"
  #dM[,"value"] <- as.numeric(dM[,"value"])
  
  return(dM)  
}



unitemats <- function(A, B, withrandomness=TRUE){
  nA <- nrow(A)
  mA <- ncol(A)
  nB <- nrow(B)
  mB <- ncol(B)
  
  if(withrandomness){
    belowA <- matrix(floor(0.5*(nB + mA)/(nB*mA) + runif(nB*mA)), nrow=nB, ncol=mA)
    aboveB <- matrix(floor(0.5*(nA + mB)/(nA*mB) + runif(nA*mB)), nrow=nA, ncol=mB)
  } else {
    belowA <- matrix(0, nrow=nB, ncol=mA)
    aboveB <- matrix(0, nrow=nA, ncol=mB)
  }
  row.names(belowA) <- row.names(B)
  colnames(belowA) <- colnames(A)
  row.names(aboveB) <- row.names(A)
  colnames(aboveB) <- colnames(B)
  
  return(cbind(rbind(A, belowA), rbind(aboveB, B)))
}

eigspaces <- function(M){
  cs <- colSums(M)
  cs[cs==0] <- 1
  Uinv <- diag(1/cs)
  row.names(Uinv) <- colnames(M)
  colnames(Uinv) <- colnames(M)
  
  rs <- rowSums(M)
  rs[rs==0] <- 1
  Dinv <- diag(1/rs)
  row.names(Dinv) <- row.names(M)
  colnames(Dinv) <- row.names(M)
  
  L <- as.matrix(M)%*%Uinv
  R <- Dinv%*%as.matrix(M)
  
  C <- L%*%t(R)
  P <- t(R)%*%L
  
  
  resC <- eigen(C)
  resP <- eigen(P)
  #Gr <- apply(resC$vectors, MARGIN=2, function(x)Re(x)*sign(cor(Re(x), rs)))
  #Gr <- apply(resC$vectors, MARGIN=2, function(x)Re(x)*sign(x[which.max(rs)]))
  Gr <- apply(resC$vectors, MARGIN=2, function(x)Re(x)*sign(Re(x[1])))
  Qr <- apply(resP$vectors, MARGIN=2, function(x)Re(x)*sign(Re(x[1])))
  
  resC <- eigen(t(C))
  resP <- eigen(t(P))
  #Gl <- apply(resC$vectors, MARGIN=2, function(x)Re(x)*sign(cor(Re(x), rs)))
  #Gl <- apply(resC$vectors, MARGIN=2, function(x)Re(x)*sign(x[which.max(rs)]))
  Gl <- apply(resC$vectors, MARGIN=2, function(x)Re(x)*sign(Re(x[1])))
  Ql <- apply(resP$vectors, MARGIN=2, function(x)Re(x)*sign(Re(x[1])))
  
  row.names(Gl) <- row.names(M)
  colnames(Gl) <- row.names(M)
  
  row.names(Gr) <- row.names(M)
  colnames(Gr) <- row.names(M)
  
  
  row.names(Ql) <- colnames(M)
  colnames(Ql) <- colnames(M)
  
  row.names(Qr) <- colnames(M)
  colnames(Qr) <- colnames(M)
  
  #rankofmatrix <- qr(M)$rank
  
  return(list(Gr=Gr, Gl=Gl, C=C, Qr=Qr, Ql=Ql, P=P, values=resC$values, R=R, L=L))
  
}

#df <- rio::import(file = "~/../Downloads/data/M1_Matrix_30x27.csv")
df <- rio::import(file = "~/../Downloads/data/M2_Matrix_30_x_36_reordered.csv")
df[1:5, 1:5]

row.names(df) <- df$Country
df <- df %>% 
  select(-Country) %>%
  as.data.frame()

Mcpfinal <- as.matrix(df)
resfinal <- eigspaces(Mcpfinal)

# ----
# comparing with random matrix
Mcprand <- randmat(Mcpfinal)
row.names(Mcprand) <- Mcprand$origin
Mcprand <- Mcprand %>% select(-origin) %>% as.matrix()
plot(rowSums(Mcprand[row.names(Mcpfinal),]), rowSums(Mcpfinal))
resrand <- eigspaces(Mcprand[row.names(Mcpfinal),])
plot(Re(resfinal$values) ~ Re(resrand$values),
     xlab = "Randomized",
     ylab = "Real")
abline(a=0, b=1)
abline(lm(log(Re(resfinal$values)) ~ log(Re(resrand$values))), col="red")
# ---


# num communities
mult.fig(1)
plot(resfinal$values, rep(1, length(resfinal$values)))
nkk <- round(sum(Re(resfinal$values)))

mycolors <- c("#e41a1c",
              "#377eb8",
              "#4daf4a",
              "#984ea3")
letcex <- 0.2

crealcxtys <- rep(0.2, nrow(Mcpfinal))

mult.fig(2, mfrow = c(1,2))
plot(resfinal$Gl[,2], resfinal$Gl[,3], 
     xlim = c(-0.4,0.4),
     ylim = c(-0.5, 0.35),
     xlab = "Economic Complexity Index",
     ylab = "3rd left-eigenvector",
     main = "Left Eigenspace",
     type="n")
grid()
points(0,0, pch=21, bg="black")
text(resfinal$Gl[,2], resfinal$Gl[,3], labels = row.names(resfinal$C), cex=0.8, col=mycolors[4])
for(i in c(1:nrow(Mcpfinal))){
  points(resfinal$Gl[i,2], resfinal$Gl[i,3], cex=(1*crealcxtys[i]/max(crealcxtys))^2)
}


plot(resfinal$Gr[,2], resfinal$Gr[,3], 
     xlim = c(-0.4,0.45),
     ylim = c(-0.6, 0.35),
     xlab = "2nd right-eigenvector",
     ylab = "3rd right-eigenvector",
     main = "Right Eigenspace",
     type="n")
grid()
points(0,0, pch=21, bg="black")
text(resfinal$Gr[,2], resfinal$Gr[,3], labels = row.names(resfinal$C), cex=0.8, col=mycolors[4])
for(i in c(1:nrow(Mcpfinal))){
  points(resfinal$Gr[i,2], resfinal$Gr[i,3], cex=(1*crealcxtys[i]/max(crealcxtys))^2)
}




complexities1 <- diag(resfinal$Gr%*%t(resfinal$Gr)) - diag(resfinal$Gl%*%t(resfinal$Gl))
complexities2 <- sqrt(diag(resfinal$Gr%*%diag(resfinal$values)%*%t(resfinal$Gr))) - sqrt(diag(resfinal$Gl%*%diag(resfinal$values)%*%t(resfinal$Gl)))
complexities3 <- sqrt(diag(resfinal$Gr%*%diag(resfinal$values)%*%t(resfinal$Gr)))
complexities4 <- sqrt(diag(resfinal$Gr[,c(1:nkk)]%*%t(resfinal$Gr[,c(1:nkk)])))
complexities4full <- sqrt(diag(resfinal$Gr%*%t(resfinal$Gr)))
complexities5 <- diag(resfinal$Gl%*%diag(1-resfinal$values)%*%t(resfinal$Gl))
complexities6 <- diag(resfinal$Gr%*%diag(resfinal$values)%*%t(resfinal$Gl))
complexities7 <- diag(resfinal$Gr%*%t(resfinal$Gl))


ctysfull <- as.matrix(sort(complexities4full, decreasing = TRUE), ncol=1)
colnames(ctysfull) <- "Andres' Collective Knowhow metric FULL"

ctys <- as.matrix(sort(complexities4, decreasing = TRUE), ncol=1)
colnames(ctys) <- "Andres' Collective Knowhow metric"

divs <- as.matrix(sort(rowSums(Mcpfinal), decreasing = TRUE), ncol=1)
colnames(divs) <- "Diversity"

eci <- as.matrix(sort(resfinal$Gl[,2], decreasing = TRUE), ncol=1)
colnames(eci) <- "eci"


cccdf <- ctysfull %>%
  merge(ctys, by=0) %>%
  merge(divs, by.x="Row.names", by.y=0) %>%
  merge(eci, by.x="Row.names", by.y=0) %>%
  mutate(Country = Row.names) %>%
  select(-Row.names) %>%
  as.data.frame()

mult.fig(1)
plot(cccdf[,colnames(cccdf)!="Row.names"])




























# -----------------------------------------------------------
#install.packages("NbClust",dependencies = TRUE)
library(NbClust)
cccdf <- as.data.table(cccdf)

nb <- NbClust(resfinal$Gl[,c(2:nkk)], diss=NULL, distance = "euclidean", 
              min.nc=2, max.nc=8, method = "kmeans", 
              index = "all", alphaBeale = 0.1)
cccdf[, nbnumclusters := length(unique(nb$Best.partition))]

cnames <- names(nb$Best.partition)
clusters_dt <- as.data.table(nb$Best.partition)
clusters_dt[, Country:=cnames]

setkey(clusters_dt, Country)
setkey(cccdf, Country)
cccdf[clusters_dt, on='Country', nbCluster := i.V1]





# ===========================================================
# FIGURES
# -----------------------------------------------------------
namefile <- paste0("Pepon_matrix.png")
folder.figs <- "C:/Users/agomez/Dropbox/Harvard/LittleProjects/StochasticPS/figures/"
png(paste0(folder.figs,namefile), 
    height=6, width=10, units="in", res=200)

mult.fig(1)
imshow(mat2im(Mcpfinal[order(-rowSums(Mcpfinal)), order(-colSums(Mcpfinal))]), 1, 
       main = "Triangularity of the matrix",
       colmap = colorRampPalette(c('gray', 'gray10'))(100))

dev.off()

# =================================================================
# PLOTTING
# mycolors <- c("#e41a1c",
#               "#377eb8",
#               "#4daf4a",
#               "#984ea3")
mycolors <- brewer.pal(n = max(cccdf$nbCluster), name = "Dark2")
mycolors <- mycolors[c(1:length(mycolors))]
names(mycolors) <- as.character(sort(unique(cccdf$nbCluster)))

letcex <- 1
alf <- 0.9

#------
namefile <- paste0("Pepon_figs.png")
png(paste0(folder.figs,namefile), 
    height=4, width=14, units="in", res=600)

mult.fig(3, mfrow = c(1,3))
plot(density(Re(resfinal$values), bw = 0.05), col="gray30", lwd=3, 
     xlab = "Eigenvalues",
     ylab = "Density",
     main = "No. communities = No. large eigenvalues")
rug(jitter(Re(resfinal$values)), lwd=4, col="dodgerblue")

plot(resfinal$Gl[,2], resfinal$Gl[,3], 
     # xlim = c(-0.4,0.4),
     # ylim = c(-0.5, 0.35),
     xlab = "Economic Complexity Index",
     ylab = "3rd left-eigenvector",
     main = "Left Eigenspace",
     type="n")
grid()
legend("topleft", title = "Communities",
       cex = 1.0,
       legend = names(mycolors),
       pch = 21, 
       col = mycolors,
       pt.bg = mycolors)
points(0,0, pch=21, bg="black")
text(resfinal$Gl[,2], resfinal$Gl[,3], 
     labels = row.names(resfinal$C), 
     cex=1, 
     col=alpha("gray50", alf))
for(cl in unique(clusters_dt$V1)){
  points(resfinal$Gl[clusters_dt[V1==cl]$Country,2], 
         resfinal$Gl[clusters_dt[V1==cl]$Country,3],
         pch = 21,
         bg = alpha(mycolors[cl], alf),
         col = mycolors[cl],
         cex=letcex)
}


plot(resfinal$Gr[,2], resfinal$Gr[,3], 
     # xlim = c(-0.4,0.45),
     # ylim = c(-0.6, 0.35),
     xlab = "2nd right-eigenvector",
     ylab = "3rd right-eigenvector",
     main = "Right Eigenspace",
     type="n")
grid()
points(0,0, pch=21, bg="black")
text(resfinal$Gr[,2], resfinal$Gr[,3], 
     labels = row.names(resfinal$C), 
     cex=1, 
     col=alpha("gray50", alf))
for(cl in unique(clusters_dt$V1)){
  points(resfinal$Gr[clusters_dt[V1==cl]$Country,2], 
         resfinal$Gr[clusters_dt[V1==cl]$Country,3],
         pch = 21,
         bg = alpha(mycolors[cl], alf),
         col = mycolors[cl],
         cex=letcex)
}

dev.off()



head(cccdf)
colnames(cccdf) <- c("FULL", "ACKM", "DIVERSITY", "ECI", "Country", "NBNUMCLUSTERS", "NBCLUSTER")
rio::export(as.data.frame(cccdf[order(-ACKM)]), paste0(folder.figs, "PeponData.csv"))
