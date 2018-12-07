rm(list = ls(all = TRUE))  # resets R to fresh
gc()

pryr::mem_used()

library(plyr)
library(dplyr)
library(tidyr)
library(data.table)
library(sfsmisc)


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

M1 = matrix(c(1,1,1,1,1,1,1,1,1,1,
              1,1,1,1,1,1,1,1,1,0,
              1,1,1,1,1,1,1,1,0,0,
              1,1,1,1,1,1,1,0,0,0,
              1,1,1,1,1,1,0,0,0,0,
              1,1,1,1,1,0,0,0,0,0,
              1,1,1,1,0,0,0,0,0,0,
              1,1,1,0,0,0,0,0,0,0,
              1,1,0,0,0,0,0,0,0,0,
              1,0,0,0,0,0,0,0,0,0), 
            byrow = TRUE,
            nrow = 10)
M2 = matrix(c(1,0,0,0,0,1,1,1,1,1,
              0,0,0,0,0,1,1,1,1,0,
              0,0,0,0,0,1,1,1,0,0,
              0,0,0,0,0,1,1,0,0,0,
              0,0,0,0,0,1,0,0,0,0,
              1,1,1,1,1,1,0,0,0,0,
              1,1,1,1,0,0,0,0,0,0,
              1,1,1,0,0,0,0,0,0,0,
              1,1,0,0,0,0,0,0,0,0,
              1,0,0,0,0,0,0,0,0,0), 
            byrow = TRUE,
            nrow = 10)

M3 = matrix(c(0,0,0,0,0,0,1,1,1,1,
              0,0,0,0,0,1,1,1,1,0,
              0,0,0,0,1,1,1,1,0,0,
              0,0,0,1,1,1,1,0,0,0,
              0,0,1,1,1,1,0,0,0,0,
              0,1,1,1,1,0,0,0,0,0,
              1,1,1,1,0,0,0,0,0,0,
              1,1,1,0,0,0,0,0,0,1,
              1,1,0,0,0,0,0,0,1,1,
              1,0,0,0,0,0,0,1,1,1), 
            byrow = TRUE,
            nrow = 10)
M4 = matrix(c(1,1,1,1,1,1,1,1,1,1,
              1,1,1,1,1,1,1,1,1,0,
              1,1,1,1,1,1,1,1,0,0,
              1,1,1,1,1,1,1,0,0,0,
              1,1,1,1,1,1,0,0,0,0,
              1,1,1,1,1,0,0,0,0,0,
              1,1,1,1,0,1,1,1,1,0,
              1,1,1,0,0,1,1,1,0,0,
              1,1,0,0,0,1,1,0,0,0,
              1,0,0,0,0,1,0,0,0,0), 
            byrow = TRUE,
            nrow = 10)

M6 <- unitemats(M1, M2)

# # Adding artificial countries with single product
# M1 <- rbind(M1, diag(ncol(M1)))
# M2 <- rbind(M2, diag(ncol(M2)))
# M3 <- rbind(M3, diag(ncol(M3)))


row.names(M1) <- paste0("C", c(1:nrow(M1)))
colnames(M1) <- paste0("P", c(1:ncol(M1)))
row.names(M2) <- paste0("C", c(1:nrow(M2)))
colnames(M2) <- paste0("P", c(1:ncol(M2)))
row.names(M3) <- paste0("C", c(1:nrow(M3)))
colnames(M3) <- paste0("P", c(1:ncol(M3)))
row.names(M4) <- paste0("C", c(1:nrow(M4)))
colnames(M4) <- paste0("P", c(1:ncol(M4)))

row.names(M6) <- paste0("C", c(1:nrow(M6)))
colnames(M6) <- paste0("P", c(1:ncol(M6)))

res <- list()
res[[1]] <- eigspaces(M1)
res[[2]] <- eigspaces(M2)
res[[3]] <- eigspaces(M3)
res[[4]] <- eigspaces(M4)
res[[5]] <- eigspaces(M6)

round(t(res[[1]]$Gl)%*%res[[1]]$Gr, 2) # This is supposed to be a diagonal matrix

(lambdas <- round(res[[1]]$values, 2))
round(t(res[[1]]$Gl)%*%res[[1]]$C%*%res[[1]]$Gr, 2) # this is supposed to be diag(lambdas)
round(solve(t(res[[1]]$Gl))%*%diag(lambdas)%*%solve(res[[1]]$Gr), 2) # this is supposed to get C back
round(res[[1]]$Gr%*%diag(lambdas)%*%t(res[[1]]$Gl), 2) # this is supposed to get C back
res[[1]]$C

mycolors <- c("#e41a1c",
              "#377eb8",
              "#4daf4a",
              "#984ea3")
letcex <- 1

mult.fig(6, mfrow = c(2,3))
plot(1,1, 
     xlim = c(-1, 1),
     ylim = c(-1, 1),
     type="n")
grid()
points(0,0, pch=21, bg="black")
for(i in c(1:length(res))){
  text(res[[i]]$Gl[,2], res[[i]]$Gl[,3], labels = row.names(res[[i]]$C), cex=letcex, col=mycolors[i])
}

plot(1,1, 
     xlim = c(-1, 1),
     ylim = c(-1, 1),
     type="n")
grid()
points(0,0, pch=21, bg="black")
for(i in c(1:length(res))){
  text(res[[i]]$Gr[,2], res[[i]]$Gr[,3], labels = row.names(res[[i]]$C), cex=letcex, col=mycolors[i])
}

plot(1,1, 
     xlim = c(-1, 1),
     ylim = c(-1, 1),
     type="n")
grid()
points(0,0, pch=21, bg="black")
for(i in c(1:length(res))){
  text(res[[i]]$Gl[,2], res[[i]]$Gr[,2], labels = row.names(res[[i]]$C), cex=letcex, col=mycolors[i])
}



plot(1,1, 
     xlim = c(-1, 1),
     ylim = c(-1, 1),
     type="n")
grid()
points(0,0, pch=21, bg="black")
for(i in c(1:length(res))){
  text(res[[i]]$Gl[,2], res[[i]]$Gl[,4], labels = row.names(res[[i]]$C), cex=letcex, col=mycolors[i])
}
legend("topleft", legend=c("Perfect triang.", "Noisy col.", "Spec.", "Two triangles"),
       col = mycolors, lty=1, seg.len = 1, cex=0.6)

plot(1,1, 
     xlim = c(-1, 1),
     ylim = c(-1, 1),
     type="n")
grid()
points(0,0, pch=21, bg="black")
for(i in c(1:length(res))){
  text(res[[i]]$Gr[,2], res[[i]]$Gr[,4], labels = row.names(res[[i]]$C), cex=letcex, col=mycolors[i])
}

plot(1,1, 
     xlim = c(-1, 1),
     ylim = c(-1, 1),
     type="n")
grid()
points(0,0, pch=21, bg="black")
for(i in c(1:length(res))){
  text(res[[i]]$Gl[,3], res[[i]]$Gr[,3], labels = row.names(res[[i]]$C), cex=letcex, col=mycolors[i])
}




##############################################

# -----------------
# RECONSTRUCTING THE TRAJECTORY OF A COUNTRY
numCs <- 100
numAs <- 1000
numPs <- 2000
frac <- 0.2
Ccamat1 <- matrix(0, nrow = numCs, ncol=numAs)
Ccamat2 <- matrix(0, nrow = numCs*frac, ncol=numAs*frac)
Papmat1 <- matrix(0, nrow = numAs, ncol=numPs)
Papmat2 <- matrix(0, nrow = numAs*frac, ncol=numPs*frac)
rCs <- 10^(-seq(0.01, 0.9, length.out = numCs))
tAs <- 10^(-seq(0.01, 1.5, length.out = numAs))
qPs <- 10^(-seq(0.1, log10(numAs/2), length.out = numPs))
for(ai in c(1:numAs)){
  for(ci in c(1:numCs)){
    Ccamat1[ci, ai] <- min(floor(runif(1) + rCs[ci]*tAs[ai]^0.5), 1)
    if(ci>numCs*(1-frac) & ai>numAs*(1-frac)){
      Ccamat2[ci-numCs*(1-frac), ai-numAs*(1-frac)] <- min(floor(runif(1) + rCs[ci]^2*tAs[ai]), 1)
    }
    #Ccamat[ci, ai] <- rCs[ci] + tAs[ai]
  }
  for(pi in c(1:numPs)){
    Papmat1[ai, pi] <- min(floor(runif(1) + qPs[pi]*tAs[ai]), 1)
    if(pi>numPs*(1-frac) & ai>numAs*(1-frac)){
      Papmat2[ai-numAs*(1-frac), pi-numPs*(1-frac)] <- min(floor(runif(1) + qPs[pi]*tAs[ai]), 1)
    }
    #Papmat[ai, pi] <- qPs[pi] + tAs[ai]
  }
}
row.names(Ccamat1) <- paste0("Ca", c(1:nrow(Ccamat1)))
colnames(Ccamat1) <- paste0("Aa", c(1:ncol(Ccamat1)))
row.names(Ccamat2) <- paste0("Cb", c(1:nrow(Ccamat2)))
colnames(Ccamat2) <- paste0("Ab", c(1:ncol(Ccamat2)))
Ccamat <- unitemats(Ccamat1, Ccamat2)

row.names(Papmat1) <- paste0("Aa", c(1:nrow(Papmat1)))
colnames(Papmat1) <- paste0("Pa", c(1:ncol(Papmat1)))
row.names(Papmat2) <- paste0("Ab", c(1:nrow(Papmat2)))
colnames(Papmat2) <- paste0("Pb", c(1:ncol(Papmat2)))
Papmat <- unitemats(Papmat1, Papmat2)



ctyxty <- rowSums(Ccamat)
prodxty <- colSums(Papmat)
prodxty[prodxty==0] <- 1
Jinv <- diag(1/prodxty)
row.names(Jinv) <- colnames(Papmat)
colnames(Jinv) <- colnames(Papmat)

Mcp <- floor(Ccamat%*%Papmat%*%Jinv)
dim(Mcp)
rs <- rowSums(Mcp)
cs <- colSums(Mcp)
Mcpfinal <- Mcp[rs!=0, cs!=0]
# row.names(Mcpfinal) <- paste0("C", c(1:nrow(Mcpfinal)))
# colnames(Mcpfinal) <- paste0("P", c(1:ncol(Mcpfinal)))

dim(Mcpfinal)
(rs <- rowSums(Mcpfinal))
(cs <- colSums(Mcpfinal))
crealcxtys <- ctyxty[row.names(Mcpfinal)]
prealcxtys <- prodxty[colnames(Mcpfinal)]

# ----------------------------
resfinal <- eigspaces(Mcpfinal)
# ----------------------------




library(pheatmap)
pheatmap(Mcpfinal, cluster_rows = FALSE, cluster_cols = FALSE)
pheatmap(Mcpfinal[order(-rs), order(-cs)], cluster_rows = FALSE, cluster_cols = FALSE)



ctytrajectorymat <- matrix(0, nrow=20, ncol=ncol(Mcpfinal))
for(rr in c(1:nrow(ctytrajectorymat))){
  ctytrajectorymat[rr,] <- c(rep(1,2*rr), rep(0,ncol(Mcpfinal)-2*rr))
}
row.names(ctytrajectorymat) <- paste0("Year", c(1:nrow(ctytrajectorymat)))
colnames(ctytrajectorymat) <- paste0("P", c(1:ncol(ctytrajectorymat)))
# embed in eigenspaces
rs <- rowSums(ctytrajectorymat)
rs[rs==0] <- 1
Dinv <- diag(1/rs)
row.names(Dinv) <- row.names(ctytrajectorymat)
colnames(Dinv) <- row.names(ctytrajectorymat)

Rnew <- Dinv%*%ctytrajectorymat
Glnew <- Rnew%*%resfinal$Ql
Glnew <- Glnew[,c(1:nrow(ctytrajectorymat))]
Grnew <- solve(t(Glnew))



mult.fig(4, mfrow = c(2,2))
plot(resfinal$Gl[,2], resfinal$Gl[,3], 
     # xlim = c(-0.5,0.5),
     # ylim = c(-0.5, 0.5),
     type="n")
grid()
points(0,0, pch=21, bg="black")
text(resfinal$Gl[,2], resfinal$Gl[,3], labels = row.names(resfinal$C), cex=0.8, col=mycolors[4])
for(i in c(1:nrow(Mcpfinal))){
  points(resfinal$Gl[i,2], resfinal$Gl[i,3], cex=(2*crealcxtys[i]/max(crealcxtys))^2)
}
text(Glnew[,2]/sqrt(Glnew[,2]%*%Glnew[,2]), Glnew[,3]/sqrt(Glnew[,3]%*%Glnew[,3]), labels = row.names(ctytrajectorymat), cex=0.8, col=mycolors[1])


plot(resfinal$Gr[,2], resfinal$Gr[,3], 
     # xlim = c(-0.8,0.8),
     # ylim = c(-0.8, 0.8),
     type="n")
grid()
points(0,0, pch=21, bg="black")
text(resfinal$Gr[,2], resfinal$Gr[,3], labels = row.names(resfinal$C), cex=0.8, col=mycolors[4])
for(i in c(1:nrow(Mcpfinal))){
  points(resfinal$Gr[i,2], resfinal$Gr[i,3], cex=(2*crealcxtys[i]/max(crealcxtys))^2)
}
text(Grnew[,2]/sqrt(Grnew[,2]%*%Grnew[,2]), Grnew[,3]/sqrt(Grnew[,3]%*%Grnew[,3]), labels = row.names(ctytrajectorymat), cex=0.8, col=mycolors[1])


plot(resfinal$Gl[,2], resfinal$Gr[,2], 
     # xlim = c(-0.4,0.4),
     # ylim = c(-0.2, 0.6),
     type="n")
grid()
abline(a=0, b=1)
points(0,0, pch=21, bg="black")
#text(resfinal$Gl[,2], resfinal$Gr[,2], labels = row.names(resfinal$C), cex=0.8, col=mycolors[4])
for(i in c(1:nrow(Mcpfinal))){
  points(resfinal$Gl[i,2], resfinal$Gr[i,2], cex=(2*crealcxtys[i]/max(crealcxtys))^2)
}

plot(resfinal$Gl[,3], resfinal$Gr[,3], 
     # xlim = c(-0.4,0.4),
     # ylim = c(-0.2, 0.6),
     type="n")
grid()
abline(a=0, b=1)
points(0,0, pch=21, bg="black")
#text(resfinal$Gl[,2], resfinal$Gr[,2], labels = row.names(resfinal$C), cex=0.8, col=mycolors[4])
for(i in c(1:nrow(Mcpfinal))){
  points(resfinal$Gl[i,3], resfinal$Gr[i,3], cex=(2*crealcxtys[i]/max(crealcxtys))^2)
}

complexities1 <- diag(resfinal$Gr%*%t(resfinal$Gr)) - diag(resfinal$Gl%*%t(resfinal$Gl))
complexities2 <- sqrt(diag(resfinal$Gr%*%diag(resfinal$values)%*%t(resfinal$Gr))) - sqrt(diag(resfinal$Gl%*%diag(resfinal$values)%*%t(resfinal$Gl)))
complexities3 <- sqrt(diag(resfinal$Gr%*%diag(resfinal$values)%*%t(resfinal$Gr)))
complexities4 <- sqrt(diag(resfinal$Gr%*%t(resfinal$Gr)))
complexities5 <- diag(resfinal$Gl%*%diag(1-resfinal$values)%*%t(resfinal$Gl))
complexities6 <- diag(resfinal$Gr%*%diag(resfinal$values)%*%t(resfinal$Gl))
complexities7 <- diag(resfinal$Gr%*%t(resfinal$Gl))


# new data point, completely diversified
Rnew <- rep(1,ncol(Mcpfinal))/ncol(Mcpfinal)
Glnew <- Rnew%*%resfinal$Ql
Grnew <- 1/Glnew[,c(1:nrow(Mcpfinal))]
ptxy <- sqrt(Grnew%*%Grnew)





################################################

library(NMF)

M <- Mcpfinal
kk = nrow(Mcpfinal) -1
nmfout <- nmf(M, kk)
nmfH <- nmfout@fit@W
nmfV <- t(nmfout@fit@H)

pcaout <- svd(M)
pcaH <- (pcaout$u %*% diag(pcaout$d))[,c(1:kk)]
pcaV <- (pcaout$v)[,c(1:kk)]

pca2out <- svd(diag(1/rowSums(Mcpfinal))%*%M%*%diag(1/colSums(Mcpfinal))*sum(Mcpfinal))
pca2H <- (pca2out$u %*% diag(pca2out$d))[,c(1:kk)]
pca2V <- (pca2out$v)[,c(1:kk)]

cout <- eigen(M%*%diag(1/colSums(Mcpfinal))%*%t(M))
pout <- eigen(t(M)%*%diag(1/rowSums(Mcpfinal))%*%M)
dout <- eigen(diag(1/rowSums(Mcpfinal))%*%M%*%t(M)%*%diag(1/rowSums(Mcpfinal)))



cH <- cout$vectors
pH <- pout$vectors
dH <- dout$vectors

mult.fig(3, mfrow=c(1,3))
plot(rowSums(nmfH), crealcxtys)
plot(sqrt(diag(pcaH%*%t(pcaH))), rowSums(Mcpfinal))
plot(sqrt(diag(pcaout$u[,c(1:kk)]%*%t(pcaout$u[,c(1:kk)]))), crealcxtys)

complexities8 <- sqrt(diag(nmfH%*%t(nmfH)))
complexities9 <- sqrt(diag(pcaH%*%t(pcaH)))
complexities10 <- sqrt(diag(pca2H%*%t(pca2H)))
round(cor(data.frame(complexities1, 
               complexities2,
               complexities3,
               complexities4,
               complexities5,
               complexities6,
               complexities7,
               complexities8,
               complexities9,
               complexities10,
               rowSums(Mcpfinal),
               crealcxtys)),4)

#vec <- rowSums(Mcpfinal)
vec <- complexities10
mult.fig(1)
plot(vec, crealcxtys,
     # xlim = c(-0.4,0.4),
     # ylim = c(-0.2, 0.6),
     type="n")
grid()
points(0,0, pch=21, bg="black")
#text(resfinal$Gl[,2], resfinal$Gr[,2], labels = row.names(resfinal$C), cex=0.8, col=mycolors[4])
for(i in c(1:nrow(Mcpfinal))){
  if(substr(names(crealcxtys)[i], 1, 2)=="Ca"){
    points(vec[i], crealcxtys[i], cex=(2*crealcxtys[i]/max(crealcxtys))^2, col="green")
  } else {
    points(vec[i], crealcxtys[i], cex=(2*crealcxtys[i]/max(crealcxtys))^2, col="blue")
  }
}
inds1 <- which(substr(names(crealcxtys), 1, 2)=="Ca")
inds2 <- which(substr(names(crealcxtys), 1, 2)!="Ca")
abline(lm(crealcxtys[inds1] ~ vec[inds1]), col="green")
abline(lm(crealcxtys[inds2] ~ vec[inds2]), col="blue")




(nmfround4M <- round(nmfH %*% t(nmfV), 2))
(pcaround4M <- round(pcaH %*% t(pcaV), 2))
(nmfround0M <- round(nmfH %*% t(nmfV), 0))
(pcaround0M <- round(pcaH %*% t(pcaV), 0))
(floorM <- floor(W %*% H))
plot(as.numeric(approxM), as.numeric(M))








M <- Mcpfinal
Pcp = M
dc = rowSums(M)
up = colSums(M)
#m <- nrow(M)*ncol(M)/(nrow(M)+ncol(M))
m <- 2*sum(M)
#m <- sum(M)
for(cc in c(1:nrow(M))){
  for(pp in c(1:ncol(M))){
    Pcp[cc,pp] <- dc[cc]*up[pp]/m
  }
}

Pcplong <- as.data.frame(Pcp) %>%
  mutate(cty = row.names(Pcp)) %>%
  gather(prod, pcp, -cty) %>%
  as.data.frame()

Pcplong$mlnp <- -log(-log(Pcplong$pcp))
Pcplong$cty <- as.factor(Pcplong$cty)
Pcplong$prod <- as.factor(Pcplong$prod)

reslm <- lm(mlnp ~ cty + prod - 1, data=Pcplong)
coefnames2 <- row.names(summary(reslm)$coef)
msanames <- coefnames2[substr(coefnames2, 1, 3)=="cty"]
actnames <- coefnames2[substr(coefnames2, 1, 4)=="prod"]

fe_msa <- as.data.table(summary(reslm)$coef[msanames, ])
fe_msa$msacode <- as.factor(substr(msanames, 4, nchar(msanames)))










library(mdatools)


mat2im <- function(M){
  dM <- as.data.frame(melt(M))
  dM <- as.matrix(dM[,"value"], ncol=1)
  attr(dM, "width") = ncol(M)
  attr(dM, "height") = nrow(M)
  attr(dM, "name") <- "Image"
  #dM[,"value"] <- as.numeric(dM[,"value"])

  return(dM)  
}

M5 = matrix(c(1,1,1,1,1,1,1,1,1,0,
              1,1,1,1,1,1,1,0,1,0,
              1,1,1,1,0,1,1,1,0,0,
              1,1,1,1,1,1,0,0,0,0,
              1,1,1,1,1,1,0,0,0,0,
              1,1,1,0,1,0,0,0,0,0,
              1,0,1,1,0,0,0,0,0,1,
              1,1,1,0,0,0,0,0,0,0,
              1,1,0,0,0,0,0,0,0,0,
              1,0,0,0,0,1,0,0,0,1), 
            byrow = TRUE,
            nrow = 10)
row.names(M5) <- paste0("C", c(1:nrow(M5)))
colnames(M5) <- paste0("P", c(1:ncol(M5)))


#M <- Mcpfinal[order(-rowSums(Mcpfinal)), order(-colSums(Mcpfinal))]
M <- M5
res <- eigspaces(M)
#res$P[res$P<median(as.numeric(res$P))] <- 0
#res$C[res$C<median(as.numeric(res$C))] <- 0
# Mnp <- floor(matrix(runif(M), nrow=nrow(M)) + round(M%*%res$P,4)^4)
# Mnc <- floor(matrix(runif(M), nrow=nrow(M)) + round(t(res$C)%*%M,4)^4)
# Mncp <- floor(matrix(runif(M), nrow=nrow(M)) + round(t(res$C)%*%M%*%res$P,4))
powerexp <- 1
Mnp <- (M%*%res$P)^powerexp
Mnc <- (t(res$C)%*%M)^powerexp

powerexp <- 2.0
Mn2pc <- ((Mnp + Mnc)/2)^powerexp
Mncp <- (t(res$C)%*%M%*%res$P)^powerexp

#persistence <- 0.9
persistence <- (rowSums(M)%o%colSums(M)/(max(rowSums(M))*max(colSums(M))))^0.1
#persistence <- rowSums(M)%o%colSums(M)/sum(M)
Mnew <- floor(matrix(runif(M), nrow=nrow(M)) + 0.001 + (1-M)*Mn2pc + persistence*M)
(Mn2pc*0.2 + persistence*0.8)
(Mn2pc*(1-M) + persistence*M)


mult.fig(6, mfrow = c(2,3))
imshow(mat2im(M), 1, 
       main = "Mcp at time t",
       colmap = colorRampPalette(c('gray90', 'gray10'))(100))
imshow(mat2im(Mnp), 1, 
       main = "P-based density",
       colmap = colorRampPalette(c('gray90', 'gray10'))(100))
imshow(mat2im(Mnc), 1, 
       main = "C-based density",
       colmap = colorRampPalette(c('gray90', 'gray10'))(100))

imshow(mat2im(persistence), 1, 
       main = "\"persistence\"",
       colmap = colorRampPalette(c('gray90', 'gray10'))(100))
imshow(mat2im(Mn2pc), 1, 
       main = "(C+P)/2 density",
       colmap = colorRampPalette(c('gray90', 'gray10'))(100))
imshow(mat2im(Mnew), 1, 
       main = "t+1",
       colmap = colorRampPalette(c('gray90', 'gray10'))(100))

sum(Mnew)/sum(M)

evolvemat <- function(M0, steps, powerexp=2.0){
  
  M <- M0
  for(ii in c(1:steps)){
    res <- eigspaces(M)
    
    Mnp <- (M%*%res$P)
    Mnc <- (t(res$C)%*%M)
    
    Mn2pc <- (0.51*Mnp + 0.51*Mnc)^powerexp
    
    #persistence <- rowSums(M)%o%colSums(M)/sum(M)
    #persistence <- 0.9
    #M <- floor(matrix(runif(M), nrow=nrow(M)) + (1-persistence)*Mn2pc + persistence*M)

        
    persistence <- (rowSums(M)%o%colSums(M)/(max(rowSums(M))*max(colSums(M))))^0.1
    #persistence <- rowSums(M)%o%colSums(M)/sum(M)
    M <- floor(matrix(runif(M), nrow=nrow(M)) + 0.0001 + (1-M)*Mn2pc + persistence*M)
    
    
  }
  
  return(M)
}

mult.fig(6, mfrow = c(2,3))
imshow(mat2im(M), 1, 
       main = "Current",
       colmap = colorRampPalette(c('gray90', 'gray10'))(100))
M2 <- evolvemat(M, 3)
imshow(mat2im(M2), 1, 
       main = "Step 1",
       colmap = colorRampPalette(c('gray90', 'gray10'))(100))
M3 <- evolvemat(M2, 3)
imshow(mat2im(M3), 1, 
       main = "Step 3",
       colmap = colorRampPalette(c('gray90', 'gray10'))(100))

M4 <- evolvemat(M3, 3)
imshow(mat2im(M4), 1, 
       main = "Step 5",
       colmap = colorRampPalette(c('gray90', 'gray10'))(100))
M5 <- evolvemat(M4, 3)
imshow(mat2im(M5), 1, 
       main = "Step 7",
       colmap = colorRampPalette(c('gray90', 'gray10'))(100))
M6 <- evolvemat(M5, 3)
imshow(mat2im(M6), 1, 
       main = "Step 9",
       colmap = colorRampPalette(c('gray90', 'gray10'))(100))


