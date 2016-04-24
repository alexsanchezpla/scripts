require(Hmisc)
require(energy)
require(FactoMineR)

vecCorrs <- function (x, y){
  corrS <- rcorr(x, y, type="spearman")
  corrP <- rcorr(x, y, type="pearson")
  return(unlist(list(rhoS=corrS$r[1,2],rhoP=corrP$r[1,2],
                     pvalS=corrS$P[1,2], pvalP=corrP$P[1,2])))
}

matCorrs <- function (X, Y){
  if ((nrow(X)!=nrow(Y))||(ncol(X)!=ncol(Y))) stop('matrices dimensions do not match')
  corrsList<- matrix(NA, nrow=nrow(X), ncol=4)
  colnames(corrsList) <- c("r (Sp)", "r (Pear)",  "p (Sp)","p (Pear)")
  for (i in 1:nrow(X)){
    corrs<- vecCorrs(X[i,],Y[i,])
    corrsList[i,] <- corrs
  }
  rownames(corrsList) <- rownames(X)
  return(corrsList)
}

distCorrs <- function (x, y){
  distCorr <- dcor(x, y, index=1.0)
  return(unlist(list(distCorr)))
}

matDistCorr <- function (X, Y){
  if ((nrow(X)!=nrow(Y))||(ncol(X)!=ncol(Y))) stop('matrices dimensions do not match')
  distCorrsList<- matrix(NA, nrow=nrow(X), ncol=1)
  colnames(distCorrsList) <- c("DistCor")
  for (i in 1:nrow(X)){
    DistCorr<- distCorrs(X[i,],Y[i,])
    distCorrsList[i,] <- DistCorr
  }
  rownames(distCorrsList) <- rownames(X)
  return(distCorrsList)
}

allCorrs  <- function (x, y){
  corrS <- rcorr(x, y, type="spearman")
  corrP <- rcorr(x, y, type="pearson")
  distCorr <- dcor(x, y, index=1.0)
  return(unlist(list(rhoS=corrS$r[1,2],  rhoP=corrP$r[1,2], distCorr,
                     pvalS=corrS$P[1,2], pvalP=corrP$P[1,2])))
}

sort1 <- function (X, col,DEC=TRUE, ...){
  return(X[sort.list(X[,col], decreasing=DEC), ])
}

matAllCorrs  <- function (X, Y){
  if ((nrow(X)!=nrow(Y))||(ncol(X)!=ncol(Y))) stop('matrices dimensions do not match')
  corrsList<- matrix(NA, nrow=nrow(X), ncol=5)
  colnames(corrsList) <- c("r (Sp)", "r (Pear)", "distCor", "p (Sp)",  "p (Pear)")
  for (i in 1:nrow(X)){
    corrs<- allCorrs(X[i,],Y[i,])
    corrsList[i,] <- corrs
  }
  corrsList<- cbind(corrsList, adj.Spear.Pval= p.adjust(corrsList[,"p (Sp)"],"fdr"))
  corrsList<- cbind(corrsList, adj.Pear.Pval = p.adjust(corrsList[,"p (Pear)"],"fdr"))
  rownames(corrsList) <- rownames(X)
  corrsList <- sort1(corrsList,4, DEC=FALSE)
  return(corrsList)
}

multivCorr <- function(X,Y){
  # stopifnot(require(energy))
  # stopifnot(require(FactoMineR))
  # RV1 <-coeffRV(X, Y)
  # cat("p-value : ", RV1$p.value,"\n")
  # dcor1<-dcor(X, Y)
  # cat("DistCorr: ", dcor1,"\n")
  coin1 <- cia(X, Y)
  cat("RV coeff: ", coin1$coinertia$RV,"\n")
}


# test
#
# (X <- round(matrix (rnorm(30)*10, ncol=6),1))
# (Y <- round(X + matrix (rnorm(30)*10, ncol=6),1))
# (m1<-matCorrs(X,Y))
# (m2<-matDistCorr(X,Y))
# (m12<- matAllCorrs (X, Y))
# sort1(m12,1)
# sort1(m12,3)
# dcor(X,Y,1)
# coeffRV(X,Y)
# multivCorr(X,Y)
## Adding an NA to X will make the preceeding fail
#  X[1,1] <-NA
# (m1<-matCorrs(X,Y))
# (m2<-matDistCorr(X,Y)) # yields an Error
# (m12<- matAllCorrs (X, Y)) # yields an Error
# sort1(m12,4, DEC=FALSE)
