
In order to avoid including genes with too many zeroes --which means suggest genes that are not expressed in most samples-- those with \emph{\textbf{more than 66\% of zeros}} will be removed from the data.


zeros <-function(x){which (x==0)}

discard <- function(x, where, howMany){
  return(length(zeros(x[where])) > howMany) 
}

discardByPercent <- function(x, percentage){ 
    maxZeros <- ceiling (length(x)*percentage)
    return (discard(x,1:length(x),maxZeros))
}

discardBy2Groups <- function(x, grp1, grp2, percent){
    d1 <- discardByPercent (x[grp1], percent)
    d2 <- discardByPercent (x[grp2], percent)
  return (d1 || d2)
}

discardBy3Groups <- function(x, grp1, grp2, grp3, percent){
    d1 <- discardByPercent (x[grp1], percent)
    d2 <- discardByPercent (x[grp2], percent)
    d3 <- discardByPercent (x[grp3], percent)
  return (d1 || d2 || d3)
}

# test1
exData <- matrix (0, nrow=10, ncol=10)
for (i in 1:10){
     for(j in 1:i)
         exData[i,j]<-1
}
exData
d1 <-apply(exData, 1, discardByPercent, 0.66); length(d1); sum(d1);show(exData1 <- exData[!d1,])
d1 <-apply(exData, 1, discardByPercent, 0.5); length(d1); sum(d1);show(exData1 <- exData[!d1,])
d1 <-apply(exData, 1, discardByPercent, 0.33); length(d1); sum(d1);show(exData1 <- exData[!d1,])
d1 <-apply(exData, 1, discardByPercent, 0.2); length(d1); sum(d1);show(exData1 <- exData[!d1,])

# test2
grup1<- 1:4; grup2 <-5:10
d <- apply (exData,1,discardBy2Groups,grup1, grup2, 0.5)
for (i in 1:nrow(exData)){
    d<-discardBy2Groups(exData[i,], grup1, grup2, 0.5)
    print (d)
}
d1 <-apply(exData, 1, ; length(d1); sum(d1);show(exData1 <- exData[!d1,])

#
# discardedA <- apply(dataRNAseq, 1, discardByPercent, 0.65)
# length(discardedA); sum(discardedA)
# dataRNAseqA <- dataRNAseq [!discardedA,]
# dim(dataRNAseqA)
# RNAseqSymbolsA <-rownames(dataRNAseqA)
#



discardedA <- apply(protData_A, 1, discardA)
#sum(discardedA)
descartaA <- rep(NA,nrow(protData_A))
for (i in 1:nrow(protData_A)) descartaA[i]<-discardA(protData_A[i,])
discardedA<- sum(descartaA)

discardB <- function(x){
  return (discard(x,1:7,3) || discard(x,8:14,3))
}
# discardedB <- apply(protData_B, 1, discardB) # PERQUE NO FUNCIONA? 
# sum(discardedB)
descartaB <- rep(NA,nrow(protData_B))
for (i in 1:nrow(protData_B)) descartaB[i]<-discardB(protData_B[i,])
discardedB<- sum(descartaB)

protDataA <- protData_A [!descartaA,]
protDataB <- protData_B [!descartaB,]


