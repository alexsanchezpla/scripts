# missing
#http://www.stefvanbuuren.nl/mi/Software.html

path<-"/media/dades/Dropbox/upc/TFM/treball/miRNA/3/data"
pA<-file.path(path, "FASE 1_miRNA FFPE_A.csv")
pB<-file.path(path, "FASE 1_miRNA FFPE_B.csv")

daA<-read.csv(pA,sep=";",header = TRUE, dec=",",
              na.strings="Undetermined",colClasses=c(rep("factor",7),rep("numeric",242)))
summary(daA)

daB<-read.csv(pB,sep=";",header = TRUE,dec=",",              ,
              na.strings="Undetermined",colClasses=c(rep("factor",7),rep("numeric",144)))
summary(daB)


daA[1:5,1:8]
daB[1:5,1:8]
table(daA[,7],daB[,7])

da<-merge(daA,daB,by=c("GROUP","GROUP.CODE","COMP.1","COMP.2","COMP.3","OPEN.ARRAY_DAY","PATIENT"))
dim(da)
dim(daA)+dim(daB)

boxplot(t(da[,1])~t(da[,-c(1:7)]))

names(t(da[,-c(1:7)]))
rownames(da)
da[,1]
?boxplot
tda<-t(da)
boxplot(t(da[,1])~t(da[,-c(1:7)]),data=da)
boxplot(da[,1]~as.matrix(da[,-c(1:7)]))
class(da[,1])
class(da[,-c(1:7)])
class(t(da[,-c(1:7)]))


summary(da[,-c(1:6)])

class(da[,-c(1:6)])
dat<-class(t(da[,-c(1:7)]))

kk<-as.data.frame(t(da[,-c(1:7)]))
boxplot(kk,las=1,names=1:36)



####
# geometric mean
?geomMean

apply(kk,2,function(x)prod(x,na.rm=TRUE)^(1/length(x)))
kk[,1]

geomMean(kk[,1])


gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

gm_mean(kk)

install.packages("psych")
library(psych)
?geometric.mean
skk<-scale(kk,scale=FALSE,center=TRUE)
mg<-geometric.mean(kk,na.rm=TRUE)

msd<-apply(kk,2,sd,na.rm=TRUE)

kkn<-(kk-mg)/msd
boxplot(kkn,las=1)


geomean <- function(x, na.rm = FALSE, trim = 0, ...)
{
  exp(mean(log(x, ...), na.rm = na.rm, trim = trim, ...))
}

geomean(kk,na.rm=TRUE)
M<-apply(kk,2,geomean,na.rm=TRUE)


geosd <- function(x, na.rm = FALSE, ...)
{
  exp(sd(log(x, ...), na.rm = na.rm, ...))
}

S<-apply(kk,2,geosd,na.rm=TRUE)
         

kkn<-(kk-M)/S
boxplot(kkn,las=1)         
         
         
kkc<-kk-M
S<-apply(kkc,2,geosd,na.rm=TRUE)
kkcn<-(kkc-M)/S
boxplot(kkcn,las=1)

kn<-scale(kk,scale=TRUE,center=TRUE)
boxplot(kn,las=1)

kkc<-kk-M
S<-apply(kkc,2,geosd,na.rm=TRUE)
kkcn<-(kkc-M)/sqrt(M)
boxplot(kkcn,las=1)


# number of missing per individual  
ind_mis<-apply(da,1,function(x){sum(is.na(x))})
names(ind_mis)<-1:36
barplot(ind_mis,las=1)
ind_mis_per<-ind_mis*100/dim(da)[2]
ind_mis_per

# number of missing depending on the variable
var_mis<-apply(da,2,function(x){sum(is.na(x))})
names(var_mis)<-names(da)
barplot(var_mis,las=1)
sum(var_mis>8)
barplot(var_mis[which(var_mis>8)],las=2)

var_mis_per<-var_mis*100/dim(da)[1]
barplot(var_mis_per[which(var_mis_per>25)],las=2,ylim=c(0,100))
barplot(var_mis_per[which(var_mis_per>24) && which(var_mis_per<25)],las=2,ylim=c(0,100))

# distribution of the Ct values
un<-unlist(da[,-c(1:7)])
plot(density(un,na.rm=TRUE))
hist(un,las=1,ylim=c(0,3000))
den<-density(un,na.rm=TRUE)
lines(den$x,den$y*10^4.5)


biocLite("HTqPCR")
library(HTqPCR)



tda<-t(da[,-c(1:7)])
head(tda)
write.table(tda,"sortida.txt",sep="\t")
?readCtData
del<-read.delim("sortida.txt")
#https://stat.ethz.ch/pipermail/bioconductor/2012-April/045182.html
pcr<-readCtData(del$File)
pq.norm <- normalizeCtData(raw.cat, norm = "quantile")

cat(paste(readLines("sortida.txt")))


raw<-readCtData("sortida.txt",path=path, header=FALSE,n.features =393,format="plain",
                column.info=list(feature))


raw <- readCtData('borig.csv', header=F, n.features = 768,
                  column.info=list(feature=6, position=1, Ct=11), path="./orig")




tda<-t(da[-27,-c(1:7)])
tda[which(is.na(tda))]<-"Undetermined"
path2<-"/media/dades/Dropbox/upc/TFM/treball/miRNA/3/HTqPCR"
for (i in 1:35) {
  #a <- data.frame(mymat[, i])
  #mytime <- format(Sys.time(), "%b_%d_%H_%M_%S_%Y")
  #myfile <- file.path(path, paste0(mytime, "_", i, ".txt"))
  a<-tda[,i]
  myfile <- file.path(path2, paste0("Sample_", i, ".txt"))
    write.table(a, file = myfile, sep = "\t", row.names = TRUE, col.names = FALSE,
              quote = FALSE, append = FALSE)
}


# conjunt de 'files'
names(da)
da$GROUP
arxius<-paste("Sample_",1:35,".txt",sep="")
# arxius<-paste(arxius,".txt",sep="")
nomdf<-data.frame(File=arxius,Treatment=da[-27,]$GROUP)
nomdf
write.table(nomdf,file.path(path2,"files.txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE,append=FALSE)



#exPath <- system.file("exData", package="HTqPCR")
exFiles <- read.delim(file.path(path2, "files.txt"))
raw <- readCtData(files=exFiles$File, path=path2,header=FALSE,
                  n.features=386,column.info=list(feature=1,Ct=2))

# Example of adding missing information (random data in this case)
experimentData(raw)<-MIAME(name="Arnau Serra-Cayuela & Alexandre Sánchez Pla",lab="VHIR")
featureClass(raw)  <- factor(rep(c("A", "B", "C"), each=384/3))
pData(raw)[,"rep"]   <- c(1,1,2,2,3,3)


colnames(exprs(raw))
pData(raw)$sampleName<-colnames(exprs(raw))
  
  c(rep("A",11),rep("B",11),rep("C",13))
 
<-c(rep("A",11),rep("B",11),rep("C",13))

raw_nomis <- qpcrImpute(raw,iterMax=1000000,tol=100)

normCt <- normalizeCtData(raw, norm = "geometric.mean")
normCt
slotNames(normCt)
slotNames(normCt)
normCt@phenoData@.__classVersion__
pData(normCt)

normCt@assayData
exprs(raw)
exprs(normCt)[1,]
exprs(raw)[1,]

getCt(raw) # same as exprs

par(mfrow=c(2,1))
#boxplot(tda,las=1)
boxplot(exprs(raw),las=1)
boxplot(exprs(normCt),las=1)


#library(HTqPCR)
library(nondetects)

###################################################
### code chunk number 3: nondetects.Rnw:98-106
###################################################
conds <- paste("Sample",1:36,sep="_")
resids <- matrix(nrow=nrow(normCt), ncol=ncol(normCt))
for(i in 1:nrow(normCt)){
  for(j in 1:ncol(normCt)){
    ind <- which(conds==conds[j])
    resids[i,j] <- exprs(normCt)[i,j]-mean(exprs(normCt)[i,ind])
  }
}


###################################################
### code chunk number 4: nondetects.Rnw:110-113
###################################################
iND <- which(featureCategory(normCt)=="Undetermined", arr.ind=TRUE)
iD <- which(featureCategory(normCt)!="Undetermined", arr.ind=TRUE)
boxes <- list("observed"=-resids[iD], "non-detect"=-resids[iND])


###################################################
### code chunk number 5: nondetects.Rnw:115-117
###################################################
boxplot(boxes, main="",ylim=c(-12,12),
        ylab=expression(paste("-",Delta,"Ct residuals",sep="")))


###################################################
### code chunk number 6: nondetects.Rnw:121-123
###################################################
oncogene2013 <- qpcrImpute(oncogene2013, 
                           groupVars=c("sampleType","treatment"))


###################################################
### code chunk number 7: nondetects.Rnw:128-130
###################################################
normCt <- normalizeCtData(oncogene2013, norm = "deltaCt", 
                          deltaCt.genes = "Becn1")


###################################################
### code chunk number 8: nondetects.Rnw:134-135
###################################################
normCt <- normCt[-which(featureNames(normCt)=="Becn1"),]


###################################################
### code chunk number 9: nondetects.Rnw:139-148
###################################################
conds <- paste(pData(normCt)$sampleType,
               pData(normCt)$treatment,sep=":")
resids <- matrix(nrow=nrow(normCt), ncol=ncol(normCt))
for(i in 1:nrow(normCt)){
  for(j in 1:ncol(normCt)){
    ind <- which(conds==conds[j])
    resids[i,j] <- exprs(normCt)[i,j]-mean(exprs(normCt)[i,ind])
  }
}


###################################################
### code chunk number 10: nondetects.Rnw:152-155
###################################################
iI <- which(featureCategory(normCt)=="Imputed", arr.ind=TRUE)
iD <- which(featureCategory(normCt)!="Imputed", arr.ind=TRUE)
boxes <- list("observed"=-resids[iD], "imputed"=-resids[iI])


###################################################
### code chunk number 11: nondetects.Rnw:157-159
###################################################
boxplot(boxes, main="",ylim=c(-12,12),
        ylab=expression(paste("-",Delta,"Ct residuals",sep="")))


###################################################
### code chunk number 12: nondetects.Rnw:175-177
###################################################
library(nondetects)
data(sagmb2011)


###################################################
### code chunk number 13: nondetects.Rnw:188-190
###################################################
library(nondetects)
data(nature2008)


###################################################
### code chunk number 14: nondetects.Rnw:199-200
###################################################
sessionInfo()
