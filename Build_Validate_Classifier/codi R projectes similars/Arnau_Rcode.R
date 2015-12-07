#######################################################
#                                                     #
# R script to analyze PCR data                        #
#                                                     #
# Author: Arnau Serra-Cayuela                         #
# Supervisor: Alexandre Sánchez Pla                   #
#                                                     #
#######################################################


####################################################################################################
## ---- wkdir ----
# set the working directory

setwd("/media/dades/Dropbox/upc/TFM/treball/miRNA/3/")



####################################################################################################
## ---- packages ----
source("./scripts/install.R") # functions to install packages (CRAN and Bioconductor)

# r-base CRAN packages
instifnot("rminer")
instifnot("e1071") # note that in Linux r-base-dev must be installed (or install r-cran-e1071)
# TAKE CARE!!! load this package before CMA, since it mask the function 'tune' 

instifnot("glmnet") # package for lasso
instifnot("randomForest") # package for Random Forest
instifnot("gplots") # for Venn diagrams

instifnot("knitr") # to create tables (ktable)
instifnot("xtable")
#instifnot("vcd") # to estimate kappa
#instifnot("pROC")

instifnot("missMDA") # for missings
instifnot("rrcov") # for ROBust PCA (detect outliers func: PcaHubert)



##
# bioconductor packages
#instifnotbio("multtest")
instifnotbio("impute") # for missings
instifnotbio("HTqPCR")
instifnotbio("limma") # if limma feature selection is gonna be used
instifnotbio("e1071")
instifnot("glmpath")
instifnotbio("CMA") # machine learning package




####################################################################################################
## ---- data ----

###
# loading data
#da<-read.csv("./data/cancer_miRNA.csv",sep=";",header=TRUE,dec=",")

path<-"/media/dades/Dropbox/upc/TFM/treball/miRNA/3/data" # where are the files
pA<-file.path(path, "FASE 1_miRNA FFPE_A.csv") # plate A
pB<-file.path(path, "FASE 1_miRNA FFPE_B.csv") # plate B

#load plate A
daA<-read.csv(pA,sep=";",header = TRUE, dec=",",
              na.strings="Undetermined",colClasses=c(rep("factor",7),rep("numeric",242)))


#load plate B
daB<-read.csv(pB,sep=";",header = TRUE,dec=",",              ,
              na.strings="Undetermined",colClasses=c(rep("factor",7),rep("numeric",144)))

# check all variables are well loaded
summary(daA) 
summary(daB)
table(daA[,7],daB[,7])

# Merge both datasets (plate A and B)
da<-merge(daA,daB,sort=FALSE,
          by=c("GROUP","GROUP.CODE","COMP.1","COMP.2","COMP.3","OPEN.ARRAY_DAY","PATIENT"))
#write.csv(da,paste(path,"/merged.csv",sep=""))

# check all is ok
dim(da)
dim(daA)+dim(daB)


####################################################################################################
## ---- missPattern----

# are missings?
nmis<-sum(is.na(da))
# load functions to plot missing pattern
source("./scripts/missings.R")


####################################################################################################
## ---- missPatternSample ----
numMiss(da,1,per=TRUE)


####################################################################################################
## ---- missPatternVariables ----
va<-numMiss(da,2,per=TRUE,plot=FALSE)
missVar(va)


####################################################################################################
## ---- data_no_27 ----
dac<-da[-27,] # data clean 



####################################################################################################
## ---- CtDist ----
# distribution of the Ct values
un<-unlist(dac[,-c(1:7)])
#plot(density(un,na.rm=TRUE))
hist(un,las=1,ylim=c(0,3000),xlab="Ct value",main="",
     col="#f4dbdb",border="#d9baba")
den<-density(un,na.rm=TRUE)
lines(den$x,den$y*10^4.5,lwd=3)



####################################################################################################
## ---- missImputation----

#impute library
impKnn<-impute.knn(t(dac[,-c(1:7)]))

#missMDA library
nb<-estim_ncpPCA(t(dac[,-c(1:7)]),method="Regularized",method.cv="Kfold",nbsim=10)
impPCA<-imputePCA(t(dac[,-c(1:7)]),ncp=nb$ncp,method="Regularized",scale=TRUE)

# imputation by 40 ==> NA<-40 (num max cycles)
imp40<-t(dac[,-c(1:7)])
imp40[which(is.na(imp40)==TRUE)]<-40
sum(is.na(imp40)) # check

tda<-t(dac[,-c(1:7)]) # matrix of the data (change row for columns)


# check for correct missing assigment
# t(da)[which(is.na(t(dac)))]
# impKnn$data[which(is.na(t(dac)))]
# which(is.na(dac))
# kk[450]
# impKnn$data[1450]

# which(is.na(tda))
# i<-49
# tda[i];impKnn$data[i];impPCA$completeObs[i]




####################################################################################################
## ---- imputationDiagnostic----

marge<-par()$mar
par(mar=c(5.1,4.1,2.1,2.1))
un<-unlist(dac[,-c(1:7)])
plot(density(un,na.rm=TRUE),lwd=15,lty=1,xlim=c(0,45),col="grey",
     las=1,xlab="Ct values",main="")

unimpPCA<-unlist(impPCA$completeObs)
den<-density(unimpPCA)
lines(den$x,den$y,lwd=2,lty=1,col="blue")

unimpKnn<-unlist(impKnn$data)
den<-density(unimpKnn)
lines(den$x,den$y,lwd=2,lty=1,col="red")

unimp40<-unlist(lapply(imp40,c))
den<-density(unimp40)
lines(den$x,den$y,lwd=2,lty=1,col="black")

legend("topright",c("With missings","PCA imputation","Knn imputation","Ct40 imputation"),
       col=c("grey","blue","red","black"),lwd=2, lty=1, pch=NA)

par(mar=marge)

####################################################################################################
## ---- imputationDiagnostic2 ----

un<-unlist(dac[,-c(1:7)])
plot(density(un,na.rm=TRUE),lwd=6,lty=1,xlim=c(0,70),col="grey",
     las=1,xlab="Imputed Ct values",main="",ylim=c(0,0.12))

den<-density(impPCA$completeObs[which(is.na(tda))])
lines(den$x,den$y,col="blue",lwd=2)
den<-density(impKnn$data[which(is.na(tda))])
lines(den$x,den$y,col="red",lwd=2)
# Other: NA changed to 40
den<-density(unimp40[which(is.na(tda))])
lines(den$x,den$y,col="black",lwd=2)
legend("topright",c("With missings","PCA imputation","Knn imputation","Ct40 imputation"),
       col=c("grey","blue","red","black"),lwd=2, lty=1, pch=NA)






####################################################################################################
## ---- outliers ----
X.sc <- scale(t(impKnn$data))
X.HubPCA <- PcaHubert(X.sc)
summary(X.HubPCA)
X.HubPCA <- PcaHubert(X.sc,k=5,alpha=0.7)
plot(X.HubPCA,las=1)

#############################################################################################
####################### NORMALIZATION
#############################################################################################

####################################################################################################
## ---- normalization----

# we just use the matrix with the Ct values imputed by Knn 
#(PCA showed similar results, and imp40 provoques a lot of bias)

Kda<-impKnn$data # create new matrix, with complete observations
sum(is.na(Kda)) # check there is no missings

source("./scripts/geometricMeanSD.R")

# Geometric mean normalization
# (Xij-???j)/sdj

M<-apply(Kda,2,geomean) # mean of each variable
S<-apply(Kda,2,geosd) # standard deviation of each variable
KdaN<-sweep(Kda,2,M) # deduc the mean of each variable to the self variable values
KdaN<-sweep(KdaN,2,S,FUN="/") # divide previous values with the sd of each variable


####################################################################################################
## ---- normalizationBoxplot----
# comparison 
par(mfrow=c(2,1))
par(mar=c(2.1,4.1,2.1,2.1),cex.axis=.8)
boxplot(Kda,las=1,ylab="Ct values")
mtext("A",side=3,cex=.8)
par(mar=c(4.1,4.1,2.1,2.1))
boxplot(KdaN,las=1,xlab="Patient",ylab="Normalized Ct values")
mtext("B",side=3,cex=.8) 



####################################################################################################
## ---- HTqPCR----

Kda<-impKnn$data

# create the table to be read by HTqPCR
# create a file for each Patient
path2<-"/media/dades/Dropbox/upc/TFM/treball/miRNA/3/HTqPCR"
dir.create(path2, showWarnings = TRUE, recursive = FALSE) # create folder
for (i in 1:35) {
  a<-Kda[,i]
  myfile <- file.path(path2, paste0("Sample_", i, ".txt"))
  write.table(a, file = myfile, sep = "\t", row.names = TRUE, col.names = FALSE,
              quote = FALSE, append = FALSE)
}

# File that indicates which files are in the folder created previously
# conjunt de 'files'
arxius<-paste("Sample_",1:35,".txt",sep="")
# arxius<-paste(arxius,".txt",sep="")
nomdf<-data.frame(File=arxius,Treatment=dac$GROUP)
write.table(nomdf,file.path(path2,"files.txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE,append=FALSE)

# Read/load the files
exFiles <- read.delim(file.path(path2, "files.txt"))
raw <- readCtData(files=exFiles$File, path=path2,header=FALSE,
                  n.features=386,column.info=list(feature=1,Ct=2))

# Example of adding missing information (random data in this case)
experimentData(raw)<-MIAME(name="Arnau Serra-Cayuela & Alexandre Sánchez Pla",lab="VHIR")
# featureClass(raw)  <- factor(rep(c("A", "B", "C"), each=384/3))
# pData(raw)[,"rep"]   <- c(1,1,2,2,3,3)

# Normalization
normCt <- normalizeCtData(raw, norm = "geometric.mean")
normQt <- normalizeCtData(raw, norm = "quantile")
normRK<- normalizeCtData(raw, norm = "norm.rankinvariant")
# is the function wrong?
# https://stat.ethz.ch/pipermail/bioconductor/2014-May/059853.html
source("./scripts/normalizeCtData2.R")
normCt2 <- normalizeCtData2(raw, norm = "geometric.mean")


####################################################################################################
## ---- normalization2----

par(mfrow=c(2,2),mar=c(2,4,2,2))
boxplot(exprs(normCt),las=1,main="HTqPCR geometric mean")
boxplot(exprs(normQt),las=1,main="HTqPCR Quantile")
boxplot(exprs(normRK),las=1,main="HTqPCR Rank invariant")
boxplot(exprs(normCt2),las=1,main="HTqPCR geometric mean modified") 


####################################################################################################
## ---- normalization3----

par(mar=c(2.1,4,2,2))
plot(density(Kda),xlim=c(-10,35),ylim=c(0,0.12),lwd=15,col="gray80",
     main="",xlab="",las=1)
den<-density(KdaN) # geometric mean
lines(den$x,den$y,col="black",lwd=4)
den<-density(exprs(normQt))   # quantile HTqPCR
lines(den$x,den$y,col="red",lwd=4)
den<-density(exprs(normRK)) # rank invariant HTqPCR
lines(den$x,den$y,col="yellow",lwd=4)
den<-density(exprs(normCt2)) # geometric mean modified HTqPCR
lines(den$x,den$y,col="green",lwd=4)
den<-density(exprs(normCt)) # geometric mean HTqPRC
lines(den$x,den$y,col="blue",lwd=4)
legend("topright",c("Without norm.","Geom.","Geom. HTqPCR","Geom. HTqPCR 2","Quantile","Rank invariant"),
       col=c("gray80","black","blue","green","red","yellow"),pch=NA,lwd=2, lty=1,)




####################################################################################################
####################################################################################################
####################################################################################################


####################################################################################################
## ---- completeDataframe ----

dim(dac)
Nda<-dac #new data.frame
Nda[,-c(1:7)]<-t(KdaN) # with the values normalized and without outliers
summary(Nda)
write.csv(Nda,paste(path,"/gNormF1.csv",sep=""),row.names=FALSE)



#%%%%%%%%%%%%%%%%%%%%%% END of PREPROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
#  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  %%%%%%%%%%%%%%%%%%  COMPARISIONS
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  




####################################################################################################
## ---- comparison----

# load the previous data (to not run all the time the past code)
path<-"/media/dades/Dropbox/upc/TFM/treball/miRNA/3/data"
Nda<-read.csv(paste(path,"/gNormF1.csv",sep=""),header=TRUE,sep=",",dec=".",
              colClasses=c(rep("factor",7),rep("numeric",386)))

# Function to select the comparison

# select comparison
# 1 Cancer vs PSA cancer
# 2 Cancer vs no cancer
# 3

comparativa<-function(dataframe,numDummy,Dummy){
  del<-1:numDummy
  del<-del[-Dummy]
  dataframe<-dataframe[,-del]
  dataframe<-dataframe[(dataframe[,1]==0 | dataframe[,1]==1),]
  dataframe<-droplevels(dataframe)
  return(dataframe)
}

# Dummy
# 3: COMP.1
# 4: COMP.2
# 5: COMP.3
da2<-comparativa(Nda,7,3)




# ####################################################################################################
## ---- iterations -----
 
iterations <-2 # number of iterations



# ####################################################################################################
## ---- comp1varSel ----

# Set the seed for reproducible prouposes
set.seed(1984)
# create object learningsets
# specify the type of resampling, as well as its characteristics.
# resampling by 5-fold cross-validation with 1,000 iterations
CV5f<-GenerateLearningsets(y=da2[,1],method="CV",fold=4,niter=iterations,strat=TRUE)

# featrure selection methods
Tvar1<-system.time(
  var_t<-GeneSelection(X=as.matrix(da2[,-1]),y=da2[,1],learningsets=CV5f,method="t.test"))
Tvar2<-system.time(
  var_lim<-GeneSelection(X=as.matrix(da2[,-1]),y=da2[,1],learningsets=CV5f,method="limma"))
Tvar3<-system.time(
  var_wil<-GeneSelection(X=as.matrix(da2[,-1]),y=da2[,1],learningsets=CV5f,method="wilcox.test"))
Tvar4<-system.time(
  var_lasso<-GeneSelection(X=as.matrix(da2[,-1]),y=da2[,1],learningsets=CV5f,method="lasso"))
Tvar5<-system.time(
  var_rf<-GeneSelection(X=as.matrix(da2[,-1]),y=da2[,1],learningsets=CV5f,method="rf"))



# ####################################################################################################
## ---- variablesRank ----
# Functions to extract the top10(or other) variables, for each feature selection method, and also for
# the combination of all the methods. A Venn diagram can be ploted
source("./scripts/variablesRank.R")


# ####################################################################################################
## ---- comp1variablesPlot ----

metodes<-list(var_lim,var_t,var_lasso,var_wil,var_rf) # genesel object class createt previously
nomMetodes<-c("Limma","t-test","Lasso","Wilcox","Random Forest") # name of the methods
#setwd("/media/dades/Dropbox/upc/TFM/treball/miRNA/figures/")
vars<-commonTopVariables(da2,CV=CV5f,metodes,nomMetodes,nVar=10) # output


# ####################################################################################################
# ##################################       Classifiers    ###########################################
# ####################################################################################################



# ####################################################################################################
## ---- numVar ---
# number of variables to build the model
numVar<-5


# ####################################################################################################
## ---- comp1LDA ----

TC1_1<-system.time(
  comp1_t_lda<-classification(X=as.matrix(da2[,-1]),y=da2[,1],learningsets=CV5f,classifier=ldaCMA,genesel=var_t,nbgene=numVar))

TC1_2<-system.time(
  comp1_lim_lda<-classification(X=as.matrix(da2[,-1]),y=da2[,1],learningsets=CV5f,classifier=ldaCMA,genesel=var_lim,nbgene=numVar))

TC1_3<-system.time(
  comp1_wil_lda<-classification(X=as.matrix(da2[,-1]),y=da2[,1],learningsets=CV5f,classifier=ldaCMA,genesel=var_wil,nbgene=numVar))

TC1_4<-system.time(
  comp1_lasso_lda<-classification(X=as.matrix(da2[,-1]),y=da2[,1],learningsets=CV5f,classifier=ldaCMA,genesel=var_lasso,nbgene=numVar))

TC1_5<-system.time(
  comp1_rf_lda<-classification(X=as.matrix(da2[,-1]),y=da2[,1],learningsets=CV5f,classifier=ldaCMA,genesel=var_rf,nbgene=numVar))


# ####################################################################################################
## ---- comp1LR ----

T1class7<-system.time(
  comp1_t_lr<-classification(X=as.matrix(da2[,-1]),y=da2[,1],learningsets=CV5f,classifier=plrCMA,lambda=0,genesel=var_t,nbgene=numVar,
                             control=glm.control(epsilon = 1e-8, maxit = 1000)))

T1class8<-system.time(
  comp1_lim_lr<-classification(X=as.matrix(da2[,-1]),y=da2[,1],learningsets=CV5f,classifier=plrCMA,lambda=0,genesel=var_lim,nbgene=numVar,
                               control=glm.control(epsilon = 1e-8, maxit = 1000)))

T1class9<-system.time(
  comp1_wil_lr<-classification(X=as.matrix(da2[,-1]),y=da2[,1],learningsets=CV5f,classifier=plrCMA,lambda=0,genesel=var_wil,nbgene=numVar,
                               control=glm.control(epsilon = 1e-8, maxit = 1000)))

T1class10<-system.time(
  comp1_lasso_lr<-classification(X=as.matrix(da2[,-1]),y=da2[,1],learningsets=CV5f,classifier=plrCMA,lambda=0,genesel=var_lasso,nbgene=numVar,
                                 control=glm.control(epsilon = 1e-8, maxit = 1000)))

T1class11<-system.time(
  comp1_rf_lr<-classification(X=as.matrix(da2[,-1]),y=da2[,1],learningsets=CV5f,classifier=plrCMA,lambda=0,genesel=var_rf,nbgene=numVar,
                              control=glm.control(epsilon = 1e-8, maxit = 1000)))


# ####################################################################################################
## ---- comp1SVM ----


tune_t_svm <- tune(X=as.matrix(da2[,-1]),y=da2[,1], learningsets=CV5f,genesel=var_t, classifier = svmCMA,
                 grids = list(cost = c(0.1, 1, 5, 10, 50, 100, 500), degree = 2:4), kernel="polynomial",
                 strat=TRUE,nbgene=numVar)
T1class12<-system.time(
  comp1_t_svm <- classification(X=as.matrix(da2[,-1]),y=da2[,1], learningsets = CV5f, classifier = svmCMA, kernel = "polynomial",tuneres=tune_t_svm,probability=TRUE,nbgene=numVar,genesel=var_t))


tune_lim_svm <- tune(X=as.matrix(da2[,-1]),y=da2[,1], learningsets=CV5f,genesel=var_lim, classifier = svmCMA,
                 grids = list(cost = c(0.1, 1, 5, 10, 50, 100, 500), degree = 2:4), kernel="polynomial",
                 strat=TRUE,nbgene=numVar)
T1class13<-system.time(
  comp1_lim_svm <- classification(X=as.matrix(da2[,-1]),y=da2[,1], learningsets = CV5f, classifier = svmCMA, kernel = "polynomial",tuneres=tune_lim_svm,probability=TRUE,nbgene=numVar,genesel=var_lim))


tune_wil_svm <- tune(X=as.matrix(da2[,-1]),y=da2[,1], learningsets=CV5f,genesel=var_wil, classifier = svmCMA,
                 grids = list(cost = c(0.1, 1, 5, 10, 50, 100, 500), degree = 2:4), kernel="polynomial",
                 strat=TRUE,nbgene=numVar)
T1class14<-system.time(
  comp1_wil_svm <- classification(X=as.matrix(da2[,-1]),y=da2[,1], learningsets = CV5f, classifier = svmCMA, kernel = "polynomial",tuneres=tune_wil_svm,probability=TRUE,nbgene=numVar,genesel=var_wil))


tune_lasso_svm <- tune(X=as.matrix(da2[,-1]),y=da2[,1], learningsets=CV5f,genesel=var_lasso, classifier = svmCMA,
                 grids = list(cost = c(0.1, 1, 5, 10, 50, 100, 500), degree = 2:4), kernel="polynomial",
                 strat=TRUE,nbgene=numVar)
T1class15<-system.time(
  comp1_lasso_svm <- classification(X=as.matrix(da2[,-1]),y=da2[,1], learningsets = CV5f, classifier = svmCMA, kernel = "polynomial",tuneres=tune_lasso_svm,probability=TRUE,nbgene=numVar,genesel=var_lasso))


tune_rf_svm <- tune(X=as.matrix(da2[,-1]),y=da2[,1], learningsets=CV5f,genesel=var_rf, classifier = svmCMA,
                 grids = list(cost = c(0.1, 1, 5, 10, 50, 100, 500), degree = 2:4), kernel="polynomial",
                 strat=TRUE,nbgene=numVar)
T1class16<-system.time(
  comp1_rf_svm <- classification(X=as.matrix(da2[,-1]),y=da2[,1], learningsets = CV5f, classifier = svmCMA, kernel = "polynomial",tuneres=tune_rf_svm,probability=TRUE,nbgene=numVar,genesel=var_rf))



# ####################################################################################################
## ---- comp1RF ----

# tune
tune_t_rf<- tune(X=as.matrix(da2[,-1]),y=da2[,1],learningsets=CV5f,classifier=rfCMA,grids=list(mtry = ceiling(c(0.1, 0.25, 0.5)*sqrt(ncol(da2[,-1]))), nodesize = c(1,2,3)),
               strat=TRUE,genesel=var_t)
T1class17<-system.time(
  comp1_t_rf<-classification(X=as.matrix(da2[,-1]),y=da2[,1],learningsets=CV5f,classifier=rfCMA, genesel=var_t,tuneres=tune_t_rf))


tune_lim_rf<- tune(X=as.matrix(da2[,-1]),y=da2[,1],learningsets=CV5f,classifier=rfCMA,grids=list(mtry = ceiling(c(0.1, 0.25, 0.5)*sqrt(ncol(da2[,-1]))), nodesize = c(1,2,3)),
               strat=TRUE,genesel=var_lim)
T1class18<-system.time(
  comp1_lim_rf<-classification(X=as.matrix(da2[,-1]),y=da2[,1],learningsets=CV5f,classifier=rfCMA, genesel=var_lim,tuneres=tune_lim_rf))

tune_wil_rf<- tune(X=as.matrix(da2[,-1]),y=da2[,1],learningsets=CV5f,classifier=rfCMA,grids=list(mtry = ceiling(c(0.1, 0.25, 0.5)*sqrt(ncol(da2[,-1]))), nodesize = c(1,2,3)),
               strat=TRUE,genesel=var_wil)
T1class19<-system.time(
  comp1_wil_rf<-classification(X=as.matrix(da2[,-1]),y=da2[,1],learningsets=CV5f,classifier=rfCMA, genesel=var_wil,tuneres=tune_wil_rf))

tune_lasso_rf<- tune(X=as.matrix(da2[,-1]),y=da2[,1],learningsets=CV5f,classifier=rfCMA,grids=list(mtry = ceiling(c(0.1, 0.25, 0.5)*sqrt(ncol(da2[,-1]))), nodesize = c(1,2,3)),
               strat=TRUE,genesel=var_lasso)
T1class20<-system.time(
  comp1_lasso_rf<-classification(X=as.matrix(da2[,-1]),y=da2[,1],learningsets=CV5f,classifier=rfCMA, genesel=var_lasso,tuneres=tune_lasso_rf))


# ####################################################################################################
## ---- comp1Others ----


tune_pls <- tune(X=as.matrix(da2[,-1]),y=da2[,1],learningsets=CV5f,classifier=pls_ldaCMA,grids=list(comp=1:5),strat=TRUE)
T1class21<-system.time(
  comp1_plsDA<-classification(X=as.matrix(da2[,-1]),y=da2[,1],learningsets=CV5f,classifier=pls_ldaCMA,tuneres=tune_pls))

tune_svm <- tune(X=as.matrix(da2[,-1]),y=da2[,1],learningsets=CV5f,classifier=svmCMA,grids=list(cost = c(0.1, 1, 5, 10, 50, 100, 500), degree = 2:4),kernel="polynomial",strat=TRUE)
T1class22<-system.time(
  comp1_svm<-classification(X=as.matrix(da2[,-1]),y=da2[,1],learningsets=CV5f,classifier=svmCMA,tuneres=tune_svm,kernel="polynomial",probability=TRUE))

tune_rf<- tune(X=as.matrix(da2[,-1]),y=da2[,1],learningsets=CV5f,classifier=rfCMA,grids=list(mtry = ceiling(c(0.1, 0.25, 0.5)*sqrt(ncol(da2[,-1]))), nodesize = c(1,2,3)),strat=TRUE)
T1class23<-system.time(
  comp1_rf<-classification(X=as.matrix(da2[,-1]),y=da2[,1],learningsets=CV5f,classifier=rfCMA,tuneres=tune_rf))


tune_lasso <- tune(X=as.matrix(da2[,-1]),y=da2[,1],learningsets=CV5f,classifier=LassoCMA,grids=list(norm.fraction = seq(from=0.1, to=0.9, length=9)),strat=TRUE)
T1class24<-system.time(
  comp1_lasso<-classification(X=as.matrix(da2[,-1]),y=da2[,1],learningsets=CV5f,classifier=LassoCMA,tuneres=tune_lasso))

tune_elastic <- tune(X=as.matrix(da2[,-1]),y=da2[,1],learningsets=CV5f,classifier=ElasticNetCMA,strat=TRUE,
                     grids=list(norm.fraction = seq(from=0.1, to=0.9, length=5)))
T1class25<-system.time(
  comp1_elastic<-classification(X=as.matrix(da2[,-1]),y=da2[,1],learningsets=CV5f,classifier=ElasticNetCMA,tuneres=tune_elastic))

tune_plr <- tune(X=as.matrix(da2[,-1]),y=da2[,1],learningsets=CV5f,classifier=plrCMA,grids=list(lambda = 2^{-4:4}),strat=TRUE)
T1class26<-system.time(
  comp1_plr<-classification(X=as.matrix(da2[,-1]),y=da2[,1],learningsets=CV5f,classifier=plrCMA,tuneres=tune_plr))




# ####################################################################################################
# ##################################       Comparisons    ###########################################
# ####################################################################################################


# ####################################################################################################
## ---- nicePlot ----
source("./scripts/nicePlot.R")


# ####################################################################################################
## ---- comp1evalCMA ----

comp1ifierlist<-list(comp1_t_lda,comp1_lim_lda,comp1_wil_lda,comp1_lasso_lda,comp1_rf_lda,
                     comp1_t_lr,comp1_lim_lr,comp1_wil_lr,comp1_lasso_lr,comp1_rf_lr,
                     comp1_t_svm,comp1_lim_svm,comp1_wil_svm,comp1_lasso_svm,comp1_rf_svm,
                     comp1_t_rf,comp1_lim_rf,comp1_wil_rf,comp1_lasso_rf,
                     comp1_plsDA,comp1_svm,comp1_rf,
                     comp1_lasso,comp1_elastic,comp1_plr)

nomClass<-c("t-LDA","Limma-LDA","Wilcox-LDA","Lasso-LDA","RF-LDA",
            "t-LR","Limma-LR","Wilcox-LR","Lasso-LR","RF-LR",
            "t-SVM","Limma-SVM","Wilcox-SVM","Lasso-SVM","RF-SVM",
            "t-RF","Limma-RF","Wilcox-RF","Lasso-RF",
            "PLS-DA","SVM","RF","Lasso","ElasticNet","PLR")


# ####################################################################################################
## ---- comp1miss ----
nicePlot(comp1ifierlist,nomClass,"misclassification")

# ####################################################################################################
## ---- comp1sen ----
nicePlot(comp1ifierlist,nomClass,"sensitivity")

# ####################################################################################################
## ---- comp1spe ----
nicePlot(comp1ifierlist,nomClass,"specificity")

# ####################################################################################################
## ---- comp1auc ----
nicePlot(comp1ifierlist,nomClass,"auc")

#nicePlot(classifierlist,nomClass,"brier score")
#nicePlot(classifierlist,nomClass,"average probability")



####################################################################################################
#################################   Validation Set  ################################################
####################################################################################################










####################################################################################################
####################################################################################################
#######################################   OTHERS  #################################################
####################################################################################################
####################################################################################################




# ####################################################################################################
## ---- extract ----
# Functions to extract the top10(or other) variables, for each feature selection method, and also for
# the combination of all the methods. A Venn diagram can be ploted
source("./scripts/extractClass.R")
classifierlist<-list(comp1_t_lda,comp1_lim_lda,comp1_wil_lda,comp1_lasso_lda,comp1_rf_lda)
extracted<-extractList(classifierlist)

# ####################################################################################################
## ---- comp1MR ----


nomClass<-c("t-LDA","Limma-LDA","Wilcox-LDA","Lasso-LDA","RF-LDA")
va<-ev(extracted,measure="mis")
colorsa<-c(seq(c("blue")))
y<-rep(seq(1,5,1),4)
colo<-c("blue","orange","green","red","yellow")     
colors<-c(colo[y],rep("grey",5),"magenta")
#colors<-c(rep("blue",4),rep("orange",4),rep("green",4),rep("red",4),rep("yellow",4),rep("grey",5),"magenta")
marge<-par()$mar
par(mar=c(4.1,4.1,5.1,8.1))
boxplot(va,las=1,names=nomClass,ylab=expression(bold("misclassification rate")),
        ylim=c(0,1),col=colors,xaxt="n")

# boxplot(va,las=1,ylab=expression(bold("misclassification rate")),
#         ylim=c(0,1),xaxt="n")

axis(1, at=seq(1,length(nomClass), by=1), labels = FALSE)
text(seq(1, length(nomClass), by=1), par("usr")[3]-0.05, labels=nomClass, srt = 45,
     xpd = TRUE,adj=c(1,0.5),cex=.8,font=2)
legend(length(nomClass)+2,0.8,xpd=TRUE,legend=c("t-test","Limma","Wilcoxon","Lasso","Random Forest","Embedded","stepwise"), col=unique(colors),pch=15,title=expression(bold("Feature Selection")),bty="n")
par(mar=marge)

adjust<- -0.1
mtext("mean",side=3,adj=adjust,line=2,font=2,cex=.7)
adjust<-0.02
for(i in seq(1,length(va),by=2)){
  mtext(sprintf("%.2f",mean(va[[i]],na.rm=TRUE)),side=3,adj=adjust,line=3,cex=.7)
  adjust<-adjust+0.08
}
adjust<-0.06
for(i in seq(2,length(va),by=2)){
  mtext(sprintf("%.2f",mean(va[[i]],na.rm=TRUE)),side=3,adj=adjust,line=1,cex=.7)
  adjust<-adjust+0.08
}

# sd
adjust<- -0.093
mtext("(sd)",side=3,adj=adjust,line=1,font=2,cex=.7)
adjust<-0.01
for(i in seq(1,length(va),by=2)){
  mtext(paste("(",sprintf("%.2f",sd(va[[i]],na.rm=TRUE)),")",sep=""),side=3,adj=adjust,line=2,cex=.7)
  adjust<-adjust+0.0812
}
adjust<-0.055
for(i in seq(2,length(va),by=2)){
  mtext(paste("(",sprintf("%.2f",sd(va[[i]],na.rm=TRUE)),")",sep=""),side=3,adj=adjust,line=0,cex=.7)
  adjust<-adjust+0.0812
}


######################################

extracted[[1]]$y[[1]]
mmetric(as.factor(extracted[[1]]$y[[1]]),as.factor(extracted[[1]]$yhat[[1]]),"CONF")

mmetric(as.factor(extracted[[1]]$y[[1]]),as.factor(extracted[[1]]$yhat[[1]]),"CE")
va<-ev2(extracted,measure="mis")
boxplot(va,las=1,ylim=c(0,1),col="blue")
va<-ev2(extracted,measure="sen")
boxplot(va,las=1,ylim=c(0,1),col="blue")



###
par(mfrow=c(1,3))
extracted<-extractList(classifierlist)
va<-ev(extracted,measure="mis")
boxplot(va,las=1,ylim=c(0,1),col="blue")
extracted<-extractList2(classifierlist)
va<-ev2(extracted,measure="mis")
boxplot(va,las=1,ylim=c(0,1),col="red")
compare(classifierlist,plot=TRUE,measure="misclassification",ylim=c(0,1))
###
par(mfrow=c(1,3))
extracted<-extractList(classifierlist)
va<-ev(extracted,measure="sen2")
boxplot(va,las=1,ylim=c(0,1),col="blue")
extracted<-extractList2(classifierlist)
va<-ev2(extracted,measure="sen")
boxplot(va,las=1,ylim=c(0,1),col="red")
compare(classifierlist,plot=TRUE,measure="sensitivity",ylim=c(0,1))
###
extracted<-extractList(classifierlist)
va<-ev(extracted,measure="spe")
boxplot(va,las=1,ylim=c(0,1),col="blue")
extracted<-extractList2(classifierlist)
va<-ev2(extracted,measure="spe")
boxplot(va,las=1,ylim=c(0,1),col="red")
compare(classifierlist,plot=TRUE,measure="specificity",ylim=c(0,1))








extracted<-extractList2(classifierlist)
va<-ev2(extracted,measure="brier")
boxplot(va,las=1,col="blue")

va<-ev(extracted,measure="brier score")
boxplot(va,las=1,ylim=c(0,1))


# http://support.minitab.com/en-us/minitab/17/topic-library/basic-statistics-and-graphs/tables/other-statistics-and-tests/using-kappa-statistics-and-kendall-s-coefficients/