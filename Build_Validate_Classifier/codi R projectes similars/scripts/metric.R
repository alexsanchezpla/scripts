y<-c(1,1,1,0,0,0,1)
yhat<-c(1,1,1,0,0,0,0)
mmetric(as.factor(y),as.factor(yhat),metric="CONF")
mmetric(as.factor(y),as.factor(yhat),metric="TPR")


### regression examples: y - desired values; x - predictions
y=c(95.01,96.1,97.2,98.0,99.3,99.7);x=95:100
print(mmetric(y,x,"ALL"))
print(mmetric(y,x,"MAE"))
m=mmetric(y,x,c("MAE","RMSE","RAE","RSE"))
print(m)
cat(names(m)[3],"=",round(m[3],digit=2),"\n",sep="")
print(mmetric(y,x,c("COR","R2","Q2")))
print(mmetric(y,x,c("TOLERANCE","NAREC"),val=0.5))
print(mmetric(y,x,"THEILSU2",val=94.1)) # val = 1-ahead random walk, c(y,94.1), same as below
print(mmetric(y,x,"THEILSU2",val=c(94.1,y[1:5]))) # val = 1-ahead random walk (previous y values)
print(mmetric(y,x,"MASE",val=c(88.1,89.9,93.2,94.1))) # val = in-samples
val=vector("list",length=4)
val[[2]]=0.5;val[[3]]=94.1;val[[4]]=c(88.1,89.9,93.2,94.1)
print(mmetric(y,x,c("MAE","NAREC","THEILSU2","MASE"),val=val))
# user defined error function example:
# myerror = number of samples with absolute error above 0.1% of y: 
myerror=function(y,x){return (sum(abs(y-x)>(0.001*y)))}
print(mmetric(y,x,metric=myerror))
# example that returns a list since "REC" is included:
print(mmetric(y,x,c("MAE","REC","TOLERANCE"),val=1))

### pure binary classification 
y=factor(c("a","a","a","a","b","b","b","b"))
x=factor(c("a","a","b","a","b","a","b","a"))
print(mmetric(y,x,"CONF")$conf)
print(mmetric(y,x,"ALL"))
print(mmetric(y,x,metric=c("ACC","TPR","ACCLASS")))

mmetric(as.factor(extracted[[1]]$y[[1]]),as.factor(extracted[[1]]$yhat[[1]]),"CONF")$conf
mmetric(as.factor(extracted[[1]]$y[[1]]),as.factor(extracted[[1]]$yhat[[1]]),"ALL")
mmetric(as.factor(extracted[[1]]$y[[1]]),as.factor(extracted[[1]]$yhat[[1]]),"ALL",TC=0)

comparison<-compare(classifierlist,plot=TRUE,measure="auc")
comparison

eva<-evaluation(classifierlist[[1]],measure="brier score")
eva
slotNames(eva)
sd(eva@score)/sqrt(length(eva@score))





aa<-list()
for(i in 1:CV5f@iter){
  aa[[i]]<-round(mmetric(as.factor(extracted[[1]]$y[[i]]),as.factor(extracted[[1]]$yhat[[i]]),"ALL",TC=1),2)
}
 
df <- data.frame(matrix(unlist(aa), nrow=length(aa), byrow=T))
apply(df,2,mean,na.rm=TRUE)

stderr <- function(x){
  sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
}

apply(df,2,stderr)
#

xtable(df)
kable(df[,c(1,4,8)])
aa<-round(mmetric(as.factor(extracted[[1]]$y[[1]]),as.factor(extracted[[1]]$yhat[[1]]),"ALL",TC=1),2)




# wrap all the metrics
aa<-list()
bb<-list()
for(j in 1:length(classifierlist)){  # for each classifier
  for(i in 1:CV5f@iter){  # for each iteration
  aa[[i]]<-mmetric(as.factor(extracted[[j]]$y[[i]]),as.factor(extracted[[j]]$yhat[[i]]),"ALL",TC=1)
  }
  bb[[j]]<-data.frame(matrix(unlist(aa), nrow=length(aa), byrow=T))
}
names(bb[[1]])<-names(aa[[1]])
 

for(i in 1:length(bb)){
  column<-dim(bb[[i]])[2]
  bb[[i]][,column+1]<-evaluation(classifierlist[[i]],measure="misclassification")@score 
  bb[[i]][,column+2]<-evaluation(classifierlist[[i]],measure="sensitivity")@score  
  bb[[i]][,column+3]<-evaluation(classifierlist[[i]],measure="specificity")@score  
  bb[[i]][,column+4]<-evaluation(classifierlist[[i]],measure="auc")@score  
  bb[[i]][,column+5]<-evaluation(classifierlist[[i]],measure="brier score")@score  
  bb[[i]][,column+6]<-evaluation(classifierlist[[i]],measure="average probability")@score  
}
names(bb[[1]])[18:23]<-c("Miscalss","sensitivity","specificity","auc","brier","average prob")
 



# mean and standard deviation
M<-list()
S<-list()
for(i in 1:length(bb)){
  M[[i]]<-apply(bb[[i]],2,mean) 
  S[[i]]<-apply(bb[[i]],2,stderr)
}


M2<-do.call("rbind",M)
S2<-do.call("rbind",S)


do.call(function(x)apply(x,2,mean) ,bb)

merge.all <- function(x, y)
  merge(x, y, all=TRUE, by="Sample")

out <- Reduce(merge.all, bb)
dd<-do.call("rbind", bb)
apply(dd,2,mean,na.rm=TRUE)
apply(dd,2,stderr)

M<-apply(dd,2,mean,na.rm=TRUE)
S<-apply(dd,2,stderr)
v1<-seq(1,2*17,2)
v2<-(1:(2*17))[-v1]

dff<-data.frame()
for(i in (17*2-1)){  
  dff[,i]<-M[]
  dff[,i+1]
}



### probabilities binary classification 
y=factor(c("a","a","a","a","b","b","b","b"))
px=matrix(nrow=8,ncol=2)
px[,1]=c(1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3)
px[,2]=1-px[,1]
print(px)
print(mmetric(y,px,"CONF")$conf)
print(mmetric(y,px,"CONF",D=0.5,TC=2)$conf)
print(mmetric(y,px,"CONF",D=0.3,TC=2)$conf)
print(mmetric(y,px,metric="ALL",D=0.3,TC=2))
print(mmetric(y,px,metric=c("ACC","AUC","AUCCLASS","BRIER","BRIERCLASS","CE"),D=0.3,TC=2))

### pure multi-class classification 
y=factor(c("a","a","b","b","c","c"))
x=factor(c("a","a","b","c","b","c"))
print(mmetric(y,x,metric="CONF")$conf)
print(mmetric(y,x,metric="CONF",TC=-1)$conf)
print(mmetric(y,x,metric="CONF",TC=1)$conf)
print(mmetric(y,x,metric="ALL"))
print(mmetric(y,x,metric=c("ACC","ACCLASS","KAPPA")))
print(mmetric(y,x,metric=c("ACC","ACCLASS","KAPPA"),TC=1))

### probabilities multi-class 
y=factor(c("a","a","b","b","c","c"))
px=matrix(nrow=6,ncol=3)
px[,1]=c(1.0,0.7,0.5,0.3,0.1,0.7)
px[,2]=c(0.0,0.2,0.4,0.7,0.3,0.2)
px[,3]=1-px[,1]-px[,2]
print(px)
print(mmetric(y,px,metric=c("AUC","AUCCLASS","NAUC"),TC=-1,val=0.1))
print(mmetric(y,px,metric=c("AUC","NAUC"),TC=3,val=0.1))
print(mmetric(y,px,metric=c("ACC","ACCLASS"),TC=-1))
print(mmetric(y,px,metric=c("CONF"),TC=3,D=0.5)$conf)
print(mmetric(y,px,metric=c("ACCLASS"),TC=3,D=0.5))
print(mmetric(y,px,metric=c("CONF"),TC=3,D=0.7)$conf)
print(mmetric(y,px,metric=c("ACCLASS"),TC=3,D=0.7))

### ordinal multi-class (example in Ricardo Sousa PhD thesis 2012)
y=ordered(c(rep("a",4),rep("b",6),rep("d",3)),levels=c("a","b","c","d"))
x=ordered(c(rep("c",(4+6)),rep("d",3)),levels=c("a","b","c","d"))
print(mmetric(y,x,metric="CONF")$conf)
print(mmetric(y,x,metric=c("CE","MAEO","MSEO","KENDALL")))
# note: only y needs to be ordered
x=factor(c(rep("b",4),rep("a",6),rep("d",3)),levels=c("a","b","c","d"))
print(mmetric(y,x,metric="CONF")$conf)
print(mmetric(y,x,metric=c("CE","MAEO","MSEO","KENDALL")))

### ranking (multi-class) 
y=matrix(nrow=1,ncol=12);x=y
# http://www.youtube.com/watch?v=D56dvoVrBBE
y[1,]=1:12
x[1,]=c(2,1,4,3,6,5,8,7,10,9,12,11)
print(mmetric(y,x,metric="KENDALL"))
print(mmetric(y,x,metric="ALL"))

y=matrix(nrow=2,ncol=7);x=y
y[1,]=c(2,6,5,4,3,7,1)
y[2,]=7:1
x[1,]=1:7
x[2,]=1:7
print(mmetric(y,x,metric="ALL"))

### mining, several runs, prob multi-class
## Not run: 
data(iris)
M=mining(Species~.,iris,model="rpart",Runs=2)
R=mmetric(M,metric="CONF",aggregate="no")
print(R[[1]]$conf)
print(R[[2]]$conf)
print(mmetric(M,metric="CONF",aggregate="mean"))
print(mmetric(M,metric="CONF",aggregate="sum"))
print(mmetric(M,metric=c("ACC","ACCLASS"),aggregate="no"))
print(mmetric(M,metric=c("ACC","ACCLASS"),aggregate="mean"))
print(mmetric(M,metric="ALL",aggregate="no"))
print(mmetric(M,metric="ALL",aggregate="mean"))

## End(Not run)

### mining, several runs, regression
## Not run: 
data(sin1reg)
S=sample(1:nrow(sin1reg),40)
M=mining(y~.,data=sin1reg[S,],model="ksvm",search=2^3,Runs=10)
R=mmetric(M,metric="MAE")
print(mmetric(M,metric="MAE",aggregate="mean"))
miR=meanint(R) # mean and t-student confidence intervals
cat("MAE=",round(miR$mean,digits=2),"+-",round(miR$int,digits=2),"\n")
print(mmetric(M,metric=c("MAE","RMSE")))
print(mmetric(M,metric=c("MAE","RMSE"),aggregate="mean"))
R=mmetric(M,metric="REC",aggregate="no")
print(R[[1]]$rec)
print(mmetric(M,metric=c("TOLERANCE","NAREC"),val=0.2))
print(mmetric(M,metric=c("TOLERANCE","NAREC"),val=0.2,aggregate="mean"))

## End(Not run)