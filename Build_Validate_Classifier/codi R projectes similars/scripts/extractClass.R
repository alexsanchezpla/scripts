
# Function to extract the prob, yhat and y for each model and iteration
# from a cloutput object (CMA):: classification

extract<-function(classificator){
  prob<-list() # probability to be assigned in group 0
  yhat<-list() # predicted value (0 or 1)
  y<-list() # real value/label (0 or 1)
  for(i in 1:length(classificator)){
    prob[[i]]<-classificator[[i]]@prob[,1]
    yhat[[i]]<-classificator[[i]]@yhat
    y[[i]]<-classificator[[i]]@y
  }  
  return(list(prob=prob,yhat=yhat,y=y))
}


###
# Function to perform extract for several classifiers

extractList<-function(llistat){
  ex<-lapply(classifierlist,extract)
  return(ex)
}




# Function to calculate misclassification
mis<-function(y, yhat){
  tab<-table(y,yhat)
  mis<-1-(sum(diag(tab))/sum(tab))
  return(mis)
}


# Function to calculate sensitivity
sen<-function(y,yhat){
  tab<-table(y,yhat)
  if(dim(tab)[2]==2){
    sen<-tab[2,2]/(tab[2,2]+tab[2,1]) 
  }else{
    sen<-NA
  }
  return(sen)
}

sen2<-function(y,yhat){
  tab<-table(y,yhat,useNA="always")
  sen<-tab[2,2]/(tab[2,2]+tab[2,1]) 
  return(sen)
}

spe<-function(y,yhat){
  tab<-table(y,yhat,useNA="always")
  spe<-tab[1,1]/(tab[1,1]+tab[1,2]) 
  return(spe)
}

kappa1<-function(y,yhat){
  tab<-table(y,yhat,useNA="always")
  kappa<-as.numeric(Kappa(tab)[[1]][1])
  return(kappa)
}


auc1<-function(y,prob){
  auc<-as.numeric(roc(y,prob)$auc)
  return(auc)
}

ev<-function(class,measure="mis"){ 
  tots<-list()
  for(j in 1:length(class)){   # classificator 
    valors<-list()
    for(i in 1:length(class[[1]][[1]])){ # iteration
      if(measure=="mis"){
        valors[[i]]<-mis(class[[j]]$y[[i]],class[[j]]$yhat[[i]])
      } else if(measure=="sen"){
        valors[[i]]<-sen(class[[j]]$y[[i]],class[[j]]$yhat[[i]])
      }   else if (measure=="sen2"){
        valors[[i]]<-sen2(class[[j]]$y[[i]],class[[j]]$yhat[[i]])
      }  else if (measure=="spe"){
        valors[[i]]<-spe(class[[j]]$y[[i]],class[[j]]$yhat[[i]])
      }  else if (measure=="kappa1"){
        if((sum(class[[j]]$y[[i]])==0) || (sum(class[[j]]$yhat[[i]]))==0){
          valors[[i]]<-NA           
        }else{
          valors[[i]]<-kappa1(class[[j]]$y[[i]],class[[j]]$yhat[[i]]) 
        }         
      }  else if (measure=="auc1"){
        valors[[i]]<-auc1(class[[j]]$y[[i]],class[[j]]$prob[[i]])
      }
    }
    tots[[j]]<-unlist(valors)
  }  
  return(tots)
}

###################################
# using mmetric (rminer) function


extract2<-function(classificator){
  prob<-list() # probability to be assigned in group 0
  yhat<-list() # predicted value (0 or 1)
  y<-list() # real value/label (0 or 1)
  for(i in 1:length(classificator)){
    prob[[i]]<-classificator[[i]]@prob
    yhat[[i]]<-classificator[[i]]@yhat
    y[[i]]<-classificator[[i]]@y
  }  
  return(list(prob=prob,yhat=yhat,y=y))
}


###
# Function to perform extract for several classifiers

extractList2<-function(llistat){
  ex<-lapply(classifierlist,extract2)
  return(ex)
}






ev2<-function(class,measure="mis"){ 
  tots<-list()
  for(j in 1:length(class)){   # classificator 
    valors<-list()
    for(i in 1:length(class[[1]][[1]])){ # iteration
      if(measure=="mis"){
        valors[[i]]<-mmetric(as.factor(extracted[[j]]$y[[i]]),as.factor(extracted[[j]]$yhat[[i]]),"CE")/100
       } else if(measure=="sen"){
        valors[[i]]<-mmetric(as.factor(extracted[[j]]$y[[i]]),as.factor(extracted[[j]]$yhat[[i]]),"TPR")[1]/100
       }   else if (measure=="spe"){
        valors[[i]]<-mmetric(as.factor(extracted[[j]]$y[[i]]),as.factor(extracted[[j]]$yhat[[i]]),"TPR")[1]/100

       }  else if (measure=="brier"){
         valors[[i]]<-mmetric(as.factor(extracted[[j]]$y[[i]]),as.matrix(extracted[[j]]$prob[[i]]),"BRIER")/100
         #       }  else if (measure=="kappa1"){
#         if((sum(class[[j]]$y[[i]])==0) || (sum(class[[j]]$yhat[[i]]))==0){
#           valors[[i]]<-NA           
#         }else{
#           valors[[i]]<-kappa1(class[[j]]$y[[i]],class[[j]]$yhat[[i]]) 
#         }         
#       }  else if (measure=="auc1"){
#         valors[[i]]<-auc1(class[[j]]$y[[i]],class[[j]]$prob[[i]])
       }
    }
    tots[[j]]<-unlist(valors)
  }  
  return(tots)
}


