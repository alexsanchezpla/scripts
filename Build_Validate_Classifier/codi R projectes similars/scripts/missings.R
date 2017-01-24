# Function to observe the missings pattern
# dades: Data frame with the data
# type: 1 to study observations; 2 for variables
# per: boolean. TRUE if the percentage is to be shown
# perVar: percentage threshold to used to plot
# plot: boolean. Whether the plot should be done

numMiss<-function(dades,type=1,per=TRUE,perVar=20,plot=TRUE){
  # 1: by rows (observations)
  # 2: by columns (variables)
  if(type==1){
    type<-1
    dim<-2
    #noms<-paste("ind",1:36,sep="-")
    noms<-1:dim(dades)[1]
  } else if (type==2){
    type<-2
    dim<-1
    vars<-which(unlist(lapply(dades,class))!="factor") # avoid non numeric columns
    dades<-dades[,vars]
    noms<-names(dades)
  } else{
    stop("type must be 1 (observations) or 2 (variables)")
  }
  mis<-apply(dades,type,function(x){sum(is.na(x))}) # number of missing by row or column (type 1 or 2)
  names(mis)<-noms
  mis_per<-mis*100/dim(da)[dim] # calculate the percentage for each row or column (type 1 or 2)
  # plot by percentage of missings
  if(plot==TRUE){
    if(per==TRUE){  # plot by percentage
      if(type==2){      # variables
        barplot(mis_per[which(mis_per>perVar)],las=1,ylim=c(0,100),"Variable",
                col="#dedbdb",border="#bebebe")      
      } else{           # observations
        barplot(mis_per,las=1,ylim=c(0,100),xlab="Observation",ylab="percentage of missings",
                col="#dedbdb",border="#bebebe")      
      } 
    # plot by number of missings  
    } else{
      barplot(mis,las=1) 
    }    
  }else{}  
  return(list(
    num=sort(mis,decreasing = TRUE),
    per=sort(mis_per,decreasing = TRUE)))    
}


missVar<-function(numMiss){
  va<-numMiss$per
  talls<-cut(va,seq(0,max(va)+4,5),include.lowest = TRUE)
  talls<-droplevels(talls)
  barplot(summary(talls),las=1,ylim=c(0,max(summary(talls)*1.1)),
          ylab="Number of variables",xlab="Percentage of missing values",
          col="#dedbdb",border="#bebebe")
}
 


############################################################################
# Antic

# # number of missing per individual  
# ind_mis<-apply(da,1,function(x){sum(is.na(x))})
# names(ind_mis)<-paste("ind",1:36,sep="-")
# barplot(ind_mis,las=1)
# ind_mis_per<-ind_mis*100/dim(da)[2]
# ind_mis_per
# barplot(ind_mis_per,las=1,ylim=c(0,100))
# 
# # number of missing depending on the variable
# var_mis<-apply(da,2,function(x){sum(is.na(x))})
# names(var_mis)<-names(da)
# barplot(var_mis,las=1)
# sum(var_mis>8)
# barplot(var_mis[which(var_mis>8)],las=2)
# 
# var_mis_per<-var_mis*100/dim(da)[1]
# barplot(var_mis_per[which(var_mis_per>25)],las=2,ylim=c(0,100))
# barplot(var_mis_per[which(var_mis_per>24) && which(var_mis_per<25)],las=2,ylim=c(0,100))
