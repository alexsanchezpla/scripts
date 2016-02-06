# Function to extend the 'toplist' CMA function, since itjust gives data from one iteration at a time
# This wrapper gives all the most used variables ALL the iterations.

#var_sel_met: feature selection output from an object 'gensel' (CMA) :: GeneSelection
# nvar: Number of variables to retain
# CV: cross-validation output from object class 'learningsets':: GenerateLearningsets

topVariables<-function(var_sel_met,nvar=10,CV){
  llista<-list()
  for(i in 1:CV@iter){
    llista[[i]]<-unlist(toplist(var_sel_met,iter=i,show=FALSE)$index)
  }
  llista<-unlist(llista)
  llista<-sort(table(llista),decreasing = TRUE)[1:nvar]
  names(llista)<-as.numeric(names(llista))+1 # +1 because in 1 there is the dummy variable
  return (llista)
}



############################################################################################
# Function to see which variables appear at the top X of ALL feature selection methods in ALL iterations
# plot Venn diagram, uo to 5 comparisons.

# data: dataframe with all the data
# CV: cross-validation output
# llista: object from the output of CMA function GeneSelection (genesel type)
# nomMethod: array with the names of the feature selection method employed
# require package gplots

commonTopVariables<-function(data,CV,llista,nomMethod=c("A","B","C","D","E"),nVar=10){
  vars<-list()
  for(i in 1:length(llista)){
    vars[[i]]<-as.numeric(names(topVariables(llista[[i]],nvar=nVar,CV=CV)))
    names(vars)[[i]]<-nomMethod[i] # for the venn diagram
  }
  # to create a Venn diagram
  par(mar=c(2,2,2,2),cex=0.8,cex.lab=1)
  q<-venn(vars,show.plot=TRUE)
  
  #function to see which variables are the same (have variables in common)
  similarTopVar<-function(llistaMetodes){
    colnum<-Reduce(intersect, llistaMetodes)
    return(colnum)
  }
  
  colnum<-similarTopVar(vars)
  varname<-names(data[,colnum])
  output<-list(colNum=colnum,varName=varname,varMethods=vars,venn=q)
  return(output)
}



