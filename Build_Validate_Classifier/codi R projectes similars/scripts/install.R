## to install R packages



instifnot<-function(paquet){
  if(!require(paquet,character.only = TRUE)){
    install.packages(paquet)
    library(paquet,character.only=TRUE)
  }
}


#source("http://bioconductor.org/biocLite.R") # load bioconductor (requires internet connection)
# once the package is download, the internet connection is not requred anymore
instifnotbio<-function(paquet){
  if(!require(paquet,character.only = TRUE)){
    biocLite(paquet,suppressUpdates=TRUE) 
  }
  
}