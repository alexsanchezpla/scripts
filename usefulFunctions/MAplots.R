#' Draw or or several MA-plot from an expression data matrix
#' @param myData
#' @param selected
#' @param paired
#' @param xcol
#' @param ycol
#' @keywords MAplot, normalization
#' examples
#' 

oneMAPlot <-  function (myData, selected=1, paired=FALSE, xcol=NULL, ycol=NULL){
  #  opt<- par(mfrow=c(1,1)) 
  if (paired){
    x <- myData[, xcol]+1
    y <- myData[, ycol] +1
  }else{
    medians <-apply(myData, 1, median, na.rm=TRUE)
    selected
    x<- myData[, selected]+1
    y <-medians
  }
  M <- log(x+1)-log(y)
  A <- 0.5*(log(x+1)+log(y))
  if (paired){
    plotTitle = paste("MAplot: ",colnames(myData)[xcol], "vs", colnames(myData)[ycol], sep=" ")
  }else{
    plotTitle = paste("MAplot: ",colnames(myData)[selected],sep=" ")
  }
  plot(A, M, main=plotTitle )  # colnames(myData)[selected])
  abline(h=0)
  #  par(opt)
}

blockMAPlots <- function (myData, rc, paired=FALSE){
  opt<- par(mfrow=c(rc,rc))
  if (paired){} # implementation pending
  medians <-apply(myData,1,median, na.rm=TRUE)
  for (i in 1:ncol(myData))
  {
    x<-myData[,i]
    #x[is.na(x)]<-0
    M <- log(x+1)-log(medians)
    A <- 0.5*(log(x+1)+log(medians))
    name = paste("MAplot: ",colnames(myData)[i],sep=" ")
    plot(A, M, main=colnames(myData)[i])
  }
  par(opt)
}
