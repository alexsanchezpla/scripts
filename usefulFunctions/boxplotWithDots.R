boxplotWithDots<- function (myExpres, lev, aTitle, groupLabels, addBoxplot=TRUE)
{
  beeswarm(myExpres~lev, 
           ylab="Expression", xlab="Groups",
           main=aTitle,
           labels=groupLabels)  
  if(addBoxplot)
    boxplot(myExpres~lev, add = T, names = c("",""), col="#0000ff22")
  # Segons un post de: https://www.r-statistics.com/2011/03/beeswarm-boxplot-and-plotting-it-with-r/
}
