#' boxplotWithDots
#'
#' This function allows plottin a dotplot and eventually overplotting a boxplot on it
#'
#' @param myExpres Numeric values
#' @param lev A factor for drawing the plot as plot(expres~lev). Pay attention to the order of values in the expres vector and levels in the factor
#' @param aTitle The title for the plot
#' @param groupLabels Names for the groups
#' @param addBoxplot Define if a boxplot has to be overdrawn (set to TRUE when the number of points is big enough)
#' @importFrom graphics boxplot
#' @importFrom beeswarm beeswarm
#' @examples
#' expres <-c(rnorm(10,5,2), rnorm(10,10,2))
#' trat <- as.factor(c(rep("CT",10), rep("TR",10)))
#' titol <- "Treatment effect"
#' groupLab<- c("Control", "Treatment")
#' boxplotWithDots(myExpres=expres, lev=trat, aTitle=titol, groupLabels=groupLab)
#' @export
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
