expres <-c(rnorm(10,5,2), rnorm(10,10,2))
trat <- as.factor(c(rep("CT",10), rep("TR",10)))
titol <- "Treatment effect"
groupLab<- c("Control", "Treatment")
require(beeswarm)
boxplotWithDots(myExpres=expres, lev=trat, aTitle=titol, groupLabels=groupLab, addBoxplot = FALSE)
boxplotWithDots(myExpres=expres, lev=trat, aTitle=titol, groupLabels=groupLab, addBoxplot = TRUE)
