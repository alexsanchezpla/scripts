nicePlot<-function(classifierlist,nomClass,measure){
  marge<-par()$mar
  par(mar=c(5.1,3.1,3.1,2.1))
  comparison<-compare(classifierlist,plot=TRUE,measure=measure,
                      names=nomClass,col="gray70",
                      ylim=c(0,1), xaxt="n")
  
  axis(1, at=seq(1,length(nomClass), by=1), labels = FALSE)
   text(seq(1, length(nomClass), by=1), par("usr")[3]-0.05, labels=nomClass, srt = 45,
       xpd = TRUE,adj=c(1,0.5),cex=.8,font=2)
   par(mar=marge)
}


##
# Millorar

# nicePlot<-function(classifierlist,nomClass,measure){
#   marge<-par()$mar
#   par(mar=c(4.1,4.1,5.1,8.1))
#   comparison<-compare(classifierlist,plot=TRUE,measure=measure,
#                       names=nomClass,col="navy",
#                       ylim=c(0,1), xaxt="n")
#   
#   axis(1, at=seq(1,length(nomClass), by=1), labels = FALSE)
#   
#   text(seq(1, length(nomClass), by=1), par("usr")[3]-0.05, labels=nomClass, srt = 45,
#        xpd = TRUE,adj=c(1,0.5),cex=.8,font=2)
#   
#   legend(length(nomClass)+1,0.8,xpd=TRUE,
#          legend=c("t-test","Limma","Wilcoxon","Lasso","Random Forest","Embedded","stepwise"),
#          col="blue",pch=15,title=expression(bold("Feature Selection")),bty="n")
#   par(mar=marge)
# }
