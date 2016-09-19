#' Representant PCAs
#'
#' @param X El nom de la matriu de la que es calculara el PCA (de t(X))
#' @param labels Vector d'etiquetes per a la visualització dels punts. Per defecte es NULL
#' @param colors Vector de colors per als punts. De la mateixa mida que el de les etiquetes.
#' @param dataDesc Breu descripcio de les dades per al titol
#' @param scale Variable logica per decidir si s'escalen les dades
#' @param transpose Variable logica per decidir si es vol fer el PCA de la matriu original o de la seva trasposta
#' @keywords Multivariate, Plots
#' @examples
#' data(USArrests)
#' plotPCA (X=USArrests, scale=TRUE, transpose=TRUE, dataDesc="Arrestos en EUA",  labels=rownames(USArrests))
#' plotPCA()

plotPCA <- function ( X, labels=NULL, colors=NULL, dataDesc="", scale=FALSE, transpose=FALSE)
{
  if (transpose) X <- t(X)
  pcX<-prcomp(X, scale=scale)
  loads<- round(pcX$sdev^2/sum(pcX$sdev^2)*100,1)
  xlab<-c(paste("PC1",loads[1],"%"))
  ylab<-c(paste("PC2",loads[2],"%"))
  if (is.null(colors)) colors=1
  plot(pcX$x[,1:2],xlab=xlab, ylab=ylab, col=colors)
  if (!(is.null(labels))) text(pcX$x[,1], pcX$x[,2], labels, pos=3, cex=0.8)
  title(paste("Plot of first 2 PCs for values in", dataDesc, sep=" "), cex=0.8)
}
