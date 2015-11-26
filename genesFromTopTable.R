#' Extract a list (a characyer vector indeed) of genes from a topTable outputted by the limma software
#' @param filename
#' @param entrezOnly
#' @param adjOrrawP
#' @param Pcutoff
#' @param FCcutoff
#' @param id2select
#' @param cols2select
#' @param return
#' @keywords
#' @seealso
#' @export
#' @examples
#'
genesFromTopTable <- function (filename, entrezOnly = TRUE, 
                               adjOrrawP = "adj",
                               Pcutoff = 0.05,
                               FCcutoff = 1,
                               id2Select = "EntrezsA", # c("EntrezsA","SYMBOL", "other")
                               cols2Select = 2){
  
  topTab <- read.csv(file.path(listsDir, filename), head=TRUE, sep=";", dec=",", row.names=1)
  dim(topTab)
  colnames(topTab)
  if (entrezOnly) {
    selectedEntrez <- !is.na(topTab[,"EntrezsA"])
    topTab <- topTab[selectedEntrez,]
  }
  dim(topTab)
  if (Pcutoff < 1){
    if (adjOrrawP=="adj"){
      selectedP <- topTab[,"adj.P.Val"] < Pcutoff
    }else{
      selectedP <- topTab[,"P.Value"] < Pcutoff
    }
    topTab<- topTab[selectedP, ]
  }
  dim(topTab)
  if (FCcutoff > 0){
    selectedFC <-(abs(topTab[,"logFC"]) > FCcutoff)
    topTab<- topTab[selectedFC, ]
    dim(topTab)
  }
  
  if (id2Select=="SYMBOL"){
    geneList <- topTab[,"SYMBOL"]
  }else{
    if (id2Select=="Entrez"){
      geneList <- topTab[,"EntrezsA"]
    }else{
      if(cols2Select> 0){
        geneList <- topTab[,cols2Select]
      }else{
        geneList <-rownames(topTab)
      }}}
  length(geneList)
  return(geneList)
}
