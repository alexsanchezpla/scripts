#' Extract a list (a character vector indeed) of genes from a topTable outputted by the limma software
#' @param aTopTable
#' @param filename
#' @param entrezOnly
#' @param uniqueIds
#' @param adjOrrawP
#' @param Pcutoff
#' @param FCcutoff
#' @param updown
#' @param id2select
#' @param cols2select
#' @keywords
#' @seealso
#' @export
#' @examples
#'

genesFromTopTable <- function (aTopTable,
                               filename=NULL, 
                               entrezOnly = TRUE, 
                               uniqueIds=TRUE, 
                               adjOrrawP = "adj",
                               Pcutoff = 0.05,
                               FCcutoff = 1,
                               updown = "both", # c("both","up","down")
                               id2Select = "ENTREZ",  # c("ENTREZ","SYMBOL", "other") 
                                                      # Alias ("Entrez", "EntrezsA", "Symbols", "SymbolsA") are accepted
                               cols2Select = 2){
  
  if (! is.null(filename)){
    topTab <- read.csv(file.path(listsDir, filename), head=TRUE, sep=";", dec=",", row.names=1)
  }else{
    topTab=aTopTable
  }
  if (entrezOnly) {
    selectedEntrez <- topTab$EntrezsA !="---" # !is.na(topTab[,"EntrezsA"])
    topTab <- topTab[selectedEntrez,]
  }
  if (Pcutoff < 1){
    if (adjOrrawP=="adj"){
      selectedP <- topTab[,"adj.P.Val"] < Pcutoff
    }else{
      selectedP <- topTab[,"P.Value"] < Pcutoff
    }
    topTab<- topTab[selectedP, ]
  }
  if (FCcutoff > 0){
    selectedFC <-(abs(topTab[,"logFC"]) > FCcutoff)
    topTab<- topTab[selectedFC, ]
  }
  if (updown!="both"){
    if (updown=="up"){
      selectedup <-topTab[,"logFC"] > 0     
      topTab<- topTab[selectedup, ]
    }else{
      if (updown=="down"){
        selecteddown <-topTab[,"logFC"] < 0
        topTab<- topTab[selecteddown, ]
      }
    }
  }
  if (toupper(substr(id2Select,1,6))=="SYMBOL"){
    geneList <- topTab[,"SymbolsA"]
  }else{
    if (toupper(substr(id2Select,1,6))=="ENTREZ"){
      geneList <- topTab[,"EntrezsA"]
    }else{
      if(sum(cols2Select)> 0){
        geneList <- topTab[,cols2Select]
      }else{
        geneList <-rownames(topTab)
      }}}
  # length(geneList)
  if (length(cols2Select)==1){
      if(uniqueIds) geneList <- unique(geneList)
      geneList <- as.character(geneList)
    }
  return(geneList)
}


