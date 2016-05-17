#' Extract a list (a character vector indeed) of genes from a topTable outputted by the limma software
#' @param aTopTable A top table file such as produced by the limma package or the UEB pipeline
#' @param filename Name of the file that contains the topTable (to be used if the topTable is not in memory yet)
#' @param entrezOnly To decide if only lines with a meaningful entrez symbol should be retained. Defaults to TRUE.
#' @param uniqueIds To decide if duplicate identifiers should be removed after the selection defined by other parameters. Defaults to TRUE.
#' @param adjOrrawP To decide if filtering is based on raw or adjusted p-values. Defaults to 'adj'.
#' @param Pcutoff The filtering cutoff for p-values (if these are raw or adjusted is set by 'adjorRaw' parameter).
#' @param FCcutoff The filtering cutoff for fold change.
#' @param updown To decide which values to retain : up-regulated (FoldC >0 ), downregulated (FoldC < 0) or both. Defaults to 'both'.
#' @param id2Select To decide which identifiers are to be selected. Defaults to 'ENTREZ'. Possible values are: ("ENTREZ","SYMBOL", NULL)  but the alias ("Entrez", "EntrezsA", "Symbols", "SymbolsA") are accepted.
#' @param cols2Select To decide which (and how many) columns are to be selected. Set 'cols2Select' to a (single or vector) value if no preferred identifier is known. It is used ony when 'id2select' is set to NULL. This parameter is dangerous because if more than 2 columns are selected the output is not a character vector and cannot be used as such.
#' @keywords genelists, filtering
#' @seealso limma
#' @export
#' @examples
#' fileName<- system.file("extdata", "ExpressAndTop_AvsB.csv2", package = "geneLists")
#' AvsB <- read.table(fileName, header=TRUE, sep=";", dec=",", row.names=1)
#' entrezs_01_up  <- genesFromTopTable (AvsB, entrezOnly = TRUE, uniqueIds=TRUE, adjOrrawP = "adj", Pcutoff = 0.01, FCcutoff=1, updown="up", id2Select = "ENTREZ", cols2Select =0)
#' length(entrezs_01_up)
#' table_01_up  <- genesFromTopTable (AvsB, entrezOnly = TRUE, uniqueIds=TRUE, adjOrrawP = "adj", Pcutoff = 0.01, FCcutoff=1, updown="up", id2Select = NULL, cols2Select =1:3)
#' dim(table_01_up)
genesFromTopTable <- function (aTopTable,
                               filename=NULL,
                               entrezOnly = TRUE,
                               uniqueIds=TRUE,
                               adjOrrawP = "adj",
                               Pcutoff = 0.05,
                               FCcutoff = 1,
                               updown = "both",
                               id2Select = "ENTREZ",
                               cols2Select = 0)
{
    if (! is.null(filename)){
      topTab <- read.csv(filename, header=TRUE, sep=";", dec=",", row.names=1)
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
  if (!is.null(id2Select)) {
    if (toupper(substr(id2Select,1,6))=="SYMBOL"){
    geneList <- topTab[,"SymbolsA"]
   }else{
    if (toupper(substr(id2Select,1,6))=="ENTREZ"){
      geneList <- topTab[,"EntrezsA"]
      }
    }
  }else{
      if(sum(cols2Select)> 0){
        geneList <- topTab[,cols2Select]
      }else{
        geneList <-rownames(topTab)
      }
  }
if (length(cols2Select)==1){
      if(uniqueIds) geneList <- unique(geneList)
      geneList <- as.character(geneList)
}
return(geneList)
}


