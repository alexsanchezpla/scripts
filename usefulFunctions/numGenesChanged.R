#' Given a limma-outputted topTable (TT) it returns the number of up or down regulated genes that woud be returned if the cutoff was set at different values.
#' @param TT A top table object such as produced by the limma package or the UEB pipeline.
#' @param cName Name to the comparison that produced the top Table. Defaults to "comparison".
#' @keywords genelists, filtering
#' @seealso limma
#' @export
#' @examples
#' AvsB <- read.table("https://raw.githubusercontent.com/alexsanchezpla/scripts/master/Gene_List_Management/dades/ExpressAndTop_AvsB.csv2", head=T, sep=";", dec=",", row.names=1)
#' genesChanged <- numGenesChanged (AvsB, "Group A vs group B")
  numGenesChanged <- function (TT, cName="comparison"){
    Bup   <- sum(TT$t>0 & TT$B >0 )
    Bdown <- sum(TT$t <=0 & TT$B >0 )
    adjP001Up <- sum(TT$t>0 & TT$adj.P.Val < 0.01 )
    adjP001Down <- sum(TT$t< 0 & TT$adj.P.Val < 0.01 )
    adjP005Up <-   sum(TT$t>0 & TT$adj.P.Val < 0.05 )
    adjP005Down <- sum(TT$t< 0 & TT$adj.P.Val < 0.05 )
    adjP025Up <-   sum(TT$t>0 & TT$adj.P.Val < 0.25 )
    adjP025Down <- sum(TT$t< 0 & TT$adj.P.Val < 0.25 )    
    P001Up <-   sum(TT$t>0 & TT$P.Value < 0.01 )
    P001Down <- sum(TT$t< 0 & TT$P.Value < 0.01 )    
    P005Up <-   sum(TT$t>0 & TT$P.Value < 0.05 )
    P005Down <- sum(TT$t< 0 & TT$P.Value < 0.05 ) 
    nGenes <- data.frame(comparisonName= c(Bup, Bdown, 
                                           adjP001Up, adjP001Down, adjP005Up, adjP005Down, adjP025Up, adjP025Down, 
                                           P001Up, P001Down,  P005Up, P005Down))
    rowNames <- c("upReg-B>0", "downReg-B>0", 
                  "upReg-Adjusted-p-val < 0.01", "downReg-Adjusted-p-val < 0.01",
                  "upReg-Adjusted-p-val < 0.05", "downReg-Adjusted-p-val < 0.05",
                  "upReg-Adjusted-p-val < 0.25", "downReg-Adjusted-p-val < 0.25",
                  "upReg-P value < 0.01 ", "downReg-P value < 0.01", 
                  "upReg-P value < 0.05", "downReg-P value < 0.05")
    rownames(nGenes)   <-rowNames
    colnames(nGenes)[1] <- cName
    return(nGenes)
  }
