quickAnalysis <- function ( expres, groupingVar, maxGenes=250, 
                            aTitle="", aFName="Results", 
                            outputDir=".", outputFType="xls",
                            pvalThreshold=0.05, useAdjP=TRUE,
                            plotSelected = TRUE,
                            groupNames =NULL)
{
  print(paste("ANALYIS:",aTitle,"...PROCESS STARTED",sep=" "))
  ### Packages needed in both analysis
  
  ### ANALYSIS USING multtest library
  
  resT <- mt.maxT(expres, classlabel = as.factor(groupingVar))
  if (useAdjP){
    resT.selected <-resT[resT$adjp < pvalThreshold,]
  }else{
    resT.selected <-resT[resT$rawp < pvalThreshold,]
  }
  gNames.multtest <-rownames(resT.selected)

  ### ANALYSIS USING limma
  
  design<- model.matrix(~ -1 + factor(groupingVar,levels=unique(groupingVar)))
  colnames(design) <- unique(groupingVar)
  
  numParameters <- ncol(design)
  parameterNames <- colnames(design)
  contrastNames <- paste(parameterNames[2],parameterNames[1],sep="-")
  contrastsMatrix <- matrix(c(-1,1),nrow=ncol(design))
  rownames(contrastsMatrix) <- parameterNames
  colnames(contrastsMatrix) <- contrastNames
  
  fit<-lmFit(expres, design)
  fit.main<-contrasts.fit(fit, contrastsMatrix)
  fit.main<-eBayes(fit.main)
  
  # top.Diff <- topTable(fit.main, n=maxGenes, adjust="fdr")
  top.Diff.all <- topTable(fit.main, n=nrow(expres), adjust="fdr")
  if (useAdjP){
    top.Diff <-top.Diff.all[top.Diff.all$adj.P.Val < pvalThreshold,]
  }else{
    top.Diff <-top.Diff.all[top.Diff.all$P.Value < pvalThreshold,]
  }
  
  gNames.limma <-rownames(top.Diff)
  selectedGenes <-union(gNames.multtest, gNames.limma)
  
  top.Diff.common <- top.Diff.all[selectedGenes, ]
  resT.common <- resT[selectedGenes, ]
  
  cat (paste("\nNumber of genes selected using permutations test (p <0.05): ", 
               nrow(resT.selected),sep=""), "\n")
  cat (paste("Number of genes selected using linear model test (p <0.05): ",
               nrow(top.Diff),sep=""), "\n")
  cat (paste ("Number of genes selected by both test (p <0.05)          : ",
                length(intersect(gNames.multtest,gNames.limma)),sep=""), "\n")
  
  cat (paste("ANALYIS:",aTitle , "... PROCESS COMPLETED",sep=" "), "\n")

  if (outputFType=="xls") {
    xlsFName<- file.path(outputDir,paste(aFName, "xls", sep="."))
    write.xlsx(resT.selected, file=xlsFName, sheetName = "multtest")
    write.xlsx(top.Diff, file=xlsFName, append=TRUE, sheetName = "limma")
  }else{
    aFName1 <- file.path(outputDir,paste(aFName, "multtest","html", sep="."))
    hwrite (resT.selected, page =aFName1)
    aFName1<- file.path(outputDir,paste(aFName, "limma","html", sep="."))
    hwrite (top.Diff, page= aFName1)
  }
  
  if (plotSelected && (length(selectedGenes) >0)){
    plotsFName<- file.path(outputDir,paste(aFName,"Plots", "pdf", sep="."))
    pdf(file.path(outputDir, plotsFName))
    myExpres <- as.matrix(expres[selectedGenes, ])
    lev <- as.factor(groupingVar)
    varName <- "Gene Expression"
    if (is.null(groupNames)) groupNames <- unique(groupingVar)
    for (i in 1:nrow(myExpres)){
      desc <- paste("logFC=", round(top.Diff.common$logFC[i],3), 
                    ", p-val=", round(top.Diff.common$P.Value[i],6),
                    ", Adj-p=", round(top.Diff.common$adj.P.Val[i],6), sep="")
      beeswarm(myExpres[i,]~lev, 
               ylab="Expression", xlab=varDescript,
               main=paste(selectedGenes[i],desc, sep="\n"),
               labels=groupNames)  
      boxplot(myExpres[i,]~lev, add = T, names = c("",""), col="#0000ff22")
      # Segons un post de: https://www.r-statistics.com/2011/03/beeswarm-boxplot-and-plotting-it-with-r/
    }
    dev.off()
  }  
  return(list(genes=selectedGenes, 
              resT = resT.selected,
              topTable = top.Diff.all 
              )
         )
}



### Exemple

testQuickAnalysis <- function(){
  fName<- "https://raw.githubusercontent.com/alexsanchezpla/scripts/master/Exemple_Analisis_BioC/results/AvsB.Expres.csv"
  expresAndTopAvsB <- read.table (fName, head=T, sep=",", dec=".", row.names=1)
  expresAvsB <- expresAndTopAvsB [,5:14] 
  dim(expresAvsB)
  head(expresAvsB)
  groupVar <- c(rep("A", 5), rep("B",5))
  
  ### Crida amb càrrega dels paquets necessaris
  stopifnot(require(limma))
  stopifnot(require(multtest))
  stopifnot(require(hwriter))
  stopifnot(require(xlsx))
  stopifnot(require(beeswarm))

  AvsB <- quickAnalysis ( expres =expresAvsB,
                      groupingVar = groupVar ,
                      maxGenes=250,
                      aTitle="",
                      aFName="Results",
                      outputDir=".",
                      pvalThreshold=0.05,
                      useAdjP = TRUE,
                      plotSelected=TRUE)
}
