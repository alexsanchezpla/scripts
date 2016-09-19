analyzeData <- function ( expres, groupingVar, maxGenes=250, 
                          aTitle="", aFName="Results", outputDir=".", 
                          pvalThreshold=0.05)
{
  print(paste("ANALYIS:",aTitle,"...PROCESS STARTED",sep=" "))
  ### Packages needed in both analysis
  stopifnot(require(limma))
  stopifnot(require(multtest))
  stopifnot(require(hwriter))
  
  ### ANALYSIS USING multtest library
  
  resT <- mt.maxT(expres, classlabel = as.factor(groupingVar))
  resT.selected <-resT[resT$rawp < pvalThreshold,]
  gNames.multtest <-rownames(resT.selected)
  aTitle1 <- paste(aTitle,"\n","method: multtest",sep="")
  aFName1 <- file.path(outputDir,paste(aFName, "multtest","html", sep="."))
  hwrite (resT.selected, page =aFName1)
  
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
  top.Diff <-top.Diff.all[top.Diff.all$P.Value < pvalThreshold,]
  gNames.limma <-rownames(top.Diff)
  
  aTitle1<-paste(aTitle,"method: limma",sep="\n")
  aFName1<- file.path(outputDir,paste(aFName, "limma","html", sep="."))
  hwrite (top.Diff, page= aFName1)
  
  selectedGenes <-union(gNames.multtest,gNames.limma)
  
  print (paste("Number of genes selected using permutations test (p <0.05): ", 
               nrow(resT.selected),sep=""))
  print (paste("Number of genes selected using linear model test (p <0.05): ",
               nrow(top.Diff),sep=""))
  print (paste ("Number of genes selected by both test (p <0.05)          : ",
                length(intersect(gNames.multtest,gNames.limma)),sep=""))
  
  print(paste("ANALYIS:",aTitle,"...PROCESS COMPLETED",sep=" "))

  return(list(genes=selectedGenes, 
              resT = resT.selected,
              topTable = top.Diff.all 
              )
         )
}

### Exemple

fName<- "https://raw.githubusercontent.com/alexsanchezpla/scripts/master/Exemple_Analisis_BioC/results/AvsB.Expres.csv"
expresAndTopAvsB <- read.table (fName, head=T, sep=",", dec=".", row.names=1)
expresAvsB <- expresAndTopAvsB [,5:14] 
dim(expresAvsB)
head(expresAvsB)
groupVar <- c(rep("A", 5), rep("B",5))

# Assignació de valors als paràmetres
expres =expresAvsB 
groupingVar = groupVar 
maxGenes=250
aTitle=""
aFName="Results"
outputDir="."
pvalThreshold=0.05
