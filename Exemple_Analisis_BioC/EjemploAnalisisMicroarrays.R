## ----librerias, eval=TRUE------------------------------------------------
installifnot <- function (pckgName){
  if(!(require(pckgName, character.only=TRUE))){
    source("http://Bioconductor.org/biocLite.R")
    biocLite(pckgName)
  }
}
installifnot("Biobase")
installifnot("hgu133a.db")
installifnot("affy")
installifnot("affyPLM")
installifnot("arrayQualityMetrics")
installifnot("genefilter")
installifnot("limma")
installifnot("hgu133a.db")
installifnot("annotate")
installifnot("annaffy")
installifnot("hwriter")
installifnot("gplots")
installifnot("GOstats")

## ----preparaDirectorios, eval=TRUE---------------------------------------
workingDir <-getwd()
dataDir <-file.path(workingDir, "data")
resultsDir <- file.path(workingDir,"results")
celfilesDir <- file.path(workingDir,"celfiles")
setwd(workingDir)

## ----tomaMuestras, eval=FALSE--------------------------------------------
## muestras <- read.csv2(file.path(dataDir, "Asignacion_de_muestras_a_grupos.csv"),
##                       head=T)
## misMuestras <- as.character (muestras$Sample)
## paraAnalisis <- c(sample(misMuestras[1:6], 5),
##                   sample(misMuestras[7:20], 5),
##                   sample(misMuestras[23:48], 5))
## alAnalisis <-muestras[muestras$Sample %in% paraAnalisis,]
## write.table(alAnalisis, file=file.path(dataDir, "targets.txt"),
##             sep="\t", row.names=FALSE, quote=FALSE)

## ----phenoData1,echo=FALSE,print=FALSE,results=tex, eval=TRUE------------
require(Biobase)
my.targets <-read.AnnotatedDataFrame(file.path(dataDir,"targets.txt"), 
  header = TRUE, row.names = 1) 
stopifnot(require(xtable))
x.big<-xtable(pData(my.targets)[,1:4],
    label='table.targets',
    caption='Archivo targets.txt con la asignación a cada muestra de su condición experimental')
print(x.big,tabular.environment='longtable',floating=FALSE)

## ----affybatch.create, eval=TRUE-----------------------------------------
require(affy)
sampleInfo <- read.AnnotatedDataFrame(file.path(dataDir,"targets.txt"), 
    header = TRUE, row.names = 1, sep="\t")
fileNames <- rownames(pData(sampleInfo))
rawData <- read.affybatch(filenames=file.path(celfilesDir,fileNames),
                          phenoData=sampleInfo)
show(rawData)

## ----preajustes, echo=F, eval=TRUE---------------------------------------
colores <- c(rep("yellow", 5), rep("blue", 5), rep("red", 5))
grupos <- pData(rawData)$Group
numSamples <- nrow(pData(rawData))
sampleNames <-paste( pData(rawData)$SampleIDs, grupos, sep=".")
colnames(exprs(rawData))<-sampleNames

## ----plotHist,echo=F,fig=T, eval=TRUE------------------------------------
hist(rawData, main="Signal distribution", col=colores, lty=1:numSamples)
legend (x="topright", legend=sampleNames , col=colores, lty=1:numSamples)

## ----boxplot, fig=T, eval=TRUE-------------------------------------------
boxplot(rawData, cex.axis=0.6, col=colores, las=2, names=sampleNames, 
        main="Signal distribution for selected chips")

## ----plotPCA, results=H, eval=TRUE---------------------------------------
plotPCA <- function ( X, labels=NULL, colors=NULL, dataDesc="", scale=FALSE)
{
  pcX<-prcomp(t(X), scale=scale) # o prcomp(t(X))
  loads<- round(pcX$sdev^2/sum(pcX$sdev^2)*100,1)
  xlab<-c(paste("PC1",loads[1],"%"))
  ylab<-c(paste("PC2",loads[2],"%"))
  if (is.null(colors)) colors=1
  plot(pcX$x[,1:2],xlab=xlab,ylab=ylab, col=colors, 
       xlim=c(min(pcX$x[,1])-10, max(pcX$x[,1])+10))
  text(pcX$x[,1],pcX$x[,2], labels, pos=3, cex=0.8)
  title(paste("Plot of first 2 PCs for expressions in", dataDesc, sep=" "), cex=0.8)
}

## ----plotPCA2D, fig=T, eval=TRUE-----------------------------------------
plotPCA(exprs(rawData), labels=sampleNames, dataDesc="selected samples")

## ----distAnalisis, fig=T, eval=TRUE--------------------------------------
  manDist <-  dist(t(exprs(rawData))) 
  heatmap (as.matrix(manDist),  col=heat.colors(16))  

## ----plotDendro, echo=F,fig=T, eval=TRUE---------------------------------
clust.euclid.average <- hclust(dist(t(exprs(rawData))),method="average")
plot(clust.euclid.average, labels=sampleNames, main="Hierarchical clustering of samples",  hang=-1)

## ----computeDeg, echo=F, eval=TRUE---------------------------------------
deg<-AffyRNAdeg(rawData, log.it=T)

## ----plotDeg,echo=F, fig=T, eval=TRUE------------------------------------
plotAffyRNAdeg(deg) 
legend (x="bottomright", legend=sampleNames, col=colores,  lty=1:numSamples, cex=0.7) 

## ----affyPLM, eval=TRUE--------------------------------------------------
stopifnot(require(affyPLM))
if (!(exists("Pset"))) 
  Pset<- fitPLM(rawData)

## ----plotRLE, echo=F, fig=T, eval=TRUE-----------------------------------
RLE(Pset, main = "Relative Log Expression", 
    names=sampleNames, las=2, col=colores, cex.axis=0.7, ylim=c(-7,7))

## ----plotNUSE, echo=F, fig=T, eval=TRUE----------------------------------
NUSE(Pset, main = "Normalized Unscaled Standard Errors", las=2, 
     names=sampleNames, las=2, col=colores, cex.axis=0.7,  ylim=c(0.5,2.5))

## ----arrayQuality, eval=TRUE---------------------------------------------
stopifnot(require(arrayQualityMetrics))
arrayQualityMetrics(rawData, intgroup="Group",
                    outdir = file.path(resultsDir, "arrayQuality"), 
                    force=TRUE)

## ----normalization.rma, eval=TRUE----------------------------------------
stopifnot(require(affy))
eset_rma <- rma(rawData)    
eset_rma
## ----normBoxPlot, fig=T, eval=TRUE---------------------------------------
boxplot(eset_rma, main="RMA", names=sampleNames, 
        cex.axis=0.7, col=colores,las=2)

## ----filtraje------------------------------------------------------------
require(genefilter)
filtered <- nsFilter(eset_rma, require.entrez=TRUE,
         remove.dupEntrez=TRUE, var.func=IQR,
         var.cutoff=0.5, var.filter=TRUE,
         filterByQuantile=TRUE, feature.exclude="^AFFX")

## ----filtrado------------------------------------------------------------
names(filtered)
class(filtered$eset)
print(filtered$filter.log)
eset_filtered <-filtered$eset

## ----saveData------------------------------------------------------------
save(eset_rma, eset_filtered, file=file.path(resultsDir, "datos.normalizados.Rda"))

## ----writeNormalized-----------------------------------------------------
write.csv2(exprs(eset_rma), file.path(resultsDir, "Datos.Normalizados.csv2"))

## ----matDesign, eval=TRUE------------------------------------------------
design<-matrix(
            c(1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,1,1,1,1,1),
            nrow=15,byrow=F)
colnames(design)<-c("A", "B", "L")
rownames(design) <-  sampleNames 
print(design)

## ----setContrasts, eval=TRUE---------------------------------------------
require(limma)
cont.matrix <- makeContrasts (
      AvsB = B-A,
      AvsL = L-A,
      BvsL = L-B,
      levels=design)
print(cont.matrix)

## ----linearmodelfit,echo=F, eval=TRUE------------------------------------
require(limma)
fit<-lmFit(eset_rma, design)
fit.main<-contrasts.fit(fit, cont.matrix)
fit.main<-eBayes(fit.main)

## ----print=FALSE, echo=TRUE, eval=TRUE-----------------------------------
topTab_AvsB <- topTable (fit.main, number=nrow(fit.main), coef="AvsB", adjust="fdr")
topTab_AvsL <- topTable (fit.main, number=nrow(fit.main), coef="AvsL", adjust="fdr")
topTab_BvsL  <- topTable (fit.main, number=nrow(fit.main) , coef="BvsL", adjust="fdr")

## ----volcanos, results=tex,echo=FALSE, eval=TRUE-------------------------
volcanoplot(fit.main, coef="AvsB", highlight=10, names=fit.main$ID)
for(i in 1:ncol(cont.matrix)){
  compName <-colnames(cont.matrix)[i]
  file=paste("volcanoPlot", compName, ".pdf", sep="")
  pdf(file=file.path(workingDir, "images", file), paper="special", width=6, height=6)
  volcanoplot(fit.main, coef=i, highlight=10, names=fit.main$ID, 
            main=paste("Differentially expressed genes",compName, sep="\n"))
  abline(v=c(-1,1))
  dev.off()
  cat("\\includegraphics{", file, "}\n\n", sep="")
}

## ----CuantosGenes, echo=F, eval=TRUE-------------------------------------
cat("Numero de genes con un p--valor inferior a 0.05 en cada comparación:\n")
  cat ("En la comparación 'A vs B': ", sum(topTab_AvsB$adj.P.Val<=0.05),"\n")
  cat ("En la comparación 'A vs L': ", sum(topTab_AvsL$adj.P.Val<=0.05),"\n")
  cat ("En la comparación 'B vs L': ", sum(topTab_BvsL$adj.P.Val<=0.05),"\n")  

cat("\nNumero de genes con un p--valor inferior a 0.01 en cada comparación:\n")
  cat ("En la comparación 'A vs B': ", sum(topTab_AvsB$adj.P.Val<=0.01),"\n")
  cat ("En la comparación 'A vs L': ", sum(topTab_AvsL$adj.P.Val<=0.01),"\n")
  cat ("En la comparación 'B vs L': ", sum(topTab_BvsL$adj.P.Val<=0.01),"\n")  

## ----topGenesAvsB,echo=FALSE,print=FALSE,results=tex, eval=TRUE----------
require(Biobase)
stopifnot(require(xtable))
AvsB10<-xtable(topTab_AvsB[1:10,1:7],
    label='topTab_AvsB',
    caption='10 genes más expresados diferencialmente en la comparación A vs B')
print(AvsB10, tabular.environment='longtable',floating=FALSE)

## ----topGenesAvsL,echo=FALSE,print=FALSE,results=tex, eval=TRUE----------
require(Biobase)
stopifnot(require(xtable))
AvsL10<-xtable(topTab_AvsL[1:10,1:7],
    label='topTab_AvsL',
    caption='10 genes más expresados diferencialmente en la comparación A vs L')
print(AvsL10, tabular.environment='longtable',floating=FALSE)

## ----topGenesBvsL,echo=FALSE,print=FALSE,results=tex, eval=TRUE----------
BvsL10<-xtable(topTab_BvsL[1:10,1:7],
    label='topTab_BvsL',
    caption='10 genes más expresados diferencialmente en la comparación B vs L')
print(BvsL10, tabular.environment='longtable',floating=FALSE)

## ----decideTests.1, echo=F, eval=TRUE------------------------------------
res<-decideTests(fit.main, method="separate", adjust.method="fdr", p.value=0.01, lfc=1)

## ----resumeDecideTests, eval=TRUE----------------------------------------
sum.res.rows<-apply(abs(res),1,sum)
res.selected<-res[sum.res.rows!=0,] 
print(summary(res))

## ----decideTests.2, echo=F, eval=TRUE------------------------------------
res<-decideTests(fit.main, method="separate", adjust.method="fdr", p.value=0.05, lfc=1)
sum.res.rows<-apply(abs(res),1,sum)
res.selected<-res[sum.res.rows!=0,] 
print(summary(res))

## ----venn1,fig=T, eval=TRUE----------------------------------------------
vennDiagram (res.selected[,1:3], main="Genes in common #1", cex=0.9)

## ----anota1, eval=TRUE---------------------------------------------------
require(hgu133a.db)
hgu133a()

## ----annaffy-------------------------------------------------------------
require(annaffy)
genesSelected <- rownames(res.selected)
at <- aafTableAnn(genesSelected, "hgu133a.db")
saveHTML (at, file.path(resultsDir, "anotations.html"), 
          "Annotations for selected genes")

## ----htmlPages-----------------------------------------------------------
require(annotate)
require(hgu133a.db)
listOfTables <- list(AvsB = topTab_AvsB, AvsL = topTab_AvsL, BvsL = topTab_BvsL) 
for (i in 1:length(listOfTables)){
  # Seleccionamos la "topTable"
  topTab <- listOfTables[[i]]
  # Escogemos los grupos de sondas a incluir en la tabla
  whichGenes<-topTab["P.Value"]<0.05
  selectedIDs <- rownames(topTab)[whichGenes]
  # Los convertimos a identificadores Entrez ("EG") y a Gene Symbols
  genes<- getEG(selectedIDs, "hgu133a.db")
  simbols <-getSYMBOL(selectedIDs, "hgu133a.db")
  # Haremos la columna de Entrez sea hiperenlazable
  paraEnlace <- list (misgenes=genes)
  # Preparamos el data.frame con el que se creará el archivo de resultados
  otherNames = data.frame(selectedIDs, simbols, topTab[whichGenes,-1])
  names(otherNames) = c("Affy ID", "Gene Symbol", colnames(topTab)[-1])
  # Invocamos la función "htmlpage"
  comparison <- names(listOfTables)[i]
  htmlpage(paraEnlace, 
           filename =file.path(resultsDir, 
           paste("Selected Genes in comparison ",comparison,".html", sep="")) , 
           title = paste("Diff. expressed genes in comparison ", comparison, sep=""), 
           othernames = otherNames, 
           table.head = c("Entrez IDs", names(otherNames)),
           table.center = TRUE, 
           repository=list("en"))
}

## ----prepareData, eval=TRUE----------------------------------------------
probeNames<-rownames(res)
probeNames.selected<-probeNames[sum.res.rows!=0]
exprs2cluster <-exprs(eset_rma)[probeNames.selected,]
colnames(exprs2cluster)<-sampleNames
color.map <- function(grupo) { 
  if (grupo=="A"){
    c<- "yellow" 
  }else{ 
    if (grupo=="B"){
      c<- "red"
    }else{
      c<- "blue"
   }
  }
return(c)}

## ----plotHeatMap1, fig=T, eval=TRUE--------------------------------------
grupColors <- unlist(lapply(pData(eset_rma)$Group, color.map))
heatmap(exprs2cluster, col=rainbow(100), ColSideColors=grupColors, cexCol=0.9)

## ----plotHeatMap2, fig=T, eval=TRUE--------------------------------------
grupColors <- unlist(lapply(pData(eset_rma)$Group, color.map))
require("gplots")
heatmap.2(exprs2cluster, 
          col=bluered(75), scale="row",
          ColSideColors=grupColors, key=TRUE, symkey=FALSE, 
          density.info="none", trace="none", cexCol=1)

## ----GOAnalysis----------------------------------------------------------
require(GOstats)
listOfTables <- list(AvsB = topTab_AvsB, AvsL = topTab_AvsL, BvsL = topTab_BvsL) 
for (i in 1:length(listOfTables)){
  # Seleccionamos la "topTable"
  topTab <- listOfTables[[i]]
  # Definimos el universo de genes: todos los que se han incluido en el análisis
  # EL programa trabaja con identificadores "entrez" y no admite duplicados
  
  entrezUniverse = unique(getEG(as.character(topTab$ID), "hgu133a.db"))
  
  # Escogemos los grupos de sondas a incluir en el análisis
  # Este análisis trabaja bien con varios centenares de genes 
  # por lo que es habitual basarse en p-valores sin ajustar para incluirlos
  
  whichGenes<-topTab["adj.P.Val"]<0.05
  geneIds <-   unique(getEG(as.character(topTab$ID[whichGenes]),"hgu133a.db"))
  
  # Creamos los "hiperparámetros" en que se basa el análisis
  GOparams = new("GOHyperGParams",
    geneIds=geneIds, universeGeneIds=entrezUniverse,
    annotation="org.Hs.eg.db", ontology="BP",
    pvalueCutoff=0.001, conditional=FALSE,
    testDirection="over")
  KEGGparams = new("KEGGHyperGParams",
    geneIds=geneIds, universeGeneIds=entrezUniverse,
    annotation="org.Hs.eg.db",  
    pvalueCutoff=0.01, testDirection="over")

  # Ejecutamos los análisis

  GOhyper = hyperGTest(GOparams)
  KEGGhyper = hyperGTest(KEGGparams)
  
# Creamos un informe html con los resultados
   comparison = names(listOfTables)[i]
   GOfilename =file.path(resultsDir, 
     paste("GOResults.",comparison,".html", sep=""))
   KEGGfilename =file.path(resultsDir, 
     paste("KEGGResults.",comparison,".html", sep=""))
  htmlReport(GOhyper, file = GOfilename, summary.args=list("htmlLinks"=TRUE))
  htmlReport(KEGGhyper, file = KEGGfilename, summary.args=list("htmlLinks"=TRUE))
}

## ----listaArchivos,echo=FALSE,print=FALSE,results=tex, eval=TRUE---------
require(gdata)
listaArchivos <-read.table(file.path(resultsDir,"listaArchivos.txt"), head=TRUE, sep="\t") 
stopifnot(require(xtable))
x.big<-xtable(listaArchivos,
    label='listaArchivos',
    caption='Lista de archivos generados en este análisis')
print(x.big,tabular.environment='longtable',floating=FALSE)

## ----listaArchivos2html,echo=FALSE,print=FALSE, eval=TRUE----------------
require(hwriter)
hwrite(listaArchivos,file.path(resultsDir, "listaArchivos.html"))

