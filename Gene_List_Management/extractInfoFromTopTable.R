extractInfo <- function (x, compName, pattern, outDir, adjOrraw="adj", 
                         pCutOff=0.01, fcCutoff=1, UpDown="both"){
  geneList <- genesFromTopTable (x, entrezOnly = TRUE, uniqueIds=TRUE, #uniqueEntrezs, 
                                 adjOrrawP = adjOrraw, Pcutoff = pCutOff, FCcutoff = fcCutoff, 
                                 updown=UpDown, id2Select = "EntrezsA" , cols2Select =2)
  symbolsList <- as.character(genesFromTopTable (x, entrezOnly = TRUE, uniqueIds=TRUE, #uniqueEntrezs, 
                                                 adjOrrawP = adjOrraw, Pcutoff = pCutOff, FCcutoff = fcCutoff, 
                                                 id2Select = "SymbolsA" , cols2Select =1))
  universeList <- genesFromTopTable (x, entrezOnly = TRUE, uniqueIds=TRUE, #uniqueEntrezs, 
                                     adjOrrawP = "raw", Pcutoff = 1, FCcutoff = 0, 
                                     id2Select = "EntrezsA" , cols2Select =2)
  geneList <- as.character(geneList); symbolsList <- as.character(symbolsList); universeList <- as.character(universeList);
  
  write.table(geneList, file=file.path(outDir, paste(compName,"EntrezSelected.csv", sep=".")), 
              row.names=FALSE, quote=FALSE, col.names=FALSE)
  write.table(symbolsList, file=file.path(outDir, paste(compName,"SymbolsSelected.csv", sep=".")), 
              row.names=FALSE, quote=FALSE, col.names=FALSE)
  write.table(universeList, file=file.path(outDir, paste(compName,"EntrezAll.csv", sep=".")), 
              row.names=FALSE, quote=FALSE, col.names=FALSE)
  
  columns2select <- grep(pattern, colnames(x))
  write.csv(x[,columns2select], file=file.path(outDir, paste(compName,"Expres.csv", sep=".")), 
            row.names=TRUE, quote=FALSE)
  xbyEnt<- aggregate(x[,columns2select],by= list(x$EntrezsA), mean); 
  rownames(xbyEnt)<-xbyEnt[,1];xbyEnt<-xbyEnt[,-1]
  write.csv(xbyEnt,
            file=file.path(outDir, paste(compName,"ExpresByEntrez.csv", sep=".")), 
            row.names=TRUE, quote=FALSE)
  return(list(geneList, universeList))
}

