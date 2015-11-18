### 2.3 Funcio per crear la topTable anotada
### Funcions per a extreure identificadors sinonims d'una llista 

extractSinonims <- function(my.strings)
{
  my.sinonims <- list()
       if (runMulticore ==1 || runMulticore ==3) { 
            my.sinonims <- mclapply(my.strings, function(x) unlist(strsplit(x, " /// ")))

         } else {
            my.sinonims <- lapply(my.strings, function(x) unlist(strsplit(x, " /// ")))
       }

  return(my.sinonims)
}

midSinonims <- function(my.IDs)
{
  my.indexes <- grep(" /// ", my.IDs)
  my.sinonimIDs <- my.IDs[my.indexes]
    
  my.IDs[my.indexes] <- extractSinonims(my.sinonimIDs)
  return(my.IDs)
}


### knowledge2Names: Per a cada identificador dels gens contrueix un link cap al fitxer html que conte les anotacions
###
### Parametres
###
###    my.genes: vector de caracters amb els identificadors dels gens
###    knowledgeFile: nom del fitxer html que conte les anotacions dels gens
###
### Exemple
###
###    > myProbes <- c("101054_at", "101878_at", "102209_at")
###    > anotFileName <- "Annotations.html"
###    > ProbesLinked <- knowledge2gNames(myProbes, anotFileName)
###    > print(ProbesLinked)
###    [1] "<A HREF=\"Annotations.html#101054_at\" TARGET=\"_blank\">101054_at</A>"
###    [2] "<A HREF=\"Annotations.html#101878_at\" TARGET=\"_blank\">101878_at</A>"
###    [3] "<A HREF=\"Annotations.html#102209_at\" TARGET=\"_blank\">102209_at</A>"

knowledge2gNames <- function(my.genes, knowledgeFile)
{
  my.genes <- as.character(my.genes)

  my.annotation <- character()
  
  if (!is.null(knowledgeFile))
  {
    for (i in 1:length(my.genes))
    {
      my.annotation[i] <- paste("<A HREF=\"", knowledgeFile,"#", my.genes[i] ,"\" TARGET=\"_blank\">", my.genes[i],"</A>", sep = "")
    }
  }
  return(my.annotation)
}


annotateTopTable2 <- function (topTab,
                               fName,
                               Title = "Genes selected", 
                               anotPackage,
                               EntrezIDs = NULL,
                               SymbolIDs = NULL,
                               #comparison = "",
                               anotFilename = "annotations"
                               ) 
{
  #gNames <- as.character(topTab$ID)
  if (!is.null(topTab$ID)){
	gNames <- as.character(as.integer(topTab$ID))
  }else{
   	gNames <- as.character(rownames(topTab))
  }
  
  linkedGeneNanes <- knowledge2gNames(gNames, knowledgeFile =  paste(anotFilename, "html", sep="."))
  
  if (is.null(EntrezIDs))
  { 
    myenvirENTREZID <- eval(parse(text = paste(anotPackage, "ENTREZID", sep = "")))
    EntrezIDs <- unlist(mget(gNames, env = myenvirENTREZID, ifnotfound = NA))
  }else{
    EntrezIDs <- midSinonims(EntrezIDs[gNames])
  }
 
  if (is.null(SymbolIDs))
  { 
    myenvirSYMBOL <- eval(parse(text = paste(anotPackage, "SYMBOL", sep = "")))
    SymbolIDs <- unlist(mget(gNames, env = myenvirSYMBOL, ifnotfound = NA))
  }else{
    SymbolIDs <- midSinonims(SymbolIDs[gNames])
  }

  aux.SymbolIDs <- as.character(SymbolIDs)
  names(aux.SymbolIDs) <- names(SymbolIDs)
       if (runMulticore ==1 || runMulticore ==3) { 
           lS <- mclapply(aux.SymbolIDs, function(l){unlist(strsplit(l, "//"))})
         } else {
           lS <- lapply(aux.SymbolIDs, function(l){unlist(strsplit(l, "//"))})
       }
   
  linkedList <- list(en = EntrezIDs)
  topTab <- topTab[, -1]

                                        # Si existeixen sinonims "otherNames" ha de ser una llista i no un data.frame com estava
  otherNames = list(affyIDs = linkedGeneNanes, GeneSymbols = lS, topTab)  
  htmlpage(linkedList, filename = fName, title = Title, othernames = otherNames,
           table.head = c("EntrezID", "affyIDs", "GeneSymbols", names(topTab)), 
           table.center = TRUE, repository = list("en"), digits=4)
}

annotateTopTable2.Old <- function (topTab, fName, Title = "Genes selected", 
                               anotPackage, EntrezIDs = NULL, SymbolIDs = NULL, anotFileName = "Annotated.Genes.html") 
 {
    
  if (!is.null(topTab$ID)){
	gNames <- as.character(as.integer(topTab$ID))
  }else{
  	gNames <- as.character(rownames(topTab))
  }


    if (is.null(EntrezIDs)){ 
        stopifnot(require(old2db (anotPackage), character.only=T))
        myenvirENTREZID <- eval(parse(text = paste(anotPackage, "ENTREZID", sep = "")))
        EntrezIDs <- unlist(mget(gNames, env = myenvirENTREZID, ifnotfound = NA))
    }else{
        EntrezIDs <- EntrezIDs[gNames]
    }
 
    if (is.null(SymbolIDs)){ 
        stopifnot(require(old2db (anotPackage), character.only=T))
        myenvirSYMBOL <- eval(parse(text = paste(anotPackage, "SYMBOL",  sep = "")))
        SymbolIDs <- unlist(mget(gNames, env = myenvirSYMBOL, ifnotfound = NA))
    }else{
        SymbolIDs <- SymbolIDs[gNames]
    }
    linkedList <- list(en = EntrezIDs)
    topTab <- topTab[, -1]
    otherNames = cbind(affyIDs = gNames, GeneSymbols = SymbolIDs, topTab)
    htmlpage(linkedList, filename = fName, title = Title, othernames = otherNames, 
        table.head = c("EntrezID", names(otherNames)), table.center = TRUE, 
        repository = list("en"), digits=4)
}

