##############################################
### Funcions per a crear l'arxiu de'enllaços
##############################################

addToLinksFile <- function(linksFile, aFName, categ = "", desc = "", subcateg = "")
{
  if (!is.null(linksFile))
  {
    write(paste(aFName, categ, subcateg, desc, sep = "\t"), file = linksFile, append = TRUE)
  }
}
 
createLinksFile <- function (linksFile)
 {
   linia <- paste("FileName", "Category", "Subcategory", "Description", sep = "\t")
   write(linia, file = linksFile, append = FALSE)
 }
 

FixLinksFile <- function(lFile)
{
  my.LnkFile <- read.table(lFile, header = TRUE, sep = "\t")

  my.newFile <- unique(my.LnkFile)
  
  write.table(my.newFile, lFile, sep = "\t", row.names = FALSE, quote = FALSE)
}


######################
### Exemples d'us
######################
###
### addToLinksFile(linksFileName,
###                MainReportName,
###                categ = 'INFO', 
###                desc = "Report with description of the main results")
### 
###


#####################################################################################################################
### Funcions per a escriure el fitxer de Results_Files
#####################################################################################################################

### printHeader: Crea la cap??alera del fitxer html de resultats
###                 
###
### Parametres:
###
###   FileName : Nom del fitxer de resultats
###
### OBERVACIONS IMPORTANTS:
###
###   

printHeader <- function(FileName = "ResultFiles")
{
  outfile <- file(paste(FileName, ".html", sep = ""), open = "wt")
  
  txt <- "<html> \n"
  txt <- paste(txt, "<head> \n", sep = "")
  txt <- paste(txt, "<meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\"/>", sep = "")
  txt <- paste(txt, "<title>List of results files</Title>", sep = "")
  txt <- paste(txt, "</head> \n", sep = "")
  txt <- paste(txt, "<style> \n", sep = "")
  txt <- paste(txt, "  <!--td {font-family: Helvetica,Arial; font-size: 14px;}--> \n", sep = "")
  txt <- paste(txt, "  <!--h1 {font-family: Helvetica,Arial; font-size: 22px;}--> \n", sep = "")
  txt <- paste(txt, "</style> \n", sep = "")
  txt <- paste(txt, "<body bgcolor=\"White\" text=\"Black\"> \n", sep = "")
  
  cat(txt, file = outfile, sep = "")
  
  close(outfile)
}


### printGroupHeader: Construeix la cap??alera de presentacio que es veura a la pagina
###                 
###
### Parametres:
###
###   FileName : Nom del fitxer de resultats
###        UEB : Per defecte pren valor TRUE. Si es FALSE crea la cap??alarea adaptada a la UB
###
### OBERVACIONS IMPORTANTS:
###
###   

printGroupHeader <- function(FileName = "ResultFiles", UEB = TRUE)
{
  outfile <- file(paste(FileName, ".html", sep = ""), open = "at")

             txt <- "<table width=\"100%\"  border=\"0\"> \n"
  txt <- paste(txt, "   <tr><td width=\"24%\"> \n", sep = "")
  txt <- paste(txt, "	    <div align=\"center\"> \n", sep = "")
  
  if(UEB)
  {
    txt <- paste(txt, "		 <a href=\"http://www.vhir.org\" target=\"z\"> \n", sep = "")
    txt <- paste(txt, "			<img src=\"images/IR.jpg\" width=\"195\" height=\"73\" border=\"0\"> \n", sep= "")
    txt <- paste(txt, "		 </a> \n", sep= "")
    txt <- paste(txt, "	    </div> \n", sep= "")
    txt <- paste(txt, "	</td> \n", sep= "") 
    txt <- paste(txt, "	<td width=\"53%\"> \n", sep= "")
    txt <- paste(txt, "	    <h1 align=\"center\">Unitat d'Estad&iacute;stica i Bioinform&agrave;tica</h1> \n", sep= "")
    txt <- paste(txt, "	    <div align=\"center\">Vall d'Hebron Institut de Recerca</div> \n", sep= "")
    txt <- paste(txt, "       </td> \n", sep= "") 
    txt <- paste(txt, "	<td width=\"23%\"> \n", sep= "")
    txt <- paste(txt, "	    <div align=\"center\"> \n", sep= "")
    txt <- paste(txt, "		 <a href=\"http://ueb.vhir.org\" target=\"z\"> \n", sep= "")
    txt <- paste(txt, "			<img src=\"images/UEBblanc.jpg\" width=\"204\" height=\"48\" border=\"0\"> \n", sep= "")
  }else{

    txt <- paste(txt, "		 <a href=\"http://www.ub.edu\" target=\"z\"> \n", sep = "")
    txt <- paste(txt, "			<img src=\"images/UB.jpg\" width=\"195\" height=\"73\" border=\"0\"> \n", sep= "")
    txt <- paste(txt, "		 </a> \n", sep= "")
    txt <- paste(txt, "	    </div> \n", sep= "")
    txt <- paste(txt, "	</td> \n", sep= "") 
    txt <- paste(txt, "	<td width=\"53%\"> \n", sep= "")
    txt <- paste(txt, "	    <h1 align=\"center\">Grup de Recerca en Estad&iacute;stica i Bioinform&agrave;tica</h1> \n", sep= "")
    txt <- paste(txt, "	    <div align=\"center\">Departament d'Estad&iacute;stica - Universitat de Barcelona</div> \n", sep= "")
    txt <- paste(txt, "       </td> \n", sep= "") 
    txt <- paste(txt, "	<td width=\"23%\"> \n", sep= "")
    txt <- paste(txt, "	    <div align=\"center\"> \n", sep= "")
    txt <- paste(txt, "		 <a href=\"http://estbioinfo.stat.ub.es\" target=\"z\"> \n", sep = "")
    txt <- paste(txt, "			<img src=\"images/Logo%20EstBioinfo.jpg\" width=\"204\" height=\"48\" border=\"0\"> \n", sep= "")
  }

  txt <- paste(txt, "		 </a> \n", sep= "")
  txt <- paste(txt, "	    </div> \n", sep= "")
  txt <- paste(txt, "	</td>  \n", sep= "")
  txt <- paste(txt, "   </tr> \n", sep= "")
  txt <- paste(txt, "</table> \n", sep= "")

  cat(txt, file = outfile, sep = "")
  
  close(outfile)
}


### printAnalysisDetails: Construeix la taula que conte la informacio general de l'estudi, com dades de contacte, 
###                       titol de l'analisi
###
### Parametres
###
###    FileName : Nom del fitxer de resultats
###   Info.list : Llista amb la informacio general 
###
### Exemple
###
###    myInfolist <- list(To = "Nom i cognoms de l'usuari",
###                       Description = "Diferentially expressed genes associated with...",
###                       Analysts = "Nom del primer analista and Alex Sanchez",
###                       Contact = "Alex Sanchez (alesanchez@ir.vhebron.net)")
###    printAnalysisDetails(FileName = "ResultFiles", Info.list =  myInfoList)
###
###    # Aixo creara un fitxer que contindra una taula de l'estil
###    <table width="100%" border="0"> 
###        <tr><td height="35" colspan="2" valign="center" bgcolor="#d0d0f0"> 
###                <h1>Statistical Analysis Report</h1>
###             </td>
###        </tr> 
###        <tr height="35"><td width="15%" align="right" bgcolor="#e0e0f0"><b>To</b></td> 
###                        <td width="85%" bgcolor="#f0f0ff"><i>Nom i cognoms de l'usuari</i></td> 
###        </tr> 
###        <tr height="35"><td align="right" bgcolor="#e0e0f0"><b>Description</b></td> 
###                        <td  bgcolor="#f0f0ff"><i>Diferentially expressed genes associated with...</i></td> 
###	                   <br> 
###	   </tr> 
###	   <tr height="35"><td align="right" bgcolor="#e0e0f0"><b>Analysts</b></td> 
###	                   <td  bgcolor="#f0f0ff"><i>Nom del primer analista and Alex Sanchez</i></td>
###                        <br> 
###	   </tr> 
###	   <tr height="35"><td align="right" bgcolor="#e0e0f0"><b>Contact</b></td> 
###	                   <td bgcolor="#f0f0ff"><i>Alex Sanchez (alesanchez@ir.vhebron.net)</i></td> 
###	      <br> 
###	   </tr> 
###    </table> 
###    <hr> 

printAnalysisDetails <- function(FileName = "ResultFiles", Info.list)
{
  outfile <- file(paste(FileName, ".html", sep = ""), open = "at")

	     txt <- "<table width=\"100%\" border=\"0\"> \n"  
  txt <- paste(txt, "     <tr><td height=\"35\" colspan=\"2\" valign=\"center\" bgcolor=\"#bb66aa\"> \n", sep = "")
  txt <- paste(txt, "                <h1>Statistical Analysis Report</h1></td></tr> \n", sep = "")
  txt <- paste(txt, "     <tr height=\"35\"><td width=\"15%\" align=\"right\" bgcolor=\"#e0e0f0\"><b>To</b></td> \n", sep = "")
  txt <- paste(txt, "         <td width=\"85%\" bgcolor=\"#f0f0ff\"><i>", Info.list$To, "</i></td> \n", sep = "")
  txt <- paste(txt, "     </tr> \n", sep = "")
  txt <- paste(txt, "     <tr height=\"35\"><td align=\"right\" bgcolor=\"#e0e0f0\"><b>Description</b></td> \n", sep = "")
  txt <- paste(txt, "         <td  bgcolor=\"#f0f0ff\"><i>", Info.list$Description, "</i></td> \n", sep = "")
  txt <- paste(txt, "	      <br> \n", sep = "")
  txt <- paste(txt, "	  </tr> \n", sep = "")
  txt <- paste(txt, "	  <tr height=\"35\"><td align=\"right\" bgcolor=\"#e0e0f0\"><b>Analysts</b></td> \n", sep = "")
  txt <- paste(txt, "	      <td  bgcolor=\"#f0f0ff\"><i>", Info.list$Analysts, "</i></td> \n", sep = "")
  txt <- paste(txt, "	      <br> \n", sep = "")
  txt <- paste(txt, "	  </tr> \n", sep = "")
  txt <- paste(txt, "	  <tr height=\"35\"><td align=\"right\" bgcolor=\"#e0e0f0\"><b>Contact</b></td> \n", sep = "")
  txt <- paste(txt, "	      <td bgcolor=\"#f0f0ff\"><i>", Info.list$Contact, "</i></td> \n", sep = "")
  txt <- paste(txt, "	      <br> \n", sep = "")
  txt <- paste(txt, "	  </tr> \n", sep = "")
  txt <- paste(txt, "</table> \n", sep = "")
  txt <- paste(txt, "<hr> \n", sep = "")
  
  cat(txt, file = outfile, sep = "")

  close(outfile)
}  


### printVHIRfooter: Inserta la imatge amb logo del VHIR al peu de la pagina de resultats
###                 
###
### Parametres:
###
###   FileName : Nom del fitxer de resultats
###
### OBERVACIONS IMPORTANTS:
###   Si UEB es FALSE no hauria d'entrar a fer aixo
###   

printVHIRfooter <- function(FileName = "ResultFiles")
{
  outfile <- file(paste(FileName, ".html", sep = ""), open = "at")

  txt <- "<table width=\"100%\"  border=\"0\">\n"
  txt <- paste(txt, "   <tr height=\"50\"><td width=\"100%\"> </td></tr> \n", sep="")
  txt <- paste(txt, "   <tr><td width=\"100%\"><div align=\"center\"> \n", sep="")
  txt <- paste(txt, "      <a href=\"http://www.vhir.org\" target=\"z\"><img width=\"90%\" src=\"images/imatgelogotip.jpg\" border=\"0\"></a> \n", sep="")
  txt <- paste(txt, "   </div></td></tr> \n", sep="")

  txt <- paste(txt, "</table> \n", sep= "")

  cat(txt, file = outfile, sep = "")
  
  close(outfile)
}


### write.ResultFiles: 
###                 
###
### Parametres:
###
###    FileName :
###   Info.list :
###
### OBERVACIONS IMPORTANTS:
###
###

write.section <- function(my.info,
                          filename,
                          my.id = "",
                          sectionTitle = "", # Era 'Files' que s'afegia darrera els titol de seccio
                          IndexDir = "") 
{
  outfile <- file(paste(filename, ".html", sep = ""), open = "at")


  txt <- paste("<h1 id=\"", my.id, "\">", sectionTitle, "</h1> \n", sep = "")

  txt <- paste(txt, "<table border=\"0\"> \n", sep = "")

  txt <- paste(txt, "    <tr>" , sep = "")
  for (j in 1:ncol(my.info))
  {
    txt <- paste(txt, "<td bgcolor=\"", c("#cc99cc", "#cc77cc")[j%%2 + 1], "\"><b>", colnames(my.info)[j], "</b></td>\n", sep = "")
  }
  txt <- paste(txt, "    </tr> \n", sep = "")

  for (i in 1:nrow(my.info))
  {
    txt <- paste(txt, "    <tr>" , sep = "")
    for (j in 1:ncol(my.info))
    {
      my.cell <- my.info[i, j]
      if(j == 1)
      {
        my.cell <- paste("<a href=\"", paste(IndexDir, my.info[i, j], sep=""),"\" target=\"z\">", my.cell, "</a>", sep = "")
      }
            
      txt <- paste(txt,
                   "<td bgcolor=\"", c("#e0e0ff", "#d0d0f0", "#f0f0ff", "#e0e0f0")[i%%2 * 2 + j%%2 + 1], "\">",
                   my.cell,
                   "</td> \n",
                   sep = "")
    }
    txt <- paste(txt, "    </tr> \n" , sep = "")
  }

  txt <- paste(txt, "</table> \n", sep = "")

  cat(txt, file = outfile, sep = "")

  close(outfile)
}


### printAnalysisDetails: 
###                 
###
### Parametres:
###
###    FileName :
###
### OBERVACIONS IMPORTANTS:
###
###

closeHtml <- function(FileName = "ResultFiles")
{
  outfile <- file(paste(FileName, ".html", sep = ""), open = "at")
  
  cat("</body>", "</html>", file = outfile, sep = "\n")

  close(outfile)
}

### SOME POSSIBLE GROPS OF CATEGORIES TO STRUCTURE THE OUTPUT FILE
### GEN (for 'General')
### Use this categories to prepare a results browser for a general analysis

#' myGENCategs
#'
#' This function returns names and descriptions for a general results file
#' @keywords categories
#' @export myGENCategs
#' @examples 
#' myCategs <-myGENCategs()
#' categorias <- myCategs[["categs"]]
#' nombres <-  myCategs[["noms"]]
myGENCategs <- function(){
  GENCategs <- c("Reports and result summaries",
                 "Data and other input",
                 "Quality Control",
                 "Pre-processing",
                 "Analysis",
                 "Annotations",
                 "Multiple Comparisons",
                 "Cluster",
                 "GO Analysis",
                 "KEGG Analysis",
                 "Ingenuity Pathways Analysis",
                 "Power Analysis")
  GENNames <- c("INFO", "DATA", "QC", "NORM", "ANALYSIS",
                "ANNOT", "MULTCOMP", "CLUSTER", "GO", "KEGG", "IPA", "POWER")
  return(list(categs=GENCategs, noms=GENNames))
}

#' myFACategs
#'
#' This function returns names and descriptions for a 'Functional Analysis' results file
#' @keywords categories
#' @export myFACategs
#' @examples 
#' myCategs <-myFACategs()
#' categorias <- myCategs[["categs"]]
#' nombres <-  myCategs[["noms"]]
myFACategs <- function(){
  FACategs <- c("Reports and result summaries",
                "Input Files for Biological Significance Analysis",
                "Gene Set Expression Analysis (GAGE) Results",
                "DAVID Analysis Results",
                "GOProfiles Analysis Results",
                "GO Analysis",
                "KEGG Analysis",
                "Ingenuity Pathways Analysis")
  FANames <- c("INFO", "INPUT", "GAGE", "DAVID", "GOPROF", "GOEA", "KEGGEA", "IPA")
  return(list(categs=FACategs, noms=FANames))
}

################


#' A Function to create HTML from file of links
#'
#' This function allows you to create an HTML from a text file with a predefined structure
#' @param lFile Input file consisting of 4 tab separated fields: Filename, Categ, Subcateg, Description
#' @param outputDir Name of the directory where the HTML is to be written
#' @param info.list List with informations to build the header: To, Description, Analysts, Contact
#' @param categs.descs Description of categories in the links file (Header of sections)
#' @param categs.names Short names for categories in the links file (to be used in the input file)
#' @param IndexDir Name of subdirectory where the files to be linked are located
#' @param UEB TRUE if UEB header is to be used
#' @param resultsFileName Name of output file
#' @keywords HTML, links
#' @export LinksFile2Html
#' @examples 
#' linksFileName <-"LinksFileName_Example.csv"
#' resultsDir <- "."
#' # UEB <- TRUE # if UEB == TRUE => UEB header is used
#' htmlInfo <- list(To = "My Client", # nomClients,
#'                  Description = "Analysis of Biological Significance for some selected gene list",
#'                  Analysts = "Ferran Brianso and Alex Sanchez",               # analistes,
#'                  Contact = "Alex Sanchez (alex.sanchez@vhir.org)")
#' GENCategs <- myGENCategs() 
#' fcategs <- GENCategs[["categs"]]
#' fNames <-  GENCategs[["noms"]]
#' # Invoke main function 
#' # Not run automatically. Check by yourself uncommenting next lines
#' # setwd("examples")
#' # LinksFile2Html(linksFileName, resultsDir, htmlInfo, 
#' #                categs.descs = fCategs, categs.names = fNames,
#' #                IndexDir = "", UEB = TRUE, resultsFileName="Resultats")

LinksFile2Html <- function(lFile,
                           outputDir,
                           info.list,
                           categs.descs=NULL,
                           categs.names=NULL,
                           IndexDir = "",
                           UEB = TRUE,
                           resultsFileName="ResultFiles")
{
  FixLinksFile(lFile)

  if (is.null(categs.descs)){
      my.names <- c("Result Summary",
                "Data",
                "Quality Control",
                "Pre-processing",
                "Analysis",
                "Annotations",
                "Multiple Comparisons",
                "Cluster",
                "GO Analysis",
                "KEGG Analysis",
                "Ingenuity Pathways Analysis",
                "Power Analysis")
  }else{
      my.names<-categs.descs
  }
  if (is.null(categs.names)){
           names(my.names) <- c("INFO", "DATA", "QC", "NORM", "ANALYSIS",
                                "ANNOT", "MULTCOMP", "CLUSTER", "GO", "KEGG", "IPA", "POWER")
  }else{
      names(my.names)<-categs.names
  }
  
  my.LnkFile <- read.table(lFile, header = TRUE, sep = "\t", stringsAsFactor = FALSE)
  my.cats <- unique(my.LnkFile[, 2])
  
  resFile <- file.path(outputDir, resultsFileName)
  
  printHeader(resFile)
  printGroupHeader(resFile, UEB = UEB)
  printAnalysisDetails(resFile, info.list)
  
  i <- 0 # Per generar el numero de categoria en el fitxer resultFiles.html de forma dinamica
  for(categ in names(my.names)[names(my.names) %in% my.cats])
  {
    i <- i+1
    my.ind <- which(my.LnkFile[, 2] == categ)
    write.section(my.info = my.LnkFile[my.ind, c(1,4)],  
                  filename = resFile,
                  my.id = categ,                  
                  sectionTitle = paste(i, ". ", my.names[categ], sep = ""),                  
                  IndexDir = IndexDir)
  }

  if (UEB) printVHIRfooter(resFile)

  closeHtml(resFile)
}

