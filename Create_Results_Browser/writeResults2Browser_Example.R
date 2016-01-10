###
### Example of using function LinksFile2Html
###


source ("writeResults2Browser.R") # To be replaced soon by 'library(writeResults2Browser.R")'

# Set values for input function

linksFileName <-"LinksFileName_Example.csv"
resultsDir <- "."
UEB <- TRUE # if UEB == TRUE => UEB header is used

htmlInfo <- list(To = "My Client", # nomClients,
                 Description = "Analysis of Biological Significance for some selected gene list",
                 Analysts = "Ferran Brianso and Alex Sanchez",               # analistes,
                 Contact = "Alex Sanchez (alex.sanchez@vhir.org)"
)

fCategs <- myGENCategs # General Categories Descriptions defined in the package. Type myGENcategs to see them
fNames <- myGENNames   # General Categories Names defined in the package. Type myGENcategs to see them

# Invoke main function 

LinksFile2Html(linksFileName,
               resultsDir,
               htmlInfo,
               categs.descs = fCategs,
               categs.names = fCategs,
               IndexDir = "",   # IndexDir : "ResultFiles/" (Paràmetre opcional)
               UEB = UEB,      
               resultsFileName="Resultats")

