# http://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/
# http://r-pkgs.had.co.nz/git.html
# http://kbroman.org/github_tutorial/pages/routine.html

if (!require("devtools")) install.packages("devtools", dep=TRUE)
if (!require("roxygen2")) install.packages("roxygen2", dep=TRUE)

require(devtools)
require(roxygen2)
SO <- version[["os"]]
if (SO=="linux-gnu")
{pkgDir <- "~/Dropbox (VHIR)/Scripts/Rpackages"
}else{
  pkgDir <- "C:/users/ALexandre/Dropbox (VHIR)/Scripts/Rpackages"
}

setwd(pkgDir)
packgName <-"geneLists"

# Per crear el paquet fer servir la instrucció de devtools
# create(packgName)

# Cal que al menys hi hagi una funció a la carpeta "R" escrita en format "roxygen". Per exemple:

#' Print hello
#'
#' This function prints hello 
#'
#' @param fname First name
#' @param lname Last name
#' @export
#' @examples
#' hello(fname="Tomas",lname="Greif")
#' hello(fname="Your",lname="Name")
hello <- function(fname, lname) {
  cat(paste("Hello",fname,lname,"!"))
}

# A partir de l'anterior es crearà la plantilla pel paquet
setwd(file.path(pkgDir, packgName))
# document()

# Opcional: Si se quiere crear una vinyeta
# use_vignette("avignette") 

# per instalar
setwd(pkgDir)
install(packgName)  # desde local
check()

install_github(paste('alexsanchezpla/scripts/Rpackages/', packgName, sep="")) # desde github
require(packgName)
?pckgName
