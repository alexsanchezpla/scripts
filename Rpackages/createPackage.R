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

# Si tenim dades i les volem documentar podem fer-ho posant en la carpeta "R" un arxiu ".R" amb un text com:
#' This is data to be included in my package
#'
#' @author My Name \email{blahblah@@roxygen.org}
#' @references \url{data_blah.com}
"data-name"
# devtools::use_data(AvsB, geneLists, overwrite = TRUE)

# Veure també: http://stackoverflow.com/questions/2310409/how-can-i-document-data-sets-with-roxygen

# A partir de l'anterior es crearà la plantilla pel paquet
setwd(file.path(pkgDir, packgName))
document()

# Opcional: Si es vol crear una vinyeta es pot generar la plantilla fent:
# Si un cop creada el tornem a invocar ens donara error, es a dir nomes s'ha de fer un cop
# use_vignette("avignette")

# per comprovar si funciona fem:
setwd(file.path(pkgDir, packgName))
check()

# per instalar
setwd(pkgDir)
install(packgName)  # desde local

# per posar-lo a github actualitzat
# des d'una consola del sistema (linux) o de git (windows)fer:
# cd Dropbox\ \(VHIR\)/Scripts
# git add . -A
# git commit -m "Missatge"
# git push -u origin master

# A partir d'això podem instal·lar-lo en qualsevol ordinador amb la instrucció 'install_github'
require(devtools)
install_github(paste('alexsanchezpla/scripts/Rpackages/', packgName, sep="")) # desde github
require(packgName, character.only=TRUE)

