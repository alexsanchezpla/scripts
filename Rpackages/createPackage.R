# http://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/
# http://r-pkgs.had.co.nz/git.html

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
# Per crear el paquet
# create(packgName)

# agregar funciones a carpeta "R"

#' A Cat Function
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' cat_function()

# procesar documents a partir de plantilles
# setwd(file.path(pkgDir, packgName))
# document()

# Opcional: Si se quiere crear una vinyeta
# use_vignette("avignette") # En linux no funciona. Pot ser pels caràcters?

# per instalar
setwd(pkgDir)
install(packgName)  # desde local


install_github(paste('alexsanchezpla/scripts/Rpackages/', packgName, sep="")) # desde github
require(packgName)
?pckgName
