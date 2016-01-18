# http://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/
# http://r-pkgs.had.co.nz/git.html

if (!require("devtools")) install.packages("devtools", dep=TRUE)
if (!require("roxygen2")) install.packages("roxygen2", dep=TRUE)

require(devtools)
require("roxygen2")
SO <- version[["os"]]
if (SO=="linux-gnu")
{pkgDir <- "~/Dropbox (VHIR)/Scripts/Rpackages"
}else{
  pkgDir <- "E:/Dropbox (VHIR)/Scripts/Rpackages"
}

setwd(pkgDir)
create("links2File")

# agregar funciones a carpeta "R" 

#' A Cat Function
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' cat_function()

# procesar documento 
setwd(file.path(pkgDir, "links2File"))
document()

# Opcional: Si se quiere crear una vinyeta
# use_vignette("avignette") En linux no funciona. Pot ser pels caràcters?


check()
# per instalar 
setwd(pkgDir)
install("links2File")  # desde local


install_github('alexsanchezpla/scripts/Rpackages/links2File') # desde github
require(links2File)
?links2File
