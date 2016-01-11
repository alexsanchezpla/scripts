# http://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/
# http://r-pkgs.had.co.nz/git.html

require("devtools")
require("roxygen2")
pkgDir <- "E:/Dropbox (VHIR)/Scripts/Rpackages"
setwd(pkgDir)
create("links2File")


# agregar funciones a carpeta "R" 

# añadir documentacion en script funcion.R 

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

check()
# per instalar 
setwd(pkgDir)
# install("links2File")  # desde local


install_github('alexsanchezpla/scripts/Rpackages/links2File') # desde github
require(links2File)
?links2File