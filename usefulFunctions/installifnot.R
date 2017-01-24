installBiocifnot <- function(pckgName){
  if (!(require(pckgName, character.only = TRUE))) {
    source("http://Bioconductor.org/biocLite.R")
    biocLite(pckgName)
    require(pckgName, character.only = TRUE)
  }
}
# Example
# installBiocifnot("CMA")

installifnot <- function(pckgName){
  if (!(require(pckgName, character.only = TRUE))) {
    install.packages(pckgName, dep = TRUE)
    require(pckgName, character.only = TRUE)
  }
}
# Example
# installifnot("xlsx")

installGitifnot <- function(pathGit, pckgName, proxy = FALSE, urlproxy = "conf_www.ir.vhebron.net", portproxy = 8081, force.install = FALSE){
  if (!(require(pckgName, character.only = TRUE)) | force.install) {
    installifnot("devtools")
    # require(curl)
    if (proxy) set_config(use_proxy(url = urlproxy, port = portproxy))
    install_github(file.path(pathGit,pckgName), force = force.install)
    require(pckgName, character.only = TRUE)
  }
}
# Example
# installGitifnot("miriammota","mmotaF",force.install = TRUE)
