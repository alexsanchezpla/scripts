
workingDir <-getwd()
dataDir <-file.path(workingDir, "datos")
resultsDir <- file.path(workingDir,"results")
targetsFile <-"targLimma.txt"
anotPackage <-"hgu95av2.db"


setwd(workingDir)

options(width=80)
options(digits=5)
memory.limit(4095)


processa <-function(Filename){
  Sweave(paste(Filename,"Rnw",sep="."))
  shell(paste("copy", paste(Filename,"tex",sep="."), ".\\informes\\*.tex", sep=" "))
  shell("if exist .\\images\\*.* copy .\\images\\*.* .\\informes\\images\\*.*")
  shell("if exist .\\images\\*.eps del .\\images\\*.eps")
  shell("if exist .\\images\\*.pdf del .\\images\\*.pdf")
  shell(paste("erase", paste(Filename,"tex",sep="."),sep=" "))
  Stangle(paste(Filename,"Rnw",sep=".")) 
}

