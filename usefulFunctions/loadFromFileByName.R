# Funcio per crear una variable amb el contingut d'un arxiu binari 
# Serveix per poder donar el nom que desitgem a la variable enlloc del que tenia quan s'ha fet el "load"

loadFromFileByName <-function (fileName, aVarName){
  tempEnv =new("environment")
  load (fileName, tempEnv)
  varNames <-ls(tempEnv)
  myVarName <- varNames[varNames %in% aVarName]
  load (fileName)
  myVar <- eval(parse(text=myVarName))
  return(myVar)
}

