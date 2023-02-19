#' Load variable from binary file by position
#'
#' This function creates a variable with the contents of a binary file. It allows to give the variable a custom name instead of the name it had when it was saved using `load()`.
#'
#' @param fileName Character string. Path to the binary file to be loaded.
#' @param pos Numeric. Position of the variable to be loaded in the list of variables stored in the binary file. Defaults to 1.
#'
#' @return The variable loaded from the binary file.
#'
#' @export
loadFromFileByPos <- function (fileName, pos=1) {
  tempEnv <- new("environment")
  load(fileName, tempEnv)
  varNames <- ls(tempEnv)
  myVarName <- varNames[pos]
  load(fileName)
  myVar <- eval(parse(text = myVarName))
  return(myVar)
}

#' Load variable from binary file by name
#'
#' This function creates a variable with the contents of a binary file. It allows to give the variable a custom name instead of the name it had when it was saved using `load()`.
#'
#' @param fileName Character string. Path to the binary file to be loaded.
#' @param aVarName Character string. Name of the variable to be loaded.
#'
#' @return The variable loaded from the binary file.
#'
#' @export
loadFromFileByName <- function (fileName, aVarName) {
  tempEnv <- new("environment")
  load(fileName, tempEnv)
  varNames <- ls(tempEnv)
  myVarName <- varNames[varNames %in% aVarName]
  load(fileName)
  myVar <- eval(parse(text = myVarName))
  return(myVar)
}

