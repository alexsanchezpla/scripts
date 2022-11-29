
###
### Build  a table with info about variables to create the excel file with group names
###

varsTable <- function(x){
  index <- 1:ncol(x)
  vars <- colnames(x)
  clase <- sapply(x, class)
  tipo <- sapply(x, function (column) ifelse (is.factor(column) | is.character(column), "n", "c"))
  df <- data.frame(index=index, vars=vars, clase=clase, tipo=tipo, Group="pendingDefinition")
}

###
### Create the groups object to be used for MFA
###

extractVarNames <- function(varsAndGroupsDF, label){
  return(with(varsAndGroupsDF, vars[Group==label]))
}

creaGrups <- function (myDf, GroupsNames){
  tipoDatos <- sapply(myDf, class)
  myVarsNames <- data.frame(idx =1:ncol(myDf), 
                            vars = colnames(myDf),
                            clase = tipoDatos)
  varsAndGroups <- left_join(myVarsNames, GroupsNames, by="vars") %>%
    select(index, vars, clase.x, Group, tipo)
  groupsOfVars <- as.data.frame.matrix(as.table(with(varsAndGroups,
                                                     table(Group, tipo)))) 
  groupTypes <- groupsOfVars$varTypes <- 
    colnames(groupsOfVars)[apply(groupsOfVars, 1, function(r) which (r > 0))]
  numericCols <- 1:(ncol(groupsOfVars)-1)
  numGroups <-  groupsOfVars$sizes <-  rowSums(groupsOfVars[,numericCols])
  groupNames <- rownames(groupsOfVars)
  groupsOfVars <- groupsOfVars |> dplyr::select(sizes, varTypes)
  listOfVars <- list()
  for (i in 1:length(groupNames)){
    listOfVars [[i]]<- list (vars =extractVarNames(varsAndGroups,
                                                   groupNames[i]), 
                             num =groupsOfVars$sizes[i] , 
                             name = groupNames[i], 
                             type = groupsOfVars$varTypes[1] )
  }
  names(listOfVars) <- groupNames
  listOfDatasets<- list()
  for (i in 1:length(listOfVars)){
    listOfDatasets[[i]] <- myDf[, listOfVars[[i]][[1]]]
  }
  names(listOfDatasets) <- names(listOfVars)
  return(list(groupsData=listOfDatasets, groupsVars=listOfVars, 
              groupsNames=groupNames, groupsSizes=numGroups, groupTypes=groupTypes))
}

### 
### checkFactorialStructure of the created groups object
###

showGroupsinList <- function(alistOfVars){
  for (i in 1:length(alistOfVars)){
    cat(names(alistOfVars[i]), "\n")
    cat("\t", "numVars = ", alistOfVars[[i]][[2]], "\n")
    cat("\t", "groupName = ", alistOfVars[[i]][[3]], "\n")
    cat("\t", "groupType = ", alistOfVars[[i]][[4]], "\n")
  }
}

###
### Create a unique data frame from all or some of the datasets on which the original data has been broken.
###

MultMerge2 <- function (lst, all.x = TRUE, all.y = TRUE, by = NULL) 
{
  # lst <- list(...) # The original version had "..." instead of "lñst" as argument
  if (length(lst) == 1) 
    return(lst[[1]])
  if (!is.null(by)) {
    for (i in seq_along(lst)) {
      rownames(lst[[i]]) <- lst[[i]][[by]]
      lst[[i]][by] <- NULL
    }
  }
  unames <- DescTools::SplitAt(make.unique(unlist(lapply(lst, colnames)), 
                                           sep = "."), cumsum(sapply(head(lst, -1), ncol)) + 1)
  for (i in seq_along(unames)) colnames(lst[[i]]) <- unames[[i]]
  res <- Reduce(function(y, z) merge(y, z, all.x = all.x, all.y = all.x), 
                lapply(lst, function(x) data.frame(x, rn = row.names(x))))
  rownames(res) <- res$rn
  res$rn <- NULL
  seq_ord <- function(xlst) {
    jj <- character(0)
    for (i in seq_along(xlst)) {
      jj <- c(jj, setdiff(xlst[[i]], jj))
    }
    return(jj)
  }
  ord <- seq_ord(lapply(lst, rownames))
  res[ord, ]
  if (!is.null(by)) {
    res <- data.frame(row.names(res), res)
    colnames(res)[1] <- by
    rownames(res) <- c()
  }
  return(res)
}

### test

# df1 <- data.frame(matrix(rnorm(20), nrow=10)); rownames(df1) <- paste("row",0:9, sep="")
# df2 <- data.frame(min=letters[1:10], may=LETTERS[1:10]); rownames(df2) <- paste("row",0:9, sep="")
# df3 <- data.frame(df1 < 0); rownames(df3) <- paste("row",0:9, sep="")
# dfList <- list(df1, df2, df3)
# library(DescTools)
# MultMerge2(list(df1, df2, df3))
# DescTools::MultMerge(df1, df2, df3)


showText <- function (aText){
  cat("\n",aText,"\n")
  cat(paste(rep("=", nchar(aText)), collapse=""),"\n")
}

checkFactorialStructure <- function (mylistOfGroups)
  #                          varsList, groupSizes, groupNames, groupTypes)
{
  mylistOfDataSets <- mylistOfGroups[1][[1]]
  uniqueDataSet <- MultMerge2 (mylistOfDataSets)
  varsList <- colnames(uniqueDataSet)
  actualTypes <- sapply(uniqueDataSet[,varsList], class)
  grupos <- character()
  tipos <- character()
  for (i in 1:length(mylistOfDataSets)){
    groupNames <-mylistOfGroups$groupsNames[i]
    groupSizes <- mylistOfGroups$groupsSizes[i]
    groupTypes <- mylistOfGroups$groupTypes[i]
    grupos <- c(grupos, rep(groupNames, groupSizes))
    assignedTypes<- tipos <- c(tipos, rep(groupTypes, groupSizes))
  }
  showText("Global dimensions")
  show(dim(uniqueDataSet))
  showText("Structure of each group")
  showGroupsinList(mylistOfGroups$groupsVars) 
  showText("Data type vs group label")
  show(table(actualTypes, assignedTypes))
  return(data.frame(Variable= varsList, TipoActual=actualTypes,
                    Grupo=grupos, TipoAsignado= tipos ))
}

