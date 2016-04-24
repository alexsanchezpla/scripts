MI <- function(x,y,h){
# x : vector of methylations
# y : vector of expressions  
# h : the std of the Gaussian kernel for density estimation 

if(length(x) != length(y))  stop("Different number of samples!")

M <- length(x) # samples
aux <- 0
two.h.square <- 2*h^2

for(i in 1:M){
  # kernel distance between the i.th data point and all other points j for each gene k
 tmpx <- x - x[i] 
 tmpx <- exp(-tmpx^2/(two.h.square))
 tmpy <- y - y[i]
 tmpy <- exp(-tmpy^2/(two.h.square))
 tmp <- tmpx*tmpy
 aux <- aux + log(M*sum(tmp)/(sum(tmpx)*sum(tmpy)))
 }
aux/M
}

cMI <- function(dataMeth, dataExp, t, h=0.3){
# dataMeth : input data methylation (a vector, not a matrix!)
# dataExp  : input data expression (a vector, not a matrix!)
# t : moves from 0 to 1
  n <- length(dataMeth)
  if(length(dataExp) != n)  stop("Different number of samples!")
  if(t < 0 | t > 1)  stop("t value is out of range")
  filter <- dataMeth < t
  ss <- sum(filter)
  if(ss != 0){
  x <- dataMeth[filter]
  y <- dataExp[filter]
  aux <- MI(x,y,h)*ss/n
  }
  else{
    aux <- 0
  }
  ss <- sum(!filter)
  if(ss != 0){
  x <- dataMeth[!filter]
  y <- dataExp[!filter]
  aux + MI(x,y,h)*ss/n
  }
  else{
    aux
  }
}

computeCMI <- function (methData, exprData){
  # check consistency
  # Provar un try/catch
  if((nrow(exprData)!=nrow(methData))||(ncol(exprData)!=ncol(methData)))
    stop("Expression amd methylation data must have same dimensions")
  tt <- seq(0,1,by=0.01)
  nt <- length(tt)
  ngenes <- nrow(exprData)
  nsamples <- ncol(exprData)
  cmi <- matrix(rep(0,nt*ngenes), ncol=nt)
  for (j in 1:nrow(methData)){
    # cat (j, " ")
    methVec <- methData[j,]
    exprVec <- exprData[j,]
    # df.xy <- na.omit(data.frame(methVec=methData[j,], exprVec=exprData[j,]))
    for(i in 1:nt){
      cmi[j ,i] <- cMI(methVec, exprVec, t=tt[i])
    }
  }
  return(cmi)
}

