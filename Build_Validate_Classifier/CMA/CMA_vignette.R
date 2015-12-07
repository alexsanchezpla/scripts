### R code from vignette source 'CMA_vignette.rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
library(CMA)


###################################################
### code chunk number 2: CMA_vignette.rnw:400-405 (eval = FALSE)
###################################################
## data(golub)
## golubY <- golub[,1]
## 
## golubX <- as.matrix(golub[,-1])
## 


###################################################
### code chunk number 3: CMA_vignette.rnw:414-420 (eval = FALSE)
###################################################
## 
## loodat <- GenerateLearningsets(y=golubY, method ="LOOCV")
## class(loodat)
## getSlots(class(loodat))
## show(loodat)
## 


###################################################
### code chunk number 4: CMA_vignette.rnw:425-430 (eval = FALSE)
###################################################
## 
## set.seed(321)
## fiveCVdat <- GenerateLearningsets(y=golubY, method = "CV", fold = 5, strat = TRUE)
## 
## 


###################################################
### code chunk number 5: CMA_vignette.rnw:437-446 (eval = FALSE)
###################################################
## 
## set.seed(456)
## MCCVdat <- GenerateLearningsets(y=golubY, method = "MCCV", niter = 3,
##                                 ntrain = floor(2/3*length(golubY)), strat = TRUE)
## 
## set.seed(651)
## bootdat <- GenerateLearningsets(y=golubY, method = "bootstrap", niter = 3, strat = TRUE)
##  
## 


###################################################
### code chunk number 6: CMA_vignette.rnw:456-462 (eval = FALSE)
###################################################
## 
## varsel_fiveCV <- GeneSelection(X = golubX, y=golubY, learningsets = fiveCVdat, method = "wilcox.test")
## varsel_MCCV <- GeneSelection(X = golubX, y=golubY, learningsets = MCCVdat, method = "wilcox.test")
## varsel_boot <- GeneSelection(X = golubX, y=golubY, learningsets = bootdat, method = "wilcox.test")
##  
## 


###################################################
### code chunk number 7: CMA_vignette.rnw:468-476 (eval = FALSE)
###################################################
## 
## show(varsel_fiveCV)
## toplist(varsel_fiveCV, iter=1)
## seliter <- numeric()
## for(i in 1:5)
## seliter <- c(seliter, toplist(varsel_fiveCV, iter = i, top = 10, show = FALSE)$index)
## sort(table(seliter), dec = TRUE)
## 


###################################################
### code chunk number 8: CMA_vignette.rnw:506-511 (eval = FALSE)
###################################################
## set.seed(351) 
## tuningstep <- CMA:::tune(X = golubX, y=golubY, learningsets = fiveCVdat,
##                       genesel = varsel_fiveCV, nbgene = 100, classifier = svmCMA,
##                       grids = list(cost = c(0.1, 1, 10, 100, 200)),probability=T)
## 


###################################################
### code chunk number 9: CMA_vignette.rnw:516-520 (eval = FALSE)
###################################################
##                       
## show(tuningstep)
## unlist(best(tuningstep))
## 


###################################################
### code chunk number 10: CMA_vignette.rnw:527-531 (eval = FALSE)
###################################################
## 
## par(mfrow = c(2,2))
## for(i in 1:4)
## plot(tuningstep, iter = i, main = paste("iteration", i))


###################################################
### code chunk number 11: CMA_vignette.rnw:537-555 (eval = FALSE)
###################################################
## 
## #class_loo <- classification(X = golubX, y=golubY, learningsets = loodat,
## #                      genesel = varsel_loo, nbgene = 100, classifier = svmCMA,
## #                      cost = 100)
##                       
## class_fiveCV <- classification(X = golubX, y=golubY, learningsets = fiveCVdat,
##                       genesel = varsel_fiveCV, nbgene = 100, classifier = svmCMA,
##                       cost = 0.1,probability=T)
##                       
## class_MCCV <- classification(X = golubX, y=golubY, learningsets = MCCVdat,
##                       genesel = varsel_MCCV, nbgene = 100, classifier = svmCMA,
##                       cost = 0.1,probability=T)
##                       
## class_boot <- classification(X = golubX, y=golubY, learningsets = bootdat,
##                       genesel = varsel_boot, nbgene = 0.1, classifier = svmCMA,
##                       cost = 100,probability=T)
##                       
## 


###################################################
### code chunk number 12: CMA_vignette.rnw:564-568 (eval = FALSE)
###################################################
## 
## resultlist <- list(class_fiveCV, class_MCCV, class_boot)
## result <- lapply(resultlist, join)
## 


###################################################
### code chunk number 13: CMA_vignette.rnw:574-582 (eval = FALSE)
###################################################
## 
## 
## schemes <- c("five-fold CV", "MCCV", "bootstrap")
## 
## par(mfrow = c(3,1))
## for(i in seq(along = result))
##  plot(result[[i]], main = schemes[i])
## 


###################################################
### code chunk number 14: CMA_vignette.rnw:588-591 (eval = FALSE)
###################################################
## 
## invisible(lapply(result, ftable))
## 


###################################################
### code chunk number 15: CMA_vignette.rnw:595-599 (eval = FALSE)
###################################################
## 
## par(mfrow = c(2,2))
## for(i in seq(along = result)) roc(result[[i]])
## 


###################################################
### code chunk number 16: CMA_vignette.rnw:605-609 (eval = FALSE)
###################################################
## 
## totalresult <- join(result)
## ftable(totalresult)
## 


###################################################
### code chunk number 17: CMA_vignette.rnw:618-624 (eval = FALSE)
###################################################
## 
## av_MCCV <- evaluation(class_MCCV, measure = "average probability")
## show(av_MCCV)
## boxplot(av_MCCV)
## summary(av_MCCV)
## 


###################################################
### code chunk number 18: CMA_vignette.rnw:639-642 (eval = FALSE)
###################################################
## av_obs_MCCV <- evaluation(class_MCCV, measure = "average probability", scheme = "obs")
## show(av_obs_MCCV)
## 


###################################################
### code chunk number 19: CMA_vignette.rnw:648-651 (eval = FALSE)
###################################################
## 
## obsinfo(av_obs_MCCV, threshold = 0.6)
## 


###################################################
### code chunk number 20: CMA_vignette.rnw:675-685
###################################################

data(khan)
khanY <- khan[,1]
khanX <- as.matrix(khan[,-1])
dim(khanX) # Trasposta respecte l'habitual!

set.seed(27611)
fiveCV5iter <- GenerateLearningsets(y=khanY, method = "CV", fold = 5, niter = 5, strat = TRUE)




###################################################
### code chunk number 21: CMA_vignette.rnw:693-698
###################################################

class_dlda <- classification(X = khanX, y=khanY, learningsets = fiveCV5iter,
                             classifier = dldaCMA)
                             



###################################################
### code chunk number 22: CMA_vignette.rnw:705-710
###################################################


genesel_da <- GeneSelection(X=khanX, y=khanY, learningsets = fiveCV5iter,
                            method = "t.test", scheme = "one-vs-all")



###################################################
### code chunk number 23: CMA_vignette.rnw:717-730
###################################################

class_lda <- classification(X = khanX, y=khanY, learningsets = fiveCV5iter,
                             classifier = ldaCMA, genesel= genesel_da,
                             nbgene = 10)
                             
class_fda <- classification(X = khanX, y=khanY, learningsets = fiveCV5iter,
                             classifier = fdaCMA, genesel = genesel_da,
                             nbgene = 10, comp = 2)
                             
class_qda <- classification(X = khanX, y=khanY, learningsets = fiveCV5iter,
                             classifier = qdaCMA, genesel = genesel_da,
                             nbgene = 1)
                             


###################################################
### code chunk number 24: CMA_vignette.rnw:737-742
###################################################

set.seed(876)
class_scda <- classification(X = khanX, y=khanY, learningsets = fiveCV5iter,
                             classifier = scdaCMA, tuninglist = list(grids = list()))



###################################################
### code chunk number 25: CMA_vignette.rnw:747-751
###################################################
require(plsgenomics)
class_plsda <- classification(X = khanX, y=khanY, learningsets = fiveCV5iter,
                             classifier = pls_ldaCMA)



###################################################
### code chunk number 26: CMA_vignette.rnw:757-763
###################################################
                             
dalike <- list(class_dlda, class_lda, class_fda, class_qda, class_scda, class_plsda)
par(mfrow = c(3,1))
comparison <- compare(dalike, plot = TRUE, measure = c("misclassification", "brier score", "average probability"))
print(comparison)
                      


