library(smacof)

fnEuclDist <- function(mX){
  
  #calculate cross product
  mTCross <- tcrossprod(mX)
  
  #calculate diagonals
  vC <- diag(mTCross)
  
  #calculate distance matrix
  mDist <- outer(vC, vC, FUN = "+") - 2 * mTCross
  
  # take square root
  mDist <- sqrt(mDist)
  return(mDist)
}

fnSMACOF<- function(mDiss, mX = NULL, dEps = 1e-6, verbose = TRUE){
  # get size
  iN <- nrow(mDiss)
  
  # check dissimilarity matrix
  
  #fill mX with random numbers
  if (is.null(mX)){
    mX = matrix(runif(2*iN)*max(mDiss), ncol = 2, nrow = iN)
  }
  
  #set counter
  iK <- 0
  
  #initialize Z and compute distances, stress
  mZ <- mX
  mDist <- fnEuclDist(mZ)
  dSumDiss <- sum(mDiss^2)
  dStress <- 0.5*sum((mDist - mDiss)^2)/dSumDiss
  dTol <- dEps + 1
  
  while (iK == 0 || dTol > dEps){
    #update k
    iK <- iK + 1
    
    #find Bx
    mF <- ifelse(mDist > 0, mDiss/mDist, 0)
    mF1 <- diag(rowSums(mF))
    mBx <- mF1 - mF
    
    #update X, z + calculate distances
    mX = (1/iN)*(mBx %*% mZ)
    mZ = mX
    mDist <- fnEuclDist(mZ)
    
    #caclulate stress and find toll
    dStressOld <- dStress
    dStress <- 0.5*sum((mDist - mDiss)^2)/dSumDiss
    dTol <- dStressOld - dStress
    
    #update
    if(verbose){
      cat("iteration ", formatC(iK), " Obj value ", formatC(dStress, digits = 5),
        " improvement ", formatC(dTol), "\n")
    }
  }
  return(list(dStress = dStress, mX = mX))
}

set.seed(1234)
load("~/UML/basket.Rdata")
mBasket <- as.matrix(basket)
mDiss <- fnEuclDist(mBasket)

lResults <- fnSMACOF(mDiss = mDiss)
lBestResults <- lResults


for (i in 1:1000){
  print(i)
  lResults <- fnSMACOF(mDiss = mDiss, verbose = FALSE)
  if(lResults$dStress<= lBestResults$dStress){
    lBestResults <- lResults
  }
}

plot(x = lBestResults$mX[,1], y = lBestResults$mX[,2], xlab = "Axis 1", ylab = "Axis 2")
text(x =  lBestResults$mX[,1], y = lBestResults$mX[,2], colnames(basket), col = "red")


#comparison
lResSmacof1 <- mds(delta = mDiss, type = "ratio", itmax = 1)

lResSmacof <- mds(delta = mDiss, type = "ratio")#, init = lResSmacof1$conf)
lResOwn <- fnSMACOF(mDiss, mX = lResSmacof$init)

plot(x = lResOwn$mX[,1], y = lResOwn$mX[,2], xlab = "Axis 1", ylab = "Axis 2")
text(x =  lResOwn$mX[,1], y = lResOwn$mX[,2], colnames(basket), col = "red")

plot(x = lResSmacof$conf[,1], y = lResOwn$conf[,2], xlab = "Axis1", ylab = "Axis 2")
text(x =lResSmacof$conf[,1], y = lResOwn$conf[,2], colnames(basket), col = "red")


