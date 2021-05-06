"UML-AssignmentI.R
Purpose:

Version [Only most relevant listed]:
0   Initialisation
1   Including figures [only ex. e left] - no comments yet.

Date:
  2021/01/18

Author:
  Stan Thijssen (TI-BDS)"

# Packages ----------------------------------------------------------------
# Cleaning
rm(list = ls())
graphics.off()

if (!require("rstudioapi"))
  install.packages("rstudioapi")
if (!require("dplyr"))
  install.packages("dplyr")
if (!require("tidyr"))
  install.packages("tidyr")
if (!require("FactoMineR"))
  install.packages("FactoMineR")
if (!require("Rfast"))
  install.packages("Rfast")
if (!require("corrplot"))
  install.packages("corrplot")
if (!require("smacof"))
  install.packages("smacof")
if (!require("PMA"))
  install.packages("PMA")
if (!require("psych"))
  install.packages("psych")
if (!require("car"))
  install.packages("car")

library(rstudioapi) # pwd
library(dplyr)
library(tidyr)
library(FactoMineR)
library(psych) 
library(Rfast)
library(corrplot)
library(PMA)
library(smacof)


# Functions --------------------------------------------------------------------
EuclDist <- function(mX){
  "
    Purpose:
      Calculate Euclidian distances between rows of n x p matrix.
      
    Inputs:
      mX:     matrix (n x p), with input vectors
      
    Return:
      mDist:  matrix (n x n), with sqrt of distance between column i and j of x 
              stored at element i,j.
  "
  iN <- nrow(mX)
  mDist <- matrix(NA, nrow=iN, ncol=iN)
  for(i in 1:iN){
    for(j in 1:iN){
      vTemp <- mX[i,] - mX[j,]
      mDist[i,j] <- sqrt(t(vTemp) %*% vTemp)
    }
  }
  return(mDist)
}

Stress <- function(mDiss, mX, mB, normalized=TRUE){
  
  iN <- nrow(mX)
  dEtaDelta2 <- (sum(lower.tri(mDiss, diag = FALSE) * (mDiss ^ 2 )))
  dEta2 <- iN * tr(t(mX) %*%  mX)
  dRho <- tr(t(mX) %*% mB %*% mX)
  dSigR <- dEtaDelta2 + dEta2 - 2 * dRho
  
  if(normalized == TRUE){
    dSig <- dSigR / dEtaDelta2
    return(dSig)
  } else{
    return(dSigR)
  }
}

SMACOF <- function(mX, mX0 = NA, typeX = 'coordinates', dEps = 10^-6, iDim = 2){
  "
    Purpose:
      Perform SMACOF using majorization algorithm.
      
    Inputs:
      mX:     matrix, either
              - n x p, with initial coordinates
              - n x n, with dissimilarities
      typeX:  coordinates (defealt) or dissimilarities
      dEps:
      
    Outputs:
      iN      number of rows in mX
      
    Return:
      mDist:  double, (normalized) Stress value of the final configuration
      ?Final: ?, final configuration itself
  "
  iN <- nrow(mX)
  iP <- ncol(mX)
  
  if(typeX == 'coordinates'){
    mDiss <- EuclDist(mX)
  } else {
    mDiss <- mX
  }
  
  if (sum(is.na(mX0)) > 0) {
    mZ <- matrix(data = rnorm(n*iDim,), nrow = n, ncol = iDim)
  } else {
    mZ <- mX0
  }
  
  mDZ <- EuclDist(mZ)
  mF <- mDiss / mDZ
  mF[!is.finite(mF)] <- 0         
  mF1 <- diag(rowSums(mF))
  mBZ <- mF1 - mF
  dSigmaOld <- Stress(mDiss = mDiss, mX = mX, mB = mBZ, normalized=TRUE)
  dSigma <- dSigmaOld
  dDiff <- dSigma - dSigmaOld
  #cat('Iteration:', 0, 'Sigma:', dSigma, 'Difference', dDiff)
  #cat("\n")
  k <- 0
  while(k == 0 || dSigmaOld - dSigma > dEps){
    k <- k + 1
    mX <- (1/iN) * mBZ %*% mZ  
    mZ <- mX
    mDZ <- EuclDist(mZ)
    mF <- mDiss / mDZ
    mF[!is.finite(mF)] <- 0         
    mF1 <- diag(rowSums(mF))
    mBZ <- mF1 - mF
    dSigmaOld <- dSigma
    dSigma <- Stress(mDiss = mDiss, mX = mX, mB = mBZ, normalized=TRUE)
    dDiff <- dSigma - dSigmaOld
    cat('Iteration:', k, 'Sigma:', dSigma, 'Difference', -dDiff)
    cat("\n")
  }
  return(list(dSigma = dSigma, mZ = mZ))
} 
# Main --------------------------------------------------------------------
# Magic Numbers
iDim <- 2

# Initialisation
setwd(dirname(getActiveDocumentContext()$path))
load("basket.RData")
mcor <- cor(basket) 
mX <- 1 - cor(basket)
mX0 <- mds(mX, ndim = iDim, verbose = TRUE, itmax = 1)$conf
lResultsCustom <- SMACOF(mX = mX, mX0 = mX0, typeX = 'dissimilarities', dEps = 10^-6, iDim =2)
lResultsPkg <- mds(mX, ndim = iDim, verbose = TRUE, eps = 10^-6, init = mX0, itmax = 250)

mds(mX, ndim = iDim, verbose = TRUE, eps = 10^-6, init = mX0, itmax = 250)

# Generating Figures ------------------------------------------------------
ggplot(as.data.frame(lResultsCustom$mZ), aes(x=D1, y=D2)) + geom_point(size=2) +
  geom_text(label=rownames(lResultsCustom$mZ), vjust = 1)

ggplot(as.data.frame(lResultsPkg$conf), aes(x=D1, y=D2)) + geom_point(size=2) +
  geom_text(label=rownames(lResultsPkg$conf), vjust = 1)