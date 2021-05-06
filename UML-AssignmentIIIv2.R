"UML-AssignmentII.R
  Purpose:
      The purpose of this R-script is to solve a Multi Arm Bandit problem.


  Version:
  0   Initiatlisation
  1   Functions implemented for UCB and TS
  2   Start with TS - should refine list with articles - PROBLEM!!

  Date:
    2021/02/05

  Author:
    Team 1"

# Install and load packages.
if (!require("pacman"))
  install.packages("pacman")
pacman::p_load(
  #base,
  #BBmisc,
  #cluster,
  #GGally,
  #ggplot2,
  #klustR,
  #latex2exp,
  #NbClust,
  readr
  #rstudioapi,
  #texreg,
  #tictoc,
  #tidytuesdayR,
  #tidyverse
)

# Clear variables before running and set working directory to current folder.
remove(list = ls())
setwd(dirname(getActiveDocumentContext()$path))

# Functions ----------------------------------------------------------------
fSelectUCB <- function(mSelectUCB, iC){
  "
    Purpose:
      Select an arm based un the Upper-Confidence Bound (UCB)
      
    Inputs:
      mSelectUCB  matrix, including the relevant information, for UCB, i.e.
                  $iArm: integer, index of the arm A
                  $iNt: integer, how often has arm A has been drawn
                  $dQ: double, sample average for arm A
      iC          integer, controls the degree of exploration
            
    Return:
      iArm:    integer, index indicating which arm has been chosen.
  "
  
  VObj <- dfSelection$dQ + iC * sqrt(log(sum(dfSelection$iNt)) / dfSelection$iNt) 
  iArm <- dfSelection$iArm[which.max(vObj)]
  
  return(iArm)
}

fSelectTS <- function(mSelectTS, iNSample = 1000){
  "
    Purpose:
      Select an arm based un Thompson Sampling (TS)
      
    Inputs:
      mSelectTS     matrix, including the relevant information, for UCB, i.e.
                    - iArm: integer, index of the arm A
                    - iNt: integer, how often has arm A has been drawn
                    - dQ: double, sample average for arm A
                    - dAlpha: double, (current) value of alpha
                    - dBeta: double, (current) value of beta 
      iNSample      sample size to base means on, by default 1000.    
    Return:
      iArticle:    integer, index indicating which arm has been chosen.
  "
  
  iNumberArms <- nrow(mSelectTS)
  mSampleResults <- matrix(NA, nrow = iNumberArms, ncol = iNSample)
  for(i in 1:iNumberArms){
    mSampleResults[i,] <- rbeta(iNSample, mSelectTS[i, 4], mSelectTS[i, 5])
  }
  vObj <- rowMeans(mSampleResults)
  iArticle <- mSelectTS[which.max(vObj),1]
  return(iArticle)
}

fData <- function(dfData, iNArticles){
  "
    Purpose:
      Cleans and samples data.
      
    Inputs:
      dfData   dataframe, with all observations
            
    Return:
      dfSample dataframe, to work with from now on. 
  "
  iNObs <- dim(dfData)[1]
  iNTimeSteps <- dim(unique(dfData[,1]))[1]
  iNArticles <- dim(unique(dfData[,2]))[1]
  mIndArticles <- matrix(0, nrow = iNArticles, ncol = iNTimeSteps)
  vArticles <- unique(dfData[,2])
  for(i in 1:iNTimeSteps){
    dfDataTemp <- subset(dfData, DateTime == unique(dfData$DateTime)[i])
    for(j in 1:iNArticles){
      if(sum(dfDataTemp$Article == as.numeric(vArticles[j,1])) > 0){
        mIndArticles[j,i] = 1
      }
    }
  }
  vFullArticles <- rowSums(mIndArticles)
  colnames(dfData)[1:3] <- c("DateTime","Article", "Result")
  
  vNPerTimeStep <- rep(NA, iNTimeSteps)
  cat('First for-loop:\n')
  for(i in 1:iNTimeSteps){
    cat('Iteration', i, '/', iNTimeSteps, '.\n')
    vNPerTimeStep[i] <- dim(subset(dfData, DateTime == unique(dfData$DateTime)[i]))[1]
  }
  iMinObsPerTimeStep <- min(vNPerTimeStep)
  dfSample <- NULL
  cat('Second for-loop:\n')
  for(i in 1:iNTimeSteps){
    cat('Iteration', i, '/', iNTimeSteps, '. \n')
    dfDataTemp <- subset(dfData, DateTime == unique(dfData$DateTime)[i])
    dfDataTemp <- dfDataTemp[sample(vNPerTimeStep[i], iMinObsPerTimeStep, replace = FALSE),]
    dfSample <- rbind(dfSample, dfDataTemp)
  }
  return(dfSample)
}

  
# Main --------------------------------------------------------------------
# Magic numbers
dAlpha0TS = 2    # Initial parameters for TS
dBeta0TS = 2     # Initial parameters for TS

# Initialisation
dfData <- read_table2("Documents/Studie/TI - BDS/Unsupervised Machine Learning/Assignment 3/ydata-fp-td-clicks-v1_0.20090505")
dfSample <- fData(dfData, iNArticles)
iNArticles <- dim(unique(dfSample[,2]))[1]
write.csv(dfSample, 'dfSample.csv')

mSelectTS <- matrix(0, nrow = iNArticles, ncol = 5)
mSelectTS[,1] <- unique(dfSample$Article)
mSelectTS[,4] <- dAlpha0TS
mSelectTS[,5] <- dBeta0TS
colnames(mSelectTS) <- c("Article", "Number pulls", "Average", "Alpha", "Beta")

vArticle = rep(NA, iNTimeSteps)
for(i in 1:iNTimeSteps){
  vArticle[i] = fSelectTS(mSelectTS)
  dfSampleTemp <- subset(dfSample, DateTime == unique(dfData$DateTime)[i])
  dfSampleTemp <- subset(dfSampleTemp, Article == vArticle[i])
  print(mean(dfSampleTemp$Result))
  #if 
  #mSelectTS
}






