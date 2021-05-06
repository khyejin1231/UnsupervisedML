"UML-AssignmentII.R
Purpose:
    The purpose of this R-script is to perform


Version:
0   Initiatlisation
1   Means per cluster

Date:
  2021/01/29

Author:
  Team 1"

# Install and load packages.
if (!require("pacman"))
  install.packages("pacman")
pacman::p_load(
  base,
  BBmisc,
  cluster,
  GGally,
  ggplot2,
  klustR,
  latex2exp,
  NbClust,
  rstudioapi,
  texreg,
  tictoc,
  tidytuesdayR,
  tidyverse
)


# Clear variables before running and set working directory to current folder.
remove(list = ls())
setwd(dirname(getActiveDocumentContext()$path))


# Functions ---------------------------------------------------------------

# Main --------------------------------------------------------------------
# Magic numbers
iSeed <- 4322
iMinC <- 2
iMaxC <- 11
iM <- 10000
sMethod <- 'kmeans'
sDistance <- 'euclidean'
sIndex <- 'ch'

lVars = c(
  'track_popularity',
  'danceability',
  'energy',
  'loudness',
  'speechiness',
  'acousticness',
  'instrumentalness',
  'liveness',
  'valence',
  'tempo',
  'duration_ms'
)

# Load and edit data
set.seed(iSeed)
tuesdata <- tidytuesdayR::tt_load('2020-01-21')
dfSpotify <- as.data.frame(tuesdata$spotify_songs)
remove(tuesdata)
dfSpotify <- subset(dfSpotify, select = lVars)
dfSpotify <- dfSpotify[sample(nrow(dfSpotify), 10000), ]
dfSpotify[,!names(dfSpotify) %in% c("track_popularity")] = normalize(
  dfSpotify[,!names(dfSpotify) %in% c("track_popularity")],
  method = "standardize",
  range = c(0, 1),
  margin = 1L,
  on.constant = "quiet"
)

# Determine optimal number of clusters
# tic()
# iOptNClusters <- as.integer(
#   NbClust(
#     dfSpotify[, !names(dfSpotify) %in% c("track_popularity")],
#     min.nc = iMinC,
#     max.nc = iMaxC,
#     method = sMethod,
#     index = sIndex,
#   )$Best.nc[1]
# )
# toc()
# 
# tic()
# lResultsKm <-
#   kmeans(dfSpotify[, !names(dfSpotify) %in% c("track_popularity")], centers = iOptNClusters, nstart = iM)
# ggpairs(dfSpotify[, !names(dfSpotify) %in% c("track_popularity")], mapping = aes(colour = factor(lResultsKm$cluster)))
# dfSpotify$cluster <- lResultsKm$cluster
# save.image('Intermediate.Rdata')
# toc()

# Part II -----------------------------------------------------------------
load('Intermediate.Rdata')
dfCluster1 <- subset(dfSpotify, cluster == 1)
cat('Mean cluster 1:', mean(dfCluster1$track_popularity))
model1 <- lm(track_popularity ~ . - cluster, data = dfCluster1)
summary(model1)

dfCluster2 <- subset(dfSpotify, cluster == 2)
cat('Mean cluster 2:', mean(dfCluster2$track_popularity))
model2 <- lm(track_popularity ~ . - cluster, data = dfCluster2)
summary(model2)

# ### Parallel coordinates plot
# pacoplot(
#   data = dfSpotify[, !names(dfSpotify) %in% c("track_popularity", "cluster")],
#   clusters = dfSpotify$cluster,
#   labelSizes = list(yaxis = 16, yticks = 12),
#   measures = list(avg = mean)
# )
#
# ggparcoord(
#   dfCluster1[, !names(dfCluster1) %in% c("track_popularity")],
#   columns = 1:10,
#   scale = "std",
#   groupColumn = "cluster",
#   alphaLines = 0.5
# ) + facet_wrap(~ cluster) + theme(axis.text.x = element_text(size = 7))
#
# ggparcoord(
#   dfCluster2[, !names(dfCluster2) %in% c("track_popularity")],
#   columns = 1:10,
#   scale = "std",
#   groupColumn = "cluster",
#   alphaLines = 0.5
# ) + facet_wrap(~ cluster) + theme(axis.text.x = element_text(size = 7))

### OLS regression results
texreg(
  list(model1, model2),
  label = "tab:OLSresults",
  caption = "OLS regression results per cluster.",
  float.pos = "h"
)