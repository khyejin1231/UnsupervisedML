library(PMA)
source("utilities.R")

data <- load("C:/Users/40532/Desktop/FIFA2017_NL.RData")
data <- eval(parse(text = data))
X <- scale(data[, 4:32])



pca <- function(X, cumulative_contribution = 0.75){
  n <- nrow(X)
  p <- ncol(X)
  k <- 0
  ratio <- 0
  res_svd <- svd(X/sqrt(n-1))
  comp_loadings <- res_svd$v
  comp_var <- res_svd$d ** 2
  factor_loadings <- matrix(0L, p, p)
  for(i in 1:p){
      factor_loadings[i, ] <- sqrt(comp_var[i]) * comp_loadings[,i] 
  }
  total_var <- sum(comp_var)
  while(k < p & ratio < cumulative_contribution){
    k <- k + 1
    ratio <- sum(comp_var[k])/total_var 
  }
  trun_comp_loadings <- comp_loadings[, 1:k]
  trun_comp_var <- comp_var[1:k]
  trun_factor_loadings <- rbind(factor_loadings[1:k, ], colSums(factor_loadings[1:k, ]))
  return(list("comp_loadings" = comp_loadings, "comp_var" = comp_var, 
              "factor_loadings" = factor_loadings, "trun_comp_loadings" = trun_comp_loadings,
              "trun_comp_var" = trun_comp_var, "trun_factor_loadings" = trun_factor_loadings))
}

l2norm <- function(a) sqrt(sum(a^2))

s_rank_1 <- function(X, c2, imax){
  lsm <- svd(X)$v
  v <- lsm[,1,drop = FALSE]
  i <- 0
  while(i < imax){
    i <- i + 1
    Xv <- X%*%v
    u <- Xv/l2norm(Xv)
    uX <- t(X)%*%u
    lambda <- binary_search(uX, c2)
    v <- soft_l2norm(uX, lambda)
  }
  sigma <- t(u)%*%X%*%v
  return(list("u"=u, "v"=v, "sigma"=sigma[,1]))
}




spca <- function(X, c2, imax, K){
  R <- X
  n <- nrow(R)
  p <- ncol(R)
  U <- matrix(0L,n,K)
  V <- matrix(0L,p,K)
  sigma <- c()
  for(k in 1:K){
    res <- s_rank_1(R, c2, imax)
    U[,k] <- res$u
    V[,k] <- res$v
    sigma[k] <- res$sigma
    R <- R - res$sigma*res$u %*% t(res$v)
  }
  return(list("U"=U, "V"=V, "sigma"=sigma))
}

##with package
X_pca_cor <- prcomp(X, tol = 0.08, rank. = 10)
summary(X_pca_cor)
print(X_pca_cor)
plot(X_pca_cor, type = "l")

cv_res <- SPC.cv(X, sumabsvs = seq(1.5, sqrt(ncol(X)), length.out = 100), nfolds = 5)
c2 <- cv_res$bestsumabsv
X_sparse <- SPC(X, sumabs = c2, K = 10, center = FALSE, trace = FALSE)


##with our own functions
pca_res <- pca(X)
spca_res <- spca(X, c2 = c2, imax = 100, K = 10)