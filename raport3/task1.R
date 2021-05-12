library(far)
library(glmnet)
library(xtable)

X <- matrix(rnorm(500*450, 0, 1/sqrt(500)), nrow=500, ncol=450)
ks <- c(5, 20, 50)

# TODO: for iii) - vii) estimate FDR and power
# TODO: estimate MSE of b and mu for all apart iii)

for(k in ks){
  b <- c(rep(10, k), rep(0, 450-k))
  for(i in 1:100){
    eps <- 2*rnorm(500)
    Y <- X %*% b + eps
    
    # i)
    # TODO: least squares regression
    
    # ii)
    # TODO: ridge regression with LASSO
    #       with tuning parameters selected by CV
    
    # iii)
    # TODO: knockoffs with ridge with LASSO keeping
    #       FDR at level 0.2
    
    # iv)
    # TODO: Adaptive LASSO I
    
    # v)
    # TODO: Adaptive LASSO II
    
    # vi)
    # TODO: Adaptive SLOPE, keep FDR at level 0.2
    
    # vii)
    # For extra points +5, tempting
    # TODO: Adaptive Bayesian SLOPE with FDR at level 0.2
    
  }
}
