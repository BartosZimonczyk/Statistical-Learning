library(far)
library(glmnet)
library(xtable)

X <- matrix(rnorm(500*450, 0, 1/sqrt(500)), nrow=500, ncol=450)
ks <- c(5, 20, 50)
task <- 'task1'
rp1 <- 100

# TODO: for iii) - vii) estimate FDR and power
# TODO: estimate MSE of b and mu for all apart iii)

final_table_fdr <- matrix(nrow=3, ncol=5)
final_table_power <- matrix(nrow=3, ncol=5)
final_table_mse_b <- matrix(nrow=3, ncol=6)
final_table_mse_mu <- matrix(nrow=3, ncol=6)


for(j in 1:3){
  k <- ks[j]
  b <- c(rep(10, k), rep(0, 450-k))
  # a)
  # columns:
  # 1. knockoffs ridge/lasso
  # 2. adaptive LASSO 1
  # 3. adaptive LASSO 2
  # 4. adaptive SLOPE
  # 5. adaptive Bayesian SLOPE
  fdr_matrix <- matrix(nrow=rp1, ncol=5)
  power_matrix <- matrix(nrow=rp1, ncol=5)
  
  # b)
  # columns:
  # 1. least squeares
  # 2. ridge/LASSO
  # 3. adaptive LASSO 1
  # 4. adaptive LASSO 2
  # 5. adaptive SLOPE
  # 6. adaptive Bayesian SLOPE
  mse_b_matrix <- matrix(nrow=rp1, ncol=6)
  mse_mu_matrix <- matrix(nrow=rp1, ncol=6)
  for(i in 1:rp1){
    eps <- 2*rnorm(500)
    Y <- X %*% b + eps
    
    # i)
    # TODO: least squares regression
    
    mse_b_matrix[i, 1] <- NA_real_
    mse_mu_matrix[i, 1] <- NA_real_
    
    # ii)
    # TODO: ridge regression with LASSO
    #       with tuning parameters selected by CV
    
    mse_b_matrix[i, 2] <- NA_real_
    mse_mu_matrix[i, 2] <- NA_real_
    
    # iii)
    # TODO: knockoffs with ridge with LASSO keeping
    #       FDR at level 0.2
    
    fdr_matrix[i, 1] <- NA_real_
    power_matrix[i, 1] <- NA_real_
    
    # iv)
    # TODO: Adaptive LASSO I
    
    fdr_matrix[i, 2] <- NA_real_
    power_matrix[i, 2] <- NA_real_
    mse_b_matrix[i, 3] <- NA_real_
    mse_mu_matrix[i, 3] <- NA_real_
    
    # v)
    # TODO: Adaptive LASSO II
    
    fdr_matrix[i, 3] <- NA_real_
    power_matrix[i, 3] <- NA_real_
    mse_b_matrix[i, 4] <- NA_real_
    mse_mu_matrix[i, 4] <- NA_real_
    
    # vi)
    # TODO: Adaptive SLOPE, keep FDR at level 0.2
    
    fdr_matrix[i, 4] <- NA_real_
    power_matrix[i, 4] <- NA_real_
    mse_b_matrix[i, 5] <- NA_real_
    mse_mu_matrix[i, 5] <- NA_real_
    
    # vii)
    # For extra points +5, tempting
    # TODO: Adaptive Bayesian SLOPE with FDR at level 0.2
    
    fdr_matrix[i, 5] <- NA_real_
    power_matrix[i, 5] <- NA_real_
    mse_b_matrix[i, 6] <- NA_real_
    mse_mu_matrix[i, 6] <- NA_real_
  }
  final_table_fdr[j, ] <- apply(fdr_matrix, 2, mean)
  final_table_power[j, ] <- apply(power_matrix, 2, mean)
  final_table_mse_b[j, ] <- apply(mse_b_matrix, 2, mean)
  final_table_mse_mu[j, ] <- apply(mse_mu_matrix, 2, mean)
}

df_fdr <- as.data.frame(
  final_table_fdr,
  row.names=c('k=5', 'k=20', 'k=50'),
  col.names=c()
)






