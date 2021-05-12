library(glmnet)
library(xtable)
library(mvtnorm)

set.seed(42)

X <- matrix(rnorm(500*450, 0, 1/sqrt(500)), nrow=500, ncol=450)

# s <- min(eigen(sigma)$values)
# s <- min(2*s,1)
# sseq <- c(rep(s,p))
# V <- 2*diag(sseq)-diag(sseq)%*%solve(sigma)%*%diag(sseq)
# mu<-X-X%*%solve(sigma)%*%diag(sseq)

# knockoffs setup, used in all knockoff scenarios
s <- 2
sseq <- rep(s, 450)
V <- 2 * diag(sseq) - diag(sseq) %*% diag(sseq)
mu_knock <- X - X%*%diag(sseq)
Xn <- mu_knock + rmvnorm(500, rep(0, 450), V) / sqrt(500)
cX <- cbind(X, Xn)

# basic setup to perform simulations
ks <- c(5, 20, 50)
task <- 'task1'
rp1 <- 100

# TODO: for iii) - vii) estimate FDR and power
# TODO: estimate MSE of b and mu for all apart iii)

final_table_fdr <- matrix(nrow=3, ncol=6)
final_table_power <- matrix(nrow=3, ncol=6)
final_table_mse_b <- matrix(nrow=3, ncol=7)
final_table_mse_mu <- matrix(nrow=3, ncol=7)

# totaly based on prof. Bogdan idea to calculate power and FDR
power_general <- function(u, b, b_hat, k, q=0.2){
  result <- sort(abs(u), decreasing = T, index.return = T)
  fd <- cumsum(u[result$ix] < 0)
  nd <- cumsum(u[result$ix] > 0)
  
  fdr <- (fd + 1)/nd
  
  fdp <- 0
  tp <- 0
  u1 <- which(fdr < q)
  if(length(u1) > 0){
    indopt <- max(u1)
    a1 <- result$ix[1:indopt]
    a2 <- which(u > 0)
    a3 <- intersect(a1, a2)
    new_b_hat <- rep(0, 450)
    new_b_hat[a3] <- b_hat[a3]
    tp <- sum(abs(b[a3]) > 0)
    fd <- length(a3) - tp
    if(length(a3) > 0){
      fdp <- fd/(fd+tp)
    }
    
  }
  return(c(tp/k, fdp))
}

find_lambda <- function(){
  result <- sort(abs(u), decreasing = T, index.return = T)
  fd <- cumsum(u[result$ix] < 0)
  nd <- cumsum(u[result$ix] > 0)
  
  fdr <- (fd + 1)/nd
  
  u1 <- which(fdr < q)
  
  return(min(fdr[u1]))
}

for(j in 1:3){
  # creating new beta for each value of k
  k <- ks[j]
  b <- c(rep(10, k), rep(0, 450-k))
  
  # a) power and FDR
  # columns:
  # 1. knockoffs ridge
  # 2. knockoffs LASSO
  # 3. adaptive LASSO 1
  # 4. adaptive LASSO 2
  # 5. adaptive SLOPE
  # 6. adaptive Bayesian SLOPE
  fdr_matrix <- matrix(nrow=rp1, ncol=6)
  power_matrix <- matrix(nrow=rp1, ncol=6)
  
  # b) MSE of b and mu
  # columns:
  # 1. least squeares
  # 2. ridge
  # 3. LASSO
  # 4. adaptive LASSO 1
  # 5. adaptive LASSO 2
  # 6. adaptive SLOPE
  # 7. adaptive Bayesian SLOPE
  mse_b_matrix <- matrix(nrow=rp1, ncol=7)
  mse_mu_matrix <- matrix(nrow=rp1, ncol=7)
  cat('Currently at k =', k, '\n')
  for(i in 1:rp1){
    cat('\nCurrent iteration:', i, '\n')
    eps <- 2*rnorm(500)
    Y <- X %*% b + eps
    
    # i)
    # TODO: least squares regression
    cat('OLS...\t\t\t\t')
    
    obj_ols <- lm(Y~X-1)
    b_hat_ols <- summary(obj_ols)$coefficients[,1]
    
    mse_b_matrix[i, 1] <- sum((b_hat_ols - b)^2)
    mse_mu_matrix[i, 1] <- sum((X %*% (b_hat_ols - b))^2)
    
    cat('DONE\n')
    
    # ii)
    # TODO: ridge with tuning parameter selected by CV
    cat('Ridge...\t\t\t')
    
    obj_ridge <- cv.glmnet(X, Y, standarize=F, intercept=F, alpha=0)
    b_hat_ridge <- coef(obj_ridge, s='lambda.min')[2:451]
    
    mse_b_matrix[i, 2] <- sum((b_hat_ridge - b)^2)
    mse_mu_matrix[i, 2] <- sum((X %*% (b_hat_ridge - b))^2)
    
    cat('DONE\n')
    
    # TODO: LASSO with tuning parameters selected by CV
    cat('LASSO...\t\t\t')
    
    obj_lasso <- cv.glmnet(X, Y, standarize=F, intercept=F)
    b_hat_lasso <- coef(obj_ridge, s='lambda.min')[2:451]
    
    mse_b_matrix[i, 3] <- sum((b_hat_lasso - b)^2)
    mse_mu_matrix[i, 3] <- sum((X %*% (b_hat_lasso - b))^2)
    
    cat('DONE\n')
    
    # iii)
    # TODO: knockoffs with ridge with FDR at level 0.2
    cat('Ridge knockoffs...\t\t')
    
    lambda_fdr_ridge <- 0
    obj_ridge_knockoffs <- glmnet(cX, Y, standarize=F, intercept=F, alpha=0, lambda=lambda_fdr_ridge)
    b_hat_ridge_knockoffs <- coef(obj_ridge_knockoffs, s='lambda.min')
    bhrk <- b_hat_ridge_knockoffs
    result_ridge <- abs(bhrk[2:451]) - abs(bhrk[452:901])
    
    res_stats_ridge <- power_general(result_ridge, b, b_hat_ridge, k, q=0.2)
    
    fdr_matrix[i, 1] <- res_stats_ridge[2]
    power_matrix[i, 1] <- res_stats_ridge[1]
    
    cat('DONE\n')
    
    # TODO: knockoffs with LASSO with FDR at level 0.2
    cat('LASSO knockoffs...\t\t')
    
    lambda_fdr_lasso <- 0
    obj_lasso_knockoffs <- glmnet(cX, Y, standarize=F, intercept=F, lambda=lambda_fdr_lasso)
    b_hat_lasso_knockoffs <- coef(obj_lasso_knockoffs, s='lambda.min')
    bhlk <- b_hat_lasso_knockoffs
    result_lasso <- abs(bhlk[2:451]) - abs(bhlk[452:901])
    
    res_stats_lasso <- power_general(result_lasso, b, b_hat_lasso, k, q=0.2)
    
    fdr_matrix[i, 2] <- res_stats_lasso[2]
    power_matrix[i, 2] <- res_stats_lasso[1]
    
    cat('DONE\n')
    
    # iv)
    # TODO: Adaptive LASSO I
    cat('Adaptive LASSO I...\t\t')
    
    fdr_matrix[i, 3] <- NA_real_
    power_matrix[i, 3] <- NA_real_
    mse_b_matrix[i, 4] <- NA_real_
    mse_mu_matrix[i, 4] <- NA_real_
    
    cat('DONE\n')
    
    # v)
    # TODO: Adaptive LASSO II
    cat('Adaptive LASSO II...\t\t')
    
    fdr_matrix[i, 4] <- NA_real_
    power_matrix[i, 4] <- NA_real_
    mse_b_matrix[i, 5] <- NA_real_
    mse_mu_matrix[i, 5] <- NA_real_
    
    cat('DONE\n')
    
    # vi)
    # TODO: Adaptive SLOPE, keep FDR at level 0.2
    cat('Adaptive SLOPE...\t\t')
    
    fdr_matrix[i, 5] <- NA_real_
    power_matrix[i, 5] <- NA_real_
    mse_b_matrix[i, 6] <- NA_real_
    mse_mu_matrix[i, 6] <- NA_real_
    
    cat('DONE\n')
    
    # vii)
    # For extra points +5, tempting
    # TODO: Adaptive Bayesian SLOPE with FDR at level 0.2
    cat('Adaptive Bayesian SLOPE...\t')
    
    fdr_matrix[i, 6] <- NA_real_
    power_matrix[i, 6] <- NA_real_
    mse_b_matrix[i, 7] <- NA_real_
    mse_mu_matrix[i, 7] <- NA_real_
    
    cat('DONE\n')
  }
  final_table_fdr[j, ] <- apply(fdr_matrix, 2, mean)
  final_table_power[j, ] <- apply(power_matrix, 2, mean)
  final_table_mse_b[j, ] <- apply(mse_b_matrix, 2, mean)
  final_table_mse_mu[j, ] <- apply(mse_mu_matrix, 2, mean)
}

df_fdr <- as.data.frame(
  final_table_fdr,
  row.names=c('k=5', 'k=20', 'k=50')
)

colnames(df_fdr) <- c(
  'ridge_knockoffs',
  'LASSO_knockoffs',
  'AdaLASSO1',
  'AdaLASSO2',
  'AdaSLOPE',
  'AdaBayesSLOPE'
)

df_power <- as.data.frame(
  final_table_power,
  row.names=c('k=5', 'k=20', 'k=50')
)

colnames(df_power) <- c(
  'ridge_knockoffs',
  'LASSO_knockoffs',
  'AdaLASSO1',
  'AdaLASSO2',
  'AdaSLOPE',
  'AdaBayesSLOPE'
)

df_mse_b <- as.data.frame(
  final_table_mse_b,
  row.names=c('k=5', 'k=20', 'k=50')
)

colnames(df_mse_b) <- c(
  'LS',
  'ridge',
  'LASSO',
  'AdaLASSO1',
  'AdaLASSO2',
  'AdaSLOPE',
  'AdaBayesSLOPE'
)

df_mse_mu <- as.data.frame(
  final_table_mse_mu,
  row.names=c('k=5', 'k=20', 'k=50')
)

colnames(df_mse_mu) <- c(
  'LS',
  'ridge',
  'LASSO',
  'AdaLASSO1',
  'AdaLASSO2',
  'AdaSLOPE',
  'AdaBayesSLOPE'
)

write.csv(df_fdr, file=paste(task, 'FDR.csv', sep=''))
write.csv(df_power, file=paste(task, 'power.csv', sep=''))
write.csv(df_mse_b, file=paste(task, 'MSE_b.csv', sep=''))
write.csv(df_mse_mu, file=paste(task, 'MSE_mu.csv', sep=''))

cat('Results wrtitten to .csv files!\n')

?glmnet
