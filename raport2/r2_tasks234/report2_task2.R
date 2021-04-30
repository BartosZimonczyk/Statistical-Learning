#############################################################
# task 2

library(mvtnorm)
library(far)
library(glmnet)
library(xtable)
library(bigstep)




ks <- c(20, 100, 200)
beta_signal <- 3.5
sigma_matrix <- matrix(rep(0, 950*950), nrow=950)
diag(sigma_matrix) <- 1
rp2 <- 100

final_table_b <- matrix(nrow=3, ncol=4)
final_table_Y <- matrix(nrow=3, ncol=4)

ptm <- proc.time()
for(i in 1:3){
  k <- ks[i]
  b <- c(rep(beta_signal, k), rep(0, 950-k))
  b_l2 <- matrix(nrow=rp2, ncol=4)
  Y_l2 <- matrix(nrow=rp2, ncol=4)
  for(j in 1:rp2){
    X <- matrix(rnorm(1000 * 950, 0, sqrt(1/1000)), nrow=1000)
    eps <- rnorm(1000)
    Y <- X %*% b + eps
    
    # lambda that minimizes SURE == lambda that minimizes MSE
    lambda_sure = 950 / sum(b^2)
    
    # ridge with lambda that minimize SURE
    obj_sure <- glmnet(X, Y, lambda=lambda_sure, alpha=0, intercept=FALSE, standardize=FALSE)
    b_hat_sure <- coef(obj_sure)[2:951]
    
    b_l2[j, 1] <- sum((b_hat_sure - b)^2)
    Y_l2[j, 1] <- sum((X%*%b_hat_sure - X%*%b)^2)
    
    # ridge for lambda chosen by 10 fold CV
    obj_cv <- cv.glmnet(X, Y, alpha=0, intercept=FALSE, standardize=FALSE)
    b_hat_cv <- coef(obj_cv, s='lambda.min')[2:951]
    
    
    b_l2[j, 2] <- sum((b_hat_cv - b)^2)
    Y_l2[j, 2] <- sum((X%*%b_hat_cv - X%*%b)^2)
    
    # OLS
    obj_ols <- lm(Y~X-1)
    b_hat_ols <- summary(obj_ols)$coefficients[,1]
    
    
    b_l2[j, 3] <- sum((b_hat_ols - b)^2)
    Y_l2[j, 3] <- sum((X%*%b_hat_ols - X%*%b)^2)
    
    # OLS with model selected by mBIC2
    d <- prepare_data(Y, X, verbose = FALSE)
    obj_mbic2 <- fast_forward(d, crit=mbic2)
    s_obj_mbic2 <- summary(obj_mbic2)
    b_hat_mbic2 <- rep(0, 950)
    if(length(obj_mbic2$model) > 0){
      b_hat_mbic2[as.numeric(obj_mbic2$model)] <- s_obj_mbic2$coefficients[-1, 1]
    }

    b_l2[j, 4] <- sum((b_hat_mbic2 - b)^2)
    Y_l2[j, 4] <- sum((X%*%b_hat_mbic2 - X%*%b)^2)
    
    cat('Current iteration:', j, 'for k =', k, '\n')
    cat('Elapsed time:', (proc.time()-ptm)[3], '\n')
  }
  final_table_b[i, ] <- apply(b_l2, 2, mean)
  final_table_Y[i, ] <- apply(Y_l2, 2, mean)
}

final_table_b <- as.data.frame(
  final_table_b,
  row.names=c('k=20', 'k=100', 'k=200'),
  col.names=c('ridge_sure', 'ridge_cv', 'ols', 'mbic2')
)

final_table_Y <- as.data.frame(
  final_table_Y,
  row.names=c('k=20', 'k=100', 'k=200'),
)

colnames(final_table_b) <- c('ridge_sure', 'ridge_cv', 'ols', 'mbic2')
colnames(final_table_Y) <- c('ridge_sure', 'ridge_cv', 'ols', 'mbic2')

write.csv(final_table_b, file='r2result_task2_b.csv')
write.csv(final_table_Y, file='r2result_task2_Y.csv')













